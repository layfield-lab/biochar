# Feature: Linked Multi-Unit Structures <!-- rq-5a19bf9b -->

`src/biochar/pipeline/multi_unit.py` builds "island" structures: several separate aromatic units,
each a PAH from `PAHAssembler`, joined in a tree by linkers, and embedded so every unit stays flat
while the links rotate. It is the architecture of Wood, Mašek & Erastova 2024's low-temperature
building blocks, where the fused-sheet path can only build one ring system at every temperature.

The capability is opt-in. `BiocharGenerator` reaches it through three hooks — skeleton, oxygen
target, geometry — and a request that does not opt in never enters any of them.

The recurring hazard is that a linked structure looks like a successful generation whether or not
it is the architecture the caller asked for: a merged ring system, a hydroxylated linker, or a
coplanar layout all export, type, and pass grompp. So the requirements pin the architecture itself
— how many ring systems, what joins them, what the file says about them — not only validity.

## Feature API <!-- rq-c23659e6 -->

- `GeneratorConfig(linked_units: bool = False, num_units: Optional[int] = None, linker_mix: Optional[Dict[str, float]] = None, ...)` <!-- rq-0a635e0c -->
  - `linked_units` switches the capability on. With it off, `num_units` and `linker_mix` raise
    rather than being ignored.
  - `num_units` is an integer of at least 1, and `target_num_carbons` must hold a benzene per unit.
  - `linker_mix` maps linker names — `biaryl`, `diaryl_ether`, `aryl_alkyl_ether`, `propyl` — to
    non-negative weights with a positive sum.

- `plan_units(config, rng) -> UnitPlan` <!-- rq-22bc4f3c -->
  - Resolves the unit count and the linker that joins each non-core unit.

- `build_linked_skeleton(config, plan, aliphatic_enabled) -> Tuple[CarbonSkeleton, LinkedSkeleton]` <!-- rq-4a66d28f -->
  - Builds and joins the units, tagging unit atoms `unit_id` and linker atoms `linker_id`.

- `embed_linked_units(mol, seed) -> Tuple[Chem.Mol, np.ndarray]` <!-- rq-0b5952bf -->
  - 3D coordinates for a fully decorated linked molecule.

- `LinkedSkeleton` <!-- rq-eb8ad54c -->
  - The census: carbons per unit (core first), linker counts by type, pendant methyls. Exposed as
    `BiocharGenerator.linked_units` and `BiocharResult.linked_units`.

## The Default Path Is Untouched <!-- rq-f8203042 -->

A request without `linked_units=True` builds exactly what it built before this capability existed.
That is checked against a fixture captured from the generator before the hooks were added, not
against a second run of the current code — two runs of changed code agree with each other.

```gherkin
Feature: Leave every existing request unchanged

  @rq-b5d5c316
  Scenario: A request that does not opt in builds what it built before
    Given a configuration without linked_units
    When the structure is generated
    Then its topology, atom types and formula equal the pre-change baseline

  @rq-52b8fb84
  Scenario: Linked-unit parameters without the opt-in are refused
    Given a configuration that sets num_units or linker_mix without linked_units
    When the configuration is created
    Then it raises rather than building a fused sheet
```

## Unit Count and Linkers Follow Temperature <!-- rq-7cee9ec6 -->

The bands are a design choice anchored on Wood's six final models, the only temperature-labelled
structures in their library: BC400 is 4 ring systems, BC600 is 3, BC800 is 1. Below 500 °C a
linked request builds 4 or 5 units, from 500 to 700 °C 3, and from 700 °C the ordinary fused sheet.
Without a temperature it builds 4. Two-unit structures, common in Wood's library at every
composition, are reachable through an explicit count.

The low-temperature linker mix is the inter-unit links of Wood blocks 0001–0008 that this module
builds (diaryl ether 15, biaryl 9, propyl 4; the two diaryl ketones are deferred). The
mid-temperature mix is a judgment, not a count, with aryl–alkyl ether standing in for the
β-O-4-type links of 0052/0053.

Splitting a sheet adds perimeter hydrogen, so k units have an H/C floor. A temperature-derived
count shrinks, with a warning, when its floor is above the target — the caller asked for a
temperature, not for that count. An explicit count is built as asked and warned about.

```gherkin
Feature: Choose units and linkers from temperature unless told otherwise

  @rq-e918b504
  Scenario: The unit count follows the temperature band
    Given linked requests at each temperature below
    When the units are planned
    Then each plans the number of units listed
      | temperature | units  |
      | 400         | 4 or 5 |
      | 600         | 3      |
      | 800         | 1      |

  @rq-d9b926a8
  Scenario: At 700 degrees C and above the fused sheet is built unchanged
    Given a linked request at 800 degrees C
    When the structure is generated
    Then it equals the structure generated without the opt-in

  @rq-1ca67e8f
  Scenario: An explicit unit count and linker mix override the band
    Given a linked request at 400 degrees C with num_units 3 and only biaryl linkers
    When the units are planned
    Then 3 units joined by 2 biaryl links are planned

  @rq-85b48c1e
  Scenario: A temperature-derived count that cannot reach the target H/C shrinks
    Given a linked request at 600 degrees C for 50 carbons of wood
    When the units are planned
    Then fewer than 3 units are planned
    And a warning says why
```

## The Structure Is the Architecture Requested <!-- rq-9b18de79 -->

A linked structure has exactly as many ring systems as units, joined by one linker per non-core
unit, and every linker is the motif its name says. Linker carbons carry only hydrogen and their
bridging bonds: the oxygen stage's spill onto sp3 carbons would otherwise turn Ar–O–CH2–Ar into a
hemiacetal and a propyl bridge into an oxygenated arm the census does not report. A bridging group
the oxygen stage places across two units would merge two ring systems, so such a placement is
redrawn.

```gherkin
Feature: Build the requested units and linkers

  @rq-2b213b2e
  Scenario: The structure has one ring system per unit
    Given a linked request for k units
    When the structure is generated
    Then it has exactly k ring systems and k - 1 linkers

  @rq-6f9371cd
  Scenario: Each linker type builds its motif
    Given a linked request for each linker type, that type only
    When the skeleton is built
    Then every inter-unit connection is the motif listed for its type
      | linker           | motif                                              |
      | biaryl           | a single bond between aromatic carbons             |
      | diaryl_ether     | an oxygen bonded to two aromatic carbons           |
      | aryl_alkyl_ether | aromatic C, O, CH2, aromatic C                     |
      | propyl           | three sp3 carbons in a chain between aromatic carbons |

  @rq-b7b7c2cd
  Scenario: Linker carbons are never oxygenated
    Given a high-O/C linked request whose linkers carry carbon
    When the structure is generated
    Then no linker carbon is bonded to an oxygen outside its bridge

  @rq-7705e753
  Scenario: An explicit bridging group does not merge units
    Given a linked request with explicit ether groups
    When the structure is generated
    Then it still has one ring system per unit
```

## Composition Is Judged as for a Single Sheet <!-- rq-9612965b -->

Total carbon stays near the target, and H/C and O/C face the same strict validation as a fused
sheet. Linker oxygen counts toward O/C: the oxygen stage sizes its count from the carbon count and
ignores oxygen already present, so the linker oxygens are taken off its target.

```gherkin
Feature: Hit the requested composition

  @rq-227a0ddd
  Scenario: A BC400-like request passes strict validation
    Given a linked request for 104 carbons at H/C 0.635 and O/C 0.19 at 400 degrees C, wood
    When the structure is generated in strict mode
    Then it passes validation with 4 or 5 ring systems

  @rq-4115991e
  Scenario: A BC600-like request passes strict validation
    Given a linked request for 272 carbons at 600 degrees C, wood
    When the structure is generated in strict mode
    Then it passes validation with 3 ring systems

  @rq-8ff54b34
  Scenario: Linker oxygen counts toward the O/C target
    Given a linked request whose every linker is a diaryl ether
    When the structure is generated
    Then the realised O/C, counting linker oxygen, is within tolerance of the target
```

## Units Stay Flat and Links Rotate <!-- rq-84076187 -->

The fused-sheet embedder lays a large structure out flat and skips the force field, which would
make every link coplanar. The linked embedder places each unit rigidly from its hex-lattice
positions, builds linker atoms from bond lengths, angles and torsions, searches the torsions
against the validator's clash floors, then relieves contacts inside units with the hex-lattice
tilt and a tethered force-field pass that is kept only if it adds no defect and every unit stays
flat.

```gherkin
Feature: Embed linked units without flattening the links

  @rq-3a47b3fb
  Scenario: A linked structure embeds without clashes or bond-length errors
    Given a generated linked structure
    When its geometry is validated
    Then there are no steric clashes and no unusual bond lengths

  @rq-6d4235f9
  Scenario: Every unit stays flat
    Given a generated linked structure
    When each unit's ring atoms are fitted to a plane
    Then every unit's out-of-plane RMS is within the lattice flatness limit

  @rq-dacb8d13
  Scenario: Biaryl links are twisted out of plane
    Given a generated linked structure with a biaryl link
    When the torsion across the link is measured
    Then the two rings are between 20 and 160 degrees apart

  @rq-a9db6c13
  Scenario: The embedding is reproducible
    Given two generations with the same seed
    When their coordinates are compared
    Then they are identical
```

## Linkers Type and Simulate as Existing OPLS-AA Chemistry <!-- rq-04f1b337 -->

```gherkin
Feature: Type linker atoms with existing force-field types

  @rq-a7ae7f39
  Scenario: Linker atoms take their OPLS-AA types
    Given a linked structure for each linker type, that type only
    When atoms are typed
    Then the linker atoms type as listed
      | linker           | types                                    |
      | biaryl           | CA for both attachment carbons           |
      | diaryl_ether     | OS for the bridging oxygen               |
      | aryl_alkyl_ether | OS for the oxygen, CT with HC for the CH2 |
      | propyl           | CT with HC for every bridge carbon       |

  @rq-e72878a2
  Scenario: grompp accepts a topology with each linker
    Given an exported linked structure for each linker type, that type only
    When grompp processes it
    Then every one is accepted without warnings
      | linker           |
      | biaryl           |
      | diaryl_ether     |
      | aryl_alkyl_ether |
      | propyl           |
```

## The Files Say What Was Built <!-- rq-df3cb704 -->

The unit sizes and linker census are on the generator and the result, and the `.top` and `.itp`
headers carry them, because a fused sheet and an island structure of the same formula are
indistinguishable from the composition alone.

```gherkin
Feature: Record the linked architecture

  @rq-73fdab80
  Scenario: The topology headers name the units and links
    Given an exported linked structure
    When the .top and .itp are read
    Then their headers state the unit sizes and the linker counts

  @rq-420e40ce
  Scenario: The result carries the census
    Given a linked structure generated through generate_biochar
    When the result is inspected
    Then result.linked_units reports the unit sizes and linker counts
```
