# Feature: 3D Geometry Embedding and Clash Detection <!-- rq-0f766d5b -->

`CoordinateGenerator` (`src/biochar/pipeline/geometry_3d.py`) turns a bonded molecular graph into 3D
coordinates, and `GeometryValidator` decides whether the result is physically acceptable.

The recurring failure mode in this area is not a wrong coordinate — it is a **validation
threshold that rejects correct structures**. Every threshold below was retuned at least once
because it flagged real chemistry as an error, and the symptom each time was the same: every
point in a parameter sweep degrading to the fallback path while the geometry itself was fine.
The requirements here therefore pin *what must not be called a clash* as tightly as they pin
what must.

## Embedding Path Selection <!-- rq-80b04b19 -->

Two embedding strategies exist, and the choice between them is made by heavy-atom count, not
by user request.

- Molecules of **80 heavy atoms or fewer** embed with RDKit ETKDGv3 (falling back to ETKDGv2)
  followed by MMFF94 relaxation.
- Molecules of **more than 80 heavy atoms** use a 2D-first path: compute the flat aromatic
  layout — the hex-lattice positions the skeleton builder tagged, or a Kamada–Kawai layout when
  any *ring* carbon is untagged — and promote it to 3D at `z = 0`. Untagged atoms outside the
  rings are placed from their parent and do not take a sheet off the lattice (see *A Decorated
  Sheet Stays on the Lattice* below). The Kamada–Kawai layout is then
  force-field minimised. The hex lattice is not minimised: it is already exact, and the only
  relief it gets is local (see *Hex-Lattice Sheets Keep Their Lattice* below).

The reason for the split is that ETKDGv3 folds or collapses large flat fused-ring systems.
A large PAH sheet is genuinely planar, and an embedding algorithm that treats it as a flexible
molecule produces a crumpled sheet with irreducible self-clashes.

Force-field iteration count also scales with size: 300 iterations at or below 100 heavy atoms,
400 above 100, and 500 above 300.

Independently of size, a skeleton carrying a **ring-substituting nitrogen** (pyridinic,
pyrrolic, or graphitic) prefers the flat hex-lattice path, because ETKDGv3 buckles the sheet
around the substituted ring atom.

```gherkin
Feature: Select an embedding path that suits the molecule

  @rq-bb1fbb12
  Scenario: A small molecule uses the distance-geometry embedder
    Given a molecule with 80 or fewer heavy atoms
    When coordinates are generated
    Then ETKDG embedding followed by MMFF94 relaxation produces the coordinates

  @rq-265ec385
  Scenario: A large flat sheet uses the 2D-first path
    Given a molecule with more than 80 heavy atoms
    When coordinates are generated
    Then the flat 2D layout is promoted to 3D rather than embedded by ETKDG
    And the resulting sheet is not folded

  @rq-c4e683c4
  Scenario: Embedding retries a bounded number of times
    Given an embedding attempt that fails to produce coordinates
    When the generator retries
    Then it makes at most 3 independent embedding attempts before giving up
```

## Clash Detection <!-- rq-c2463a78 -->

A clash is atom overlap. The generic floor for a contact between atoms `i` and `j` is
**0.75 × (vdW radius i + vdW radius j)** — 2.55 Å for a carbon–carbon pair. A contact is only
reported when it is *meaningfully* below that floor; see the tolerance section.

### Hydrogen bonds are not clashes <!-- rq-00cacb46 -->

The generic 0.75 × vdW floor evaluates to **2.04 Å** for an O/H pair, which sits squarely
inside the physical H···A range of roughly 1.6–2.2 Å. Every intramolecular hydrogen bond
between adjacent hydroxyl groups therefore trips the generic floor. On high-oxygen chars
(O/C ≥ ~0.2, e.g. 400 °C softwood) such pairs are unavoidable, so strict validation failed on
every seed and every sweep point fell back.

A polar hydrogen pointing at an N or O acceptor is held to a reduced floor of
**`HBOND_MIN_H_ACCEPTOR_DISTANCE` = 1.5 Å** instead. The pair qualifies when the donor and
acceptor are both in `HBOND_DONOR_ACCEPTOR_ELEMENTS` (N, O) and the D–H···A angle is at least
**`HBOND_MIN_DHA_ANGLE_DEG` = 90°**.

The 90° gate is deliberately permissive. Real hydrogen bonds are near-linear (typically
> 120°); requiring only 90° asks that the hydrogen point *toward* the acceptor rather than the
acceptor being jammed into the side of the D–H bond. Hydrogen-bonded pairs are **not** excluded
outright — a contact shorter than 1.5 Å is too close even for a low-barrier hydrogen bond and
is still a clash.

```gherkin
Feature: Distinguish hydrogen bonds from steric clashes

  @rq-6b70c04c
  Scenario: An intramolecular hydrogen bond is not reported as a clash
    Given two adjacent hydroxyl groups whose H sits 2.0 Angstrom from the neighbouring O
    And the D-H...A angle is at least 90 degrees
    When geometry validation runs
    Then the contact is not reported as a steric clash

  @rq-6c0f608d
  Scenario: An overlap too short even for a hydrogen bond is still a clash
    Given a donor-acceptor pair whose H...A distance is below 1.5 Angstrom
    When geometry validation runs
    Then the contact is reported as a steric clash

  @rq-a63f775d
  Scenario: A non-polar contact does not get the hydrogen-bond floor
    Given a hydrogen and a carbon at 2.0 Angstrom with no N or O involved
    When geometry validation runs
    Then the contact is judged against the generic vdW floor, not the hydrogen-bond floor

  @rq-51d0fb61
  Scenario: A high-oxygen structure validates without falling back
    Given a structure generated at an O/C ratio of 0.2 or above
    When strict validation runs
    Then hydrogen bonds between adjacent hydroxyls do not force the fallback path
```

### A clash needs to be deeper than embedding noise <!-- rq-56e62d3f -->

Neither floor is meaningful to 0.01 Å. ETKDG embedding and MMFF relaxation leave contacts
scattered by more than that, so a hard `distance < floor` comparison is a knife edge: a contact
a hundredth of an ångström inside the floor gets reported as a clash.

A contact is a clash only when it is more than **`CLASH_SEVERITY_TOLERANCE` = 0.05 Å** below
its applicable floor. That value is well below genuine overlap — the real clashes this system
exists to catch (two functional groups embedded on top of each other) sit 0.2–0.5 Å inside the
floor. Anything within 0.05 Å is embedding noise, and GROMACS energy minimisation would move
the atom in its first step.

```gherkin
Feature: Tolerate knife-edge contacts that are embedding noise

  @rq-21834b0c
  Scenario: A contact marginally inside the floor is not a clash
    Given a contact 0.01 Angstrom below its applicable clash floor
    When geometry validation runs
    Then it is not reported as a clash

  @rq-b6b192ce
  Scenario: A genuine overlap is still reported
    Given a contact 0.3 Angstrom below its applicable clash floor
    When geometry validation runs
    Then it is reported as a clash with its severity

  @rq-92fec87b
  Scenario: Severity is reported so a marginal contact can be told from a real one
    Given any reported clash
    When the validation message is produced
    Then it states how far below the floor the contact sits
```

## Hex-Lattice Sheets Keep Their Lattice, Except Where a Flat Lattice Is Wrong <!-- rq-8e5f5dc1 -->

When `CoordinateGenerator.used_hex_lattice` is true, the generic clash resolver does not run. It
displaces atoms pairwise with no notion of a ring, so on a fused sheet it would try to move hundreds
of atoms and shatter the lattice — replacing a correct structure with a broken one.

A flat lattice is correct almost everywhere. It is wrong at **helicene-type edges**, and this used
to be misdescribed as "peri contacts that are real geometry":

- at a **cove** — two edge carbons 2.46 Å apart and four bonds apart, the [4]helicene motif —
  both substituents point into the same vacant lattice site and land 0.1–0.6 Å apart, and the cove
  carbons themselves sit at 2.46 Å, inside the 2.50 Å C···C floor;
- at a **fjord** — five bonds apart, [5]helicene — the edge carbons are 1.42 Å apart unbonded, and
  a substituent can land on the other arm's ring carbon.

Those are not physical features: real helicenes twist out of plane, and GROMACS energy minimisation
of an exported sheet opens a cove from 2.46 Å to ~2.95 Å and buckles its carbons up to ~1.0 Å out
of plane. Left unrelieved, they failed strict validation on nearly every structure above 80 heavy
atoms.

So the hex-lattice path runs its own relief, `CoordinateGenerator.relieve_lattice_crowding`, after
embedding and planarity enforcement:

1. substituents that clash are turned rigidly out of plane about their ring atom, the two sides of a
   crowded site in opposite directions;
2. a UFF relaxation then lets each cove open, with every ring atom tethered within
   `LATTICE_TETHER_RADIUS` of its lattice position and every aromatic bond held to
   `LATTICE_BOND_RANGE`.

The relaxed result is kept only if it has no more bond-length errors and no more clashes than the
tilted one, and its ring atoms stay within `LATTICE_MAX_SHEET_RMS` of a plane. Otherwise the tilted
coordinates are kept. A structure with no clash is returned unchanged.

Stage 1 is not optional. The overlaps it removes are close enough to coincident that a force field
started from them throws atoms apart: ring bonds stretched to 2–28 Å.

The relaxation also needs care in two places. UFF cannot parametrise an aromatic sheet it cannot
kekulise, and the de-aromatised fallback types ring carbons sp³ and folds the sheet, so the force
field sees a copy with aromatic atoms typed sp². And a perfectly flat start is a saddle point, so
ring atoms get a small seeded out-of-plane offset.

Two groups that crowd each other must tilt opposite ways. Directions are handed out greedily,
deepest contact first, so a pair can arrive already sharing one: at a fjord each arm's hydroxyl
oxygen lands on the *other* arm's carbon, those two contacts are the deepest, and each group gets a
direction on its own before their O···O contact is reached. Such a pair also tries turning its
second group the other way. Without this the outcome depended on where every hydroxyl hydrogen
happened to point, since that decides which contacts come first: rotating those hydrogens
correctly (rq-7515b128) took one sheet from 0 to 115 clashes and another from 74 to 0 (both
at 150 C, H/C 0.35, seed 4, measured before #74 stopped the skeleton builder overshooting).

The relaxation is allowed to converge. `LATTICE_RELAX_MAX_ITERATIONS` is a ceiling the minimiser
usually stops well short of; at 400 it cut descents off with contacts still open (100 C,
H/C 0.6, seed 1 kept 19 clashes).

A skeleton that is itself malformed — the elongated builder can overshoot its carbon target and
leave unbonded ring carbons at bond distance — is not repaired here. Relief makes it no worse, but
may leave it failing.

```gherkin
Feature: Keep the hex lattice, relieving only what a flat lattice gets wrong

  @rq-49e2ed82
  Scenario: Hex-lattice structures keep their lattice
    Given a structure whose coordinates came from the hex-lattice path
    When the generator finishes geometry
    Then the generic clash resolver does not run
    And every aromatic carbon-carbon bond stays at the lattice length

  @rq-3c400653
  Scenario: Small-molecule structures still get clash resolution
    Given a structure embedded by ETKDG rather than the hex lattice
    When the generator finishes geometry
    Then clash resolution runs

  @rq-57d3fc46
  Scenario: A cove's substituents are pulled apart
    Given a hex-lattice sheet with a cove whose two hydrogens were placed on top of each other
    When lattice crowding is relieved
    Then no contact in the structure is a steric clash

  @rq-da5ef2fc
  Scenario: Relief opens the cove rather than only moving its hydrogens
    Given a hex-lattice sheet with a cove
    When lattice crowding is relieved
    Then the two cove carbons are further apart than the flat lattice's 2.46 Angstrom

  @rq-13532295
  Scenario: A relieved sheet stays flat
    Given a hex-lattice sheet with crowded edges
    When lattice crowding is relieved
    Then its ring atoms deviate from their best-fit plane by no more than LATTICE_MAX_SHEET_RMS

  @rq-ac64435a
  Scenario: A relaxation that makes things worse is discarded
    Given a relaxation that would leave more clashes or bond-length errors than it started with
    When lattice crowding is relieved
    Then the tilted, unrelaxed coordinates are returned instead

  @rq-37a3d53e
  Scenario: A sheet with nothing to relieve is left alone
    Given a hex-lattice structure with no steric clash
    When lattice crowding is relieved
    Then its coordinates are unchanged

  @rq-9a59d62c
  Scenario: Large hex-lattice structures pass strict validation
    Given a strict configuration of 100 carbons whose composition target is reachable
    When structures are generated across several seeds
    Then none fails validation on a steric clash

  @rq-7cab77bd
  Scenario: Two groups crowding each other tilt opposite ways
    Given two substituents at a fjord, each crowding the other arm's carbon and each other
    When lattice crowding is relieved
    Then they turn to opposite sides of the sheet and no longer overlap

  @rq-47354e4f
  Scenario: The relaxation is not cut off before it converges
    Given a decorated 100-carbon sheet whose relaxation needs more than 400 iterations
    When the generator finishes geometry
    Then no contact in the structure is a steric clash
```

### A Decorated Sheet Stays on the Lattice <!-- rq-d6926052 -->

Aliphatic decoration — pendant methyls, and the `-CH2-OH` an O/C target can spill onto them —
adds carbons the skeleton builder never tagged, and so does a carboxyl group. The lattice path
once required *every* carbon to be tagged, so a single methyl sent the whole sheet to the
Kamada–Kawai layout and past this relief: at 150 C and H/C 0.35, seeds 2, 3 and 5 failed strict
validation on overlaps up to 0.40 Å inside the floor (issue #66).

The lattice is now kept whenever every **ring** carbon is tagged and every untagged heavy atom
hangs, through other untagged atoms, off a tagged one. Untagged atoms are placed breadth-first
from their parent: a ring atom's substituent points radially out, as an oxygen always has; an
off-lattice atom's children and hydrogens take its tetrahedral (four bonds), trigonal (three) or
bent (a two-bond oxygen) directions, with the in-plane direction on the side with more room.
Bond lengths are the bond-order-scaled covalent radii sums the validator expects, except a
phenolic C–O, which keeps 1.36 Å. Each such group is a substituent the relief can tilt.

The relief cannot clear everything a decorated sheet carries. Methyls on *both* carbons of a cove
land ~0.2 Å apart, and with `-CH2-OH` arms a rigid tilt leaves overlaps a force field then blows
apart. A free embedding — Kamada–Kawai plus the clash resolver and force field, every atom
movable — often clears those. So when the relief leaves a clash on a sheet with off-lattice
carbons, the generator also embeds it freely and keeps whichever result has fewer geometry
errors, the lattice one on a tie. Neither path is better everywhere: across 100–150 C, H/C
0.35–0.6, O/C 0–0.3 and five seeds, the lattice won at light decoration and the free embedding at
heavy. Keeping the better of the two took the clash-free count from 73 to 77 of 90 and made
14 points better. One point (150 C, 0.35/0.3, seed 5) went from 26 to 30 clashes: the free
embedding had fewer clashes but seven ring bonds stretched to ~2 Å, which is why the comparison
counts every geometry error, not just clashes. A sheet with no off-lattice carbon never embeds
twice.

```gherkin
Feature: Keep a decorated sheet on the lattice

  @rq-06dfeecf
  Scenario: A sheet with pendant methyls takes the hex-lattice path
    Given a hex-lattice sheet of more than 80 heavy atoms carrying untagged pendant methyls
    When coordinates are generated
    Then the hex-lattice path is used
    And every tagged atom sits at its lattice position

  @rq-91985117
  Scenario: An off-lattice carbon is tetrahedral and can be tilted
    Given a hex-lattice sheet carrying pendant methyls
    When coordinates are generated
    Then every methyl has three hydrogens at tetrahedral angles and no bond-length error
    And each methyl and its hydrogens form a substituent the relief can tilt

  @rq-d0935f4f
  Scenario: Lightly decorated sheets pass strict validation
    Given a strict configuration of 150 carbons at H/C 0.35 that reaches its target by methyl decoration
    When structures are generated for seeds 2, 3 and 5
    Then none fails validation on a steric clash

  @rq-ab5203c4
  Scenario: A decorated sheet the relief cannot clear keeps the better embedding
    Given a hex-lattice sheet with off-lattice carbons whose relief leaves a steric clash
    When the generator finishes geometry
    Then it also embeds the sheet freely
    And keeps whichever embedding has fewer geometry errors

  @rq-0cacf78d
  Scenario: An undecorated sheet is embedded once
    Given a hex-lattice sheet with no off-lattice carbon whose relief leaves a steric clash
    When the generator finishes geometry
    Then no free embedding is attempted
```

### A hydroxyl hydrogen keeps its bond angle <!-- rq-7515b128 -->

On the flat path each hydrogen was placed "opposite the parent's other neighbours". For a
hydroxyl or carboxyl OH that direction is the C–O axis itself, so the X–O–H angle came out at
180° on the roughly 30% of hydroxyls no clash drew attention to. That is wrong geometry, and it is
a singular torsion for any force field later handed the structure — the relief pass's UFF energy
started at 10¹⁷ until it was fixed.

A hydrogen on a heteroatom with a single heavy neighbour is placed at
`HYDROXYL_H_PLACEMENT_ANGLE_DEG` (108.5°) to that bond, in the sheet plane, on the side with more
room.

A crowded hydroxyl hydrogen is then moved by turning it about the C–O bond, which keeps the angle.
It used to be translated off that cone — a z offset of up to 0.9 Å, then the bond rescaled — which
bent C–O–H to 95–135°, and to 178° where the offset cancelled the hydrogen's own out-of-plane
component and left it on the C–O axis. Across 90 structures (100–150 C, H/C 0.35–0.6, O/C 0–0.3,
five seeds) 104 of 1,522 hydroxyls were outside 100–120°.

```gherkin
Feature: Place hydroxyl hydrogens at a real bond angle

  @rq-653a5fd2
  Scenario: A hydroxyl on a large sheet is bent, not linear
    Given a structure on the hex-lattice path carrying hydroxyl groups
    When coordinates are generated
    Then every C-O-H angle is between 100 and 120 degrees

  @rq-3e250d24
  Scenario: A crowded hydroxyl keeps its angle when it is moved
    Given a decorated hex-lattice sheet whose phenolic and -CH2-OH hydrogens are crowded
    When coordinates are generated
    Then every C-O-H angle is between 100 and 120 degrees
```

### The validator and the relief share one definition of a clash <!-- rq-5079b4df -->

`_clash_pairs` is the single implementation of "a judged pair more than
`CLASH_SEVERITY_TOLERANCE` inside its floor". `GeometryValidator._check_steric_clashes` reports what
it returns, and the relief pass resolves what it returns. A second implementation would let the
relief call a structure clean that strict validation then rejects.

```gherkin
Feature: One definition of a clash

  @rq-4eb90e63
  Scenario: The validator reports exactly the pairs the clash finder returns
    Given any structure with steric clashes
    When geometry validation runs
    Then the clashes it reports are the pairs the shared clash finder returns, in the same order
```

## Force-Field Refinement Is Not Gated on Clashes <!-- rq-e05e7026 -->

On the non-hex-lattice path every structure goes through force-field refinement, whether or not
clash detection reported anything. Clash *resolution* is conditional — there is nothing to
displace when no contact is flagged — but refinement is unconditional.

A clean clash report does not mean the embedding is strain-free. The two rules look at different
things: clash detection measures non-bonded contacts, and ETKDG can satisfy every one of those
while leaving a compressed aromatic bond of 1.16 Å where 1.40 Å is expected. Only the force-field
pass relaxes that.

Gating refinement on the clash report would make it depend on whichever contacts the clash rules
happen to flag. Such a dependency is invisible while it holds — refinement runs by side effect
for as long as some contact is always reported — and silently stops the moment the clash rules
become more precise. The two passes are therefore kept independent.

```gherkin
Feature: Refine geometry regardless of what clash detection found

  @rq-acf97ed2
  Scenario: A structure with no reported clashes is still refined
    Given a structure embedded off the hex-lattice path
    And geometry validation reports no steric clash
    When the generator finishes geometry
    Then force-field refinement runs

  @rq-5658c4b6
  Scenario: Clash resolution does not run when there is nothing to resolve
    Given a structure embedded off the hex-lattice path
    And geometry validation reports no steric clash
    When the generator finishes geometry
    Then the clash resolver is not invoked
```

## Bond-Length Validation <!-- rq-ff2c049b -->

`COVALENT_RADII` are single-bond radii, so their sum predicts a single bond only. The expected
length is scaled by bond order through `BOND_ORDER_LENGTH_FACTORS`: single 1.00, aromatic 0.92,
double 0.87, triple 0.79. An aromatic C–C is 1.40 Å, not the 1.52 Å the unscaled radii sum
implies; a C=O is 1.23 Å, not 1.42 Å.

Before this scaling every out-of-range message stated the wrong expectation, even when the bond
really was out of range — a diagnostic that misleads is worse than no diagnostic.

The accepted band is **`BOND_LENGTH_MIN_FACTOR` = 0.85** to **`BOND_LENGTH_MAX_FACTOR` = 1.40**
of the scaled expectation. These are tighter than the previous 0.8/1.5 band on purpose:
correcting the expectation downward for aromatic and multiple bonds would otherwise lower the
absolute floor and let a genuinely compressed bond through.

```gherkin
Feature: Judge bond length against the bond order actually present

  @rq-9a9488ab
  Scenario: An aromatic C-C is measured against the aromatic expectation
    Given an aromatic carbon-carbon bond of 1.40 Angstrom
    When bond-length validation runs
    Then the bond is accepted
    And the expected length reported is the aromatic value, not the single-bond radii sum

  @rq-64623ade
  Scenario: A carbonyl is measured against the double-bond expectation
    Given a C=O bond of 1.23 Angstrom
    When bond-length validation runs
    Then the bond is accepted

  @rq-50d3a0c4
  Scenario: A compressed bond is still rejected after order scaling
    Given a bond shorter than 0.85 times its order-scaled expected length
    When bond-length validation runs
    Then the bond is reported as out of range
```

## Every Geometry Error Is Reported <!-- rq-2ceb3ae6 -->

A structure report names every geometry error found, not a sample of them.

Completeness is what makes the report a diagnostic. Reporting a fixed-size prefix of the error
list biases it toward whichever check runs first: a sheet carrying one out-of-range bond and a
long tail of clash contacts shows the contacts and drops the bond, and the one error worth acting
on is the one that never appears. A caller that wants a summary can truncate a complete list;
a caller given a truncated list cannot recover what was dropped.

```gherkin
Feature: Report every geometry error, not a sample

  @rq-1f2ce9fc
  Scenario: A structure with more than three geometry errors reports all of them
    Given a structure whose geometry validation finds more than three errors
    When the structure report is produced
    Then the report contains every one of those errors
```

## Validation Does Not Rewrite What It Validates <!-- rq-0066e1ac -->

Structure validation is a read. The molecule that comes out of it is the molecule that went in.

The check that made this worth stating is `ChemicalFeasibilityValidator`'s connectivity test,
which asks RDKit to sanitise the molecule and records a warning if that raises. Sanitisation is
not a read: it kekulises, and kekulisation of a fused-ring sheet that has no kekulé structure
rewrites aromatic bonds to single on its way to raising. Run against the caller's own molecule,
the question answered itself and destroyed the subject in the same call.

What came back was worse than either state — atoms still flagged aromatic, every bond single —
and nothing downstream re-perceives the molecule, so that inconsistency is what was written to
the topology. The valence check in the same validator runs *before* the sanitisation, so the
report said the structure was fine while handing back one that was not.

Sanitising a copy answers the same question about the same graph.

```gherkin
Feature: Leave the molecule as it was found

  @rq-75c9724b
  Scenario: Validating a molecule that will not kekulise leaves its bonds alone
    Given a fused-ring molecule RDKit cannot kekulise
    When chemical feasibility validation runs on it
    Then its bond types and aromatic flags are unchanged

  @rq-8ebaed02
  Scenario: A molecule RDKit rejects is still reported
    Given a molecule whose sanitisation raises
    When chemical feasibility validation runs on it
    Then the report carries a sanitisation warning
```

## A Molecule Leaves Embedding Preparation Knowing Its Rings <!-- rq-07eab6b8 -->

The step that kekulises a molecule for embedding — or strips its aromaticity when the graph will
not kekulise — returns it with its ring information available, on both paths.

The de-aromatising path sanitises, and sanitisation is allowed to fail: a biochar graph carrying
one bad valence is still a graph the force field can be handed, so the failure is swallowed and
the geometry goes ahead. But `SanitizeMol` clears the molecule's computed properties, ring
information among them, before it does any of its work. A sanitisation that raises partway
therefore hands back a molecule that knows *less* about itself than the one that went in, and
every later question of the form "is this atom in a ring?" raises `RingInfo not initialized`
rather than answering.

Nothing near the failure notices. The molecule embeds, coordinates come out, and the generation
dies three stages downstream in atom typing — which asks that question of every atom — reporting
a missing-ring-info precondition rather than the valence that actually went wrong. Perceiving the
rings before returning costs nothing and does not sanitise, so it cannot disturb the bond types the
heteroatom stages set.

```gherkin
Feature: Do not hand on a molecule whose rings have been forgotten

  @rq-a157df17
  Scenario: A molecule whose sanitisation failed still knows its rings
    Given a molecule that cannot be kekulised and fails sanitisation
    When it is prepared for embedding
    Then its ring information is available

  @rq-ac433b55
  Scenario: A molecule that kekulises cleanly also knows its rings
    Given a molecule that kekulises
    When it is prepared for embedding
    Then its ring information is available
```

## Embedding Returns the Molecule It Was Given, Carrying the Coordinates It Found <!-- rq-2a552a1a -->

`generate_3d_coordinates` returns the caller's own molecule with a conformer on it. The
bond-order-rewritten copy that embedding and the force fields worked on does not escape.

The copy exists because RDKit's embedders and force fields need integer bond orders, and a large
fused-ring biochar sheet frequently will not kekulise. `_kekulize_or_dearomatize` then strips
**every** aromatic bond to single and clears every aromatic flag — a legitimate scaffold for
computing coordinates, and a molecule that no longer claims a single aromatic ring.

Returning that copy handed the rest of the pipeline a sheet whose chemistry was gone. Valence
validation reported `Valence 3 below minimum 4` for every ring carbon in the molecule at once,
atom typing fell back to its ring-and-degree proxy for aromatic carbon rather than reading the
flag, and aromatic-planarity enforcement found no aromatic ring to enforce anything on. Only the
coordinates were a result; the chemistry was scaffolding.

The refinement pass that runs afterwards, `validate_and_relax`, therefore builds its own
force-field-safe copy instead of parametrising what it is handed. It returns coordinates and
nothing else, so nothing about the caller's molecule needs to survive into it — and MMFF and UFF
would both refuse an aromatic sheet, quietly, leaving the coordinates unrefined and
refined-looking.

```gherkin
Feature: Do not return the scaffold in place of the molecule

  @rq-d5b18b3a
  Scenario: A structure that will not kekulise still comes back aromatic
    Given a fused-ring molecule RDKit cannot kekulise
    When coordinates are generated for it
    Then the returned molecule still has its aromatic bonds and flags

  @rq-79908020
  Scenario: The returned molecule carries the coordinates that were computed
    Given a molecule that embeds successfully
    When coordinates are generated for it
    Then the returned molecule's conformer holds those coordinates

  @rq-592b8c66
  Scenario: Refinement relaxes an aromatic molecule rather than declining it
    Given an aromatic molecule and a set of coordinates
    When force-field refinement runs
    Then the coordinates are changed by the refinement
```

## Cross-references <!-- rq-13a85be1 -->

- Clash warnings on hex-lattice structures interact with strict-mode seed retry in the sweep
  driver; see `parameter-sweep.md` when that document exists.
- The oxygen placement that produces the adjacent hydroxyls behind most hydrogen-bond contacts
  is specified in `heteroatom-assignment.md`.
- Background and the full reasoning for the clash-threshold changes live in
  `docs/solutions/bugs/physical-features-misread-as-geometry-errors.md`.
- Atom typing keeps its own ring-perception backstop for molecules arriving unperceived from
  anywhere else; see `rq-c6ab7cbe` in `opls-typing.md`. The failure that made it necessary is
  written up in `docs/solutions/bugs/failed-sanitisation-forgets-the-rings.md`.
