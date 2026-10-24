# Feature: Carbon Skeleton Growth <!-- rq-f2d55972 -->

`PAHAssembler` (`src/biochar/pipeline/carbon_skeleton.py`) builds the polycyclic aromatic
skeleton every structure starts from. It is stage one of the pipeline: everything downstream —
heteroatoms, geometry, typing, export — decorates what this produces.

A skeleton is assembled from **whole rings**. That single fact shapes every requirement here.
Carbon counts are reachable only in the increments ring fusion allows, a "10% pentagon" request is
a probability applied per ring rather than a share of the finished structure, and aromaticity is
not an input at all — it is whatever the ring topology turns out to be.

So most of what a caller passes is an approximation, and the module's job is to be honest about
which. The failure mode is a request that is quietly not honoured: a parameter accepted and
ignored, a count met by a different number, or a fallback that returns something unrelated to
what was asked for.

## Feature API <!-- rq-f61f4d50 -->

- `PAHAssembler(seed: int = None)` <!-- rq-5840d1b3 -->
  - `generate(target_num_carbons: int, prefer_larger_pahs: bool = True, defect_fraction: float = 0.0, target_h_c: Optional[float] = None, heptagon_fraction: float = 0.0) -> CarbonSkeleton`
    builds a skeleton of approximately `target_num_carbons` carbons.
  - Every parameter it accepts changes the result. A knob that cannot affect the output is not
    part of the signature.
  - Targets within the library's range are satisfied exactly from a pre-validated structure;
    larger targets take the closest library seed and grow by ring fusion.

- `CarbonSkeleton` <!-- rq-419ffa85 -->
  - The assembled `mol` plus its ring composition. The realised carbon count is read from the
    molecule, not assumed from the request.

- `SkeletonValidator.validate(skeleton) -> Tuple[bool, List[str]]` <!-- rq-83521d1c -->
  - Structural checks on an assembled skeleton.

## A Target Within the Library Is Met Exactly <!-- rq-f7dda304 -->

`PAH_LIBRARY` holds pre-validated structures from benzene at 6 carbons up to 40. A target matching
one of those sizes is satisfied from the library rather than grown, so it is met exactly and the
result is a structure whose aromaticity and ring closure are already known good.

```gherkin
Feature: Use the validated library where it covers the request

  @rq-4af71330
  Scenario: A target matching a library size is met exactly
    Given a target carbon count that a library structure provides
    When a skeleton is generated
    Then the skeleton has exactly that many carbons
```

## A Larger Target Is Reached From Above, Never From Below <!-- rq-3f012547 -->

Above the library, growth proceeds by fusing whole rings, so only certain carbon counts are
reachable. When the requested count is not one of them, growth continues to the next reachable
count rather than stopping short: **the realised count is greater than or equal to the request,
never less.**

The direction is the requirement. A skeleton one ring short of the request has fewer edge sites
than the caller planned for, which propagates into the achievable H/C and O/C; overshooting by a
ring leaves the composition targets reachable. Deviations of up to about five carbons occur at
awkward targets.

The realised count is a property of the returned molecule, so a caller who needs the exact figure
measures it rather than assuming the request was met.

```gherkin
Feature: Reach an unachievable count from above

  @rq-59e23c2e
  Scenario: A target above the library is met or exceeded
    Given a target carbon count larger than any library structure
    When a skeleton is generated
    Then the skeleton has at least that many carbons

  @rq-5435c9d8
  Scenario: The overshoot stays within a ring or two
    Given a range of awkward target carbon counts
    When skeletons are generated for each
    Then none exceeds its target by more than eight carbons

  @rq-e8417061
  Scenario: A sheet that closes with an odd count does not keep growing
    Given large targets and seeds whose growth closes the sheet with an odd number of carbons
    When skeletons are generated for each
    Then none exceeds its target by more than eight carbons
```

The bound held only for the small targets it was first tested on. A grown sheet needs an even
carbon count to have a Kekulé structure, and the parity step used to fuse at the top-ranked edge
whatever that edge added. A hexagon on open edge adds four carbons and cannot change parity, so an
odd sheet grew hexagon by hexagon until a cove (three new carbons) happened to rank first: requests
for 150, 400 and 628 carbons returned 228, 936 and 1680, without complaint. The parity step now
fuses only where the addition is odd, or leaves the sheet odd for a retry.

The same runaway was behind issue #65. Elongated growth (`target_h_c` above the compact ceiling)
returned sheets up to 55% over the request. The extra hexagons wrapped around until the sheet met
itself, leaving pairs of ring carbons at bond distance on the lattice with no bond between them,
which no later relaxation can pull apart. Both are now held as invariants, for the elongated build
as much as the compact one.

```gherkin
Feature: Elongated growth honors the size and the lattice

  @rq-5e74710c
  Scenario: An elongated skeleton stays within a ring or two of the request
    Given targets from 100 to 250 carbons and target H/C values above the compact ceiling
    When skeletons are generated across several seeds
    Then none exceeds its target by more than eight carbons

  @rq-20d8d5ea
  Scenario: No two lattice carbons sit at bond distance without a bond
    Given targets from 100 to 250 carbons, with and without a target H/C
    When skeletons are generated across several seeds
    Then every pair of carbons closer than 1.6 Å on the lattice is bonded
```

## Defect Fractions Are Per-Ring Probabilities <!-- rq-9c8fe79c -->

`defect_fraction` and `heptagon_fraction` are the probability that **any one ring added during
growth** is a pentagon or a heptagon. They are not the fraction of rings in the finished
structure.

The realised share lands below the requested probability, because the seed a skeleton grows from
is all hexagons and only the added rings are subject to the draw. Asking for 0.3 yields roughly
0.22 of the finished rings.

Above about ninety carbons with heptagons requested, pentagons are placed two to a heptagon (see
"A Curved Sheet Is Chemically Possible and Cheap to Perceive"), so their realized share follows the
heptagons rather than `defect_fraction` alone, and a pentagon draw with no heptagon to sit on is
lost.

Naming them as probabilities rather than fractions matters because the parameter is how the Wood
et al. curvature model is reached, and a caller reproducing published ring statistics is comparing
against realised fractions.

```gherkin
Feature: Apply defect probabilities per ring addition

  @rq-c6b675bb
  Scenario: A zero defect fraction produces a pure hexagonal skeleton
    Given a defect fraction and a heptagon fraction of zero
    When a skeleton is generated
    Then every ring is a hexagon

  @rq-0c91f5d5
  Scenario: A positive defect fraction produces pentagons
    Given a positive defect fraction
    When a large skeleton is generated
    Then the skeleton contains at least one pentagon
    And the realised pentagon share does not exceed the requested probability
```

## Aromaticity Is an Outcome, Not a Request <!-- rq-cf3d7a56 -->

Aromaticity is decided by the ring topology the assembler builds. There is no way to ask for a
skeleton of a given aromaticity and no mechanism that could honour such a request, so the
assembler does not accept one.

A parameter that is accepted and documented as unused is worse than no parameter: it reads at the
call site as a request being made, and a caller who sets it has no way to discover that nothing
changed. The composition-level aromaticity target is handled in `generation-config.md`, which
clamps it to what the builder can realise — that is where the concept belongs.

```gherkin
Feature: Do not accept a request that cannot be honoured

  @rq-b6eec273
  Scenario: The assembler takes no aromaticity argument
    Given the skeleton generation entry point
    When its accepted parameters are inspected
    Then none of them is an aromaticity target
```

## A Skeleton Is Never Returned at an Unrelated Size <!-- rq-fb1a5f14 -->

If growth cannot produce a skeleton, that is a failure of the request and is raised. Returning a
fixed small structure instead — the library's pyrene, at 16 carbons — answers a request for 200
carbons with something an eighth of the size, and every downstream stage then works on it as
though it were what was asked for.

A caller can recover from an exception. A caller cannot recover from a plausible-looking molecule
of the wrong size, because nothing later in the pipeline re-checks the request.

```gherkin
Feature: Fail rather than substitute an unrelated skeleton

  @rq-f953b9b0
  Scenario: A skeleton that cannot be grown raises
    Given a request that graph growth cannot satisfy
    When a skeleton is generated
    Then the failure is raised rather than answered with a fixed structure
    And the message says another seed may succeed
```

## Growth Depends on the Seed, and a Failed Growth Is Retried <!-- rq-4b588285 -->

Growth is not deterministic in the request alone. The seed breaks ties between equally ranked
fusion edges and draws ring sizes, so one seed can close a sheet with an odd count or unbalanced
sublattices — no Kekulé structure — where the next does not. The failure message once said the
build "does not depend on the seed, so retrying will not help"; seeds 1 and 3 then succeeded where
2 failed, on identical arguments.

So a failed growth is retried on sub-seeds derived from the caller's seed before anything is
raised. Attempt zero is the caller's seed itself, so a build that succeeds first time is unchanged.
Later attempts are derived from the pair (seed, attempt), not from `seed + attempt`: with the
additive form, seeds 1 and 2 shared four of their five retries, so a "different" seed mostly
replayed the same growth.

```gherkin
Feature: Retry a seed-dependent failure

  @rq-a3c395ed
  Scenario: A growth that fails on the caller's seed is retried
    Given a request whose growth fails on the caller's seed alone
    When a skeleton is generated
    Then a skeleton of at least the requested size is returned

  @rq-7504e809
  Scenario: Neighboring seeds do not share their retries
    Given two consecutive seeds
    When the sub-seeds each would retry with are listed
    Then the two lists have no sub-seed in common
```

## A Curved Sheet Is Chemically Possible and Cheap to Perceive <!-- rq-2d15a0c0 -->

Pentagons and heptagons do not tile the plane, so once two growth fronts meet around a defect their
predicted atoms land *near* existing atoms rather than on them. The lattice rules for flat growth —
merge within 0.05 Å, bond anything under 1.45 Å — bonded those near-misses, giving carbons of degree
four to ten. No sanitization accepts that, and every sheet at the Wood et al. ring fractions
(`defect_fraction` 0.154, `heptagon_fraction` 0.077) above about 200 carbons failed on every seed.
Curved growth therefore plans each fusion first and refuses one that would put an atom in an
impossible place.

A sheet that is chemically valid can still be unusable. RDKit perceives aromaticity ring by ring,
then searches subsets of up to six fused rings until every ring **bond** is aromatic, and that
search is combinatorial in the number of rings. A pentagon or heptagon is never aromatic alone; a
5-7 pair is (azulene, ten electrons), but only the pair's outer bonds are marked, so its shared bond
needs a second complementary ring. One bare bond in a sheet of seventy rings makes every
sanitization — the pipeline does many — run for minutes. Wood et al.'s own curved blocks avoid this
by clustering: their heptagons carry two or three fused pentagons.

So above roughly ninety carbons (thirty rings, where the exhaustive search stops being cheap), a
heptagon is placed only as a 5-7-5 cluster, and a pentagon only on a heptagon. A finished sheet with
a bare bond is a failed attempt. Pentagons without heptagons cannot be paired at all, so a
pentagon-only sheet above that size is usually refused — the error says so and names the fix.

```gherkin
Feature: Build curved sheets that later stages can process

  @rq-22536e48
  Scenario: Every carbon in a curved sheet has at most three bonds
    Given the Wood et al. ring fractions
    When skeletons of a few hundred carbons are generated
    Then no carbon has more than three carbon neighbors

  @rq-8f5204cd
  Scenario: Wood-sized curved sheets build on every seed
    Given the Wood et al. ring fractions
    And targets of 102, 256 and 608 carbons
    When skeletons are generated across several seeds
    Then each is returned with pentagons and heptagons

  @rq-f1e2eb6a
  Scenario: Aromaticity perception on a large curved sheet terminates quickly
    Given a curved sheet of several hundred carbons
    When RDKit perceives its aromaticity
    Then every ring bond is aromatic and perception takes well under a second

  @rq-5af71ee2
  Scenario: A large pentagon-only request is refused with the reason
    Given a pentagon-only request of several hundred carbons
    When a skeleton is generated
    Then the failure names the missing heptagons as the cause
```

### Known limit: curved sheets carry more edge than flat ones <!-- rq-37961357 -->

A flat projection cannot close around a 5-7-5 cluster, whose net disclination is a cone. Fusions
refused for collision leave ragged fjords, and an occasional hole where a front could not close.
At 608 carbons a curved sheet has 0.25–0.34 of its carbons on an edge, where a flat flake has
0.14, and every edge carbon carries hydrogen or a group. A curved request at a low H/C (Wood's BC600c
at 0.234, BC800c at 0.105) is therefore built, but fails strict validation on H/C. Closing that gap
needs growth on a curved surface rather than a plane, and is not done here.

## A Library Entry That Cannot Be Loaded Is Reported <!-- rq-e7260a44 -->

The library is parsed once at construction. An entry whose SMILES will not parse, or will not
sanitise, is dropped from the working set.

Dropping it silently shrinks the library without changing anything a caller can see: targets that
were met exactly begin to be grown instead, and the structure that was pre-validated is quietly
replaced by one that was not. `PAH_LIBRARY` states that all its entries are validated; a dropped
entry means that statement is no longer true, and the discrepancy is worth surfacing rather than
absorbing.

```gherkin
Feature: Surface a library entry that does not load

  @rq-4dc5ce35
  Scenario: An unparseable library entry is reported
    Given a library containing an entry whose SMILES cannot be parsed
    When the assembler is constructed
    Then the entry is reported rather than silently dropped
```

## Cross-references <!-- rq-a60c7e3b -->

- The composition targets that drive `target_num_carbons` and `target_h_c`, and the aromaticity
  floor applied above this layer, are specified in `generation-config.md`.
- The hex-lattice coordinates this stage attaches, and why large sheets skip clash resolution,
  are specified in `geometry-embedding.md`.
- Functional-group placement onto the edge sites this stage produces is specified in
  `heteroatom-assignment.md`.
