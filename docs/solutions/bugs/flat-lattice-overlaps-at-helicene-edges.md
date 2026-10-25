---
title: Overlaps on a flat lattice were dismissed as artefacts, and strict mode rejected nearly every large sheet
date: 2026-09-28
category: bugs
module: geometry_3d
problem_type: bug
component: geometry
severity: high
applies_when:
  - "Strict validation fails on most structures above ~80 heavy atoms"
  - "A clash on the hex-lattice path is about to be written off as an artefact"
  - "Handing an RDKit molecule that already has a conformer to a force field"
  - "Placing a substituent by a radial or 'opposite the neighbours' rule"
  - "A sheet carrying methyls, -CH2-OH or carboxyl groups is not taking the hex-lattice path"
  - "Choosing between two embeddings by a count of validation errors"
related_components:
  - validation
  - biochar_generator
  - carbon_skeleton
tags:
  - steric-clash
  - hex-lattice
  - helicene
  - cove
  - uff
  - conformer
  - masked-bug
  - aliphatic-decoration
  - kamada-kawai
---

# Overlaps on a flat lattice were dismissed as artefacts, and strict mode rejected nearly every large sheet

Issue #62.

## Context

On the hex-lattice path (>80 heavy atoms), ring carbons sit on an exact flat lattice and every
substituent is placed radially outward from its parent. The generic clash resolver is skipped
there, and the codebase's guidance said to ignore clash warnings on these structures: *"peri
contacts are real geometry; GROMACS energy minimisation resolves them."*

## Symptom

With the default `strict=True`, 100-carbon structures failed on **every** seed, whatever the
composition, even with no oxygen at all (`H_C_ratio=0.35, O_C_ratio=0.0`: 6/6 failed).
Temperature mode at 500 °C passed 1 in 6 at 100 C and above.

## What was actually there

Cataloguing every reported contact across 69 structures, and grouping each by how far apart its
two ring anchors are in space and through the graph, gave one family. It is **helicene-type
edges**, which a flat lattice cannot represent:

| Motif | Anchors | What the flat placement produced |
|---|---|---|
| bay (phenanthrene) | 2.84 Å, 3 bonds | H···O ~1.6 Å |
| **cove** ([4]helicene) | 2.46 Å, 4 bonds | substituents 0.1–0.6 Å apart; ring C···C 2.46 Å, inside the 2.50 Å floor |
| fjord ([5]helicene) | 1.42 Å, 5 bonds | substituent 0.1–0.3 Å from the other arm's ring carbon |

The "peri contacts at ~2.44 Å" in the old comment were real, but they weren't what was failing.
The contacts left behind were 6–7× shorter.

GROMACS energy minimisation settled what was physical. It opened each cove from 2.46 Å to
~2.95 Å, took cove H···H from 0.57 Å to 1.95 Å, and buckled ring carbons up to 1.0 Å out of
plane, which is exactly how a real [4]helicene twists. So the ring carbons had to move. Moving
only H and O, as first proposed in #62, could never have passed strict mode.

> This is the inverse of `physical-features-misread-as-geometry-errors.md`. There, a threshold
> called real chemistry a clash. Here, the prose called real overlaps chemistry. Both are settled
> the same way: take the structure to something that knows the physics (GROMACS EM) and see what
> it does.

## Fix

`CoordinateGenerator.relieve_lattice_crowding`, run in `_generate_geometry` on the hex-lattice
branch after embedding and planarity enforcement. It does three things:

1. **Tilt.** Each crowded substituent group turns rigidly out of plane about its ring atom, with
   the two sides of a site going opposite ways.
2. **Tethered relaxation.** A UFF pass on an sp2-typed copy. Ring atoms have flat-bottomed tethers
   (0.3 Å), and aromatic bonds are held to 1.40–1.42 Å.
3. **Guard.** The relaxed result is kept only if it has no more bond errors and no more clashes
   than the tilted one, and it's still flat (ring RMS ≤ 0.75 Å). Otherwise the tilted coordinates
   are kept.

Separately, a polar H is now placed at 108.5° to its bond, not collinear with it.

Results:

- **Temperature grid (500–800 °C, 60–200 C, 6 seeds):** clash failures went from 5–6 of 6 in every
  cell at 100 C and above to **none**.
- **Plain 100 C (H/C 0.35, O/C 0–0.1):** 18/18 pass, against 0/18 before.
- **The same example through GROMACS energy minimisation:** 1,068 steps to 444 kJ/mol, against
  1,337 steps to 620 kJ/mol before.

## Three traps, each of which produced a convincing wrong answer

**1. The stale conformer.** The molecule reaching the relief already carries its embedded
conformer. `AddConformer` appends a second one, and `UFFGetMoleculeForceField` uses the first. The
force field therefore minimised the *raw overlapping* coordinates every time and ignored the input
entirely.

The tell was an energy (1.3 × 10¹⁷) that did not change when atoms were moved by 3 Å. An energy
that doesn't respond to the coordinates you changed is minimising different coordinates.
`RemoveAllConformers()` first.

**2. The collinear hydroxyl.** "Place H opposite the parent's other neighbours" is the C–O axis
itself when the parent has one heavy neighbour. 30% of hydroxyls sat at a 180° C–O–H, the ones no
clash drew attention to, and `_optimize_h_positions` cannot rotate a collinear H. It made UFF's
torsions singular.

Separately, RDKit marks a phenolic O conjugated (sp2), which UFF reads as a carbonyl `O_2`, so the
copy retypes it sp3.

**3. The NaN that validated.** `validate_geometry` returns early on NaN coordinates, reporting only
"NaN values found". A metric that counted *clash* errors therefore scored a NaN structure as
clean. An early version of the tilt produced NaN (a degenerate axis once a group had turned 90°)
and looked like a 35/45 success.

Count every error kind, or check finiteness first.

## Also learned

- **UFF can't take a non-kekulisable aromatic sheet**, and `_kekulize_or_dearomatize` types ring
  carbons sp3, which folds the sheet. Single bonds plus explicit sp2 hybridisation and
  `SetNoImplicit(True)` work. The single-bond rest length (~1.47 Å) then needs a *stiff* distance
  constraint: at k = 500 the bonds settled at 1.46 Å. A relaxed bond sits at the constraint's upper
  edge, so the range should top out at the target (1.40–1.42 Å), not straddle it.
- **A perfectly flat start is a saddle point**, and a minimiser never leaves it. Seed a small
  out-of-plane offset.
- **Force fields don't start from near-coincident atoms.** Without the tilt stage, UFF stretched
  ring bonds to 2–28 Å.
- **Tuning on cached inputs can mislead.** The prototype cleaned several oversized skeletons from a
  cache built before the hydroxyl fix. The same skeletons from the real pipeline, with correct H
  placement, did not relieve. When a result depends on input details, re-capture the input from
  the real pipeline before tuning.

## Follow-up: decorated sheets fell off the lattice (issue #66)

The relief above only ran on sheets whose *every* carbon carried a lattice tag. Aliphatic
decoration adds pendant methyl carbons that the skeleton builder never tagged, so one methyl sent
the whole sheet to the Kamada–Kawai layout and past the relief. Nothing logged the switch: the
symptom was only strict failures at 150 C, H/C 0.35 (seeds 2, 3 and 5, overlaps up to 0.40 Å). A
coverage predicate that is stricter than the code needs routes work silently. When a path
"should" apply and the symptom says it didn't, check the gate before the stage.

Fixed by requiring only *ring* carbons to be tagged (`_lattice_covers_skeleton`,
`src/biochar/pipeline/geometry_3d.py:211`) and placing the rest outward with tetrahedral or
trigonal bond slots (`_offlattice_layout`, `geometry_3d.py:277`). Three things were learned on the
way:

- **Bond chemistry must come from the original molecule.** The embedding's working copy is
  de-aromatised, so a phenolic C–O read from it looks aliphatic and was placed at 1.42 Å instead
  of 1.36 Å. The layout takes `chemistry=mol` (`geometry_3d.py:1057`). This is the same family as
  `the-working-copy-that-became-the-result.md`.
- **Neither embedding dominates.** On the lattice, light decoration clears. Heavy decoration
  (H/C 0.6) puts methyls on both carbons of a cove about 0.2 Å apart. A rigid tilt cannot clear
  that, and the tethered UFF then blows apart. The free Kamada–Kawai + force-field path clears it
  instead. The generator now runs the free embedding only when the relief leaves a clash on a
  decorated sheet, and keeps the better of the two (`biochar_generator.py:953`). Across a 90-point
  grid (100–150 C, H/C 0.35–0.6, O/C 0–0.3, five seeds), clash-free results went from 73 to 77.
- **Compare total errors, not clashes.** Choosing by clash count alone picked a free embedding
  with 26 clashes and seven ring bonds stretched to about 2 Å over an intact lattice sheet with 30.
  `validate_geometry` also stops at the first NaN with a single error, so the comparison rejects
  non-finite coordinates first.

## Follow-up: crowded hydroxyls bent their angle, and the relief hinged on it

`_optimize_h_positions` moved a crowded O–H hydrogen by adding a z offset (up to 0.9 Å) to a
dihedral candidate and rescaling only the bond, so the C–O–H angle drifted to 95–135°. Where the
offset cancelled the hydrogen's own out-of-plane component it landed on the C–O axis: 178°, a
singular torsion for any later force field. On the 90-point grid (100/150 C × H/C 0.35/0.5/0.6 ×
O/C 0/0.1/0.3 × seeds 1–5), 104 of 1,522 hydroxyls were outside 100–120°, at 9 points.

Keeping only the dihedral sweep (every candidate a rotation about C–O, so the angle holds) fixed
all 104. It was tried once before and reverted, because 150 C, H/C 0.35, O/C 0.3, seed 4 went from
0 to 115 clashes. That was not caused by the hydrogens. It was exposed by them:

- **The tilt can hand two crowding groups the same direction.** Directions go out greedily,
  deepest contact first. At a fjord (two unbonded ring carbons 1.42 Å apart) each arm's OH oxygen
  sits 0.06 Å from the *other* arm's carbon. Those two contacts are the deepest, so each group is
  given +1 on its own, and when their O···O contact comes up both already tilt the same way: O···O
  0.72 Å, H···O 0.48 Å. UFF starts at 3.7 × 10⁸ kcal/mol, throws the oxygens 3–7 Å off their
  carbons, and the relax is rightly rejected. That leaves the tilted sheet with 100+ clashes.
  Which contacts come first depends on every hydroxyl hydrogen's position, so any change to them
  reshuffles which sheets fail. The old placement had the same problem at a different seed
  (O/C 0.1, seed 4: 74 clashes). Fix: a pair sharing a direction also tries flipping its second
  group (`_tilt_crowded_substituents`). The fixture `_fjord_hydroxyls` in
  `tests/test_lattice_relief.py` reproduces it on a [5]helicene.
- **The tethered UFF was cut off mid-descent.** With the tilt fixed, seed 4 still kept one fjord
  C···C contact 0.06 Å inside the floor at 400 iterations. At 1,000 it cleared from every z-noise
  seed tried. `LATTICE_RELAX_MAX_ITERATIONS` is now 2,000. The minimiser stops when it converges,
  so grid wall time went from 42 s to 48 s (one run each).

All three changes are needed. Ablations on the same grid, against 77 clash-free points and 554
total clashes before:

| Change | Clash-free | Total clashes | Bad angles | Points worse |
|---|---|---|---|---|
| Iterations only | 80 | 428 | 103 | 0 |
| Cone + iterations | 81 | 430 | 0 | 1 (O/C 0.3 seed 4: 0 → 1) |
| Cone + flip (400 its) | 82 | 299 | 0 | 1 (same point: 0 → 1) |
| **Cone + flip + iterations** | **85** | **206** | **0** | **0** |

**Re-measured after #74.** #74 stopped the elongated skeleton builder from overshooting, which
removed most of the fjords behind the seed-4 failures above. Those points then pass on main
without this change, so their regression tests were replaced. On the new main, the same 90-point
grid goes from 82 clash-free points and 342 clashes to 82 and 284, with bad angles 6 → 0. A second
grid of decorated, oxygenated sheets (100/150 C × H/C 0.5/0.6 × O/C 0.1/0.2/0.3 × seeds 6–15)
goes from 119 to 120 clash-free, with bad angles 110 → 0 and no point worse. Both remaining
changes still matter there. Without the flip, four 150 C, H/C 0.6, O/C 0 sheets keep 8–17 more
clashes. At 400 iterations, 100 C, H/C 0.6, O/C 0, seed 1 keeps 19 clashes; at 2,000 it keeps none
on the conda-forge RDKit (macOS) used for the grid. CI's pip-installed RDKit 2026.9.1 (Linux) keeps all 19 at 2,000, so a
clash-free assertion at that seed failed there. The regression test checks that the minimizer
converges on a fixture instead, since a seed's outcome is build-dependent and convergence is not.

One point is worse, and it is a known exception: 100 C, H/C 0.6, O/C 0.3, seed 3 goes from 0
clashes to 1. The lattice relief fails on this sheet with or without the change, so the generator
keeps the free embedding. That embedding then ends with one O···H contact 1.43 Å apart, 0.07 Å
inside the 1.50 Å hydrogen-bond floor. Turning the flip off restores 0 at this point but costs
the four sheets above. Nothing in the free path's force field is wrong. Hydroxyl start positions
just select a different local minimum. It fails strict validation by that one contact.

The lesson generalises the one above about the gate: when a local change to the input of a
multi-stage pipeline moves the outcome a lot, look for a greedy, order-dependent step downstream
before blaming the change. A greedy choice that is right for most inputs makes the result depend
on input ordering, so small changes look like they cause large regressions.

## What this does not fix

- **Oversized skeletons.** Some elongated builds overshoot their carbon target (asked for 200, got
  310) and leave unbonded ring carbons at bond distance (fjords). Relief makes them no worse but
  can't rescue them. That's a skeleton-builder defect.
- **The Kamada–Kawai path.** When H/C-reaching aliphatic carbons were untagged, the sheet took the
  Kamada–Kawai layout, not the hex lattice. This bullet first said its residual contacts were at
  most ~0.1 Å inside the floor. That understated them: at 150 C, H/C 0.35, seeds 2, 3 and 5 failed
  strict on 11–17 clashes up to **0.40 Å** inside it. Fixed in #66: decorated sheets now stay on
  the lattice, with a free embedding kept only when it leaves fewer errors (rq-d6926052).
- **Heavily decorated sheets.** After #74, seven grid points at H/C 0.6 still keep 16–56 clashes
  after both embeddings: 100 C, O/C 0, seed 5; 100 C, O/C 0.3, seed 2; and 150 C, O/C 0, seeds 1–5.
  The fjord contacts at 150 C, H/C 0.35, seed 4 that this bullet used to describe cleared once the
  relaxation was allowed to converge (see above).
- **Composition.** 650–800 °C targets below a skeleton's H/C floor are #63.
