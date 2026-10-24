---
title: Curved sheets failed on every seed, and the obvious fix would have hung instead
date: 2026-10-23
category: bugs
module: carbon_skeleton
problem_type: bug
component: skeleton
severity: high
applies_when:
  - "A skeleton build with defect_fraction or heptagon_fraction fails, or hangs in RDKit"
  - "A SkeletonError appears on one seed and not its neighbors"
  - "A grown sheet comes back far larger than requested"
  - "Anything calls SanitizeMol or SetAromaticity on a large sheet containing pentagons or heptagons"
related_components:
  - biochar_generator
  - workflows.sweep
tags:
  - pentagon
  - heptagon
  - curvature
  - rdkit
  - aromaticity
  - kekulization
  - parity
  - seed
---

# Curved sheets failed on every seed, and the obvious fix would have hung instead

Found reproducing Wood, Masek & Erastova 2024's six models with `generate_biochar`. Three
separate defects lived behind one error message.

## 1. "Retrying will not help" — and the next seed worked

`SkeletonError` said the compact hex build "does not depend on the seed". It does: the RNG breaks
ties between equally ranked fusion edges. At 104 C, seed 2 closed the sheet with an odd count; an
odd bipartite graph has no Kekulé structure. The parity step then fused at the top-ranked edge
regardless of what it added. A hexagon on open edge adds four atoms and cannot fix parity, so the
sheet grew hexagon by hexagon until a cove (three atoms) ranked first — to 463 C, still odd. Worse,
the same loop *succeeded* silently elsewhere: requests for 150, 400 and 628 C returned 228, 936 and
1680. The overshoot requirement (≤ 8 C) was only ever tested on targets under 100.

**Fix:** parity fuses only where the addition is odd; a failed build retries on sub-seeds derived
from `(seed, attempt)`. Not `seed + attempt` — that made seeds 1 and 2 share four of five retries.

## 2. Near-miss atoms bonded into degree-10 carbons

Odd rings do not tile the plane. When two fronts met around a defect, predicted atoms landed
0.1–1.4 Å from existing ones instead of on them, and the lattice rule "bond anything under 1.45 Å"
bonded them all. Every Wood-ratio sheet over ~200 C had carbons with up to ten neighbors.

**Fix:** curved growth plans each fusion against contact bands (merge < 0.35 Å, bond 1.25–1.60 Å
at degree ≤ 3, otherwise clear by 2.0 Å) and refuses fusions that break them.

## 3. The valid sheets then hung RDKit for minutes

With valid topology, `SetAromaticity` ran 0.5 s at 104 C, 2.5 s at 130 C, and never finished at
200 C. RDKit (`applyHuckelToFused`, 2025.09) stops only when every ring **bond** is aromatic,
otherwise enumerating every subset of up to six fused rings. A pentagon or heptagon is never
aromatic alone. A 5-7 pair is (10 π), but RDKit marks only the pair's *envelope*, so the shared
bond needs a second complementary ring. One bare bond in a 70-ring sheet means C(70, 6) subsets,
on every sanitization.

Two wrong theories cost time here and are worth not repeating: "every atom aromatic" (RDKit counts
bonds, not atoms) and "a 5-7 pair suffices" (the shared bond stays bare). Small extracted fragments
also misled, because cut bonds left sp3 atoms. What settled it was reading RDKit's source and then
measuring RDKit itself on accepted sheets.

**Fix:** above ~30 rings, heptagons are placed as 5-7-5 clusters (Wood's own motif — their
heptagons carry two or three pentagons), and a finished sheet with a bare bond is a failed attempt,
checked by `_aromaticity_search_is_bounded` *before* RDKit sees it. Converting first would hang,
not fail. Curved sheets are kekulized by maximum matching, since RDKit's backtracking kekulizer is
the other place this blows up.

## Still open

- Curved sheets carry 0.25–0.34 of their carbons on an edge at 608 C, against 0.14 flat, so low
  H/C targets fail strict validation (rq-37961357). Needs growth on a curved surface.
- RDKit's `Kekulize` gives up on any sheet of ~600 C, **flat ones included**, and on Wood's own
  0053/0056. The pipeline evidently never re-kekulizes these, but nothing guarantees that.
