---
title: Verify OPLS type mappings against a real oplsaa.ff, at all three depths
date: 2026-07-16
category: conventions
module: constants
problem_type: convention
component: testing_framework
severity: high
applies_when:
  - "Adding or changing an entry in GROMACS_OPLS_TYPE_MAP"
  - "Adding a functional group that introduces a bond or angle never emitted before"
  - "Adding any force-field constant to biochar/constants.py"
  - "Reviewing a claim that force-field data has been verified"
related_components:
  - opls_typing
  - gromacs_export
  - constants
tags:
  - opls-aa
  - force-field
  - gromacs
  - ffbonded
  - silent-failure
  - verification
---

# Verify OPLS type mappings against a real oplsaa.ff, at all three depths

## Context

An exported `.top` carries **no force-field parameters of its own**. `TOPFileWriter.write`
emits `[ bonds ]` and `[ angles ]` with a `funct` code and no `c0`/`c1`:

```python
f.write(f"{i:5d} {j:5d} 1\n")                    # bond: indices + funct, no params
f.write(f"{n1+1:5d} {j+1:5d} {n2+1:5d} 1\n")     # angle: same
```

So 100% of the physics comes from the `#include "oplsaa.ff/forcefield.itp"`, resolved by
the `opls_XXX` names written into `[ atoms ]` via `GROMACS_OPLS_TYPE_MAP`. Those names are
the entire contract.

That contract has now been checked at three successive depths, and **each depth was
believed sufficient until the next defect proved it wasn't**:

| Depth | Check | Defeated by |
|-------|-------|-------------|
| 1 | The mapped type **exists** | Every wrong mapping named a real type |
| 2 | The mapped type has the right **element/mass** | `CA-S-CA` — right element, right mass, still broken |
| 3 | Every emitted **bond/angle combination resolves** | now enforced — see Examples |

**Depth 1** (fixed in c26bacd): `SS` (thioether sulfur) mapped to `opls_209`, which is a
*carbon* — bonded type `CT`, mass 12.011. Pyridinic and pyrrolic N mapped to carbon types
too. `grompp` accepts any name that exists, so nothing errored; the simulation ran the
wrong chemistry.

**Depth 2**: the fix added element/mass consistency checks
(`TestTypesAreTheRightElement` in `tests/test_opls_type_map.py`). Believed sufficient at
the time — a prior session recorded that *"element consistency is the check that catches
this class."* Depth 3 falsified that.

**Depth 3**: `SS -> opls_222` is genuinely sulfur (mass 32.06, bonded type `S`), so
depth 2 passes cleanly — yet a thioether bridges two aromatic carbons and emits a
`CA-S-CA` angle, and `ffbonded.itp` has **no such angletype**. Only `CA S CT` (thioanisole)
and `CA S CM` exist, both `104.200` deg / `518.816` kJ/mol/rad². `grompp` rejected the
topology with "No default Angle types". Reproduced end-to-end with
`functional_groups={"thioether": 2}`.

Depth 3 is now checked for every functional group and every N-doping mode, and the
`CA-S-CA` angle is supplied by `SUPPLEMENTARY_ANGLE_PARAMS`. Turning the check on found
**two further gaps the same day**, which is the strongest evidence that the depth is the
right one — see "What the check found immediately" below.

## Guidance

When a generator maps internal types onto an external force field and relies on that force
field to supply the parameters, verification must reach all three depths. They catch
**non-overlapping** failure classes; passing depth 2 says nothing about depth 3.

The depth-3 lookup needs an indirection that is easy to get wrong: `ffbonded.itp` is keyed
by **bonded type** (column 2 of `ffnonbonded.itp`), *not* by the `opls_XXX` name. Several
opls names collapse onto one bonded type — `opls_145` and `opls_521` are both `CA`. The
chain is:

```
internal type -> GROMACS_OPLS_TYPE_MAP -> opls_XXX
              -> ffnonbonded.itp column 2 -> bonded type
              -> ffbonded.itp [bondtypes] / [angletypes]
```

Corollary, and the rule is narrower than "never hand-copy a value": **a local table may
hold only what the force field does not define.** Three tables that ignored this
(`OPLS_LJ_PARAMS`, `OPLS_BOND_PARAMS`, `OPLS_ANGLE_PARAMS`) drifted into an AMBER/OPLS
chimera with no single provenance and were deleted. The distinction is what makes
`SUPPLEMENTARY_ANGLE_PARAMS` safe where they were not:

- A value that **also exists** in the force field is duplication. Two sources of truth
  diverge — that is the entire failure history above.
- A value that exists **nowhere else** cannot drift. There is nothing to diverge from.

So a supplement is legitimate, and its boundary is mechanically enforceable rather than a
matter of discipline: `TestSupplementDoesNotShadowForcefield` fails if any entry duplicates
a stock angletype. Every entry also carries a provenance comment naming why the force field
lacks it and where the number came from.

## Why This Matters

The two failure modes are asymmetric, and the dangerous one is silent:

- **Wrong-but-existing type name** — runs to completion, emits a `.top` that looks
  perfectly normal, and simulates the wrong molecule. Nothing flags it: not `grompp`, not
  the suite, not reading the output. It surfaces only as physically implausible results,
  if anyone notices. This shipped, undetected, for an unknown period.
- **Missing bonded combination** — `grompp` refuses outright. Loud, but only at simulation
  time, which the test suite never reaches. A contributor can add a functional group, watch
  every test pass, and hand a user a topology that cannot start.

Both trace to the same design: the `.top` asserts names and delegates all physics. The
names are unverified data pretending to be code.

## When to Apply

- Adding or changing any `GROMACS_OPLS_TYPE_MAP` entry.
- Adding a functional group or heteroatom that can introduce a bond/angle combination never
  emitted before — the new combination has never been checked by anything.
- Any time "verified against the force field" is claimed: ask *at what depth*. The phrase
  has now been true-but-insufficient twice.
- Generally: when a data-driven mapping feeds a downstream lookup table, the *combinations*
  are a separate verification surface from the *entries*.

## Examples

### The depth-3 check

It lives in `tests/test_opls_type_map.py::TestBondedResolution`, parametrised over every
functional group and every N-doping mode. It is deliberately **not** reproduced here: a
copy of it in this file would be a second source of truth that drifts out of sync with the
real one — the exact failure this document is about.

What matters conceptually is the lookup chain, because the indirection in the middle is the
part that is easy to get wrong:

```python
# internal type -> opls name -> BONDED type -> ffbonded combination
def to_bonded(internal):
    return bonded_type_of[GROMACS_OPLS_TYPE_MAP[internal]]   # col 2 of ffnonbonded.itp

triple = (min(i, k), j, max(i, k))          # angles fix the middle atom
ok = triple in angletypes or internal_triple in SUPPLEMENTARY_ANGLE_PARAMS
```

The second clause is the whole design: a combination is acceptable if the force field
resolves it **or** we supply it deliberately. Anything else is a gap.

The check has teeth — deleting the `CA-S-CA` entry from `SUPPLEMENTARY_ANGLE_PARAMS` makes
it fail again:

```
AssertionError: 'thioether' emits terms no forcefield entry and no
SUPPLEMENTARY_ANGLE_PARAMS covers: ['angle CA-SS-CA (CA-S-CA)']
```

### What the check found immediately

Turning it on surfaced two gaps nobody was tracking. Both were held by `xfail(strict=True)`
until their fixes landed — so the fix could never silently leave a stale exemption behind —
and both are now closed. They are worth studying because they are **different kinds of
problem that look identical from the test output**:

- **`carboxyl` emits `O_2-C-OH` — a typing bug, not a missing parameter.** `AtomTyper`
  never assigns the carboxylic-acid types at all: its `OH2` branch sits behind `else` after
  `num_bonds == 1` and `== 2`, so it needs an oxygen with 3+ bonds, which cannot occur.
  `O`/`OH2`/`HO2` are dead code and a `-COOH` types as ketone O + phenol OH. OPLS *does*
  define the real angle — `O_3 C OH`, 121.000/669.440, "RCOOH". **Supplementing here would
  be the wrong fix**: it would cement the wrong chemistry and hide the typing bug behind a
  plausible number. The fix is correct typing.
- **`num_pyridinic` emitted `NC-CA-OH` — a genuine force-field gap, now supplied.** A
  hydroxyl on a ring carbon adjacent to a pyridinic N, reachable at the *default* O/C ratio.
  OPLS simply omits it (3-hydroxypyridine is plausible chemistry). This one *did* belong in
  the supplement, and the value is the phenol angle `CA-CA-OH` (120.000 / 585.760)
  transcribed verbatim — same hydroxyl-on-aromatic geometry, differing only in the far ring
  atom. It is the worked example of the supplement's intended use.

The lesson in the pair: when the check fires, the first question is not "what value do I
add?" but **"is the force field missing this, or are we asking it the wrong question?"**
`carboxyl` was the wrong question (a typing bug); `num_pyridinic` was a real gap. Only the
second justifies a supplement.

### What this check still does not cover

Know the edges, or it will be trusted further than it earns:

- **CI runs it now, but only because the force field is installed.** CI installs
  `gromacs-data` (the force-field files, not the binaries) and points `GMXDATA` at it, and a
  verify step fails loudly if `oplsaa.ff` is ever not discoverable — so a green badge does
  mean depth 3 actually ran. This was not always true: for a long time CI installed no
  GROMACS and every `requires_oplsaa` test silently skipped, the exact failure these tests
  exist to catch. If the `gromacs-data` install is ever dropped, the guarantee goes with it.
- **It checks tables, not `grompp`.** It proves a combination has an entry; it does not run
  GROMACS. An end-to-end `grompp` smoke test over one molecule per group would be strictly
  stronger.
- **It does not check dihedrals**, which the exporter also emits with `funct` only. The same
  class of gap could exist there and would not be caught.
- **Resolvable is not correct.** Depth 3 says a parameter *exists*, not that it describes
  the chemistry. `SS -> opls_222` and `NGR -> opls_379` remain deliberate approximations, and
  the supplemented `CA-S-CA` and `NPY-CA-OH` are nearest-analog values, not QM-validated.

### Known limitation: the analog values are not QM-validated (and, for now, that is fine)

The supplemented angles (`CA-S-CA`, `NPY-CA-OH`) and the approximate type mappings
(`SS -> opls_222`, `NGR -> opls_379`) are nearest-analog choices with traceable provenance
but no quantum-chemistry backing. Deriving them properly means, per parameter: build the real
model compound (3-hydroxypyridine, diphenyl sulfide), optimize + run a frequency calculation
at an OPLS-consistent level of theory (HF/6-31G\* to match the original fit), extract the
force constant via the Seminario method, and reconcile the GROMACS ½-factor. That is a
parameterization sub-project, not a one-off number, and it multiplies with every new gap.

**Decision (2026-07-17): deferred, deliberately.** No current use of this package depends
tightly enough on the geometry or dynamics *at these specific sites* to justify the work. For
bulk structural, surface-area, or adsorption-capacity results, a slightly-off aryl-S-aryl or
pyridinol bending constant is deep in the noise, and "chemically reasonable nearest-analog" is
the right amount of rigor.

**Revisit when** a specific result actually hinges on one of these sites — e.g. binding or
dynamics studied *at* a pyridinic-N-hydroxyl pocket or a thioether bridge, or any move toward
publication-grade parameters for these environments. The natural first target is `NPY-CA-OH`:
3-hydroxypyridine is a tiny, clean model, and the calculation would tell you directly whether
the phenol borrow holds. Charges are a separate axis and would need their own treatment (the
ionized types were derived by analogy too); QM-validating a *bend* does not touch the charge.

### Two traps that look like the answer

Both are tempting, confident, and wrong. Recorded because the next reader will reach for
them:

- **"The HA/HC epsilon gap is a clean factor of 2 — it's a kcal/kJ conversion bug."**
  It isn't. Every epsilon in the deleted table decodes as a round number in kcal/mol, so
  the conversion was fine. `HA`'s `(0.2600, 0.0630)` is *exactly* AMBER parm94's HA;
  `NA` is exactly AMBER's N; `CA`/`CT`/`C` are genuine OPLS. The table was a chimera, not a
  mis-conversion — there was no single error to undo, which is precisely why it was deleted
  rather than corrected.
- **"A wrong force constant means the table's convention is wrong."**
  Also no. The bond table's convention was *sound* — OPLS literature `kr` = half the
  GROMACS `kb` — and `CA-HA` (367), `CT-CT` (268), `CA-OH` (450) matched the stock file
  exactly. Individual values were still wrong (`CA-S` 340 vs 250; `CA-CA` 770 vs 469).
  Internal self-consistency across most rows is not evidence any row is right.

## Related

- PR [#20](https://github.com/layfield-lab/biochar/pull/20) — depth-1 fix
  (wrong-element mappings); added `tests/test_opls_type_map.py`.
- PR [#21](https://github.com/layfield-lab/biochar/pull/21) — deleted the three
  hand-copied tables and the unreachable `[ atomtypes ]` fallback.
- PR [#23](https://github.com/layfield-lab/biochar/pull/23) — supplied the
  `CA-S-CA` angle and landed the depth-3 check described here.
- PR [#25](https://github.com/layfield-lab/biochar/pull/25) — reconciled the
  `feat/ph-protonation` lineage: fixed the carboxyl typing (retiring that gap) and settled
  `NGR -> opls_379`.
- PR [#26](https://github.com/layfield-lab/biochar/pull/26) — gave CI `gromacs-data`,
  so the depth-3 check runs on every PR instead of skipping.
- Issue #3 — S-doping (thiol/thioether); origin of the types behind the `CA-S-CA` gap.
- `biochar/constants.py` — `GROMACS_OPLS_TYPE_MAP`, the `SS -> opls_222` note, and
  `SUPPLEMENTARY_ANGLE_PARAMS` with the rule that bounds it.
- `tests/test_opls_type_map.py` — the three depths against an installed `oplsaa.ff`. Its
  companion `tests/test_constants_ff.py` checks the map against a committed ground-truth
  table so element consistency runs even when no force field is present (e.g. locally
  without `GMXDATA`); the two are complementary, not duplicates.
- `docs/api/constants.rst` — states that LJ/bonded parameters live in `oplsaa.ff`, not in
  this module.
