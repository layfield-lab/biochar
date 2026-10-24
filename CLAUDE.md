# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

**Read `AGENTS.md` first.** It holds the working rules that apply to any coding agent — how requirements
move with the code, which green runs are not evidence, and the force-field rules. This file adds the
Claude-Code-specific configuration and the architecture reference; it does not repeat `AGENTS.md`.

## Project knowledge

- `rqm/ARCHITECTURE.md` — system design overview and the index of the requirements set. **Read this
  first** when orienting to the codebase.
- `rqm/*.md` — feature requirements with Gherkin scenarios. Each heading and scenario carries a stable
  `rq-XXXXXXXX` ID; tests reference those IDs in comments, and `rqm/registry.json` maps them both ways.
  These are **prospective** — the behaviour we promise. When changing behaviour in a specified area,
  update the requirement in the same change, then `bash tools/rqm/rq index`.
- `docs/solutions/` — documented solutions to past problems (bugs, conventions, patterns), organized by
  category with YAML frontmatter (`module`, `tags`, `problem_type`). These are **retrospective** — what
  we learned. Relevant when implementing or debugging in a documented area.
- `docs/reviews/` — dated, commit-pinned assessments of the codebase. **Point-in-time**: what was
  true at one commit, not a plan and not scheduled. Re-verify before acting on a finding.
- `CONCEPTS.md` — shared domain vocabulary. Relevant when orienting to the codebase or discussing domain
  concepts.

`rqm/` and `docs/solutions/` are complements, not duplicates: one states the invariant, the other records
how it was learned. A fix usually earns an entry in both.

```bash
bash tools/rqm/rq check          # which requirements no test defends
bash tools/rqm/rq show rq-xxxxxxxx   # one requirement and everything referencing it
```

Not every module is specified yet — `rqm/ARCHITECTURE.md` lists what is covered and what is not. An
unspecified module is not unverified; it means its tests answer to no written requirement.

## Commands

```bash
# Install from conda-forge (recommended)
conda install -c conda-forge biochar

# or from PyPI
pip install biochar

# Development install (editable, with extras)
pip install -e ".[dev,viz]"

# Install the pre-commit gate (once per clone)
bash tools/hooks/install.sh

# Fast tier -- ~30s. Always pass -n auto; serially this is minutes.
python -m pytest -m "not slow" -n auto

# Full suite, including the slow regressions -- ~1m25s under -n auto.
python -m pytest -n auto

# Run a single test
python -m pytest tests/test_generator.py::TestCarbonSkeleton::test_pah_assembler_benzene -v

# PAH quality suite (not pytest — run directly)
python tests/test_pah_quality.py

# Build docs
cd docs && make html
```

The test discovery root is `tests/` (configured in `pyproject.toml`), which holds 40 files.

**The suite has two tiers.** 15 of 1064 tests account for most of the runtime — end-to-end regressions
that generate real structures across several seeds. They carry the `slow` marker so the tier you run
constantly stays fast. Measured 2026-08-04 at `27b8244` on 14 workers, three runs each within 0.5s: the
fast tier is **1049 tests in ~30s**, the full suite **1064 in ~1m25s**. When adding a test that takes
more than a few seconds, mark the individual case `slow`, not its whole class — a class-level marker
exiles its fast siblings from the pre-commit tier.

**Both figures assume `-n auto`, and the difference is the whole point.** Serial numbers for this suite
have been quoted before without saying they were serial — an earlier "~22 minutes for everything" was
one, and it is roughly fifteen times the parallel figure. It made the full suite look like something to
route around when it costs a minute and a half. Always pass `-n auto`, and when recording a timing here,
say which one it is.

**The full-suite figure depends on how xdist distributes, not only on how much work there is.**
`pyproject.toml` sets `addopts = "--dist worksteal"`; without it the same suite takes 3m30s. The default
`--dist load` hands each worker a chunk up front, so the 43-second regressions bunch onto a few workers
and serialise there while the rest idle. The arithmetic is the tell: the fast tier alone is ~30s and the
slow tier alone is ~49s, and any full-suite figure much above their sum is a scheduling artefact rather
than work. If the full suite ever reads slower than about a minute and a half again, check that the flag
is still in effect before assuming a test got heavier.

**A run that skips more than once verified less than it appears to.** With GROMACS and MOPAC on `PATH`
the expected result is exactly **1 skip** — `test_grompp_smoke.py`'s deliberate case, which requires
`gmx` to be *absent*. Anything higher means the forcefield or QM-charge tests opted out silently. The
usual cause is `export GMXDATA="$CONDA_PREFIX/share/gromacs"` in a non-interactive shell, where
`CONDA_PREFIX` is empty, so it expands to `/share/gromacs` and does nothing. Write the path literally:

```bash
export GMXDATA="/opt/miniconda3/envs/biochar/share/gromacs"
export PATH="/opt/miniconda3/envs/biochar/bin:$PATH"
```

**The commit gate is `tools/hooks/pre-commit`**, installed by `bash tools/hooks/install.sh`. It runs
secrets, ruff, the fast tier, and `rq index`. Bypass with `git commit --no-verify`; CI still runs
everything.

`.claude/settings.json` carries permissions only. It used to hold a `PreToolUse` hook that ran the
full suite before commits; it was removed because it did not do that. Its `matcher` was `"Bash"` —
every command — and the `"if": "Bash(git commit*)"` that was meant to narrow it is not part of the
hook schema, so it was ignored. The result was the whole suite, unfiltered by the `slow` marker,
in front of arbitrary shell commands, against a 600s timeout it could not meet. If you want a
gate, `tools/hooks/pre-commit` above is it.

## Architecture

The package lives entirely in `src/biochar/` (src layout — the package is importable only when
installed, so a test run cannot silently exercise the working tree instead of the installed package).
After moving branches or a fresh clone, run `pip install -e ".[dev,ml]"` or imports will fail.

**`pyproject.toml` is the source of truth for what is public.** `[project.scripts]` declares the four console
entry points; `src/biochar/__init__.py`'s `__all__` declares the public Python API. Do not infer the public surface
from this file — check those two.

| Console script | Module | Purpose |
|---|---|---|
| `biochar-gen` | `cli/cli.py` | Generate a single structure or series |
| `biochar-sweep` | `cli/sweep_cli.py` | Run a declarative factorial parameter sweep |
| `biochar-md-setup` | `cli/md_setup_cli.py` | Generate GROMACS run directories from a sweep manifest |
| `biochar-condense` | `cli/condensation_cli.py` | Set up a Wood et al. 2024 condensation-annealing run |

### Module map

The package is organised into subpackages (`pipeline/`, `charges/`, `export/`, `workflows/`,
`models/`, `cli/`) that follow the import graph, plus `constants.py` kept at the package root
since every subpackage imports it. Full layering detail lives in `rqm/ARCHITECTURE.md`; this is
a quick index.

`pipeline/biochar_generator.py` orchestrates the single-molecule pipeline (five steps, below).
Beyond that pipeline:

- **`workflows/sweep.py`** — declarative parameter-sweep driver; `run_sweep`, `expand_grid`,
  `build_point`, and the sweep manifest that `export/md_setup.py` consumes.
- **`export/md_setup.py`** — writes GROMACS run inputs (`.mdp` templates plus a driver script) per
  structure; `setup_md_from_manifest`, ion profiles, pre-solvation stages.
- **`workflows/condensation.py`** — Wood et al. 2024 condensation-annealing setup (parallel
  construction mode).
- **`models/temperature_model.py`** — data-driven temperature × feedstock property model;
  `TemperatureModel`, `properties`, `VALID_FEEDSTOCKS`.
- **`charges/qm_charges.py`** — LigParGen-style QM partial charges (1.14*CM1A) via an external
  MOPAC binary. Requires a MOPAC install; raises `QMChargeError` when unavailable.
- **`charges/ml_charges.py`** — ML-based partial charge refinement.
- **`pipeline/multi_unit.py`** — opt-in linked multi-unit ("island") structures
  (`GeneratorConfig.linked_units`): several aromatic units joined by biaryl, ether and propyl
  linkers, unit count by temperature, and a rigid-fragment embedder. Specified in
  `rqm/linked-units.md`.
- **`pipeline/valence.py`** — valence validation system (see `docs/guides/VALENCE_SYSTEM.md`).
- **`cli/cli.py`**, **`cli/sweep_cli.py`**, **`cli/md_setup_cli.py`**, **`cli/condensation_cli.py`**
  — argument parsing plus the three decisions only this layer can make: precedence between a
  loaded config and the command line, the exit status, and consistency between flags that are
  individually valid but jointly impossible. Everything else lives in the module each one wraps.
  Specified in `rqm/cli-arguments.md`.

### Generation pipeline (single molecule)

`BiocharGenerator.generate()` runs five sequential steps, each in its own module under `pipeline/`:

1. **`carbon_skeleton.py`** — `PAHAssembler` builds a PAH graph. For targets ≤40 C it picks an exact SMILES from `PAH_LIBRARY` (in `constants.py`). For larger targets it selects the closest library seed then iteratively fuses rings (hexagons or pentagons controlled by `defect_fraction`) until the carbon count is within tolerance. A custom hex-lattice position tracker keeps all coordinates consistent so `geometry_3d` receives a flat, strain-free sheet.

2. **`heteroatom_assignment.py`** — `OxygenAssigner` places functional groups (phenolic, carboxyl, ether, etc.) either from an explicit dict or derived from `O_C_ratio`. `HydrogenAssigner` then fills remaining free valences. A critical function `_fix_heteroatom_bond_types` (also called in `biochar_generator.py` after geometry and after validation) corrects RDKit's tendency to mark ether C–O bonds as `AROMATIC` during sanitisation.

3. **`geometry_3d.py`** — `CoordinateGenerator` embeds in 3D. Molecules ≤80 heavy atoms use RDKit ETKDGv3/v2 + MMFF94; larger molecules use the pre-computed hex-lattice 2D coords. When `generator.used_hex_lattice` is `True`, the generic clash resolver is **skipped** — it would shatter the ring lattice — and `relieve_lattice_crowding` runs instead: it tilts substituents crowded at coves/fjords out of plane, then lets the coves open under a tethered UFF pass that keeps bond lengths and overall flatness (rq-8e5f5dc1).

4. **`opls_typing.py`** — `AtomTyper` maps each atom to an internal OPLS-AA type (CA, HA, OH, OS, …). At export time `GROMACS_OPLS_TYPE_MAP` in `constants.py` translates these to GROMACS `opls_XXX` names.

5. **`validation.py`** — `ValidationEngine` checks composition ratios and geometry. A steric clash surviving the hex-lattice relief is a real defect, not an artefact.

### Surface pipeline

`generate_surface()` → `SurfaceBuilder` (`workflows/surface_builder.py`). Generates each sheet through the single-molecule pipeline, flattens to the xy plane (SVD best-fit rotation), then stacks along z with spacing `pore_diameter + 3.4 Å` (graphene interlayer distance). Identical sheets are generated once and deep-copied. The combined `.gro` has all sheets as separate residues; `.top` references one `.itp` (identical sheets) or one `.itp` per unique sheet type.

### GROMACS export

`export/gromacs_export.py` contains three writers: `GROFileWriter`, `ITPFileWriter`, `TOPFileWriter`, orchestrated by `GromacsExporter`. Coordinates are converted from RDKit Å → GROMACS nm (× 0.1). Residue names are hard-limited to 5 characters by GROMACS `.gro` format.

### Key constraints

- `molecule_name` / residue name: **≤5 characters** (GROMACS hard limit).
- Carbonyl, quinone, and lactone functional groups fall back silently to phenolic — pure aromatic PAH edges have no free valence for these.
- `_fix_heteroatom_bond_types` must be called after any RDKit `SanitizeMol` pass that touches a molecule containing ether oxygens.
- The hex-lattice path (large molecules) keeps an exact lattice except at helicene-type edges (coves, fjords), where a flat lattice overlaps atoms; `relieve_lattice_crowding` relieves those. Do not dismiss a clash on these structures as an artefact — see rq-8e5f5dc1 and `docs/solutions/bugs/flat-lattice-overlaps-at-helicene-edges.md`.
