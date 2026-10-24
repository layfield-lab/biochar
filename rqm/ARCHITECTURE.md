# Architecture <!-- rq-7efb46e4 -->

System design overview for `biochar`. This is the entry point for the requirements set in
`rqm/`; individual feature documents specify behaviour and carry the traceability IDs that link
requirements to the tests defending them.

Domain vocabulary lives in `CONCEPTS.md` at the repository root and is not repeated here. The
distinction between an **internal atom type**, an **OPLS type name**, and a **bonded type** is
load-bearing throughout the export path — read it there before working on typing or export.

## What this project is <!-- rq-bb785e6d -->

A Python library that generates biochar polycyclic aromatic hydrocarbon structures with
realistic composition and exports them as GROMACS topologies parameterised against OPLS-AA.
Inputs are composition targets (carbon count, O/C and H/C ratios, aromaticity, functional group
distribution, pH) or a pyrolysis temperature and feedstock from which those targets are derived.
Outputs are `.gro`, `.itp`, and `.top` files ready for molecular dynamics.

## Public surface <!-- rq-12b97d18 -->

`pyproject.toml` is the source of truth for what is public — `[project.scripts]` declares the
console entry points and `src/biochar/__init__.py`'s `__all__` declares the Python API. Do not infer
the public surface from this document.

| Console script | Module | Purpose |
|---|---|---|
| `biochar-gen` | `cli/cli.py` | Generate a single structure or series |
| `biochar-sweep` | `cli/sweep_cli.py` | Run a declarative factorial parameter sweep |
| `biochar-md-setup` | `cli/md_setup_cli.py` | Generate GROMACS run directories from a sweep manifest |
| `biochar-condense` | `cli/condensation_cli.py` | Set up a Wood et al. 2024 condensation-annealing run |

The four CLI modules parse arguments only; logic lives in the module each one wraps.

## Layering <!-- rq-18307b81 -->

The package is organised into subpackages under `src/biochar/` that follow the import graph, plus
one module — `constants.py` — kept at the package root because it is genuinely shared: every
subpackage below imports it, so nesting it inside any one of them would force the others to reach
across subpackage boundaries into an internal of a sibling.

The src layout is deliberate: because there is no importable `biochar/` at the repository root, the
package must be installed to be imported. A test run therefore cannot pass by exercising the working
tree while the installed package is broken — the failure mode that lets a packaging bug reach users
with CI green.

- **`constants.py`** (package root) — no intra-package imports; imported by every subpackage below
- **`pipeline/`** — the single-molecule generation pipeline and its supporting chemistry:
  `biochar_generator` (orchestrator), `carbon_skeleton`, `heteroatom_assignment`, `geometry_3d`,
  `opls_typing`, `valence`, `validation`, `protonation`. Depends only on `constants` and modules
  within the same subpackage.
- **`charges/`** — partial-charge backends: `qm_charges` (no intra-package deps), `ml_charges`
  (depends on `pipeline.opls_typing`)
- **`export/`** — GROMACS file writers: `gromacs_export` (depends on `constants`,
  `pipeline.opls_typing`), `md_setup` (depends on `workflows.condensation` for
  `anneal_spec_for_htt`)
- **`workflows/`** — higher-level orchestration built on `pipeline/`: `sweep`, `condensation`,
  `surface_builder` (also depends on `export.gromacs_export` and `pipeline.heteroatom_assignment`)
- **`models/`** — `temperature_model`, the data-driven temperature × feedstock composition model;
  no intra-package deps
- **`cli/`** — the four console-entry-point modules; parse arguments only

`export.md_setup` → `workflows.condensation` runs against that order, and is deliberate. Both
modules render the Wood annealing schedule, and two implementations of one published protocol
drift apart silently because both keep producing runnable input — they already had, disagreeing on
peak temperature and on which axis a semi-isotropic barostat holds. One implementation is worth
more here than layer purity, and no cycle is possible: `condensation` imports nothing from within
the package. If it ever needs to, move `AnnealSpec` and `anneal_spec_for_htt` down to
`constants.py` rather than duplicating them.

The dependency order is otherwise `constants` → `pipeline` → `charges`/`export` → `workflows` →
`cli`, with
`models` sitting alongside `pipeline` (both depend only on `constants`, and `pipeline` reaches into
`models` for the temperature/feedstock defaults). `__init__.py`'s `__all__` and `pyproject.toml`'s
`[project.scripts]` are what actually declare the public surface (see above); the subpackage
boundaries are an internal organisation, not a second public contract — a name's dotted path inside
`biochar.*` is free to move as long as it stays reachable from the top-level package.

This was a flat module directory until the subpackage split, and the one real cycle it had is
already gone: `generate_surface` used to live in `biochar_generator.py` despite being a
`workflows`-level entry point, which made `surface_builder` import `biochar_generator` at module
scope while `biochar_generator` imported `SurfaceBuilder` inside a function to dodge the resulting
circular import. `generate_surface` now lives in `workflows/surface_builder.py` alongside the class
it wraps, and the deferred import is gone.

## Single-molecule pipeline <!-- rq-60a0144a -->

`BiocharGenerator.generate()` runs five sequential stages, each in its own module.

1. **Carbon skeleton** (`pipeline/carbon_skeleton.py`) — `PAHAssembler` builds the PAH graph.
   Targets of 40 carbons or fewer take an exact SMILES from `PAH_LIBRARY`; larger targets select
   the closest seed and fuse rings until the count is within tolerance, with pentagon frequency
   controlled by `defect_fraction`. A hex-lattice position tracker keeps coordinates consistent
   so the geometry stage receives a flat, strain-free sheet.

2. **Heteroatom assignment** (`pipeline/heteroatom_assignment.py`) — `OxygenAssigner` places
   functional groups; `HydrogenAssigner` fills the remaining valences. Specified in
   `heteroatom-assignment.md`.

3. **Geometry** (`pipeline/geometry_3d.py`) — `CoordinateGenerator` embeds in 3D;
   `GeometryValidator` judges the result. Specified in `geometry-embedding.md`.

4. **Typing** (`pipeline/opls_typing.py`) — `AtomTyper` assigns internal types;
   `GROMACS_OPLS_TYPE_MAP` translates them at export. Specified in `opls-typing.md`.

5. **Validation** (`pipeline/validation.py`) — `ValidationEngine` checks composition ratios and
   geometry.

With `linked_units=True` the skeleton step builds several separate aromatic units joined in a tree by
linkers instead of one fused sheet (`pipeline/multi_unit.py`), the oxygen target is taken net of
linker oxygen, and a dedicated embedder keeps each unit flat while the links rotate. Off by default;
specified in `linked-units.md`.

## Surface pipeline <!-- rq-e13b2348 -->

`generate_surface()` delegates to `SurfaceBuilder`. Each sheet goes through the single-molecule
pipeline, is flattened to the xy plane by an SVD best-fit rotation, then stacked along z with
spacing `pore_diameter + 3.4 Å` — the graphene interlayer distance. Identical sheets are
generated once and deep-copied. The combined `.gro` carries each sheet as a separate residue;
the `.top` references one `.itp` for identical sheets, or one per unique sheet type.

## Export <!-- rq-92c69795 -->

`export/gromacs_export.py` holds `GROFileWriter`, `ITPFileWriter`, and `TOPFileWriter`,
orchestrated by `GromacsExporter`. Coordinates convert from RDKit ångström to GROMACS nanometres.
Residue names are hard-limited to 5 characters by the `.gro` format, which constrains
`molecule_name`.

## Standing constraints <!-- rq-faf5d18d -->

These hold across the whole system and are easy to violate locally.

- `molecule_name` and residue name must be **5 characters or fewer** — a GROMACS format limit.
- `_fix_heteroatom_bond_types` must be called after **any** `SanitizeMol` pass touching a
  molecule with ether oxygens.
- A clash on a hex-lattice structure is a defect, not an artefact: the lattice is exact, but its
  flat placement overlaps atoms at coves and fjords, which the geometry stage relieves
  (rq-8e5f5dc1).
- Requested composition is a target, not a promise. Measure the realised structure.
- Only the `opls_XXX` name reaches GROMACS. A name identifying the wrong element silently
  simulates the wrong chemistry.

## Requirements set <!-- rq-4f998061 -->

| Document | Covers |
|---|---|
| `geometry-embedding.md` | Embedding path selection, clash detection, bond-length validation, refinement, error reporting |
| `heteroatom-assignment.md` | Functional group placement, fallbacks, oxygen spill, hydrogen saturation |
| `opls-typing.md` | Internal typing, GROMACS mapping, forcefield verification depths |
| `gromacs-export.md` | Coordinate units, atom naming, exclusions and 1–4 pairs, impropers, charge reporting |
| `md-setup.md` | Stage ordering, ion placement, warning suppression, annealing schedule, run provenance |
| `generation-config.md` | Composition resolution, the aromatic floor, refused requests, strictness |
| `temperature-model.md` | Feedstock curves, per-property support, evidence behind an estimate |
| `surface-stacking.md` | Sheet separation, geometry labelling, copied-sheet state, topology agreement |
| `carbon-skeleton.md` | Library targets, growth from above, defect probabilities, honest parameters |
| `parameter-sweep.md` | Grid expansion, point naming, seed retry, failure modes, the manifest |
| `ph-protonation.md` | Henderson-Hasselbalch direction, sampling, the census, recorded provenance |
| `valence-validation.md` | Aromatic ring members, anion ranges, doping sites, safe bond addition |
| `condensation-annealing.md` | Wood anchors, the schedule, packing verification, run provenance |
| `charge-backends.md` | Backend refusal, charge conservation, the 1.14 factor, model provenance |
| `cli-arguments.md` | Config precedence, exit status, cross-argument consistency, stream separation |
| `linked-units.md` | Opt-in island structures: unit count by temperature, linker motifs, flat units, composition, records |

Every module the 2026-07-30 review listed is now specified, and so are the four console entry
points. Those were previously described here as holding no logic of their own; specifying them
showed otherwise. Precedence between a loaded config and the command line, the exit status a
pipeline branches on, and consistency between two individually-valid flags are all decided at
that layer and nowhere else — and each had a defect.

What remains unspecified is `constants.py`, whose tables are checked by
`tests/test_constants_ff.py` against the forcefield rather than against a promise this repository
makes.

An absent document is not an absent requirement — it means the behaviour is currently defended
by tests that no specification points at. `tools/rqm/rq check` reports the converse: specified
requirements that no test defends.
