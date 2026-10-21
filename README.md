# Biochar Simulator — Structure Generator

A Python package for generating realistic biochar molecular structures for GROMACS molecular dynamics simulations. Supports single molecules, temperature/composition series, and **porous slit-pore surfaces**.

---

## Overview

The Biochar Simulator builds polycyclic aromatic hydrocarbon (PAH) structures with user-specified compositional parameters and exports GROMACS-ready force field files. Supports both pure hexagonal aromatic skeletons and topologically disordered structures with pentagon ring defects.

**Key capabilities:**

- **Size**: 6 to 200+ carbons (exact count for most targets)
- **Composition**: H/C and O/C ratios with configurable tolerances
- **Aromaticity**: 100% aromatic skeleton (graphene-nanoflake topology, or defective with pentagons)
- **Ring Defects**: Optional pentagon insertion during graph growth (`defect_fraction` parameter)
- **Functional Groups**: Exact counts via dict API — phenolic, carboxyl, ether, carbonyl, quinone, lactone, hydroxyl
- **Porous Surfaces**: Slit-pore systems of stacked PAH sheets with user-controlled pore diameter
- **Output**: GROMACS `.gro`, `.top`, `.itp` files with OPLS-AA force field

> **Related project:** [**biochar-pfas**](https://github.com/layfield-lab/biochar-pfas)
> drives this package to set up PFAS-sorption simulations (PFOA/PFOS/PFBS),
> inserting the ligands through the generic
> [pre-solvation seam](#pre-solvation-molecule-insertion-extension-seam).

---

## Installation

### conda (recommended)

```bash
conda install -c conda-forge biochar
```

### PyPI

```bash
pip install biochar
```

### Requirements

- Python 3.9+
- RDKit ≥ 2023.9
- NumPy ≥ 1.24, SciPy ≥ 1.10, NetworkX ≥ 3.1

---

## Quick Start

### Single molecule

```python
from biochar.pipeline.biochar_generator import generate_biochar

mol, coords, gro_path, top_path, itp_path = generate_biochar(
    target_num_carbons=100,
    H_C_ratio=0.5,
    O_C_ratio=0.1,
    output_directory="output",
    basename="biochar_100C",
    seed=42,
)
```

### With specific functional groups

Use a dict to place **exact counts** of each group type:

```python
mol, coords, gro, top, itp = generate_biochar(
    target_num_carbons=50,
    functional_groups={"phenolic": 3, "carboxyl": 1, "ether": 2},
    output_directory="output",
    basename="biochar_fg",
    seed=42,
)
```

When `functional_groups` is `None` (default), total oxygen is controlled by `O_C_ratio` and placed as phenolic groups.

### With pentagon ring defects

Add topological disorder by inserting 5-membered rings during growth:

```python
# ~15% of rings will be pentagons instead of hexagons
mol, coords, gro, top, itp = generate_biochar(
    target_num_carbons=60,
    defect_fraction=0.15,  # probability per ring addition
    output_directory="output",
    basename="biochar_defects",
    seed=42,
)
```

`defect_fraction` ranges from 0.0 (pure hexagonal PAH) to 1.0 (all pentagons). Typical values: 0.1–0.2.

#### Understanding `defect_fraction`

`defect_fraction` is a *request probability*, not a guaranteed pentagon fraction. During skeleton growth each ring-addition event has a `defect_fraction` chance of *attempting* a pentagon; parity and adjacency constraints reject many of these attempts. Observed pentagon counts are lower than the requested fraction:

| `defect_fraction` | ~50 C skeleton | ~100 C skeleton |
|:-----------------:|:--------------:|:---------------:|
| 0.05              | 0–1 pentagons  | 1–3 pentagons   |
| 0.15              | 1–3 pentagons  | 2–6 pentagons   |
| 0.30              | 2–5 pentagons  | 4–10 pentagons  |

Use `result.ring_composition` to inspect the actual counts after generation:

```python
result = generate_biochar(target_num_carbons=60, defect_fraction=0.15, seed=42)
print(result.ring_composition)  # e.g. {"hexagons": 12, "pentagons": 2}
```

### Slit-pore surface

```python
from biochar.workflows.surface_builder import generate_surface

# Two identical sheets, 10 Å pore
sheets, gro, top, itps = generate_surface(
    target_num_carbons=50,
    functional_groups={"phenolic": 2, "ether": 1},
    pore_diameter=10.0,
    output_directory="output",
    basename="slit_pore",
    seed=42,
)

# Asymmetric pore — different chemistry on each wall
sheets, gro, top, itps = generate_surface(
    pore_diameter=8.0,
    sheet_overrides=[
        {"functional_groups": {"phenolic": 3}, "target_num_carbons": 40},
        {"functional_groups": {"carboxyl": 2}, "target_num_carbons": 50},
    ],
    output_directory="output",
    basename="asymmetric_pore",
)
```

### Batch generation (temperature/composition series)

```python
from biochar.pipeline.biochar_generator import generate_biochar_series

configs = [
    {"molecule_name": "BC400", "target_num_carbons": 80,  "H_C_ratio": 0.65, "O_C_ratio": 0.20, "seed": 1},
    {"molecule_name": "BC600", "target_num_carbons": 100, "H_C_ratio": 0.55, "O_C_ratio": 0.12, "seed": 2},
    {"molecule_name": "BC800", "target_num_carbons": 120, "H_C_ratio": 0.40, "O_C_ratio": 0.05, "seed": 3},
]

results = generate_biochar_series(
    configurations=configs,
    output_directory="output/temperature_series",
    create_combined_top=True,   # writes combined.top
)
```

Then run in GROMACS:

```bash
cd output/temperature_series
gmx grompp -f md.mdp -p combined.top -o topol.tpr
gmx mdrun -deffnm topol
```

---

## Parameter sweeps (declarative grids)

`generate_biochar_series` takes a hand-written list of configs. For a full
**factorial grid** — every combination of pyrolysis temperature × feedstock ×
oxygen content — use the sweep driver instead. You describe the axes once in a
YAML/JSON file and the pipeline expands the grid, generates every structure,
and writes a manifest table for downstream analysis.

### CLI

```bash
# Emit a starter config you can edit
biochar-sweep template > my_sweep.yaml

# Run a sweep
biochar-sweep run examples/sweeps/temperature_grid.yaml
```

### Config format

```yaml
name: temperature_grid
output_directory: sweep_out/temperature_grid
seed: 1000              # base RNG seed; each grid point gets a distinct offset
max_retries: 8          # strict-validation seed retries before fallback
on_validation_fail: fallback   # fallback | skip | strict
name_template: "T{temperature}_{feedstock}"

axes:                   # cartesian product of every list below
  temperature: [300, 400, 500, 600, 700]
  feedstock:   [softwood, hardwood]

fixed:                  # applied to every grid point
  target_num_carbons: 100
```

Temperature and feedstock are mapped to `H_C_ratio` / `O_C_ratio` targets via
`biochar.models.temperature_model`, so a single temperature axis drives the oxygen
chemistry automatically. The example
[`examples/sweeps/oxygen_group_grid.yaml`](examples/sweeps/oxygen_group_grid.yaml)
sweeps explicit `functional_groups` counts instead.

### Python API

```python
from biochar.workflows.sweep import load_sweep_config, run_sweep

cfg = load_sweep_config("examples/sweeps/temperature_grid.yaml")
summary = run_sweep(cfg, quiet=True)

print(summary["n_built"], "structures")   # -> 10
import pandas as pd
df = pd.read_csv(summary["manifest_csv"])
```

### Output

```
<output_directory>/
├── manifest.csv          # one row per grid point (see below)
├── manifest.json         # same data + run metadata
└── structures/
    ├── 000_T300_softwood/   # .gro / .top / .itp for this point
    ├── 001_T300_hardwood/
    └── ...
```

The **manifest** carries one row per structure with the axis values
(`axis_temperature`, `axis_feedstock`, …), the achieved composition
(`molecular_formula`, `H_C_ratio`, `O_C_ratio`, `functional_groups`), the
build `status`, `seed_used`, `n_attempts`, validation counts, and the paths to
the three GROMACS files — ready to join against downstream analysis.

### Validation and the `status` column

Each grid point is generated **strict-first**: the driver retries with
successive seeds (`max_retries`) seeking a structure that passes full
composition + geometry validation. If none passes, `on_validation_fail`
decides what happens:

| `status` | Meaning |
|---|---|
| `strict_pass` | Passed strict validation on one of the retry seeds |
| `fallback` | Strict never passed; built with `strict=False` (files still written) |
| `skipped` | Strict never passed and `on_validation_fail: skip` — no files written |
| `failed` | A non-validation error occurred (captured in the `error` column) |

`on_validation_fail: strict` writes no row at all: the sweep is refused, naming
the point and its validation errors. Choose it when the grid is only meaningful
if every point is in tolerance.

This table is the same set as `biochar.workflows.sweep.POINT_STATUSES`, and a
test fails if the two drift apart.

The seed-retry → fallback design exists so a sweep completes deterministically
rather than stalling when a particular grid point cannot pass strict validation
(e.g. minor flat-sheet steric clashes that resolve under GROMACS energy
minimisation) — the `.gro/.top/.itp` files are still well-formed.

---

## Configuration Parameters

### Single molecule (`GeneratorConfig` / `generate_biochar`)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `target_num_carbons` | int | 50 | Target carbon count |
| `H_C_ratio` | float | 0.5 | Target H/C molar ratio |
| `H_C_tolerance` | float | 0.10 | Allowed H/C error (fraction) |
| `O_C_ratio` | float | 0.1 | Target O/C molar ratio |
| `O_C_tolerance` | float | 0.10 | Allowed O/C error (fraction) |
| `aromaticity_percent` | float | 90.0 | Target % aromatic carbons |
| `functional_groups` | dict\|None | `None` | Exact group counts, e.g. `{"phenolic": 2}` |
| `defect_fraction` | float | 0.0 | Probability [0, 1) each ring is a pentagon |
| `allow_aliphatic` | bool | True | Allow pendant sp3 (aliphatic) carbon to reach H/C above the pure-aromatic ceiling (see below). Set False to force a purely aromatic structure. |
| `charge_method` | str | `"opls"` | Partial charge source: `"opls"`, `"ml"`, or `"qm"` |
| `molecule_name` | str | `"BC"` | Residue name (≤5 chars for GROMACS) |
| `periodic_box` | bool | False | Add periodic box vectors to `.gro` |
| `seed` | int\|None | None | RNG seed for reproducibility |

### Slit-pore surface (`SurfaceConfig` / `generate_surface`)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `target_num_carbons` | int | 50 | Carbons per sheet |
| `H_C_ratio` | float | 0.3 | Target H/C per sheet |
| `O_C_ratio` | float | 0.05 | Target O/C per sheet |
| `functional_groups` | dict\|None | `None` | Groups applied to all sheets |
| `defect_fraction` | float | 0.0 | Pentagon probability per ring per sheet |
| `pore_diameter` | float | 10.0 | Gap between sheet surfaces (Å) |
| `num_sheets` | int | 2 | Number of parallel sheets |
| `sheet_overrides` | list\|None | `None` | Per-sheet config dicts (length = `num_sheets`) |
| `box_padding_xy` | float | 1.0 | Box padding in x/y (nm) |
| `box_padding_z` | float | 1.0 | Box padding in z (nm) |
| `system_name` | str | `"SLIT"` | Name in `.top [ system ]` section |
| `sheet_base_name` | str | `"SHT"` | Residue name base (≤3 chars) |
| `seed` | int\|None | None | RNG seed |

### Available functional groups

| Group | O added | Notes |
|---|---|---|
| `"phenolic"` | 1 | Ar–OH; always works |
| `"hydroxyl"` | 1 | Same as phenolic for pure PAH |
| `"carboxyl"` | 2 | Ar–C(=O)(OH); adds extra C |
| `"ether"` | 1 | Ar–O–Ar bridge across two edge sites |
| `"carbonyl"` | 1 | Falls back to phenolic with warning |
| `"quinone"` | 2 | Falls back to phenolic with warning |
| `"lactone"` | 2 | Falls back to phenolic with warning |

Carbonyl, quinone, and lactone require ≥2 free valence on one carbon, which is unavailable on pure aromatic PAH edge sites — they warn and substitute phenolic automatically.

### Partial charge methods (`charge_method`)

| `charge_method` | Description | Requirement |
|---|---|---|
| `"opls"` (default) | Static OPLS-AA lookup table | None |
| `"ml"` | GP model trained on OPLS reference charges | `pip install biochar[ml]` (scikit-learn) |
| `"qm"` | 1.14×CM1A via MOPAC AM1 semiempirical | `conda install -c conda-forge mopac` |

The `"qm"` backend reproduces the LigParGen 1.14*CM1A methodology using an external MOPAC binary. It runs a single-point AM1 calculation on the generated 3D geometry and maps the resulting Mulliken charges through the CM1A correction. See `docs/qm-charge-backend.md` for full details. Raises `biochar.QMChargeError` if MOPAC is not on PATH or exits with an error.

---

## H/C ratio control

A hydrogen only attaches at a carbon with a free valence. In a fully aromatic
(sp2) flake those are exactly the **perimeter** carbons, so the achievable H/C
is bounded by the perimeter/area ratio — roughly 0.44 at 50 C and 0.36 at 100 C
for a maximally condensed sheet, falling as the sheet grows. Requesting a
higher H/C than that ceiling used to silently produce a hydrogen-deficient
structure. The generator now reaches the requested H/C in three stages:

1. **Honest reporting.** If a target is above what the built structure can
   carry, `result.composition.h_c_ceiling` and `.h_c_target_unreachable` are set
   and a warning is logged, instead of the shortfall surfacing only later in
   validation.
2. **H/C-aware skeleton growth.** When a higher H/C is requested, the skeleton
   is grown into a less-condensed (more elongated, higher-perimeter) shape,
   raising the pure-aromatic ceiling toward ~0.5.
3. **Aliphatic decoration.** Above the pure-aromatic ceiling, the generator
   builds a smaller aromatic core and attaches pendant sp3 **methyl** groups so
   the total carbon count still matches `target_num_carbons` while H/C reaches
   the target (0.6–0.8). Aromaticity drops accordingly — matching the aliphatic
   content of low-temperature biochar — and OPLS typing maps the added carbons
   to `opls_135` (CT) with `opls_140` (HC) hydrogens.

To force a purely aromatic structure (H/C then capped at the aromatic ceiling,
with a warning), set `allow_aliphatic=False` or request
`aromaticity_percent ≥ 99`.

---

## Generation Pipeline

### Single molecule

```
target_num_carbons
      │
      ▼
Carbon skeleton (PAH seed + ring-growth)
      │  Parity-aware hex-lattice builder; 100% aromatic for any size
      ▼
Oxygen assignment (functional groups dict or O/C-ratio-driven)
      │
      ▼
Hydrogen assignment (fill valences → target H/C ratio)
      │
      ▼
3D geometry
      │  ≤80 heavy atoms: ETKDGv3 / ETKDGv2 embedding + MMFF94
      │  >80 heavy atoms: 2D-first embedding (flat graphene sheet) + FF minimisation
      ▼
OPLS-AA atom typing & partial charges
      │
      ▼
Validation (composition, geometry, steric clashes)
      │
      ▼
GROMACS export (.gro / .top / .itp)
```

### Slit-pore surface

```
SurfaceConfig
      │
      ▼
Generate N sheets (each via single-molecule pipeline above)
      │  Identical sheets: generate once, deep-copy remainder
      │  Distinct sheets:  generate each independently
      ▼
Flatten each sheet to xy plane (SVD best-fit plane rotation)
      │
      ▼
Stack along z: sheet_i centroid at z = i × (pore_diameter + 3.4 Å)
      │
      ▼
Compute periodic box (bounding box + padding)
      │
      ▼
Centre system in box
      │
      ▼
GROMACS export
      │  .gro  — all N sheets as separate residues, single file
      │  .itp  — one file (identical sheets) or one per sheet (distinct)
      └─ .top  — includes forcefield + itp(s); [ molecules ] count = N
```

---

## Supported Sizes

| Range | Strategy | Aromaticity | Count accuracy |
|---|---|---|---|
| 6–40 C | Exact PAH library match | 100% | Exact |
| 41–200+ C | Library seed + 4-node ring growth | 100% | ≤5% error |

### PAH library (18 validated entries)

| Molecule | Carbons | Type |
|---|---|---|
| benzene | 6 | Classic |
| naphthalene | 10 | Classic |
| anthracene / phenanthrene | 14 | Linear / angular |
| pyrene | 16 | Pericondensed |
| chrysene / tetracene / triphenylene | 18 | Various |
| pentacene / picene / hex_lattice_22 | 22 | Various |
| coronene | 24 | 7-ring pericondensed |
| hexacene / dibenzo_bc_ef_coronene | 26 | Various |
| hex_lattice_28 | 28 | Compact nanoflake |
| hex_lattice_30 | 30 | Compact nanoflake |
| hex_lattice_38 | 38 | Compact nanoflake |
| hex_lattice_40 | 40 | Compact nanoflake |

---

## Typical H/C ratios by pyrolysis temperature

| Temperature | H/C | O/C | Notes |
|---|---|---|---|
| 300–400 °C | 0.6–0.8 | 0.15–0.25 | Partially carbonised, high O |
| 500–600 °C | 0.4–0.6 | 0.08–0.15 | Moderately graphitic |
| 700–800 °C | 0.2–0.4 | 0.02–0.08 | Highly graphitic, low O |

---

## Output Files

| File | Format | Contents |
|---|---|---|
| `.gro` | GROMACS structure | Atom positions in **nm**, box vectors |
| `.top` | GROMACS topology | Force field include, molecule definitions |
| `.itp` | Include topology | Atoms, bonds, angles, dihedrals for one molecule type |

Coordinates are in **nanometers** (GROMACS convention; RDKit Å × 0.1).

For surfaces, a single `.gro` contains all sheets as separate residues, and the `.top` references one `.itp` with a molecule count (identical sheets) or one `.itp` per unique sheet type.

---

## Source Modules

| Module | Responsibility |
|---|---|
| `carbon_skeleton.py` | PAH library lookup and ring-growth engine |
| `heteroatom_assignment.py` | Oxygen (functional groups) and hydrogen placement |
| `geometry_3d.py` | 3D coordinate generation and clash resolution |
| `opls_typing.py` | OPLS-AA atom types and partial charges |
| `gromacs_export.py` | `.gro` / `.top` / `.itp` writers (single and multi-sheet) |
| `surface_builder.py` | Slit-pore surface assembly (`SurfaceBuilder`, `SurfaceConfig`) |
| `validation.py` | Composition, chemistry, and geometry checks |
| `constants.py` | OPLS-AA parameters, PAH library, VdW radii |
| `biochar_generator.py` | Public API: `generate_biochar`, `generate_surface`, `generate_biochar_series` |

---

## GROMACS run setup

`biochar-sweep` builds structures; `biochar-md-setup` turns each row of a
sweep manifest into a ready-to-submit GROMACS run directory, following the
dry-anneal / solvate / wet-equilibrate protocol (Wood, Mašek & Erastova
annealing schedule). It writes files only — no `gmx` binary is invoked — so it
runs anywhere and the resulting directories are submitted on the user's own
workstation or HPC cluster.

### CLI

```bash
biochar-md-setup sweep_out/temperature_grid/manifest.csv \
    --output-root sweep_out/temperature_grid/md_runs \
    --ion-profile mn_calcareous_default
```

### Python API

```python
from biochar.export.md_setup import setup_md_from_manifest

results = setup_md_from_manifest(
    "sweep_out/temperature_grid/manifest.csv",
    output_root="sweep_out/temperature_grid/md_runs",
    ion_profile="mn_calcareous_default",
)
for r in results:
    print(r["label"], r["status"], r["run_dir"])
```

Only manifest rows with `status` in `{strict_pass, fallback}` are processed;
`skipped`/`failed` rows are reported back (with a `skipped_reason`) rather
than silently dropped.

### Per-structure pipeline

Each run directory gets its own `.mdp` files and a single executable
`run_pipeline.sh` that chains the stages with `gmx grompp`/`mdrun`:

1. **Dry EM** — steepest descent, `emtol = 500`
2. **Anneal NVT** — 300 K pre-equilibration with velocity generation
3. **Anneal NPT** — simulated annealing, Berendsen 100 bar, 300→1000→300 K over 5.5 ns
4. **Final dry NPT** — Berendsen 1 bar, 300 K, 2 ns
5. **Box + solvate** — `editconf -d <pad> -bt cubic` (box size computed from
   the actual final dry structure, not hand-tuned per run) + `solvate`
6. **Ion addition** — one `genion` call per non-zero cation in the chosen
   ion profile (correct valence per species: Ca²⁺/Mg²⁺ → `-pq 2`,
   Na⁺/K⁺ → `-pq 1`), each with `-conc` and `-neutral`
7. **Wet EM** — `emtol = 1000`
8. **Wet NVT** then **wet NPT** — semiisotropic pressure coupling, C-rescale

### Ion profiles (Minnesota water chemistry)

`IonProfile` is a small dataclass (`ca_mM`, `mg_mM`, `na_mM`, `k_mM`,
`counter_ion`) selected by name via `--ion-profile` / `ion_profile=`:

| Profile | Description |
|---|---|
| `pure_water` | No background electrolyte — solvation-only control |
| `mn_calcareous_default` | Illustrative Ca-HCO₃-type Minnesota glacial-till groundwater (Ca²⁺/Mg²⁺-dominated, low Na⁺/K⁺) |
| `na_dominated` | Sodium-dominated scenario, e.g. softened or road-salt-impacted water |

> **These are illustrative starting points, not measured values.** Minnesota
> groundwater ion concentrations vary widely by aquifer and county. For a
> specific site, override the profile with monitoring-well data — see the
> MN DNR (*Bulletin 26, Natural Quality of Minnesota's Ground Water*) and MPCA
> ambient groundwater monitoring reports — by constructing a custom
> `IonProfile(...)` and passing it directly instead of a preset name.

### Pre-solvation molecule insertion (extension seam)

`MDSetupConfig.pre_solvation_stage` is a domain-neutral hook for inserting extra
molecules into the equilibrated slab's box **after dry-anneal, before
solvation** — the pipeline anneals the bare surface, inserts the molecules
(`gmx insert-molecules`), then runs box+solvate+ions and the wet stages against
a topology that already accounts for them:

```python
from biochar.export.md_setup import MDSetupConfig, PreSolvationStage, MoleculeInsertion

stage = PreSolvationStage(
    name="Insert sorbate molecules",
    insertions=[MoleculeInsertion("MOL.gro", n_copies=4, n_try=500)],
    solvation_top="merged.top",     # topology used from solvation onward
    extra_files=["merged.top", "MOL.gro"],  # copied into the run dir
)
setup_md_from_manifest(manifest, output_root, config=MDSetupConfig(pre_solvation_stage=stage))
```

`md_setup` stays agnostic about *what* the molecules are — a downstream workflow
builds the coordinates and merged topology and hands them in through this seam.
[**biochar-pfas**](https://github.com/layfield-lab/biochar-pfas) is one such
workflow: it parametrizes PFAS anions (PFOA/PFOS/PFBS), merges their topology
with a biochar structure, and drives this package through the seam to set up
PFAS-sorption simulations.

---

## Wood-style condensation (bulk & surface models)

`biochar.workflows.condensation` is a **parallel construction mode** that reproduces the
Wood, Mašek & Erastova (2024, *Cell Rep. Phys. Sci.* 5, 102037) protocol:
instead of a single flat sheet, it packs **many copies of one biochar molecule**
into a box and condenses them into an amorphous **bulk** via HTT-scaled simulated
annealing, then expands that bulk into an exposed **surface**. It is **setup-only**
— it writes the box, the exact Wood `.mdp` files, and run scripts; you run the
~45 ns × 3 annealing on your own GROMACS build.

### CLI

```bash
# Generate one molecule fresh, then condense 50 copies (HTT sets the anneal temperature)
biochar-condense generate --output-dir cond_600 --copies 50 --htt 600 \
    --carbons 100 --hc 0.23 --oc 0.07 --name BC600

# ...or condense copies of an existing molecule
biochar-condense from-files --output-dir cond --gro mol.gro --itp mol.itp --copies 50 --htt 600
```

Each writes a directory with `run_condensation.sh` (pack → EM → NVT → NPT anneal
→ NPT final, ×3 repeats), `run_surface.sh` (z-expand ~10 nm + semi-isotropic NPT),
and `analyze.sh` (true density, SASA, convergence, TEM).

### Python API

```python
from biochar.workflows.condensation import generate_and_condense, add_surface_and_validation

run = generate_and_condense("cond_600", n_copies=50, htt_c=600)  # pack + anneal setup
add_surface_and_validation(run)                                   # + surface + analysis
# then, on your GROMACS build:  ./run_condensation.sh -> ./run_surface.sh -> ./analyze.sh
```

### HTT-scaled annealing (Wood Tables 6 & 7)

| HTT | anneal peak T | timestep |
|---|---|---|
| 400 °C | 1000 K | 1 fs |
| 600 °C | 2000 K | 0.5 fs |
| 800 °C | 3000 K | 0.5 fs |

Per model: EM → NVT 10 ns @ peak T → NPT anneal 25 ns (hold peak → cool to 300 K
→ hold, Berendsen 100 bar) → NPT final 10 ns (300 K, 1 bar). The dry bulk is then
z-expanded and relaxed under semi-isotropic pressure coupling to form the surface.
Density/roughness are validated with `gmx freevolume` (He probe) and `gmx sasa`
(N₂ probe). This does not replace the single-molecule generator — it's an
additional mode.

---

## Testing

```bash
# Unit tests (constants, skeleton, heteroatoms, geometry, OPLS, validation, generator)
python3 -m pytest tests/test_generator.py -v

# Surface builder tests (config, geometry, GROMACS export, convenience function)
python3 -m pytest tests/test_surface_builder.py -v

# PAH quality suite (sizes 6–200C, compositions, seeds)
python3 tests/test_pah_quality.py
```

The PAH quality suite reports:
- PAH library SMILES validity and kekulizability
- Atom count accuracy (target vs actual)
- Bond errors and steric clashes per size
- Ring planarity (Å deviation from best-fit plane)
- H/C and O/C ratio accuracy across compositions

---

## Known Limitations

| Issue | Notes |
|---|---|
| H/C has a lower floor for small structures | A maximally condensed flake still caps H on every edge carbon, so very small sheets cannot go *below* ~0.44 (50 C) / ~0.36 (100 C). Requests below that floor return the most-condensed structure; use a larger sheet for lower H/C. |
| Steric clash count increases with size | Large flat aromatics have H···H near-contacts; use GROMACS energy minimisation after generation for production runs |
| Geometry validation thresholds | The built-in validator uses strict VdW radii; some reported "clashes" are artefacts of the flat starting structure and resolve under MD |
| Amorphous packing is a search and can fail | `pore_type="amorphous"` places randomly rotated sheets until none violate `min_separation`, which does not always converge for many sheets in a small box. It raises by default; set `amorphous_fallback="slit"` to degrade to a slit stack instead, which is reported and recorded in the `.gro` title. |

---

## References

- Jorgensen, W. L. et al. "Development and Testing of the OPLS All-Atom Force Field." *J. Am. Chem. Soc.* 118.45 (1996): 11225–11236.
- RDKit: Open-source cheminformatics. https://www.rdkit.org
