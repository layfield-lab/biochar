"""
Default-path baseline: linked multi-unit structures (rq-PENDING) are opt-in, and
a request that does not opt in must build exactly what it built before.

The fixture was captured from the generator *before* the linked-unit hooks
existed (main at 02b2574). The comparison is of graph topology in atom-index
order, OPLS types and composition, not of canonical SMILES or coordinates: CI
runs several RDKit releases, and both of those can move between releases on
their own, which would turn this into a test of RDKit rather than of the hooks.
A hook that consumed the RNG differently, reordered atoms, or placed one more
group changes the topology, so it still fails here.

Regenerate only when a default-path change is intended:
    PYTHONPATH=src python tests/test_default_path_baseline.py --write
"""

import json
import sys
from pathlib import Path

import pytest

from biochar.pipeline.biochar_generator import BiocharGenerator, GeneratorConfig

FIXTURE = Path(__file__).parent / "fixtures" / "default_path_baseline.json"

# Both sides of the 80-heavy-atom embedding split, temperature mode on and off,
# explicit functional groups on and off. strict=False: the baseline records what
# the default path builds, including structures strict validation would reject.
CASES = {
    "default_50": dict(target_num_carbons=50),
    "hc045_100": dict(target_num_carbons=100, H_C_ratio=0.45, O_C_ratio=0.05),
    "wood_400_80": dict(target_num_carbons=80, temperature=400, feedstock="wood"),
    "groups_40": dict(target_num_carbons=40, functional_groups={"phenolic": 2, "carboxyl": 1}),
}
SEEDS = (1, 2)


def _signature(kwargs: dict, seed: int) -> dict:
    gen = BiocharGenerator(GeneratorConfig(seed=seed, strict=False, **kwargs))
    gen.generate()
    mol = gen.mol
    return {
        "atoms": [[a.GetAtomicNum(), a.GetFormalCharge()] for a in mol.GetAtoms()],
        "bonds": sorted(
            [min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
             max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
             b.GetBondTypeAsDouble()]
            for b in mol.GetBonds()
        ),
        "types": [gen.atom_types[i] for i in range(mol.GetNumAtoms())],
        "formula": gen.composition.molecular_formula,
    }


def _capture() -> dict:
    return {f"{name}/seed{seed}": _signature(kw, seed) for name, kw in CASES.items() for seed in SEEDS}


@pytest.fixture(scope="module")
def baseline():
    return json.loads(FIXTURE.read_text())


@pytest.mark.parametrize("name", sorted(CASES))
@pytest.mark.parametrize("seed", SEEDS)
def test_default_path_matches_pre_change_baseline(baseline, name, seed):
    expected = baseline[f"{name}/seed{seed}"]
    actual = _signature(CASES[name], seed)
    assert actual["formula"] == expected["formula"]
    assert actual["atoms"] == expected["atoms"]
    assert actual["bonds"] == expected["bonds"]
    assert actual["types"] == expected["types"]


if __name__ == "__main__" and "--write" in sys.argv:
    FIXTURE.write_text(json.dumps(_capture(), separators=(",", ":")) + "\n")
    print(f"wrote {FIXTURE}")
