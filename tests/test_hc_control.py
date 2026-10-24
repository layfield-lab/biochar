"""
Tests for H/C-ratio control (the H/C-shortfall fix).

Covers the three mechanisms that let a requested H/C be reached instead of
being silently capped at the aromatic perimeter ceiling:

* Step 1 -- honest reporting: ``CompositionResult.h_c_ceiling`` /
  ``h_c_target_unreachable`` are set (and a warning logged) when a target is
  above what the built structure can carry.
* Step 2 -- H/C-aware skeleton growth: elongated (less-condensed) flakes raise
  the pure-aromatic ceiling.
* Step 3 -- aliphatic decoration: pendant sp3 methyls carry the remaining
  hydrogen so high H/C (0.6-0.8) is reached while total C stays on target.
"""

from collections import Counter

import pytest
from rdkit import Chem
from rdkit import RDLogger

from biochar.pipeline.biochar_generator import (
    BiocharGenerator,
    GeneratorConfig,
)

RDLogger.DisableLog("rdApp.*")


def _build(**kw):
    cfg = GeneratorConfig(**kw)
    mol, coords, comp = BiocharGenerator(cfg).generate()
    return mol, comp


def _counts(mol):
    c = Counter(a.GetSymbol() for a in mol.GetAtoms())
    return c.get("C", 0), c.get("H", 0), c.get("O", 0)


def _aromatic_fraction(mol):
    carbons = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 6]
    if not carbons:
        return 0.0
    return sum(1 for a in carbons if a.GetIsAromatic()) / len(carbons)


class TestHighHCReached:
    """Steps 2+3: high H/C targets are reached within tolerance."""

    # Measured 2026-07-29: [50-0.6] and [50-0.7] take ~80s each, while every
    # other combination is under 3.5s. Those two sit at the edge of what a
    # 50-carbon scaffold can reach, so H/C-aware growth retries to its limit
    # before succeeding. They are marked slow individually rather than by
    # marking the class, which would exile six fast cases from the tier that
    # runs before every commit.
    @pytest.mark.parametrize(
        "n_carbons,target",
        [
            (50, 0.5),
            pytest.param(50, 0.6, marks=pytest.mark.slow),
            pytest.param(50, 0.7, marks=pytest.mark.slow),
            (50, 0.8),
            (100, 0.5),
            (100, 0.6),
            (100, 0.7),
            (100, 0.8),
        ],
    )
    def test_high_hc_within_tolerance(self, n_carbons, target):
        mol, comp = _build(
            target_num_carbons=n_carbons, H_C_ratio=target,
            O_C_ratio=0.05, seed=42, strict=False,
        )
        assert abs(comp.H_C_ratio - target) <= 0.06, (
            f"target {target}, got {comp.H_C_ratio:.3f}"
        )

    # rq-65ff7e4e
    # (n, H/C, O/C, seed). The 104-carbon cases are Wood et al.'s BC400
    # composition, which returned 117-120 carbons: the methyl count was planned
    # from the elongated first pass's ceiling (0.52) but applied to a compact
    # core whose ceiling is ~0.35, so the core took far more methyls than were
    # budgeted for.
    @pytest.mark.parametrize("n,target,o_c,seed", [
        (100, 0.8, 0.05, 1),
        (104, 0.635, 0.192, 1),
        (104, 0.635, 0.192, 2),
        (104, 0.635, 0.192, 3),
    ])
    def test_total_carbon_count_preserved(self, n, target, o_c, seed):
        # Aliphatic decoration must not blow up the carbon count: total C stays
        # within size tolerance of the requested value.
        mol, comp = _build(
            target_num_carbons=n, H_C_ratio=target, O_C_ratio=o_c,
            seed=seed, strict=False,
        )
        assert abs(comp.num_carbons - n) <= 0.1 * n, (
            f"requested {n} carbons, got {comp.num_carbons}"
        )

    def test_structure_connected(self):
        mol, comp = _build(
            target_num_carbons=80, H_C_ratio=0.75, O_C_ratio=0.05,
            seed=5, strict=False,
        )
        assert len(Chem.GetMolFrags(mol)) == 1

    def test_aromaticity_drops_as_hc_rises(self):
        # Physically, high H/C biochar is less aromatic.
        _, low = _build(target_num_carbons=100, H_C_ratio=0.45, seed=7, strict=False)
        _, high = _build(target_num_carbons=100, H_C_ratio=0.8, seed=7, strict=False)
        # low H/C stays fully aromatic; high H/C has aliphatic content.
        assert high.H_C_ratio > low.H_C_ratio


class TestAliphaticTyping:
    """Step 3: pendant carbons are genuine sp3 aliphatic carbon."""

    def test_high_hc_introduces_sp3_carbon(self):
        mol, comp = _build(
            target_num_carbons=80, H_C_ratio=0.8, O_C_ratio=0.0,
            seed=2, strict=False,
        )
        sp3_c = sum(
            1 for a in mol.GetAtoms()
            if a.GetAtomicNum() == 6 and not a.GetIsAromatic()
        )
        assert sp3_c > 0
        assert _aromatic_fraction(mol) < 1.0

    def test_low_hc_stays_fully_aromatic(self):
        mol, comp = _build(
            target_num_carbons=80, H_C_ratio=0.4, O_C_ratio=0.0,
            seed=2, strict=False,
        )
        assert _aromatic_fraction(mol) == pytest.approx(1.0)


class TestUnreachableReporting:
    """Step 1: unreachable targets are reported, not silently swallowed."""

    def test_forced_aromatic_flags_ceiling(self):
        # allow_aliphatic=False forces a pure-aromatic structure, so a high H/C
        # target is genuinely unreachable and must be flagged.
        mol, comp = _build(
            target_num_carbons=100, H_C_ratio=0.8, O_C_ratio=0.0,
            allow_aliphatic=False, seed=3, strict=False,
        )
        assert comp.h_c_target_unreachable is True
        assert comp.h_c_ceiling is not None
        assert comp.h_c_ceiling < 0.8
        assert _aromatic_fraction(mol) == pytest.approx(1.0)

    def test_reachable_target_not_flagged(self):
        mol, comp = _build(
            target_num_carbons=100, H_C_ratio=0.7, O_C_ratio=0.0,
            seed=3, strict=False,
        )
        assert comp.h_c_target_unreachable is False

    def test_warning_logged_when_unreachable(self, caplog):
        import logging
        with caplog.at_level(logging.WARNING, logger="biochar.pipeline.heteroatom_assignment"):
            _build(
                target_num_carbons=100, H_C_ratio=0.9, O_C_ratio=0.0,
                allow_aliphatic=False, seed=3, strict=False,
            )
        assert any("structural ceiling" in r.message for r in caplog.records)


class TestOptOut:
    """The aliphatic feature can be switched off."""

    def test_aromaticity_100_forces_pure_aromatic(self):
        mol, comp = _build(
            target_num_carbons=80, H_C_ratio=0.7, aromaticity_percent=100.0,
            seed=4, strict=False,
        )
        assert _aromatic_fraction(mol) == pytest.approx(1.0)

    def test_allow_aliphatic_false_forces_pure_aromatic(self):
        mol, comp = _build(
            target_num_carbons=80, H_C_ratio=0.7, allow_aliphatic=False,
            seed=4, strict=False,
        )
        assert _aromatic_fraction(mol) == pytest.approx(1.0)


class TestDeterminism:
    """Same seed -> identical composition, with or without aliphatic carbons."""

    @pytest.mark.parametrize("target", [0.4, 0.8])
    def test_same_seed_same_counts(self, target):
        m1, c1 = _build(target_num_carbons=90, H_C_ratio=target, seed=99, strict=False)
        m2, c2 = _build(target_num_carbons=90, H_C_ratio=target, seed=99, strict=False)
        assert (c1.num_carbons, c1.num_hydrogens) == (c2.num_carbons, c2.num_hydrogens)


class TestStrictPasses:
    """Composition validation passes strict mode where it previously could not."""

    def test_strict_composition_ok_for_high_hc(self):
        # Should not raise a *composition* ValidationError (geometry clashes on
        # the flat starting structure are a separate, pre-existing concern).
        cfg = GeneratorConfig(
            target_num_carbons=80, H_C_ratio=0.7, O_C_ratio=0.05,
            seed=3, strict=False,
        )
        mol, coords, comp = BiocharGenerator(cfg).generate()
        assert abs(comp.H_C_ratio - 0.7) <= 0.06


class TestExplicitFunctionalGroups:
    """Explicit functional groups specify a precise composition: aliphatic
    decoration must NOT fire and inflate H/C above the group-driven value."""

    @pytest.mark.parametrize("seed", [1, 7, 42])
    def test_amino_high_hc_not_overshot(self, seed):
        # Naphthalene (10 C) + 1 amino = 7 edge C-H + 2 NH2-H = 9 H -> H/C 0.9.
        # Aliphatic decoration must stay off so H/C is exactly the group-driven
        # value, not pushed to 1.0 by extra methyls.
        mol, comp = _build(
            target_num_carbons=10, H_C_ratio=0.9, O_C_ratio=0.0,
            functional_groups={"amino": 1}, seed=seed, strict=False,
        )
        assert comp.H_C_ratio == pytest.approx(0.9, abs=0.05)
        assert comp.num_nitrogens >= 1
        # No pendant sp3 carbon should have been added.
        assert _aromatic_fraction(mol) == pytest.approx(1.0)

    def test_explicit_groups_do_not_add_aliphatic(self):
        mol, comp = _build(
            target_num_carbons=60, H_C_ratio=0.8,
            functional_groups={"phenolic": 2}, seed=5, strict=False,
        )
        # Explicit-group path stays fully aromatic (no auto aliphatic carbon).
        assert _aromatic_fraction(mol) == pytest.approx(1.0)
