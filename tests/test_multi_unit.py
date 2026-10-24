"""
Linked multi-unit ("island") structures -- rqm/linked-units.md.

The architecture is what these tests pin, not only validity: a merged ring
system, an oxygenated linker or a coplanar layout all export and pass grompp,
so each is checked directly on the generated molecule.
"""

import logging
import math
import random

import numpy as np
import pytest
from rdkit import Chem

from biochar.constants import LATTICE_MAX_SHEET_RMS
from biochar.pipeline.biochar_generator import (
    BiocharGenerator,
    GeneratorConfig,
    generate_biochar,
)
from biochar.pipeline.geometry_3d import GeometryValidator
from biochar.pipeline.multi_unit import LINKER_TYPES, plan_units

BC400 = dict(target_num_carbons=104, H_C_ratio=0.635, O_C_ratio=0.19,
             temperature=400, feedstock="wood")


def _generate(strict=False, seed=1, **kw):
    gen = BiocharGenerator(GeneratorConfig(seed=seed, strict=strict, linked_units=True, **kw))
    gen.generate()
    return gen


def _ring_systems(mol):
    """Rings merged when they share an atom -- the measure used on Wood's blocks."""
    mol = Chem.Mol(mol)
    Chem.FastFindRings(mol)
    systems = []
    for ring in map(set, mol.GetRingInfo().AtomRings()):
        merged = [s for s in systems if s & ring]
        for s in merged:
            systems.remove(s)
            ring |= s
        systems.append(ring)
    return systems


def _unit_atoms(mol, uid):
    return [a.GetIdx() for a in mol.GetAtoms() if a.HasProp("unit_id") and a.GetIntProp("unit_id") == uid]


def _plane_normal(coords):
    centred = coords - coords.mean(axis=0)
    return np.linalg.svd(centred, full_matrices=False)[2][-1], centred


def _single_linker(linker, seed=1, **kw):
    params = dict(target_num_carbons=60, H_C_ratio=0.6, O_C_ratio=0.1,
                  num_units=3, linker_mix={linker: 1.0})
    params.update(kw)
    return _generate(seed=seed, **params)


# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------

class TestConfiguration:
    def test_off_by_default_and_round_trips(self):
        cfg = GeneratorConfig()
        assert cfg.linked_units is False and cfg.num_units is None and cfg.linker_mix is None
        cfg = GeneratorConfig(linked_units=True, num_units=3, linker_mix={"biaryl": 1.0})
        assert GeneratorConfig.from_dict(cfg.to_dict()) == cfg

    # rq-52b8fb84
    @pytest.mark.parametrize("kw", [dict(num_units=3), dict(linker_mix={"biaryl": 1.0})])
    def test_linked_parameters_without_the_opt_in_are_refused(self, kw):
        with pytest.raises(ValueError, match="linked_units=True"):
            GeneratorConfig(**kw)

    @pytest.mark.parametrize("kw", [
        dict(num_units=0),
        dict(num_units=2.5),
        dict(linker_mix={"ketone": 1.0}),
        dict(linker_mix={"biaryl": 0.0}),
        dict(linker_mix={"biaryl": -1.0, "propyl": 2.0}),
        dict(num_units=6, target_num_carbons=30),
    ])
    def test_unusable_linked_parameters_raise(self, kw):
        with pytest.raises(ValueError):
            GeneratorConfig(linked_units=True, **kw)


# --------------------------------------------------------------------------
# Planning
# --------------------------------------------------------------------------

class TestPlanning:
    # rq-e918b504
    @pytest.mark.parametrize("temperature,allowed", [
        (400, {4, 5}), (600, {3}), (800, {1}), (None, {4}),
    ])
    def test_unit_count_follows_the_temperature_band(self, temperature, allowed):
        for seed in range(5):
            cfg = GeneratorConfig(linked_units=True, target_num_carbons=300,
                                  temperature=temperature, feedstock="wood" if temperature else None,
                                  seed=seed)
            plan = plan_units(cfg, random.Random(seed))
            assert plan.num_units in allowed
            assert len(plan.linkers) == plan.num_units - 1
            assert set(plan.linkers) <= set(LINKER_TYPES)

    # rq-1ca67e8f
    def test_explicit_count_and_mix_override_the_band(self):
        cfg = GeneratorConfig(linked_units=True, target_num_carbons=104, temperature=400,
                              num_units=3, linker_mix={"biaryl": 1.0})
        plan = plan_units(cfg, random.Random(0))
        assert plan.num_units == 3
        assert plan.linkers == ["biaryl", "biaryl"]

    # rq-85b48c1e
    def test_a_count_that_cannot_reach_the_target_hc_shrinks(self, caplog):
        cfg = GeneratorConfig(linked_units=True, target_num_carbons=50,
                              temperature=600, feedstock="wood")
        with caplog.at_level(logging.WARNING, logger="biochar.pipeline.multi_unit"):
            plan = plan_units(cfg, random.Random(0))
        assert plan.num_units < 3
        assert "cannot reach H/C" in caplog.text
        # The same band at Wood's BC600 size can reach it and keeps 3.
        big = GeneratorConfig(linked_units=True, target_num_carbons=272,
                              temperature=600, feedstock="wood")
        assert plan_units(big, random.Random(0)).num_units == 3

    # rq-d9b926a8
    def test_high_temperature_builds_the_fused_sheet_unchanged(self):
        kw = dict(target_num_carbons=60, temperature=800, feedstock="wood", seed=2, strict=False)
        plain = BiocharGenerator(GeneratorConfig(**kw))
        plain.generate()
        linked = BiocharGenerator(GeneratorConfig(linked_units=True, **kw))
        linked.generate()
        assert linked.linked_units is None
        assert Chem.MolToSmiles(linked.mol) == Chem.MolToSmiles(plain.mol)
        assert np.allclose(linked.coords, plain.coords)


# --------------------------------------------------------------------------
# Architecture
# --------------------------------------------------------------------------

class TestArchitecture:
    # rq-2b213b2e
    @pytest.mark.parametrize("k,n", [(3, 60), (4, 80), (5, 104)])
    def test_one_ring_system_per_unit(self, k, n):
        gen = _generate(target_num_carbons=n, H_C_ratio=0.6, O_C_ratio=0.1, num_units=k)
        assert len(_ring_systems(gen.mol)) == k
        assert gen.linked_units.num_units == k
        assert sum(gen.linked_units.linker_counts.values()) == k - 1
        sizes = gen.linked_units.unit_carbons
        assert min(sizes) >= 6 and sizes[0] == max(sizes)

    # rq-6f9371cd
    @pytest.mark.parametrize("linker", LINKER_TYPES)
    def test_each_linker_builds_its_motif(self, linker):
        mol = _single_linker(linker).mol
        units = {a.GetIdx(): a.GetIntProp("unit_id") for a in mol.GetAtoms() if a.HasProp("unit_id")}
        bridges = {}
        for a in mol.GetAtoms():
            if a.HasProp("linker_id"):
                bridges.setdefault(a.GetIntProp("linker_id"), []).append(a)
        direct = [
            b for b in mol.GetBonds()
            if b.GetBeginAtomIdx() in units and b.GetEndAtomIdx() in units
            and units[b.GetBeginAtomIdx()] != units[b.GetEndAtomIdx()]
        ]
        if linker == "biaryl":
            assert len(direct) == 2 and not bridges
            for b in direct:
                assert b.GetBondType() == Chem.BondType.SINGLE
                assert b.GetBeginAtom().GetIsAromatic() and b.GetEndAtom().GetIsAromatic()
            return
        assert not direct and len(bridges) == 2
        expected = {"diaryl_ether": [8], "aryl_alkyl_ether": [6, 8], "propyl": [6, 6, 6]}[linker]
        for atoms in bridges.values():
            assert sorted(a.GetAtomicNum() for a in atoms) == expected
            ends = [n for a in atoms for n in a.GetNeighbors() if n.GetIdx() in units]
            assert len(ends) == 2 and all(n.GetIsAromatic() for n in ends)
            assert units[ends[0].GetIdx()] != units[ends[1].GetIdx()]
            for a in atoms:
                if a.GetAtomicNum() == 6:
                    assert not a.GetIsAromatic()
                    assert sum(1 for n in a.GetNeighbors() if n.GetAtomicNum() == 1) == 2

    # rq-b7b7c2cd
    @pytest.mark.parametrize("linker", ["aryl_alkyl_ether", "propyl"])
    def test_linker_carbons_are_never_oxygenated(self, linker):
        # High O/C on small units: the aromatic edges cannot hold it all, so
        # the oxygen stage spills onto sp3 carbons wherever it can.
        gen = _single_linker(linker, target_num_carbons=40, O_C_ratio=0.35, H_C_ratio=0.8,
                             num_units=4)
        mol = gen.mol
        for a in mol.GetAtoms():
            if a.HasProp("linker_id") and a.GetAtomicNum() == 6:
                oxygens = [n for n in a.GetNeighbors() if n.GetAtomicNum() == 8]
                assert all(n.HasProp("linker_id") for n in oxygens), (
                    f"linker carbon {a.GetIdx()} carries a non-bridge oxygen"
                )

    # rq-2b213b2e -- the smallest accepted budget: a benzene core runs out of
    # attachment sites, and the build must move on rather than crash.
    @pytest.mark.parametrize("seed", range(10))
    def test_five_units_at_the_smallest_budget_build(self, seed):
        gen = _generate(seed=seed, target_num_carbons=30, num_units=5, H_C_ratio=0.8,
                        O_C_ratio=0.05, linker_mix={"biaryl": 1.0})
        assert len(_ring_systems(gen.mol)) == 5

    # rq-a9db6c13 -- seed=None is a fresh structure per call, as on the fused path.
    def test_unseeded_linked_builds_vary(self):
        kw = dict(target_num_carbons=60, H_C_ratio=0.45, O_C_ratio=0.1, strict=False,
                  linked_units=True, num_units=3, linker_mix={"diaryl_ether": 1.0})
        coords = []
        for _ in range(3):
            gen = BiocharGenerator(GeneratorConfig(seed=None, **kw))
            gen.generate()
            coords.append(gen.coords)
        assert not all(c.shape == coords[0].shape and np.allclose(c, coords[0]) for c in coords[1:])

    # rq-7705e753
    @pytest.mark.parametrize("seed", [1, 2, 3])
    def test_explicit_bridging_groups_do_not_merge_units(self, seed):
        gen = _generate(seed=seed, target_num_carbons=80, num_units=3,
                        linker_mix={"biaryl": 1.0},
                        functional_groups={"ether": 3, "phenolic": 2})
        assert len(_ring_systems(gen.mol)) == 3


# --------------------------------------------------------------------------
# Composition
# --------------------------------------------------------------------------

class TestComposition:
    # rq-227a0ddd
    @pytest.mark.parametrize("seed", [1, 2, 3])
    def test_bc400_like_request_passes_strict(self, seed):
        gen = _generate(strict=True, seed=seed, **BC400)
        assert len(_ring_systems(gen.mol)) in (4, 5)
        assert gen.validation_report[0]

    # rq-4115991e
    @pytest.mark.slow
    @pytest.mark.parametrize("seed", [1, 2, 3])
    def test_bc600_like_request_passes_strict(self, seed):
        gen = _generate(strict=True, seed=seed, target_num_carbons=272,
                        temperature=600, feedstock="wood")
        assert len(_ring_systems(gen.mol)) == 3

    # rq-8ff54b34
    def test_linker_oxygen_counts_toward_oc(self):
        gen = _single_linker("diaryl_ether", target_num_carbons=80, num_units=5,
                             O_C_ratio=0.1, H_C_ratio=0.55)
        comp = gen.composition
        assert gen.linked_units.num_linker_oxygens == 4
        assert abs(comp.O_C_ratio - 0.1) <= 0.1 * 0.1 + 1e-9
        assert comp.num_oxygens == sum(1 for a in gen.mol.GetAtoms() if a.GetAtomicNum() == 8)


# --------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------

class TestGeometry:
    @pytest.fixture(scope="class")
    def bc400(self):
        return [_generate(seed=s, **BC400) for s in (1, 2, 3)]

    # rq-3a47b3fb
    def test_no_clashes_or_bond_length_errors(self, bc400):
        for gen in bc400:
            valid, errors = GeometryValidator.validate_geometry(gen.mol, gen.coords)
            assert valid, errors[:5]

    # rq-6d4235f9
    def test_every_unit_stays_flat(self, bc400):
        for gen in bc400:
            for uid in range(gen.linked_units.num_units):
                normal, centred = _plane_normal(gen.coords[_unit_atoms(gen.mol, uid)])
                rms = float(np.sqrt(np.mean((centred @ normal) ** 2)))
                assert rms <= LATTICE_MAX_SHEET_RMS, (uid, rms)

    # rq-dacb8d13
    @pytest.mark.parametrize("seed", [1, 2, 3])
    def test_biaryl_links_are_twisted(self, seed):
        gen = _single_linker("biaryl", seed=seed)
        mol, coords = gen.mol, gen.coords
        for b in mol.GetBonds():
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            ai, aj = mol.GetAtomWithIdx(i), mol.GetAtomWithIdx(j)
            if not (ai.HasProp("unit_id") and aj.HasProp("unit_id")):
                continue
            ui, uj = ai.GetIntProp("unit_id"), aj.GetIntProp("unit_id")
            if ui == uj:
                continue
            ni, _ = _plane_normal(coords[_unit_atoms(mol, ui)])
            nj, _ = _plane_normal(coords[_unit_atoms(mol, uj)])
            # Planes have no direction, so the angle between them is 0-90°;
            # "between 20 and 160 degrees apart" is at least 20 here.
            angle = math.degrees(math.acos(min(1.0, abs(float(ni @ nj)))))
            assert angle >= 20.0, f"biaryl rings only {angle:.1f} degrees apart"

    # rq-a9db6c13
    def test_embedding_is_reproducible(self):
        a = _generate(seed=4, **BC400)
        b = _generate(seed=4, **BC400)
        assert Chem.MolToSmiles(a.mol) == Chem.MolToSmiles(b.mol)
        assert np.array_equal(a.coords, b.coords)


# --------------------------------------------------------------------------
# Typing and provenance
# --------------------------------------------------------------------------

class TestTypingAndRecords:
    # rq-a7ae7f39
    @pytest.mark.parametrize("linker", LINKER_TYPES)
    def test_linker_atoms_take_their_opls_types(self, linker):
        gen = _single_linker(linker)
        mol, types = gen.mol, gen.atom_types
        for a in mol.GetAtoms():
            if a.HasProp("linker_id"):
                if a.GetAtomicNum() == 8:
                    assert types[a.GetIdx()] == "OS"
                else:
                    assert types[a.GetIdx()] == "CT"
                    hs = [n.GetIdx() for n in a.GetNeighbors() if n.GetAtomicNum() == 1]
                    assert hs and all(types[h] == "HC" for h in hs)
            if a.HasProp("unit_id") and any(
                n.HasProp("linker_id") or (n.HasProp("unit_id") and n.GetIntProp("unit_id") != a.GetIntProp("unit_id"))
                for n in a.GetNeighbors()
            ):
                assert types[a.GetIdx()] == "CA"

    # rq-73fdab80
    def test_topology_headers_name_units_and_links(self, tmp_path):
        gen = _generate(**BC400)
        _gro, top, itp = gen.export_gromacs(output_directory=str(tmp_path), basename="isl")
        census = gen.linked_units.describe()
        for path in (top, itp):
            assert census in open(path).read()

    # rq-420e40ce
    def test_result_carries_the_census(self):
        result = generate_biochar(seed=1, write_files=False, linked_units=True, **BC400)
        assert result.linked_units is not None
        assert result.linked_units.num_units == len(_ring_systems(result.mol))
        mol, coords, gro, top, itp = result  # positional unpacking unchanged
        assert gro is None
