"""
Tests for hex-lattice crowding relief (rqm/geometry-embedding.md, rq-8e5f5dc1).

The hex-lattice path places ring carbons on an exact flat lattice and every
substituent radially. At a cove -- two edge carbons 2.46 A apart, four bonds
apart -- both hydrogens point into the same vacant lattice site. The fixture is
benzo[c]phenanthrene, the smallest PAH with a cove, built exactly the way the
flat path builds a sheet: carbons on the lattice, each H 1.09 A out along the
bisector of its carbon's two ring neighbours.
"""

import itertools
import random
import re

import numpy as np
import pytest

from rdkit import Chem

import biochar.pipeline.geometry_3d as geometry_3d
from biochar.constants import LATTICE_MAX_SHEET_RMS
from biochar.pipeline.biochar_generator import (
    BiocharGenerator,
    GeneratorConfig,
    ValidationError,
)
from biochar.pipeline.carbon_skeleton import PAHAssembler
from biochar.pipeline.geometry_3d import (
    CoordinateGenerator,
    GeometryValidator,
    _bent_hydrogen_position,
    _clash_pairs,
    _substituent_groups,
    _tilt_crowded_substituents,
)
from biochar.pipeline.heteroatom_assignment import attach_aliphatic_carbons

LATTICE_BOND = 1.42


def _flat_pah_on_lattice(ring_centres):
    """A PAH whose rings are hexagons centred at *ring_centres* (2D), flat at
    z = 0, with radially placed hydrogens -- the flat path's placement rule."""
    vertices = []

    def vertex_index(p):
        for k, q in enumerate(vertices):
            if np.linalg.norm(p - q) < 1e-3:
                return k
        vertices.append(p)
        return len(vertices) - 1

    rings = []
    for cx, cy in ring_centres:
        ring = []
        for k in range(6):
            a = np.radians(30.0 + 60.0 * k)
            ring.append(vertex_index(np.array([cx, cy]) + LATTICE_BOND * np.array([np.cos(a), np.sin(a)])))
        rings.append(ring)

    mol = Chem.RWMol()
    for _ in vertices:
        atom = Chem.Atom(6)
        atom.SetIsAromatic(True)
        mol.AddAtom(atom)
    bonded = set()
    for ring in rings:
        for k in range(6):
            i, j = sorted((ring[k], ring[(k + 1) % 6]))
            if (i, j) not in bonded:
                bonded.add((i, j))
                mol.AddBond(i, j, Chem.BondType.AROMATIC)

    coords = [np.array([p[0], p[1], 0.0]) for p in vertices]
    for c in range(len(vertices)):
        nbrs = [n.GetIdx() for n in mol.GetAtomWithIdx(c).GetNeighbors()]
        if len(nbrs) != 2:
            continue
        h = mol.AddAtom(Chem.Atom(1))
        mol.AddBond(c, h, Chem.BondType.SINGLE)
        out = coords[c] - np.mean([coords[n] for n in nbrs], axis=0)
        coords.append(coords[c] + 1.09 * out / np.linalg.norm(out))

    mol = mol.GetMol()
    Chem.SanitizeMol(mol)
    coords = np.array(coords)
    conformer = Chem.Conformer(mol.GetNumAtoms())
    for i, xyz in enumerate(coords):
        conformer.SetAtomPosition(i, xyz.tolist())
    mol.AddConformer(conformer, assignId=True)
    return mol, coords


def _benzo_c_phenanthrene():
    """[4]helicene: four rings, each fused one step further round the same way."""
    step = LATTICE_BOND * np.sqrt(3.0)
    centres, p = [np.zeros(2)], np.zeros(2)
    for turn in (0.0, 60.0, 120.0):
        p = p + step * np.array([np.cos(np.radians(turn)), np.sin(np.radians(turn))])
        centres.append(p)
    return _flat_pah_on_lattice(centres)


def _cove_carbons(mol, coords):
    """The two edge carbons 2.46 A apart that are four bonds apart."""
    dm = Chem.GetDistanceMatrix(mol)
    carbons = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6]
    for i in carbons:
        for j in carbons:
            if i < j and dm[i][j] == 4 and abs(np.linalg.norm(coords[i] - coords[j]) - 2.46) < 0.02:
                return i, j
    raise AssertionError("fixture has no cove")


def _sheet_rms(mol, coords):
    ring = [a.GetIdx() for a in mol.GetAtoms() if a.IsInRing()]
    centred = coords[ring] - coords[ring].mean(axis=0)
    normal = np.linalg.svd(centred, full_matrices=False)[2][-1]
    return float(np.sqrt(np.mean((centred @ normal) ** 2)))


@pytest.fixture
def cove():
    mol, coords = _benzo_c_phenanthrene()
    assert mol.GetNumAtoms() == 30  # C18H12
    assert _clash_pairs(mol, coords), "the flat placement must crowd the cove"
    return mol, coords


class TestCoveRelief:
    # rq-57d3fc46
    def test_cove_substituents_are_pulled_apart(self, cove):
        mol, coords = cove
        relieved = CoordinateGenerator(seed=0).relieve_lattice_crowding(mol, coords)
        assert _clash_pairs(mol, relieved) == []

    # rq-da5ef2fc
    def test_cove_opens_rather_than_only_moving_its_hydrogens(self, cove):
        mol, coords = cove
        i, j = _cove_carbons(mol, coords)
        relieved = CoordinateGenerator(seed=0).relieve_lattice_crowding(mol, coords)
        assert np.linalg.norm(relieved[i] - relieved[j]) > 2.46 + 0.1

    # rq-13532295
    def test_relieved_sheet_stays_flat_and_keeps_its_bonds(self, cove):
        mol, coords = cove
        relieved = CoordinateGenerator(seed=0).relieve_lattice_crowding(mol, coords)
        assert _sheet_rms(mol, relieved) <= LATTICE_MAX_SHEET_RMS
        for bond in mol.GetBonds():
            if bond.GetIsAromatic():
                length = np.linalg.norm(
                    relieved[bond.GetBeginAtomIdx()] - relieved[bond.GetEndAtomIdx()]
                )
                assert length == pytest.approx(LATTICE_BOND, abs=0.05)

    def test_relief_is_deterministic(self, cove):
        mol, coords = cove
        a = CoordinateGenerator(seed=3).relieve_lattice_crowding(Chem.Mol(mol), coords)
        b = CoordinateGenerator(seed=3).relieve_lattice_crowding(Chem.Mol(mol), coords)
        assert np.array_equal(a, b)

    def test_relieved_coordinates_are_carried_by_the_molecule(self, cove):
        mol, coords = cove
        relieved = CoordinateGenerator(seed=0).relieve_lattice_crowding(mol, coords)
        assert np.allclose(mol.GetConformer().GetPositions(), relieved)


class TestReliefNeverMakesThingsWorse:
    # rq-ac64435a
    def test_a_worse_relaxation_is_discarded(self, cove, monkeypatch):
        mol, coords = cove

        def wrecked(mol_, coords_, seed=None):
            broken = coords_.copy()
            broken[0] += np.array([3.0, 0.0, 0.0])  # tear a ring bond
            return broken

        monkeypatch.setattr(geometry_3d, "_tethered_lattice_relax", wrecked)
        relieved = CoordinateGenerator(seed=0).relieve_lattice_crowding(mol, coords)
        assert np.array_equal(relieved, _tilt_crowded_substituents(mol, coords))

    def test_a_relaxation_that_cannot_be_set_up_falls_back_to_the_tilt(self, cove, monkeypatch):
        mol, coords = cove
        monkeypatch.setattr(geometry_3d, "_tethered_lattice_relax", lambda *a, **k: None)
        relieved = CoordinateGenerator(seed=0).relieve_lattice_crowding(mol, coords)
        assert np.array_equal(relieved, _tilt_crowded_substituents(mol, coords))

    # rq-37a3d53e
    def test_a_sheet_with_nothing_to_relieve_is_left_alone(self):
        # Naphthalene has no cove: its flat lattice placement is already clean.
        mol, coords = _flat_pah_on_lattice([(0.0, 0.0), (LATTICE_BOND * np.sqrt(3.0), 0.0)])
        assert _clash_pairs(mol, coords) == []
        relieved = CoordinateGenerator(seed=0).relieve_lattice_crowding(mol, coords)
        assert np.array_equal(relieved, coords)


class TestOneDefinitionOfAClash:
    # rq-5079b4df -- the validator and the relief share one clash definition
    # rq-4eb90e63
    def test_validator_reports_exactly_the_clash_finders_pairs(self, cove):
        mol, coords = cove
        reported = [
            tuple(int(x) for x in re.search(r"atoms (\d+) and (\d+)", e).groups())
            for e in GeometryValidator._check_steric_clashes(mol, coords)
        ]
        assert reported == [(i, j) for i, j, _d, _f in _clash_pairs(mol, coords)]
        assert reported, "fixture must have clashes for this to mean anything"


def _generate_on_hex_lattice(monkeypatch, **config):
    """Generate a structure, asserting it really took the relief path."""
    calls = {"n": 0}
    original = CoordinateGenerator.relieve_lattice_crowding

    def spy(self, mol, coords):
        calls["n"] += 1
        return original(self, mol, coords)

    monkeypatch.setattr(CoordinateGenerator, "relieve_lattice_crowding", spy)
    result = BiocharGenerator(GeneratorConfig(**config)).generate()
    assert calls["n"] == 1, "expected the hex-lattice relief path"
    return result


class TestHydroxylPlacement:
    # rq-7515b128 -- a hydroxyl hydrogen keeps its bond angle
    # rq-653a5fd2
    def test_hydroxyls_on_a_large_sheet_are_bent(self, monkeypatch):
        mol, coords, _comp = _generate_on_hex_lattice(
            monkeypatch, target_num_carbons=100, H_C_ratio=0.35, O_C_ratio=0.1,
            seed=2, strict=False,
        )
        angles = []
        for o in mol.GetAtoms():
            if o.GetAtomicNum() != 8:
                continue
            hs = [n.GetIdx() for n in o.GetNeighbors() if n.GetAtomicNum() == 1]
            heavy = [n.GetIdx() for n in o.GetNeighbors() if n.GetAtomicNum() != 1]
            if len(hs) == 1 and len(heavy) == 1:
                u = coords[hs[0]] - coords[o.GetIdx()]
                v = coords[heavy[0]] - coords[o.GetIdx()]
                cos = u @ v / (np.linalg.norm(u) * np.linalg.norm(v))
                angles.append(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))
        assert angles, "expected hydroxyl groups on this structure"
        assert all(100.0 <= a <= 120.0 for a in angles), sorted(angles)


class TestStrictValidationOnLargeSheets:
    # rq-9a59d62c
    @pytest.mark.slow
    @pytest.mark.parametrize("o_c", [0.0, 0.1, 0.2])
    def test_no_steric_clash_failures_at_100_carbons(self, monkeypatch, o_c):
        # H/C 0.35 is reachable at 100 C, so composition cannot be what fails;
        # before relief every one of these seeds failed on a clash.
        for seed in range(1, 6):
            try:
                _generate_on_hex_lattice(
                    monkeypatch, target_num_carbons=100, H_C_ratio=0.35,
                    O_C_ratio=o_c, seed=seed,
                )
            except ValidationError as e:
                assert "Steric clash" not in str(e), f"seed {seed}: {str(e)[:300]}"


def _methyl_decorated_sheet(n_methyls=6, seed=0):
    """A >80-carbon hex-lattice sheet carrying pendant methyls, built the way
    aliphatic decoration builds one: the methyl carbons have no lattice tag."""
    core = PAHAssembler(seed=seed).generate(target_num_carbons=100).mol
    return attach_aliphatic_carbons(core, n_methyls, random.Random(seed))


def _angle(a, centre, b):
    u, v = a - centre, b - centre
    cos = u @ v / (np.linalg.norm(u) * np.linalg.norm(v))
    return float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))


class TestDecoratedSheetsKeepTheLattice:
    # rq-d6926052 -- a decorated sheet stays on the lattice
    # rq-06dfeecf
    def test_methyl_decorated_sheet_takes_the_hex_lattice_path(self):
        mol = _methyl_decorated_sheet()
        untagged = [a.GetIdx() for a in mol.GetAtoms()
                    if a.GetAtomicNum() == 6 and not a.HasProp("init_x")]
        assert len(untagged) == 6, "fixture must carry untagged methyl carbons"

        gen = CoordinateGenerator(seed=0)
        out, coords = gen.generate_3d_coordinates(mol)

        assert gen.used_hex_lattice
        for atom in out.GetAtoms():
            if atom.HasProp("init_x"):
                lattice = (atom.GetDoubleProp("init_x"), atom.GetDoubleProp("init_y"))
                assert coords[atom.GetIdx()][:2] == pytest.approx(lattice, abs=1e-6)

    # rq-91985117
    def test_methyls_are_placed_tetrahedrally(self):
        out, coords = CoordinateGenerator(seed=0).generate_3d_coordinates(
            _methyl_decorated_sheet()
        )
        assert GeometryValidator._validate_bond_lengths(out, coords) == []
        methyls = [a for a in out.GetAtoms()
                   if a.GetAtomicNum() == 6 and not a.HasProp("init_x")]
        for c in methyls:
            ci = c.GetIdx()
            nbrs = [n.GetIdx() for n in c.GetNeighbors()]
            hs = [n for n in nbrs if out.GetAtomWithIdx(n).GetAtomicNum() == 1]
            assert len(hs) == 3
            for k, i in enumerate(nbrs):
                for j in nbrs[k + 1:]:
                    assert _angle(coords[i], coords[ci], coords[j]) == pytest.approx(109.47, abs=2.0)

    # rq-d6926052 -- a phenolic C-O keeps 1.36 Angstrom; the embedding's
    # de-aromatised working copy must not decide it is aliphatic.
    def test_phenolic_oxygen_keeps_its_bond_length(self):
        sheet = Chem.RWMol(_methyl_decorated_sheet())
        edge = [a.GetIdx() for a in sheet.GetAtoms() if a.GetIsAromatic() and a.GetDegree() < 3]
        for site in edge[:3]:
            sheet.AddBond(site, sheet.AddAtom(Chem.Atom(8)), Chem.BondType.SINGLE)
        mol = sheet.GetMol()
        mol.UpdatePropertyCache(strict=False)

        gen = CoordinateGenerator(seed=0)
        out, coords = gen.generate_3d_coordinates(mol)

        assert gen.used_hex_lattice
        lengths = [
            np.linalg.norm(coords[b.GetBeginAtomIdx()] - coords[b.GetEndAtomIdx()])
            for b in out.GetBonds()
            if {b.GetBeginAtom().GetAtomicNum(), b.GetEndAtom().GetAtomicNum()} == {6, 8}
        ]
        assert len(lengths) == 3
        assert lengths == pytest.approx([1.36] * 3, abs=1e-6)

    # rq-91985117 -- deeper off-lattice groups: -CH2-OH and a carboxyl
    def test_hydroxymethyl_and_carboxyl_keep_their_angles(self):
        sheet = Chem.RWMol(_methyl_decorated_sheet())
        methyl = next(a.GetIdx() for a in sheet.GetAtoms()
                      if a.GetAtomicNum() == 6 and not a.HasProp("init_x"))
        hydroxyl_o = sheet.AddAtom(Chem.Atom(8))
        sheet.AddBond(methyl, hydroxyl_o, Chem.BondType.SINGLE)
        site = next(a.GetIdx() for a in sheet.GetAtoms() if a.GetIsAromatic() and a.GetDegree() < 3)
        acid_c = sheet.AddAtom(Chem.Atom(6))
        sheet.AddBond(site, acid_c, Chem.BondType.SINGLE)
        acid_o = [sheet.AddAtom(Chem.Atom(8)) for _ in range(2)]
        sheet.AddBond(acid_c, acid_o[0], Chem.BondType.DOUBLE)
        sheet.AddBond(acid_c, acid_o[1], Chem.BondType.SINGLE)
        mol = sheet.GetMol()
        mol.UpdatePropertyCache(strict=False)

        gen = CoordinateGenerator(seed=0)
        out, coords = gen.generate_3d_coordinates(mol)

        assert gen.used_hex_lattice
        assert GeometryValidator._validate_bond_lengths(out, coords) == []
        ring_c = next(n.GetIdx() for n in out.GetAtomWithIdx(methyl).GetNeighbors() if n.IsInRing())
        assert _angle(coords[ring_c], coords[methyl], coords[hydroxyl_o]) == pytest.approx(109.47, abs=2.0)
        assert _angle(coords[acid_o[0]], coords[acid_c], coords[acid_o[1]]) == pytest.approx(120.0, abs=2.0)
        assert _angle(coords[site], coords[acid_c], coords[acid_o[0]]) == pytest.approx(120.0, abs=2.0)
        # The C-O-H angle of a crowded -CH2-OH is TestCrowdedHydroxyls' job.

    # rq-91985117
    def test_a_methyl_is_a_substituent_the_relief_can_tilt(self):
        out, _coords = CoordinateGenerator(seed=0).generate_3d_coordinates(
            _methyl_decorated_sheet()
        )
        owned = {a for members in _substituent_groups(out).values() for a in members}
        for c in out.GetAtoms():
            if c.GetAtomicNum() == 6 and not c.HasProp("init_x"):
                assert c.GetIdx() in owned
                assert all(n.GetIdx() in owned for n in c.GetNeighbors()
                           if n.GetAtomicNum() == 1)


def _hydroxylated_decorated_sheet():
    """The methyl-decorated sheet with every methyl turned into -CH2-OH and an
    OH on every free edge carbon: crowded enough that some of each kind of
    hydroxyl hydrogen has to be moved."""
    sheet = Chem.RWMol(_methyl_decorated_sheet())
    methyls = [a.GetIdx() for a in sheet.GetAtoms()
               if a.GetAtomicNum() == 6 and not a.HasProp("init_x")]
    edge = [a.GetIdx() for a in sheet.GetAtoms() if a.GetIsAromatic() and a.GetDegree() < 3]
    for site in methyls + edge:
        sheet.AddBond(site, sheet.AddAtom(Chem.Atom(8)), Chem.BondType.SINGLE)
    mol = sheet.GetMol()
    mol.UpdatePropertyCache(strict=False)
    return Chem.AddHs(mol, addCoords=False)


def _hydroxyl_angles(mol, coords):
    """{hydroxyl H index: (C-O-H angle, whether the carbon is aromatic)}."""
    angles = {}
    for o in mol.GetAtoms():
        hs = [n.GetIdx() for n in o.GetNeighbors() if n.GetAtomicNum() == 1]
        heavy = [n for n in o.GetNeighbors() if n.GetAtomicNum() != 1]
        if o.GetAtomicNum() == 8 and len(hs) == 1 and len(heavy) == 1:
            c = heavy[0]
            angles[hs[0]] = (
                _angle(coords[c.GetIdx()], coords[o.GetIdx()], coords[hs[0]]),
                c.GetIsAromatic(),
            )
    return angles


class TestCrowdedHydroxyls:
    # rq-7515b128 -- a crowded hydroxyl hydrogen is turned about its C-O bond,
    # not translated off it. Its out-of-plane candidates used to add a z offset
    # and rescale only the bond, which bent C-O-H anywhere from 95 to 178
    # degrees.
    # rq-3e250d24
    def test_crowded_hydroxyls_keep_their_angle(self, monkeypatch):
        mol = _hydroxylated_decorated_sheet()
        seen = {}
        original = geometry_3d._optimize_h_positions

        def spy(m, coords, **kwargs):
            seen["before"] = coords.copy()
            out = original(m, coords, **kwargs)
            seen["after"] = out.copy()
            return out

        monkeypatch.setattr(geometry_3d, "_optimize_h_positions", spy)
        gen = CoordinateGenerator(seed=0)
        out, coords = gen.generate_3d_coordinates(mol)
        assert gen.used_hex_lattice

        moved = {
            h for h in np.flatnonzero(
                np.linalg.norm(seen["after"] - seen["before"], axis=1) > 1e-6
            )
        }
        angles = _hydroxyl_angles(mol, coords)
        crowded = {h: angles[h] for h in moved if h in angles}
        # The premise: both kinds of hydroxyl were crowded enough to be moved.
        assert any(aromatic for _a, aromatic in crowded.values())
        assert any(not aromatic for _a, aromatic in crowded.values()), "no -CH2-OH was moved"
        assert all(100.0 <= a <= 120.0 for a, _ in angles.values()), sorted(angles.values())

    # rq-3e250d24
    @pytest.mark.slow
    def test_hydroxyls_on_a_150_carbon_sheet_are_not_linear(self, monkeypatch):
        # Before the fix, seven of this structure's fifteen hydroxyls came out
        # at 178 degrees: a z offset that cancelled the H's out-of-plane
        # component left it on the C-O axis.
        mol, coords, _comp = _generate_on_hex_lattice(
            monkeypatch, target_num_carbons=150, H_C_ratio=0.5, O_C_ratio=0.1,
            seed=5, strict=False,
        )
        angles = [a for a, _ in _hydroxyl_angles(mol, coords).values()]
        assert angles, "expected hydroxyl groups on this structure"
        assert all(100.0 <= a <= 120.0 for a in angles), sorted(angles)


def _fjord_hydroxyls():
    """[5]helicene on the flat lattice with an OH on each fjord carbon.

    The two fjord carbons are unbonded and 1.42 A apart, so each radially
    placed oxygen lands 0.06 A from the other arm's carbon -- what the
    elongated skeleton builder produced at 150 C, H/C 0.35, seed 4."""
    step = LATTICE_BOND * np.sqrt(3.0)
    centres, p = [np.zeros(2)], np.zeros(2)
    for turn in (0.0, 60.0, 120.0, 180.0):
        p = p + step * np.array([np.cos(np.radians(turn)), np.sin(np.radians(turn))])
        centres.append(p)
    mol, coords = _flat_pah_on_lattice(centres)
    carbons = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6]
    fjord = next(
        (i, j) for i, j in itertools.combinations(carbons, 2)
        if mol.GetBondBetweenAtoms(i, j) is None
        and np.linalg.norm(coords[i] - coords[j]) < 1.5
    )
    rw, coords = Chem.RWMol(mol), list(coords)
    for c in fjord:
        o = next(n.GetIdx() for n in rw.GetAtomWithIdx(c).GetNeighbors() if n.GetAtomicNum() == 1)
        rw.GetAtomWithIdx(o).SetAtomicNum(8)
        out = coords[o] - coords[c]
        coords[o] = coords[c] + 1.36 * out / np.linalg.norm(out)
    heavy = [x for k, x in enumerate(coords) if rw.GetAtomWithIdx(k).GetAtomicNum() != 1]
    for c in fjord:
        o = next(n.GetIdx() for n in rw.GetAtomWithIdx(c).GetNeighbors() if n.GetAtomicNum() == 8)
        h = rw.AddAtom(Chem.Atom(1))
        rw.AddBond(o, h, Chem.BondType.SINGLE)
        coords.append(_bent_hydrogen_position(coords[o], coords[c], 0.96, heavy))
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    return mol, np.array(coords)


class TestReliefDoesNotDependOnHydroxylHydrogens:
    # rq-8e5f5dc1 -- relief must not hinge on where hydroxyl hydrogens point.
    # rq-340c44e4 -- seed 4 at O/C 0.3 is the point the first attempt at
    # rq-3e250d24 broke (0 -> 115 clashes): two unbonded fjord carbons 1.42 A
    # apart each carry an OH aimed at the other, and the tilt turned both the
    # same way. Seed 4 at O/C 0.1 failed the same way before the change.
    @pytest.mark.slow
    @pytest.mark.parametrize("o_c", [0.1, 0.3])
    def test_fjord_hydroxyls_are_relieved(self, monkeypatch, o_c):
        mol, coords, _comp = _generate_on_hex_lattice(
            monkeypatch, target_num_carbons=150, H_C_ratio=0.35, O_C_ratio=o_c,
            seed=4, strict=False,
        )
        assert _clash_pairs(mol, coords) == []

    # rq-7cab77bd
    def test_groups_sharing_a_tilt_direction_still_separate(self):
        # Two hydroxyls on the arms of a fjord, aimed at each other: each O
        # sits almost on the other arm's carbon. Their first contacts are with
        # those carbons, so each group is handed a direction on its own --
        # both +1 -- before the O...O contact between them is reached.
        mol, coords = _fjord_hydroxyls()
        tilted = _tilt_crowded_substituents(mol, coords)
        o1, o2 = (a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 8)
        assert np.linalg.norm(tilted[o1] - tilted[o2]) > 2.0
        assert np.sign(tilted[o1][2]) != np.sign(tilted[o2][2])


class TestStrictValidationOnDecoratedSheets:
    # rq-d0935f4f -- issue #66: these seeds reach H/C 0.35 at 150 C through
    # methyl decoration, and used to fall off the lattice onto Kamada-Kawai.
    @pytest.mark.parametrize("seed", [2, 3, 5])
    def test_no_steric_clash_failures_at_150_carbons(self, monkeypatch, seed):
        try:
            _generate_on_hex_lattice(
                monkeypatch, target_num_carbons=150, H_C_ratio=0.35,
                O_C_ratio=0.0, seed=seed,
            )
        except ValidationError as e:
            assert "Steric clash" not in str(e), f"seed {seed}: {str(e)[:300]}"


def _spy_embeddings(monkeypatch):
    """Record the use_hex_lattice flag of every embedding the generator asks for."""
    flags = []
    original = CoordinateGenerator.generate_3d_coordinates

    def spy(self, mol, *args, **kwargs):
        flags.append(kwargs.get("use_hex_lattice", True))
        return original(self, mol, *args, **kwargs)

    monkeypatch.setattr(CoordinateGenerator, "generate_3d_coordinates", spy)
    return flags


class TestDecoratedSheetFallback:
    # rq-ab5203c4
    @pytest.mark.slow
    def test_an_unclearable_decorated_sheet_keeps_the_better_embedding(self, monkeypatch):
        # 150 C at H/C 0.6 is heavily methylated; the relief leaves clashes
        # at doubly decorated coves.
        flags = _spy_embeddings(monkeypatch)
        relieved, free = {}, {}
        relieve = CoordinateGenerator.relieve_lattice_crowding
        refine = BiocharGenerator._refine_free_embedding

        def relieve_spy(self, mol, coords):
            out = relieve(self, mol, coords)
            relieved["errors"] = GeometryValidator.validate_geometry(mol, out)[1]
            return out

        def refine_spy(self, generator, mol, coords):
            out = refine(self, generator, mol, coords)
            free["errors"] = out[2]
            return out

        monkeypatch.setattr(CoordinateGenerator, "relieve_lattice_crowding", relieve_spy)
        monkeypatch.setattr(BiocharGenerator, "_refine_free_embedding", refine_spy)
        mol, coords, _ = BiocharGenerator(GeneratorConfig(
            target_num_carbons=150, H_C_ratio=0.6, O_C_ratio=0.0, seed=4, strict=False,
        )).generate()

        assert flags == [True, False]
        assert any("Steric clash" in e for e in relieved["errors"])
        final = GeometryValidator.validate_geometry(mol, coords)[1]
        assert len(final) == min(len(relieved["errors"]), len(free["errors"]))

    # rq-0cacf78d
    def test_an_undecorated_sheet_is_embedded_once(self, monkeypatch):
        # Relief is stubbed out so clashes are certain to remain: whether the
        # real relief clears a given seed differs between RDKit/NumPy builds.
        flags = _spy_embeddings(monkeypatch)
        seen = {}

        def no_relief(self, mol, coords):
            seen["offlattice"] = self.offlattice_carbons
            seen["clashes"] = len(_clash_pairs(mol, coords))
            return coords

        monkeypatch.setattr(CoordinateGenerator, "relieve_lattice_crowding", no_relief)
        BiocharGenerator(GeneratorConfig(
            target_num_carbons=150, H_C_ratio=0.35, O_C_ratio=0.0, seed=4, strict=False,
        )).generate()
        # The premise, or this test passes without reaching the gated branch.
        assert seen["offlattice"] == 0 and seen["clashes"] > 0
        assert flags == [True]
