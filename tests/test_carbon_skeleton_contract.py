"""
Tests for the rqm/carbon-skeleton.md scenarios -- what a caller asks the
assembler for, and how faithfully the request is answered.

Distinct from the size and quality checks elsewhere: these are about fidelity.
A skeleton is assembled from whole rings, so several requests can only be
approximated, and the point is which approximations are honest. A parameter
accepted and ignored, a count met by a different number, or a fallback returning
a fixed 16-carbon structure are all the same failure wearing different clothes.

Sizes are read from PAH_LIBRARY rather than hard-coded, so extending the library
moves the tests with it.

Every scenario here passes. The three that carried xfail(strict=True) were
retired by XPASS when their fixes landed.
"""

import inspect

import pytest

from biochar.constants import PAH_LIBRARY
from biochar.pipeline.carbon_skeleton import PAHAssembler


def _library_sizes():
    return sorted({int(v["num_atoms"]) for v in PAH_LIBRARY.values()})


def _generate(target, **kw):
    return PAHAssembler(seed=1).generate(target_num_carbons=target, **kw)


def _carbons(skeleton):
    return sum(1 for a in skeleton.mol.GetAtoms() if a.GetAtomicNum() == 6)


def _ring_sizes(skeleton):
    return [len(r) for r in skeleton.mol.GetRingInfo().AtomRings()]


class TestLibraryTargets:
    # rq-4af71330
    @pytest.mark.parametrize("target", _library_sizes())
    def test_a_library_size_is_met_exactly(self, target):
        assert _carbons(_generate(target)) == target


class TestGrowthReachesFromAbove:
    # Targets chosen to sit between reachable ring counts.
    AWKWARD = [41, 43, 45, 47, 49, 53, 57, 61, 71, 83, 97]

    # rq-59e23c2e
    @pytest.mark.parametrize("target", AWKWARD)
    def test_a_grown_skeleton_is_never_smaller_than_requested(self, target):
        actual = _carbons(_generate(target))
        assert actual >= target, (
            f"requested {target} carbons and got {actual}; a skeleton short of "
            f"the request has fewer edge sites than the caller planned for"
        )

    # rq-5435c9d8
    def test_the_overshoot_stays_within_a_ring_or_two(self):
        deltas = {t: _carbons(_generate(t)) - t for t in self.AWKWARD}
        worst = max(deltas.values())
        assert worst <= 8, (
            f"overshoot of {worst} carbons is more than a ring or two: {deltas}"
        )


class TestOddSheetsDoNotRunAway:
    # rq-e8417061
    # (target, seed, target_h_c): each grew to 228, 936 and 1680 carbons before
    # the parity step was made to fuse only where the addition is odd.
    RUNAWAYS = [(150, 3, 0.5), (400, 6, None), (628, 5, None)]

    @pytest.mark.parametrize("target,seed,target_h_c", RUNAWAYS)
    def test_an_odd_sheet_stays_within_a_ring_or_two(self, target, seed, target_h_c):
        skeleton = PAHAssembler(seed=seed).generate(target, target_h_c=target_h_c)
        over = _carbons(skeleton) - target
        assert 0 <= over <= 8, (
            f"requested {target} carbons, got {_carbons(skeleton)}; the parity "
            f"step kept adding hexagons, which cannot change parity"
        )


class TestElongatedGrowth:
    """Issue #65: elongated growth returned sheets up to 55% oversized that had
    wrapped onto themselves, leaving unbonded lattice carbons at bond distance.
    On main this grid failed 18 of 90 cases."""

    GRID = [(n, h) for n in (100, 150, 200, 250) for h in (0.35, 0.45, 0.6)]

    # rq-5e74710c
    @pytest.mark.parametrize("target,target_h_c", GRID)
    def test_an_elongated_skeleton_stays_within_a_ring_or_two(self, target, target_h_c):
        for seed in range(1, 7):
            skeleton = PAHAssembler(seed=seed).generate(target, target_h_c=target_h_c)
            over = _carbons(skeleton) - target
            assert 0 <= over <= 8, (
                f"seed {seed}: requested {target} carbons at H/C {target_h_c}, "
                f"got {_carbons(skeleton)}"
            )

    # rq-20d8d5ea
    @pytest.mark.parametrize("target,target_h_c", GRID + [(200, None), (250, None)])
    def test_no_lattice_carbons_at_bond_distance_without_a_bond(self, target, target_h_c):
        import numpy as np

        for seed in range(1, 7):
            mol = PAHAssembler(seed=seed).generate(target, target_h_c=target_h_c).mol
            xy = np.array([[a.GetDoubleProp("init_x"), a.GetDoubleProp("init_y")]
                           for a in mol.GetAtoms()])
            d = np.linalg.norm(xy[:, None, :] - xy[None, :, :], axis=-1)
            bonded = {frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
                      for b in mol.GetBonds()}
            close = [(i, j) for i, j in zip(*np.nonzero(np.triu(d < 1.6, k=1)))
                     if frozenset((int(i), int(j))) not in bonded]
            assert not close, (
                f"seed {seed}, {target} C, H/C {target_h_c}: {len(close)} unbonded "
                f"carbon pairs under 1.6 Å -- the sheet wrapped onto itself"
            )


class TestSeedDependentFailuresAreRetried:
    # rq-a3c395ed
    def test_a_growth_that_fails_on_the_callers_seed_is_retried(self):
        # Seed 2 at 104 carbons closed the compact sheet odd on its own seed
        # and raised "retrying will not help"; seeds 1 and 3 built fine.
        skeleton = PAHAssembler(seed=2).generate(104, target_h_c=0.635)
        assert _carbons(skeleton) >= 104

    # rq-7504e809
    def test_neighboring_seeds_do_not_share_retries(self):
        from biochar.pipeline.carbon_skeleton import _MAX_CURVED_BUILD_ATTEMPTS, _sub_seed

        n = _MAX_CURVED_BUILD_ATTEMPTS
        first = {_sub_seed(1, k) for k in range(n)}
        second = {_sub_seed(2, k) for k in range(n)}
        assert not first & second, (
            f"seeds 1 and 2 retry with shared sub-seeds {first & second}, so a "
            f"different seed replays the same growth"
        )


WOOD_FRACTIONS = dict(defect_fraction=0.154, heptagon_fraction=0.077)


class TestCurvedSheets:
    # rq-22536e48
    @pytest.mark.parametrize("seed", [1, 2, 3])
    def test_no_carbon_has_more_than_three_bonds(self, seed):
        skeleton = PAHAssembler(seed=seed).generate(256, **WOOD_FRACTIONS)
        worst = max(a.GetDegree() for a in skeleton.mol.GetAtoms())
        assert worst <= 3, f"a carbon has {worst} carbon neighbors"

    # rq-8f5204cd
    @pytest.mark.parametrize("target", [102, 256, 608])
    def test_wood_sized_curved_sheets_build_on_every_seed(self, target):
        for seed in (1, 2, 3):
            skeleton = PAHAssembler(seed=seed).generate(target, **WOOD_FRACTIONS)
            rings = skeleton.ring_composition
            assert _carbons(skeleton) >= target
            assert rings["pentagons"] > 0 and rings["heptagons"] > 0, (
                f"seed {seed}: a curved request returned rings {rings}"
            )

    # rq-f1e2eb6a
    def test_aromaticity_perception_on_a_large_curved_sheet_is_quick(self):
        """Run in a subprocess: on a bare-bond sheet RDKit does not fail, it
        runs for minutes, and a test that hangs reports nothing useful.

        The Kekule form comes from a matching, not ``Chem.Kekulize``: RDKit's
        own kekulization gives up on any sheet this size, flat ones included.
        """
        import subprocess
        import sys

        code = (
            "import time\n"
            "import networkx as nx\n"
            "from rdkit import Chem\n"
            "from biochar.pipeline.carbon_skeleton import PAHAssembler, _mol_to_graph\n"
            "sk = PAHAssembler(seed=1).generate(400, defect_fraction=0.154,"
            " heptagon_fraction=0.077).mol\n"
            "G = _mol_to_graph(sk)\n"
            "double = {frozenset(e) for e in nx.max_weight_matching(G, maxcardinality=True)}\n"
            "rw = Chem.RWMol()\n"
            "for _ in range(sk.GetNumAtoms()): rw.AddAtom(Chem.Atom(6))\n"
            "for b in sk.GetBonds():\n"
            "    i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()\n"
            "    rw.AddBond(i, j, Chem.BondType.DOUBLE if frozenset((i, j)) in double"
            " else Chem.BondType.SINGLE)\n"
            "m = rw.GetMol()\n"
            "f = Chem.SanitizeFlags\n"
            "Chem.SanitizeMol(m, f.SANITIZE_ALL ^ f.SANITIZE_SETAROMATICITY)\n"
            "t = time.time()\n"
            "Chem.SetAromaticity(m, Chem.AromaticityModel.AROMATICITY_RDKIT)\n"
            "bare = sum(1 for b in m.GetBonds() if b.IsInRing() and not b.GetIsAromatic())\n"
            "print(time.time() - t, bare)\n"
        )
        try:
            out = subprocess.run(
                [sys.executable, "-c", code], capture_output=True, text=True, timeout=60
            )
        except subprocess.TimeoutExpired as exc:
            raise AssertionError(
                "RDKit aromaticity perception on a 400-C curved sheet ran past 60 s"
            ) from exc
        assert out.returncode == 0, out.stderr[-2000:]
        elapsed, bare = out.stdout.split()
        assert int(bare) == 0, f"{bare} ring bonds were left non-aromatic"
        assert float(elapsed) < 1.0, f"perception took {float(elapsed):.1f} s"

    # rq-5af71ee2
    def test_a_large_pentagon_only_request_names_the_missing_heptagons(self):
        from biochar.pipeline.carbon_skeleton import SkeletonError

        with pytest.raises(SkeletonError) as exc:
            PAHAssembler(seed=1).generate(300, defect_fraction=0.1)
        assert "heptagon" in str(exc.value), str(exc.value)


class TestDefectFractions:
    # rq-c6b675bb
    def test_zero_defect_fraction_gives_a_pure_hexagonal_skeleton(self):
        sizes = _ring_sizes(_generate(80, defect_fraction=0.0, heptagon_fraction=0.0))
        assert sizes, "skeleton has no rings"
        assert set(sizes) == {6}, f"expected only hexagons, got sizes {sorted(set(sizes))}"

    # rq-0c91f5d5
    def test_a_positive_defect_fraction_produces_pentagons_below_the_probability(self):
        requested = 0.3
        sizes = _ring_sizes(_generate(80, defect_fraction=requested))
        pentagons = sum(1 for s in sizes if s == 5)

        assert pentagons > 0, "a positive defect fraction produced no pentagons"
        realised = pentagons / len(sizes)
        # Only rings added during growth are subject to the draw; the seed is
        # all hexagons, so the realised share sits below the probability.
        assert realised <= requested, (
            f"realised pentagon share {realised:.3f} exceeds the requested "
            f"per-ring probability {requested}"
        )


class TestAromaticityIsNotARequest:
    # rq-b6eec273
    def test_the_assembler_takes_no_aromaticity_argument(self):
        params = inspect.signature(PAHAssembler.generate).parameters
        aromatic = [p for p in params if "aromatic" in p.lower()]
        assert not aromatic, (
            f"generate() accepts {aromatic}, which it cannot honour; aromaticity "
            f"is decided by the ring topology"
        )


class TestNoUnrelatedFallback:
    # rq-f953b9b0
    def test_a_skeleton_that_cannot_be_grown_raises(self, monkeypatch):
        from biochar.pipeline.carbon_skeleton import SkeletonError

        monkeypatch.setattr(PAHAssembler, "_build_from_seed", lambda *a, **k: None)

        with pytest.raises(SkeletonError) as exc:
            _generate(200)

        message = str(exc.value)
        assert "200" in message, f"the failure does not state the request: {message}"
        # Growth is seed-dependent (rq-4b588285); the message once claimed the
        # opposite while seeds 1 and 3 succeeded where 2 failed.
        assert "another seed may succeed" in message, (
            f"the failure should say a different seed can succeed: {message}"
        )

    def test_the_failure_is_catchable_as_a_runtime_error(self, monkeypatch):
        """Subclassing RuntimeError keeps existing handlers working.

        workflows.sweep catches bare Exception for non-validation failures and
        records the point as failed.
        """
        from biochar.pipeline.carbon_skeleton import SkeletonError

        assert issubclass(SkeletonError, RuntimeError)
        monkeypatch.setattr(PAHAssembler, "_build_from_seed", lambda *a, **k: None)
        with pytest.raises(RuntimeError):
            _generate(200)


class TestLibraryLoadIsReported:
    # rq-4dc5ce35
    def test_an_unparseable_library_entry_is_reported(self, monkeypatch, caplog):
        import logging

        broken = dict(PAH_LIBRARY)
        broken["not_a_molecule"] = {
            "smiles": "this is not smiles",
            "num_atoms": 99,
            "num_aromatic": 99,
            "molecular_formula": "C99",
            "references": "fixture",
        }
        monkeypatch.setattr(
            "biochar.pipeline.carbon_skeleton.PAH_LIBRARY", broken, raising=False
        )
        monkeypatch.setattr(PAHAssembler, "_SIZE_INDEX", None)

        with caplog.at_level(logging.WARNING, logger="biochar"):
            assembler = PAHAssembler(seed=1)

        assert "not_a_molecule" not in assembler.pahs, "fixture entry unexpectedly loaded"
        assert caplog.records, (
            "a library entry was dropped from the working set with nothing "
            "reported; PAH_LIBRARY claims every entry is validated"
        )
