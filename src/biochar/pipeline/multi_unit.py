"""
Linked multi-unit ("island") structures.

Low-temperature biochar is not one fused sheet. Wood, Mašek & Erastova 2024
(Cell Rep. Phys. Sci., DOI 10.1016/j.xcrp.2024.102036) build their BC400 model
from a block of four separate aromatic ring systems, BC600 from three, and only
BC800 from a single sheet; the smaller units hang off a core through biaryl
bonds, ether bridges and short alkyl chains. This module builds that
architecture: several aromatic units, each from :class:`PAHAssembler`, joined
in a tree by linkers, and a 3D embedding that keeps every unit flat while the
links rotate.

It is opt-in (``GeneratorConfig.linked_units``) and runs beside the fused-sheet
path rather than inside it; :class:`BiocharGenerator` calls it through three
small hooks (skeleton, oxygen target, geometry). Specified in
``rqm/linked-units.md``.

Unit and linker membership travel on atom properties -- ``unit_id`` on every
atom of an aromatic unit, ``linker_id`` on every heavy atom of a linker -- set
here at skeleton time. Atom properties survive the oxygen, nitrogen and
hydrogen stages (the hex-lattice ``init_x``/``init_y`` already rely on that),
while atom indices do not: hydrogen trimming removes atoms.
"""

from __future__ import annotations

import logging
import math
import random
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from rdkit import Chem

from ..constants import LATTICE_MAX_SHEET_RMS
from .carbon_skeleton import CarbonSkeleton, PAHAssembler, SkeletonError
from .geometry_3d import (
    GeometryValidator,
    _clash_floor_matrix,
    _clash_pairs,
    _get_excluded_pairs,
    _tethered_lattice_relax,
    _tilt_crowded_substituents,
)
from .heteroatom_assignment import (
    _fix_heteroatom_bond_types,
    _safe_sanitize,
    attach_aliphatic_carbons,
)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Linkers
# ---------------------------------------------------------------------------

#: Every linker this module builds, with the atoms each puts between the two
#: aromatic carbons it joins. Wood's low-temperature blocks also link units
#: through diaryl ketones (0004, 0005) and glycerol-type beta-O-4 arms (0045);
#: the first needs a ketone carbonyl type the typer does not have yet, and the
#: second is side-chain chemistry (#68, #69).
LINKER_BRIDGES: Dict[str, Tuple[int, ...]] = {
    "biaryl": (),                    # Ar-Ar
    "diaryl_ether": (8,),            # Ar-O-Ar
    "aryl_alkyl_ether": (8, 6),      # Ar-O-CH2-Ar
    "propyl": (6, 6, 6),             # Ar-CH2CH2CH2-Ar
}
LINKER_TYPES: Tuple[str, ...] = tuple(LINKER_BRIDGES)


def linker_carbons(name: str) -> int:
    return sum(1 for z in LINKER_BRIDGES[name] if z == 6)


def linker_oxygens(name: str) -> int:
    return sum(1 for z in LINKER_BRIDGES[name] if z == 8)


# Temperature bands. A design choice anchored on Wood's six final models -- the
# only temperature-labelled structures in their library: BC400/BC400c (blocks
# 0045/0046) have 4 ring systems, BC600/BC600c (0052/0053) have 3, BC800/BC800c
# have 1. Edges sit at the midpoints. Not a fit: the 57-block library carries
# no temperature labels, and two-unit blocks (0017-0019, 0024, 0043, 0049,
# 0051) occur at every composition, so k=2 is left to an explicit num_units.
LINKED_LOW_T_MAX = 500.0     # below: 4-5 units, low-T linker mix
LINKED_MID_T_MAX = 700.0     # below: 3 units, mid-T mix; at or above: 1 unit
LINKED_LOW_T_UNITS = (4, 5)
LINKED_MID_T_UNITS = 3
LINKED_DEFAULT_UNITS = 4     # linked_units=True with no temperature

# Low-T mix: the inter-unit links of Wood blocks 0001-0008 this module can build
# (30 links: 15 diaryl ether, 9 biaryl, 4 propyl, 2 diaryl ketone; ketones
# deferred, so weights are over the remaining 28).
LINKER_MIX_LOW_T: Dict[str, float] = {
    "diaryl_ether": 15 / 28,
    "biaryl": 9 / 28,
    "propyl": 4 / 28,
}
# Mid-T mix: a judgment from 0009 (2 propyl), 0013-0016 (2 diaryl ether each)
# and 0052/0053 (diaryl ether plus Ar-C-C-O-Ar beta-O-4-type links), not a
# count. aryl_alkyl_ether stands in for the beta-O-4 links until #69.
LINKER_MIX_MID_T: Dict[str, float] = {
    "diaryl_ether": 0.65,
    "aryl_alkyl_ether": 0.20,
    "propyl": 0.15,
}

# Smallest unit: benzene.
MIN_UNIT_CARBONS = 6
# Arm sizes PAHAssembler builds exactly (its library); only the core is grown
# to an arbitrary count, so it is built last and absorbs the remainder.
_ARM_SIZES = (6, 10, 14, 16, 18, 22, 24, 30, 40)
# Core share of the aromatic carbon budget tried when sizing for H/C.
_CORE_FRACTIONS = (0.3, 0.45, 0.6, 0.75)
# Probability a new unit hangs off the core rather than an earlier arm. Wood's
# blocks show both shapes: arms on a core (0045, 0052) and chains (0001, 0007).
_CORE_PARENT_BIAS = 0.6


@dataclass
class UnitPlan:
    """How many units to build and which linker joins each non-core unit."""

    num_units: int
    linkers: List[str]


@dataclass
class LinkedSkeleton:
    """Census of a linked structure, reported on the generator and the result."""

    unit_carbons: List[int]
    linker_counts: Dict[str, int]
    num_methyls: int = 0
    parents: List[int] = field(default_factory=list)

    @property
    def num_units(self) -> int:
        return len(self.unit_carbons)

    @property
    def num_linker_oxygens(self) -> int:
        return sum(linker_oxygens(t) * n for t, n in self.linker_counts.items())

    def describe(self) -> str:
        """One line for file headers, e.g. ``4 units (64, 10, 6, 6 C); links: ...``."""
        sizes = ", ".join(str(n) for n in self.unit_carbons)
        links = ", ".join(
            f"{n} {t.replace('_', ' ')}" for t, n in sorted(self.linker_counts.items()) if n
        )
        text = f"{self.num_units} linked aromatic units ({sizes} C); links: {links}"
        if self.num_methyls:
            text += f"; {self.num_methyls} pendant methyls"
        return text


# ---------------------------------------------------------------------------
# Planning
# ---------------------------------------------------------------------------

def _band(temperature: Optional[float]) -> str:
    if temperature is None:
        return "default"
    if temperature < LINKED_LOW_T_MAX:
        return "low"
    if temperature < LINKED_MID_T_MAX:
        return "mid"
    return "high"


def _hc_floor(num_units: int, total_carbons: int, linkers: Sequence[str]) -> float:
    """
    The lowest saturated H/C *num_units* units can reach at *total_carbons*.

    Splitting a sheet adds perimeter, so H/C has a floor that rises with the
    unit count. The least-hydrogen split is a compact core plus benzene arms: a
    compact PAH of n carbons has about sqrt(6n) edge hydrogens (coronene 24 ->
    12, circumcoronene 54 -> 18), each arm carries 5 (one site taken by its
    link), each link costs the core or parent one more, and every linker CH2
    carries two.
    """
    link_c = sum(linker_carbons(t) for t in linkers)
    arms = num_units - 1
    core = total_carbons - link_c - MIN_UNIT_CARBONS * arms
    if core < MIN_UNIT_CARBONS:
        return float("inf")
    hydrogens = math.sqrt(6.0 * core) - arms + 5 * arms + 2 * link_c
    return hydrogens / total_carbons


def plan_units(config, rng: random.Random) -> UnitPlan:
    """
    Resolve the unit count and linker sequence for an opted-in *config*.

    An explicit ``num_units`` is honoured (with a warning when its H/C floor is
    above the target); a temperature-derived count shrinks, with a warning,
    until it fits the carbon budget and can reach the target H/C -- the caller
    asked for a temperature, not for that count.
    """
    band = _band(config.temperature)
    if config.num_units is not None:
        k = config.num_units
    elif band == "low":
        k = rng.choice(LINKED_LOW_T_UNITS)
    elif band == "mid":
        k = LINKED_MID_T_UNITS
    elif band == "high":
        k = 1
    else:
        k = LINKED_DEFAULT_UNITS

    if config.linker_mix is not None:
        mix = config.linker_mix
    else:
        mix = LINKER_MIX_MID_T if band == "mid" else LINKER_MIX_LOW_T
    names = [t for t in LINKER_TYPES if mix.get(t, 0.0) > 0]
    weights = [mix[t] for t in names]
    linkers = rng.choices(names, weights=weights, k=max(k - 1, 0))

    N = config.target_num_carbons
    explicit = config.num_units is not None
    hc_target = None if config.functional_groups else config.H_C_ratio

    def fits(kk: int, lk: Sequence[str]) -> bool:
        if N - sum(linker_carbons(t) for t in lk) < MIN_UNIT_CARBONS * kk:
            return False
        if hc_target is None or kk == 1:
            return True
        return _hc_floor(kk, N, lk) <= hc_target * (1.0 + config.H_C_tolerance)

    if not explicit:
        k0 = k
        while k > 1 and not fits(k, linkers[: k - 1]):
            k -= 1
        linkers = linkers[: k - 1]
        if k < k0:
            logger.warning(
                "Linked units: %d units cannot reach H/C %.3f at %d carbons (each "
                "split adds perimeter hydrogen); building %d instead.",
                k0, hc_target if hc_target is not None else float("nan"), N, k,
            )
    else:
        # An explicit count the carbon budget cannot hold even with
        # carbon-free linkers is refused at config time; here the drawn
        # linkers may still be too carbon-hungry, so trade them down.
        while N - sum(linker_carbons(t) for t in linkers) < MIN_UNIT_CARBONS * k:
            heavy = max(range(len(linkers)), key=lambda i: linker_carbons(linkers[i]))
            linkers[heavy] = "biaryl"
        if hc_target is not None and k > 1 and _hc_floor(k, N, linkers) > hc_target * (
            1.0 + config.H_C_tolerance
        ):
            logger.warning(
                "Linked units: num_units=%d at %d carbons has an H/C floor of about "
                "%.3f, above the target %.3f; validation will report the miss. "
                "Fewer units or more carbons would reach it.",
                k, N, _hc_floor(k, N, linkers), hc_target,
            )
    return UnitPlan(num_units=k, linkers=list(linkers))


# ---------------------------------------------------------------------------
# Skeleton
# ---------------------------------------------------------------------------

def _edge_sites(mol: Chem.Mol, atoms: Sequence[int], used: set) -> List[int]:
    """Aromatic carbons of a unit that can take a link, least crowded first.

    A site whose ring neighbours both carry hydrogen sits on an open edge; one
    next to a fusion carbon sits in a bay, where an attached unit is crowded by
    the peri hydrogen. Sites bonded to an existing attachment are skipped.
    """
    blocked = set(used)
    for u in used:
        blocked.update(n.GetIdx() for n in mol.GetAtomWithIdx(u).GetNeighbors())
    open_edge, bay = [], []
    for i in atoms:
        a = mol.GetAtomWithIdx(i)
        if a.GetAtomicNum() != 6 or not a.GetIsAromatic() or a.GetDegree() >= 3 or i in blocked:
            continue
        if all(n.GetDegree() <= 2 for n in a.GetNeighbors()):
            open_edge.append(i)
        else:
            bay.append(i)
    return open_edge + bay


def _saturated_hc(mol: Chem.Mol) -> float:
    """H/C once every free valence carries hydrogen (explicit H counted)."""
    carbons = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 6]
    h = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 1)
    h += sum(1 for a in carbons if a.GetIsAromatic() and a.GetDegree() < 3)
    return h / len(carbons) if carbons else 0.0


def _arm_sizes(arm_budget: int, num_arms: int, core_size: int) -> List[int]:
    """Split *arm_budget* into *num_arms* library sizes, none above the core."""
    each = arm_budget / num_arms
    allowed = [s for s in _ARM_SIZES if s <= max(core_size, MIN_UNIT_CARBONS)]
    return [min(allowed, key=lambda s: abs(s - each))] * num_arms


class _Assembly:
    """One candidate linked skeleton under construction."""

    def __init__(self):
        self.rw = Chem.RWMol()
        self.units: List[List[int]] = []

    def add_unit(self, unit: Chem.Mol) -> int:
        uid = len(self.units)
        index = {}
        for a in unit.GetAtoms():
            new = Chem.Atom(a.GetAtomicNum())
            new.SetIsAromatic(a.GetIsAromatic())
            idx = self.rw.AddAtom(new)
            atom = self.rw.GetAtomWithIdx(idx)
            for prop in ("init_x", "init_y"):
                if a.HasProp(prop):
                    atom.SetDoubleProp(prop, a.GetDoubleProp(prop))
            atom.SetIntProp("unit_id", uid)
            index[a.GetIdx()] = idx
        for b in unit.GetBonds():
            self.rw.AddBond(index[b.GetBeginAtomIdx()], index[b.GetEndAtomIdx()], b.GetBondType())
            nb = self.rw.GetBondBetweenAtoms(index[b.GetBeginAtomIdx()], index[b.GetEndAtomIdx()])
            nb.SetIsAromatic(b.GetIsAromatic())
        self.units.append(list(index.values()))
        return uid

    def link(self, linker: str, a: int, b: int, linker_id: int) -> None:
        """Join aromatic carbons *a* and *b* through *linker*; saturate its carbons."""
        chain = [a]
        for z in LINKER_BRIDGES[linker]:
            idx = self.rw.AddAtom(Chem.Atom(z))
            self.rw.GetAtomWithIdx(idx).SetIntProp("linker_id", linker_id)
            chain.append(idx)
        chain.append(b)
        for u, v in zip(chain, chain[1:]):
            self.rw.AddBond(u, v, Chem.BondType.SINGLE)
        # A linker carbon with a free valence is a site for the oxygen stage's
        # aliphatic spill, which would turn Ar-O-CH2-Ar into a hemiacetal and a
        # propyl bridge into an oxygenated arm the census does not report.
        # Explicit hydrogens leave it none.
        for idx in chain[1:-1]:
            if self.rw.GetAtomWithIdx(idx).GetAtomicNum() == 6:
                for _ in range(2):
                    h = self.rw.AddAtom(Chem.Atom(1))
                    self.rw.AddBond(idx, h, Chem.BondType.SINGLE)

    def mol(self, sanitize: bool = True) -> Chem.Mol:
        out = self.rw.GetMol()
        if not sanitize:
            # Enough for degree-based counting while candidates are compared;
            # sanitising every candidate costs more than building it.
            out.UpdatePropertyCache(strict=False)
            return out
        _safe_sanitize(out)
        return _fix_heteroatom_bond_types(out)


class _UnitCache:
    """Units built once per (size, role) for all candidates of one request."""

    def __init__(self, config, seed_base: int):
        self.config = config
        self.seed_base = seed_base
        self._units: Dict[Tuple[int, bool], Chem.Mol] = {}

    def get(self, size: int, core: bool) -> Chem.Mol:
        key = (size, core)
        if key not in self._units:
            # Arms are pure-hex library PAHs; ring defects apply to the core
            # only -- Wood's arms are benzene to pyrene size, and small grown
            # units with pentagons are unrepresentative and fail growth more.
            self._units[key] = PAHAssembler(seed=self.seed_base).generate(
                size,
                defect_fraction=self.config.defect_fraction if core else 0.0,
                heptagon_fraction=self.config.heptagon_fraction if core else 0.0,
            ).mol
        return self._units[key]


def _build_candidate(
    plan: UnitPlan, core_n: int, arm_sizes: Sequence[int], units: _UnitCache,
    sanitize: bool = True,
) -> Tuple[Chem.Mol, List[int]]:
    """Build and join the units of one candidate; return (mol, parents)."""
    rng = random.Random(units.seed_base)
    asm = _Assembly()
    asm.add_unit(units.get(core_n, core=True))
    for size in arm_sizes:
        asm.add_unit(units.get(size, core=False))

    used: set = set()
    parents = [-1]
    for i, linker in enumerate(plan.linkers, start=1):
        if i == 1 or rng.random() < _CORE_PARENT_BIAS:
            parent = 0
        else:
            parent = rng.randrange(1, i)
        tmp = asm.rw.GetMol()
        tmp.UpdatePropertyCache(strict=False)
        p_sites = _edge_sites(tmp, asm.units[parent], used)
        c_sites = _edge_sites(tmp, asm.units[i], used)
        if not p_sites:
            # A small parent can run out of sites; the core rarely does.
            parent = 0
            p_sites = _edge_sites(tmp, asm.units[0], used)
        if not p_sites or not c_sites:
            raise RuntimeError("no free edge site to attach a linked unit")
        # Choose among the least crowded few so layouts vary with the seed.
        a = rng.choice(p_sites[: max(1, len(p_sites) // 2)])
        b = rng.choice(c_sites[: max(1, len(c_sites) // 2)])
        asm.link(linker, a, b, i)
        used.update((a, b))
        parents.append(parent)
    return asm.mol(sanitize), parents


def build_linked_skeleton(
    config, plan: UnitPlan, aliphatic_enabled: bool
) -> Tuple[CarbonSkeleton, LinkedSkeleton]:
    """
    Build the linked carbon/oxygen skeleton for *plan*.

    Linker carbons are drawn first; the aromatic budget is split into a core and
    library-size arms over several core fractions, and the candidate whose
    saturated H/C is closest to the target wins. When even that falls short by
    more than 0.02, pendant methyls make up the difference exactly as the
    fused path does, with the core shrunk so the total carbon count holds.
    """
    N = config.target_num_carbons
    k = plan.num_units
    budget = N - sum(linker_carbons(t) for t in plan.linkers)
    target = None if config.functional_groups else config.H_C_ratio
    seed_base = 0 if config.seed is None else int(config.seed)

    candidates = []
    seen = set()
    for frac in _CORE_FRACTIONS:
        arms = _arm_sizes(budget * (1.0 - frac), k - 1, int(budget * frac))
        core_n = budget - sum(arms)
        if core_n < max(arms) or core_n < MIN_UNIT_CARBONS:
            continue
        key = (core_n, tuple(arms))
        if key in seen:
            continue
        seen.add(key)
        candidates.append((core_n, arms))
    if not candidates:
        arms = [MIN_UNIT_CARBONS] * (k - 1)
        candidates.append((budget - sum(arms), arms))

    units = _UnitCache(config, seed_base)
    best = None
    failures = []
    for core_n, arms in candidates:
        try:
            mol, _ = _build_candidate(plan, core_n, arms, units, sanitize=False)
        except SkeletonError as exc:
            # Growth fails at some sizes regardless of seed; another core
            # fraction is another size.
            failures.append(str(exc))
            continue
        hc = _saturated_hc(mol)
        score = abs(hc - target) if target is not None else abs(core_n / budget - 0.5)
        if best is None or score < best[0]:
            best = (score, hc, core_n, arms)
    if best is None:
        raise SkeletonError(
            f"could not build a core for any linked-unit split of "
            f"{N} carbons into {k} units: " + "; ".join(failures[:2])
        )
    _, hc, core_n, arms = best
    mol, parents = _build_candidate(plan, core_n, arms, units)

    methyls = 0
    if target is not None and aliphatic_enabled and target > hc + 0.02:
        # Each methyl turns one edge C-H into C-CH3: +1 C, net +2 H. Shrink the
        # core by the estimate so total C stays ~N, then size m exactly on the
        # rebuilt structure: H/C = (H + 2m)/(C + m) = target.
        n_est = round(N * (target - hc) / (2.0 - hc))
        n_est = max(1, min(n_est, core_n - max(arms)))
        for shift in (0, 1, -1, 2, -2):
            try:
                mol, parents = _build_candidate(plan, core_n - n_est - shift, arms, units)
                n_est += shift
                break
            except SkeletonError:
                continue
        else:
            raise SkeletonError(f"could not build a {core_n - n_est}-carbon core for methyl decoration")
        carbons = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
        h_now = _saturated_hc(mol) * carbons
        edge = sum(
            1 for a in mol.GetAtoms()
            if a.GetAtomicNum() == 6 and a.GetIsAromatic() and a.GetDegree() < 3
        )
        methyls = max(0, min(round((target * carbons - h_now) / (2.0 - target)), edge))
        if methyls:
            mol = attach_aliphatic_carbons(mol, methyls, random.Random(config.seed))
        core_n -= n_est

    unit_carbons = [
        sum(1 for a in mol.GetAtoms() if a.HasProp("unit_id") and a.GetIntProp("unit_id") == u)
        for u in range(k)
    ]
    counts: Dict[str, int] = {}
    for t in plan.linkers:
        counts[t] = counts.get(t, 0) + 1
    census = LinkedSkeleton(
        unit_carbons=unit_carbons, linker_counts=counts, num_methyls=methyls, parents=parents
    )
    skeleton = PAHAssembler._make_skeleton(mol)
    total = skeleton.num_carbons
    if abs(total - N) > config.size_tolerance * N:
        logger.warning(
            "Linked units: built %d carbons for a %d-carbon target (outside the "
            "%.0f%% size tolerance).", total, N, config.size_tolerance * 100,
        )
    return skeleton, census


def unit_bridging_atoms(mol: Chem.Mol) -> List[int]:
    """Non-linker heteroatoms bonded to two different units.

    The oxygen stage can place a bridging group (an explicit ``ether`` at the
    ortho carbons of a biaryl closes a furan) across two units. That merges two
    ring systems and closes a cycle in the unit tree, which neither the census
    nor the embedder can represent.
    """
    out = []
    for a in mol.GetAtoms():
        if a.GetAtomicNum() in (1, 6) or a.HasProp("linker_id") or a.HasProp("unit_id"):
            continue
        units = {
            n.GetIntProp("unit_id") for n in a.GetNeighbors() if n.HasProp("unit_id")
        }
        if len(units) > 1:
            out.append(a.GetIdx())
    return out


# ---------------------------------------------------------------------------
# 3D embedding
# ---------------------------------------------------------------------------

# Bond lengths (Å) for atoms placed off the lattice, by element pair; anything
# else falls back to a covalent-radius sum.
_BOND_LENGTHS = {
    frozenset((6, 1)): 1.09,
    frozenset((8, 1)): 0.96,
    frozenset((7, 1)): 1.01,
    frozenset((16, 1)): 1.34,
    frozenset((6, 7)): 1.40,
    frozenset((6, 16)): 1.77,
}
_COVALENT = {1: 0.31, 6: 0.76, 7: 0.71, 8: 0.66, 16: 1.05}


def _bond_length(mol: Chem.Mol, i: int, j: int) -> float:
    ai, aj = mol.GetAtomWithIdx(i), mol.GetAtomWithIdx(j)
    zi, zj = ai.GetAtomicNum(), aj.GetAtomicNum()
    bond = mol.GetBondBetweenAtoms(i, j)
    pair = frozenset((zi, zj))
    if pair == frozenset((6, 8)):
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            return 1.23
        carbon = ai if zi == 6 else aj
        return 1.37 if carbon.GetIsAromatic() else 1.43
    if pair == frozenset((6,)):
        if ai.GetIsAromatic() and aj.GetIsAromatic():
            return 1.48  # biaryl
        if ai.GetIsAromatic() or aj.GetIsAromatic():
            return 1.51
        return 1.53
    if pair in _BOND_LENGTHS:
        return _BOND_LENGTHS[pair]
    return _COVALENT.get(zi, 0.76) + _COVALENT.get(zj, 0.76)


def _ideal_angle(mol: Chem.Mol, idx: int) -> float:
    """Bond angle (radians) about atom *idx* for placing its neighbours."""
    a = mol.GetAtomWithIdx(idx)
    z, deg = a.GetAtomicNum(), a.GetDegree()
    trigonal = a.GetIsAromatic() or any(
        b.GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.AROMATIC) for b in a.GetBonds()
    )
    if z == 8 and deg == 2:
        heavy = sum(1 for n in a.GetNeighbors() if n.GetAtomicNum() != 1)
        return math.radians(118.0 if heavy == 2 else 109.0)
    if z == 16 and deg == 2:
        return math.radians(100.0)
    if trigonal or deg == 3:
        return math.radians(120.0)
    return math.radians(109.47)


def _unit(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    return v / n if n > 1e-9 else v


def _perpendicular(u: np.ndarray, ref: np.ndarray) -> np.ndarray:
    """A unit vector perpendicular to *u*, as close to *ref* as possible."""
    p = ref - (ref @ u) * u
    if np.linalg.norm(p) < 1e-6:
        p = np.cross(u, [1.0, 0.0, 0.0])
        if np.linalg.norm(p) < 1e-6:
            p = np.cross(u, [0.0, 1.0, 0.0])
    return _unit(p)


def _new_directions(placed: List[np.ndarray], n_new: int, angle: float, total: int,
                    normal: np.ndarray) -> List[np.ndarray]:
    """
    Directions for *n_new* further neighbours of an atom whose placed
    neighbours lie along *placed* (unit vectors), for an atom with *total*
    neighbours in all. *normal* is the local plane normal: trigonal and bent
    centres keep their new neighbours in that plane, the way an aromatic ring
    keeps its hydrogens and a phenol its O-H.
    """
    if n_new <= 0:
        return []
    if not placed:
        return [_unit(np.array([1.0, 0.0, 0.0]))] + _new_directions(
            [np.array([1.0, 0.0, 0.0])], n_new - 1, angle, total, normal
        )
    if len(placed) == 1:
        u = placed[0]
        q = _unit(np.cross(normal, u)) if np.linalg.norm(np.cross(normal, u)) > 1e-6 \
            else _perpendicular(u, np.array([0.0, 0.0, 1.0]))
        if total >= 4:
            # Tetrahedral: one neighbour in plane, two either side of it.
            n = _perpendicular(u, normal)
            return [
                _unit(math.cos(angle) * u + math.sin(angle) * (math.cos(phi) * q + math.sin(phi) * n))
                for phi in (math.radians(90.0), math.radians(210.0), math.radians(330.0))
            ][:n_new]
        dirs = [
            _unit(math.cos(angle) * u + math.sin(angle) * q),
            _unit(math.cos(angle) * u - math.sin(angle) * q),
        ]
        return dirs[:n_new]
    if len(placed) == 2:
        b = _unit(-(placed[0] + placed[1]))
        if total >= 4:
            w = _unit(np.cross(placed[0], placed[1]))
            half = math.radians(54.75)
            return [
                _unit(math.cos(half) * b + math.sin(half) * w),
                _unit(math.cos(half) * b - math.sin(half) * w),
            ][:n_new]
        return [b][:n_new]
    return [_unit(-np.sum(placed, axis=0))][:n_new]


def _place_outward(mol: Chem.Mol, coords: np.ndarray, placed: np.ndarray,
                   start: Sequence[int], blocked: set, normal: np.ndarray) -> None:
    """Breadth-first placement of unplaced atoms reachable from *start*, not
    crossing *blocked*. Writes into *coords* and *placed* in place."""
    queue = list(start)
    while queue:
        x = queue.pop(0)
        atom = mol.GetAtomWithIdx(x)
        nbrs = [n.GetIdx() for n in atom.GetNeighbors()]
        done = [j for j in nbrs if placed[j]]
        todo = [j for j in nbrs if not placed[j] and j not in blocked]
        if not todo:
            continue
        # Heavy atoms first so hydrogens fill the remaining directions.
        todo.sort(key=lambda j: mol.GetAtomWithIdx(j).GetAtomicNum() == 1)
        dirs = _new_directions(
            [_unit(coords[j] - coords[x]) for j in done], len(todo),
            _ideal_angle(mol, x), len(nbrs), normal,
        )
        for j, d in zip(todo, dirs):
            coords[j] = coords[x] + d * _bond_length(mol, x, j)
            placed[j] = True
            queue.append(j)


def _nerf(a: np.ndarray, b: np.ndarray, c: np.ndarray, length: float,
          angle: float, torsion: float) -> np.ndarray:
    """Place D bonded to C with |CD| = length, angle BCD, torsion ABCD."""
    bc = _unit(c - b)
    n = _unit(np.cross(b - a, bc))
    if np.linalg.norm(n) < 1e-6:
        n = _perpendicular(bc, np.array([0.0, 0.0, 1.0]))
    m = np.cross(n, bc)
    d2 = np.array([
        -length * math.cos(angle),
        length * math.sin(angle) * math.cos(torsion),
        length * math.sin(angle) * math.sin(torsion),
    ])
    return c + d2[0] * bc + d2[1] * m + d2[2] * n


def _rotation_between(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Rotation matrix taking unit vector *u* onto unit vector *v*."""
    w = np.cross(u, v)
    s, c = float(np.linalg.norm(w)), float(u @ v)
    if s < 1e-9:
        if c > 0:
            return np.eye(3)
        axis = _perpendicular(u, np.array([0.0, 0.0, 1.0]))
        return 2.0 * np.outer(axis, axis) - np.eye(3)
    k = w / s
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + s * K + (1 - c) * (K @ K)


def _axis_rotation(axis: np.ndarray, angle: float) -> np.ndarray:
    k = _unit(axis)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + math.sin(angle) * K + (1 - math.cos(angle)) * (K @ K)


@dataclass
class _Link:
    """One linker as the embedder sees it: parent and child attachment carbons,
    the bridge atoms between them, and the child unit."""

    parent_atom: int
    child_atom: int
    bridge: List[int]
    child_unit: int


class _LinkedLayout:
    """Rigid units plus the links between them, posable by link torsions."""

    def __init__(self, mol: Chem.Mol):
        self.mol = mol
        n = mol.GetNumAtoms()
        self.unit_of = np.full(n, -1)
        for a in mol.GetAtoms():
            if a.HasProp("unit_id"):
                self.unit_of[a.GetIdx()] = a.GetIntProp("unit_id")
        self.num_units = int(self.unit_of.max()) + 1
        self.linker_atoms = {a.GetIdx() for a in mol.GetAtoms() if a.HasProp("linker_id")}
        self.links = self._find_links()
        self.body = self._bodies()
        self.local = self._local_coordinates()

    def _find_links(self) -> List[_Link]:
        mol = self.mol
        by_id: Dict[int, List[int]] = {}
        for i in self.linker_atoms:
            by_id.setdefault(mol.GetAtomWithIdx(i).GetIntProp("linker_id"), []).append(i)
        links = []
        # Biaryl links have no bridge atoms: find them as inter-unit C-C bonds.
        for b in mol.GetBonds():
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            ui, uj = self.unit_of[i], self.unit_of[j]
            if ui >= 0 and uj >= 0 and ui != uj:
                links.append((min(ui, uj), max(ui, uj), i if ui < uj else j, j if ui < uj else i, []))
        for _lid, atoms in by_id.items():
            ends = []
            for i in atoms:
                for nb in mol.GetAtomWithIdx(i).GetNeighbors():
                    if self.unit_of[nb.GetIdx()] >= 0:
                        ends.append((int(self.unit_of[nb.GetIdx()]), nb.GetIdx(), i))
            (u1, a1, b1), (u2, a2, b2) = sorted(ends)[:2]
            # Order bridge atoms from the lower unit's end to the higher's.
            path = [b1]
            while path[-1] != b2:
                nxt = [
                    nb.GetIdx() for nb in mol.GetAtomWithIdx(path[-1]).GetNeighbors()
                    if nb.GetIdx() in atoms and nb.GetIdx() not in path
                ]
                path.append(nxt[0])
            links.append((u1, u2, a1, a2, path))
        # Orient as a tree from the core (unit 0) outward.
        placed = {0}
        ordered: List[_Link] = []
        pending = list(links)
        while pending:
            progressed = False
            for lk in list(pending):
                u1, u2, a1, a2, path = lk
                if u1 in placed and u2 not in placed:
                    ordered.append(_Link(a1, a2, list(path), u2))
                    placed.add(u2)
                elif u2 in placed and u1 not in placed:
                    ordered.append(_Link(a2, a1, list(reversed(path)), u1))
                    placed.add(u1)
                else:
                    continue
                pending.remove(lk)
                progressed = True
            if not progressed:
                raise RuntimeError("linked units do not form a tree")
        return ordered

    def _bodies(self) -> List[List[int]]:
        """Each unit's atoms plus everything hanging off it (not linkers)."""
        mol = self.mol
        owner = self.unit_of.copy()
        queue = [i for i in range(mol.GetNumAtoms()) if owner[i] >= 0]
        while queue:
            x = queue.pop(0)
            for nb in mol.GetAtomWithIdx(x).GetNeighbors():
                j = nb.GetIdx()
                if owner[j] < 0 and j not in self.linker_atoms:
                    # A hydrogen on a linker carbon belongs to the linker.
                    if nb.GetAtomicNum() == 1 and x in self.linker_atoms:
                        continue
                    owner[j] = owner[x]
                    queue.append(j)
        self.owner = owner
        return [[int(i) for i in np.nonzero(owner == u)[0]] for u in range(self.num_units)]

    def _local_coordinates(self) -> np.ndarray:
        """Every body atom in its own unit's frame: lattice at z = 0, then
        substituents and hydrogens placed outward."""
        mol = self.mol
        n = mol.GetNumAtoms()
        coords = np.zeros((n, 3))
        placed = np.zeros(n, dtype=bool)
        normal = np.array([0.0, 0.0, 1.0])
        for u in range(self.num_units):
            ring = [i for i in self.body[u] if self.unit_of[i] == u]
            for i in ring:
                a = mol.GetAtomWithIdx(i)
                if a.HasProp("init_x") and a.HasProp("init_y"):
                    coords[i] = (a.GetDoubleProp("init_x"), a.GetDoubleProp("init_y"), 0.0)
                    placed[i] = True
            missing = [i for i in ring if not placed[i]]
            if missing:
                raise RuntimeError(f"unit {u} has atoms without lattice positions: {missing[:5]}")
            centre = coords[ring].mean(axis=0)
            coords[ring] -= centre
            # A group bonded to two atoms of this unit (an in-unit ether
            # closing a furan at a bay) sits between its anchors, offset
            # outward so both bonds have their length; placing it from one
            # anchor alone leaves the other bond 2-5 Å long.
            ring_set = set(ring)
            for i in self.body[u]:
                if placed[i] or i in ring_set:
                    continue
                anchors = [nb.GetIdx() for nb in mol.GetAtomWithIdx(i).GetNeighbors()
                           if nb.GetIdx() in ring_set]
                if len(anchors) < 2:
                    continue
                p1, p2 = coords[anchors[0]], coords[anchors[1]]
                mid = 0.5 * (p1 + p2)
                half = 0.5 * float(np.linalg.norm(p2 - p1))
                length = _bond_length(mol, anchors[0], i)
                outward = _unit(mid - coords[ring].mean(axis=0))
                along = _unit(p2 - p1)
                outward = _unit(outward - (outward @ along) * along)
                coords[i] = mid + outward * math.sqrt(max(length ** 2 - half ** 2, 0.0))
                placed[i] = True
            blocked = self.linker_atoms | {
                i for i in range(n) if self.owner[i] != u and self.owner[i] >= 0
            }
            # A substituent with two ring anchors (an in-unit ether bridge)
            # is placed from whichever anchor reaches it first.
            _place_outward(mol, coords, placed, ring, blocked, normal)
        self.placed_local = placed
        return coords

    def pose(self, torsions: Dict[int, List[float]]) -> np.ndarray:
        """Global coordinates for one choice of link torsions."""
        mol = self.mol
        n = mol.GetNumAtoms()
        out = np.zeros((n, 3))
        done = np.zeros(n, dtype=bool)
        z = np.array([0.0, 0.0, 1.0])
        frames = {0: (np.eye(3), np.zeros(3))}  # unit -> (R, t): global = R @ local + t
        out[self.body[0]] = self.local[self.body[0]]
        done[self.body[0]] = True
        for li, link in enumerate(self.links):
            tors = torsions[li]
            a = link.parent_atom
            ring_nbrs = [
                nb.GetIdx() for nb in mol.GetAtomWithIdx(a).GetNeighbors()
                if self.unit_of[nb.GetIdx()] == self.unit_of[a]
            ]
            outward = _unit(out[a] - out[ring_nbrs].mean(axis=0))
            chain = [ring_nbrs[0], a]
            pos = {ring_nbrs[0]: out[ring_nbrs[0]], a: out[a]}
            atoms = link.bridge + [link.child_atom]
            for step, x in enumerate(atoms):
                length = _bond_length(mol, chain[-1], x)
                if step == 0:
                    p = out[a] + outward * length
                else:
                    p = _nerf(pos[chain[-3]], pos[chain[-2]], pos[chain[-1]], length,
                              _ideal_angle(mol, chain[-1]), tors[step - 1])
                pos[x] = p
                chain.append(x)
            for x in link.bridge:
                out[x] = pos[x]
                done[x] = True
            # Rigidly place the child body: its attachment carbon at the end of
            # the chain, its own outward direction pointing back along the last
            # bond, then turned about that bond by the last torsion.
            b = link.child_atom
            u = link.child_unit
            c_nbrs = [
                nb.GetIdx() for nb in mol.GetAtomWithIdx(b).GetNeighbors()
                if self.unit_of[nb.GetIdx()] == u
            ]
            local_out = _unit(self.local[b] - self.local[c_nbrs].mean(axis=0))
            back = _unit(pos[chain[-2]] - pos[b])
            R1 = _rotation_between(local_out, back)
            # Start the child plane coplanar with the previous bond's plane so
            # the twist torsion is measured from a defined reference.
            local_normal = R1 @ z
            ref = _perpendicular(back, np.cross(back, pos[chain[-3]] - pos[chain[-2]])
                                 if len(chain) >= 3 else z)
            ref_now = _perpendicular(back, local_normal)
            align = math.atan2(np.cross(ref_now, ref) @ back, ref_now @ ref)
            R = _axis_rotation(back, align + tors[-1]) @ R1
            t = pos[b] - R @ self.local[b]
            frames[u] = (R, t)
            idx = self.body[u]
            out[idx] = self.local[idx] @ R.T + t
            done[idx] = True
        # Linker hydrogens, now that both ends of every bridge are fixed.
        bridge_atoms = [x for link in self.links for x in link.bridge]
        _place_outward(mol, out, done, bridge_atoms, set(), z)
        return out

    def torsion_counts(self) -> List[int]:
        # One torsion per bridge atom after the first, plus the child twist.
        return [max(len(link.bridge), 0) + 1 for link in self.links]

    def unit_plane_rms(self, coords: np.ndarray) -> float:
        worst = 0.0
        for u in range(self.num_units):
            ring = [i for i in range(len(coords)) if self.unit_of[i] == u]
            centred = coords[ring] - coords[ring].mean(axis=0)
            normal = np.linalg.svd(centred, full_matrices=False)[2][-1]
            worst = max(worst, float(np.sqrt(np.mean((centred @ normal) ** 2))))
        return worst


# Biphenyl twists ~44° in the gas phase; a coplanar start has the ortho
# hydrogens 2.0 Å apart, outside the clash floor but not where a real biaryl
# sits, so biaryl twists are drawn away from coplanar.
_BIARYL_TWIST_RANGE = (math.radians(30.0), math.radians(150.0))
_EMBED_TRIALS = 60
_EMBED_SWEEPS = 2
_SWEEP_ANGLES = np.radians(np.arange(0.0, 360.0, 30.0))


def _sample_torsions(layout: _LinkedLayout, rng: np.random.Generator) -> Dict[int, List[float]]:
    out = {}
    for li, link in enumerate(layout.links):
        count = len(link.bridge) + 1
        tors = list(rng.uniform(0.0, 2 * math.pi, count))
        # Alkyl torsions start near anti or gauche, where a real chain sits.
        for j in range(count - 1):
            tors[j] = rng.choice([math.pi, math.pi / 3, -math.pi / 3]) + rng.normal(0.0, 0.15)
        if not link.bridge:
            lo, hi = _BIARYL_TWIST_RANGE
            tors[-1] = rng.uniform(lo, hi) * rng.choice([1.0, -1.0])
        out[li] = tors
    return out


def _overlap_score(coords: np.ndarray, floor: np.ndarray, judged: np.ndarray) -> float:
    d = np.sqrt(((coords[:, None, :] - coords[None, :, :]) ** 2).sum(-1))
    return float(np.where(judged, np.clip(floor - d, 0.0, None), 0.0).sum())


def embed_linked_units(mol: Chem.Mol, seed: Optional[int] = None) -> Tuple[Chem.Mol, np.ndarray]:
    """
    3D coordinates for a fully decorated linked molecule.

    Units are placed rigidly from their own hex-lattice positions, so each stays
    exactly flat; linker atoms are built from bond lengths, angles and link
    torsions; the torsions are searched to minimise overlap against the clash
    floors the validator uses. Torsions cannot move atoms inside a unit, so the
    two hex-lattice relief stages follow: crowded substituents tilt out of
    plane, then a tethered UFF pass is kept only if it adds no defect and every
    unit stays flat. ``relieve_lattice_crowding`` itself is not used: it judges
    flatness by one plane through every ring atom, which rejects any structure
    whose units are rotated apart.
    """
    layout = _LinkedLayout(mol)
    rng = np.random.default_rng(0 if seed is None else seed)
    excluded = _get_excluded_pairs(mol)
    floor, judged = _clash_floor_matrix(mol, set(), excluded)

    best_t, best_s = None, float("inf")
    for _ in range(_EMBED_TRIALS):
        t = _sample_torsions(layout, rng)
        s = _overlap_score(layout.pose(t), floor, judged)
        if s < best_s:
            best_t, best_s = t, s
        if s == 0.0:
            break
    # Coordinate descent over each torsion, keeping biaryl twists off-planar.
    for _ in range(_EMBED_SWEEPS):
        if best_s == 0.0:
            break
        for li, link in enumerate(layout.links):
            for j in range(len(best_t[li])):
                for ang in _SWEEP_ANGLES:
                    if not link.bridge:
                        lo, hi = _BIARYL_TWIST_RANGE
                        if not lo <= abs(((ang + math.pi) % (2 * math.pi)) - math.pi) <= hi:
                            continue
                    trial = {k: list(v) for k, v in best_t.items()}
                    trial[li][j] = ang
                    s = _overlap_score(layout.pose(trial), floor, judged)
                    if s < best_s:
                        best_t, best_s = trial, s
    coords = layout.pose(best_t)

    if _clash_pairs(mol, coords, excluded):
        coords = _tilt_crowded_substituents(mol, coords)
    if _clash_pairs(mol, coords, excluded):
        relaxed = _tethered_lattice_relax(mol, coords, seed)
        if relaxed is not None:
            def defects(c):
                return (
                    len(GeometryValidator._validate_bond_lengths(mol, c)),
                    len(_clash_pairs(mol, c, excluded)),
                )
            if (
                defects(relaxed) <= defects(coords)
                and layout.unit_plane_rms(relaxed) <= LATTICE_MAX_SHEET_RMS
            ):
                coords = relaxed

    out = Chem.Mol(mol)
    conf = Chem.Conformer(out.GetNumAtoms())
    for i, xyz in enumerate(coords):
        conf.SetAtomPosition(i, xyz.tolist())
    out.RemoveAllConformers()
    out.AddConformer(conf, assignId=True)
    return out, coords
