"""
3D Geometry Generation and Validation

Generate 3D coordinates for biochar molecules and validate geometry.
"""

import logging
from typing import Tuple, List, Optional
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial.distance import pdist, squareform

from ..constants import (
    VDW_RADII,
    COVALENT_RADII,
    HBOND_MIN_H_ACCEPTOR_DISTANCE,
    HBOND_MIN_DHA_ANGLE_DEG,
    HBOND_DONOR_ACCEPTOR_ELEMENTS,
    CLASH_SEVERITY_TOLERANCE,
    BOND_ORDER_LENGTH_FACTORS,
    BOND_LENGTH_MIN_FACTOR,
    BOND_LENGTH_MAX_FACTOR,
    HYDROXYL_H_PLACEMENT_ANGLE_DEG,
    LATTICE_TILT_MAX_DEG,
    LATTICE_TILT_STEP_DEG,
    LATTICE_TETHER_RADIUS,
    LATTICE_TETHER_FORCE_CONSTANT,
    LATTICE_BOND_RANGE,
    LATTICE_BOND_FORCE_CONSTANT,
    LATTICE_RELAX_MAX_ITERATIONS,
    LATTICE_MAX_SHEET_RMS,
)

logger = logging.getLogger(__name__)


def _get_excluded_pairs(mol: Chem.Mol) -> set:
    """
    Return the set of (min_idx, max_idx) atom-index pairs that should be
    **excluded** from non-bonded clash detection.

    Standard molecular-mechanics exclusion rules:
    * 1-2 pairs (directly bonded)
    * 1-3 pairs (angle neighbours, sharing one atom)

    In a flat graphene-like sheet the 1-3 distance is
    √3 × L ≈ 2.46 Å (L = 1.42 Å aromatic bond), which is below the
    0.75 × vdW-sum threshold for C-C (2.55 Å).  Without this exclusion
    every second-neighbour carbon pair is flagged as a clash, producing
    hundreds of false positives that the resolver can never fix.
    """
    excluded: set = set()
    for atom in mol.GetAtoms():
        i = atom.GetIdx()
        for nbr in atom.GetNeighbors():
            j = nbr.GetIdx()
            # 1-2 pair
            excluded.add((min(i, j), max(i, j)))
            # 1-3 pairs: all neighbours of j except i
            for nbr2 in nbr.GetNeighbors():
                k = nbr2.GetIdx()
                if k != i:
                    excluded.add((min(i, k), max(i, k)))
    return excluded


def _get_hbond_pairs(mol: Chem.Mol, coords: np.ndarray) -> set:
    """
    Return the set of (min_idx, max_idx) pairs that are hydrogen bonds rather
    than steric clashes.

    A pair qualifies when a polar hydrogen (one covalently bound to N or O)
    points at an N/O acceptor: the H...A distance is within the physical H-bond
    range and the D-H...A angle is at least
    :data:`~biochar.constants.HBOND_MIN_DHA_ANGLE_DEG`.

    Callers hold these pairs to
    :data:`~biochar.constants.HBOND_MIN_H_ACCEPTOR_DISTANCE` instead of the
    generic 0.75 × vdW-sum floor.  They are *not* excluded outright, so a real
    overlap between a donor H and an acceptor is still caught.

    This mirrors the 1-2/1-3 exclusion in :func:`_get_excluded_pairs`: both
    strip contacts that are ordinary molecular geometry rather than atom
    overlap.  The acceptor's own D...A heavy-atom distance (~2.6-2.8 Å for
    O...O) already clears the 2.28 Å O/O floor, so only the H...A pair needs
    special handling.
    """
    # Map each polar hydrogen to the donor heavy atom it is bonded to.
    polar_h: dict = {}
    acceptors: set = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in HBOND_DONOR_ACCEPTOR_ELEMENTS:
            acceptors.add(atom.GetIdx())
        elif atom.GetAtomicNum() == 1:
            nbrs = atom.GetNeighbors()
            if len(nbrs) == 1 and nbrs[0].GetAtomicNum() in HBOND_DONOR_ACCEPTOR_ELEMENTS:
                polar_h[atom.GetIdx()] = nbrs[0].GetIdx()

    if not polar_h or not acceptors:
        return set()

    hbonds: set = set()
    for h_idx, donor_idx in polar_h.items():
        h_pos = coords[h_idx]
        # Vector H->D, used for the D-H...A angle.
        v_hd = coords[donor_idx] - h_pos
        n_hd = np.linalg.norm(v_hd)
        if n_hd < 1e-6:
            continue
        for a_idx in acceptors:
            # The donor itself is a 1-2 pair; both are already excluded upstream,
            # but skip explicitly so this helper stands on its own.
            if a_idx == donor_idx:
                continue
            v_ha = coords[a_idx] - h_pos
            n_ha = np.linalg.norm(v_ha)
            if n_ha < 1e-6:
                continue
            cos_theta = float(np.dot(v_hd, v_ha) / (n_hd * n_ha))
            angle = np.degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0)))
            if angle >= HBOND_MIN_DHA_ANGLE_DEG:
                hbonds.add((min(h_idx, a_idx), max(h_idx, a_idx)))
    return hbonds


def _clash_floor(
    mol: Chem.Mol, i: int, j: int, hbond_pairs: set
) -> float:
    """
    Minimum allowed non-bonded distance (Å) for the atom pair *(i, j)*.

    Defaults to 0.75 × vdW-sum — a real clash is atom overlap, so this keeps
    normal non-bonded proximity (e.g. peri H-H ≈ 1.95 Å in PAHs) from being
    flagged.  Hydrogen-bonded donor-acceptor pairs get the reduced H-bond
    floor instead.
    """
    if (min(i, j), max(i, j)) in hbond_pairs:
        return HBOND_MIN_H_ACCEPTOR_DISTANCE
    r_vdw_i = VDW_RADII.get(mol.GetAtomWithIdx(i).GetSymbol(), 1.70)
    r_vdw_j = VDW_RADII.get(mol.GetAtomWithIdx(j).GetSymbol(), 1.70)
    return 0.75 * (r_vdw_i + r_vdw_j)


def _clash_floor_matrix(
    mol: Chem.Mol, hbond_pairs: set, excluded: set
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Every pair's :func:`_clash_floor` at once, and a mask of the pairs that are
    judged at all (not self, not 1-2 or 1-3).
    """
    n = mol.GetNumAtoms()
    radii = np.array([VDW_RADII.get(a.GetSymbol(), 1.70) for a in mol.GetAtoms()])
    floor = 0.75 * (radii[:, None] + radii[None, :])
    for i, j in hbond_pairs:
        floor[i, j] = floor[j, i] = HBOND_MIN_H_ACCEPTOR_DISTANCE
    judged = np.ones((n, n), dtype=bool)
    np.fill_diagonal(judged, False)
    for i, j in excluded:
        judged[i, j] = judged[j, i] = False
    return floor, judged


def _clash_pairs(
    mol: Chem.Mol, coords: np.ndarray, excluded: Optional[set] = None
) -> List[Tuple[int, int, float, float]]:
    """
    Every steric clash as ``(i, j, distance, floor)`` with ``i < j``, in
    ascending ``(i, j)`` order.

    The one definition of a clash: a judged pair more than
    :data:`~biochar.constants.CLASH_SEVERITY_TOLERANCE` inside its floor.
    :meth:`GeometryValidator._check_steric_clashes` reports these, and the
    hex-lattice relief pass resolves them, so the two cannot disagree about
    what counts.
    """
    if excluded is None:
        excluded = _get_excluded_pairs(mol)
    floor, judged = _clash_floor_matrix(mol, _get_hbond_pairs(mol, coords), excluded)
    distances = squareform(pdist(coords))
    clashing = np.triu(judged & (distances < floor - CLASH_SEVERITY_TOLERANCE))
    return [
        (int(i), int(j), float(distances[i, j]), float(floor[i, j]))
        for i, j in zip(*np.nonzero(clashing))
    ]


def _read_skeleton_positions(mol: Chem.Mol) -> dict:
    """
    Collect hex-lattice positions tagged on the molecule by
    :mod:`carbon_skeleton` (atom properties ``init_x`` / ``init_y``).

    Returns ``{atom_idx: np.array([x, y, 0])}`` for heavy atoms that have
    both props.  Empty dict if no atom is tagged.
    """
    positions: dict = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        if atom.HasProp("init_x") and atom.HasProp("init_y"):
            try:
                x = atom.GetDoubleProp("init_x")
                y = atom.GetDoubleProp("init_y")
            except Exception:
                continue
            positions[atom.GetIdx()] = np.array([x, y, 0.0])
    return positions


def _lattice_covers_skeleton(mol: Chem.Mol, skel_pos: dict) -> bool:
    """
    Whether the hex-lattice tags are enough to lay *mol* out flat: every ring
    carbon is tagged, and every untagged heavy atom hangs, through other
    untagged atoms, off a tagged one.

    Untagged atoms outside the rings -- pendant methyls from aliphatic
    decoration, a carboxyl carbon, the oxygens -- are placed from their parent
    by :func:`_offlattice_layout`. Requiring *every* carbon to be tagged sent a
    whole sheet to the Kamada-Kawai fallback over one methyl (issue #66).
    """
    if not skel_pos:
        return False
    # *mol* is the throwaway heavy-atom copy, whose ring information is never
    # initialised; perceiving it here touches nothing the caller holds.
    Chem.FastFindRings(mol)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.IsInRing() and atom.GetIdx() not in skel_pos:
            return False
    reached, stack = set(skel_pos), list(skel_pos)
    while stack:
        for nbr in mol.GetAtomWithIdx(stack.pop()).GetNeighbors():
            if nbr.GetAtomicNum() != 1 and nbr.GetIdx() not in reached:
                reached.add(nbr.GetIdx())
                stack.append(nbr.GetIdx())
    return all(a.GetIdx() in reached for a in mol.GetAtoms() if a.GetAtomicNum() != 1)


def _bond_slots(axis: np.ndarray, side: np.ndarray, n_bonds: int, element: int) -> List[np.ndarray]:
    """
    Unit directions of the bonds an off-lattice atom makes besides the one to
    its parent, which lies along -*axis*. In-plane directions come first, so
    heavy children take them and hydrogens take what is left.

    Four bonds is tetrahedral: one bond in the sheet plane on *side*, two
    symmetric above and below it. Three is trigonal in plane. Two is linear for
    carbon and bent (tetrahedral angle, in plane) for anything else -- an ether
    or hydroxyl oxygen.
    """
    z = np.cross(axis, side)
    along, across = 1.0 / 3.0, np.sqrt(8.0) / 3.0  # cos/sin of 180 - 109.47
    if n_bonds >= 4:
        return [
            along * axis + across * (np.cos(phi) * side + np.sin(phi) * z)
            for phi in np.radians([0.0, 120.0, 240.0])
        ]
    if n_bonds == 3:
        return [0.5 * axis + np.sqrt(3.0) / 2.0 * sign * side for sign in (1.0, -1.0)]
    if n_bonds == 2:
        if element == 6:
            return [axis]
        return [along * axis + across * side]
    return []


def _offlattice_bond_length(parent: Chem.Atom, atom: Chem.Atom, bond: Chem.Bond) -> float:
    """Bond length for an off-lattice atom: phenolic C-O keeps its 1.36 Å;
    anything else is the bond-order-scaled covalent radii sum, which is also
    what :meth:`GeometryValidator._validate_bond_lengths` expects."""
    pair = {parent.GetAtomicNum(), atom.GetAtomicNum()}
    if pair == {6, 8} and (parent.GetIsAromatic() or atom.GetIsAromatic()):
        return 1.36
    radii = COVALENT_RADII.get(parent.GetSymbol(), 0.76) + COVALENT_RADII.get(atom.GetSymbol(), 0.76)
    return radii * BOND_ORDER_LENGTH_FACTORS.get(str(bond.GetBondType()), 1.00)


def _offlattice_layout(
    mol: Chem.Mol, placed: dict, chemistry: Chem.Mol
) -> Tuple[dict, dict]:
    """
    Place every untagged heavy atom of *mol* (all atoms explicit, hydrogens
    included) outward from the tagged ones, and return ``(positions, slots)``.

    Bond lengths are read from *chemistry*, the same molecule before it was
    kekulised or de-aromatised for embedding: the working copy has lost the
    aromatic flags that tell a phenolic C-O from an aliphatic one.

    A child of a ring atom points radially out of the sheet, away from the
    ring atom's other placed neighbours. A child of an off-lattice atom
    takes one of that atom's :func:`_bond_slots`, so a methyl is tetrahedral and
    a carboxyl trigonal instead of every substituent stacking on one line.
    *slots* maps each off-lattice atom to the directions still free for its
    hydrogens. Atoms are placed breadth-first in index order, so the layout is
    deterministic.
    """
    positions = dict(placed)
    slots: dict = {}
    frontier = sorted(placed)
    while frontier:
        next_frontier = []
        for p_idx in frontier:
            parent = mol.GetAtomWithIdx(p_idx)
            for atom in sorted(parent.GetNeighbors(), key=lambda a: a.GetIdx()):
                idx = atom.GetIdx()
                if atom.GetAtomicNum() == 1 or idx in positions:
                    continue
                bl = _offlattice_bond_length(
                    chemistry.GetAtomWithIdx(p_idx), chemistry.GetAtomWithIdx(idx),
                    chemistry.GetBondBetweenAtoms(p_idx, idx),
                )
                if p_idx in slots:
                    direction = slots[p_idx].pop(0)
                else:
                    others = [
                        positions[n.GetIdx()] for n in parent.GetNeighbors()
                        if n.GetIdx() != idx and n.GetIdx() in positions
                    ]
                    direction = positions[p_idx] - np.mean(others, axis=0) if others else np.array([1.0, 0.0, 0.0])
                    norm = float(np.linalg.norm(direction))
                    direction = direction / norm if norm > 1e-8 else np.array([1.0, 0.0, 0.0])
                pos = positions[p_idx] + direction * bl
                positions[idx] = pos

                # The in-plane slot goes on whichever side has more room.
                side = np.cross(np.array([0.0, 0.0, 1.0]), direction)
                if np.linalg.norm(side) < 1e-8:
                    side = _perpendicular_unit(direction)
                side = side / np.linalg.norm(side)
                crowd = [q for k, q in positions.items() if k not in (idx, p_idx)]
                if crowd:
                    def clearance(s):
                        probe = pos + s * bl
                        return min(float(np.linalg.norm(probe - q)) for q in crowd)
                    if clearance(-side) > clearance(side):
                        side = -side
                slots[idx] = _bond_slots(direction, side, atom.GetDegree(), atom.GetAtomicNum())
                next_frontier.append(idx)
        frontier = next_frontier
    return positions, slots


def _bent_hydrogen_position(
    parent_pos: np.ndarray,
    grandparent_pos: np.ndarray,
    bond_length: float,
    occupied: List[np.ndarray],
) -> np.ndarray:
    """
    Place a hydrogen on *parent* at :data:`HYDROXYL_H_PLACEMENT_ANGLE_DEG` to
    the grandparent bond, in the z = 0 plane of the flat layout, on whichever
    side has more room from the *occupied* heavy-atom positions.
    """
    axis = parent_pos - grandparent_pos
    axis = axis / np.linalg.norm(axis)
    side = np.cross(np.array([0.0, 0.0, 1.0]), axis)
    if np.linalg.norm(side) < 1e-8:
        side = _perpendicular_unit(axis)
    side = side / np.linalg.norm(side)
    bend = np.radians(180.0 - HYDROXYL_H_PLACEMENT_ANGLE_DEG)
    others = [p for p in occupied
              if np.linalg.norm(p - parent_pos) > 1e-6
              and np.linalg.norm(p - grandparent_pos) > 1e-6]

    def clearance(pos):
        return min((float(np.linalg.norm(pos - p)) for p in others), default=np.inf)

    candidates = [
        parent_pos + bond_length * (np.cos(bend) * axis + sign * np.sin(bend) * side)
        for sign in (1.0, -1.0)
    ]
    return max(candidates, key=clearance)


def _heavy_atom_flat_layout(mol: Chem.Mol, seed: Optional[int] = None) -> np.ndarray:
    """
    Produce a flat (z=0) 2D layout for the heavy-atom scaffold of a PAH.

    Uses networkx.kamada_kawai_layout on the heavy-atom bond graph, then
    scales so the mean bond length is ~1.4 Å.  This is dramatically more
    robust on large fused polycyclic biochar graphs than either
    AllChem.Compute2DCoords (which collapses ~25 pairs to distance 0) or
    ETKDGv3 3D embedding + flatten (which folds then loses separation when
    projected).  Biochar graphs are planar by construction, so a planar
    spring-like layout converges cleanly.

    Args:
        mol: molecule with heavy atoms (no explicit hydrogens).
        seed: optional RNG seed for reproducibility of the layout.

    Returns:
        (N, 3) ndarray of 3D coords in Å, with z = 0 for every atom.
    """
    import networkx as nx

    G = nx.Graph()
    for a in mol.GetAtoms():
        G.add_node(a.GetIdx())
    for b in mol.GetBonds():
        G.add_edge(b.GetBeginAtomIdx(), b.GetEndAtomIdx())

    # Kamada-Kawai needs a connected graph; handle components independently.
    n = mol.GetNumAtoms()
    coords = np.zeros((n, 3))

    components = list(nx.connected_components(G))
    x_offset = 0.0
    for comp in components:
        subG = G.subgraph(comp).copy()
        if subG.number_of_nodes() == 1:
            only_node = next(iter(subG.nodes()))
            coords[only_node] = [x_offset, 0.0, 0.0]
            x_offset += 2.0
            continue
        pos = nx.kamada_kawai_layout(subG)
        sub_coords = np.array([pos[i] for i in sorted(pos.keys())])

        # Scale so the mean bond length is ~1.4 Å (aromatic C-C baseline).
        bls = [
            np.linalg.norm(np.array(pos[u]) - np.array(pos[v]))
            for u, v in subG.edges()
        ]
        if bls:
            scale = 1.4 / np.mean(bls)
            sub_coords = sub_coords * scale
        # Shift so component's min-x sits past previous component's max-x.
        min_x = sub_coords[:, 0].min()
        sub_coords[:, 0] += x_offset - min_x
        x_offset = sub_coords[:, 0].max() + 2.0

        for sub_idx, orig_idx in enumerate(sorted(pos.keys())):
            coords[orig_idx] = [sub_coords[sub_idx, 0], sub_coords[sub_idx, 1], 0.0]

    return coords


def _with_rings_perceived(mol: Chem.Mol) -> Chem.Mol:
    """Return *mol* with its ring information available, perceiving it if not.

    `FastFindRings` is a symmetrized-SSSR-free ring perception: it answers ring
    membership without sanitising, so it cannot disturb the bond types the
    heteroatom stages set. Idempotent -- a molecule that already knows its
    rings is returned untouched.
    """
    try:
        mol.GetRingInfo().NumRings()
    except RuntimeError:
        Chem.FastFindRings(mol)
    return mol


def _mmff_force_field(mol: Chem.Mol):
    """An MMFF94 force field for *mol*, or None when MMFF cannot type it.

    `MMFFGetMoleculeForceField` needs the properties object, and asking for a
    force field without one raises rather than returning None.
    """
    props = AllChem.MMFFGetMoleculeProperties(mol)
    if props is None:
        return None
    return AllChem.MMFFGetMoleculeForceField(mol, props)


def _adopt_conformer(mol: Chem.Mol, source: Chem.Mol) -> Chem.Mol:
    """Give *mol* the coordinates *source* carries, and return *mol*.

    The molecule embedding works on is not the molecule embedding was asked
    about: `_kekulize_or_dearomatize` hands back a copy whose bond orders were
    rewritten to something RDKit's embedders and force fields accept. Only its
    coordinates are a result — its chemistry is scaffolding, and on the
    de-aromatising branch it is a molecule with no aromatic ring in it at all.

    Falls back to *source* if the two have drifted apart in atom count, since
    then there is no index-for-index mapping to carry coordinates over.
    """
    if source.GetNumAtoms() != mol.GetNumAtoms():
        return source
    try:
        embedded = source.GetConformer(0)
    except (ValueError, RuntimeError):
        return source
    conformer = Chem.Conformer(mol.GetNumAtoms())
    for i in range(mol.GetNumAtoms()):
        conformer.SetAtomPosition(i, embedded.GetAtomPosition(i))
    mol.RemoveAllConformers()
    mol.AddConformer(conformer, assignId=True)
    return mol


def _kekulize_or_dearomatize(mol: Chem.Mol) -> Chem.Mol:
    """
    Ensure the returned mol is safe for FF/embedding by either kekulizing it
    or, if that fails (non-kekulizable graph — common for large biochar PAHs
    with odd-parity subgraphs), converting every aromatic bond to SINGLE and
    clearing all aromatic flags.

    What comes back is a *working copy*, not a result. The de-aromatising
    branch throws away every aromatic bond and flag in the molecule, so a
    caller that lets this copy escape hands its own caller a fused-ring sheet
    that no longer claims a single aromatic ring — see `_adopt_conformer`.

    After this call, `AllChem.Compute2DCoords`, `AllChem.EmbedMolecule`,
    `AllChem.UFFGetMoleculeForceField`, and `AllChem.MMFFGetMoleculeProperties`
    all work deterministically.

    The molecule also comes back with its rings perceived, whichever branch it
    took. `SanitizeMol` clears the computed properties -- ring information
    among them -- before it does anything else, so a sanitisation that raises
    partway leaves the molecule knowing less about itself than when it went in:
    every later "is this atom in a ring?" raises `RingInfo not initialized`
    rather than answering. Swallowing the sanitisation failure is deliberate
    (the geometry is still usable), but swallowing it silently is what turned
    one bad valence into a crash three stages downstream, in atom typing.
    """
    m = Chem.RWMol(mol)
    try:
        Chem.Kekulize(m, clearAromaticFlags=False)
        return _with_rings_perceived(m)
    except Exception:
        pass

    # Fallback: strip aromaticity entirely. Downstream FF treats the ring as
    # sp3/sp2 single-bonded; MMFF may still refuse without conjugation info,
    # but UFF works and produces reasonable bond lengths.
    for bond in m.GetBonds():
        if bond.GetBondType() == Chem.BondType.AROMATIC:
            bond.SetBondType(Chem.BondType.SINGLE)
            bond.SetIsAromatic(False)
    for atom in m.GetAtoms():
        atom.SetIsAromatic(False)
    try:
        Chem.SanitizeMol(
            m,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
            ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY,
        )
    except Exception:
        pass
    return _with_rings_perceived(m)


def _perpendicular_unit(v: np.ndarray) -> np.ndarray:
    """Return an arbitrary unit vector perpendicular to *v* (v must be non-zero)."""
    v = np.asarray(v, dtype=float)
    # Pick the component with the smallest absolute value to avoid a near-parallel
    # cross product, which would give a numerically unstable result.
    if abs(v[0]) <= abs(v[1]) and abs(v[0]) <= abs(v[2]):
        perp = np.array([0.0, -v[2], v[1]])
    elif abs(v[1]) <= abs(v[2]):
        perp = np.array([-v[2], 0.0, v[0]])
    else:
        perp = np.array([-v[1], v[0], 0.0])
    return perp / np.linalg.norm(perp)


def _optimize_h_positions(
    mol: Chem.Mol,
    all_coords: np.ndarray,
    threshold: float = 1.8,
    num_angles: int = 12,
    allow_out_of_plane: bool = True,
    num_passes: int = 2,
) -> np.ndarray:
    """
    Rotate hydroxyl (O-H) hydrogens about their X-O bond to minimise steric
    clashes, optionally turning them out of the aromatic plane.

    **Scope**: Only H atoms bonded to a heteroatom (O, N, S …) are touched.
    Aromatic C-H bonds are geometrically locked in the sp2 ring plane —
    rotating them out of plane produces an unphysical sp3-like geometry and
    is the wrong fix.  Peri C-H contacts that survive this step are small
    enough (<2 Å) for GROMACS energy minimisation to resolve.

    For each eligible H whose nearest non-bonded neighbour is closer than
    *threshold* Å, the H is swept through *num_angles* equally-spaced dihedral
    angles around the grandparent→parent bond axis (e.g. C–O for phenolic OH),
    which keeps both the bond length and the X–O–H angle.  When
    *allow_out_of_plane* is False only the two in-plane dihedrals are tried.
    The position that maximises the minimum non-bonded distance is kept.

    *num_passes* sweeps are performed (default 2) so that moving H_i updates
    the environment seen by H_j before H_j is evaluated.

    Args:
        mol:               RDKit molecule (connectivity only).
        all_coords:        (N, 3) coordinate array in Å — **modified in place**.
        threshold:         Nearest non-bonded contact (Å) that triggers
                           optimisation.  Default 1.8 Å.
        num_angles:        Dihedral samples (default 12 → every 30°).
        allow_out_of_plane: Sample dihedrals that leave the plane (default True).
        num_passes:        Full sweeps over eligible H atoms (default 2).

    Returns:
        The (possibly modified) coordinate array (same object as *all_coords*).
    """
    excluded = _get_excluded_pairs(mol)

    def _min_nb_dist(h_idx: int, h_pos: np.ndarray) -> float:
        """Minimum distance from *h_pos* to any non-excluded, non-self atom."""
        min_d = float("inf")
        for j in range(mol.GetNumAtoms()):
            if j == h_idx:
                continue
            pair = (min(h_idx, j), max(h_idx, j))
            if pair in excluded:
                continue
            d = float(np.linalg.norm(h_pos - all_coords[j]))
            if d < min_d:
                min_d = d
        return min_d

    angles = np.linspace(0.0, 2.0 * np.pi, num_angles, endpoint=False)
    if not allow_out_of_plane:
        # On a flat layout the in-plane positions are the two dihedrals 0 and π.
        angles = np.array([0.0, np.pi])

    for _pass in range(num_passes):
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 1:
                continue
            h_idx = atom.GetIdx()
            h_pos = all_coords[h_idx].copy()

            # Only optimise H atoms that currently have a short contact
            if _min_nb_dist(h_idx, h_pos) >= threshold:
                continue

            # Parent atom (O, N, C, …)
            nbrs = atom.GetNeighbors()
            if not nbrs:
                continue
            parent = nbrs[0]

            # ── Key constraint ──────────────────────────────────────────────
            # Skip C-H bonds.  In a flat aromatic (sp2) system the C-H bond
            # direction is uniquely determined by the ring geometry; rotating
            # it produces an unphysical geometry and would push the H out of
            # the molecular plane.  Only O-H (and N-H, S-H …) bonds have
            # genuine dihedral freedom around the heteroatom–carbon bond.
            if parent.GetAtomicNum() == 6:
                continue
            # ────────────────────────────────────────────────────────────────

            par_idx = parent.GetIdx()
            par_pos = all_coords[par_idx].copy()

            # Grandparent: first heavy-atom neighbour of parent (excludes self)
            grand_nbrs = [
                n for n in parent.GetNeighbors()
                if n.GetAtomicNum() != 1 and n.GetIdx() != h_idx
            ]
            if not grand_nbrs:
                continue
            grand_pos = all_coords[grand_nbrs[0].GetIdx()].copy()

            # Rotation axis = grandparent → parent (normalised)
            axis = par_pos - grand_pos
            a_norm = float(np.linalg.norm(axis))
            if a_norm < 1e-8:
                continue
            axis = axis / a_norm

            # Decompose h_vec into components along / perpendicular to the axis
            h_vec = h_pos - par_pos
            bl = float(np.linalg.norm(h_vec))
            if bl < 1e-8:
                continue

            h_par = np.dot(h_vec, axis) * axis    # component along rotation axis
            h_perp = h_vec - h_par                 # component in the sweep plane
            perp_norm = float(np.linalg.norm(h_perp))

            # Orthonormal basis in the sweep plane
            u_hat = (h_perp / perp_norm if perp_norm > 1e-8
                     else _perpendicular_unit(axis))
            v_hat = np.cross(axis, u_hat)

            best_pos = h_pos.copy()
            best_d = _min_nb_dist(h_idx, h_pos)

            # Every candidate is a rotation about the axis, so it keeps the
            # bond length and the X-O-H angle (rq-7515b128). Translating a
            # candidate off this cone and rescaling only the bond -- what a
            # z offset did -- bent the angle anywhere from 95 to 178 degrees.
            for theta in angles:
                candidate = par_pos + h_par + perp_norm * (
                    np.cos(theta) * u_hat + np.sin(theta) * v_hat
                )
                d = _min_nb_dist(h_idx, candidate)
                if d > best_d:
                    best_d = d
                    best_pos = candidate

            all_coords[h_idx] = best_pos

    return all_coords


def _substituent_groups(mol: Chem.Mol) -> dict:
    """
    Map each ring atom to the non-ring atoms hanging off it.

    A group is a connected tree of non-ring atoms -- an H, a hydroxyl, a
    carboxyl -- attached to exactly one ring atom, so it can be turned about
    that atom as a rigid body. A chain bridging two ring atoms has no single
    pivot and is left alone.
    """
    in_ring = {a.GetIdx() for a in mol.GetAtoms() if a.IsInRing()}
    seen: set = set()
    groups: dict = {}
    for atom in mol.GetAtoms():
        start = atom.GetIdx()
        if start in in_ring or start in seen:
            continue
        members, stack, anchors = [], [start], set()
        seen.add(start)
        while stack:
            u = stack.pop()
            members.append(u)
            for nbr in mol.GetAtomWithIdx(u).GetNeighbors():
                v = nbr.GetIdx()
                if v in in_ring:
                    anchors.add(v)
                elif v not in seen:
                    seen.add(v)
                    stack.append(v)
        if len(anchors) == 1:
            groups.setdefault(anchors.pop(), []).extend(members)
    return groups


def _ring_normal_at(mol: Chem.Mol, coords: np.ndarray, idx: int) -> Optional[np.ndarray]:
    """Unit normal of the sheet at ring atom *idx*, from two of its ring neighbours."""
    nbrs = [n.GetIdx() for n in mol.GetAtomWithIdx(idx).GetNeighbors() if n.IsInRing()]
    if len(nbrs) < 2:
        return None
    normal = np.cross(coords[nbrs[0]] - coords[idx], coords[nbrs[1]] - coords[idx])
    norm = float(np.linalg.norm(normal))
    return normal / norm if norm > 1e-6 else None


def _rotate_about(points: np.ndarray, origin: np.ndarray, axis: np.ndarray, angle: float) -> np.ndarray:
    """Rotate *points* by *angle* (radians) about the line through *origin* along unit *axis*."""
    v = points - origin
    c, s = np.cos(angle), np.sin(angle)
    return origin + v * c + np.cross(axis, v) * s + np.outer(v @ axis, axis) * (1.0 - c)


def _tilt_crowded_substituents(mol: Chem.Mol, coords: np.ndarray) -> np.ndarray:
    """
    Turn substituents that clash out of the sheet plane, the two sides of a
    crowded site in opposite directions.

    At a cove both substituents point into the same vacant lattice site; tilting
    one up and the other down is the direction a real [4]helicene twists. Each
    group turns rigidly about its ring atom, so its bond lengths and internal
    angles are untouched, by the smallest tilt that clears its contacts (or the
    one that leaves the least overlap, if none does).

    This pass is what makes the force field safe to run afterwards: the
    overlaps it removes are close enough to near-coincident that a force field
    started from them throws atoms apart.
    """
    coords = coords.copy()
    excluded = _get_excluded_pairs(mol)
    floor, judged = _clash_floor_matrix(mol, _get_hbond_pairs(mol, coords), excluded)
    floor = floor - CLASH_SEVERITY_TOLERANCE

    groups = _substituent_groups(mol)
    owner = {a: anchor for anchor, members in groups.items() for a in members}
    origin_pos, tilt_axis = {}, {}
    for anchor, members in groups.items():
        normal = _ring_normal_at(mol, coords, anchor)
        root = [a for a in members if mol.GetBondBetweenAtoms(a, anchor) is not None]
        if normal is None or not root:
            continue
        radial = coords[root[0]] - coords[anchor]
        radial = radial - (radial @ normal) * normal
        axis = np.cross(normal, radial)
        if np.linalg.norm(axis) < 1e-6:
            continue
        origin_pos[anchor] = coords[members].copy()
        tilt_axis[anchor] = axis / np.linalg.norm(axis)

    direction: dict = {}
    angle: dict = {}

    def worst_overlap(atoms):
        d = np.linalg.norm(coords[atoms][:, None, :] - coords[None, :, :], axis=2)
        return float(np.where(judged[atoms], floor[atoms] - d, -np.inf).max())

    def set_tilt(anchor, a):
        coords[groups[anchor]] = _rotate_about(
            origin_pos[anchor], coords[anchor], tilt_axis[anchor], a
        )
        angle[anchor] = a

    steps = np.radians(np.arange(
        LATTICE_TILT_STEP_DEG, LATTICE_TILT_MAX_DEG + 1e-9, LATTICE_TILT_STEP_DEG
    ))
    # A tilt that clears one contact can create another, so sweep a few times.
    for _ in range(4):
        pending = [c for c in _clash_pairs(mol, coords, excluded)
                   if c[0] in owner or c[1] in owner]
        if not pending:
            break
        pending.sort(key=lambda c: c[2] - c[3])  # deepest overlap first
        for i, j, _d, _f in pending:
            movers = list(dict.fromkeys(
                owner[k] for k in (i, j) if k in owner and owner[k] in origin_pos
            ))
            if not movers:
                continue
            atoms = np.array([a for anchor in movers for a in groups[anchor]])
            if worst_overlap(atoms) < 0:
                continue  # already cleared by an earlier move this sweep
            for k, anchor in enumerate(movers):
                if anchor not in direction:
                    partner = [direction[m] for m in movers if m in direction]
                    direction[anchor] = -partner[0] if partner else (1.0 if k == 0 else -1.0)
            # Directions are handed out greedily, deepest contact first, so two
            # groups crowding each other can arrive already sharing one -- both
            # tilt the same way and never separate. Which pairs get there first
            # depends on where every hydroxyl hydrogen happens to point, which
            # made the whole relief hostage to them. Such a pair also tries
            # turning its second group the other way.
            options = [dict(direction)]
            if len(movers) == 2 and direction[movers[0]] == direction[movers[1]]:
                options.append({**direction, movers[1]: -direction[movers[1]]})
            start = {anchor: angle.get(anchor, 0.0) for anchor in movers}
            best = (worst_overlap(atoms), dict(start), direction)
            for option in options:
                for step in steps:
                    for anchor in movers:
                        set_tilt(anchor, option[anchor] * max(step, abs(start[anchor])))
                    overlap = worst_overlap(atoms)
                    if overlap < best[0]:
                        best = (overlap, {anchor: angle[anchor] for anchor in movers}, option)
                    if overlap < 0:
                        break
                if best[0] < 0:
                    break
            direction = best[2]
            for anchor, a in best[1].items():
                set_tilt(anchor, a)
    return coords


def _lattice_force_field_copy(mol: Chem.Mol) -> Chem.Mol:
    """
    A copy of *mol* UFF can parametrise without kekulising: aromatic bonds as
    single bonds, aromatic atoms typed sp2.

    UFF refuses an aromatic sheet it cannot kekulise, and the usual fallback --
    :func:`_kekulize_or_dearomatize` -- leaves ring carbons typed sp3, which
    buckles a flat sheet into a blob. Typing them sp2 keeps UFF's trigonal
    angles and planar inversion terms; the single-bond rest length (~1.46 Å)
    is overridden by distance constraints in the caller.
    """
    working = Chem.RWMol(mol)
    for bond in working.GetBonds():
        if bond.GetIsAromatic() or bond.GetBondType() == Chem.BondType.AROMATIC:
            bond.SetBondType(Chem.BondType.SINGLE)
            bond.SetIsAromatic(False)
    for atom in working.GetAtoms():
        if atom.GetIsAromatic():
            atom.SetIsAromatic(False)
            atom.SetHybridization(Chem.HybridizationType.SP2)
        elif atom.GetAtomicNum() == 8 and all(
            b.GetBondType() == Chem.BondType.SINGLE for b in atom.GetBonds()
        ):
            # RDKit marks a phenolic O conjugated (sp2), which UFF reads as a
            # carbonyl-type O_2 and opens C-O-H to 120°. It is an sp3 oxygen.
            atom.SetHybridization(Chem.HybridizationType.SP3)
        # Every hydrogen is already an explicit atom. Without this, a ring
        # carbon's two single ring bonds read as room for an implicit H, which
        # UFF then warns about.
        atom.SetNoImplicit(True)
    working.UpdatePropertyCache(strict=False)
    return working.GetMol()


def _tethered_lattice_relax(
    mol: Chem.Mol, coords: np.ndarray, seed: Optional[int] = None
) -> Optional[np.ndarray]:
    """
    UFF-minimise *coords* with every ring atom tethered near its position and
    every aromatic bond held to the lattice length. Returns None when UFF
    cannot be set up.

    A perfectly flat start is a saddle a minimiser never leaves -- every
    out-of-plane force is zero -- so ring atoms get a small seeded z offset.
    """
    working = _lattice_force_field_copy(mol)
    start = coords.copy()
    ring = [a.GetIdx() for a in mol.GetAtoms() if a.IsInRing()]
    rng = np.random.default_rng(0 if seed is None else seed)
    start[ring, 2] += rng.normal(0.0, 0.05, len(ring))

    conformer = Chem.Conformer(working.GetNumAtoms())
    for i, xyz in enumerate(start):
        conformer.SetAtomPosition(i, xyz.tolist())
    # The molecule arrives carrying its embedded conformer; UFF uses the
    # first one, so it must go or the pass relaxes the wrong coordinates.
    working.RemoveAllConformers()
    working.AddConformer(conformer, assignId=True)
    try:
        ff = AllChem.UFFGetMoleculeForceField(working)
    except Exception:
        return None
    if ff is None:
        return None
    for i in ring:
        ff.UFFAddPositionConstraint(i, LATTICE_TETHER_RADIUS, LATTICE_TETHER_FORCE_CONSTANT)
    low, high = LATTICE_BOND_RANGE
    for bond in mol.GetBonds():
        if bond.GetIsAromatic():
            ff.UFFAddDistanceConstraint(
                bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), False,
                low, high, LATTICE_BOND_FORCE_CONSTANT,
            )
    ff.Initialize()
    ff.Minimize(maxIts=LATTICE_RELAX_MAX_ITERATIONS)
    relaxed = np.array(working.GetConformer().GetPositions())
    return relaxed if np.all(np.isfinite(relaxed)) else None


class CoordinateGenerator:
    """Generate valid 3D coordinates for molecules."""

    def __init__(self, seed: int = None):
        self.seed = seed
        # Set to True after generate_3d_coordinates() if hex-lattice positions
        # (from carbon_skeleton init_x / init_y props) were used.  When True,
        # callers should skip any subsequent force-field minimization because
        # running MMFF/UFF on the dearomatized mol (all-SINGLE bonds) would
        # incorrectly treat sp2 carbons as sp3 and buckle the flat sheet.
        self.used_hex_lattice: bool = False
        # Carbons the lattice path placed off the lattice (pendant methyls,
        # -CH2- of an aliphatic hydroxyl, carboxyl carbons).  Zero unless
        # used_hex_lattice is True.
        self.offlattice_carbons: int = 0
        if seed is not None:
            np.random.seed(seed)

    def generate_3d_coordinates(
        self,
        mol: Chem.Mol,
        force_aromatic_planarity: bool = True,
        num_conformers: int = 1,
        use_distance_geometry: bool = True,
        max_embedding_attempts: int = 3,
        max_ff_iterations: int = None,
        use_hex_lattice: bool = True,
    ) -> Tuple[Chem.Mol, np.ndarray]:
        """
        Generate 3D coordinates for molecule with enhanced relaxation.

        Args:
            mol: RDKit molecule
            force_aromatic_planarity: If True, enforce planarity on aromatic rings
            num_conformers: Number of conformers to generate
            use_distance_geometry: Use distance geometry for coordinate generation
            max_embedding_attempts: Number of independent embedding attempts (default 3)
            max_ff_iterations: Max force field iterations (default: 500 for large molecules, 300 for small)
            use_hex_lattice: If False, a large molecule takes the Kamada-Kawai
                layout even when it carries hex-lattice positions

        Returns:
            (mol_with_coords, coords_array)
        """
        # Ensure aromaticity is properly perceived before embedding.
        # Do NOT call Kekulize here: it fails for large graphene-like systems and
        # corrupts the aromatic bond types.  ETKDGv3 can embed aromatic bond types
        # natively; UFF/MMFF operate on the sanitized mol after embedding.
        try:
            Chem.SetAromaticity(mol, Chem.AromaticityModel.AROMATICITY_MDL)
        except Exception:
            pass

        # Add hydrogens if not already present
        if mol.GetNumHeavyAtoms() == mol.GetNumAtoms():
            mol = Chem.AddHs(mol)

        # Adaptive force field iterations based on molecule size
        if max_ff_iterations is None:
            num_heavy_atoms = mol.GetNumHeavyAtoms()
            if num_heavy_atoms > 300:
                max_ff_iterations = 500
            elif num_heavy_atoms > 100:
                max_ff_iterations = 400
            else:
                max_ff_iterations = 300

        # For large flat aromatic systems (> 80 heavy atoms), ETKDGv3 frequently
        # folds or collapses the structure.  Use 2D-first embedding instead:
        #   1. Compute 2D layout (flat aromatic ring system)
        #   2. Promote to 3D with z = 0
        #   3. Perturb non-ring atoms slightly in z for FF convergence
        #   4. Run force field minimization
        use_2d_first = mol.GetNumHeavyAtoms() > 80

        # Ring-substituting nitrogen (pyridinic/pyrrolic/graphitic) is best
        # embedded via the flat hex-lattice path: ETKDGv3 tends to buckle the
        # sheet around the substituted ring atom, producing irreducible steric
        # clashes.  When such an N is present AND the skeleton carries
        # hex-lattice positions, use the 2D-first (flat) path regardless of size.
        ring_info = mol.GetRingInfo()
        has_ring_nitrogen = any(
            a.GetAtomicNum() == 7
            and ring_info.NumAtomRings(a.GetIdx()) > 0
            and sum(1 for n in a.GetNeighbors() if n.GetAtomicNum() != 1) >= 2
            for a in mol.GetAtoms()
        )
        if has_ring_nitrogen:
            use_2d_first = True

        best_mol = None
        best_energy = float('inf')
        self.used_hex_lattice = False  # reset each call
        self.offlattice_carbons = 0

        if use_2d_first:
            # Produce a kekulized (or de-aromatized) working copy so that
            # Compute2DCoords, FF setup, etc. don't raise KekulizeException
            # on non-kekulizable biochar graphs (all 100+ atoms flagged).
            mol_copy = Chem.RWMol(_kekulize_or_dearomatize(mol))
            try:
                # 2D-first embedding for large flat aromatics:
                #  1. Compute 2D layout on the heavy-atom-only scaffold
                #  2. Place heavy atoms at z=0 (flat PAH)
                #  3. Compute H positions geometrically (not from 2D conventions)
                #  4. FF minimization to relax H and functional groups
                heavy_atom_indices = [
                    a.GetIdx() for a in mol_copy.GetAtoms() if a.GetAtomicNum() != 1
                ]
                mol_no_h = Chem.RemoveAllHs(mol_copy.GetMol(), sanitize=False)

                # --- Layout strategy (priority order) ---
                # 1. Hex-lattice positions threaded from carbon_skeleton.py
                #    via atom props init_x / init_y.  Gives perfect 1.42 Å
                #    aromatic bonds and planar hex geometry for all sizes.
                # 2. Kamada-Kawai spring layout (fallback when props absent).
                #    Works well up to ~100 atoms; degrades at 200-500 atoms.
                skel_pos = _read_skeleton_positions(mol_no_h)
                use_skel = use_hex_lattice and _lattice_covers_skeleton(mol_no_h, skel_pos)
                if use_skel:
                    self.used_hex_lattice = True
                    self.offlattice_carbons = sum(
                        1 for a in mol_no_h.GetAtoms()
                        if a.GetAtomicNum() == 6 and a.GetIdx() not in skel_pos
                    )

                # Off-lattice atoms' remaining bond directions, keyed by
                # mol_copy index; their hydrogens are placed along these.
                h_slots: dict = {}
                if use_skel:
                    # Ring atoms come directly from hex-lattice geometry.
                    # Everything else -- oxygens, pendant methyls, carboxyl
                    # carbons -- is placed outward from its tagged parent.
                    heavy_pos, h_slots = _offlattice_layout(
                        mol_copy,
                        {heavy_atom_indices[nh]: pos for nh, pos in skel_pos.items()},
                        chemistry=mol,
                    )
                else:
                    # Fallback: kamada-kawai spring layout
                    flat_heavy = _heavy_atom_flat_layout(mol_no_h, seed=self.seed)
                    # Map mol_no_h atom indices (0..Nheavy-1) back to mol_copy
                    # heavy-atom indices (which can interleave with H indices).
                    heavy_pos = {}  # mol_copy_idx -> np.array([x, y, z])
                    for nh_idx, orig_idx in enumerate(heavy_atom_indices):
                        heavy_pos[orig_idx] = flat_heavy[nh_idx].copy()

                # Build full coordinate array; start with heavy atoms
                all_coords = np.zeros((mol_copy.GetNumAtoms(), 3))
                for idx, coord in heavy_pos.items():
                    all_coords[idx] = coord

                # Place H atoms: radiate from parent at correct bond length
                for atom in mol_copy.GetAtoms():
                    if atom.GetAtomicNum() != 1:
                        continue
                    idx = atom.GetIdx()
                    parent = atom.GetNeighbors()[0]
                    par_idx = parent.GetIdx()
                    par_pos = heavy_pos.get(par_idx, all_coords[par_idx])

                    # An off-lattice carbon (a methyl, a -CH2-) is sp3: its
                    # hydrogens take the tetrahedral directions its heavy
                    # neighbours left free.  Radially they would all coincide.
                    if parent.GetAtomicNum() == 6 and h_slots.get(par_idx):
                        all_coords[idx] = par_pos + h_slots[par_idx].pop(0) * 1.09
                        continue

                    # Direction: opposite to average of other neighbor directions
                    other_nbrs = [
                        n.GetIdx() for n in parent.GetNeighbors() if n.GetIdx() != idx
                    ]
                    if other_nbrs:
                        nbr_avg = np.mean(
                            [heavy_pos.get(n, all_coords[n]) for n in other_nbrs], axis=0
                        )
                        h_dir = par_pos - nbr_avg
                    else:
                        h_dir = np.array([1.0, 0.0, 0.0])
                    norm = np.linalg.norm(h_dir)
                    if norm < 1e-8:
                        h_dir = np.array([1.0, 0.0, 0.0])
                    else:
                        h_dir = h_dir / norm

                    # Choose bond length by parent-H element
                    par_elem = parent.GetAtomicNum()
                    bl = 1.09 if par_elem == 6 else 0.96  # C-H or O-H

                    # A hydrogen on a heteroatom with a single heavy neighbour
                    # (hydroxyl, carboxyl OH) has a bond angle to keep, not a
                    # radial direction: "opposite the other neighbours" is the
                    # C-O axis itself, and a collinear X-O-H is both wrong and
                    # a singular torsion for any later force field.
                    if par_elem != 6 and len(other_nbrs) == 1:
                        grand_pos = heavy_pos.get(other_nbrs[0], all_coords[other_nbrs[0]])
                        all_coords[idx] = _bent_hydrogen_position(
                            par_pos, grand_pos, bl, list(heavy_pos.values())
                        )
                        continue
                    all_coords[idx] = par_pos + h_dir * bl

                # Rotate clashing H atoms to minimise non-bonded contacts.
                # Called for both the hex-lattice (use_skel=True) and the
                # Kamada-Kawai fallback paths: both produce a flat z=0 layout,
                # so out-of-plane sampling (z ≠ 0) is safe and effective.
                _optimize_h_positions(mol_copy, all_coords)

                # Create 3D conformer
                conf3d = Chem.Conformer(mol_copy.GetNumAtoms())
                for i, xyz in enumerate(all_coords):
                    conf3d.SetAtomPosition(i, xyz.tolist())
                mol_copy.RemoveAllConformers()
                mol_copy.AddConformer(conf3d, assignId=True)
                best_mol = mol_copy

                # FF minimization:
                # When hex-lattice positions are used (use_skel=True) we SKIP the
                # force field entirely.  The hex lattice gives perfect 1.42 Å CC
                # bonds; running MMFF/UFF on the dearomatized mol (all-SINGLE bonds)
                # would incorrectly treat carbons as sp3 and buckle the sheet into a
                # 3D blob.  H and O positions from geometric placement are good enough
                # for GROMACS input — the production MD run will relax them further.
                #
                # For the KK fallback (use_skel=False) the FF is still useful to
                # refine bond lengths from the spring layout.
                if not use_skel:
                    try:
                        mol_for_ff = mol_copy.GetMol()
                        mmff_props = AllChem.MMFFGetMoleculeProperties(mol_for_ff)
                        ff = (
                            AllChem.MMFFGetMoleculeForceField(mol_for_ff, mmff_props)
                            if mmff_props is not None
                            else None
                        )
                        if ff is None:
                            ff = AllChem.UFFGetMoleculeForceField(mol_for_ff)
                        if ff is not None:
                            ff.Minimize(maxIts=max_ff_iterations, forceTol=1e-4)
                            best_energy = ff.CalcEnergy()
                            # Pull minimized coordinates back into mol_copy
                            conf_min = mol_for_ff.GetConformer(0)
                            conf_target = mol_copy.GetConformer(0)
                            for i in range(mol_for_ff.GetNumAtoms()):
                                conf_target.SetAtomPosition(i, conf_min.GetAtomPosition(i))
                    except Exception:
                        pass
            except Exception:
                best_mol = None
        else:
            # Kekulize (or de-aromatize on failure) once outside the retry loop
            # so EmbedMolecule / MMFF don't throw on the same issue repeatedly.
            mol_safe = _kekulize_or_dearomatize(mol)

            # Embedders are tried best-first: ETKDGv3's torsion knowledge gives
            # the best geometry when it applies, plain distance geometry the
            # worst. But on a substituted fused-PAH sheet the knowledge-based
            # terms are frequently unsatisfiable, and a failure is not cheap --
            # each one burns its full maxIterations budget before returning -1.
            #
            # Measured on 400 C softwood (73 heavy atoms), before this loop
            # remembered anything:
            #
            #   ETKDGv3          18.1s  -> failed
            #   ETKDGv2          18.7s  -> failed
            #   EmbedParameters   0.06s -> succeeded
            #
            # and that ran once per attempt, so generate() spent 110 of its
            # 111 seconds re-discovering the same two failures three times.
            #
            # Only the random seed changes between attempts; the molecule does
            # not. An embedder that exhausted 200 iterations on this graph is
            # not going to succeed on the next seed, so record the failure and
            # stop paying for it. Ordering and first-attempt behaviour are
            # unchanged -- a molecule ETKDGv3 can embed still gets ETKDGv3.
            embedders = [
                ("ETKDGv3", AllChem.ETKDGv3),
                ("ETKDGv2", AllChem.ETKDGv2),
                ("EmbedParameters", AllChem.EmbedParameters),
            ]
            exhausted: set = set()

            for attempt in range(max_embedding_attempts):
                seed = (self.seed + attempt) if self.seed is not None else attempt
                mol_copy = Chem.RWMol(mol_safe)

                result = -1
                for name, params_fn in embedders:
                    if name in exhausted:
                        continue
                    params = params_fn()
                    params.randomSeed = seed
                    params.maxIterations = 200
                    result = AllChem.EmbedMolecule(mol_copy, params)
                    if result == 0:
                        break
                    exhausted.add(name)
                    logger.debug(
                        "%s could not embed this molecule (attempt %d); "
                        "skipping it for the remaining attempts",
                        name, attempt,
                    )

                if result != 0:
                    # Every embedder is exhausted; further attempts would only
                    # re-run the ones that already failed.
                    break

                ff_energy = float('inf')
                # MMFF first — requires MMFFMolProperties, else call throws.
                try:
                    mmff_props = AllChem.MMFFGetMoleculeProperties(mol_copy)
                    if mmff_props is not None:
                        ff = AllChem.MMFFGetMoleculeForceField(mol_copy, mmff_props)
                        if ff is not None:
                            ff.Minimize(maxIts=max_ff_iterations, forceTol=1e-6)
                            ff_energy = ff.CalcEnergy()
                except Exception:
                    pass

                if ff_energy == float('inf'):
                    try:
                        ff = AllChem.UFFGetMoleculeForceField(mol_copy)
                        if ff is not None:
                            ff.Minimize(maxIts=max_ff_iterations, forceTol=1e-6)
                            ff_energy = ff.CalcEnergy()
                    except Exception:
                        pass

                if ff_energy < best_energy:
                    best_energy = ff_energy
                    best_mol = mol_copy

        # Use best embedding, or fall back to 2D if all attempts failed.
        # The embedding is coordinates; the molecule they belong to is the one
        # that came in, not the bond-order-rewritten copy they were computed on.
        if best_mol is not None:
            mol = _adopt_conformer(mol, best_mol)
        else:
            # Final fallback: use 2D coordinates with z=0
            AllChem.Compute2DCoords(mol)

        # Extract coordinates
        try:
            conf = mol.GetConformer(0)
            coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
        except Exception:
            coords = np.zeros((mol.GetNumAtoms(), 3))

        # Apply planarity constraints to aromatic rings if requested
        if force_aromatic_planarity:
            coords = self._enforce_aromatic_planarity(mol, coords)

        return mol, coords

    def _enforce_aromatic_planarity(
        self, mol: Chem.Mol, coords: np.ndarray, tolerance: float = 0.5
    ) -> np.ndarray:
        """
        Enforce planarity on aromatic ring systems.

        Args:
            mol: RDKit molecule
            coords: Current 3D coordinates
            tolerance: Maximum deviation from plane (Angstrom)

        Returns:
            Adjusted coordinates
        """
        # Find aromatic rings
        ring_info = mol.GetRingInfo()
        aromatic_rings = [
            ring for ring in ring_info.AtomRings()
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
        ]

        adjusted_coords = coords.copy()

        for ring_atoms in aromatic_rings:
            if len(ring_atoms) < 3:
                continue

            # Get ring atom coordinates
            ring_coords = coords[list(ring_atoms)]

            # Fit plane to ring atoms
            plane_normal, plane_point = self._fit_plane(ring_coords)

            # Project ring atoms onto plane
            for i, atom_idx in enumerate(ring_atoms):
                point = adjusted_coords[atom_idx]
                projected = point - np.dot(point - plane_point, plane_normal) * plane_normal
                adjusted_coords[atom_idx] = projected

        return adjusted_coords

    @staticmethod
    def _fit_plane(points: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Fit a plane to 3D points using PCA.

        Returns:
            (plane_normal, plane_point)
        """
        center = points.mean(axis=0)
        centered_points = points - center

        # SVD
        u, s, vt = np.linalg.svd(centered_points)
        normal = vt[2]  # Normal is the smallest singular vector

        return normal, center

    def validate_and_relax(
        self, mol: Chem.Mol, coords: np.ndarray, max_iterations: int = 100
    ) -> Tuple[np.ndarray, bool]:
        """
        Validate and relax coordinates using force field.

        Args:
            mol: RDKit molecule
            coords: 3D coordinates
            max_iterations: Maximum relaxation iterations

        Returns:
            (relaxed_coords, converged)

        Only coordinates come back, so the refinement runs on its own
        force-field-safe copy rather than on *mol*. MMFF and UFF both need
        integer bond orders, and the molecule reaching here is the caller's own
        — aromatic bonds included, since embedding no longer hands its
        kekulised working copy back. Parametrising the caller's molecule
        directly would fail on exactly the fused-ring sheets this pass exists
        to relax, and fail quietly: both attempts are guarded, so the
        coordinates would come back untouched and refined-looking.
        """
        working = _kekulize_or_dearomatize(mol)

        conf = Chem.Conformer(working.GetNumAtoms())
        for i, coord in enumerate(coords):
            conf.SetAtomPosition(i, tuple(coord))
        working.RemoveAllConformers()
        working.AddConformer(conf)

        converged = False
        # MMFF94 first. NOTE: MMFFGetMoleculeForceField requires
        # MMFFMolProperties; UFF is the fallback when MMFF has no parameters.
        for build_ff in (_mmff_force_field, AllChem.UFFGetMoleculeForceField):
            try:
                ff = build_ff(working)
            except Exception:
                ff = None
            if ff is None:
                continue
            try:
                converged = ff.Minimize(maxIts=max_iterations) == 0
            except Exception:
                continue
            relaxed = working.GetConformer(0)
            coords = np.array([
                list(relaxed.GetAtomPosition(i))
                for i in range(working.GetNumAtoms())
            ])
            break

        # The caller's molecule carries the coordinates this pass settled on,
        # as it did when the refinement ran on the molecule itself.
        _adopt_conformer(mol, working)
        return coords, converged

    def relieve_lattice_crowding(
        self, mol: Chem.Mol, coords: np.ndarray
    ) -> np.ndarray:
        """
        Relieve the overlaps an exact flat lattice leaves at helicene-type
        edges, and return the new coordinates.

        Two stages. Crowded substituents are first tilted out of plane
        (:func:`_tilt_crowded_substituents`), which removes the near-coincident
        atoms; a tethered UFF relaxation (:func:`_tethered_lattice_relax`) then
        lets each cove open and buckle locally, as it does in a real molecule
        and in GROMACS energy minimisation. The relaxed result is kept only
        when it has no more bond-length errors and no more clashes than the
        tilted one, and the sheet is still flat
        (:data:`~biochar.constants.LATTICE_MAX_SHEET_RMS`); otherwise the
        tilted coordinates are returned. Neither stage can make the geometry
        worse than it arrived.
        """
        coords = self._relieved_coordinates(mol, coords)
        # The caller's molecule carries the coordinates this pass settled on,
        # as it does after validate_and_relax.
        if mol.GetNumConformers():
            conformer = mol.GetConformer(0)
            for i, xyz in enumerate(coords):
                conformer.SetAtomPosition(i, xyz.tolist())
        return coords

    def _relieved_coordinates(self, mol: Chem.Mol, coords: np.ndarray) -> np.ndarray:
        excluded = _get_excluded_pairs(mol)
        if not _clash_pairs(mol, coords, excluded):
            return coords
        tilted = _tilt_crowded_substituents(mol, coords)
        if not _clash_pairs(mol, tilted, excluded):
            return tilted
        relaxed = _tethered_lattice_relax(mol, tilted, self.seed)
        if relaxed is None:
            return tilted

        def defects(c):
            return (
                len(GeometryValidator._validate_bond_lengths(mol, c)),
                len(_clash_pairs(mol, c, excluded)),
            )

        ring = [a.GetIdx() for a in mol.GetAtoms() if a.IsInRing()]
        centred = relaxed[ring] - relaxed[ring].mean(axis=0)
        normal = np.linalg.svd(centred, full_matrices=False)[2][-1]
        sheet_rms = float(np.sqrt(np.mean((centred @ normal) ** 2)))

        bonds_after, clashes_after = defects(relaxed)
        bonds_before, clashes_before = defects(tilted)
        if (
            bonds_after <= bonds_before
            and clashes_after <= clashes_before
            and sheet_rms <= LATTICE_MAX_SHEET_RMS
        ):
            return relaxed
        return tilted


class ClashResolver:
    """Resolve steric clashes in molecular geometries."""

    @staticmethod
    def resolve_clashes(
        mol: Chem.Mol,
        coords: np.ndarray,
        max_iterations: int = 10,
        displacement_step: float = 0.1,
        use_vdw_radii: bool = True,
    ) -> np.ndarray:
        """
        Iteratively resolve steric clashes by displacing atoms away from each other.

        Algorithm:
        - Detect pairs of atoms that are too close (violating Van der Waals radii)
        - Displace the lighter atom (usually H) away from its neighbor
        - Repeat until no clashes remain or max iterations reached

        Args:
            mol: RDKit molecule
            coords: 3D coordinates
            max_iterations: Maximum number of clash resolution iterations
            displacement_step: Step size for atom displacement (Angstrom)
            use_vdw_radii: If True, use Van der Waals radii; else use fixed threshold

        Returns:
            Adjusted coordinates with reduced clashes
        """
        adjusted_coords = coords.copy()

        for iteration in range(max_iterations):
            clashes = ClashResolver._detect_clashes(mol, adjusted_coords, use_vdw_radii)

            if not clashes:
                # No more clashes
                break

            # Displace clashing atoms
            for atom_i, atom_j in clashes:
                # Prefer to move the lighter atom (usually H)
                mass_i = mol.GetAtomWithIdx(atom_i).GetMass()
                mass_j = mol.GetAtomWithIdx(atom_j).GetMass()

                if mass_i < mass_j:
                    moving_atom = atom_i
                else:
                    moving_atom = atom_j

                # Calculate displacement direction (away from the other atom)
                if moving_atom == atom_i:
                    other_atom = atom_j
                else:
                    other_atom = atom_i

                direction = adjusted_coords[moving_atom] - adjusted_coords[other_atom]
                distance = np.linalg.norm(direction)

                if distance > 1e-6:  # Avoid division by zero
                    direction_normalized = direction / distance
                    adjusted_coords[moving_atom] += displacement_step * direction_normalized

        return adjusted_coords

    @staticmethod
    def _detect_clashes(
        mol: Chem.Mol,
        coords: np.ndarray,
        use_vdw_radii: bool = True,
    ) -> List[Tuple[int, int]]:
        """
        Detect steric clashes between atoms.

        1-2 (bonded) and 1-3 (angle) atom pairs are excluded because in flat
        graphene-like sheets the 1-3 C-C distance is ≈2.46 Å — below the
        0.75 × vdW threshold (2.55 Å) — so they would produce hundreds of
        false positives.

        Args:
            mol: RDKit molecule
            coords: 3D coordinates
            use_vdw_radii: If True, use Van der Waals radii; else use fixed 2.0 Å threshold

        Returns:
            List of (atom_i, atom_j) pairs with clashes
        """
        # Pre-compute 1-2 and 1-3 exclusions (see _get_excluded_pairs).
        excluded = _get_excluded_pairs(mol)
        # Hydrogen bonds get a reduced floor here too.  Without this the
        # resolver pushes apart the very H-bonds the MMFF pass re-forms on the
        # next iteration, burning iterations and distorting -OH geometry.
        hbond_pairs = _get_hbond_pairs(mol, coords) if use_vdw_radii else set()
        clashes = []
        distances = squareform(pdist(coords))

        for i in range(mol.GetNumAtoms()):
            for j in range(i + 1, mol.GetNumAtoms()):
                if (i, j) in excluded:
                    continue

                if use_vdw_radii:
                    # A "clash" is real overlap — use 0.75 × vdW-sum as the hard
                    # threshold. Adding +0.2 to the vdW sum (as the old code did)
                    # flagged every peri-H in a PAH as a clash, and the resolver
                    # then destroyed bond geometry trying to fix "clashes" that
                    # were just normal aromatic-ring proximity.
                    min_distance = _clash_floor(mol, i, j, hbond_pairs)
                else:
                    # Use fixed threshold
                    min_distance = 1.5

                # Same tolerance as the validator (see _check_steric_clashes):
                # a contact within CLASH_SEVERITY_TOLERANCE of the floor is not a
                # clash, so the resolver leaves it alone rather than nudging an
                # atom the validator would already accept.
                if distances[i, j] < min_distance - CLASH_SEVERITY_TOLERANCE:
                    clashes.append((i, j))

        return clashes


class GeometryValidator:
    """Validate 3D molecular geometry."""

    @staticmethod
    def validate_geometry(
        mol: Chem.Mol, coords: np.ndarray
    ) -> Tuple[bool, List[str]]:
        """
        Validate 3D geometry for issues.

        Args:
            mol: RDKit molecule
            coords: 3D coordinates

        Returns:
            (is_valid, error_messages)
        """
        errors = []

        # Check for NaN coordinates
        if np.any(np.isnan(coords)):
            errors.append("NaN values found in coordinates")
            return False, errors

        # Check for steric clashes
        clash_errors = GeometryValidator._check_steric_clashes(mol, coords)
        errors.extend(clash_errors)

        # Check bond lengths
        bond_errors = GeometryValidator._validate_bond_lengths(mol, coords)
        errors.extend(bond_errors)

        # Check for unphysical geometries
        if np.max(coords) > 1000 or np.min(coords) < -1000:
            errors.append("Coordinates are unphysically large")

        return len(errors) == 0, errors

    @staticmethod
    def _check_steric_clashes(
        mol: Chem.Mol, coords: np.ndarray, clash_threshold: float = 2.0
    ) -> List[str]:
        """
        Check for steric clashes between atoms with enhanced reporting.

        Args:
            mol: RDKit molecule
            coords: 3D coordinates
            clash_threshold: Minimum allowed distance (Angstrom)

        Returns:
            List of error messages with clash severity
        """
        # Pre-compute 1-2 and 1-3 exclusions: in a flat graphene sheet the
        # 1-3 C-C distance (≈2.46 Å) falls below the 0.75×vdW threshold
        # (2.55 Å), so without this exclusion every angle neighbour pair is
        # reported as a clash.
        excluded = _get_excluded_pairs(mol)
        # Intramolecular O-H...O / N-H...O contacts are hydrogen bonds, not
        # clashes; they get a reduced floor (see _get_hbond_pairs).  Only an
        # overlap deeper than CLASH_SEVERITY_TOLERANCE is reported: a contact
        # within it is embedding noise the FF/EM step relaxes, not a clash.
        hbond_pairs = _get_hbond_pairs(mol, coords)

        errors = []
        for i, j, distance, min_distance in _clash_pairs(mol, coords, excluded):
            symbol_i = mol.GetAtomWithIdx(i).GetSymbol()
            symbol_j = mol.GetAtomWithIdx(j).GetSymbol()
            severity = min_distance - distance  # How far below minimum
            # An H-bonded pair that still violates the reduced floor is
            # genuinely too close — label it so the report is diagnostic.
            clash_type = "H-bond" if (i, j) in hbond_pairs else \
                         "H-H" if symbol_i == "H" and symbol_j == "H" else \
                         "H-C" if symbol_i in ["H", "C"] and symbol_j in ["H", "C"] else \
                         "Other"
            errors.append(
                f"Steric clash: atoms {i} and {j} "
                f"(distance: {distance:.2f}, min: {min_distance:.2f}, "
                f"severity: {severity:.2f} Å, type: {clash_type})"
            )

        return errors

    @staticmethod
    def _validate_bond_lengths(
        mol: Chem.Mol, coords: np.ndarray
    ) -> List[str]:
        """
        Validate bond lengths are reasonable.

        Args:
            mol: RDKit molecule
            coords: 3D coordinates

        Returns:
            List of error messages
        """
        errors = []

        for bond in mol.GetBonds():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()

            atom_i = mol.GetAtomWithIdx(i)
            atom_j = mol.GetAtomWithIdx(j)

            symbol_i = atom_i.GetSymbol()
            symbol_j = atom_j.GetSymbol()

            # Calculate distance
            distance = np.linalg.norm(coords[i] - coords[j])

            # Expected distance: sum of covalent radii, scaled by bond order.
            # The radii are single-bond values, so an unscaled sum over-predicts
            # every aromatic and multiple bond (aromatic C-C would be reported
            # as "expected 1.52" when it is really 1.40).
            r_cov_i = COVALENT_RADII.get(symbol_i, 0.76)
            r_cov_j = COVALENT_RADII.get(symbol_j, 0.76)
            factor = BOND_ORDER_LENGTH_FACTORS.get(str(bond.GetBondType()), 1.00)
            expected_distance = (r_cov_i + r_cov_j) * factor

            min_dist = expected_distance * BOND_LENGTH_MIN_FACTOR
            max_dist = expected_distance * BOND_LENGTH_MAX_FACTOR

            if distance < min_dist or distance > max_dist:
                errors.append(
                    f"Unusual bond length: atoms {i}-{j} "
                    f"(distance: {distance:.2f}, expected: {expected_distance:.2f})"
                )

        return errors

    @staticmethod
    def check_aromaticity(mol: Chem.Mol) -> Tuple[bool, int]:
        """
        Check if aromaticity is preserved.

        Returns:
            (is_valid, num_aromatic_atoms)
        """
        num_aromatic = sum(
            1 for atom in mol.GetAtoms() if atom.GetIsAromatic()
        )
        return num_aromatic > 0, num_aromatic

    @staticmethod
    def measure_ring_planarity(mol: Chem.Mol, coords: np.ndarray) -> Tuple[float, str]:
        """
        Measure planarity of aromatic rings.

        Args:
            mol: RDKit molecule
            coords: 3D coordinates

        Returns:
            (avg_deviation, assessment)
        """
        ring_info = mol.GetRingInfo()
        aromatic_rings = [
            ring for ring in ring_info.AtomRings()
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
        ]

        if not aromatic_rings:
            return 0.0, "No aromatic rings"

        deviations = []

        for ring_atoms in aromatic_rings:
            if len(ring_atoms) < 3:
                continue

            ring_coords = coords[list(ring_atoms)]
            normal, point = CoordinateGenerator._fit_plane(ring_coords)

            # Calculate deviations from plane
            for coord in ring_coords:
                deviation = abs(np.dot(coord - point, normal))
                deviations.append(deviation)

        avg_deviation = np.mean(deviations) if deviations else 0.0
        if avg_deviation < 0.1:
            assessment = "Good planarity"
        elif avg_deviation < 0.3:
            assessment = "Acceptable planarity"
        else:
            assessment = "Poor planarity"

        return avg_deviation, assessment
