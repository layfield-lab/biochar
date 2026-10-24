"""
Carbon Skeleton Generation

Generates aromatic carbon frameworks for biochar structures using PAHs or random graphs.
"""

import math
import random
from typing import List, Tuple, Set, Optional, Dict
from dataclasses import dataclass

import networkx as nx
import logging

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, RWMol

from ..constants import PAH_LIBRARY


# Target aromatic C–C bond length (Å) for hex-lattice layout.  Used to
# rescale seed 2D coords from RDKit (which returns ~1.5 Å bonds) and as the
# edge length for all ring-fusion geometry.
_AROMATIC_CC_BOND = 1.42

# When position-tracking is active, two computed ring-vertex positions that
# lie within this distance (Å) of each other are treated as the same atom.
# This handles shared corners between adjacent hexagons in a fused lattice:
# both growth fronts independently predict the same Cartesian coordinate for
# what is physically one carbon atom.
_POSITION_MERGE_TOLERANCE = 0.05  # Å

# Contact rules for growth with pentagons/heptagons. Odd rings do not tile the
# plane, so once two growth fronts meet around a defect their predicted vertices
# land *near* existing atoms rather than on them. The pure-hexagon rules (merge
# within 0.05 Å, bond anything under 1.45 Å) then bonded those near-misses,
# giving carbons of degree 4-10 that no sanitization accepts -- every Wood-ratio
# sheet of a few hundred carbons failed that way. With defects a fusion is
# therefore planned before it is applied and refused when it would put a new
# atom in an impossible place. Bands, in Å from an existing atom:
#   < merge tolerance       the same atom (wider than the lattice's 0.05, since
#                           off-lattice fronts meet with a few tenths of drift)
#   [_DEFECT_BOND_MIN, _DEFECT_BOND_MAX]  a bond, provided both atoms stay
#                           at degree <= 3
#   anything else under _DEFECT_CLEARANCE  a collision: refuse the fusion
# 2.0 Å sits below the shortest legitimate 1,3 distance a mixed 5/6/7 sheet
# produces (~2.2 Å, two heptagons sharing a vertex) and well above a bond.
_DEFECT_MERGE_TOLERANCE = 0.35
_DEFECT_BOND_MIN = 1.25
_DEFECT_BOND_MAX = 1.60
_DEFECT_CLEARANCE = 2.0

# Attempts per build. A failed growth is seed-dependent -- the RNG breaks ties
# between equally ranked fusion edges and draws ring sizes -- so a fresh
# sub-seed is a real second chance, not a repeat.
_MAX_BUILD_ATTEMPTS = 8

# Curved sheets reject more attempts (a bare bond, a parity that will not fix)
# and each costs ~0.2 s at 600 carbons, so they get more of them.
_MAX_CURVED_BUILD_ATTEMPTS = 24

# Ring count up to which RDKit's exhaustive aromaticity search (every subset of
# up to six fused rings) is cheap enough to allow. Measured on RDKit 2025.09:
# 25 rings 0.06 s, 35 rings 0.5 s, 45 rings 2.5 s, 68 rings did not finish in
# 25 s -- per sanitization, and the pipeline sanitizes many times. 30 rings
# (~100 carbons) keeps the worst case near 0.2 s.
_RDKIT_EXHAUSTIVE_RING_LIMIT = 30

# Best-ranked edges a 5-7-5 cluster is tried on before the heptagon draw is
# given up as a hexagon. Each try copies the sheet for rollback, so the bound
# keeps a blocked draw from costing a pass over the whole perimeter.
_CLUSTER_SITES_TRIED = 6

# Atoms a 5-7-5 cluster adds on open edge (heptagon 5, two pentagons 3 each).
_CLUSTER_ATOMS = 11

# Hexagons a curved sheet may be padded with while looking for a fusion that
# fixes its parity. Each adds four carbons, so this bounds the overshoot.
_PARITY_PADDING_TRIES = 1

logger = logging.getLogger(__name__)


class SkeletonError(RuntimeError):
    """Raised when a carbon skeleton cannot be built for the request.

    Subclasses RuntimeError so callers already catching that keep working.
    Growth depends on the seed (it breaks ties between fusion edges and draws
    ring sizes), and the assembler has already tried several sub-seeds derived
    from it before raising. A different seed can still succeed, so this is not
    a statement that the request is impossible.
    """


def _sub_seed(seed: Optional[int], attempt: int):
    """Seed for retry *attempt* of a build seeded with *seed*.

    Attempt 0 is the caller's seed itself, so a build that succeeds first time
    is unchanged. Later attempts are derived from (seed, attempt) rather than
    ``seed + attempt``: with the additive form, seeds 1 and 2 shared four of
    their five retries, so a "different" seed mostly replayed the same growth.
    ``None`` stays ``None`` (an unseeded build is unseeded on every attempt).
    """
    if seed is None or attempt == 0:
        return seed
    return f"{seed}:{attempt}"


@dataclass
class CarbonSkeleton:
    """Represents a carbon skeleton/framework."""

    mol: Chem.Mol
    smiles: str
    num_carbons: int
    num_aromatic_carbons: int
    aromaticity_percent: float
    ring_composition: Dict[str, int] = None  # {"hexagons": N, "pentagons": M}


# ---------------------------------------------------------------------------
# Ring-fusion geometry helpers (propagate hex-lattice positions during growth)
# ---------------------------------------------------------------------------

def _compute_fused_ring_positions(
    pu: Tuple[float, float],
    pv: Tuple[float, float],
    outward: Tuple[float, float],
    n_sides: int,
    num_new: int,
) -> List[Tuple[float, float]]:
    """
    Compute Cartesian positions for *num_new* new vertices of a regular
    *n_sides*-gon sharing edge (pu, pv), placed on the *outward* side.

    Returns positions in the order matching ``_fuse_hexagon`` /
    ``_fuse_pentagon``:  v → new[0] → new[1] → … → new[num_new-1] → u.
    """
    pu_a = np.array(pu, dtype=float)
    pv_a = np.array(pv, dtype=float)
    edge = pv_a - pu_a
    L = float(np.linalg.norm(edge))
    if L < 1e-8:
        # Degenerate; fall back to a small lattice step along +y
        return [(pu_a[0], pu_a[1] + (i + 1) * _AROMATIC_CC_BOND) for i in range(num_new)]

    e_hat = edge / L
    # Perpendicular (right-handed in 2D)
    n_right = np.array([-e_hat[1], e_hat[0]])
    out = np.array(outward, dtype=float)
    if np.linalg.norm(out) < 1e-8:
        out = n_right
    n_hat = n_right if float(np.dot(n_right, out)) >= 0.0 else -n_right

    apothem = L / (2.0 * math.tan(math.pi / n_sides))
    R = L / (2.0 * math.sin(math.pi / n_sides))
    midpoint = (pu_a + pv_a) / 2.0
    center = midpoint + n_hat * apothem

    theta_v = math.atan2(pv_a[1] - center[1], pv_a[0] - center[0])
    # Choose angular direction (CCW or CW) so traversal wraps AWAY from the
    # shared edge (through the new vertices) rather than across it to u.
    cross = e_hat[0] * n_hat[1] - e_hat[1] * n_hat[0]
    direction = 1.0 if cross > 0 else -1.0

    step = 2.0 * math.pi / n_sides
    positions = []
    for k in range(num_new):
        theta = theta_v + direction * (k + 1) * step
        px = center[0] + R * math.cos(theta)
        py = center[1] + R * math.sin(theta)
        positions.append((float(px), float(py)))
    return positions


def _outward_direction(
    G: nx.Graph,
    positions: Dict[int, Tuple[float, float]],
    u: int,
    v: int,
) -> Tuple[float, float]:
    """
    Estimate the outward normal direction at fusable edge (u, v).

    A fusable edge has both endpoints at degree 2.  Each endpoint has
    exactly one "other neighbour" in the existing graph; those neighbours
    sit on the interior side.  The outward direction points from their
    centroid toward the edge midpoint.
    """
    mid = ((positions[u][0] + positions[v][0]) / 2.0,
           (positions[u][1] + positions[v][1]) / 2.0)
    interior = []
    for nbr in G.neighbors(u):
        if nbr != v and nbr in positions:
            interior.append(positions[nbr])
    for nbr in G.neighbors(v):
        if nbr != u and nbr in positions:
            interior.append(positions[nbr])
    if not interior:
        return (0.0, 1.0)
    cx = sum(p[0] for p in interior) / len(interior)
    cy = sum(p[1] for p in interior) / len(interior)
    return (mid[0] - cx, mid[1] - cy)


# ---------------------------------------------------------------------------
# Utility: convert RDKit Mol (carbon-only) to NetworkX graph
# ---------------------------------------------------------------------------

def _mol_to_graph(mol: Chem.Mol) -> nx.Graph:
    """Convert an RDKit carbon-only Mol to a NetworkX graph."""
    G = nx.Graph()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            G.add_node(atom.GetIdx())
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if i in G.nodes and j in G.nodes:
            G.add_edge(i, j)
    return G


def _graph_to_mol(
    G: nx.Graph,
    positions: Optional[Dict[int, Tuple[float, float]]] = None,
    kekulize_by_matching: bool = False,
) -> Optional[Chem.Mol]:
    """
    Convert a NetworkX carbon graph to an RDKit Mol.

    Uses aromatic bond type so RDKit handles kekulization internally,
    which correctly perceives all fused-ring aromatic carbons.
    Falls back to max-matching Kekule assignment if sanitization fails.

    If *positions* is provided, each carbon atom is tagged with double
    properties ``init_x`` / ``init_y`` so downstream geometry generation
    can use the hex-lattice layout directly (bypassing Compute2DCoords /
    kamada-kawai on large fused graphs).

    *kekulize_by_matching* skips the all-aromatic attempt. RDKit kekulizes by
    backtracking, which on a large sheet with odd rings can run for minutes;
    a maximum matching is polynomial and gives the same Kekule structure class.
    """
    if G.number_of_nodes() == 0:
        return None

    node_list = sorted(G.nodes())
    node_to_idx = {n: i for i, n in enumerate(node_list)}

    def _attach_positions(m: Chem.Mol) -> None:
        if positions is None:
            return
        for orig_node, idx in node_to_idx.items():
            if orig_node in positions:
                x, y = positions[orig_node]
                atom = m.GetAtomWithIdx(idx)
                atom.SetDoubleProp("init_x", float(x))
                atom.SetDoubleProp("init_y", float(y))

    if kekulize_by_matching:
        return _graph_to_mol_by_matching(G, node_list, node_to_idx, _attach_positions)

    # Primary: all-aromatic bonds, let SanitizeMol kekulize
    rwmol = RWMol()
    for _ in node_list:
        a = Chem.Atom(6)
        a.SetNoImplicit(False)
        rwmol.AddAtom(a)
    for u, v in G.edges():
        rwmol.AddBond(node_to_idx[u], node_to_idx[v], Chem.BondType.AROMATIC)
    mol = rwmol.GetMol()
    try:
        Chem.SanitizeMol(mol)
        Chem.SetAromaticity(mol, Chem.AromaticityModel.AROMATICITY_MDL)
        _attach_positions(mol)
        return mol
    except Exception:
        pass

    # Fallback: explicit Kekule via max matching (requires even node count)
    return _graph_to_mol_by_matching(G, node_list, node_to_idx, _attach_positions)


def _graph_to_mol_by_matching(G, node_list, node_to_idx, attach_positions):
    """Kekule structure from a perfect matching; ``None`` if there is none."""
    if G.number_of_nodes() % 2 != 0:
        return None

    matching = nx.max_weight_matching(G, maxcardinality=True)
    double_bonds: Set[Tuple[int, int]] = {
        (min(u, v), max(u, v)) for u, v in matching
    }
    rwmol2 = RWMol()
    for _ in node_list:
        rwmol2.AddAtom(Chem.Atom(6))
    for u, v in G.edges():
        iu, iv = node_to_idx[u], node_to_idx[v]
        key = (min(u, v), max(u, v))
        btype = Chem.BondType.DOUBLE if key in double_bonds else Chem.BondType.SINGLE
        rwmol2.AddBond(iu, iv, btype)
    mol2 = rwmol2.GetMol()
    try:
        if len(matching) * 2 != G.number_of_nodes():
            return None  # no perfect matching: some carbon would be a radical
        Chem.SanitizeMol(mol2)
        Chem.SetAromaticity(mol2, Chem.AromaticityModel.AROMATICITY_MDL)
        attach_positions(mol2)
        return mol2
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Hex-lattice seed builder
# ---------------------------------------------------------------------------

def _build_hex_lattice_graph(num_rings: int) -> nx.Graph:
    """
    Build a compact graphene-like PAH graph from a hexagonal lattice.

    Uses axial hex coordinates. Starts with the centre cell and spirals
    outward.  Each hex cell is a 6-membered ring; shared edges between
    adjacent cells are automatically deduplicated.

    The resulting graph is always bipartite (A/B sublattice of the
    honeycomb), which guarantees a perfect matching exists (by König's
    theorem) when the graph has an even number of nodes.

    Args:
        num_rings: Number of hexagonal cells to place (1 = benzene, 7 = coronene)

    Returns:
        NetworkX graph of the carbon backbone.
    """
    # Generate hex cell centres in axial coordinates using a spiral
    cells = _hex_spiral(num_rings)

    # For each cell, compute its 6 corner positions (using "pointy-top" hex).
    # Corners are shared between adjacent cells.
    # We identify corners by rounding their (x, y) to avoid float duplicates.
    corner_id_map: dict = {}  # (rounded_x, rounded_y) -> node_id
    cell_corners: list = []   # list of lists of node_ids per cell
    next_id = 0

    for q, r in cells:
        # Centre of hex cell in Cartesian (pointy-top orientation)
        cx = math.sqrt(3) * (q + r / 2.0)
        cy = 1.5 * r

        corners = []
        for k in range(6):
            angle = math.pi / 3.0 * k + math.pi / 6.0
            px = cx + math.cos(angle)
            py = cy + math.sin(angle)
            key = (round(px, 4), round(py, 4))
            if key not in corner_id_map:
                corner_id_map[key] = next_id
                next_id += 1
            corners.append(corner_id_map[key])
        cell_corners.append(corners)

    # Build the graph: nodes are corner IDs, edges are ring edges
    G = nx.Graph()
    G.add_nodes_from(range(next_id))
    for corners in cell_corners:
        for i in range(6):
            G.add_edge(corners[i], corners[(i + 1) % 6])

    return G


def _hex_spiral(n: int) -> List[Tuple[int, int]]:
    """Return *n* hex-cell axial coordinates in a compact spiral.

    Cells are emitted from the centre outwards, ring by ring.  Each ring
    is traversed using the correct ring-edge direction vectors so that
    consecutive cells on the same ring are always adjacent (share an
    edge), producing a compact 2-D cluster rather than a strip.
    """
    if n <= 0:
        return []
    cells = [(0, 0)]
    if n == 1:
        return cells

    # Direction vectors for traversing the *perimeter* of ring r.
    # Starting at (r, 0) and going counter-clockwise in axial coords.
    ring_dirs = [(-1, 1), (-1, 0), (0, -1), (1, -1), (1, 0), (0, 1)]
    radius = 1
    while len(cells) < n:
        q, r = radius, 0  # Start of ring `radius`
        for d in range(6):
            for _ in range(radius):
                if len(cells) >= n:
                    return cells
                cells.append((q, r))
                q += ring_dirs[d][0]
                r += ring_dirs[d][1]
        radius += 1
    return cells[:n]


# ---------------------------------------------------------------------------
# Ring-growth engine (shared by PAHAssembler and RandomGraphGenerator)
# ---------------------------------------------------------------------------

def _find_node_at_position(
    positions: Dict[int, Tuple[float, float]],
    target: Tuple[float, float],
    tolerance: float = _POSITION_MERGE_TOLERANCE,
) -> Optional[int]:
    """
    Return the id of a node already at *target* (within *tolerance* Å),
    or ``None`` if no such node exists.

    Called during ring fusion to detect shared corners in the hex lattice:
    when a newly-computed vertex position coincides with an existing node the
    two should be the same atom, not independent duplicates.
    """
    tx, ty = target
    for node_id, (px, py) in positions.items():
        if math.hypot(px - tx, py - ty) < tolerance:
            return node_id
    return None


def _fuse_hexagon(
    G: nx.Graph,
    u: int,
    v: int,
    next_node: int,
    positions: Optional[Dict[int, Tuple[float, float]]] = None,
) -> int:
    """
    Fuse a 6-membered ring onto edge (u, v).

    When *positions* is provided the four ring vertices are located
    analytically **before** any graph modification.  If a computed position
    coincides with an existing node (within ``_POSITION_MERGE_TOLERANCE``)
    that node is **reused** rather than duplicated — correctly modelling the
    shared-corner topology of a graphene lattice when two growth fronts meet.

    Returns the next free node id.
    """
    if positions is not None and u in positions and v in positions:
        # 1. Compute the four vertex positions of the new ring FIRST.
        outward = _outward_direction(G, positions, u, v)
        new_pos_list = _compute_fused_ring_positions(
            positions[u], positions[v], outward, n_sides=6, num_new=4
        )

        # 2. For each vertex: reuse an existing node if the position is
        #    already occupied, otherwise allocate a fresh node.
        new_nodes: List[int] = []
        for pos in new_pos_list:
            existing = _find_node_at_position(positions, pos)
            if existing is not None:
                new_nodes.append(existing)
            else:
                G.add_node(next_node)
                positions[next_node] = pos
                new_nodes.append(next_node)
                next_node += 1

        # 3. Add the ring edges (nx.Graph silently ignores duplicate edges).
        ring_seq = [v, new_nodes[0], new_nodes[1], new_nodes[2], new_nodes[3], u]
        for a, b in zip(ring_seq, ring_seq[1:]):
            G.add_edge(a, b)

        return next_node

    else:
        # No position tracking: original behaviour (always allocate 4 new nodes).
        new_nodes = list(range(next_node, next_node + 4))
        for n in new_nodes:
            G.add_node(n)
        G.add_edge(v, new_nodes[0])
        G.add_edge(new_nodes[0], new_nodes[1])
        G.add_edge(new_nodes[1], new_nodes[2])
        G.add_edge(new_nodes[2], new_nodes[3])
        G.add_edge(new_nodes[3], u)
        return next_node + 4


def _fuse_pentagon(
    G: nx.Graph,
    u: int,
    v: int,
    next_node: int,
    positions: Optional[Dict[int, Tuple[float, float]]] = None,
) -> int:
    """
    Fuse a 5-membered ring onto edge (u, v).

    Applies the same position-merge logic as ``_fuse_hexagon``: computed
    vertices that coincide with existing nodes are reused rather than
    duplicated.

    Returns the next free node id.
    """
    if positions is not None and u in positions and v in positions:
        outward = _outward_direction(G, positions, u, v)
        new_pos_list = _compute_fused_ring_positions(
            positions[u], positions[v], outward, n_sides=5, num_new=3
        )

        new_nodes: List[int] = []
        for pos in new_pos_list:
            existing = _find_node_at_position(positions, pos)
            if existing is not None:
                new_nodes.append(existing)
            else:
                G.add_node(next_node)
                positions[next_node] = pos
                new_nodes.append(next_node)
                next_node += 1

        ring_seq = [v, new_nodes[0], new_nodes[1], new_nodes[2], u]
        for a, b in zip(ring_seq, ring_seq[1:]):
            G.add_edge(a, b)

        return next_node

    else:
        new_nodes = list(range(next_node, next_node + 3))
        for n in new_nodes:
            G.add_node(n)
        G.add_edge(v, new_nodes[0])
        G.add_edge(new_nodes[0], new_nodes[1])
        G.add_edge(new_nodes[1], new_nodes[2])
        G.add_edge(new_nodes[2], u)
        return next_node + 3


def _fuse_heptagon(
    G: nx.Graph,
    u: int,
    v: int,
    next_node: int,
    positions: Optional[Dict[int, Tuple[float, float]]] = None,
) -> int:
    """
    Fuse a 7-membered ring onto edge (u, v) — adds five new vertices.

    Heptagons introduce the negative (saddle) curvature that, together with
    pentagons (positive curvature), disrupts graphitic stacking in
    non-graphitizing carbons (Wood et al. 2024). Same position-merge logic as
    ``_fuse_hexagon`` / ``_fuse_pentagon``.

    Returns the next free node id.
    """
    if positions is not None and u in positions and v in positions:
        outward = _outward_direction(G, positions, u, v)
        new_pos_list = _compute_fused_ring_positions(
            positions[u], positions[v], outward, n_sides=7, num_new=5
        )

        new_nodes: List[int] = []
        for pos in new_pos_list:
            existing = _find_node_at_position(positions, pos)
            if existing is not None:
                new_nodes.append(existing)
            else:
                G.add_node(next_node)
                positions[next_node] = pos
                new_nodes.append(next_node)
                next_node += 1

        ring_seq = [v] + new_nodes + [u]
        for a, b in zip(ring_seq, ring_seq[1:]):
            G.add_edge(a, b)

        return next_node

    else:
        new_nodes = list(range(next_node, next_node + 5))
        for n in new_nodes:
            G.add_node(n)
        ring_seq = [v] + new_nodes + [u]
        for a, b in zip(ring_seq, ring_seq[1:]):
            G.add_edge(a, b)
        return next_node + 5


def _plan_defect_fusion(
    G: nx.Graph,
    positions: Dict[int, Tuple[float, float]],
    u: int,
    v: int,
    n_sides: int,
) -> Optional[Tuple[List[Tuple[Optional[int], Tuple[float, float]]], List[Tuple[int, int]]]]:
    """Plan fusing an *n_sides* ring onto (u, v) under the defect contact rules.

    Returns ``(vertices, contacts)`` -- one ``(existing_node_or_None, position)``
    per new ring vertex, and the extra bonds the ring closes onto atoms it
    meets -- or ``None`` when the ring would collide with the sheet or push any
    carbon past three bonds. Nothing is modified.
    """
    outward = _outward_direction(G, positions, u, v)
    new_pos = _compute_fused_ring_positions(
        positions[u], positions[v], outward, n_sides=n_sides, num_new=n_sides - 2
    )

    vertices: List[Tuple[Optional[int], Tuple[float, float]]] = []
    pending_contacts: List[Tuple[int, int]] = []  # (vertex index, existing node)
    for i, (px, py) in enumerate(new_pos):
        merged = None
        bonded: List[int] = []
        for node, (qx, qy) in positions.items():
            d = math.hypot(px - qx, py - qy)
            if d < _DEFECT_MERGE_TOLERANCE:
                if merged is not None:
                    return None
                merged = node
            elif d < _DEFECT_BOND_MIN:
                return None
            elif d <= _DEFECT_BOND_MAX:
                bonded.append(node)
            elif d < _DEFECT_CLEARANCE:
                return None
        if merged is not None:
            if merged in (u, v) or any(m == merged for m, _ in vertices):
                return None
            # A merged corner's neighborhood is already part of the sheet;
            # the distances above were measured from a point a few tenths away
            # from it, so its own bonds and neighbors are not new contacts.
            vertices.append((merged, (px, py)))
        else:
            vertices.append((None, (px, py)))
            pending_contacts.extend((i, b) for b in bonded)

    # A contact closes a second ring through the sheet; refuse it when that
    # ring would be smaller than a pentagon. Its size is one (the contact) plus
    # the way round the new ring from vertex i to an existing ring atom x plus
    # the path from x back to the contacted atom b.
    n = n_sides
    existing_at = {0: v, n - 1: u}  # ring order: v, new..., u
    for k, (m, _) in enumerate(vertices):
        if m is not None:
            existing_at[k + 1] = m
    ring_atoms = set(existing_at.values())
    for i, b in pending_contacts:
        if b in ring_atoms:
            continue  # a bond of the new ring itself, not a contact
        dist_b = nx.single_source_shortest_path_length(G, b, cutoff=3)
        for k, x in existing_at.items():
            around = abs((i + 1) - k)
            around = min(around, n - around)
            if x in dist_b and 1 + around + dist_b[x] < 5:
                return None

    # Degree check over every edge the fusion would add. Placeholder ids
    # ("new", i) stand in for vertices that do not exist yet.
    ids = [m if m is not None else ("new", i) for i, (m, _) in enumerate(vertices)]
    ring_seq = [v] + ids + [u]
    ring_edges = {frozenset(p) for p in zip(ring_seq, ring_seq[1:])}
    contact_edges = {frozenset((ids[i], b)) for i, b in pending_contacts} - ring_edges
    added = [
        e for e in ring_edges | contact_edges
        if not all(isinstance(n, int) for n in e) or not G.has_edge(*e)
    ]
    extra: Dict = {}
    for e in added:
        for n in e:
            extra[n] = extra.get(n, 0) + 1
    for n, k in extra.items():
        base = G.degree(n) if isinstance(n, int) else 0
        if base + k > 3:
            return None

    contacts = [(i, b) for i, b in pending_contacts
                if frozenset((ids[i], b)) in contact_edges]
    return vertices, contacts


def _apply_defect_fusion(
    G: nx.Graph,
    positions: Dict[int, Tuple[float, float]],
    u: int,
    v: int,
    plan,
    next_node: int,
) -> Tuple[int, Tuple[int, ...]]:
    """Apply a plan from ``_plan_defect_fusion``.

    Returns the next free node id and the new ring's atoms.

    Contact bonds are added immediately rather than in a pass at the end, so
    the degrees the next fusion is planned against are the real ones.
    """
    vertices, contacts = plan
    ids: List[int] = []
    for existing, pos in vertices:
        if existing is not None:
            ids.append(existing)
        else:
            G.add_node(next_node)
            positions[next_node] = pos
            ids.append(next_node)
            next_node += 1
    ring_seq = [v] + ids + [u]
    for a, b in zip(ring_seq, ring_seq[1:]):
        G.add_edge(a, b)
    for i, b in contacts:
        G.add_edge(ids[i], b)
    return next_node, tuple(ring_seq)


def _first_clear_fusion(
    G: nx.Graph,
    positions: Dict[int, Tuple[float, float]],
    ranked_edges: List[Tuple[int, int]],
    ring_preference: Tuple[int, ...],
    odd_only: bool = False,
):
    """First (u, v, ring, plan) that passes the defect contact rules.

    Ring sizes are tried in *ring_preference* order, each across all edges in
    rank order, so the drawn size is kept wherever it fits before falling back.
    With *odd_only*, only fusions adding an odd number of atoms qualify (the
    parity step).
    """
    seen = set()
    for ring in ring_preference:
        if ring in seen:
            continue
        seen.add(ring)
        for u, v in ranked_edges:
            plan = _plan_defect_fusion(G, positions, u, v, ring)
            if plan is None:
                continue
            if odd_only and sum(1 for m, _ in plan[0] if m is None) % 2 == 0:
                continue
            return u, v, ring, plan
    return None


def _hexagon_new_atom_count(
    G: nx.Graph,
    positions: Dict[int, Tuple[float, float]],
    u: int,
    v: int,
) -> int:
    """How many atoms ``_fuse_hexagon`` would add at (u, v); nothing modified."""
    outward = _outward_direction(G, positions, u, v)
    new_pos = _compute_fused_ring_positions(
        positions[u], positions[v], outward, n_sides=6, num_new=4
    )
    return sum(1 for p in new_pos if _find_node_at_position(positions, p) is None)


def _choose_ring_size(
    rng: random.Random,
    pentagon_fraction: float,
    heptagon_fraction: float,
    nodes_remaining: int,
) -> int:
    """Pick a ring size (5, 6, or 7) for the next fusion.

    Pentagons/heptagons are drawn with their requested probabilities; hexagons
    are the remainder. A ring is only eligible if it does not overshoot the
    remaining node budget by more than one (pentagon +3, hexagon +4,
    heptagon +5), so growth lands near the target size.
    """
    penta_ok = nodes_remaining >= 3
    hexa_ok = nodes_remaining >= 3          # hexagon is the fallback; allow near the end
    hepta_ok = nodes_remaining >= 5

    r = rng.random()
    if penta_ok and r < pentagon_fraction:
        return 5
    if hepta_ok and r < pentagon_fraction + heptagon_fraction:
        return 7
    if hexa_ok:
        return 6
    return 5  # tiny budget: a pentagon overshoots least


def _edge_carbon_fraction(mol: Chem.Mol) -> float:
    """
    Fraction of carbons on the perimeter (heavy-atom degree < 3).

    On a fully H-saturated flake every such carbon carries exactly one
    hydrogen, so this fraction *is* the maximum achievable H/C for the
    structure (its "H/C ceiling").  Computed on the hydrogen-free skeleton.
    """
    carbons = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 6]
    if not carbons:
        return 0.0
    edge = sum(1 for a in carbons if a.GetDegree() < 3)
    return edge / len(carbons)


def _attach_seed_positions(
    mol: Chem.Mol, positions: Optional[Dict[int, Tuple[float, float]]]
) -> None:
    """Attach hex-lattice ``init_x``/``init_y`` props to *mol* in place."""
    if not positions:
        return
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx in positions:
            x, y = positions[idx]
            atom.SetDoubleProp("init_x", x)
            atom.SetDoubleProp("init_y", y)


def _compactness_for_hc(target_h_c: Optional[float], compact_ceiling: float) -> float:
    """
    Map a requested H/C ratio to a growth ``compactness`` in [0, 1], given the
    *measured* H/C ceiling of the maximally condensed flake for this size.

    Returns 1.0 (fully compact) when the target is at or below the compact
    ceiling.  Above it, compactness drops toward 0 (elongated, higher-perimeter
    growth) in proportion to how far above the ceiling the target sits.
    """
    if target_h_c is None or compact_ceiling <= 0:
        return 1.0
    if target_h_c <= compact_ceiling:
        return 1.0
    excess = (target_h_c - compact_ceiling) / compact_ceiling
    return max(0.0, min(1.0, 1.0 - excess))


def _grow_graph(
    G: nx.Graph,
    target_nodes: int,
    seed: Optional[int] = None,
    defect_fraction: float = 0.0,
    positions: Optional[Dict[int, Tuple[float, float]]] = None,
    compactness: float = 1.0,
    heptagon_fraction: float = 0.0,
) -> nx.Graph:
    """
    Grow a carbon graph by iteratively fusing rings.

    By default only 6-membered rings are added (+4 nodes each).  When
    *defect_fraction* > 0 the algorithm occasionally fuses a 5-membered
    ring (+3 nodes) instead, introducing pentagonal defects analogous to
    those observed in disordered graphitic carbon.

    Each new ring shares one edge with the existing graph.  The algorithm
    prefers edges whose neighbours have higher total degree, promoting
    compact 2D growth.

    A parity guard ensures the final graph has an even number of nodes
    (required for a valid Kekulé / perfect-matching assignment):
      - Pure hexagon mode (defect_fraction == 0): adds one extra hexagon.
      - Defect mode (defect_fraction >  0): adds one extra pentagon
        because odd + 3 = even, whereas odd + 4 remains odd.

    Args:
        G: Seed graph (will NOT be mutated; a copy is made).
        target_nodes: Desired number of nodes.
        seed: Optional RNG seed for tie-breaking and pentagon selection.
        defect_fraction: Probability [0, 1) that any given ring addition
            is a 5-membered (pentagon) ring rather than a 6-membered ring.

    Returns:
        Grown graph (new object).
    """
    G = G.copy()
    rng = random.Random(seed)
    next_node = max(G.nodes()) + 1 if G.nodes() else 0
    # Defect growth with positions is planned against the contact rules (see
    # _DEFECT_CLEARANCE). Pure-hexagon growth keeps its historical path: on an
    # exact lattice every predicted vertex is a lattice site, so there are no
    # near-misses to guard against.
    guarded = positions is not None and (defect_fraction > 0 or heptagon_fraction > 0)
    if guarded:
        return _grow_curved(
            G, target_nodes, rng, defect_fraction, heptagon_fraction,
            positions, next_node,
        )

    # When position-tracking is active some ring fusions add zero new nodes
    # (shared corners are merged into existing nodes).  We must therefore
    # continue the loop until the *unique* node count reaches target_nodes,
    # not merely until we've attempted target/min_step iterations.
    #
    # For the common no-merging case the two conditions are equivalent:
    # the seed is chosen so (target - seed) % 4 == 0, each step adds
    # exactly 4 nodes, and we reach target_nodes on the last iteration
    # regardless of whether the guard reads  "< target"  or
    # "+ 4 <= target".

    while G.number_of_nodes() < target_nodes:
        fusable = _get_fusable_edges(G)
        if not fusable:
            break

        # Edge selection controls the perimeter/area ratio of the flake, which
        # in turn sets the maximum achievable H/C (edge carbons carry the H).
        #   compactness == 1.0 -> always fuse at the highest-degree (interior)
        #     edge: maximally condensed flake, lowest H/C ceiling (default,
        #     byte-identical to the historical behaviour).
        #   compactness  < 1.0 -> with probability (1 - compactness) fuse at a
        #     low-degree (peripheral / tip) edge instead, elongating growth,
        #     raising the perimeter and thus the reachable H/C.
        if compactness >= 1.0:
            fusable.sort(key=lambda e: (-_edge_neighbor_degree(G, e[0], e[1]), rng.random()))
        elif rng.random() < (1.0 - compactness):
            fusable.sort(key=lambda e: (_edge_neighbor_degree(G, e[0], e[1]), rng.random()))
        else:
            fusable.sort(key=lambda e: (-_edge_neighbor_degree(G, e[0], e[1]), rng.random()))
        u, v = fusable[0]

        nodes_remaining = target_nodes - G.number_of_nodes()

        # Decide ring size for this step. Pentagon (+3) / heptagon (+5) are drawn
        # with their requested fractions; hexagon (+4) is the remainder. When no
        # defects are requested this is always a hexagon.
        if defect_fraction > 0 or heptagon_fraction > 0:
            ring = _choose_ring_size(rng, defect_fraction, heptagon_fraction, nodes_remaining)
        else:
            ring = 6

        if ring == 5:
            next_node = _fuse_pentagon(G, u, v, next_node, positions=positions)
        elif ring == 7:
            next_node = _fuse_heptagon(G, u, v, next_node, positions=positions)
        else:
            next_node = _fuse_hexagon(G, u, v, next_node, positions=positions)

    # When position-tracking was active, some ring fusions may have created
    # non-adjacent pairs at bond distance (1.42 Å) that were never directly
    # bonded (see ``_add_missing_lattice_bonds`` docstring).  Fix these before
    # building the mol so the graph has the correct topology.  The 1.45 Å
    # tolerance safely distinguishes genuine bonds from all non-bonded C-C
    # pairs (min non-bonded = 2.30 Å in pentagons, 2.46 Å in hexagons).
    if positions is not None:
        _add_missing_lattice_bonds(G, positions)

    if positions is not None and defect_fraction <= 0 and heptagon_fraction <= 0:
        # Pure-hexagon parity. A hexagon adds four atoms on open edge and so
        # cannot change parity; only one that closes onto an existing corner
        # (a cove: three new atoms) can. The generic loop below fused the
        # top-ranked edge whatever it added, so an odd sheet kept growing
        # hexagon by hexagon until a cove happened to come up first -- a 104-C
        # request grew to 463 C and was then unkekulizable anyway. Fuse at the
        # best-ranked edge that adds an odd count, or leave the sheet odd for
        # the caller's retry if none does.
        if G.number_of_nodes() % 2 != 0:
            fusable = _get_fusable_edges(G)
            fusable.sort(key=lambda e: (-_edge_neighbor_degree(G, e[0], e[1]), rng.random()))
            for u, v in fusable:
                if _hexagon_new_atom_count(G, positions, u, v) % 2 == 1:
                    _fuse_hexagon(G, u, v, next_node, positions=positions)
                    _add_missing_lattice_bonds(G, positions)
                    break
        return G

    # Parity guard without positions: ensure even node count for Kekulé
    # assignment. Without position merging every fusion adds its full count,
    # so a pentagon or heptagon (odd) restores parity in one step.
    _parity_limit = max(10, G.number_of_nodes())
    for _ in range(_parity_limit):
        if G.number_of_nodes() % 2 == 0:
            break
        fusable = _get_fusable_edges(G)
        if not fusable:
            break
        fusable.sort(key=lambda e: (-_edge_neighbor_degree(G, e[0], e[1]), rng.random()))
        u, v = fusable[0]
        if defect_fraction > 0:
            # Pentagon adds odd → odd + odd = even  ✓
            next_node = _fuse_pentagon(G, u, v, next_node, positions=positions)
        elif heptagon_fraction > 0:
            # Heptagon also adds odd → restores parity when only heptagons are on
            next_node = _fuse_heptagon(G, u, v, next_node, positions=positions)
        else:
            # Hexagon adds even → pure-hexagon safety net
            next_node = _fuse_hexagon(G, u, v, next_node, positions=positions)

    return G


def _grow_curved(
    G: nx.Graph,
    target_nodes: int,
    rng: random.Random,
    defect_fraction: float,
    heptagon_fraction: float,
    positions: Dict[int, Tuple[float, float]],
    next_node: int,
) -> nx.Graph:
    """Ring-fusion growth with pentagons and heptagons (see ``_grow_graph``).

    Specified in rqm/carbon-skeleton.md (rq-2d15a0c0), including its known
    limit: a flat projection cannot close around a cluster's disclination, so
    curved sheets carry about twice the edge of flat ones (rq-37961357).

    Two things differ from flat hexagon growth.

    Every fusion is planned against the contact rules first (``_DEFECT_*``), so
    fronts that meet around a defect never bond atoms that are not really
    neighbors.

    Defects are placed in fused 5-7 clusters where possible. RDKit perceives
    aromaticity ring by ring, then searches *subsets* of up to six fused rings
    until every ring bond is aromatic -- and that search is combinatorial in
    the number of rings in the sheet. Pentagons and heptagons are aromatic only
    in 5-7 pairs, and a pair leaves its own shared bond bare
    (``_uncovered_defect_bonds``), so a defect needs two complementary
    neighbors: a heptagon carrying two pentagons. A sheet with a bare bond
    makes every sanitization search all subsets -- the pipeline sanitizes many
    times -- which took 0.5 s at 104 C and did not finish at 200 C. Wood et
    al.'s curved sheets cluster their defects this way (their heptagons carry
    two or three fused pentagons), which is why their 256-carbon block
    perceives in under half a second.

    Growth cannot guarantee coverage -- holes and contact-closed rings appear
    that it never planned -- so callers measure the finished sheet with
    ``_aromaticity_search_is_bounded`` and treat a failure as a failed attempt.
    """
    rings: List[Tuple[int, ...]] = []  # defect rings added, as ordered atom tuples
    # Bonds inside some hexagon. Library seeds are all-hexagon, so every seed
    # bond starts covered.
    hex_bonds = {frozenset(e) for e in G.edges()}

    def fuse(choice):
        nonlocal next_node
        u, v, ring, plan = choice
        next_node, atoms = _apply_defect_fusion(G, positions, u, v, plan, next_node)
        if ring == 6:
            hex_bonds.update(_ring_bonds(atoms))
        else:
            rings.append(atoms)
        return atoms

    def ranked(edges):
        return sorted(edges, key=lambda e: (-_edge_neighbor_degree(G, e[0], e[1]), rng.random()))

    def on_ring(ring_atoms):
        rs = set(ring_atoms)
        return [e for e in _get_fusable_edges(G) if e[0] in rs and e[1] in rs]

    def attach_pentagon(heptagons, odd_only=False):
        """Fuse a pentagon onto the first heptagon in *heptagons* that takes one."""
        for h in heptagons:
            choice = _first_clear_fusion(G, positions, ranked(on_ring(h)), (5,), odd_only=odd_only)
            if choice is not None:
                fuse(choice)
                return True
        return False

    def place_cluster(edges):
        """Fuse a heptagon carrying two pentagons; roll back if it cannot.

        A heptagon is placed only with two pentagons on it, which covers every
        bond of the cluster; with one, the shared 5-7 bond stays bare. So the
        cluster is built as a unit and undone if the second pentagon has
        nowhere to go.
        """
        nonlocal next_node
        for u, v in edges[:_CLUSTER_SITES_TRIED]:
            plan = _plan_defect_fusion(G, positions, u, v, 7)
            if plan is None:
                continue
            saved = (G.copy(), dict(positions), next_node, len(rings), set(hex_bonds))
            heptagon = fuse((u, v, 7, plan))
            if attach_pentagon([heptagon]) and attach_pentagon([heptagon]):
                return True
            G_saved, pos_saved, next_node, n_rings, hex_saved = saved
            G.clear()
            G.update(G_saved)
            positions.clear()
            positions.update(pos_saved)
            del rings[n_rings:]
            hex_bonds.clear()
            hex_bonds.update(hex_saved)
        return False

    # Pentagon draws waiting for a heptagon to sit on. With heptagons requested,
    # a pentagon is placed only on a heptagon: alone it can never be covered.
    # Only on sheets too big for RDKit's exhaustive search, though -- a compact
    # flake has about one ring per three carbons -- since below that a lone
    # defect costs nothing and clusters (11 atoms each) would crowd a small
    # sheet out of heptagons altogether.
    pending_pentagons = 0
    clustering = (heptagon_fraction > 0
                  and target_nodes > 3 * _RDKIT_EXHAUSTIVE_RING_LIMIT)

    while G.number_of_nodes() < target_nodes:
        fusable = _get_fusable_edges(G)
        if not fusable:
            break
        fusable = ranked(fusable)
        nodes_remaining = target_nodes - G.number_of_nodes()

        # A cluster adds 11 atoms and so flips parity; hexagons do not. Fixing
        # an odd count after the target is reached costs the most overshoot,
        # so an odd sheet takes its last cluster while it still fits.
        if (clustering and G.number_of_nodes() % 2 == 1
                and _CLUSTER_ATOMS <= nodes_remaining < _CLUSTER_ATOMS + 4
                and place_cluster(fusable)):
            pending_pentagons = max(0, pending_pentagons - 2)
            continue

        ring = _choose_ring_size(rng, defect_fraction, heptagon_fraction, nodes_remaining)

        if clustering and ring == 5:
            # On an existing heptagon if one has room (Wood's heptagons carry
            # two or three); otherwise owed to the next heptagon.
            if attach_pentagon([r for r in rings if len(r) == 7]):
                continue
            pending_pentagons += 1
            ring = 6

        if clustering and ring == 7:
            # A cluster adds 11 atoms; near the target it would overshoot by
            # more than a ring, so the draw becomes a hexagon there.
            if nodes_remaining >= _CLUSTER_ATOMS and place_cluster(fusable):
                pending_pentagons = max(0, pending_pentagons - 2)
                continue
            ring = 6

        # A hexagon, or a defect placed alone (small sheets, or pentagon-only
        # growth). Take the best-ranked edge that can carry it without a
        # collision, falling back to a hexagon. A drawn defect that fits nowhere is a lost draw, not a
        # stall -- the realized share already sits below the probability.
        choice = _first_clear_fusion(G, positions, fusable, (ring, 6))
        if choice is None:
            break  # every front is blocked
        fuse(choice)

    # Parity needs a fusion adding an odd number of atoms. Cheapest first: a
    # hexagon closing a cove (+3, no defect), then a pentagon on an existing
    # heptagon (+3, stays covered) -- or, on a sheet small enough not to
    # cluster, a lone pentagon or heptagon. With none available, pad with one
    # hexagon -- which opens new coves -- and look again. A whole cluster
    # (+11) would also do, but overshoots by more than a ring or two; a sheet
    # still odd is left for the caller to retry with another sub-seed.
    for attempt in range(_PARITY_PADDING_TRIES + 1):
        if G.number_of_nodes() % 2 == 0:
            break
        edges = ranked(_get_fusable_edges(G))
        choice = _first_clear_fusion(G, positions, edges, (6,), odd_only=True)
        if choice is None and not clustering:
            choice = _first_clear_fusion(G, positions, edges, (5, 7), odd_only=True)
        if choice is not None:
            fuse(choice)
            break
        if clustering and attach_pentagon([r for r in rings if len(r) == 7], odd_only=True):
            break
        if attempt == _PARITY_PADDING_TRIES:
            break  # still odd: the caller retries rather than overshoot further
        choice = _first_clear_fusion(G, positions, edges, (6,))
        if choice is None:
            break
        fuse(choice)

    return G


def _ring_bonds(ring) -> Set[frozenset]:
    """Bonds of a ring given as an ordered atom tuple."""
    return {frozenset((ring[i], ring[(i + 1) % len(ring)])) for i in range(len(ring))}


def _uncovered_defect_bonds(defect_rings, hex_bonds: Set[frozenset]) -> Set[frozenset]:
    """Defect-ring bonds RDKit's cheap aromaticity passes will not mark.

    Mirrors RDKit's ``applyHuckelToFused`` (2025.09), which stops searching as
    soon as every ring bond is aromatic. A hexagon of sp2 carbons is aromatic
    alone and marks all six bonds. A pentagon (5 pi) or heptagon (7) is never
    aromatic alone; a 5-7 pair sharing one bond is (10, azulene), but RDKit
    marks only the pair's *envelope* -- the bonds in exactly one of the two
    rings -- so the shared bond needs a second pair to cover it: a heptagon
    carrying two pentagons, or a pentagon between two heptagons. No other pair
    of 5/6/7 rings is 4n+2, so anything left over is reached only by subsets of
    three to six rings, enumerated over every ring in the sheet.
    """
    covered = set(hex_bonds)
    bonds = [_ring_bonds(r) for r in defect_rings]
    for ri, rb in zip(defect_rings, bonds):
        # Any all-sp2 ring of 4n+2 atoms is aromatic alone (a 10- or 14-atom
        # hole), up to RDKit's 24-atom ring limit.
        if len(ri) % 4 == 2 and len(ri) <= 24:
            covered |= rb
    for i, ri in enumerate(defect_rings):
        for j in range(i + 1, len(defect_rings)):
            if len(ri) + len(defect_rings[j]) != 12 or len(ri) == 6:
                continue
            shared = bonds[i] & bonds[j]
            if len(shared) == 1:  # RDKit does not pair rings sharing more
                covered |= (bonds[i] | bonds[j]) - shared
    return set().union(*bonds) - covered if bonds else set()


def _needy_rings(defect_rings, hex_bonds) -> List[Tuple[int, ...]]:
    """Defect rings carrying a bond ``_uncovered_defect_bonds`` leaves bare."""
    bare = _uncovered_defect_bonds(defect_rings, hex_bonds)
    return [r for r in defect_rings if _ring_bonds(r) & bare]


def _sheet_rings(G: nx.Graph) -> List[Tuple[int, ...]]:
    """SSSR rings of a carbon graph, each as atoms in ring order.

    Measured on the finished graph rather than trusted from growth bookkeeping,
    because a contact bond closes rings growth never planned. No aromaticity
    perception is run, so this is safe on any sheet.
    """
    rw = RWMol()
    order = sorted(G.nodes())
    idx = {n: i for i, n in enumerate(order)}
    for _ in order:
        rw.AddAtom(Chem.Atom(6))
    for a, b in G.edges():
        rw.AddBond(idx[a], idx[b], Chem.BondType.SINGLE)
    m = rw.GetMol()
    m.UpdatePropertyCache(strict=False)
    return [tuple(r) for r in Chem.GetSymmSSSR(m)]


def _aromaticity_search_is_bounded(rings: List[Tuple[int, ...]]) -> bool:
    """Whether RDKit can perceive aromaticity on a sheet with these rings.

    RDKit stops searching ring subsets once every ring bond is aromatic, and
    otherwise enumerates every subset of up to six rings. A sheet passes when
    it has no bare bond (``_uncovered_defect_bonds``), or when it has few
    enough rings that the full enumeration is cheap anyway -- which is what
    keeps small pentagon-only sheets (no heptagon to pair with) buildable.
    """
    if len(rings) <= _RDKIT_EXHAUSTIVE_RING_LIMIT:
        return True
    hex_bonds = set().union(*[_ring_bonds(r) for r in rings if len(r) == 6])
    # Rings of other sizes are holes where a front could not close (a flat
    # sheet cannot close around a 5-7-5 cluster's net disclination). Their
    # bonds also border real rings, so they cost RDKit nothing as long as
    # those are covered.
    others = [r for r in rings if len(r) != 6]
    return not _needy_rings(others, hex_bonds)


def _add_missing_lattice_bonds(
    G: nx.Graph,
    positions: Dict[int, Tuple[float, float]],
    tolerance: float = 1.45,
) -> None:
    """
    Add any bond-distance C-C bonds that are missing from the graph.

    This can occur when position-merging places two existing nodes at adjacent
    hex-lattice sites (1.42 Å) as non-adjacent members of the same fused ring
    (e.g. positions new[0] and new[3]).  The ring-sequence edges connect them
    through intermediates but not directly; the direct edge must be added here.

    The tolerance (1.45 Å < 2.46 Å = min non-bonded distance in graphene)
    guarantees that only genuine lattice bonds are added, never spurious ones.

    Operates **in place** on G.
    """
    pos_list = [(n, positions[n]) for n in G.nodes() if n in positions]
    for i, (na, (ax, ay)) in enumerate(pos_list):
        for nb, (bx, by) in pos_list[i + 1:]:
            if G.has_edge(na, nb):
                continue
            if math.hypot(ax - bx, ay - by) < tolerance:
                G.add_edge(na, nb)


def _get_fusable_edges(G: nx.Graph) -> List[Tuple[int, int]]:
    """Return edges where both endpoints have degree <= 2 (open for ring fusion)."""
    return [(u, v) for u, v in G.edges() if G.degree(u) == 2 and G.degree(v) == 2]


def _edge_neighbor_degree(G: nx.Graph, a: int, b: int) -> int:
    """Sum of degrees of neighbours of a and b (excluding each other)."""
    na = {n for n in G.neighbors(a) if n != b}
    nb = {n for n in G.neighbors(b) if n != a}
    return sum(G.degree(n) for n in na | nb)


# ---------------------------------------------------------------------------
# PAH Assembler (main public API)
# ---------------------------------------------------------------------------

class PAHAssembler:
    """
    Build aromatic carbon skeletons from validated PAH building blocks.

    Strategy:
      within PAH_LIBRARY's range (6-40C): use a pre-validated PAH, exactly when
        the target matches a library size and the nearest otherwise.
      beyond it: take the closest library structure as a seed, then grow by ring
        fusion until the count is reached or passed.

    Growth proceeds in whole rings, so only certain carbon counts are reachable
    and the realised count is >= the request, never below (rqm/carbon-skeleton.md).
    """

    # Sorted list of (num_carbons, pah_name) for quick lookup
    _SIZE_INDEX: Optional[List[Tuple[int, str]]] = None

    def __init__(self, seed: int = None):
        self.seed = seed
        # Instance-level RNG — does not touch the global random state.
        # _grow_graph creates its own random.Random(seed) internally, so this
        # is only needed for any future per-instance randomness.
        self._rng = random.Random(seed)

        # Preload PAH library
        self.pahs = {}
        # A dropped entry shrinks the working library without changing anything a
        # caller can see: a target that was met exactly from a pre-validated
        # structure starts being grown instead. PAH_LIBRARY states that every
        # entry is validated, so a drop means that claim is no longer true and is
        # worth saying out loud rather than absorbing.
        for name, data in PAH_LIBRARY.items():
            mol = Chem.MolFromSmiles(data["smiles"])
            if mol is None:
                logger.warning(
                    "PAH_LIBRARY entry %r has a SMILES RDKit cannot parse; it is "
                    "unavailable, so targets near %s carbons will be grown rather "
                    "than taken from the library.",
                    name, data.get("num_atoms", "?"),
                )
                continue
            try:
                Chem.SanitizeMol(mol)
            except Exception as exc:
                logger.warning(
                    "PAH_LIBRARY entry %r failed sanitisation (%s: %s); it is "
                    "unavailable, so targets near %s carbons will be grown rather "
                    "than taken from the library.",
                    name, type(exc).__name__, exc, data.get("num_atoms", "?"),
                )
                continue
            nC = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
            self.pahs[name] = {"mol": mol, "num_carbons": nC, "data": data}

        # Build size index (sorted ascending by carbon count)
        if PAHAssembler._SIZE_INDEX is None:
            PAHAssembler._SIZE_INDEX = sorted(
                [(info["num_carbons"], name) for name, info in self.pahs.items()]
            )

    @staticmethod
    def _compute_seed_positions(
        seed_mol: Chem.Mol,
    ) -> Optional[Dict[int, Tuple[float, float]]]:
        """
        Generate 2D positions for a library PAH seed and rescale so the
        mean C–C bond length equals ``_AROMATIC_CC_BOND``.

        Works reliably on small library PAHs (benzene, naphthalene,
        pyrene, coronene, …) where RDKit's Compute2DCoords produces a
        clean hexagonal layout.  Returns ``None`` if the layout can't be
        computed.
        """
        try:
            mol = Chem.Mol(seed_mol)
            AllChem.Compute2DCoords(mol)
            conf = mol.GetConformer()
            positions: Dict[int, Tuple[float, float]] = {}
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() != 6:
                    continue
                p = conf.GetAtomPosition(atom.GetIdx())
                positions[atom.GetIdx()] = (float(p.x), float(p.y))
            if not positions:
                return None

            # Rescale so mean bond length ≈ _AROMATIC_CC_BOND (RDKit's 2D
            # coords come out with ~1.5 Å bonds; we want 1.42 Å).
            bond_lengths = []
            for bond in mol.GetBonds():
                a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                if a in positions and b in positions:
                    dx = positions[a][0] - positions[b][0]
                    dy = positions[a][1] - positions[b][1]
                    bond_lengths.append(math.hypot(dx, dy))
            if bond_lengths:
                mean_bl = sum(bond_lengths) / len(bond_lengths)
                if mean_bl > 1e-6:
                    scale = _AROMATIC_CC_BOND / mean_bl
                    positions = {
                        k: (v[0] * scale, v[1] * scale) for k, v in positions.items()
                    }
            return positions
        except Exception:
            return None

    def generate(
        self,
        target_num_carbons: int,
        prefer_larger_pahs: bool = True,
        defect_fraction: float = 0.0,
        target_h_c: Optional[float] = None,
        heptagon_fraction: float = 0.0,
    ) -> CarbonSkeleton:
        """
        Generate a carbon skeleton of approximately *target_num_carbons*.

        Args:
            target_num_carbons: Desired carbon count. Reachable only in the
                increments ring fusion allows, so the realised count is >= this
                and never below it. Measure the returned molecule for the exact
                figure.

            There is deliberately no aromaticity argument. Aromaticity is decided
            by the ring topology this builds, and a parameter that cannot be
            honoured reads at the call site as a request being made. The
            composition-level target lives on GeneratorConfig, which clamps it to
            what this builder can realise (rqm/generation-config.md).
            defect_fraction: Probability [0, 1) that any ring addition
                during graph growth is a 5-membered (pentagon) ring.
                0.0 = pure hexagonal PAH (default).
                0.1 = roughly 10% pentagons.
            heptagon_fraction: Probability [0, 1) that any ring addition is a
                7-membered (heptagon) ring. Together with pentagons this adds the
                curvature of non-graphitizing carbons (Wood et al. 2024).
            target_h_c: Optional target H/C ratio.  When provided and above the
                compact ceiling for this size, growth is steered toward a
                less-condensed (higher-perimeter) flake so more edge carbons are
                available to carry hydrogen.  ``None`` keeps fully compact growth.
        """
        mol = None

        if defect_fraction <= 0.0 and heptagon_fraction <= 0.0:
            # --- Route 1: exact library match (only for pure hexagon mode) ---
            for name, info in self.pahs.items():
                if info["num_carbons"] == target_num_carbons:
                    mol = Chem.Mol(info["mol"])
                    # Attach hex-lattice positions from Compute2DCoords so
                    # downstream geometry generation preserves planarity and
                    # uniform bond lengths.
                    seed_positions = self._compute_seed_positions(mol)
                    if seed_positions:
                        for atom in mol.GetAtoms():
                            idx = atom.GetIdx()
                            if idx in seed_positions:
                                x, y = seed_positions[idx]
                                atom.SetDoubleProp("init_x", x)
                                atom.SetDoubleProp("init_y", y)
                    break

        # --- Route 2: library seed + graph growth (with optional pentagons) ---
        if mol is None:
            mol = self._build_from_seed(
                target_num_carbons, defect_fraction, target_h_c, heptagon_fraction
            )

        # No fallback structure. Growth failing is a failure of the request, and
        # answering a 200-carbon request with the library's 16-carbon pyrene
        # produces a molecule every later stage decorates, embeds, types and
        # exports as though it were what was asked for -- nothing downstream
        # re-checks the count. A caller can recover from an exception; it cannot
        # recover from a plausible molecule of the wrong size.
        if mol is None:
            curved = defect_fraction > 0.0 or heptagon_fraction > 0.0
            attempts = _MAX_CURVED_BUILD_ATTEMPTS if curved else _MAX_BUILD_ATTEMPTS
            hint = ""
            if defect_fraction > 0.0 and heptagon_fraction <= 0.0:
                hint = (
                    " Pentagons without heptagons cannot be paired, and a sheet "
                    f"of more than {_RDKIT_EXHAUSTIVE_RING_LIMIT} rings (about "
                    "100 carbons) with an unpaired pentagon on its edge makes "
                    "RDKit's aromaticity perception run for minutes, so such "
                    "sheets are refused; add heptagons (Wood et al. 2024: "
                    "heptagon_fraction about half of defect_fraction)."
                )
            raise SkeletonError(
                f"could not build a carbon skeleton for "
                f"target_num_carbons={target_num_carbons} "
                f"(defect_fraction={defect_fraction}, "
                f"heptagon_fraction={heptagon_fraction}) after "
                f"{attempts} growth attempts derived from "
                f"seed={self.seed}. Growth depends on the seed, so retrying with "
                f"another seed may succeed; otherwise reduce the target or relax "
                f"the ring fractions.{hint}"
            )

        Chem.SetAromaticity(mol, Chem.AromaticityModel.AROMATICITY_MDL)
        return self._make_skeleton(mol)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _build_from_seed(
        self, target: int, defect_fraction: float = 0.0,
        target_h_c: Optional[float] = None,
        heptagon_fraction: float = 0.0,
    ) -> Optional[Chem.Mol]:
        """
        Build a carbon skeleton of *target* carbons using a library seed
        plus ring-growth.

        Pure hexagon mode (defect_fraction == 0):
            Each ring addition adds exactly 4 nodes, so seed selection uses
            the parity constraint (target - seed_size) % 4 == 0.

        Defect mode (defect_fraction > 0):
            Pentagon additions add 3 nodes, so the parity constraint is
            relaxed — any even-count seed works because the growth loop
            handles parity itself.  Up to *_MAX_CURVED_BUILD_ATTEMPTS* attempts
            are made with different sub-seeds if growth stalls short of the
            target or kekulization fails.
        """
        _MAX_ELONGATE_RETRIES = 4

        if defect_fraction <= 0.0 and heptagon_fraction <= 0.0:
            # --- Pure hexagon route ---
            # Always build the reliable, fully-compact skeleton first (largest
            # fitting seed; deterministic and historically always kekulizable).
            # This is byte-identical to the previous behaviour AND measures the
            # true compact H/C ceiling for this size.
            compact_mol, adj_target = self._build_pure_hex_compact(target)
            if compact_mol is None:
                return None
            compact_ceiling = _edge_carbon_fraction(compact_mol)

            # No target, or the target is already reachable by the compact flake
            # (within a small margin): keep it -- no elongation, no extra builds.
            if target_h_c is None or target_h_c <= compact_ceiling + 0.02:
                return compact_mol

            # Target exceeds the compact ceiling: grow a less-condensed
            # (higher-perimeter) flake from a small seed.  Keep whichever
            # candidate's edge-carbon fraction is *closest* to the target -- so
            # an over-elongated flake that shoots past the target is rejected in
            # favour of the compact flake (or a milder elongation) rather than
            # returning an H/C that is now too high.
            comp = _compactness_for_hc(target_h_c, compact_ceiling)
            base_seed = self.seed if self.seed is not None else 0
            best_mol = compact_mol
            best_err = abs(compact_ceiling - target_h_c)
            for attempt in range(_MAX_ELONGATE_RETRIES):
                cand = self._build_pure_hex_elongated(
                    adj_target, comp, base_seed + attempt
                )
                if cand is None:
                    continue
                frac = _edge_carbon_fraction(cand)
                err = abs(frac - target_h_c)
                if err < best_err:
                    best_mol, best_err = cand, err
                if frac >= target_h_c:
                    # Reached/overshot the target; further tries only move away.
                    break
            return best_mol

        else:
            # --- Defect mode: relax parity, retry on kekulization failure ---
            # Use the largest even-count library seed that fits
            seed_mol = None
            seed_carbons = 0
            for nC, name in reversed(self._SIZE_INDEX):
                if nC <= target and name in self.pahs:
                    seed_mol = self.pahs[name]["mol"]
                    seed_carbons = nC
                    break

            if seed_mol is None:
                seed_mol = Chem.MolFromSmiles(PAH_LIBRARY["benzene"]["smiles"])
                seed_carbons = 6

            seed_positions = self._compute_seed_positions(seed_mol)

            if seed_carbons >= target:
                out = Chem.Mol(seed_mol)
                if seed_positions:
                    for atom in out.GetAtoms():
                        idx = atom.GetIdx()
                        if idx in seed_positions:
                            x, y = seed_positions[idx]
                            atom.SetDoubleProp("init_x", x)
                            atom.SetDoubleProp("init_y", y)
                return out

            base_seed = self.seed if self.seed is not None else 0
            for attempt in range(_MAX_CURVED_BUILD_ATTEMPTS):
                G = _mol_to_graph(seed_mol)
                grown_positions = dict(seed_positions) if seed_positions else None
                # Defect mode stays compact: pentagon insertion already perturbs
                # topology and H count, and stacking elongation on top raises the
                # kekulization-failure rate.  H/C-aware elongation is a
                # pure-hexagon feature (see the defect_fraction == 0 branch).
                G = _grow_graph(
                    G, target,
                    seed=_sub_seed(base_seed, attempt),
                    defect_fraction=defect_fraction,
                    positions=grown_positions, compactness=1.0,
                    heptagon_fraction=heptagon_fraction,
                )
                # Every front blocked before the target: a sheet below the
                # request breaks rq-59e23c2e, so it is a failed attempt.
                if G.number_of_nodes() < target:
                    continue
                # Checked before RDKit sees the sheet: an isolated defect makes
                # its aromaticity perception exponential (see _grow_curved), so
                # converting first would hang rather than fail.
                rings = _sheet_rings(G)
                if not _aromaticity_search_is_bounded(rings):
                    continue
                # Above the exhaustive limit the check above passes a sheet
                # whose defect draws all happened to be lost, and retrying
                # until one does would quietly answer a defect request with a
                # flat sheet. Such an attempt is a failed one.
                if len(rings) > _RDKIT_EXHAUSTIVE_RING_LIMIT:
                    sizes = {len(r) for r in rings}
                    if ((defect_fraction > 0 and 5 not in sizes)
                            or (heptagon_fraction > 0 and 7 not in sizes)):
                        continue
                mol = _graph_to_mol(G, positions=grown_positions, kekulize_by_matching=True)
                if mol is not None:
                    return mol
                # kekulization failed — try again with a different sub-seed

            return None  # all retries exhausted

    def _build_pure_hex_compact(
        self, target: int
    ) -> Tuple[Optional[Chem.Mol], int]:
        """
        Build a maximally condensed pure-hexagon skeleton from the *largest*
        fitting library seed (fast, deterministic, always kekulizable).

        Returns ``(mol, adjusted_target)``; *adjusted_target* may be bumped up to
        satisfy the ``(target - 6) % 4 == 0`` parity constraint when falling back
        to the benzene seed.
        """
        seed_mol = None
        seed_carbons = 0
        for nC, name in reversed(self._SIZE_INDEX):
            if nC <= target and (target - nC) % 4 == 0 and name in self.pahs:
                seed_mol = self.pahs[name]["mol"]
                seed_carbons = nC
                break

        if seed_mol is None:
            rem = (target - 6) % 4
            if rem != 0:
                target = target + (4 - rem)
            seed_mol = Chem.MolFromSmiles(PAH_LIBRARY["benzene"]["smiles"])
            seed_carbons = 6

        seed_positions = self._compute_seed_positions(seed_mol)

        if seed_carbons >= target:
            out = Chem.Mol(seed_mol)
            _attach_seed_positions(out, seed_positions)
            return out, target

        # Seed-dependent: the RNG breaks ties between equally ranked edges, so
        # one seed can close the sheet with an odd count or unbalanced
        # sublattices (no Kekule structure) where the next does not. Attempt 0
        # is self.seed, so a build that succeeds first time is unchanged.
        for attempt in range(_MAX_BUILD_ATTEMPTS):
            G = _mol_to_graph(seed_mol)
            grown_positions = dict(seed_positions) if seed_positions else None
            G = _grow_graph(
                G, target, seed=_sub_seed(self.seed, attempt), defect_fraction=0.0,
                positions=grown_positions, compactness=1.0,
            )
            mol = _graph_to_mol(G, positions=grown_positions)
            if mol is not None:
                return mol, target
        return None, target

    def _build_pure_hex_elongated(
        self, target: int, compactness: float, seed: int
    ) -> Optional[Chem.Mol]:
        """
        Build a less-condensed pure-hexagon skeleton from the *smallest* fitting
        seed, so growth (and its ``compactness`` control) shapes most of the
        flake.  Returns the mol, or ``None`` on kekulization failure.
        """
        seed_mol = None
        seed_carbons = 0
        for nC, name in self._SIZE_INDEX:
            if nC <= target and (target - nC) % 4 == 0 and name in self.pahs:
                seed_mol = self.pahs[name]["mol"]
                seed_carbons = nC
                break
        if seed_mol is None:
            seed_mol = Chem.MolFromSmiles(PAH_LIBRARY["benzene"]["smiles"])
            seed_carbons = 6

        if seed_carbons >= target:
            return None  # seed already fills the target; nothing to elongate

        seed_positions = self._compute_seed_positions(seed_mol)
        G = _mol_to_graph(seed_mol)
        grown_positions = dict(seed_positions) if seed_positions else None
        G = _grow_graph(
            G, target, seed=seed, defect_fraction=0.0,
            positions=grown_positions, compactness=compactness,
        )
        return _graph_to_mol(G, positions=grown_positions)

    @staticmethod
    def _make_skeleton(mol: Chem.Mol) -> CarbonSkeleton:
        num_carbons = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
        num_aromatic = sum(
            1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and a.GetIsAromatic()
        )
        aromaticity = (num_aromatic / num_carbons * 100) if num_carbons > 0 else 0
        try:
            smiles = Chem.MolToSmiles(mol)
        except Exception:
            smiles = f"C{num_carbons}_PAH"
        rings = mol.GetRingInfo().AtomRings()
        ring_comp = {
            "hexagons": sum(1 for r in rings if len(r) == 6),
            "pentagons": sum(1 for r in rings if len(r) == 5),
            "heptagons": sum(1 for r in rings if len(r) == 7),
        }
        return CarbonSkeleton(
            mol=mol,
            smiles=smiles,
            num_carbons=num_carbons,
            num_aromatic_carbons=num_aromatic,
            aromaticity_percent=aromaticity,
            ring_composition=ring_comp,
        )


# ---------------------------------------------------------------------------
# RandomGraphGenerator (kept for backward compatibility)
# ---------------------------------------------------------------------------

class RandomGraphGenerator:
    """
    Generate random aromatic carbon frameworks.

    Delegates to hex-lattice builder + ring growth for all sizes.
    """

    def __init__(self, seed: int = None):
        self.seed = seed
        self._rng = random.Random(seed)

    def generate(
        self, target_num_carbons: int, target_aromaticity: float = 100.0
    ) -> CarbonSkeleton:
        # Delegate to PAHAssembler for consistent parity-aware seed selection
        asm = PAHAssembler(seed=self.seed)
        return asm.generate(target_num_carbons, target_aromaticity)

    def _create_pah_template(self, num_carbons: int) -> str:
        if num_carbons <= 6:
            return PAH_LIBRARY["benzene"]["smiles"]
        elif num_carbons <= 10:
            return PAH_LIBRARY["naphthalene"]["smiles"]
        elif num_carbons <= 14:
            return PAH_LIBRARY["anthracene"]["smiles"]
        else:
            return PAH_LIBRARY["pyrene"]["smiles"]


# ---------------------------------------------------------------------------
# Skeleton Validator
# ---------------------------------------------------------------------------

class SkeletonValidator:
    """Validate carbon skeletons for chemical feasibility."""

    @staticmethod
    def validate(skeleton: CarbonSkeleton) -> Tuple[bool, List[str]]:
        errors = []

        if skeleton.mol is None:
            errors.append("Molecule object is None")
            return False, errors

        if skeleton.num_carbons == 0:
            errors.append("No carbon atoms in skeleton")
            return False, errors

        if skeleton.num_aromatic_carbons > skeleton.num_carbons:
            errors.append("More aromatic carbons than total carbons")
            return False, errors

        if skeleton.aromaticity_percent < 0 or skeleton.aromaticity_percent > 100:
            errors.append(f"Invalid aromaticity: {skeleton.aromaticity_percent}%")
            return False, errors

        mol_graph = nx.Graph()
        mol_graph.add_nodes_from(range(skeleton.mol.GetNumAtoms()))
        for bond in skeleton.mol.GetBonds():
            mol_graph.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

        if not nx.is_connected(mol_graph):
            errors.append("Molecular graph is not connected")
            return False, errors

        return len(errors) == 0, errors
