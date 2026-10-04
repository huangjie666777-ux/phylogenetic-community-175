"""Community mass construction and p=1 Kantorovich-Rubinstein tree distance.

A placement result distributes one read over every candidate edge according
to the edge's like_weight_ratio (not just the best edge, and not moved to the
tip).  Each share is mass at the grafting point encoded by the candidate's
distal_length on that edge; pendant_length does not enter the geometry.

For two normalised mass distributions mu, nu on the same tree, the p=1
Kantorovich-Rubinstein (Wasserstein-1 / earth-mover) distance is

    KR(mu, nu) = sum over edges  integral_edge |cumulative distal mass| dx,

where the integrand at a point is the signed net mass in the subtree on the
distal side of that point.  Contributions are non-negative and sum exactly
to the distance.
"""

from __future__ import annotations

from dataclasses import dataclass

from .jplacetree import parse_jplace_tree
from .tree import TreeGraph


@dataclass
class Atom:
    """Mass located on edge eid at coordinate x from the parent node u."""

    eid: int
    x: float
    mass: float


@dataclass
class SampleMass:
    sample_id: str
    atoms: list[Atom]
    effective_total: float          # sum of placed reads (normalisation basis)
    excluded: list[dict]            # [{id, count, reason}]


def placement_index(job: dict) -> dict[str, dict[int, dict]]:
    """Map query id -> edge_num -> placement row from a stored job's jplace."""
    jplace = job["jplace"]
    fields = list(jplace["fields"])
    index: dict[str, dict[int, dict]] = {}
    for entry in jplace["placements"]:
        qid = entry["n"][0]
        rows = {}
        for row in entry["p"]:
            rec = dict(zip(fields, row))
            rows[int(rec["edge_num"])] = rec
        index[qid] = rows
    return index


def build_sample_mass(sample_id: str, counts: dict[str, int],
                      job: dict) -> SampleMass:
    """Distribute reads over candidate insertion points.

    Read ids listed in the job's 'unplaceable' map (or absent from the
    placement table) are excluded with their reason; zero-count ids are
    silently dropped.  Like_weight_ratio is used across *all* candidate
    edges, including candidates with small support.
    """
    index = placement_index(job)
    atoms: list[Atom] = []
    excluded: list[dict] = []
    effective_total = 0
    for qid, count in counts.items():
        if count == 0:
            continue
        reason = job["unplaceable"].get(qid)
        rows = index.get(qid)
        if reason is not None or not rows:
            excluded.append({
                "id": qid,
                "count": count,
                "reason": reason or ("no placement record for this id in the "
                                     "referenced job"),
            })
            continue
        effective_total += count
        edge_count = job["edge_count"]
        partial: list[Atom] = []
        for edge_num in range(edge_count):
            rec = rows[edge_num]
            # placement optimisation returned a valid interior point
            distal = float(rec["distal_length"])
            pendant = float(rec["pendant_length"])
            if not (distal >= 0.0 and pendant >= 0.0):
                excluded.append({
                    "id": qid,
                    "count": count,
                    "reason": ("candidate on edge " + str(edge_num) +
                               " has a negative branch length"),
                })
                break
            weight = float(rec["like_weight_ratio"])
            if weight <= 0.0:
                continue
            partial.append(Atom(eid=edge_num, x=distal,
                                mass=count * weight))
        else:
            atoms.extend(partial)
            continue
        effective_total -= count
    return SampleMass(sample_id=sample_id, atoms=atoms,
                      effective_total=float(effective_total),
                      excluded=excluded)


def geometry_of(job: dict) -> TreeGraph:
    """Recover oriented tree geometry from the job's annotated jplace tree."""
    return parse_jplace_tree(job["jplace"]["tree"])


def pairwise_distances(samples: list[SampleMass],
                       g: TreeGraph) -> tuple[list[list[float]], list[dict]]:
    """Compute the symmetric KR distance matrix and per-pair contributions.

    Returns (matrix, pair_records); pair_records is in input order over
    unordered pairs i<j, each with every edge's non-negative contribution.
    """
    n = len(samples)
    # Normalise each sample by its effective (placed) read total.
    normed: list[dict[int, list[tuple[float, float]]]] = []
    for smp in samples:
        total = smp.effective_total
        edge_atoms: dict[int, list[tuple[float, float]]] = {}
        for atom in smp.atoms:
            # x stored as distal (from v); convert to u-coordinate
            edge = g.edges[atom.eid]
            x = edge.length - atom.x
            edge_atoms.setdefault(atom.eid, []).append((x, atom.mass / total))
        normed.append(edge_atoms)

    # Net distal-subtree mass D(v) for every pair: D_i(v) - D_j(v).
    # Root the orientation at the jplace root (the node with no parent edge).
    parent_edge: dict[int, int] = {}
    children: dict[int, list[int]] = {}
    for eid, e in enumerate(g.edges):
        parent_edge[e.v] = eid
        children.setdefault(e.u, []).append(e.v)
    root = next(n for n in g.adj if n not in parent_edge)

    order: list[int] = []
    stack = [root]
    while stack:
        u = stack.pop()
        order.append(u)
        stack.extend(children.get(u, ()))

    def net_subtree(i: int, j: int) -> dict[int, float]:
        edge_atom_net: dict[int, list[tuple[float, float]]] = {}
        for eid in range(len(g.edges)):
            mi = {x: m for x, m in normed[i].get(eid, ())}
            mj = {x: m for x, m in normed[j].get(eid, ())}
            xs = sorted(set(mi) | set(mj))
            edge_atom_net[eid] = [(x, mi.get(x, 0.0) - mj.get(x, 0.0))
                                  for x in xs]
        # subtree(v) includes all atoms on edges in v's descendant tree.
        net = {n: 0.0 for n in g.adj}
        for u in reversed(order):
            acc = 0.0
            for v in children.get(u, ()):
                eid = parent_edge[v]
                acc += net[v]
                for _, dm in edge_atom_net[eid]:
                    acc += dm
            net[u] = acc
        return net, edge_atom_net

    matrix = [[0.0] * n for _ in range(n)]
    pair_records: list[dict] = []
    for i in range(n):
        for j in range(i + 1, n):
            net, edge_atom_net = net_subtree(i, j)
            contributions = []
            distance = 0.0
            for eid, e in enumerate(g.edges):
                c = _edge_integral(e.length, net[e.v], edge_atom_net[eid])
                if -1e-12 < c < 0.0:
                    c = 0.0
                contributions.append({"edge_num": eid, "contribution": c})
                distance += c
            matrix[i][j] = matrix[j][i] = distance
            pair_records.append({
                "sample_i": samples[i].sample_id,
                "sample_j": samples[j].sample_id,
                "distance": distance,
                "edge_contributions": contributions,
            })
    return matrix, pair_records


def _edge_integral(length: float, distal_base: float,
                   atoms: list[tuple[float, float]]) -> float:
    """Integrate |net distal mass| along one oriented edge.

    distal_base is D(v), the net mass strictly below the child node.  An atom
    at u-coordinate x adds its signed mass to the distal side for points
    strictly beyond x; at an exact atom location the left/right constant
    pieces are each integrated over their own intervals, so endpoint atoms
    (x=0 or x=length) fall out naturally.
    """
    # D(v) counts all atoms on this edge too; subtract them to get the base.
    base = distal_base - sum(dm for _, dm in atoms)
    prev_x = 0.0
    level = base
    total = 0.0
    for x, dm in sorted(atoms):
        total += abs(level) * (x - prev_x)
        level += dm
        prev_x = x
    total += abs(level) * (length - prev_x)
    return total
