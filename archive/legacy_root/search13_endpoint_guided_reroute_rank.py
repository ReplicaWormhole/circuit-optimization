"""Exact residual cut ranks for endpoint-guided prefixes and tail reroutes.

For a fixed phase-informed five-CNOT prefix P5 and exact14 diagonalizer
U14, the residual target U14 P5^-1 is computed over Q(zeta_96).  Cut-rank
compatibility is necessary for implementing this fixed target with a given
eight-CNOT tail, but not for other valid output eigenspace gauges or a
jointly varied prefix.
"""

import json
from fractions import Fraction

from search13_firstpair_mask_relocation import ROOT
from search13_secondpair_mask_transport import CyclotomicPair
from search13_crosspair_single_reroute_rank import CUTS, exact_rank, realignment
from search13_joint_gauge_phased_exact import base_matrices


BASE = ROOT / "topology14_exact_matchgate_rational.json"
SCREEN = ROOT / "search13_endpoint_guided_screen_result.json"
SELECTED = (("2->0", 0), ("2->0", 3), ("1->3", 0), ("1->3", 1))
TAIL = ((0, 2), (0, 2), (1, 3), (1, 3),
        (0, 1), (0, 1), (2, 3), (2, 3))
ROUTES = ((0, 3), (3, 0), (1, 2), (2, 1))
PHASES = {7: Fraction(-1, 4), 6: Fraction(-1, 2),
          14: Fraction(1, 8), 11: Fraction(-1, 8)}


def phase_prefix(path):
    gates = [("rz", 1, Fraction(5, 8)), ("rz", 2, Fraction(3, 4)),
             ("rz", 3, Fraction(3, 8))]
    rows = [8, 4, 2, 1]
    seen = set()
    for a, b in path:
        gates.append(("cx", a, b))
        rows[b] ^= rows[a]
        mask = rows[b]
        if mask in PHASES and mask not in seen:
            gates.append(("rz", b, PHASES[mask]))
            seen.add(mask)
    return gates


def apply(gates, exact):
    matrix = exact.identity(16)
    for gate in gates:
        if gate[0] == "cx":
            _, control, target = gate
            cmask, tmask = 8 >> control, 8 >> target
            out = [None]*16
            for row in range(16):
                new = row ^ tmask if row & cmask else row
                out[new] = matrix[row]
            matrix = out
        else:
            _, wire, angle = gate
            exponent = 24*angle
            assert exponent.denominator == 1
            positive = exact.z**int(exponent)
            negative = positive**-1
            mask = 8 >> wire
            matrix = [[(positive if row & mask else negative)*value
                       for value in matrix[row]] for row in range(16)]
    return matrix


def inverse(gates):
    return [(name, a, b) if name == "cx" else (name, a, -b)
            for name, a, b in reversed(gates)]


def rank_map(matrix, exact):
    return {"".join(map(str, cut)):
            exact_rank(realignment(matrix, cut, exact), exact) for cut in CUTS}


def reroutes():
    for position in range(8):
        for replacement in ROUTES:
            edges = list(TAIL)
            edges[position] = replacement
            cap = {"".join(map(str, cut)): 2**sum(
                (a in cut) != (b in cut) for a, b in edges) for cut in CUTS}
            yield {"position": position, "edge": replacement,
                   "edges": edges, "cut_capacities": cap}


def main():
    exact = CyclotomicPair()
    full_tail, _ = base_matrices(exact)
    source = json.loads(BASE.read_text())["gates"][:13]
    prefix6 = [(g["gate"], g["control"], g["target"])
               if g["gate"] == "cx" else
               (g["gate"], g["qubit"],
                Fraction(g["theta"].replace("*pi", "")) if "*pi" in g["theta"]
                else Fraction(g["theta"].replace("pi", "1")))
               for g in source]
    p6 = apply(prefix6, exact)
    routes = list(reroutes())
    screen = json.loads(SCREEN.read_text())
    rows = []
    for orientation, index in SELECTED:
        item = screen["selected"][orientation][index]
        path = tuple(tuple(edge) for edge in item["path"])
        p5gates = phase_prefix(path)
        p5inv = apply(inverse(p5gates), exact)
        residual = exact.matmul(exact.matmul(full_tail, p6), p5inv)
        ranks = rank_map(residual, exact)
        compatible = [i for i, route in enumerate(routes)
                      if all(ranks[key] <= cap for key, cap in
                             route["cut_capacities"].items())]
        row = {"orientation": orientation, "selected_index": index,
               "prefix": path, "endpoint": item["endpoint"],
               "visited_masks": item["covered_required_masks"],
               "residual_ranks": ranks,
               "compatible_route_indices": compatible}
        rows.append(row)
        print(json.dumps(row, sort_keys=True), flush=True)
    result = {"field": "Q(zeta_96)", "target": "fixed U14 P5^-1",
              "tail_edges": TAIL, "routes": routes, "rows": rows}
    (ROOT / "search13_endpoint_guided_reroute_rank_result.json").write_text(
        json.dumps(result, indent=2)+"\n")


if __name__ == "__main__":
    main()
