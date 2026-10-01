"""Bounded three-output-CNOT Givens gauge synthesis.

The output CNOT word is a computational-basis permutation.  Exact output
labels select only Givens rotations within equal-eigenvalue subspaces.  A
Schmidt-rank integer program gives a necessary CNOT count for each modified
tail.  Qiskit is then run on a limited selection; its count is an upper bound.
"""

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np
from qiskit import transpile
from qiskit.quantum_info import Operator
from scipy.optimize import Bounds, LinearConstraint, milp

from check_circuit import evaluate
from qiskit_crosscheck import from_qiskit, to_qiskit
from topology16_exact_block import exact_block
from topology16_15_exact import reverse_block


ROOT = Path(__file__).resolve().parent
PAIRS = list(itertools.combinations(range(4), 2))
CUTS = [(0,), (1,), (2,), (3,), (0, 1), (0, 2), (0, 3)]


def schmidt_rank(unitary, cut):
    other = [q for q in range(4) if q not in cut]
    axes = list(cut) + [q + 4 for q in cut] + other + [q + 4 for q in other]
    matrix = unitary.reshape([2] * 8).transpose(axes).reshape(
        4 ** len(cut), 4 ** len(other))
    return int(np.count_nonzero(np.linalg.svd(matrix, compute_uv=False) > 1e-9))


def cnot_rank_lower_bound(unitary):
    rows = []
    requirements = []
    for cut in CUTS:
        rank = schmidt_rank(unitary, cut)
        requirements.append(int(math.ceil(math.log2(rank))))
        rows.append([int((a in cut) != (b in cut)) for a, b in PAIRS])
    result = milp(np.ones(len(PAIRS)), integrality=np.ones(len(PAIRS)),
                  bounds=Bounds(np.zeros(len(PAIRS)), np.full(len(PAIRS), 20)),
                  constraints=LinearConstraint(np.asarray(rows),
                                               np.asarray(requirements), np.inf))
    if not result.success:
        raise RuntimeError(result.message)
    return int(round(result.fun)), requirements


def gauge_circuit(tail, record, angle):
    extras = [{"gate": "cx", "control": c, "target": t}
              for c, t in record["output_cx"]]
    circuit = to_qiskit({"n": 4, "gates": tail + extras})
    q, r = record["wire_pair"]
    circuit.rxx(angle, 3 - q, 3 - r)
    circuit.ryy(angle, 3 - q, 3 - r)
    return circuit


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--angle-over-pi", type=float, default=0.5)
    p.add_argument("--max-trials", type=int, default=20)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--result", type=Path, required=True)
    args = p.parse_args()
    source = json.loads((ROOT / "topology16_15_exact_candidate.json").read_text())
    split = 13 + len(exact_block()) + len(reverse_block())
    prefix, tail = source["gates"][:split], source["gates"][split:]
    records = json.loads((ROOT / "delete14_gauge_pair_screen3_result.json").read_text())["records"]
    records = [record for record in records if len(record["output_cx"]) == 3]
    angle = args.angle_over_pi * math.pi
    unique = {}
    for record in records:
        key = (tuple(record["relabeled_labels"]), tuple(record["wire_pair"]))
        unique.setdefault(key, record)
    ranked = []
    for record in unique.values():
        circuit = gauge_circuit(tail, record, angle)
        lower, cuts = cnot_rank_lower_bound(Operator(circuit).data)
        ranked.append({**record, "tail_rank_lower_bound": lower,
                       "cut_rank_requirements": cuts})
    # Cover each novel wire pair, then select the lowest necessary CNOT counts.
    ranked.sort(key=lambda row: (row["tail_rank_lower_bound"],
                                 row["wire_pair"], row["output_cx"]))
    selected = []
    pair_counts = {pair: 0 for pair in PAIRS}
    for row in ranked:
        pair = tuple(row["wire_pair"])
        if pair_counts[pair] >= 4 or len(selected) >= args.max_trials:
            continue
        selected.append(row)
        pair_counts[pair] += 1
    results = []
    best = None
    best_candidate = None
    for row in selected:
        circuit = gauge_circuit(tail, row, angle)
        optimized = transpile(circuit, basis_gates=["u", "cx"],
                              optimization_level=3, seed_transpiler=args.seed)
        candidate = {"n": 4, "gates": prefix + from_qiskit(optimized)["gates"]}
        check = evaluate(candidate)
        result = {"output_cx": row["output_cx"],
                  "wire_pair": row["wire_pair"],
                  "tail_rank_lower_bound": row["tail_rank_lower_bound"],
                  "tail_cnot_count": optimized.count_ops().get("cx", 0),
                  "total_cnot_count": check["cnot_count"],
                  "off_diagonal_error": check["off_diagonal_error"],
                  "valid": check["valid_diagonalizer"]}
        results.append(result)
        if best is None or (result["total_cnot_count"],
                            result["off_diagonal_error"]) < (
                                best["total_cnot_count"], best["off_diagonal_error"]):
            best, best_candidate = result, candidate
    args.output.write_text(json.dumps(best_candidate, indent=2) + "\n")
    summary = {"angle_over_pi": args.angle_over_pi,
               "unique_motifs": len(unique), "ranked_count": len(ranked),
               "selected_count": len(selected),
               "rank_bound_histogram": {str(k): sum(r["tail_rank_lower_bound"] == k
                                                    for r in ranked)
                                        for k in sorted(set(r["tail_rank_lower_bound"]
                                                            for r in ranked))},
               "best": best, "results": results}
    args.result.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "results"}, indent=2))


if __name__ == "__main__":
    main()
