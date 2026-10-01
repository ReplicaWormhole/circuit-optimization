"""Test commuting Fourier-butterfly reorderings on five-CNOT prefixes.

The fixed-order scans in runs 61 and 68 found prefix5+tail10. Their Qiskit
tail synthesis may depend on the source gate order despite the first two and
last two butterflies acting on disjoint wire pairs. This search tests the
three other legal orders for each endpoint with inverse distance <=4.
"""

import json
from pathlib import Path

import numpy as np
from qiskit import transpile
from qiskit.quantum_info import Operator

from algebraic14_alt_pair_sweep import prefix_gates
from algebraic14_five_prefix import inverse_paths
from algebraic14_fivemask_prefix import tours
from check_circuit import evaluate
from qiskit_crosscheck import from_qiskit, to_qiskit


ROOT = Path(__file__).resolve().parent


def tail_variants(gates):
    first = (gates[18:28], gates[28:38])
    middle = gates[38:41]
    # Gate 41 is Rz on wire 3 and commutes with the wire-(0,1) block.
    last = ([gates[42], *gates[43:52]], [gates[41], *gates[52:62]])
    return [(a, b, first[a] + first[1 - a] + middle
             + last[b] + last[1 - b])
            for a in range(2) for b in range(2)]


def main():
    baseline = json.loads((ROOT / "baseline_18.json").read_text())
    original = to_qiskit({"n": 4, "gates": baseline["gates"][18:]})
    original_u = Operator(original).data
    variants = []
    for a, b, gates in tail_variants(baseline["gates"]):
        candidate_u = Operator(to_qiskit({"n": 4, "gates": gates})).data
        error = float(np.max(np.abs(candidate_u - original_u)))
        assert error < 1e-12, (a, b, error)
        variants.append((a, b, gates))

    models = json.loads((ROOT / "algebraic14_fivemask_models.json").read_text())
    inverses = inverse_paths(4)
    trials = 0
    best = None
    best_candidate = None
    counts = {}
    for model_index, model in enumerate(models):
        endpoints, _ = tours(model["support"])
        for rows, edges in endpoints.items():
            inverse = inverses.get(rows)
            if inverse is None:
                continue
            inverse_gates = [{"gate": "cx", "control": c, "target": t}
                             for c, t in inverse]
            for a, b, tail in variants:
                if a == 0 and b == 0:
                    continue  # Fixed order already tested in run 68.
                optimized = transpile(to_qiskit({"n": 4, "gates":
                                     inverse_gates + tail}),
                                     basis_gates=["u", "cx"],
                                     optimization_level=3,
                                     seed_transpiler=42)
                tail_count = optimized.count_ops().get("cx", 0)
                trials += 1
                counts[str(tail_count)] = counts.get(str(tail_count), 0) + 1
                record = {"model_index": model_index, "endpoint_rows": rows,
                          "inverse_distance": len(inverse),
                          "first_order": a, "last_order": b,
                          "tail_cnot_count": tail_count,
                          "total_cnot_count": 5 + tail_count}
                if best is None or (record["total_cnot_count"],
                                    record["inverse_distance"]) < (
                                        best["total_cnot_count"],
                                        best["inverse_distance"]):
                    best = record
                    best_candidate = {"n": 4, "gates":
                                      prefix_gates(edges,
                                        model["coefficients_pi_over_8"])
                                      + from_qiskit(optimized)["gates"]}
    if best_candidate is not None:
        path = ROOT / "fourier_reorder14_best.json"
        path.write_text(json.dumps(best_candidate, indent=2) + "\n")
        check = evaluate(best_candidate)
        best["candidate"] = path.name
        best["off_diagonal_error"] = check["off_diagonal_error"]
        best["valid"] = check["valid_diagonalizer"]
    result = {"unitary_orders": 4, "trials": trials,
              "tail_cnot_histogram": counts, "best": best}
    (ROOT / "fourier_reorder14_result.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
