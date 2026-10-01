"""Search exact five-mask/five-CNOT parity tours and Fourier-tail absorption.

Each CNOT must expose a new required nonlocal mask, so a constrained DFS
enumerates all possible five-step tours for each of the 54 phase models from
run 65. Endpoints within inverse CNOT distance <=4 are tested with Qiskit
level-3 tail synthesis, seed 42. This is a restricted topology search.
"""

import json
from pathlib import Path

from qiskit import transpile

from algebraic14_alt_pair_sweep import prefix_gates
from algebraic14_five_prefix import BASIS, inverse_paths
from check_circuit import evaluate
from qiskit_crosscheck import from_qiskit, to_qiskit


ROOT = Path(__file__).parent


def tours(required):
    required = set(required)
    ends = {}
    paths = 0

    def visit(rows, seen, edges):
        nonlocal paths
        if len(edges) == 5:
            if seen == required:
                paths += 1
                ends.setdefault(rows, edges)
            return
        for control in range(4):
            for target in range(4):
                if control == target:
                    continue
                mask = rows[target] ^ rows[control]
                if mask not in required or mask in seen:
                    continue
                new_rows = list(rows)
                new_rows[target] = mask
                visit(tuple(new_rows), seen | {mask},
                      edges + ((control, target),))

    visit(BASIS, set(), ())
    return ends, paths


def main():
    models = json.loads((ROOT / "algebraic14_fivemask_models.json").read_text())
    baseline = json.loads((ROOT / "baseline_18.json").read_text())
    tail = baseline["gates"][18:]
    inverses = inverse_paths(4)
    summary = []
    best = None
    best_candidate = None
    total_evaluated = 0
    for index, model in enumerate(models):
        support = tuple(model["support"])
        coefficients = model["coefficients_pi_over_8"]
        endpoints, paths = tours(support)
        counts = {}
        local_best = None
        for rows, edges in endpoints.items():
            inverse = inverses.get(rows)
            distance = len(inverse) if inverse is not None else ">4"
            counts[str(distance)] = counts.get(str(distance), 0) + 1
            if inverse is None:
                continue
            tail_input = {"n": 4, "gates": [
                *({"gate": "cx", "control": c, "target": t}
                  for c, t in inverse), *tail]}
            optimized = transpile(to_qiskit(tail_input),
                                  basis_gates=["u", "cx"],
                                  optimization_level=3, seed_transpiler=42)
            candidate = {"n": 4, "gates": prefix_gates(edges, coefficients)
                         + from_qiskit(optimized)["gates"]}
            check = evaluate(candidate)
            total_evaluated += 1
            record = {"model_index": index, "support": support,
                      "endpoint_rows": rows, "inverse_distance": len(inverse),
                      "tail_cnot_count": optimized.count_ops().get("cx", 0),
                      "total_cnot_count": check["cnot_count"],
                      "off_diagonal_error": check["off_diagonal_error"],
                      "valid": check["valid_diagonalizer"]}
            if local_best is None or (record["total_cnot_count"],
                                      record["off_diagonal_error"]) < (
                                          local_best["total_cnot_count"],
                                          local_best["off_diagonal_error"]):
                local_best = record
            if best is None or (record["total_cnot_count"],
                                record["off_diagonal_error"]) < (
                                    best["total_cnot_count"],
                                    best["off_diagonal_error"]):
                best, best_candidate = record, candidate
        summary.append({"model_index": index, "support": support,
                        "five_step_paths": paths,
                        "distinct_endpoints": len(endpoints),
                        "inverse_distance_histogram": counts,
                        "best_total_cnot_count": local_best["total_cnot_count"]
                        if local_best else None})
    if best_candidate is not None:
        path = ROOT / "algebraic14_fivemask_best.json"
        path.write_text(json.dumps(best_candidate, indent=2) + "\n")
        best["candidate"] = path.name
    print(json.dumps({"models_tested": len(models),
                      "models_with_paths": sum(item["five_step_paths"] > 0
                                               for item in summary),
                      "total_tail_syntheses": total_evaluated,
                      "summary": summary, "best": best}, indent=2))


if __name__ == "__main__":
    main()
