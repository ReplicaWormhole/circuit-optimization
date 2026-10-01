"""Six-CNOT open-prefix search across three nontrivial triple-orbit gauges.

Build on the run-32 BFS implementation, replacing the required triple masks
and their angles by shift-invariant gauge choices that cancel each original
triple phase in turn. For all endpoints at CNOT distance <=3 from identity,
resynthesize the inverse endpoint plus fixed Fourier tail with Qiskit.
"""

import json
from pathlib import Path

from qiskit import transpile

import algebraic16_to15_prefix6 as base
from check_circuit import evaluate
from qiskit_crosscheck import from_qiskit, to_qiskit


ROOT = Path(__file__).parent
TRIPLES = (11, 13, 14, 7)
COEFFICIENTS = (1, 2, 3, 0)
PAIR = 6


def main():
    baseline = json.loads((ROOT / "baseline_18.json").read_text())
    tail = baseline["gates"][18:]
    inverses = base.inverse_paths(3)
    results = []
    overall_best = None
    overall_candidate = None
    for cancelled in range(3):
        shift = -COEFFICIENTS[cancelled]
        phases = {mask: f"{coefficient + shift}*pi/8"
                  for mask, coefficient in zip(TRIPLES, COEFFICIENTS)
                  if coefficient + shift}
        phases[PAIR] = "-pi/2"
        base.REQUIRED = tuple(phases)
        base.FLAGS = {mask: 1 << i for i, mask in enumerate(base.REQUIRED)}
        base.ALL_FLAGS = (1 << len(base.REQUIRED)) - 1
        base.ANGLE = phases
        endpoints, states = base.endpoints()
        evaluated = 0
        best = None
        best_candidate = None
        for rows, edges in endpoints:
            inverse = inverses.get(rows)
            if inverse is None:
                continue
            tail_input = {"n": 4, "gates": [
                *({"gate": "cx", "control": c, "target": t}
                  for c, t in inverse), *tail]}
            optimized = transpile(to_qiskit(tail_input),
                                  basis_gates=["u", "cx"],
                                  optimization_level=3, seed_transpiler=42)
            candidate = {"n": 4, "gates": base.prefix_gates(edges)
                         + from_qiskit(optimized)["gates"]}
            check = evaluate(candidate)
            evaluated += 1
            record = {"endpoint_rows": rows,
                      "inverse_cnot_count": len(inverse),
                      "tail_cnot_count": optimized.count_ops().get("cx", 0),
                      "total_cnot_count": check["cnot_count"],
                      "off_diagonal_error": check["off_diagonal_error"],
                      "valid": check["valid_diagonalizer"]}
            if best is None or (record["total_cnot_count"],
                                record["off_diagonal_error"]) < (
                                    best["total_cnot_count"],
                                    best["off_diagonal_error"]):
                best, best_candidate = record, candidate
        if best_candidate:
            path = ROOT / f"algebraic16_to15_gauge_{cancelled}_best.json"
            path.write_text(json.dumps(best_candidate, indent=2) + "\n")
            best["candidate"] = path.name
            if overall_best is None or (best["total_cnot_count"],
                                        best["off_diagonal_error"]) < (
                                            overall_best["total_cnot_count"],
                                            overall_best["off_diagonal_error"]):
                overall_best, overall_candidate = best, path.name
        results.append({"cancelled_triple_mask": TRIPLES[cancelled],
                        "required_masks": list(phases),
                        "states_explored": states,
                        "six_cnot_endpoints": len(endpoints),
                        "tails_evaluated": evaluated,
                        "best": best})
    print(json.dumps({"results": results,
                      "overall_best": overall_best,
                      "overall_candidate": overall_candidate}, indent=2))


if __name__ == "__main__":
    main()
