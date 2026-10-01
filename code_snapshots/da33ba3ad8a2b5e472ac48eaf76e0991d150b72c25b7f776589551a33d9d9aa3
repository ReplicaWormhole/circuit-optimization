"""Search five-CNOT prefixes for the alternate CZ pair-mask gauge.

The modular phase classification found four support sets with pair mask 9
instead of the baseline's mask 6. For each, enumerate five-CNOT open parity
prefixes and synthesize the inverse endpoint plus Fourier tail with Qiskit.
"""

import json
from pathlib import Path

from qiskit import transpile

from algebraic14_five_prefix import BASIS, enumerate_prefixes, inverse_paths
from check_circuit import evaluate
from qiskit_crosscheck import from_qiskit, to_qiskit


ROOT = Path(__file__).parent


def prefix_gates(edges, coefficients):
    gates = []
    for qubit in range(4):
        coefficient = coefficients.get(str(8 >> qubit), 0)
        if coefficient:
            gates.append({"gate": "rz", "qubit": qubit,
                          "theta": f"{coefficient}*pi/8"})
    rows = list(BASIS)
    seen = set()
    required = {int(mask) for mask in coefficients
                if int(mask).bit_count() >= 2}
    for control, target in edges:
        gates.append({"gate": "cx", "control": control, "target": target})
        rows[target] ^= rows[control]
        mask = rows[target]
        if mask in required and mask not in seen:
            gates.append({"gate": "rz", "qubit": target,
                          "theta": f"{coefficients[str(mask)]}*pi/8"})
            seen.add(mask)
    assert seen == required
    return gates


def main():
    records = json.loads((ROOT / "algebraic14_fourmask_models.json").read_text())
    models = [record for record in records if 9 in record["support"]]
    assert len(models) == 4
    baseline = json.loads((ROOT / "baseline_18.json").read_text())
    tail = baseline["gates"][18:]
    inverses = inverse_paths(4)
    results = []
    overall = None
    overall_path = None
    for index, model in enumerate(models):
        support = tuple(model["support"])
        coefficients = model["coefficients_pi_over_8"]
        ends, states = enumerate_prefixes(support)
        evaluated = 0
        best = None
        best_candidate = None
        distance_histogram = {}
        for rows, edges in ends:
            inverse = inverses.get(rows)
            distance = len(inverse) if inverse is not None else ">4"
            distance_histogram[str(distance)] = (
                distance_histogram.get(str(distance), 0) + 1)
            if inverse is None:
                continue
            modified_tail = {"n": 4, "gates": [
                *({"gate": "cx", "control": c, "target": t}
                  for c, t in inverse), *tail]}
            optimized = transpile(to_qiskit(modified_tail),
                                  basis_gates=["u", "cx"],
                                  optimization_level=3, seed_transpiler=42)
            candidate = {"n": 4, "gates": prefix_gates(edges, coefficients)
                         + from_qiskit(optimized)["gates"]}
            check = evaluate(candidate)
            evaluated += 1
            result = {"endpoint_rows": rows,
                      "inverse_distance": len(inverse),
                      "tail_cnot_count": optimized.count_ops().get("cx", 0),
                      "total_cnot_count": check["cnot_count"],
                      "off_diagonal_error": check["off_diagonal_error"],
                      "valid": check["valid_diagonalizer"]}
            if best is None or (result["total_cnot_count"],
                                result["off_diagonal_error"]) < (
                                    best["total_cnot_count"],
                                    best["off_diagonal_error"]):
                best, best_candidate = result, candidate
        if best_candidate:
            path = ROOT / f"algebraic14_alt_pair_{index}_best.json"
            path.write_text(json.dumps(best_candidate, indent=2) + "\n")
            best["candidate"] = path.name
            if overall is None or (best["total_cnot_count"],
                                   best["off_diagonal_error"]) < (
                                       overall["total_cnot_count"],
                                       overall["off_diagonal_error"]):
                overall, overall_path = best, path.name
        results.append({"support": support,
                        "coefficients_pi_over_8": coefficients,
                        "states_explored": states,
                        "five_cnot_endpoints": len(ends),
                        "inverse_distance_histogram": distance_histogram,
                        "tails_evaluated": evaluated,
                        "best": best})
    print(json.dumps({"results": results,
                      "overall_best": overall,
                      "overall_candidate": overall_path}, indent=2))


if __name__ == "__main__":
    main()
