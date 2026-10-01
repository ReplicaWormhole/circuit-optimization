"""Three PyZX optimizer modes applied directly to the exact14 gate list."""

import json
from pathlib import Path

import pyzx as zx

from check_circuit import evaluate
from pyzx_search14 import from_pyzx, reduced_graph, to_pyzx


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "topology14_exact_matchgate_rational.json"


def main():
    source = json.loads(SOURCE.read_text())
    circuit = to_pyzx(source)
    modes = {
        "basic": lambda: zx.optimize.basic_optimization(circuit.copy()),
        "full": lambda: zx.optimize.full_optimize(circuit.copy()),
        "zx": lambda: zx.extract_circuit(reduced_graph(circuit),
                                         optimize_cnots=2),
    }
    rows = []
    best = None
    for name, build in modes.items():
        try:
            candidate = from_pyzx(build())
            check = evaluate(candidate)
            row = {"mode": name, "cnot_count": check["cnot_count"],
                   "off_diagonal_error": check["off_diagonal_error"],
                   "valid_diagonalizer": check["valid_diagonalizer"]}
            if check["valid_diagonalizer"] and (
                    best is None or row["cnot_count"] < best[0]["cnot_count"]):
                best = row, candidate
        except Exception as error:
            row = {"mode": name, "error": repr(error)}
        rows.append(row)
    if best:
        (ROOT / "search13_pyzx_exact14_best.json").write_text(
            json.dumps(best[1], indent=2) + "\n")
    print(json.dumps({"source": SOURCE.name, "rows": rows,
                      "best": best[0] if best else None}, indent=2))


if __name__ == "__main__":
    main()
