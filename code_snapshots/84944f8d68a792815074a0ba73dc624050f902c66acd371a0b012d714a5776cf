"""Bounded PyZX reductions of the exact 14-CNOT cycle diagonalizer.

These rewrites preserve the full unitary up to global phase. Each extracted
candidate is independently checked against the right shift; a shorter valid
candidate would be checked exactly before any mathematical claim.
"""

import json
from pathlib import Path

import pyzx as zx

from check_circuit import evaluate
from pyzx_search14 import from_pyzx, to_pyzx


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "topology14_exact_matchgate_rational.json"
OUT = ROOT / "pyzx_exact14_scan_result.json"


def graph_extract(circuit, optimize_cnots):
    graph = circuit.to_graph()
    zx.simplify.full_reduce(graph, quiet=True)
    return zx.extract_circuit(graph, optimize_cnots=optimize_cnots)


def main():
    source = json.loads(SOURCE.read_text())
    circuit = to_pyzx(source)
    rows = []
    methods = [
        ("basic", lambda: zx.optimize.basic_optimization(circuit.copy())),
        ("full", lambda: zx.optimize.full_optimize(circuit.copy())),
        *[(f"extract_{level}", lambda level=level: graph_extract(circuit, level))
          for level in range(4)],
    ]
    for name, method in methods:
        try:
            candidate = from_pyzx(method())
            check = evaluate(candidate)
            row = {"method": name, "cnot_count": check["cnot_count"],
                   "gate_count": len(candidate["gates"]),
                   "valid": bool(check["valid_diagonalizer"]),
                   "max_offdiagonal": float(check["off_diagonal_error"])}
            if row["cnot_count"] <= 14:
                path = ROOT / f"pyzx_exact14_{name}.json"
                path.write_text(json.dumps(candidate, indent=2) + "\n")
                row["candidate"] = path.name
            rows.append(row)
        except Exception as error:
            rows.append({"method": name, "error": repr(error)})
    result = {"source": SOURCE.name, "pyzx_version": zx.__version__,
              "methods": rows,
              "scope": "heuristic exact-unitary rewrites; no lower bound"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
