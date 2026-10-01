"""Build a 15-CNOT candidate with two exact modified Fourier butterflies.

The two-qubit identity used twice is symbolically certified by
algebraic16_exact.py. The special one-qubit gate has the exact matrix
[[1, sqrt(2)+i], [-sqrt(2)+i, 1]]/2; a companion script checks the full
four-qubit cycle identity over Q(zeta_32) from the exact-spec JSON.
"""

import json
from pathlib import Path

from algebraic16_exact import PHI, REPLACEMENT, certify_replacement
from check_circuit import evaluate


ROOT = Path(__file__).parent
PREFIX = [
    {"gate": "rz", "qubit": 1, "theta": "5*pi/8"},
    {"gate": "rz", "qubit": 2, "theta": "3*pi/4"},
    {"gate": "rz", "qubit": 3, "theta": "3*pi/8"},
    {"gate": "cx", "control": 3, "target": 1},
    {"gate": "cx", "control": 1, "target": 2},
    {"gate": "rz", "qubit": 2, "theta": "-pi/4"},
    {"gate": "cx", "control": 3, "target": 2},
    {"gate": "rz", "qubit": 2, "theta": "-pi/2"},
    {"gate": "cx", "control": 0, "target": 2},
    {"gate": "rz", "qubit": 2, "theta": "pi/8"},
    {"gate": "cx", "control": 1, "target": 2},
    {"gate": "rz", "qubit": 2, "theta": "-pi/8"},
    {"gate": "cx", "control": 3, "target": 2},
]

SPECIAL_MATRIX = [["1/2", "(sqrt(2)+i)/2"],
                  ["(-sqrt(2)+i)/2", "1/2"]]


def relabel(gates, wire_map):
    result = []
    for gate in gates:
        item = gate.copy()
        if item["gate"] == "cx":
            item["control"] = wire_map[item["control"]]
            item["target"] = wire_map[item["target"]]
        else:
            item["qubit"] = wire_map[item["qubit"]]
        result.append(item)
    return result


def numeric_candidate():
    baseline = json.loads((ROOT / "baseline_18.json").read_text())
    gates = PREFIX + REPLACEMENT + relabel(REPLACEMENT, {0: 1, 2: 3})
    gates += baseline["gates"][38:]
    return {"n": 4, "gates": gates}


def exact_spec(candidate):
    gates = []
    special_count = 0
    for gate in candidate["gates"]:
        if gate["gate"] == "u3" and gate.get("phi") == PHI:
            gates.append({"gate": "u3_algebraic", "qubit": gate["qubit"],
                          "theta": "2*pi/3",
                          "phi": "pi-atan(1/sqrt(2))",
                          "lam": "-pi+atan(1/sqrt(2))",
                          "matrix": SPECIAL_MATRIX})
            special_count += 1
        else:
            gates.append({key: ("0*pi" if key in ("theta", "phi", "lam")
                                 and value == 0 else value)
                          for key, value in gate.items()})
    assert special_count == 2
    return {"n": 4, "gates": gates}


def main():
    numeric = numeric_candidate()
    num_path = ROOT / "algebraic16_to15_exact_candidate.json"
    num_path.write_text(json.dumps(numeric, indent=2) + "\n")
    spec = exact_spec(numeric)
    spec_path = ROOT / "algebraic16_to15_exact_spec.json"
    spec_path.write_text(json.dumps(spec, indent=2) + "\n")
    block_valid, global_phase = certify_replacement()
    print(json.dumps({"numeric_candidate": num_path.name,
                      "exact_spec": spec_path.name,
                      "two_qubit_identity_exact": block_valid,
                      "two_qubit_global_phase": str(global_phase),
                      "numeric_check": {k: v for k, v in evaluate(numeric).items()
                                        if k != "diagonal_eigenvalues"}}, indent=2))
    if not block_valid:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
