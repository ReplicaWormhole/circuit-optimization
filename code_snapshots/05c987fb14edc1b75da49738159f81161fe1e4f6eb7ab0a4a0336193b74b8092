"""Convert exact 16-CNOT circuit to RX/RZ/H for independent ring audit.

Uses RY(theta) = RZ(pi/2) RX(theta) RZ(-pi/2); the gate list is chronological.
"""

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def main():
    candidate = json.loads((ROOT / "topology16_exact_angle_candidate.json").read_text())
    gates = []
    for gate in candidate["gates"]:
        if gate["gate"] == "ry":
            q = gate["qubit"]
            gates.extend([{"gate": "rz", "qubit": q, "theta": "-pi/2"},
                          {"gate": "rx", "qubit": q, "theta": gate["theta"]},
                          {"gate": "rz", "qubit": q, "theta": "pi/2"}])
        else:
            gates.append(gate)
    output = ROOT / "topology16_rxrz_candidate.json"
    output.write_text(json.dumps({"n": 4, "gates": gates}, indent=2) + "\n")
    print(json.dumps({"output": str(output), "gate_count": len(gates),
                      "cnot_count": sum(g["gate"] == "cx" for g in gates)}))


if __name__ == "__main__":
    main()
