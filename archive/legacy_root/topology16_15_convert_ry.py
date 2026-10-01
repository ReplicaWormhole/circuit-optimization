"""Convert exact 15-CNOT candidate from RY to RX/RZ for integer-ring audit."""

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def main():
    source = json.loads((ROOT / "topology16_15_exact_candidate.json").read_text())
    gates = []
    for gate in source["gates"]:
        if gate["gate"] == "ry":
            q = gate["qubit"]
            gates.extend([{"gate": "rz", "qubit": q, "theta": "-pi/2"},
                          {"gate": "rx", "qubit": q, "theta": gate["theta"]},
                          {"gate": "rz", "qubit": q, "theta": "pi/2"}])
        else:
            gates.append(gate)
    output = ROOT / "topology16_15_rxrz_candidate.json"
    output.write_text(json.dumps({"n": 4, "gates": gates}, indent=2) + "\n")
    print(json.dumps({"output": str(output), "gate_count": len(gates),
                      "cnot_count": sum(g["gate"] == "cx" for g in gates)}))


if __name__ == "__main__":
    main()
