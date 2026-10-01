"""Bounded PyZX simplification of the exact 15-CNOT V4 circuit.

All outputs are expanded to one-qubit rotations and CNOTs, then checked
against V4. This tests a different search family from parity-prefix BFS.
"""

import json
from fractions import Fraction
from pathlib import Path

import pyzx as zx
from pyzx.circuit import gates

from check_circuit import evaluate
from exact_check import rational_pi


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "topology16_15_exact_candidate.json"


def to_pyzx(candidate):
    circuit = zx.Circuit(4)
    for gate in candidate["gates"]:
        name = gate["gate"].lower()
        if name == "cx":
            circuit.add_gate(gates.CNOT(gate["control"], gate["target"]))
        elif name == "h":
            circuit.add_gate(gates.HAD(gate["qubit"]))
        elif name in ("rx", "ry", "rz"):
            q = gate["qubit"]
            angle = rational_pi(gate["theta"])
            if name == "ry":
                circuit.add_gate(gates.ZPhase(q, Fraction(-1, 2)))
                circuit.add_gate(gates.HAD(q))
                circuit.add_gate(gates.ZPhase(q, angle))
                circuit.add_gate(gates.HAD(q))
                circuit.add_gate(gates.ZPhase(q, Fraction(1, 2)))
            elif name == "rx":
                circuit.add_gate(gates.HAD(q))
                circuit.add_gate(gates.ZPhase(q, angle))
                circuit.add_gate(gates.HAD(q))
            else:
                circuit.add_gate(gates.ZPhase(q, angle))
        else:
            raise ValueError(f"unsupported input gate: {name}")
    return circuit


def angle(phase):
    value = Fraction(phase)
    numerator, denominator = value.numerator, value.denominator
    if numerator == 0:
        return "0*pi"
    if denominator == 1:
        return f"{numerator}*pi"
    return f"{numerator}*pi/{denominator}"


def from_pyzx(circuit):
    output = []
    for gate in circuit.gates:
        name = gate.name
        if name == "CNOT":
            output.append({"gate": "cx", "control": gate.control,
                           "target": gate.target})
        elif name == "CZ":
            output.extend([{"gate": "h", "qubit": gate.target},
                           {"gate": "cx", "control": gate.control,
                            "target": gate.target},
                           {"gate": "h", "qubit": gate.target}])
        elif name == "SWAP":
            a, b = gate.control, gate.target
            output.extend([{"gate": "cx", "control": a, "target": b},
                           {"gate": "cx", "control": b, "target": a},
                           {"gate": "cx", "control": a, "target": b}])
        elif name == "HAD":
            output.append({"gate": "h", "qubit": gate.target})
        elif name == "ZPhase":
            output.append({"gate": "rz", "qubit": gate.target,
                           "theta": angle(gate.phase)})
        elif name == "XPhase":
            output.append({"gate": "rx", "qubit": gate.target,
                           "theta": angle(gate.phase)})
        else:
            raise ValueError(f"unsupported output gate: {name}")
    return {"n": 4, "gates": output}


def main():
    source = json.loads(SOURCE.read_text())
    circuit = to_pyzx(source)
    trials = [
        ("basic", lambda: zx.optimize.basic_optimization(circuit.copy())),
        ("full", lambda: zx.optimize.full_optimize(circuit.copy())),
        ("zx", lambda: zx.extract_circuit(
            reduced_graph(circuit), optimize_cnots=2)),
    ]
    results = []
    for label, method in trials:
        try:
            candidate = from_pyzx(method())
            check = evaluate(candidate)
            path = ROOT / f"pyzx14_{label}.json"
            path.write_text(json.dumps(candidate, indent=2) + "\n")
            results.append({"method": label, "candidate": path.name,
                            "cnot_count": check["cnot_count"],
                            "off_diagonal_error": check["off_diagonal_error"],
                            "valid": check["valid_diagonalizer"]})
        except Exception as error:
            results.append({"method": label, "error": repr(error)})
    print(json.dumps(results, indent=2))


def reduced_graph(circuit):
    graph = circuit.to_graph()
    zx.simplify.full_reduce(graph, quiet=True)
    return graph


if __name__ == "__main__":
    main()
