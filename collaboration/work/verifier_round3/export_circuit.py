"""Render the accepted exact gate list without numerical reconstruction."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT / "collaboration/work/algebraic13_new/canonical_algebraic_ansatz.json"
OUTPUT = Path(__file__).with_name("EXPLICIT_CIRCUIT.md")

def render():
    raw = SOURCE.read_bytes()
    data = json.loads(raw)
    acceptance = json.loads((ROOT / "collaboration/EXACT13_ACCEPTANCE.json").read_text())
    digest = hashlib.sha256(raw).hexdigest()
    assert digest == acceptance["sha256"]
    gates = data["gates"]
    assert len(gates) == 73
    assert sum(g["gate"] == "cx" for g in gates) == 13
    assert sum(g["gate"] in ("rx", "ry", "rz") for g in gates) == 60
    lines = ["# Explicit accepted 13-CNOT diagonalizer", "",
        "Current incumbent, not a proven optimum. The retained bound is 6 <= C_min <= 13.", "",
        "Gate model: ancilla-free four qubits, arbitrary exact single-qubit rotations and fixed directed CNOT on any ordered pair. Other native entanglers have separate costs.", "",
        "Qubit 0 is most significant. Execute rows in ascending order; U = G_73 ... G_1. The target is U V4 = D U, where V4|x0 x1 x2 x3> = |x3 x0 x1 x2>. Output eigenvalue order is arbitrary.", "",
        "R_a(theta) = cos(theta/2) I - i sin(theta/2) sigma_a. The signed half-angle columns specify the exact rotations, including their branches; theta/pi is shown only when supplied exactly. Every radical sqrt denotes the nonnegative real root. Full-angle columns preserve the source metadata. Zero rotations are retained.", "",
        "Source: `collaboration/work/algebraic13_new/canonical_algebraic_ansatz.json`", "",
        "SHA256: `" + digest + "`", "",
        "Existing acceptance: `collaboration/EXACT13_ACCEPTANCE.json`, certificate runs 370 and 371 using the same full-matrix implementation, with separate relation and scalar audits. This export checks transcription only; it does not repeat or add an exact matrix certificate.", "",
        "| Step | Gate | Wire(s) | theta/pi | cos(theta/2) | sin(theta/2) | cos(theta) | sin(theta) | Coordinate |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for step, g in enumerate(gates, 1):
        if g["gate"] == "cx":
            values = [str(step), "CX", str(g["control"]) + " -> " + str(g["target"])] + ["—"] * 6
        else:
            values = [str(step), g["gate"].upper(), str(g["qubit"]), g.get("rational_pi", "—"), g["cos_half_exact"], g["sin_half_exact"], g["cos_exact"], g["sin_exact"], str(g["coordinate"])]
        lines.append("| " + " | ".join(values) + " |")
    lines += ["", "D output labels (basis states 0000 through 1111, in binary order). Labels 0, 1, 2, 3 mean 1, i, -1, -i, respectively; D_rr = i**labels[r]:", "", "`" + json.dumps(data["labels"]) + "`", "", "Reproduce: `python3 collaboration/work/verifier_round3/export_circuit.py`", ""]
    return "\n".join(lines)

if __name__ == "__main__":
    OUTPUT.write_text(render())
    print(OUTPUT.relative_to(ROOT))
