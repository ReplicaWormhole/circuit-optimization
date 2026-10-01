import unittest
import json
from pathlib import Path

from check_circuit import angle, evaluate


class CircuitCheckerTests(unittest.TestCase):
    def test_known_one_cnot_two_qubit_diagonalizer(self):
        result = evaluate({"n": 2, "gates": [
            {"gate": "cx", "control": 0, "target": 1},
            {"gate": "h", "qubit": 0},
        ]})
        self.assertEqual(result["cnot_count"], 1)
        self.assertTrue(result["valid_diagonalizer"])

    def test_identity_does_not_diagonalize_four_cycle(self):
        result = evaluate({"n": 4, "gates": []})
        self.assertFalse(result["valid_diagonalizer"])
        self.assertFalse(result["beats_18_cnot_reference"])

    def test_reconstructed_paper_baseline(self):
        baseline = json.loads((Path(__file__).parent / "baseline_18.json").read_text())
        result = evaluate(baseline)
        self.assertEqual(result["cnot_count"], 18)
        self.assertTrue(result["valid_diagonalizer"])
        self.assertLess(result["off_diagonal_error"], 1e-12)
        self.assertFalse(result["diagonalizes_total_spin"])

    def test_angle_and_gate_validation(self):
        self.assertAlmostEqual(angle("-3*pi/4"), -3.0 * angle("pi") / 4)
        with self.assertRaises(ValueError):
            evaluate({"n": 4, "gates": [{"gate": "cx", "control": 0, "target": 0}]})
        with self.assertRaises(ValueError):
            evaluate({"n": 4, "gates": [{"gate": "swap", "control": 0, "target": 1}]})


if __name__ == "__main__":
    unittest.main()
