import json
import unittest
from pathlib import Path

from exact_check import exact_eigenvalue_labels


ROOT = Path(__file__).parent


class ExactCheckTests(unittest.TestCase):
    def test_paper_baseline(self):
        candidate = json.loads((ROOT / "baseline_18.json").read_text())
        result = exact_eigenvalue_labels(candidate)
        self.assertTrue(result["exact_diagonalizer"])
        self.assertEqual(result["cnot_count"], 18)
        self.assertEqual(result["eigenvalue_multiplicities"], [6, 3, 4, 3])

    def test_seventeen_cnot_candidate(self):
        candidate = json.loads((ROOT / "algebraic_boundary_candidate.json").read_text())
        result = exact_eigenvalue_labels(candidate)
        self.assertTrue(result["exact_diagonalizer"])
        self.assertEqual(result["cnot_count"], 17)
        self.assertEqual(result["eigenvalue_multiplicities"], [6, 3, 4, 3])

    def test_identity_is_not_diagonalizer(self):
        result = exact_eigenvalue_labels({"n": 4, "gates": []})
        self.assertFalse(result["exact_diagonalizer"])

    def test_decimal_angles_are_not_claimed_exact(self):
        with self.assertRaisesRegex(ValueError, "requires angle strings"):
            exact_eigenvalue_labels({"n": 4, "gates": [
                {"gate": "rz", "qubit": 0, "theta": 0.1}]})


if __name__ == "__main__":
    unittest.main()
