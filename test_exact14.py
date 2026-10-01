"""Regression checks for the exact four-qubit 14-CNOT certificate."""

import copy
import json
import unittest
from pathlib import Path

from delete14_integer_audit import check as integer_audit
from exact_check import exact_eigenvalue_labels


ROOT = Path(__file__).resolve().parent
CANDIDATE = ROOT / "topology14_exact_matchgate_rational.json"


class ExactFourteenCnotTest(unittest.TestCase):
    def test_two_exact_arithmetics_agree(self):
        candidate = json.loads(CANDIDATE.read_text())
        field = exact_eigenvalue_labels(candidate)
        integers = integer_audit(candidate)
        self.assertTrue(field["exact_diagonalizer"])
        self.assertTrue(integers["exact_uv_equals_du"])
        self.assertEqual(field["cnot_count"], 14)
        self.assertEqual(integers["cnot_count"], 14)
        self.assertEqual(field["cyclotomic_conductor"], 96)
        names = ("1", "i", "-1", "-i")
        self.assertEqual([names[i] for i in field["output_labels"]],
                         integers["eigenvalue_labels"])

    def test_angle_change_breaks_certificate(self):
        candidate = copy.deepcopy(json.loads(CANDIDATE.read_text()))
        candidate["gates"][0]["theta"] = "3*pi/8"
        self.assertFalse(exact_eigenvalue_labels(candidate)["exact_diagonalizer"])
        self.assertFalse(integer_audit(candidate)["exact_uv_equals_du"])


if __name__ == "__main__":
    unittest.main()
