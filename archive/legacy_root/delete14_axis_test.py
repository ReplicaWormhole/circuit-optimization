"""Regression checks for the axis-angle SU2 search serialization."""

import unittest

import numpy as np
import torch

from check_circuit import one_qubit_matrix
import delete14_search
import delete14_tail_search


class AxisAngleTest(unittest.TestCase):
    def test_zero_axis_serializes_to_identity(self):
        for module in (delete14_search, delete14_tail_search):
            with self.subTest(module=module.__name__):
                gate = {"gate": "u3", **module.axis_to_u3(0.0, 0.0, 0.0)}
                self.assertTrue(np.allclose(one_qubit_matrix(gate), np.eye(2),
                                            atol=1e-14))

    def test_zero_axis_gradient_is_finite(self):
        for module in (delete14_search, delete14_tail_search):
            with self.subTest(module=module.__name__):
                vector = torch.zeros(3, dtype=torch.float64, requires_grad=True)
                matrix = module.axis_rotation(*vector)
                matrix[0, 1].imag.backward()
                self.assertTrue(torch.isfinite(vector.grad).all())

    def test_nonzero_axis_matches_u3_up_to_global_phase(self):
        for module in (delete14_search, delete14_tail_search):
            for vector in ((0.1, -0.2, 0.3), (-2.1, 0.7, 1.3),
                           (0.0, 0.0, 0.4)):
                with self.subTest(module=module.__name__, vector=vector):
                    x, y, z = vector
                    length = np.linalg.norm(vector)
                    scale = np.sin(length / 2) / length
                    expected = np.array([
                        [np.cos(length / 2) - 1j * z * scale,
                         (-1j * x - y) * scale],
                        [(-1j * x + y) * scale,
                         np.cos(length / 2) + 1j * z * scale]])
                    actual = one_qubit_matrix(
                        {"gate": "u3", **module.axis_to_u3(*vector)})
                    phase = expected[0, 0] / actual[0, 0]
                    self.assertTrue(np.allclose(phase * actual, expected,
                                                atol=1e-13))


if __name__ == "__main__":
    unittest.main()
