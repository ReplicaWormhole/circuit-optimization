"""Focused native matrix, conventions, and baseline regression tests."""
import copy
import json
import hashlib
import math
from pathlib import Path
import unittest

import numpy as np

from native_gate_check import circuit_matrix, embedded_two_qubit, evaluate, native_matrix, phase_error

ROOT = Path(__file__).resolve().parent


def compiled_f(a, b, first=0, second=1, n=2):
    return {"n": n, "gates": [
        {"gate": "rx", "qubit": first, "theta": -math.pi / 2},
        {"gate": "rx", "qubit": second, "theta": -math.pi / 2},
        {"gate": "cx", "control": first, "target": second},
        {"gate": "rx", "qubit": first, "theta": -2 * a},
        {"gate": "rz", "qubit": second, "theta": -2 * b},
        {"gate": "cx", "control": first, "target": second},
        {"gate": "rx", "qubit": first, "theta": math.pi / 2},
        {"gate": "rx", "qubit": second, "theta": math.pi / 2}]}


class NativeGateTest(unittest.TestCase):
    def test_fixed_exchange_matrices(self):
        iswap = native_matrix({"gate": "iswap"})
        expected = np.array([[1, 0, 0, 0], [0, 0, 1j, 0], [0, 1j, 0, 0], [0, 0, 0, 1]])
        np.testing.assert_allclose(iswap, expected, atol=1e-15)
        root = native_matrix({"gate": "sqrt_iswap"})
        np.testing.assert_allclose(root @ root, expected, atol=1e-15)

    def test_two_cnot_identity(self):
        for a,b in [(math.pi/4, math.pi/8), (.23,-.71), (0,0)]:
            for first,second in [(0,2),(2,0),(1,2)]:
                with self.subTest(a=a,b=b,first=first,second=second):
                    native = {"n": 3, "gates": [{"gate":"xx_yy","qubits":[first,second],"a":a,"b":b}]}
                    self.assertLess(phase_error(circuit_matrix(native), circuit_matrix(compiled_f(a,b,first,second,3))), 1e-14)

    def test_embedding_msb_and_order(self):
        # A deliberately asymmetric local matrix tests ordered wire semantics.
        local_cx = np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
        forward = embedded_two_qubit(local_cx, 0, 2, 3)
        reverse = embedded_two_qubit(local_cx, 2, 0, 3)
        self.assertEqual(np.argmax(forward[:,4]),5) # |100> -> |101>
        self.assertEqual(np.argmax(reverse[:,1]),5) # |001> -> |101>
        self.assertEqual(np.argmax(forward[:,1]),1)

    def test_cp_compilation_and_phase(self):
        theta=.73
        native={"n":2,"gates":[{"gate":"cp","qubits":[0,1],"theta":theta}]}
        compiled={"n":2,"gates":[{"gate":"rz","qubit":0,"theta":theta/2},{"gate":"rz","qubit":1,"theta":theta/2},{"gate":"cx","control":0,"target":1},{"gate":"rz","qubit":1,"theta":-theta/2},{"gate":"cx","control":0,"target":1}]}
        self.assertLess(phase_error(circuit_matrix(native),circuit_matrix(compiled)),1e-14)

    def test_invalid_gates_and_angles(self):
        good={"gate":"xx_yy","qubits":[0,1],"a":0,"b":0}
        bad=[dict(good,qubits=[0,0]),dict(good,qubits=[0,2]),dict(good,qubits=[True,1]),dict(good,qubits=[0]),dict(good,a=float('nan')),dict(good,b=float('inf')),dict(good,a='pi/0'),dict(good,a=True),{"gate":"cp","qubits":[0,1],"theta":'bad'}, {"gate":"rz","qubit":0,"theta":float('inf')}, {"gate":"unknown"}, None]
        for gate in bad:
            with self.subTest(gate=gate), self.assertRaises(ValueError):
                circuit_matrix({"n":2,"gates":[gate]})
        for n in [True,1,9,2.0]:
            with self.assertRaises(ValueError): circuit_matrix({"n":n,"gates":[]})
        for tol in [0,-1,float('nan'),float('inf'),True]:
            with self.assertRaises(ValueError): evaluate({"n":2,"gates":[]},tol)

    def test_native_baseline(self):
        baseline=json.loads((ROOT/'collaboration/native_baseline.json').read_text())
        source=ROOT/'topology14_exact_matchgate_rational.json'
        reference=json.loads(source.read_text())
        self.assertEqual(baseline['provenance']['compiled_sha256'], hashlib.sha256(source.read_bytes()).hexdigest())
        result=evaluate(baseline,reference=reference)
        self.assertTrue(result['valid_diagonalizer'])
        self.assertEqual(result['two_qubit_count'],10)
        self.assertEqual(result['native_family_gate_count'],4)
        self.assertEqual(result['cnot_compilation_upper_bound'],14)
        self.assertIsNone(result['minimum_cnot_cost'])
        self.assertLess(result['reference_matrix_error_up_to_phase'],1e-14)
        altered=copy.deepcopy(baseline)
        altered['gates'][23]['a']='pi/8'
        self.assertFalse(evaluate(altered)['valid_diagonalizer'])

    def test_parallel_depth(self):
        candidate={"n":4,"gates":[{"gate":"iswap","qubits":[0,1]}, {"gate":"iswap","qubits":[2,3]}, {"gate":"cp","qubits":[1,2],"theta":.2}]}
        self.assertEqual(evaluate(candidate)['two_qubit_depth'],2)


if __name__ == '__main__': unittest.main()
