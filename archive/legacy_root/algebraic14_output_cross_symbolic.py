"""Exact all-angle rank certificate for cross-pair output RZZ gauges.

Let A be the last two-CNOT butterfly on wires (0,1), and B its copy on
(2,3). For q in one pair and r in the other, the four row blocks of
RZZ(q,r;delta)(A tensor B), realigned across q|rest, have the form

    A[a,b] tensor D[a] B,  a,b in {0,1},

where D[a] is an invertible diagonal gate on r. Exact rank 4 of the 4x4
realignment of A proves each two-row group is independent. If D[0] and
D[1] are independent, the groups are independent by their second tensor
factor; if proportional, rank 4 of A proves the same. Thus each one-wire
Schmidt rank is 4 for *every* real delta, with no angle grid assumption.
Four CNOTs are consequently necessary for this fixed four-wire block.
"""

import itertools
import json
from pathlib import Path

import sympy as sp

from exact_check import rational_pi


ROOT = Path(__file__).parent
FIELD = sp.QQ.cyclotomic_field(32)
ZETA = FIELD.unit
ZERO, ONE = FIELD.zero, FIELD.one
HALF = FIELD.from_sympy(sp.Rational(1, 2))
IMAGINARY = ZETA ** 8


def phase(angle):
    fraction = rational_pi(angle)
    exponent = fraction.numerator * 32 // (4 * fraction.denominator)
    assert exponent * 4 * fraction.denominator == 32 * fraction.numerator
    return ZETA ** exponent


def one_wire(gate):
    positive = phase(gate["theta"])
    negative = positive ** -1
    cosine = (positive + negative) * HALF
    sine = (positive - negative) / (2 * IMAGINARY)
    if gate["gate"] == "rx":
        return ((cosine, -IMAGINARY * sine),
                (-IMAGINARY * sine, cosine))
    if gate["gate"] == "rz":
        return ((negative, ZERO), (ZERO, positive))
    raise ValueError(gate)


def butterfly():
    gates = json.loads((ROOT / "baseline_18.json").read_text())["gates"][42:52]
    matrix = [[ONE if row == column else ZERO for column in range(4)]
              for row in range(4)]
    for gate in gates:
        if gate["gate"] == "cx":
            assert (gate["control"], gate["target"]) == (0, 1)
            matrix[2], matrix[3] = matrix[3], matrix[2]
        else:
            local = one_wire(gate)
            mask = 2 >> gate["qubit"]
            for row in range(4):
                if row & mask:
                    continue
                partner = row | mask
                first, second = matrix[row], matrix[partner]
                matrix[row] = [local[0][0] * x + local[0][1] * y
                               for x, y in zip(first, second)]
                matrix[partner] = [local[1][0] * x + local[1][1] * y
                                   for x, y in zip(first, second)]
    return matrix


def realignment(matrix, qubit):
    result = []
    for output_bit in range(2):
        for input_bit in range(2):
            row = []
            for output_other in range(2):
                for input_other in range(2):
                    output = (output_bit << (1 - qubit)) | (output_other << qubit)
                    source = (input_bit << (1 - qubit)) | (input_other << qubit)
                    row.append(matrix[output][source])
            result.append(row)
    return result


def determinant(matrix):
    total = ZERO
    for permutation in itertools.permutations(range(4)):
        inversions = sum(permutation[i] > permutation[j]
                         for i in range(4) for j in range(i + 1, 4))
        term = ONE
        for row, column in enumerate(permutation):
            term *= matrix[row][column]
        total += -term if inversions % 2 else term
    return total


def main():
    matrix = butterfly()
    cases = []
    for qubit in range(2):
        value = determinant(realignment(matrix, qubit))
        cases.append({"qubit_in_butterfly": qubit,
                      "determinant": str(FIELD.to_sympy(value)),
                      "nonzero_exact": value != ZERO})
    result = {"cyclotomic_conductor": 32,
              "butterfly_schmidt_determinants": cases,
              "all_nonzero": all(case["nonzero_exact"] for case in cases),
              "all_angle_conclusion":
                  "Every cross-pair RZZ output gauge has rank 4 across each one-wire cut, hence at least 4 CNOTs in the fixed final block."}
    print(json.dumps(result, indent=2))
    if not result["all_nonzero"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
