"""Exact selected-axis witnesses for the two gate-3 shapes among 40 survivors.

These are only partial necessary-condition witnesses, not diagonalizers.
All CNOT/Pauli/shift matrices have integer entries. Spectral projector
numerators have exact Gaussian-integer entries represented in complex128.
"""

from collections import Counter
from itertools import permutations, product
import json

import numpy as np

from global_lower_bound6_terminal_filter import (
    ADJACENT, EDGES, degrees, full_forward_support, terminal_obstructions,
)
from global_lower_bound6_two_terminal_axes import terminal_opposite_degree_two_wires
from global_lower_bound6_chain_certificate import chain_match
from global_lower_bound6_adjacent_v_certificate import adjacent_v_match
from global_lower_bound6_opposite_chain_certificate import opposite_chain_match


def cnot(control, target):
    gate = np.zeros((16, 16), dtype=np.int64)
    target_bit = 1 << (3 - target)
    control_bit = 1 << (3 - control)
    for state in range(16):
        out = state ^ target_bit if state & control_bit else state
        gate[out, state] = 1
    return gate


def pauli(letter, wire):
    one = {
        'I': np.eye(2, dtype=np.int64),
        'X': np.array([[0, 1], [1, 0]], dtype=np.int64),
        'Z': np.diag([1, -1]).astype(np.int64),
    }
    out = np.array([[1]], dtype=np.int64)
    for j in range(4):
        out = np.kron(out, one[letter] if j == wire else one['I'])
    return out


def compose_circuit(oriented):
    out = np.eye(16, dtype=np.int64)
    for control, target in oriented:
        out = cnot(control, target) @ out
    return out


def shift_matrix():
    shift = np.zeros((16, 16), dtype=np.int64)
    for state in range(16):
        shift[((state & 1) << 3) | (state >> 1), state] = 1
    return shift


def spectral_defect_zero(axis, shift):
    powers = [np.linalg.matrix_power(shift, r).astype(np.complex128) for r in range(4)]
    numerators = [sum(((-1j) ** (k * r) * powers[r] for r in range(4)),
                      np.zeros((16, 16), dtype=np.complex128)) for k in range(4)]
    return all(np.max(np.abs(numerators[a] @ axis @ numerators[b] @ axis @ numerators[c])) == 0
               for a, b, c in permutations(range(4), 3))


def classify_remaining():
    count = Counter()
    examples = {}
    for schedule in product(EDGES, repeat=5):
        deg = degrees(schedule)
        if min(deg) < 2 or not full_forward_support(schedule):
            continue
        if any(terminal_obstructions(schedule, deg)):
            continue
        if len(terminal_opposite_degree_two_wires(schedule, deg)) >= 2:
            continue
        if chain_match(schedule, deg) or adjacent_v_match(schedule, deg):
            continue
        if opposite_chain_match(schedule, deg):
            continue
        lows = [j for j in range(4) if deg[j] == 2]
        early = [j for j in lows if max(t for t, edge in enumerate(schedule) if j in edge) == 2]
        assert len(early) == 1
        j = early[0]
        k = next(wire for wire in schedule[2] if wire != j)
        shape = 'adjacent_star' if tuple(sorted((j, k))) in ADJACENT else 'opposite_chain'
        count[shape] += 1
        examples.setdefault(shape, [list(edge) for edge in schedule])
    assert count == Counter({'adjacent_star': 8, 'opposite_chain': 32})
    return count, examples


def verify_witness(name, oriented, selected_z, selected_x):
    unitary = compose_circuit(oriented)
    assert np.array_equal(unitary.T @ unitary, np.eye(16, dtype=np.int64))
    shift = shift_matrix()
    full_z = np.linalg.multi_dot([pauli('Z', j) for j in range(4)])
    opposite_x = pauli('X', 1) @ pauli('X', 3)
    image_z = unitary @ pauli('Z', selected_z) @ unitary.T
    image_x = unitary @ pauli('X', selected_x) @ unitary.T
    assert np.array_equal(image_z, full_z)
    assert np.array_equal(image_x, opposite_x)
    assert np.array_equal(image_z @ shift, shift @ image_z)
    assert np.array_equal(image_x @ shift @ shift, shift @ shift @ image_x)
    assert spectral_defect_zero(image_z, shift)
    assert spectral_defect_zero(image_x, shift)
    transformed_shift = unitary.T @ shift @ unitary
    offdiag_nonzero = int(np.count_nonzero(transformed_shift - np.diag(np.diag(transformed_shift))))
    assert offdiag_nonzero > 0
    return {
        'name': name,
        'oriented_cnot_schedule': [list(edge) for edge in oriented],
        'selected_input_z_wire': selected_z,
        'selected_input_x_wire': selected_x,
        'selected_z_image': 'Z0 Z1 Z2 Z3',
        'selected_x_image': 'X1 X3',
        'spectral_path_zero_for_both': True,
        'transformed_shift_offdiagonal_nonzero_entries': offdiag_nonzero,
    }


def main():
    count, examples = classify_remaining()
    witnesses = [
        verify_witness('adjacent_star', [(0, 1), (2, 0), (3, 2), (0, 3), (1, 3)], 2, 1),
        verify_witness('opposite_chain', [(0, 1), (2, 3), (0, 2), (1, 0), (3, 1)], 2, 3),
    ]
    for witness in witnesses:
        unoriented = [sorted(edge) for edge in witness['oriented_cnot_schedule']]
        assert unoriented == examples[witness['name']]
    result = {
        'remaining_schedule_shapes': dict(sorted(count.items())),
        'unoriented_example_schedules': examples,
        'partial_axis_witnesses': witnesses,
        'new_exclusions': 0,
    }
    with open('global_lower_bound6_gate3_selected_witness_result.json', 'w', encoding='utf-8') as out:
        json.dump(result, out, indent=2, sort_keys=True)
        out.write('\n')
    print(json.dumps({'remaining': sum(count.values()), 'shapes': count,
                      'new_exclusions': 0}, sort_keys=True))


if __name__ == '__main__':
    main()
