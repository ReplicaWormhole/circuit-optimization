"""Exact trace certificate and exhaustive five-CNOT pair-schedule closure.

The mathematical lemmas are in global_lower_bound6_complete_filter.md.
The adjacent-chain rank-one certificates are separately recomputed by
global_lower_bound6_chain_certificate.py and
global_lower_bound6_adjacent_v_certificate.py; the opposite-chain symmetry
certificate is global_lower_bound6_opposite_chain_certificate.py.
"""

from collections import Counter
from itertools import product
import json

import sympy as sp


EDGES = tuple((a, b) for a in range(4) for b in range(a + 1, 4))
ADJACENT = {(0, 1), (0, 3), (1, 2), (2, 3)}
OPPOSITE = {(0, 2), (1, 3)}


def exact_trace_tensors():
    one = {
        'I': sp.eye(2),
        'X': sp.Matrix([[0, 1], [1, 0]]),
        'Y': sp.Matrix([[0, -sp.I], [sp.I, 0]]),
        'Z': sp.diag(1, -1),
    }
    shift = sp.zeros(16)
    for state in range(16):
        shift[((state & 1) << 3) | (state >> 1), state] = 1
    squared = shift * shift
    output = {}
    for pair in sorted(OPPOSITE):
        matrix = []
        for a in 'XYZ':
            row = []
            for b in 'XYZ':
                factors = [one[a] if j == pair[0] else
                           one[b] if j == pair[1] else one['I'] for j in range(4)]
                operator = sp.kronecker_product(*factors)
                row.append(int(sp.trace(squared * operator)))
            matrix.append(row)
        assert matrix == [[4, 0, 0], [0, 4, 0], [0, 0, 4]]
        output[''.join(map(str, pair))] = matrix
    assert sp.trace(squared) == 4
    return output


def degree(schedule):
    out = [0] * 4
    for a, b in schedule:
        out[a] += 1
        out[b] += 1
    return out


def full_forward_support(schedule):
    for input_wire in range(4):
        support = {input_wire}
        for edge in schedule:
            if support.intersection(edge):
                support.update(edge)
        if len(support) != 4:
            return False
    return True


def terminal_degree_two_edges(schedule, deg):
    out = []
    for j in range(4):
        if deg[j] != 2:
            continue
        second = [t for t, edge in enumerate(schedule) if j in edge][1]
        edge = schedule[second]
        if all(not set(edge).intersection(later) for later in schedule[second + 1:]):
            out.append((j, edge))
    return out


def penultimate_chain_type(schedule, deg):
    found = []
    for j in range(4):
        if deg[j] != 2 or j not in schedule[3] or j in schedule[4]:
            continue
        k = next(wire for wire in schedule[3] if wire != j)
        if k not in schedule[4]:
            continue
        ell = next(wire for wire in schedule[4] if wire != k)
        first = tuple(sorted((j, k)))
        endpoint = tuple(sorted((j, ell)))
        if first in OPPOSITE:
            kind = 'opposite_first_trace'
        elif endpoint in OPPOSITE:
            kind = 'adjacent_chain_opposite_endpoint'
        else:
            kind = 'adjacent_v_both_endpoints_adjacent'
        found.append((kind, (j, k, ell)))
    return found


def main():
    traces = exact_trace_tensors()
    counts = Counter()
    examples = {}
    for schedule in product(EDGES, repeat=5):
        deg = degree(schedule)
        if min(deg) < 2:
            counts['degree_rejected'] += 1
            continue
        if not full_forward_support(schedule):
            counts['forward_support_rejected'] += 1
            continue
        terminal = terminal_degree_two_edges(schedule, deg)
        if any(edge in ADJACENT for _, edge in terminal):
            counts['terminal_adjacent_gram'] += 1
            continue
        if terminal:
            assert all(edge in OPPOSITE for _, edge in terminal)
            counts['terminal_opposite_trace'] += 1
            continue
        chains = penultimate_chain_type(schedule, deg)
        assert chains, f'unclassified schedule {schedule}'
        kinds = {kind for kind, _ in chains}
        assert len(kinds) == 1, f'conflicting chain types {schedule}: {chains}'
        kind = next(iter(kinds))
        counts[kind] += 1
        examples.setdefault(kind, {'edges': [list(edge) for edge in schedule],
                                   'selected_chain': list(chains[0][1])})
    assert sum(counts.values()) == 6 ** 5
    assert counts == Counter({
        'degree_rejected': 5196,
        'forward_support_rejected': 1572,
        'terminal_adjacent_gram': 576,
        'terminal_opposite_trace': 288,
        'adjacent_chain_opposite_endpoint': 48,
        'adjacent_v_both_endpoints_adjacent': 48,
        'opposite_first_trace': 48,
    })
    result = {
        'target': 'four-qubit shift with arbitrary local gates, all CNOT pairs and directions',
        'exact_half_shift_opposite_pair_pauli_trace_tensors': traces,
        'ordered_unoriented_pair_schedules': 6 ** 5,
        'exclusive_rejection_counts': dict(sorted(counts.items())),
        'unclassified_schedule_count': 0,
        'examples_for_penultimate_classes': examples,
    }
    with open('global_lower_bound6_complete_filter_result.json', 'w', encoding='utf-8') as out:
        json.dump(result, out, indent=2, sort_keys=True)
        out.write('\n')
    print(json.dumps(result['exclusive_rejection_counts'], sort_keys=True))


if __name__ == '__main__':
    main()
