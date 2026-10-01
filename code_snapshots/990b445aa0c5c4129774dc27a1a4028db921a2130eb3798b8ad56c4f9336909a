"""Exact rational nullspaces and square-coordinate checks for 24D Grams."""

import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
SHAPES = ('adjacent_adjacent', 'adjacent_opposite', 'opposite_adjacent')


def main():
    result = {}
    for shape in SHAPES:
        data = json.loads((ROOT / f'global_lower_bound7_degree3_centered_gram_{shape}_result.json').read_text())
        gram = sp.Matrix(data['gram'])
        null = gram.nullspace()
        monomials = data['monomials']
        forced = [data['basis_labels'][i] for i in range(24)
                  if all(v[monomials.index([i,i])] == 0 for v in null)]
        rank = gram.cols - len(null)
        assert rank == data['rank_mod_prime_1000003']
        result[shape] = {'exact_rank':rank, 'exact_nullity':len(null),
                         'forced_zero_coefficient_squares':forced,
                         'ordered_triple_j_k_l':data['ordered_triple_j_k_l']}
    (ROOT / 'global_lower_bound7_degree3_centered_nullspace_result.json').write_text(
        json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps(result, sort_keys=True))


if __name__ == '__main__':
    main()
