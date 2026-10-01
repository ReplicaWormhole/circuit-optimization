"""Exact 127 mixed-cut ranks of the ungauged exact14 eight-CNOT tail.

Compute each rank directly over Q(zeta_96); reductions modulo 97 and 193
are independent lower-bound cross-checks. Output is an exact target profile
for subsequent necessary seven-CNOT topology screens.
"""

from collections import Counter
import json
from pathlib import Path

import numpy as np
from sympy import primitive_root

from search13_crosspair_single_reroute_rank import exact_rank
from search13_fulltail_invariant_exact_ansatz import tail_matrix
from search13_rank8_neighborhood_mixedcut import field_mod, rank_mod, realignment
from search13_sparse_gauge_clifford_pairs_certificate import mixed_realign
from search13_secondpair_mask_transport import CyclotomicPair


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "search13_identity_tail_allcuts_result.json"
MASKS = tuple(range(1, 255, 2))
PRIMES = (97, 193)


def main():
    assert not OUT.exists()
    exact = CyclotomicPair()
    tail = tail_matrix(exact)
    exact_ranks = {}
    for index, mask in enumerate(MASKS):
        exact_ranks[mask] = exact_rank(mixed_realign(tail, mask, exact), exact)
        if index % 16 == 15:
            print(json.dumps({"completed_masks": index + 1}), flush=True)
    assert len(exact_ranks) == 127 and exact_ranks[83] == exact_ranks[163] == 4
    modular = {}
    for prime in PRIMES:
        root = pow(primitive_root(prime), (prime - 1) // 96, prime)
        matrix = np.asarray([[field_mod(v, root, prime) for v in row]
                             for row in tail], dtype=np.int64)
        modular[prime] = {mask: rank_mod(realignment(matrix, mask), prime)
                          for mask in MASKS}
    assert all(modular[p][m] <= exact_ranks[m] for p in PRIMES for m in MASKS)
    discrepancies = {str(p): [m for m in MASKS if modular[p][m] != exact_ranks[m]]
                     for p in PRIMES}
    result = {"source": "topology14_exact_matchgate_rational.json",
              "output_gauge": "identity", "field": "Q(zeta_96)",
              "mask_convention": "8 bits: output qubits 0..3 then input qubits 0..3; odd masks represent complementary cut once",
              "exact_mixed_cut_ranks": {str(k): v for k, v in exact_ranks.items()},
              "rank_histogram": dict(Counter(exact_ranks.values())),
              "modular_rank_discrepancies": discrepancies,
              "scope": "exact 127-cut rank profile of the original eight-CNOT tail; no seven-CNOT synthesis"}
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"histogram": result["rank_histogram"],
                      "modular_discrepancies": {k: len(v) for k, v in discrepancies.items()}}))


if __name__ == "__main__":
    main()
