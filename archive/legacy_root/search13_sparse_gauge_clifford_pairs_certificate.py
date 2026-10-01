"""Exact rank certificate for the minimum of the sparse Clifford-pair family."""

import json
from pathlib import Path

from search13_crosspair_single_reroute_rank import exact_rank
from search13_fulltail_invariant_exact_ansatz import LABELS, tail_matrix
from search13_fulltail_invariant_exact_witness import gauge_matrix
from search13_secondpair_mask_transport import CyclotomicPair


ROOT = Path(__file__).resolve().parent
SCREEN = ROOT / "search13_sparse_gauge_clifford_pairs_result.json"
OUT = ROOT / "search13_sparse_gauge_clifford_pairs_certificate_result.json"
MASKS = (83, 163)


def mixed_realign(matrix, mask, exact):
    left_axes = tuple(j for j in range(8) if mask & (1 << j))
    right_axes = tuple(j for j in range(8) if not mask & (1 << j))
    nleft, nright = 1 << len(left_axes), 1 << len(right_axes)
    out = [[exact.zero for _ in range(nright)] for _ in range(nleft)]
    for y in range(16):
        for x in range(16):
            bits = tuple((y >> (3 - j)) & 1 for j in range(4)) + \
                   tuple((x >> (3 - j)) & 1 for j in range(4))
            l = sum(bits[j] << (len(left_axes) - 1 - k)
                    for k, j in enumerate(left_axes))
            r = sum(bits[j] << (len(right_axes) - 1 - k)
                    for k, j in enumerate(right_axes))
            out[l][r] = matrix[y][x]
    return out


def main():
    assert not OUT.exists()
    scan = json.loads(SCREEN.read_text())
    for mask in MASKS:
        assert scan["minimum_mod97_ranks_by_mask"][str(mask)]["rank"] == 12
        example = scan["minimum_mod97_ranks_by_mask"][str(mask)]["example"]
        assert example["row_pairs"] == [[1, 14]] and example["phase_indices"] == [[0, 0]]
    exact = CyclotomicPair()
    g, _ = gauge_matrix(exact)
    h = (exact.z**12 + exact.z**-12) * exact.half
    assert 2*h*h == exact.one and LABELS[1] == LABELS[14]
    old_a, old_b = g[1][:], g[14][:]
    g[1] = [h*(a+b) for a, b in zip(old_a, old_b)]
    g[14] = [h*(a-b) for a, b in zip(old_a, old_b)]
    target = exact.matmul(g, tail_matrix(exact))
    ranks = {str(mask): exact_rank(mixed_realign(target, mask, exact), exact)
             for mask in MASKS}
    assert ranks == {"83": 12, "163": 12}
    OUT.write_text(json.dumps({
        "field": "Q(zeta_96)", "family_screen": SCREEN.name,
        "attaining_gauge": "left-compose exact G by real Hadamard on same-eigenvalue output rows 1 and 14",
        "hadamard_entry": "(zeta_96^12 + zeta_96^-12)/2 = 1/sqrt(2)",
        "unitarity_reason": "the pair block is [[h,h],[h,-h]] with 2h^2=1; all other rows unchanged",
        "eigenvalue_preservation_reason": "output labels of rows 1 and 14 agree",
        "exact_mixed_cut_ranks": ranks,
        "minimum_over_enumerated_family": "both ranks >=12 by exhaustive mod97 screen, and this exact gauge attains 12 at both masks",
        "scope": "exact minimum for these two masks and this finite sparse family only; no global output-gauge obstruction",
    }, indent=2, sort_keys=True) + "\n")
    print(json.dumps(ranks))


if __name__ == "__main__":
    main()
