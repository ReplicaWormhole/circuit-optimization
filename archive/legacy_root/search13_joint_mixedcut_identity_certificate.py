"""Exact mixed-cut ranks of the ungauged exact14 eight-CNOT tail."""

import json
from pathlib import Path

from exact_check import exact_eigenvalue_labels
from search13_crosspair_single_reroute_rank import exact_rank
from search13_fulltail_invariant_exact_ansatz import tail_matrix
from search13_sparse_gauge_clifford_pairs_certificate import mixed_realign
from search13_secondpair_mask_transport import CyclotomicPair


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "topology14_exact_matchgate_rational.json"
OUT = ROOT / "search13_joint_mixedcut_identity_certificate_result.json"
MASKS = (83, 163)


def main():
    assert not OUT.exists()
    source = json.loads(SOURCE.read_text())
    assert exact_eigenvalue_labels(source)["exact_diagonalizer"]
    exact = CyclotomicPair()
    tail = tail_matrix(exact)
    ranks = {str(mask): exact_rank(mixed_realign(tail, mask, exact), exact)
             for mask in MASKS}
    assert ranks == {"83": 4, "163": 4}
    OUT.write_text(json.dumps({
        "source": SOURCE.name, "output_gauge": "identity",
        "field": "Q(zeta_96)", "prefix_cnot_count": 6,
        "tail_cnot_count": 8, "exact_V4_diagonalizer": True,
        "exact_mixed_cut_ranks": ranks,
        "scope": "exact ranks at two mixed cuts for the original eight-CNOT tail; seven-CNOT synthesis remains open",
    }, indent=2, sort_keys=True) + "\n")
    print(json.dumps(ranks))


if __name__ == "__main__":
    main()
