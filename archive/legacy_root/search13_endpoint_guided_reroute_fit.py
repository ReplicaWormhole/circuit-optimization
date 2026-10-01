"""Joint free-local fit for endpoint-guided prefixes and cross-pair tail reroutes.

Selected routes replace one 13 or 01 CNOT near the first/final-pair
boundary with directed 12 or 03.  Each fit starts from an endpoint-guided
13-CNOT candidate with the original tail and perturbs all local SU(2)
angles.  The objective is label-free cycle diagonalization.
"""

import argparse
import json

import numpy as np

from check_circuit import evaluate
from delete14_search import u3_to_axis
from search13_joint_gauge_five_endpoint import ROOT, candidate, fit
from search13_endpoint_guided_reroute_rank import TAIL


CASES = (
    (1, "warm", 3, (1, 2)),  # misses 7,11; reroute final first-pair CX
    (1, "warm", 4, (0, 3)),  # misses 7,11; reroute initial final-pair CX
    (2, "phase", 3, (1, 2)), # misses 6,7
    (0, "warm", 4, (0, 3)),  # misses 11,14
)
SIGMAS = (0.02, 0.15)


def decode_layers(gates):
    chunks = [[]]
    for gate in gates:
        if gate["gate"] == "cx":
            chunks.append([])
        else:
            chunks[-1].append(gate)
    assert len(chunks) == 14
    local = []
    for chunk in chunks:
        layer = chunk[-4:]
        assert [g["gate"] for g in layer] == ["u3"]*4
        assert [g["qubit"] for g in layer] == list(range(4))
        local.append([u3_to_axis(g) for g in layer])
    return np.asarray(local)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=25060)
    parser.add_argument("--maxiter", type=int, default=300)
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    rows = []
    for case, (source_case, source_start, position, replacement) in enumerate(CASES):
        source_path = ROOT / f"search13_endpoint_guided_case{source_case}_{source_start}.json"
        source = json.loads(source_path.read_text())
        old_edges = [(g["control"], g["target"]) for g in source["gates"]
                     if g["gate"] == "cx"]
        assert len(old_edges) == 13
        assert tuple(old_edges[5:]) == TAIL
        edges = old_edges.copy()
        assert edges[5+position] != replacement
        edges[5+position] = replacement
        warm = decode_layers(source["gates"])
        for start, sigma in enumerate(SIGMAS):
            initial = warm+rng.normal(0, sigma, warm.shape)
            angles, record = fit(edges, initial, args.maxiter)
            circuit = candidate(edges, angles)
            path = ROOT / f"search13_endpoint_guided_reroute_case{case}_s{start}.json"
            path.write_text(json.dumps(circuit, indent=2)+"\n")
            checked = evaluate(circuit)
            record.update(case=case, start=start, sigma=sigma,
                          source=source_path.name, prefix_edges=edges[:5],
                          rerouted_tail_position=position,
                          replacement_edge=replacement,
                          tail_edges=edges[5:], candidate=path.name,
                          off_diagonal_error=checked["off_diagonal_error"],
                          valid_diagonalizer=checked["valid_diagonalizer"])
            rows.append(record)
            print(json.dumps(record, sort_keys=True), flush=True)
    result = {"cases": CASES, "sigmas": SIGMAS, "seed": args.seed,
              "maxiter_per_fit": args.maxiter, "rows": rows,
              "best_loss": min(rows, key=lambda r: r["loss"]),
              "best_off_diagonal": min(rows, key=lambda r: r["off_diagonal_error"])}
    (ROOT / "search13_endpoint_guided_reroute_fit_result.json").write_text(
        json.dumps(result, indent=2)+"\n")


if __name__ == "__main__":
    main()
