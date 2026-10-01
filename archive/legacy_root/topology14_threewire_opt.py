"""Test a 4-CNOT tail with a two-CNOT three-wire boundary factorization.

After the exact ten-CNOT prefix, retain the two-CNOT F2 interaction on
wires (0,1), but use only cross edges (2,1) and (1,3) to diagonalize the
three-wire boundary on (1,2,3). All local SU2 gates are optimized and the
eigenbasis is unconstrained.
"""

import argparse
import json
from pathlib import Path

import torch

from topology14_fourcx_opt import prefix_data, loss, candidate


torch.set_num_threads(1)
EDGES = ((2, 1), (0, 1), (0, 1), (1, 3))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--steps", type=int, required=True)
    parser.add_argument("--lr", type=float, default=0.04)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    torch.manual_seed(args.seed)
    prefix_gates, prefix = prefix_data()
    angles = torch.nn.Parameter(0.3 * torch.randn(5, 4, 3, dtype=torch.float64))
    optimizer = torch.optim.Adam([angles], lr=args.lr)
    best = float("inf")
    best_angles = None
    for step in range(args.steps):
        optimizer.zero_grad()
        objective = loss(angles, EDGES, prefix)
        value = float(objective.detach())
        if value < best:
            best, best_angles = value, angles.detach().clone().numpy()
        objective.backward()
        optimizer.step()
        if step % 200 == 0:
            print(json.dumps({"step": step, "loss": value, "best": best}),
                  flush=True)
    args.output.write_text(json.dumps(candidate(prefix_gates, best_angles, EDGES),
                                      indent=2) + "\n")
    print(json.dumps({"topology": EDGES, "best_loss": best,
                      "output": str(args.output)}))


if __name__ == "__main__":
    main()
