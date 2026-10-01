"""Optimize an eight-CNOT tail with one butterfly CNOT moved to a cross edge.

The six-CNOT gauged parity prefix is fixed. The Fourier tail starts from
the exact 15-CNOT circuit without its CZ, but its final CNOT (2,3) is
replaced by (2,1). All interleaved one-qubit SU2 layers are then free.
"""

import argparse
import json
from pathlib import Path

import torch

from topology14_eightcx_opt import source_data, layers_from_suffix, emit
from topology14_fourcx_opt import loss


torch.set_num_threads(1)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--steps", type=int, required=True)
    parser.add_argument("--lr", type=float, default=0.02)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    torch.manual_seed(args.seed)
    prefix, prefix_u, suffix = source_data()
    initial, edges = layers_from_suffix(suffix)
    assert edges[-1] == (2, 3) and len(edges) == 8
    edges[-1] = (2, 1)
    parameters = torch.nn.Parameter(torch.tensor(initial, dtype=torch.float64) +
                                    0.03 * torch.randn(initial.shape, dtype=torch.float64))
    optimizer = torch.optim.Adam([parameters], lr=args.lr)
    initial_loss = float(loss(parameters, edges, prefix_u).detach())
    best = float("inf")
    best_parameters = None
    for step in range(args.steps):
        optimizer.zero_grad()
        objective = loss(parameters, edges, prefix_u)
        value = float(objective.detach())
        if value < best:
            best = value
            best_parameters = parameters.detach().clone().numpy()
        objective.backward()
        optimizer.step()
        if step % 200 == 0:
            print(json.dumps({"step": step, "loss": value, "best": best}),
                  flush=True)
    args.output.write_text(json.dumps(emit(prefix, best_parameters, edges),
                                      indent=2) + "\n")
    print(json.dumps({"initial_loss": initial_loss, "best_loss": best,
                      "topology": edges, "output": str(args.output)}))


if __name__ == "__main__":
    main()
