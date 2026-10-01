# Joint prefix mutation with continuous output-gauge continuation

Ledger run 325 tested a 13-CNOT whole-circuit search derived from the exact
14-CNOT diagonalizer. In each proposed topology, one prefix CNOT is deleted,
a different retained prefix CNOT is rerouted, and both CNOTs of one suffix
matchgate block are rerouted. The two explicit topologies are in
`search13_joint_prefix_orbit_fit.py`. Every one-qubit layer throughout the
resulting circuit is fitted jointly: 14 layers, four wires, three axis-angle
parameters per wire, or 168 free real parameters.

The continuation first fits the full circuit to the orbit of the exact
diagonalizer under the complete commuting output gauge with sector sizes
`(6,4,3,3)`. Each of two outer stages computes the best gauge by a numerical
block Procrustes SVD, then fits all local layers by process overlap. A final
stage minimizes the label-free cycle off-diagonal Frobenius loss. This differs
from run 249's endpoint-selected five-CNOT prefixes with unchanged tail and
run 257's single boundary tail reroutes. It also differs from the adaptive
topology search: proposals are simultaneous changes to the exact circuit's
prefix and an entire suffix block, and the first objective is process
continuation over the continuous output eigenspace gauge.

Both zero and seeded random starts were tried for each topology, with two
70-iteration process stages and one 180-iteration direct stage. All eight
saved process and final circuits were checked with `check_circuit.py`; each
has 13 CNOTs and fails to diagonalize the shift. The best final candidate is
`search13_joint_prefix_orbit_case1_random_direct.json`, with maximum cycle
off-diagonal error `0.3866559358256528` and direct Frobenius loss
`0.07873301259438109`. Its optimizer stopped at the iteration limit, so
further local refinement is possible. These four bounded fits have no no-go
implication for other starts, topologies, or 13-CNOT circuits.

Ledger run 327 continued that best candidate for up to 800 L-BFGS-B
iterations from the saved parameters and perturbations of sizes 0.02 and
0.08. The candidate's saved parameters reproduce the run-325 loss to
better than $10^{-8}$ before refinement. All three continuations ended
at loss about `0.078733012516` and failed the independent checker,
with maximum off-diagonal errors about `0.386658`. This supports the
interpretation of a local basin; it does not prove that the topology
has no exact solution.

Reproduce:

```bash
python3 search13_joint_prefix_orbit_fit.py --seed 32700 --outer 2 --process-maxiter 70 --direct-maxiter 180
python3 check_circuit.py search13_joint_prefix_orbit_case1_random_direct.json
python3 experiment_log.py show 325
python3 experiment_log.py audit
```
