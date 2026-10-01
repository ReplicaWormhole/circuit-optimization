# One bounded full-local single-deletion refit from exact12

Board132, numerical_sol. Source: run406 exact12 candidate, independently
certified by run407; SHA25658a862a4acef4e768c4613e9b142504592f5237a68539941148653a7f083eea2.
This tests deletions of the new12CX source, distinct from historical searches
deleting gates from an older13CX circuit. No matrix or optimizer computation
is authorized until coordinator reservation, snapshot and dispatch.

Delete each of the twelve source CX in chronological ordinal order0..11,
once. Each resulting word has exactly11CX and12 free local layers, four
qubits each. All144 real SU2 rotation-vector coordinates are free; no source
local gate or wire is frozen. SU2 restriction discards only global phases.
Source one-qubit gates are merged per wire in each layer after that specific
CX deletion, converted to principal rotation vectors, and retained as base
vectors, base matrices and complete base/deleted-source gate lists.

One NumPy Generator(PCG64(10012101)) supplies twelve independent chronological
144-coordinate normal draws, sigma0.025 radians, to perturb all coordinates.
Each optimizer starts from its perturbed vector; base, noise and initial
vectors and complete gate lists are retained. There are no restarts.

Method SciPy L-BFGS-B with analytic Torch autograd gradients, max80iterations
per deletion, max4000objective/gradient calls per deletion, max20line-search
steps, ftol1e-15, gtol1e-10. One thread,110seconds internal and120seconds
external; optimization has at most8seconds per case, leaving validation and
output time inside the total bound. No case is deepened after its first fit.
Every initial/best serialized word has59 gates (48U3 plus11CX), actual counts,
shared numerical checker metrics, loss, stop reason and vector history.
All cases remain in the output even if budget prevents a fit. No exact
certificate is run by this diagnostic.

Root must reserve the exact config object using its exact method string,
create run workspace using experiment_log.py workspace, copy family.py,
source_candidate.json,topologies.json,config.json andPLAN.md alongside the
reserved fit.py, then link board132. The script guards all frozen keys,
candidate/topology/family/shared checker/runtime/dependency hashes and
matching active ledger/board ownership before matrices or random draws.
Exclusive execution marker prevents silent rerun. Reserved independent
preflight will check collapse, serialization, Torch/NumPy matrix agreement and
gradients before coordinator dispatch of the fit.

Command after authorization:
`timeout 120s python3 experiments/runs/ID/fit.py`.

Algebraic_sol's static router audit finds CX01 commutes with the merged block
and can be moved adjacent to CH. Literal router deletion fails, but it does
not exclude local-basis refitting. This run does not silently reorder the
source word, and any failure is bounded numerical evidence, not an11CX lower
bound, topology exclusion or stationary-point proof.
