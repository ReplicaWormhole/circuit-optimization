# One bounded full-local single-deletion refit from exact11

Board140, owner numerical10. Source: accepted run411 exact11 candidate,
independently certified by run412; SHA256
557c55cddefb697bcd4e038db300117d0aaaf6f6af97e851e0831baf257dbe4a.
This is distinct from run409's exact12 source deletion. Preparation performs
no matrices, RNG, AD or fitting. Root reserves, snapshots, links and executes
only after independently reserved preflight passes. Preserve all preparation
unexecuted if algebraic work obtains exact10 before fit dispatch.

Delete each of the eleven source CX in chronological ordinal order0..10 once.
Retain every other CX's order and direction; no other gate is deleted or
reordered. Source CX positions/edges are 0/02,2/13,4/01,5/23,7/02,10/32,
13/02,16/32,21/10,25/12,28/32. Each resulting topology has exactly10CX with
11 full four-wire local SU2 layers before/between/after the retained CXs.
All132 real rotation-vector coordinates are free, layer-major then wire0..3
then axesx,y,z. Qubit0 is MSB; gates are chronological; target is UV4=DU.

For each deletion, merge source one-qubit gates per wire until the next
retained CX, multiplying chronologically on the left. Convert merged U2 gates
to principal SU2 rotation vectors by dropping determinant phase and choosing
nonnegative trace; these phases are global in the complete circuit. Retain
the literal source-deleted gate list, base132 vector, base local matrices,
base complete serialized gate list, noise132 and initial132 vectors and full
initial gate list before fitting. Every serialized word has54 elementary
gates, namely44U3 and10CX; retain best132 vectors, complete best gate lists,
objective/gradient-norm history, optimizer counts, stop reasons and checker
metrics for all11 cases, including skipped or failed cases.

One NumPy Generator(PCG64(10012201)) allocates eleven consecutive
normal(0,0.025,size=132) draws at startup in deletion ordinal order0..10.
All draws are saved before any fit. Each fit starts once from its base plus
assigned perturbation. No restart, topology variant or deepening is allowed.

SciPy L-BFGS-B uses Torch complex128/float64 analytic gradients of
sum(abs(offdiag(U V4 Udag))^2)/16. Maximum80iterations,4000objective/gradient
calls and20line-search steps per case; ftol1e-15,gtol1e-10. A loss1e-22 trigger
stops that case only, without claiming exactness. OneCPUthread; at most
8seconds per case, optimizer deadline at96seconds from main start,
internal105seconds and external115seconds. Each case uses the minimum of
its8second deadline and the global96second deadline. Guard/import setup and
initialization consume the same internal/global budget. Internal deadline
preserves current best vectors/lists/history and records failure; external
kill is distinguishable from a completed result. Exclusive execution marker
rejects a silent rerun. No new numerical or symbolic search follows failure.

Root must reserve target cycle, method exactly as config, agent numerical10,
parent411, exact canonical config bytes, and fit.py's frozen code hash. Then
create experiments/runs/ID with experiment_log.py workspace ID --code PATH,
copy family.py,source_candidate.json,topologies.json,config.json,PLAN.md and
dependencies.json alongside fit.py, and link board140 to ID. Only root runs
`timeout 115s python3 experiments/runs/ID/fit.py` after preflight passes.
Guard checks canonical rawconfig equality with ledger, fixed values/types,
source/family/topology/PLAN/dependency/checker/acceptance hashes, runtime
versions and binary dependency hashes, retained accepted run411 source,
reserved code snapshot, active ledger status/agent/parent/target/method and
board status/ownership/run link before matrices or RNG. Internal import of
family.py constructs no circuit matrix/tensor/random draw.

Independent preflight should check every collapse/topology, base versus
literal deletion up to global phase, complete base/initial serializer,
NumPy/Torch matrix/loss agreement and directional gradients from one frozen
separate audit direction per case. Fit output must undergo complete endpoint
audit before any exactness investigation or incumbent decision.

Any failure is bounded numerical evidence, never a topology exclusion,
local-minimum proof, exact certificate, incumbent promotion or ten-CX lower
bound. The established interval remains6<=C_min<=11 absent independent exact
certification of a new lower-cost gate list. Source and acceptance artifacts
remain untouched. No shared checker modification or dependency installation.
