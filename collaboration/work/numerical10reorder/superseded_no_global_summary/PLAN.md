# Round23 two-start reordered ten-CX diagnostic preparation

Board148 owner and ledger agent numerical10reorder, targetV4,parent416.
Read-only source is immutable416/cases/07_delete_cx07/result.json best132;
raw SHA853d53241573363a5baf9d86f833cdf9ec5c6304cc0b2b2fa168466015864ef7.
Source oldbest gate list SHA26b9f34ab733deb77127dcd666929ae8c8a515746892f5ddced57fdf4d36985b.
Source416 was independently endpoint-validated by417 and closed failed.
This new schedule changes noncommuting CX chronology, so it is a distinct
finite diagnostic rather than repeating416's literal deletion topology.

Old CX sequence02,13,01,23,02,32,02,10,12,32; new CX sequence
02,13,01,23,02,32,10,02,12,32. Swap old CX ordinals6 and7 only; retain all
other edges/directions and all132 local coordinates at their original layer
positions. This is NOT a commutation or equivalent-circuit claim. New topology
is exactly[[0,2],[1,3],[0,1],[2,3],[0,2],[3,2],[1,0],[0,2],[1,2],[3,2]].
Qubit0 is MSB; gates chronological; diagonalization targetUV4=DU. Eleven
full4wireSU2 layers before/between/after the ten fixed CXs,132 coordinates
(layer then wire0..3 then axesx,y,z). All coordinates remain free; no control
sector, CH wrapper, final output basis, target axis or phase anchor is frozen.
No sourcecollapse is used: base132 is transplanted verbatim from saved best132.
Literal newbase preserves every saved old U3 local angle and replaces only
CX chronology. Retain old/new topology, raw source JSON and oldbest full word,
base132 JSON and literal54gate newbase word as preparation evidence.

One NumPy Generator(PCG64(10012301)) supplies exactly two consecutive
normal(0,sigma,size=132) draws, firstsigma0.025 then0.25, allocated and saved
upfront before any fit. Case0 and case1 both start from unchangedbase+assigned
noise. Two starts only; no other RNG, restart, topology variant, deepening or
ranking between starts. Per-case optimizer best retention is allowed; no
winning start or globalbest is selected or written.

OneCPUthread; SciPy L-BFGS-B with Torch analyticfloat64/complex128 gradients
of sum(abs(offdiag(UV4Udag))^2)/16; max250iterations,max10000calls,maxls20,
ftol1e-15,gtol1e-10. At most15seconds percase and globaloptimizerdeadline35s
from main start; internal45s/external60s. Setup/initialization consumes the same
internal/global budget. Loss1e-22 stops only thatcase and does not certify
exactness. Internal deadline salvages current best vectors/lists/histories;
external kill yields an incomplete failure requiring truthful closeout.

Save each case's54gate/10CX base, literalbase, initial and best complete words;
base/noise/initial/best132 vectors; base/initial/best numerical checker metrics
and losses; analytic gradient at retained best, objective/gradient-norm history,
iterations/calls/stopreason/timing. Both case records remain present if budget
skips a fit. No new computation follows failure. Family.py is an unchanged copy
of numerical10delete/family.py; collapse helper is unused in this diagnostic.

Root reserves exactcanonicalraw config bytes (no trailingnewline), exactmethod,
agentnumerical10reorder,targetV4,parent416 and frozen fit.py. Create workspace
using experiment_log.py workspace ID --code PATH. Copy fit.py,family.py,
source_case_result.json,source_best_candidate.json,base_vector.json,
base_candidate.json,topologies.json,input_manifest.json,config.json,PLAN.md,
dependencies.json alongside reserved script; link board148 with actualowner.
Root executes `timeout 60s python3 experiments/runs/ID/fit.py` only after
by-hand algebraic review and separately reserved independent preflight pass.

Guard precedes matrices/RNG/AD: rawcanonicalconfig equality with activeledger;
fixedkeys/value/types; sameownerboardlink; ledger agent/parent/target/method;
code snapshot; script/family/checker/source/topology/base/PLAN/inputmanifest/
dependency snapshots; runtime versions/binary dependencies; immutable origin
file hashes and source416closedfailed provenance; exactunchangedsourcebase
coordinates; complete literalnewbase gate list; exclusiveexecutionmarker.
Manifest origins include416case7 result/word/globalresult/config/topologies,
417independentresult and accepted411/acceptance. No sharedfile edits or installs.

By-hand review has identified only a restricted classical-sector/computational-x
obstruction. Full local bases after earlyu-to-x interaction and sector mixing
can escape it, so this family must not be narrowed on that basis. Formal review
and dispatch decisions remain root-owned. Finite numerical results establish
no exactcertificate, topology exclusion, localminimum or ten-CXlowerbound;
accepted11 incumbent remains unchanged absent independent exactcertification.
