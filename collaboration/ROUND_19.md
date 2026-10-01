# Round19: reordered full98 Bell-decoder attempt

2026-10-01. The exact12-CNOT success goal remains unmet and active. The
accepted13 incumbent is unchanged. This is a different chronology from the
old full98 family, not an exclusion of either family.

## Hypothesis and frozen fit

Move F02/F13 Bell decoders before the fan-in:
`L0 F02 L1 F13 L2 CX01 A1 CX21 B1 CX31 L3 CX12 L4 F01 L5 F23 L6`.
The algebraic reviewer checked the initializer signs/rotvec by hand, including
control vector2*pi/(3*sqrt3)*(-1,+1,-1), target(-pi/2,0,0), postX control(pi,0,0),
and first two F(-pi/4,0). Both pair input locals are in L0; onlyX0 in L1,
onlyX1 in L2. The full circuit still has its fixed CNOT suffix.

Fit402 froze seed9261801, independent default_rng Gaussian noise sigma0.1 on
all98 coordinates, one start, zero restarts. All98 coordinates were free.
Objective ||offdiag(UV4U-dagger)||_F^2/16, SciPy L-BFGS-B with Torch AD gradient;
max1000iterations/20000optimizer calls,maxls30,ftol1e-15,gtol1e-10,
110seconds internal/120seconds external,one thread, plus at mostone explicitly
counted best-loss recomputation. Code SHA
3f4777c2492e34e47fb5655df4c288f432f8dbad0e74230afeff872c257448ac;
family SHA9936da7408d5202a00e784daa4b33d950efa68229bec06c006457862411303c8;
config SHA1fb2864e4f9750324080d95f8193d46dce75c3939373036b4187b8e080c4616a.
Reservation/workspace/board106link preceded computation.

## Preflight and guarded failures

Fit401 was reserved but unexecuted: static source/config review found missing
top-level noise_sigma. Preserve inputs/snapshot; ledgerfailed/board105inconclusive.
Corrected402 added that field under a fresh reservation. No scientific evaluation
occurred in401.

Independent preflight403 executed once but stopped at exact module-version guard:
Torch distribution metadata2.11.0 differed from module2.11.0+cu130. Zero
coordinate/RNG/matrix evaluations; ledgerfailed/board103inconclusive.
Refrozen404 bound actual module version under board107 and a new ledger row,
retaining old inputs. It completed once exit0 in1.729seconds:3fixedpoints,
6Torch/NumPyfull-family builds,6native/compiledlist builds,2prefix/direct
products,14top-level matrix builds,external30seconds/one thread/noAD/targetloss.
All checks passed: maximumerror1.4433e-15, decoderprefix directcomparison
1.86794e-16, native4CX4F and compiled12CX. Fullbase matrix was not mistaken
for the decoder prefix. This is finite preflight evidence, not a certificate.

## Numerical result

After404 passed, the single402command exited0. Runtime1.364887827seconds,
113iterations127objective-gradient calls plus1finalrecomputation. Bestloss
0.3750000000000021, gradientnorm3.92257e-8. Best compiled12-CNOT list has
maxoff0.9999999999999969 and is invalid. The optimizer's convergence flag is
not diagonalization. Ledger402failed/board106inconclusive; all inputs, complete
initial/best lists, bestvector, history, versions, outputs and conclusion kept.
Coordinator independently checked saved hashes, best-selection and budget
arithmetic without further matrices. Independent endpoint audit108 passed (POSTFIT402_AUDIT.json): all four full
lists match coordinates/compilation upglobalphase below8e-16, unitarity below
1.2e-15. Independent initial/best losses0.8181143009791337/0.375000000000002
and maxoff0.9018603550676626/0.9999999999999974 confirm invalidity. Loss
reduction did not reduce the largest residual. Initialvector matches404hash;
history/best-selection/iteration and evaluation limits independently checked.
This finite validation is not an exact certificate.

## Interpretation and next hypothesis

One seeded small perturbation of the Bell-derived prefix did not find a valid
member. It neither proves a local minimum nor excludes this family. The algebraic
change of basis and prefix encodability are not a factorization of the required
suffix. Best next hypothesis: construct an initialization that resolves the
2x2 weighted Bell-label exchange blocks of (CZ_A tensor I_B)SWAP_AB, rather than
initialize only the decoder and perturb an otherwise unsolved suffix. Read-only board109 derives three H and three H S-dagger two-state blocks,
plus four singleton states, and proposes a reversible label router followed
by conditional block rotations. Its actual gate cost and realizability in
the remaining fourCX+twoF suffix are unresolved. See
[the suffix derivation](work/numerical98belllabelsuffix/REPORT.md). Any search
requires a distinct board task and frozen finite budget.
The exact12 target still requires full-list independent validation and an
exact certificate. No accepted-circuit/baseline changes occurred.

## Checks and review

Ledger/board audits passed; saved input/hash/best-selection/budget assertions
passed; all reserved scripts parse. Historical ledger rows, completed earlier
board hypotheses/events, accepted13 and exact14 bytes were preserved. Sixteen
focused ledger/board/checker tests passed in this continuation (existing SQLite
ResourceWarnings remain). Independent finite preflight and endpoint checks
are the relevant matrix validation for this new experimental family.
Adversarial review found/corrected the native-checker mismatch and caught the
missing config field before fit execution; exact runtime-version guard caught
the distribution/module mismatch before matrices. A suffix count wording error
was corrected to six native entanglers (fourCX+twoF). No unresolved blocking
artifact finding; numerical evidence remains insufficient for exact synthesis,
local minimality, or family exclusion. All research scripts/results/databases
are intentional kept artifacts; Python caches/SQLite journals remain ignored.
