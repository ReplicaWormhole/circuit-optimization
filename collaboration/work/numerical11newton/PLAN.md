# Fresh bounded Newton refinement of deletion9 only

Board136, numerical_sol. This is a NEW budget after unsuccessful run409,
not a restart/deepening within that run. It uses only the source409 best
deletion9 point, removing source406 gate29 (CX3->2) and retaining11CX with
all144 local coordinates free. Deletion7 is not included.

Input candidate is source409/best_candidate.json,
SHA256f67d8d40a7458f7813b3cba5e16a2e49f54a16c94ed916072e1e0bd584f1e429.
Immutable family.py is byte-identical to409 and independentpreflight408.
source_point.json records the exact144vector, topology, source result hash
and deletion identity. No matrix or optimizer is run during preparation.

Construct a480real residual by concatenating real then imaginary parts of
the240 off-diagonal entries of U V4 Udag. Compute its480x144 Jacobian by
Torch autograd. Before Newton, one deterministic normalized cosine direction
checks Jd against central differences at h1e-6; componentwise tolerance is
1e-6absolute+1e-5relative. No RNG is used (seed0 recorded).

SVD J=L diag(s) Rt; retain s>1e-8*s_max. Each minimum-norm step is
delta=-R diag(1/s) Lt residual on retained singular directions. Up to8trial
alphas2^-j,j0..7 seek a strict squared-residual decrease. At most20Newton
steps, no restart, one thread. Targetmaxoff1e-13, sharedcheckertolerance1e-9.
Internal55seconds and external60seconds bound the entire new invocation.
Save singularvalues/ranks, stepnorms, everytrial, acceptedalphas, FDcheck,
stopreason, completeinitial/bestvectors and literal11CX gate lists. A numerical
pass still needs later exactrecognition/certification.

Root reserves only after independentpostfit validation of409 passes. Guard
requires direct matching config, activeboard136ownedbynumerical_sol and run
link, parent_id409 and closedfailed409 status. Source/candidate/family/runtime/
dependency hashes and completeconfigkeys are checked. Exclusive execution
marker rejects reruns. On terminalerrors, savebestvectors/word without another
matrixcall. Root owns reservation/dispatch/ledgerfinish/boardclosure/commit.

Copy refine.py,family.py,source_point.json,initial_candidate.json,config.json,
PLAN.md to the reserved workspace. Authorized command:
`timeout 60s python3 experiments/runs/ID/refine.py`.

SVD retained rank is a numerical regularization diagnostic. Small residuals
or successfulNewtonsteps establish no exact identity, localminimum, exclusion
or lowerbound. No shared checker is edited and no dependency is installed.
