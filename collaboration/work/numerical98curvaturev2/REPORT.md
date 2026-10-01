# Static re-review of repaired run 398 curvature diagnostic

Board 99. Read-only source/config review; no matrix, objective, assignment,
gradient, Hessian, optimizer, or certificate calculation was run.

Frozen artifacts:

- Script: `experiments/runs/398/diagnostic.py`
- Config: `experiments/runs/398/config.json`
- Code SHA256: `5e67087e57f9b7d3a8915ac41339db7877b22df958a1c60e67896456cac9b086`
- Config SHA256: `a5e6ecbf165e548321c5b178edb0c4d2941cfa72dfbce2179fd89b95acaeae13`

Both hashes match the coordinator's frozen values, and a syntax-only AST parse
passed. The diagnostic loads the exact run396 best vector for seed 9261703 and
uses it for the original off-diagonal objective and the fixed optimal-label
branch. It performs no fitting, perturbation, or random draw. The stated limits
match the code: at most two first gradients, 196 second-derivative rows total,
and 17 Hungarian assignments. It sets Torch and BLAS thread limits to one.

The prior timeout issue is repaired. The assignment-gap and both Hessian-row
loops are inside one `BudgetStop` handler. Assignment alternatives and the gap
are initialized before work; each objective record is appended before computing
its second-derivative rows, and completed rows plus their count are retained.
After an internal timeout the base native/compiled lists and partial result are
still serialized. The constrained assignments ban every duplicate slot of the
selected root in one row, preserving capacities elsewhere. The raw row-residual
cost is the complete labeled residual. Hessians are symmetrized, eigenvectors
are stored by columns with canonical signs, and the output scope avoids claims
of a local minimum, exclusion, or certificate.

No blocking static issue found. The internal deadline is cooperative at
assignment and second-derivative row boundaries; eigendecomposition and
serialization/checking are atomic phases bounded by the external 120-second
timeout, as the plan states. Findings were sent directly to root,
verifier98fresh, and numerical98. This review provides no evidence about the
computed curvature or the circuit family.
