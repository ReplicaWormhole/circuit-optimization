# Run 362: fresh π/24 Euler gauge attempt

Reset from original run359 numerical candidate; the previous run361 π/8 constraints were discarded. Same thirteen chronological CNOTs, four qubits without ancillas, q0 most significant, UV4=DU with fixed original output eigenvalue labels and free sector bases.

Command: `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 collaboration/work/algebraic13_new/fine24_gauge_snap.py`.

Completed the authorized twelve stages, at most nine prescribed coordinates per stage, at most400 residual evaluations per stage and one thread. Actual stages used25–41 evaluations. All twelve fits passed.108 Euler coordinates were prescribed nearest π/24 values; a total115 coordinates recognize rational multiples of π with denominator at most96 within1e-10, leaving53 generic. This tolerance-based recognition is numerical evidence.

`fine24_gauge_candidate.json`: checker passes thirteen CNOTs, maximum off-diagonal error8.081e-14, eigenvalue-root error1.155e-14. Source run359 remains the more accurate numerical candidate.

`fine24_rounded_candidate.json`: maximum off-diagonal error0.2548 and exact checker returns `exact_diagonalizer:false`, conductor96. Exact rejection command: `python3 exact_check.py collaboration/work/algebraic13_new/fine24_rounded_candidate.json`.

`fine24_result.json` records per-stage selected coordinates, Jacobian ranks, residuals, evaluations and low-denominator recognition. The Jacobian rank was56 at all stage starts. Some resulting Euler theta coordinates are0 orπ, creating chart degeneracies in phi/lambda; generic coordinate values do not directly demonstrate necessary noncyclotomic gates. Neither rational-angle existence nor a thirteen-CNOT impossibility is established.

Adversarial review: fresh source reset and denominator replacements inspected; original and old run361 artifacts preserved; no shared checker, native code, status documents or dependencies changed. Numerical and exact checks separated. Python syntax compilation passes. All fine24 artifacts are intentional retained research evidence. Run362 and hypothesis8 closed inconclusive; no exact certification or incumbent promotion. Coordinator owns commits.
