# Run 389: Full coupled 98-coordinate Hessian

Ledger record: `python3 experiment_log.py show 389`. Board hypothesis: 68. Parent ledger run: 387 (frozen point derived from run385 seed 9261101).

## Hypothesis and scope

Evaluate one full automatic-differentiation gradient and symmetric 98x98 Hessian for the equivalent family with 4 native CX and 4 native XX+YY gates (compiled cost 12 CX). This probes all local, retained target-interior, and F-angle mixed directions at the frozen run385 seed9261101 point. The frozen negative-curvature trigger is lambda_min < -1e-6. No optimization, perturbed fit, restart, or conditional run is included.

## Frozen inputs and limits

- Source point: `run387_frozen_parent.json`, ledger run387 `config.point116`, tracing to source run385 seed9261101.
- Coordinates: 84 rotation-vector values in local-layer order [L0prime,L1prime,L2,L3,L4,L5,L6] (q-major within each layer), A1 on q1 (3), B1 on q1 (3), then F(a,b) for pairs [(0,2),(1,3),(0,1),(2,3)] (8).
- Chronology: L0prime, CX01, A1, CX21, B1, CX31, L1prime, F02, L2, F13, L3, CX12, L4, F01, L5, F23, L6.
- Exact boundary map from [L0,A,B,L1,L2..L6]: L0prime(2)=A2 L0(2), L0prime(3)=B3 A3 L0(3), L1prime(0)=L1(0) B0 A0, L1prime(2)=L1(2)B2. A1/B1 remain interior. For this source seed A=B=I.
- Seed 9261101; one gradient; one 98x98 Hessian; one thread; hard external timeout 240 seconds.
- Frozen code SHA-256 `739c8585bcb5b0146a3d735dd06dc2c0be65a7eaf9c923ffffa766303e741957`; frozen config SHA-256 `7698988d1a8c6d8eee0ff63c0d98dbf9345217872a7a27340cdae8d7dbe6d4e6`.
- Independent prefit approval: verifier98 reports source remap error 4.61e-16, full98 direct-vs-serialized error 5.24e-16, native-vs-compiled error 3.86e-16, 4 native CX + 4 F, compiled 12 CX.

## Reproduction command

Reproduce only in a fresh copy of this workspace, preserving this completed output directory: `timeout --signal=TERM --kill-after=5s 240s env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 hessian98.py`. The script refuses to overwrite an existing `result.json`.

## Outputs and independent checks

`result.json` saves input point, full gradient, full Hessian, eigenvalues/eigenvectors, residuals, elapsed time, serialized base gate lists and checker results. `base_native.json` and `base_compiled.json` are included even though the point is not a valid diagonalizer. The assigned verifier independently audits the source embedding and gate-list matrices.

## Conclusion and limitations

Close the run according to the frozen threshold. A failed trigger means only that this selected floating-point Hessian did not authorize follow-up fits; it is not a local-minimum proof or exclusion of 12-CX solutions. Any passing gate list requires independent exact certification before changing the 13-CX incumbent.
