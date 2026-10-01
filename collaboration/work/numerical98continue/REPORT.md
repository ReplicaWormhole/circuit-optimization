# Full98 continuation strategy review

Board hypothesis 79. This assignment is read-only review of the saved run393 failure and recommendation of one distinct follow-up after the coordinator's separately frozen 32-start batch. No fit, derivative, matrix evaluation, symbolic calculation, or exact certificate was run here.

## Evidence from run393

Run393 used one fresh all-98 seed and completed 284 L-BFGS-B iterations / 301 objective-gradient evaluations in 3.439 seconds, with no restarts. Its normalized Frobenius objective improved from 0.908501360688 to 0.454898570731. Termination was by relative objective reduction; the final gradient norm, 3.8973e-8, remains above the requested 1e-10 tolerance.

The maximum off-diagonal entry moved in the opposite direction: 0.529568978351 initially and 0.674363073634 at the best Frobenius point. Both saved circuits are invalid. The independent audit confirmed all four coordinate/native/compiled lists and the 4CX+4XX/YY to 12CX count. This is a bounded single-start result, not a family exclusion, local-minimum result, basin claim, or exact certificate.

## Run394 endpoint distribution

Run394 finished all 32 starts in 117.3 seconds. Its best normalized Frobenius objective is approximately `0.125`, selected from seed 9261616 (sigma 0.3); the saved checker reports that endpoint invalid with maximum off-diagonal entry `0.70710678`. Several other starts reached the same `0.125` plateau and the same maximum-entry scale. A distinct endpoint from seed 9261620 (sigma 0.3) has loss `0.1357233` but a substantially smaller maximum off-diagonal entry `0.35355339`. Thus ranking only by the current Frobenius loss selects a visibly worse maximum-entry residual, repeating the metric mismatch seen in run393. Verifier post-fit audit should be consulted before any endpoint is used as a polish source.

## One recommended next diagnostic

After independent review of run394, use the best audited low-maximum-entry endpoint (currently seed 9261620) for one separately frozen full98 polish with a dynamic eigenlabel objective:

`min_D ||U V4 - D U||_F^2 / 16`,

where `D` is diagonal with the exact cycle spectrum and multiplicities `+1` (6), `-1` (4), `+i` (3), `-i` (3). At each evaluation, choose the minimum-cost assignment of these repeated eigenvalues to diagonal entries by Hungarian assignment; form the gradient with the chosen `D` held fixed for that evaluation. Fix a deterministic tie rule and record assignments because label changes make the objective piecewise smooth. Keep all 98 coordinates free and compare final candidates with the original off-diagonal objective and independent checker; any passing candidate still needs exact certification.

This is a distinct diagnostic motivated by the observed spread: the Frobenius-best endpoints stall at a `0.125` plateau with maximum error `sqrt(1/2)`, while another endpoint has a higher Frobenius loss and half the maximum error. A minimax-only polish was the initial suggestion before examining the batch; the known-spectrum assignment objective is the recommended next test because it directly includes both diagonal spectral mismatch and off-diagonal residual while preserving arbitrary eigenvalue order and degenerate multiplicities. This remains an unexecuted hypothesis, not a claim that it will escape or reach a solution.

## Batch review status

The frozen batch artifacts are at `experiments/runs/394/{batch.py,config.json,PLAN.md}`. Static review found the parameter assembly (90 local/interior normal coordinates plus 8 uniform F angles), seed/sigma sequence, per-start iteration/evaluation caps, shared monotonic deadline, pre-call maxfun guard, best-sample retention, per-start checkpointing, and complete native/compiled endpoint saving consistent with the stated budget. A start interrupted by the internal deadline retains its best evaluated point before the batch stops. No script was run or modified.

The independent verifier's post-fit report is retained under the run394 workspace. Numerical98 also supplied a distinct structural direction: regroup wires as `A=(q0,q2)` and `B=(q1,q3)`, express the cycle as `(S_A tensor I) SWAP_AB`, diagonalize the within-pair swap in a Bell basis, then synthesize the resulting unordered-label 2x2 eigentransforms. The match between that target basis and the fixed local+F schedule is not established; this is a research direction, not a construction claim.
