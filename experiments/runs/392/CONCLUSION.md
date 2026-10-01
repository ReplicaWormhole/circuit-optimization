# Run 392 conclusion

The corrected exact-column directional finite-difference diagnostic completed with five NumPy/SciPy objective evaluations in `0.02846 s`, one thread, under the externally enforced 120-second timeout. The centered directional curvatures were `2.0650e-8` at `h=1e-3` and `8.8818e-10` at `h=5e-4`; the base loss was `0.5922239996302496`. Both estimates are roundoff-scale and far above the frozen `-1e-6` escape trigger, so they do not show a meaningful negative direction. Their signs are unresolved at this scale.

The input was exactly the run389 point and saved least eigenvector column `Q[:,0]`, with hashes and source checks recorded in `FROZEN_INPUT.json` and the ledger. The run found no candidate, made no optimization steps, and gives no full-Hessian positivity, local or global minimum, 12-CNOT exclusion, or exact certificate. No follow-up escape fit was run.
