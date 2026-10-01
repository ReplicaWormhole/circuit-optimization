# Run 398 conclusion

Both full98 Hessians completed in 4.839 seconds on one thread at exactly run
396 seed 9261703's best vector. No optimizer, perturbation or random draw was
used. The distinct-root assignment gap is 0.2133883472436626 in normalized
objective units, ignoring duplicate slots of the same eigenvalue.

Original off-diagonal loss: 0.1357233047033634; gradient norm 1.77612e-8;
least Hessian eigenvalue -1.4609760402e-9. Fixed-label branch loss:
0.1464466094067266; gradient norm 1.91640e-8; least eigenvalue -1.5584236709e-9.
Neither crosses the frozen -1e-6 negative-curvature trigger.

The independent read-only audit validates source hashes, point, complete base
gate lists, matrix/compiled equivalence, unitarity, losses, feasible row-root
exclusions and saved gap arithmetic. It also checks Hessian symmetry, all
eigenpair residuals, orthonormal eigenvector columns and canonical signs.
See POSTFIT398_AUDIT.json. That audit did not independently recalculate Hessian
entries or solve the constrained assignments. A separate finite-difference run
checks the two least directions. Run 399 stopped at a stale source guard with
zero evaluations; replacement run 400 completed all ten frozen independent
matrices. Its finite differences found no decrease at either signed step and
did not cross the trigger. Finite-step effects and cancellation leave the tiny
AD eigenvalue signs unresolved; this is not a full Hessian validation.

The base candidate, base_compiled.json, has 12 CNOTs and remains invalid with
maximum off-diagonal entry 0.3535533934. Ledger 398 is terminal inconclusive;
no trigger-based escape fit follows. Neither absence of the trigger nor a
positive numerical assignment gap proves local minimality or family exclusion.
No exact certificate or incumbent change occurred.
