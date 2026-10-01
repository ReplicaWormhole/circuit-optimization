# Run 394 conclusion

All 32 frozen starts completed in 117.314 seconds within the single-threaded
285-second internal and 300-second external limits. Inputs, seed rule,
optimizer settings, topology, code hash and dependency hashes are recorded in
PLAN.md, config.json, the reserved ledger row and the immutable code snapshot.
The result includes every initial and best 98-coordinate vector, evaluation
counts, termination, complete gate lists and environment versions.

Best seed: 9261616. Saved loss: 0.12499999999999986. Independent loss:
0.12500000000000003. Maximum off-diagonal entry: 0.7071067812.
The best candidate is best_compiled.json (12 chronological CNOTs).
The ledger records this candidate and terminal status `failed`.

The independent NumPy/SciPy audit checked all 64 endpoint native/compiled
lists and both global best lists, frozen RNG vectors, coordinate agreement,
compiler agreement, unitarity and objective. See VERIFIER_POSTFIT_AUDIT.json
and VERIFIER_REPORT.md. None passes diagonalization. No exact certificate
was produced and the exact 13-CNOT incumbent is unchanged.

This bounded batch does not exclude the full98 family or establish a local
minimum, topology exclusion or lower bound. No extra fits were performed.
The next proposed hypothesis changes the objective to an eigenvalue-assignment
residual while preserving all eigenvalue orders and degenerate bases. It is
unexecuted and requires a separate frozen budget and ledger reservation.
