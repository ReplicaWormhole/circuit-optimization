# Run401 unexecuted static failure

The reserved runner accesses cfg["noise_sigma"], but frozen config contains
that value only under initialization. Root and independent verifier identified
the defect before execution. No RNG start, matrix, objective, derivative or
optimizer was evaluated. Preserve all frozen inputs, snapshot and ledger row.
Ledger401 is failed; board105 inconclusive. A corrected config is a new attempt,
not an edit or retry under this ID. No candidate or scientific conclusion follows.
