# Run 390 conclusion and correction

No objective evaluation occurred. The frozen input's `least_eigenvector` was mistakenly copied from the first ROW of the serialized `eigenvectors_columns` matrix. The script guard correctly compared against the first COLUMN (`Q[:, 0]`) and refused to run. The external 120-second timeout and single-thread environment were active; elapsed time was below one second. This is a failed preflight, not numerical evidence. The script and frozen input are preserved byte-for-byte.

Correction to the initial ledger note: the guard was correct; the input vector orientation was wrong. The error is documented here and in the board handoff. No ledger row was rewritten. Run391 then mistakenly changed the guard to accept the row and evaluated that wrong direction; its five curvature evaluations do not validate the AD least eigenvalue and are recorded separately. A proposed corrected input now uses exact `Q[:,0]`; it has not been reserved or evaluated pending coordinator review.

Metadata note: the first run390 ledger config contains empty `point_sha256` and `direction_sha256` fields, but it records the hash of the complete frozen input file containing those values. The corrected run392 ledger records both hashes directly. The run390 ledger row was preserved without rewriting.
