# Run 391 conclusion

The five centered finite-difference evaluations completed under the external 120-second timeout and one-thread environment, but the frozen direction was the first ROW of the serialized `eigenvectors_columns` matrix. NumPy eigensolver vectors are columns, so this was not the least-eigenvalue direction. The checked direction had curvature about `0.3430388` at both steps; this result does not validate the AD least eigenvalue and must not be used as such.

Run389's `eigenvectors_columns` is a JSON encoding of a 98-by-98 matrix. The proper vector is `np.asarray(eigenvectors_columns)[:, 0]`. The exact point, mistaken direction, output, hashes, five evaluation count, and limitation are preserved. No fit or optimizer was run. A corrected test must use a new board hypothesis and ledger reservation.

Metadata note: the initial run391 ledger config likewise left the point and direction hash fields empty while recording the hash of the full frozen input. Run392 records both actual hashes directly; no prior ledger row was rewritten.
