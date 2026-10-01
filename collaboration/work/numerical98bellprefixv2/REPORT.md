# Corrected config snapshot for the reordered Bell-decoder full98 fit

Board 106. This separate directory preserves the failed-to-preflight 401
snapshot while correcting its config schema. The source and family are byte
identical to the prior snapshot; only the copied config adds the required
top-level `noise_sigma: 0.1` accessed by `search.py`. The nested initialization
record already had the same value. No numerical vector was generated and no
matrix, objective, optimizer, or fit was run.

Artifacts and SHA256:

- `search.py`: `3f4777c2492e34e47fb5655df4c288f432f8dbad0e74230afeff872c257448ac`
- `family.py`: `9936da7408d5202a00e784daa4b33d950efa68229bec06c006457862411303c8`
- corrected `config.json`: `1fb2864e4f9750324080d95f8193d46dce75c3939373036b4187b8e080c4616a`

Static validation parsed both Python files and JSON, checked every literal
top-level `cfg[...]` key used by the runner exists, asserted top-level and
nested sigma are both 0.1, verified code/family hashes, and verified each
external dependency hash. All passed. The single start and limits are
unchanged: seed 9261801, Bell-decoder base plus normal noise on all 98
coordinates, one L-BFGS-B run capped at 1000 iterations/20000 evaluations,
110-second internal/120-second external limit, one thread.

Root owns the new ledger reservation, verifier preflight, and eventual
execution. The corrected snapshot is ready to copy into a new reserved numeric
workspace; it must not be substituted into or written over run 401.
