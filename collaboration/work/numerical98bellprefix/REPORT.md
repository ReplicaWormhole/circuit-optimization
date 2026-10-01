# Reordered Bell-decoder full98 preparation

Board 105. This directory contains a frozen one-start search proposal; root
owns its ledger reservation, run workspace, verifier gate, and execution. No
matrix/objective evaluation, initialization-vector draw, optimizer, or fit was
run during preparation.

The changed entangler chronology is
`L0 F02 L1 F13 L2 CX01 A1 CX21 B1 CX31 L3 CX12 L4 F01 L5 F23 L6`.
The 98-coordinate layout remains 84 local rotation-vector coordinates, six
free A1/B1 rotation coordinates, and eight F angles. `family.py` implements
the differentiable matrix and matching native serializer. `search.py` makes
one deterministic start: the Bell-decoder prefix described in
`algebraic98bellprefix/REPORT.md` plus `default_rng(9261801).normal(0,0.1,98)`
on every coordinate. All 98 coordinates are passed freely to one Torch-jacobian
SciPy L-BFGS-B fit.

Frozen artifacts and hashes:

- `search.py`: `3f4777c2492e34e47fb5655df4c288f432f8dbad0e74230afeff872c257448ac`
- `family.py`: `9936da7408d5202a00e784daa4b33d950efa68229bec06c006457862411303c8`
- `config.json`: `e2b8583b9a53f947055d60520d8282f72b96ce66ce81eb3b1f1b6180ebbbbce6`

The script checks the frozen code, family, algebra-report, and dependency
hashes before work. It refuses a preexisting result or execution marker, saves
initial native/compiled lists before optimization, and retains the best
evaluated point including line-search samples. The result includes the initial
and best vectors, optimizer status/counts/history, separately counted optional
single best-loss recomputation, gate-list checker outputs, elapsed time, and
runtime versions. A failed fit is only bounded evidence for this one reordered
family and start; it cannot establish a lower bound, family exclusion, or an
exact 12-CX candidate.

Only syntax and JSON parsing were checked here. Independent finite preflight
and root execution remain pending; those steps must verify the Bell-prefix
coordinates, chronology, and both serializers before the one reserved fit.
