# Explicit 19-CX parity-difference baseline preparation

Board 116. This report records a finite candidate preparation, not a
verification result. I have not evaluated a circuit matrix or run the
numerical or exact checker. The checker is guarded to run only from a
matching active reserved workspace.

## Gate construction

The candidate combines the Bell-coordinate decoder, the parity relabeling,
and the common-phase branch diagonalizer derived in
`../verifier98difference/REPORT.md`. Qubit 0 is the most significant bit and
all listed gates are chronological.

1. Bell decoder on pairs (0,2) and (1,3): `CX02, H0, CX13, H1`.
2. Difference register `(a,b) -> (a,d=a xor b)`: `CX01, CX23`.
3. `CS†` on (0,2), with diagonal `diag(1,1,1,-i)`, using local phase
   gates `P(-pi/4)` on both wires, then `CX02, P(+pi/4)_2, CX02`.
4. Apply `CX(0->2)` controlled by `d1=q3`, using the standard exact
   six-CX Toffoli decomposition with controls 0,3 and target 2.
5. Apply `H0` controlled by `d0=q1`. The chronological decomposition is
   `Ry0(-pi/4), H0, CX(1->0), H0, Ry0(+pi/4)`: the active control branch is
   `Ry(pi/4) Z Ry(-pi/4)=H`, while the inactive branch is identity.
6. Apply `H2` controlled by `d1=1,d0=0`. Conjugate `CCX(1,3->2)` by
   `V=Ry2(-pi/4)`, where `V X V†=H`; wrap it in `X1` before and after to
   invert the q1 control.

The CX count is `2 + 2 + 2 + 6 + 1 + 6 = 19`. The list uses only the exact
checker primitives `cx`, `h`, `x`, `ry`, and rational-pi `u3` phase gates.
All selector gates are controls on the parity-register bits produced by
step 2. The route is an explicit exact comparator hypothesis, not a 12-CX
candidate.

## Planned checks and limits

After root reserves the check, `check_baseline.py` will validate source and
runtime hashes, require an active ledger row whose code hash matches the
script, save the constructed candidate, and run each checker once with one
thread: the NumPy native-gate checker has a 30-second timeout and the exact
cyclotomic checker a 60-second timeout. Timeout outputs and return codes are
retained. No optimizer or search is present.

Passing both checks would establish an exact 19-CX example for this route,
not an incumbent change, a 12-CX construction, or an optimality statement.
