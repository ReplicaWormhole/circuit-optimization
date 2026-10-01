# Exact 14-CNOT right-cycle diagonalizer

The full rational-angle gate list is
[`topology14_exact_matchgate_rational.json`](topology14_exact_matchgate_rational.json).
It applies the first 13 gates of `topology16_15_exact_candidate.json` (six
CNOTs), then the following eight-CNOT tail. Wires are numbered `0,1,2,3`
from most significant to least significant bit; gates are chronological.

For each wire, `ZYZ(theta,phi,lambda)` means chronological
`Rz(lambda), Ry(theta), Rz(phi)`. The input local layer is

| wire | theta | phi | lambda |
|---|---:|---:|---:|
| 0 | 3π/4 | 0 | 0 |
| 1 | π/2 | −π/2 | −π/3 |
| 2 | π/2 | −3π/2 | π/6 |
| 3 | 3π/4 | −π | π |

Apply the disjoint two-qubit blocks `F(0,2;π/4,π/8)` and
`F(1,3;π/4,π/8)`, where

`F(c,t;a,b) = exp(i a X_c X_t + i b Y_c Y_t)`.

The center local layer is

| wire | theta | phi | lambda |
|---|---:|---:|---:|
| 0 | π/2 | 0 | −π |
| 1 | π | 0 | 3π/2 |
| 2 | π | 0 | 3π/2 |
| 3 | π/2 | −3π/2 | 0 |

Finally apply `F(0,1;π/8,π/8)` and `F(2,3;π/8,π/8)`. No output
single-qubit gates are needed.

## Two-CNOT identity

Let `B = Rx_c(π/2) Rx_t(π/2)`. Each `F(c,t;a,b)` is implemented by the
chronological sequence

1. `Rx_c(−π/2), Rx_t(−π/2)`;
2. `CX(c,t)`;
3. `Rx_c(−2a), Rz_t(−2b)`;
4. `CX(c,t)`;
5. `Rx_c(π/2), Rx_t(π/2)`.

Indeed, `CX X_c CX = X_c X_t` and `CX Z_t CX = Z_c Z_t`. Conjugation by
`B` fixes `X_c X_t` and maps `Z_c Z_t` to `Y_c Y_t`; the two Pauli
products commute. Thus each block costs exactly two CNOTs, giving
`6 + 4×2 = 14` CNOTs.

## Verification and status

`topology14_exact_matchgate_candidate.py` generates the rational gate list.
`python3 exact_check.py topology14_exact_matchgate_rational.json` proves
`U V4 = D U` over `Q(ζ_96)` with output eigenvalue labels
`[0,3,1,1,2,0,2,2,0,2,0,0,3,1,3,0]`; it reports 14 CNOTs. This is
an exact algebraic certificate of right-cycle diagonalization, independent
of the numerical search that suggested the angles. The independent
`delete14_integer_audit.py` checks the same identity using integer
polynomials modulo $z^{32}-z^{16}+1$ and no floating-point arithmetic.
`algebraic14_raw_matchgate_certificate.py` separately checks the four
XX/YY blocks and their equality to the 14-CNOT gate list up to global
phase. Qiskit independently simulates the same list with input
off-diagonal error `4.83e−16` and transpiles it to 14 CNOTs. The local
checker reports `6.70e−16`.

This circuit diagonalizes the right-cycle `V4`. It does not diagonalize
total spin; the local checker reports total-spin off-diagonal error
`2.6457513110645867`. No lower bound or optimality claim is made.
