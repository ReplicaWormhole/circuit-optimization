# Parity-difference controlled-basis diagonalizer

Board hypothesis 112. By-hand derivation only; no matrices, symbolic
computation, search, objective evaluation, or gate-list test was performed.

## Router and branch action

Let `A=(q0,q2)` carry `a=(a0,a1)` and `B=(q1,q3)` carry `b=(b0,b1)` in the
Bell-label coordinates. Apply the disjoint router
`R=CX(q0->q1) CX(q2->q3)`, which changes coordinates from `(a,b)` to
`(a,d)` with `d=a xor b` stored in `(q1,q3)`. From the Bell-label shift
`T|a,b> = (-1)^(b0 b1)|b,a>`, the routed operator obeys

`R T R† |a,d> = (-1)^((a0 xor d0)(a1 xor d1)) |a xor d,d>`.

Consequently each two-bit control value `d` is invariant and its target
operator on `A` is exactly `T_d = CZ_A X_A^d`, with operator order as written
(the `X` acts first). In particular this is not simply `CZ_A` with a phase
on the input: the phase uses the post-flip bits.

## Explicit branch diagonalizer with shared phase correction

Let `B=CS†_(q0->q2)` be the diagonal two-qubit phase
`B|a0,a1>=i^(-a0 a1)|a0,a1>`. By evaluating its action on a basis label and
its flipped label, the phase conjugation gives:

- `B T_00 B† = CZ_A`, still diagonal.
- `B T_10 B† = i^(a1) X_0`.
- `B T_01 B† = i^(a0) X_2`.
- `B T_11 B† = E X_0 X_2`, where `E` is diagonal with phase `i` on even
  parity `a0 xor a1=0` and phase `1` on odd parity. Equivalently
  `E=e^(i pi/4) e^(i pi Z_0 Z_2/4)`.

Thus explicit fixed-branch diagonalizers are especially simple:
`W_00=I`, `W_10=H_0`, `W_01=H_2`, and
`W_11=H_0 CX(q0->q2)`. For the last branch, the CNOT sends
`X_0 X_2` to `X_0` and `Z_0 Z_2` to `Z_2`; the final Hadamard sends `X_0`
to `Z_0`, so both factors become diagonal. The factors `i^(a1)` and
`i^(a0)` in branches 10 and 01 are already diagonal and commute with the
respective bit flips.

The full basis change after routing is the controlled unitary
`C = sum_d |d><d|_(q1,q3) tensor W_d`. There is a useful shared-selector
factorization. Let `P` be `CX(q0->q2)` controlled only on `d1=q3=1`. On
branch `d=01`, this leaves `i^(a0) X_2` unchanged because its diagonal phase
depends only on the control q0 and `X_2` itself is invariant under the CNOT
conjugation rule. On `d=11`, it maps `E X_0 X_2` to
`i^(1-a1) X_0`. Thus after `P`, apply `H_0` controlled on `d0=q1=1`, and
`H_2` controlled on `d=01` (q1=0,q3=1). These three operations implement
all four branch diagonalizers while sharing the `d1` selector. Since
arbitrary output-label order is allowed, `C B R` is sufficient after the
Bell prefix; there is no need to append `R†`. The branch spectrum has
multiplicities `+1,-1,+i,-i = 6,4,3,3`, consistent with the Bell-label
derivation.

## Cost and fit to the reserved suffix

The unconditional phase correction `B=CS†` has an exact two-CX realization.
The router costs two CX. The specialized exact Bell decoder for each pair has
matrix `H_c CX(c->t)` (chronological gates: first CX, then H), so the
two-pair prefix costs two CX total. Together the prefix, router and `B` use
six CX, leaving six for `C` under a twelve-CX target. This specialized-prefix
accounting is distinct from compiling the two free `F` gates of a generic
variational member at two CX each.

The shared-selector `P` is a Toffoli with controls q3 and q0 and target q2.
The standard exact no-ancilla Toffoli decomposition uses six CX by itself;
the two controlled Hadamards still need implementation. So that literal
exact-Toffoli-plus-separate-H construction does not fit the six-CX remainder.
This is not an impossibility proof: relative-phase Toffoli behavior may be
usable because diagonal output phases are free, and a merged decomposition may
share gates with the controlled Hadamards. No such <=6-CX decomposition has
been supplied or certified here.

## Actionable handoffs

The equations, branch ordering and honest cost scope were sent to the
coordinator and the numerical/verifier peers. The useful next synthesis
question is whether the `d1`-controlled Toffoli and the two controlled
Hadamards admit a phase-aware merged decomposition within six CX. A candidate
must be supplied as a complete gate list and independently checked; numerical
or structural derivation alone is not an exact certificate.
