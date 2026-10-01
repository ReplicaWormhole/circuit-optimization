# Fixed parity-gauged Fourier boundary: four crossing CNOTs required

The four wires are ordered `0,1,2,3` from most to least significant. Consider
the fixed boundary operator used in `algebraic14_full_parity_boundary.py`:

1. Apply `RZZ(1,2,δ)`.
2. Apply `CX(0,2)` and `CX(3,1)`.
3. Apply the two disjoint Fourier butterflies on pairs `(0,2)` and `(1,3)`.

This is an exact commuting global-parity gauge of the known 15-CNOT
diagonalizer. The question is whether *this operator* can be synthesized with
three CNOTs so that its four-CNOT implementation can be reduced by one.

Across the bipartition `{0,1}|{2,3}`, realign its matrix as

`M[(o0o1,i0i1),(o2o3,i2i3)] = B[o0o1o2o3,i0i1i2i3]`.

With `t=exp(-iδ/2)`, the exact normalized two-qubit butterfly matrices `F`
and `R` in `topology14_boundary_schmidt_symbolic.py` give

`M[(o0o1,i0i1),(o2o3,i2i3)] = F[o0o2,i0i2] R[o1o3,i1i3] * (t if i1=i2 else t^-1)`.

An irrelevant global phase from the butterfly normalization does not affect
Schmidt rank. Symbolic evaluation gives

`det M = 1/256`,

independent of nonzero `t`. Hence the operator-Schmidt rank across this cut is
16 for every real `δ`.

A CNOT crossing this cut has operator-Schmidt rank 2; a local gate or CNOT
within either side has rank 1. Operator-Schmidt rank is submultiplicative under
operator composition. Any CNOT circuit for this fixed boundary therefore needs
at least `ceil(log2(16)) = 4` CNOTs crossing `{0,1}|{2,3}`. In particular, the
proposed three-CNOT boundary is impossible for all `δ`.

This does **not** lower-bound every diagonalizer of `V4`. A different Fourier
boundary, circuit split, or rotation inside degenerate eigenspaces changes the
target operator and may escape this obstruction.

Checks: `python3 topology14_boundary_schmidt_symbolic.py` gives exact determinant
`1/256` (ledger run 89). `python3 topology14_boundary_schmidt.py` independently
finds rank 16 on `δ=kπ/16`, `k=0,...,16` (run 86). The root agent independently
matched the symbolic operator against the actual circuit matrix at `δ=0` to
`4.3e-16` up to global phase. Runs 93 and 94 sampled one- and two-pair
degenerate-eigenspace rotations; their induced boundaries had minimum ranks 14
and 13, respectively, across this cut. These sampled results are numerical
evidence only and do not exclude other eigenbasis choices.
