# Static review of the reordered Bell-decoder initializer

Hypothesis104. This is a by-hand review of the existing Bell-prefix derivation
and proposed full98 initializer. No matrix multiplication/evaluation, symbolic
or numerical search, optimizer, derivative, or certificate was performed.

## Conventions

Use pair `A=(q0,q2)` with control `q0`, target `q2`, and pair
`B=(q1,q3)` with control `q1`, target `q3`. Chronological gates are written
left-to-right in the schedule, so their total matrix product acts with the
rightmost factor first. The native interaction convention is
`F_ct(a,b)=exp(i a X_c X_t+i b Y_c Y_t)`. Single-qubit rotation vectors
represent `exp[-i(xX+yY+zZ)/2]`, up to scalar phase when comparing to U(2)
gates.

The existing encoder factorization is

`E_ct = K_ct F_ct(pi/4,0) P_ct`,

`P_ct=R_y,c(pi/2)H_c = X_c`,
`K_ct=(S_c R_y,c(-pi/2)) tensor R_x,t(pi/2)`.

Since the Bell states are columns of `E_pair=E02 tensor E13`, the useful
prefix for target conjugation is its decoder. Indeed, with `U=W E_pair†`,
`U V4 U†=W(E_pair† V4 E_pair)W†`. Taking the adjoint in matrix-product
order gives

`E_ct†=(X_c tensor I_t) F_ct(-pi/4,0)
       ((R_y,c(pi/2)S_c†) tensor R_x,t(-pi/2))`.

Thus the sign of both prefix interactions is **negative**: `F02` and `F13`
are each `F(-pi/4,0)`. The post-entangler Pauli is `X` on the pair's control.

## Rotation-vector check by hand

The target input local is `R_x(-pi/2)`, with SO(3) rotation vector
`(-pi/2,0,0)`. The post local `X=exp(-i pi X/2)` uses `(pi,0,0)` on the
control.

For the control input local,

`S† = e^(-i pi/4) R_z(-pi/2)`,

so up to global phase

`R_y(pi/2)S† ~ R_y(pi/2)R_z(-pi/2)`.

The SU(2) quaternion product, in scalar/vector convention
`q0 I - i(qx X+qy Y+qz Z)`, is

`(1/sqrt(2),0,1/sqrt(2),0) *
 (1/sqrt(2),0,0,-1/sqrt(2))
 = (1/2,-1/2,1/2,-1/2)`.

This is the axis-angle rotation with angle `2pi/3` and unit axis
`(-1,1,-1)/sqrt(3)`, hence rotation vector

`(2pi/(3sqrt(3))) * (-1,+1,-1)`.

This agrees with the proposed control vector up to the stated global phase.
The local in the decoder formula is the matrix product `R_y(pi/2)S†`; if
serializing as successive gates, apply `S†` first and then `R_y(pi/2)`.

## Exact schedule placement

A nonduplicating realization of the two decoder prefixes in
`L0 F02 L1 F13 L2` is:

- `L0`: on both controls `q0,q1`, use rotation vector
  `(2pi/(3sqrt(3)))(-1,+1,-1)`; on both targets `q2,q3`, use
  `(-pi/2,0,0)`; identities on any other local slots (there are no other
  wires here).
- `F02=F(-pi/4,0)`.
- `L1`: apply only `X_0`, rotation vector `(pi,0,0)` on `q0`; identity on
  `q1,q2,q3`.
- `F13=F(-pi/4,0)`.
- `L2`: apply only `X_1`, rotation vector `(pi,0,0)` on `q1`; identity on
  `q0,q2,q3`.

The pair-B input locals can be in `L0` even though `F13` comes later, because
`F02` acts on disjoint wires and commutes with all pair-B locals. Likewise
`X_0` in `L1` commutes with the disjoint `F13`. The product therefore groups
into `E02† tensor E13†` without local duplication. The claimed schedule
placement and signs follow from the algebraic factorization; this review did
not instantiate or numerically compare gate matrices.

This fixes an initializer only. The remaining local coordinates, `A1/B1`, the
three-CX star, and tail gates must remain free for a full coupled98 fit. It is
a new reordered family, not a candidate, and it does not prove the suffix can
diagonalize the Bell-label exchange operator. No existing family is excluded.
