# Pre-star opposite-pair Bell encoders

Hypothesis100; read-only structural proposal for the reordered full98 word

`L0 F02 L1 F13 L2 CX01 A1 CX21 B1 CX31 L3 CX12 L4 F01 L5 F23 L6`.

It retains four directed CX, four native `F=exp(i a XX+i b YY)` gates, seven
four-qubit local layers, and the six interior `A1/B1` coordinates, for the same
98-coordinate count and nominal 12-CX compilation. This is a new ordered
family; it neither searches nor excludes the old word. No numerical or symbolic
matrix evaluation, optimizer, or exact certificate was performed.

## Exact local-plus-F Bell encoder

For an ordered pair `(c,t)` with control/first wire `c` and target/second wire
`t`, define

`F_ct(a,b)=exp(i a X_c X_t+i b Y_c Y_t)`,

`R_y(theta)=exp(-i theta Y/2)`, `R_x(theta)=exp(-i theta X/2)`, and
`S=diag(1,i)`. Then, up to an irrelevant overall phase,

`CX_(c->t) = (S_c tensor R_x,t(pi/2))
              (R_y,c(-pi/2) tensor I)
              F_ct(pi/4,0)
              (R_y,c(pi/2) tensor I)`.

Indeed, conjugating `XX` by `R_y,c(-pi/2)` gives `ZX`, so the middle
entangling part and target rotation combine to

`exp(-i pi X_t/4) exp(i pi Z_c X_t/4)
 = exp[-i pi (I-Z_c)X_t/4]`.

On control zero this is identity; on control one it is `-i X_t`. Multiplying
by `S_c` changes the latter to `X_t`, yielding CX. Consequently the encoder unitary
that maps computational pair labels to the Bell ordering

`00 -> (|00>+|11>)/sqrt(2)`,
`01 -> (|01>+|10>)/sqrt(2)`,
`10 -> (|00>-|11>)/sqrt(2)`,
`11 -> (|01>-|10>)/sqrt(2)`

is `E_ct=CX_(c->t)(H_c tensor I_t)`, and one explicit factorization is

`E_ct = (S_c tensor R_x,t(pi/2))
        (R_y,c(-pi/2) tensor I)
        F_ct(pi/4,0)
        (R_y,c(pi/2) H_c tensor I_t)`.

Here `H` is the Hadamard. The displayed order is matrix-product order, so the
rightmost local factor acts first. The `S` and `H` are single-qubit gates;
scalar phase is immaterial to the diagonalization equation. This gives an explicit pair encoder factorization `E_ct=K_ct F_ct(pi/4,0) P_ct`,
where `P_ct=(R_y,c(pi/2)H_c) tensor I_t = (X_c tensor I_t)` and
`K_ct=(S_c R_y,c(-pi/2)) tensor R_x,t(pi/2)`. These two encoders act on
disjoint registers. The target-conjugation orientation is important: if `E`
has the Bell states as its columns, then the Bell-coordinate target is
`V_B=E_pair† V4 E_pair`, not `E_pair V4 E_pair†`. Thus for a total
`U=W E_pair†`, applied chronologically as the decoder prefix and then `W`,
`U V4 U†=W V_B W†`. The decoder is

`E_ct† = (X_c tensor I_t) F_ct(-pi/4,0)
         ((R_y,c(pi/2) S_c†) tensor R_x,t(-pi/2))`.

In chronological form on each pair, apply the input locals
`R_y,c(pi/2)S_c†=U3(pi/2,0,-pi/2)` and `R_x,t(-pi/2)`, then
`F_ct(-pi/4,0)`, then `X_c`. One precise schedule for both decoders in the
ordered word `L0 F02 L1 F13 L2` is: put both pairs' input locals in `L0`, set
`L1` to `X_0` only, and set `L2` to `X_1` only; all other local factors in
these three layers are identity. The `(q1,q3)` input locals may be placed in
`L0` because they commute with the disjoint earlier `F02`. Use
`F02=F13=F(-pi/4,0)`. This realizes `E02†` followed by `E13†` at the prefix
level without duplicating any local factor.

## Why moving them before the CX fan-in could help

Group registers as `A=(q0,q2)` and `B=(q1,q3)`. The prior exact algebra gives
`V4=(S_A tensor I_B) SWAP_AB`, where `S` swaps the two bits inside a register.
For the Bell encoder `E_pair` whose columns are the Bell states, conjugating
into Bell coordinates gives `E_pair† V4 E_pair`; this diagonalizes the
internal bit swap and yields

`(E_A tensor E_B)† V4 (E_A tensor E_B) = (CZ_A tensor I_B) SWAP_AB`.

Thus in that basis the target does not require an arbitrary dense eigenvector
map. It consists of a label exchange with a known phase sign. Equal-sign
label pairs can be resolved into the `+1/-1` combinations and opposite-sign
pairs into `+i/-i` combinations; the previous derivation gives the exact
blocks and multiplicities `(6,4,3,3)`.

The proposed reordering places the pair Bell decoders `E_pair†` at the start,
before the three-CX star `CX01, CX21, CX31`. In the total matrix `U=W E_pair†`,
this conjugates the target to its Bell-coordinate form before `W` acts. This
gives the star a plausible role in routing/mixing the Bell labels into the
two-dimensional exchange blocks,
while the remaining `CX12`, `F01`, and `F23` tail can plausibly complete the
label-pair combinations and phase cleanup. In the old order, the first three
CX act before the opposite-pair Bell decoders, so this direct “change to Bell
coordinates, then exchange labels” factorization is obscured. That is the motivation for testing
a distinct word, not evidence that the tail already implements the needed
transform.

## What remains unproved

- The Bell identity for `V4` and local-plus-`F` encoder/decoder formulas are
  exact algebraic gate identities. The two-decoder prefix fits the stated
  chronology using `L0,F02,L1,F13,L2`.
- No factorization has been found for the remaining label-exchange eigenbasis
  in the exact `CX01 A1 CX21 B1 CX31 ... CX12 F01 F23` suffix. In particular,
  the six `A1/B1` rotations and intervening local layers must realize the
  controlled permutation/sign structure without undoing the Bell encoding.
- It is not shown that the fixed 98-coordinate chart can absorb all required
  local factors with the stated chronology, nor that the remaining tail can
  synthesize the exchange eigenvectors at four native CX plus two native F.
  Prefix encodability alone is not a complete diagonalizer.
- This does not prove a valid 12-CX circuit, an exact certificate, or a lower
  bound. Failure of this reordered family would not exclude the previous
  ordering or other 12-CX schedules.

Actionable next step after a separately frozen budget: initialize the
reordered full98 family from this exact decoder prefix and keep all 98
coordinates free during the bounded fit. The prefix is an initialization, not
a permanently constrained prefix or a replacement for the full coupled-family
investigation. Any full candidate still needs an independent gate-list check
and exact certificate.
