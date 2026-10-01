# Opposite-pair Bell encoding for the full98 family

Hypothesis81; read-only structural derivation. No optimizer, numeric matrix
evaluation, symbolic search, circuit fitting, or exact certificate was run.
The target is the right cyclic shift
`V4 |x0,x1,x2,x3> = |x3,x0,x1,x2>` with q0 most significant.

## Register factorization

Regroup the four wires into two ordered two-qubit registers
`A=(q0,q2)` and `B=(q1,q3)`. For computational strings
`A=(a0,a1)` and `B=(b0,b1)`, the right shift gives

`(A,B) -> ((b1,b0), A)`.

Let `S` swap the two bits within a register, and let `SWAP_AB` exchange the
two registers. Then, exactly in this regrouped tensor order,

`V4 = (S_A tensor I_B) SWAP_AB`.

There is no fermionic sign in this identity: it is just the permutation of the
four qubit tensor factors.

Use the Bell basis on each register:

| label | Bell state | eigenvalue of S |
| --- | --- | --- |
| 00 | `(|00>+|11>)/sqrt(2)` | +1 |
| 01 | `(|01>+|10>)/sqrt(2)` | +1 |
| 10 | `(|00>-|11>)/sqrt(2)` | +1 |
| 11 | `(|01>-|10>)/sqrt(2)` | -1 |

With this ordering, the pair-swap matrix is `CZ=diag(1,1,1,-1)`. Since the
same Bell transform is applied to both registers, it commutes with their
register exchange. In Bell-label coordinates `a,b in {00,01,10,11}`, the
conjugated right shift is therefore

`V' = (CZ_A tensor I_B) SWAP_AB`,

and its action is

`V'|a,b> = sigma_b |b,a>`, where `sigma_a=(-1)^(a0*a1)`.

## Explicit eigenbasis and degeneracies

For `a=b`, `|a,a>` is already an eigenvector with eigenvalue `sigma_a`. For
each unordered pair `a != b`, the span of `|a,b>` and `|b,a>` is invariant,
with matrix

`[[0, sigma_a], [sigma_b, 0]]`

in that ordered basis. Its eigenvalues obey `lambda^2=sigma_a sigma_b`, and
normalized eigenvectors are

`(|a,b> + (lambda/sigma_a)|b,a>)/sqrt(2)`.

Thus equal-sign labels give eigenvalues `+1,-1`, while opposite-sign labels
give `+i,-i`. The diagonal label states contribute three `+1` and one `-1`;
the three unordered pairs among the three `+1` labels contribute three more
of each `+1,-1`; the three pairs between a `+1` and `-1` label contribute
three each of `+i,-i`. Total multiplicities are `(6,4,3,3)`, matching the
right-cycle spectrum. This gives a concrete eigenbasis with explicit freedom
to choose order and phase within degenerate eigenspaces.

## Relation to the fixed 4CX+4F word

The full98 word contains `F02` and `F13`, one native XX/YY interaction on
each opposite-pair register, followed by `F01` and `F23` across the two
registers. This graph matches the pair regrouping: the first two F gates can
carry the entangling part of the two Bell-basis changes, and the later
cross-register gates are where a label-exchange eigenbasis could be factored.
For example, `F(pi/4,0)` is locally equivalent to a CNOT-class entangler and
can supply the entangling part of a Bell transform.

This is a synthesis hypothesis, not a decomposition of the frozen word. The
word has three fan-in CNOTs and the retained `A1/B1` rotations before `F02` and
`F13`; those gates must be incorporated into the reverse factorization. The
actionable next construction is to factor the explicit target eigenbasis above
against the exact frozen chronology

`L0 CX01 A1 CX21 B1 CX31 L1 F02 L2 F13 L3 CX12 L4 F01 L5 F23 L6`,

absorbing only permitted adjacent local factors, and solve the remaining
two-register exchange transform in the `CX12/F01/F23` tail. Any resulting
member retains the family cost of four native CX plus four F, compiling to
12 CNOTs. A match would still require complete serialized gate lists,
independent numerical checks, and an exact certificate before promotion.

The construction is related to the previously documented orbit-Fourier
families, but uses the tensor factorization into opposite-pair registers and
the Bell-label exchange blocks as its basis-selection rule. It does not
replace their orbit encoders or prove a smaller gate count. It makes no
lower-bound, total-spin, or strong-Schur claim. If the exact word cannot
realize this target basis, that does not exclude other bases in the 98-family
or other 12-CNOT circuits.

Independent review by verifier98fresh confirms the register order, Bell-label
signs, `lambda/sigma_a` eigenvector coefficient, and multiplicities. The
review also agrees that matching this basis to the fixed word remains a
hypothesis until a factorization is actually derived.
