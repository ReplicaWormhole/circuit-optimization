# Exact realness, unitarity, and nonzero projective factors

This proves local-gate properties of `canonical_algebraic_ansatz.json`; it does not prove the four-qubit UV4=DU identity.

For the five exceptional squared sine values, use the exact inequalities `1 < sqrt(2) < 3/2`:

- `(252−45sqrt(2))/367` lies between `369/734` and `207/367`, hence strictly between0 and1.
- `(3−sqrt(2))/7` lies between `3/14` and `2/7`, hence strictly between0 and1.
- `5/12−5sqrt(2)/48` lies between `25/96` and `5/16`, hence strictly between0 and1.
- `(184+72sqrt(2))/367` lies between `256/367` and `292/367`, hence strictly between0 and1.
- `(19+3sqrt(2))/28` lies between `11/14` and `47/56`, hence strictly between0 and1.

The other non-rational-angle squared sine values are1/6,1/3,2/5,3/5,2/3,5/6, also strictly between0 and1. Thus every square root used for their signed sine/cosine is real, and `s²+c²=1` by construction. Every rational-angle sine/cosine is real and obeys the same identity. The exact program additionally verifies this polynomial identity for all60 rotations without floating-point zero tests.

For any exact real `s,c` with `s²+c²=1`, define the normalized local matrix `R=cos(theta/2)I−i sin(theta/2)P`, with the signed half-angle radicals in the ansatz. The half-angle identities give `cos²(theta/2)+sin²(theta/2)=1` and `2sin(theta/2)cos(theta/2)=s`, with signs selected from the prescribed angle branch. Hence `R†R=I`, because `P†=P`, `P²=I`.

When `c != −1`, `A=(1+c)I−i sP` satisfies `A†A=2(1+c)I`. Since real `c∈[−1,1]` and `c != −1`, its scalar norm is strictly positive, so `A` is a nonzero scalar multiple of the corresponding normalized unitary rotation. For irrational-angle gates the strict sine-square bounds above imply `c∈(−1,1)`. For rational-angle gates `c=−1` is recognized exactly and the program uses `P` as projective representative, which is itself unitary and a nonzero scalar multiple of the exact π rotation. Rational angles with `c=1` have `A=2I`.

Multiplication by the extra symbolic common denominator in the verifier also preserves nonzero projective equivalence: all displayed radical denominators have positive radicands and positive integer denominators. The circuit contains exactly13 chronological CNOTs; its local scalar factors do not change that count. If the separate exact verifier proves every entry of `A_total V4−D A_total` is zero, division by the nonzero accumulated scalar yields the same identity for the normalized unitary gate list. No numerical tolerance can establish this final global premise.
