# Independent analytical11CX certificate: PASS

Board138, scientific run412. Root executed the frozen standalone checker
once under the reserved90-second budget:
`timeout 90s python3 experiments/runs/412/certify.py`.
It completed in0.176396088 seconds with one thread, seed0 and no RNG.
The source is the authoritative30-gate11CX candidate also checked in411.

Independent exact Fraction arithmetic in Q[z]/(z^16+1), z=exp(i*pi/16),
establishes all256 entries of U V4=D U and all256 entries of U Udag=I.
Both residual/failure lists are empty, and the literal gate count is11CX.
The checker independently discovers row-root labels in order(1,i,-1,-i):

`[0,0,0,2,0,1,3,0,0,1,2,3,2,3,1,2]`.

Multiplicities are(6,3,4,3). The labels agree with separate run411 and the
by-hand sector prediction; they were not supplied as checker inputs.
Thus this is an exact ancillary-free four-qubit right-shift diagonalizer
in the all-to-all CX plus arbitrary one-qubit gate model. It gives an11CX
upper bound, not optimality, a lower bound, total-spin diagonalization,
a strong Schur transform, or exclusion of other circuits.

## Independent implementation and source audit

The verifier is deliberately adapted from independent run407, with every
conductor-sensitive quantity reviewed:16rational coefficients, reduction
z^16=-1, exponent period32, conjugationz->z^-1, i=z^8,
1/sqrt2=(z^4+z^-4)/2, phase exponent16f and rotationhalfangle exponent8f.
The polynomial is Phi32, giving a unique exact coefficient representation
at the stated embedding. Unsupported nonintegral exponents are rejected.
Primitive rotations/U3 retain all phases, chronological multiplication uses
independent row updates, CX uses independent MSB row permutations and the
right shift is constructed directly from the bit labels.

Only standard-library modules are imported; no NumPy, SymPy, shared circuit
checker or shared unitary-construction code is used. The source audits cover
all256 row-root comparisons, all256 row inner products, actual CX count,
complete config-key/value checking, immutable inputs, runtime hashes, active
ledger/board ownership and one-shot guard. Syntax compilation passed.

## Saved artifact audit

After execution, read-only hashing verifies equality of prep/run/result
script, config and candidate hashes; candidate bytes also equal411's input.
The saved exact matrix hash equals the recorded result. Its shape is16x16,
every entry has16rational coefficients and all stored denominators are positive.
No matrix or coefficient recomputation was performed during this audit.

| Artifact | SHA256 |
| --- | --- |
| candidate.json | 557c55cddefb697bcd4e038db300117d0aaaf6f6af97e851e0831baf257dbe4a |
| certify.py | f0c9c05f872139937e1f9fd5ad22006033ec1d3f7e3467788f4d1199b4b6eb09 |
| config.json | dd00ad8759bf193ed07fb92ce0b6c47674dae1f3c40d1563c2e43d4a5aadc7b5 |
| run412 unitary_exact.json | fbb007935825221f250efc5d86fb0dd3db1e438f92706ef7eac0ae457a81e965 |

No confirmed issue was found in adversarial source/output review. Remaining
scope limitation is optimality; coordinator acceptance/ledger closure and
the independent review of the11CX candidate remain separate provenance steps.
The unexecuted Newton comparator is preserved and its board136 closed without
a fictitious run. This packet and run412 outputs are intentional research
evidence; Python caches are ignored. Root owns the focused local commit and
incumbent promotion. Shared checker code was not edited; no remote action
or extra numerical/exact run was performed by this collaborator.
