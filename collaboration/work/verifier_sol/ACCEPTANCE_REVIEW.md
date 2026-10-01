# Independent saved-list and certificate audit of run406

Board127; reviewer verifier_sol. This review is read-only source/log/hash
inspection plus by-hand algebra. No additional circuit matrix, objective,
search, optimizer or exact-check invocation was performed. Root retains
incumbent/acceptance and commit ownership. The separate Fraction/Phi16
independent certifier is pending; this report does not claim it ran.

## Findings

No blocking defect found. The complete saved12CX candidate and run406's exact
certificate agree with the independently derived four-sector diagonal labels.
The result proves existence of this ancilla-free12CX diagonalizer within the
specified gate model, subject to the ordinary trust in the audited exact
checker/runtime. It establishes neither optimality nor a global lower bound.

## Frozen provenance and complete gate word

All five run406 inputs (check_merged.py,candidate.json,config.json,MANIFEST.json,
PLAN.md) are byte-identical to the authoritative candidate12 pack. Script,
config and candidate byte hashes match result.json and the reservation:

- script: bdccafc08432d507a766e94c509e463b09294b5666577742b9c95747ca181533
- config: 98245eac648421862b5d3ee34dd54d1f826b8ebca18af39db61875de3bfa6b15
- candidate: 58a862a4acef4e768c4613e9b142504592f5237a68539941148653a7f083eea2

Checker source hashes equal config. The exclusive execution_started marker
binds run406 to the same script/candidate. Ledger406 is complete and retains
the same frozen_config/config-byte-hash/code binding. Its canonical candidate
hash2c4f47648585599b52d2324c11ae6442f342ad43ac782581a4af925ef96fab71
differs from the saved candidate byte hash only by JSON formatting: the entire
parsed candidate, including metadata and gate word, is equal.

AST-selected literal builders reconstruct the full saved candidate exactly
without constructing gate matrices. There are39 chronological gates:
12CX,6H,4Rz,8Ry,5U3,2X,2Rx. All qubits lie in0..3; no ancilla or native
uncompiled gate is used. The complete chronological CNOT subsequence is:

CX02,CX13,CX01,CX23,CX02,CX32,CX02,CX32,CX10,CX32,CX12,CX32.

The gate blocks, using zero-based inclusive positions, are Bell0..3,
router4..5, mergedBP6..18, exactCH19..23, finalCCH24..38. Their CNOT counts
are2+2+4+1+3=12. Every elementary two-qubit gate is a CNOT, so the cost is
an actual gate-list count rather than an unknown compilation estimate.

The N3 reflection is the requested target conjugation
Rz2(-pi),CX02,Rz2(pi), exactly controlled(-X), with no dropped active minus
sign. N1,N2,N4 target conjugations and the output T0 agree with the four
block products in MERGE4_STATIC.md. The final wrapper chronology implements
Vprime=Ry(-3pi/4)Rx(-pi/2), with negative q1 control and positive q3 control.
It maps activeY to H and inactiveZ to Y confined to d00. This agrees with
REPAIRED13_STATIC.md, using the firstP control order(x=q0,v=q3) implied by
the merged block. Qubit0 is most significant and gates are chronological.

## Numerical result and exact output

The numerical checker exited0 in0.07517625301261432seconds, without timeout.
It reports12CNOTs, two-qubit depth9, valid_diagonalizer=true, tolerance1e-9,
unitarity error1.5543122344752192e-15, maximum off-diagonal entry
5.673584001061413e-16, and root error7.105427357601002e-15. These are
numerical evidence; they are not the exact certificate.

The exact checker exited0 in0.38095643199631013seconds without timeout. It
reports exact_diagonalizer=true,12CNOTs and cyclotomic conductor16. Its16
output labels are exactly the independently predicted list:

[0,0,0,2,0,1,3,0,2,1,0,3,2,3,1,2]

in root order(1,i,-1,-i), with multiplicities(6,3,4,3). Equivalently the
explicit diagonal is
[1,1,1,-1,1,i,-i,1,-1,i,1,-i,-1,-i,i,-1], with counts(6,4,3,3) in order
(1,-1,i,-i). Wrapper result.json agrees with both checker logs and the
frozen expected-label guard. Wrapper/numerical/exact stderr files are empty.
Both invocations stayed below their30/60second limits. The saved command
uses one-thread environment and external100second cap; no optimization or
random draws occurred.

## Exact proof coverage from source inspection

The frozen exact_check.py builds U exactly over Q(zeta_16), with rational-pi
angles parsed as Fractions. Its H/Rx/Ry/Rz/U3 primitives match the stated
gate conventions. Phase-only U3(0,0,lambda) is exactly diag(1,exp(i lambda)).
All listed primitive gates are analytically unitary, and CNOTs are
permutations, so the exact product U is unitary without a floating-point
assumption. The checker does not separately test UdagU, but its supported
primitive definitions establish it for this word.

Chronological local updates multiply rows on the left; CNOT row permutations
use masks1<<(3-qubit), matching MSB-zero convention. For every output row
and all16 input columns the checker tests

row[right_shift(basis)] = lambda * row[basis],

where right_shift(basis)=((basis&1)<<3)|(basis>>1). This is precisely the
full row relation U V4=D U. A row must match exactly one root; zero rows
would match multiple roots and fail. Thus no columns, degenerate sectors,
inactive selector branches or output rows are omitted. Unitarity gives the
equivalent U V4 Udag=D statement. Root labels are allowed in arbitrary
order; the fixed label check is an extra consistency check for this word.

The exact matrix construction shares only the angle-regex parser with the
floating checker, while the field arithmetic and matrix updates are separate.
Independent Fraction/Phi16 arithmetic would further reduce implementation
trust and is being prepared separately; it is not included as completed
verification here. Neither checker certifies total-spin diagonalization or
a strong Schur transform, and this review makes neither claim.

## Adversarial closeout and cleanup

Review followed ~/.codex/adversarial_review.md. It checked full word/control
chronology, active/inactive phases, N3 sign, T0 sector phase, exact field
coverage, eigenvalue order, all hashes versus semantic canonical candidate
identity, timeout/result handling and claim strength. No confirmed issue
remains in this audit. Pending work: separate independent Fraction certificate
and coordinator-owned acceptance/commit. Own report and preparation directories
are intentional evidence; existing concurrent dirty/untracked work was preserved.
