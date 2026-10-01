# Independent acceptance audit of the exact11CX construction

Board137; verifier11router. Read-only saved-word/source/log/hash review, with
by-hand phase derivation. No additional numerical or exact circuit matrix,
objective, optimizer or checker invocation was performed for this audit.
Root owns formal acceptance, documentation and commits.

## Findings and scope

No blocking defect found. Run411's exact cyclotomic certificate and independent
run412's Fraction/Phi32 certificate validate the same complete30gate,11CX,
ancilla-free four-qubit diagonalizer U V4=D U. All256 diagonalization relations
and all256 exact unitarity relations are covered by412; discovered row labels
agree with411 and the independent by-hand prediction. This supports acceptance
of C_min<=11 in the arbitrary-one-qubit/all-to-all-directed-CNOT model.
It proves no optimality, new lower bound, total-spin diagonalization, strong
Schur transform or membership in an older fixed-coordinate family.

## Full saved word and provenance

All six411 inputs (check_analytic.py,candidate.json,config.json,MANIFEST.json,
PLAN.md,SECTOR_REVIEW.md) are byte-identical to the frozen candidate11 pack.
The entire parsed candidate reconstructs from AST-selected literal builders
without matrices. Its30chronological gates are11CX,6H,4Rz,8Ry,1U3, all on
qubits0..3. Every elementary two-qubit gate is a CNOT; there is no omitted
native-gate compilation cost or ancilla.

The chronological CX subsequence is
CX02,CX13,CX01,CX23,CX02,CX32,CX02,CX32,CX10,CX12,CX32.
Zero-based blocks are Bell0..3,router4..5,modifiedmerge6..18,CH19..23,
finaltwo-selectors24..29, with counts2+2+4+1+2=11.

The modified final merged reflection is chronological
Ry2(7pi/8),CX32,Ry2(-7pi/8), implementing
N4prime=Ry(-pi/4)N4=-cos(pi/8)X+sin(pi/8)Z exactly. The T0 phase remains
present. The final block is H2,CX12,H2 then
Ry2(3pi/8),CX32,Ry2(-3pi/8), implementing CZ(u,y) then controlled
Nv=sin(pi/8)X+cos(pi/8)Z with positive u=q1,v=q3 controls.
No obsolete negative-control wrapper or target phase was silently retained.

411 hashes match result/config/ledger and frozen pack:

- script: dc3686e6ee071789c489c7b8f419c2fcbd65a6e5e97f7642b7445e232e1315ea
- config: d021c7a2c74289f7b2e48aa81ec730b5290fb624a71492c5b849233430a48c63
- candidate: 557c55cddefb697bcd4e038db300117d0aaaf6f6af97e851e0831baf257dbe4a

The current numerical/exact checker source hashes equal the frozen config.
Ledger411 is complete with its matching frozen_config/config-hash/code binding.
Its canonical candidate8ce7cca9a51201abfda0c3b6963cd5c5c717edb69f8134eb0acaa0faca3634d7
is semantically identical, including metadata and full word; bytehash difference
is JSON formatting only. No candidate input changed after either reservation.

## Numerical and first exact certificate

Run411 numerical checker exited0 in0.07517395302420482seconds, no timeout,
with actual11CX,depth8,valid=true,tolerance1e-9. Maximum offdiagonal entry is
5.039344865125116e-16, unitarity error1.7763568394002505e-15 and root error
8.43769498715119e-15. These are numerical evidence, not the exact certificate.

The exact checker exited0 in0.44670272001530975seconds without timeout and
reports conductor32,11CX,exact_diagonalizer=true. All16output labels match
[0,0,0,2,0,1,3,0,0,1,2,3,2,3,1,2] in roots(1,i,-1,-i), counts(6,3,4,3).
The explicit D is[1,1,1,-1,1,i,-i,1,1,i,-1,-i,-1,-i,i,-1].
All wrapper/checker stderr files are empty. Both calls remain below60/90second
budgets. This checker exactly tests each row against all16right-shift columns,
with chronological left multiplication and MSB-zero qubit masks. Its supported
rational-pi primitive gates are analytically unitary;412 additionally checks
unitarity explicitly in a separate arithmetic implementation.

The full by-hand phase audit in candidate11/SECTOR_REVIEW.md agrees: sector
outputs are CZ_xy,(-i)^yZ_x,i^xZ_y,f_i(y)Z_x. In d11 the scalar identity term
and relative phase in f_i=diag(i,1) are retained while its Pauli axis rotates
M->H->Z with positive sign. d00 remains CZ instead of receiving the12word's
Y permutation, explaining the allowed change in output eigenvalue order.

## Independent Fraction certificate and coverage

Run412 copied certify.py/config/candidate/PLAN/MANIFEST unchanged from
numerical11exact. Its candidate bytehash equals411 exactly. Saved hashes:

- independent script: f0c9c05f872139937e1f9fd5ad22006033ec1d3f7e3467788f4d1199b4b6eb09
- independent config: dd00ad8759bf193ed07fb92ce0b6c47674dae1f3c40d1563c2e43d4a5aadc7b5
- exact unitary artifact: fbb007935825221f250efc5d86fb0dd3db1e438f92706ef7eac0ae457a81e965

Result/candidate/script/config/unitary hashes and raw-config ledger binding
agree; ledger412 is complete. The saved exact matrix has16rows/16columns,
each entry containing16rational coefficients. No matrix entries were
recomputed in this review. Runtime is0.1763960880052764seconds, stderrempty;
all thread environment variables are1. The result reports exact_valid=true,
exact_unitarity=true, exact_row_eigenvector_identity=true,11CX,30gates,
empty row/unitarity failure lists, and the same16discovered roots/counts.

Source inspection confirms genuine separate arithmetic and parsing: imports
include standard Fraction/dataclass/re/math but no NumPy,SymPy or repository
checker. Coefficients are rational in Q[z]/(z^16+1), z=exp(i*pi/16).
Multiplication degrees0..30 reduce with z^16=-1; complex conjugation uses
z->z^-1 modulo32. Phi32(z)=z^16+1 is the irreducible cyclotomic polynomial,
so coefficient equality is exact field equality. Phase exponents16*f and
half-angle exponents8*f correctly represent all angle strings; i=z^8 and
H's1/sqrt2=(z^4+z^-4)/2. Rx/Ry/Rz/U3 matrices match the declared conventions.

Independent bit-list construction gives the right-shift permutation
[0,8,1,9,2,10,3,11,4,12,5,13,6,14,7,15]. For every row and column, the
certifier discovers its matching root from z^(8k) and checks U V4=D U.
It separately evaluates every row inner product with every conjugate row,
checking all256 entries of U Udag=I exactly. For square U this proves full
unitarity and the equivalent U V4 Udag=D. No inactive branches, columns or
degenerate sectors are omitted. Roots are discovered rather than assumed
from the frozen prediction; their agreement is an additional cross-check.

## Closeout

Adversarial review followed ~/.codex/adversarial_review.md: actual word/cost,
reflection/selector signs, scalar phases, MSB/chronology, exact ring scaling,
all-row/all-column coverage, independent unitarity, hashes and claim strength.
No unresolved defect found. Remaining step is coordinator-owned formal
EXACT11 acceptance and commit; this reviewer changed no incumbent. Existing
12/13/14 certificates and concurrent work remain preserved. Own reports and
packs are intentional evidence. Optimality and smaller CNOT counts remain open.
