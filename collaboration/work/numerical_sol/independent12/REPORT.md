# Independent exact certificate: PASS

Board128, scientific run407. Root executed the frozen checker exactly once:
`timeout 90s python3 experiments/runs/407/certify.py`.
Runtime was0.198070682 seconds under the90-second limit, with seed0,
no randomness and all four thread environment variables set to1.

The independently implemented checker uses only exact Fraction coefficients
in Q[z]/(z^8+1), z=exp(i*pi/8). This polynomial is the cyclotomic polynomial
Phi16, so the eight coefficient representation is exact at the stated
embedding. Multiplication reduces z^8=-1; complex conjugation maps z to z^-1.
No SymPy, NumPy, shared exact checker or shared circuit-construction helper
is imported. Primitive gate matrices, chronological composition and the MSB
right-shift permutation are independently constructed.

For the complete39-gate literal candidate, run407 establishes:

- exactly12 chronological CX gates;
- every one of256 row-eigenvector equations UV4=DU holds exactly;
- every one of256 unitarity equations UUdag=I holds exactly;
- root labels in order(1,i,-1,-i) are
  `[0,0,0,2,0,1,3,0,2,1,0,3,2,3,1,2]`;
- root multiplicities are(6,3,4,3).

All row-eigenvector and unitarity failure lists are empty. The checker
discovers the labels independently rather than receiving them as frozen
expected labels. Its discovered labels agree with the separate run406
certificate. The direct certificate means U V4=D U for the independently
constructed right shift |x0x1x2x3> -> |x3x0x1x2>; qubit0 is most significant.

## Saved-input and output audit

Read-only hashing after execution confirmed equality between this preparation
packet, the run407 inputs, and every corresponding result hash:

| Artifact | SHA256 |
| --- | --- |
| candidate.json | 58a862a4acef4e768c4613e9b142504592f5237a68539941148653a7f083eea2 |
| certify.py | dcd3956ba7194ffc869dc03cee1f92db46ac249c0f57029398e1ad0a41094e6f |
| config.json | 044ddc496be587abd767b1f8a14b3914f08441f1bb13410363f14509601334d2 |
| run407 unitary_exact.json | 2db88f4e63708ac41b0f22468216bd9bacedec2600b0010fe61cafbc5641dbff |

The saved matrix has16 rows and16 columns, with each entry represented by
eight rational coefficients. Execution marker, Python3.13.9 runtime versions,
Python/Fraction/math-extension hashes and single-thread environment are
recorded in run407. The current checker source/config exactly match their
frozen values. No new matrix execution occurred during this output audit.

## Adversarial review

No confirmed issue found in the independent certificate or saved-artifact
audit. Source review covered polynomial reduction, exact conjugation, all
rotation and U3 phase conventions, CX row permutation, chronological order,
MSB bit labels, independent right-shift construction, unique row-root
matching, all row inner products, complete config-key auditing, hash/runtime
guards, active ledger/board guards and exclusive execution marker.
The verifier and coordinator also reviewed the independent source before
execution. Syntax check `python3 -m py_compile .../independent12/certify.py`
passed. Board audit checks hashes/JSON only and is not a mathematical proof.

This certifies an exact ancilla-free12CX diagonalizer in the stated gate set.
It supplies an upper bound, not optimality, a new lower bound, total-spin
diagonalization, a strong Schur transform, or exclusion of other families.
Coordinator acceptance and ledger finish are separate provenance steps.

Prepared source, input packet, this report and run407 outputs are retained
research artifacts for the coordinator's focused local commit. Python cache
is ignored. Shared checkers and old numerical98phase12screen files were not
modified; no push or other remote state change was performed.
