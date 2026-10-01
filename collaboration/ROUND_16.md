# Round 16: broad full98 initialization batch

The continuing goal is a valid 12-CNOT diagonalizer of V4 with independent
gate-list validation and an exact certificate. That goal remains open.
The three collaborators used GPT-6 Luna: algebraic structure (board81),
numerical strategy review (board79), and independent verification (board78).
The coordinator owned batch board80 and scientific run394.

Run394 freed all 98 coordinates: 84 local rotation-vector coordinates, six
interior target-wire coordinates and eight XX/YY angles. The fixed chronological
word is L0 CX01 A1 CX21 B1 CX31 L1 F02 L2 F13 L3 CX12 L4 F01 L5 F23 L6.
Four native CX and four F blocks compile to 12 CX. Seeds9261600..9261631
used local scales cycling through 0.3,1,2,3 and independent uniform F angles.
Each start had one L-BFGS-B fit, at most800 iterations and20000 evaluations,
within a single-threaded285-second internal/300-second external batch cap.
The board claim, ledger reservation, workspace, hashes and independent preflight
preceded optimization. No unrecorded restart occurred.

All32 starts completed in117.314 seconds. The best seed9261616 has loss0.125
and maximum off-diagonal entry0.7071067812. Independent NumPy/SciPy evaluation
checked all64 endpoint lists and both best lists, coordinate agreement,
compilation, counts, unitarity, RNG rule and hashes. No circuit passes and no
near-zero endpoint occurred. Run394 is terminal failed; board80 inconclusive.
This is bounded numerical evidence, not family exclusion or a lower bound.

The algebraic collaborator derived the opposite-pair factorization
V4=(S_A tensor I_B) SWAP_AB, A=(q0,q2), B=(q1,q3), and the Bell-label
eigenblocks with multiplicities(6,4,3,3). The verifier independently checked
the signs and eigenvectors. Matching the exact eigenbasis to the full98 word
remains an unproved synthesis hypothesis.

The next proposed diagnostic minimizes the eigenvalue-assignment residual
min_D ||UV4-DU||F^2/16, allowing every eigenvalue order with the required
multiplicities. It would change the optimization path without fixing a basis
inside degenerate eigenspaces. This proposal has not been executed or reserved.

Evidence: [run394](../experiments/runs/394/CONCLUSION.md),
[algebra](work/algebraic98continue/REPORT.md),
[numerical review](work/numerical98continue/REPORT.md), and
[independent audit](work/verifier98continue/REPORT.md).
The accepted13-CNOT incumbent and exact14 baseline remain unchanged.

Closeout checks: 16 focused unit tests passed; ledger and board audits passed.
Adversarial review found no unresolved substantive defect. Remaining risk is
mathematical: neither numerical endpoint validation nor these failures certify
exactness or exclude the family. All run artifacts are retained evidence;
Python caches remain ignored. No dependencies or remote state were changed.
