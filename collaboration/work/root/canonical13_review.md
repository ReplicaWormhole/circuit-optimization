# Canonicalization review

Identity: for a control, U=Rz(alpha)Ry(beta)Rz(gamma); for a target, U=Rx(alpha)Ry(beta)Rx(gamma). The first factor commutes with the CNOT, so it is multiplied into the next local layer. Idle-wire gates commute and are moved entirely. Numerical Euler reconstruction is checked at every twoqubit layer. SU2 normalization of U3 gates introduces only global phase. Final Rz output phases commute with the computationally diagonal target and may be omitted.

Result has13 CNOTs and64 rotation angles, or60 after four finalZ removals. Root independently compared the full matrix upglobalphase(1.12e-15), reconstructed the stripped gate list from final-layer metadata, and checked13 CNOTs and numerical diagonalization. Qiskit in the existing venv reports2.27e-15 offdiagonal error. A default-python Qiskit import failed because it is not installed there; no installation occurred, and the venv run supplies actual evidence.

Adversarial review: no confirmed defect. Euler branches and degenerate charts can affect angle recognition; no exact angle or lower-bound claim is made. Generic angles persist. Files are kept research artifacts. README distinguishes numerical13 success from exact14 incumbent. Root independently rechecked all16 run363 saved matrices/counts/schedules/ranking and confirmed failed twelve-CNOT local fits do not exclude architectures.
