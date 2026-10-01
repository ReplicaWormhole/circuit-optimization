# Numerical residual correction — run 352

Command: `python3 search13_fresh_chain_residual_generator.py`.

For run339 source U and D=U V U†, nearest-root output masks Q have multiplicities(6,3,4,3). Spectral projectors P are computed as fourth-root polynomials in D. The unitary polar factor W of sum QP intertwines D with the prescribed diagonal matrix to about1e-14; the minimum overlap singular value is0.99949835. This is a numerical correction, not a native gate circuit.

The full generator H=i log(W) is expanded in all256 Pauli strings, qubit0 the leftmost factor. Dominant terms are YXZI(-0.0135109), ZYZI(-0.0123368), YXIZ(-0.00730890) and YXZZ(+0.00730890). Squared coefficient masses by weights1–4 are about8.70e-7,6.47e-5,4.20e-4,7.02e-5. This motivates investigating distributed interactions rather than assuming a missing one-qubit correction. It is not a proof that another correction requires these supports; spectral-sector gauge freedom remains.

Independent check reconstructed H from the saved256 real coefficients, formed exp(-iH), and obtained intertwining error1.02e-14. Reconstruction, Hermiticity and unitarity checks pass. Adversarial review: no confirmed defect. No CNOT cost is assigned to W, and no upper bound, lower bound, novelty or optimality follows. Stored source hash, coefficients and numerical errors are kept artifacts.
