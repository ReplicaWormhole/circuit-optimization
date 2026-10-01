# Fresh chain curvature escapes — run 350

Command: `python3 search13_fresh_chain_curvature_escape.py --maxiter 400`.
Source is run339 trial1; all local parameters remain free and every candidate has 13 native CNOTs.

Full automatic Hessian at ordinary loss 0.00121583422757 has minimum eigenvalue -8.09e-9 and gradient norm 3.89e-8. These floating-point figures do not certify positive semidefiniteness or a minimum. Hessian symmetry and extreme-mode finite differences pass. Three lowest positive modes above 1e-4, indices112–114, were probed by signed pi steps before local optimization. Four fits return to the source basin; two reach worse loss0.005955870961. All six fail independent matrix checks.

Root independently reconstructed all serialized gates with NumPy, verified saved losses within1e-11, maximum off-diagonal entries within1e-12 and exactly13 CNOTs. Adversarial review found no confirmed defect; remaining risk is incomplete basin coverage, coordinate-dependent curvature and numerical precision. No circuit or global exclusion results.
