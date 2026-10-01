# Hermitian spectral reweighting from the slot8 basin

Run333 uses H_t = (V+V†)/2 + t(V-V†)/(2i) at t = 1/4, 4, -1/4, -4. Its eigenvalues on the four V sectors are 1,t,-1,-t, all distinct. Thus an exact diagonalizer of H_t is an exact diagonalizer of V and conversely, with full freedom within each degenerate sector. Explicitly V = a H_t + b H_t^3, where b=(i-t)/(t^3-t) and a=(t^3-i)/(t^3-t). The denominator is nonzero for the chosen t. This is an operator identity on all four sectors; it does not assert that a nonzero H_t residual gives a good V residual.

Four bounded fits reweighted the same run332 13-CNOT topology and source local angles, up to 200 iterations on normalized H_t offdiagonal Frobenius loss and then 300 on the original cycle loss. All 168 local parameters varied. All four weighted candidates and all four polished candidates independently fail the original cycle check. All polished candidates return to the same loss about 0.060192795058 and maximum offdiagonal error about 0.25491354. Changing sector weights alone did not escape this basin in these trials.

Checks: source topology and loss roundtrip; H_t Hermiticity and distinct eigenvalues; seeded directional finite-difference gradient check for every t; both stage losses reconstructed from serialized gates with the independent NumPy matrix builder (agreement <1e-11); CNOT counts and validity reevaluated. The general polynomial identity follows by substituting 1,t,-1,-t; exact SymPy interpolation was also checked for each chosen rational t. Adversarial review found no confirmed defect. All fit outcomes are numerical evidence, not exclusions or exact circuit certificates.

Reproduce in a separate copy without the existing result file:

```bash
python3 search13_hermitian_reweight.py --maxiter 200 --direct-maxiter 300
python3 experiment_log.py show 333
```

Scripts, candidates, results, and this summary are kept research artifacts.
