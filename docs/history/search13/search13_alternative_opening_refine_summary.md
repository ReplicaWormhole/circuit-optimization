# Alternative opening refinement — run 348

Command: `python3 search13_alternative_opening_refine.py --seed 35000 --maxiter 650`.

Source: run 347 proposal 0, a 13-CNOT pure local-layer circuit with first interaction 02. Four starts use perturbation scales 0, 0.05, 0.2, 0.6 and seeds 35000–35003. All 168 local parameters are free.

All four serialized candidates fail the independent checker. The first three converge to ordinary off-diagonal Frobenius squared loss divided by 16 of 0.001455266622; the fourth reaches a worse basin. Best maximum off-diagonal entry is 0.03519943712. This is smaller than the original chain source maximum entry but has worse ordinary loss. Neither is a valid diagonalizer.

Root reconstructed all four matrices with the independent NumPy builder, checked saved objective values within 1e-11, maximum entries within 1e-12, and exactly 13 CNOTs. No exact upper bound, topology exclusion, or optimality claim follows.

Adversarial review: checked objective distinction, arbitrary local layers, native CNOT count, initialization, and numerical scope. No confirmed defect. Remaining risk: bounded local optimization can miss valid circuits on these schedules.
