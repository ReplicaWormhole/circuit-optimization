# Round 17: dynamic spectrum-assignment objective

The 12-CNOT goal remains active and unmet. Three GPT-6 Luna collaborators
provided algebraic derivations (boards 85/91), static numerical reviews (83/90),
and independent verification (82/89). The coordinator owned runs 395/396 and
boards 84/87. Shared checkers, dependencies, baselines and the incumbent are unchanged.

The independently checked objective is min_D ||UV4-DU||F²/16, where diagonal D
contains 6(+1), 4(-1), 3(+i), 3(-i), allowing all eigenvalue orders and degenerate
bases. For unitary U this equals the original off-diagonal loss plus minimum
diagonal spectral mismatch, so it has the same exact zeros. Raw row-cost
Hungarian assignment includes both terms. Selected labels are held fixed for
each gradient. Distinct-root ties can create kinks; every selection is recorded.

Run 395 used exactly run 394 seed 9261620's audited best vector, with all 98
coordinates free. One L-BFGS-B fit was limited to 1,000 iterations and 20,000
objective-gradient calls, 110 seconds internally and 120 seconds externally,
on one thread, without perturbations or restarts. It stopped after one iteration
and three evaluations in 0.06673 seconds by relative reduction tolerance.
Best labeled loss: 0.1464466094; original loss: 0.1357233047; maximum off-diagonal
entry: 0.3535533925. No meaningful escape occurred. Independent verification
checked all five endpoint vectors and their ten complete gate lists.

Run 396 then tested eight fresh full98 initializations, seeds 9261700–9261707:
90 normal coordinates at scales 0.3, 1, 2, 3 in cyclic order and eight uniform
F angles. It used the same objective, one fit per start, at most 1,000 iterations
and 20,000 optimizer calls per start, and a total one-thread budget of 110 seconds
internally and 120 seconds externally. All eight completed in 42.767 seconds.
Best seed 9261703 returned to the same plateau: labeled loss 0.1464466094,
original loss 0.1357233047, maximum off-diagonal entry 0.3535533934.
The independent audit checked all sixteen endpoint lists and the global best,
matrices, 12-CX count, unitarity, hashes, RNG rules, both losses and history
multiplicities. None passes. Each retained endpoint also had one scalar original
loss diagnostic and gate-list checks outside the optimizer call counter; these
were documented before execution and remained within the same time cap.

Both runs are terminal failed. Inputs, code hashes, budgets and commands were
reserved before computation; independent preflight approvals preceded execution.
Run 395's finish-time ledger note recorded that its post-fit audit was pending;
the completed audit is retained without rewriting that historical row.
Board proposals 86/88 remain unexecuted duplicates whose dependencies prevented
claims. Recorded handoffs identify replacements 87/89. No computation or ledger
row belongs to either duplicate.

The algebraic report conditionally explains the loss split using
 a=(2+sqrt(2))/4 at eight attenuated diagonal entries: original=(1-a²)/2 and
labeled=1-a. Eight exact unit-magnitude diagonal entries would isolate an 8D
spectator block and an 8D active block. Exact block structure and its realization
by this gate word remain unproved. This is a synthesis clue, not an obstruction.

The best next hypothesis tests the repeated low-loss point for negative
curvature, with a newly frozen diagnostic budget and an independently checked
direction. Optimizer termination proves neither local minimality nor family
exclusion. Run 389's earlier Hessian sampled a different, higher-loss point and
objective.

Evidence: [run 395](../experiments/runs/395/CONCLUSION.md),
[run 396](../experiments/runs/396/CONCLUSION.md),
[objective derivation](work/algebraic98label/REPORT.md),
[plateau derivation](work/algebraic98label8/REPORT.md),
[static review](work/numerical98label8/REPORT.md), and
[independent audit](work/verifier98label8/REPORT.md).

Verification: 16 focused unit tests passed, with existing SQLite resource warnings.
Ledger/board audits and py_compile passed. Historical rows and accepted 13-CX /
exact 14-CX artifacts are unchanged. Adversarial review found no unresolved
substantive defect; accounting and deadline limitations are documented. Research
outputs are retained evidence, Python caches ignored. No remote state changed.
