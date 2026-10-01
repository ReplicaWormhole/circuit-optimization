# Static review of full98 dual-objective curvature diagnostic

Board hypothesis 92. This review inspected source/configuration only. No matrix,
objective, gradient, assignment, Hessian, optimizer, or certificate calculation
was run by this reviewer.

The reserved copy is `experiments/runs/397/diagnostic.py`, with
`experiments/runs/397/config.json` and `PLAN.md`. The code SHA256 is
`2e518feeaa4a57bf6a3497a17696643049bf8123750f81628606478ac9f0ccda`, and the
config SHA256 is
`f70c638de7ac8fecd9ffc23750ae4b5b1ceecd92236855ecc3f31052f9933d1a`; both
match the coordinator's frozen hashes.

The diagnostic binds the exact run396 best vector at seed 9261703 and uses the
same 98-coordinate family chart for the original off-diagonal objective and
the fixed-label residual. There is no fit or perturbation. It computes at most
two gradients and two Hessians. The fixed labels are taken from the raw-cost
Hungarian assignment and held constant in the branch derivative. Symmetrizing
the Hessian before diagonalization, saving eigenvectors as columns, and fixing
each vector's sign by its largest-magnitude component are implemented
consistently.

The 16 assignment-gap alternatives correctly ban the selected root value in
one selected row across *all* duplicate slots (`slots == labels[i]`), while
leaving the remaining root capacities unchanged. This rules out zero-cost
same-root column swaps; the raw row-residual costs already represent the full
labeled residual. The one-thread controls, 110-second internal deadline,
120-second external plan, and limited interpretation are recorded.

One execution-path issue was sent to root before the diagnostic ran:
`BudgetStop` can be raised in the assignment loop before the current
`try/except`, which starts only around the Hessians. If the internal deadline
expires during that loop, the process exits without writing partial assignment
state, base gate lists, or `result.json`. Root should put assignment and Hessian
work under a common timeout handler and serialize partial state, or otherwise
ensure the assignment section cannot raise uncaught `BudgetStop`. No source
changes were made by this reviewer. This is a pre-execution review finding,
not evidence about curvature or the circuit family.
