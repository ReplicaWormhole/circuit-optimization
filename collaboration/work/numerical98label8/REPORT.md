# Static review of fresh8 full98 eigenlabel batch

Board hypothesis 90. This is a read-only review of reserved run 396. No
objective, matrix, optimizer, derivative, or certificate computation was run.

Frozen files are `experiments/runs/396/batch.py` and `config.json`. Their
SHA-256 hashes match the plan: `c0e5eeed37086b68f588d291c9fdfed69a71043e8fb52646c264bbc7b24c0c47`
and `775abe1ab572124817e907db842df1c938a5f4dbb0f2998c0c1aa34fb3a9194c`.
Static AST parsing passed.

## Checks

- The eight seeds are 9261700 through 9261707. Each draws 90 local/interior
  coordinates with the cyclic sigmas `[0.3, 1, 2, 3]`, then eight F angles
  uniformly from `[-pi/2, pi/2]`. This matches 84 local + 6 A1/B1 + 8 F
  coordinates, with all 98 free.
- One fit runs per seed, with no restart loop. The pre-evaluation guards stop
  at the shared 110-second deadline or before evaluation 20,001. The external
  120-second cap leaves ten seconds for endpoint serialization and the final
  checkpoint.
- The assignment routine and exact repeated spectrum come from frozen run395
  code with the SciPy version pinned. Each fit recomputes the assignment,
  holds the selected labels fixed for that gradient, and records the selected
  columns for every evaluation. The raw row-cost sum is already the complete
  labeled objective; no separate off-diagonal term is added.
- The best point is updated on each successful evaluation, including line
  search points. When the internal budget interrupts an active fit, the best
  evaluated point is serialized before that row is checkpointed. Full native
  and compiled endpoint lists are saved.

## Actionable accounting notes

1. If the total 20,000-call cap is intended to include every objective
   calculation, the post-fit `family.objective(...)` call used to store
   `original_loss` is one additional uncounted scalar objective evaluation
   per completed seed. If `maxfun` means optimizer calls only, record this
   diagnostic separately in the result accounting.
2. If the shared deadline expires between starts, the next start currently
   enters `minimize`, immediately raises `BudgetStop` before any evaluation,
   and records `stop='no_evaluated_point'`. Previously completed rows remain
   checkpointed, but the final stop reason is less informative than
   `overall_deadline_before_start`; a pre-start deadline check would clarify
   it.

Neither note drops an evaluated best point. The independent verifier's
prefit remains required before execution.
