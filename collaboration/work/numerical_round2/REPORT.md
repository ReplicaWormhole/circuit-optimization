# Changed twelve-CNOT architecture — hypothesis19, run372

Researcher numerical_round2 proposed the frozen plan, warm-start conversion,
optimizer/serialization interfaces and independent NumPy validation recipe.
Coordinator implemented and executed run.py because the child command launcher
failed. This was one shared run. Findings were exchanged directly among all
three researchers and persisted on the board.

Delete original slot6 (01), rewire original7 (32) to02 and8 (12) to03.
Keep fourteen arbitrary local layers and168 parameters. The native ordered list
is 01,32,12,10,23,21,02,03,10,23,21,01. Its minimum unordered pair-order
Hamming distance from the13 single-deletion canonical schedules is2; novelty
is not exhaustive over all historical configurations or circuit equivalences.
Source, code, dependency hashes and complete budget are in plan.json.
Ledger reservation372 precedes optimization and archives code bytes.

Warm preparation merges the exact source's approximate local matrices into
fourteen U3 layers. The unmutated reconstructed source passes shared checking
and axis-objective loss<1e-20; this validates a numerical warm start, not exact
radical serialization. No dependencies were added.

| Seed | Perturbation std | Iterations | Best ordinary loss | Maximum offdiagonal |
|---|---|---|---|---|
|926201|0|200 (iteration limit)|0.3352359084631581|0.7561025755556912|
|926202|0.15|110 (optimizer convergence)|0.3125000000000019|0.999999999999997|

Both native12 gate lists fail diagonalization. Separate NumPy Kronecker gates
and MSB bit permutations reconstruct each saved gate list and agree with the
shared checker. Unitarity errors<9e-16. Best by ordinary loss is case1; case0
has smaller maximum residual. Both are retained; case1 is content-addressed by
the ledger. Optimizer convergence is not circuit validity.

Budget: two starts max200 iterations each,60 seconds each,120 seconds total
compute, one thread. Actual elapsed8.052 seconds; no timeout. Run and hypothesis
finished inconclusive. Failure excludes neither the topology nor twelve-CNOT
solutions. No passing candidate was sent for exact certification.

Recommendation: derive an analytically decorated central three-wire rewrite
before selecting a different bounded topology. Do not deepen this unchanged
schedule without a new mechanism. No follow-up dispatched.
