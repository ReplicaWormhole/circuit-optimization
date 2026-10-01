# One-step topology escape from the changed-prefix near-hit

Ledger run 328 starts from `search13_joint_prefix_orbit_refine_0.json`, the
best refined 13-CNOT candidate of run 327. Its six-CNOT exact prefix had
already been changed by a deletion and another reroute, and a complete
suffix pair block had also been rerouted in run 325. Its label-free cycle
off-diagonal Frobenius loss is `0.0787330125161619` and it is invalid.

At each of all 13 CNOT slots, the script substitutes each of the other 11
directed CNOT edges and scores the resulting circuit without changing its
168 warm-start local angles: 143 raw tests in total. The six lowest-loss
changes on distinct slots that also change the unordered interacting pair
receive full 250-iteration L-BFGS-B local-layer fits. Pure direction
reversals are screened but not chosen for fitting because they preserve the
interaction graph.

All six saved candidates have 13 CNOTs and fail independent
`check_circuit.py` validation. The best maximum cycle off-diagonal error is
`0.38641938214453964` for slot 5 changed from `0->2` to `1->2`, in
`search13_joint_prefix_one_step_escape_rank1.json`. Its Frobenius loss
`0.07879651107629726` is slightly above the source's. The mutation with
the lowest **raw** loss, at slot 12, fitted to loss `0.31830084371767436`
and maximum entry `0.9600890925539337`. Raw warm-start score was therefore
a poor predictor of fitted quality in this small sample.

This is a local one-step search from a specific changed-prefix basin. It
differs from the earlier adaptive scan seeded by a simple prefix deletion,
which sampled random mutations and used a different near-hit. The numerical
failures do not exclude other mutations, initializations, or 13-CNOT
circuits. A next test could fit a deliberately diverse set of proposals
across **all** slots instead of selecting by raw loss alone, or jointly
change two interacting slots to cross the basin boundary.

Reproduce:

```bash
python3 search13_joint_prefix_one_step_escape.py --fits 6 --maxiter 250
python3 check_circuit.py search13_joint_prefix_one_step_escape_rank1.json
python3 experiment_log.py show 328
python3 experiment_log.py audit
```
