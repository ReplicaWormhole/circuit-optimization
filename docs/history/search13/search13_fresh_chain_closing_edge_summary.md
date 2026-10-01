# Closing-edge replacements of the fresh chain13 lead

Ledger run 340 starts from `search13_fresh_chain_refine_trial1.json`, whose
13-CNOT schedule uses only the chain edges `01`, `12`, and `23`. At each
of all 13 slots, the experiment replaces that CNOT by `0->3` at even slots
or `3->0` at odd slots, adding closing-edge connectivity. Every fit uses
the same source local angles and freely varies all 168 local SU(2)
parameters for at most 300 L-BFGS-B iterations.

All 14 fixed local matrices are identities and all serialization chunks
are empty. No fixed chunks from the earlier exact-14-derived construction
are reused. The source layout is checked to contain exactly four `u3`
gates per local layer and 13 CNOTs. Its reconstructed loss
`0.0012158342275705673` differs from the saved source loss by
`3.0357660829594124e-18`.

All 13 saved candidates independently failed `check_circuit.py`, and their
saved-gate NumPy Frobenius losses reproduce optimizer records within
`1e-11`. The best mutation is slot 6, with loss `0.08193722893084353` and
maximum cycle off-diagonal entry `0.3483137071611512`. Both are worse than
the source's loss `0.00121583422757` and maximum entry `0.04303458`.
No further refinement was launched because no better basin appeared.

This is a failure of one warm start on each of 13 specified schedules. It
does not exclude those topologies, different initializations or a general
13-CNOT diagonalizer. The validated incumbent remains the exact 14-CNOT
circuit.

Recorded command:

```bash
python3 search13_fresh_chain_closing_edge.py --maxiter 300
python3 check_circuit.py search13_fresh_chain_closing_edge_slot6.json
python3 experiment_log.py show 340
python3 experiment_log.py audit
```

Preserve the existing result when reproducing in a separate workspace copy.
