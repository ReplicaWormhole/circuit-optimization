# Coordinated two-slot topology fits from run 332

Run 334 stopped before synthesis: its duplicate guard detected that a
planned proposal recreated a previous fitted directed schedule. The source
loss roundtrip passed. No candidate or topology exclusion resulted.

Corrected run 335 fitted 20 deterministic two-slot mutations from
`search13_joint_prefix_slot8_refine_trial1.json`. The slot-pair families
were ten prefix/tail pairs, seven adjacent tail pairs, and three separated
tail pairs. Replacements prefer the complementary disjoint wire pair, with
mixed directions. A deterministic fallback uses other supports when needed
to avoid duplicates. Proposal selection did not use numerical raw loss.
Every fit varied all 168 local SU(2) parameters, from the same source warm
start, for at most 300 L-BFGS-B iterations.

The source optimizer loss `0.06019279505817157` reconstructed as
`0.06019279505817168`, an error of `1.1102230246251565e-16`. All 20 saved
13-CNOT candidates independently failed `check_circuit.py`. Case 15 had
the best maximum cycle off-diagonal error `0.2549137120576426`, with loss
`0.06019279508950225`; it returned to the old basin and did not improve it.
Case 9 had the second lowest fitted loss, `0.10070215086680244`.

The exclusion scope was exactly the contemporaneous root-level files
matching `search13*.json` whose JSON objects had `n=4` and a `gates` field:
171 directed schedules and 158 unordered interaction schedules. This is
not coverage of every historical topology. Later files can change this
proposal-generation set. Replay should use the complete frozen schedules
in `search13_joint_prefix_two_slot_v2_frozen_proposals.json`; each also
reconstructs from the result's `source_topology` plus its row's `changes`.

These are bounded numerical failures, not exclusions of the tested
topologies or of circuits below 14 CNOTs. The next small test uses large
seeded local-angle perturbations on the two lowest final-loss topologies,
to test the common warm-start limitation.

Recorded command:

```bash
python3 search13_joint_prefix_two_slot_v2.py --maxiter 300
python3 check_circuit.py search13_joint_prefix_two_slot_v2_case15.json
python3 experiment_log.py show 335
python3 experiment_log.py audit
```

Preserve the saved result and frozen proposals when reproducing in a
separate workspace copy.

Editorial correction after the run: the v2 script's opening docstring now
mentions the support fallback instead of claiming every replacement was
disjoint. Executable code is unchanged; the ledger retains the original
immutable code snapshot. Root independently reconstructed every fitted
loss and checked all twenty candidate schedules and validity results.
