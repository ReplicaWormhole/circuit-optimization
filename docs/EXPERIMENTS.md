# Scientific run workspaces

Reserve every numerical or symbolic search in `experiments.sqlite3` before
execution. The database remains the authoritative attempt index. A run
workspace groups the human-readable plan, source script, outputs, and review
notes; it does not replace the ledger's code and candidate snapshots.

## New run

1. Claim a bounded hypothesis on the [collaboration board](../collaboration/README.md)
   when working in a coordinated round.
2. Write the search script and freeze the seed, complete topology rule,
   optimizer, iteration or time limit, and compute budget.
3. Reserve the run with `python3 experiment_log.py reserve ... --code PATH`.
   Use the ID printed by that command.
4. Run `python3 experiment_log.py workspace ID --code PATH` to create
   `experiments/runs/<ID>/` and copy the reserved script after checking its
   hash. Fill in `PLAN.md`, place run-specific inputs and outputs there, and
   link the board hypothesis. Omit `--code` only when no script was reserved.
5. Execute the bounded search with one compute thread. Save the complete
   [gate list](CANDIDATES.md) and independent numerical checks. Finish the
   ledger row even when the run fails or is inconclusive, then write the
   workspace conclusion.

The `workspace` command is safe to repeat: it never overwrites an existing
plan or a changed script. Keep generated candidates in the run directory as
readable evidence;
`finish --candidate` separately archives a content-addressed copy under
`candidates/`. Use `python3 experiment_log.py audit` to check ledger snapshots.

Historical scripts and outputs are in `archive/legacy_root/`. Their flat
filenames and companion paths are preserved; see the [archive guide](../archive/README.md).
Do not rewrite old ledger rows or archived snapshots to fit the new layout.
