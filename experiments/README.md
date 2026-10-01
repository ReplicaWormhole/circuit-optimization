# Experiment workspaces

Each new reserved scientific run gets `runs/<ledger-id>/` for its plan, code,
outputs, checks, and conclusion. Create it with
`python3 experiment_log.py workspace ID --code PATH` after reserving the ledger row.
See the [run protocol](../docs/EXPERIMENTS.md).

The authoritative database is still `../experiments.sqlite3`. Its immutable
candidate and code snapshots remain in `../candidates/` and
`../code_snapshots/`; this directory is for readable run-specific work.
