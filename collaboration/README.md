# Coordinated circuit research

The coordinator connects three researchers: algebraic, numerical, and native
gate sets. Each works independently on a bounded hypothesis, exchanges useful
findings directly, and leaves a durable handoff here. This board complements
the existing experiment ledger; it does not replace its duplicate-run guard,
code snapshots, or candidate checks.

## Interface

All commands run from the repository root. SQLite transactions prevent two
agents claiming the same hypothesis. Dependencies must be completed before
dependent work can be claimed. Failed and inconclusive experiments get new
follow-up hypotheses rather than being silently retried.
For a new scientific run, reserve its ledger ID, then create
`experiments/runs/<ID>/` with `python3 experiment_log.py workspace ID --code PATH`.
The [run protocol](../docs/EXPERIMENTS.md) explains the workspace contents.

```bash
python3 collaboration_board.py propose --agent root --title 'New architecture' \
  --gate-set cx --proposal 'State the distinct hypothesis and prior exclusions' \
  --budget '2 starts, 80 iterations each, one CPU thread'
python3 collaboration_board.py claim 1 --agent numerical
# Reserve a scientific run with experiment_log.py; use its actual returned ID.
python3 collaboration_board.py link-run 1 --agent numerical --run 355
python3 collaboration_board.py handoff 1 --agent numerical --to algebraic \
  --message 'Candidate path, residuals, exactness question and next action'
python3 collaboration_board.py submit 1 --agent numerical --candidate candidate.json \
  --costs '{"native_two_qubit":13,"compiled_cnot":13}' \
  --evidence numerical --notes 'Bounded numerical result; not exact'
python3 collaboration_board.py close 1 --agent numerical --status inconclusive \
  --outcome 'No passing candidate; state scope and distinct next hypothesis'
python3 collaboration_board.py list
python3 collaboration_board.py events
python3 collaboration_board.py submissions
python3 collaboration_board.py audit
```

IDs above illustrate syntax; never assume a run or hypothesis ID. A handoff
records a message, not an ownership transfer. The coordinator assigns a new
hypothesis to its recipient. Live agent messages provide timely notification;
the persisted event log provides restart/handoff context.

Submissions are content-addressed JSON copies. Evidence is `numerical` or
`exact_claim`; neither makes a submission a verified incumbent. Cost fields
are submitting-agent claims pending independent checks. Coordinator acceptance
requires a separate review record with exact gate specification, reproducible
certificate commands, independent results and scope. Keep incumbents separate
by elementary gate set; never compare an arbitrary two-qubit unitary with a
fixed entangler merely by calling each one gate.

The accepted exact13 CNOT incumbent is recorded in `EXACT13_ACCEPTANCE.json`.
The exact14 certificates remain preserved baselines. First-round native work
establishes a separate gate-model baseline; its ten entanglers are not ten CNOTs.
There is no automatic search scheduler: the active coordinator enforces budgets
and dispatches follow-ups. Saved board records do not prove that processes are
alive. Check process/session handles before recovering interrupted runs.

## First round

- Algebraic: identify constructions outside exhausted fixed-prefix families;
  share an analytic insight with the native-gate researcher.
- Numerical: one new bounded architecture experiment with frozen proposal
  provenance, full local freedom, saved gates, independent numerical checks.
- Native gates: explicit matrix conventions, checker tests, and a native
  regrouping of the existing exact circuit checked against its compiled list.
- Coordinator: board, ownership/dependency tests, incumbent verification,
  independent review and final shared report.

`board.sqlite3`, submissions, work directories and round reports are kept
research artifacts; SQLite journals and Python caches are ignored.

## Current results

See `incumbents.json` for separate gate-set records, `ROUND_13.md` for the latest
bounded work, `CHATS.md` for actual recorded exchanges, and `ROUND_1.md` for the
first team's historical outcomes. These
files are reviewed records, not automatic promotions from submissions.
