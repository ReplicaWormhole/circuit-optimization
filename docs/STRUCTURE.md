# Repository layout

Use the repository root for the entry points and state that existing commands
and collaboration code import directly. Keep new work in the directory for its
purpose rather than adding another root-level search file.

| Material | Location |
| --- | --- |
| Current result, instructions, and navigation | `README.md`, `AGENTS.md`, `docs/STATUS.md`, `docs/INDEX.md` |
| Reusable checkers and command-line tools | Existing root `*.py` files; keep their import paths stable |
| Preserved circuit inputs, seed files, and scientific ledger | Existing root JSON files and `experiments.sqlite3`; acceptance records in `collaboration/` |
| New bounded numerical or symbolic run | `experiments/runs/<ID>/`, after ledger reservation; keep plan, code, output, and conclusion together |
| Content-addressed candidates and code snapshots | `candidates/` and `code_snapshots/`, managed by the scientific ledger |
| Agent hypotheses, handoffs, and submissions | `collaboration/` and agent-owned `collaboration/work/<agent>/` |
| Mathematical notes | `docs/research/lower_bounds/` or `docs/research/constructions/`; label what is proved, conditional, or exploratory |
| Search chronology and summaries | `docs/history/`, grouped by search family |
| Historical scripts and saved outputs | `archive/legacy_root/`; preserve their bytes and use `archive/verify.py` |
| Manuscripts | `writing/drafts/`; local source papers and notes in ignored `writing/sources/` |
| Tests | Existing root `test_*.py` files until the executable import surface is migrated |

The root still contains compatibility entry points and a few JSON inputs
because current commands and recorded artifacts use those paths. Do not
move one without checking imports, ledger references, and the acceptance
records. The accepted 14-CNOT derivation has a root link because an immutable
submission names that path.

The archive is a frozen research record, not the destination for new results.
The scientific ledger tracks attempts; the collaboration board tracks
ownership and handoffs. Neither certifies a mathematical claim by itself.
