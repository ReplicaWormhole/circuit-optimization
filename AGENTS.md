# Circuit optimization working agreements

## Repository purpose

Find exact ancilla-free diagonalizers of the four-qubit right shift. Main
artifacts: code, math, computation, source-backed notes, workflow.
Read `README.md`, `docs/STATUS.md`, and `collaboration/README.md` first.
This directory is now an independent nested Git repository. Parent history
and tracking remain intact; run Git commands from this directory.

## Important paths

- Verification: `check_circuit.py`, `exact_check.py`, `delete14_integer_audit.py`.
- CX incumbent: `collaboration/work/algebraic13_new/canonical_algebraic_ansatz.json`,
  accepted in `collaboration/EXACT13_ACCEPTANCE.json`. Preserve the exact14
  baseline `topology14_exact_matchgate_rational.json` and its derivation.
- Scientific attempt ledger: `experiments.sqlite3`, managed by `experiment_log.py`.
- Collaboration: `collaboration_board.py`, `collaboration/board.sqlite3`.
- Saved evidence: `candidates/`, `code_snapshots/`, `collaboration/submissions/`.
- Historical flat search namespace: `archive/legacy_root/`.
- Manuscripts: `writing/drafts/`. Local source materials in `writing/sources/`
  are ignored and must not be added to Git.
- Historical search chronology: `docs/history/SEARCH_14_STATUS.md`.
- Preserve archived search files, immutable snapshots, and previous ledger rows.

## Commands

Existing Python environment supplies NumPy, SciPy, Torch and SymPy. Qiskit is
optional and is not installed in the system `python3` used for this setup;
install no dependencies without explaining the need. No global formatter or
build is required.

```bash
python3 -m unittest test_archive.py test_collaboration_board.py test_check_circuit.py test_exact_check.py test_exact14.py test_experiment_log.py test_native_gate_check.py
python3 archive/verify.py
python3 check_circuit.py topology14_exact_matchgate_rational.json
python3 exact_check.py topology14_exact_matchgate_rational.json
python3 delete14_integer_audit.py topology14_exact_matchgate_rational.json
python3 experiment_log.py summary
python3 experiment_log.py audit
python3 collaboration_board.py list
python3 collaboration_board.py events
python3 collaboration_board.py audit
```

## Collaborative research contract

Use a coordinator and three collaborators: algebraic construction, numerical
architecture search, and alternative gate sets/independent verification.
They must share findings directly and record actionable handoffs on the board.
Delegate bounded independent hypotheses, not duplicate unchanged searches.
Claim a board hypothesis before work; reserve every numerical/symbolic search
in the scientific ledger before execution, then link its run to the board.
Record seed, complete topology rule, optimizer, iteration/time limit, code hash,
best candidate and limitations. Create `experiments/runs/<ID>/` with
`python3 experiment_log.py workspace ID --code PATH` after reservation; keep the run plan,
script, outputs and conclusion together. Finish unsuccessful runs too. Structural
read-only investigations need board outcomes but no fictitious search run.

Use per-agent directories under `collaboration/work/<agent>/`; never overwrite
another agent's files. Coordinator owns shared board code, policies and commits.
Only the assigned verifier changes new native-gate code. Coordinate shared
checker edits before making them. Limit each optimizer process to one compute
thread; first-round budget is one bounded run per numerical researcher.
Do not start unlimited follow-up loops without a newly stated resource budget.

## Math and verification contract

Qubit zero is most significant. Gates are chronological; the target is
`U V4 = D U`, with arbitrary eigenvalue order and arbitrary orthonormal bases
inside degenerate eigenspaces. Exact diagonalization is distinct from total-spin
diagonalization and a strong Schur transform.
Keep numerical evidence, exact claims, and independently certified results
separate. Optimizer failures are not topology exclusions or global lower bounds.
Compare native entangler counts and depth only within precisely defined gate
sets; retain a separate CNOT compiled cost where known. Unknown cost is not zero.
The six-CNOT bound does not transfer automatically to other gate sets.
Never promote a submitted candidate to incumbent without independent gate-list
validation and exact certification. The board does not certify proofs.

## Cleanliness and closeout

Preserve user work. Ignore Python caches and SQLite journals; retain research
JSON, notes, scripts and databases as intentional evidence. No destructive
cleanup. Review diffs using `~/.codex/adversarial_review.md`, run focused tests,
and make a local commit after verification. No push or remote configuration
without explicit authorization. Report remaining changes and unresolved risks.
For this collaboration, stage only its assigned files and evidence; leave
changes produced concurrently by other research tasks unstaged. The shared
experiment ledger may be committed as a whole snapshot; do not rewrite its rows
to isolate one team's work.
