# Artifact index

Commands below run from the repository root. Links are relative to this index.
Start with the current status and accepted circuit; use the research index for
derivations and bounded lower-bound investigations.

| Purpose | Starting point |
| --- | --- |
| Current result and next question | [Current status](STATUS.md) |
| Where new material belongs | [Repository layout](STRUCTURE.md) |
| Accepted exact 13-CNOT circuit | [Gate list](../collaboration/work/algebraic13_new/canonical_algebraic_ansatz.json), [acceptance](../collaboration/EXACT13_ACCEPTANCE.json), [explicit prescription](../collaboration/work/verifier_round3/EXPLICIT_CIRCUIT.md) |
| Preserved exact 14-CNOT circuit | [Gate list](../topology14_exact_matchgate_rational.json), [derivation](research/constructions/topology14_exact_matchgate_derivation.md) |
| Lower bounds and constructions | [Research index](research/README.md) |
| Other gate sets | [Incumbent index](../collaboration/incumbents.json), [native-gate conventions](../collaboration/native_gates.md) |
| Numerical and exact checkers | `check_circuit.py`, `exact_check.py`, `native_gate_check.py` |
| Candidate format | [Gate-list reference](CANDIDATES.md) |
| Scientific attempts | `experiments.sqlite3`, `experiment_log.py`, [run protocol](EXPERIMENTS.md) |
| Search snapshots | `candidates/`, `code_snapshots/`, [legacy scripts and outputs](../archive/README.md), [archive integrity check](../archive/verify.py) |
| Collaboration | [Board protocol](../collaboration/README.md), `collaboration_board.py`, [round 13](../collaboration/ROUND_13.md) |
| Manuscripts and local sources | [Writing index](../writing/README.md) |
| Historical narrative | [Search chronology](history/SEARCH_14_STATUS.md), [12-CNOT summaries](history/search12/README.md), [13-CNOT search notes](history/search13/README.md), [former full README](history/README_before_reorganization.md) |

The ledger, acceptance record, and original source files provide different
types of evidence. Use the acceptance record to identify certified circuits;
use the ledger for bounded attempts; use the board for ownership and handoffs.
