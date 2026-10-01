# Bell-label difference construction and twelve-CNOT certification

This round follows the full coupled 98-coordinate investigations through run
404. The accepted starting incumbent was the exact 13-CNOT circuit; the
exact14 baseline and all previous experiments are preserved.

## Strategy and collaboration

The structural hypothesis was to decode the opposite-pair Bell labels, store
their XOR in a difference register, and diagonalize the resulting four small
control sectors. This addresses the missing exchange-eigenbasis suffix of the
earlier full98 searches. It does not assume that an independently synthesized
circuit belongs to their fixed chronological word.

The three collaborators exchanged actionable findings directly and on the
board. At the user's request, the active collaborators were switched to
GPT-6.1 Sol: algebraic_sol derived phase repairs and a four-CNOT block merge;
numerical_sol independently audited all four block products; verifier_sol
checked the construction and prepared the literal gate-list certificates.
Interrupted cheaper-agent preparation is retained and marked unexecuted.

Run405 established an over-budget comparator: its 19-CNOT primitive circuit
passed numerical and exact cyclotomic checks. The subsequent by-hand work
replaced exact Toffoli selectors with phase-repaired relative-phase selectors,
then merged the controlled phase and first selector from five CNOTs to four.
The resulting proposed cost is 2 Bell-decoder + 2 difference-router + 4 merge
+ 1 controlled-H + 3 final selector = 12 CNOTs, with four qubits and no ancilla.

The derivation and independent branch reviews are in
`work/algebraic_sol/MERGED12.md`, `work/numerical_sol/MERGE4_STATIC.md`, and
`work/verifier_sol/MERGE4_STATIC.md`. Every phase and inactive selector branch
is retained. These structural notes alone do not promote an incumbent.

## Certification record

Run406 reserved board127 and the complete 39-gate candidate before computation.
Seed0, no randomness, one thread, one numerical invocation30s and one exact
invocation60s, external100s. Frozen hashes are in its MANIFEST.json and ledger.
Numerical0.0752s and exact0.3810s both pass, with numerical maxoff5.67e-16.
The exact cyclotomic checker verifies all 16 rows and independently predicted
root labels. Board127's original proposal says38 gates; the frozen manifest and
actual word correctly specify39. No gate changed after reservation.

Run407 reserved board128, seed0, no randomness, one thread and one90s exact
invocation. Standalone Fraction arithmetic in Q[z]/(z^8+1), z=exp(i*pi/8),
importing no existing checker, NumPy or SymPy, independently proves all 256
UV4-DU entries and all 256 UUdag-I entries zero. It discovers the same labels
and roots(1,i,-1,-i) with multiplicities(6,3,4,3); measured0.1981s. Its full
exact matrix, source/config/dependency hashes and outputs are preserved.
Both runs are closed complete/symbolic, without optimizer or retries. Each
workspace keeps the sole/best candidate, inputs, outputs, command and limits.
Standard ledger finish/audit additionally recomputes routine floating candidate
metrics during ingestion/audit, as bookkeeping validation rather than searches.

Independent source, full gate-list and certificate review passed. Acceptance is
recorded in EXACT12_ACCEPTANCE.json and incumbents.json. The13- and14-CNOT
circuits and historical acceptance records remain byte-identical. Established
interval is now6<=Cmin<=12; the lower bound is inherited, not re-proved here.

Best next hypothesis: absorb a difference-router CNOT into the merged controlled
block using output computational-permutation freedom, seeking 11 CNOTs. This is
untested and no further computation is dispatched.

## Scope

The original full98 family's exact membership and minimum remain open.
Numerical failures from earlier rounds exclude neither that family nor a
global 12-CNOT solution. Any accepted circuit here must be justified by its
complete gate list and exact certificate. Optimality, a strong Schur transform,
total-spin diagonalization, and novelty are not claimed.

## Review, verification and retained artifacts

Adversarial review followed ~/.codex/adversarial_review.md; no blocking issue
remains after independent exact certification. Twenty focused board, ledger,
numerical/exact and native checker tests passed. Existing SQLite connection
ResourceWarnings in the test suite did not fail tests. Ledger/board audits
passed; historical ledger rows through404 and board/events at initial HEAD
were unchanged. Own scripts, JSON, logs, matrices, databases, reports and
unexecuted preparations are intentional research evidence; Python caches and
journals remain ignored. No dependencies installed or remote changes.
