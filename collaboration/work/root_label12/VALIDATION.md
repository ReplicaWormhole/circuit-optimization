# Validation and adversarial review

Artifact types: math, computation, research scripts, workflow.

The focused repository command
`python3 -m unittest test_collaboration_board.py test_check_circuit.py test_exact_check.py test_exact14.py test_experiment_log.py test_native_gate_check.py`
passed all22 tests in2.522 seconds. The existing suite emitted ResourceWarnings
about unclosed SQLite connections; no unrelated repair was attempted. These
tests check infrastructure and selected checker behavior, not mathematical
optimality or new candidate exactness.

Independent verification checks all four complete gate lists, source optimizer
roundtrip (1.45e-15 up to phase), selected pairs/labels, fixed-target objectives,
counts and hashes. All candidates fail. No new exact proof is claimed.

Final `experiment_log.py audit`:375 runs,190 archived candidates, no problems.
`collaboration_board.py audit`:four submission hashes/JSON valid; this does not
certify mathematics. Ledger summary:zero running. Both SQLite integrity checks
return ok. Helper hashes and the immutable search-code archive match. Process
inspection found no new optimizer or existing exact13 verifier process.

Adversarial review used `~/.codex/adversarial_review.md`, examined the numerical
code, source layout, independent audit and all shared document changes. No
confirmed defect remains. Initial alternative swap ranking is explicitly
superseded, fixed-target losses are distinguished from label-free offdiagonal
loss, and the timing guard is described as cooperative. The reduced-target
partial-trace invariant for the recommended joint swap was independently checked
with exact Gaussian-integer sums; it proves only inequivalence under free
product-local output conjugation, not success or a new bound. Existing exact13
certification is not re-executed. The mathematical optimum remains unresolved.

`git diff --check` passes. Generated research JSONs, scripts, reports, chat logs,
scientific code/candidate archives and database snapshots are retained evidence.
No pre-existing changes were present or overwritten; no remote operation.
