# Persistent objective and completion requirements

Objective: reduce the CNOT count, coordinate agents to find the optimum and
construct the associated circuit explicitly, and show their actual chats.

Previous goal-related turn: progress. It accepted the completed exact13
certificate, synchronized the authoritative incumbent, and completed a new
bounded12CX experiment. It did not prove the optimum.

Current starting revision61efeb1 is clean. All372 ledger runs and19 hypotheses
are terminal; completed round2 handles are not reused as running jobs.
No live optimizer/certificate process was found. New round3 hypotheses20–22
have distinct owners, hypotheses and finite budgets.

Requirements for full completion:

1. An explicit chronological native-CX circuit attaining the proved minimum,
   improving the incumbent wherever a smaller circuit exists, independently
   validated and exactly certified. If the incumbent is actually optimal,
   proving that fact is necessary; numerical loss alone is insufficient.
2. A global lower bound matching the construction in the ancilla-free,
   arbitrary-one-qubit, all-to-all directed-CX model, with arbitrary eigenvalue
   order and basis freedom in degenerate eigenspaces. A fixed-block obstruction
   or failed optimizer cannot prove optimality.
3. Readable exact gate prescription for that circuit, source hash, conventions,
   verifier commands and evidence with proof limitations.
4. Actual operational researcher exchanges, with author/recipient/provenance,
   visible in chat artifacts rather than an invented conversation.
5. Preserved prior work, terminal ledger rows for finished runs, scoped
   adversarial review, focused validation and local commit; no push.

At the start, requirements1–2 are incomplete: current exact upper13, lower6.
The accepted exact source is explicit JSON but its readable presentation and
current-round chats are being improved. The persistent goal stays active until
the entire objective is verified; finishing a bounded round is not completion.

Closeout evidence: runs373–374 completed two distinct bounded mechanisms;
all saved new native gate lists fail independent validation. Requirements1–2
remain incomplete. Requirements3–4 now have concrete artifacts: source-bound
EXPLICIT_CIRCUIT.md and verbatim board/peer CHATS.md. Neither proves optimum.
The independent joint-axis lemma lacks a topology-level forcing argument.
Goal status remains active; no completion, pause or blocked update is justified.
