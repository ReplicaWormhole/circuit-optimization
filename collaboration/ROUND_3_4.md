# Continued search, explicit circuit and researcher chats

The persistent objective is to reduce CNOT count, prove the optimum, construct
the corresponding circuit explicitly, and show actual researcher exchanges.
The goal remains active. The optimum has not been established: the exact CX
incumbent is 13 and the preserved global lower bound is six.

Starting revision:61efeb1, clean worktree. All 372 prior runs and 19 hypotheses
were terminal. Historical/completed agent handles were not treated as active
jobs. Three new owned researchers claimed hypotheses20–22, then reused their
own handles for distinct hypotheses23–25 with a separately stated budget.
No existing exact13 certificate assignment or execution was duplicated.

## Visible chats and explicit circuit

`CHATS.md` exports verbatim board messages with event IDs, timestamps, authors,
recipients and message kinds. Each researcher also retains live sent/received
messages in its own CHATS.md, including corrections. This is operational team
communication, not an invented reconstruction or a private-reasoning dump.
The coordinator requested the chat and circuit files in Codex panels; the app
returned queued. Both artifacts are available directly in the workspace.

`work/verifier_round3/EXPLICIT_CIRCUIT.md` is a readable full chronological
prescription:73 gates,13 CX,60 local rotations. Signed exact half-angle radicals
preserve rotation branches; full-angle and rational-pi metadata, output labels
and their root mapping are supplied. Exporter verifies the accepted source hash.
Independent table parsing matched every row and metadata field. This is the
current exact circuit, not a proved optimal circuit; transcription checking
does not add another full-matrix certificate.

## Structural findings and independent review

For C=CX(a,b), E=CX(b,c), F=CX(a,c), and tensor-product local layers L,M
commuting with C, chronological C,L,E,M,C equals L,E,F,M. The bare identity
follows from binary updates. For product local gates, the C-centralizer is
control-Rz and target-Rx up to scalar phases, with arbitrary spectators.
Independent block/partial-trace arguments prove this class; failed fits do not
exclude it. The accepted circuit lacks an immediately applicable aba motif.

The first architecture motivation incorrectly described a one-edge mutation.
The actual precursor requires two mutations [01,32,12] to [01,12,01]. This was
corrected in reports and subsequent messages; the original message remains in
the chat log. Direct twelve-CX schedules were already frozen and unaffected.

A joint halfshift lemma was independently reviewed. With UV4U*=D and
A_j=U* p_j U, two Pauli axes on distinct input wires, each having nonzero
off-diagonal entries, cannot both have images commuting with V4 squared.
Each commutation forces D squared labels equal across one bit flip; two flips
would make its plus/minus multiplicities divisible by four. V4 squared has
multiplicities10 and6, contradiction. No support restriction is required.
No schedule-specific forcing argument has been supplied, so this lemma does
not establish a seven-CX bound or exclude any of the3024 surviving schedules.

## Run373: direct changed twelve-CX motifs

One reserved run, two new schedules, one seeded start each, full local freedom,
max300 iterations and120 seconds per case,240 seconds total, one thread.
Slots0..2 or6..8 of the canonical schedule become [None,12,02]. The full native
lists are saved. They differ from unchanged single-deletion schedules and
run372; historical equivalence-class novelty is not claimed.

| Seed | Iterations | Ordinary loss | Maximum offdiagonal |
|---|---|---|---|
|926301|130|0.26146233225017707|0.4487654196973814|
|926302|157|0.25857489404993894|0.6728896343792273|

Elapsed5.864 seconds. Both saved12-CX gate lists are invalid, independently
confirmed by the verifier and a separate NumPy reconstruction. Run373 is
terminal failed; hypothesis21 is inconclusive. Optimizer convergence is not
diagonalization or a topology exclusion.

## Run374: explicit thirteen-to-twelve bridge

A distinct new run fits each true [01,12,01] thirteen-CX precursor, imposes the
two interior local commutations, collapses it using the reviewed identity, then
polishes the twelve-CX endpoint. Per motif:100 unrestricted13 iterations,
100 hard-constrained13 iterations,200 unrestricted12 iterations;120 seconds
each,240 total, one thread. The middle local layer becomes identity and its
matrix moves into the next layer by the chronological product L3@L2.

Both unrestricted13 precursors already fail, maxoff0.6281 and0.6651. Thus no
passing circuit was lost merely through the commutation restriction. Independent
saved-gate matrix comparisons validate constrained13 to collapsed12 up to global
phase at5.36e-16 and3.56e-16. All eight saved gate lists have the correct13/12
counts and fail diagonalization.

| Motif | Seed | Final12 maximum offdiagonal |
|---|---|---|
|0..2|926401|0.6104871878877347|
|6..8|926402|0.6674475006279953|

Elapsed10.649 seconds. Run374 is terminal failed; hypothesis23 is inconclusive.
The collapse implementation is validated; no smaller diagonalizer, exact12
certificate, topology exclusion or optimality proof follows. No additional
unchanged fits were launched.

## Next action and remaining requirements

The bridge approach needs a new initialization mechanism, since its free13
precursors failed. A next bounded hypothesis should use spectral-label changes
or a fresh analytically structured basis, rather than deepen these same local
fits. Any new budget and complete proposal must be stated, claimed and reserved
before execution. A new passing gate list needs independent exact certification
before promotion. The joint lemma needs an actual schedule-level implication
before it can strengthen the global lower bound. Until construction and global
lower bound match, the persistent optimum goal is incomplete.

All six hypotheses and both new scientific runs are terminal. Research JSON,
scripts, chats, reports, code snapshots and databases are kept artifacts;
caches/journals are ignored. No dependency install, destructive cleanup or
remote action. Coordinator performs final audits/review and a local commit.

Final verification:22 focused tests pass, with SQLite ResourceWarnings from
the existing test suite recorded in the coordinator review. Ledger audit374
runs/189 candidates has no problems; board auditfour submissions passes content
hash/JSON checks. Independent rendered-circuit fidelity and bridge matrix
checks pass. See `work/root_round3/ADVERSARIAL_REVIEW.md` for corrections and
residual risks. Compilation/testing is not an optimality proof.
