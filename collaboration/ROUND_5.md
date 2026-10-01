# Bounded spectral-label round

Read-only reconciliation found a clean tree at fb16033, three completed owned
round3 workers, zero scientific processes and 374 terminal ledger runs. Agent
registry and process inspection were independent of ledger status. No exact
certificate assignment was duplicated. Details: `work/root_label12/RECONCILIATION.md`.

Three distinct board assignments26–28 covered algebraic output-label freedom,
one numerical experiment and independent gate-list verification. They exchanged
findings directly and persisted handoffs; `CHATS.md` exports actual board
messages, while each researcher retains sent/received text in its own directory.

The algebraic derivation uses traces(V^k)=(16,2,4,2) to establish eigenvalue
multiplicities(6,3,4,3) for roots(1,i,-1,-i). For a fixed unitary U, minimizing
||UV-DU||F over admissible diagonal D reduces to a capacity-constrained assignment
on diagonal entries of UVU†. Output order and all orthonormal bases within each
degenerate eigenspace remain mathematically free; these basis changes are not
automatically free in the circuit gate model. The initial suggestion to rank
swaps by diagonal target-cost increase was superseded before execution by the
frozen offdiagonal-coupling rule, without an extra numerical search.

Run375 starts from the invalid twelve-CX deletion lead
`work/numerical12/deletion_slot1_deep.json`, source SHA256
9a9dd0e0b430e935f99f4d998fb25d5e1b48160ad16a6c46a27e2574c59c2398.
Its independent maxoff is0.2588190369994409. The fixed topology replaces slot1
of its thirteen-slot schedule by identity; all168 local parameters are free.
Four different-label pair swaps are ranked by descending
|B_rs|²+|B_sr|² with lexicographic pair tie breaks, B=UVU†. Each frozen target
receives at most100L-BFGS iterations and three minimum-norm Newton steps with
SVD cutoff1e-6 and at most eight backtracking trials. Total computation budget
is120seconds including screening, one BLAS/Torch thread, no random seeds.
Code/config/helper/source hashes and complete schedule are retained before
execution. Runtime checks occur between operations, not by preempting a
single Jacobian/SVD call. Actual total14.848544887seconds is below the cap.

| Swapped rows | Saved twelve-CX maxoff | Valid diagonalizer |
| --- | ---: | --- |
| 2,6 | 0.2588190451027492 | No |
| 3,7 | 0.2588190451147769 | No |
| 8,13 | 0.2794011825428628 | No |
| 9,12 | 0.2794011768358336 | No |

All four fits reach the100-iteration limit; Newton corrections do not produce
a passing list. Full chronological candidates and frozen target labels are
saved under `work/numerical_label12/`. Independent matrices use explicit
computational-basis CX columns and Kronecker U3 operators, with no shared
simulation import. See `work/verifier_label12/REPORT.md` and independent JSON.
No exact-certification target and no incumbent promotion follows. A failed
finite target search excludes neither other labels, eigenspace bases, initial
points nor the topology; it proves no new lower bound.

Costs remain separate. In the ancilla-free all-to-all directed-CX plus arbitrary
one-qubit model, the accepted exact circuit has13CX and global interval
6 <= C_min <= 13. Its full signed-radical gate prescription remains
`work/verifier_round3/EXPLICIT_CIRCUIT.md`. In the separate native set
{CX,F(a,b)=exp(i a XX+i b YY), arbitrary one-qubit gates}, independently tunable
real a,b give the preserved baseline6CX+4F, native depth8, compiled-CX upper14.
Its native JSON has numerical checking and analytic provenance, not a direct
exact native certificate. The six-CX bound does not transfer to that model.

Recommended next round: compose swaps(2,6)(3,7) on this source, alongside one
analytically chosen four-row target, rather than repeat single swaps. The first
composition flips wire1 conditioned on wire0=0 and wire2=1, ignoring wire3;
it is not a whole-wire localX. Partial trace onto wire1 has eigenvalues
{1+i,1-i} for the base target and {3-i,-1+i} for this composition. Product-local
output conjugation preserves these eigenvalues, so this target is not merely
a free local reordering of the base. This distinguishes a target, not an
optimizer basin or a circuit construction. See the algebraic report for the
comparison to both component single swaps. Freeze at most two targets,200
iterations and five SVD corrections each,120seconds total, one thread; reserve
a new ledger run first. No follow-up is dispatched. A passing gate list needs
a separate exact-recognition/certification assignment; global optimality still
requires a matching lower bound.

Repository focused tests pass22/22, with pre-existing SQLite connection
ResourceWarnings recorded separately. Candidate/provenance review and final
audits are in `work/root_label12/VALIDATION.md`. Research scripts, JSON lists,
board/ledger snapshots, source archive and chat records are intentional kept
artifacts. No new dependency, checker change or remote operation is required.
