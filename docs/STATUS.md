# Current research status

Reviewed on 2026-10-01 through runs 387–388 and collaboration round 13. Check the live
[scientific ledger](../experiments.sqlite3) and
[board](../collaboration/board.sqlite3) before starting new work; this page is
a reviewed summary, not a process monitor.

## Established result

In the ancilla-free, all-to-all model with arbitrary one-qubit gates and
directed CNOTs, **6 <= C_min <= 13**. The
[13-CNOT gate list](../collaboration/work/algebraic13_new/canonical_algebraic_ansatz.json)
has exact signed-radical metadata and an
[acceptance record](../collaboration/EXACT13_ACCEPTANCE.json). The
[14-CNOT rational-pi circuit](../topology14_exact_matchgate_rational.json)
remains an independently checked baseline. The optimum is unknown. Cycle
diagonalization does not establish total-spin diagonalization or a strong
Schur transform.

## Latest reviewed work

The projected 24-coordinate Hessian and independent finite-difference checks
in runs 387–388 did not meet the frozen negative-curvature trigger; no escape
fit was run. Exact control-wire absorption gives an equivalent 98-coordinate
description of the relaxed family with boundary remapping, without reducing
the CNOT count. A native 4-CX plus 4-F construction compiles to 12 CNOTs,
but no passing member is known. These bounded results exclude neither the
family nor a 12-CNOT solution. See [round 13](../collaboration/ROUND_13.md).

The next suggested investigation is a bounded full coupled 98-coordinate
diagnostic. It has **not** been dispatched; it needs a frozen hypothesis and
resource budget before a run is reserved.

## Where to inspect evidence

- [Accepted and baseline circuits](../collaboration/incumbents.json)
- [Round reports](../collaboration/README.md) and [actual recorded exchanges](../collaboration/CHATS.md)
- [Historical search chronology](history/SEARCH_14_STATUS.md)
- [Handoff notes](../collaboration/HANDOFF.md)

Use `python3 experiment_log.py summary` and `python3 collaboration_board.py list`
for live local status. A numerical fit or a board submission is not an exact
certificate or a global lower bound.
