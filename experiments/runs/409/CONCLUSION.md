# Run409 bounded deletion refit

All twelve single-CX deletions of the exact12 source were fitted once, with144
free real SU2 local coordinates and PCG64seed10012101 noise sigma0.025. One
thread; max80iterations/max4000calls,8s/case,102s optimization,110s internal/
120s external limits. Observed7.108seconds; shell exit0. Scientific status is
failed because no endpoint meets tolerance1e-9; process completion is separate.

Best: ordinal9, sourcegate29 CX(3,2), first CNOT of final selector. Eleven CNOTs,
59 total gates. Loss2.279470641322442e-15, maxoff5.058049457903065e-8,
unitarityerror1.11e-15, gradientnorm9.2465e-8;80iterations89evals, iterationcap.
Second near-hit ordinal7: maxoff8.9871e-8, loss7.03908e-15,80iterations90evals.
Every other fit also fails; ordinal4 reached71iterations beforeftol stop.
The initial vectors, noise, base/deleted/initial/best complete gate lists,
histories, metrics and stopreasons for all12 cases are preserved in cases/.
No restart, extraiteration, rounding or exactcheck happened in this run.

Independent postfit verification is separately frozen/reserved beforeexecution.
The near-hit motivates a NEW bounded numerical precision-refinement hypothesis,
with its own budget/board/ledger; it is not a valid circuit or exact certificate.
Failure proves no topology/family exclusion, local minimum or global lowerbound.
Accepted exact12 remainsunchanged. Routineledgeringestion recomputesfloating
candidate metrics asbookkeeping, withoutsearch or exactclaim.
