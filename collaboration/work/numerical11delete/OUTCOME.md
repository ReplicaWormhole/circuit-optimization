# Run409: two near-diagonal11CX candidates, neither accepted

The root executed the frozen one-shot12case diagnostic after independent
preflight408 passed. Run409 completed in7.108003442 seconds:1067objective/
gradient calls,951total iterations,12cases, all11CX. No case hit its8-second
or the total110-second cap; the slowest optimizer case took1.01684seconds.
Eleven cases reached the80-iteration limit; deletion4 stopped after71iterations
on SciPy's relative function-reduction criterion. No case passes the frozen
shared-checker tolerance1e-9. The exact12 certificate is unchanged.

## Best candidates and actionable stop reasons

| Deleted CX ordinal | Source gate position / edge | Best loss | Max off-diagonal | Best gradient norm | Stop |
| --- | --- | --- | --- | --- | --- |
| 9 | 29 / 3->2, first CX of final RCCX | 2.279470641322442e-15 | 5.058049457903065e-8 | 9.246505792997596e-8 | 80 iterations |
| 7 | 16 / 3->2, last CX of merged B/P | 7.039082434994391e-15 | 8.987102717635442e-8 | 1.321867847733767e-7 | 80 iterations |
| 11 | 33 / 3->2, last CX of final RCCX | 0.06600168357225329 | 0.34435171940094095 | 4.3945895151319397e-7 | 80 iterations |

The other best losses were approximately0.535169,0.551990,0.375,0.125,
0.125,0.125,0.125001,0.25 and0.125000017. Small gradients at nonzero loss
are numerical endpoint diagnostics; no Hessian or local-minimum proof was
performed. For the two near-zero-loss candidates, the80-iteration cap rather
than a verified precision floor terminated the fit. Their residuals exceed
1e-9, so they remain unsuccessful numerical candidates, not exact11CX results.

The global best saved gate list is
`experiments/runs/409/best_candidate.json`, SHA256
`f67d8d40a7458f7813b3cba5e16a2e49f54a16c94ed916072e1e0bd584f1e429`.
Per-case complete base, initial and best lists, source deletion lists, local
base matrices and vectors are retained in409/cases/. Initial/best serialized
lists have59gates:48U3 and11CX. All144 coordinates were free. The initial
vectors equal the saved source-collapse base plus the corresponding fixed
PCG64 perturbation from the twelve-draw seed10012101 allocation, sigma0.025.

## Read-only provenance and output audit

The saved run script, family, input candidate, config, global noise file and
best word hashes match the frozen preparation packet and result hashes.
All12 per-case candidate hashes and result.json copies match their global
records. Pure JSON audit confirms11CX in each deleted-source/base/initial/best
word and144-dimensional base/noise/initial vectors. Every saved initial vector
is exactly the recorded componentwise base-plus-noise sum. No circuit matrix,
RNG or optimizer was executed during this outcome audit.

Source candidate remains the run406/407 exact12 word,
SHA25658a862a4acef4e768c4613e9b142504592f5237a68539941148653a7f083eea2.
Main script SHA256f90819b38d17156c1f20bbf4b21274f4023a4f39829f52e524d8a7426cc5593c;
family SHA256b7361e44000a2a0221a61fdef1173c6d68f509995aa279c3b9446f9b3dd04e0b;
config SHA256277e8509a1eb57531fe800f79c646b1c1e6564fbc3c66c6641576cf3cb7f9748.
Runtime versions and single-thread settings match the frozen config.
Independent preflight408 checked132 matrices,12AD directions and24finite-
difference losses across all12 base/initial cases. Its maximum phase-adjusted
matrix error was1.78e-15 and maximum directional gradient discrepancy1.34e-9.
Separate postfit numerical validation is assigned to verifier_sol and must
pass before interpreting the near-hit architecture as a refinement source.

## Next bounded hypothesis

The coordinator explicitly authorized preparation of one NEW SVD/minimum-norm
Newton refinement of deletion9 only, keeping all144 free local coordinates.
It has a separate20step/8backtracking/60-second/one-thread seed0 budget with
no RNG or restart, and must receive a fresh board/ledger reservation after
independent postfit validation. This is not continuation inside run409 and
does not silently refine deletion7. Its full real/imaginary off-diagonal
residual Jacobian and one deterministic finite-difference sanity direction
will test whether the remaining error can be reduced below tolerance.
Exact recognition/certification remains a further separate requirement.

Deletion7 points to a three-CX merged-phase block, while deletion9 points to
a two-CX final conditional basis block when all local gates change. Both
findings were shared with algebraic_sol and verifier_sol. They are structural
leads rather than proven syntheses or topology exclusions. The alternative
router rewrites from algebraic_sol remain distinct untested hypotheses.

Adversarial output review found no confirmed source/input/hash/count/stop-
reason defect. Numerical failure under finite caps establishes no11CX lower
bound, no family exclusion and no local minimum. All outputs are retained
research evidence; caches are ignored. Root owns ledger finish, board close,
independent postfit acceptance and the focused local commit.
