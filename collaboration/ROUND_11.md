# Phase-aware rewrite and restricted 12-CNOT search

The accepted exact 13-CNOT circuit remains the incumbent in the ancilla-free,
all-to-all CX plus arbitrary local-unitary model. The proved interval remains
6 <= C_min <= 13; this round neither proves optimality nor excludes 12 CNOTs.
Initial clean commit6a84fce was reconciled using the live agent registry,
process inspection, board and ledger. Completed handles were reused only after
verification. No existing incumbent certificate assignment was duplicated.

## Exact-preserving orbit rewrite: runs383–384

Numerical board51 reserved a deterministic commuting/cancellation pass before
execution:100sweeps,50000rule checks,120seconds, one thread. Actual1225checks
and0.00089seconds merged five rational-pi same-axis pairs. The complete circuit
changed149 to143 gates, with67 CX unchanged and82 to76 local gates. Independent
board52 replayed the trace and checked matrix equality. Root board57 reserved
one separate180second exact-check invocation on the changed gate list;
run384 passed in approximately0.47seconds, conductor64, all entries of UV4-DU
zero and eigenvalue multiplicities6,3,4,3. This is a certified higher-cost
construction, not an incumbent improvement.

Changed full list: `work/numerical_label12/joint_simplify/candidate.json`,
raw SHA25699742270d05745155e74379602293bfb7bf2c19782d57def349521a570a000f5.
Certificate: `work/root_phase/exact_output.json`. Archive normalized hash:
d1ed2d2cccbcce23cfc3f74c0c6890db4817ff0a1dbd862d84e33b39e4b8bc4b.

Algebraic board53 derives the phase condition for a fixed per-orbit Fourier
basis: r(j+1) conjugate(r(j)) must be constant on each orbit, so
r(j)=r(0) rho^j with rho^m=1. Orbit-constant phases are sufficient but not
necessary when eigenvalue labels can move. This condition is restricted to
the fixed orbit basis; it does not capture all degenerate-eigenspace freedom.
Relative-phase reversible replacements require phase verification.

## Competitive restricted family: run385

After the rewrite produced no CX saving, a new finite budget was explicitly
assigned to boards54–56. Chronological sequence:
L0; CX01,CX21,CX31; L1; F02; L2; F13; L3; CX12;
L4; F01; L5; F23; L6, with F(a,b)=exp(i a XX+i b YY).
Seven arbitrary four-wire local layers contribute84 parameters and four F
blocks contribute8. The three-CX prefix has no interior local layers.
Arbitrary L0 means this is parity of the current basis bits, without an
asserted original Hamming-weight tag or preserved downstream invariant.

In the precisely defined hybrid {CX,F,arbitrary locals} model the family
uses4 CX+4 F, eight native entanglers. Each F has a two-CX compiler, yielding
a12-CX ceiling. These are resource counts for a tested family, not an upper
bound on diagonalization: no passing member was found. The assigned verifier
owned the new native helper and approved actual source/config/helper hashes
before reservation. Two frozen Gaussian0.3 starts, seeds9261101/9261102,
500L-BFGS iterations each,240seconds total, one thread, no extra polishing.
Actual run385 took2.228609seconds; the initially rejected assumed384 board
link was corrected immediately to the actual reserved385 before fitting.

Both starts failed. Seed1:94iterations, off-diagonal squared loss/16
0.592223999630, max off-diagonal0.943861941802. Seed2:128iterations,
loss0.500000000000, max off-diagonal0.707106786628. Independent SciPy
matrix construction verifies saved92parameters, seed initialization and
chronology; native/compiled matrix differences are below4.18e-16. Focused
kernel gradient and15 embedding/compiler cases pass. Full gate lists,
checkpoints, configs and actual direct exchanges are retained in the three
`parity12/` directories and central `CHATS.md`. No exact certification is
warranted for these failing lists; no failure is a family exclusion.

## Next-round recommendation (not dispatched)

Return to the lower-loss run377 12-CX checkpoint and test a distinct, frozen
family with additional local freedom inside the parity prefix, rather than
repeating these two failed restricted starts. Specify and independently validate the mapping of the168-coordinate run377
checkpoint into the new family; if no such mapping exists, use explicitly
frozen fresh seeds instead of claiming a warm start. First compare unordered CX-pair
sequences and parameter rules against the ledger. Propose at most two warm
starts,200 iterations each,120seconds total, one thread, with full seed/source
hashes before reservation. Independent native/compiled gate-list checks and
exact recognition/certification remain mandatory for any passing candidate.

For a lower-bound route, include actual six-CX gate-template reachability
constraints before attempting symbolic elimination. Full-support commuting
observables alone are vacuous as a strengthening once all light cones are
full. The existing3024 surviving six-CX skeletons remain unchanged.
Do not dispatch another search without a newly frozen hypothesis and budget.
