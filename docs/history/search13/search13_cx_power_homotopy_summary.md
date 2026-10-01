# CNOT-power continuation toward a 13-CNOT diagonalizer

The exact 14-CNOT circuit uses the four-wire right shift convention of
the repository. For one selected CNOT with control $c$ and target $t$,
run 214 replaced it by

\[
C_{c,t}(\alpha)=H_t\,\operatorname{CP}_{c,t}(\pi\alpha)\,H_t,
\qquad 0\leq\alpha\leq1.
\]

Here $\operatorname{CP}(\phi)$ multiplies only the computational state with
both selected bits equal to one by $e^{i\phi}$. Consequently
$C_{c,t}(1)=\operatorname{CX}_{c,t}$ and $C_{c,t}(0)=I$; a direct matrix
check gave maximum error $2.31\times10^{-16}$ at $\alpha=1$. The gate is
locally equivalent to an Ising $ZZ$ interaction at intermediate angles.
All one-qubit layers between successive entanglers were freely varied
within $\mathrm{SU}(2)$. The label-free objective was
$\|\operatorname{offdiag}(UV_4U^\dagger)\|_F^2/16$.

Run 214 tested CNOT positions 3, 6, and 12, each with two independent
local-layer starts. The prescribed stages were
$\alpha=(0.875,0.75,0.5,0.25,0)$, with L-BFGS-B capped at 100 iterations
per stage. All six $\alpha=0$ gate lists have exactly 13 CNOTs and failed
the independent cycle checker. The best endpoint, position 6 seed 14202,
had objective 0.0688143997 and maximum off-diagonal entry 0.3040636428.
Direct NumPy reconstruction of its saved gates matched the optimizer loss
within $8.75\times10^{-16}$.

For position 6, the two starts reached losses around $10^{-7}$ or less at
$\alpha=0.5$ but around 0.0636 at $\alpha=0.25$. Run 217 replayed seed
14202 deterministically to $\alpha=0.5$ and saved its full 15-by-4-by-3
local-angle array. The replayed losses matched run 214 exactly; the
half-angle snapshot SHA-256 is
2d31e0afa6b90902cbf28697b91eaf260b8feff1499265369f45415ce3f8b8be.
It continued through ten smaller steps from 0.45 to 0, with a
250-iteration cap and two small perturbations at each of 0.35 and 0.25.
The endpoint objective was 0.0688142950, maximum off-diagonal entry
0.3037625286, and fourth-root eigenvalue error 0.3496376831. The saved
gate list was independently checked and still fails to diagonalize $V_4$.
All 11 parameter snapshot hashes were verified.

Run 221 tested a second movable entangler. It drove the power of CNOT 6
through $\alpha=(0.75,0.5,0.35,0.2,0)$ while also optimizing the power
$\beta$ of either CNOT 7 or CNOT 8, all local $\mathrm{SU}(2)$ layers, and
a small growing quadratic penalty for $\beta-1$. At the final stage
$\beta$ was fixed exactly to one, leaving a 13-CNOT gate list. There
were two starts per partner and 150 iterations per stage. In the four
tested paths, $\beta$ stayed within 0.00073 of one at every intermediate
stage, so this penalized experiment barely moved the partner gate.
The best endpoint had objective 0.0688142958 and maximum off-diagonal
entry 0.3038136554. It failed the independent checker. All 20 local
parameter and partner-power snapshots were hash-verified.

These are bounded numerical paths. Several stages reached their
iteration caps, and even a converged local minimum would not give a
global lower bound. The next genuinely distinct continuation would
allow a wider or unpenalized partner-power excursion before imposing
the exact-CNOT endpoint, or move an entangler to another edge. A
fixed-family analytic obstruction would require a separate proof.
