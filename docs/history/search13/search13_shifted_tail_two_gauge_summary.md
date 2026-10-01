# Shifted eight-CNOT tail: two-mask gauges and seven-CNOT fit

The exact 14-CNOT circuit has a six-CNOT prefix in gates 0–12. The
remaining eight-CNOT tail begins at gate 13 and contains both first-layer
matchgates, the center one-qubit layer, and the two terminal matchgates.
Replacing this tail with seven CNOTs would give a 13-CNOT cycle
diagonalizer, provided the replacement implements the same tail up to an
allowed output diagonal gauge.

Run 209 used exact modular specializations in $\mathbb F_{193}$ and
$\mathbb F_{769}$ to screen 36 unordered pairs of cross-terminal-pair
parity masks. Each angle numerator was chosen independently from
$\{-4,-2,2,4\}$ in units of $\pi/16$, giving 576 fixed gauged tails.
The four single-wire operator-Schmidt ranks were 4 in every case, and
the ranks across the two other balanced cuts were both 16. Across the
terminal-pair cut, 468 cases had rank 16 and 108 had rank 14. Every case
passed a necessary seven-CNOT cut-graph screen; the minimum graph count
allowed by these ranks was six. Independent complex matrix ranks
matched all 576 modular profiles.

Run 211 selected one rank-14 case: output gauge
$\exp(-i\pi Z_5/8)\exp(i\pi Z_6/8)$, where masks 5 and 6 denote the
binary parity masks on the four output wires. Two distinct seven-CNOT
star-plus-chord schedules satisfy all seven cut inequalities for this
target. Each schedule received two arbitrary local $\mathrm{SU}(2)$
initializations, and each L-BFGS-B fit was capped at 300 iterations.
The loss was phase-insensitive process infidelity to the fixed gauged
tail. The best loss was 0.6222087667; direct NumPy reconstruction of the
saved 13-CNOT gate list reproduced it. Its maximum cycle
off-diagonal entry was 0.8661734447. The smallest cycle off-diagonal
entry among all four fits was 0.7715541754. None passed the independent
cycle checker.

The rank screen is necessary, and the four failed numerical fits exclude
neither seven-CNOT synthesis of another gauged tail nor a 13-CNOT
diagonalizer. A next bounded search could enumerate more chronological
seven-edge orders passing the exact cut constraints and optimize the
label-free cycle objective, allowing the output eigenbasis to change
beyond these fixed two-mask gauges.
