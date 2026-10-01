# Explicit orbit-address Fourier family

Hypothesis42; structural derivation, no optimization, symbolic screening or
exact13 certificate reproduction. Read existing algebraic16_to15_derivation.md,
topology14_exact_matchgate_derivation.md and docs/history/search13/search13_joint_prefix_orbit_summary.md.
The proposed construction is an explicit abstract exact unitary family, not a
serialized new CX circuit and not an improvement over the accepted exact13.
Ancilla-free four wires, q0 MSB, arbitrary locals and all-to-all directed CX.

## Constructive lemma

The right cycle partitions computational strings into these disjoint orbits:

| Orbit in right-cycle order | New addresses |
|---|---|
| (0) | (0) |
| (15) | (1) |
| (5,10) | (2,3) |
| (1,8,4,2) | (4,5,6,7) |
| (3,9,12,6) | (8,9,10,11) |
| (7,11,13,14) | (12,13,14,15) |

Let P be exactly this computational-basis permutation. Then P V P† shifts
each quartet address by one modulo4 and the pair address by one modulo2.
Apply F4_{k,j}=i^(k j)/2 independently on address sectors q0q1=01,10,11.
Apply H2 on q3 when q0q1q2=001; retain singleton rows0,1. Since
F4 S4=diag(1,i,-1,-i) F4 and H2 S2=diag(1,-1) H2, the explicitly defined
U=blockdiag(I2,H2,F4,F4,F4) P obeys U V=D U exactly. Eigenvalue capacities
are6,3,4,3. No total-spin or strong-Schur claim follows.

The family consists of1536 orbit-address permutations: permute the three
quartet sectors(3!), choose each quartet's cyclic starting position(4³),
choose the pair's start(2), and order the two singleton addresses(2).
Every such choice has the same abstract Fourier proof. This is an explicit
finite construction family, not an already-executed enumeration. Continuous
eigenspace freedom can subsequently be used, but is not free in CX cost.

## Transparent but deliberately loose CX compilation bound

For k conditioning bits, a controlled Rz or Ry rotation has generator
Pi_control times the target Pauli, where Pi_control=product(I±Z)/2. Expansion
gives2^k commuting Pauli strings. A weight w Pauli rotation is implemented
by basis changes and a parity ladder using2(w-1) CX. Summing over control
subsets yields k*2^k CX per controlled rotation; signs for negative controls
only change angles and require no CX. A phase exp(i theta Pi_control) costs
at most sum_S 2max(|S|-1,0)=(k-2)2^k+2 CX; the empty term is global phase.
These identities specify exact rational-pi rotations for the construction.

Using R_a(theta)=exp(-i theta a/2), H=i Ry(pi/2) Rz(pi). Chronologically,
controlledH therefore applies controlledRz(pi), controlledRy(pi/2), and the
control-only phase exp(i(pi/2)Pi_control). Bounds are18CX for k=2 and58CX
for k=3. This controlled phase must not be discarded as a global scalar.

For each quartet sector, conditional positive QFT4 needs two conditionalH
gates and controlledS(pi/2) on the two low wires, under the two outer sector
conditions. The controlledS is the four-bit projector phase
exp(i(pi/2) Pi_sector Pi_q2=1 Pi_q3=1), costing34CX. Chronological H_q2,
controlledS, H_q3 gives F4 with output low-bit reversal. Output reversal still
diagonalizes and merely reorders root labels, so no SWAP is needed. Three
sector blocks therefore cost at most3*(18+34+18)=210CX; pairH costs58.

Any computational-basis permutation on16rows uses at most15 arbitrary row
transpositions. Along a Gray path of Hamming distance h<=4, an arbitrary
endpoint swap uses2h-1<=7 adjacent row transpositions. Each adjacent swap is
an X on one bit conditioned on the other three bits. Conjugating a four-bit
projector pi-phase by local H on that target implements it, at cost34CX by
the same Pauli expansion. Thus P has generic upper bound15*7*34=3570CX and
the complete abstract construction has compiled bound3838CX. This bound is
intentionally loose and establishes feasibility in the precise gate model;
it is not an estimate near the incumbent and not a3838-gate serialized list.
It uses no ancillas, native multi-control primitives counted as free, or
unknown costs treated as zero.

## Actionable bounded next proposal

The useful research question is whether the explicit Boolean orbit-address
structure and sector controls can be synthesized jointly with cancellations,
instead of generic transposition compilation. A future algebraic round could
choose one fixed P above, derive its algebraic-normal-form Boolean mapping,
then synthesize a finite specified reversible template together with its
conditional Fourier gates. It must give native chronological gates and
compiled CX counts before comparing with13. This route differs from the two
current whole-circuit topology mutations and prior fixed-prefix eigenspace
Procrustes fits; it does not assert those older orbit studies were exhaustive.
The coarse bound supplies no credible expectation of beating13 by itself.

Adversarial review: gate model and abstract-vs-native distinction explicit;
right-cycle direction checked from binary strings; Fourier sign checked by
column shift; conditionalH branch phase retained; negative controls handled
through projector signs. No confirmed defect. Remaining challenge is actual
efficient synthesis, not the abstract diagonalization identity. Report/chat
files are intentional kept evidence; coordinator owns shared commit.

Verifier independently confirms each cost sum: four-bit phase34, conditionalH
18/58 for two/three controls, quartet70, threequartets210 and full loose3838.
Coordinator's separate Clifford obstruction agrees with this construction's
non-Clifford controlledS/nonlinear-encoder content; no Clifford-only completion
or new CX lower bound is claimed. The special00address sector explicitly uses
I2 andH2, not uniformF4 on allfour sectors.
