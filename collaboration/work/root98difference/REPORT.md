# Difference-register route and phase-sharing hypothesis

Board111. By-hand coordinator derivation; no matrix evaluation or search.
Target conventions q0MSB, chronological gates, UV4=DU, arbitrary output order.
Bell coordinates A=(q0,q2),B=(q1,q3) satisfy
T|a,b>=(-1)^(b0 b1)|b,a>. R=CX01 CX23 maps to (a,d=a xor b), giving
R T R†|a,d>=(-1)^((a0 xor d0)(a1 xor d1))|a xor d,d>.
Thus T_d=CZ_A X0^d0 X2^d1, with X acting first.

A final inverse router is unnecessary: if C R T R† C† is diagonal, C R
already diagonalizes. Output permutations preserve diagonality and are free
choices of eigenvalue order, not free circuit gates.

The exact Bell prefix E_pair† uses CX02,H0,CX13,H1 chronologically: two
CNOTs. Its generic F compiler allocation is four; that allocation is not a
minimal cost for these specialized Clifford decoder parameters. Prior all98
search counts stay unchanged because generic F angles were free.

Let B=CS† on A, beta(a)=i^(-a0 a1). Direct amplitude multiplication gives
B T_d B† phase beta(a xor d)*sigma(a xor d)*conj(beta(a)).
For d10 this is i^a1 times X0; for d01 i^a0 times X2; for d11 it is i on
even A parity and1 on odd parity times X0X2; for d00 it remains CZ_A.
Let P=CX0->2 controlled on d1=q3. P leaves d01 unchanged; d11 becomes
i^(1-a1)X0. Then H0 conditioned on d0=q1=1 and H2 conditioned on
(q1,q3)=(0,1) diagonalize every branch. The selectors preserve d.

The literal primitive construction costs19CX: Bellprefix2 + R2 + B2 + P6
(Toffoli) + CH1 + doubly-controlledH6. It is a structural comparator,
not a new incumbent. B uses phase(-pi/4) on each A wire, CX02,
phase(+pi/4) on q2, CX02. CH is oneCX conjugated by target locals;
doubly-controlledH is Toffoli conjugated by target Ry(±pi/4).
All formulas were independently reviewed by algebraic/verifier collaborators.

Next hypothesis: relative-phase three-CX versions of both six-CX selectors
could yield13 before correction, but arbitrary relative phases can break
basis mixing. Characterize phases that commute to the output or merge withB;
one further CNOT saving would be needed for12. No such complete gate list
has been established. A separate ledger-frozen check is required before
computing even the literal19CX comparator. This does not exclude full98,
other ordered families, or any12CX diagonalizer.
