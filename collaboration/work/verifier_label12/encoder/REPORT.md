# Independent exact encoder verification

Hypothesis47. No independent search. Qubit0 MSB, chronological reversible
permutation P|x>=|p(x)>. Direct orbit-address derivation gives table
[0,4,7,8,6,2,11,12,5,9,3,13,10,14,15,1], independently matched to frozen
numerical config. Before reservation/search, sourcehash, configSHA
5e28c0afff49bac0ac84a6adc8b9aa85172cb662c78e87c6acc9a14c48d28afc,
all108catalog16basis actions and conservative CXbounds0/1/6/34 were verified
and directly sent to numerical researcher with board handoff.

Algebraic44three-control-X sequence independently evaluates exactly to fixedP
on16inputs; all64 outputbits match stated ANF masks. Cost44*34=1496CXupper.
No optimized synthesis claim. Check_encoder.py records hash-bound evidence.

Run381 result independently composes all saved native actions on every16basis
input, reproducing fixedP exactly with zero mismatch:10steps7CX+3CCX,
compiledCXupper7+3*6=25, noCCCX; negativecontrols are exact localX
conjugations. Search stopped exactfound at122364expansions/depth10 under frozen
beam128/depth12/200000expansion/120second/one-thread limits. This is a finite
heuristic discovery, no minimum proof. Full result hash and composed action
are retained in beam_check.json. Native basis-state integer verification is
exact for these permutation gates; it does not certify an unspecified
compiledCX list or a complete cycle diagonalizer.

Adding previously reviewed generic conditionalFourier upper268 gives a
constructive total293CXupper, which is far above accepted13. Encoder improvement
is separate from diagonalizer improvement. No nativeF/CX cost conflation,
unknowncost treatedzero, exact13 recertification or promotion.

Adversarial review found no confirmed issue in bitordering, chronology,
negativecontrols, ANF or cost upperbound semantics. Independent check scripts
and integer audits pass; limited beam gives no lowerbound/impossibility.
Future work should simplify the explicit encoder together with Fourier
controls under a newly stated budget, rather than claim25CX is wholecircuit.
Intentional scripts/JSON/report/chats retained; coordinator owns commit,
no push.

Board49 assigned focused compiler followup: implemented phase_compiler.py.
Projector phase coefficients are exact Fraction rational-pi; cyclic Gray
parities group by highest target, visit all nonempty subsets, restore target
and cost2^m−2CX. Every emitted phase polynomial is checked exactly on all
conditioned-bit basis words against theta*(projector−1/2^m). Dropped term is
input-independent global phase across all branches. Negative controls change
coefficient signs only. ControlledH chronologicalRy(−pi/4), conditionedZ,
Ry(pi/4) uses H=Ry(pi/4)ZRy(−pi/4); unrelated branches cancel rotations.

Complete orbit_encoder67.json compiles381encoder25CX, global bit-reversedF4
2CX, upper00 sameF4 conditionalinverse26CX, and001pairH14CX: total67CX.
All local rotation fields exact rational-pi strings. No cancellation pass
performed. Full numeric checker passes maxoff5.79e-16, unitarity1.55e-15;
focused exactphase checks are not a full exactUV certificate. That execution
is explicitly deferred to separately reserved rootbudget. Sourcehash bound,
no67CX incumbent promotion since accepted13 already better. Frozen381catalog
CCCX34 unchanged; newlyprovedcompiler14 available only future definitions.
