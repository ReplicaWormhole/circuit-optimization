# Endpoint-guided five-CNOT prefix search

Runs 247 and 249 sought a 13-CNOT diagonalizer of the four-qubit right
shift. The target exact 14-CNOT circuit has a six-CNOT prefix with binary
endpoint rows `(8,5,10,1)`, followed by eight CNOTs in two first-pair and
two final-pair matchgate blocks. A five-CNOT endpoint differing by one CNOT
on `02` or `13` may allow that correction to be absorbed into the first
pair operator. The existing exact Cartan screen finds that orientations
`2->0` and `1->3` have the appropriate isolated two-CNOT first-pair
invariant; this is a necessary structural clue, not a whole-circuit
construction.

Run 247 enumerated all `12^5 = 248,832` directed five-CNOT paths in binary
row arithmetic. Each of the two chosen endpoint classes has 209 paths.
After excluding paths within one edge change of a simple deletion of the
exact six-CNOT prefix, 414 paths remain. Of those, 293 visit none, 106
visit one, and 15 visit two of the exact prefix's nonlocal phase masks
`{6,7,11,14}`. None visits three or four. The endpoint counts agree with
the older independent enumeration in
`search13_firstpair_prefix_screen_result.json`. This finite count does not
exclude general local gates, changed tail interactions, or other endpoint
gauges.

Run 249 selected four paths with two visited masks, spanning both endpoint
classes and different visited-mask pairs. Their full local `SU(2)` layers
were optimized across both the five-CNOT prefix and the original
eight-CNOT tail interaction pattern. Each path used a phase-informed seed
with the exact tail's local frames and a transfer seed from run 230's
near-hit. All eight L-BFGS-B fits were capped at 350 iterations and failed
`check_circuit.py`.

The smallest off-diagonal Frobenius loss was `0.125`, but that candidate's
maximum off-diagonal matrix entry was `0.70710678`. The candidate with the
smallest maximum entry had loss `0.1417625585` and independently checked
maximum entry `0.2665569719`, essentially the same basin as run 230's
prefix-deletion near-hit. Neither is a valid 13-CNOT diagonalizer. No exact
gate check was warranted because no candidate passed numerical validation.

The failure is limited to the four selected paths, two initializations
each, the original eight-CNOT tail interactions, and the iteration cap.
Every selected path omits two of the exact nonlocal masks. A revised
hypothesis is to reroute one **cross-pair** CNOT inside the eight-CNOT tail
while retaining an endpoint-guided five-CNOT prefix, then optimize all
local frames together. This lets the changed tail transport an omitted
phase across the first/final-pair boundary, which the fixed tail of run
249 could not do by its CNOT pattern alone. The revised search should be
screened against cut-rank capacities before fitting.

Reproduce:

```bash
python3 search13_endpoint_guided_screen.py
python3 search13_endpoint_guided_fit.py --seed 25050 --maxiter 350 --phase-sigma 0.04 --warm-sigma 0.08
python3 check_circuit.py search13_endpoint_guided_case1_warm.json
python3 experiment_log.py show 247
python3 experiment_log.py show 249
```

## Cross-pair tail reroutes

Run 255 built the fixed residual target `U14 P5^{-1}` for the four
phase-informed five-CNOT prefixes of run 249 and computed operator
Schmidt ranks exactly over `Q(zeta_96)`. On cuts `(0,1,2,3,01,02,03)`,
every residual has ranks `(4,4,4,4,16,16,16)`. All 32 tail schedules
formed by replacing one of the original eight CNOTs with directed `03`
or `12` have sufficient cut capacities. Thus this exact screen excludes
none of them. It concerns the fixed target `U14 P5^{-1}` only; a
different output eigenbasis or jointly varied prefix changes the
operator.

Run 257 fitted four boundary-adjacent reroutes using endpoint-guided
five-CNOT prefixes. The replacement was the final first-pair CNOT
`13 -> 12` or the initial final-pair CNOT `01 -> 03`. Two warm-start
perturbations were used per route. All single-qubit gates across the
13-CNOT circuits varied under the label-free cycle loss, with a cap of
300 L-BFGS-B iterations per fit. None of the eight candidates passed
`check_circuit.py`. The lowest Frobenius loss was `0.125` with maximum
off-diagonal entry `0.7071`; the smallest independently checked
maximum entry was `0.4715229142` with loss `0.2355171485`. This is worse
than the unchanged-tail endpoint near-hit at `0.2665569719`. No exact
circuit check was run on these invalid numerical candidates.

The exact whole-tail cut ranks are saturated and cannot distinguish the
routes tested here. A revised test should place a synthesis boundary
*after* the altered first-pair operation and compute exact residual
four-qubit ranks or two-qubit Cartan data there. That would screen whether
the omitted prefix phase can be absorbed without leaving an over-complex
final layer before spending more optimizer time. Any such certificate
would apply to its stated fixed intermediate operator, not to arbitrary
jointly synthesized circuits.

```bash
python3 search13_endpoint_guided_reroute_rank.py
python3 search13_endpoint_guided_reroute_fit.py --seed 25060 --maxiter 300
python3 check_circuit.py search13_endpoint_guided_reroute_case0_s0.json
python3 experiment_log.py show 255
python3 experiment_log.py show 257
```
