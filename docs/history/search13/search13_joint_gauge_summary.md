# Paired eigenspace gauges with a jointly fitted seven-CNOT tail

Runs 225 and 227 tested a restricted 13-CNOT construction. The first six
CNOTs and the first 13 gates of the exact 14-CNOT diagonalizer were fixed.
The remaining target was the entire eight-CNOT tail, so the numerical fit
crossed both first-pair and final-pair blocks.

Run 225 rotated output rows `(8,11)` and `(1,14)` by real Givens angles
`pi/4`. Each pair lies in a single degenerate eigenspace, of eigenvalue
`+1` and `-i`, respectively. The gauged full tail had numerical balanced-cut
Schmidt ranks `(16,15,16)` for cuts `01|23`, `02|13`, and `03|12`, using a
singular-value tolerance of `1e-10`. Both tested interleaved seven-CNOT
schedules passed this necessary rank screen. Each schedule was fitted with
arbitrary local `SU(2)` layers before, between, and after its CNOTs, using
one specified random seed and at most 350 L-BFGS-B iterations. The better
process-overlap loss was `0.4682348383`, but `check_circuit.py` measured
maximum cycle off-diagonal error `0.6556081352`. The candidate is invalid.

Run 227 continued that better seven-CNOT schedule from its fitted local
gates, changing the objective to the label-free cycle off-diagonal
Frobenius loss. It converged after 75 iterations to loss `0.515625`, with
`check_circuit.py` maximum off-diagonal error `0.7071067879`. The second
candidate is also invalid. These are numerical failures for the stated
gauges, schedules, and initializations, not a lower bound on 13 CNOTs.

Code, complete gate lists, and machine-readable records are in
`search13_joint_gauge_fit.py`, `search13_joint_gauge_continuation.py`, and
their corresponding `search13_joint_gauge_*.json` files. The ledger archived
the better candidates under runs 225 and 227. No exact check was run on
either candidate because neither passes the numerical diagonalization check.

Reproduction commands, from this directory:

```bash
python3 search13_joint_gauge_fit.py --maxiter 350 --seed 25010
python3 check_circuit.py search13_joint_gauge_candidate_1.json
python3 search13_joint_gauge_continuation.py --maxiter 500
python3 check_circuit.py search13_joint_gauge_continuation_candidate.json
python3 experiment_log.py show 225
python3 experiment_log.py show 227
```

A distinct next hypothesis is to include complex phases in the two
eigenspace rotations and screen whether a pair lowers two balanced ranks at
once. Such a gauge would change the compatible seven-CNOT interaction
multisets before numerical synthesis. The present real rotations lowered
only one rank and both interleaved fits remained far from the target.

## Complex-phased follow-up

Run 232 evaluated 17,472 ordered specifications of two disjoint
same-eigenvalue row rotations. Each angle was `pi/8` or `pi/4`; each phase
was one of `-pi/2`, `-pi/4`, `pi/4`, and `pi/2`. It used numerical SVD with
tolerance `1e-9`. No sampled gauge passed the necessary cut-rank graph
screen for a three-CNOT **fixed final layer**. This finite-grid numerical
screen is not a theorem about all phases or changed synthesis boundaries.

The most useful full-tail gauge rotates pairs `(4,7)` and `(8,11)` by
`pi/4`, both with phase `-pi/2`. Run 233 confirmed over `Q(zeta_96)` that
its full-tail balanced-cut Schmidt ranks are exactly `(14,12,14)` on cuts
`01|23`, `02|13`, and `03|12`. Its fixed final-layer ranks are exactly
`(2,12,12)`. The first is compatible with some seven-CNOT interaction
graphs. The second cannot have a three-CNOT realization because a
three-CNOT circuit has rank at most `8` across any cut. This exact
obstruction applies to this fixed gauged operator only.

Run 233 also certified a nearby phased gauge with phases `(-pi/4,pi/4)`:
its full-tail ranks are `(16,11,16)` and its fixed-final ranks are
`(3,11,11)`. Its fixed final layer likewise cannot use three CNOTs.

Run 235 fitted the first gauge to the seven-CNOT chronological tail
`(02,01,03,01,12,03,01)`, with arbitrary local `SU(2)` gates. The fixed
target process fit reached loss `0.6183999`; an independent
`check_circuit.py` invocation measured maximum cycle off-diagonal error
`0.6027011`, so its 13-CNOT circuit is invalid. Gauge-free continuation
converged to off-diagonal Frobenius loss `0.5544733`, but maximum entry
error rose to `0.8300947`; that circuit is invalid too. These are bounded
numerical failures, not a lower bound.

Follow-up reproduction commands:

```bash
python3 search13_joint_gauge_phased_rank.py --top 20
python3 search13_joint_gauge_phased_exact.py
python3 search13_joint_gauge_phased_fit.py --seed 25020 --process-maxiter 350 --direct-maxiter 500
python3 check_circuit.py search13_joint_gauge_phased_process_candidate.json
python3 check_circuit.py search13_joint_gauge_phased_direct_candidate.json
```

## Gauge-free interleaved tails

Run 241 changed the topology rather than prescribing an output gauge. It
tested three chronological seven-CNOT full tails after the exact six-CNOT
prefix. Every schedule used at least five undirected interaction pairs,
interleaved repeated pairs, and passed the necessary single-wire and
balanced-cut capacity test for the ungauged full tail, whose ranks are
`(4,4,4,4,16,16,16)` on cuts `(0,1,2,3,01,02,03)`.

Each schedule was optimized directly for the label-free cycle
off-diagonal Frobenius loss from two starts: a transfer warm start from
run 235's seven-CNOT local frames and a seeded random start. The best
local frames were transferred between successive schedules. Six fits used
at most 350 L-BFGS-B iterations each. The strongest result used the
chronological tail `(02,01,23,02,12,01,03)` and reached loss
`0.2235680754`; independent `check_circuit.py` evaluation found maximum
off-diagonal entry `0.4658954362`. It is an invalid 13-CNOT candidate.
The other two schedules reached best losses `0.3034` and `0.4375`.
These bounded numerical failures do not exclude any schedule globally,
or any other 13-CNOT topology.

```bash
python3 search13_joint_gauge_interleaved_direct.py --seed 25030 --maxiter 350 --sigma 0.4
python3 check_circuit.py search13_joint_gauge_interleaved_s0_transferred.json
python3 experiment_log.py show 241
```

## Earlier boundary: five-CNOT full-parity endpoints

Run 244 replaced the exact six-CNOT prefix by two wholesale five-CNOT
linear networks, each at least two edge changes from every simple
one-CNOT deletion of the exact prefix. Their endpoints were
`(15,5,14,9)` and `(15,12,11,5)`, respectively. Thus each exposes the
full-parity mask `15` on one wire, where a local phase can act as a
nonlocal input gauge. The first kept the original eight-CNOT tail
interaction list; the second interleaved its disjoint first and final
pairs. All local `SU(2)` layers across both prefix and tail were fitted
jointly with the label-free cycle loss.

One warm start transferred the complete local frames of run 230's
near-hit, and one start used seeded random angles, per topology. With at
most 300 L-BFGS-B iterations per fit, all four candidates failed the
diagonalization checker. The best Frobenius loss was `0.375`; the smallest
independently checked maximum off-diagonal entry was `0.463751985` on a
different candidate. That latter optimization stopped at its iteration
limit. Full-parity exposure alone did not improve the run-230 near-hit,
whose maximum off-diagonal error was about `0.267`.

These endpoints are also far from the exact six-CNOT prefix endpoint
`(8,5,10,1)`: the first preserves only one row and the second preserves
none. A revised search should enumerate five-CNOT paths whose endpoint is
one CNOT from `(8,5,10,1)` while visiting a useful nonlocal phase mask,
then jointly fit the compensating first-pair operator and eight-CNOT
tail. This tests a structured phase-transport hypothesis instead of
selecting full-parity exposure alone. Run 244 is a bounded numerical
failure, not an exact obstruction.

```bash
python3 search13_joint_gauge_five_endpoint.py --seed 25040 --maxiter 300 --warm-sigma 0.05 --random-sigma 0.4
python3 check_circuit.py search13_joint_gauge_five_endpoint_1_random.json
python3 experiment_log.py show 244
```
