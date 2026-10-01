# Five-CNOT parity prefixes at the first matchgate pair

The exact 14-CNOT circuit has a six-CNOT phase-polynomial prefix with final
linear rows `(8,5,10,1)`, followed by the pair blocks `02,13,01,23`.
The first two blocks have Cartan operator

`F = exp(i*pi*XX/4 + i*pi*YY/8)`.

This screen asks whether a five-CNOT parity prefix can replace the six-CNOT
prefix while leaving only one CNOT on `02` or `13` to absorb into the first
two-qubit block. For an output CNOT `c->t`, its required final rows are the
exact prefix rows with row `t` XORed by row `c`:

| CNOT | Required rows |
|---|---|
| `0->2` | `(8,5,2,1)` |
| `2->0` | `(2,5,10,1)` |
| `1->3` | `(8,5,10,4)` |
| `3->1` | `(8,4,10,1)` |

Runs 194 and 196 exhaust all `12^5 = 248832` directed five-CNOT schedules.
For each path, a newly produced nonlocal row can receive a phase rotation.
The path must visit every nonlocal mask in one of the repository's 54
five-mask or eight four-mask rational-`pi/8` parity models. None of the paths
reaching the four required endpoints visits the masks of any of these 62
models. The five-mask scan has 3,584 support-complete paths and 686 distinct
endpoints, but zero at the four required endpoints. The four-mask scan also
has zero such paths. These are exact finite enumeration results within the
catalogued phase-polynomial ansatz, not a lower bound on general 13-CNOT
diagonalizers.

The first-block absorption is possible in two of the hypothetical cases.
Let `L02` and `L13` be the exact input `ZYZ` frames from
`topology14_exact_matchgate_derivation.md`. For each pair, form
`B = F L_pair CX`, where `CX` is in the indicated orientation. In the magic
basis `Q`, normalize `B` to determinant one and set
`M = (Q^dagger B Q)^T (Q^dagger B Q)`. Exact arithmetic in
`Q(zeta_96)` gives:

| CNOT | `tr(M)` | Result |
|---|---|---|
| `0->2` | `-2*zeta_96^16` | nonreal; outside the two-CNOT Cartan face |
| `2->0` | `0` | `tr(M^2)=tr(M^3)=0`; Cartan `(pi/4,pi/8,0)` |
| `1->3` | `0` | `tr(M^2)=tr(M^3)=0`; Cartan `(pi/4,pi/8,0)` |
| `3->1` | `2*zeta_96^8-2*zeta_96^24` | nonreal; outside the two-CNOT Cartan face |

The determinant is one, so three vanishing power traces imply characteristic
polynomial `lambda^4+1`, exactly that of `F`. For a two-CNOT block the third
Cartan coordinate is zero, making `tr(M)` real. Thus the two nonreal cases
cannot be implemented by two CNOTs as isolated first blocks. The two
absorbable orientations would work at the two-qubit level, but the enumerated
five-CNOT parity prefixes cannot reach their required endpoints.

No numerical joint fit was run because this screen found no eligible prefix
in the stated ansatz. A distinct next test is to allow a phase mask to be
implemented inside the first two-qubit block, jointly fitting its local
frames and the five-CNOT prefix. That changes the phase-polynomial support
condition and should be logged as a separate bounded search.

Reproduce with:

```bash
python3 search13_firstpair_prefix_screen.py
python3 search13_firstpair_fourmask_screen.py
python3 experiment_log.py show 194
python3 experiment_log.py show 196
```
