# Two-shear orbit coordinates and branch Fourier bases

Write a computational basis state as `x0x1x2x3`, with `x0` most significant.
The right cycle `V4` maps `(x0,x1,x2,x3)` to `(x3,x0,x1,x2)`.

Run 108's best reversible shears are

`S1(x) = x ⊕ 1111 · (x2⊕x3)(x1⊕x3)` and

`S2(x) = x ⊕ 0011 · x1(x0⊕x2⊕x3)`.

The products and sums here are in `F2`; the masks in the run108 file are
`S1=(15,3,5,0,0)` and `S2=(3,4,11,0,0)`. After `S2 S1`, the three
length-four orbits are cosets of `H={0000,0110,1010,1100}`. The remaining
four states are the fourth coset. `H` is characterized by `x3=0` and
`x0⊕x1⊕x2=0`.

Apply `CX(0,2), CX(1,2)` after the shears. The resulting bits `y` have
position register `p=(y0,y1)` and orbit-label register `ℓ=(y2,y3)`. The
conjugated right cycle preserves `ℓ` and acts on the four position values
`00,01,10,11` as follows:

| Label `ℓ` | Position permutation | Interpretation | Exact diagonalizing basis |
| --- | --- | --- | --- |
| `00` | `[0,2,1,3]` | `SWAP(p0,p1)` | chronological `CX(0,1), H(0)` |
| `01` | `[1,2,3,0]` | increment modulo 4 | chronological `H(0), CP(+π/2), H(1)` |
| `10` | `[3,0,1,2]` | decrement modulo 4 | same Fourier basis |
| `11` | `[1,3,0,2]` | Gray-code increment | chronological `CX(0,1)` then Fourier basis |

Here `CP(+π/2)` is exact with two CNOTs and rational-π `Rz` gates.
One chronological sequence, up to global phase, is
`CX(0,1), Rz_1(-π/4), CX(0,1), Rz_0(π/4), Rz_1(π/4)`.
The Fourier basis gives eigenvalues `(1,-1,i,-i)` for increment and
`(1,-1,-i,i)` for decrement. The Bell basis gives `(1,1,1,-1)` for SWAP.
These identities were checked exactly over `Q(i,√2)` in run 116.

A label-preserving direct sum of these four branch bases diagonalizes the
conjugated `V4` exactly (run 119). Its operator-Schmidt ranks are `4,4,2,2`
for one-wire cuts `0,1,2,3` and `3,7,7` for two-wire cuts `01,02,03`.
All-cut edge enumeration gives a four-CNOT lower bound for **this chosen
block**. This bound does not account for a different eigenbasis or a joint
synthesis with the nonlinear coordinate map.

The conventional exact synthesis of the two shears alone took 21–22 CNOTs
in the algebraic agent's run 115, before the two linear-postmap CNOTs.
Consequently that direct implementation cannot give a 14-CNOT full circuit.
For separate synthesis using the chosen block, a different preprocessor
would need at most ten CNOTs including the linear postmap. The algebraic
agent's run 118 tested relative-phase versions of all 128 aligned two-shear
maps. Its best compiled preprocessors still required at least 16 CNOTs,
and none had input phase defects constant on every `V4` orbit. This closes
the tested separate two-shear architecture at the 14-CNOT target. A joint
synthesis or a different coordinate map is outside that conclusion.

Reproduce the finite-field action with `python3 archive/legacy_root/topology14_orbit_action.py`,
the exact branch identities with `python3 archive/legacy_root/topology14_orbit_branch_basis.py`,
and the direct-sum rank certificate with
`python3 archive/legacy_root/topology14_orbit_block_rank.py` (ledger runs 112, 116, 119).
