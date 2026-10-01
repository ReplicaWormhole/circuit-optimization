# A phase-compatible 13-CX selector construction

Board 123; by-hand structural derivation only, pending independent complete
gate-list numerical and exact certification. No matrices or search ran.
Use the notation and target conventions in `LITERAL_OBSTRUCTIONS.md`.

## Corrected literal 13-CX obstruction

Changing the anchor to either `B=CSdg_xy` or `B=CS_xy` fixes the `d=10`
quarter-turn error but does not fix the unconjugated RCCX P selector.
After B, the `d=11` operator is

`K_11=f_r(x xor y) X_x X_y`, where `f_r(0)=r`, `f_r(1)=1`,
and `r=i` for CSdg or `r=-i` for CS.

On `v=1`, the P target blocks are `K0=Z,K1=epsilon Y` for controls `(v,x)`,
or `K0=I,K1=epsilon Y` for `(x,v)`, where `epsilon=+/-1` covers reversing
the target T signs. Direct block multiplication gives

`P XX Pdagger = -epsilon Y_x` or `-epsilon Y_x Z_y`, respectively.

P routes the original parity to y, so `f_r(x xor y)` becomes `f_r(y)`.
The subsequent H_x preserves Y_x up to sign. On d=11 the final CCH is
identity or a target-only operation on y, hence cannot remove the Y_x
off-diagonal factor. All these corrected literal 13-CX variants fail too.
The independent verifier confirmed the same factor obstruction by hand.

## Repair the first selector

Conjugate the P RCCX target by Sdg: chronological `S_y,RCCX,Sdg_y`.
This maps its active Y block to X and leaves the inactive Z block unchanged.
Use positive target T signs and controls `(v,x)` for definiteness. Its v=0
branch is identity, while its v=1 branch has target blocks `Z` at x=0 and
`X` at x=1. After B and this P, the four sector operators are

| d | operator |
|---|---|
| 00 | CZ_xy |
| 10 | r^y X_x |
| 01 | -(-r)^x X_y |
| 11 | f_r(y) Z_y X_x |

For d=01, Z conjugates X to -X while X conjugates X to itself.
For d=11, direct target blocks give `P XX Pdagger=X_x Z_y` and the routed
parity phase is `f_r(y)`. These checks determine the displayed signs.
Apply the exact controlled H_x when u=1. The d=10 and d=11 operators now
contain Z_x and are diagonal; d=01 still needs H_y.

## Repair the final inactive branch

The earlier conjugation mapped inactive Z to `(Z-X)/sqrt(2)`, which mixes
computational states in an inactive sector. Instead choose the target unitary

`Vprime=Ry(-3*pi/4) Rx(-pi/2)`.

It satisfies `Vprime Y Vprime_dagger=H=(X+Z)/sqrt(2)` and
`Vprime Z Vprime_dagger=Y`. First Rx maps Z to Y and Y to -Z;
then Ry(-3*pi/4) maps -Z to H and leaves Y fixed.

Use final RCCX controls `(c,v)` with `c=1-u`, target y, wrapped by Vprime.
Implement the negative control with X_u before and after the wrapper. The
target wrapper chronology is

`Ry_y(3*pi/4),Rx_y(pi/2),RCCX(c,v;y),Rx_y(-pi/2),Ry_y(-3*pi/4)`.

Its sectors are Y_y on d=00, H_y on d=01, and identity on d=10 and d=11.
On d=00, Y_y is a computational permutation with phases, so it preserves
diagonality: `Y_y CZ_xy Y_y = Z_x CZ_xy`. On d=01, H_y maps X_y to Z_y.
Thus every displayed sector becomes diagonal. The final inactive error is
harmless because it is a permutation in an already diagonal sector, rather
than because all inactive phases are automatically harmless.

## Complete block-level proposal and costs

Chronological sequence:

1. Exact Bell decoder `CX02,H0,CX13,H1`: 2 CX.
2. Difference router `CX01,CX23`: 2 CX.
3. B=CSdg_xy: 2 CX. One explicit chronology is phase(-pi/4) on x and y,
   CXxy, phase(+pi/4) on y, CXxy.
4. `S_y,RCCX(v,x;y),Sdg_y`: 3 CX, with the literal primitive word above.
5. Exact CH(u;x): 1 CX. Chronology is
   `Ry_x(-pi/4),H_x,CX(u,x),H_x,Ry_x(pi/4)`.
6. Negative-control Vprime-wrapped RCCX(c,v;y): 3 CX, as above.

Total: `2+2+2+3+1+3=13 CX`. Every angle is rational pi; all phases can use
`u3(0,0,theta)` in the repository JSON format. This is a structurally derived
13-CX comparator, not an accepted incumbent or an improvement over 13.

Swapping P controls to `(x,v)` also appears valid by the same sector method:
v=0 is CZ_xy and v=1 is the exact CXxy. Its d=10 coefficient changes from
r^y to (-r)^y, and d=11 loses Z_y; all sectors remain diagonal. Use the
definiteness choice above for a first bounded checker invocation.

## Remaining 12-CX issue and review

The d=10 sector still needs a quarter-turn anchor: replacing B by CZ again
fails there, since the corrected P is identity on v=0. A further CX saving
must come from a different factorization or an actual cross-block identity,
not from the previously failed CZ replacement. No such identity is proved.

Adversarial review corrected an intermediate parity-sign slip using direct
2x2 target blocks. The final table uses those direct blocks. Chronological
Vprime inverses, negative control, inactive branch, and honest 13-CX cost
were checked by hand. Independent full-list numerical and exact checks are
still required; no symbolic or numerical matrix was evaluated here.
