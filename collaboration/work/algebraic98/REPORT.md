# Exact 116-to-98 prefix absorption map

Board hypothesis 66 (claimed by `algebraic98`): derive the exact boundary
remapping for the relaxed prefix and state precisely what a full 98-coordinate
diagnostic covers. This was a read-only structural derivation; no search,
numerical evaluation, certificate attempt, or gate-list edit was performed.

## Family and coordinate count

The relaxed chronology is

`L0 CX01 A CX21 B CX31 L1 F02 L2 F13 L3 CX12 L4 F01 L5 F23 L6`,

where each `Lk`, `A`, and `B` is a four-qubit product of arbitrary SU(2)
gates, and `Fct(a,b)=exp(i a XcXt + i b YcYt)`. Qubit zero is most
significant and the gates are chronological. The nine four-qubit local layers
contribute 108 real rotation-vector coordinates, and the four F blocks
contribute eight, for 116 total. For each single-qubit triple `(x,y,z)`, the
optimizer uses `exp[-i(x X + y Y + z Z)/2]`. The vector is listed in local
layer order `L0,A,B,L1,...,L6`, with qubits 0,1,2,3 inside each layer, followed
by the four F angle pairs. When serialized as circuit JSON, each SU(2) matrix
is converted to a `u3(theta,phi,lam)` representation up to scalar phase; those
Euler angles are not the optimizer's coordinate chart.

Only `A` and `B` factors on the common target q1 need to stay in the interior.
Every other interior factor is on a control wire of the three-CNOT group and
can be moved to one of the free boundary layers without crossing a CNOT on
which that wire is a control. The required same-wire multiplication order is:

- q0: move `A0` and `B0` forward across `CX21,CX31`; absorb as
  `L1'(0) = L1(0) B0 A0`.
- q2: move `A2` backward across `CX01` and `B2` forward across `CX31`;
  absorb as `L0'(2) = A2 L0(2)` and `L1'(2) = L1(2) B2`.
- q3: move `A3` and `B3` backward across `CX01,CX21`; absorb, preserving
  chronological order, as `L0'(3) = B3 A3 L0(3)`.
- q1: retain `A1` and `B1` in the interior, each with three SU(2) coordinates.

All other factors in `L0,L1` are unchanged. The later layers `L2,...,L6` and
all F angles are unchanged. Thus the transformed chronology is

`L0' CX01 A1 CX21 B1 CX31 L1' F02 L2 F13 L3 CX12 L4 F01 L5 F23 L6`.

It has seven arbitrary local layers (84 coordinates), six target-wire
interior coordinates, and eight F coordinates: 98 total. This is equality of
realizable circuit families: the stated commutations map every 116-coordinate
instance to this form; conversely, each 98-coordinate instance is obtained in
the original form by setting the control-wire components of `A,B` to identity.
The map uses arbitrary free boundary SU(2) factors. Any scalar phase carried
by a U(2) gate convention contributes only an overall circuit phase and does
not change the diagonalization objective. This is not a reduction in entangler
count or proof of parameter-manifold dimension.

## Embedding and diagnostic scope

For the 385/386 warm start, `A=B=I` embeds without changing any boundary
matrices; the six retained q1 interior coordinates are zero. A diagnostic
must still include all 98 coordinates to probe coupled boundary, target
interior, later-local, and F directions, including mixed second derivatives.
The 116-to-98 map is not generally obtained by deleting 18 coordinates while
holding `L0,L1` fixed: for nonidentity controls, the four displayed boundary
products must first be formed.

The independent verifier reports a generic nonidentity random-matrix check of
the full remapping with maximum matrix discrepancy `5.36e-16` up to global
phase, and an exact match for the frozen zero-A/B embedding. This supports the
implemented map and that particular embedding; it is not a certificate for a
diagonalizing circuit.

The products above are matrix products. Rotation vectors must not be added
componentwise either; if a 98-vector is required after a nonidentity
remapping, compose the SU(2) matrices and take coordinates in the explicitly
chosen chart (for example, a rotation-vector logarithm), recording branch and
scalar-phase handling. A full Hessian in one specified 98-coordinate
chart tests curvature only at that represented point in that chart. A negative
eigenvalue supplies a local descent direction; a nonnegative Hessian would
not certify global family optimality or a 12-CNOT obstruction. The prior
fixed-boundary 24-coordinate test did not include these coupled directions.
Run 389's full chart has gradient norm `1.60e-8` and least computed eigenvalue
`-2.70e-9`, far above its frozen negative-curvature trigger `-1e-6`; 57
eigenvalues lie within `1e-6` of zero. These floating values do not resolve
curvature in the flat cluster or prove a local minimum or saddle. The inherited
seed gate list remains invalid (maximum off-diagonal entry `0.9439`). No
perturbation fit was authorized by that diagnostic.

## Recommended next distinct hypothesis (not dispatched)

Test whether a fresh basin can reach a diagonalizer in the same coupled 98
family: freeze one reproducible, genuinely 98-coordinate initialization with
nonzero `A1` and `B1` (not an embedding of run 385), then run one separately
reserved, single-threaded, externally capped refinement with no restarts. A
material loss reduction or a passing candidate would show the inherited
high-loss point was not representative of the family; failure would remain a
bounded optimizer result, not a family exclusion. Independently validate any
complete gate list, and require the exact certificate process before changing
the 13-CNOT incumbent. The run budget, seed, chart, code hash, and stop rule
must be frozen in the board and ledger before execution. Do not repeat the
same-point 24- or 98-coordinate Hessian.

## Actionable handoffs

The exact map, coordinate ordering, and zero-A/B embedding were sent to
`numerical98` and `verifier98` via direct messages and persisted board
handoffs. Independent full-matrix remapping and seed-embedding checks are
reported above. No exact diagonalization certificate was attempted here.
