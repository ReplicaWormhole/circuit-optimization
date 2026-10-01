# Three-edge tree reachability and one six-CX survivor

Hypothesis64; structural derivation and read-only inspection, no enumeration,
symbolic elimination, numerical screening or certification. Four wires, q0
MSB, arbitrary local gates. Operator evolution direction is explicit: W Zj W†,
where W=U† for diagonalizer convention U V U†=D. Historical skeleton notes
use W's forward chronological direction; do not silently reverse that schedule.

## Fresh-leaf tree lemma

Consider exactly three distinct CX edges forming a four-wire tree, with
arbitrary local gates between them, and one initially single-wire Pauli axis.
For a final full-support image, all three edges must be active and each add
one fresh wire to the support. Otherwise fewer than four wires can enter
the causal cone. Thus this is a rooted growing-tree operator path, not an
assumption that every ordering gives full support.

At each active gate the observable, if product up to that point, has a Pauli
axis on the active endpoint and I on the fresh endpoint. A control axis
aX+bY+cZ becomes (aX+bY) tensor X+cZ tensor I. It has product form iff
c=0 (spreading XY plane) or a=b=0 (commuting Z, no fresh support). A target
axis aX+bY+cZ becomes I tensor aX+Z tensor(bY+cZ), with product form iff
a=0 (spreading YZ plane) or b=c=0 (commuting X, no fresh support). The
full-support requirement rules out the commuting case at each tree edge.

If both components are nonzero, operator Schmidt rank across the cut formed
by deleting that tree edge is2. Every later tree edge lies entirely on one
side of that cut, so subsequent gates and local rotations cannot reduce the
rank. Therefore final full product and full support are equivalent to choosing
the spreading-plane condition at every fresh-leaf crossing, with the actual
preceding local rotations included. This is an exact template condition,
not merely existence of an abstract commuting operator on full support.

Only axes originating on the first edge's endpoints can have full support
after three gates; a later starting wire misses an already used edge.
Hence at most two one-axis-per-input-wire images can be full-support products.
Moreover, after the first edge both candidate images have product support
on that edge's two old wires. Their factors need not match there: controlX
and targetZ, for example, give XcXt and ZcZt. The SECOND tree edge must
attach a common fresh third wire to that same old component for both images
to reach full support. At this second edge, the new wire receives the same
Pauli X (old endpoint control) or Z (old endpoint target) in both images,
independent of the old spreading-plane axis. Common local rotations preserve
proportionality of these new-wire factors. If the third edge touches that
wire, it attaches the final fresh wire; product/full-support assumptions force
both proportional axes to spread, preserving proportionality on the old
wire as well. If it does not touch the shared wire, only local rotations
act there. Thus both final images share a proportional factor. If both
products commute with V, the independently reviewed cyclic product lemma
makes them ±Q^tensor4 and ±R^tensor4. Sharing that fresh-wire axis forces
Q,R parallel, making the images dependent, impossible for two distinct input
Pauli axes. Thus at most one independent full-support cycle-commuting product
image is attainable in this three-edge tree template. One is attainable by
spreading a root axis through a tree and matching final local axes. This
does not say every tree image is product or that later gates preserve this cap.

## One read-only degree4422 example

The first degree4422 representative inspected in
global_lower_bound7_dihedral_orbits_result.json is
[(0,1),(0,2),(2,3),(1,3),(0,1),(0,1)], with degrees(4,4,2,2).
This is a representative from the archived144class, whose count remains144
in the current necessary screen; no new schedule enumeration was performed.
Its first three edges form a growing tree from wire0, so the above lemma
applies to that prefix. The last three edges include a second route and
repeated01, allowing new entanglement or rank reduction across cuts previously
protected in the pure tree. Therefore the prefix lemma cannot force final
six-gate product images. No survivor is excluded.

For degree-two wire2, an axis commuting through its first02 gate is still
one-body before23 after appropriate local propagation. Its remaining active
path is23,13,01,01, with the last pair repeated. It is not a no-revisit tree.
This identifies an actual small template to study next, not a commutant-only
support system. Any elimination needs shared gate/local parameters and its
own ledger reservation. The global interval remains6<=Cmin<=13.

## Projected AB curvature and prefix gauge review

A24-coordinate A/B Hessian with the remaining92coordinates fixed tests only
that subspace. A sufficiently negative direction checked by independent
finite differences gives a numerical descent direction there; it cannot
establish global optimum or a topology exclusion. A PSD projected Hessian
does not establish full116PSD or a local minimum. Run387's reported least
value near-8.55e-17 failed its frozen-1e-6 threshold; no escape fits triggered.

Independently reviewed root PREFIX_GAUGES.md: A/B factors on control wires
0,2,3 move to free boundary layers without crossing their own controlled
gate, preserving chronological same-wire multiplication. Only two targetq1
gates need remain internally, six rather than24 extra family coordinates.
This is exact family equivalence with boundaries free, not permission to
drop24-subspace coordinates while freezing those boundaries. No gauge nullity
or Hessian dimension is inferred from the family representation alone.

Adversarial review: growing-tree and full-support premises stated; protecting
cut argument uses unique tree edge; six-template revisit explicitly prevents
extrapolation; projected curvature and gauge scope separated. Remaining risk
is absence of a six-gate forcing theorem. No six-CX exclusion, incumbent
promotion or duplicated certificate. Reports/chats retained evidence.

Verifier independently accepted the stronger at-most-one tree lemma under
its full-support/product/cycle-commuting premises. Attainability is an example
with a suitable fanout tree and locals, not a claim for every fixed ordered
tree. No six-CX forcing or exclusion was obtained. The free-boundary prefix
representation has98coordinates (original92 plus six target-wire angles);
the executed fixed-boundary24audit remains unchanged.

Closeout correction: an earlier wording incorrectly located the common
factor on the first tree edge. That is false for controlX/targetZ. The proof
above now locates it on the fresh wire introduced by the second edge. The
lemma statement is unchanged; no search or exclusion followed the defect.
