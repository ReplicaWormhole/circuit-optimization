# Restricted parity-prefix family concurrence

Hypothesis56; structural math/workflow only, no fitting, symbolic screening
or certificate execution. Gate conventions q0 MSB, chronological, four wires,
ancilla-free arbitrary locals plus directed all-to-all CX. The exact13 incumbent
persists; this is an untested family until the numerical researcher's reserved
round is executed and independently checked.

The exact chronological primitive rule is

L0; CX0→1,CX2→1,CX3→1; L1; F02; L2; F13; L3; CX1→2;
L4; F01; L5; F23; L6,

where each L is a product of four independently parameterized SU2 gates and
Fct(a,b)=exp(i a XcXt+i b YcYt). There are no free local gates inside the
three-CX prefix. Seven local layers give84real coordinates; four F blocks
give eight angles, total92. Each F has the exact two-CX compilation from
topology14_exact_matchgate_derivation.md, so the compiled count is12 before
any demonstrated cancellation. Native gate model {CX,F(a,b),locals} has
four CX and four F, eight entanglers; this is not an eight-CX construction.

On computational-basis inputs to the prefix, q1 becomes q0 xor q1 xor q2 xor
q3 and the other bits stay unchanged. This is an exact property of the
three-gate group. An arbitrary preceding L0 means this tag need not represent
the original input Hamming parity. Subsequent F13 and arbitrary local layers
need not preserve Z1, and they are not constrained to do so. Thus a proof
assuming parity remains a final output label would not apply to this family.
Although Fct commutes with ZcZt and with full Z parity, arbitrary local layers
and the bridge preclude importing that restricted invariant as a no-go result.
No diagonalizing member or invariant obstruction has been proved here.

Frozen numerical budget is two seeds, at most500iterations each,240seconds
total, one thread. Literal parameterization, gate list, seeds, code/source
hashes and complete checks belong to the numerical run and verifier preflight.
Historical duplication review should compare the full primitive sequence,
restricted layer semantics and initialization rule. Reversing a CX alone is
not a new full-local architecture because local H conjugation changes its
direction; this restricted family additionally has fixed block-internal
relations, which must not be silently replaced by fully free local layers.
Existing exact14's six-CX correction prefix differs from the three-CX grouped
prefix, but that observation alone is not an exhaustive historical audit.

The prior proposed commutant/support-only lower-bound system is corrected in
joint_simplify/REPORT.md. With full supports it contains known exact13
observables independently of six-CX reachability and is vacuous. Future
lower-bound work needs actual shared-template reachability or independently
derived axis constraints before reserving any symbolic elimination.

Adversarial review: literal chronology, prefix scope,92parameter count and
separate native/compiled counts agree. No confirmed defect. Remaining issue
is existence of any passing member, not cost accounting. Structural artifacts
and actual exchanges retained; no extra search or certification duplicated.
