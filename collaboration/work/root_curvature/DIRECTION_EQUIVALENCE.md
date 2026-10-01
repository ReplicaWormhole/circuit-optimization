# Direction-only changes are not new unrestricted-local families

Board39, structural derivation only. For an ordered pair c,t, write
C_c,t=(I+Z_c)/2 tensor I_t +(I-Z_c)/2 tensor X_t. Conjugating by H_c H_t
replaces Z_c by X_c and X_t by Z_t, giving
(I_c tensor I_t + X_c tensor I_t + I_c tensor Z_t - X_c tensor Z_t)/2.
Reordering terms equals I_c tensor(I+Z_t)/2 + X_c tensor(I-Z_t)/2=C_t,c.
Thus (H_c H_t) C_c,t (H_c H_t)=C_t,c, exactly.

In a chronological CX schedule with arbitrary independent one-qubit layers
between each entangler and at both ends, these local Hadamards can be absorbed
into the neighboring local layers. This proves equality of reachable unitary
sets for every orientation of the SAME unordered chronological pair sequence,
at the same CNOT count. SU2-only layers may realize H up to phase, which is
irrelevant to UV=DU. The statement does not apply to restricted local ansatzes,
fixed angles or absent local layers. It does not say two different unordered
pair sequences are equivalent, nor that identical initial coordinates behave
the same under optimization.

Historical ordered topology hashes and failed runs are preserved. No count of
distinct canonical families is inferred without an actual audited grouping.
Future architecture proposals should change unordered pair order/incidence,
not rely on direction reversal alone. New fits still require a source,
bounded budget and reservation even if intended as alternate initializations.

Independent verifier confirms the projector identity, absorption of consecutive
local factors and the SU2 phase caveat. It distinguishes reachable-family
equality from optimizer behavior. Direct peer findings and board handoff saved.
