# Independent reordered-suffix review

Board150, verifier10reorder. Conventions: chronological gates, q0 MSB, x=q0,u=q1,y=q2,v=q3. This is a static by-hand review; no angle, matrix, RNG, AD or optimizer evaluation.

Accepted exact11 is prefix4 plus merge4 plus CH1 plus selector2. The proposed ten-CX chronology is02,13,01,23|02,32,10,02,12,32, with prefix4 and suffix6. Saved416ordinal7 chronology is02,13,01,23|02,32,02,10,12,32. Swapping the noncommuting02/10 changes the bare action. Full arbitrary-local-gate families may nevertheless be equivalent; no equivalence or inequivalence proof is available.

Bare triangle is independently valid including phases: chronological CX02,CX10,CX12 sends(x,u,y) to(x xor u,u,y xor x xor u). Chronological CX10,CX02 gives exactly the same computational-basis map. CNOT entries are unit positive, so equality carries no hidden scalar phase.

This does not prove the corresponding wrapped controlled-H/reflection identity. Local x rotations inside CH do not commute with x-controlled target operations; replacing the bare triangle while ignoring those wrappers may alter eigenvalue phases and bases. Need complete chronological primitive word and all16 sector eigenidentities, or exact full-list certification, before claiming success.

The previous corner-product obstruction fixes computational u/v sectors and earlier x/y prefix/CH. Moving CH into the merge alters the x basis before the final x-controlled reflection, outside that theorem. Unrestricted local refitting also frees prefix locals, so numerical architecture should not be described as holding the exact prefix fixed.

Transplanting saved416best7 raw132vector verbatim is a deterministic initialization on a new topology; it asserts neither preserved unitary nor decreased loss. Its source endpoint is invalid maxoff.20445939300014332, not near-exactness. Root-authorized separately frozen preflight must independently validate serializer/NumPy/Torch/gradient consistency before fitting. No10CX solution or exclusion follows from the structural proposal.

New restricted early-CH obstruction independently reviewed by hand: uv10 with M=diag_x(I,S), D=Sdg_y Zx gives W=M†(Sdg_y Xx)M=[[0,I],[Z,0]], since Sdg*S=I and Sdg²=Z. For Q=diag(Q0,Q1), E=Q1Q0† and P=Q0ZQ0†, QWQ† has upperE† and lowerEP. If E is phase times a reflection, E² is scalar. x-only K diagonalization requires proportional offdiagonal target blocks; EP=rhoE† then forces P scalar, impossible for conjugatedZ. Later computational-x-blockdiagonal gates preserve nonzero xoffdiagonal blocks. This assumes no later nondiagonal x locals and frozen computational sector/control basis; the full132 free family violates those restrictions and is not excluded.
