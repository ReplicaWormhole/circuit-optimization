# Algebraic round 3: decorated rewrite and joint halfshift obstruction

Artifact types: math and workflow. Hypothesis20, owner algebraic_round3.
One bounded structural investigation; no optimizer, symbolic search, new ledger run, or duplicate exact13 certification. The accepted exact13 incumbent is preserved. All-to-all arbitrary local gates plus directed CNOT, ancilla-free; qubit0 most significant and gate lists chronological.

## Decorated three-wire identity

Let a,b,c be distinct and C_ab be CNOT(a,b). For local layers L,M each commuting with C_ab, the chronological lists

    C_ab, L, C_bc, M, C_ab
    L, C_bc, C_ac, M

implement the same unitary. In matrix multiplication order the first is
C_ab M C_bc L C_ab = M C_ab C_bc C_ab L = M C_ac C_bc L.
The binary action of the bare chronological triple is
(a,b,c)->(a,b xor a,c)->(a,b xor a,c xor b xor a)->(a,b,c xor b xor a),
which equals C_bc followed by C_ac. A sufficient allowed local family is
Rz(a), Rx(b), and arbitrary one-qubit gates on spectator c. The commutation requirement is sufficient, not necessary or a full classifier.

The canonical source has chronological CX pairs
01,32,12,10,23,21,01,32,12,10,23,21,01.
It has no contiguous aba motif to which this rule immediately applies.
Two new12CX architecture hypotheses replace slots0..2 or6..8 of this list
([01,32,12]) by [12,02]. They are obtained by TWO architecture mutations: [01,32,12] becomes [01,12,01] (32->12 and last12->01); the bare aba reduction then gives [12,02]. They are NOT exact reductions of the accepted source: local decorations and the original32 edge change. Numerical researcher froze fixed13-slot templates with [None,12,02], full local freedom and canonical warm starts plus noise. Their performance must be independently measured.

## Joint halfshift condition

Use W=U dagger so W dagger V4 W=D is diagonal. Let p_j and p_k be unit single-qubit Pauli axes on distinct wires, each having nonzero off-diagonal entries in the computational basis. Define A=W p_j W dagger and B=W p_k W dagger.

Claim: A and B cannot BOTH commute with V4 squared.

Proof: Pulling both commutators back gives [p_j,D squared]=[p_k,D squared]=0. The off-diagonal entry of p_j forces the corresponding diagonal D squared labels to agree on every pair of strings differing in bit j, and similarly for k. Thus D squared is constant on each four-string square indexed by the other two bits. Its +1 and -1 multiplicities are both divisible by4. But V4 squared is SWAP02 SWAP13. Each two-qubit SWAP has multiplicities3 and1; their tensor product has +1 multiplicity3*3+1*1=10 and -1 multiplicity3*1+1*3=6. Unitary conjugation preserves these multiplicities, a contradiction.

This condition imposes no support restriction. The earlier proper-support formulation is redundant: minimal support of an operator commuting with V4 squared is invariant under swapping opposite pairs, so a nonscalar proper-support operator must be confined to an opposite pair, already excluded for input-axis images. The support-free joint condition may constrain full-support selected axes in a future joint classifier.

More generally, m distinct input axes with nonzero off-diagonal entries commuting with D squared force both multiplicities divisible by2^m, hence m<=1. No existing sixCX schedule is excluded by this note: a schedule-specific argument must force two such commutations and establish their transversality. Surviving3024 schedules remain necessary possibilities, not circuit witnesses. This is a structural necessary condition, not a sevenCX proof or optimum certificate.

## Review and handoff

Direct algebraic/numerical/verifier exchanges are preserved verbatim in CHATS.md and actionable messages on board20. Independent verifier confirmed both the decorated rewrite and the support-free joint condition; see CHATS.md. An adversarial review corrected the initial mistaken one-mutation architecture motivation; the frozen12CX schedules and their results are unchanged. Self adversarial review checked chronological identity, spectator scope, support-vacuity issue, both multiplicities, and distinction between a topology hypothesis and a source-equivalent rewrite. No additional computation was needed for the definition-level derivations. Coordinator owns final review, shared records and commit.
