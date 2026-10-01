# Verifier round 3

Artifact types: math, computation, workflow, rendering code.

## Explicit circuit

`EXPLICIT_CIRCUIT.md` lists every chronological gate in the accepted 13-CNOT source, including zero rotations. `export_circuit.py` preserves signed half-angle radicals, full-angle metadata, rational multiples of pi where supplied, coordinate identifiers, and output labels. The source hash is checked against `EXACT13_ACCEPTANCE.json` before export. Gate set: arbitrary exact single-qubit gates and directed CX on any pair, four qubits, no ancillas. Current interval remains 6 <= C_min <= 13; optimality is unproved.

Verification: `python3 collaboration/work/verifier_round3/export_circuit.py` passed; separate parsing of Markdown rows compared every cell with canonical JSON and passed (73 rows, 13 CX, 60 rotations, all exact metadata and output labels). A first ad hoc parser had a leading-delimiter bug and failed before comparison; corrected strip-and-split parser passed. No full-matrix certificate rerun and no new exact acceptance.

## Independent structural review

For distinct wires a,b,c, let C=CX(a,b), E=CX(b,c), F=CX(a,c). Bare identity C E C = E F = F E follows on bits: (a,b,c) maps to (a,b,c+b+a). If [L,C]=[M,C]=0 then chronological C,L,E,M,C has matrix C M E L C = M C E C L = M F E L, hence chronological L,E,F,M. Sufficient local decorations include Rz(a), Rx(b), and arbitrary rotations on spectator c. Without both commutations this rewrite is not justified. The accepted schedule has no adjacent aba entangler motif, so the numerical architecture mutation is a new fitting hypothesis, not an equivalent replacement.

Joint obstruction independently reviewed: if U V4 U†=D and two distinct input-wire Pauli axes p_j,p_k each have a nonzero off-diagonal component, then [D²,p_j]=[D²,p_k]=0 is impossible. Each commutation forces D² labels equal across its bit flip; combined, all four vertices of each two-bit square have equal labels. Thus the +1 multiplicity is divisible by four. But V4² fixes precisely the four binary strings with x0=x2 and x1=x3, so trace=4 and multiplicities are (16+4)/2=10 and6. Contradiction. Equivalently images U†pU commuting with V4² yield the same obstruction; conjugation orientation must be stated consistently. This is not a circuit-count lower bound until a given schedule forces both commutations.

## Adversarial review

No confirmed artifact defects after checking chronological order, rotation signs, source hash, label order, and claim scope. The full-angle columns do not replace signed half-angle prescription. Matrix certification is inherited from accepted record, not reproduced here. Native gate counts remain separate from CX counts. Remaining risk: optimum is unknown, new numerical results need independent checks and exact certification before promotion. All files here are kept research evidence; coordinator owns commits.

## Numerical independent check

`python3 check_circuit.py collaboration/work/numerical_round3/case0.json` and corresponding `case1.json` both return expected exit1 for invalid diagonalizers. Counts are 12CX each; maximum off-diagonal residuals are 0.4487654196973814 and 0.6728896343792273. Unitarity errors are 8.881784197001252e-16 and 5.551115123125783e-16. Thus neither submission improves the incumbent and no exact certificate run is warranted. These two failed fits do not exclude their topologies or prove a lower bound.
