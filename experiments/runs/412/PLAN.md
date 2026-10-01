# Independent exact analytical11CX validation in Phi32

Board138, numerical_sol. Source is the authoritative verifier11router/
candidate11/candidate.json, with30 chronological gates and11CX. No circuit
matrix or coefficient computation is allowed before coordinator reservation.

This deliberately adapts standalone run407's independent Fraction verifier,
using16coefficients modulo z^16+1 (Phi32) at z=exp(i*pi/16). The quotient is
exact at this embedding. Multiplication reduces every degree>=16 by z^16=-1;
z powers have period32; conjugation mapsz->z^-1; i=z^8 and1/sqrt2=(z^4+z^-4)/2.
Phase exp(i*f*pi) uses exponent16f, rotationhalfangle uses exponent8f.
Nonintegral exponents are rejected. All input phases are retained exactly.

Independent primitive matrices, chronological row updates and CX permutations
construct U, with qubit0 most significant. Independent bit-label permutation
constructs V4|x0x1x2x3>=|x3x0x1x2>. Discover one root in(1,i,-1,-i) per row,
verify all256 entries of UV4=DU and all256 of UUdag=I, count11CX and root
multiplicities(6,3,4,3), retain every exact16-coefficient matrix entry. No
NumPy,SymPy, shared checker or existing unitary/gate implementation is imported.
Expectedlabels from other validators are not supplied to this certificate.

One deterministic90-second invocation, one thread, seed0/noRNG, nooptimizer.
Freeze directconfig and source/candidate/Python/Fraction/math-extension hashes;
root reserves exactmethod "exact Fraction polynomial arithmetic in Q[z]/(z^16+1)"
and matching codehash, createsworkspace, copiescandidate/config/PLAN and
linksboard138. Guard checksactiveledger/config/method/boardowner/runlink and
exclusiveexecutionmarker beforematrices. Rootownsdispatch, terminalledger,
acceptance, closure andcommit. Command:
`timeout 90s python3 experiments/runs/ID/certify.py`.

This certifies only the complete saved literalancilla-free11CX gate list.
No optimum, newlowerbound, totalspindiagonalization, strongSchur or topology
exclusion follows. Unsupported angles failratherthanbeingrounded.
