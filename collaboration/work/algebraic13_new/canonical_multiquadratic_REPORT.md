# Exact13 CNOT radical certificate: run370

Command: `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 collaboration/work/algebraic13_new/canonical_multiquadratic_verify.py`.

The one authorized300-second exact attempt completed in14.32 seconds. Every one of256 entries of `U_projective V4−D U_projective` reduced to an empty coefficient dictionary by exact Python `Fraction` arithmetic. Output: `exact_diagonalizer:true`,13 CNOTs,60 locally normalized rotations. `canonical_multiquadratic_residuals.json` is the empty list. Source SHA256: `d500ff2a79c4717ff38d935d71b67bd798b06a7c6899a275fec25570b86b8d13`, corresponding to `canonical_algebraic_ansatz.json`.

## Exact algebra

The coefficient field is Gaussian `Q(sqrt(2))`: four rational components per coefficient. Each positive radical is represented by a positive `Q(sqrt(2))` multiplier and a product of selected positive generator roots. The five generator radicands are:

1. `1/2−sqrt(2)/4`;
2. `2/3`;
3. `3/5`;
4. `(252−45sqrt(2))/367`;
5. `(3−sqrt(2))/7`.

Each representation was found by an exact subset square test and checked by the identity `factor² × product(radicands)=desired_square`, with factor positivity determined by rational comparisons. Signed sine/cosine branches were resolved exactly from the input expressions. Generator relations and input hash are saved in `canonical_multiquadratic_basis.json`.

Elements use at most32 bitmask monomials. Multiplication uses mask XOR and the product of radicands for shared generator bits. No floating-point arithmetic, numerical zero tolerance, primitive-element construction or per-entry symbolic simplification is used in global assembly or zero testing. SymPy is used only to parse initial exact expressions, extract their squared `Q(sqrt(2))` values, and establish exact signs.

The output field called `basis_rank` means the number of selected generators, **not a proved algebraic independence rank or field degree**. Independence is unnecessary: the formal algebra maps to the actual positive real roots, so an empty dictionary implies exact zero. A nonempty dictionary would not establish nonzero if generators were dependent; this verifier makes no such negative inference.

## Normalized gates

Every local gate exactly satisfies `sin²(theta)+cos²(theta)=1` in the dictionary algebra. The projective representation is `(1+cos(theta))I−i sin(theta)P`, except `cos=-1`, where pureP is used. The independent coordinator local review (`collaboration/work/root/radical13_local_exact_review.json`) and `canonical_local_unitarity_proof.md` establish exact realness and nonzero projective factors. Hence dividing by the nonzero total scalar transfers the exact global identity to the normalized unitary gate prescription. Gates are chronological, q0 is most significant, target is the right shift, with output labels inherited from run359.

## Review and state

Adversarial review checked Gaussian-field multiplication, generator-intersection reduction, positive-root branch selection, left multiplication by local gates, CNOT row permutation and right-shift column permutation. No confirmed implementation defect found. Remaining requirement before incumbent promotion: coordinator independent reproduction and review of this certificate. The earlier300-second SymPy expansion timeout remains preserved as run368; this distinct sparse arithmetic method is run370. No original artifacts or shared checker/native code/status documents were changed; no collaborator commit. The exact claim is submission4 on hypothesis14, closed completed. All evidence artifacts retained intentionally.
