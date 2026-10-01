# Independent spectral-label review

Hypothesis27 covers math, computation and workflow. No optimization and no new
exact-certification run is performed here. The independent checker imports no
shared circuit simulation. It assembles U3 matrices by Kronecker products,
directed CX by computational-basis bit columns, and the right cycle by
`x → (x >> 1) | ((x & 1) << 3)`. Qubit0 is most significant and chronological
gates multiply on the left. The target is UV=DU.

The source hash is
`9a9dd0e0b430e935f99f4d998fb25d5e1b48160ad16a6c46a27e2574c59c2398`.
The source independently has12CX, maxoff0.2588190369994409 and normalized
offdiagonal squared loss0.050240473580837317. It is invalid. Nearest labels
already have counts6,3,4,3; source-check JSON records every label and residual.

Character traces(V^k)=(16,2,4,2) give root multiplicities(6,3,4,3) in root order
(1,i,-1,-i), by the four-term discrete Fourier transform. For unitary U,
right multiplication by U† gives ||UV-DU||F=||UVU†-D||F. Consequently capacity
assignment on diagonal entries selects closest admissible labels at this U.
A finite collection of swapped targets and local fits does not quantify over
all output orders, all degenerate eigenspace bases, or all twelve-CX circuits.

The reviewed frozen plan has four unequal-label swaps ranked by descending
offdiagonal coupling and index tie breaks;100L-BFGS iterations and at most
three truncated-SVD steps per target, with120seconds overall and one thread.
The timing starts before matrix target screening. Hungarian repeated-root
assignment is used; no lexicographic minimum label-vector tie guarantee is
claimed by its config. No randomness is used. Source and helper hashes are
recorded in config. Numerical passing, if any, still requires separately
planned independent exact certification before promotion.

Costs stay separate: fixed directed all-to-all CX with arbitrary locals has
an accepted exact13 incumbent and these submitted lists have12CX. A failing
list supplies no12CX upper bound. The different native set containing CX and
F(a,b)=exp(i a XX+i b YY), with independently tunable real a,b and arbitrary
locals, has a10-entangler baseline (6CX+4F), depth8, compiled-CX upper14.
Its direct JSON checker is numerical with analytic regrouping provenance.
This does not give a10CX circuit or transfer the six-CX lower bound. Unknown
compiled cost would remain unknown, never zero.

Saved-output verification and final adversarial review follow below.

Independent checks of all four saved lists found12CX each, unitary errors at
most1.56e-15, and maxoff respectively0.2588190451027492,
0.25881904511477694,0.2794011825428628,0.2794011768358336. All are invalid.
The final selected-target losses are0.05111126056639771,
0.05111126056641396,0.2883009523150334,0.28830095231503283. For the latter two,
nearest-root labels differ from selected targets; the two residuals must not
be confused. Source layout, including adjacent local layers around identity
slot1, passes a sequential gate parser; optimizer/source matrices agree up to
phase within1.45e-15. All four target pair scores and labels agree independently.

Read-only exact arithmetic review of an unexecuted possible joint swap
(2,6)(3,7): this flips q1 precisely when q0=0 and q2=1, regardless of q3.
Reduced q1 diagonal spectra for base/swap26/swap37/joint are respectively
{1+i,1-i}, {3+i,-1-i}, {1-i,1+i}, {3-i,-1+i}. Gaussian-integer sums reproduce
these values independently. Distinct reduced spectra prove that joint target
is not equivalent to those three by product-local output conjugation. This
is a potentially distinct initialization target, not a tested improvement.

Commands: `python3 collaboration/work/verifier_label12/verify.py` with output
JSON path followed by saved input paths; `OPENBLAS_NUM_THREADS=1
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3
collaboration/work/verifier_label12/audit.py`. Both pass. Script timing guard
is checked between evaluations/steps; a single Jacobian/SVD operation is not
preemptively interrupted. Observed14.848545seconds is safely below120seconds.

Adversarial review: no confirmed defect in gate conventions, chronology,
normalization, source binding, target scope or reported numerical costs.
One reporting hazard (optimizer loss before Newton versus final target loss)
was resolved by independent final matrix residuals. No exact success or
optimality claim is supported. Only assigned per-agent files are written;
all JSON/script/report/chat artifacts are intentional evidence for the
coordinator's local commit. No push and no incumbent changes.
