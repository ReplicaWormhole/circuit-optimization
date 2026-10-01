# Two composed spectral-label targets

Hypothesis29: one structural report; no numerical screening, optimization or
exact-certification work. Source is the invalid twelve-CX saved candidate
`collaboration/work/numerical12/deletion_slot1_deep.json`; fixed base labels
are those independently reviewed for run375 in numerical_label12/targets.json.
Labels0,1,2,3 mean roots1,i,-1,-i. Wires are ordered q0q1q2q3 in binary rows.
The gate model is arbitrary one-qubit gates plus directed all-to-all CX,
four wires and no ancillas. Incumbent13CX and global interval6<=Cmin<=13 persist.

TargetA composes disjoint swaps(2,6)(3,7). TargetB composes disjoint
swaps(8,13)(9,12). Every swap only permutes output labels, so both preserve
the necessary exact cycle multiplicities(6,3,4,3). These label changes specify
objectives for fresh fits; they do not append free gates to the old circuit.

TargetA's permutation flips q1 iff q0=0 and q2=1, independently of q3.
TargetB's permutation flips q1 and q3 simultaneously iff q0=1 and q2=0:
1000↔1101 and1001↔1100, with every other row fixed. Neither permutation
is a tensor product of individual-wire permutations.

For a diagonal root target D, reduced matrices obtained by tracing other
wires are diagonal root sums. Under any product-local output conjugation W,
the reduced matrix on qj is conjugated by Wj, preserving its eigenvalues.
Thus unequal reduced spectra prove that two targets cannot be related by
free product-local output gates, even allowing those gates to be nondiagonal.

| Target | Reduced q1 spectrum | Reduced q3 spectrum |
|---|---|---|
| Base | {1+i,1-i} | {0,2} |
| A | {3-i,-1+i} | {0,2} |
| B | {1+i,1-i} | {-2-2i,4+2i} |

For q1 in A, rows2,3 change from(-1,i) to(1,-i), adding2-2i to
the q1=0 sum, with the opposite change in q1=1. For q3 in B, rows8,12
change respectively1→-i and i→-1, adding-2-2i to the even-row sum0;
the odd-row sum2 gains2+2i. The q1 sum for B and q3 sum for A are unchanged.
The two B singles individually have q3 spectra{-1-i,3+i}. Therefore A,B
are mutually distinct by this invariant; A differs from base and its two
singles, and B differs from base and its own two singles. This does not rule
out every other gate-cost-preserving symmetry or establish different basins.

The numerical budget is two targets, each200L-BFGS iterations and at most
five SVD corrections,120seconds total and one compute thread. The numerical
researcher owns freezing exact label vectors and source/code hashes, reserving
the experiment before execution, fitting all168 local coordinates, and saving
complete chronological lists. The verifier independently checks these lists.

Repeated run375 residuals near sin(pi/12)=(sqrt6-sqrt2)/4 suggest a special
angle might describe a stationary basin. This exact elementary identity follows
from sin(pi/4-pi/6); the identification with a circuit basin remains conjectural.
A residual may arise as a sine entry of a two-dimensional rotation block, but
no such invariant block or forced angle has been proved for this candidate.
Finite fitting failure and apparent angle recognition imply no topology or
twelve-CX impossibility. No additional fit was made for this structural report.

Adversarial review: exact changes in Gaussian-integer partial-trace sums agree
with total trace2; permutations and label capacities are preserved. No
confirmed defect. Remaining risks are unknown fitting behavior and incomplete
classification of all circuit symmetries. Report/chat files are intentional
kept evidence; coordinator owns shared commit.

Verifier independently reproduced every listed reduced spectrum by direct
Gaussian-integer sums and confirmed B's conditional bit action. Numerical
researcher's direct run376 report gives A maxoff0.1517277915509333 and selected
target loss0.009459905009667346, B maxoff0.7071067918307067 and target
loss0.46343390751450736; both remain invalid. This reported improvement in A
already disproves treating the prior0.258819 plateau as a universal numerical
floor. Independent list checks belong to the verifier's separate report.

Coordinator's unscheduled next proposal is a label-free Hessian audit of the
improved A source, followed only if justified by two signed curvature starts.
First compute the label-free gradient: stopping a selected-D Gauss-Newton fit
does not establish stationarity of the offdiagonal objective. Retain gradient
norm, Hessian eigenpair residual and directional curvature. Coordinate/gauge
redundancies can produce near-zero curvature; negative curvature below a
stated numerical threshold motivates bounded perturbation fits but gives no
exact saddle certificate. Positive curvature or failed signed fits supplies
neither global optimality nor a topology exclusion. No audit was executed here.
