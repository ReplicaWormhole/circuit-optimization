# Exact joint trace obstruction in one fixed first-axis branch

This note concerns the six-CNOT ordered pair schedule

\[
(01),(23),(12),(01),(23),(03),
\]

and the **fixed** first selected-axis image

\[
A=(Z_0X_1+X_0Z_3+Y_0X_1Z_3)/\sqrt3.
\]

The first axis is fixed to this operator, not classified over its full
24-dimensional image span. The second selected axis on input wire 3 has
the sound 48-dimensional output-span overapproximation described in
`global_lower_bound7_joint_aa_pair.md`. Requiring it to commute with
this A gives a 12-dimensional real Pauli span. All results below apply
to **every** B in that span that is a Hermitian involution and obeys the
single-axis spectral-path identities. The span is an overapproximation
of realizable circuit images, so excluding this entire B family is
sound for the fixed-A branch.

## Exact algebraic family

The twelve centralizer labels are
`IIIZ, IIXZ, IIYZ, IIZZ, ZIIX, ZIIY, ZIXX, ZIXY, ZIYX, ZIYY, ZIZX, ZIZY`.
Form the 78 degree-two monomials in their real coefficients. The
spectral-path Gram has integer entries and exact rational rank 69.
For any real zero, the four coefficients of
`IIIZ, IIZZ, ZIZX, ZIZY` vanish. The remaining eight coefficients
`x0,...,x7` correspond, in order, to
`IIXZ, IIYZ, ZIIX, ZIIY, ZIXX, ZIXY, ZIYX, ZIYY`.
The exact rational Groebner basis for the spectral-path equations
plus \(B^2=I\) is saved in
`global_lower_bound7_joint_aa_centralizer_certificate_result.json`.

Here is a parameterization of every real solution, derived from that
basis. Put \(u=x_2=x_0\), \(v=x_3=x_1\), \(a=x_5\),
\(b=x_6\), \(c=x_7=-x_4\), and \(s=b-a\). Its equations give

\[
ab+c^2=0,\qquad a^2+b^2+2c^2=1/3,
\]

so \(s^2=1/3\). Also \(u^2+v^2=1/3\), \(v^2=bs\), and
\(uv+sc=0\). Thus, for \(\epsilon\in\{\pm1\}\) and real
\(|t|\leq1\), every solution obeys

\[
s=\epsilon/\sqrt3,\quad u=t/\sqrt3,\quad
a=-\epsilon t^2/\sqrt3,\quad b=\epsilon(1-t^2)/\sqrt3,
\quad v=\eta\sqrt{1-t^2}/\sqrt3,\quad
c=-\epsilon\eta t\sqrt{1-t^2}/\sqrt3,
\qquad \eta\in\{\pm1\}
\]

At \(|t|=1\), both choices of \(\eta\) give the same point. The
conclusion below uses only \(s,u,a\), so \(\eta\) is immaterial.
Direct Pauli trace rows then give

\[
\begin{aligned}
\operatorname{tr}(V_4 A)&=2i/\sqrt3,&
\operatorname{tr}(V_4^2 A)&=0,\\
\operatorname{tr}(V_4 B)&=2i\epsilon/\sqrt3,&
\operatorname{tr}(V_4^2 B)&=0,\\
\operatorname{tr}(V_4 AB)&=(8t-2\epsilon)/3,&
\operatorname{tr}(V_4^2 AB)&=(8t-4\epsilon t^2)/3.
\end{aligned}
\]

## Finite diagonal-eigenvalue assignment

If A and B came from local input axes \(p_0,p_3\) of one
diagonalizer U, write \(D=U^\dagger V_4U\) and let \(z_0,z_3\)
be their respective input Bloch Z components. For each positive
integer power r, cyclicity of trace and the computational-diagonal
form of D force

\[
\operatorname{tr}(V_4^r A)=z_0\operatorname{tr}(D^rZ_0),\quad
\operatorname{tr}(V_4^r B)=z_3\operatorname{tr}(D^rZ_3),\quad
\operatorname{tr}(V_4^r AB)=z_0z_3\operatorname{tr}(D^rZ_0Z_3).
\]

The eigenvalues of D are \(1,i,-1,-i\) with multiplicities
\((6,3,4,3)\). Group its 16 diagonal entries by the input bits
\((x_0,x_3)\). Each of the four blocks contains four entries.
There are exactly **7,168** possible four-block count arrays with
these multiplicities. Exact enumeration gives 76 satisfying both
single-axis moment requirements for real \(|z_0|,|z_3|<1\), 52 of
these with a real first AB moment, and **zero** satisfying both AB
moments for any \(\epsilon,t\). The script uses integer Gaussian
moments and rational arithmetic for t; it applies no floating-point
acceptance threshold. The strict Bloch bound is harmless here:
the first-axis moment forces \(z_0=2/(\sqrt3 m_0)\) with nonzero
integer \(m_0\), so \(|z_0|=1\) is impossible, and similarly for B.

For the earlier explicit B witness, \(\epsilon=-1,t=1\), and its
input Bloch components are \(z_0=z_3=1/\sqrt3\). The first moments
would force the \((0,0)\) block sum of eigenvalues to be 3. The
second moments would force all four squared eigenvalues in that
block to be +1, hence all four eigenvalues are ±1 and their sum is
even. This is a short independent contradiction for that single pair.

Reproduce the exact computations with:

```bash
python3 global_lower_bound7_joint_aa_centralizer_gram.py
python3 global_lower_bound7_joint_aa_centralizer_certificate.py
python3 global_lower_bound7_joint_aa_trace_moments.py
python3 global_lower_bound7_joint_aa_family_trace.py
```

Ledger runs 275, 277, 272, and 281 record these computations. The
two temporal cut screens (runs 268 and 271) excluded zero of the 401
dihedral classes and are recorded in their separate scripts/results.

**Scope:** This is an exact necessary-condition contradiction for the
fixed canonical A branch and every B in its sound 12-dimensional
centralizer. Other A values in the first selected-axis image span
remain unclassified. No schedule or dihedral orbit is excluded, and
the established six-CNOT survivor count remains **3,208** (401
dihedral classes). This does not establish a global seven-CNOT lower
bound.
