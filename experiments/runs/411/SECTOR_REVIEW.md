# Independent complete phase check of analytic11CX proposal

Board137; by-hand only, no matrices, symbolic search, optimizer or checkers.
Let x=q0,y=q2,u=q1,v=q3. Begin with certified12 merged M using the firstP
order(x,v), and write N4=(-X+Z)/sqrt2. Replace its last controlled reflection
by N4prime=Ry(-pi/4)N4=-cos(pi/8)X+sin(pi/8)Z. The last T_x commutes with
the controlled rotation, so Mnew=C_v Ry_y(-pi/4) M exactly, including phase.
N4prime=Ry(-7pi/8) X Ry(7pi/8): chronological Ry_y(7pi/8),CX(v,y),Ry_y(-7pi/8).
The other three reflection blocks and T_x remain unchanged; merged cost4CX.

After Mnew and unchanged CH(u;x), the four sector targets are:

- d00: CZ_xy.
- d10: (-i)^y Z_x, already diagonal.
- d01: i^x H_y because Ry(-pi/4) X Ry(pi/4)=H.
- d11: Z_x[(1+i)I_y/2+(i-1)M_y/2], M_y=(Z-X)/sqrt2.

The last line retains the scalar/relative phase of f_i(y)=diag(i,1), rather
than checking only its eigenvector axis. Ry(-pi/4) Z Ry(pi/4)=M.

Replace the old3CX final selector by chronological CZ(u,y) then controlled
Nv(v,y), Nv=sin(pi/8)X+cos(pi/8)Z. The first gate applies Z_y iff u=1;
Z M Z=H. The second reflection maps H to Z exactly by Nv H Nv=Z.
Thus the final sector targets, including every scalar, are CZ_xy,
(-i)^y Z_x,i^x Z_y,f_i(y) Z_x. Both d00,d10 stay diagonal. d11's scalar
identity coefficient is unchanged throughout and its M coefficient rotates
to H then Z with positive sign, preserving the complete f_i factor.

Nv=Ry(-3pi/8) X Ry(3pi/8), so its literal oneCX chronology is
Ry_y(3pi/8),CX(v,y),Ry_y(-3pi/8). CZ is H_y,CX(u,y),H_y. The final block
costs2CX. Total2Bell+2router+4merge+1CH+2final=11CX, no ancilla.

In MSB computational order q0q1q2q3 the predicted exact D is
[1,1,1,-1,1,i,-i,1,1,i,-1,-i,-1,-i,i,-1],
labels[0,0,0,2,0,1,3,0,0,1,2,3,2,3,1,2] in roots(1,i,-1,-i).
Compared with the12 word, d00 now remains CZ instead of being permuted by
Y_y; this is allowed arbitrary output eigenvalue order. Multiplicities are
6/3/4/3 in the stated root order. Angle denominators8 require conductor32.

Adversarial review checked multiplication order, N4prime and Nv axes/signs,
Tx commutation, both diagonal and mixing sectors, all scalar phases, final
positive controls, honest11CX count and label prediction. No unresolved
issue found in this by-hand derivation. Complete saved-list numerical/exact
checks plus independent certification remain required before acceptance.
