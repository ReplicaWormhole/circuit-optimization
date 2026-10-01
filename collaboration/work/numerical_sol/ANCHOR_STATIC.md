# Static anchor steering and remaining CX saving

Board124 initially described a16-word CSsign/control/Tsign anchor. The
coordinator stopped that literal proposal before source/config generation,
reservation or matrices after algebraic_sol found a d11 obstruction. That
proposal used active Y in first P and consequently retained Y_x after CH.
No numerical exclusion is claimed; no experiment was needed for this step.

The first P can instead be target-conjugated by S: chronological
S2,RCCX(v=q3,x=q0;t=q2),Sdg2. ActiveY mapsX and inactiveZ remainsZ. In the
v=1 sector its output-side error relative to ideal CX(x->y) is
E(x,y)=(-1)^((1-x)y), y=q2. In d11, exact P maps the exchange to a scalar
depending on y times X_x. Conjugation by E replaces X_x by(-1)^y X_x;
CH H_x therefore gives a diagonal operator with changed root ordering.
This is a by-hand sector argument, not a complete gate-list certificate.

Algebraic_sol proposed final target V'=Ry(-3pi/4)Rx(-pi/2), with activeY
mapping to H and inactiveZ mapping to Y. Literal wrapper chronology is
X1,Ry2(3pi/4),Rx2(pi/2),RCCX(1,3;2),Rx2(-pi/2),Ry2(-3pi/4),X1.
The inactive Y lies only in d00 for controls(1-u,v); Y is a computational
basis permutation up to phases, preserving diagonality of that sector CZ.
This needs independent complete matrix/exact checking before a valid anchor
claim. With exact CSdg B, nominal count is13CX, not the requested12.

No shared checker edits, new dependency, random draw, gate-matrix evaluation
or optimizer was performed. The exact12 objective remains open.
