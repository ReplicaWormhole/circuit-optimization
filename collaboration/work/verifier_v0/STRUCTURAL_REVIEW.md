# Independent routed target and old v0 source review

Board157/verifier_v0. By-hand only; no matrix,RNG,AD,symbolicsearch or angle recognition. Q0 is MSB; gates chronological. Let Bell(a,b)=(|0,b>+(-1)^a|1,1 xor b>)/sqrt2. Swapping its physical wires multiplies by(-1)^(ab).

Accepted411 prefix P has chronologyCX02,H0,CX13,H1,CX01,CX23. Rightshift exchanges opposite Bellpairs and reverses the old second pair: E†V4E sends |a,c,b,d> to(-1)^(cd)|c,a,d,b>. Routerlabels x=a,u=a xor c,y=b,v=b xor d. Thus P V4 P† maps(x,u,y,v) to(x xor u,u,y xor v,v) with phase(-1)^((x xor u)*(y xor v)). This is exactly W=CZ_xy Xx^u Xy^v with flips acting first.

V0 restriction uses n3wireorderx,u,y and operator W0=CZ_xy CX(u->x). It is not a threewire cyclic shift. u0branchCZ_xy;u1branch inxblocks [[0,I_y],[Z_y,0]], targety1phaseoperatoriY_x. Computationalu is conserved; eigensectorphases matter.

Algebraic_v0/old_source.json has15gates4CX02,02,10,12 and no u locals. Its first two controlledreflection axes areN1=(-X+Y)/sqrt2 andN3=-X. Their product N3N1=(I-iZ)/sqrt2=Rz_y(pi/2). Following exactT_x activephase e^iπ/4 givesS_y exactly, so mergedv0blockM=diag_x(I,S_y). The CH wrapper chronologicalRy_x(-π/4),H_x,CXux,H_x,Ry_x(π/4) has inactiveI andactiveH. FinalCZuy usesH_y,CXuy,H_y.

Foru0 allactionscommute withCZ_xy, which remains diagonal. Foru1 M W10 M† hasupperSdg andlowerS Z=Sdg; it isSdg_y X_x. CH turnsX_x intoZ_x; finalZ_y commutes withtargetphase. Exactdiagonal inxuyMSBorder is[1,1,1,-i,1,-1,-1,i]. No usectorscalar was discarded. Collapsing x/y locals intoSU2 maychange a single globalphase, but mustnot introduce orremove relativeusectorphase.

Tilting x control and adding later x readout escapes priorrestricted earlyCH proof, which assumed firstcomputationalxcontrol and everylatergate xblockdiagonal. The30coordinatefamily has only x/y localfreedom,sectorpreservingu; unrestrictedusectormixing is outside it. Fullcomparatorprefix4CX+reduced4CX+two barev->yCX costs10, but reducedv0validation imposes no v1solution. ExactcompletefullV4certification plus independentvalidation requiredbeforepromotion.

New reduced tilted-axis construction independentlyreviewed byhand: K=(X+Y)/sqrt2,J=(Y+Z)/sqrt2,D=(K+Z)/sqrt2,R=ZJZ,K0=ZDZ. DZD=K; xblocks afterD,CZxy,K0,CRux,CZxy areZK,I,ZJK,J for(u,y)=00,01,10,11. KXK=Y andJYJ=Z give the requiredI,Z,Z,iZ. FinalCZuy scalar- in11preservesconjugationroots. PredD=[1,1,1,i,1,-1,-1,-i].

Complete algebraic_v0/new_source.json24gates4CX chronology02,10,02,12 reviewed: exactU3(0,0,π)=Z preservesD/K0 phases; D polarπ4azimuthπ4,K0polarπ4azimuth-3π4. CR usesRxπ4 ZRx-π4=(Z-Y)/sqrt2. AdjacentH2 pair retainedcancels exactly. No matrices or exactarithmeticchecker ran in thisreview.

Particularnewliteral bare-v extension is staticallyinvalid: Ured commutesZy because onlynetyoperations areCZxy/CZuy. InsertingbareCXvy immediatelyafterfirstCXxy insideHy,CXxy,Hy producesZy^v afterthatCZ,whichcommutes throughremainingUred. EndbareCXvy givesFv1=Xy Zy Ured. RoutedWv1 mapsy->1-xory andhasinvertibleoffdiagonalxblocks; conjugatingbyy-monomialFv1 cannot remove them. This provesonlythisliteralbareextension failsv1,notany30coordinatefreefamily orphase-awarevcontrolledextension. Reducedcertificate is aseparatefact fromfullV4candidateexistence.
