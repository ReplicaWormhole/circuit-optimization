"""Exact rational sign checks for real, nonzero projective local rotations.

This verifies local factors only, not full diagonalization.
"""
import json
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[3]
SOURCE=ROOT/'collaboration/work/algebraic13_new/canonical_algebraic_ansatz.json'
OUT=Path(__file__).with_name('radical13_local_exact_review.json')
def sign_qsqrt2(expr):
 rt=s.sqrt(2);e=s.expand(s.simplify(expr));b=e.coeff(rt);a=s.simplify(e-b*rt)
 assert a.is_Rational and b.is_Rational,(e,a,b)
 if not b:return s.sign(a)
 if a>=0 and b>0:return 1
 if a<=0 and b<0:return -1
 if a>0:return s.sign(a*a-2*b*b)
 return s.sign(2*b*b-a*a)
def main():
 c=json.loads(SOURCE.read_text());proof=[]
 for g in c['gates']:
  if g['gate']=='cx':continue
  si=s.sympify(g['sin_exact']);co=s.sympify(g['cos_exact']);radicands=[]
  for e in (si,co):
   for node in s.preorder_traversal(e):
    if isinstance(node,s.Pow) and node.exp==s.Rational(1,2):
     assert sign_qsqrt2(node.base)>0;radicands.append(str(node.base))
  assert s.simplify(si**2+co**2-1)==0
  if s.simplify(1+co)==0:a,b=0,1;reason='purePauli'
  elif si==0:assert co==1;a,b=2,0;reason='factor2'
  else:assert sign_qsqrt2(si**2)>0;a,b=1+co,-s.I*si;reason='real normalized cosine strictly between-1 and1'
  if g['gate']=='rx':m=((a,b),(b,a))
  elif g['gate']=='ry':m=((a,-s.I*b),(s.I*b,a))
  elif g['gate']=='rz':m=((a+b,0),(0,a-b))
  else:raise ValueError(g)
  den=s.lcm([s.denom(s.radsimp(s.sympify(z))) for row in m for z in row]);assert den.is_Integer and den>0,den
  proof.append({'coordinate':g['coordinate'],'positive_radicands':radicands,'projective_nonzero_reason':reason,'cleared_denominator':int(den)})
 assert len(proof)==60
 OUT.write_text(json.dumps({'source':str(SOURCE.relative_to(ROOT)),'rotations_checked':len(proof),'sign_method':'Exact rational comparisons for a+b sqrt2, squaring only after signs fixed; no floating sign tests','checks':proof,'scope':'Local realness normalization andnonzero projective scaling only'},indent=2)+'\n')
 print('All60local exact conditions verified')
if __name__=='__main__':main()
