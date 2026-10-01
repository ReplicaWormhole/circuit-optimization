"""Independent exact positive-root relation audit over Q(sqrt2)."""
import json,hashlib,math
from fractions import Fraction as F
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[3];SOURCE=ROOT/'collaboration/work/algebraic13_new/canonical_algebraic_ansatz.json';OUT=Path(__file__).with_name('radical13_tower_relations_result.json')
def mul(x,y):return x[0]*y[0]+2*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def div(x,y):
 norm=y[0]*y[0]-2*y[1]*y[1];assert norm
 return mul(x,(y[0]/norm,-y[1]/norm))
def sqrt_rat(x):
 if x<0:return None
 a,b=math.isqrt(x.numerator),math.isqrt(x.denominator)
 return F(a,b) if a*a==x.numerator and b*b==x.denominator else None
def sign(x):
 a,b=x
 if not b:return (a>0)-(a<0)
 if a>=0 and b>0:return 1
 if a<=0 and b<0:return -1
 z=a*a-2*b*b if a>0 else 2*b*b-a*a
 return (z>0)-(z<0)
def sqrt_k(x):
 a,b=x
 if not b:
  r=sqrt_rat(a)
  if r is not None:return r,F(0)
  r=sqrt_rat(a/2)
  return (F(0),r) if r is not None else None
 n=sqrt_rat(a*a-2*b*b)
 if n is None:return None
 for t in (n,-n):
  u=sqrt_rat((a+t)/2)
  if u:
   result=(u,b/(2*u));assert mul(result,result)==x
   return result if sign(result)>0 else (-result[0],-result[1])
 return None
def pair(expr):
 e=s.expand(s.simplify(expr));b=e.coeff(s.sqrt(2));a=s.simplify(e-b*s.sqrt(2));assert a.is_Rational and b.is_Rational
 return F(int(a.p),int(a.q)),F(int(b.p),int(b.q))
def encoded(x):return [str(t) for t in x]
def main():
 assert not OUT.exists();c=json.loads(SOURCE.read_text());basis=[];products=[(F(1),F(0))];rows=[];cache={}
 for gate in c['gates']:
  if gate['gate']=='cx':continue
  for trig in ('sin_exact','cos_exact'):
   q=pair(s.sympify(gate[trig])**2)
   if q not in cache:
    found=None
    for mask,product in enumerate(products):
     k=sqrt_k(div(q,product))
     if k is not None:found=(mask,k);break
    if found is None:
     assert sign(q)>0;mask=1<<len(basis);basis.append(q);products.extend([mul(p,q) for p in products.copy()]);found=mask,(F(1),F(0))
    mask,k=found;assert mul(mul(k,k),products[mask])==q
    assert sign(k)>0 or q==(0,0);cache[q]=found
   mask,k=cache[q];rows.append({'coordinate':gate['coordinate'],'trig':trig,'squared_value':encoded(q),'mask':mask,'positive_coefficient':encoded(k)})
 OUT.write_text(json.dumps({'source':str(SOURCE.relative_to(ROOT)),'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'basis_radicands':[encoded(q) for q in basis],'generator_count':len(basis),'distinct_squared_values':len(cache),'relations':rows,'scope':'Exact positive-root relations; independence of generators not asserted; no fullcircuit proof'},indent=2)+'\n');print('Generator count',len(basis),'distinct squared values',len(cache),'basis',[encoded(x) for x in basis])
if __name__=='__main__':main()
