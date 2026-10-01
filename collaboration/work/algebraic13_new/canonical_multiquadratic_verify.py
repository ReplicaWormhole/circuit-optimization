"""Exact Fraction multiquadratic algebra; dictionary zero is a sufficient certificate."""
import json,time,signal,sys,hashlib,math
from fractions import Fraction as F
from pathlib import Path
import sympy as sp
O=Path(__file__).resolve().parent;S=O/'canonical_algebraic_ansatz.json';start=time.monotonic();state={'source_sha256':hashlib.sha256(S.read_bytes()).hexdigest(),'exact_diagonalizer':None};basis=[];products=[(F(1),F(0))];proofs=[]
def save(**kw):
 state.update(kw);state['elapsed']=time.monotonic()-start;state['basis_rank']=len(basis);(O/'canonical_multiquadratic_result.json').write_text(json.dumps(state,indent=2)+'\n');print(json.dumps(state),flush=True)
def alarm(sig,frame):save(stage='timeout',limitation='Incomplete computation; no certificate');sys.exit(2)
signal.signal(signal.SIGALRM,alarm);signal.alarm(300)
def ka(x,y):return x[0]+y[0],x[1]+y[1]
def kn(x):return -x[0],-x[1]
def km(x,y):return x[0]*y[0]+2*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def kd(x,y):
 norm=y[0]*y[0]-2*y[1]*y[1];assert norm
 return km(x,(y[0]/norm,-y[1]/norm))
def positive(x):
 a,b=x
 if b==0:return a>0
 if a==0:return b>0
 if a>0 and b>0:return True
 if a<0 and b<0:return False
 return (a*a>2*b*b) if a>0 else (2*b*b>a*a)
def rsqrt(x):
 if x<0:return None
 a,b=math.isqrt(x.numerator),math.isqrt(x.denominator)
 return F(a,b) if a*a==x.numerator and b*b==x.denominator else None
def ksqrt(x):
 a,b=x
 if b==0:
  u=rsqrt(a)
  if u is not None:return(u,F(0))
  v=rsqrt(a/2)
  if v is not None:return(F(0),v)
  return None
 rho=rsqrt(a*a-2*b*b)
 if rho is None:return None
 for z in (rho,-rho):
  u=rsqrt((a+z)/2)
  if u is None or u==0:continue
  t=(u,b/(2*u))
  if km(t,t)==x:return t if positive(t) else kn(t)
 return None
def js(x):return [str(z) for z in x]
ZERO=(F(0),)*4;ONE=(F(1),F(0),F(0),F(0));IMAG=(F(0),F(0),F(1),F(0))
def ca(x,y):return tuple(a+b for a,b in zip(x,y))
def cn(x):return tuple(-z for z in x)
def cm(x,y):
 r=ka(km(x[:2],y[:2]),kn(km(x[2:],y[2:])));i=ka(km(x[:2],y[2:]),km(x[2:],y[:2]));return r+i
def ck(x,k):return km(x[:2],k)+km(x[2:],k)
def add(x,y):
 out=dict(x)
 for mask,c in y.items():
  z=ca(out.get(mask,ZERO),c)
  if z==ZERO:out.pop(mask,None)
  else:out[mask]=z
 return out
def neg(x):return {m:cn(c) for m,c in x.items()}
def times(x,y):
 out={}
 for m,c in x.items():
  for n,d in y.items():
   p=m^n;v=ck(cm(c,d),products[m&n]);z=ca(out.get(p,ZERO),v)
   if z==ZERO:out.pop(p,None)
   else:out[p]=z
 return out
def scalar(c):return {} if c==ZERO else {0:c}
def phase(x):return {m:cm(c,IMAG) for m,c in x.items()}
def trig(expr,coord,name):
 e=sp.sympify(expr)
 if e==0:proofs.append({'coordinate':coord,'name':name,'zero':True});return {}
 sq=sp.expand(sp.simplify(e**2));a=sp.simplify(sq.subs(sp.sqrt(2),0));b=sp.simplify((sq-a)/sp.sqrt(2));assert a.is_Rational and b.is_Rational,(coord,expr,sq,a,b)
 q=(F(str(a)),F(str(b)));assert positive(q);sgn=1 if e.is_positive is True else -1 if e.is_negative is True else None;assert sgn is not None,(coord,e)
 for mask,p in enumerate(products):
  root=ksqrt(kd(q,p))
  if root is not None:
   assert km(km(root,root),p)==q
   proofs.append({'coordinate':coord,'name':name,'square':js(q),'sign':sgn,'mask':mask,'K_factor':js(root),'relation_verified':True});return{mask:(sgn*root[0],sgn*root[1],F(0),F(0))}
 mask=1<<len(basis);basis.append(q);products.extend([km(p,q) for p in list(products)]);proofs.append({'coordinate':coord,'name':name,'square':js(q),'sign':sgn,'mask':mask,'K_factor':['1','0'],'new_generator':True});save(stage='basis_construction',last_coordinate=coord);return{mask:(F(sgn),F(0),F(0),F(0))}
c=json.loads(S.read_text());gates=[];unit={0:ONE}
for g in c['gates']:
 if g['gate']=='cx':gates.append((g,None));continue
 si=trig(g['sin_exact'],g['coordinate'],'sin');co=trig(g['cos_exact'],g['coordinate'],'cos');assert add(times(si,si),times(co,co))==unit
 if add(unit,co)=={}:a={};b=unit
 else:a=add(unit,co);b=neg(phase(si))
 if g['gate']=='rx':m=((a,b),(b,a))
 elif g['gate']=='ry':m=((a,neg(phase(b))),(phase(b),a))
 else:m=((add(a,b),{}),({},add(a,neg(b))))
 gates.append((g,m))
(O/'canonical_multiquadratic_basis.json').write_text(json.dumps({'basis_radicands':[js(q) for q in basis],'representation_proofs':proofs,'source_sha256':state['source_sha256'],'semantics':'Positive formal square roots; even if dependent, a zero dictionary guarantees exact zero.'},indent=2)+'\n')
save(stage='basis_complete',basis_dimension=2**len(basis),cnot_count=sum(g['gate']=='cx' for g in c['gates']),local_normalizations_exact=True)
u=[[unit.copy() if r==j else {} for j in range(16)] for r in range(16)]
for pos,(g,m) in enumerate(gates):
 if g['gate']=='cx':
  cmask,tmask=1<<(3-g['control']),1<<(3-g['target']);v=[None]*16
  for r in range(16):v[r^tmask if r&cmask else r]=u[r]
  u=v
 else:
  mask=1<<(3-g['qubit'])
  for r in range(16):
   if r&mask:continue
   t=r^mask;x,y=u[r],u[t];u[r]=[add(times(m[0][0],a),times(m[0][1],b)) for a,b in zip(x,y)];u[t]=[add(times(m[1][0],a),times(m[1][1],b)) for a,b in zip(x,y)]
 if pos%5==0:save(stage='assembly',gate_position=pos,max_terms=max(len(z) for row in u for z in row))
roots=[ONE,IMAG,cn(ONE),cn(IMAG)];bad=[]
for r,row in enumerate(u):
 for j in range(16):
  z=add(row[((j&1)<<3)|(j>>1)],neg(times(scalar(roots[c['labels'][r]]),row[j])))
  if z:bad.append({'row':r,'column':j,'terms':{str(m):[str(v) for v in a] for m,a in z.items()}})
(O/'canonical_multiquadratic_residuals.json').write_text(json.dumps(bad,indent=2)+'\n')
save(stage='certificate' if not bad else 'unreduced_residuals',exact_diagonalizer=True if not bad else None,entries_checked=256,nonempty_residual_dictionaries=len(bad),limitation=None if not bad else 'Nonzero dictionaries can represent zero if generators dependent; no disproof.')
