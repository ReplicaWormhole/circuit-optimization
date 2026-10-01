"""Bounded standalone exact radical projective verification; no floating-point zero tests."""
import json,time,signal,sys
from pathlib import Path
import sympy as s
O=Path(__file__).resolve().parent;start=time.monotonic();progress={'stage':'initializing','elapsed':0,'exact_diagonalizer':None}
def save(**kw):
 progress.update(kw);progress['elapsed']=time.monotonic()-start;(O/'canonical_projective_progress.json').write_text(json.dumps(progress,indent=2)+'\n');print(json.dumps(progress),flush=True)
def timeout(sig,frame):save(stage='timeout',limitation='Incomplete exact verification; no certificate');sys.exit(2)
signal.signal(signal.SIGALRM,timeout);signal.alarm(300)
c=json.loads((O/'canonical_algebraic_ansatz.json').read_text());I=s.I
# A= (1+cos(theta))I-i sin(theta)P = 2cos(theta/2)R(theta).
# When cos=-1, use P instead (R is a nonzero scalar multiple).
# All supplied real sin/cos satisfy sin²+cos²=1 exactly; normalization separately checked.
mat=[]
for g in c['gates']:
 if g['gate']=='cx':mat.append(None);continue
 si=s.sympify(g['sin_exact']);co=s.sympify(g['cos_exact']);assert s.simplify(si**2+co**2-1)==0
 if s.simplify(1+co)==0:a,b=s.Integer(0),s.Integer(1)
 else:a,b=1+co,-I*si
 if g['gate']=='rx':m=((a,b),(b,a))
 elif g['gate']=='ry':m=((a,-I*b),(I*b,a))
 else:m=((a+b,s.Integer(0)),(s.Integer(0),a-b))
 # Remove common scalar denominators to avoid growing rational denominators.
 den=s.lcm([s.denom(s.radsimp(z)) for row in m for z in row]);m=tuple(tuple(s.expand(s.radsimp(z*den)) for z in row) for row in m);mat.append(m)
save(stage='local_unitarity_checked',rotations=sum(m is not None for m in mat),cnot_count=sum(g['gate']=='cx' for g in c['gates']))
u=[[s.Integer(r==j) for j in range(16)] for r in range(16)]
def reduce(e):return s.sqrtdenest(s.radsimp(s.expand(e)))
for pos,(g,m) in enumerate(zip(c['gates'],mat)):
 if g['gate']=='cx':
  cm,tm=1<<(3-g['control']),1<<(3-g['target']);v=[None]*16
  for r in range(16):v[r^tm if r&cm else r]=u[r]
  u=v
 else:
  mask=1<<(3-g['qubit'])
  for r in range(16):
   if r&mask:continue
   t=r^mask;x,y=u[r],u[t];u[r]=[reduce(m[0][0]*a+m[0][1]*b) for a,b in zip(x,y)];u[t]=[reduce(m[1][0]*a+m[1][1]*b) for a,b in zip(x,y)]
 if pos%5==0:save(stage='assembling',gate_position=pos,max_expression_ops=max(s.count_ops(z) for row in u for z in row))
save(stage='checking_equation')
roots=[1,I,-1,-I];checked=0
for r,row in enumerate(u):
 for j in range(16):
  e=row[((j&1)<<3)|(j>>1)]-roots[c['labels'][r]]*row[j];z=s.simplify(reduce(e))
  if z!=0:save(stage='nonzero_or_unsimplified_residual',row=r,column=j,residual=str(z),exact_diagonalizer=False,limitation='A nonzero symbolic form may require further reduction');sys.exit(1)
  checked+=1
 save(stage='checking_equation',rows_checked=r+1,entries_checked=checked)
save(stage='certified',exact_diagonalizer=True,entries_checked=checked,interpretation='Projective UV4=DU exactly; each local normalized rotation is unitary because exact sin²+cos²=1 and prescribed real branches. All projective scalar factors nonzero by cos != -1 or explicit P branch.')
