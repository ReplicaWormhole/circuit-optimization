"""Bounded Euler-coordinate gauge freezing of source run 359; numerical only until exact check."""
import sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[3]))
import numpy as np,torch
from scipy.optimize import least_squares
from threadpoolctl import threadpool_limits
from check_circuit import right_shift,evaluate
from search13_adaptive_topology import matrix_for_schedule
from search13_fulltail_invariant_gauge_rank import circuit_matrix
from delete14_search import kron_gate
ROOT=Path(__file__).resolve().parents[3]; OUT=Path(__file__).resolve().parent
torch.set_num_threads(1)
source=json.loads((ROOT/'collaboration/work/root/labelswap_svd_refine_candidate.json').read_text())
x=np.array([[g[k] for k in ('theta','phi','lam')] for g in source['gates'] if g['gate']=='u3']).ravel()
schedule=tuple((g['control'],g['target']) for g in source['gates'] if g['gate']=='cx');cxs=matrix_for_schedule(schedule)
v=torch.tensor(right_shift(4),dtype=torch.complex128);u0=circuit_matrix(source['gates']);roots=np.array([1,1j,-1,-1j]);labels=np.argmin(abs(np.diag(u0@v.numpy()@u0.conj().T)[:,None]-roots),axis=1);d=torch.tensor(roots[labels])
def rt(t):
 u=torch.eye(16,dtype=torch.complex128)
 for slot in range(14):
  for wire in range(4):
   th,ph,la=t.reshape(14,4,3)[slot,wire];c=torch.cos(th/2);s=torch.sin(th/2)
   g=torch.stack((torch.stack((c,-torch.exp(1j*la)*s)),torch.stack((torch.exp(1j*ph)*s,torch.exp(1j*(ph+la))*c))))
   u=kron_gate(g,wire)@u
  if slot<13:u=cxs[slot]@u
 r=(u@v-d[:,None]*u)/4
 return torch.cat((r.real.ravel(),r.imag.ravel()))
def residual(z):return rt(torch.tensor(z)).detach().numpy()
def jac(z):return torch.autograd.functional.jacobian(rt,torch.tensor(z,requires_grad=True),vectorize=True).detach().numpy()
def candidate(z,exact=False):
 out=[];p=0
 for g in source['gates']:
  h=dict(g)
  if h['gate']=='u3':
   for k in ('theta','phi','lam'):
    h[k]= f'{int(round(z[p]/np.pi*24))}*pi/24' if exact else float(z[p]);p+=1
  out.append(h)
 return {'n':4,'gates':out}
fixed={};rows=[]
with threadpool_limits(limits=1):
 for stage in range(12):
  free=np.array([i for i in range(168) if i not in fixed]);j=jac(x)[:,free];_,sv,vh=np.linalg.svd(j,full_matrices=True);rank=int(np.sum(sv>1e-8));null=vh[rank:].T;P=null@null.T;selected=[]
  for _ in range(min(9,len(free)-rank)):
   delta=np.round(x[free]/np.pi*24)*np.pi/24-x[free];sens=np.diag(P);score=abs(delta)/np.sqrt(np.maximum(sens,1e-30));score[sens<1e-5]=np.inf
   if not np.isfinite(score).any():break
   q=int(np.argmin(score));selected.append(int(free[q]));col=P[:,q].copy();P-=np.outer(col,col)/P[q,q]
  if not selected:break
  trial=dict(fixed);trial.update({i:float(round(x[i]/np.pi*24)*np.pi/24) for i in selected});remain=np.array([i for i in range(168) if i not in trial]);base=x.copy()
  for i,val in trial.items():base[i]=val
  def expand(y):z=base.copy();z[remain]=y;return z
  fit=least_squares(lambda y:residual(expand(y)),base[remain],jac=lambda y:jac(expand(y))[:,remain],max_nfev=400,ftol=1e-13,xtol=1e-13,gtol=1e-13)
  loss=float(2*fit.cost);ok=loss<1e-20
  rows.append({'stage':stage,'rank':rank,'added':selected,'fixed':len(trial),'nfev':fit.nfev,'loss':loss,'accepted':ok});print(json.dumps(rows[-1]),flush=True)
  if not ok:break
  x=expand(fit.x);fixed=trial
  (OUT/f'fine24_stage{stage}.json').write_text(json.dumps(candidate(x),indent=2)+'\n')
(OUT/'fine24_gauge_candidate.json').write_text(json.dumps(candidate(x),indent=2)+'\n');rounded=candidate(x,True);(OUT/'fine24_rounded_candidate.json').write_text(json.dumps(rounded,indent=2)+'\n')
result={'rows':rows,'fixed_count':len(fixed),'labels':labels.tolist(),'numeric_check':evaluate(candidate(x)),'rounded_check':evaluate(rounded),'limitations':'Gauge snap search; numerical fit does not certify exactness.'};(OUT/'fine24_result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
