"""Canonical run365 constants then pair endpoint, max8 bounded constrained fits."""
import sys,json,hashlib
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[3]))
import numpy as np,torch
from scipy.optimize import least_squares
from threadpoolctl import threadpool_limits
from check_circuit import right_shift,evaluate
from search13_fulltail_invariant_gauge_rank import circuit_matrix
from search13_adaptive_topology import matrix_for_schedule
from delete14_search import kron_gate
R=Path(__file__).resolve().parents[3];O=Path(__file__).resolve().parent;S=R/'collaboration/work/root/collective/canonical13_stripped_candidate.json';src=json.loads(S.read_text());x=np.array([g['theta'] for g in src['gates'] if g['gate']!='cx']);schedule=tuple((g['control'],g['target']) for g in src['gates'] if g['gate']=='cx');cs=matrix_for_schedule(schedule);v=torch.tensor(right_shift(4),dtype=torch.complex128);u=circuit_matrix(src['gates']);roots=np.array([1,1j,-1,-1j]);labels=np.argmin(abs(np.diag(u@v.numpy()@u.conj().T)[:,None]-roots),axis=1);d=torch.tensor(roots[labels]);torch.set_num_threads(1)
I=torch.eye(2,dtype=torch.complex128);P={'rx':torch.tensor([[0,1],[1,0]],dtype=torch.complex128),'ry':torch.tensor([[0,-1j],[1j,0]],dtype=torch.complex128),'rz':torch.diag(torch.tensor([1,-1],dtype=torch.complex128))}
def rt(t):
 u=torch.eye(16,dtype=torch.complex128);a=b=0
 for g in src['gates']:
  if g['gate']=='cx':u=cs[b]@u;b+=1
  else:
   m=torch.cos(t[a]/2)*I-1j*torch.sin(t[a]/2)*P[g['gate']];u=kron_gate(m,g['qubit'])@u;a+=1
 r=(u@v-d[:,None]*u)/4
 return torch.cat((r.real.ravel(),r.imag.ravel()))
def res(t):return rt(torch.tensor(t)).detach().numpy()
def jac(t):return torch.autograd.functional.jacobian(rt,torch.tensor(t,requires_grad=True),vectorize=True).detach().numpy()
def cand(t,rounded=False):
 out=[];a=0
 for g in src['gates']:
  h=dict(g)
  if h['gate']!='cx':h['theta']=f'{int(round(t[a]/np.pi*24))}*pi/24' if rounded else float(t[a]);a+=1
  out.append(h)
 return {'n':4,'gates':out}
fixed={i:float(round(t/np.pi*24)*np.pi/24) for i,t in enumerate(x) if abs(t-round(t/np.pi*24)*np.pi/24)<1e-10};rows=[];mapping=[{'coordinate':i,'gate':g['gate'],'qubit':g['qubit']} for i,g in enumerate([g for g in src['gates'] if g['gate']!='cx'])];snapshot={'source':str(S.relative_to(R)),'source_sha256':hashlib.sha256(S.read_bytes()).hexdigest(),'schedule':schedule,'mapping':mapping,'initial_fixed':fixed,'source_parameters':x.tolist(),'labels':labels.tolist()};(O/'canonical_source_map.json').write_text(json.dumps(snapshot,indent=2)+'\n')
with threadpool_limits(limits=1):
 for stage in range(8):
  trial=dict(fixed)
  if stage==0:purpose='recognized_constants'
  elif stage==1:trial[3]=0.;trial[5]=np.pi;purpose='paired_initial_endpoint_ry_q1_and_q2'
  else:
   free=np.array([i for i in range(len(x)) if i not in fixed]);j=jac(x)[:,free];_,s,vh=np.linalg.svd(j,full_matrices=True);rank=sum(s>1e-8);N=vh[rank:].T;diag=np.sum(N*N,axis=1);delta=np.round(x[free]/np.pi*24)*np.pi/24-x[free];score=abs(delta)/np.sqrt(np.maximum(diag,1e-30));score[diag<1e-7]=np.inf
   if not np.isfinite(score).any():break
   q=int(free[np.argmin(score)]);trial[q]=float(round(x[q]/np.pi*24)*np.pi/24);purpose=f'nullspace_snap_{q}'
  free=np.array([i for i in range(len(x)) if i not in trial]);base=x.copy()
  for i,val in trial.items():base[i]=val
  def expand(y):z=base.copy();z[free]=y;return z
  fit=least_squares(lambda y:res(expand(y)),base[free],jac=lambda y:jac(expand(y))[:,free],max_nfev=400,ftol=1e-13,xtol=1e-13,gtol=1e-13);loss=float(2*fit.cost);ok=loss<1e-20;row={'stage':stage,'purpose':purpose,'fixed_count':len(trial),'nfev':fit.nfev,'loss':loss,'accepted':ok};rows.append(row);print(json.dumps(row),flush=True)
  if ok:x=expand(fit.x);fixed=trial;(O/f'canonical_stage{stage}.json').write_text(json.dumps(cand(x),indent=2)+'\n')
  elif stage!=1:break
(O/'canonical_gauge_candidate.json').write_text(json.dumps(cand(x),indent=2)+'\n');(O/'canonical_rounded_candidate.json').write_text(json.dumps(cand(x,True),indent=2)+'\n');result={'source':snapshot,'rows':rows,'fixed':fixed,'numeric_check':evaluate(cand(x)),'rounded_check':evaluate(cand(x,True)),'limitation':'Bounded recognition only; no exact claim.'};(O/'canonical_result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result['numeric_check']),flush=True);print(json.dumps(result['rounded_check']),flush=True)
