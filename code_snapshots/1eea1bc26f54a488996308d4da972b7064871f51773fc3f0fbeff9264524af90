"""Bounded truncated-SVD Newton refinement; no exactness claim."""
import sys,json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
import numpy as np
import torch
from threadpoolctl import threadpool_limits
from check_circuit import evaluate,right_shift
from search13_adaptive_topology import matrix_for_schedule,serialize
from search13_joint_prefix_orbit_fit import unitary
from search13_prefix_perturb import decode_layers
from search13_fulltail_invariant_gauge_rank import circuit_matrix
HERE=Path(__file__).resolve().parent
SOURCE=ROOT/'search13_fresh_labelswap_refine_trial0.json'
def main():
 out=HERE/'labelswap_svd_refine_result.json';assert not out.exists();torch.set_num_threads(1)
 c=json.loads(SOURCE.read_text());x=decode_layers(c).ravel();schedule=tuple((g['control'],g['target']) for g in c['gates'] if g['gate']=='cx');fixed=[torch.eye(16,dtype=torch.complex128) for _ in range(14)];cx=matrix_for_schedule(schedule);vn=right_shift(4);v=torch.as_tensor(vn,dtype=torch.complex128);roots=np.array([1,1j,-1,-1j]);u=circuit_matrix(c['gates']);labels=np.argmin(abs(np.diag(u@vn@u.conj().T)[:,None]-roots[None,:]),axis=1);d=torch.as_tensor(roots[labels],dtype=torch.complex128)
 def residual(z):
  u=unitary(z.reshape(14,4,3),fixed,cx);r=(u@v-d[:,None]*u)/4;return torch.cat((r.real.ravel(),r.imag.ravel()))
 def array(z):return residual(torch.as_tensor(z,dtype=torch.float64)).detach().numpy()
 rows=[]
 with threadpool_limits(limits=1):
  for iteration in range(20):
   r=array(x);loss=float(r@r)
   if loss<1e-28:break
   t=torch.tensor(x,dtype=torch.float64,requires_grad=True);j=torch.autograd.functional.jacobian(residual,t,vectorize=True).detach().numpy()
   step,_,rank,singular=np.linalg.lstsq(j,-r,rcond=1e-6)
   accepted=False
   for k in range(8):
    alpha=2.**(-k);new=x+alpha*step;rr=array(new);newloss=float(rr@rr)
    if newloss<loss:x=new;accepted=True;break
   row={'iteration':iteration,'before':loss,'after':newloss,'rank':int(rank),'minimum_retained_singular':float(singular[rank-1]),'step_norm':float(np.linalg.norm(step)),'alpha':alpha,'accepted':accepted};rows.append(row);print(json.dumps(row),flush=True)
   if not accepted:break
 c=serialize([[] for _ in range(14)],schedule,x.reshape(14,4,3));p=HERE/'labelswap_svd_refine_candidate.json';p.write_text(json.dumps(c,indent=2)+'\n');q=evaluate(c);assert q['cnot_count']==13
 out.write_text(json.dumps({'source':str(SOURCE.relative_to(ROOT)),'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'labels':labels.tolist(),'schedule':schedule,'rcond':1e-6,'max_steps':20,'rows':rows,'candidate':str(p.relative_to(ROOT)),'checker':q,'scope':'Numerical13CNOT candidate; requires exact certificate'},indent=2)+'\n');print(json.dumps(q),flush=True)
if __name__=='__main__':main()
