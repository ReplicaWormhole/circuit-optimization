"""Change output eigenlabels by targeted pair swaps, then free-label polish.

All labels retain multiplicities. Fixed-label failures exclude no eigenbasis
or topology; only native13-CNOT saved gate lists are candidates.
"""
import hashlib,json
from pathlib import Path
import numpy as np
import torch
from scipy.optimize import minimize
from check_circuit import evaluate,right_shift
from search13_adaptive_topology import fit,matrix_for_schedule,serialize
from search13_joint_prefix_orbit_fit import unitary
from search13_prefix_perturb import decode_layers
from search13_fulltail_invariant_gauge_rank import circuit_matrix
ROOT=Path(__file__).resolve().parent
SOURCE=ROOT/'search13_fresh_chain_refine_trial1.json'
STEM='search13_fresh_chain_label_swaps'

def main():
 out=ROOT/f'{STEM}_result.json';plan=ROOT/f'{STEM}_plan.json';assert not out.exists() and not plan.exists()
 c=json.loads(SOURCE.read_text());x0=decode_layers(c);schedule=tuple((g['control'],g['target']) for g in c['gates'] if g['gate']=='cx');assert len(schedule)==13
 vnp=right_shift(4);u0=circuit_matrix(c['gates']);d0=u0@vnp@u0.conj().T;roots=np.array([1,1j,-1,-1j]);labels=np.argmin(abs(np.diag(d0)[:,None]-roots[None,:]),axis=1);assert np.bincount(labels,minlength=4).tolist()==[6,3,4,3]
 pairs=sorted([(max(abs(d0[a,b]),abs(d0[b,a])),a,b) for a in range(16) for b in range(a+1,16) if labels[a]!=labels[b]],key=lambda t:(-t[0],t[1],t[2]))[:6]
 proposals=[]
 for score,a,b in pairs:
  target=labels.copy();target[a],target[b]=target[b],target[a];proposals.append({'swapped_rows':[a,b],'source_max_transition':float(score),'labels':target.tolist()})
 plan.write_text(json.dumps({'source':SOURCE.name,'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'schedule':schedule,'nearest_labels':labels.tolist(),'proposals':proposals,'fixed_maxiter':400,'polish_maxiter':300},indent=2)+'\n')
 torch.set_num_threads(1);fixed=[torch.eye(16,dtype=torch.complex128) for _ in range(14)];cxs=matrix_for_schedule(schedule);v=torch.as_tensor(vnp,dtype=torch.complex128);rows=[]
 for case,p in enumerate(proposals):
  target=torch.as_tensor(roots[p['labels']],dtype=torch.complex128)
  def objective(flat):
   x=torch.tensor(np.asarray(flat).reshape(14,4,3),dtype=torch.float64,requires_grad=True);u=unitary(x,fixed,cxs);loss=(u@v-target[:,None]*u).abs().square().sum().real/16;grad=torch.autograd.grad(loss,x)[0];return float(loss.detach()),grad.detach().numpy().ravel().copy()
  initial,_=objective(x0.ravel());opt=minimize(objective,x0.ravel(),jac=True,method='L-BFGS-B',options={'maxiter':400,'ftol':1e-15,'gtol':1e-11,'maxls':30});x=opt.x.reshape(14,4,3)
  for stage in ('fixed','polish'):
   if stage=='polish':x,record=fit(schedule,x,fixed,v,300)
   else:record={'initial_loss':initial,'loss':float(opt.fun),'nit':int(opt.nit),'nfev':int(opt.nfev),'optimizer_success':bool(opt.success),'optimizer_message':str(opt.message)}
   candidate=serialize([[] for _ in range(14)],schedule,x);path=ROOT/f'{STEM}_case{case}_{stage}.json';assert not path.exists();path.write_text(json.dumps(candidate,indent=2)+'\n');q=evaluate(candidate);assert q['cnot_count']==13
   u=circuit_matrix(candidate['gates']);d=u@vnp@u.conj().T;ordinary=float(np.sum(abs(d-np.diag(np.diag(d)))**2)/16);label_loss=float(np.sum(abs(u@vnp-roots[p['labels'],None]*u)**2)/16)
   assert abs(record['loss']-(label_loss if stage=='fixed' else ordinary))<1e-11
   row={'case':case,'stage':stage,**p,**record,'ordinary_loss':ordinary,'fixed_label_loss':label_loss,'candidate':path.name,'cnot_count':q['cnot_count'],'off_diagonal_error':q['off_diagonal_error'],'valid_diagonalizer':q['valid_diagonalizer']};rows.append(row);print(json.dumps(row),flush=True)
 best=min(rows,key=lambda r:r['off_diagonal_error']);out.write_text(json.dumps({'plan':plan.name,'rows':rows,'best_candidate':best['candidate'],'scope':'six targeted label swaps at one13CNOT order, fixedlabel warmfits thenfree-labelpolish; no topology exclusion'},indent=2)+'\n')
if __name__=='__main__':main()
