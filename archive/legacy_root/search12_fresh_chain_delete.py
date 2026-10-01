"""Two deletions from fresh chain13 lead, warm and perturbed full-local fits.

Keep 14 local layers (two become adjacent at the removed interaction).
Their product is arbitrary local, so the extra layer is redundant and
adds no entangling resource. This is a 12-CNOT ansatz, not a global screen.
"""
import argparse,json
from pathlib import Path
import numpy as np
import torch
from scipy.optimize import minimize
from check_circuit import evaluate,right_shift
from search13_adaptive_topology import matrix_for_schedule,objective
from search13_prefix_perturb import decode_layers
from delete14_search import axis_to_u3
ROOT=Path(__file__).resolve().parent
SOURCE=ROOT/'search13_fresh_chain_refine_trial1.json'
torch.set_num_threads(1)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--maxiter',type=int,default=400);p.add_argument('--seed',type=int,default=34200);args=p.parse_args()
    out=ROOT/'search12_fresh_chain_delete_result.json';assert not out.exists(),out
    source=json.loads(SOURCE.read_text());x0=decode_layers(source);schedule=[(g['control'],g['target']) for g in source['gates'] if g['gate']=='cx'];assert len(schedule)==13
    fixed=[torch.eye(16,dtype=torch.complex128) for _ in range(14)];base=matrix_for_schedule(schedule);v=torch.as_tensor(right_shift(4),dtype=torch.complex128)
    assert abs(objective(x0.ravel(),fixed,base,v,False)-.0012158342275705642)<1e-10
    rows=[]
    for deletion in (0,12):
        cxs=list(base);cxs[deletion]=torch.eye(16,dtype=torch.complex128)
        for start,sigma in enumerate((0.,.3)):
            seed=args.seed+10*deletion+start
            initial=x0+np.random.default_rng(seed).normal(0,sigma,x0.shape)
            initial_loss=objective(initial.ravel(),fixed,cxs,v,False)
            result=minimize(lambda x:objective(x,fixed,cxs,v,True),initial.ravel(),jac=True,method='L-BFGS-B',options={'maxiter':args.maxiter,'ftol':1e-15,'gtol':1e-11,'maxls':30})
            angles=result.x.reshape(14,4,3);gates=[]
            for slot in range(14):
                for wire in range(4):gates.append({'gate':'u3','qubit':wire,**axis_to_u3(*angles[slot,wire])})
                if slot<13 and slot!=deletion:gates.append({'gate':'cx','control':schedule[slot][0],'target':schedule[slot][1]})
            candidate={'n':4,'gates':gates};path=ROOT/f'search12_fresh_chain_delete_slot{deletion}_start{start}.json';path.write_text(json.dumps(candidate,indent=2)+'\n')
            check=evaluate(candidate);assert check['cnot_count']==12
            assert abs(float(result.fun)-objective(angles.ravel(),fixed,cxs,v,False))<1e-10
            row={'deletion':deletion,'start':start,'sigma':sigma,'seed':seed,'initial_loss':initial_loss,'loss':float(result.fun),'nit':int(result.nit),'nfev':int(result.nfev),'optimizer_success':bool(result.success),'candidate':path.name,'cnot_count':check['cnot_count'],'off_diagonal_error':check['off_diagonal_error'],'valid_diagonalizer':check['valid_diagonalizer']};rows.append(row);print(json.dumps(row),flush=True)
    best=min(rows,key=lambda r:r['off_diagonal_error'])
    out.write_text(json.dumps({'source':SOURCE.name,'source_schedule':schedule,'seed':args.seed,'maxiter':args.maxiter,'rows':rows,'best_candidate':best['candidate'],'scope':'Two specified12-CNOT deletions, two numerical starts each; arbitrary local gates and one redundant adjacent layer; no global exclusion'},indent=2)+'\n')
if __name__=='__main__':main()
