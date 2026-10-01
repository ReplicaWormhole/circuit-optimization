"""Continuously relocate one CNOT, count only native endpoints.

C(a)=(I+C)/2+exp(i*pi*a)(I-C)/2. At a=1 it is CNOT, at a=0 identity.
The intermediate product involves non-native powered entanglers and is
not a13CNOT circuit. Endpoint candidates contain exactly13 ordinary CNOTs.
"""
import argparse,json
from pathlib import Path
import numpy as np
import torch
from check_circuit import evaluate,right_shift
from search13_adaptive_topology import matrix_for_schedule,objective,serialize,fit as direct_fit
from search13_joint_prefix_orbit_fit import fit
from search13_prefix_perturb import decode_layers
ROOT=Path(__file__).resolve().parent
SOURCE=ROOT/'search13_fresh_chain_refine_trial1.json'
torch.set_num_threads(1)

def power(c,a):
    identity=torch.eye(16,dtype=torch.complex128)
    return (identity+c)/2+np.exp(1j*np.pi*a)*(identity-c)/2

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--maxiter',type=int,default=150);p.add_argument('--polish-maxiter',type=int,default=300);args=p.parse_args()
    out=ROOT/'search13_fresh_chain_edge_homotopy_result.json';assert not out.exists(),out
    source=json.loads(SOURCE.read_text());x0=decode_layers(source);schedule=tuple((g['control'],g['target']) for g in source['gates'] if g['gate']=='cx');assert len(schedule)==13
    base=matrix_for_schedule(schedule);fixed=[torch.eye(16,dtype=torch.complex128) for _ in range(14)];chunks=[[] for _ in range(14)];v=torch.as_tensor(right_shift(4),dtype=torch.complex128)
    initial_loss=objective(x0.ravel(),fixed,base,v,False);assert abs(initial_loss-.0012158342275705642)<1e-10
    rows=[]
    for case,(slot,new_edge) in enumerate(((6,(0,3)),(8,(0,2)))):
        new=matrix_for_schedule([new_edge])[0];old=base[slot]
        assert torch.max(abs(power(old,1)-old))<1e-14 and torch.max(abs(power(old,0)-torch.eye(16)))<1e-14
        angles=x0.copy();stages=[]
        for alpha in (1.,.875,.75,.5,.25,.125,0.):
            interpolation=power(new,1-alpha)@power(old,alpha)
            assert torch.max(abs(interpolation.conj().T@interpolation-torch.eye(16)))<1e-13
            cxs=list(base);cxs[slot]=interpolation
            before=objective(angles.ravel(),fixed,cxs,v,False)
            angles,record=fit(lambda x:objective(x,fixed,cxs,v,True),angles,args.maxiter)
            stage={'alpha':alpha,'initial_loss':before,**record,'native_single_cnot_endpoint':alpha in (0.,1.)};stages.append(stage)
            print(json.dumps({'case':case,'slot':slot,'alpha':alpha,'loss':record['loss'],'native_endpoint':alpha in (0.,1.)}),flush=True)
        endpoint=list(schedule);endpoint[slot]=new_edge
        # At alpha0 the old powered interaction is identity and the new one a CNOT.
        assert torch.max(abs(cxs[slot]-new))<1e-14
        candidate=serialize(chunks,endpoint,angles);path=ROOT/f'search13_fresh_chain_edge_homotopy_case{case}_endpoint.json';path.write_text(json.dumps(candidate,indent=2)+'\n');q=evaluate(candidate);assert q['cnot_count']==13
        polished,record=direct_fit(endpoint,angles,fixed,v,args.polish_maxiter)
        candidate2=serialize(chunks,endpoint,polished);path2=ROOT/f'search13_fresh_chain_edge_homotopy_case{case}_polish.json';path2.write_text(json.dumps(candidate2,indent=2)+'\n');check=evaluate(candidate2);assert check['cnot_count']==13
        row={'case':case,'slot':slot,'old_edge':schedule[slot],'new_edge':new_edge,'stages':stages,'endpoint_candidate':path.name,'endpoint_off_diagonal_error':q['off_diagonal_error'],'endpoint_valid':q['valid_diagonalizer'],'polish':record,'candidate':path2.name,'off_diagonal_error':check['off_diagonal_error'],'valid_diagonalizer':check['valid_diagonalizer'],'cnot_count':check['cnot_count']};rows.append(row)
        print(json.dumps({'case':case,'polish_loss':record['loss'],'maxoff':check['off_diagonal_error'],'valid':check['valid_diagonalizer']}),flush=True)
    best=min(rows,key=lambda r:r['off_diagonal_error'])
    out.write_text(json.dumps({'source':SOURCE.name,'source_loss':initial_loss,'source_schedule':schedule,'maxiter':args.maxiter,'polish_maxiter':args.polish_maxiter,'rows':rows,'best_candidate':best['candidate'],'scope':'Two specified powered-entangler homotopies; only native endpoint/polish gate lists count as13CNOT candidates; bounded failures exclude nothing globally'},indent=2)+'\n')
if __name__=='__main__':main()
