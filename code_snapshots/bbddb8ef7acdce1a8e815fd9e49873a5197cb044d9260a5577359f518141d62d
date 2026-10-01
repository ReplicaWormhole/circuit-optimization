"""Multistarts from fitted tournament's alternative opening02 topology.

This source has smaller max entry but worse ordinary loss than chain13;
all figures remain numerical. No exact or topology-exclusion claim.
"""
import argparse,json
from pathlib import Path
import numpy as np
import torch
from check_circuit import evaluate,right_shift
from search13_adaptive_topology import fit,matrix_for_schedule,objective,serialize
from search13_prefix_perturb import decode_layers
ROOT=Path(__file__).resolve().parent
SOURCE=ROOT/'search13_fresh_chain_fitted_tournament_p0_deep.json'
torch.set_num_threads(1)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--seed',type=int,default=35000);p.add_argument('--maxiter',type=int,default=650);args=p.parse_args()
    out=ROOT/'search13_alternative_opening_refine_result.json';assert not out.exists(),out
    source=json.loads(SOURCE.read_text());x0=decode_layers(source);schedule=tuple((g['control'],g['target']) for g in source['gates'] if g['gate']=='cx');assert len(schedule)==13 and set(schedule[0])=={0,2}
    prior=json.load(open(ROOT/'search13_fresh_chain_fitted_tournament_result.json'));expected=next(r['loss'] for r in prior['deep_rows'] if r['proposal']==0)
    fixed=[torch.eye(16,dtype=torch.complex128) for _ in range(14)];chunks=[[] for _ in range(14)];v=torch.as_tensor(right_shift(4),dtype=torch.complex128)
    loss=objective(x0.ravel(),fixed,matrix_for_schedule(schedule),v,False);assert abs(loss-expected)<1e-10
    rows=[]
    for trial,sigma in enumerate((0.,.05,.2,.6)):
        seed=args.seed+trial;initial=x0+np.random.default_rng(seed).normal(0,sigma,x0.shape)
        angles,record=fit(schedule,initial,fixed,v,args.maxiter)
        candidate=serialize(chunks,schedule,angles);path=ROOT/f'search13_alternative_opening_refine_trial{trial}.json';path.write_text(json.dumps(candidate,indent=2)+'\n')
        check=evaluate(candidate);assert check['cnot_count']==13
        assert abs(record['loss']-objective(angles.ravel(),fixed,matrix_for_schedule(schedule),v,False))<1e-10
        row={'trial':trial,'sigma':sigma,'seed':seed,**record,'candidate':path.name,'cnot_count':check['cnot_count'],'off_diagonal_error':check['off_diagonal_error'],'valid_diagonalizer':check['valid_diagonalizer']};rows.append(row);print(json.dumps(row),flush=True)
    best=min(rows,key=lambda r:r['off_diagonal_error'])
    out.write_text(json.dumps({'source':SOURCE.name,'source_roundtrip_loss':loss,'schedule':schedule,'seed':args.seed,'maxiter':args.maxiter,'rows':rows,'best_candidate':best['candidate'],'scope':'Four local fits from one alternativeopening13CNOT topology, no global exclusion'},indent=2)+'\n')
if __name__=='__main__':main()
