"""Deepen fresh thirteen-CNOT chain lead, without fixed inherited chunks."""
import argparse
import json
from pathlib import Path
import numpy as np
import torch
from check_circuit import evaluate,right_shift
from search13_adaptive_topology import fit,matrix_for_schedule,objective,serialize
from search13_prefix_perturb import decode_layers
ROOT=Path(__file__).resolve().parent
SOURCE=ROOT/'search13_fresh_schedules_case2_start1.json'
torch.set_num_threads(1)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--maxiter',type=int,default=1200);p.add_argument('--seed',type=int,default=34000);args=p.parse_args()
    out=ROOT/'search13_fresh_chain_refine_result.json';assert not out.exists(),out
    source=json.loads(SOURCE.read_text());initial=decode_layers(source)
    schedule=tuple((g['control'],g['target']) for g in source['gates'] if g['gate']=='cx')
    assert len(schedule)==13
    fixed=[torch.eye(16,dtype=torch.complex128) for _ in range(14)];chunks=[[] for _ in range(14)];v=torch.as_tensor(right_shift(4),dtype=torch.complex128)
    loss=objective(initial.ravel(),fixed,matrix_for_schedule(schedule),v,False)
    assert abs(loss-.0012158348375901386)<1e-10,loss
    rows=[]
    for trial,sigma in enumerate((0.,.02,.1)):
        start=initial+np.random.default_rng(args.seed+trial).normal(0,sigma,initial.shape)
        angles,record=fit(schedule,start,fixed,v,args.maxiter)
        saved=serialize(chunks,schedule,angles);path=ROOT/f'search13_fresh_chain_refine_trial{trial}.json';path.write_text(json.dumps(saved,indent=2)+'\n')
        check=evaluate(saved);assert check['cnot_count']==13
        row={'trial':trial,'sigma':sigma,'seed':args.seed+trial,**record,'candidate':path.name,'cnot_count':check['cnot_count'],'off_diagonal_error':check['off_diagonal_error'],'valid_diagonalizer':check['valid_diagonalizer']};rows.append(row);print(json.dumps(row),flush=True)
    best=min(rows,key=lambda r:r['off_diagonal_error'])
    out.write_text(json.dumps({'source':SOURCE.name,'source_roundtrip_loss':loss,'schedule':schedule,'seed':args.seed,'maxiter':args.maxiter,'rows':rows,'best_candidate':best['candidate'],'scope':'Three bounded local continuations on one fresh chain13 topology; numerical results require exact confirmation if valid'},indent=2)+'\n')
if __name__=='__main__':main()
