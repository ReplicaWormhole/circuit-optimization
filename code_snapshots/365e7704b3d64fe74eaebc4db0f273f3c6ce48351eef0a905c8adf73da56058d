"""Whole-circuit 13-CNOT fits without inherited prefix or local chunks.

Two unordered pattern families with two direction variants and Gaussian local starts
per order. This is bounded numerical exploration, not global exclusion.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import torch
from check_circuit import evaluate, right_shift
from search13_adaptive_topology import admissible, fit, matrix_for_schedule, objective, serialize
ROOT=Path(__file__).resolve().parent
torch.set_num_threads(1)


def topology(candidate):
    return tuple((g['control'],g['target']) for g in candidate['gates'] if g['gate']=='cx')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed',type=int,default=33800)
    parser.add_argument('--maxiter',type=int,default=300)
    parser.add_argument('--family',choices=['roundrobin-chain','ring-star'],default='roundrobin-chain')
    args=parser.parse_args()
    stem='search13_fresh_schedules' if args.family=='roundrobin-chain' else 'search13_fresh_ring_star'
    out=ROOT/f'{stem}_result.json'
    assert not out.exists(),out
    # This limited duplicate guard covers retained root-level search13 candidates.
    seen=set()
    for path in ROOT.glob('search13*.json'):
        data=json.loads(path.read_text())
        if isinstance(data,dict) and data.get('n')==4 and 'gates' in data:seen.add(topology(data))
    roundrobin=((0,1),(2,3),(0,2),(1,3),(0,3),(1,2))*2+((0,1),)
    chain=((0,1),(2,3),(1,2))*4+((0,1),)
    ring=((0,1),(1,2),(2,3),(0,3))*3+((0,1),)
    star=((0,1),(0,2),(0,3))*4+((0,1),)
    patterns=(roundrobin,chain) if args.family=='roundrobin-chain' else (ring,star)
    schedules=[]
    for pattern in patterns:
        for variant in range(2):
            # No projection to an exact target circuit or degenerate eigenspace gauge.
            schedule=tuple(edge[::-1] if (slot+variant)%2 else edge for slot,edge in enumerate(pattern))
            assert len(schedule)==13 and admissible(schedule)
            assert schedule not in seen
            assert schedule not in schedules
            schedules.append(schedule)
    fixed=[torch.eye(16,dtype=torch.complex128) for _ in range(14)]
    chunks=[[] for _ in range(14)]
    v=torch.as_tensor(right_shift(4),dtype=torch.complex128)
    rows=[]
    for case,schedule in enumerate(schedules):
        for start,sigma in enumerate((.8,2.0)):
            seed=args.seed+10*case+start
            initial=np.random.default_rng(seed).normal(0,sigma,(14,4,3))
            angles,record=fit(schedule,initial,fixed,v,args.maxiter)
            saved=serialize(chunks,schedule,angles)
            path=ROOT/f'{stem}_case{case}_start{start}.json'
            path.write_text(json.dumps(saved,indent=2)+'\n')
            check=evaluate(saved)
            assert topology(saved)==schedule and check['cnot_count']==13
            assert abs(record['loss']-objective(angles.ravel(),fixed,matrix_for_schedule(schedule),v,False))<1e-10
            row={'case':case,'start':start,'seed':seed,'normal_std':sigma,**record,'candidate':path.name,'cnot_count':check['cnot_count'],'off_diagonal_error':check['off_diagonal_error'],'valid_diagonalizer':check['valid_diagonalizer']}
            rows.append(row);print(json.dumps(row),flush=True)
    best=min(rows,key=lambda r:r['off_diagonal_error'])
    out.write_text(json.dumps({'family':args.family,'seed':args.seed,'maxiter':args.maxiter,'schedules':schedules,'excluded_prior_directed_count':len(seen),'initialization':'independent Gaussian SU2-axis coordinates, no inherited fixed gates; not Haar-distributed','rows':rows,'best_candidate':best['candidate'],'scope':'Four specified full 13-CNOT schedules, two independent starts each, bounded local fits; no topology no-go'},indent=2)+'\n')
if __name__=='__main__':main()
