"""Append output-label CNOT gauge, delete earlier CNOT, fit native13.

Computational-basis permutations preserve exact diagonalization. The
source is approximate; deletion/fits remain numerical, no exact claim.
"""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
import torch
from check_circuit import evaluate,right_shift
from delete14_search import fixed_matrix
from search13_adaptive_topology import fit,matrix_for_schedule,objective,serialize
ROOT=Path(__file__).resolve().parent
SOURCE=ROOT/'search13_fresh_chain_refine_trial1.json'
torch.set_num_threads(1)

def split(gates,deleted):
    chunks=[[]];schedule=[];index=0
    for gate in gates:
        if gate['gate']=='cx':
            if index!=deleted:
                schedule.append((gate['control'],gate['target']));chunks.append([])
            index+=1
        else:chunks[-1].append(gate)
    assert index==14 and len(schedule)==13 and len(chunks)==14
    matrices=[]
    for chunk in chunks:
        m=torch.eye(16,dtype=torch.complex128)
        for gate in chunk:m=torch.as_tensor(fixed_matrix(gate),dtype=torch.complex128)@m
        matrices.append(m)
    return chunks,schedule,matrices

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--maxiter',type=int,default=250);args=p.parse_args()
    out=ROOT/'search13_fresh_chain_output_gauge_delete_result.json';assert not out.exists(),out
    source=json.loads(SOURCE.read_text());source_check=evaluate(source);assert source_check['cnot_count']==13
    v=torch.as_tensor(right_shift(4),dtype=torch.complex128);zero=np.zeros((14,4,3));rows=[]
    for gauge_id,edge in enumerate(((0,2),(1,3),(0,3))):
        gauged={'n':4,'gates':source['gates']+[{'gate':'cx','control':edge[0],'target':edge[1]}]}
        check=evaluate(gauged);assert check['cnot_count']==14
        assert abs(check['off_diagonal_error']-source_check['off_diagonal_error'])<1e-12
        for deleted in (0,4,8):
            chunks,schedule,fixed=split(gauged['gates'],deleted)
            # Zero free layers reproduce the deleted source gauge circuit.
            warm=serialize(chunks,schedule,zero);wu=evaluate(warm);assert wu['cnot_count']==13
            angles,record=fit(schedule,zero,fixed,v,args.maxiter)
            candidate=serialize(chunks,schedule,angles);path=ROOT/f'search13_fresh_chain_output_gauge_delete_g{gauge_id}_d{deleted}.json';path.write_text(json.dumps(candidate,indent=2)+'\n')
            q=evaluate(candidate);assert q['cnot_count']==13
            assert abs(record['loss']-objective(angles.ravel(),fixed,matrix_for_schedule(schedule),v,False))<1e-10
            row={'gauge':edge,'gauge_id':gauge_id,'deleted':deleted,'schedule':schedule,'angles':angles.tolist(),**record,'candidate':path.name,'cnot_count':q['cnot_count'],'off_diagonal_error':q['off_diagonal_error'],'valid_diagonalizer':q['valid_diagonalizer']};rows.append(row);print(json.dumps({k:row[k] for k in ['gauge','deleted','loss','off_diagonal_error','valid_diagonalizer']}),flush=True)
    best=min(rows,key=lambda r:r['off_diagonal_error'])
    out.write_text(json.dumps({'source':SOURCE.name,'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'maxiter':args.maxiter,'rows':rows,'best_candidate':best['candidate'],'scope':'Nine specified output-permutation/delete13CNOT fits with arbitrary free local layers, warmstart one per topology, no exact exclusion'},indent=2)+'\n')
if __name__=='__main__':main()
