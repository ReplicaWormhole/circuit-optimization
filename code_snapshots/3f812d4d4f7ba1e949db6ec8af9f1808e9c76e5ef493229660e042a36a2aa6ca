"""Change the fit landscape while preserving all cycle eigenspaces.

H_t=Re(V)+t Im(V) has eigenvalues 1,t,-1,-t. For real t not in
{0,1,-1}, its exact diagonalizers are exactly those of V. Numerical failure
of a bounded fit excludes nothing globally.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import torch
from check_circuit import evaluate, right_shift
from search13_adaptive_topology import fit, matrix_for_schedule, objective, serialize
from search13_joint_prefix_orbit_fit import CASES, unitary
from search13_multi_pair_block_rewire import prepare
from search13_prefix_perturb import decode_layers
from search13_joint_prefix_orbit_fit import fit as generic_fit
ROOT=Path(__file__).resolve().parent
SOURCE=ROOT/'search13_joint_prefix_slot8_refine_trial1.json'
torch.set_num_threads(1)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--maxiter',type=int,default=200)
    parser.add_argument('--direct-maxiter',type=int,default=300)
    args=parser.parse_args()
    out=ROOT/'search13_hermitian_reweight_result.json'
    assert not out.exists(),out
    source=json.loads(SOURCE.read_text());x0=decode_layers(source)
    deleted,rewire,blocks=CASES[1]
    chunks,cxs,fixed,_=prepare(deleted,blocks,prefix_rewire=rewire)
    schedule=[(g['control'],g['target']) for g in cxs];schedule[8]=(1,2)
    assert tuple(schedule)==tuple((g['control'],g['target']) for g in source['gates'] if g['gate']=='cx')
    cxm=matrix_for_schedule(schedule)
    v=torch.as_tensor(right_shift(4),dtype=torch.complex128)
    baseline=objective(x0.ravel(),fixed,cxm,v,False)
    assert abs(baseline-0.06019279505817157)<1e-9,baseline
    real=(v+v.conj().T)/2;imag=(v-v.conj().T)/(2j)
    rows=[]
    for trial,t in enumerate((.25,4.,-.25,-4.)):
        assert t not in (0.,1.,-1.)
        h=real+t*imag
        assert torch.max(torch.abs(h-h.conj().T))<1e-14
        spectrum=np.linalg.eigvalsh(h.numpy())
        assert np.allclose(np.unique(np.round(spectrum,10)), sorted((1.,t,-1.,-t)))
        normal=float(h.abs().square().sum())
        def target_loss(flat):
            x=torch.tensor(np.asarray(flat).reshape(14,4,3),dtype=torch.float64,requires_grad=True)
            u=unitary(x,fixed,cxm);d=u@h@u.conj().T
            off=d-torch.diag(torch.diagonal(d))
            loss=off.abs().square().sum().real/normal
            grad=torch.autograd.grad(loss,x)[0].detach().numpy().ravel().copy()
            return float(loss.detach()),grad
        # Check gradient in a deterministic direction before optimizing.
        rng=np.random.default_rng(33500+trial);direction=rng.normal(size=x0.size);direction/=np.linalg.norm(direction)
        value,grad=target_loss(x0.ravel());eps=1e-6
        fd=(target_loss(x0.ravel()+eps*direction)[0]-target_loss(x0.ravel()-eps*direction)[0])/(2*eps)
        assert abs(fd-np.dot(grad,direction))<1e-7,(fd,np.dot(grad,direction))
        angles,weighted_record=generic_fit(target_loss,x0,args.maxiter)
        weighted_cycle_loss=objective(angles.ravel(),fixed,cxm,v,False)
        weighted=serialize(chunks,schedule,angles)
        path=ROOT/f'search13_hermitian_reweight_trial{trial}_weighted.json'
        path.write_text(json.dumps(weighted,indent=2)+'\n')
        weighted_check=evaluate(weighted)
        final,record=fit(schedule,angles,fixed,v,args.direct_maxiter)
        path2=ROOT/f'search13_hermitian_reweight_trial{trial}_direct.json'
        candidate=serialize(chunks,schedule,final);path2.write_text(json.dumps(candidate,indent=2)+'\n')
        check=evaluate(candidate)
        assert check['cnot_count']==weighted_check['cnot_count']==13
        assert abs(record['loss']-objective(final.ravel(),fixed,cxm,v,False))<1e-10
        row={'trial':trial,'t':t,'hermitian_eigenvalues':spectrum.tolist(),'initial_weighted_loss':value,'gradient_direction_check_error':float(abs(fd-np.dot(grad,direction))),'weighted':weighted_record,'weighted_cycle_loss':weighted_cycle_loss,'weighted_candidate':path.name,'weighted_off_diagonal_error':weighted_check['off_diagonal_error'],'weighted_valid':weighted_check['valid_diagonalizer'],'direct':record,'candidate':path2.name,'cnot_count':check['cnot_count'],'off_diagonal_error':check['off_diagonal_error'],'valid_diagonalizer':check['valid_diagonalizer']}
        rows.append(row);print(json.dumps({'trial':trial,'t':t,'weighted_cycle_loss':weighted_cycle_loss,'direct_loss':record['loss'],'maxoff':check['off_diagonal_error'],'valid':check['valid_diagonalizer']}),flush=True)
    best=min(rows,key=lambda r:r['off_diagonal_error'])
    out.write_text(json.dumps({'source':SOURCE.name,'baseline_loss':baseline,'maxiter':args.maxiter,'direct_maxiter':args.direct_maxiter,'rows':rows,'best_candidate':best['candidate'],'scope':'Four Hermitian spectral reweightings and direct-polish fits at one 13-CNOT topology; no exact exclusion'},indent=2)+'\n')
if __name__=='__main__':main()
