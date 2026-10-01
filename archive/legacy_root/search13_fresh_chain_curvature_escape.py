"""Numerical curvature and directed escapes of the fresh chain local basin.

Full 168-parameter Hessian of label-free cycle diagonalization loss. No
positivity result here is a certified lower bound or global exclusion.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import torch
from check_circuit import evaluate, right_shift
from search13_adaptive_topology import fit, matrix_for_schedule, objective, serialize
from search13_joint_prefix_orbit_fit import unitary
from search13_prefix_perturb import decode_layers

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / 'search13_fresh_chain_refine_trial1.json'
torch.set_num_threads(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--maxiter', type=int, default=250)
    args = parser.parse_args()
    out = ROOT / 'search13_fresh_chain_curvature_escape_result.json'
    assert not out.exists(), out
    candidate = json.loads(SOURCE.read_text())
    x0 = decode_layers(candidate).ravel()
    chunks = [[] for _ in range(14)]
    fixed = [torch.eye(16, dtype=torch.complex128) for _ in range(14)]
    schedule = tuple((g['control'], g['target']) for g in candidate['gates'] if g['gate']=='cx')
    assert schedule == tuple((g['control'], g['target']) for g in candidate['gates'] if g['gate']=='cx')
    cxm = matrix_for_schedule(schedule)
    v = torch.as_tensor(right_shift(4), dtype=torch.complex128)
    def scalar(x):
        u = unitary(x.reshape(14,4,3), fixed, cxm)
        d = u @ v @ u.conj().T
        off = d - torch.diag(torch.diagonal(d))
        return off.abs().square().sum().real / 16
    x = torch.tensor(x0, dtype=torch.float64, requires_grad=True)
    loss = scalar(x)
    grad = torch.autograd.grad(loss, x)[0].detach().numpy()
    baseline = objective(x0, fixed, cxm, v, False)
    assert abs(float(loss.detach())-baseline)<1e-12
    assert abs(baseline-0.0012158342275705642)<1e-8
    h = torch.autograd.functional.hessian(scalar, x).detach().numpy()
    asym = float(np.max(np.abs(h-h.T)))
    assert asym<1e-10,asym
    vals, vecs = np.linalg.eigh((h+h.T)/2)
    checks=[]
    for index in (0,167):
        direction=vecs[:,index];eps=1e-4
        fd=(objective(x0+eps*direction,fixed,cxm,v,False)+objective(x0-eps*direction,fixed,cxm,v,False)-2*baseline)/eps**2
        assert abs(fd-vals[index])<2e-5,(index,fd,vals[index])
        checks.append({'index':index,'eigenvalue':float(vals[index]),'finite_difference':float(fd)})
    print(json.dumps({'baseline':baseline,'gradient_norm':float(np.linalg.norm(grad)), 'min_eigenvalue':float(vals[0]),'max_eigenvalue':float(vals[-1]),'negative_below_1e8':int(np.sum(vals < -1e-8))}),flush=True)
    # Probe three lowest-curvature independent directions at both signs.
    # Large steps can cross a basin even when the Hessian has no negative mode.
    negative=np.flatnonzero(vals < -1e-6)
    indices=negative[:3].tolist() if len(negative) else np.flatnonzero(vals > 1e-4)[:3].tolist()
    assert len(indices)==3, indices
    rows=[]
    for index in indices:
        for sign in (-1,1):
            step=sign*np.pi
            initial=x0+step*vecs[:,index]
            fitted,record=fit(schedule,initial.reshape(14,4,3),fixed,v,args.maxiter)
            saved=serialize(chunks,schedule,fitted)
            path=ROOT/f'search13_fresh_chain_curvature_escape_mode{index}_sign{sign}.json'
            path.write_text(json.dumps(saved,indent=2)+'\n')
            check=evaluate(saved)
            assert check['cnot_count']==13
            assert abs(record['loss']-objective(fitted.ravel(), fixed, cxm, v, False))<1e-10
            row={'mode':index,'step':step,**record,'candidate':path.name,'cnot_count':check['cnot_count'],'off_diagonal_error':check['off_diagonal_error'],'valid_diagonalizer':check['valid_diagonalizer']}
            rows.append(row);print(json.dumps(row),flush=True)
    best=min(rows,key=lambda r:r['off_diagonal_error'])
    out.write_text(json.dumps({'source':SOURCE.name,'baseline':baseline,'maxiter':args.maxiter,'gradient_norm':float(np.linalg.norm(grad)),'hessian_asymmetry':asym,'hessian_eigenvalues':vals.tolist(),'finite_difference_checks':checks,'rows':rows,'best_candidate':best['candidate'],'scope':'Numerical Hessian at one coordinate point and six bounded fits at one topology; no global exclusion'},indent=2)+'\n')

if __name__=='__main__':main()
