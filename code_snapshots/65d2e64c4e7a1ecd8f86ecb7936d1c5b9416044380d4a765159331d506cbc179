"""Exact-eigenlabel least-squares refinement of the fresh chain13 near-hit.

Labels inferred once from the source and multiplicities checked. The
UV=DU equation leaves arbitrary basis freedom inside each labeled sector.
A failed local solve excludes neither other labels nor the topology.
"""
import argparse,json
from pathlib import Path
import numpy as np
import torch
from scipy.optimize import least_squares
from threadpoolctl import threadpool_limits
from check_circuit import evaluate,right_shift
from search13_adaptive_topology import matrix_for_schedule,serialize
from search13_joint_prefix_orbit_fit import unitary
from search13_prefix_perturb import decode_layers
from search13_fulltail_invariant_gauge_rank import circuit_matrix
ROOT=Path(__file__).resolve().parent
SOURCE=ROOT/'search13_fresh_chain_label_swaps_case4_polish.json'
torch.set_num_threads(1)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--max-nfev',type=int,default=10);p.add_argument('--seed',type=int,default=35800);args=p.parse_args()
    out=ROOT/'search13_fresh_labelswap_refine_result.json';assert not out.exists(),out
    source=json.loads(SOURCE.read_text());x0=decode_layers(source).ravel()
    schedule=tuple((g['control'],g['target']) for g in source['gates'] if g['gate']=='cx');assert len(schedule)==13
    fixed=[torch.eye(16,dtype=torch.complex128) for _ in range(14)];cxs=matrix_for_schedule(schedule);chunks=[[] for _ in range(14)]
    v_np=right_shift(4);u0=circuit_matrix(source['gates']);diagonal=np.diag(u0@v_np@u0.conj().T)
    roots=np.array([1,1j,-1,-1j]);labels=np.argmin(abs(diagonal[:,None]-roots[None,:]),axis=1)
    counts=np.bincount(labels,minlength=4);assert counts.tolist()==[6,3,4,3],counts
    distances=np.sort(abs(diagonal[:,None]-roots[None,:]),axis=1)
    assert np.min(distances[:,1]-distances[:,0])>1.,distances
    v=torch.as_tensor(v_np,dtype=torch.complex128);d=torch.as_tensor(roots[labels],dtype=torch.complex128)
    def residual_tensor(x):
        u=unitary(x.reshape(14,4,3),fixed,cxs)
        r=(u@v-d[:,None]*u)/4
        return torch.cat((r.real.ravel(),r.imag.ravel()))
    def residual(x):return residual_tensor(torch.as_tensor(x,dtype=torch.float64)).detach().numpy().copy()
    def jacobian(x):
        t=torch.tensor(x,dtype=torch.float64,requires_grad=True)
        return torch.autograd.functional.jacobian(residual_tensor,t,vectorize=True).detach().numpy().copy()
    # Independent identity check against the serialized source.
    expected=(u0@v_np-roots[labels,None]*u0)/4
    assert abs(np.dot(residual(x0),residual(x0))-np.sum(abs(expected)**2))<1e-12
    j0=jacobian(x0);assert j0.shape==(512,168)
    direction=np.random.default_rng(args.seed).normal(size=168);direction/=np.linalg.norm(direction);eps=1e-6
    fd=(residual(x0+eps*direction)-residual(x0-eps*direction))/(2*eps)
    jerror=float(np.max(abs(fd-j0@direction)));assert jerror<1e-7,jerror
    print(json.dumps({'labels':labels.tolist(),'counts':counts.tolist(),'initial_squared_residual':float(np.dot(residual(x0),residual(x0))),'jacobian_direction_error':jerror}),flush=True)
    rows=[]
    for trial,sigma in enumerate((0.,)):
        start=x0+np.random.default_rng(args.seed+trial).normal(0,sigma,x0.shape)
        with threadpool_limits(limits=1):
            result=least_squares(residual,start,jac=jacobian,method='trf',max_nfev=args.max_nfev,ftol=1e-13,xtol=1e-13,gtol=1e-11)
        candidate=serialize(chunks,schedule,result.x.reshape(14,4,3));path=ROOT/f'search13_fresh_labelswap_refine_trial{trial}.json';path.write_text(json.dumps(candidate,indent=2)+'\n')
        check=evaluate(candidate);assert check['cnot_count']==13
        u=circuit_matrix(candidate['gates']);independent=float(np.sum(abs(u@v_np-roots[labels,None]*u)**2)/16)
        assert abs(independent-2*result.cost)<1e-10
        row={'trial':trial,'sigma':sigma,'seed':args.seed+trial,'squared_residual':float(2*result.cost),'nfev':int(result.nfev),'njev':int(result.njev),'optimality':float(result.optimality),'optimizer_success':bool(result.success),'message':str(result.message),'candidate':path.name,'cnot_count':check['cnot_count'],'off_diagonal_error':check['off_diagonal_error'],'valid_diagonalizer':check['valid_diagonalizer']};rows.append(row);print(json.dumps(row),flush=True)
    best=min(rows,key=lambda r:r['off_diagonal_error'])
    out.write_text(json.dumps({'source':SOURCE.name,'schedule':schedule,'labels':labels.tolist(),'multiplicities':counts.tolist(),'minimum_nearest_label_margin':float(np.min(distances[:,1]-distances[:,0])),'jacobian_direction_error':jerror,'seed':args.seed,'max_nfev':args.max_nfev,'rows':rows,'best_candidate':best['candidate'],'scope':'Fixed-label Gauss-Newton refinement of numerical13CNOT nearhit; no exact certificate'},indent=2)+'\n')
if __name__=='__main__':main()
