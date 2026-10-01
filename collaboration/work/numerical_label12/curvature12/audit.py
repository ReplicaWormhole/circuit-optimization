"""One reserved projected24-dimensional automatic-differentiation Hessian."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import json,hashlib,time,importlib.util
from pathlib import Path
import numpy as np
import torch
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[4];HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('relaxed_source',ROOT/'collaboration/work/numerical_label12/relaxed12/search.py');family=importlib.util.module_from_spec(spec);spec.loader.exec_module(family)
from check_circuit import right_shift
def main():
    cfg=json.loads((HERE/'config.json').read_text());assert not (HERE/'audit_result.json').exists()
    for name,digest in cfg['dependencies'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest
    assert hashlib.sha256((ROOT/cfg['source']).read_bytes()).hexdigest()==cfg['source_sha256']
    point=torch.tensor(cfg['point116'],dtype=torch.float64);v=torch.tensor(right_shift(4),dtype=torch.complex128);torch.set_num_threads(1);begun=time.monotonic()
    def loss(ab):
        z=torch.cat((point[:12],ab,point[36:]));u=family.matrix(z);b=u@v@u.conj().T;off=b-torch.diag(torch.diagonal(b));return (off.abs()**2).sum().real/16
    with threadpool_limits(limits=1):
        ab=point[12:36].clone().requires_grad_(True);value=loss(ab);gradient=torch.autograd.grad(value,ab)[0].detach().numpy();h=torch.autograd.functional.hessian(loss,ab,vectorize=True).detach().numpy();sym=(h+h.T)/2;values,vectors=np.linalg.eigh(sym);direction=vectors[:,0]
        if direction[np.argmax(abs(direction))]<0:direction=-direction
        result={'point116':cfg['point116'],'projected_slice':[12,36],'objective':float(value.detach()),'gradient24':gradient.tolist(),'gradient_norm':float(np.linalg.norm(gradient)),'hessian24':h.tolist(),'symmetry_max_error':float(np.max(abs(h-h.T))),'eigenvalues':values.tolist(),'least_eigenvector24':direction.tolist(),'least_eigenvalue':float(values[0]),'eigenpair_residual':float(np.linalg.norm(sym@direction-values[0]*direction)),'vector_norm':float(np.linalg.norm(direction)),'elapsed_seconds':time.monotonic()-begun,'threshold':-1e-6,'below_threshold':bool(values[0]<-1e-6),'scope':'FloatingAD projectedHessian, no exactstationarity/minimum claim; conditionalfitsrequireverifierFDreview/newreservation'}
    (HERE/'audit_result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:w for k,w in result.items() if k not in ['point116','hessian24','least_eigenvector24','gradient24']}))
if __name__=='__main__':main()
