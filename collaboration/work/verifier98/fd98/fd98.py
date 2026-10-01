"""Independent centered finite differences along the frozen full98 least-Hessian direction."""
import os
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key]='1'
import hashlib, json, math, time
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from threadpoolctl import threadpool_limits
HERE=Path(__file__).resolve().parent
ROOT=Path(__file__).resolve().parents[3]
CFG=json.loads((HERE/'FROZEN_INPUT.json').read_text())

def cx(c,t):
    M=np.zeros((16,16),complex)
    for j in range(16): M[j^(1<<(3-t)) if ((j>>(3-c))&1) else j,j]=1
    return M

def oneq_rotation(v):
    x,y,z=map(float,v);length=math.sqrt(x*x+y*y+z*z)
    c=math.cos(length/2);s=0.5*np.sinc(length/(2*math.pi))
    return np.array([[c-1j*z*s,(-1j*x-y)*s],[(-1j*x+y)*s,c+1j*z*s]],complex)

def kron_local(vectors):
    M=np.array([[1]],complex)
    for v in vectors:M=np.kron(M,oneq_rotation(v))
    return M

def embed(P,q):
    M=np.array([[1]],complex)
    for k in range(4):M=np.kron(M,P if k==q else np.eye(2))
    return M

def fgate(a,b,c,t):
    X=np.array([[0,1],[1,0]],complex);Y=np.array([[0,-1j],[1j,0]],complex)
    return expm(1j*(a*embed(X,c)@embed(X,t)+b*embed(Y,c)@embed(Y,t)))

def matrix(z):
    local=np.asarray(z[:84]).reshape(7,4,3);ab=np.asarray(z[84:90]).reshape(2,3);f=np.asarray(z[90:]).reshape(4,2)
    pairs=((0,2),(1,3),(0,1),(2,3));U=np.eye(16,dtype=complex)
    U=kron_local(local[0])@U;U=cx(0,1)@U
    U=embed(oneq_rotation(ab[0]),1)@U;U=cx(2,1)@U
    U=embed(oneq_rotation(ab[1]),1)@U;U=cx(3,1)@U
    U=kron_local(local[1])@U;U=fgate(*f[0],*pairs[0])@U
    U=kron_local(local[2])@U;U=fgate(*f[1],*pairs[1])@U
    U=kron_local(local[3])@U;U=cx(1,2)@U
    U=kron_local(local[4])@U;U=fgate(*f[2],*pairs[2])@U
    U=kron_local(local[5])@U;U=fgate(*f[3],*pairs[3])@U
    U=kron_local(local[6])@U
    return U

def objective(z):
    V=np.zeros((16,16),complex)
    for j in range(16):V[((j&1)<<3)|(j>>1),j]=1
    U=matrix(z);B=U@V@U.conj().T;off=B-np.diag(np.diag(B))
    return float(np.vdot(off,off).real/16)

def main():
    script_hash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if script_hash!=CFG['code_sha256']: raise RuntimeError('frozen code hash mismatch')
    result_path=ROOT/'experiments/runs/389/result.json'
    if hashlib.sha256(result_path.read_bytes()).hexdigest()!=CFG['run389_result_sha256']:
        raise RuntimeError('run389 result changed after direction freeze')
    cfg=json.loads(result_path.read_text())
    x=np.asarray(CFG['point98'],float);v=np.asarray(CFG['least_eigenvector'],float)
    if x.shape!=(98,) or v.shape!=(98,) or abs(np.linalg.norm(v)-1)>2e-12:raise RuntimeError('bad frozen input shape/norm')
    if not np.array_equal(x,np.asarray(cfg['point98'],float)):raise RuntimeError('frozen point mismatch')
    if not np.array_equal(v,np.asarray(cfg['eigenvectors_columns'],float)[0]):raise RuntimeError('frozen least direction mismatch')
    started=time.monotonic()
    with threadpool_limits(limits=1):
        f0=objective(x)
        rows=[]
        for h in (1e-3,5e-4):
            fp=objective(x+h*v);fm=objective(x-h*v)
            curvature=(fp-2*f0+fm)/(h*h)
            rows.append({'h':h,'f_plus':fp,'f_minus':fm,'f0':f0,'centered_directional_curvature':curvature})
    out={'source_run':389,'source_result_sha256':CFG['run389_result_sha256'],'source_script_sha256':CFG['hessian_script_sha256'],
         'fd_script_sha256':script_hash,'point_sha256':CFG['point98_sha256'],'direction_sha256':CFG['least_eigenvector_sha256'],
         'objective_evaluations':5,'optimizer_starts':0,'optimizer_iterations':0,'threads':1,
         'elapsed_seconds':time.monotonic()-started,'base_loss':f0,'saved_AD_least_eigenvalue':float(cfg['least_eigenvalue']),
         'steps':rows,'interpretation':'Centered finite differences independently check one AD Hessian direction. This does not establish a full Hessian, minimum, diagonalizer, or topology exclusion.'}
    (HERE/'fd_result.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
