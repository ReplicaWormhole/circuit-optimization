"""Independent NumPy objective and frozen directional finite differences."""
import json
import hashlib
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from verify import cycle, matrix


def axis_matrix(point, schedule):
    u=np.eye(16,dtype=complex)
    paulis=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],complex)
    for layer, angles in enumerate(np.asarray(point).reshape(14,4,3)):
        op=np.array([[1]],complex)
        for axis in angles:
            radius=np.linalg.norm(axis)
            local=np.cos(radius/2)*np.eye(2)-1j*(0.5*np.sinc(radius/(2*np.pi)))*np.einsum('i,ijk->jk',axis,paulis)
            op=np.kron(op,local)
        u=op@u
        if layer<13 and schedule[layer] is not None:
            c,t=schedule[layer]
            op=np.zeros((16,16),complex)
            for x in range(16):
                y=x^(1<<(3-t)) if ((x>>(3-c))&1) else x
                op[y,x]=1
            u=op@u
    return u


def loss(point,schedule):
    u=axis_matrix(point,schedule); b=u@cycle()@u.conj().T
    off=b-np.diag(np.diag(b))
    return float(np.sum(abs(off)**2)/16)


def main():
    cfgpath, auditpath, outpath=map(Path,sys.argv[1:4])
    cfg=json.loads(cfgpath.read_text()); audit=json.loads(auditpath.read_text())
    point=np.asarray(audit['base_axis']); direction=np.asarray(audit['least_eigenvector'])
    assert point.shape==(168,) and direction.shape==(168,)
    assert abs(np.linalg.norm(direction)-1)<1e-10
    schedule=cfg['optimizer_schedule']; base=loss(point,schedule)
    checks=[]
    for h in (1e-3,5e-4):
        minus=loss(point-h*direction,schedule);plus=loss(point+h*direction,schedule)
        checks.append({'h':h,'minus_loss':minus,'base_loss':base,'plus_loss':plus,
                       'central_directional_gradient':(plus-minus)/(2*h),
                       'central_directional_curvature':(plus-2*base+minus)/(h*h)})
    u=axis_matrix(point,schedule)
    hessian=np.asarray(audit['hessian']); hs=(hessian+hessian.T)/2
    gradient=np.asarray(audit['gradient']); eigenvalue=audit['eigenvalues'][0]
    output={'mapping':'C-order14layers,4wires,3axiscomponents x,y,z; exp(-i axis.sigma/2)',
            'point_sha256':hashlib.sha256(auditpath.read_bytes()).hexdigest(),
            'eigenvector_norm':float(np.linalg.norm(direction)),
            'unitarity_error':float(np.max(abs(u.conj().T@u-np.eye(16)))),
            'independent_loss':base,'saved_eigenvalue':eigenvalue,
            'hessian_symmetry_max_error':float(np.max(abs(hessian-hessian.T))),
            'eigenpair_residual':float(np.linalg.norm(hs@direction-eigenvalue*direction)),
            'gradient_norm':float(np.linalg.norm(gradient)),
            'saved_gradient_directional':float(gradient@direction),
            'checks':checks}
    Path(outpath).write_text(json.dumps(output,indent=2)+'\n');print(json.dumps(output,indent=2))


if __name__=='__main__':main()
