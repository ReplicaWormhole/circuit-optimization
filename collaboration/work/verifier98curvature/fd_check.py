"""Reserved finite-difference validation; do not execute without run399 reservation.
Exactly 10 independent NumPy matrix evaluations; no autodiff, Hessian, fitting,
or Hungarian assignment at perturbed points.
"""
import hashlib, importlib.util, json, os
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
BASE=ROOT/'collaboration/work/verifier98fresh/independent_audit.py'
spec=importlib.util.spec_from_file_location('fd_independent_matrix',BASE)
base=importlib.util.module_from_spec(spec); spec.loader.exec_module(base)
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def label_value(u, labels, v):
    return float(np.linalg.norm(u@v-np.diag(labels)@u,'fro')**2/16)
def off_value(u,v):
    b=u@v@u.conj().T
    off=b-np.diag(np.diag(b))
    return float(np.linalg.norm(off,'fro')**2/16)
def main():
    manifest_path=HERE/'fd_manifest.json'
    frozen_path=HERE/'fd_config.json'
    if (HERE/'fd_result.json').exists(): raise RuntimeError('refusing to repeat reserved FD evaluation')
    frozen=json.loads(frozen_path.read_text())
    m=json.loads(manifest_path.read_text())
    if sha(Path(__file__))!=frozen['script_sha256'] or sha(manifest_path)!=frozen['manifest_sha256']:
        raise SystemExit('frozen FD script or manifest hash mismatch')
    if frozen.get('external_seconds')!=30 or frozen.get('threads')!=1 or frozen.get('max_matrix_evaluations')!=10:
        raise SystemExit('FD runtime budget differs from freeze')
    if m.get('schema')!='full98-independent-finite-difference-v1' or m.get('total_independent_matrix_evaluations')!=frozen['max_matrix_evaluations']:
        raise SystemExit('manifest schema/evaluation count mismatch')
    for rel,h in m['source_hashes'].items():
        if sha(ROOT/rel)!=h: raise SystemExit(f'frozen FD source changed: {rel}')
    parent=json.loads((ROOT/'experiments/runs/398/result.json').read_text())
    if sha(ROOT/'experiments/runs/398/result.json')!=m['source_result_sha256']:
        raise SystemExit('run398 source result hash mismatch')
    point=np.asarray(m['point'],dtype=float)
    if point.shape!=(98,) or hashlib.sha256(np.asarray(point,dtype='<f8').tobytes()).hexdigest()!=m['point_sha256_f64le']:
        raise SystemExit('frozen point malformed or hash mismatch')
    if not np.array_equal(point,np.asarray(parent['point98'],dtype=float)) or parent['seed']!=m['seed']:
        raise SystemExit('frozen point/seed differ from source diagnostic')
    records={r['name']:r for r in parent['records']}
    labels=np.asarray([complex(r,i) for r,i in m['fixed_labels_real_imag']],dtype=complex)
    source_labels=np.asarray([complex(r,i) for r,i in parent['labels_real_imag']],dtype=complex)
    if not np.array_equal(labels,source_labels):
        raise SystemExit('fixed FD labels differ from run398 frozen base assignment')
    for d in m['directions']:
        v=np.asarray(d['direction'],dtype=float); src=records[d['objective']]
        stored=np.asarray(src['eigenvectors_columns'],dtype=float)[:,0]
        if v.shape!=(98,) or not np.array_equal(v,stored): raise SystemExit('direction is not frozen source eigenvector column 0')
        if hashlib.sha256(np.asarray(v,dtype='<f8').tobytes()).hexdigest()!=d['direction_sha256_f64le']:
            raise SystemExit('direction hash mismatch')
        if abs(float(d['eigenvalue'])-float(src['least_eigenvalue']))>1e-15: raise SystemExit('direction eigenvalue binding mismatch')
    v4=base.right_shift()
    calls=[]
    def evaluate(point_,tag):
        if len(calls)>=10: raise RuntimeError('matrix-evaluation cap exceeded')
        u=base.coordinate_matrix(np.asarray(point_,dtype=float))
        unit=float(np.max(np.abs(u.conj().T@u-np.eye(16))))
        if unit>2e-12: raise SystemExit(f'nonunitary independent point {tag}: {unit}')
        calls.append({'index':len(calls)+1,'tag':tag,'matrix_unitarity_max_error':unit})
        return u
    report_centers={}
    with threadpool_limits(limits=1):
        # Separate center calls keep the frozen total at exactly ten matrix evaluations.
        uc_off=evaluate(point,'center_offdiagonal')
        report_centers['offdiagonal']=off_value(uc_off,v4)
        uc_label=evaluate(point,'center_fixed_label')
        report_centers['fixed_label_branch']=label_value(uc_label,labels,v4)
        samples={}
        for d in m['directions']:
            name=d['objective']; direction=np.asarray(d['direction'],dtype=float)
            for h in m['steps']:
                values={}
                for sign in (-1,1):
                    tag=f'{name}_h{h:g}_'+('minus' if sign<0 else 'plus')
                    u=evaluate(point+sign*float(h)*direction,tag)
                    values[sign]={name: off_value(u,v4) if name=='offdiagonal' else label_value(u,labels,v4)}
                objective=name
                f0=report_centers[objective]
                fm,fp=values[-1][objective],values[1][objective]
                fd1=(fp-fm)/(2*float(h))
                fd2=(fp-2*f0+fm)/(float(h)**2)
                g=np.asarray(records[objective]['gradient'],dtype=float)
                predicted=float(np.dot(g,direction))
                target=float(records[objective]['least_eigenvalue'])
                samples[f'{name}|h={h:g}']={
                  'direction_source':name,'objective':objective,'step':float(h),
                  'minus_value':fm,'center_value':f0,'plus_value':fp,
                  'central_first_difference':fd1,'autodiff_gradient_dot_direction':predicted,
                  'first_difference_abs_error':abs(fd1-predicted),
                  'central_second_difference':fd2,'source_least_eigenvalue':target,
                  'second_difference_abs_error_vs_source_eigenvalue':abs(fd2-target)}
    if len(calls)!=10: raise SystemExit(f'expected exactly10 independent matrix evaluations, got {len(calls)}')
    out={'schema':m['schema'],'run_id':int(HERE.name),'frozen_config_sha256':sha(frozen_path),
         'manifest_sha256':sha(manifest_path),'fd_script_sha256':sha(Path(__file__)),'source_run398_sha256':sha(ROOT/'experiments/runs/398/result.json'),
         'source_point_sha256_f64le':m['point_sha256_f64le'],'matrix_evaluation_count':len(calls),
         'matrix_evaluations':calls,'center_values':report_centers,'finite_difference_samples':samples,
         'runtime_budget':{'external_seconds':frozen['external_seconds'],'threads':frozen['threads'],'max_matrix_evaluations':frozen['max_matrix_evaluations']},
         'scope':'Independent finite differences along two frozen least-eigenvector columns for their matching objectives only; no AD/full Hessian/optimizer/assignment solve'}
    (HERE/'fd_result.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
