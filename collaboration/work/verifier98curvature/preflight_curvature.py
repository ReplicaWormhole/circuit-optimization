"""Finite independent source/point/circuit preflight for run396 seed9261703.
No optimization, autodiff, Hessian, derivative, or finite differences.
"""
from __future__ import annotations
import hashlib, importlib.util, json
from pathlib import Path
import numpy as np, scipy
ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'collaboration/work/verifier98fresh/independent_audit.py'
spec=importlib.util.spec_from_file_location('independent_curve_base',BASE)
base=importlib.util.module_from_spec(spec); spec.loader.exec_module(base)
FIX=ROOT/'collaboration/work/verifier98label/objective_fixtures.py'
spec=importlib.util.spec_from_file_location('independent_curve_objective',FIX)
fixtures=importlib.util.module_from_spec(spec); spec.loader.exec_module(fixtures)

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def load_family(path):
    s=importlib.util.spec_from_file_location('curve_family',path)
    m=importlib.util.module_from_spec(s); s.loader.exec_module(m); return m

def main():
    run=ROOT/'experiments/runs/396'; cfg=json.loads((run/'config.json').read_text())
    result_path=run/'result.json'; result=json.loads(result_path.read_text())
    if sha(run/'batch.py') != cfg['code_sha256'] or result['config'] != cfg or result['run_id'] != 396:
        raise SystemExit('run396 frozen source/config/result binding mismatch')
    for rel,h in cfg['dependencies'].items():
        if sha(ROOT/rel)!=h: raise SystemExit(f'dependency changed: {rel}')
    audit=json.loads((run/'POSTFIT_AUDIT.json').read_text())
    if audit.get('result_sha256') != sha(result_path): raise SystemExit('prior independent audit is not bound to current run396 result')
    row=result['rows'][result['best_index']]
    if row['seed'] != 9261703 or cfg['seeds'][result['best_index']] != 9261703:
        raise SystemExit('global best is not reserved seed9261703')
    x=np.asarray(row['best'],dtype=float)
    if x.shape!=(98,) or not np.isfinite(x).all(): raise SystemExit('best point is not a finite 98-vector')
    # The frozen implementation draws both arrays sequentially from one generator.
    rng=np.random.default_rng(9261703)
    expected_initial=np.r_[rng.normal(0,3.0,90),rng.uniform(-np.pi/2,np.pi/2,8)]
    if not np.array_equal(np.asarray(row['initial'],dtype=float),expected_initial):
        raise SystemExit('seed9261703 initial vector does not match frozen RNG stream')
    native_path=run/row['native']; compiled_path=run/row['compiled']
    if sha(native_path)!=row['native_sha256'] or sha(compiled_path)!=row['compiled_sha256']:
        raise SystemExit('saved endpoint list hash mismatch')
    family=load_family(ROOT/cfg['family_module'])
    u=base.coordinate_matrix(x)
    native,ng=base.gate_list(native_path); compiled,cg=base.gate_list(compiled_path)
    gen_native=family.serialize(x); gen_compiled=family.helper.compile_native(gen_native)
    regenerated,rg=base.gate_list(gen_native); regenerated_compiled,rcg=base.gate_list(gen_compiled)
    errs={'coordinates_vs_saved_native':base.phase_error(u,native),
          'saved_native_vs_compiled':base.phase_error(native,compiled),
          'coordinates_vs_regenerated_native':base.phase_error(u,regenerated),
          'regenerated_native_vs_compiled':base.phase_error(regenerated,regenerated_compiled)}
    unitarity={'native':float(np.max(np.abs(native.conj().T@native-np.eye(16)))),
               'compiled':float(np.max(np.abs(compiled.conj().T@compiled-np.eye(16))))}
    costs={'saved_native':(sum(g['gate']=='cx' for g in ng),sum(g['gate']=='xx_yy' for g in ng)),
           'saved_compiled_cx':sum(g['gate']=='cx' for g in cg),
           'regenerated_native':(sum(g['gate']=='cx' for g in rg),sum(g['gate']=='xx_yy' for g in rg)),
           'regenerated_compiled_cx':sum(g['gate']=='cx' for g in rcg)}
    obj=fixtures.objective_components(u)
    if max(errs.values())>2e-12 or max(unitarity.values())>2e-12:
        raise SystemExit('independent coordinate/native/compiled or unitarity check failed')
    if costs['saved_native']!=(4,4) or costs['saved_compiled_cx']!=12:
        raise SystemExit(f'unexpected saved gate costs {costs}')
    if abs(obj['label_residual_direct']-row['loss'])>2e-12 or abs(obj['original_offdiagonal_loss']-row['original_loss'])>2e-12:
        raise SystemExit('independent objectives differ from run396 row')
    out={'run_id':396,'seed':9261703,'auditor_script_sha256':sha(Path(__file__)),'result_sha256':sha(result_path),
         'config_sha256':sha(run/'config.json'),'batch_sha256':sha(run/'batch.py'),
         'assignment_source_sha256':sha(ROOT/cfg['assignment_module']),
         'prior_postfit_audit_sha256':sha(run/'POSTFIT_AUDIT.json'),
         'vector_sha256_f64le':hashlib.sha256(np.asarray(x,dtype='<f8').tobytes()).hexdigest(),
         'native_sha256':sha(native_path),'compiled_sha256':sha(compiled_path),
         'initial_rng_exact':True,'coordinate_matrix_vs_saved_circuit_phase_errors':errs,
         'unitarity_max_errors':unitarity,'gate_counts':costs,
         'independent_objective':obj,'recorded_label_loss':row['loss'],
         'recorded_original_offdiagonal_loss':row['original_loss'],
         'checker_valid_diagonalizer':row['checker']['valid_diagonalizer'],
         'scope':'finite point/hash/circuit/objective preflight only; no AD/full-Hessian computations'}
    outpath=ROOT/'collaboration/work/verifier98curvature/PREFIT_POINT_AUDIT.json'
    outpath.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__': main()
