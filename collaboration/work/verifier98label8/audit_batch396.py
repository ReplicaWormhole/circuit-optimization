"""Independent finite preflight and endpoint audit for reserved run396.

No optimizer, derivative, finite difference, Hessian, or parameter search.
"""
from __future__ import annotations
import hashlib
import importlib.util
import json
import sys
from pathlib import Path
import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
BASE_PATH = ROOT / 'collaboration/work/verifier98fresh/independent_audit.py'
BASE_SPEC = importlib.util.spec_from_file_location('verifier98_base396', BASE_PATH)
base = importlib.util.module_from_spec(BASE_SPEC)
BASE_SPEC.loader.exec_module(base)
FIX_PATH = ROOT / 'collaboration/work/verifier98label/objective_fixtures.py'
FIX_SPEC = importlib.util.spec_from_file_location('label_fixtures396', FIX_PATH)
fixtures = importlib.util.module_from_spec(FIX_SPEC)
FIX_SPEC.loader.exec_module(fixtures)
LABELS = fixtures.LABELS


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def family_for(cfg):
    spec = importlib.util.spec_from_file_location('family396_audit', ROOT / cfg['family_module'])
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def check_frozen(run: Path, cfg: dict):
    if digest(run / 'batch.py') != cfg['code_sha256']:
        raise SystemExit('batch script differs from frozen code hash')
    if scipy.__version__ != cfg['scipy_version']:
        raise SystemExit('SciPy version differs from frozen config')
    for rel, expected in cfg['dependencies'].items():
        if digest(ROOT / rel) != expected:
            raise SystemExit(f'frozen dependency hash mismatch: {rel}')
    if len(cfg['seeds']) != 8 or cfg['seeds'] != list(range(9261700, 9261708)):
        raise SystemExit('seed list differs from reservation')
    if cfg['local_sigmas'] != [0.3, 1.0, 2.0, 3.0]:
        raise SystemExit('local sigma cycle differs from reservation')
    expected_caps = {'maxiter':1000,'maxfun':20000,'maxls':40,
                     'internal_seconds':110,'external_seconds':120,'threads':1,
                     'nearhit_loss':1e-18}
    for key, val in expected_caps.items():
        if cfg.get(key) != val:
            raise SystemExit(f'budget differs for {key}')
    if cfg.get('assignment_module') != 'experiments/runs/395/fit.py':
        raise SystemExit('frozen assignment module path changed')
    if cfg['dependencies'].get(cfg['assignment_module']) != digest(ROOT / cfg['assignment_module']):
        raise SystemExit('run395 assignment function hash mismatch')
    if cfg.get('parameter_order','').count('98') == 0 or 'all coordinates free' not in cfg.get('parameter_order',''):
        raise SystemExit('98-free parameter statement missing')
    # Ensure the batch uses the frozen evaluator and no coordinate mask.
    source = (run / 'batch.py').read_text()
    required = ('labels,columns,cost=label_module.assignment(u.detach().numpy(),v_np)',
                'd=torch.tensor(labels,dtype=torch.complex128)',
                'residual=u@torch.tensor(v_np,dtype=torch.complex128)-d[:,None]*u',
                'torch.autograd.grad(value,x)[0]',
                "'maxiter':cfg['maxiter']", "'maxfun':cfg['maxfun']")
    for fragment in required:
        if fragment not in source:
            raise SystemExit(f'batch implementation does not contain expected fragment: {fragment}')
    if 'bounds=' in source or 'jacobian mask' in source.lower():
        raise SystemExit('unexpected coordinate restriction in fit source')
    if (run / 'result.json').exists():
        raise SystemExit('result already exists; prefit mode refuses to proceed')
    return expected_caps


def rng_starts(cfg):
    records=[]
    for i, seed in enumerate(cfg['seeds']):
        rng=np.random.default_rng(seed)
        sigma=cfg['local_sigmas'][i % len(cfg['local_sigmas'])]
        x=np.r_[rng.normal(0,sigma,90),rng.uniform(-np.pi/2,np.pi/2,8)]
        if x.shape != (98,) or not np.isfinite(x).all():
            raise SystemExit(f'seed {seed} did not produce a finite 98-vector')
        records.append({'seed':seed,'sigma':sigma,'vector':x,
                        'vector_sha256_f64le':hashlib.sha256(np.asarray(x,dtype='<f8').tobytes()).hexdigest(),
                        'normal90_sha256_f64le':hashlib.sha256(np.asarray(x[:90],dtype='<f8').tobytes()).hexdigest(),
                        'uniform8_sha256_f64le':hashlib.sha256(np.asarray(x[90:],dtype='<f8').tobytes()).hexdigest()})
    return records


def preflight(run: Path, cfg: dict):
    caps=check_frozen(run,cfg)
    starts=rng_starts(cfg)
    family=family_for(cfg)
    reps=[]
    for rec in starts[:4]:
        x=rec['vector']
        native=family.serialize(x)
        compiled=family.helper.compile_native(native)
        un, ng=base.gate_list(native)
        uc, cg=base.gate_list(compiled)
        um=base.coordinate_matrix(x)
        phase_coord=base.phase_error(um,un)
        phase_comp=base.phase_error(un,uc)
        unitarity_n=float(np.max(np.abs(un.conj().T@un-np.eye(16))))
        unitarity_c=float(np.max(np.abs(uc.conj().T@uc-np.eye(16))))
        counts=(sum(g['gate']=='cx' for g in ng),sum(g['gate']=='xx_yy' for g in ng),sum(g['gate']=='cx' for g in cg))
        objective=fixtures.objective_components(um)
        if counts != (4,4,12) or max(phase_coord,phase_comp,unitarity_n,unitarity_c)>2e-12:
            raise SystemExit(f"representative seed {rec['seed']} matrix/compiler check failed")
        if max(objective['identity_abs_difference'],objective['cost_abs_difference'])>2e-12:
            raise SystemExit(f"representative seed {rec['seed']} raw-label objective identity failed")
        multiplicities={}
        for val in (1+0j,-1+0j,1j,-1j):
            multiplicities[str(val)]=sum(abs(z-val)<1e-12 for z in fixtures.assign_labels(um)[0])
        if [multiplicities[str(z)] for z in (1+0j,-1+0j,1j,-1j)] != [6,4,3,3]:
            raise SystemExit('assigned spectrum multiplicities differ')
        reps.append({'seed':rec['seed'],'sigma':rec['sigma'],
                     'coordinate_vs_native_error_up_to_phase':phase_coord,
                     'native_vs_compiled_error_up_to_phase':phase_comp,
                     'native_unitarity_max_error':unitarity_n,
                     'compiled_unitarity_max_error':unitarity_c,
                     'native_cx_count':counts[0],'native_xx_yy_count':counts[1],
                     'compiled_cx_count':counts[2],
                     'label_loss':objective['label_residual_direct'],
                     'assigned_cost_loss':objective['raw_assignment_cost_normalized'],
                     'original_offdiagonal_loss':objective['original_offdiagonal_loss'],
                     'original_offdiagonal_max':objective['original_offdiagonal_max'],
                     'assignment_columns':objective['assignment_columns'],
                     'assignment_cost_identity_error':objective['cost_abs_difference'],
                     'direct_conjugated_identity_error':objective['identity_abs_difference']})
    summary={'run_id':396,'board_hypothesis_id':89,
             'batch_sha256':digest(run/'batch.py'),'config_sha256':digest(run/'config.json'),
             'assignment_source_sha256':digest(ROOT/cfg['assignment_module']),
             'scipy_version':scipy.__version__,'caps':caps,
             'seeds':[{'seed':r['seed'],'sigma':r['sigma'],'vector_sha256_f64le':r['vector_sha256_f64le'],
                       'normal90_sha256_f64le':r['normal90_sha256_f64le'],
                       'uniform8_sha256_f64le':r['uniform8_sha256_f64le']} for r in starts],
             'representative_checks':reps,
             'scope':'RNG provenance all8; matrix/compiler/objective fixtures at four sigma representatives; no optimizer or derivatives'}
    return summary


def postfit(run: Path, cfg: dict):
    if digest(run/'batch.py') != cfg['code_sha256']:
        raise SystemExit('batch source hash changed after reservation')
    for rel, expected in cfg['dependencies'].items():
        if digest(ROOT/rel) != expected:
            raise SystemExit(f'postfit dependency hash changed: {rel}')
    result_path=run/'result.json'
    result=json.loads(result_path.read_text())
    if result.get('run_id') != 396 or result.get('config') != cfg:
        raise SystemExit('result run/config mismatch')
    expected_starts=rng_starts(cfg)
    rows=result.get('rows',[])
    if len(rows)>8 or len(rows)==0:
        raise SystemExit('batch row count outside reserved 1..8')
    family=family_for(cfg)
    audited=[]
    used=set()
    for i,row in enumerate(rows):
        start=expected_starts[i]
        seed=start['seed']
        x0=np.asarray(row['initial'],dtype=float)
        x=np.asarray(row['best'],dtype=float)
        if row.get('seed') != seed or row.get('sigma') != start['sigma'] or not np.array_equal(x0,start['vector']):
            raise SystemExit(f'seed/sigma/initial RNG mismatch at row {i}')
        if x.shape!=(98,) or not np.isfinite(x).all():
            raise SystemExit(f'best vector for seed {seed} is invalid')
        npath=run/row['native']; cpath=run/row['compiled']
        used.update((npath.name,cpath.name))
        if digest(npath)!=row['native_sha256'] or digest(cpath)!=row['compiled_sha256']:
            raise SystemExit(f'gate-list hash mismatch for seed {seed}')
        u=base.coordinate_matrix(x)
        un,ng=base.gate_list(npath); uc,cg=base.gate_list(cpath)
        c1=base.phase_error(u,un); c2=base.phase_error(un,uc)
        unitarity_n=float(np.max(np.abs(un.conj().T@un-np.eye(16))))
        unitarity_c=float(np.max(np.abs(uc.conj().T@uc-np.eye(16))))
        counts=(sum(g['gate']=='cx' for g in ng),sum(g['gate']=='xx_yy' for g in ng),sum(g['gate']=='cx' for g in cg))
        if counts!=(4,4,12) or max(c1,c2,unitarity_n,unitarity_c)>2e-12:
            raise SystemExit(f'matrix/count/unity check failed for seed {seed}')
        obj=fixtures.objective_components(u)
        if abs(float(row['loss'])-obj['label_residual_direct'])>2e-12:
            raise SystemExit(f'label objective differs for seed {seed}')
        if abs(float(row['original_loss'])-obj['original_offdiagonal_loss'])>2e-12:
            raise SystemExit(f'original objective differs for seed {seed}')
        if row.get('history') is None or len(row['history'])!=row['evaluations']:
            raise SystemExit(f'history/evaluation count mismatch for seed {seed}')
        hist_changes=0; prev=None
        for j,h in enumerate(row['history'],1):
            cols=h.get('columns')
            if h.get('evaluation') != j or len(cols or [])!=16:
                raise SystemExit(f'malformed assignment history for seed {seed}')
            counts_hist=np.bincount(np.asarray(cols,dtype=int),minlength=16)
            label_mult=[sum(int(col in np.where(LABELS==z)[0]) for col in cols) for z in (1,-1,1j,-1j)]
            if label_mult != [6,4,3,3]:
                raise SystemExit(f'history assignment multiplicities differ for seed {seed}')
            if prev is not None and cols != prev:
                hist_changes+=1
            prev=cols
        hist_losses=[float(h['loss']) for h in row['history']]
        best_hist_index=int(np.argmin(hist_losses)) if hist_losses else None
        if best_hist_index is None or abs(hist_losses[best_hist_index]-float(row['loss']))>2e-12:
            raise SystemExit(f'best loss is absent from evaluation history for seed {seed}')
        hist_labels=LABELS[np.asarray(row['history'][best_hist_index]['columns'],dtype=int)]
        independent_labels=np.asarray([complex(re,im) for re,im in obj['labels_by_row_real_imag']],dtype=complex)
        if not np.array_equal(hist_labels,independent_labels):
            raise SystemExit(f'best evaluation label vector differs under independent SciPy for seed {seed}')
        checker=row.get('checker',{})
        audited.append({'seed':seed,'sigma':start['sigma'],'evaluations':row['evaluations'],
                        'iterations':row['nit'],'termination':row['termination'],
                        'label_loss_independent':obj['label_residual_direct'],
                        'assigned_cost_loss':obj['raw_assignment_cost_normalized'],
                        'original_offdiagonal_loss_independent':obj['original_offdiagonal_loss'],
                        'original_offdiagonal_max':obj['original_offdiagonal_max'],
                        'assignment_columns':obj['assignment_columns'],
                        'history_assignment_changes':hist_changes,
                        'coordinate_vs_native_error_up_to_phase':c1,
                        'native_vs_compiled_error_up_to_phase':c2,
                        'native_unitarity_max_error':unitarity_n,'compiled_unitarity_max_error':unitarity_c,
                        'native_cx_count':counts[0],'native_xx_yy_count':counts[1],'compiled_cx_count':counts[2],
                        'native_sha256':digest(npath),'compiled_sha256':digest(cpath),
                        'checker_valid':checker.get('valid_diagonalizer'),
                        'checker_offdiagonal_error':checker.get('off_diagonal_error')})
    if result.get('best_index') != min(range(len(rows)),key=lambda i: rows[i]['loss']):
        raise SystemExit('global best index is not minimum saved label loss')
    best=rows[result['best_index']]
    for name in ('best_native.json','best_compiled.json'):
        path=run/name
        if not path.exists():
            raise SystemExit(f'missing global {name}')
    if (run/'best_native.json').read_bytes() != (run/best['native']).read_bytes() or (run/'best_compiled.json').read_bytes() != (run/best['compiled']).read_bytes():
        raise SystemExit('global best gate lists do not match best row')
    used.update(('best_native.json','best_compiled.json'))
    gate_files={p.name for p in run.glob('endpoint_*_native.json')}|{p.name for p in run.glob('endpoint_*_compiled.json')}|{p.name for p in run.glob('best_*json')}
    if gate_files!=used:
        raise SystemExit(f'extra/missing saved gate lists: {gate_files ^ used}')
    return {'run_id':396,'board_hypothesis_id':89,'batch_sha256':digest(run/'batch.py'),
            'config_sha256':digest(run/'config.json'),'result_sha256':digest(result_path),
            'assignment_source_sha256':digest(ROOT/cfg['assignment_module']),
            'rows_audited':len(rows),'elapsed_seconds':result.get('elapsed_seconds'),
            'stop_reason':result.get('stop_reason'),'best_index':result['best_index'],
            'best_seed':best['seed'],'best_label_loss':audited[result['best_index']]['label_loss_independent'],
            'endpoints':audited,
            'scope':'independent RNG, complete gate-list, dynamic-label and original residual audit; no exact certificate'}


def main():
    if len(sys.argv)!=3 or sys.argv[2] not in ('prefit','postfit'):
        raise SystemExit('usage: audit_batch396.py RUN_DIR prefit|postfit')
    run=Path(sys.argv[1]); run=run if run.is_absolute() else ROOT/run
    cfg=json.loads((run/'config.json').read_text())
    report=preflight(run,cfg) if sys.argv[2]=='prefit' else postfit(run,cfg)
    out=HERE/('PREFIT_APPROVAL.json' if sys.argv[2]=='prefit' else 'POSTFIT_AUDIT.json')
    out.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__': main()
