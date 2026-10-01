"""Read-only comparison of explicit chronological unordered CX pair sequences."""
import json,sqlite3,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4];HERE=Path(__file__).resolve().parent
cfg=json.loads((ROOT/'collaboration/work/numerical_label12/curvature/config.json').read_text())
base=cfg['optimizer_schedule'];proposals=[]
for slot,edge in [(6,[0,3]),(8,[0,2])]:
    s=[None if e is None else list(e) for e in base];s[slot]=edge
    proposals.append({'changed_optimizer_slot':slot,'edge':edge,'schedule':s,'unordered_native':[sorted(e) for e in s if e is not None]})
matches=[[] for _ in proposals];scanned=0
def walk(obj,origin):
    global scanned
    if isinstance(obj,list):
        if 6<=len(obj)<=18 and all(e is None or isinstance(e,(list,tuple)) and len(e)==2 and all(type(q) is int and 0<=q<4 for q in e) and e[0]!=e[1] for e in obj):
            seq=[sorted(e) for e in obj if e is not None];scanned+=1
            for k,p in enumerate(proposals):
                if seq==p['unordered_native']:matches[k].append(origin)
        else:
            for i,v in enumerate(obj):walk(v,origin+f'[{i}]')
    elif isinstance(obj,dict):
        if isinstance(obj.get('gates'),list):
            edges=[[g['control'],g['target']] for g in obj['gates'] if isinstance(g,dict) and g.get('gate')=='cx'];walk(edges,origin+'.gates')
        for key,v in obj.items():
            if key!='gates':walk(v,origin+'.'+key)
con=sqlite3.connect(ROOT/'experiments.sqlite3');rows=list(con.execute('select id,config_json from attempts'))
for run,text in rows:walk(json.loads(text),f'ledger:{run}')
paths={p for p in ROOT.rglob('*.json') if '.git' not in p.parts and '__pycache__' not in p.parts}
for p in sorted(paths):
    if HERE in p.parents or ROOT/'collaboration/work/verifier_label12/architecture' in p.parents:continue
    try:walk(json.loads(p.read_text()),str(p.relative_to(ROOT)))
    except (ValueError,KeyError,TypeError):pass
for p,m in zip(proposals,matches):p['matches']=sorted(set(m))
out={'ledger_rows':len(rows),'json_files':len(paths),'explicit_sequence_observations':scanned,'proposals':proposals,'scope':'Explicit JSON schedule/gate arrays in all ledgerconfigs, all repository JSON artifacts. Does not infer implicit topology-generating rules or prove universal historical novelty.'}
(HERE/'duplicate_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
