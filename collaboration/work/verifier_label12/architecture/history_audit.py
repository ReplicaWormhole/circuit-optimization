"""Read-only explicit prior-config/archived-CX-sequence audit, not a search."""
import json
import sqlite3
from pathlib import Path

ROOT=Path(__file__).resolve().parents[4]
BASE=[(0,1),None,(1,2),(1,0),(2,3),(2,1),(0,1),(3,2),(1,2),(1,0),(2,3),(2,1),(0,1)]

def canonical(sequence):
    return tuple(tuple(sorted(pair)) for pair in sequence if pair is not None)

def sequences(obj,path='root'):
    found=[]
    if isinstance(obj,list):
        if len(obj) in (12,13) and all(x is None or (isinstance(x,(list,tuple)) and len(x)==2 and all(type(q) is int and 0<=q<4 for q in x) and x[0]!=x[1]) for x in obj):
            if sum(x is not None for x in obj)==12:found.append((path,canonical(obj)))
        for k,v in enumerate(obj):found.extend(sequences(v,f'{path}[{k}]'))
    elif isinstance(obj,dict):
        for k,v in obj.items():found.extend(sequences(v,f'{path}.{k}'))
    return found

wanted=[]
for slot,pair in ((6,(0,3)),(8,(0,2))):
    schedule=BASE.copy();schedule[slot]=pair
    wanted.append({'slot':slot,'pair':list(pair),'unordered':canonical(schedule),'matches':[]})
with sqlite3.connect(ROOT/'experiments.sqlite3') as db:
    rows=list(db.execute('select id,config_json,candidate_path from attempts where id<=378'))
for run,config,candidate in rows:
    for field,seq in sequences(json.loads(config)):
        for target in wanted:
            if seq==target['unordered']:target['matches'].append({'run':run,'kind':'explicit_config_sequence','field':field})
    if candidate:
        path=ROOT/candidate
        if path.exists():
            data=json.loads(path.read_text())
            cx=[(g['control'],g['target']) for g in data.get('gates',[]) if g.get('gate')=='cx']
            if len(cx)==12:
                for target in wanted:
                    if canonical(cx)==target['unordered']:target['matches'].append({'run':run,'kind':'archived_candidate','path':candidate})
report={'scope':'All378prior ledger configs explicit12CX pair sequences and archived candidate paths; dynamic generated unarchived schedules and encoded rules not exhaustively reconstructed.', 'prior_runs':len(rows),'targets':wanted}
Path(__file__).with_name('history_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
