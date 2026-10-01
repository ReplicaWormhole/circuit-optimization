"""Bounded exact-preserving chronological gate-list rewriting."""
import json,time,hashlib,re
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4];HERE=Path(__file__).resolve().parent
def support(g):return {g['control'],g['target']} if g['gate']=='cx' else {g['qubit']}
def commute(a,b):
    if not support(a)&support(b):return True
    if a['gate']==b['gate']=='cx':return a['control']!=b['target'] and b['control']!=a['target']
    if a['gate']=='cx' and b['gate']!='cx':return (b['gate']=='rz' and b['qubit']==a['control']) or (b['gate']=='rx' and b['qubit']==a['target'])
    if b['gate']=='cx':return commute(b,a)
    return a['gate']==b['gate'] and a['gate'] in ('rx','ry','rz') and a['qubit']==b['qubit']
def angle(g):
    s=g.get('theta');assert isinstance(s,str)
    m=re.fullmatch(r'([+-]?\d+)\*pi(?:/(\d+))?',s);assert m,s
    return Fraction(int(m[1]),int(m[2] or 1))
def main():
    cfg=json.loads((HERE/'config.json').read_text());assert not (HERE/'result.json').exists()
    p=ROOT/cfg['source'];assert hashlib.sha256(p.read_bytes()).hexdigest()==cfg['source_sha256'];source=json.loads(p.read_text());gates=source['gates'].copy();trace=[];checks=0;begun=time.monotonic();stop='fixpoint'
    for sweep in range(100):
        changed=False
        for i,a in enumerate(gates):
            for j in range(i+1,len(gates)):
                checks+=1
                if checks>50000 or time.monotonic()-begun>120:stop='rulecap' if checks>50000 else 'wallcap';break
                b=gates[j];replacement=None;rule=None
                if a==b and a['gate'] in ('cx','h'):replacement=[];rule='involution cancellation'
                elif a['gate']==b['gate'] and a['gate'] in ('rx','ry','rz') and a['qubit']==b['qubit']:
                    value=angle(a)+angle(b);replacement=[] if not value else [{**a,'theta':f'{value.numerator}*pi/{value.denominator}'}];rule='exact sameaxis merge'
                if rule:
                    trace.append({'sweep':sweep,'i':i,'j':j,'rule':rule,'a':a,'b':b,'crossed':gates[i+1:j],'replacement':replacement})
                    gates=gates[:i]+replacement+gates[i+1:j]+gates[j+1:];changed=True;break
                if not commute(a,b):break
            if changed or stop!='fixpoint':break
        if stop!='fixpoint' or not changed:break
    else:stop='sweepcap'
    candidate={'n':4,'gates':gates,'provenance':{'source':cfg['source'],'source_sha256':cfg['source_sha256'],'scope':'Exact-preserving reviewed local rewrite; independentnewgate certificate pending'}}
    p=HERE/'candidate.json';p.write_text(json.dumps(candidate,indent=2)+'\n')
    result={'config':cfg,'elapsed_seconds':time.monotonic()-begun,'rule_checks':checks,'stop':stop,'trace':trace,'before_gates':len(source['gates']),'after_gates':len(gates),'before_cx':sum(g['gate']=='cx' for g in source['gates']),'after_cx':sum(g['gate']=='cx' for g in gates),'candidate':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'scope':'Exact algebraic rewrite rules; new saved gate list needs independent full exact certificate; no optimalityclaim'}
    (HERE/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['config','trace']}))
if __name__=='__main__':main()
