"""Select a six-CNOT survivor with two short selected tails, one AA."""

from collections import Counter
from itertools import product
import json
from pathlib import Path

import global_lower_bound7_schedule_filter as prior
import global_lower_bound7_degree3_tail_screen as tails
import global_lower_bound7_adjacent_pair_support_filter as pair_filter

ROOT = Path(__file__).resolve().parent


def surviving(schedule):
    degree = prior.degrees(schedule)
    if min(degree)<2 or not prior.all_forward_cones_full(schedule):
        return False
    if prior.tail_rejection(schedule,degree) is not None:
        return False
    return not any(len(cone)==2 and tuple(sorted(cone)) not in prior.OPPOSITE
                   for j in range(4)
                   for cone,_ in [pair_filter.first_commuting_cone(schedule,j)])


def tail_geometry(schedule,j):
    active,_=tails.first_commuting_tail(schedule,j)
    if len(active)!=2:
        return None
    k=next(w for w in active[0][1] if w!=j)
    ell=next(w for w in active[1][1] if w!=j)
    return (k,ell,tails.geometry(j,active))


def main():
    counts=Counter()
    examples={}
    for schedule in product(prior.EDGES,repeat=6):
        if not surviving(schedule) or prior.degrees(schedule)!=(3,3,3,3):
            continue
        aa=tail_geometry(schedule,0)
        if aa is None or aa[:2]!=(1,3):
            continue
        for h in (1,2,3):
            other=tail_geometry(schedule,h)
            if other is None:
                continue
            key=other[2]
            counts[key]+=1
            examples.setdefault(key,{'schedule':[list(e) for e in schedule],
                                     'aa_wire':0,'aa_ordered_partners':[1,3],
                                     'second_wire':h,'second_ordered_partners':list(other[:2])})
    result={'degree_3333_aa0_13_other_short_tail_occurrences':dict(sorted(counts.items())),
            'examples_by_second_geometry':examples}
    (ROOT/'global_lower_bound7_joint_aa_schedule_select_result.json').write_text(
        json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result,sort_keys=True))


if __name__=='__main__':
    main()
