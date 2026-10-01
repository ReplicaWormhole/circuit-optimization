"""CNOT-factorized tensor-network min-cut screen for 401 D4 classes.

Each CNOT has an exact operator-Schmidt decomposition with two terms.
Represent it by two wire-local tensors joined by one dimension-two bond.
All worldline segments are dimension two; fixed product input legs are
contracted.  A cut of c bonds between the two output sides bounds every
circuit-column Schmidt rank by 2**c for every direction and local gate.
"""

from collections import Counter
import json
from pathlib import Path

import networkx as nx

ROOT=Path(__file__).resolve().parent
CUTS=((0,1),(0,3))
INF=1000


def bond(graph,a,b):
    for u,v in ((a,b),(b,a)):
        capacity=graph[u][v]['capacity']+1 if graph.has_edge(u,v) else 1
        graph.add_edge(u,v,capacity=capacity)


def mincut(schedule,left):
    graph=nx.DiGraph()
    previous=[None]*4
    for t,(a,b) in enumerate(schedule):
        va=('g',t,a)
        vb=('g',t,b)
        bond(graph,va,vb)
        for q,v in ((a,va),(b,vb)):
            if previous[q] is not None:
                bond(graph,previous[q],v)
            previous[q]=v
    source=('source',)
    sink=('sink',)
    graph.add_nodes_from((source,sink))
    for q in range(4):
        output=('out',q)
        graph.add_node(output)
        if previous[q] is not None:
            bond(graph,previous[q],output)
        if q in left:
            graph.add_edge(source,output,capacity=INF)
        else:
            graph.add_edge(output,sink,capacity=INF)
    value,(source_side,sink_side)=nx.minimum_cut(graph,source,sink,capacity='capacity')
    assert value<=2
    return int(value)


def main():
    data=json.loads((ROOT/'global_lower_bound7_dihedral_orbits_result.json').read_text())
    assert data['class_count']==401
    classes=[]
    partition=Counter()
    by_degree=Counter()
    examples={}
    for row in data['representatives']:
        schedule=tuple(tuple(edge) for edge in row['representative'])
        values={''.join(map(str,left)):mincut(schedule,left) for left in CUTS}
        category='excluded_rank_at_most_two' if min(values.values())<=1 else 'survives_cnot_bond_cut'
        partition[category]+=1
        by_degree[(category,row['degree_sequence'])]+=1
        examples.setdefault(category,{'schedule':row['representative'],'cut_values':values})
        classes.append({'representative':row['representative'],'degree_sequence':row['degree_sequence'],
                        'cut_values':values,'category':category})
    assert sum(partition.values())==401
    result={'source_classes':401,'source_schedules':3208,
            'exclusive_class_partition':dict(sorted(partition.items())),
            'exclusive_schedule_partition':{k:8*v for k,v in sorted(partition.items())},
            'by_degree_sequence':{category:{degree:by_degree[(category,degree)]
                                            for degree in ('3333','4332','4422','5322')
                                            if by_degree[(category,degree)]}
                                  for category in sorted(partition)},
            'examples':examples,'classes':classes}
    (ROOT/'global_lower_bound7_cnot_bond_mincut_result.json').write_text(
        json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:result[k] for k in
                      ('exclusive_class_partition','exclusive_schedule_partition','by_degree_sequence')},
                     sort_keys=True))


if __name__=='__main__':main()
