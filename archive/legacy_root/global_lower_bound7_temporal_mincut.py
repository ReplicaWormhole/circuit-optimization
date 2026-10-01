"""Tensor-network min-cut Schmidt-rank screen for 401 six-CNOT D4 classes.

Each CNOT is overapproximated by an arbitrary two-qubit tensor.  Contract
fixed product input legs into their first gate tensors; every remaining
wire segment has dimension two.  A cut of c segments between output sides
bounds every circuit-column Schmidt rank by 2**c, for all local gates and
CNOT orientations.  An adjacent balanced output cut with c <= 1 conflicts
with the i-eigenspace theorem in global_lower_bound7_cut_schmidt.md.
"""

from collections import Counter
import json
from pathlib import Path

import networkx as nx

ROOT=Path(__file__).resolve().parent
CUTS=((0,1),(0,3))
INF=1000


def add_undirected_bond(graph,a,b):
    for u,v in ((a,b),(b,a)):
        capacity=graph[u][v]['capacity']+1 if graph.has_edge(u,v) else 1
        graph.add_edge(u,v,capacity=capacity)


def cut_bound(schedule,left):
    graph=nx.DiGraph()
    previous=[None]*4
    for t,edge in enumerate(schedule):
        gate=('g',t)
        graph.add_node(gate)
        for q in edge:
            if previous[q] is not None:
                add_undirected_bond(graph,previous[q],gate)
            previous[q]=gate
    source=('source',)
    sink=('sink',)
    graph.add_nodes_from((source,sink))
    for q in range(4):
        output=('out',q)
        graph.add_node(output)
        if previous[q] is not None:
            add_undirected_bond(graph,previous[q],output)
        if q in left:
            graph.add_edge(source,output,capacity=INF)
        else:
            graph.add_edge(output,sink,capacity=INF)
    value,(source_side,sink_side)=nx.minimum_cut(graph,source,sink,capacity='capacity')
    assert value<=2
    return int(value),sorted(node[1] for node in source_side
                             if len(node)==2 and node[0]=='g')


def main():
    data=json.loads((ROOT/'global_lower_bound7_dihedral_orbits_result.json').read_text())
    assert data['class_count']==401 and data['source_survivors']==3208
    partition=Counter()
    by_degree=Counter()
    examples={}
    per_class=[]
    for row in data['representatives']:
        schedule=tuple(tuple(edge) for edge in row['representative'])
        bounds={''.join(map(str,left)):cut_bound(schedule,left) for left in CUTS}
        minimum=min(value for value,_ in bounds.values())
        category='excluded_rank_at_most_two' if minimum<=1 else 'survives_temporal_cut'
        partition[category]+=1
        by_degree[(category,row['degree_sequence'])]+=1
        examples.setdefault(category,{'schedule':row['representative'],
                                      'bounds':{k:v[0] for k,v in bounds.items()}})
        per_class.append({'representative':row['representative'],
                          'degree_sequence':row['degree_sequence'],
                          'cut_min_edges':{k:v[0] for k,v in bounds.items()},
                          'category':category})
    assert sum(partition.values())==401
    result={'source_classes':401,'source_schedules':3208,
            'exclusive_class_partition':dict(sorted(partition.items())),
            'exclusive_schedule_partition':{k:8*v for k,v in sorted(partition.items())},
            'by_degree_sequence':{category:{degree:by_degree[(category,degree)]
                                            for degree in ('3333','4332','4422','5322')
                                            if by_degree[(category,degree)]}
                                  for category in sorted(partition)},
            'examples':examples,'classes':per_class}
    (ROOT/'global_lower_bound7_temporal_mincut_result.json').write_text(
        json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:result[k] for k in
                      ('exclusive_class_partition','exclusive_schedule_partition','by_degree_sequence')},
                     sort_keys=True))


if __name__=='__main__':main()
