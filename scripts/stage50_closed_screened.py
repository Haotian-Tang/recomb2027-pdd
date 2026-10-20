#!/usr/bin/env python3
"""Exact original-edge closed oracle with generic convex-range screening.

Screening uses only observed coordinates and coordinate differences. For any
normalized nonnegative mixture, a linear observable lies between its smallest
and largest selected endpoint value. This excludes no feasible closed support.
No reduction metadata or tree-specific hardness lemma is used.
"""
from fractions import Fraction as F
from itertools import combinations
from stage50_closed_oracle import nonnegative_point, serial


def exact_closed_minimum(case):
    E=[(u,v,F(ell)) for u,v,ell in case['edges']]
    adj={}
    for u,v,ell in E:
        assert u!=v and ell>0 and v not in adj.setdefault(u,{})
        adj[u][v]=ell
        adj.setdefault(v,{})[u]=ell
    assert len(E)==len(adj)-1
    obs={v:F(d) for v,d in case['observed'].items()}
    y=F(case['pendant']); assert y>0
    paths={}
    for r in obs:
        assert len(adj[r])==1
        d={r:F(0)};todo=[r]
        for u in todo:
            for v,ell in adj[u].items():
                if v not in d:d[v]=d[u]+ell;todo.append(v)
        assert len(d)==len(adj)
        paths[r]=d
    # Contrasts first: inexpensive checks often expose incompatible supports.
    observables=[({v:paths[r][v]-paths[s][v] for v in adj},obs[r]-obs[s])
                 for r,s in combinations(obs,2)]
    observables += [(paths[r],obs[r]-y) for r in obs]
    ranges=[([(min(values[u],values[v]),max(values[u],values[v])) for u,v,_ in E],target)
            for values,target in observables]
    stats=dict(supports_examined=0,range_rejections=0,lp_supports=0,bases_examined=0)
    for k in range(1,min(case['cap'],len(E))+1):
        for support in combinations(range(len(E)),k):
            stats['supports_examined']+=1
            if any(target<min(intervals[i][0] for i in support) or target>max(intervals[i][1] for i in support)
                   for intervals,target in ranges):
                stats['range_rejections']+=1
                continue
            stats['lp_supports']+=1
            # Duplicate endpoint columns can be merged in the closed model.
            vertices=list(dict.fromkeys(v for i in support for v in E[i][:2]))
            equations=[[F(1)]*len(vertices)+[F(1)]]
            equations += [[paths[r][v] for v in vertices]+[obs[r]-y] for r in obs]
            point=nonnegative_point(equations,len(vertices),stats)
            if point is None:continue
            allocations={i:[F(0),F(0)] for i in support}
            for vertex,mass in zip(vertices,point):
                if mass:
                    index=next(i for i in support if vertex in E[i][:2])
                    allocations[index][0 if E[index][0]==vertex else 1]+=mass
            sources=[dict(edge=i,p=a+b,x=b/(a+b)) for i,(a,b) in allocations.items() if a+b]
            assert sum(s['p'] for s in sources)==1
            for r,expected in obs.items():
                value=y
                for s in sources:
                    u,v,ell=E[s['edge']]
                    assert s['p']>0 and 0<=s['x']<=1
                    value+=s['p']*min(paths[r][u]+s['x']*ell,paths[r][v]+(1-s['x'])*ell)
                assert value==expected
            return dict(minimum=len(sources),witness=sources,statistics=stats)
    return dict(minimum=None,witness=None,statistics=stats)
