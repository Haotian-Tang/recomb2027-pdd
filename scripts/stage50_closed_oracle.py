#!/usr/bin/env python3
"""Independent exact closed-edge oracle for small original-tree cases.

No project imports. Enumerates original-edge supports and nonnegative basic
solutions to normalization and observed-distance equations. Vertices of the
original tree are allowed. Call in a subprocess with an external deadline.
"""
from fractions import Fraction as F
from itertools import combinations


def eliminate(matrix, width):
    a=[list(map(F,row)) for row in matrix]
    pivots=[]
    for column in range(width):
        found=next((i for i in range(len(pivots),len(a)) if a[i][column]),None)
        if found is None:
            continue
        rank=len(pivots)
        a[rank],a[found]=a[found],a[rank]
        scale=a[rank][column]
        a[rank]=[x/scale for x in a[rank]]
        for i in range(len(a)):
            if i!=rank and a[i][column]:
                scale=a[i][column]
                a[i]=[x-scale*y for x,y in zip(a[i],a[rank])]
        pivots.append(column)
    if any(not any(row[:width]) and row[width] for row in a):
        return None
    return a[:len(pivots)],pivots


def nonnegative_point(equations,n,statistics):
    reduced=eliminate(equations,n)
    if reduced is None:
        return None
    rows,_=reduced
    rank=len(rows)
    # sum(z)=1 makes the feasible nonnegative polyhedron compact.
    for columns in combinations(range(n),rank):
        statistics['bases_examined']+=1
        basic=eliminate([[row[j] for j in columns]+[row[-1]] for row in rows],rank)
        if basic is None or len(basic[1])!=rank:
            continue
        values=[row[-1] for row in basic[0]]
        if any(x<0 for x in values):
            continue
        point=[F(0)]*n
        for j,x in zip(columns,values):
            point[j]=x
        assert all(sum(a*x for a,x in zip(row[:-1],point))==row[-1] for row in equations)
        return point
    return None


def exact_closed_minimum(case):
    E=[(u,v,F(ell)) for u,v,ell in case['edges']]
    adj={}
    for u,v,ell in E:
        assert u!=v and ell>0
        assert v not in adj.setdefault(u,{})
        adj[u][v]=ell
        adj.setdefault(v,{})[u]=ell
    assert len(E)==len(adj)-1
    obs={v:F(d) for v,d in case['observed'].items()}
    y=F(case['pendant'])
    assert y>0
    paths={}
    for r in obs:
        assert len(adj[r])==1
        d={r:F(0)};todo=[r]
        for u in todo:
            for v,ell in adj[u].items():
                if v not in d:
                    d[v]=d[u]+ell;todo.append(v)
        assert len(d)==len(adj)
        paths[r]=d
    stats=dict(supports_examined=0,bases_examined=0)
    for k in range(1,min(case['cap'],len(E))+1):
        for support in combinations(range(len(E)),k):
            stats['supports_examined']+=1
            endpoints=[v for i in support for v in E[i][:2]]
            equations=[[F(1)]*(2*k)+[F(1)]]
            equations += [[paths[r][v] for v in endpoints]+[obs[r]-y] for r in obs]
            point=nonnegative_point(equations,2*k,stats)
            if point is None:
                continue
            sources=[]
            for j,index in enumerate(support):
                a,b=point[2*j:2*j+2]
                if a+b:
                    sources.append(dict(edge=index,p=a+b,x=b/(a+b)))
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


def serial(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {str(k):serial(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [serial(v) for v in x]
    return x
