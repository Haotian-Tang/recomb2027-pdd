#!/usr/bin/env python3
"""Exact controls for the <=2 strict-source envelope lemma.

Standard library / Fraction only; no frozen code imported. The finite controls
are not a proof of the general statement and do not establish novelty.
"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
from datetime import datetime, timezone
import hashlib,json,time

ROOT=Path(__file__).resolve().parents[1]


def geom(edges):
    V=sorted({v for e in edges for v in e[:2]})
    adj={v:{} for v in V}
    for u,v,l in edges:
        assert l>0
        adj[u][v]=adj[v][u]=F(l)
    R=[v for v in V if len(adj[v])==1]
    D={}
    for r in V:
        D[r]={r:F(0)};stack=[r]
        for v in stack:
            for u,l in adj[v].items():
                if u not in D[r]:D[r][u]=D[r][v]+l;stack.append(u)
    return V,adj,R,D


def evaluate(edges,sources,y,D):
    assert y>0 and sum(s['p'] for s in sources)==1
    mass={v:F(0) for v in D}
    for s in sources:
        u,v,l=edges[s['edge']];p,x=s['p'],s['x']
        assert p>0 and 0<x<1
        mass[u]+=p*(1-x);mass[v]+=p*x
    f={r:y+sum(mass[v]*D[r][v] for v in D) for r in D}
    for r in D:
        direct=y+sum(s['p']*min(D[r][edges[s['edge']][0]]+s['x']*edges[s['edge']][2],
                              D[r][edges[s['edge']][1]]+(1-s['x'])*edges[s['edge']][2])
                     for s in sources)
        assert direct==f[r]
    return mass,f


def two_pass_envelope(adj,observed):
    """max_r d_r-dist(v,r) using two linear-work passes; None means -infinity."""
    root=min(adj);parent={root:None};order=[root]
    for v in order:
        for u in adj[v]:
            if u!=parent[v]:parent[u]=v;order.append(u)
    down={}
    for v in reversed(order):
        candidates=[observed[v]] if v in observed else []
        candidates.extend(down[u]-l for u,l in adj[v].items()
                          if parent.get(u)==v and down[u] is not None)
        down[v]=max(candidates) if candidates else None
    up={root:None};out={};operations=0
    for v in order:
        items=[]
        if v in observed:items.append((observed[v],'self'))
        if up[v] is not None:items.append((up[v],'up'))
        for u,l in adj[v].items():
            if parent.get(u)==v and down[u] is not None:items.append((down[u]-l,u))
        # Retain the two largest entries without sorting a high-degree list.
        best=[]
        for item in items:
            operations+=1
            if not best or item[0]>best[0][0]:best.insert(0,item)
            elif len(best)==1 or item[0]>best[1][0]:best.insert(1,item)
            if len(best)>2:best.pop()
        out[v]=best[0][0]
        for u,l in adj[v].items():
            if parent.get(u)!=v:continue
            choice=next((value for value,key in best if key!=u),None)
            up[u]=choice-l if choice is not None else None
    return out,operations


def inverse(adj,f):
    m={v:(2-len(adj[v])-sum((f[v]-f[u])/l for u,l in adj[v].items()))/2 for v in adj}
    y=(sum((2-len(adj[v]))*f[v] for v in adj)
       -sum(l for v in adj for u,l in adj[v].items() if v<u))/2
    return m,y


def covers(edges,m):
    P={v for v in m if m[v]>0}
    allowed=[i for i,(u,v,l) in enumerate(edges) if u in P and v in P]
    return [list(Fset) for k in (1,2) for Fset in combinations(allowed,k)
            if {v for i in Fset for v in edges[i][:2]}==P]


def source_split(edges,m,selected):
    degree={v:sum(v in edges[i][:2] for i in selected) for v in m}
    out=[]
    for i in selected:
        u,v,l=edges[i];a,b=m[u]/degree[u],m[v]/degree[v]
        out.append({'edge':i,'p':a+b,'x':b/(a+b)})
    return out


def solve_at_most_two(edges,observed):
    # The solver itself does not construct the all-pairs distance matrix used
    # by the independent controls. Input is a valid positive-length tree.
    adj={}
    for u,v,l in edges:
        assert l>0
        adj.setdefault(u,{})[v]=F(l);adj.setdefault(v,{})[u]=F(l)
    V=list(adj);R=[v for v in V if len(adj[v])==1]
    assert set(observed)==set(R) and all(len(adj[v])!=2 for v in V)
    f,op=two_pass_envelope(adj,observed)
    if any(f[r]!=observed[r] for r in R):return {'feasible':False,'reason':'envelope_changes_observed_leaf'}
    m,y=inverse(adj,f)
    if min(m.values())<0 or y<=0:return {'feasible':False,'reason':'nonphysical_envelope','mass':m,'y':y}
    if sum(x>0 for x in m.values())>4:return {'feasible':False,'reason':'more_than_four_positive_vertices'}
    cs=covers(edges,m)
    if not cs:return {'feasible':False,'reason':'no_strict_at_most_two_cover','mass':m,'y':y}
    assert len(cs)==1
    sources=source_split(edges,m,cs[0])
    m2={v:F(0) for v in V}
    for s in sources:
        u,v,l=edges[s['edge']]
        m2[u]+=s['p']*(1-s['x']);m2[v]+=s['p']*s['x']
    f2={v:y for v in V}
    # At most four positive endpoints, hence at most four linear tree scans.
    for r in V:
        if m2[r]==0:continue
        distance={r:F(0)};queue=[r]
        for v in queue:
            for u,l in adj[v].items():
                if u not in distance:distance[u]=distance[v]+l;queue.append(u)
        for v in V:f2[v]+=m2[r]*distance[v]
    assert m2==m and f2==f
    return {'feasible':True,'sources':sources,'source_edges':cs[0],'mass':m,'y':y,
            'all_vertex_potential':f,'forward_work_items':op}


def positive_control(name,edges,selected,p,x):
    V,adj,R,D=geom(edges)
    sources=[{'edge':i,'p':a,'x':b} for i,a,b in zip(selected,p,x)]
    y=F(1,10**9)
    m,f=evaluate(edges,sources,y,D);observed={r:f[r] for r in R}
    direct={v:max(observed[r]-D[v][r] for r in R) for v in V}
    fast,ops=two_pass_envelope(adj,observed)
    assert direct==fast==f
    m2,y2=inverse(adj,fast);assert m2==m and y2==y
    cs=covers(edges,m);assert cs==[list(selected)]
    sol=solve_at_most_two(edges,observed);assert sol['source_edges']==list(selected)
    return {'name':name,'edges':edges,'sources':sources,'pendant':y,'observed':observed,
            'true_endpoint_mass':m,'true_potential':f,'envelope':fast,
            'strict_covers_up_to_two':cs,'reconstructed':sol,'passed':True}


def named_boundaries():
    star=[(0,3,F(1)),(1,3,F(1)),(2,3,F(1))]
    V,adj,R,D=geom(star)
    # x is from leaf towards the center for these edge orientations.
    left=[{'edge':0,'p':F(1,2),'x':F(1,2)}, {'edge':1,'p':F(1,2),'x':F(1,2)}]
    right=[{'edge':0,'p':F(1,3),'x':F(1,4)}, {'edge':1,'p':F(2,3),'x':F(5,8)}]
    m0,f0=evaluate(star,left,F(1),D);m1,f1=evaluate(star,right,F(1),D)
    assert m0==m1 and f0==f1
    adjacent={'edges':star,'worlds':[left,right],'pendant':F(1),
              'observed':{r:f0[r] for r in R},'endpoint_mass':m0,
              'scope':'Supports unique among <=2 strict sources, but weights and locations differ even with the same mean pendant.'}
    # Root-supplied three-source quartet control, checked from original distances.
    quartet=[(4,0,F(1)),(4,1,F(1)),(4,5,F(1)),(5,2,F(1)),(5,3,F(1))]
    Vq,aq,Rq,Dq=geom(quartet)
    w0=[{'edge':0,'p':F(1,6),'x':F(3,4)}, {'edge':1,'p':F(1,6),'x':F(3,4)},
        {'edge':2,'p':F(2,3),'x':F(15,16)}]
    w1=[{'edge':3,'p':F(1,6),'x':F(3,4)}, {'edge':4,'p':F(1,6),'x':F(3,4)},
        {'edge':2,'p':F(2,3),'x':F(1,16)}]
    mq0,fq0=evaluate(quartet,w0,F(11,8),Dq);mq1,fq1=evaluate(quartet,w1,F(11,8),Dq)
    dq={r:F(3) for r in Rq}
    assert {r:fq0[r] for r in Rq}=={r:fq1[r] for r in Rq}==dq
    fqmax,_=two_pass_envelope(aq,dq)
    assert fq0!=fq1 and fqmax!=fq0 and fqmax!=fq1
    three={'edges':quartet,'worlds':[w0,w1],'pendant':F(11,8),'observed':dq,
           'endpoint_masses':[mq0,mq1],'true_potentials':[fq0,fq1],'envelope':fqmax,
           'scope':'k=3 has two different strict supports even with the same weight multiset and positive mean pendant; no >2 identifiability claim.'}
    # Missing reference: one source on unobserved C arm versus two on A/B.
    c=[{'edge':2,'p':F(1),'x':F(1,2)}]
    mc,fc=evaluate(star,c,F(1,2),D)
    assert {r:fc[r] for r in (0,1)}=={r:f0[r] for r in (0,1)}=={0:F(2),1:F(2)}
    assert fc[2]!=f0[2]
    missing={'edges':star,'references':[0,1],'worlds':[c,left], 'pendants':[F(1,2),F(1)],
             'observations':{0:F(2),1:F(2)},'omitted_leaf_observations':[fc[2],f0[2]],
             'scope':'The all-leaf assumption is necessary for this conclusion.'}
    # Degree-two path: two equal-weight placements on different pairs of edges.
    path=[(0,1,F(1)),(1,2,F(1)),(2,3,F(1)),(3,4,F(1))]
    Vp,ap,Rp,Dp=geom(path)
    p0=[{'edge':0,'p':F(1,2),'x':F(1,2)}, {'edge':3,'p':F(1,2),'x':F(1,2)}]
    p1=[{'edge':1,'p':F(1,2),'x':F(1,2)}, {'edge':2,'p':F(1,2),'x':F(1,2)}]
    _,fp0=evaluate(path,p0,F(1),Dp);_,fp1=evaluate(path,p1,F(1),Dp)
    assert {r:fp0[r] for r in Rp}=={r:fp1[r] for r in Rp}=={0:F(3),4:F(3)}
    subdivision={'edges':path,'worlds':[p0,p1],'pendant':F(1),'observed':{0:F(3),4:F(3)},
                 'scope':'Degree-two vertices can destroy the edge-support conclusion; the reference is not reduced.'}
    inconsistent=solve_at_most_two(star,{0:F(100),1:F(1),2:F(1)})
    assert not inconsistent['feasible'] and inconsistent['reason']=='envelope_changes_observed_leaf'
    equal_three={r:F(2) for r in R};two_no=solve_at_most_two(star,equal_three)
    assert not two_no['feasible']
    return {'adjacent_nonidentifiability':adjacent,'three_source_nonidentifiability':three,
            'missing_leaf_failure':missing,'degree_two_failure':subdivision,
            'invalid_data_rejection':inconsistent,'strict_star_at_most_two_rejection':two_no}


def main():
    start=time.perf_counter()
    trees={
        'single_edge':[(0,1,F(2))],
        'three_star':[(0,3,F(1)),(1,3,F(2)),(2,3,F(3))],
        'quartet':[(0,4,F(1,3)),(1,4,F(2)),(4,5,F(3,2)),(2,5,F(5,3)),(3,5,F(7,4))],
        'five_caterpillar':[(0,5,F(1)),(1,5,F(2)),(5,6,F(1,7)),(2,6,F(4,3)),
                            (6,7,F(7,3)),(3,7,F(2,5)),(4,7,F(5))],
        'four_star':[(0,4,F(1,1000)),(1,4,F(2)),(2,4,F(17)),(3,4,F(5,2))]}
    cases=[]
    settings=[([F(1,2),F(1,2)],[F(1,3),F(3,4)]),
              ([F(1,3),F(2,3)],[F(2,5),F(1,7)]),
              ([F(1,101),F(100,101)],[F(1,1000),F(999,1000)])]
    for name,edges in trees.items():
        for i in range(len(edges)):
            cases.append(positive_control(name+'_one_'+str(i),edges,[i],[F(1)],[F(2,5)]))
        for pair in combinations(range(len(edges)),2):
            for j,(p,x) in enumerate(settings):
                cases.append(positive_control(name+'_pair_'+str(pair)+'_'+str(j),edges,list(pair),p,x))
    output={'passed':True,'created_utc':datetime.now(timezone.utc).isoformat(),
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'scope':'Exact finite controls for a separately stated analytic <=2-source lemma. '
                    'Not a prior/novelty certificate, biological validation, or general-k recovery result.',
            'case_count':len(cases),'tree_count':len(trees),'cases':cases,
            'boundary_controls':named_boundaries(),'elapsed_seconds':time.perf_counter()-start}
    (ROOT/'results/stage23_two_source_envelope.json').write_text(json.dumps(output,indent=2,default=str)+'\n')
    print(json.dumps({k:v for k,v in output.items() if k not in ('cases','boundary_controls')},indent=2))
    return output


if __name__=='__main__':main()
