#!/usr/bin/env python3
"""Exact <=3-source solver by at most one latent potential deficit.

The solver uses no optimizer and no all-pairs matrix. Its finite control harness
uses frozen geometry plus an independent old full endpoint-polytope oracle.
"""
from fractions import Fraction as F
from itertools import combinations
from collections import deque, Counter
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,time,random
from stage23_two_source_envelope import two_pass_envelope, inverse, source_split

ROOT=Path(__file__).resolve().parents[1]

def adjacency(edges):
    a={}
    for u,v,l in edges:
        assert u!=v and l>0
        assert v not in a.get(u,{})
        a.setdefault(u,{})[v]=F(l);a.setdefault(v,{})[u]=F(l)
    seen={next(iter(a))};todo=list(seen)
    for v in todo:
        for u in a[v]:
            if u not in seen:seen.add(u);todo.append(u)
    assert len(seen)==len(a) and len(edges)==len(a)-1
    assert all(len(x)!=2 for x in a.values())
    return a

def strict_cover(edges,mass):
    """Minimum edge cover of the positive-induced forest, or None for isolates."""
    positive={v for v,m in mass.items() if m>0}
    if not positive:return None
    neighbors={v:{} for v in positive}
    for i,(u,v,l) in enumerate(edges):
        if u in positive and v in positive:
            neighbors[u][v]=i;neighbors[v][u]=i
    if any(not ns for ns in neighbors.values()):return None
    active=set(positive);degree={v:len(ns) for v,ns in neighbors.items()}
    queue=deque(v for v in positive if degree[v]<=1)
    matched=set();chosen=[]
    while queue:
        v=queue.popleft()
        if v not in active:continue
        ns=[u for u in neighbors[v] if u in active]
        if not ns:
            active.remove(v);continue
        assert len(ns)==1
        u=ns[0];chosen.append(neighbors[v][u]);matched.update((u,v))
        active.remove(v);active.remove(u)
        for x in (u,v):
            for w in neighbors[x]:
                if w in active:
                    degree[w]-=1
                    if degree[w]<=1:queue.append(w)
    assert not active
    matching_count=len(chosen)
    for v in positive-matched:chosen.append(next(iter(neighbors[v].values())))
    chosen=list(dict.fromkeys(chosen))
    assert len(chosen)==len(positive)-matching_count
    assert {v for i in chosen for v in edges[i][:2]}==positive
    return chosen

def solve_three(edges,observed):
    """Rational input, all leaves, positive reduced tree. Minimum count if <=3."""
    adj=adjacency(edges);leaves={v for v in adj if len(adj[v])==1}
    assert set(observed)==leaves
    obs={v:F(d) for v,d in observed.items()}
    h,_=two_pass_envelope(adj,obs)
    if any(h[r]!=obs[r] for r in leaves):return {'feasible':False,'reason':'inconsistent_leaf_envelope'}
    base,y0=inverse(adj,h)
    if y0<=0:return {'feasible':False,'reason':'envelope_pendant_nonpositive'}
    examined=0
    def attempt(v,t):
        nonlocal examined
        examined+=1
        m=base.copy();y=y0
        if v is not None:
            m[v]-=sum(1/l for l in adj[v].values())*t/2
            for u,l in adj[v].items():m[u]+=t/(2*l)
            y-=t/2
        if y<=0 or any(a<0 for a in m.values()):return None
        support=strict_cover(edges,m)
        if support is None or len(support)>3:return None
        support.sort()  # At most three entries, constant work.
        f=h.copy()
        if v is not None:f[v]+=t
        return {'feasible':True,'minimum_sources_within_cap':len(support),'support':support,
                'mass':m,'pendant':y,'potential':f,'envelope':h,
                'deficit_vertex':v,'deficit':t,'sources':source_split(edges,m,support),
                'candidate_potentials_examined':examined}
    sol=attempt(None,F(0))
    if sol:return sol
    for v in adj:
        if len(adj[v])!=3 or base[v]<=0:continue
        upper=min(2*base[v]/sum(1/l for l in adj[v].values()),2*y0)
        cuts={F(0),upper}
        for u,l in adj[v].items():
            crossing=-2*l*base[u]
            if 0<crossing<upper:cuts.add(crossing)
        cuts=sorted(cuts)
        # Physical/sign status is constant on each open interval. Its midpoint
        # suffices, and nonnegative-mass endpoints are checked separately.
        candidates=sorted(set(cuts[1:]+[(a+b)/2 for a,b in zip(cuts,cuts[1:])]))
        for t in candidates:
            sol=attempt(v,t)
            if sol:return sol
    return {'feasible':False,'reason':'no_at_most_three_strict_representation',
            'candidate_potentials_examined':examined}

def check_witness(edges,obs,sol):
    from stage23_two_source_envelope import geom,evaluate
    V,a,R,D=geom(edges)
    if not sol['feasible']:return
    m,f=evaluate(edges,sol['sources'],sol['pendant'],D)
    assert {r:f[r] for r in R}==obs and m==sol['mass'] and f==sol['potential']
    assert len(sol['sources'])==sol['minimum_sources_within_cap']<=3

def geometry_bound(edges,sources,y):
    from stage23_two_source_envelope import geom,evaluate
    V,a,R,D=geom(edges);m,f=evaluate(edges,sources,y,D)
    obs={r:f[r] for r in R};h={v:max(obs[r]-D[v][r] for r in R) for v in V}
    deficits=[v for v in V if f[v]>h[v]]
    occupied={}
    for v in V:
        if len(a[v])==1:continue
        directions=set()
        for s in sources:
            u,w,l=edges[s['edge']]
            if v==u:directions.add(w);continue
            if v==w:directions.add(u);continue
            choices=[z for z,ll in a[v].items() if D[v][u]==ll+D[z][u]]
            assert len(choices)==1;directions.add(choices[0])
        occupied[v]=len(directions)
        assert (v in deficits)==(len(directions)==len(a[v]))
    budget=sum(len(a[v])-2 for v in deficits)
    assert budget<=max(len(sources)-2,0)
    return obs,{'deficit_vertices':deficits,'branch_budget':budget,'source_count':len(sources)}

def full_endpoint_oracle(edges,obs,cap=3):
    """Independent old coordinate model: enumerate full closed polytope vertices.

    A support is strict feasible iff its face has a point with every selected
    endpoint and y positive. Averaging all vertices of that face is complete.
    """
    import stage17_strict_bounds as old
    V,a,R,D=old.geom(edges);idx={v:i for i,v in enumerate(V)}
    A=[[D[r][v] for v in V]+[F(1)] for r in R]+[[F(1)]*len(V)+[F(0)]]
    vertices=old.bases(A,[obs[r] for r in R]+[F(1)])
    for k in range(1,cap+1):
        for selected in combinations(range(len(edges)),k):
            U={v for i in selected for v in edges[i][:2]}
            face=[x for x in vertices if all(x[idx[v]]==0 for v in V if v not in U)]
            if face and any(x[-1]>0 for x in face) and all(any(x[idx[v]]>0 for x in face) for v in U):return k
    return None

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    from stage23_two_source_envelope import geom,evaluate
    start=time.perf_counter();rng=random.Random(2413)
    trees={
      'star3':[(0,3,F(1)),(1,3,F(2)),(2,3,F(3))],
      'star4':[(0,4,F(1)),(1,4,F(2)),(2,4,F(3)),(3,4,F(4))],
      'quartet':[(0,4,F(1)),(1,4,F(2)),(4,5,F(3,2)),(2,5,F(4)),(3,5,F(5))],
      'caterpillar5':[(0,5,F(1)),(1,5,F(2)),(5,6,F(1,7)),(2,6,F(4,3)),(6,7,F(3,2)),(3,7,F(5)),(4,7,F(1,100))]}
    cases=[]
    for name,E in trees.items():
        for selected in combinations(range(len(E)),3):
            for config in (0,1):
                weights=[F(1,3)]*3 if config==0 else [F(1,10),F(3,10),F(3,5)]
                positions=[F(1,2)]*3 if config==0 else [F(1,1000),F(2,7),F(999,1000)]
                src=[{'edge':e,'p':p,'x':x} for e,p,x in zip(selected,weights,positions)]
                y=F(1,10**9);obs,bound=geometry_bound(E,src,y);sol=solve_three(E,obs)
                assert sol['feasible'];check_witness(E,obs,sol)
                cases.append({'tree':name,'edges':E,'sources':src,'pendant':y,'observed':obs,'bound':bound,'solution':sol})
    previous=json.loads((ROOT/'results/stage23_two_source_envelope.json').read_text())
    # The prior positive controls' key is obtained explicitly below.
    low=[]
    for c in previous['cases']:
        E=[(u,v,F(l)) for u,v,l in c['edges']];obs={int(r):F(d) for r,d in c['observed'].items()}
        sol=solve_three(E,obs);assert sol['feasible'] and sol['support']==[s['edge'] for s in c['sources']]
        check_witness(E,obs,sol);low.append(sol['minimum_sources_within_cap'])
    higher=[]
    for name,E in trees.items():
        for k in range(1,min(len(E),6)+1):
            for rep in range(4):
                sel=rng.sample(range(len(E)),k);ww=[rng.randint(1,9) for _ in sel]
                src=[{'edge':e,'p':F(w,sum(ww)),'x':F(rng.randint(1,99),100)} for e,w in zip(sel,ww)]
                obs,bound=geometry_bound(E,src,F(1,11))
                higher.append({'tree':name,'edges':E,'sources':src,'pendant':F(1,11),'observed':obs,'bound':bound})
    oracle=[]
    for name in ('star3','star4','quartet'):
        E=trees[name];V,a,R,D=geom(E)
        specs=[]
        for c in [c for c in cases if c['tree']==name][::max(1,len([c for c in cases if c['tree']==name])//4)]:specs.append(c['observed'])
        specs.extend([{r:F(3) for r in R},{r:F(1) for r in R}])
        for rep in range(8):specs.append({r:F(rng.randint(1,18),3) for r in R})
        for obs in specs:
            sol=solve_three(E,obs);ans=full_endpoint_oracle(E,obs)
            assert sol['feasible']==(ans is not None),(name,obs,sol,ans)
            if ans is not None:assert sol['minimum_sources_within_cap']==ans;check_witness(E,obs,sol)
            oracle.append({'tree':name,'edges':E,'observed':obs,'exact_oracle_min_up_to_three':ans,'solution':sol})
    # Equal unit-star data genuinely need all arms in the strict model.
    E=[(i,4,F(1)) for i in range(4)];obs={i:F(2) for i in range(4)}
    sol=solve_three(E,obs);assert not sol['feasible'] and full_endpoint_oracle(E,obs) is None
    oracle.append({'tree':'unit_star4_needs4','edges':E,'observed':obs,'exact_oracle_min_up_to_three':None,'solution':sol})
    result={'created_utc':datetime.now(timezone.utc).isoformat(),'passed':True,'script_sha256':sha(Path(__file__)),
       'dependencies':{p:sha(ROOT/p) for p in ['scripts/stage23_two_source_envelope.py','scripts/stage17_strict_bounds.py','results/stage23_two_source_envelope.json']},
       'three_source_cases':cases,'three_source_case_count':len(cases),'previous_low_source_cases':len(low),
       'higher_source_geometry_cases':higher,'higher_source_geometry_count':len(higher),
       'oracle_comparisons':oracle,'oracle_comparison_count':len(oracle),
       'three_source_output_counts':dict(Counter(c['solution']['minimum_sources_within_cap'] for c in cases)),
       'shifted_witnesses':sum(c['solution']['deficit']>0 for c in cases),
       'elapsed_seconds':time.perf_counter()-start,
       'scope':'Finite exact controls for a proved candidate structural characterization; no prior clearance or empirical speed comparison.'}
    (ROOT/'results/stage24_three_source.json').write_text(json.dumps(result,indent=2,default=str)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['three_source_cases','higher_source_geometry_cases','oracle_comparisons']},indent=2))

if __name__=='__main__':main()
