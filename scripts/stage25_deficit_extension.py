#!/usr/bin/env python3
"""Complete exact <=4 strict-source solver by local deficit polygons.

Rational arithmetic; no optimizer or dense distance matrix in solve_four.
The finite harness independently uses original distances and an older exact
endpoint-polytope oracle. Stage01--24 files are read-only dependencies.
"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
from datetime import datetime, timezone
import hashlib,json,time,random
from stage23_two_source_envelope import two_pass_envelope,inverse
ROOT=Path(__file__).resolve().parents[1]

def tree(edges):
    adj={}; lookup={}
    for i,(u,v,ll) in enumerate(edges):
        l=F(ll); assert l>0 and u!=v and (min(u,v),max(u,v)) not in lookup
        adj.setdefault(u,{})[v]=l;adj.setdefault(v,{})[u]=l
        lookup[min(u,v),max(u,v)]=i
    assert len(edges)==len(adj)-1 and all(len(ns)!=2 for ns in adj.values())
    seen={next(iter(adj))}; order=list(seen)
    for v in order:
        for u in adj[v]:
            if u not in seen:seen.add(u);order.append(u)
    assert len(seen)==len(adj)
    return adj,lookup

def polygon_representatives(constraints,dim):
    """One relative-interior representative of every face, possibly extras.

    Q={x: a.x<=b}; bounded by nonnegative deficits and nonnegative pendant.
    All intersections of dim independent active rows yield the vertices.
    In dimension two, vertices, all pair midpoints, and their global average
    cover every face, including when Q is a segment or a singleton.
    """
    vertices=set()
    for active in combinations(constraints,dim):
        if dim==1:
            (a,b),=active
            if not a[0]:continue
            x=(b/a[0],)
        else:
            (a,b),(c,d)=active; det=a[0]*c[1]-a[1]*c[0]
            if not det:continue
            x=((b*c[1]-a[1]*d)/det,(a[0]*d-b*c[0])/det)
        if all(sum(ai*xi for ai,xi in zip(a,x))<=b for a,b in constraints):vertices.add(x)
    vv=sorted(vertices)
    if not vv:return [],0
    points=set(vv)
    for a,b in combinations(vv,2):points.add(tuple((x+y)/2 for x,y in zip(a,b)))
    points.add(tuple(sum(v[j] for v in vv)/len(vv) for j in range(dim)))
    return sorted(points),len(vv)

def small_cover(positive,lookup,edges,cap=4):
    """Only the <=8 positive vertices are inspected, never all tree edges."""
    if not positive or len(positive)>2*cap:return None
    possible=[]
    for u,v in combinations(sorted(positive),2):
        i=lookup.get((u,v))
        if i is not None:possible.append(i)
    for size in range(1,min(cap,len(possible))+1):
        for selected in combinations(possible,size):
            if {v for i in selected for v in edges[i][:2]}==set(positive):return list(selected)
    return None

def solve_four(edges,observed):
    adj,lookup=tree(edges); leaves={v for v in adj if len(adj[v])==1}
    assert set(observed)==leaves
    obs={v:F(d) for v,d in observed.items()}
    h,_=two_pass_envelope(adj,obs)
    if any(h[v]!=obs[v] for v in leaves):return {'feasible':False,'reason':'inconsistent_leaf_envelope'}
    baseline,Y=inverse(adj,h)
    if Y<=0:return {'feasible':False,'reason':'nonpositive_baseline_pendant'}
    positive={v for v,b in baseline.items() if b>0};negative={v for v,b in baseline.items() if b<0}
    # A two-degree-three deficit set affects at most eight endpoint masses.
    if len(positive)>16 or len(negative)>8:
        return {'feasible':False,'reason':'too_many_unchangeable_baseline_masses',
                'baseline_positive_count':len(positive),'baseline_negative_count':len(negative)}
    best=None;stats={'deficit_sets_considered':0,'local_polytopes_constructed':0,
                     'polytope_vertices':0,'candidate_points':0,'strict_physical_points':0}
    def inspect(S,coeff,x):
        nonlocal best
        stats['candidate_points']+=1
        if any(t<=0 for t in x):return
        y=Y-sum(F(len(adj[v])-2,2)*t for v,t in zip(S,x))
        if y<=0:return
        affected=set(coeff)
        m={v:baseline[v] for v in positive-affected}
        for v,row in coeff.items():
            mass=baseline[v]+sum(a*t for a,t in zip(row,x))
            if mass<0:return
            if mass>0:m[v]=mass
        assert not (negative-affected)
        assert sum(m.values())==1
        stats['strict_physical_points']+=1
        selected=small_cover(set(m),lookup,edges)
        if selected is None:return
        if best is None or len(selected)<len(best['support']):
            best={'support':selected,'positive_mass':m,'pendant':y,'deficit_vertices':list(S),'deficits':list(x)}
    if not negative:inspect((),{},())
    triples=[v for v in adj if len(adj[v])==3]
    sets=((v,) for v in adj if 3<=len(adj[v])<=4)
    # Iterators avoid storing quadratically many candidate sets.
    def all_sets():
        yield from sets
        yield from combinations(triples,2)
    for S in all_sets():
        stats['deficit_sets_considered']+=1
        affected=set(S)
        for v in S:affected.update(adj[v])
        if negative-affected or len(positive-affected)>8:continue
        dim=len(S); coeff={v:[F(0)]*dim for v in affected}
        for j,v in enumerate(S):
            coeff[v][j]=-sum(1/l for l in adj[v].values())/2
            for u,l in adj[v].items():coeff[u][j]+=1/(2*l)
        constraints=[]
        for j in range(dim):constraints.append((tuple(-F(i==j) for i in range(dim)),F(0)))
        constraints.append((tuple(F(len(adj[v])-2,2) for v in S),Y))
        for v in sorted(affected):constraints.append((tuple(-a for a in coeff[v]),baseline[v]))
        stats['local_polytopes_constructed']+=1
        reps,nv=polygon_representatives(constraints,dim);stats['polytope_vertices']+=nv
        for x in reps:inspect(S,coeff,x)
    if best is None:return {'feasible':False,'reason':'no_at_most_four_strict_representation','statistics':stats}
    degree={v:sum(v in edges[i][:2] for i in best['support']) for v in best['positive_mass']}
    sources=[]
    for i in best['support']:
        u,v,l=edges[i];a=best['positive_mass'][u]/degree[u];b=best['positive_mass'][v]/degree[v]
        sources.append({'edge':i,'p':a+b,'x':b/(a+b)})
    f=h.copy()
    for v,t in zip(best['deficit_vertices'],best['deficits']):f[v]+=t
    best.update({'feasible':True,'minimum_sources_within_cap':len(sources),'sources':sources,
                 'potential':f,'envelope':h,'statistics':stats})
    return best

def direct_geometry(edges):
    adj,_=tree(edges);D={}
    for r in adj:
        D[r]={r:F(0)};order=[r]
        for v in order:
            for u,l in adj[v].items():
                if u not in D[r]:D[r][u]=D[r][v]+l;order.append(u)
    return adj,D,[v for v in adj if len(adj[v])==1]

def direct_observe(edges,sources,y):
    adj,D,R=direct_geometry(edges)
    assert y>0 and sum(s['p'] for s in sources)==1
    f={r:y for r in adj};m={v:F(0) for v in adj}
    for s in sources:
        u,v,l=edges[s['edge']];x,p=s['x'],s['p'];assert 0<x<1 and p>0
        m[u]+=p*(1-x);m[v]+=p*x
        for r in adj:f[r]+=p*min(D[r][u]+x*F(l),D[r][v]+(1-x)*F(l))
    return {r:f[r] for r in R},f,m

def witness_check(edges,obs,result):
    if not result['feasible']:return
    obs2,f,m=direct_observe(edges,result['sources'],result['pendant'])
    assert obs2==obs and f==result['potential']
    assert {v:a for v,a in m.items() if a>0}==result['positive_mass']
    assert len(result['sources'])==result['minimum_sources_within_cap']<=4

def full_endpoint_oracle(edges,obs,cap=4):
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
    start=time.perf_counter(); rng=random.Random(2513)
    star4=[(i,4,F(1)) for i in range(4)]
    star5=[(i,5,F(i+1,3)) for i in range(5)]
    cat5=[(0,5,F(1)),(1,5,F(1)),(5,6,F(1)),(2,6,F(1)),(6,7,F(1)),(3,7,F(1)),(4,7,F(1))]
    weighted=[(u,v,F(i+1,(i%3)+1)) for i,(u,v,l) in enumerate(cat5)]
    cases=[]
    for name,E in [('star4',star4),('weighted_star5',star5),('caterpillar5',cat5),('weighted_caterpillar5',weighted)]:
        sels=list(combinations(range(len(E)),4))
        if len(sels)>10:sels=[sels[i] for i in sorted(rng.sample(range(len(sels)),10))]
        for selected in sels:
            for rep in (0,1):
                pp=[F(1,4)]*4 if rep==0 else [F(1,10),F(1,5),F(3,10),F(2,5)]
                xx=[F(1,2)]*4 if rep==0 else [F(1,1000),F(2,7),F(3,5),F(999,1000)]
                sources=[{'edge':e,'p':p,'x':x} for e,p,x in zip(selected,pp,xx)]
                obs,f,m=direct_observe(E,sources,F(1,10**9)); sol=solve_four(E,obs)
                assert sol['feasible'];witness_check(E,obs,sol)
                cases.append({'tree':name,'edges':E,'sources':sources,'pendant':F(1,10**9),'observed':obs,'solution':sol})
    special=[]
    def special_case(name,E,obs,expected):
        sol=solve_four(E,obs);oracle=full_endpoint_oracle(E,obs)
        assert oracle==expected and sol['feasible']==(expected is not None)
        if expected is not None:assert sol['minimum_sources_within_cap']==expected;witness_check(E,obs,sol)
        special.append({'name':name,'edges':E,'observed':obs,'exact_oracle_min_up_to_four':oracle,'solution':sol})
    special_case('degree_four_deficit_is_essential',star4,{i:F(2) for i in range(4)},4)
    special_case('two_nonadjacent_degree_three_deficits_are_essential',cat5,{0:F(13,4),1:F(13,4),2:F(7,2),3:F(13,4),4:F(13,4)},4)
    unitstar5=[(i,5,F(1)) for i in range(5)]
    special_case('five_arm_star_requires_more_than_four',unitstar5,{i:F(2) for i in range(5)},None)
    # A valid closed solution at y=0 need not admit a positive-pendant solution.
    special_case('zero_pendant_only',star4,{i:F(1) for i in range(4)},None)
    special_case('invalid_leaf_lipschitz',star4,{0:F(100),1:F(1),2:F(1),3:F(1)},None)
    for c in [c for c in cases if c['tree'] in ('caterpillar5','weighted_caterpillar5')][::6]:
        expected=full_endpoint_oracle(c['edges'],c['observed'])
        assert expected==c['solution']['minimum_sources_within_cap']
        special.append({'name':'generated_exact_oracle','edges':c['edges'],'observed':c['observed'],
                        'exact_oracle_min_up_to_four':expected,'solution':c['solution']})
    # Lower-dimensional polytope controls independent of the tree application.
    geometry=[]
    systems=[('singleton',[( (F(1),F(0)),F(1)),((-F(1),F(0)),-F(1)),((F(0),F(1)),F(2)),((F(0),-F(1)),-F(2))],{(F(1),F(2))}),
      ('segment',[((F(1),F(0)),F(1)),((-F(1),F(0)),F(0)),((F(0),F(1)),F(0)),((F(0),-F(1)),F(0))],{(F(0),F(0)),(F(1,2),F(0)),(F(1),F(0))}),
      ('triangle',[((-F(1),F(0)),F(0)),((F(0),-F(1)),F(0)),((F(1),F(1)),F(1))],{(F(0),F(0)),(F(1),F(0)),(F(0),F(1)),(F(1,2),F(0)),(F(0),F(1,2)),(F(1,2),F(1,2)),(F(1,3),F(1,3))})]
    for name,cs,required in systems:
        reps,nv=polygon_representatives(cs,2);assert required<=set(reps)
        geometry.append({'name':name,'constraints':cs,'representatives':reps,'vertex_count':nv})
    old=json.loads((ROOT/'results/stage24_three_source.json').read_text());previous=[]
    for c in old['oracle_comparisons']:
        E=[(u,v,F(l)) for u,v,l in c['edges']];obs={int(v):F(d) for v,d in c['observed'].items()}
        sol=solve_four(E,obs);expected=full_endpoint_oracle(E,obs)
        assert sol['feasible']==(expected is not None)
        if expected is not None:assert sol['minimum_sources_within_cap']==expected;witness_check(E,obs,sol)
        previous.append({'edges':E,'observed':obs,'exact_oracle_min_up_to_four':expected,'solution':sol})
    result={'created_utc':datetime.now(timezone.utc).isoformat(),'passed':True,'script_sha256':sha(Path(__file__)),
      'dependencies':{p:sha(ROOT/p) for p in ['scripts/stage23_two_source_envelope.py','scripts/stage17_strict_bounds.py','results/stage24_three_source.json']},
      'generated_case_count':len(cases),'generated_cases':cases,'special_oracle_count':len(special),'special_oracle_cases':special,
      'previous_oracle_count':len(previous),'previous_oracle_cases':previous,'polytope_geometry_controls':geometry,
      'elapsed_seconds':time.perf_counter()-start,
      'scope':'Exact <=4-source characterization controls, not novelty clearance or speed superiority. Old full-polytope oracle is only used in the finite harness.'}
    (ROOT/'results/stage25_deficit_extension.json').write_text(json.dumps(result,indent=2,default=str)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['generated_cases','special_oracle_cases','previous_oracle_cases','polytope_geometry_controls']},indent=2))

if __name__=='__main__':main()
