#!/usr/bin/env python3
"""General-K exact PDD anchored-growth inference with explicit resource unknowns.

All-leaf, positive reduced tree, strict original-edge sources, free positive
mean pendant, rational exact observations. Main writes new Stage27 only.
"""
import os
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[_key]='1'
from fractions import Fraction as F
from itertools import combinations
from collections import deque
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,time
from stage23_two_source_envelope import two_pass_envelope,inverse
from stage25_deficit_extension import tree
ROOT=Path(__file__).resolve().parents[1]


class ResourceLimit(Exception):pass


class Budget:
    def __init__(self,seconds=None,max_states=None,max_bases=None,max_faces=None):
        self.started=time.perf_counter();self.deadline=None if seconds is None else self.started+seconds
        self.limits={'growth_states':max_states,'vertex_bases':max_bases,'mass_faces':max_faces}
        self.counts={key:0 for key in self.limits}
    def check(self):
        if self.deadline is not None and time.perf_counter()>=self.deadline:raise ResourceLimit('wall_time_limit')
    def touch(self,key):
        self.check();self.counts[key]+=1
        limit=self.limits.get(key)
        if limit is not None and self.counts[key]>limit:raise ResourceLimit(key+'_limit')


def dot(a,b):return sum((x*y for x,y in zip(a,b)),F(0))


def square_solve(A,b):
    n=len(b);M=[[F(x) for x in row]+[F(v)] for row,v in zip(A,b)]
    for j in range(n):
        p=next((i for i in range(j,n) if M[i][j]),None)
        if p is None:return None
        M[j],M[p]=M[p],M[j];c=M[j][j];M[j]=[x/c for x in M[j]]
        for i in range(n):
            if i!=j and M[i][j]:
                c=M[i][j];M[i]=[x-c*y for x,y in zip(M[i],M[j])]
    return tuple(row[-1] for row in M)


def polytope_vertices(constraints,dim,budget=None):
    """All vertices of bounded {x:a*x<=b}; lower-dimensional cases included."""
    budget=budget or Budget();vertices=set()
    for I in combinations(range(len(constraints)),dim):
        budget.touch('vertex_bases')
        x=square_solve([constraints[i][0] for i in I],[constraints[i][1] for i in I])
        if x is not None and all(dot(a,x)<=b for a,b in constraints):vertices.add(x)
    return sorted(vertices)


def mass_face_representatives(vertices,mass_rows,budget=None):
    """Average all vertices of every intersection of mass-zero faces.

    mass_rows=(constant, coefficients), nonnegative throughout the polytope.
    Vertex-index bitsets deduplicate equal faces. This includes the whole Q.
    """
    if not vertices:return
    budget=budget or Budget();dim=len(vertices[0]);allbits=(1<<len(vertices))-1
    zero_sets=[]
    for constant,row in mass_rows:
        bits=0
        for i,x in enumerate(vertices):
            value=constant+dot(row,x);assert value>=0
            if value==0:bits|=1<<i
        if bits and bits!=allbits:zero_sets.append(bits)
    zero_sets=sorted(set(zero_sets));seen={allbits};todo=[allbits]
    while todo:
        budget.touch('mass_faces');bits=todo.pop();selected=[];rest=bits
        while rest:
            low=rest & -rest;selected.append(low.bit_length()-1);rest-=low
        point=tuple(sum((vertices[i][j] for i in selected),F(0))/len(selected) for j in range(dim))
        yield point
        for z in zero_sets:
            intersection=bits & z
            if intersection and intersection not in seen:seen.add(intersection);todo.append(intersection)


def minimum_forest_cover(positive,lookup,edges,cap):
    """Maximum forest matching by leaf removal, followed by unmatched edges."""
    P=set(positive)
    if not P or len(P)>2*cap:return None
    adj={v:set() for v in P}
    for u,v in combinations(sorted(P),2):
        if (u,v) in lookup:adj[u].add(v);adj[v].add(u)
    if any(not adj[v] for v in P):return None
    active=set(P);degree={v:len(adj[v]) for v in P};queue=deque(sorted(v for v in P if degree[v]==1));matching=[]
    while queue:
        u=queue.popleft()
        if u not in active or degree[u]!=1:continue
        v=next(v for v in adj[u] if v in active);matching.append(lookup[min(u,v),max(u,v)])
        for a in (u,v):
            active.remove(a)
        for a in (u,v):
            for w in adj[a]:
                if w in active:
                    degree[w]-=1
                    if degree[w]==1:queue.append(w)
    covered={v for i in matching for v in edges[i][:2]};selected=set(matching)
    for u in sorted(P-covered):
        v=min(adj[u]);selected.add(lookup[min(u,v),max(u,v)])
    assert len(selected)==len(P)-len(matching)
    assert {v for i in selected for v in edges[i][:2]}==P
    return sorted(selected) if len(selected)<=cap else None


def components(S,adj):
    left=set(S);answer=[]
    while left:
        seed=min(left);left.remove(seed);C={seed};todo=[seed]
        for v in todo:
            for u in adj[v]:
                if u in left:left.remove(u);C.add(u);todo.append(u)
        answer.append(C)
    return answer


def solve_general(edges,observed,cap=3,time_limit_seconds=30.0,max_states=200000,
                  max_bases=2000000,max_faces=500000,record_feasible_sets=False):
    """Complete exact result or explicit unknown; no optimizer is used.

    On a resource stop an incumbent certifies feasibility but NOT its minimality.
    No full-tree scan occurs per grown candidate; all local sets are O(K).
    """
    if not isinstance(cap,int) or isinstance(cap,bool) or cap<1:raise ValueError('cap must be a positive integer')
    budget=Budget(time_limit_seconds,max_states,max_bases,max_faces)
    stats={'cap':cap,'candidate_vertex_polytopes':0,'empty_polytopes':0,'local_mass_entries':0,
           'flux_rejections':0,'unchangeable_mass_rejections':0,'strict_face_points':0,
           'states_with_nonpositive_selected_baseline':0,'disconnected_growth_states':0,
           'maximum_deficit_dimension':0,'maximum_local_vertex_count':0,'seed_masks':0,'frontier_states':0}
    best=None;min_deficits=None;feasible_sets=set();adj=None;h=None
    def summarize():return dict(stats,**budget.counts,elapsed_seconds=time.perf_counter()-budget.started)
    def materialize(raw,minimum=False):
        if raw is None:return None
        mass=raw['positive_mass'];S=raw['support'];degree={v:sum(v in edges[i][:2] for i in S) for v in mass};sources=[]
        for i in S:
            u,v,l=edges[i];a=mass[u]/degree[u];b=mass[v]/degree[v]
            sources.append({'edge':i,'p':a+b,'x':b/(a+b)})
        potential=h.copy()
        for v,t in zip(raw['deficit_vertices'],raw['deficits']):potential[v]+=t
        result=dict(raw,sources=sources,potential=potential)
        result['minimum_sources_within_cap' if minimum else 'source_count']=len(S)
        return result
    def complete_reject(reason):return {'status':'complete','feasible':False,'reason':reason,'statistics':summarize()}
    try:
        budget.check();adj,lookup=tree(edges);leaves={v for v in adj if len(adj[v])==1}
        if set(observed)!=leaves:raise ValueError('Requires every original leaf and no extra observations')
        obs={v:F(d) for v,d in observed.items()};h,_=two_pass_envelope(adj,obs)
        if any(h[r]!=obs[r] for r in leaves):return complete_reject('inconsistent_leaf_envelope')
        b,Y=inverse(adj,h);stats['vertex_count']=len(adj)
        if Y<=0:return complete_reject('nonpositive_baseline_pendant')
        positive={v for v,x in b.items() if x>0};negative={v for v,x in b.items() if x<0}
        stats.update(baseline_positive_count=len(positive),baseline_negative_count=len(negative))
        posbound=2*cap if cap<=2 else 3*cap-2;negbound=0 if cap<=2 else 4*cap-8
        if len(positive)>posbound or len(negative)>negbound:return complete_reject('baseline_mass_count_obstruction')
        # Eligibility depends on degree, NOT on baseline sign.
        eligible={v for v in adj if 3<=len(adj[v])<=cap};anchors=positive & eligible
        stats['positive_eligible_anchors']=len(anchors);maxweight=max(cap-2,0)
        def inspect(S,coeff,point):
            nonlocal best,min_deficits
            budget.check()
            if any(t<=0 for t in point):return
            y=Y-sum(F(len(adj[v])-2,2)*t for v,t in zip(S,point))
            if y<=0:return
            affected=set(coeff);mass={v:b[v] for v in positive-affected}
            for v,row in coeff.items():
                stats['local_mass_entries']+=1;a=b[v]+dot(row,point)
                if a<0:return
                if a>0:mass[v]=a
            assert not (negative-affected)
            if len(mass)>2*cap:return
            assert sum(mass.values())==1
            stats['strict_face_points']+=1;selected=minimum_forest_cover(set(mass),lookup,edges,cap)
            if selected is None:return
            min_deficits=len(S) if min_deficits is None else min(min_deficits,len(S))
            if record_feasible_sets:feasible_sets.add(tuple(S))
            if best is None or len(selected)<len(best['support']):
                best={'support':selected,'positive_mass':mass,'pendant':y,'deficit_vertices':list(S),'deficits':list(point)}
        def evaluate(S):
            stats['maximum_deficit_dimension']=max(stats['maximum_deficit_dimension'],len(S))
            if not S:
                if not negative:inspect((),{},())
                return
            CC=components(S,adj)
            if len(CC)>1:stats['disconnected_growth_states']+=1
            if any(b[v]<=0 for v in S):stats['states_with_nonpositive_selected_baseline']+=1
            if any(sum((b[v] for v in C),F(0))<=0 for C in CC):
                stats['flux_rejections']+=1;return
            affected=set(S)
            for v in S:affected.update(adj[v])
            if negative-affected or len(positive-affected)>2*cap:
                stats['unchangeable_mass_rejections']+=1;return
            coeff={v:[F(0)]*len(S) for v in affected}
            for j,v in enumerate(S):
                coeff[v][j]=-sum(1/l for l in adj[v].values())/2
                for u,l in adj[v].items():coeff[u][j]+=1/(2*l)
            constraints=[]
            for j in range(len(S)):constraints.append((tuple(-F(i==j) for i in range(len(S))),F(0)))
            constraints.append((tuple(F(len(adj[v])-2,2) for v in S),Y))
            for v in sorted(affected):constraints.append((tuple(-x for x in coeff[v]),b[v]))
            stats['candidate_vertex_polytopes']+=1
            vertices=polytope_vertices(constraints,len(S),budget)
            stats['maximum_local_vertex_count']=max(stats['maximum_local_vertex_count'],len(vertices))
            if not vertices:stats['empty_polytopes']+=1;return
            rows=[(b[v],coeff[v]) for v in sorted(affected)]
            for point in mass_face_representatives(vertices,rows,budget):inspect(S,coeff,point)
        # Canonical seed mask R=S intersect anchors. Every admissible smaller
        # prefix appears as its own terminal branch, including the empty set.
        anchor_list=sorted(anchors)
        for nr in range(min(len(anchor_list),maxweight)+1):
            for Rtuple in combinations(anchor_list,nr):
                budget.check();R=frozenset(Rtuple);weight=sum(len(adj[v])-2 for v in R)
                if weight>maxweight:continue
                stats['seed_masks']+=1
                rejected=frozenset(anchors-R);frontier=set()
                for v in R:frontier.update(u for u in adj[v] if u in eligible and u not in R and u not in rejected)
                stack=[(R,rejected,frozenset(frontier),weight)]
                while stack:
                    budget.check();stats['frontier_states']+=1
                    current,rejected,frontier,weight=stack.pop()
                    if not frontier:
                        budget.touch('growth_states');evaluate(tuple(sorted(current)));continue
                    v=min(frontier);remaining=frontier-{v}
                    stack.append((current,rejected|{v},remaining,weight))
                    newweight=weight+len(adj[v])-2
                    if newweight<=maxweight:
                        included=current|{v};expanded=set(remaining)
                        expanded.update(u for u in adj[v] if u in eligible and u not in included and u not in rejected)
                        stack.append((included,rejected,frozenset(expanded),newweight))
        if best is None:return complete_reject('no_strict_representation_within_budget')
        answer=materialize(best,True)
        return dict(answer,status='complete',feasible=True,statistics=summarize(),
                    minimum_positive_deficit_vertices_within_cap=min_deficits,
                    feasible_deficit_sets=[list(S) for S in sorted(feasible_sets)] if record_feasible_sets else None)
    except ResourceLimit as exc:
        return {'status':'unknown','feasibility_established':best is not None,'incumbent':materialize(best,False),
                'reason':str(exc),'statistics':summarize(),
                'scope':'Resource truncation does not prove infeasibility or minimum source count'}


def validate_original_paths(edges,observed,solution):
    if solution['status']=='unknown':solution=solution.get('incumbent')
    elif not solution['feasible']:return
    if solution is None:return
    adj={v:{} for e in edges for v in e[:2]}
    for u,v,l in edges:adj[u][v]=adj[v][u]=F(l)
    D={};src=solution['sources'];y=F(solution['pendant'])
    assert y>0 and src and sum(F(s['p']) for s in src)==1
    assert len({s['edge'] for s in src})==len(src)
    for s in src:
        assert F(s['p'])>0 and 0<F(s['x'])<1
        for root in edges[s['edge']][:2]:
            if root in D:continue
            ds={root:F(0)};todo=[root]
            for v in todo:
                for u,l in adj[v].items():
                    if u not in ds:ds[u]=ds[v]+l;todo.append(u)
            D[root]=ds
    for r,d in observed.items():
        q=y
        for s in src:
            u,v,l=edges[s['edge']];p,x=F(s['p']),F(s['x']);l=F(l)
            q+=p*min(D[u][r]+x*l,D[v][r]+(1-x)*l)
        assert q==F(d)


def alternating_fixture(cap):
    """K midpoint source leaves, K-3 empty leaves; alternating deficits."""
    s=cap-2;inner=list(range(100,100+2*s-1));pairs=[];source_edge_ids=[];leaf=0
    for j,v in enumerate(inner):
        count=2 if j in (0,len(inner)-1) else 1
        for _ in range(count):
            if j%2==0:source_edge_ids.append(len(pairs))
            pairs.append((v,leaf,F(1)));leaf+=1
    pairs.extend((u,v,F(1)) for u,v in zip(inner,inner[1:]))
    sources=[{'edge':i,'p':F(1,cap),'x':F(1,2)} for i in source_edge_ids]
    assert len(sources)==cap
    from stage26_scale import observations
    observed=observations(pairs,sources,F(1))
    return {'name':'alternating_'+str(cap),'edges':pairs,'observed':observed,'generating_sources':sources,
            'generating_pendant':F(1),'known_generating_deficit_vertices':inner[::2],'known_generating_deficits':[F(1,cap)]*s}


def serial(value):
    if isinstance(value,F):return str(value)
    if isinstance(value,dict):return {str(k):serial(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [serial(v) for v in value]
    return value


def decode_case(c):
    return [(u,v,F(l)) for u,v,l in c['edges']],{int(v):F(d) for v,d in c['observed'].items()}


def geometry_controls():
    answer=[]
    # A lower-dimensional triangle embedded in R^3, with seven mass faces.
    Q=[(tuple(-F(i==j) for i in range(3)),F(0)) for j in range(3)]
    Q += [(tuple([F(1)]*3),F(1)),(tuple([F(-1)]*3),F(-1))]
    rows=[(F(0),tuple(F(i==j) for i in range(3))) for j in range(3)]
    V=polytope_vertices(Q,3);P=list(mass_face_representatives(V,rows))
    assert len(V)==3 and len(P)==7 and all(sum(x)==1 for x in P)
    answer.append({'name':'embedded_triangle_dimension3','constraints':Q,'mass_rows':rows,'vertices':V,'face_representatives':P})
    # A segment embedded in R^4: x0=x1, x2=x3=0.
    Q=[((F(-1),0,0,0),F(0)),((F(1),0,0,0),F(1)),((1,-1,0,0),F(0)),((-1,1,0,0),F(0)),
       ((0,0,1,0),F(0)),((0,0,-1,0),F(0)),((0,0,0,1),F(0)),((0,0,0,-1),F(0))]
    rows=[(F(0),(F(1),0,0,0)),(F(1),(F(-1),0,0,0))]
    V=polytope_vertices(Q,4);P=list(mass_face_representatives(V,rows))
    assert len(V)==2 and len(P)==3 and all(x[0]==x[1] and x[2]==x[3]==0 for x in P)
    answer.append({'name':'embedded_segment_dimension4','constraints':Q,'mass_rows':rows,'vertices':V,'face_representatives':P})
    Q=[]
    for j in range(4):
        a=tuple(F(i==j) for i in range(4));Q += [(a,F(1,4)),(tuple(-x for x in a),F(-1,4))]
    V=polytope_vertices(Q,4);assert V==[(F(1,4),)*4]
    answer.append({'name':'singleton_dimension4','constraints':Q,'vertices':V})
    return answer


def main():
    started=time.process_time();deadline=time.perf_counter()+145;records=[];regression=[];scales=[]
    oracle_path=ROOT/'results/stage27_general_oracle.json';oracle=json.loads(oracle_path.read_text())
    for c in oracle['cases']:
        E,d=decode_case(c);remaining=deadline-time.perf_counter()
        r=solve_general(E,d,c['cap'],time_limit_seconds=max(0,min(15,remaining)),record_feasible_sets=True)
        validate_original_paths(E,d,r);expected=c['oracle']['minimum']
        agrees=None if r['status']=='unknown' else (r.get('minimum_sources_within_cap') if r['feasible'] else None)==expected
        assert agrees is not False,(c['name'],expected,r)
        records.append({'name':c['name'],'group':c['group'],'edges':E,'observed':d,'cap':c['cap'],
                        'expected_minimum_from_independent_original_support_oracle':expected,'solution':r,
                        'agrees_with_oracle':agrees,'original_path_witness_checked':bool(r.get('feasible') or r.get('feasibility_established'))})
    # Existing exact regression inputs: no old writer or solver main is called.
    old=json.loads((ROOT/'results/stage26_adversarial.json').read_text())
    old2=json.loads((ROOT/'results/stage25_deficit_extension.json').read_text())
    cases=[('stage26_adversarial',i,c,c['oracle_minimum_up_to_four']) for i,c in enumerate(old['cases'])]
    for key in ('generated_cases','special_oracle_cases','previous_oracle_cases'):
        for i,c in enumerate(old2[key]):
            expected=c.get('exact_oracle_min_up_to_four',c['solution'].get('minimum_sources_within_cap'))
            cases.append(('stage25_deficit_extension/'+key,i,c,expected))
    for origin,index,c,expected in cases:
        E,d=decode_case(c);r=solve_general(E,d,4,time_limit_seconds=max(0,min(5,deadline-time.perf_counter())))
        validate_original_paths(E,d,r)
        agrees=None if r['status']=='unknown' else (r.get('minimum_sources_within_cap') if r['feasible'] else None)==expected
        assert agrees is not False,(origin,index,expected,r)
        regression.append({'origin':origin,'index':index,'expected_minimum':expected,'solution':r,'agrees':agrees})
    resource=[]
    c=alternating_fixture(5)
    for kwargs in ({'time_limit_seconds':0},{'max_states':0},{'max_bases':0},{'max_faces':0}):
        r=solve_general(c['edges'],c['observed'],5,**kwargs);assert r['status']=='unknown'
        assert 'minimum_sources_within_cap' not in r;validate_original_paths(c['edges'],c['observed'],r)
        resource.append({'name':'hard5_limit_'+next(iter(kwargs)),'edges':c['edges'],'observed':c['observed'],'cap':5,'kwargs':kwargs,'solution':r})
    from stage26_scale import caterpillar,observations,center_profile
    E=caterpillar(4,False);src=[{'edge':0,'p':F(1),'x':F(1,2)}];d=observations(E,src,F(1))
    r=solve_general(E,d,5,max_bases=0);assert r['status']=='unknown' and r['feasibility_established']
    assert 'minimum_sources_within_cap' not in r['incumbent'];validate_original_paths(E,d,r)
    resource.append({'name':'resource_limit_keeps_only_feasibility_incumbent','edges':E,'observed':d,'cap':5,'kwargs':{'max_bases':0},'solution':r})
    for n,k in ((100,5),(1000,5),(1000,6)):
        c=center_profile(n,'four');E,d=c['edges'],c['observed'];runs=[]
        for repeat in range(3):
            r=solve_general(E,d,k,time_limit_seconds=max(0,min(12,deadline-time.perf_counter())))
            validate_original_paths(E,d,r)
            if r['status']=='complete':assert r['feasible'] and r['minimum_sources_within_cap']==3
            runs.append(r)
        scales.append({'name':'balanced_center_'+str(n)+'_cap'+str(k),'edges':E,'observed':d,'cap':k,'expected_minimum':3,'runs':runs})
    geometry=geometry_controls()
    dependencies={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in (
      'scripts/stage23_two_source_envelope.py','scripts/stage25_deficit_extension.py','scripts/stage26_scale.py',
      'review/stage26_anchor_localization.md','review/stage27_parameter_bound.md',
      'results/stage26_adversarial.json','results/stage25_deficit_extension.json')}
    allsol=[c['solution'] for c in records]+[c['solution'] for c in regression]+[r for c in scales for r in c['runs']]
    out={'created_utc':datetime.now(timezone.utc).isoformat(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'dependencies':dependencies,'primary_cases':records,'regression':regression,'resource_controls':resource,
         'geometry_controls':geometry,'scale_cases':scales,'summary':{'primary_cases':len(records),'regression_cases':len(regression),
          'completed':sum(r['status']=='complete' for r in allsol),'unknown':sum(r['status']=='unknown' for r in allsol),
          'resource_controls':len(resource),'geometry_controls':len(geometry),'scale_calls':sum(len(c['runs']) for c in scales),
          'cpu_seconds':time.process_time()-started},'passed':True,
         'scope':'Exact rational synthetic/oracle controls, not biological observations or noisy inference. Resource limits return unknown; finite tests do not prove the theorem. No frozen experiment main executed.'}
    (ROOT/'results/stage27_general_fpt.json').write_text(json.dumps(serial(out),indent=2)+'\n')
    print(json.dumps(out['summary'],indent=2));return out


if __name__=='__main__':main()
