#!/usr/bin/env python3
"""Exact <=3 strict-source decision/minimum in linear arithmetic work.

Uses constant-size local mass checks after one global envelope/inversion pass.
Only main writes the new Stage25 result; imported helpers/writers are not run.
"""
from fractions import Fraction as F
from pathlib import Path
from datetime import datetime, timezone
import hashlib,json,time
from stage23_two_source_envelope import two_pass_envelope,inverse
from stage25_deficit_extension import tree,polygon_representatives,small_cover
ROOT=Path(__file__).resolve().parents[1]

def solve_three_linear(edges,observed):
    adj,lookup=tree(edges);leaves={v for v in adj if len(adj[v])==1}
    assert set(observed)==leaves
    obs={v:F(d) for v,d in observed.items()}
    h,_=two_pass_envelope(adj,obs)
    if any(h[v]!=obs[v] for v in leaves):return {'feasible':False,'reason':'inconsistent_leaf_envelope'}
    b,Y=inverse(adj,h)
    if Y<=0:return {'feasible':False,'reason':'nonpositive_baseline_pendant'}
    positive={v for v,m in b.items() if m>0};negative={v for v,m in b.items() if m<0}
    stats={'vertex_count':len(adj),'baseline_positive_count':len(positive),'baseline_negative_count':len(negative),
           'degree_three_vertices_considered':0,'local_polytopes_constructed':0,
           'candidate_points':0,'local_mass_entries_processed':0}
    # Any nonzero deficit changes only v plus its three neighbors; a strict
    # three-edge support has at most six positive endpoints.
    if len(positive)>10 or len(negative)>4:
        return {'feasible':False,'reason':'too_many_unchangeable_baseline_masses','statistics':stats}
    def inspect(v,coeff,t):
        stats['candidate_points']+=1
        if v is not None and t<=0:return None
        y=Y-t/2
        if y<=0:return None
        affected=set(coeff)
        m={u:b[u] for u in positive-affected}
        for u,a in coeff.items():
            stats['local_mass_entries_processed']+=1
            mass=b[u]+a*t
            if mass<0:return None
            if mass>0:m[u]=mass
        assert not (negative-affected)
        if len(m)>6:return None
        selected=small_cover(set(m),lookup,edges,cap=3)
        if selected is None:return None
        assert sum(m.values())==1
        if v is not None:assert len(selected)==3
        degree={u:sum(u in edges[i][:2] for i in selected) for u in m}
        sources=[]
        for i in selected:
            u,w,l=edges[i];a=m[u]/degree[u];c=m[w]/degree[w]
            sources.append({'edge':i,'p':a+c,'x':c/(a+c)})
        f=h.copy()
        if v is not None:f[v]+=t
        return {'feasible':True,'minimum_sources_within_cap':len(selected),'support':selected,
                'sources':sources,'pendant':y,'positive_mass':m,'potential':f,'envelope':h,
                'deficit_vertex':v,'deficit':t,'statistics':dict(stats)}
    if not negative:
        solution=inspect(None,{},F(0))
        if solution is not None:return solution
    for v in adj:
        if len(adj[v])!=3:continue
        stats['degree_three_vertices_considered']+=1
        affected={v,*adj[v]}
        if negative-affected or len(positive-affected)>6:continue
        coeff={v:-sum(1/l for l in adj[v].values())/2}
        coeff.update({u:1/(2*l) for u,l in adj[v].items()})
        constraints=[((-F(1),),F(0)),((F(1,2),),Y)]
        constraints.extend(((-a,),b[u]) for u,a in coeff.items())
        stats['local_polytopes_constructed']+=1
        reps,_=polygon_representatives(constraints,1)
        for x in reps:
            solution=inspect(v,coeff,x[0])
            if solution is not None:return solution
    return {'feasible':False,'reason':'no_at_most_three_strict_representation','statistics':stats}

def check_witness(edges,observed,sol):
    """Original-path evaluation from <=6 endpoints, not an all-pairs matrix."""
    if not sol['feasible']:return
    adj,_=tree(edges);D={};mass={v:F(0) for v in adj}
    for s in sol['sources']:
        u,v,l=edges[s['edge']];p,x=s['p'],s['x']
        assert p>0 and 0<x<1
        mass[u]+=p*(1-x);mass[v]+=p*x
        for r in (u,v):
            if r in D:continue
            D[r]={r:F(0)};todo=[r]
            for a in todo:
                for b,ll in adj[a].items():
                    if b not in D[r]:D[r][b]=D[r][a]+ll;todo.append(b)
    assert sum(mass.values())==1 and sol['pendant']>0
    f={r:sol['pendant'] for r in adj}
    for s in sol['sources']:
        u,v,l=edges[s['edge']];p,x=s['p'],s['x']
        for r in adj:f[r]+=p*min(D[u][r]+F(l)*x,D[v][r]+F(l)*(1-x))
    assert f==sol['potential'] and {r:f[r] for r in observed}==observed
    assert {v:m for v,m in mass.items() if m>0}==sol['positive_mass']
    assert len(sol['sources'])==sol['minimum_sources_within_cap']<=3
    assert [s['edge'] for s in sol['sources']]==sol['support']

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    start=time.perf_counter();old_path=ROOT/'results/stage24_three_source.json';bench_path=ROOT/'results/stage25_benchmark.json'
    old=json.loads(old_path.read_text());bench=json.loads(bench_path.read_text());cases=[]
    def run(name,kind,c,expected):
        edges=[(u,v,F(l)) for u,v,l in c['edges']];obs={int(v):F(d) for v,d in c['observed'].items()}
        solution=solve_three_linear(edges,obs)
        assert solution['feasible']==expected['feasible'],(name,solution,expected)
        if expected['feasible']:
            assert solution['minimum_sources_within_cap']==expected['minimum_sources_within_cap']
            check_witness(edges,obs,solution)
        cases.append({'name':name,'kind':kind,'edges':edges,'observed':obs,'expected_cap_three':expected,'solution':solution})
    for i,c in enumerate(old['three_source_cases']):
        run('stage24_generated_'+str(i),'frozen_exact_generated',c,
            {k:c['solution'][k] for k in ['feasible','minimum_sources_within_cap']})
    for i,c in enumerate(old['oracle_comparisons']):
        k=c['exact_oracle_min_up_to_three'];expected={'feasible':k is not None}
        if k is not None:expected['minimum_sources_within_cap']=k
        run('stage24_oracle_'+str(i),'frozen_exact_endpoint_oracle',c,expected)
    for c in bench['cases']:
        expected={'feasible':c['structural_solution']['feasible']}
        if expected['feasible']:expected['minimum_sources_within_cap']=c['structural_solution']['minimum_sources_within_cap']
        run(c['name'],c['kind'],c,expected)
    edge=[(0,1,F(2))]
    run('single_edge_strict','analytic_boundary',{'edges':edge,'observed':{0:F(11,12),1:F(19,12)}},{'feasible':True,'minimum_sources_within_cap':1})
    run('single_edge_zero_pendant','analytic_boundary',{'edges':edge,'observed':{0:F(1),1:F(1)}},{'feasible':False})
    run('single_edge_endpoint','analytic_boundary',{'edges':edge,'observed':{0:F(1),1:F(3)}},{'feasible':False})
    dependencies=['scripts/stage23_two_source_envelope.py','scripts/stage25_deficit_extension.py',
                  'results/stage24_three_source.json','results/stage25_benchmark.json']
    result={'passed':True,'created_utc':datetime.now(timezone.utc).isoformat(),'script_sha256':sha(Path(__file__)),
            'dependencies':{p:sha(ROOT/p) for p in dependencies},'case_count':len(cases),
            'case_groups':{'stage24_generated':100,'stage24_endpoint_oracle':41,'stage25_benchmark':16,'single_edge_boundaries':3},
            'public_generating_geometry_count':sum(c['kind']=='frozen_public_simulator_generating_geometry' for c in cases),
            'cases':cases,'elapsed_seconds':time.perf_counter()-start,
            'scope':'Pure solver calls on frozen input tables and analytic boundaries. No frozen writer or timed baseline rerun. Runtime is not a paired performance claim.'}
    (ROOT/'results/stage25_three_source_linear.json').write_text(json.dumps(result,indent=2,default=str)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='cases'},indent=2))
    return result

if __name__=='__main__':main()
