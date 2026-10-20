#!/usr/bin/env python3
"""Exact missing-leaf projection, lifting and two-source fiber controls.

Classical tree projection is explicit; no all-leaf identifiability theorem is
silently applied to unsuppressed partial observations. All writers <=25 frozen.
"""
from fractions import Fraction as F
from collections import deque
from itertools import combinations
from pathlib import Path
from datetime import datetime,timezone
from hashlib import sha256
import json,time,sys
ROOT=Path(__file__).resolve().parents[1]


def geometry(edges):
    adj={v:{} for e in edges for v in e[:2]};lookup={}
    for i,(u,v,ll) in enumerate(edges):
        l=F(ll)
        assert l>0 and u!=v and v not in adj[u]
        adj[u][v]=adj[v][u]=l;lookup[frozenset((u,v))]=i
    assert len(edges)==len(adj)-1 and all(len(ns)!=2 for ns in adj.values())
    seen={next(iter(adj))};todo=list(seen)
    for v in todo:
        for u in adj[v]:
            if u not in seen:seen.add(u);todo.append(u)
    assert len(seen)==len(adj)
    return adj,lookup


def distances(adj,root):
    d={root:F(0)};todo=[root]
    for v in todo:
        for u,l in adj[v].items():
            if u not in d:d[u]=d[v]+l;todo.append(u)
    return d


def observe(edges,sources,y,observed=None):
    adj,_=geometry(edges);y=F(y)
    assert y>0 and sources and sum(F(s['p']) for s in sources)==1
    assert len({s['edge'] for s in sources})==len(sources)
    R=sorted(observed if observed is not None else [v for v in adj if len(adj[v])==1])
    D={v:distances(adj,v) for s in sources for v in edges[s['edge']][:2]}
    result={r:y for r in R}
    for s in sources:
        u,v,l=edges[s['edge']];l=F(l);p,x=F(s['p']),F(s['x'])
        assert p>0 and 0<x<1
        for r in R:result[r]+=p*min(D[u][r]+x*l,D[v][r]+(1-x)*l)
    return result


def observed_backbone(edges,observed_vertices):
    """Return exact observed hull, reduced edges, chains, and omitted gates.

    Original tree must be reduced; observed_vertices must be >=2 original leaves.
    Original vertex IDs are retained at reduced-backbone vertices.
    """
    E=[(u,v,F(l)) for u,v,l in edges];adj,lookup=geometry(E);R=set(observed_vertices)
    assert len(R)>=2 and R<={v for v in adj if len(adj[v])==1}
    active=set(adj);degree={v:len(adj[v]) for v in adj}
    queue=deque(v for v in adj if degree[v]==1 and v not in R)
    while queue:
        v=queue.popleft()
        if v not in active:continue
        active.remove(v)
        for u in adj[v]:
            if u in active:
                degree[u]-=1
                if degree[u]==1 and u not in R:queue.append(u)
    H={v:{u:l for u,l in adj[v].items() if u in active} for v in active}
    kept={v for v in H if len(H[v])!=2};assert {v for v in H if len(H[v])==1}==R
    used=set();B=[];chains=[];original_to_backbone={};interior_vertex={}
    for start in sorted(kept):
        for nxt in sorted(H[start]):
            first=lookup[frozenset((start,nxt))]
            if first in used:continue
            u,v=start,nxt;offset=F(0);chain=[];j=len(B)
            while True:
                i=lookup[frozenset((u,v))];used.add(i);l=H[u][v]
                part={'edge':i,'u':u,'v':v,'start':offset,'end':offset+l,'length':l}
                chain.append(part);original_to_backbone[i]={'backbone_edge':j,**part}
                offset+=l
                if v in kept:break
                interior_vertex[v]={'edge':j,'coordinate':offset}
                w=next(w for w in H[v] if w!=u);u,v=v,w
            B.append((start,v,offset));chains.append(chain)
    gate={v:v for v in H};off_distance={v:F(0) for v in H};todo=list(H)
    for v in todo:
        for u,l in adj[v].items():
            if u not in gate:
                gate[u]=gate[v];off_distance[u]=off_distance[v]+l;todo.append(u)
    gate_edges={v:sorted(lookup[frozenset((v,u))] for u in adj[v] if u not in H) for v in H}
    omitted=[i for i in range(len(E)) if i not in used]
    violations=[{'vertex':v,'hull_degree':len(H[v]),'original_degree':len(adj[v])}
                for v in sorted(H) if len(H[v])>=3 and len(H[v])!=len(adj[v])]
    return {'edges':B,'original_edges':E,'observed_vertices':sorted(R),'hull_vertices':sorted(H),
            'hull_edges':sorted(used),'chains':chains,'original_to_backbone':original_to_backbone,
            'interior_vertex':interior_vertex,'gate_edges':gate_edges,'gate':gate,
            'off_distance':off_distance,'omitted_edges':omitted,
            'saturated_branching':not violations,'saturation_violations':violations}


def project_sources(backbone,sources,y):
    E=backbone['original_edges'];B=backbone['edges'];atoms=[];new_y=F(y)
    assert F(y)>0 and sum(F(s['p']) for s in sources)==1
    for s in sources:
        i,p,x=s['edge'],F(s['p']),F(s['x']);u,v,l=E[i]
        assert p>0 and 0<x<1
        if i in backbone['original_to_backbone']:
            part=backbone['original_to_backbone'][i]
            j=part['backbone_edge']
            coord=part['start']+(x if u==part['u'] else 1-x)*l
            atom={'kind':'edge','edge':j,'x':coord/B[j][2],'p':p,'off_distance':F(0),'origin_edge':i}
        else:
            g=backbone['gate'][u];assert backbone['gate'][v]==g
            a=min(backbone['off_distance'][u]+x*l,backbone['off_distance'][v]+(1-x)*l)
            new_y+=p*a
            if g in backbone['interior_vertex']:
                item=backbone['interior_vertex'][g];j=item['edge']
                atom={'kind':'edge','edge':j,'x':item['coordinate']/B[j][2],'p':p,'off_distance':a,'origin_edge':i,'gate':g}
            else:atom={'kind':'vertex','vertex':g,'p':p,'off_distance':a,'origin_edge':i,'gate':g}
        atoms.append(atom)
    strict=None
    if all(a['kind']=='edge' for a in atoms):
        strict=[]
        for j in sorted({a['edge'] for a in atoms}):
            selected=[a for a in atoms if a['edge']==j];p=sum(a['p'] for a in selected)
            strict.append({'edge':j,'p':p,'x':sum(a['p']*a['x'] for a in selected)/p})
    return {'atoms':atoms,'pendant':new_y,'strict_backbone_sources':strict}


def projected_observe(backbone,projection):
    adj,_=geometry(backbone['edges']);D={r:distances(adj,r) for r in backbone['observed_vertices']}
    out={r:projection['pendant'] for r in D}
    for a in projection['atoms']:
        if a['kind']=='vertex':
            for r in D:out[r]+=a['p']*D[r][a['vertex']]
        else:
            u,v,l=backbone['edges'][a['edge']]
            for r in D:out[r]+=a['p']*min(D[r][u]+a['x']*l,D[r][v]+(1-a['x'])*l)
    return out


def location_fiber(backbone,j,x):
    """Full original-edge fiber of a strictly interior backbone point."""
    coord=F(x)*backbone['edges'][j][2];E=backbone['original_edges']
    assert 0<coord<backbone['edges'][j][2]
    for part in backbone['chains'][j]:
        if part['start']<coord<part['end']:
            u,v,l=E[part['edge']];xx=(coord-part['start'])/l
            if u!=part['u']:xx=1-xx
            return {'kind':'retained_interior','backbone_edge':j,'coordinate':coord,
                    'candidates':[{'edge':part['edge'],'minimum_off_distance':F(0),'fixed_x':xx}]}
        if coord==part['end']:
            g=part['v'];assert g in backbone['interior_vertex']
            candidates=[]
            for i in backbone['omitted_edges']:
                u,v,l=E[i]
                if backbone['gate'][u]==g:
                    assert backbone['gate'][v]==g
                    du,dv=backbone['off_distance'][u],backbone['off_distance'][v]
                    candidates.append({'edge':i,'minimum_off_distance':min(du,dv),'near_endpoint':u if du<dv else v})
            assert candidates
            return {'kind':'omitted_gate','gate':g,'backbone_edge':j,'coordinate':coord,'candidates':candidates}
    raise AssertionError('Point location failed')


def lift_backbone_sources(backbone,sources,y):
    y=F(y);assert y>0 and sum(F(s['p']) for s in sources)==1
    out=[];cost=F(0);E=backbone['original_edges']
    for s in sources:
        p=F(s['p']);fib=location_fiber(backbone,s['edge'],F(s['x']))
        if fib['kind']=='retained_interior':
            c=fib['candidates'][0];out.append({'edge':c['edge'],'p':p,'x':c['fixed_x']})
        else:
            c=next(c for c in fib['candidates'] if c['minimum_off_distance']==0)
            u,v,l=E[c['edge']];a=min(l/2,y/2)
            x=a/l if u==c['near_endpoint'] else 1-a/l
            out.append({'edge':c['edge'],'p':p,'x':x});cost+=p*a
    assert len({s['edge'] for s in out})==len(out) and y-cost>0
    return {'sources':out,'pendant':y-cost,'added_off_distance_cost':cost}


def disjoint_support_certificate(edges,observed):
    """Exact original <=2 supports when reduced mass has four positive endpoints.

    Return not_applicable otherwise, with no assertion of ambiguity/infeasibility.
    """
    B=observed_backbone(edges,observed);adj,_=geometry(B['edges']);D={v:distances(adj,v) for v in adj}
    obs={v:F(d) for v,d in observed.items()};h={v:max(obs[r]-D[v][r] for r in obs) for v in adj}
    m={v:(2-len(adj[v])-sum((h[v]-h[u])/l for u,l in adj[v].items()))/2 for v in adj}
    Y=(sum((2-len(adj[v]))*h[v] for v in adj)-sum(l for u,v,l in B['edges']))/2
    P={v for v,a in m.items() if a>0}
    if any(h[r]!=obs[r] for r in obs) or any(a<0 for a in m.values()) or Y<=0 or len(P)!=4:
        return {'status':'not_applicable','positive_vertices':sorted(P),'backbone':B,'mass':m,'pendant':Y}
    pairs=[list(c) for c in combinations(range(len(B['edges'])),2)
           if {v for i in c for v in B['edges'][i][:2]}==P]
    if len(pairs)!=1:return {'status':'not_applicable','reason':'No unique disjoint pair','backbone':B}
    S=pairs[0];projected=[];fibers=[]
    for j in S:
        u,v,l=B['edges'][j];p=m[u]+m[v];x=m[v]/p
        projected.append({'edge':j,'p':p,'x':x});fibers.append(location_fiber(B,j,x))
    records=[]
    for a in fibers[0]['candidates']:
        for b in fibers[1]['candidates']:
            cs=[a,b];assert a['edge']!=b['edge']
            minimum=sum(p['p']*c['minimum_off_distance'] for p,c in zip(projected,cs))
            if minimum>=Y:continue
            slack=Y-minimum;off_count=sum('fixed_x' not in c for c in cs);sources=[];cost=F(0)
            for p,c in zip(projected,cs):
                weight=p['p'];u,v,l=B['original_edges'][c['edge']]
                if 'fixed_x' in c:x=c['fixed_x']
                else:
                    extra=min(l/2,slack/(2*off_count*weight))
                    x=extra/l if u==c['near_endpoint'] else 1-extra/l
                    cost+=weight*(c['minimum_off_distance']+extra)
                sources.append({'edge':c['edge'],'p':weight,'x':x})
            witness={'sources':sources,'pendant':Y-cost}
            assert observe(edges,sources,Y-cost,obs)==obs
            records.append({'support':sorted(s['edge'] for s in sources),'minimum_off_cost':minimum,'witness':witness})
    return {'status':'certified_unique_original_support' if len(records)==1 else 'certified_original_support_ambiguity' if records else 'certified_infeasible',
            'backbone':B,'mass':m,'pendant':Y,'projected_sources':projected,'fibers':fibers,'original_supports':records,
            'parameters_unique':len(records)==1 and all(f['kind']=='retained_interior' for f in fibers),
            'scope':'All strict original-edge explanations with at most two sources; positive free pendant'}


def source(edge,p,x):return {'edge':edge,'p':F(p),'x':F(x)}


def main():
    start=time.perf_counter();records=[];E=[(4,0,F(1)),(4,1,F(1)),(4,5,F(1)),(5,2,F(1)),(5,3,F(1))];R=[0,1,2]
    worlds=[{'sources':[source(0,F(2,5),F(1,2)),source(2,F(3,5),F(2,3))],'pendant':F(1)},
            {'sources':[source(0,F(7,10),F(2,7)),source(3,F(3,10),F(1,3))],'pendant':F(1)}]
    B=observed_backbone(E,R);obs=[observe(E,w['sources'],w['pendant'],R) for w in worlds];assert obs[0]==obs[1]
    missing=[observe(E,w['sources'],w['pendant'])[3] for w in worlds];assert missing==[F(14,5),F(3)]
    projections=[project_sources(B,w['sources'],w['pendant']) for w in worlds]
    assert all(projected_observe(B,p)==obs[0] for p in projections)
    assert disjoint_support_certificate(E,obs[0])['status']=='not_applicable'
    records.append({'name':'one_missing_leaf_two_distinct_two_edge_supports','edges':E,'observed':obs[0],'worlds':worlds,
                    'missing_leaf_values':missing,'backbone':B,'projections':projections})
    worlds=[{'sources':[source(2,1,F(1,2))],'pendant':F(1)},
            {'sources':[source(2,F(4,5),F(1,4)),source(3,F(1,5),F(1,2))],'pendant':F(1)}]
    assert observe(E,worlds[0]['sources'],1,R)==observe(E,worlds[1]['sources'],1,R)
    ps=[project_sources(B,w['sources'],w['pendant']) for w in worlds]
    assert ps[0]['strict_backbone_sources']==ps[1]['strict_backbone_sources']
    records.append({'name':'suppression_merges_distinct_original_source_edges','edges':E,'observed':observe(E,worlds[0]['sources'],1,R),
                    'worlds':worlds,'backbone':B,'projections':ps})
    # Unsaturated center: original minimum1, strict backbone minimum3.
    star=[(4,i,F(1)) for i in range(4)];BB=observed_backbone(star,[0,1,2]);original={'sources':[source(3,1,F(1,2))],'pendant':F(1)}
    pr=project_sources(BB,original['sources'],original['pendant']);dd=observe(star,original['sources'],1,[0,1,2])
    assert not BB['saturated_branching'] and pr['atoms'][0]['kind']=='vertex' and pr['strict_backbone_sources'] is None
    assert projected_observe(BB,pr)==dd
    from stage25_three_source_linear import solve_three_linear
    sol=solve_three_linear(BB['edges'],dd);assert sol['minimum_sources_within_cap']==3
    records.append({'name':'unsaturated_nonbinary_one_becomes_three','edges':star,'original':original,'observed':dd,'backbone':BB,'projection':pr,'strict_backbone_solution':sol})
    # Binary and saturated nonbinary forward/reverse controls, including gate hits.
    binary=[(6,0,F(1)),(6,1,F(1)),(6,7,F(1)),(7,2,F(1)),(7,8,F(1)),(8,3,F(1)),(8,9,F(1)),(9,4,F(1)),(9,5,F(1))]
    sat_nonbinary=[(8,0,F(1)),(8,1,F(2)),(8,2,F(3)),(8,9,F(1)),(9,3,F(2)),(9,4,F(1))]
    roundtrips=[]
    for name,EE,RR in [('quartet',E,R),('binary6',binary,[0,1,3,4,5]),('saturated_nonbinary',sat_nonbinary,[0,1,2,3])]:
        bb=observed_backbone(EE,RR);assert bb['saturated_branching']
        selections=[(i,) for i in range(len(EE))]+list(combinations(range(len(EE)),2))
        for S in selections:
            ss=[source(i,F(1,len(S)),F(j+1,len(S)+2)) for j,i in enumerate(S)];yy=F(1,1000)
            dd=observe(EE,ss,yy,RR);pp=project_sources(bb,ss,yy);assert pp['strict_backbone_sources'] is not None
            assert len(pp['strict_backbone_sources'])<=len(ss)
            assert observe(bb['edges'],pp['strict_backbone_sources'],pp['pendant'])==dd==projected_observe(bb,pp)
            lift=lift_backbone_sources(bb,pp['strict_backbone_sources'],pp['pendant'])
            assert observe(EE,lift['sources'],lift['pendant'],RR)==dd
            roundtrips.append({'tree':name,'edges':EE,'observed_vertices':RR,'original_sources':ss,'original_pendant':yy,
                               'observed':dd,'projection':pp,'lift':lift})
    # Four positive nonadjacent endpoints: direct edge, single omitted edge,
    # and deeper fiber candidates controlled by the pendant budget.
    certs=[]
    cat=[(5,0,F(1)),(5,1,F(1)),(5,6,F(1)),(6,2,F(1)),(6,7,F(1)),(7,3,F(1)),(7,4,F(1))]
    configs=[('nongate_exact_original_parameters',cat,[0,1,3,4],[source(0,F(1,3),F(2,5)),source(5,F(2,3),F(3,7))],F(1),1,True),
             ('single_missing_edge_gate_support_only',binary,[0,1,3,4,5],[source(3,F(1,3),F(1,2)),source(7,F(2,3),F(2,5))],F(1),1,False)]
    deep=binary[:3]+[(7,10,F(2)),(10,2,F(1)),(10,11,F(1))]+binary[4:]
    for Y in [F(1,2),F(2,3),F(1)]:
        configs.append(('fiber_budget_'+str(Y),deep,[0,1,3,4,5],[source(3,F(1,3),F(1,4)),source(9,F(2,3),F(2,5))],Y-F(1,6),1 if Y<=F(2,3) else 3,False))
    for name,EE,RR,ss,yy,count,unique_parameters in configs:
        dd=observe(EE,ss,yy,RR);cert=disjoint_support_certificate(EE,dd)
        assert len(cert['original_supports'])==count and cert['parameters_unique']==unique_parameters,(name,cert)
        certs.append({'name':name,'edges':EE,'observed':dd,'generating_sources':ss,'generating_pendant':yy,'certificate':cert})
    # Independent fixed-support LP checks use missing observation rows, no old main.
    from stage25_bounded_error import solve_fixed_support
    lp=[]
    for c in certs:
        certificate=c['certificate'];allowed={tuple(x['support']) for x in certificate['original_supports']}
        allfiber=list(__import__('itertools').product(*[f['candidates'] for f in certificate['fibers']]))
        for cs in allfiber:
            S=sorted(cand['edge'] for cand in cs)
            out=solve_fixed_support(c['edges'],c['observed'],F(0),S)
            assert out['status']!='unknown'
            assert (out['status']=='strict_feasible')==(tuple(S) in allowed)
            lp.append({'case':c['name'],'support':S,'solution':out})
    for i,w in enumerate(records[0]['worlds']):
        out=solve_fixed_support(E,records[0]['observed'],F(0),[s['edge'] for s in w['sources']]);assert out['status']=='strict_feasible'
        lp.append({'case':'one_missing_leaf_counterexample_'+str(i),'support':[s['edge'] for s in w['sources']],'solution':out})
    # Nonbinary closed projected center passes no <=2 strict support at zero error.
    from stage25_bounded_error import solve_all_supports
    allres=solve_all_supports(BB['edges'],records[2]['observed'],F(0));assert allres['status']=='strict_infeasible'
    output={'passed':True,'created_utc':datetime.now(timezone.utc).isoformat(),'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
            'counterexamples':records,'saturated_roundtrip_count':len(roundtrips),'saturated_roundtrips':roundtrips,
            'fiber_certificate_cases':certs,'fixed_support_lp_checks':lp,'nonbinary_backbone_at_most_two':allres,
            'dependencies':{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ['scripts/stage25_bounded_error.py','scripts/stage25_three_source_linear.py','scripts/stage25_deficit_extension.py','scripts/stage23_two_source_envelope.py']},
            'elapsed_seconds':time.perf_counter()-start,'scope':'Exact synthetic missing-leaf geometry and standard LP controls; no noise, biological performance, or novelty claim'}
    (ROOT/'results/stage26_missing_leaf.json').write_text(json.dumps(output,indent=2,default=str)+'\n')
    print(json.dumps({k:v for k,v in output.items() if k not in ['counterexamples','saturated_roundtrips','fiber_certificate_cases','fixed_support_lp_checks','nonbinary_backbone_at_most_two']},indent=2))
    return output

if __name__=='__main__':main()
