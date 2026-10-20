#!/usr/bin/env python3
"""Bounded single-thread scale controls for exact low-source tree solvers.

Three repetitions per completed solver/case cell. Child-process time limits
produce censored unknown measurements, never mathematical infeasibility.
No prior writer is called; production inputs contain exact rational numbers.
"""
import os
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[_key]='1'
from fractions import Fraction as F
from pathlib import Path
from datetime import datetime,timezone
from statistics import median
import multiprocessing as mp
import resource,signal,math,time,json,hashlib,random,platform
from stage24_three_source import solve_three as old_three
from stage25_three_source_linear import solve_three_linear
from stage25_deficit_extension import solve_four
ROOT=Path(__file__).resolve().parents[1]
SOLVERS={'new_three_linear':solve_three_linear,'old_three_quadratic':old_three,'four_local':solve_four}
LIMITS={'new_three_linear':5.0,'old_three_quadratic':3.0,'four_local':6.0}
TOTAL_CPU_BUDGET=175.0


def adjacency(edges):
    a={}
    for u,v,l in edges:
        l=F(l);assert l>0
        a.setdefault(u,{})[v]=l;a.setdefault(v,{})[u]=l
    assert len(edges)==len(a)-1 and all(len(ns)!=2 for ns in a.values())
    return a


def caterpillar(n,weighted=True):
    inner=list(range(n,2*n-2))
    pairs=[(0,inner[0]),(1,inner[0])]
    pairs += [(i+2,v) for i,v in enumerate(inner[1:-1])]
    pairs += [(n-2,inner[-1]),(n-1,inner[-1])]
    pairs += list(zip(inner,inner[1:]))
    lengths=[F(1,3),F(2,5),F(7,11),F(13,7),F(3),F(1,97)]
    return [(u,v,lengths[i%len(lengths)] if weighted else F(1)) for i,(u,v) in enumerate(pairs)]


def binary(n,kind,seed):
    rng=random.Random(seed);pairs=[(0,3),(1,3),(2,3)];depth={0:1,1:1,2:1};next_id=4
    while len(depth)<n:
        leaf=min(depth,key=lambda v:(depth[v],v)) if kind=='balanced' else rng.choice(sorted(depth))
        i=next(i for i,e in enumerate(pairs) if leaf in e);u,v=pairs.pop(i);p=v if u==leaf else u
        new,next_leaf=next_id,next_id+1;next_id+=2
        pairs.extend([(p,new),(leaf,new),(next_leaf,new)])
        d=depth[leaf]+1;depth[leaf]=depth[next_leaf]=d
    pool=[F(1,7),F(2,3),F(3,5),F(11,4),F(1,101),F(7)]
    return [(u,v,pool[rng.randrange(len(pool))]) for u,v in pairs]


def observations(edges,sources,y):
    """Linear reroot generation; separate validator uses original path scans."""
    adj=adjacency(edges);mass={v:F(0) for v in adj}
    for s in sources:
        u,v,l=edges[s['edge']];p,x=s['p'],s['x'];assert p>0 and 0<x<1
        mass[u]+=p*(1-x);mass[v]+=p*x
    assert sum(mass.values())==1 and y>0
    root=next(iter(adj));parent={root:None};order=[root];depth={root:F(0)}
    for v in order:
        for u,l in adj[v].items():
            if u!=parent[v]:parent[u]=v;depth[u]=depth[v]+l;order.append(u)
    subtree=mass.copy()
    for v in reversed(order[1:]):subtree[parent[v]]+=subtree[v]
    pred={root:y+sum(mass[v]*depth[v] for v in adj)}
    for v in order[1:]:pred[v]=pred[parent[v]]+adj[v][parent[v]]*(1-2*subtree[v])
    return {v:pred[v] for v in adj if len(adj[v])==1}


def boundary_sources(edges,k):
    ids=[len(edges)*i//k for i in range(k)]
    raw=[1]+[10**8+i for i in range(1,k)];total=sum(raw)
    xs=[F(1,10**8),1-F(1,10**8),F(2,7),F(5,13)]
    return [{'edge':e,'p':F(w,total),'x':xs[i]} for i,(e,w) in enumerate(zip(ids,raw))]


def case(name,shape,edges,sources,y,models,expected_min=None,expected_reject_cap=None,extra=None):
    adj=adjacency(edges)
    out={'name':name,'shape':shape,'edges':edges,'observed':observations(edges,sources,y),
         'leaf_count':sum(len(a)==1 for a in adj.values()),'vertex_count':len(adj),
         'generating_source_count':len(sources),'generating_pendant':y,
         'generating_sources':sources,'models':models,
         'expected_minimum_from_analytic_argument':expected_min,'analytic_reject_up_to':expected_reject_cap}
    if extra:out.update(extra)
    return out


def relabel_reorder(original,name,seed):
    rng=random.Random(seed);E=original['edges'];vertices=sorted({v for e in E for v in e[:2]})
    labels=[100003+37*i for i in range(len(vertices))];rng.shuffle(labels);mapping=dict(zip(vertices,labels))
    order=list(range(len(E)));rng.shuffle(order);new=[];edge_map={}
    flips={}
    for j,i in enumerate(order):
        u,v,l=E[i];flip=rng.choice([False,True]);flips[i]=flip
        new.append((mapping[v],mapping[u],l) if flip else (mapping[u],mapping[v],l));edge_map[i]=j
    src=[{'edge':edge_map[s['edge']],'p':s['p'],'x':1-s['x'] if flips[s['edge']] else s['x']} for s in original['generating_sources']]
    out=case(name,original['shape']+'_relabeled_reordered',new,src,original['generating_pendant'],original['models'],
             original['expected_minimum_from_analytic_argument'],original['analytic_reject_up_to'],
             {'invariance_pair':original['name'],'vertex_mapping':mapping,'old_to_new_edge':edge_map})
    assert out['observed']=={mapping[v]:d for v,d in original['observed'].items()}
    return out


def center_profile(n,model,late=True):
    E=binary(n,'balanced',26000+n);adj=adjacency(E)
    # A degree-three vertex with three internal neighbors where available.
    internal=[v for v in adj if len(adj[v])==3]
    v=max(internal,key=lambda v:sum(len(adj[u])>1 for u in adj[v]))
    incident=[e for e in E if v in e[:2]];other=[e for e in E if v not in e[:2]]
    E=other+incident if late else incident+other
    radius=min(adj[v].values())/4
    src=[]
    for i,(u,w,l) in enumerate(E):
        if v in (u,w):src.append({'edge':i,'p':F(1,3),'x':radius/l if u==v else 1-radius/l})
    name=('late' if late else 'early')+'_center_'+str(n)+'_'+model
    models=['new_three_linear','old_three_quadratic'] if model=='three' and n<=256 else ['new_three_linear'] if model=='three' else ['four_local']
    return case(name,'balanced_center_profile',E,src,F(1),models,3,extra={
        'stress':'isolated envelope baseline center; only that positive deficit can yield <=3; center ordered '+('last' if late else 'first')+' among eligible internal vertices',
        'center_vertex':v,'physical_radius':radius})


def build_cases():
    C=[]
    for n in (64,1000):
        E=[(i,n,F(1+(i%5),1+(i%3))) for i in range(n)]
        C.append(case('star_'+str(n)+'_near_two','unequal_star',E,boundary_sources(E,2),F(1,10**12),
                      ['new_three_linear','old_three_quadratic'] if n==64 else ['new_three_linear']))
    E=[(i,1000,F(1)) for i in range(1000)]
    C.append(case('star_1000_requires_1000','unit_star',E,[{'edge':i,'p':F(1,1000),'x':F(1,2)} for i in range(1000)],F(1),['new_three_linear'],1000,3))
    for n in (32,128,1000):
        E=caterpillar(n)
        C.append(case('caterpillar_'+str(n)+'_near_three','unequal_caterpillar',E,boundary_sources(E,3),F(1,10**12),
                      ['new_three_linear','old_three_quadratic'] if n<=128 else ['new_three_linear']))
    for shape,n in [('balanced',128),('balanced',1000),('random',64),('random',1000)]:
        E=binary(n,shape,26100+n)
        C.append(case(shape+'_'+str(n)+'_near_three',shape+'_unequal_binary',E,boundary_sources(E,3),F(1,10**12),
                      ['new_three_linear','old_three_quadratic'] if n<=128 else ['new_three_linear']))
    original=next(c for c in C if c['name']=='random_64_near_three')
    C.append(relabel_reorder(original,'random_64_relabel_reorder',26213))
    for n in (64,256,1000):C.append(center_profile(n,'three',True))
    C.append(center_profile(64,'three',False));C[-1]['invariance_pair']='late_center_64_three'
    for n in (32,128,512,1000):
        E=caterpillar(n,False);adj=adjacency(E);le=[i for i,(u,v,l) in enumerate(E) if len(adj[u])==1 or len(adj[v])==1]
        src=[{'edge':i,'p':F(1,n),'x':F(1,2)} for i in le]
        C.append(case('dense_caterpillar_'+str(n),'unit_caterpillar_dense_leaf_sources',E,src,F(1),
                      ['new_three_linear','old_three_quadratic'] if n<=512 else ['new_three_linear'],expected_reject_cap=3,
                      extra={'stress':'baseline mass is leaf-degree/n on every internal vertex, zero on leaves; >10 positives excludes <=3 independently of generating source count'}))
    # Four-source cells: sparse two-deficit geometry, high-degree stars, and
    # a center-profile quadratic scan stress despite a known three-source fit.
    E=[(i,200,F(1+(i%5),1+(i%3))) for i in range(200)]
    C.append(case('star_200_near_four','unequal_star',E,boundary_sources(E,4),F(1,10**12),['four_local']))
    for n in (16,64,200):
        E=caterpillar(n,False);adj=adjacency(E);le=[i for i,(u,v,l) in enumerate(E) if len(adj[u])==1 or len(adj[v])==1]
        selected=le[:2]+le[-2:];src=[{'edge':i,'p':F(1,4),'x':F(1,2)} for i in selected]
        C.append(case('two_cherries_'+str(n)+'_four','unit_caterpillar_two_cherries',E,src,F(1),['four_local'],4,
                      extra={'stress':'two separated cherry deficits; distinct negative neighbor masses require both deficits'}))
    for n in (64,200):
        E=binary(n,'random',26300+n)
        C.append(case('random_'+str(n)+'_near_four','random_unequal_binary',E,boundary_sources(E,4),F(1,10**12),['four_local']))
    for n in (16,64,200):C.append(center_profile(n,'four',True))
    E=[(i,200,F(1)) for i in range(200)]
    C.append(case('star_200_requires_200','unit_star',E,[{'edge':i,'p':F(1,200),'x':F(1,2)} for i in range(200)],F(1),['four_local'],200,4))
    return C


def compact(sol):
    keys=['feasible','minimum_sources_within_cap','support','sources','pendant','deficit_vertex','deficit',
          'deficit_vertices','deficits','statistics','candidate_potentials_examined','reason',
          'baseline_positive_count','baseline_negative_count']
    return {k:sol[k] for k in keys if k in sol}


def worker(conn,model,edges,observed,limit):
    try:
        resource.setrlimit(resource.RLIMIT_CPU,(math.ceil(limit)+1,math.ceil(limit)+2))
        wall=time.perf_counter();cpu=time.process_time()
        sol=SOLVERS[model](edges,observed)
        elapsed=time.perf_counter()-wall;cp=time.process_time()-cpu
        conn.send({'status':'complete','wall_seconds':elapsed,'cpu_seconds':cp,'solution':compact(sol)})
    except BaseException as exc:
        conn.send({'status':'error','error_type':type(exc).__name__,'message':str(exc)})
    finally:conn.close()


def timed_call(model,edges,obs):
    ctx=mp.get_context('fork');parent,child=ctx.Pipe(False);limit=LIMITS[model]
    process=ctx.Process(target=worker,args=(child,model,edges,obs,limit));process.start();child.close()
    try:
        if parent.poll(limit+0.5):
            try:out=parent.recv()
            except EOFError:out={'status':'censored','reason':'child_exited_without_result','exitcode':process.exitcode}
        else:out={'status':'censored','reason':'wall_time_limit','limit_seconds':limit+0.5}
    finally:
        if process.is_alive():
            process.join(0.1)
            if process.is_alive():process.terminate();process.join(0.2)
            if process.is_alive():process.kill();process.join()
        else:process.join()
        parent.close()
    return out


def verify(edges,observed,solution,cap):
    if not solution['feasible']:return
    adj=adjacency(edges);D={};src=solution['sources'];y=F(solution['pendant'])
    assert 1<=len(src)<=cap and len(src)==solution['minimum_sources_within_cap']
    assert y>0 and sum(s['p'] for s in src)==1 and len({s['edge'] for s in src})==len(src)
    for s in src:
        assert s['p']>0 and 0<s['x']<1
        u,v,l=edges[s['edge']]
        for root in (u,v):
            if root in D:continue
            D[root]={root:F(0)};todo=[root]
            for a in todo:
                for b,ll in adj[a].items():
                    if b not in D[root]:D[root][b]=D[root][a]+ll;todo.append(b)
    for r,d in observed.items():
        pred=y
        for s in src:
            u,v,l=edges[s['edge']];x,p=s['x'],s['p']
            pred+=p*min(D[u][r]+x*l,D[v][r]+(1-x)*l)
        assert pred==d
    assert sorted(s['edge'] for s in src)==sorted(solution['support'])


def cpu_total():
    own=resource.getrusage(resource.RUSAGE_SELF);children=resource.getrusage(resource.RUSAGE_CHILDREN)
    return own.ru_utime+own.ru_stime+children.ru_utime+children.ru_stime


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    started=time.perf_counter();cpu0=cpu_total();cases=build_cases();records=[]
    deps=['scripts/stage24_three_source.py','scripts/stage23_two_source_envelope.py',
          'scripts/stage25_three_source_linear.py','scripts/stage25_deficit_extension.py']
    out={'created_utc':datetime.now(timezone.utc).isoformat(),'script_sha256':sha(Path(__file__)),
         'dependencies':{p:sha(ROOT/p) for p in deps},'python':platform.python_version(),'platform':platform.platform(),
         'plan':{'repetitions':3,'single_thread':True,'models':list(SOLVERS),'per_call_wall_limit_seconds':{m:v+0.5 for m,v in LIMITS.items()},
                 'total_cpu_budget_seconds':TOTAL_CPU_BUDGET,'fixed_case_order':[c['name'] for c in cases],
                 'timing':'solver-only CPU and wall time in isolated fork child; inputs generated and witnesses verified outside timed region; alternate old/new order by repetition'},
         'cases':records,'validation_passed':False,'completed_planned_suite':False}
    path=ROOT/'results/stage26_scale.json'
    def save():
        out['elapsed_seconds']=time.perf_counter()-started;out['cpu_seconds_including_children']=cpu_total()-cpu0
        path.write_text(json.dumps(out,indent=2,default=str)+'\n')
    seen={}
    for c in cases:
        entry=dict(c);entry['runs']={m:[] for m in c['models']};entry['summaries']={};blocked=set();records.append(entry)
        for rep in range(3):
            order=c['models'][:] if rep%2==0 else c['models'][::-1]
            for model in order:
                if model in blocked:continue
                if cpu_total()-cpu0>=TOTAL_CPU_BUDGET:
                    blocked.add(model);entry['runs'][model].append({'status':'not_run','reason':'total_cpu_budget','repetition':rep});continue
                run=timed_call(model,c['edges'],c['observed']);run['repetition']=rep
                entry['runs'][model].append(run)
                if run['status']!='complete':
                    assert run['status']!='error',(c['name'],model,run)
                    blocked.add(model);continue
                sol=run['solution'];cap=4 if model=='four_local' else 3
                verify(c['edges'],c['observed'],sol,cap)
                expected=c['expected_minimum_from_analytic_argument']
                if expected is not None:
                    assert sol['feasible']==(expected<=cap),(c['name'],sol,expected)
                    if sol['feasible']:assert sol['minimum_sources_within_cap']==expected
                elif c['analytic_reject_up_to'] is not None and cap<=c['analytic_reject_up_to']:
                    assert not sol['feasible']
                elif c['generating_source_count']<=cap:assert sol['feasible']
        for model,runs in entry['runs'].items():
            complete=[r for r in runs if r['status']=='complete']
            if complete:
                answers={(r['solution']['feasible'],r['solution'].get('minimum_sources_within_cap')) for r in complete};assert len(answers)==1
                answer=next(iter(answers));key=(c.get('invariance_pair',c['name']),model)
                if key in seen:assert seen[key]==answer
                seen[key]=answer
                entry['summaries'][model]={'completed_repetitions':len(complete),'all_three_complete':len(complete)==3,
                    'answer':{'feasible':answer[0],'minimum_sources_within_cap':answer[1]},
                    'median_wall_seconds':median(r['wall_seconds'] for r in complete) if len(complete)==3 else None,
                    'median_cpu_seconds':median(r['cpu_seconds'] for r in complete) if len(complete)==3 else None,
                    'wall_range_seconds':[min(r['wall_seconds'] for r in complete),max(r['wall_seconds'] for r in complete)],
                    'statistics':complete[0]['solution'].get('statistics',{}),
                    'candidate_potentials_examined':complete[0]['solution'].get('candidate_potentials_examined')}
            else:entry['summaries'][model]={'completed_repetitions':0,'all_three_complete':False,'answer':'unknown_due_to_limit'}
        summaries=entry['summaries']
        if all(m in summaries and summaries[m]['completed_repetitions'] for m in ('new_three_linear','old_three_quadratic')):
            assert summaries['new_three_linear']['answer']==summaries['old_three_quadratic']['answer'];entry['paired_answers_agree']=True
        # An untimed cap-three cross-check separates a four-source witness from
        # a known minimum, without using generating source count as the oracle.
        if 'four_local' in summaries and summaries['four_local']['completed_repetitions']:
            small=solve_three_linear(c['edges'],c['observed']);verify(c['edges'],c['observed'],small,3)
            ans=summaries['four_local']['answer']
            if small['feasible']:assert ans['feasible'] and ans['minimum_sources_within_cap']==small['minimum_sources_within_cap']
            elif ans['feasible']:assert ans['minimum_sources_within_cap']==4
            entry['untimed_cap_three_crosscheck']=compact(small)
        save()
        print(json.dumps({'case':c['name'],'summaries':{m:{k:v for k,v in s.items() if k in ['all_three_complete','answer','median_wall_seconds']} for m,s in summaries.items()},'cpu_total':cpu_total()-cpu0}),flush=True)
    assert all(sha(ROOT/p)==h for p,h in out['dependencies'].items())
    out['validation_passed']=True;out['completed_planned_suite']=all(any(r['status']!='not_run' for runs in c['runs'].values() for r in runs) for c in records)
    out['scope']='Selected exact-rational generated cases and transformations; censored runs are unknown. Witness validity and independent cap-three comparisons do not identify biological truth. Three-repeat timings are local measurements, not universal software performance claims.'
    save()
    print(json.dumps({'validation_passed':True,'case_count':len(records),'elapsed_seconds':out['elapsed_seconds'],'cpu_seconds':out['cpu_seconds_including_children'],'completed_planned_suite':out['completed_planned_suite']},indent=2),flush=True)
    return out

if __name__=='__main__':main()
