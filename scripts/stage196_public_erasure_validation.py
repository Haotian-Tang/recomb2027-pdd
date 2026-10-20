#!/usr/bin/env python3
"""Stage196 bounded, preregistered erasure validation; preparation is metadata only."""
import sys
sys.dont_write_bytecode = True
import argparse
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import resource
import subprocess
import time
import traceback
from collections import Counter
from fractions import Fraction as Q

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = 'scripts/stage196_public_erasure_validation.py'
DESIGN = 'results/stage196_root_erasure_validation_design.json'
DATA = 'data/stage196_public_erasure_cases.json'
PROTOCOL = 'results/stage196_public_erasure_protocol.json'
PREFIX = 'results/stage196_public_erasure_validation'
STATES = ('D', 'K00', 'K01', 'K11')
LIMITS = dict(DP_runs=7, DP_candidate_transitions=2000000, separator_subsets=3000,
              wall_seconds=180, address_space_bytes=1073741824, physical_worlds=6,
              original_vertex_physical_evaluations=60, flow_calls=0, rank_calls=0,
              LP_calls=0, retries=0)


def need(ok, message):
    if not ok:
        raise AssertionError(message)


def encode(x):
    if isinstance(x, Q):
        return str(x)
    if isinstance(x, dict):
        return {str(k): encode(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [encode(v) for v in x]
    return x


def read(path):
    return json.loads((ROOT / path).read_text())


def sha(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def save(path, obj):
    with (ROOT / path).open('x') as f:
        json.dump(encode(obj), f, indent=2)
        f.write('\n')
    return dict(path=path, sha256=sha(path))


def event(obj):
    with (ROOT / (PREFIX + '_events.jsonl')).open('a') as f:
        f.write(json.dumps(encode(obj)) + '\n')
        f.flush()


class Limit(Exception):
    pass


class Budget:
    def __init__(self):
        self.deadline = time.monotonic() + LIMITS['wall_seconds']
        self.counts = Counter()

    def check(self):
        if time.monotonic() >= self.deadline:
            raise Limit('whole_worker_wall_budget')

    def charge(self, key, maximum=None):
        self.check()
        if maximum is not None and self.counts[key] >= maximum:
            raise Limit(key + '_budget')
        self.counts[key] += 1


def topology(incidence, support, budget):
    labels = list(dict.fromkeys(v for e in incidence for v in e))
    need(all(isinstance(v, (str, int)) and not isinstance(v, bool) for v in labels), 'labels')
    index = {v: i for i, v in enumerate(labels)}
    edges = []
    for e in incidence:
        budget.check()
        need(len(e) == 2 and e[0] != e[1], 'simple edge')
        edges.append([index[e[0]], index[e[1]]])
    n = len(labels)
    need(n >= 2 and len(edges) == n - 1 and len({frozenset(e) for e in edges}) == n - 1, 'tree size')
    adj = [[] for _ in labels]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    parent = [-2] * n
    parent[0] = -1
    order = [0]
    for u in order:
        budget.check()
        for v in adj[u]:
            if v == parent[u]:
                continue
            need(parent[v] == -2, 'cycle')
            parent[v] = u
            order.append(v)
    need(len(order) == n, 'connected')
    need(support and len(set(support)) == len(support) and all(type(e) is int and 0 <= e < len(edges) for e in support), 'support')
    children = [[v for v in adj[u] if parent[v] == u] for u in range(n)]
    P = {v for e in support for v in edges[e]}
    B = {v for v in range(n) if len(adj[v]) == 1}
    return dict(labels=labels, index=index, edges=edges, adjacency=adj, parent=parent,
                order=order, children=children, P=P, B=B, H=len(support), support=list(support))


def reachable_leaves(model, C, budget):
    C = set(C)
    reached = model['P'] - C
    stack = list(reached)
    while stack:
        budget.check()
        u = stack.pop()
        for v in model['adjacency'][u]:
            if v not in C and v not in reached:
                reached.add(v)
                stack.append(v)
    return sorted(reached & model['B']), sorted(reached)


def verify_rank_certificate(model, certificate, budget):
    paths = [[model['index'][v] for v in path] for path in certificate['paths']]
    C = [model['index'][v] for v in certificate['separator']]
    occupied = set()
    E = {frozenset(e) for e in model['edges']}
    for path in paths:
        budget.check()
        need(path and path[0] in model['P'] and path[-1] in model['B'], 'rank path terminals')
        need(len(path) == len(set(path)) and not occupied.intersection(path), 'disjoint rank paths')
        need(all(frozenset((u, v)) in E for u, v in zip(path, path[1:])), 'original rank path')
        occupied.update(path)
    s = len(paths)
    need(s >= 1 and len(set(C)) == len(C) == s, 'equal nonzero rank certificates')
    D, reached = reachable_leaves(model, C, budget)
    need(not D, 'rank separator disconnects all P from all B')
    budget.counts['rank_certificate_paths'] += s
    return dict(status='PASS_ORIGINAL_PATH_SEPARATOR',s=s, paths=paths, separator=C,
                reached_after_separator=reached, canonical_dimension=2*model['H']-s,
                canonical_unique=s == 2*model['H'], no_new_flow_or_rank=True)


def merge_state(left, right):
    if left == 'D':
        return 'D' if right in ('D', 'K00', 'K11') else None
    if right == 'D':
        return left
    h, a = int(left[1]), int(left[2])
    j, c = int(right[1]), int(right[2])
    return 'K' + str(max(h, j)) + str(a) if a == c else None


def rows_from(values, argmins, K):
    return [dict(budget=b, state=state, cost=values[(b, state)], argmin=argmins.get((b, state)))
            for b in range(K+1) for state in STATES]


def row_map(table):
    return {(r['budget'], r['state']): r for r in table['rows']}


def solve_dp(model, s, budget):
    """Receives only checked topology/support/rank, no expected delta or metric."""
    budget.charge('DP_runs', 7)
    K = s - 1
    raw = dict(status='UNKNOWN', s=s, K=K, tables=[], node_final_tables={},
               root=0, vertex_order=model['labels'], rooted_parent=model['parent'],
               children=model['children'], traceback=[])
    finals = {}
    pending = None
    try:
        for v in reversed(model['order']):
            budget.check()
            values = {(b, state): None for b in range(K+1) for state in STATES}
            args = {}
            if K >= 1:
                values[(1, 'D')] = 0
            h = int(v in model['P'])
            for a in range(h, 2):
                values[(0, 'K'+str(h)+str(a))] = a*int(v in model['B'])
            tid = len(raw['tables'])
            raw['tables'].append(dict(id=tid, node=v, processed_children=0, kind='initialize',
                                      rows=rows_from(values, args, K)))
            for step, child in enumerate(model['children'][v], 1):
                child_id = finals[child]
                right = row_map(raw['tables'][child_id])
                new = {(b, state): None for b in range(K+1) for state in STATES}
                args = {}
                transitions = []
                next_id = len(raw['tables'])
                pending = dict(id=next_id, node=v, processed_children=step, kind='merge',
                    child=child, left_table=tid, right_table=child_id, completion='IN_PROGRESS',
                    rows=rows_from(new,args,K), transitions=transitions)
                raw['tables'].append(pending)
                for b1 in range(K+1):
                    for state1 in STATES:
                        for b2 in range(K-b1+1):
                            for state2 in STATES:
                                target = merge_state(state1, state2)
                                if target is None:
                                    continue
                                budget.charge('DP_candidate_transitions', 2000000)
                                lcost = values[(b1, state1)]
                                rcost = right[(b2, state2)]['cost']
                                cost = None if lcost is None or rcost is None else lcost+rcost
                                key = (b1+b2, target)
                                j = len(transitions)
                                transitions.append(dict(left=[b1,state1], right=[b2,state2], target=list(key), cost=cost))
                                if cost is not None and (new[key] is None or cost < new[key]):
                                    new[key] = cost
                                    args[key] = j
                pending.update(rows=rows_from(new,args,K),completion='COMPLETE')
                pending = None
                tid = next_id
                values = new
            finals[v] = tid
            raw['node_final_tables'][str(v)] = tid
        root_table = finals[0]
        accepted = []
        chosen = None
        for row in raw['tables'][root_table]['rows']:
            if row['state'] not in ('D','K00','K11'):
                continue
            accepted.append(dict(budget=row['budget'],state=row['state'],cost=row['cost']))
            if row['cost'] is not None and (chosen is None or row['cost'] < chosen['cost']):
                chosen = row
        need(chosen is not None, 'root has a finite closure')
        C = set()
        stack = [(root_table, (chosen['budget'], chosen['state']))]
        while stack:
            budget.check()
            table_id, key = stack.pop()
            table = raw['tables'][table_id]
            row = row_map(table)[key]
            need(row['cost'] is not None, 'finite traceback')
            raw['traceback'].append(dict(table=table_id,key=list(key),cost=row['cost'],argmin=row['argmin']))
            if table['kind'] == 'initialize':
                if key[1] == 'D':
                    C.add(table['node'])
            else:
                tr = table['transitions'][row['argmin']]
                stack.append((table['right_table'],tuple(tr['right'])))
                stack.append((table['left_table'],tuple(tr['left'])))
        D, reached = reachable_leaves(model, C, budget)
        raw.update(status='DP_CANDIDATE',delta=chosen['cost'],separator_C=sorted(C),erased_leaves_D=D,
                   reached_vertices=reached,root_table=root_table,accepted_root_states=accepted,
                   chosen_root_key=[chosen['budget'],chosen['state']], canonical_unique=s==2*model['H'])
    except Limit as exc:
        if pending is not None:
            pending.update(rows=rows_from(new,args,K),completion='PARTIAL')
        raw.update(status='UNKNOWN_RESOURCE_LIMIT',reason=str(exc))
    except Exception as exc:
        if pending is not None:
            pending.update(rows=rows_from(new,args,K),completion='PARTIAL')
        raw.update(status='FAIL_DP',error_type=type(exc).__name__,error=str(exc),traceback=traceback.format_exc())
    return raw


def verify_bellman(model, raw, budget):
    """Verify every Bellman inequality and equality witness; do not reoptimize."""
    need(raw['status']=='DP_CANDIDATE','complete DP receipt')
    K=raw['K']; tables=raw['tables']; byid={t['id']:t for t in tables}
    need(list(byid)==list(range(len(tables))), 'complete ordered table IDs')
    expected_keys={(b,state) for b in range(K+1) for state in STATES}
    finals={}; state_count=0; transition_count=0; init_count=0
    for table in tables:
        budget.check(); v=table['node']; rows=row_map(table)
        need(len(rows)==len(table['rows'])==len(expected_keys) and set(rows)==expected_keys,'finite and INF state coverage')
        state_count+=len(rows)
        for row in rows.values():
            need(row['cost'] is None or type(row['cost']) is int and 0<=row['cost']<=len(model['B']),'integer state cost')
        if table['kind']=='initialize':
            init_count+=1
            need(table['processed_children']==0 and v not in finals,'one initialization per vertex')
            for (b,state),row in rows.items():
                value=None
                if state=='D' and b==1: value=0
                if b==0 and state!='D' and int(state[1])==int(v in model['P']): value=int(state[2])*int(v in model['B'])
                need(row['cost']==value and row['argmin'] is None,'initial state exactness')
        else:
            need(table['kind']=='merge','table kind'); left=byid[table['left_table']]; right=byid[table['right_table']]
            need(table['completion']=='COMPLETE','only complete merges can certify optimum')
            need(left['id']<table['id'] and right['id']<table['id'] and left['node']==v,'acyclic Bellman dependencies')
            j=table['processed_children']; child=model['children'][v][j-1]
            need(j==left['processed_children']+1 and right['node']==child==table['child'],'child order')
            need(right['processed_children']==len(model['children'][child]),'complete child table')
            L=row_map(left); R=row_map(right); expected=set()
            for b1,st1 in expected_keys:
                for b2,st2 in expected_keys:
                    if b1+b2>K: continue
                    if st1=='D': target='D' if st2!='K01' else None
                    elif st2=='D': target=st1
                    elif st1[-1]!=st2[-1]: target=None
                    else: target='K'+str(int(st1[1]=='1' or st2[1]=='1'))+st1[-1]
                    if target is not None: expected.add((b1,st1,b2,st2,b1+b2,target))
            seen=set()
            for tr in table['transitions']:
                budget.check(); transition_count+=1
                key=tuple(tr['left']+tr['right']+tr['target']); need(key in expected and key not in seen,'all legal predecessor transitions once'); seen.add(key)
                a=L[tuple(tr['left'])]['cost']; b=R[tuple(tr['right'])]['cost']; cost=None if a is None or b is None else a+b
                need(tr['cost']==cost,'candidate sum'); value=rows[tuple(tr['target'])]['cost']
                if cost is not None: need(value is not None and value<=cost,'every predecessor lower-bound inequality')
            need(seen==expected,'no missing legal predecessor, including INF')
            for key,row in rows.items():
                if row['cost'] is None: need(row['argmin'] is None,'INF has no argmin')
                else:
                    a=row['argmin']; need(type(a) is int and 0<=a<len(table['transitions']),'finite argmin index')
                    tr=table['transitions'][a]; need(tuple(tr['target'])==key and tr['cost']==row['cost'],'finite equality predecessor')
        finals[v]=table['id']
    need(init_count==len(model['labels']) and set(finals)==set(range(len(model['labels']))),'all vertices represented')
    for v,t in finals.items(): need(byid[t]['processed_children']==len(model['children'][v]),'all children merged')
    need(raw['node_final_tables']=={str(v):t for v,t in finals.items()},'final table mapping')
    root=byid[finals[0]]; R=row_map(root); chosen=tuple(raw['chosen_root_key']); cost=raw['delta']
    need(raw['root_table']==root['id'] and chosen[1] in ('D','K00','K11') and R[chosen]['cost']==cost,'root feasible equality')
    expected_accepted=[dict(budget=r['budget'],state=r['state'],cost=r['cost']) for r in root['rows'] if r['state']!='K01']
    need(raw['accepted_root_states']==expected_accepted,'all closed root alternatives')
    for row in expected_accepted:
        if row['cost'] is not None: need(cost<=row['cost'],'root lower bound against every finite alternative')
    C=set(); trace=[]; stack=[(root['id'],chosen)]
    while stack:
        budget.check(); tid,key=stack.pop(); t=byid[tid]; r=row_map(t)[key]
        trace.append(dict(table=tid,key=list(key),cost=r['cost'],argmin=r['argmin']))
        if t['kind']=='initialize':
            if key[1]=='D': C.add(t['node'])
        else:
            tr=t['transitions'][r['argmin']]; stack.extend([(t['right_table'],tuple(tr['right'])),(t['left_table'],tuple(tr['left']))])
    need(trace==raw['traceback'] and sorted(C)==raw['separator_C'],'traceback replay')
    need(len(C)==chosen[0]<=K,'separator budget')
    budget.counts['Bellman_states_checked']+=state_count; budget.counts['Bellman_transitions_checked']+=transition_count
    return dict(status='PASS_COMPLETE_BELLMAN_CERTIFICATE',tables=len(tables),states=state_count,transitions=transition_count,
                initialization_tables=init_count,delta=cost,new_DP_optimizations_in_checker=0)


def components_objective(model, C, budget):
    """Independent deletion-components implementation, not source-rooted reachability."""
    remaining=set(range(len(model['labels'])))-set(C); components=[]; D=[]
    while remaining:
        start=min(remaining); remaining.remove(start); todo=[start]; comp=[]
        while todo:
            budget.check(); u=todo.pop(); comp.append(u)
            for v in model['adjacency'][u]:
                if v in remaining: remaining.remove(v); todo.append(v)
        active=bool(set(comp)&model['P']); leaves=sorted(set(comp)&model['B'])
        if active: D.extend(leaves)
        components.append(dict(vertices=sorted(comp),contains_P=active,original_leaves=leaves))
    return dict(cost=len(D),D=sorted(D),components=components)


def enumerate_separators(model, s, budget):
    out=dict(status='UNKNOWN',records=[],minimum=None,minimizers=[])
    try:
        for size in range(s):
            for C in itertools.combinations(range(len(model['labels'])),size):
                budget.charge('separator_subsets',3000); value=components_objective(model,C,budget)
                rec=dict(C=list(C),**value); out['records'].append(rec)
                if out['minimum'] is None or value['cost']<out['minimum']:
                    out['minimum']=value['cost']; out['minimizers']=[list(C)]
                elif value['cost']==out['minimum']: out['minimizers'].append(list(C))
        out['status']='COMPLETE_EXHAUSTIVE_SEPARATOR_CHECK'
    except Limit as exc: out.update(status='UNKNOWN_RESOURCE_LIMIT',reason=str(exc))
    except Exception as exc: out.update(status='FAIL_ENUMERATION',error=str(exc),traceback=traceback.format_exc())
    return out


def generate_physical_world(model, specification, budget, raw):
    """Primary physical calculation inserts the prescribed strict sources."""
    raw.update(status='STARTED_PHYSICAL_GENERATION',specification=specification,sources=[],original_vertex_records=[])
    budget.charge('physical_worlds',6)
    q=list(map(Q,specification['allocations'])); y=Q(specification['mean']); H=model['H']
    need(len(q)==2*H and all(a>0 for a in q) and sum(q)==1 and y>0,'strict prescribed normalized world')
    n=len(model['labels']); augmented=[[] for _ in range(n+H)]; selected={e:i for i,e in enumerate(model['support'])}; sources=raw['sources']
    def add(u,v,length): augmented[u].append((v,length)); augmented[v].append((u,length))
    for e,(u,v) in enumerate(model['edges']):
        if e not in selected: add(u,v,Q(1)); continue
        i=selected[e]; a,b=q[2*i:2*i+2]; weight=a+b; x=b/weight; point=n+i
        add(u,point,x); add(point,v,1-x)
        sources.append(dict(edge=e,point_vertex=point,weight=weight,fraction_from_original_first=x,
                            allocation_first=a,allocation_second=b))
    sources.sort(key=lambda z:model['support'].index(z['edge']))
    records=raw['original_vertex_records']
    for r in range(n):
        budget.charge('original_vertex_physical_evaluations',60)
        distances={r:Q(0)}; stack=[r]
        while stack:
            budget.check(); u=stack.pop()
            for v,length in augmented[u]:
                if v not in distances: distances[v]=distances[u]+length; stack.append(v)
        need(len(distances)==n+H,'source insertion remains connected')
        values=[distances[source['point_vertex']] for source in sources]
        potential=y+sum((source['weight']*d for source,d in zip(sources,values)),Q(0))
        records.append(dict(vertex=r,source_distances=values,potential=potential))
    raw.update(status='GENERATED_PRESCRIBED_STRICT_WORLD',mean=y,all_original_lengths='1',expanded_vertex_count=n+H)
    return raw


def verify_physical_pair(model, specifications, worlds, delta, budget):
    n=len(model['labels']); all_dist=[]
    for r in range(n):
        dist={r:Q(0)}; todo=[r]
        while todo:
            budget.check(); u=todo.pop()
            for v in model['adjacency'][u]:
                if v not in dist: dist[v]=dist[u]+1; todo.append(v)
        all_dist.append(dist)
    world_potentials=[]; source_terms=0; leaf_equations=0
    for spec,world in zip(specifications,worlds):
        q=list(map(Q,spec['allocations'])); y=Q(spec['mean']); need(world['mean']==y,'world mean')
        need([x['edge'] for x in world['sources']]==model['support'],'same prescribed support')
        need(len(world['original_vertex_records'])==n,'all original vertices')
        vals=[]
        for r,record in enumerate(world['original_vertex_records']):
            need(record['vertex']==r,'original vertex ordering'); prediction=y
            for i,e in enumerate(model['support']):
                u,v=model['edges'][e]; a,b=q[2*i:2*i+2]; p=a+b; x=b/p
                distance=min(all_dist[r][u]+x,all_dist[r][v]+1-x)
                need(record['source_distances'][i]==distance,'independent original endpoint minimum')
                s=world['sources'][i]
                need(s['weight']==p and s['fraction_from_original_first']==x and s['allocation_first']==a and s['allocation_second']==b,'allocation/source conversion')
                prediction+=p*distance; source_terms+=1
            need(record['potential']==prediction,'independent original-vertex physical equation')
            if r in model['B']: leaf_equations+=1
            vals.append(prediction)
        world_potentials.append(vals)
    difference=[b-a for a,b in zip(*world_potentials)]
    expected_vertices={model['index'][v] for v in specifications[1]['expected_changed_vertices']}
    expected_amount=Q(specifications[1]['expected_difference'])
    need(all(x==(expected_amount if v in expected_vertices else 0) for v,x in enumerate(difference)),'prescribed full potential difference')
    D=[v for v in sorted(model['B']) if difference[v]]
    need(len(D)==delta and set(D)==expected_vertices&model['B'],'minimum changed-leaf count')
    part=set(D[:len(D)//2]); spliced={v:(world_potentials[1][v] if v in part else world_potentials[0][v]) for v in sorted(model['B'])}
    errors=[[v for v in sorted(model['B']) if spliced[v]!=world_potentials[i][v]] for i in range(2)]
    need(set(errors[0])==part and set(errors[1])==set(D)-part,'two error sets partition all differences')
    need(all(len(e)<=(delta+1)//2 for e in errors),'classical correction lower-bound control')
    budget.counts['physical_minpath_terms']+=source_terms; budget.counts['physical_leaf_equations']+=leaf_equations
    return dict(status='PASS_INDEPENDENT_PHYSICAL_PAIR',all_vertex_difference=difference,erased_leaf_indices=D,
        erased_leaf_labels=[model['labels'][v] for v in D],remaining_leaf_reports_equal=True,
        spliced_leaf_reports=[[model['labels'][v],value] for v,value in spliced.items()],error_sets=errors,
        source_minpath_terms=source_terms,original_vertex_equations=2*n,leaf_equations=leaf_equations,
        canonical_erasure_tolerance=delta-1,canonical_corruption_tolerance=(delta-1)//2,
        scope='Six prescribed strict worlds only; no corrupted-report decoder or noise bound.')


def prepare():
    design=read(DESIGN)
    for f,h in design['dependencies_sha256'].items(): need(sha(f)==h,'design dependency '+f)
    root_review=read('results/stage195_root_review.json'); need('PASS' in root_review['status'],'accepted195')
    old={c['id']:c for c in read('data/stage191_prior_graphic_cases.json')['cases']}
    records={c['id']:c for c in read('results/stage191_prior_graphic_validation.json')['records']}
    cases=[]; extra=[]; subset_counts=[]
    for item in design['cases']:
        if 'source' in item:
            archive=records[item['id']]; need(sha(archive['path'])==archive['sha256'],'archive record')
            rec=read(archive['path']); need(sha(rec['raw_path'])==rec['raw_sha256'],'archive raw')
            top=read(rec['raw_path'])['topology']; common=[top['labels'][v] for v in top['P'] if v in top['B']]
            cert=rec['independent_certificate']; paths=[x['vertices'] for x in cert['path_checks']]+[[v] for v in common]
            separator=cert['separator']+common; graph=old[item['id']]['graph_inputs']
            incidence=graph['incidence']; support=graph['support']; extra.extend([archive['path'],rec['raw_path']])
            n=len(top['labels'])
        else:
            incidence=item['edges']; support=item['support']; paths=item['rank_paths']; separator=item['rank_separator']
            n=len(incidence)+1
        s=len(paths); enumerated=item['id']!='public_tree1322'; subsets=sum(math.comb(n,k) for k in range(s)) if enumerated else 0
        subset_counts.append(subsets)
        checker=dict(exhaustive_separator_check=enumerated,expected_separator_subsets=subsets)
        if 'checker_only_expected_delta' in item: checker['expected_delta']=item['checker_only_expected_delta']
        if not enumerated: checker.update(expected_delta=1,method='fullrank_selected_original_leaf_guard')
        physical=None
        if item['id'] in ('internal_H1_quartet','internal_H1_eight_leaf'):
            changed=['u','a','b'] if item['id']=='internal_H1_quartet' else ['u','q1','q2','q1a','q1b','q2a','q2b']
            physical=[dict(name='base',allocations=['1/2','1/2'],mean='1'),
                dict(name='alternative',allocations=['1/4','3/4'],mean='5/4',expected_changed_vertices=changed,expected_difference='1/2')]
        elif item['id']=='internal_H2_caterpillar':
            physical=[dict(name='base',allocations=['1/4']*4,mean='1'),
                dict(name='alternative',allocations=['5/16','1/8','5/16','1/4'],mean='1',expected_changed_vertices=['v','b'],expected_difference='1/8')]
        cases.append(dict(id=item['id'],graph_inputs=dict(incidence=incidence,support=support),
            rank_certificate=dict(paths=paths,separator=separator),checker_only=checker,physical_checker_only=physical))
    need(subset_counts==[0,794,11,7,176,15,3] and sum(subset_counts)==1006,'predeclared metadata subset counts')
    save(DATA,dict(stage=196,status='FROZEN_NOT_EXECUTED',cases=cases,expected_primary_counts=dict(DP_runs=7,
        exhaustive_small_cases=6,separator_subsets=1006,physical_worlds=6,original_vertex_physical_evaluations=60,
        physical_minpath_terms=80,physical_leaf_equations=36),input_boundary='DP sees incidence/support and independently certified rank s only.'))
    deps=list(dict.fromkeys([SCRIPT,DATA,DESIGN,'results/stage195_root_review.json','results/stage195_structure_symbolic_review.json',
        *design['dependencies_sha256'],*extra]))
    for f,h in root_review['dependencies_sha256'].items(): need(sha(f)==h,'195 root dependency '+f)
    protocol=dict(stage=196,status='FROZEN_PREPARED_ONLY_PENDING_NONAUTHOR_REVIEW_AND_ROOT_RUN_APPROVAL',
        dependencies_sha256={f:sha(f) for f in deps},case_order=[c['id'] for c in cases],data_path=DATA,limits=LIMITS,
        expected_primary_counts=read(DATA)['expected_primary_counts'],
        phases=['Verify saved/supplied full leaf rank using disjoint original P-to-B paths and equal-size separator; no new flow.',
            'Exactly one anticipated-activity DP percase with K=s-1. Save all complete/partial tables and traceback before checking.',
            'Verify every state including INF, every compatible predecessor, finite equality argmin and closed root alternatives. No second DP optimization.',
            'For six small inputs enumerate all separators of size<=s-1; save all1006records/allminimizers. Public1322 uses fullrank selected-leaf guard.',
            'On three fixed fixtures only, generate two prescribed strict worlds after DP save. Expanded-source traversal checked by independent original endpoint minpaths; splice the two saved leaf vectors.'],
        state_order=list(STATES),vertex_order='first appearance in original incidence; root0; original adjacency child order',
        tie_rule='First finite minimum in budget/state/predecessor order; equalcost never replaces argmin.',
        receipt='Each table has all4*(K+1) states including INF/null. Each merge saves all compatible budget/state pairs, including INF candidates. No omitted legal predecessors or collapsed INF states. Bellman inequalities and selected equalities prove lower/upper bounds without rerunning optimization.',
        process='One fresh serial child, external180s watchdog,CPU180s,1GiB; no retry. Started and every artifact write-exclusive.',
        failures='Raw DP/enumeration/world outputs saved before comparisons. All7 case slots and all declared physical slots retained. Interrupted work UNKNOWN; invalid certificates FAIL; no unreported completed-outcome filtering.',
        counter_scope='DP_runs and physical_worlds count started computations; completed status is reported separately. Physical vertex evaluation counts are attempted traversals if interrupted. Incomplete merge predecessor prefixes and physical-record prefixes are retained.',
        scope='Topology-only information erasure radius. Canonical thresholds only if verified s=2H. Public data is fixed topology metadata, not biological robustness evidence.',
        preparation_accounting=dict(DP_runs=0,flow_calls=0,rank_calls=0,LP_calls=0,new_physical_worlds=0),
        execution_authorized=False)
    save(PROTOCOL,protocol)
    print(json.dumps(dict(status='PREPARED_ONLY',source_sha256=sha(SCRIPT),data_sha256=sha(DATA),protocol_sha256=sha(PROTOCOL))))


def frozen(token):
    need(token and sha(PROTOCOL)==token,'protocol SHA')
    protocol=read(PROTOCOL)
    for f,h in protocol['dependencies_sha256'].items(): need(sha(f)==h,'frozen dependency '+f)
    return protocol


def worker(token):
    frozen(token); started=read(PREFIX+'_started.json')
    need(str(os.getppid())==os.environ.get('STAGE196_PARENT_PID') and started['parent_pid']==os.getppid() and started['protocol_sha256']==token,'approved parent')
    resource.setrlimit(resource.RLIMIT_AS,(1073741824,1073741824)); resource.setrlimit(resource.RLIMIT_CPU,(180,180))
    budget=Budget(); outcomes=[]
    for case in read(DATA)['cases']:
        cid=case['id']; before=dict(budget.counts); result=dict(id=cid,status='STARTED',physical_records=[])
        model=None; raw=None; rank=None
        try:
            model=topology(case['graph_inputs']['incidence'],case['graph_inputs']['support'],budget)
            rank=verify_rank_certificate(model,case['rank_certificate'],budget); result['rank_certificate']=rank
            raw=solve_dp(model,rank['s'],budget); result['DP_raw']=save(PREFIX+'_'+cid+'_DP_raw.json',raw)
            if raw['status'].startswith('UNKNOWN'): raise Limit(raw.get('reason','DP incomplete'))
            need(raw['status']=='DP_CANDIDATE','DP returned failed certificate')
            result['Bellman_check']=verify_bellman(model,raw,budget)
            objective=components_objective(model,raw['separator_C'],budget)
            need(objective['cost']==raw['delta'] and objective['D']==raw['erased_leaves_D'],'independent component cost/erasure witness')
            result['separator_objective_check']=objective
            if case['checker_only']['exhaustive_separator_check']:
                enumeration=enumerate_separators(model,rank['s'],budget)
                result['enumeration_raw']=save(PREFIX+'_'+cid+'_enumeration_raw.json',enumeration)
                if enumeration['status'].startswith('UNKNOWN'): raise Limit(enumeration.get('reason','enumeration incomplete'))
                need(enumeration['status']=='COMPLETE_EXHAUSTIVE_SEPARATOR_CHECK','enumeration failed')
                need(len(enumeration['records'])==case['checker_only']['expected_separator_subsets'],'all predetermined separators')
                need(enumeration['minimum']==raw['delta'],'independent exhaustive minimum')
            else:
                need(rank['canonical_unique'],'public fullrank guard'); leaves=sorted(model['P']&model['B']); need(leaves,'selected original leaf')
                leaf=leaves[0]; C=model['P']-{leaf}; guard=components_objective(model,C,budget)
                need(len(C)==rank['s']-1 and guard['D']==[leaf] and guard['cost']==raw['delta']==1,'selected-leaf exact threshold')
                result['selected_leaf_guard']=dict(leaf=leaf,C=sorted(C),**guard)
            if 'expected_delta' in case['checker_only']: need(raw['delta']==case['checker_only']['expected_delta'],'predeclared analytic delta')
            result['thresholds']=dict(information_erasure_tolerance=raw['delta']-1,
                canonical_erasure_tolerance=raw['delta']-1 if rank['canonical_unique'] else None,
                canonical_corruption_tolerance=(raw['delta']-1)//2 if rank['canonical_unique'] else None,
                canonical_unique_before_erasure=rank['canonical_unique'],full_leaf_information_rank=rank['s'])
            result['status']='PASS'
        except Limit as exc: result.update(status='UNKNOWN',reason=str(exc))
        except Exception as exc: result.update(status='FAIL',error_type=type(exc).__name__,error=str(exc),traceback=traceback.format_exc())
        specs=case['physical_checker_only']
        if specs is not None:
            worlds=[]
            for i,spec in enumerate(specs):
                world={}; physical_record={}
                try:
                    need(model is not None,'physical geometry unavailable')
                    generate_physical_world(model,spec,budget,world)
                    worlds.append(world); physical_record['status']='SAVED'
                except Limit as exc:
                    world.update(status='UNKNOWN_RESOURCE_LIMIT',reason=str(exc))
                    physical_record.update(status='UNKNOWN',reason=str(exc))
                    if result['status']!='FAIL': result['status']='UNKNOWN'
                except Exception as exc:
                    world.update(status='FAIL_PHYSICAL_GENERATION',error=str(exc),traceback=traceback.format_exc())
                    physical_record.update(status='FAIL',error=str(exc)); result['status']='FAIL'
                physical_record.update(save(PREFIX+'_'+cid+'_world_'+str(i)+'.json',world))
                result['physical_records'].append(physical_record)
            if len(worlds)==2:
                try:
                    need(raw is not None and raw['status']=='DP_CANDIDATE' and rank['canonical_unique'],'certified canonical threshold required for pair')
                    result['physical_pair_check']=verify_physical_pair(model,specs,worlds,raw['delta'],budget)
                except Limit as exc:
                    result['physical_pair_check']=dict(status='UNKNOWN',reason=str(exc))
                    if result['status']!='FAIL': result['status']='UNKNOWN'
                except Exception as exc:
                    result['physical_pair_check']=dict(status='FAIL',error=str(exc),traceback=traceback.format_exc()); result['status']='FAIL'
        result['counters']={k:v-before.get(k,0) for k,v in budget.counts.items()}
        artifact=save(PREFIX+'_'+cid+'.json',result); outcomes.append(dict(id=cid,status=result['status'],**artifact))
        event(dict(event='CASE_COMPLETED',id=cid,status=result['status'],counters=result['counters']))
    status='FAIL' if any(r['status']=='FAIL' for r in outcomes) else 'PARTIAL_WITH_UNKNOWN' if any(r['status']!='PASS' for r in outcomes) else 'PASS'
    if status=='PASS':
        for key,value in read(DATA)['expected_primary_counts'].items():
            if key=='exhaustive_small_cases': continue
            need(budget.counts[key]==value,'final fixed count '+key)
    save(PREFIX+'_worker.json',dict(stage=196,status=status,records=outcomes,counters=dict(budget.counts),
        flow_calls=0,rank_calls=0,LP_calls=0,retries=0,all_fixed_cases_retained=len(outcomes)==7))


def run_approved(token):
    protocol=frozen(token)
    save(PREFIX+'_started.json',dict(stage=196,status='STARTED_NOT_A_RESULT',protocol_sha256=token,parent_pid=os.getpid(),retries=0))
    logs=[PREFIX+'_stdout.log',PREFIX+'_stderr.log']; begin=time.monotonic(); timeout=False
    with (ROOT/logs[0]).open('x') as out,(ROOT/logs[1]).open('x') as err:
        child=subprocess.Popen([sys.executable,'-B',str(ROOT/SCRIPT),'--worker','--protocol-sha',token],cwd=ROOT,stdout=out,stderr=err,
            env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','STAGE196_PARENT_PID':str(os.getpid())})
        try: code=child.wait(timeout=180)
        except subprocess.TimeoutExpired: timeout=True; child.kill(); code=child.wait()
    errors=[]; records=[]; summary=None
    for cid in protocol['case_order']:
        path=PREFIX+'_'+cid+'.json'
        try: records.append(dict(id=cid,record=read(path),path=path,sha256=sha(path)) if (ROOT/path).exists() else dict(id=cid,status='UNKNOWN',reason='no completed record'))
        except Exception as exc: errors.append(str(exc)); records.append(dict(id=cid,status='UNKNOWN',reason='unreadable record'))
    if (ROOT/(PREFIX+'_worker.json')).exists():
        try: summary=read(PREFIX+'_worker.json')
        except Exception as exc: errors.append(str(exc))
    status=summary['status'] if code==0 and not timeout and not errors and summary else 'UNKNOWN_OR_FAILURE'
    output=dict(stage=196,status=status,protocol_sha256=token,dependencies_sha256=protocol['dependencies_sha256'],
        records=records,worker_summary=summary,hard_timeout=timeout,returncode=code,elapsed_seconds=time.monotonic()-begin,
        parsing_errors=errors,raw_logs=logs,unexposed_inflight_work_possible=timeout or code!=0,retries=0,
        flow_calls=0,rank_calls=0,LP_calls=0)
    save(PREFIX+'.json',output); print(json.dumps(dict(status=status,result_sha256=sha(PREFIX+'.json'))))


if __name__=='__main__':
    parser=argparse.ArgumentParser(); actions=parser.add_mutually_exclusive_group(required=True)
    actions.add_argument('--prepare',action='store_true'); actions.add_argument('--worker',action='store_true'); actions.add_argument('--run-approved',action='store_true')
    parser.add_argument('--protocol-sha'); args=parser.parse_args()
    if args.prepare: prepare()
    elif args.worker: worker(args.protocol_sha)
    else: run_approved(args.protocol_sha)
