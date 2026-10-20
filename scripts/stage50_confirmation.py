#!/usr/bin/env python3
"""Frozen confirmation: independent closed lower bound plus strict witness.

Closed NO implies strict NO. If the closed minimum is H and a separately
replayed strict H witness exists, the strict minimum is H too. For a NO input,
a strict H+1 promise witness certifies its full strict minimum as H+1.
"""
from datetime import datetime,timezone
from fractions import Fraction as F
from pathlib import Path
from math import lcm
import copy
import hashlib
import json
import subprocess
import sys
import time

from stage50_binary_reduction import generate_case,promised_witness,witness_for_solution,sum_solutions

ROOT=Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replay(case,sources,strict=True):
    """Independent original paths; no reduction metadata is read."""
    E=[(u,v,F(ell)) for u,v,ell in case['edges']]
    adj={}
    for u,v,ell in E:
        assert ell>0 and u!=v and v not in adj.setdefault(u,{})
        adj[u][v]=ell;adj.setdefault(v,{})[u]=ell
    assert len(E)==len(adj)-1
    assert all(len(ns) in (1,3) for ns in adj.values())
    assert sum(F(s['p']) for s in sources)==1
    assert len({s['edge'] for s in sources})==len(sources)
    for r,expected in case['observed'].items():
        assert len(adj[r])==1
        ds={r:F(0)};todo=[r]
        for u in todo:
            for v,ell in adj[u].items():
                if v not in ds:ds[v]=ds[u]+ell;todo.append(v)
        assert len(ds)==len(adj)
        value=F(case['pendant'])
        for s in sources:
            p,x=F(s['p']),F(s['x'])
            assert p>0 and ((0<x<1) if strict else (0<=x<=1))
            u,v,ell=E[s['edge']]
            value+=p*min(ds[u]+x*ell,ds[v]+(1-x)*ell)
        assert value==F(expected)
    return True


def integer_scaled(case):
    scale=1
    for ell in [F(e[2]) for e in case['edges']]+[F(d)-1 for d in case['observed'].values()]:
        scale=lcm(scale,ell.denominator)
    changed=copy.deepcopy(case)
    changed['name']+='_integer_scaled'
    changed['edges']=[[u,v,str(F(ell)*scale)] for u,v,ell in case['edges']]
    changed['observed']={r:str(1+scale*(F(d)-1)) for r,d in case['observed'].items()}
    changed['rescaling']=dict(scale=str(scale),rule='Scale all edges by L, set d_new=1+L*(d_old-1), keep pendant1. Use unscaled construction for identical weights and edge fractions.')
    assert all(F(e[2]).denominator==1 for e in changed['edges'])
    assert all(F(d).denominator==1 for d in changed['observed'].values())
    return changed


def main():
    pp=ROOT/'results/stage50_confirmation_protocol.json'
    rp=ROOT/'results/stage50_confirmation.json'
    if pp.exists() or rp.exists():raise RuntimeError('Refusing to replace frozen confirmation')
    M=10**18
    specs=[
      ('two_a_yes',[[2,7],[4,10]],11,'caterpillar'),
      ('two_a_no',[[2,7],[4,10]],13,'balanced'),
      ('two_b_yes',[[5,8],[6,12]],17,'balanced'),
      ('two_b_no',[[5,8],[6,12]],18,'caterpillar'),
      ('three_a_yes',[[2],[3],[4,9]],14,'caterpillar'),
      ('three_a_no',[[2],[3],[4,9]],11,'balanced'),
      ('three_b_yes',[[2,6],[5],[7]],18,'balanced'),
      ('three_b_no',[[2,6],[5],[7]],16,'caterpillar'),
      ('four_yes',[[2],[4],[6],[9,13]],25,'balanced'),
      ('four_no',[[2],[4],[6],[9,13]],23,'caterpillar'),
      ('large_yes',[[M+2,M+6],[2*M+3,2*M+9]],3*M+11,'caterpillar'),
      ('large_no',[[M+2,M+6],[2*M+3,2*M+9]],3*M+10,'balanced'),
    ]
    bases=[generate_case(g,t,n,s) for n,g,t,s in specs]
    cases=bases+[integer_scaled(bases[i]) for i in (0,1)]
    base_indices=list(range(len(bases)))+[0,1]
    names=['stage50_confirmation.py','stage50_binary_reduction.py','stage50_closed_screened.py','stage50_closed_oracle.py']
    hashes={'scripts/'+n:sha(ROOT/'scripts'/n) for n in names}
    protocol=dict(created_utc=datetime.now(timezone.utc).isoformat(),base_cases=bases,cases=cases,
                  base_indices=base_indices,source_sha256=hashes,hard_seconds_per_case=30,
                  design='12 new base inputs plus 2 deterministic integer rescalings; no case replacement or timeout filtering.',
                  inference='Independent closed source minimum plus independently replayed strict witnesses bounds the strict source minimum.',
                  limitations='Finite correctness checks, not proof of hardness or originality; no performance comparison.')
    pp.write_text(json.dumps(protocol,indent=2)+'\n')
    worker='import json,sys; from stage50_closed_screened import exact_closed_minimum,serial; print(json.dumps(serial(exact_closed_minimum(json.load(sys.stdin)))))'
    rows=[]
    for case,index in zip(cases,base_indices):
        base=bases[index];info=base['reduction']
        solutions=sum_solutions(info['groups'],info['target'])
        expected=case['cap'] if solutions else None
        promise=promised_witness(base)
        assert replay(case,promise['sources'])
        strict=None
        if solutions:
            strict=witness_for_solution(base,solutions[0])
            assert len(strict)==case['cap'] and replay(case,strict)
        start=time.monotonic()
        try:
            proc=subprocess.run([sys.executable,'-c',worker],input=json.dumps(case),capture_output=True,text=True,cwd=ROOT/'scripts',timeout=30)
            oracle=dict(status='complete',answer=json.loads(proc.stdout)) if proc.returncode==0 else dict(status='error',stderr=proc.stderr)
        except subprocess.TimeoutExpired:
            oracle=dict(status='unknown',reason='hard_wall_timeout')
        agreement=None;minimum=None
        if oracle['status']=='complete':
            answer=oracle['answer'];agreement=answer['minimum']==expected
            if answer['minimum'] is not None:assert replay(case,answer['witness'],strict=False)
            if agreement:
                minimum=case['cap'] if solutions else case['cap']+1
                if not solutions:assert promise['source_count']==case['cap']+1
        rows.append(dict(name=case['name'],base_index=index,group_solutions=solutions,expected_closed_minimum_within_cap=expected,
                         oracle=oracle,agreement=agreement,certified_strict_minimum=minimum,
                         promise=promise,promise_replayed=True,strict_yes_witness=strict,oracle_seconds=time.monotonic()-start))
        rp.write_text(json.dumps(dict(status='running',rows=rows),indent=2)+'\n')
        print(case['name'],oracle['status'],agreement,flush=True)
    assert all(sha(ROOT/n)==h for n,h in hashes.items())
    summary=dict(base_inputs=12,rescalings=2,total=len(rows),complete=sum(r['oracle']['status']=='complete' for r in rows),
                 unknown=sum(r['oracle']['status']=='unknown' for r in rows),errors=sum(r['oracle']['status']=='error' for r in rows),
                 agreement=sum(r['agreement'] is True for r in rows),disagreement=sum(r['agreement'] is False for r in rows),
                 yes=sum(bool(r['group_solutions']) for r in rows),no=sum(not r['group_solutions'] for r in rows),
                 strict_promises_replayed=len(rows),strict_yes_replayed=sum(r['strict_yes_witness'] is not None for r in rows))
    result=dict(status='PASS' if summary['agreement']==len(rows) else 'INCOMPLETE_OR_FAILURE',protocol_sha256=sha(pp),
                source_hashes_unchanged=True,summary=summary,rows=rows)
    rp.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(summary))


if __name__=='__main__':main()
