#!/usr/bin/env python3
"""Certified rational fixed-support bounded-error PDD feasibility.

SciPy proposes bases; Fraction primal/dual identities certify conclusions.
This is a standard small LP formulation, not a novel optimization algorithm.
No frozen experimental module is imported. main writes only Stage25 output.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json, time, os, warnings
ROOT=Path(__file__).resolve().parents[1]
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[_key]='1'


def geometry(edges):
    adj={v:{} for e in edges for v in e[:2]}
    for u,v,l in edges:
        l=F(l)
        if l<=0 or u==v or v in adj[u]:raise ValueError('Invalid positive tree edge')
        adj[u][v]=adj[v][u]=l
    if len(edges)!=len(adj)-1:raise ValueError('Input is not a tree')
    D={}
    for root in adj:
        dd={root:F(0)};todo=[root]
        for v in todo:
            for u,l in adj[v].items():
                if u not in dd:dd[u]=dd[v]+l;todo.append(u)
        if len(dd)!=len(adj):raise ValueError('Disconnected tree')
        D[root]=dd
    return adj,D


def square_solve(A,b):
    n=len(b);M=[[F(v) for v in row]+[F(v)] for row,v in zip(A,b)]
    if len(A)!=n or any(len(row)!=n for row in A):raise ValueError('Not square')
    for j in range(n):
        p=next((i for i in range(j,n) if M[i][j]),None)
        if p is None:return None
        M[j],M[p]=M[p],M[j];a=M[j][j];M[j]=[v/a for v in M[j]]
        for i in range(n):
            if i!=j and M[i][j]:
                a=M[i][j];M[i]=[v-a*w for v,w in zip(M[i],M[j])]
    return [row[-1] for row in M]


def rank(rows):
    if not rows:return 0
    M=[[F(v) for v in row] for row in rows];i=0
    for j in range(len(M[0])):
        p=next((p for p in range(i,len(M)) if M[p][j]),None)
        if p is None:continue
        M[i],M[p]=M[p],M[i];a=M[i][j];M[i]=[v/a for v in M[i]]
        for p in range(i+1,len(M)):
            if M[p][j]:
                a=M[p][j];M[p]=[v-a*w for v,w in zip(M[p],M[i])]
        i+=1
        if i==len(M):break
    return i


def dot(a,b):return sum((x*y for x,y in zip(a,b)),F(0))


def verify_lp_certificate(A,b,c,cert):
    x=[F(v) for v in cert['primal']]
    lam=[F(0)]*len(A)
    for i,v in cert['dual_nonzero']:lam[i]=F(v)
    assert all(v>=0 for v in lam)
    assert all(dot(row,x)<=rhs for row,rhs in zip(A,b))
    assert all(sum((lam[i]*A[i][j] for i in range(len(A))),F(0))==c[j] for j in range(len(c)))
    assert dot(b,lam)==dot(c,x)==F(cert['objective'])
    return True


def certified_lp(A,b,c,max_exhaustive_bases=100000,proposal=True):
    """Max c*x subject to A*x<=b. Returns only exact certified optimum or unknown.

    For our two models, whenever feasible, a finite optimum has a vertex.
    Enumeration is a termination fallback, not used as an infeasibility oracle.
    """
    A=[[F(v) for v in row] for row in A];b=list(map(F,b));c=list(map(F,c))
    n=len(c);tested=set();proposals=[];proposal_status=None
    def test(I,method):
        I=tuple(sorted(I))
        if I in tested:return None
        tested.add(I);B=[A[i] for i in I]
        x=square_solve(B,[b[i] for i in I])
        if x is None or any(dot(row,x)>v for row,v in zip(A,b)):return None
        lam=square_solve([[B[i][j] for i in range(n)] for j in range(n)],c)
        if lam is None or any(v<0 for v in lam):return None
        out={'status':'optimal','primal':x,'objective':dot(c,x),'basis':list(I),
             'dual_nonzero':[[i,v] for i,v in zip(I,lam) if v],
             'certificate_method':method,'bases_tested':len(tested)}
        verify_lp_certificate(A,b,c,out);return out
    if proposal:
        try:
            from scipy.optimize import linprog
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                res=linprog([-float(v) for v in c],A_ub=[[float(v) for v in row] for row in A],
                    b_ub=[float(v) for v in b],bounds=[(None,None)]*n,
                    method='highs',options={'threads':1,'parallel':False,'presolve':False})
            proposal_status=int(res.status)
            if res.success:
                order=sorted(range(len(A)),key=lambda i:(abs(float(res.ineqlin.residual[i])),
                    -abs(float(res.ineqlin.marginals[i]))))
                active=[i for i in order if abs(float(res.ineqlin.residual[i]))<=1e-7*(1+abs(float(b[i])))]
                positive=[i for i in active if -float(res.ineqlin.marginals[i])>1e-9]
                # One greedily independent basis, prioritizing proposed dual support.
                chosen=[]
                for i in positive+[i for i in active if i not in positive]:
                    if rank([A[j] for j in chosen+[i]])>len(chosen):chosen.append(i)
                    if len(chosen)==n:break
                if len(chosen)==n:
                    cert=test(chosen,'numeric_proposal_exact_primal_dual')
                    if cert:return cert
                if len(active)>=n and comb(len(active),n)<=max_exhaustive_bases:
                    for I in combinations(active,n):
                        cert=test(I,'active_set_enumeration_exact_primal_dual')
                        if cert:return cert
        except (ImportError,ValueError,OverflowError) as exc:
            proposal_status=type(exc).__name__
    total=comb(len(A),n)
    if total<=max_exhaustive_bases:
        for I in combinations(range(len(A)),n):
            cert=test(I,'exhaustive_exact_primal_dual')
            if cert:return cert
    return {'status':'unknown','reason':'No exact primal-dual certificate reconstructed within basis limit',
            'bases_tested':len(tested),'all_bases_count':total,'proposal_status':proposal_status,
            'max_exhaustive_bases':max_exhaustive_bases}


def build_model(edges,observed,epsilon,support,phase='margin'):
    """Uses unique endpoints; m_last=1-sum(other masses), then y, then t/z."""
    adj,D=geometry(edges);eps=F(epsilon);obs={int(r):F(v) for r,v in observed.items()}
    support=sorted(support)
    if not 1<=len(support)<=2 or len(set(support))!=len(support):raise ValueError('One or two distinct edges required')
    if not obs or not set(obs)<=set(adj) or eps<0:raise ValueError('Invalid observed subset or epsilon')
    if any(i<0 or i>=len(edges) for i in support):raise ValueError('Invalid support')
    U=sorted({v for i in support for v in edges[i][:2]});u=len(U);dim=u+1
    A=[];b=[];labels=[]
    def add(row,rhs,label):A.append(list(map(F,row)));b.append(F(rhs));labels.append(label)
    for i in range(u-1):
        row=[F(0)]*dim;row[i]=-1
        if phase=='margin':row[-1]=1
        add(row,0,'endpoint_mass_'+str(U[i]))
    row=[F(1)]*(u-1)+[F(0),F(1) if phase=='margin' else F(0)]
    add(row,1,'endpoint_mass_'+str(U[-1]))
    row=[F(0)]*dim;row[-2]=-1
    if phase=='margin':row[-1]=1
    add(row,0,'pendant')
    row=[F(0)]*dim;row[-1]=-1;add(row,0,'aux_nonnegative')
    if phase=='margin':
        row=[F(0)]*dim;row[-1]=1;add(row,1,'margin_cap')
    for r,d in sorted(obs.items()):
        coeff=[D[r][v]-D[r][U[-1]] for v in U[:-1]]+[F(1),F(0)]
        upper=coeff[:];lower=[-x for x in coeff]
        if phase=='slack':upper[-1]=lower[-1]=-1
        add(upper,d+eps-D[r][U[-1]],'observation_upper_'+str(r))
        add(lower,D[r][U[-1]]-d+eps,'observation_lower_'+str(r))
    c=[F(0)]*dim;c[-1]=1 if phase=='margin' else -1
    return {'A':A,'b':b,'c':c,'labels':labels,'vertices':U,'phase':phase,
            'selected_support':support,'observed':obs,'epsilon':eps}


def witness_from_primal(edges,model,primal):
    U=model['vertices'];u=len(U);x=list(map(F,primal))
    m=dict(zip(U,x[:u-1]+[1-sum(x[:u-1])]))
    y=x[-2];support=model['selected_support'];adj,D=geometry(edges)
    degree={v:sum(v in edges[i][:2] for i in support) for v in U}
    sources=[]
    for i in support:
        a,b,_=edges[i];aa=m[a]/degree[a];bb=m[b]/degree[b]
        p=aa+bb
        sources.append({'edge':i,'p':p,'x':bb/p if p else None,
                        'endpoint_contributions':[aa,bb]})
    pred={r:y+sum((mass*D[r][v] for v,mass in m.items()),F(0)) for r in model['observed']}
    strict=min(m.values())>0 and y>0
    return {'endpoint_mass':m,'pendant':y,'sources':sources,'predicted':pred,
            'strict':strict,'max_observation_error':max(abs(pred[r]-model['observed'][r]) for r in pred)}


def solve_fixed_support(edges,observed,epsilon,support,max_exhaustive_bases=100000,proposal=True):
    """Fixed 1/2-edge strict feasibility. Observations may be a subset of vertices.

    Unknown is an explicit numerical/certificate-reconstruction outcome, never infeasible.
    """
    model=build_model(edges,observed,epsilon,support,'margin')
    opt=certified_lp(model['A'],model['b'],model['c'],max_exhaustive_bases,proposal)
    base={'selected_support':model['selected_support'],'unique_endpoints':model['vertices'],
          'epsilon':F(epsilon),'observation_count':len(observed),
          'independent_variables':len(model['c']),'calibration':'user_supplied_deterministic_budget'}
    if opt['status']=='optimal':
        w=witness_from_primal(edges,model,opt['primal'])
        assert w['max_observation_error']<=F(epsilon)
        status='strict_feasible' if opt['objective']>0 else 'closure_only'
        assert w['strict']==(status=='strict_feasible')
        return dict(base,status=status,margin_certificate=opt,margin_labels=model['labels'],witness=w)
    slack=build_model(edges,observed,epsilon,support,'slack')
    sc=certified_lp(slack['A'],slack['b'],slack['c'],max_exhaustive_bases,proposal)
    if sc['status']=='optimal' and sc['objective']<0:
        return dict(base,status='closure_infeasible',slack_certificate=sc,slack_labels=slack['labels'],
                    minimum_extra_uniform_error=-sc['objective'],margin_attempt=opt)
    return dict(base,status='unknown',margin_attempt=opt,slack_attempt=sc,
                reason='Exact certificate reconstruction did not settle strict feasibility')


def solve_all_supports(edges,observed,epsilon,**kwargs):
    results=[solve_fixed_support(edges,observed,epsilon,list(S),**kwargs)
             for k in (1,2) for S in combinations(range(len(edges)),k)]
    feasible=[r['selected_support'] for r in results if r['status']=='strict_feasible']
    unknown=[r['selected_support'] for r in results if r['status']=='unknown']
    return {'status':('strict_feasible' if feasible else 'unknown' if unknown else 'strict_infeasible'),
            'strict_feasible_supports':feasible,'unknown_supports':unknown,
            'support_identified':len(feasible)==1 and not unknown,'results':results,
            'scope':'All original edge subsets of size one or two; no greater-source competitors'}


def serialize(obj):
    if isinstance(obj,F):return str(obj)
    raise TypeError(type(obj).__name__)


def main():
    start=time.perf_counter();cases=[]
    edge=[(0,1,F(2))]
    fixtures=[
      ('exact_strict_single',edge,{0:F(11,12),1:F(19,12)},F(0),[0],'strict_feasible'),
      ('zero_pendant_boundary',edge,{0:F(1),1:F(1)},F(0),[0],'closure_only'),
      ('zero_endpoint_boundary',edge,{0:F(1),1:F(3)},F(0),[0],'closure_only'),
      ('zero_endpoint_positive_band',edge,{0:F(1),1:F(3)},F(1,100),[0],'strict_feasible'),
      ('negative_observations',edge,{0:F(-1),1:F(-1)},F(0),[0],'closure_infeasible'),
      ('lipschitz_violation',edge,{0:F(1),1:F(4)},F(0),[0],'closure_infeasible'),
      ('exact_error_threshold',edge,{0:F(1),1:F(4)},F(1,2),[0],'closure_only'),
      ('above_error_threshold',edge,{0:F(1),1:F(4)},F(5001,10000),[0],'strict_feasible'),
      ('tiny_positive_margin',edge,{0:1+F(2,10**12),1:3-F(2,10**12)},F(0),[0],'strict_feasible'),
    ]
    for name,edges,obs,eps,support,expected in fixtures:
        sol=solve_fixed_support(edges,obs,eps,support)
        assert sol['status']==expected,(name,sol)
        # An independent exact basis enumeration makes no numeric proposal.
        exact=solve_fixed_support(edges,obs,eps,support,proposal=False)
        assert exact['status']==sol['status']
        if 'margin_certificate' in sol:
            assert exact['margin_certificate']['objective']==sol['margin_certificate']['objective']
        else:assert exact['minimum_extra_uniform_error']==sol['minimum_extra_uniform_error']
        cases.append({'name':name,'edges':edges,'observed':obs,'epsilon':eps,
                      'support':support,'expected':expected,'solution':sol,'exhaustive_solution':exact})
    frozen=ROOT/'results/stage24_robust_two_source.json';raw=json.loads(frozen.read_text())
    inputs=[]
    for rec in raw['stable_noise_cases']:
        inputs.append(('stage24_'+rec['name'],rec,rec['truth']))
    rec=raw['negative_nominal_pendant'];inputs.append(('negative_nominal_pendant',rec,rec['truth']))
    # Three distinct failure mechanisms, one exact midpoint box each, both truths retained.
    for rec in raw['instability_sequence_cases'][:3]:
        for i,w in enumerate(rec['worlds']):inputs.append((rec['kind']+'_'+str(i),rec,w))
    recovered=[]
    for name,r,w in inputs:
        edges=[(a,b,F(l)) for a,b,l in r['edges']];obs={int(k):F(v) for k,v in r['observed'].items()}
        eps=F(r['epsilon']);support=sorted(s['edge'] for s in w['sources'])
        sol=solve_fixed_support(edges,obs,eps,support);assert sol['status']=='strict_feasible',(name,sol)
        recovered.append({'name':name,'edges':edges,'observed':obs,'epsilon':eps,'support':support,'solution':sol})
    # Actual all-support conclusions, including nonempty unique and ambiguous boxes.
    chosen=[raw['stable_noise_cases'][0],raw['instability_sequence_cases'][0],
            raw['instability_sequence_cases'][1],raw['negative_nominal_pendant']]
    allcases=[]
    for i,r in enumerate(chosen):
        edges=[(a,b,F(l)) for a,b,l in r['edges']];obs={int(k):F(v) for k,v in r['observed'].items()};eps=F(r['epsilon'])
        out=solve_all_supports(edges,obs,eps)
        assert out['status']=='strict_feasible' and not out['unknown_supports']
        if i in (0,3):assert out['support_identified']
        else:assert len(out['strict_feasible_supports'])>=2
        allcases.append({'name':['stable_unique','adjacent_ambiguous','single_ambiguous','negative_nominal_y_unique'][i],
                         'edges':edges,'observed':obs,'epsilon':eps,'solution':out})
    # Missing observations are permitted by this LP but invalidate all-leaf envelope claims.
    r=chosen[0];edges=[(a,b,F(l)) for a,b,l in r['edges']];obs={0:F(10)}
    out=solve_all_supports(edges,obs,F(0));assert len(out['strict_feasible_supports'])==15
    allcases.append({'name':'one_observed_leaf_all_supports_feasible','edges':edges,'observed':obs,'epsilon':F(0),'solution':out})
    # Exact >2-source-compatible prediction rejects every <=2 strict support.
    edges=[(3,0,F(1)),(3,1,F(1)),(3,2,F(1))];obs={0:F(2),1:F(2),2:F(2)}
    out=solve_all_supports(edges,obs,F(0));assert out['status']=='strict_infeasible'
    assert all(r['status']=='closure_only' for r in out['results'])
    allcases.append({'name':'strict_three_source_only_vs_at_most_two','edges':edges,'observed':obs,'epsilon':F(0),'solution':out,
                     'known_three_source_witness':{'sources':[{'edge':i,'p':F(1,3),'x':F(1,2)} for i in range(3)],'pendant':F(5,6)}})
    counts={'analytic_cases':len(cases),'recovered_stage24_inputs':len(recovered),'all_support_datasets':len(allcases),
            'all_support_lp_cases':sum(len(r['solution']['results']) for r in allcases)}
    result={'passed':True,'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
            'created_utc':datetime.now(timezone.utc).isoformat(),
            'frozen_input':{'path':str(frozen.relative_to(ROOT)),'sha256':sha256(frozen.read_bytes()).hexdigest()},
            'counts':counts,'analytic_cases':cases,'stage24_recovered_cases':recovered,'all_support_cases':allcases,
            'elapsed_seconds':time.perf_counter()-start,
            'scope':'Standard LP completion with exact rational certificates; deterministic user-supplied bands, no statistical or novelty claim'}
    path=ROOT/'results/stage25_bounded_error.json';path.write_text(json.dumps(result,indent=2,default=serialize)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('analytic_cases','stage24_recovered_cases','all_support_cases')},indent=2))
    return result

if __name__=='__main__':main()
