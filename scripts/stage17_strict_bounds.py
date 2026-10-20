#!/usr/bin/env python3
"""Exact fixed-pendant sharp lower controls and finite-point strict compression.
Standard library only. No imports from previous experimental helpers.
"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,time
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rr(A):
 A=[list(map(F,r)) for r in A];r=0;p=[]
 for c in range(len(A[0]) if A else 0):
  j=next((j for j in range(r,len(A)) if A[j][c]),None)
  if j is None:continue
  A[r],A[j]=A[j],A[r];t=A[r][c];A[r]=[v/t for v in A[r]]
  for j in range(len(A)):
   if j!=r and A[j][c]:
    t=A[j][c];A[j]=[v-t*w for v,w in zip(A[j],A[r])]
  p.append(c);r+=1
  if r==len(A):break
 return A,p
def solve(A,b):
 n=len(A[0]);R,p=rr([row+[v] for row,v in zip(A,b)])
 if n in p or len(p)!=n:return None
 x=[F(0)]*n
 for i,j in enumerate(p):x[j]=R[i][-1]
 assert all(sum(u*v for u,v in zip(row,x))==t for row,t in zip(A,b))
 return x

def geom(E):
 V=sorted({v for u,v,l in E for v in [u,v]});adj={v:{} for v in V}
 for u,v,l in E:adj[u][v]=F(l);adj[v][u]=F(l)
 leaves=[v for v in V if len(adj[v])==1];D={}
 for root in V:
  d={root:F(0)};todo=[root]
  for v in todo:
   for w,l in adj[v].items():
    if w not in d:d[w]=d[v]+l;todo.append(w)
  D[root]=d
 return V,adj,leaves,D

def bases(A,b):
 rank=len(rr(A)[1]);N=len(A[0]);out={}
 for cols in combinations(range(N),rank):
  p=solve([[row[i] for i in cols] for row in A],b)
  if p is None or min(p)<0:continue
  x=[F(0)]*N
  for i,v in zip(cols,p):x[i]=v
  out[tuple(x)]=x
 return list(out.values())

def analyze(name,E):
 V,adj,R,D=geom(E);n=len(R);I=[v for v in V if v not in R];assert n>=3
 DL=[[D[r][s] for s in R] for r in R]
 lam=solve(DL,[F(1)]*n);assert lam is not None and sum(lam)>0
 g={v:sum(lam[i]*D[r][v] for i,r in enumerate(R)) for v in V}
 assert all(g[r]==1 for r in R) and all(g[v]<1 for v in I)
 BL=[[D[r][s]-D[R[0]][s] for s in R] for r in R[1:]]+[[F(1)]*n]
 H={v:solve(BL,[D[r][v]-D[R[0]][v] for r in R[1:]]+[F(1)]) for v in V}
 assert all(sum(H[v])==1 and min(H[v])>=0 for v in V)
 assert all(min(H[v])>0 for v in I)
 for v in I:
  assert all(sum((H[v][i]-H[u][i])/l for u,l in adj[v].items())==0 for i in range(n))
 w=[F(i+1,n*(n+1)//2) for i in range(n)];wmin=min(w);kappa=min(1-g[v] for v in I)
 parent=[next(iter(adj[r])) for r in R]
 loss=sum(w[i]*(1-g[parent[i]]) for i in range(n));theta=wmin/(2*(wmin+loss/kappa));assert 0<theta<1
 fixed_y=F(3,5);a={v:F(0) for v in V}
 for i,r in enumerate(R):a[r]+=(1-theta)*w[i];a[parent[i]]+=theta*w[i]
 d=[sum(D[r][v]*a[v] for v in V)+fixed_y for r in R]
 P=[[D[r][v] for v in V] for r in R]+[[F(1)]*len(V)];target=[v-fixed_y for v in d]+[F(1)]
 vertices=bases(P,target);assert vertices
 quantitative=(1-theta)*wmin-theta*loss/kappa;assert quantitative==wmin/2
 for x in vertices:
  assert sum((1-g[v])*x[j] for j,v in enumerate(V))==theta*loss
  assert sum(x[V.index(v)] for v in I)<=theta*loss/kappa
  assert all(x[V.index(r)]>=quantitative for r in R)
 # Fixed original strict points: one point per edge. Compress weights, retaining positive y.
 points=[]
 for i,(u,v,l) in enumerate(E):
  t=F(2+(i%3),7)
  points.append(dict(edge=[u,v],position=t,distances=[(1-t)*D[r][u]+t*D[r][v] for r in R]))
 m=len(points);weights=[F((i+1)**2,sum((j+1)**2 for j in range(m))) for i in range(m)]
 old_y=F(1);obs=[sum(weights[i]*points[i]['distances'][j] for i in range(m))+old_y for j in range(n)]
 B=[[p['distances'][j]-p['distances'][0] for p in points] for j in range(1,n)]+[[F(1)]*m]
 b=[obs[j]-obs[0] for j in range(1,n)]+[F(1)]
 bv=bases(B,b);assert bv
 costs=[sum(x[i]*points[i]['distances'][0] for i in range(m)) for x in bv]
 best=bv[min(range(len(bv)),key=lambda i:costs[i])];new_y=obs[0]-min(costs)
 assert new_y>=old_y and sum(v>0 for v in best)<=n
 assert all(sum(best[i]*points[i]['distances'][j] for i in range(m))+new_y==obs[j] for j in range(n))
 # Preserve y instead: an extra moment row. This is only a finite dictionary diagnostic.
 raw=[[p['distances'][j] for p in points] for j in range(n)]+[[F(1)]*m]
 raw_target=[v-old_y for v in obs]+[F(1)];fixed_bases=bases(raw,raw_target);assert fixed_bases
 fixed_min=min(sum(v>0 for v in x) for x in fixed_bases)
 assert fixed_min<=n+1
 return dict(name=name,edges=E,reference_nodes=R,leaf_exposing_lambda=lam,exposing_values=g,
  harmonic_coordinates=H,leaf_weights=w,kappa=kappa,loss=loss,inward_fraction=theta,fixed_y=fixed_y,
  observations=d,strict_n_source_witness=[dict(edge=[r,parent[i]],position=theta,weight=w[i]) for i,r in enumerate(R)],
  uniform_leaf_mass_lower_bound=quantitative,fixed_y_endpoint_vertices=vertices,
  actual_minimum_leaf_mass_over_vertices=min(x[V.index(r)] for x in vertices for r in R),
  strict_fixed_y_minimum_source_edges=n,
  finite_point_compression=dict(points=points,input_weights=weights,input_y=old_y,input_observations=obs,
   free_y_optimal_weights=best,free_y_optimal_pendant=new_y,free_y_source_count=sum(v>0 for v in best),
   free_y_basic_vertex_count=len(bv),fixed_y_finite_dictionary_minimum_sources=fixed_min,
   fixed_y_basic_vertices=fixed_bases,
   scope='The fixed-y finite-dictionary minimum is NOT a lower bound for continuous source positions.'),passed=True)

def main():
 start=time.perf_counter();spec=[('three_leaf_star',[(0,3,1),(1,3,2),(2,3,4)]),
 ('four_leaf_binary',[(0,4,1),(1,4,2),(2,5,3),(3,5,5),(4,5,2)]),
 ('five_leaf_caterpillar',[(0,5,1),(1,5,2),(2,6,3),(3,7,4),(4,7,5),(5,6,2),(6,7,3)]),
 ('six_leaf_three_cherries',[(0,6,1),(1,6,2),(2,7,3),(3,7,4),(4,8,5),(5,8,6),(6,9,2),(7,9,3),(8,9,5)]),
 ('six_leaf_caterpillar',[(0,6,1),(1,6,2),(2,7,3),(3,8,4),(4,9,5),(5,9,6),(6,7,2),(7,8,3),(8,9,4)])]
 cases=[analyze(name,E) for name,E in spec]
 # A path is the necessary exception to a claimed n-source sharp bound.
 L=F(7);t=[F(1,5),F(2,3),F(4,5)];p=[F(1,4),F(1,3),F(5,12)];mean=sum(a*b for a,b in zip(t,p))
 assert 0<mean<1
 path_control=dict(length=L,source_positions=t,weights=p,fixed_y=F(2),single_source_position=mean,
  observations=[L*mean+2,L*(1-mean)+2],minimum_strict_sources=1)
 result=dict(created_utc=datetime.now(timezone.utc).isoformat(),passed=True,cases=cases,path_exception=path_control,
  script_sha256=sha(Path(__file__)),case_count=len(cases),fixed_y_vertex_count=sum(len(c['fixed_y_endpoint_vertices']) for c in cases),
  elapsed_seconds=time.perf_counter()-start,scope=['Sharp fixed-y n lower bounds are for n>=3; a path needs only one source.',
  'Free positive pendant strict compression uses a classical finite moment LP and retains original strict points.',
  'The general continuous fixed-y strict upper bound established here is n+1, not n.',
  'Finite point dictionary obstructions are not claimed to establish continuous strict n+1 lower bounds.'])
 (ROOT/'results/stage17_strict_bounds.json').write_text(json.dumps(result,indent=2,default=str)+'\n')
 print(json.dumps({k:v for k,v in result.items() if k not in ['cases','path_exception']},indent=2))
 print([(c['name'],len(c['fixed_y_endpoint_vertices']),c['finite_point_compression']['free_y_source_count'],c['finite_point_compression']['fixed_y_finite_dictionary_minimum_sources']) for c in cases])
if __name__=='__main__':main()
