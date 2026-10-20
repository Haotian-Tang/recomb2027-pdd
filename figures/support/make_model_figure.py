"""Draw a symbolic model schematic; no simulation or scientific fitting."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10.5,
                     'pdf.fonttype':42,'ps.fonttype':42})
fig, axes = plt.subplots(1, 2, figsize=(6.5, 2.15), gridspec_kw={'width_ratios':[1,1.3]})
navy, teal, orange, gray = '#17324d', '#007d80', '#c56616', '#a2a8ae'
for ax in axes:
    ax.set_axis_off()
def edge(ax, a, b, color=navy, ls='-', lw=1.6):
    ax.plot([a[0],b[0]],[a[1],b[1]],color=color,lw=lw,ls=ls,zorder=1)
def leaf(ax, p, label, offset=(0,0), muted=False):
    ax.scatter(*p,s=24,facecolor='white',edgecolor=gray if muted else navy,zorder=3)
    ax.text(p[0]+offset[0],p[1]+offset[1],label,ha='center',va='center',color=gray if muted else navy)
ax=axes[0]
u,v=(.9,.8),(1.9,.8)
nodes={'a':(.15,1.45),'b':(.15,.15),'c':(2.65,1.45),'d':(2.65,.15)}
edge(ax,u,v)
for k in ['a','b']:edge(ax,u,nodes[k])
for k in ['c','d']:edge(ax,v,nodes[k])
ax.scatter([u[0],v[0]],[u[1],v[1]],s=14,color=navy,zorder=3)
for k,p in nodes.items():leaf(ax,p,'$'+k+'$',(-.13 if k in 'ab' else .13,0))
q1,q2=(.50,1.147), (2.29,.462)
ax.scatter(*q1,s=55,color=teal,zorder=4)
ax.scatter(*q2,s=55,color=orange,zorder=4)
ax.text(.7,1.53,r'$q_1,\ p_1$',color=teal,ha='left')
ax.text(1.4,.10,r'$q_2,\ p_2$',color=orange,ha='left')
ax.text(1.4,-.42,r'$d_r=y+\sum_i p_i\,d_T(r,q_i)$',ha='center')
ax.set_xlim(-.15,3.0);ax.set_ylim(-.56,1.93)
ax.text(-.1,1.78,'A   Unknown source edges',fontweight='bold',ha='left')

ax=axes[1]
g=(.8,.8)
a,b,c,z=(.1,1.4),(.1,.2),(1.45,1.4),(1.45,.2)
for p in [a,b,c]:edge(ax,g,p)
edge(ax,g,z,gray,'--')
for p,label,ofs in [(a,'$a$',(-.12,.02)),(b,'$b$',(-.12,-.02)),(c,'$c$',(.1,.02))]:leaf(ax,p,label,ofs)
leaf(ax,z,'$z$',(.1,-.03),True)
q=(1.15,.477)
ax.scatter(*q,s=55,color=orange,zorder=4)
ax.text(1.45,.58,'$q$',color=orange)
ax.scatter(*g,s=16,color=navy,zorder=3)
ax.text(.74,1.04,'$g$',color=navy)
ax.annotate('',xy=(2.13,.80),xytext=(1.69,.80),arrowprops={'arrowstyle':'->','color':navy,'lw':1.5})
g2=(2.82,.8)
for p,label,ofs in [((2.18,1.4),'$a$',(-.08,.06)),((2.18,.2),'$b$',(-.08,-.04)),((3.48,1.4),'$c$',(.08,.06))]:
    edge(ax,g2,p);leaf(ax,p,label,ofs)
ax.scatter(*g2,s=60,facecolor=orange,edgecolor=navy,zorder=4)
ax.text(2.95,.55,'marked $g$',ha='left',color=orange)
ax.text(1.85,-.42,r"$y'=y+p\,d_T(q,g)$",ha='center')
ax.set_xlim(-.13,3.85);ax.set_ylim(-.56,1.93)
ax.text(-.08,1.78,'B   Projection through a missing branch',fontweight='bold',ha='left',fontsize=10.2)
fig.subplots_adjust(left=.015,right=.99,bottom=.035,top=.99,wspace=.10)
for ext in ['pdf','png']:
    fig.savefig(ROOT/('model_projection.'+ext),dpi=220)
plt.close(fig)
