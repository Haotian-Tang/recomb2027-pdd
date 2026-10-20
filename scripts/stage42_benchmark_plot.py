#!/usr/bin/env python3
"""Plot immutable Stage42 benchmark results without rerunning any solver."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
records = json.loads((ROOT / 'results/stage42_benchmark.json').read_text())['records']
style = {
    'binary_generated': ('#2378ad', 'Binary generated'),
    'essential_two_deficits': ('#9854a4', 'Essential two-deficit'),
    'known_strict_budget_boundary': ('#ba7817', 'Strict budget controls'),
    'known_closure_boundary': ('#737373', 'Closure-only controls'),
    'multifurcating_marked_gates': ('#25875a', 'Marked gates'),
}
plt.rcParams.update({'font.size': 10, 'axes.labelsize': 10, 'xtick.labelsize': 10,
                    'ytick.labelsize': 10, 'legend.fontsize': 10,
                    'pdf.fonttype': 42, 'ps.fonttype': 42})
fig, ax = plt.subplots(figsize=(6.5, 5.2))
for group, (color, label) in style.items():
    selected = [r for r in records if r['group'] == group and r['both_complete']]
    ax.scatter([r['runs']['original_support']['child_cpu_seconds'] for r in selected],
               [r['runs']['marked_deficit']['child_cpu_seconds'] for r in selected],
               color=color, s=30, edgecolor='white', linewidth=0.5, label=label, zorder=3)
for r in records:
    if r['both_complete']:
        continue
    baseline, proposed = [r['runs'][m] for m in ('original_support', 'marked_deficit')]
    assert baseline['status'] == 'unknown' and proposed['status'] == 'complete'
    fig.text(0.5, 0.87,
             'One additional case: baseline UNKNOWN at the 5 s wall cap;\n'
             f"structural solver completed in {proposed['child_cpu_seconds']:.4f} CPU s.",
             ha='center', va='center', fontsize=10,
             bbox={'boxstyle': 'round,pad=0.4', 'facecolor': '#f6f6f6', 'edgecolor': '#aaaaaa', 'lw': 0.6})
ax.plot([0.002, 5], [0.002, 5], '--', color='#999999', lw=0.8, zorder=1)
ax.set(xscale='log', yscale='log', xlim=(0.002, 5), ylim=(0.002, 5),
       xlabel='Original-support baseline: CPU seconds',
       ylabel='Marked-deficit prototype: CPU seconds')
ax.text(0.04, 0.94, 'Baseline faster', transform=ax.transAxes, va='top', fontsize=10)
ax.text(0.96, 0.04, 'Structural solver faster', transform=ax.transAxes,
        ha='right', fontsize=10)
ax.grid(True, which='major', color='#e9e9e9', linewidth=0.6)
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, 0.07),
           frameon=False, fontsize=10, ncol=2, columnspacing=1.0)
fig.suptitle('Exact minimum-source tasks: 29 jointly completed cases', fontsize=11, y=0.98)
fig.text(0.5, 0.015, 'One timed trial per method/case.\nA wall-time limit is not a CPU completion time.',
         ha='center', fontsize=10)
fig.tight_layout(rect=(0, 0.24, 1, 0.82))
(ROOT / 'exports').mkdir(exist_ok=True)
for suffix in ('pdf', 'png'):
    fig.savefig(ROOT / f'exports/stage42_benchmark.{suffix}', dpi=200)
from matplotlib.text import Text
text_sizes = [item.get_fontsize() for item in fig.findobj(Text) if item.get_text().strip()]
assert text_sizes and min(text_sizes) >= 10, text_sizes
print(f'Font audit: minimum nominal text size {min(text_sizes):g} pt at 6.5-inch figure width')
print('Saved exports/stage42_benchmark.pdf and .png')
