#!/usr/bin/env python3
"""Redraw main-text Figure 2, "What incomplete references cost", from its payload.

The three panels are:
  (a) median fall in Theta_1 against the share of references removed, exact arm, by removal regime;
  (b) the same draws on the real krepp ANI estimates of the same trees and the identical leaf sets;
  (c) references erasable before the minimum-size explanation stops being unique, against panel
      size, pooled over K = 2-10, adversarial (nearest to sources) against uniformly at random.

The payload `figures/data/fig_incompleteness_data.json` was derived from two registered records:
  results/stage299/stage299_regime_ratios.json  -> the `fall` block (median fall by arm, regime and
                                                    level, each with its draw count)
  results/stage297/stage297_analysis.json       -> the `threshold_vs_panel` block (one
                                                    [panel_size, threshold, K] triple per case)

Run from the repository root:  python3 figures/make_fig_incompleteness.py
Writes figures/incompleteness.pdf and figures/incompleteness.png. Needs numpy and matplotlib.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FixedFormatter, FixedLocator, NullLocator

ROOT = Path(__file__).resolve().parents[1]
PAYLOAD = ROOT / "figures" / "data" / "fig_incompleteness_data.json"

LEVELS = ["one", "0.10", "0.25", "0.50", "0.75"]
LEVEL_LABELS = ["1", "10%", "25%", "50%", "75%"]
COLOUR = {"nearest": "#b2182b", "clade": "#ef8a62", "taxon": "#2166ac", "random": "#92c5de"}
LABEL = {"nearest": "nearest to sources", "clade": "whole clade",
         "taxon": "whole taxon", "random": "random"}


def style():
    plt.rcParams.update({
        "font.size": 10, "axes.labelsize": 10, "axes.titlesize": 10,
        "xtick.labelsize": 10, "ytick.labelsize": 10, "legend.fontsize": 10,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.dpi": 150, "savefig.bbox": None, "pdf.fonttype": 42,
    })


def series(fall, arm, regime):
    block = fall.get(arm, {}).get(regime, {})
    return np.array([block[l][0] if l in block else np.nan for l in LEVELS], dtype=float)


def spearman(x, y):
    # Rank correlation with midranks, so the annotation tracks the data.
    x, y = np.asarray(x, float), np.asarray(y, float)

    def ranks(v):
        order = np.argsort(v, kind="mergesort")
        r = np.empty(len(v), float)
        i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
                j += 1
            r[order[i:j + 1]] = (i + j) / 2.0 + 1.0
            i = j + 1
        return r

    a, b = ranks(x), ranks(y)
    return float(np.corrcoef(a, b)[0, 1])


def main():
    style()
    payload = json.loads(PAYLOAD.read_text())
    fall, points = payload["fall"], payload["threshold_vs_panel"]
    x = np.arange(len(LEVELS))

    fig, axes = plt.subplots(1, 3, figsize=(6.5, 2.14))

    for ax, (arm, title) in zip(axes[:2], [("exact_partner", "Exact distances"),
                                           ("krepp_ani", "Real ANI distances")]):
        for regime in ("nearest", "clade", "taxon", "random"):
            ax.plot(x, series(fall, arm, regime), marker="o", ms=3.0, lw=1.4,
                    color=COLOUR[regime], zorder=3 if regime != "random" else 2)
        ax.set_xticks(x)
        ax.set_xticklabels([t if i % 2 == 0 else "" for i, t in enumerate(LEVEL_LABELS)])
        ax.set_ylim(-0.04, 1.0)
        ax.set_xlim(-0.3, 4.5)
        ax.set_xlabel("references removed")
        ax.set_title(title)
    axes[0].set_ylabel(r"median fall in $\Theta_1$")
    axes[1].set_yticklabels([])
    for regime, xy, align in [("clade", (0.12, 0.98), "left"),
                              ("nearest", (0.12, 0.83), "left"),
                              ("taxon", (0.12, 0.68), "left"),
                              ("random", (0.12, 0.53), "left")]:
        axes[1].annotate(LABEL[regime], xy, ha=align, va="center",
                         color=COLOUR[regime])

    adversarial = np.asarray(points["adversarial_nearest"], dtype=float)
    random_arm = np.asarray(points["random"], dtype=float)
    ax = axes[2]
    ax.scatter(random_arm[:, 0], random_arm[:, 1], s=2.6, color=COLOUR["random"],
               alpha=0.30, lw=0, zorder=1)
    ax.scatter(adversarial[:, 0], adversarial[:, 1], s=2.6, color=COLOUR["nearest"],
               alpha=0.40, lw=0, zorder=2)
    bins = np.logspace(np.log10(19), np.log10(365), 7)
    for arr, colour in ((random_arm, COLOUR["random"]), (adversarial, COLOUR["nearest"])):
        xs, ys = [], []
        for lo, hi in zip(bins[:-1], bins[1:]):
            mask = (arr[:, 0] >= lo) & (arr[:, 0] < hi)
            if mask.sum() >= 8:
                xs.append(np.sqrt(lo * hi))
                ys.append(np.median(arr[mask, 1]))
        ax.plot(xs, ys, color=colour, lw=1.7, marker="s", ms=3.2, zorder=4,
                mec="white", mew=0.6)
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=1, linscale=0.45)
    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_locator(FixedLocator([20, 50, 100, 361]))
    ax.xaxis.set_major_formatter(FixedFormatter(["20", "50", "100", "361"]))
    ax.yaxis.set_major_locator(FixedLocator([0, 1, 10, 100]))
    ax.yaxis.set_major_formatter(FixedFormatter(["0", "1", "10", "100"]))
    ax.set_ylim(-0.3, 1400)
    ax.set_xlim(17, 430)
    ax.set_xlabel("panel size")
    ax.set_ylabel("references erasable")
    ax.set_title("Tolerance by\npanel size")
    rho_random = spearman(random_arm[:, 0], random_arm[:, 1])
    rho_adv = spearman(adversarial[:, 0], adversarial[:, 1])
    ax.annotate(r"random $\rho=%.2f$" % rho_random, (21, 640), ha="left", va="center",
                color="#4b93b8")
    ax.annotate("nearest\n" + r"$\rho=%.2f$" % rho_adv, (21, 170), ha="left",
                va="center", color=COLOUR["nearest"], linespacing=0.9)

    for letter, ax in zip("abc", axes):
        ax.margins(0.04)
        ax.text(-0.20, 1.20, letter, transform=ax.transAxes, fontsize=10,
                fontweight="bold", va="bottom", ha="left")
    # Fixed canvas, every ordinary label at 10 pt: reserve explicit space rather
    # than let tight layout expand the PDF and shrink the text on inclusion.
    fig.subplots_adjust(left=0.115, right=0.985, bottom=0.220, top=0.790, wspace=0.52)

    out = ROOT / "figures"
    fig.savefig(out / "incompleteness.pdf", bbox_inches=None)
    fig.savefig(out / "incompleteness.png", dpi=300, bbox_inches=None)
    print("wrote", out / "incompleteness.pdf", "and", out / "incompleteness.png")
    print("  panel (c): %d adversarial (rho=%.4f) and %d random (rho=%.4f) measurements"
          % (len(adversarial), rho_adv, len(random_arm), rho_random))


if __name__ == "__main__":
    main()
