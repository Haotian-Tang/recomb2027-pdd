#!/usr/bin/env python3
"""Rebuild Figure 3, "The tolerance axis".

The original ad hoc plotting code for this figure was not retained.  This script
replaces it and reads every plotted value directly from the registered record
`results/stage296/stage296_payload.json`, so there is no intermediate payload
that can drift away from the measurements.

Two deliberate differences from the figure it replaces:

  * every error row now shows the WORST-LEAF residual, the same statistic the
    first-drop tolerance and Theta_k are defined with.  The old figure mixed
    median-leaf residuals for the simulated-noise rows with worst-leaf for the
    real ANI row, which made the row-to-row comparison ill-defined.
  * the panel (a) title no longer carries a numeric multiplier.  The old
    "10-40x" could not be reproduced from any statistic in the record.  The
    supported statement, the paired per-case median ratio of 14.5, is in the
    caption where it can be checked.

Run from the repository root:  python3 figures/make_fig_tolerance_axis.py
Writes figures/tolerance_axis.pdf and figures/tolerance_axis.png.
Needs numpy and matplotlib.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullLocator

ROOT = Path(__file__).resolve().parents[1]
RECORD = ROOT / "results" / "stage296" / "stage296_payload.json"

RED, ORANGE, BLUE, PALE = "#b2182b", "#ef8a62", "#4b93b8", "#92c5de"

# label, record accessor, colour.  Bottom row first.
# The median-leaf values are retained only for the shaded reference band.
# Row labels are deliberately terse: RECOMB requires a 10 pt minimum for
# ordinary figure text, and at 10 pt the descriptive forms do not fit beside a
# panel this wide.  The caption carries the full definitions.
ROWS = [
    ("first drop", "first_drop", "0.25"),
    ("noise 2, 10 kb", "published_noise2_l10000", PALE),
    ("noise 2, 100 bp", "published_noise2_l100", PALE),
    ("noise 1", "noise1_pooled", BLUE),
    ("real ANI", "krepp_max", RED),
]


def style():
    plt.rcParams.update({
        "font.size": 10, "axes.labelsize": 10, "axes.titlesize": 10,
        "xtick.labelsize": 10, "ytick.labelsize": 10, "legend.fontsize": 10,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6,
        "ytick.major.width": 0.6, "figure.dpi": 150,
        "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.bbox": None,
    })


def collect(payload):
    cal, full = payload["calibration_placed"], payload["full_panel"]
    out = {}
    out["first_drop"] = [c["first_drop_eps_rel"] for c in full
                         if c.get("first_drop_eps_rel")]
    out["krepp_max"] = [c["krepp_ani"]["residual_max_rel"] for c in cal
                        if "krepp_ani" in c]
    out["krepp_median"] = [c["krepp_ani"]["residual_median_rel"] for c in cal
                           if "krepp_ani" in c]
    out["noise1_pooled"] = [c[k]["residual_max_rel"] for c in cal
                            for k in ("published_noise1_l100",
                                      "published_noise1_l1000",
                                      "published_noise1_l10000") if k in c]
    for k in ("published_noise2_l100", "published_noise2_l10000"):
        out[k] = [c[k]["residual_max_rel"] for c in cal if k in c]
    return {k: sorted(v for v in vals if v and v > 0) for k, vals in out.items()}


def panel_a(ax, data, rng):
    lo = float(np.median(data["krepp_median"]))
    hi = float(np.median(data["krepp_max"]))
    ax.axvspan(lo, hi, color=RED, alpha=0.08, lw=0, zorder=0)
    for y, (label, key, colour) in enumerate(ROWS):
        v = np.asarray(data[key], float)
        jitter = rng.uniform(-0.17, 0.17, len(v))
        ax.scatter(v, y + jitter, s=4.5, color=colour, alpha=0.55, lw=0, zorder=2)
        med = float(np.median(v))
        ax.plot([med, med], [y - 0.3, y + 0.3], color=colour, lw=1.8, zorder=3)
        label_on_right = key == "published_noise2_l10000"
        ax.annotate(("%.2g" % med), (med, y + 0.17),
                    xytext=(3 if label_on_right else -3, 2),
                    textcoords="offset points",
                    ha="left" if label_on_right else "right", va="bottom",
                    color=colour)
    ax.set_yticks(range(len(ROWS)))
    ax.set_yticklabels([r[0] for r in ROWS])
    ax.set_ylim(-0.7, len(ROWS) - 0.25)
    ax.set_xscale("log")
    ax.set_xlim(3e-6, 3e0)
    ax.set_xticks([1e-5, 1e-3, 1e-1])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlabel(r"tolerance $\varepsilon$ (relative)", labelpad=2)
    ax.set_title("Error exceeds the\nfirst-drop tolerance", loc="left")
    ax.tick_params(axis="y", length=0)


def panel_b(ax, payload, data):
    # the table also carries two rows keyed by a measured error level rather
    # than a numeric tolerance ("krepp_median", "krepp_max"); the curve uses the
    # numeric grid only, and those two levels appear as the shaded band.
    def numeric(v):
        try:
            return float(v)
        except ValueError:
            return None

    rows = [r for r in payload["within_eps_table"]
            if numeric(r["eps_rel"]) is not None and numeric(r["eps_rel"]) > 0]
    x = np.array([float(r["eps_rel"]) for r in rows])
    rows = [r for _, r in sorted(zip(x, rows), key=lambda z: z[0])]
    x = np.sort(x)
    ax.fill_between(x, [r["share_m_full_is_1"] for r in rows], color="0.85",
                    lw=0, zorder=0)
    lo = float(np.median(data["krepp_median"]))
    hi = float(np.median(data["krepp_max"]))
    ax.axvspan(lo, hi, color=RED, alpha=0.08, lw=0, zorder=1)
    for key, colour, label in (("fall_nearest", RED, "nearest"),
                               ("fall_clade", ORANGE, "clade"),
                               ("fall_random", BLUE, "random")):
        ax.plot(x, [r[key] for r in rows], color=colour, lw=1.5, marker="o",
                ms=2.8, zorder=3, label=label)
    ax.set_xscale("log")
    ax.set_xlim(3e-6, 3e0)
    ax.set_ylim(0, None)
    ax.set_xticks([1e-5, 1e-3, 1e-1])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlabel(r"tolerance $\varepsilon$ (relative)", labelpad=2)
    ax.set_ylabel("fraction of draws\nbelow full-panel count", labelpad=3)
    ax.set_title("Structured loss hurts\nmore than random", loc="left")
    ax.legend(frameon=False, loc="upper right", handlelength=1.1,
              borderaxespad=0.1, labelspacing=0.2)


def main():
    style()
    payload = json.loads(RECORD.read_text())
    data = collect(payload)
    rng = np.random.default_rng(0)
    fig, axes = plt.subplots(1, 2, figsize=(6.5, 2.50),
                             gridspec_kw=dict(width_ratios=[1.12, 1.0]))
    panel_a(axes[0], data, rng)
    panel_b(axes[1], payload, data)
    for letter, ax in zip("ab", axes):
        ax.text(-0.05, 1.04, letter, transform=ax.transAxes, fontsize=10,
                fontweight="bold", va="bottom", ha="right")
    # Explicit margins: the canvas is 6.5 by 2.50 inches and every ordinary
    # label is 10 pt, so the panels are placed rather than tight-layouted.
    axes[0].set_position([0.195, 0.207, 0.295, 0.613])
    axes[1].set_position([0.685, 0.207, 0.300, 0.613])
    out = ROOT / "figures"
    fig.savefig(out / "tolerance_axis.pdf", bbox_inches=None)
    fig.savefig(out / "tolerance_axis.png", dpi=300, bbox_inches=None)
    print("wrote", out / "tolerance_axis.pdf", "and", out / "tolerance_axis.png")
    for label, key, _ in ROWS:
        v = data[key]
        print("  %-34s n=%4d  median=%.5g" % (label, len(v), float(np.median(v))))


if __name__ == "__main__":
    main()
