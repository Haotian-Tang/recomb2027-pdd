#!/usr/bin/env python3
"""Render archived measurements only; never import or rerun an inference method.

Run from any directory: python3 -B <this file> (requires matplotlib).
Writes only this manuscript's performance.{pdf,png} and provenance JSON.
Timeouts have separate categorical axes, never an imputed CPU runtime.
"""
from __future__ import annotations

from collections import Counter
import hashlib
import json
import os
from pathlib import Path
from statistics import median

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
from matplotlib.text import Text

MANUSCRIPT = Path(__file__).resolve().parents[1]
PACKAGED_ARCHIVE = MANUSCRIPT / 'support' / 'figure_archive'
ROOT = PACKAGED_ARCHIVE if PACKAGED_ARCHIVE.is_dir() else MANUSCRIPT.parents[1]

FILES = {
    "development_protocol": "results/stage42_benchmark_protocol.json",
    "development_results": "results/stage42_benchmark.json",
    "revised_development_protocol": "results/stage43_benchmark_protocol.json",
    "revised_development_results": "results/stage43_benchmark.json",
    "local_lp_validation": "results/stage43_solver_validation.json",
    "fixed_one_validation": "results/stage43_fixed_one_audit.json",
    "fixed_one_extra_validation": "results/stage43_fixed_one_extra_audit.json",
    "confirmation_protocol": "results/stage44_confirmation_protocol.json",
    "confirmation_results": "results/stage44_confirmation.json",
    "confirmation_integrity": "results/stage44_confirmation_integrity.json",
    "confirmation_replay": "results/stage44_confirmation_replay.json",
    "posthoc_both_protocol": "results/stage44_confirmation_posthoc_protocol.json",
    "posthoc_both_results": "results/stage44_confirmation_posthoc.json",
    "posthoc_remaining_protocol": "results/stage45_remaining_protocol.json",
    "posthoc_remaining_results": "results/stage45_remaining_checks.json",
    "posthoc_remaining_integrity": "results/stage45_remaining_integrity.json",
    "binary_confirmation_protocol": "results/stage50_confirmation_protocol.json",
    "binary_confirmation_results": "results/stage50_confirmation.json",
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_archives():
    # Parse every full input/result, not only the stored headline summaries.
    documents = {key: json.loads((ROOT / path).read_text()) for key, path in FILES.items()}
    protocol = documents["confirmation_protocol"]
    result = documents["confirmation_results"]
    assert result["protocol_sha256"] == sha(ROOT / FILES["confirmation_protocol"])
    assert protocol["repetitions"] == 2 and len(protocol["cases"]) == 48
    assert protocol["per_job_hard_wall_seconds"] == 2.0
    assert protocol["per_job_address_space_bytes"] == 2**30
    stream_path = ROOT / "results/stage44_confirmation_runs.jsonl"
    stream = [json.loads(line) for line in stream_path.read_text().splitlines()]
    lookup = {(r["name"], r["repetition"], r["method"]): r["result"] for r in stream}
    assert len(lookup) == len(stream) == 256
    for case, record in zip(protocol["cases"], result["records"], strict=True):
        assert case["name"] == record["name"]
        for repetition in record["repetitions"]:
            for method, run in repetition["runs"].items():
                assert lookup[case["name"], repetition["repetition"], method] == run
    return documents, stream_path


def paired_points(documents, task, enumeration, structural):
    points = []
    for record in documents["confirmation_results"]["records"]:
        if record["task"] != task:
            continue
        row = {"case_id": record["name"], "task": task,
               "original_leaves": record["original_leaves"], "cap": record["cap"]}
        for label, method in (("enumeration", enumeration), ("structural", structural)):
            runs = [r["runs"][method] for r in record["repetitions"]]
            assert len(runs) == 2
            assert len({(r["status"], r.get("feasible"), r.get("minimum")) for r in runs}) == 1
            complete = all(r["status"] == "complete" for r in runs)
            row[label] = {
                "archive_method": method,
                "statuses": [r["status"] for r in runs],
                "child_cpu_seconds": [r["child_cpu_seconds"] for r in runs],
                "process_wall_seconds": [r["process_wall_seconds"] for r in runs],
                "complete_both": complete,
                "median_complete_cpu_seconds": median(r["child_cpu_seconds"] for r in runs) if complete else None,
            }
        points.append(row)
    paired = [p for p in points if p["enumeration"]["complete_both"] and p["structural"]["complete_both"]]
    ratios = [p["enumeration"]["median_complete_cpu_seconds"] / p["structural"]["median_complete_cpu_seconds"] for p in paired]
    wins = Counter()
    for p in paired:
        x, y = (p[k]["median_complete_cpu_seconds"] for k in ("enumeration", "structural"))
        wins["structural" if y < x else "enumeration" if x < y else "ties"] += 1
    statuses = Counter((p["enumeration"]["complete_both"], p["structural"]["complete_both"]) for p in points)
    stats = {
        "inputs": len(points), "joint_complete_inputs": len(paired),
        "enumeration_complete": sum(p["enumeration"]["complete_both"] for p in points),
        "structural_complete": sum(p["structural"]["complete_both"] for p in points),
        "enumeration_faster": wins["enumeration"], "structural_faster": wins["structural"], "ties": wins["ties"],
        "median_casewise_cpu_ratio_enumeration_over_structural": median(ratios),
        "only_enumeration_completes": statuses[True, False],
        "only_structural_completes": statuses[False, True], "neither_completes": statuses[False, False],
    }
    saved = documents["confirmation_results"]["summary"]["tasks"][task]["pairwise"][f"{enumeration} vs {structural}"]
    assert stats["joint_complete_inputs"] == saved["joint_complete_both_repetitions"]
    assert stats["median_casewise_cpu_ratio_enumeration_over_structural"] == saved["median_case_ratio_first_over_second"]
    assert wins["enumeration"] == saved["first_lower_cpu_cases"]
    assert wins["structural"] == saved["second_lower_cpu_cases"]
    return points, stats


def panel(fig, spec, title, points, stats, letter):
    gs = spec.subgridspec(2, 2, height_ratios=[0.15, 1], width_ratios=[1, 0.15], hspace=0.06, wspace=0.06)
    ax = fig.add_subplot(gs[1, 0])
    top = fig.add_subplot(gs[0, 0], sharex=ax)
    right = fig.add_subplot(gs[1, 1], sharey=ax)
    corner = fig.add_subplot(gs[0, 1])
    low, high = 0.0025, 2.0
    teal, red = "#147D92", "#B54937"
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(low, high); ax.set_ylim(low, high)
    ax.plot([low, high], [low, high], color="#858B91", lw=0.9, ls=(0, (3, 3)), zorder=1)
    ax.grid(which="major", color="#E6E8EB", linewidth=0.6)
    ticks = [0.01, 0.1, 1]
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(FixedLocator(ticks))
        axis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
        axis.set_minor_locator(NullLocator())
    top.set_ylim(0, 1); right.set_xlim(0, 1)
    for margin in (top, right, corner):
        margin.set_facecolor("#F0F2F4")
        margin.tick_params(which="both", left=False, bottom=False, labelleft=False, labelbottom=False)
        for spine in margin.spines.values():
            spine.set_visible(False)
    top.text(0.02, 0.5, "S timeout", transform=top.transAxes, fontsize=10, va="center", color="#555B62")
    both = 0
    for p in points:
        x, y = (p[k]["median_complete_cpu_seconds"] for k in ("enumeration", "structural"))
        if x is not None and y is not None:
            ax.scatter(x, y, s=22, c=teal if y < x else red if x < y else "#555555", edgecolors="white", linewidths=0.35, zorder=3)
        elif x is not None:
            top.scatter(x, 0.5, marker="^", s=31, color=red, edgecolors="white", linewidths=0.4, zorder=4)
        elif y is not None:
            right.scatter(0.5, y, marker=">", s=31, color=teal, edgecolors="white", linewidths=0.4, zorder=4)
        else:
            both += 1
    corner.text(0.5, 0.5, f"{both}\nboth", ha="center", va="center", fontsize=10, color="#555B62")
    ax.set_xlabel("Enumeration CPU time (s)", labelpad=3)
    if letter == "a":
        ax.set_ylabel("Structural CPU time (s)", labelpad=3)
    top.text(0, 1.32, f"{letter}  {title}", transform=top.transAxes, fontsize=10, weight="bold")
    ax.text(0.025, 0.965,
            f"E {stats['enumeration_complete']}/{stats['inputs']}; S {stats['structural_complete']}/{stats['inputs']} complete",
            transform=ax.transAxes, fontsize=10, va="top")
    for spine in ax.spines.values():
        spine.set_color("#9CA3AA"); spine.set_linewidth(0.65)


def main():
    documents, stream_path = load_archives()
    specifications = [
        ("free_pendant", "original_support", "stage43_equality_reduced", "Free mean pendant", "a"),
        ("fixed_one_missing", "original_fixed_oracle", "stage43_fixed_one", "Fixed mean; one missing leaf", "b"),
    ]
    extracted = {task: paired_points(documents, task, enum, struct) for task, enum, struct, _, _ in specifications}
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.titlesize": 11, "axes.labelsize": 10,
                         "pdf.fonttype": 42, "ps.fonttype": 42,
                         "savefig.facecolor": "white"})
    fig = plt.figure(figsize=(6.5, 3.0))
    outer = fig.add_gridspec(1, 2, left=0.14, right=0.98, bottom=0.19, top=0.83, wspace=0.23)
    for index, (task, _, _, title, letter) in enumerate(specifications):
        panel(fig, outer[index], title, *extracted[task], letter)
    # Check all visible labels against the fixed PDF media box, not a tight crop.
    # Padding is kept inside the 6.5-inch figure so TeX cannot place ink outside
    # the manuscript text block when including this PDF at full width.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    text_bounds = []
    for text in fig.findobj(match=Text):
        if text.get_visible() and text.get_text().strip():
            box = text.get_window_extent(renderer=renderer)
            points = [value * 72 / fig.dpi for value in (box.x0, box.y0, box.x1, box.y1)]
            text_bounds.append({"text": text.get_text(), "bbox_points": points, "font_points": text.get_fontsize()})
    width_pt, height_pt = (v * 72 for v in fig.get_size_inches())
    for record in text_bounds:
        x0, y0, x1, y1 = record["bbox_points"]
        assert 2 <= x0 <= x1 <= width_pt - 2 and 2 <= y0 <= y1 <= height_pt - 2, record
        assert record["font_points"] >= 10, record
    layout_check = {"media_box_points": [width_pt, height_pt], "minimum_internal_margin_points": 2,
                    "all_visible_text_inside_media_box": True, "minimum_font_points": 10,
                    "text_bounds": text_bounds}
    figure_dir = MANUSCRIPT
    figure_dir.mkdir(exist_ok=True)
    fig.savefig(figure_dir / "performance.pdf", metadata={"Title": "Frozen exact-source confirmation", "Creator": "Archived-data matplotlib rendering", "CreationDate": None, "ModDate": None})
    fig.savefig(figure_dir / "performance.png", dpi=300)
    plt.close(fig)

    binary = documents["binary_confirmation_results"]
    binary_protocol = documents["binary_confirmation_protocol"]
    assert binary["protocol_sha256"] == sha(ROOT / FILES["binary_confirmation_protocol"])
    assert len(binary["rows"]) == len(binary_protocol["cases"]) == 14
    assert binary_protocol["base_indices"] == list(range(12)) + [0, 1]
    assert all(r["agreement"] is True for r in binary["rows"])
    provenance = {
        "purpose": "Manuscript extraction and vector rendering of existing same-task measurements; no experiment rerun.",
        "archive_root_relative_to_manuscript": os.path.relpath(ROOT, MANUSCRIPT),
        "input_files": {key: {"path": path, "sha256": sha(ROOT / path)} for key, path in FILES.items()},
        "raw_primary_stream": {"path": str(stream_path.relative_to(ROOT)), "sha256": sha(stream_path), "jobs": 256},
        "historical_source_hashes": documents["confirmation_protocol"]["script_hashes"],
        "machine_as_registered": documents["confirmation_protocol"]["machine"],
        "hardware_limitations": "The frozen protocol records x86_64 Linux and 62 visible CPUs, not processor model, clock, RAM capacity, affinity, or competing load. No current-machine specification is substituted.",
        "primary_protocol": {
            "input_count": 48, "repetitions": 2, "jobs": 256, "seed": 4410142026,
            "wall_cap_seconds": 2, "address_space_cap_bytes": 2**30,
            "cpu_unit": "child user plus system seconds from wait4",
            "isolation": documents["confirmation_protocol"]["process_isolation"],
            "timing_scope": documents["confirmation_protocol"]["timing"],
            "method_order": documents["confirmation_protocol"]["execution_order"],
            "thread_environment": documents["confirmation_protocol"]["thread_environment"],
        },
        "plot_methods": "Current structural prototype versus original-edge support enumeration. The older structural prototype remains in the appendix completion table and all raw measurements.",
        "plot_timeout_policy": "A categorical top/right margin preserves every timeout; no CPU-time value is assigned to a timeout and no follow-up time enters the figure.",
        "plot_statistics": {task: stats for task, (_, stats) in extracted.items()},
        "plot_points": {task: points for task, (points, _) in extracted.items()},
        "all_methods_primary_summary": documents["confirmation_results"]["summary"],
        "development_not_confirmation": {
            "initial": documents["development_results"]["summary"],
            "same_30_reused_after_preprocessing": documents["revised_development_results"]["summary"],
            "scope": "These 30 inputs informed equality-reduction engineering; they are not additional out-of-sample confirmation inputs.",
        },
        "posthoc_correctness_only": {
            "wall_cap_per_job_seconds": 30, "first_jobs": 2, "remaining_jobs": 5,
            "first_new_witnesses": 2, "remaining_new_witnesses": 3,
            "completed_agreement_inputs": 48, "primary_timeouts_retained": True,
            "summary": documents["posthoc_remaining_results"]["summary"],
            "independent_replay_scope": documents["posthoc_remaining_integrity"]["scope"],
        },
        "binary_reduction_controls": {
            "summary": binary["summary"], "base_indices": binary_protocol["base_indices"],
            "wall_cap_seconds": binary_protocol["hard_seconds_per_case"],
            "logic": binary_protocol["inference"],
            "cases": [{"id": r["name"], "base_index": r["base_index"],
                       "oracle_status": r["oracle"]["status"], "agreement": r["agreement"],
                       "strict_minimum": r["certified_strict_minimum"]} for r in binary["rows"]],
            "limits": "12 base inputs plus two deterministic rescalings, not 14 independent inputs. Closed lower bounds are combined with separate strict witnesses; no timing comparison or finite proof of hardness.",
        },
        "scientific_scope": "Exact minimum strict source count within a budget on small synthetic trees; no biological recovery, empirical coverage, published-software comparison, or universal speed claim.",
        "new_work": {"solver_calls": 0, "new_inputs": 0, "new_measurements": 0,
                     "work_performed": "Read saved JSON/JSONL, reconcile records, aggregate recorded CPU values, and render graphics."},
        "outputs": {name: sha(figure_dir / name) for name in ["performance.pdf", "performance.png"]},
        "layout_check": layout_check,
        "renderer_sha256": sha(Path(__file__)),
    }
    (MANUSCRIPT / "support" / "experiment_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps({"figure": "figures/performance.pdf", "statistics": provenance["plot_statistics"], "new_solver_calls": 0}, indent=2))


if __name__ == "__main__":
    main()
