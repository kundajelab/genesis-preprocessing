"""Vector figures with explicit Arial identity and no silent font substitution."""

from __future__ import annotations

import datetime
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from evidence import Evidence, digest, tsv
from fontTools.ttLib import TTFont
from matplotlib import font_manager
from matplotlib.text import Text

COLORS = ["#2a78d6", "#eb6834", "#1baf7a"]
LABELS = {
    "baseline": "BWA-MEM2",
    "bwa": "BWA",
    "bowtie2": "Bowtie2",
    "cutadapt": "Cutadapt",
    "fastp": "fastp",
    "picard": "Picard",
    "macs2": "MACS2",
    "genrich": "Genrich",
}


def configure(font_dir: Path) -> list[dict[str, Any]]:
    fonts = []
    for path in sorted(font_dir.iterdir()):
        if path.suffix.lower() not in (".ttf", ".otf"):
            continue
        font = TTFont(path)
        family = font["name"].getDebugName(1)
        if family != "Arial":
            continue
        flags = font["OS/2"].fsType
        if flags & (2 | 256 | 512):
            raise ValueError("Arial font does not permit this subset embedding")
        font_manager.fontManager.addfont(str(path))
        fonts.append(
            {
                "filename": path.name,
                "family": family,
                "style": font["name"].getDebugName(2),
                "sha256": digest(path),
                "embedding_flags": flags,
            }
        )
    if not fonts:
        raise ValueError(
            "Actual Arial font files are required; fallback is not allowed"
        )
    selected = Path(font_manager.findfont("Arial", fallback_to_default=False))
    if selected.resolve().parent != font_dir.resolve():
        raise ValueError("Arial resolved outside the declared font directory")
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial"],
            "pdf.fonttype": 42,
            "axes.unicode_minus": False,
            "mathtext.fontset": "dejavusans",
            "font.size": 8.5,
            "axes.linewidth": 0.8,
            "axes.grid": False,
            "legend.frameon": False,
            "axes.spines.top": True,
            "axes.spines.right": True,
        }
    )
    return fonts


def save(fig: Any, out: Path, name: str, rows: list[dict[str, Any]]) -> None:
    fig.canvas.draw()
    canvas = fig.bbox
    renderer = fig.canvas.get_renderer()
    hidden_ticks = set()
    for ax in fig.axes:
        for axis, limits in ((ax.xaxis, ax.get_xlim()), (ax.yaxis, ax.get_ylim())):
            for tick in axis.get_major_ticks() + axis.get_minor_ticks():
                if not min(limits) <= tick.get_loc() <= max(limits):
                    hidden_ticks.update((tick.label1, tick.label2))
    for text in fig.findobj(Text):
        if text in hidden_ticks or not text.get_visible() or not text.get_text():
            continue
        bounds = text.get_window_extent(renderer)
        if (
            bounds.x0 < -1
            or bounds.y0 < -1
            or bounds.x1 > canvas.x1 + 1
            or bounds.y1 > canvas.y1 + 1
        ):
            raise ValueError(f"Clipped figure text in {name}: {text.get_text()}")
    timestamp = datetime.datetime(2026, 10, 16, tzinfo=datetime.timezone.utc)
    fig.savefig(
        out / f"{name}.pdf", metadata={"CreationDate": timestamp, "ModDate": timestamp}
    )
    plt.close(fig)
    tsv(out / f"{name}.tsv", rows)


def quality_figure(
    e: Evidence,
    out: Path,
    name: str,
    cohort: str,
    methods: list[str],
    libraries: list[str],
    metrics: list[tuple[str, str]],
    title: str,
    labels: list[str] | None = None,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(7, 2.6))
    rows = []
    for panel, (metric, ylabel) in enumerate(metrics):
        ax = axes[panel]
        for j, method in enumerate(methods):
            values = [e.q(cohort, method, lib, metric) for lib in libraries]
            if any(v is None for v in values):
                raise ValueError("Figure requires explicit handling of missing values")
            x = [i + (j - (len(methods) - 1) / 2) * 0.2 for i in range(len(libraries))]
            ax.plot(
                x,
                values,
                "o",
                color=COLORS[j],
                label=(labels or [LABELS[m] for m in methods])[j],
                ms=4,
            )
            for lib, value in zip(libraries, values):
                rows.append(
                    {
                        "cohort": cohort,
                        "method": method,
                        "library": lib,
                        "metric": metric,
                        "value": value,
                    }
                )
        ax.set_xticks(
            range(len(libraries)), [lib.replace("SRR", "") for lib in libraries]
        )
        ax.set_xlabel("Library accession (SRR prefix)")
        ax.set_ylabel(ylabel)
        ax.set_title(f"{'ab'[panel]}  {ylabel}")
        ax.margins(x=0.18, y=0.22)
        ax.ticklabel_format(axis="y", style="plain", useOffset=False)
        if metric == "peak_jaccard" and name == "duplicates-quality":
            ax.set_ylim(0.99, 1.001)
        ax.text(
            0.98,
            0.07 if name == "duplicates-quality" else 0.97,
            f"n = {len(libraries)} libraries",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=7,
        )
    axes[0].legend(loc="best", fontsize=7)
    fig.suptitle(title, fontsize=9)
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.21, top=0.77, wspace=0.34)
    save(fig, out, name, rows)


def cost_figure(
    e: Evidence,
    out: Path,
    name: str,
    cohort: str,
    methods: list[str],
    labels: list[str],
) -> None:
    fig, ax = plt.subplots(figsize=(4.6, 2.6))
    rows = [e.c(cohort, m) for m in methods]
    for i, row in enumerate(rows):
        mid, lo, hi = (float(row[k]) for k in ("median", "minimum", "maximum"))
        ax.errorbar(
            mid,
            i,
            xerr=[[mid - lo], [hi - mid]],
            fmt="o",
            color=COLORS[i % 3],
            ecolor="grey",
            capsize=3,
        )
        ax.annotate(
            f"{mid:.1f} [{lo:.1f}, {hi:.1f}] s",
            (mid, i),
            xytext=(0, 9),
            textcoords="offset points",
            ha="center",
            fontsize=7,
        )
    ax.set_yticks(range(len(methods)), labels)
    ax.set_ylim(-0.6, len(methods) - 0.3)
    ax.set_xlim(0, max(float(r["maximum"]) for r in rows) * 1.3)
    ax.set_xlabel("Cohort processing seconds; median and range, n = 3")
    fig.suptitle(
        "Processing cost includes all libraries, excluding indexing and evaluation",
        fontsize=9,
    )
    fig.subplots_adjust(left=0.19, right=0.98, bottom=0.23, top=0.77)
    save(fig, out, name, rows)


def figures(e: Evidence, out: Path) -> None:
    quality_figure(
        e,
        out,
        "alignment-quality",
        "dap-pe-1m",
        ["baseline", "bwa", "bowtie2"],
        ["SRR31238589", "SRR31239194", "SRR31238722", "SRR31238919"],
        [
            ("peak_jaccard", "Peak base-pair Jaccard"),
            ("fixed_count_spearman", "Fixed-region Spearman"),
        ],
        "Bowtie2 differs most in two paired-end DAP libraries",
    )
    for name, methods, labels in [
        (
            "trimming",
            ["baseline", "cutadapt", "fastp"],
            ["Untrimmed", "Cutadapt", "fastp"],
        ),
        ("duplicates", ["baseline", "picard"], ["samtools", "Picard"]),
    ]:
        quality_figure(
            e,
            out,
            name + "-quality",
            "atac-arabidopsis-1m",
            methods,
            ["SRR4000468", "SRR4000469"],
            [
                ("fragment_retention", "Retained fragments / raw pairs"),
                ("peak_jaccard", "Peak base-pair Jaccard"),
            ],
            "Downstream measurements remain close in two short-read ATAC libraries",
            labels,
        )
        cost_figure(e, out, name + "-cost", "atac-arabidopsis-1m", methods, labels)
    cost_figure(
        e,
        out,
        "alignment-cost",
        "dap-pe-1m",
        ["baseline", "bwa", "bowtie2"],
        ["BWA-MEM2", "BWA", "Bowtie2"],
    )
    fig, axes = plt.subplots(1, 2, figsize=(7, 2.6))
    rows = []
    for i, metric in enumerate(["peak_count", "heldout_fragments_FRiP"]):
        for j, lib in enumerate(["SRR4000468", "SRR4000469"]):
            vals = [
                e.q(f"atac-arabidopsis-{d}m", "baseline", lib, metric) for d in (1, 5)
            ]
            axes[i].plot([1, 5], vals, "o-", color=COLORS[j], label=lib)
            for d, v in zip((1, 5), vals):
                rows.append(
                    {
                        "depth_million_pairs": d,
                        "library": lib,
                        "metric": metric,
                        "value": v,
                    }
                )
                axes[i].annotate(
                    f"{v:,.0f}" if i == 0 else f"{v:.3f}",
                    (d, v),
                    xytext=(0, 7 if j == 0 else -13),
                    textcoords="offset points",
                    ha="center",
                    fontsize=7,
                )
        axes[i].set_xticks([1, 5])
        axes[i].set_xlim(0.3, 5.7)
        axes[i].margins(y=0.3)
        axes[i].set_xlabel("Input pairs (millions)")
    axes[0].set_title("a  Peak discovery")
    axes[0].set_ylabel("Peaks")
    axes[0].legend(fontsize=7)
    axes[1].set_title("b  Other-replicate peak overlap")
    axes[1].set_ylabel("Held-out fragment FRiP")
    fig.suptitle(
        "Greater depth increases both peak discovery and held-out overlap", fontsize=9
    )
    fig.subplots_adjust(left=0.09, right=0.98, bottom=0.21, top=0.77, wspace=0.35)
    save(fig, out, "depth", rows)
    fig, axes = plt.subplots(1, 2, figsize=(7, 2.6))
    rows = []
    for j, cpu in enumerate([1, 2, 8]):
        subset = [r for r in e.macs if int(r["cpus"]) == cpu]
        speed = [float(r["speedup"]) for r in subset]
        ratios = [
            float(r["candidate_sampled_group_rss_max_mib"])
            / float(r["baseline_sampled_group_rss_max_mib"])
            for r in subset
        ]
        for ax, vals in zip(axes, [speed, ratios]):
            ax.scatter(
                [cpu + (i - 3.5) * 0.04 for i in range(len(vals))],
                vals,
                color=COLORS[j],
                s=14,
            )
        axes[0].annotate(
            f"{min(speed):.1f} to {max(speed):.1f}",
            (cpu, max(speed)),
            xytext=(0, 6),
            textcoords="offset points",
            ha="center",
            fontsize=7,
        )
        rows.extend(subset)
    for ax in axes:
        ax.axhline(1, color="black", linestyle=":", linewidth=0.8)
        ax.set_xticks([1, 2, 8])
        ax.set_xlim(0.3, 8.7)
        ax.set_xlabel("Allocated CPUs")
        ax.margins(y=0.22)
        ax.set_ylim(bottom=0)
    axes[0].set_title("a  Wall-time speedup")
    axes[0].set_ylabel("Released / optimized wall time")
    axes[1].set_title("b  Sampled process-group memory")
    axes[1].set_ylabel("Optimized / released RSS")
    fig.suptitle(
        "All 24 paired comparisons match outputs; more CPUs increase memory", fontsize=9
    )
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.21, top=0.77, wspace=0.35)
    save(fig, out, "macs-speed", rows)
    fig, axes = plt.subplots(1, 2, figsize=(7, 2.6))
    rows = []
    for j, lib in enumerate(["SRR4000468", "SRR4000469"]):
        for i, key in enumerate(
            ["reciprocal_peak_base_jaccard", "common_region_fragment_spearman"]
        ):
            vals = []
            for depth, source in [(1, e.reference), (5, e.local["reference_5m"])]:
                value = source["libraries"][lib]["metrics"][key]["value"]
                vals.append(value)
                rows.append(
                    {
                        "library": lib,
                        "depth_million_pairs": depth,
                        "metric": key,
                        "value": value,
                    }
                )
            axes[i].plot([1, 5], vals, "o-", color=COLORS[j], label=lib)
            for d, v in zip([1, 5], vals):
                axes[i].annotate(
                    f"{v:.4f}",
                    (d, v),
                    xytext=(
                        0,
                        -12
                        if (i == 0 and d == 5 and j == 0) or (i == 1 and j == 1)
                        else 7,
                    ),
                    textcoords="offset points",
                    ha="center",
                    fontsize=7,
                )
    for ax in axes:
        ax.set_xticks([1, 5])
        ax.set_xlim(0.3, 5.7)
        ax.margins(y=0.45)
        ax.set_xlabel("Input pairs (millions)")
    axes[0].set_title("a  Projected peak Jaccard")
    axes[1].set_title("b  Matched-region Spearman")
    axes[0].legend(fontsize=7, loc="lower right")
    fig.suptitle(
        "Cross-assembly agreement is high among reciprocally mapped regions", fontsize=9
    )
    fig.subplots_adjust(left=0.09, right=0.98, bottom=0.21, top=0.77, wspace=0.35)
    save(fig, out, "references", rows)
    fig, ax = plt.subplots(figsize=(4.6, 2.6))
    rows = []
    for j, lib in enumerate(["SRR4000468", "SRR4000469"]):
        values = []
        for annot in [
            "araport11 lifted to Col-CEN v1.2",
            "tair10 lifted to Col-CEN v1.2",
        ]:
            for scope in ["all_transcript_starts", "common_single_start_genes"]:
                item = e.annotation["libraries"][lib]["annotations"][annot][scope]
                value = item["metric"]["value"]
                values.append(value)
                rows.append(
                    {
                        "library": lib,
                        "annotation": annot,
                        "scope": scope,
                        "tss": item["usable_TSS"],
                        "value": value,
                    }
                )
        ax.plot([0, 1, 2, 3], values, "o", color=COLORS[j], label=lib)
    ax.set_xticks(
        [0, 1, 2, 3], ["Araport\nall", "Araport\ncommon", "TAIR\nall", "TAIR\ncommon"]
    )
    ax.set_xlim(-0.4, 3.4)
    ax.set_ylabel("TSS center / flank ratio")
    ax.legend(fontsize=7)
    ax.margins(y=0.3)
    fig.suptitle(
        "The same aligned reads yield different annotation-based scores", fontsize=9
    )
    fig.subplots_adjust(left=0.14, right=0.98, bottom=0.23, top=0.8)
    save(fig, out, "annotations", rows)
