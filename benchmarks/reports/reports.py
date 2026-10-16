"""Short reports assembled from verified measurements; no pipeline execution."""

from __future__ import annotations

from typing import Any

from evidence import COHORTS, Evidence


def escape(value: Any) -> str:
    mapping = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(mapping.get(c, c) for c in str(value))


def fmt(value: float | None, digits: int = 4) -> str:
    return "nd" if value is None else f"{value:,.{digits}f}"


def table(
    headers: list[str], rows: list[list[Any]], caption: str, widths: str | None = None
) -> str:
    layout = widths or ("l" + "r" * (len(headers) - 1))
    return (
        r"\begin{table}[H]\centering\small\caption{"
        + escape(caption)
        + "}\n"
        + r"\begin{tabular}{"
        + layout
        + r"}\toprule"
        + "\n"
        + " & ".join(escape(v) for v in headers)
        + r" \\ \midrule"
        + "\n"
        + "\n".join(" & ".join(escape(v) for v in row) + r" \\" for row in rows)
        + "\n"
        + r"\bottomrule\end{tabular}\end{table}"
        + "\n"
    )


def figure(name: str, caption: str, wide: bool = True) -> str:
    width = "7in" if wide else "4.6in"
    return (
        r"\begin{figure}[H]\centering\makebox[\linewidth][c]{\includegraphics[width="
        + width
        + "]{figures/"
        + name
        + ".pdf}}\n"
        + r"\caption{"
        + escape(caption)
        + r"}\end{figure}"
        + "\n"
    )


def heading(title: str, numbered: bool = True) -> str:
    return "\\section" + ("" if numbered else "*") + "{" + escape(title) + "}\n"


def para(text: str) -> str:
    return escape(text) + "\n\n"


def start(title: str, subtitle: str, abstract: str) -> str:
    return (
        "\\reporttitle{"
        + escape(title)
        + "}{"
        + escape(subtitle)
        + "}\n"
        + r"\begin{abstract}"
        + escape(abstract)
        + r"\end{abstract}"
        + "\n"
        + heading("Setup", False)
    )


def setup(label: str, text: str) -> str:
    return r"\textbf{" + escape(label) + ".} " + para(text)


METHODS = "We summarize recorded results without rerunning biological processing. The main Arabidopsis matrix uses seeded nested raw-template subsets (seed 60141), three fresh technical repetitions at the 1-million cap and one at the 5-million cap. Actual counts replace requested caps for short libraries. DAP retains primary mapped reads including MAPQ 0; ATAC selects proper primary pairs with both mates MAPQ at least 30, excludes MAPQ 255, organelles and marked duplicates. Biological identity is preserved. Source files are SHA-256 checked against the accompanying catalog; evidence snapshot revision ed1c50b8e4b5 includes the October 13 to 16 records. Completed local reference assessments are separately identified in that catalog."


def end(
    issues: list[str], recommendation: str, steps: list[str], methods: str = METHODS
) -> str:
    return (
        heading("Open issues", False)
        + r"\begin{itemize}"
        + "\n"
        + "".join(r"\item " + escape(s) + "\n" for s in issues)
        + r"\end{itemize}"
        + "\n"
        + heading("Recommendation", False)
        + para(recommendation)
        + r"\begin{enumerate}"
        + "\n"
        + "".join(r"\item " + escape(s) + "\n" for s in steps)
        + r"\end{enumerate}"
        + "\n"
        + heading("Appendix: Methods", False)
        + para(methods)
    )


def samples(e: Evidence, cohort: str) -> str:
    return table(
        ["Library", "Role / control", "Raw templates"],
        [
            [
                r["id"],
                r.get("control")
                or (
                    "biological replicate " + str(r["biological_replicate"])
                    if r["assay_subtype"] == "bulk-ATAC"
                    else "Control"
                ),
                f"{r['raw_templates']:,}",
            ]
            for r in e.samples
            if r["cohort"] == cohort
        ],
        f"Input identities at the {cohort} cap.",
    )


def quality_table(
    e: Evidence,
    cohort: str,
    methods: list[str],
    libraries: list[str],
    metrics: list[tuple[str, str, int]],
) -> str:
    rows = []
    for method in methods:
        for lib in libraries:
            rows.append(
                [method, lib]
                + [
                    fmt(e.q(cohort, method, lib, key), digits)
                    for key, _, digits in metrics
                ]
            )
    return table(
        ["Method", "Library"] + [label for _, label, _ in metrics],
        rows,
        "Recorded median measurements; nd means not determined.",
    )


def reports(e: Evidence) -> dict[str, str]:
    result = {}
    libs = ["SRR4000468", "SRR4000469"]
    cases = len(e.cases)
    # Overview: two pages.
    s = start(
        "Plant preprocessing benchmark evidence",
        "Arabidopsis DAP-seq and bulk ATAC; matched-stage comparisons",
        f"We reconcile {cases} completed method/repetition cases across 6 public-data subsets. In a separate MACS3 test, all 24 CPU-specific comparisons preserved the evaluated peak outputs, with 1-CPU speedups of 4.20 to 6.32 times. Most package alternatives lack evidence of a consistent biological improvement, and the broader species panel remains incomplete. We would retain scientific defaults while evaluating the pinned optimized caller within its tested scope.",
    )
    s += setup(
        "Data",
        "Two SE DAP treatments, four PE DAP treatments and two biological ATAC replicates were sampled at 1-million and 5-million template caps. DAP treatments are not extra biological replicates; repeated execution is not a new experiment.",
    )
    s += setup(
        "Scoring",
        "Peak Jaccard is intersecting divided by union peak bases. Held-out FRiP is the fraction of retained fragments overlapping the other biological replicate's peaks. Fixed-region Spearman measures count agreement, not truth. Independent PBM motif support uses matched windows; cost is secondary to these scientific comparisons.",
    )
    s += table(
        ["Data", "Treatment / replicate identity", "Cases"],
        [
            ["SE DAP", "FUS3, ABI5; assigned control", 16],
            ["PE DAP", "Four protocol-optimization treatments", 16],
            ["Bulk ATAC", "SRR4000468 and SRR4000469", 32],
        ],
        "Completed method/repetition cases, not biological sample counts.",
        "lp{3.4in}r",
    )
    s += heading("Depth changes discovery more than these small package differences")
    s += para(
        "For the 2 ATAC biological replicates, peak Jaccard increases from 0.2404 at 1 million pairs to 0.4858 at 5 million pairs. Held-out FRiP rises from 0.1433 to 0.4601 and from 0.1633 to 0.4023, while peak coverage also increases."
    )
    s += figure("depth", "TAIR10 baseline peak discovery and cross-replicate overlap.")
    s += r"\newpage" + "\n" + heading("The current evidence supports scoped decisions")
    s += para(
        f"The {cases} cases reconcile to the historical matrix count, but do not represent {cases} independent biological experiments. A separate local audit found 8 completed frozen reference baselines and 3 completed assessments; these are not added to the method-matrix denominator."
    )
    s += table(
        ["Step", "Observed evidence", "Decision"],
        [
            [
                "Alignment",
                "BWA is close; Bowtie2 changes PE DAP peaks",
                "Retain BWA-MEM2",
            ],
            [
                "Trimming",
                "Two short-read ATAC replicates change little",
                "Defer general change",
            ],
            ["Duplicates", "Picard is faster in one ATAC cohort", "Broaden validation"],
            ["Peaks", "24 exact-output optimized MACS3 comparisons", "Scoped opt-in"],
            [
                "Counting",
                "Independent read-end arithmetic agrees",
                "Retain implementation",
            ],
            [
                "References",
                "High eligible-region concordance; exclusions remain",
                "Choose explicitly",
            ],
        ],
        "Decisions supported by the reports in this series.",
        "p{.85in}p{3.1in}p{1.45in}",
    )
    s += table(
        ["Step / claim", "Remaining evidence"],
        [
            [
                "QC and signal-generation packages",
                "No completed matched-package ranking in this series",
            ],
            [
                "Other plant species",
                "Acquisition or mapping checks do not establish method performance",
            ],
            [
                "Model training or sequence design",
                "No benchmark here establishes predictive or experimental benefit",
            ],
        ],
        "Coverage gaps are distinct from failed comparisons.",
        "p{2in}p{3.55in}",
    )
    s += end(
        [
            "Subsets are not full-study validation; short libraries may reach full yield.",
            "Baseline agreement, own FRiP and peak counts are not ground truth.",
            "The nonexclusive host and CPU allocations limit timing transfer.",
        ],
        "We would use these reports to select the next acceptance tests, without changing scientific defaults.",
        [
            "Evaluate candidate changes on longer reads and additional plant genomes.",
            "Require independent support and biological reproducibility before selecting a new default.",
        ],
    )
    result["overview"] = s
    # Alignment.
    s = start(
        "Alignment package comparison",
        "Arabidopsis SE/PE DAP and bulk ATAC; BWA-MEM2, BWA and Bowtie2",
        "At the 1-million-pair cap, Bowtie2 gives peak-base Jaccard values of 0.5611 to 0.8493 against BWA-MEM2 across four PE DAP treatments. At the deeper SE DAP cap, FUS3 motif support is 466/998 baseline windows versus 465/998 Bowtie2 windows. These measurements describe agreement and independent motif support, not a validated accuracy ranking. We would retain BWA-MEM2 pending broader evidence.",
    )
    s += setup(
        "Comparisons",
        "BWA-MEM2 2.3, BWA 0.7.19 and Bowtie2 2.5.5 use the same sampled reads and declared reference. The BWA-MEM2 package reports an older banner; recorded binary identity establishes the package version. Bowtie2 uses very-sensitive-local and a 2,000-base maximum PE insert; ATAC also allows dovetail and disables mixed/discordant pairs.",
    )
    s += setup(
        "Scoring",
        "Peak Jaccard is base-pair intersection/union against baseline peaks. Fixed-region Spearman compares fragment counts on identical baseline regions. PBM support counts motif-positive 200-base summit windows with same-chromosome, GC-matched controls. These operating policies are not a pure binary substitution.",
    )
    s += samples(e, "dap-pe-1m")
    s += heading("Bowtie2 changes peak agreement unevenly across treatments")
    s += para(
        "Across 4 PE treatments, fixed-region Spearman spans 0.8793 to 0.9269 for Bowtie2. The 2 lower-peak-count libraries have Jaccard below 0.58, compared with more than 0.82 for the other 2 treatments."
    )
    s += figure(
        "alignment-quality",
        "Peak and fixed-region agreement at the 1-million-pair cap.",
    )
    s += (
        r"\newpage"
        + "\n"
        + heading("Independent motif support does not favor a general aligner change")
    )
    motifrows = [
        r
        for r in e.motifs
        if r["selection"] == "top1000" and r["arm"] in ("baseline", "bwa", "bowtie2")
    ]
    s += para(
        "The deeper SE DAP comparison covers 2 transcription factors, not 1,998 independent experiments. FUS3 has 2 unmatched control windows, so its top-ranked denominator is 998 rather than 1,000."
    )
    s += table(
        ["Library", "Method", "Motif hits / windows", "Control hits"],
        [
            [
                r["library"],
                r["arm"],
                r["peak_hits"] + "/" + r["paired_windows"],
                r["control_hits"],
            ]
            for r in motifrows
        ],
        "Independent PBM support in matched top-ranked windows.",
    )
    s += figure(
        "alignment-cost",
        "Fresh PE DAP cohort timings at the 1-million-pair cap.",
        False,
    )
    s += end(
        [
            "One species and selected TFs cannot establish a universal aligner ranking.",
            "Bowtie2 local alignment and pairing policies differ from BWA policies.",
            "Timings cover all cohort libraries, exclude indexing and evaluation, and allow 14 tool CPUs.",
        ],
        "We would retain BWA-MEM2 and investigate disagreement on held-out biological evidence.",
        [
            "Expand TF and read-length coverage.",
            "Inspect alternative-specific peaks before accepting greater retention as improvement.",
        ],
    )
    result["alignment"] = s
    # Trimming.
    s = start(
        "Adapter trimming comparison",
        "Two Arabidopsis bulk ATAC biological replicates; paired input",
        "At 1 million pairs, untrimmed input retains 458,813 and 565,032 eligible fragments. Cutadapt and fastp alter these totals by at most 62 fragments each, and peak Jaccard remains above 0.9995. The inspected archived reads are only 31 bases long, limiting generalization to longer libraries. We would retain explicit protocol-specific trimming decisions.",
    )
    s += setup(
        "Data",
        "SRR4000468 and SRR4000469 are biological replicates from seven-day Arabidopsis seedlings. The archived first 100 records per mate are 31 bases long; this is an inspection result, not a full-library length census.",
    )
    s += setup(
        "Comparisons",
        "No trimming, Cutadapt 5.2 and fastp 1.4.0 use the explicit Nextera adapter CTGTCTCTTATACACATCT, with unrelated fastp filters disabled. All arms keep the same reference, alignment, duplicate policy and MACS3 settings.",
    )
    s += setup(
        "Scoring",
        "Retention is eligible nuclear fragments divided by raw input pairs, not the fraction exiting the trimmer. Peak Jaccard is base-pair intersection/union; fixed-region Spearman uses identical baseline regions. Output-pair retention from trimming alone is not inferred from downstream retention.",
    )
    s += samples(e, "atac-arabidopsis-1m")
    s += heading("Trimming produces small changes in these short-read libraries")
    s += para(
        "The 2 Cutadapt fixed-region correlations exceed 0.99998. Fastp reports 4,995 and 3,042 peaks, compared with 4,997 and 3,043 for untrimmed input."
    )
    s += figure(
        "trimming-quality", "Downstream fragment retention and peak-base agreement."
    )
    s += r"\newpage" + "\n" + heading("The processing times are close on this cohort")
    rows = [e.c("atac-arabidopsis-1m", m) for m in ("baseline", "cutadapt", "fastp")]
    s += para(
        "Across 3 fresh technical repetitions per arm, cohort timing includes both biological libraries and subsequent processing. These are 14-CPU-budget runs; they do not isolate trimmer execution time or establish a throughput ranking for long reads."
    )
    s += figure(
        "trimming-cost",
        "Cohort wall time by trimming policy; ranges are observed.",
        False,
    )
    s += quality_table(
        e,
        "atac-arabidopsis-1m",
        ["baseline", "cutadapt", "fastp"],
        libs,
        [("fixed_count_spearman", "Spearman", 6)],
    )
    s += end(
        [
            "The first-100-read inspection cannot prove that every read is adapter-free.",
            "Long-read and adapter-rich libraries need separate comparisons.",
            "Standalone trimming pair-loss and phase-specific speed are not inferred from downstream totals.",
        ],
        "We would defer a universal trimming change and retain declared adapter policies.",
        [
            "Compare longer paired reads with measured adapter burden.",
            "Report trimmer output pairs separately from alignment and fragment selection.",
        ],
    )
    result["trimming"] = s
    # Duplicates.
    s = start(
        "Duplicate marker comparison",
        "Arabidopsis bulk ATAC; samtools and Picard",
        "At 1 million pairs per library, samtools and Picard retain identical counts in one replicate and differ by one fragment in the other. Both produce peak Jaccard of 1.0 against the samtools baseline. Median cohort time is 289.3 seconds with samtools and 243.2 seconds with Picard, but this is one short-read cohort. We would broaden validation before changing the default marker.",
    )
    s += setup(
        "Comparisons",
        "samtools 1.24 and Picard 3.5.0 mark the same aligned inputs, followed by the same proper-pair, MAPQ and organellar selection. Duplicate-marked fragments are excluded in both arms; DAP duplicate retention is a separate scientific policy.",
    )
    s += setup(
        "Scoring",
        "Retention is selected fragments/raw pairs. Jaccard uses peak-base intersection/union. Fixed-region Spearman compares counts on baseline regions. Post-filter duplicate flags cannot estimate pre-filter PCR duplication; pre-exclusion coordinate complexity is a separate measurement.",
    )
    s += samples(e, "atac-arabidopsis-1m")
    s += heading("The two markers preserve nearly identical selected data")
    s += para(
        "Both arms retain 458,813 fragments in SRR4000468; SRR4000469 retains 565,032 with samtools and 565,033 with Picard. Fixed-region Spearman is 1.0 and 0.999998 for the 2 Picard comparisons."
    )
    s += figure(
        "duplicates-quality",
        "Retained fragments and peak agreement after duplicate exclusion.",
    )
    s += (
        r"\newpage"
        + "\n"
        + heading("Picard lowers observed cohort time in this dataset")
    )
    s += para(
        "The median cohort reduction is 46.0 seconds, or 15.9%, over 3 fresh repetitions per arm. The observed samtools range is 273.2 to 293.7 seconds, versus 238.1 to 248.2 seconds for Picard."
    )
    s += figure(
        "duplicates-cost",
        "Whole-cohort processing time, excluding indexing and evaluation.",
        False,
    )
    s += table(
        ["Measurement", "1M pairs: replicate 1 / 2", "5M pairs: replicate 1 / 2"],
        [
            ["Baseline NRF", "0.99639 / 0.99520", "0.98273 / 0.97605"],
            ["Raw pairs", "1,000,000 / 1,000,000", "5,000,000 / 5,000,000"],
        ],
        "Pre-exclusion NRF uses distinct outer-coordinate groups / eligible fragments.",
        "lp{2in}p{2in}",
    )
    s += end(
        [
            "Cohort time is not the isolated marking phase and includes common processing.",
            "The larger maize-genome comparison is not established by these records.",
            "Outer-fragment complexity differs from unclipped-start definitions used elsewhere.",
        ],
        "We would treat Picard as a performance candidate, with the duplicate-exclusion policy unchanged.",
        [
            "Compare longer reads and larger genomes under matched resource limits.",
            "Inspect tie-breaking and optical-duplicate behavior before a default change.",
        ],
    )
    result["duplicates"] = s
    # Peak calling: three pages.
    s = start(
        "Peak calling package comparison",
        "Arabidopsis DAP and ATAC; caller policies and optimized MACS3",
        "MACS2 and MACS3 give similar independent motif support in two SE DAP treatments, while Genrich produces zero peaks at the tested q = 0.01 operating point. A separate optimized-MACS3 test preserves outputs in all 24 comparisons and reduces wall time by 4.20 to 42.79 times across CPU settings. Genrich thresholds are not statistically interchangeable with MACS thresholds, and parallel MACS3 uses more sampled memory. We would retain caller policies and support the pinned optimization within the tested scope.",
    )
    s += setup(
        "Data",
        "The method matrix uses 1-million and 5-million raw-template caps. The separate optimization uses 8 existing treatment/replicate input cases across SE DAP, PE DAP and ATAC, with fixed controls and input hashes.",
    )
    s += setup(
        "Scoring",
        "Peak discovery is reported with union coverage and biological replicate overlap, not treated as accuracy. Held-out FRiP uses the other replicate's peaks. Exact optimized-output equivalence compares narrowPeak/summit bytes and peak-table content excluding header comments.",
    )
    s += samples(e, "dap-se-5m")
    s += heading("Greater depth improves measured ATAC replicate agreement")
    s += para(
        "Across 2 biological replicates, TAIR10 peak-base Jaccard rises from 0.2404 to 0.4858 between 1 million and 5 million pairs. Peak discovery rises from 4,997 to 14,805 and from 3,043 to 11,551, so greater held-out FRiP is interpreted alongside expanded peak coverage."
    )
    s += figure(
        "depth", "Depth-dependent peak discovery and held-out fragment overlap."
    )
    s += (
        r"\newpage"
        + "\n"
        + heading("Caller operating points prevent a simple accuracy ranking")
    )
    s += para(
        "At q = 0.01, Genrich reports 0 peaks for both ATAC replicates at both tested depths. Diagnostic p = 0.01 calls are a different operating point, not a rescue at a matched false-discovery rate."
    )
    s += quality_table(
        e,
        "atac-arabidopsis-1m",
        ["baseline", "macs2", "genrich"],
        libs,
        [
            ("peak_count", "Peaks", 0),
            ("peak_bases", "Union bp", 0),
            ("fragments_FRiP_own", "Own FRiP", 4),
        ],
    )
    s += table(
        ["Library", "Caller", "PBM hits / matched windows"],
        [
            [r["library"], r["arm"], r["peak_hits"] + "/" + r["paired_windows"]]
            for r in e.motifs
            if r["selection"] == "top1000" and r["arm"] in ("baseline", "macs2")
        ],
        "Independent motif support at the deeper SE DAP cap.",
    )
    s += setup(
        "Policies",
        "MACS2 2.2.9.1 and MACS3 3.0.5 use BAM for SE and BAMPE for PE, with assigned DAP controls and keep-dup = 1; ATAC retains all already-filtered fragments. Genrich 0.6.2 uses ATAC cut-site mode, 100-base expansion, AUC 200 and gap 100. Its background includes header contigs, whereas the ATAC MACS effective size uses nuclear ACGT bases.",
    )
    s += setup(
        "Reproducibility",
        "Exploratory pairwise IDR uses permissive MACS p = 0.01 calls, signal-value ranking and seed 60141. There are no pooled/pseudoreplicate rescue ratios or approved plant exclusion thresholds; this is not an ENCODE acceptance claim.",
    )
    s += heading("The optimization preserves the evaluated peak outputs")
    s += para(
        "All 24 CPU-specific comparisons pass output equivalence across 8 inputs. The 192 invocations comprise 48 warmups and 144 measured calls, with 3 measured repetitions per build and CPU setting."
    )
    s += table(
        ["Allocated CPUs", "Wall-time speedup range", "RSS direction"],
        [
            [
                str(cpu),
                f"{min(float(r['speedup']) for r in e.macs if int(r['cpus']) == cpu):.2f} to {max(float(r['speedup']) for r in e.macs if int(r['cpus']) == cpu):.2f}",
                "Lower in all 8"
                if cpu == 1
                else ("Higher in all 8" if cpu == 8 else "Case-dependent"),
            ]
            for cpu in (1, 2, 8)
        ],
        "Paired optimized/released comparisons on identical input files.",
    )
    s += (
        r"\newpage"
        + "\n"
        + heading("Parallel speed gains require additional sampled memory")
    )
    s += para(
        "For the largest tested PE DAP case at 8 CPUs, median wall time falls from 84.29 to 1.97 seconds while sampled group RSS rises from 639 to 1,330 MiB. At 1 CPU, all 8 inputs instead show lower sampled RSS with the optimized build."
    )
    s += figure(
        "macs-speed", "Paired speed and sampled-memory ratios across eight input cases."
    )
    s += setup(
        "Build identity",
        "The released and optimized binaries both report 3.0.5. The optimized source is pinned to 22eb40f7b51f5b20b635614ca0b1759dc5398c5c, with matched Python and numerical dependencies; the release number alone cannot identify this optimization. Measurements use a nonexclusive x86-64 Xeon host, CPU affinity, one warmup and three measured runs.",
    )
    s += end(
        [
            "Sampled process-group RSS can double-count shared pages and is not host peak RAM.",
            "The PE DAP 5M treatment cap contains 1,657,793 actual pairs; its control contains 5M.",
            "The SE DAP 5M control contains 4,241,389 reads, not 5M.",
            "The upstream suite has one call-variant failure also seen in the release; it is not fully green.",
        ],
        "We would support the pinned optimized caller as a scoped opt-in and retain scientific peak-calling policies.",
        [
            "Validate the published deployment image and representative full-depth inputs before promotion.",
            "Compare callers using independent support and biological reproducibility at declared operating points.",
        ],
    )
    result["peak-calling"] = s
    # Counting two pages, no manufactured runtime comparison.
    ok = [r for r in e.counts if r["status"] == "OK"]
    total = sum(int(r["denominator"]) for r in ok)
    s = start(
        "Interval counting verification",
        "Genesis/Python read-end overlap counts and bedtools coverage",
        f"Across {len(ok)} evaluated library/method/repetition records, independent read-end counts have zero mismatching regions. The records contain {total:,} region evaluations, including repeated regions across methods and technical runs. This verifies arithmetic under the declared overlap convention, not biological accuracy or relative speed. We would retain the current counter and its independent checks.",
    )
    s += setup(
        "Comparisons",
        "The benchmark counts retained primary read-end spans in zero-based half-open peak intervals using Python and bedtools 2.31.1 coverage. Paired ends count independently. Fragment-level FRiP is a different quantity and is not used as a counting-equivalence target.",
    )
    s += setup(
        "Scoring",
        "A mismatch is a region whose integer overlap counts differ. Any shared base counts as overlap; no fractional NH weighting is introduced. Denominators count evaluated region instances, not distinct regulatory elements or biological replicates.",
    )
    s += table(
        ["Cohort", "Records", "Region evaluations", "Mismatches"],
        [
            [
                c,
                len([r for r in ok if r["cohort"] == c]),
                f"{sum(int(r['denominator']) for r in ok if r['cohort'] == c):,}",
                f"{sum(int(r['value']) for r in ok if r['cohort'] == c):,}",
            ]
            for c in COHORTS
        ],
        "Read-end counting comparisons; repeated region evaluations are retained.",
    )
    s += heading("The evaluated integer overlap counts agree exactly")
    s += para(
        f"The {len(ok)} evaluated records report 0 mismatches over {total:,} region instances. This denominator includes the 1-million and 5-million caps, multiple method arms, and repeated execution; none increases the number of biological samples."
    )
    s += table(
        ["Quantity", "Unit / denominator", "Equivalent?"],
        [
            ["Read-end overlap count", "Integer aligned ends per region", "Tested"],
            [
                "Fragment FRiP",
                "Fragments overlapping peaks / retained fragments",
                "Different target",
            ],
            [
                "Peak RPM",
                "Read-end count / retained ends, scaled by 1 million",
                "Normalization separate",
            ],
            [
                "Mean coverage RPKM",
                "Length-weighted coverage, including zero bases",
                "Different target",
            ],
        ],
        "Measurement conventions distinguish otherwise similar interval files.",
        "p{1.25in}p{3.1in}p{1in}",
    )
    s += (
        r"\newpage"
        + "\n"
        + heading("The evidence does not establish a counting speed advantage")
    )
    s += para(
        "The independent bedtools evaluation runs after processing, whereas the 64 method-case processing timings exclude common evaluation. There are 0 matched standalone Python-versus-bedtools timing series in this report, so relative speed and memory are nd."
    )
    s += table(
        ["Measurement", "Python", "bedtools"],
        [
            ["Matching read-end arithmetic", "PASS", "PASS"],
            ["Matched standalone runtime", "nd", "nd"],
            ["Matched standalone peak memory", "nd", "nd"],
        ],
        "nd means not determined in a matched counting-only comparison.",
    )
    s += setup(
        "Missing values",
        "A missing or unevaluated count is not zero. Any unavailable records remain in the exported counting table with their status and reason; the denominator above includes only records explicitly marked OK.",
    )
    s += setup(
        "Publication scope",
        "These checks cover the benchmark's independent integer read-end overlap calculation. They do not certify all production quantification modes, NH fractional weighting, spliced-block conventions, or coverage integration by transitivity.",
    )
    s += end(
        [
            "Repeated evaluation of a fixed peak set is not independent biological validation.",
            "Span-based overlap should not be generalized to spliced exon-block counting.",
            "No matched standalone runtime advantage is established.",
        ],
        "We would retain exact arithmetic tests and avoid choosing a counter by unmatched timings.",
        [
            "Define a counting-only test with identical interval ordering, BAM cache conditions and units.",
            "Measure time and process memory separately from conversion and file staging.",
        ],
    )
    result["interval-counting"] = s
    # References three pages.
    s = start(
        "Reference and annotation comparison",
        "Arabidopsis TAIR10 and Col-CEN; fixed inputs and explicit coordinate mapping",
        "At 1 million ATAC pairs, reciprocal peak Jaccard is 0.8984 and 0.9639 between TAIR10 and Col-CEN for two replicates. A completed local 5-million-pair assessment gives 0.9547 and 0.9571 among eligible mappings. Annotation alone changes TSS scores on identical aligned reads, and mapping exclusions limit every cross-assembly claim. We would choose reference and annotation explicitly rather than rank them by peak count or TSS enrichment.",
    )
    s += setup(
        "Data",
        "The same paired ATAC subsets and both-mate MAPQ30 selection are used for each assembly, with duplicate and organellar exclusion. Each reference supplies its own effective genome size; peak calling therefore also reflects a changed background model.",
    )
    s += setup(
        "Scoring",
        "Reciprocal conversion requires at least 95% aligned bases in each direction, unique placement and at least 95% round-trip span overlap. Jaccard uses eligible peak-base unions in the destination assembly. Spearman compares fragment counts on baseline intervals and their eligible projections; widths may differ.",
    )
    s += samples(e, "atac-arabidopsis-1m")
    s += heading("Eligible regions agree strongly across the two assemblies")
    s += para(
        "At 1 million pairs, matched-region Spearman is 0.9979 and 0.9946 for 4,992 and 3,040 eligible TAIR10 peaks. At 5 million pairs, it is 0.9990 and 0.9977; these are existing recorded results, not newly rerun analyses."
    )
    s += figure("references", "Reciprocal-coordinate agreement at two read depths.")
    s += r"\newpage" + "\n" + heading("Coordinate exclusions remain part of the result")
    rows = []
    for depth, source in [(1, e.reference), (5, e.local["reference_5m"])]:
        for lib, row in source["libraries"].items():
            rows.append(
                [depth, lib]
                + [
                    f"{row[side]['counts']['UNIQUE_RECIPROCAL']:,} / {row[side]['input_peaks']:,}"
                    for side in ["left_mapping", "right_mapping"]
                ]
            )
    s += table(
        [
            "Million pairs",
            "Library",
            "TAIR10 accepted / total",
            "Col-CEN accepted / total",
        ],
        rows,
        "Peak eligibility before cross-assembly comparison.",
    )
    s += para(
        "At the 1-million-pair cap, 8 TAIR10 and 19 Col-CEN peaks across the 2 libraries fail the reciprocal criteria. Excluded, ambiguous or unmapped intervals are not assigned to a convenient first placement."
    )
    s += setup(
        "Evidence revision",
        "The public October 14 snapshot described the 5-million-pair mapping as pending. Its later completed assessment now has a matching SHA-256 completion receipt. The accompanying compact audit preserves the result, mapping, script and experiment hashes; it does not overwrite the historical report or distribute the local BAMs.",
    )
    s += heading("Reference distribution names do not guarantee new sequence")
    s += para(
        "The A. lyrata audit reports identical case-insensitive sequence multisets for 695 contigs totaling 206,667,935 bases in the Ensembl v1.0 and Phytozome v2.1 distributions. All 695 contigs have unique full-sequence matches, while names and order differ."
    )
    s += table(
        ["Comparison", "Completed evidence", "Not established"],
        [
            [
                "A. lyrata distributions",
                "695 exact contig matches",
                "Newer sequence improves peaks",
            ],
            [
                "Maize AGPv4 / NAM5",
                "1,003 forward / 2,413 reverse chains checked",
                "Sequencing method superiority",
            ],
            [
                "Sorghum assemblies",
                "Reference acquisition records",
                "Validated peak concordance",
            ],
        ],
        "Reference preparation is distinct from a sequencing comparison.",
        "p{1.35in}p{2.2in}p{1.85in}",
    )
    s += setup(
        "Limits",
        "Uniqueness is bounded by the supplied assembly alignments, not a proof that every repetitive placement is represented. Comparison within retained coordinates does not establish overall assembly accuracy, and neither reference has independent biological ground truth here.",
    )
    s += (
        r"\newpage"
        + "\n"
        + heading("Annotation changes TSS enrichment without changing reads")
    )
    s += para(
        "On identical Col-CEN BAMs, all-start TSS sets contain 46,558 Araport and 32,986 TAIR positions. For replicate SRR4000469, the score is 2.7517 with Araport versus 3.2096 with TAIR; the common-gene comparison retains 18,428 eligible starts per annotation."
    )
    s += figure(
        "annotations",
        "TSS enrichment with all-start and common-gene annotation sets.",
        False,
    )
    s += setup(
        "Metric",
        "Tn5 cut positions use zero-based start+4 and end-6. The enrichment ratio divides the center bin by mean outer 100-base flanks across full 2,000-base windows on each side. Of 18,430 common single-start genes, 2 edge positions lack full windows; common genes can have different TSS positions across annotations.",
    )
    s += end(
        [
            "Higher TSS enrichment does not prove better annotation or library quality.",
            "Matched-coordinate agreement excludes precisely some difficult assembly regions.",
            "The broader species sequencing panel is not established by reference acquisition alone.",
        ],
        "We would record sequence identity and annotation suitability as separate decisions.",
        [
            "Choose the assembly and gene models required by the intended downstream task.",
            "Retain coordinate exclusions and compare independent regulatory evidence.",
        ],
    )
    result["references"] = s
    return result
