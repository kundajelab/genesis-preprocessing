# Report validation

Snapshot date: October 16, 2026. Scientific evidence baseline:
`ed1c50b8e4b55eb4a87e58f27120da62368c5509`.

This is a report-generation check, not a rerun of the sequencing experiments.
The [build manifest](build-manifest.json) records Python/tool versions, font
identities, recipe hashes, source-catalog hash, and individual artifact hashes.
[The visual-review receipt](review.json) binds the inspected pages and published PDFs
to their hashes. [The source catalog](data/sources.json) records input hashes and source revisions.

| Check | Result | Scope |
| --- | --- | --- |
| Input provenance | PASS | Frozen source files and six archive publication hashes; checked before and after generation |
| Case reconciliation | PASS | 64 complete, non-resumed method/repetition records; duplicates rejected |
| Focused report tests | PASS | 10 tests, including negative provenance/font/metric tests |
| Typesetting | PASS | Seven PDFs, each compiled twice with pdflatex; 16 total pages |
| TeX warnings | PASS | No overfull boxes, undefined references, missing glyphs or font warnings |
| Documentation | PASS | 144 CLI examples, 28 Pixi recipes and links/anchors in 85 documents |
| Lint and formatting | PASS | Ruff checks for the reporting code |
| Figure bounds | PASS | Rendered figure text checked against the canvas bounds |
| Fonts | PASS | Arial embedded in figures; Times-family text through the requested LaTeX `times` package |
| Visual review | PASS | All 16 rendered pages inspected at 60 dpi; revised overlapping reference labels reviewed again |
| Scientific pipeline rerun | NOT RUN | No new alignment, peak calling, model training or full-data processing |
| Full Genesis regression | NOT RUN | Reporting-only change; processing code and defaults untouched |

Commands, executed from the checkout root (local path arguments are represented
by descriptive placeholders):

```bash
pixi run --manifest-path benchmarks/reports/pixi.toml build \
  --out NEW_LOCAL_OUTPUT --font-dir LOCAL_ARIAL_DIRECTORY --tex-bin LOCAL_TEX_BIN
pixi run --manifest-path benchmarks/reports/pixi.toml test
ruff check benchmarks/reports
ruff format --check benchmarks/reports
PYTHONPATH=genesis_tools/src python tests/verify_readme.py
```

The successful build returned exit 0. It checks actual input counts, preserves
missing metrics, and rejects altered source hashes. Initial development builds
exposed an incomplete Conda TeX installation, page overflow, and crowded figure
labels; those builds were not published. Final PDFs use the tested local TeX
Live toolchain described in the build instructions. Fonts and toolchain files
remain outside Git.

The compact later reference assessment was verified against its local completion
receipt. Its full original result and BAMs are not published in this collection;
its hashes and coordinate exclusions remain available for a separate audit.
The source evidence is recorded evidence, not an independent replication.
