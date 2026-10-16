"""Build seven reports from frozen evidence into a new output directory."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

from evidence import HERE, Evidence, digest, tsv, verify_sources


def run(command: list[str], cwd: Path, log: Path) -> str:
    result = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    log.write_text(result.stdout)
    if result.returncode:
        raise RuntimeError(f"{command[0]} exited {result.returncode}; see {log.name}")
    return result.stdout


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        type=Path,
        required=True,
        help="New directory; existing paths are never overwritten",
    )
    parser.add_argument(
        "--font-dir",
        type=Path,
        required=True,
        help="Directory containing genuine Arial TTF/OTF files",
    )
    parser.add_argument(
        "--tex-bin",
        type=Path,
        help="Optional directory containing the pinned pdflatex toolchain",
    )
    args = parser.parse_args()
    if args.tex_bin:
        os.environ["PATH"] = (
            str(args.tex_bin.resolve()) + os.pathsep + os.environ["PATH"]
        )
    out = args.out.resolve()
    if out.exists():
        parser.error("--out must be a new directory")
    e = Evidence()  # Validate evidence before writing outputs.
    out.mkdir(parents=True)
    for name in ("figures", "data", "logs", "rendered", "cache", "tmp"):
        (out / name).mkdir()
    os.environ.update(
        {
            "MPLCONFIGDIR": str(out / "cache/matplotlib"),
            "TEXMFVAR": str(out / "cache/texmf"),
            "TEXMFCONFIG": str(out / "cache/texconfig"),
            "TMPDIR": str(out / "tmp"),
            "SOURCE_DATE_EPOCH": "1792108800",
            "FORCE_SOURCE_DATE": "1",
        }
    )
    from plots import configure, figures
    from reports import reports

    fonts = configure(args.font_dir.resolve())
    figures(e, out / "figures")
    tsv(out / "data/samples.tsv", e.samples)
    tsv(out / "data/cases.tsv", e.cases)
    tsv(out / "data/counting.tsv", e.counts)
    tsv(out / "data/resources.tsv", e.resources)
    shutil.copyfile(HERE / "sources.json", out / "data/sources.json")
    shutil.copyfile(HERE / "local-audit.json", out / "data/local-audit.json")
    style = (HERE / "style.tex").read_text()
    checks = []
    for name, body in reports(e).items():
        tex = style + body + "\\end{document}\n"
        if any(s in tex for s in ("\u2013", "\u2014")):
            raise ValueError("Use plain wording instead of typographic dashes")
        (out / f"{name}.tex").write_text(tex)
        for iteration in (1, 2):
            log = run(
                [
                    "pdflatex",
                    "-no-shell-escape",
                    "-interaction=nonstopmode",
                    "-halt-on-error",
                    name + ".tex",
                ],
                out,
                out / "logs" / f"{name}-compile-{iteration}.txt",
            )
        if re.search(
            r"Overfull|undefined references|Missing character|Font Warning", log
        ):
            raise ValueError(f"TeX layout/font warning in {name}")
        info = run(["pdfinfo", name + ".pdf"], out, out / "logs" / f"{name}-info.txt")
        pages = int(re.search(r"Pages:\s+(\d+)", info).group(1))
        expected = 3 if name in ("peak-calling", "references") else 2
        if pages != expected:
            raise ValueError(f"{name}: {pages} pages, expected {expected}")
        font_info = run(
            ["pdffonts", name + ".pdf"], out, out / "logs" / f"{name}-fonts.txt"
        )
        if name != "interval-counting" and "Arial" not in font_info:
            raise ValueError("Arial not embedded")
        run(
            ["pdftotext", "-layout", name + ".pdf", name + ".txt"],
            out,
            out / "logs" / f"{name}-text.txt",
        )
        run(
            [
                "gs",
                "-q",
                "-dSAFER",
                "-dBATCH",
                "-dNOPAUSE",
                "-sDEVICE=png16m",
                "-r60",
                f"-sOutputFile=rendered/{name}-%02d.png",
                name + ".pdf",
            ],
            out,
            out / "logs" / f"{name}-render.txt",
        )
        checks.append(
            {
                "report": name,
                "pages": pages,
                "compile_passes": 2,
                "status": "PASS",
                "visual_review": "NOT_ASSESSED; inspect rendered pages separately",
            }
        )
    artifact_paths = [
        p
        for p in out.rglob("*")
        if p.is_file()
        and not any(s in p.relative_to(out).parts for s in ("cache", "tmp"))
    ]
    tools = {
        tool: subprocess.run(
            [tool, "--version"] if tool not in ("pdfinfo",) else [tool, "-v"],
            text=True,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        ).stdout.splitlines()[0]
        for tool in ["pdflatex", "pdfinfo", "gs"]
    }
    verify_sources()
    manifest = {
        "snapshot_date": "2026-10-16",
        "python": sys.version.split()[0],
        "tools": tools,
        "fonts": fonts,
        "checks": checks,
        "source_catalog_sha256": digest(HERE / "sources.json"),
        "source_files_unchanged": True,
        "recipe": {
            p.name: digest(p)
            for p in sorted(HERE.iterdir())
            if p.suffix in (".py", ".tex", ".toml", ".lock")
        },
        "artifacts": {
            str(p.relative_to(out)): digest(p) for p in sorted(artifact_paths)
        },
        "scope": "Reports only; existing scientific evidence not rerun",
        "command": "pixi run --manifest-path benchmarks/reports/pixi.toml build --out NEW_DIRECTORY --font-dir ARIAL_DIRECTORY",
    }
    (out / "build-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(
        json.dumps(
            {
                "status": "PASS",
                "reports": len(checks),
                "pages": sum(r["pages"] for r in checks),
            }
        )
    )


if __name__ == "__main__":
    main()
