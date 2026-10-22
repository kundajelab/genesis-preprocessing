"""Parse published CLI recipes without executing network, pipelines or model deployment."""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import io
import re
import shlex
import subprocess
import sys
import tomllib
from collections.abc import Sequence
from pathlib import Path
from unittest.mock import patch
from urllib.parse import unquote

from genesis_tools.atac import run as atac

ROOT = Path(__file__).resolve().parents[1]
DOCS = [
    ROOT / "README.md",
    ROOT / "docs/README.md",
    ROOT / "docs/scientific-overview.md",
    ROOT / "docs/development.md",
    ROOT / "docs/nextflow-style.md",
    ROOT / "docs/diagrams/README.md",
    ROOT / "genesis_tools/README.md",
    *sorted((ROOT / "docs/usage").glob("*.md")),
]
ARCHIVE = ROOT / "docs/history/dapseq-readme-fb60096.md"
ARCHIVED_REVISION = "fb600968b3abf5d79756aaf7e6d1e4ff43a1b41d"
DAP_SECTIONS = (
    "install-and-run",
    "sample-sheet",
    "analysis-choices",
    "outputs-and-resources",
    "environments-and-docker-builds",
    "verification",
)


def anchors(text: str) -> set[str]:
    """Collect GitHub-style heading fragments and explicit IDs outside fenced code."""
    prose = re.sub(r"^```.*?^```[^\n]*$", "", text, flags=re.M | re.S)
    result = set(re.findall(r'<a\s+id="([^"]+)"', prose))
    for heading in re.findall(r"^#{1,6}\s+(.+?)\s*#*\s*$", prose, re.M):
        heading = re.sub(r"\[([^]]+)\]\([^)]*\)", r"\1", heading)
        base = re.sub(r"[^\w\-\s]", "", heading.lower()).replace(" ", "-")
        slug = base
        duplicate = 0
        while slug in result:
            duplicate += 1
            slug = f"{base}-{duplicate}"
        result.add(slug)
    return result


# Moved to genesis-curation in October 2026. Dated validation reports keep their original
# links; current documentation links to the new repository instead.
MOVED = (
    "benchmarks/",
    "docs/curation/",
    "docs/usage/ai-providers.md",
    "docs/usage/curation.md",
    "docs/usage/execution.md",
    "docs/usage/llm-deployment.md",
    "docs/usage/scatac.md",
    "docs/usage/scatac-model-handoff.md",
    "docs/usage/studies.md",
    "tests/study_pipeline_fixture.py",
    "tests/verify_study.py",
)


def check_links(document: Path) -> None:
    prose = re.sub(r"^```.*?^```[^\n]*$", "", document.read_text(), flags=re.M | re.S)
    historical = document.resolve().is_relative_to(ROOT / "validation")
    for target in re.findall(r"\]\(([^)]+)\)", prose):
        if "://" in target or target.startswith("mailto:"):
            continue
        relative, separator, fragment = unquote(target).partition("#")
        path = (document.parent / relative).resolve() if relative else document
        if not path.exists() and historical and path.is_relative_to(ROOT):
            if str(path.relative_to(ROOT)).startswith(MOVED):
                continue
        if not path.exists():
            raise AssertionError(f"Broken documentation link in {document.name}: {target}")
        if separator and path.suffix == ".md" and fragment not in anchors(path.read_text()):
            raise AssertionError(f"Broken documentation anchor in {document.name}: {target}")


def check_preserved_dap_readme() -> None:
    """Protect the full original text, including historically reported evidence."""
    body = ARCHIVE.read_text().split("<!-- BEGIN ORIGINAL README -->\n", 1)[1]
    body = body.split("<!-- END ORIGINAL README -->", 1)[0]
    for target in ("vendor/phantompeakqualtools/README.md", "docs/nextflow-style.md"):
        pinned = (
            f"https://github.com/kundajelab/genesis-preprocessing/blob/{ARCHIVED_REVISION}/{target}"
        )
        body = body.replace(f"]({pinned})", f"]({target})")
    assert hashlib.sha256(body.encode()).hexdigest() == (
        "3ce68770e50460992a2768f58b2c6ed9a69985f80b3ea4dcbcf4f02ae9f830ff"
    ), "Historical README text changed"
    readme = (ROOT / "README.md").read_text()
    index = (ROOT / "docs/README.md").read_text()
    guide = (ROOT / "docs/usage/dapseq-reference.md").read_text()
    assert "](docs/usage/dapseq-reference.md)" in readme, "README must link to the DAP guide"
    assert "](usage/dapseq-reference.md)" in index, "Documentation index must link to the DAP guide"
    assert "](history/dapseq-readme-fb60096.md)" in index, (
        "Historical DAP documentation must be linked"
    )
    for section in DAP_SECTIONS:
        assert section in anchors(readme), f"Original README bookmark lost: {section}"
        assert section in anchors(guide), f"DAP guide section lost: {section}"
    for path in (ROOT / "README.md", ROOT / "docs/development.md"):
        assert "https://github.com/kundajelab/genesis-preprocessing.git" in path.read_text()


class Parsed(Exception):
    pass


ORIGINAL_PARSE = argparse.ArgumentParser.parse_args


def parse_only(
    self: argparse.ArgumentParser,
    args: Sequence[str] | None = None,
    namespace: argparse.Namespace | None = None,
) -> argparse.Namespace:
    ORIGINAL_PARSE(self, args, namespace)
    raise Parsed


def main() -> None:
    check_preserved_dap_readme()
    # Historical commands describe the archived revision; check its links only.
    all_documents = sorted(
        {
            ROOT / name
            for name in subprocess.check_output(
                [
                    "git",
                    "ls-files",
                    "--cached",
                    "--others",
                    "--exclude-standard",
                    "--",
                    "*.md",
                ],
                cwd=ROOT,
                text=True,
            ).splitlines()
        }
    )
    for document in all_documents:
        check_links(document)
    tasks = tomllib.loads((ROOT / "pixi.toml").read_text())["tasks"]
    count = 0
    task_count = 0
    for document in DOCS:
        text = document.read_text()
        for block in re.findall(r"```bash\n(.*?)```", text, re.S):
            for line in block.replace("\\\n", " ").splitlines():
                if line.strip().startswith("pixi run "):
                    command = shlex.split(line, comments=True)[2]
                    # uv is an installed dependency, rather than a named Pixi task.
                    assert command in tasks or command == "uv", (document, line)
                    task_count += 1
                if not line.strip().startswith("pixi run pipeline-atac "):
                    continue
                # Parse documented ATAC options with the runner's own parser.
                entry, argv = (
                    atac.main,
                    ["pixi run pipeline-atac", *shlex.split(line, comments=True)[3:]],
                )
                with (
                    patch.object(sys, "argv", argv),
                    patch.object(argparse.ArgumentParser, "parse_args", parse_only),
                    contextlib.redirect_stdout(io.StringIO()),
                ):
                    try:
                        entry()
                    except Parsed:
                        count += 1
                    except SystemExit as error:
                        if error.code not in (None, 0):
                            raise AssertionError(f"Invalid documented command: {line}") from error
                        count += 1
                    else:
                        raise AssertionError("Documentation parser reached an executable action")
    print(
        f"PASS {count} documented CLI commands parse; {task_count} Pixi recipes resolve; "
        f"links/anchors in {len(all_documents)} documents resolve; original DAP README preserved"
    )


if __name__ == "__main__":
    main()
