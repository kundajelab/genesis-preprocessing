"""Focused, dependency-free checks for the Nextflow/environment style refactor."""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pysam
from genesis_tools.quantification import quantify, quantify_peaks

ROOT = Path(__file__).resolve().parents[1]

FAKE_COMMAND = r"""
import json
import os
from pathlib import Path
import sys

name = Path(sys.argv[0]).name
args = sys.argv[1:]
with Path(os.environ["BUILD_LOG"]).open("a") as stream:
    stream.write(json.dumps([name, args]) + "\n")
if name == "git":
    print("0.1.0")
elif name == "pixi":
    if args[0] == "init":
        folder = Path(args[-1])
        folder.mkdir(exist_ok=True)
        (folder / "pixi.toml").write_text("mock manifest\n")
    elif args[0] == "lock":
        manifest = Path(args[args.index("--manifest-path") + 1])
        if manifest.is_dir():
            manifest = manifest / "pixi.toml"
        manifest.with_name("pixi.lock").write_text("mock lock\n")
elif name == "docker":
    if args[0] == "build":
        if os.environ.get("FAIL_BUILD") and "ENV_NAME=GENESIS_ALIGNMENT" in args:
            sys.exit(7)
    elif args[0] == "images":
        print("REPOSITORY TAG DIGEST IMAGE_ID CREATED SIZE")
        print(args[-1].split(":")[0] + " test sha256:" + "a" * 64 + " fake now 1MB")
"""


def write_bam(path: Path, records: str) -> None:
    """Write SAM records under a fixed header to a BAM, as the aligner's output would be."""
    sam = path.with_suffix(".sam")
    header = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chrUn\tLN:500\n"
    sam.write_text(header + records)
    with (
        pysam.AlignmentFile(str(sam)) as source,
        pysam.AlignmentFile(str(path), "wb", template=source) as destination,
    ):
        for alignment in source:
            destination.write(alignment)


def scores(path: Path) -> list[float]:
    return [float(line.rsplit("\t", 1)[1]) for line in path.read_text().splitlines()]


def check_public_configuration(work: Path) -> None:
    """Keep private configuration untracked and the public image lock limited to pins."""
    tracked = subprocess.check_output(["git", "ls-files", "-z"], cwd=ROOT, text=True).split("\0")
    for name in filter(None, tracked):
        path = Path(name)
        private = (
            bool({".claude", ".codex"} & set(path.parts))
            or path.name in {".env", ".envrc"}
            or path.name.startswith(".env.")
        )
        assert not private, f"Local configuration must not be tracked: {name}"

    pins: dict[str, str] = {}
    for line in (ROOT / "conf/genesis-images.lock").read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        key, separator, value = line.partition("=")
        assert separator and re.fullmatch(r"[A-Z][A-Z0-9_]*_IMAGE", key), (
            "Image lock accepts only IMAGE keys"
        )
        assert key not in pins, "Duplicate image pin"
        assert re.fullmatch(r"[A-Za-z0-9._:/-]+@sha256:[a-f0-9]{64}", value), (
            "Image lock requires a registry image with an immutable SHA256 digest"
        )
        pins[key] = value
    sources = [
        *sorted((ROOT / "workflows").glob("*.nf")),
        *sorted((ROOT / "modules").rglob("*.nf")),
    ]
    required = {
        key for path in sources for key in re.findall(r"dotenv\('([^']+)'\)", path.read_text())
    }
    assert pins.keys() == required, "Image pins must cover exactly the workflows' images"

    work.mkdir()
    subprocess.run(["git", "init", "--quiet", str(work)], check=True)
    shutil.copyfile(ROOT / ".gitignore", work / ".gitignore")
    private_paths = [
        prefix + name
        for prefix in ("", "nested/")
        for name in (
            ".env",
            ".env.local",
            ".env.production",
            ".envrc",
            ".claude/settings.json",
            ".claude/settings.local.json",
            ".codex/config.toml",
            ".codex/rules/default.rules",
        )
    ]
    public_paths = ["conf/genesis-images.lock", "AGENTS.md"]
    ignored = subprocess.check_output(
        ["git", "check-ignore", "--no-index", "--stdin"],
        cwd=work,
        input="\n".join([*private_paths, *public_paths]) + "\n",
        text=True,
    ).splitlines()
    assert ignored == private_paths, "Ignore rules must cover local config, not public image locks"


def check_quantification(work: Path) -> None:
    """Verify BAM counting retains MAPQ 0, NH fractions, overlap rules, and original rows."""
    sizes = work / "genome.chrom.sizes"
    sizes.write_text("chr1\t1000\n")
    bam = work / "reads.bam"
    write_bam(
        bam,
        "unique\t0\tchr1\t1\t0\t50M\t*\t0\t0\t*\t*\tNH:i:1\n"
        "multi\t0\tchr1\t21\t0\t50M\t*\t0\t0\t*\t*\tNH:i:2\n"
        "multi\t256\tchr1\t221\t0\t50M\t*\t0\t0\t*\t*\tNH:i:2\n",
    )
    bed = work / "peaks.bed"
    bed.write_text("chr1\t0\t100\tfirst\nchr1\t0\t100\tduplicate\nchr1\t200\t300\tsecond\n")
    coverage = work / "coverage.bedGraph"
    coverage.write_text("chr1\t0\t10\t2\nchr1\t10\t40\t4\nchr1\t40\t1000\t0\n")
    rpm, rpkm = work / "rpm.tsv", work / "rpkm.tsv"
    quantify(
        bed=bed,
        bam=bam,
        coverage=coverage,
        chrom_sizes=sizes,
        rpm_output=rpm,
        rpkm_output=rpkm,
        weighting="NH",
    )
    # NH: the unique read and half the multimapper fall in the first peak; the denominator
    # counts each read end once (2), so 1.5 / 2 and 0.5 / 2 per million.
    assert scores(rpm) == [750000, 750000, 250000]
    assert scores(rpkm) == [1.4, 1.4, 0]
    assert [line.rsplit("\t", 1)[0] for line in rpm.read_text().splitlines()] == (
        bed.read_text().splitlines()
    )
    # Primary: the secondary alignment is dropped; MAPQ 0 primaries count fully.
    quantify(
        bed=bed,
        bam=bam,
        coverage=coverage,
        chrom_sizes=sizes,
        rpm_output=rpm,
        rpkm_output=rpkm,
        weighting="primary",
    )
    assert scores(rpm) == [1000000, 1000000, 0]
    # Overlap edge cases, matching bedtools intersect on half-open intervals.
    edges = work / "edges.bam"
    write_bam(
        edges,
        # spans [95, 105): one base in each adjacent peak
        "both\t0\tchr1\t96\t0\t10M\t*\t0\t0\t*\t*\n"
        # spans [190, 200): ends exactly where the next peak starts
        "touching\t0\tchr1\t191\t0\t10M\t*\t0\t0\t*\t*\n"
        # deletion widens the reference span to [500, 520), inside the long peak
        "deletion\t0\tchr1\t501\t0\t5M10D5M\t*\t0\t0\t*\t*\n"
        "supplementary\t2048\tchr1\t1\t0\t10M\t*\t0\t0\t*\t*\n"
        "unmapped\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n"
        "elsewhere\t0\tchrUn\t1\t0\t10M\t*\t0\t0\t*\t*\n",
    )
    regions = work / "edges.bed"
    regions.write_text("chr1\t0\t100\nchr1\t100\t200\nchr1\t200\t300\nchr1\t300\t900\n")
    output = work / "edges.tsv"
    quantify_peaks(bed=regions, bam=edges, chrom_sizes=sizes, output=output, weighting="primary")
    # Three retained read ends; chrUn is outside chrom sizes and is not counted anywhere.
    # Scores are written with 12 significant digits.
    expected = [float(f"{value:.12g}") for value in (1e6 / 3, 2e6 / 3, 0, 1e6 / 3)]
    assert scores(output) == expected, scores(output)


def check_builds(work: Path, bash: str) -> None:
    """Verify the reference scripts build both architectures and honor --no-push."""
    repo = work / "repository with spaces"
    scripts = repo / "scripts"
    scripts.mkdir(parents=True)
    for name in (
        "find-projects.sh",
        "build-dockers.sh",
        "build-project-docker.sh",
        "build-yaml-docker.sh",
    ):
        destination = scripts / name
        shutil.copyfile(ROOT / "scripts" / name, destination)
        destination.chmod(0o755)
    (repo / "dockers").mkdir()
    (repo / "dockers" / ".dockerignore").write_text(".venv/\n")
    (repo / "environments").mkdir()
    (repo / "environments" / "GENESIS_ALIGNMENT.yaml").write_text("name: genesis_alignment\n")
    (repo / "environments" / "GENESIS_TOOLS.yaml").write_text("name: genesis_tools\n")
    project = repo / "genesis_tools"
    project.mkdir()
    for name in ("pyproject.toml", "uv.lock"):
        (project / name).write_text("mock\n")
    (repo / "conf").mkdir()
    image_lock = repo / "conf/genesis-images.lock"
    unrelated = "OTHER_IMAGE=example/other@sha256:" + "b" * 64
    image_lock.write_text(unrelated + "\nGENESIS_ALIGNMENT_IMAGE=old\n")
    private_env = repo / ".env"
    private_env.write_text("LOCAL_SETTING=private-test-fixture\n")
    binaries = work / "bin"
    binaries.mkdir()
    for name in ("docker", "pixi", "uv", "git"):
        executable = binaries / name
        executable.write_text(f"#!{sys.executable}\n" + FAKE_COMMAND)
        executable.chmod(0o755)
    log = work / "commands.jsonl"
    environment = {
        **os.environ,
        "PATH": f"{binaries}:{os.environ['PATH']}",
        "BUILD_LOG": str(log),
    }

    def run(arguments: list[str], *, fail: bool = False) -> subprocess.CompletedProcess[str]:
        env = {**environment, **({"FAIL_BUILD": "1"} if fail else {})}
        return subprocess.run(
            [bash, str(scripts / "build-dockers.sh"), *arguments],
            env=env,
            text=True,
            capture_output=True,
            check=False,
            timeout=20,
        )

    result = run(["--tag", "test", "--no-push"])
    assert result.returncode == 0, result.stderr
    contents = image_lock.read_text()
    assert unrelated in contents
    assert "LOCAL_SETTING" not in contents
    assert private_env.read_text() == "LOCAL_SETTING=private-test-fixture\n"
    assert "GENESIS_ALIGNMENT_IMAGE=kundajelab/genesis_alignment@sha256:" in contents
    assert "GENESIS_TOOLS_IMAGE=kundajelab/genesis_tools@sha256:" in contents
    commands = [json.loads(line) for line in log.read_text().splitlines()]
    assert not any(command[0] == "docker" and command[1][0] == "push" for command in commands)
    builds = [args for name, args in commands if name == "docker" and args[0] == "build"]
    assert len(builds) == 2
    for args in builds:
        assert args[args.index("--platform") + 1] == "linux/amd64,linux/arm64"
    assert (repo / "environments" / "GENESIS_ALIGNMENT.lock").read_text() == "mock lock\n"
    assert not (project / ".dockerignore").exists()
    log.write_text("")
    result = run(["--tag", "published", "--push", "--", "genesis_tools"])
    assert result.returncode == 0, result.stderr
    commands = [json.loads(line) for line in log.read_text().splitlines()]
    assert sum(name == "docker" and args[0] == "push" for name, args in commands) == 1
    assert sum(name == "docker" and args[0] == "build" for name, args in commands) == 1
    assert "GENESIS_ALIGNMENT_IMAGE=kundajelab/genesis_alignment@sha256:" in image_lock.read_text()
    log.write_text("")
    result = run(["--no-push", "--", "genesis_tools"])
    assert result.returncode == 0, result.stderr
    commands = [json.loads(line) for line in log.read_text().splitlines()]
    build = next(args for name, args in commands if name == "docker" and args[0] == "build")
    assert "genesis_tools:0.1.0" in build
    assert not any(name == "docker" and args[0] == "push" for name, args in commands)
    unchanged = image_lock.read_text()
    result = run(["--tag", "broken", "--no-push"], fail=True)
    assert result.returncode != 0
    assert image_lock.read_text() == unchanged
    assert private_env.read_text() == "LOCAL_SETTING=private-test-fixture\n"


def config_params(path: Path) -> set[str]:
    """Names assigned in a config file's params block."""
    block = re.search(r"^params \{\n(.*?)^\}", path.read_text(), re.M | re.S)
    return set(re.findall(r"^    (\w+) = ", block.group(1), re.M)) if block else set()


def check_workflow_configs() -> None:
    """Each workflow has its own config and parameter schema; the shared sections stay identical."""
    workflows = sorted(path.stem for path in (ROOT / "workflows").glob("*.nf"))
    assert workflows == ["atac-seq", "dap-seq", "rna-seq", "scatac-reads", "scatac-seq"], workflows
    shared = config_params(ROOT / "nextflow.config") | {"run_dir"}
    schemas = {}
    for name in workflows:
        config = ROOT / "conf" / f"{name}.config"
        text = config.read_text()
        assert "includeConfig '../nextflow.config'" in text, config
        schema_path = f"assets/{name}.params.schema.json"
        assert f'parametersSchema = "${{projectDir.parent}}/{schema_path}"' in text, config
        schema = json.loads((ROOT / schema_path).read_text())
        declared = {key for group in schema["$defs"].values() for key in group["properties"]}
        expected = shared | config_params(config)
        assert declared == expected, (name, declared ^ expected)
        assert schema["allOf"] == [{"$ref": f"#/$defs/{key}"} for key in schema["$defs"]], name
        schemas[name] = schema["$defs"]
    # nf-schema cannot $ref another file, so schemas repeat shared sections verbatim.
    common = set.intersection(*(set(defs) for defs in schemas.values()))
    sections = {"input_output_options", "resource_options", "internal_options"}
    assert sections <= common, sections - common
    for key in set().union(*schemas.values()):
        copies = {json.dumps(defs[key], sort_keys=True) for defs in schemas.values() if key in defs}
        assert len(copies) == 1, key


def main() -> None:
    """Run focused checks without external services, Docker, or genomic data."""
    assert not (ROOT / "lib" / "WorkflowSupport.groovy").exists()
    assert not (ROOT / "custom_tools").exists()
    modules = list((ROOT / "modules").rglob("*.nf"))
    assert modules
    for path in modules:
        source = path.read_text()
        assert "WorkflowSupport" not in source, path
        assert "with-environment.sh" not in source and "uv run" not in source, path
        if re.search(r"^process \w+", source, re.M) is None:
            # Subworkflows only compose process modules, which declare tools and resources.
            assert re.search(r"^workflow \w+", source, re.M), path
            assert "    conda " not in source and "container " not in source, path
            continue
        assert "include { dotenv } from 'plugin/nf-dotenv'" in source, path
        assert "    cpus " in source and "    memory " in source and "    time " in source, path
        assert source.count("    conda ") == 1, path
    check_workflow_configs()
    interpreters = list(dict.fromkeys(["/bin/bash", shutil.which("bash") or "/bin/bash"]))
    with tempfile.TemporaryDirectory(prefix="genesis-style-") as temporary:
        work = Path(temporary)
        check_public_configuration(work / "configuration")
        quant = work / "quantification"
        quant.mkdir()
        check_quantification(quant)
        for index, bash in enumerate(interpreters):
            folder = work / f"builds-{index}"
            folder.mkdir()
            check_builds(folder, bash)
    print("Public configuration, Nextflow style, BAM quantification and mocked builds passed.")


if __name__ == "__main__":
    main()
