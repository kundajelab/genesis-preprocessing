"""Parse every deployment profile, check the default profile per host, and the podman-hpc shim."""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


# Only system folders, so the host's own scheduler or container tools do not decide the result.
BASE_PATH = "/usr/bin:/bin"


def default_profile(folder: Path, tools: list[str], nersc_host: str = "") -> str:
    """The default profile on a host with these tools, as stand-ins, on a minimal PATH."""
    folder.mkdir(parents=True)
    for name in tools:
        script = folder / name
        script.write_text("#!/bin/sh\nexit 0\n")
        script.chmod(0o755)
    return subprocess.check_output(
        [str(ROOT / "scripts/get-default-profile.sh")],
        env={**os.environ, "PATH": f"{folder}:{BASE_PATH}", "NERSC_HOST": nersc_host},
        stderr=subprocess.DEVNULL,
        text=True,
    ).strip()


def check_default_profiles(work: Path) -> None:
    """Perlmutter is NERSC, as is Slurm with podman-hpc; Slurm with Apptainer is Sherlock.

    Anything else, including another Slurm cluster, is local Docker.
    """
    if any(shutil.which(tool, path=BASE_PATH) for tool in ("sbatch", "podman-hpc", "apptainer")):
        print("SKIP: default profiles; this host's system folders have a scheduler or runtime")
        return
    for name, tools, nersc_host, expected in [
        ("local", ["docker"], "", "local,docker"),
        ("perlmutter", [], "perlmutter", "nersc"),
        ("nersc", ["sbatch", "podman-hpc"], "", "nersc"),
        ("sherlock", ["sbatch", "apptainer"], "", "sherlock"),
        ("other-slurm", ["sbatch"], "", "local,docker"),
    ]:
        assert default_profile(work / name, tools, nersc_host) == expected, name


def stand_in(folder: Path, name: str, text: str) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    (folder / name).write_text(text)
    (folder / name).chmod(0o755)


def host_environment(
    tools: Path, names: dict[str, str], nersc_host: str = "", **variables: str
) -> dict[str, str]:
    """A simulated host: stand-in tools first on a minimal PATH, and its environment variables."""
    for name, text in names.items():
        stand_in(tools, name, text)
    tools.mkdir(parents=True, exist_ok=True)
    environment = {**os.environ, "PATH": f"{tools}:{BASE_PATH}", "NERSC_HOST": nersc_host}
    for name in ("CFS", "OAK"):
        environment.pop(name, None)
    return {**environment, **variables}


def account(name: str) -> str:
    return f"#!/bin/sh\nprintf '%s\\n' '{name}'\n"


TRUE = "#!/bin/sh\nexit 0\n"


def check_default_references(work: Path) -> None:
    """Project storage for profiles naming nersc or sherlock, the checkout's references/ else."""
    script = str(ROOT / "scripts/get-default-references.sh")
    user = os.environ.get("USER", "user")
    cfs, oak = str(work / "cfs"), str(work / "oak")

    def default(profile: str, environment: dict[str, str]) -> subprocess.CompletedProcess[str]:
        return subprocess.run([script, profile], env=environment, capture_output=True, text=True)

    perlmutter = default(
        "nersc", host_environment(work / "nersc", {"sacctmgr": account("gen0035")}, CFS=cfs)
    )
    assert perlmutter.stdout.strip() == f"{cfs}/gen0035/{user}/references", perlmutter
    # Any profile list naming nersc; a GPU allocation (PROJECT_g) uses the project's storage.
    gpu = default(
        "nersc,test", host_environment(work / "gpu", {"sacctmgr": account("m1234_g")}, CFS=cfs)
    )
    assert gpu.stdout.strip() == f"{cfs}/m1234/{user}/references", gpu
    missing = default("nersc", host_environment(work / "none", {"sacctmgr": account("")}))
    assert missing.returncode == 2 and "pass -R/--references" in missing.stderr, missing
    sherlock = default("sherlock,apptainer", host_environment(work / "sherlock", {}, OAK=oak))
    assert sherlock.stdout.strip() == f"{oak}/{user}/genesis/references", sherlock
    no_oak = default("sherlock", host_environment(work / "no-oak", {}))
    assert no_oak.returncode == 2 and "pass -R/--references" in no_oak.stderr, no_oak
    local = default("local,docker", host_environment(work / "local", {}))
    assert local.stdout.strip() == str(ROOT / "references"), local


def check_pipeline_references(work: Path) -> None:
    """-R/--references, or the host's default, is created and passed to Nextflow once."""
    log = work / "pixi-arguments.txt"
    pixi = f"#!/bin/sh\nprintf '%s\\n' \"$@\" > '{log}'\n"
    sheet = work / "samples.tsv"
    work.mkdir(parents=True, exist_ok=True)
    sheet.write_text("sample_id\n")
    user = os.environ.get("USER", "user")

    def launch(environment: dict[str, str], *arguments: str) -> list[str]:
        subprocess.run(
            [str(ROOT / "scripts/run-pipeline.sh"), "-w", str(work / "workspace"), *arguments],
            env=environment,
            cwd=work,
            check=True,
            capture_output=True,
            text=True,
        )
        return log.read_text().splitlines()

    def references(arguments: list[str]) -> list[str]:
        return [arguments[i + 1] for i, a in enumerate(arguments) if a == "--references"]

    perlmutter = host_environment(
        work / "nersc",
        {"pixi": pixi, "podman-hpc": TRUE, "sacctmgr": account("gen0035")},
        "perlmutter",
        CFS=str(work / "cfs"),
    )
    default = work / "cfs" / "gen0035" / user / "references"
    assert references(launch(perlmutter, str(sheet))) == [str(default)]
    assert default.is_dir()
    # The default follows -p: here a profile list naming nersc.
    assert references(launch(perlmutter, "-p", "nersc,test", str(sheet))) == [str(default)]
    # An explicit folder wins, even on a host with no default (Sherlock without $OAK).
    no_default = host_environment(
        work / "no-oak", {"pixi": pixi, "sbatch": TRUE, "apptainer": TRUE}
    )
    given = work / "given" / "references"
    assert references(launch(no_default, "-R", str(given), str(sheet))) == [str(given)]
    assert given.is_dir()
    # A relative folder is relative to where the command runs, like the sheet and workspace
    # (compared resolved: macOS temporary folders are reached through a symlink).
    (relative,) = references(launch(no_default, "--references", "relative", str(sheet)))
    assert Path(relative).resolve() == (work / "relative").resolve(), relative
    # Nextflow arguments pass through untouched, even one named --references.
    arguments = launch(no_default, "-R", str(given), str(sheet), "-resume", "--references", "x")
    assert arguments[-3:] == ["-resume", "--references", "x"], arguments


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="genesis-profile-") as temporary:
        folder = Path(temporary)
        check_default_profiles(folder / "defaults")
        if not any(shutil.which(tool, path=BASE_PATH) for tool in ("sbatch", "apptainer")):
            check_default_references(folder / "references")
            check_pipeline_references(folder / "pipeline")
        for profile, executor, runtime in [
            ("local,docker", "local", "docker"),
            ("sherlock", "slurm", "apptainer"),
            ("sherlock,apptainer", "slurm", "apptainer"),
            ("nersc", "slurm", "podman"),
            ("slurm,apptainer", "slurm", "apptainer"),
            ("pbspro,apptainer", "pbspro", "apptainer"),
            ("lsf,apptainer", "lsf", "apptainer"),
            ("sge,apptainer", "sge", "apptainer"),
        ]:
            output = subprocess.check_output(
                [
                    "nextflow",
                    "-C",
                    str(ROOT / "nextflow.config"),
                    "config",
                    str(ROOT / "workflows"),
                    "-profile",
                    profile,
                    "-flat",
                ],
                text=True,
                timeout=30,
            )
            assert f"process.executor = '{executor}'" in output
            assert f"{runtime}.enabled = true" in output
            if profile == "nersc":
                assert "executor.queueSize = 15" in output
                assert "executor.pollInterval = '5 min'" in output
                assert "--constraint=cpu --qos=${params.queue ?: 'shared'}" in output
                assert "toMega() + 4096}M" in output
                assert "podman.runOptions = '--userns=keep-id --user $(id -u):$(id -g)'" in output
                assert "podman.registry = 'docker.io'" in output
                assert "podman.temp" not in output
                assert "process.resourceLimits.cpus = 128" in output
                assert "process.resourceLimits.memory = '247 GB'" in output
                assert "process.resourceLimits.time = '2d'" in output
        # On NERSC, `podman` is the shim, which hands every argument to podman-hpc unchanged,
        # with PATH minus its own folder: podman-hpc runs podman by name itself.
        fake = folder / "podman-hpc"
        fake.write_text('#!/bin/sh\nprintf "%s\\n" "$PATH" "$@"\n')
        fake.chmod(0o755)
        shim = ROOT / "scripts/podman-hpc-bin"
        output = subprocess.check_output(
            [str(shim / "podman"), "info", "two words"],
            env={**os.environ, "PATH": f"{shim}:{folder}:{os.environ['PATH']}"},
            text=True,
        )
        path, *arguments = output.splitlines()
        assert arguments == ["info", "two words"], arguments
        assert str(shim) not in path.split(":") and str(folder) in path.split(":"), path
    print("PASS: default profile and references per host, seven parsed profiles, podman-hpc shim")


if __name__ == "__main__":
    main()
