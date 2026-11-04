"""Parse every deployment profile, check the default profile per host, and the podman-hpc shim."""

from __future__ import annotations

import os
import shutil
import signal
import subprocess
import tempfile
import time
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


def check_submit_pipeline(work: Path) -> None:
    """submit-pipeline.sh queues run-pipeline.sh, its arguments unchanged, with site settings."""
    work.mkdir(parents=True)
    sbatch_arguments, script, pulls = work / "sbatch.txt", work / "script.sh", work / "pulls.txt"
    recorded = work / "pixi.txt"
    sbatch = f"#!/bin/sh\nprintf '%s\\n' \"$@\" > '{sbatch_arguments}'\ncat > '{script}'\necho 42\n"
    podman = f"#!/bin/sh\necho \"$@\" >> '{pulls}'\n"
    # pixi records what it would run and the environment it would run it in.
    pixi = (
        f"#!/bin/sh\n{{ printf '%s\\n' \"$@\"; env | grep -E '^(SLURM_|NXF_|GENESIS_)' | sort; }}"
        f" > '{recorded}'\n"
    )
    sheet = work / "my sheet.tsv"
    sheet.write_text("sample_id\n")
    images = sum("_IMAGE=" in line for line in (ROOT / "conf/genesis-images.lock").open())

    def submit(environment: dict[str, str], *arguments: str) -> subprocess.CompletedProcess[str]:
        command = [str(ROOT / "scripts/submit-pipeline.sh"), "-L", str(work / "logs"), *arguments]
        # Run folders go in this test's $SCRATCH/workspace, never the checkout's.
        environment = {**environment, "SCRATCH": str(work / "scratch")}
        return subprocess.run(command, env=environment, cwd=work, capture_output=True, text=True)

    def options() -> list[str]:
        return sbatch_arguments.read_text().splitlines()

    tools = {"sbatch": sbatch, "podman-hpc": podman, "pixi": pixi, "sacctmgr": account("m1")}
    nersc = host_environment(work / "nersc", tools, "perlmutter")
    run = ["-p", "nersc,test", "-R", str(work / "refs"), str(sheet), "-resume"]
    result = submit(nersc, "-T", "01:00:00", "-M", "4G", *run)
    assert result.returncode == 0, result
    assert "job 42" in result.stdout, result.stdout
    for option in ["--qos=shared", "--constraint=cpu", "--time=01:00:00", "--mem=4G"]:
        assert option in options(), (option, options())
    for option in ["--no-requeue", "--signal=B:TERM@120", "--job-name=genesis-my sheet"]:
        assert option in options(), (option, options())
    # The job runs where the command ran (compared resolved: macOS temporary folders are symlinks).
    (chdir,) = [option for option in options() if option.startswith("--chdir=")]
    assert Path(chdir.removeprefix("--chdir=")).resolve() == work.resolve(), chdir
    # Every pinned image is pulled before submitting, on the login node, and not again in the job.
    assert len(pulls.read_text().splitlines()) == images, pulls.read_text()
    text = script.read_text()
    assert "export GENESIS_IMAGES_PULLED=1" in text and "-Xmx3g" in text, text
    # The job: run-pipeline.sh gets the arguments unchanged, and Slurm settings of this job, except
    # SLURM_CONF, do not reach the tasks Nextflow submits.
    slurm = {"SLURM_JOB_ID": "42", "SLURM_MEM_PER_NODE": "4096", "SLURM_CONF": "/etc/slurm.conf"}
    job_environment = {**nersc, **slurm, "SCRATCH": str(work / "scratch")}
    job = subprocess.run(
        ["bash", str(script)], env=job_environment, cwd=work, capture_output=True, text=True
    )
    assert job.returncode == 0, job
    lines = recorded.read_text().splitlines()
    assert lines[lines.index("--input") + 1] == str(sheet), lines
    assert "-resume" in lines and "nersc,test" in lines, lines
    assert "SLURM_CONF=/etc/slurm.conf" in lines, lines
    assert not any(line.startswith(("SLURM_JOB_ID=", "SLURM_MEM_PER_NODE=")) for line in lines)
    # The heap comes first, so NXF_OPTS already set (e.g. by the checks) can override it.
    assert "NXF_ANSI_LOG=false" in lines, lines
    assert any(line.startswith("NXF_OPTS=-Xmx3g") for line in lines), lines
    assert len(pulls.read_text().splitlines()) == images, "the job pulled images again"
    # NERSC's workflow QOS runs on its own nodes.
    assert submit(nersc, "-Q", "workflow", *run).returncode == 0
    assert "--qos=workflow" in options() and "--constraint=cron" in options(), options()

    # Sherlock: a partition, normal unless given; no images to pull.
    sherlock_tools = {"sbatch": sbatch, "pixi": pixi, "apptainer": TRUE}
    sherlock = host_environment(work / "sherlock", sherlock_tools)
    pulls.unlink()
    assert submit(sherlock, str(sheet)).returncode == 0
    assert "--partition=normal" in options() and not pulls.exists(), options()
    assert submit(sherlock, "-Q", "owners", "--", str(sheet)).returncode == 0
    assert "--partition=owners" in options(), options()

    # Without a Slurm profile there is nothing to submit to.
    local = submit(host_environment(work / "local", {"docker": TRUE}), str(sheet))
    assert local.returncode == 2 and "needs a Slurm profile" in local.stderr, local

    # Slurm's SIGTERM before the walltime reaches Nextflow, and the job waits for its cleanup.
    marker = work / "cleaned-up"
    slow_pixi = (
        f"#!/bin/bash\ntrap 'sleep 1; touch \"{marker}\"; exit 143' TERM\nsleep 30 &\nwait\n"
    )
    stand_in(work / "sherlock", "pixi", slow_pixi)
    assert submit(sherlock, "-R", str(work / "refs"), str(sheet)).returncode == 0
    job_process = subprocess.Popen(
        ["bash", str(script)],
        env={**sherlock, "SCRATCH": str(work / "scratch")},
        cwd=work,
        stdout=subprocess.DEVNULL,
    )
    time.sleep(3)
    job_process.send_signal(signal.SIGTERM)
    status = job_process.wait(timeout=30)
    assert status == 143 and marker.exists(), (status, marker.exists())


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="genesis-profile-") as temporary:
        folder = Path(temporary)
        check_default_profiles(folder / "defaults")
        if not any(shutil.which(tool, path=BASE_PATH) for tool in ("sbatch", "apptainer")):
            check_default_references(folder / "references")
            check_pipeline_references(folder / "pipeline")
            check_submit_pipeline(folder / "submit")
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
    print(
        "PASS: default profile and references per host, pipeline submission, seven parsed "
        "profiles, podman-hpc shim"
    )


if __name__ == "__main__":
    main()
