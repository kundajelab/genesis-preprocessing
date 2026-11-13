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
    environment = {
        **os.environ,
        "PATH": f"{tools}:{BASE_PATH}",
        "NERSC_HOST": nersc_host,
    }
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
        "nersc",
        host_environment(work / "nersc", {"sacctmgr": account("gen0035")}, CFS=cfs),
    )
    assert perlmutter.stdout.strip() == f"{cfs}/gen0035/{user}/references", perlmutter
    # Any profile list naming nersc; a GPU allocation (PROJECT_g) uses the project's storage.
    gpu = default(
        "nersc,test",
        host_environment(work / "gpu", {"sacctmgr": account("m1234_g")}, CFS=cfs),
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

    # The persist folder, runs/ beside these references/: none for local profiles.
    script = str(ROOT / "scripts/get-default-persist-dir.sh")
    perlmutter = default(
        "nersc",
        host_environment(work / "nersc", {"sacctmgr": account("gen0035")}, CFS=cfs),
    )
    assert perlmutter.stdout.strip() == f"{cfs}/gen0035/{user}/runs", perlmutter
    gpu = default(
        "nersc,test",
        host_environment(work / "gpu", {"sacctmgr": account("m1234_g")}, CFS=cfs),
    )
    assert gpu.stdout.strip() == f"{cfs}/m1234/{user}/runs", gpu
    sherlock = default("sherlock", host_environment(work / "sherlock", {}, OAK=oak))
    assert sherlock.stdout.strip() == f"{oak}/{user}/genesis/runs", sherlock
    no_oak = default("sherlock", host_environment(work / "no-oak", {}))
    assert no_oak.returncode == 2 and "OAK is not set" in no_oak.stderr, no_oak
    local = default("local,docker", host_environment(work / "local", {}))
    assert local.returncode == 0 and local.stdout == "", local


def check_pipeline_references(work: Path) -> None:
    """-R/--references and -P/--persist-dir, or the host's defaults, are created and passed on."""
    log = work / "pixi-arguments.txt"
    pixi = f"#!/bin/sh\nprintf '%s\\n' \"$@\" > '{log}'\n"
    sheet = work / "samples.tsv"
    work.mkdir(parents=True, exist_ok=True)
    sheet.write_text("sample_id\n")
    user = os.environ.get("USER", "user")

    warnings: list[str] = []

    def launch(environment: dict[str, str], *arguments: str) -> list[str]:
        result = subprocess.run(
            [
                str(ROOT / "scripts/run-pipeline.sh"),
                "-w",
                str(work / "workspace"),
                *arguments,
            ],
            env=environment,
            cwd=work,
            check=True,
            capture_output=True,
            text=True,
        )
        warnings[:] = [result.stderr]
        return log.read_text().splitlines()

    def values(arguments: list[str], option: str) -> list[str]:
        return [arguments[i + 1] for i, a in enumerate(arguments) if a == option]

    def references(arguments: list[str]) -> list[str]:
        return values(arguments, "--references")

    def persist(arguments: list[str]) -> list[str]:
        return values(arguments, "--persist_dir")

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

    # The persist folder: the host's default, created, or -P/--persist-dir, which may be relative.
    kept = work / "cfs" / "gen0035" / user / "runs"
    assert persist(launch(perlmutter, str(sheet))) == [str(kept)] and kept.is_dir()
    (relative,) = persist(launch(perlmutter, "-P", "keep", str(sheet)))
    assert Path(relative).resolve() == (work / "keep").resolve(), relative
    assert (work / "keep").is_dir()
    # Without a default (Sherlock without $OAK), outputs stay in the workspace, with a warning.
    assert persist(launch(no_default, "-R", str(given), str(sheet))) == []
    assert "pass -P/--persist-dir" in warnings[0], warnings
    given_persist = work / "given" / "persist"
    persisted = launch(
        no_default, "-R", str(given), "--persist-dir", str(given_persist), str(sheet)
    )
    assert persist(persisted) == [str(given_persist)], persisted
    # Locally there is no default and no warning: outputs stay in the workspace.
    local = host_environment(work / "local", {"pixi": pixi, "docker": TRUE})
    assert persist(launch(local, "-R", str(given), str(sheet))) == [] and not warnings[0]


def check_submit_pipeline(work: Path) -> None:
    """submit-pipeline.sh queues run-pipeline.sh, its arguments unchanged, with site settings."""
    work.mkdir(parents=True)
    sbatch_arguments, script, pulls = (
        work / "sbatch.txt",
        work / "script.sh",
        work / "pulls.txt",
    )
    recorded = work / "pixi.txt"
    sbatch = f"#!/bin/sh\nprintf '%s\\n' \"$@\" > '{sbatch_arguments}'\ncat > '{script}'\necho 42\n"
    podman = f"#!/bin/sh\necho \"$@\" >> '{pulls}'\n"
    # pixi records what it would run and the environment it would run it in.
    pixi = (
        "#!/bin/sh\n{ printf '%s\\n' \"$@\"; env | "
        f"grep -E '^(SLURM_|NXF_|GENESIS_)' | sort; }} > '{recorded}'\n"
    )
    sheet = work / "my sheet.tsv"
    sheet.write_text("sample_id\n")
    images = sum("_IMAGE=" in line for line in (ROOT / "conf/genesis-images.lock").open())

    def submit(environment: dict[str, str], *arguments: str) -> subprocess.CompletedProcess[str]:
        command = [
            str(ROOT / "scripts/submit-pipeline.sh"),
            "-L",
            str(work / "logs"),
            *arguments,
        ]
        # Run folders go in this test's $SCRATCH/workspace, never the checkout's.
        environment = {**environment, "SCRATCH": str(work / "scratch")}
        return subprocess.run(command, env=environment, cwd=work, capture_output=True, text=True)

    def options() -> list[str]:
        return sbatch_arguments.read_text().splitlines()

    tools = {
        "sbatch": sbatch,
        "podman-hpc": podman,
        "pixi": pixi,
        "sacctmgr": account("m1"),
    }
    nersc = host_environment(work / "nersc", tools, "perlmutter", CFS=str(work / "cfs"))
    run = ["-p", "nersc,test", "-R", str(work / "refs"), str(sheet), "-resume"]
    result = submit(nersc, "-T", "01:00:00", "-M", "4G", *run)
    assert result.returncode == 0, result
    assert "job 42" in result.stdout, result.stdout
    for option in ["--qos=shared", "--constraint=cpu", "--time=01:00:00", "--mem=4G"]:
        assert option in options(), (option, options())
    for option in [
        "--no-requeue",
        "--signal=B:TERM@120",
        "--job-name=genesis-my sheet",
    ]:
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
    slurm = {
        "SLURM_JOB_ID": "42",
        "SLURM_MEM_PER_NODE": "4096",
        "SLURM_CONF": "/etc/slurm.conf",
    }
    job_environment = {**nersc, **slurm, "SCRATCH": str(work / "scratch")}
    job = subprocess.run(
        ["bash", str(script)],
        env=job_environment,
        cwd=work,
        capture_output=True,
        text=True,
    )
    assert job.returncode == 0, job
    lines = recorded.read_text().splitlines()
    assert lines[lines.index("--input") + 1] == str(sheet), lines
    assert "-resume" in lines and "nersc,test" in lines, lines
    # The default persist folder, from the account sacctmgr reports, reaches Nextflow too.
    user = os.environ.get("USER", "user")
    assert lines[lines.index("--persist_dir") + 1] == str(work / "cfs" / "m1" / user / "runs"), (
        lines
    )
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

    # A finished run is summarized from the folders Nextflow reported in the job's log; a failed
    # one is not.
    publish, task = work / "runs/my sheet", work / "scratch/my sheet/work/ab/cdef1234"
    reporting_pixi = f"""#!/bin/sh
mkdir -p '{publish}/trace' '{task}'
trace='{publish}/trace/execution_trace_1.txt'
printf 'hash\\tnative_id\\tstatus\\trealtime\\t%%cpu\\tpeak_rss\\twchar\\n' > "$trace"
printf 'ab/cdef12\\t1\\tCOMPLETED\\t1h\\t100%%\\t1 GB\\t1 GB\\n' >> "$trace"
printf '#SBATCH -c 1\\n#SBATCH --mem 2048M\\n' > '{task}/.command.run'
echo '  publish_base   : {publish}'
echo '  workDir        : {task.parent.parent}'
exit "${{EXIT:-0}}"
"""
    reporting = host_environment(
        work / "reporting", {"sbatch": sbatch, "pixi": reporting_pixi, "apptainer": TRUE}
    )
    assert submit(reporting, "-R", str(work / "refs"), str(sheet)).returncode == 0
    log = work / "logs/genesis-my sheet-7.out"
    for status in (3, 0):
        with log.open("w") as stream:
            job = subprocess.run(
                ["bash", str(script)],
                env={**reporting, "SLURM_JOB_ID": "7", "EXIT": str(status)},
                cwd=work, stdout=stream, stderr=subprocess.PIPE, text=True,
            )  # fmt: skip
        assert job.returncode == status, job.stderr
        summary = publish / "resources.tsv"
        assert summary.exists() == (status == 0), (status, log.read_text())
    summary_lines = (publish / "resources.tsv").read_text().splitlines()
    assert "cpu-hours\t1.000" in summary_lines and "node-hours\tnan" in summary_lines

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
        stderr=subprocess.DEVNULL,
    )
    time.sleep(3)
    job_process.send_signal(signal.SIGTERM)
    status = job_process.wait(timeout=30)
    assert status == 143 and marker.exists(), (status, marker.exists())


def check_shim_descriptors(work: Path, shim: Path) -> None:
    """A daemon podman-hpc leaves running cannot hold the job script's output pipes open.

    Nextflow's job script starts the container under `3>&1 1>&2 2>&3 | tee .command.err`,
    leaving fd 3 open on that pipe. podman-hpc's squashfuse daemon detaches with fds 0-2 on
    /dev/null but kept fd 3, so the tee waited, and the finished job with it, until its time
    limit. The shim closes such descriptors.
    """
    work.mkdir(parents=True)
    # A stand-in podman-hpc that leaves a detached daemon, like squashfuse_ll, then exits.
    stand_in(
        work / "bin",
        "podman-hpc",
        "#!/bin/sh\nnohup sleep 30 < /dev/null > /dev/null 2>&1 &\n"
        "echo container output\necho container error >&2\n",
    )
    out, err = work / ".command.out", work / ".command.err"
    launch = (
        f"(set -o pipefail; (podman run image | tee '{out}') 3>&1 1>&2 2>&3 | tee '{err}')"
        " 3>&1 1>&2 2>&3 & pid=$!; wait $pid"
    )
    start = time.monotonic()
    subprocess.run(
        ["bash", "-c", launch],
        env={**os.environ, "PATH": f"{shim}:{work / 'bin'}:{BASE_PATH}"},
        capture_output=True,
        timeout=60,
        check=True,
    )
    elapsed = time.monotonic() - start
    assert elapsed < 15, f"the job script waited {elapsed:.0f} s for the daemon"
    assert out.read_text() == "container output\n" and err.read_text() == "container error\n"


def check_resource_summary(work: Path) -> None:
    """scripts/summarize-resources.sh totals each task attempt that ran, from its trace and work.

    Requests come from Slurm's directives or else the container's limits; cached rows and attempts
    with no realtime are skipped; node-hours add the manager job and are nan without Slurm.
    """
    persist, scratch, home = work / "runs/demo", work / "workspace/demo/work", work / "home"
    (persist / "trace").mkdir(parents=True)
    header = "task_id\thash\tnative_id\tname\tstatus\texit\tsubmit\tduration\trealtime\t%cpu"
    header += "\tpeak_rss\tpeak_vmem\trchar\twchar\n"
    (persist / "trace/execution_trace_2026-11-01_00-00-00.txt").write_text(
        header
        + "1\tab/cdef12\t101\tA (s1)\tCOMPLETED\t0\t-\t1h 5m\t1h\t150.0%\t1 GB\t2 GB\t0\t2 GB\n"
        + "2\t12/345678\t102\tB (s1)\tFAILED\t137\t-\t31m\t30m\t-\t-\t-\t0\t1 MB\n"
    )
    (persist / "trace/execution_trace_2026-11-02_00-00-00.txt").write_text(
        header
        + "1\tab/cdef12\t101\tA (s1)\tCACHED\t0\t-\t1h 5m\t1h\t150.0%\t1 GB\t2 GB\t0\t2 GB\n"
        + "3\tff/000000\t103\tC (s1)\tABORTED\t-\t-\t-\t-\t-\t-\t-\t-\t-\n"
    )
    for task, script in {
        "ab/cdef1234": "#SBATCH -c 2\n#SBATCH --mem 4096M\n#SBATCH --qos=shared --mem=8192M\n",
        "12/3456789a": "podman run -i --cpu-shares 4096 --memory 2048m image bash\n",
    }.items():
        (scratch / task).mkdir(parents=True)
        (scratch / task / ".command.run").write_text(script)
    (home / "logs/genesis").mkdir(parents=True)
    (home / "logs/genesis/genesis-demo-900.out").write_text(
        f"Nextflow manager: job 900 on nid000001\n  publish_base   : {persist}\n"
    )
    # Each job, the two tasks and the manager, ran an hour billed 64 of a node's 256 threads.
    sacct = '#!/bin/sh\nwhile [ "$1" != -j ]; do shift; done\n'
    sacct += "for job in $(echo \"$2\" | tr , ' '); do echo '3600|billing=64,cpu=64'; done\n"
    script = [str(ROOT / "scripts/summarize-resources.sh"), str(persist), str(scratch)]

    def summarize(tools: dict[str, str]) -> dict[str, str]:
        environment = host_environment(work / f"tools-{len(tools)}", tools, HOME=str(home))
        run = subprocess.run(script, env=environment, capture_output=True, text=True)
        assert run.returncode == 0, run.stderr
        lines = (persist / "resources.tsv").read_text().splitlines()
        pairs = [line.split("\t") for line in lines]
        assert all(len(pair) == 2 for pair in pairs), lines
        assert [name for name, _ in pairs] == [
            "cpu-hours", "mem-gb-hours", "node-hours", "scratch-gb", "persistent-disk-gb",
            "gb-written", "avg-%cpu", "avg-%mem",
        ], lines  # fmt: skip
        return dict(pairs)

    values = summarize({"sacct": sacct})
    # A: 2 CPUs and 4 GB for an hour at 150% and 1 GB; B: 4 CPUs and 2 GB for half an hour, with
    # no %cpu or peak_rss, so it counts in the totals but not the averages. They wrote 2 GB + 1 MB.
    disk = {name: values.pop(name) for name in ("scratch-gb", "persistent-disk-gb")}
    assert values == {
        "cpu-hours": "4.000", "mem-gb-hours": "5.000", "node-hours": "0.750",
        "gb-written": "2.001", "avg-%cpu": "75.000", "avg-%mem": "25.000",
    }, values  # fmt: skip
    assert all(float(value) >= 0 for value in disk.values()), disk
    assert summarize({})["node-hours"] == "nan"
    # A task that ran needs its work folder for its requests.
    shutil.rmtree(scratch / "12")
    missing = subprocess.run(
        script, env=host_environment(work / "tools-missing", {}, HOME=str(home)),
        capture_output=True, text=True,
    )  # fmt: skip
    assert missing.returncode != 0 and "No work folder for 12/345678" in missing.stderr, missing


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="genesis-profile-") as temporary:
        folder = Path(temporary)
        check_default_profiles(folder / "defaults")
        check_resource_summary(folder / "resources")
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
                # The strict parser flattens maps (resourceLimits.cpus = ...); the legacy one,
                # which a shell may select, prints them whole.
                env={**os.environ, "NXF_SYNTAX_PARSER": "v2"},
                text=True,
                timeout=30,
            )
            assert f"process.executor = '{executor}'" in output
            assert f"{runtime}.enabled = true" in output
            if profile == "nersc":
                assert "executor.queueSize = 50" in output
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
        check_shim_descriptors(folder / "descriptors", shim)
    print(
        "PASS: default profile and references per host, pipeline submission, seven parsed"
        " profiles, podman-hpc shim, resource summary"
    )


if __name__ == "__main__":
    main()
