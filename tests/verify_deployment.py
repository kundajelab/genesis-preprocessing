"""Parse every deployment profile and verify local defaults regardless of host tools."""

from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="genesis-profile-") as temporary:
        folder = Path(temporary)
        for name in ("sbatch", "qsub", "docker"):
            script = folder / name
            script.write_text("#!/bin/sh\nexit 0\n")
            script.chmod(0o755)
        default = subprocess.check_output(
            [str(ROOT / "scripts/get-default-profile.sh")],
            env={**os.environ, "PATH": f"{folder}:{os.environ['PATH']}", "NERSC_HOST": ""},
            text=True,
        ).strip()
        assert default == "local,docker"
        perlmutter = subprocess.check_output(
            [str(ROOT / "scripts/get-default-profile.sh")],
            env={**os.environ, "NERSC_HOST": "perlmutter"},
            text=True,
        ).strip()
        assert perlmutter == "nersc"
        for profile, executor, runtime in [
            ("local,docker", "local", "docker"),
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
                assert "--constraint=cpu --qos=${params.nersc_qos}" in output
                assert "podman.registry = 'docker.io'" in output
                assert "podman.runOptions = '--userns=keep-id --user $(id -u):$(id -g)'" in output
        wrapper = ROOT / "scripts/nersc-podman.sh"
        fake = folder / "podman-hpc"
        fake.write_text('#!/bin/sh\nprintf "%s\\n" "$@"\n')
        fake.chmod(0o755)
        output = subprocess.check_output(
            [str(wrapper), "podman", "info", "two words"],
            env={
                **os.environ,
                "PATH": f"{folder}:{os.environ['PATH']}",
                "GENESIS_RUNTIME_DIR": str(folder / "shared-runtime"),
            },
            text=True,
        )
        assert output.splitlines() == ["info", "two words"]
        assert not (folder / "podman").exists()
        shared = folder / "shared-runtime"
        failed = subprocess.run(
            [str(wrapper), "false"],
            env={
                **os.environ,
                "PATH": f"{folder}:{os.environ['PATH']}",
                "GENESIS_RUNTIME_DIR": str(shared),
            },
            capture_output=True,
            text=True,
        )
        assert failed.returncode == 1 and "retained for recovery" in failed.stderr
        assert len(list(shared.glob("podman.*/podman"))) == 1
    print("PASS: local default, seven parsed deployment profiles and scoped Podman-HPC adapter")


if __name__ == "__main__":
    main()
