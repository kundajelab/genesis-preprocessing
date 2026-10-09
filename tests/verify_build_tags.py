"""Check build-script tag handling without invoking Pixi, Docker, or a registry."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="build-tags-") as temporary:
        root = Path(temporary)
        shutil.copytree(ROOT / "scripts", root / "scripts")
        (root / "dockers").mkdir()
        (root / "dockers/.dockerignore").touch()
        (root / "environments").mkdir()
        (root / "environments/EXAMPLE.yaml").write_text("name: example\n")
        (root / "project").mkdir()
        (root / "project/pixi.lock").touch()
        binaries = root / "bin"
        binaries.mkdir()
        mock = binaries / "mock"
        mock.write_text(
            f"#!{sys.executable}\n"
            "import json, os, sys\nfrom pathlib import Path\n"
            "tool = Path(sys.argv[0]).name\n"
            "with open(os.environ['BUILD_LOG'], 'a') as f:\n"
            "    f.write(json.dumps([tool, *sys.argv[1:]]) + '\\n')\n"
            "if tool == 'git':\n"
            "    tag = os.environ.get('TEST_GIT_TAG')\n"
            "    if not tag: sys.exit(128)\n"
            "    print(tag)\n"
            "elif tool == 'pixi':\n"
            "    Path('pixi.toml').touch(); Path('pixi.lock').touch()\n"
            "elif tool == 'docker' and sys.argv[1] == 'images':\n"
            "    print('REPOSITORY TAG DIGEST\\nexample test sha256:fixture')\n"
        )
        mock.chmod(0o755)
        for tool in ("git", "pixi", "docker"):
            (binaries / tool).symlink_to(mock)
        log = root / "calls.jsonl"
        env = {
            **os.environ,
            "PATH": f"{binaries}:{os.environ['PATH']}",
            "BUILD_LOG": str(log),
        }
        env.pop("TEST_GIT_TAG", None)

        def check(script: str, args: list[str], expected: int, tag: str | None = None) -> None:
            log.write_text("")
            dotenv = root / ".env"
            dotenv.write_text("UNCHANGED=1\n")
            result = subprocess.run(
                [str(root / "scripts" / script), *args],
                cwd=root,
                env={**env, **({"TEST_GIT_TAG": tag} if tag else {})},
                capture_output=True,
                text=True,
                check=False,
                timeout=15,
            )
            assert result.returncode == expected, (script, args, result.stderr)
            calls = [json.loads(line) for line in log.read_text().splitlines()]
            assert not any(call[:2] == ["docker", "push"] for call in calls), calls
            if expected:
                assert "--tag" in result.stderr or "explicit tag" in result.stderr
            if expected or args == ["--help"]:
                assert not any(call[0] in ("docker", "pixi") for call in calls), calls
                assert dotenv.read_text() == "UNCHANGED=1\n"
            else:
                builds = [call for call in calls if call[:2] == ["docker", "build"]]
                assert builds and any(
                    value.endswith(":" + (tag or "explicit")) for value in builds[0]
                )
            if tag is None and expected == 0:
                assert not any(call[0] == "git" for call in calls), calls

        check("build-dockers.sh", ["--help"], 0)
        for args in (["--tag"], ["--tag", ""], ["--tag", "--no-push"]):
            check("build-dockers.sh", args, 2)
        check("build-dockers.sh", ["--no-push", "--", "example"], 2)
        check("build-dockers.sh", ["--tag", "explicit", "--no-push", "--", "example"], 0)
        check("build-dockers.sh", ["--no-push", "--", "example"], 0, "v1.2.3")
        for script, project in (
            ("build-project-docker.sh", "project"),
            ("build-yaml-docker.sh", "example"),
        ):
            check(script, ["--no-push", project], 2)
            check(script, ["--no-push", project, "explicit"], 0)
            check(script, ["--no-push", project], 0, "v1.2.3")


if __name__ == "__main__":
    main()
