"""Known-answer checks for the experimental ATAC prototype.

The real-tool prototype regression ran through the benchmark runner, which now lives in
genesis-curation. The supported bulk ATAC workflow is validated by validate_supported_atac.py.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def command(*args: str) -> None:
    subprocess.run(args, cwd=ROOT, check=True)


def main() -> None:
    command(str(ROOT / "scripts/run-atac-prototype.sh"), "--help")
    for name in ("verify_metrics.py", "verify_parameters.py"):
        command(sys.executable, str(ROOT / "experimental/atac" / name))


if __name__ == "__main__":
    main()
