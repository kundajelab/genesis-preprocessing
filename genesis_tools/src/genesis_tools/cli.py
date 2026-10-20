"""Expose typed functions as command-line tools with defopt."""

from __future__ import annotations

import sys

import defopt

from .metadata import infer_metadata
from .nh import add_nh
from .quantification import mean_peak_rpkm, quantify, quantify_peaks


def main() -> None:
    """Dispatch commands and report actionable validation errors."""
    try:
        defopt.run(
            {
                "metadata": infer_metadata,
                "peak-rpm": quantify_peaks,
                "peak-rpkm": mean_peak_rpkm,
                "add-nh": add_nh,
                "quantify": quantify,
            },
            cli_options="all",
        )
    except (OSError, ValueError, EOFError) as error:
        print(f"Error: {error}", file=sys.stderr)
        raise SystemExit(1) from error


if __name__ == "__main__":
    main()
