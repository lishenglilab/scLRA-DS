"""Run legacy-named implementations through the provider-neutral CLI paths."""

from __future__ import annotations

import runpy
import sys
from pathlib import Path


def run_implementation(filename: str) -> None:
    bin_dir = Path(__file__).resolve().parents[1]
    implementation_dir = bin_dir / "deepseek"
    sys.path.insert(0, str(bin_dir))
    sys.path.insert(0, str(implementation_dir))
    runpy.run_path(str(implementation_dir / filename), run_name="__main__")
