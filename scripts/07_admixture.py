#!/usr/bin/env python3
#SBATCH --job-name=admixture-1.3.0-case
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=30:00:00
"""Run ADMIXTURE K=2-30 under SLURM."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(line_buffering=True)
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(line_buffering=True)


def load_module(module_name: str) -> None:
    """Load an environment module into the current Python process if possible."""
    modulecmd = shutil.which("modulecmd")
    if not modulecmd:
        return

    completed = subprocess.run(
        [modulecmd, "python", "load", module_name],
        check=True,
        capture_output=True,
        text=True,
    )
    exec_globals = {"__builtins__": __builtins__, "os": os, "sys": sys}
    exec(completed.stdout, exec_globals, {})


def timestamp() -> str:
    """Return a batch-log-friendly timestamp."""
    return datetime.now().strftime("%a %b %d %H:%M:%S %Y")


def main() -> int:
    """Run ADMIXTURE for K=2..30 and tee each log file."""
    load_module("admixture/1.3.0")

    print(f"[{timestamp()}] Running ADMIXTURE K=2-30, 10-fold CV, 16 threads, seed=42")
    for k_value in range(2, 31):
        print(f"[{timestamp()}] K={k_value}")
        with Path(f"log{k_value}.out").open("w") as handle:
            process = subprocess.Popen(
                [
                    "admixture",
                    "--cv=10",
                    "-j16",
                    "-s",
                    "42",
                    "graega_ldpruned.bed",
                    str(k_value),
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )
            assert process.stdout is not None
            for line in process.stdout:
                handle.write(line)
                sys.stdout.write(line)
            return_code = process.wait()
            if return_code != 0:
                raise subprocess.CalledProcessError(return_code, process.args)
            sys.stdout.flush()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
