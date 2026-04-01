#!/usr/bin/env python3
"""Master runner for the full QC + PCA pipeline."""

from __future__ import annotations

import os
import sys
from datetime import datetime

from _common import load_context, run_command


STEP_SCRIPTS = [
    "01_pre_qc_filter.py",
    "02_snp_qc.py",
    "03_sample_qc.py",
    "04_relatedness.py",
    "05_qc_summary.py",
    "06_pca.py",
]


def main() -> int:
    """Execute all QC/PCA steps and generate the reports."""
    os.environ["RUN_STAMP"] = datetime.now().strftime("%Y%m%d_%H%M%S")
    ctx = load_context()

    print("=========================================")
    print(" enipro Phase 1: QC + PCA Pipeline")
    print(f" Run: {ctx.run_stamp}")
    print(f" Output: {ctx.run_dir}")
    print("=========================================")

    for step_script in STEP_SCRIPTS:
        run_command(
            [sys.executable, str(ctx.script_dir / step_script)],
            env=ctx.env,
        )
        print()

    print("=========================================")
    print(" Pipeline complete. Generating reports...")
    print("=========================================")

    run_command(
        [sys.executable, str(ctx.repo_dir / "src" / "qc_report.py"), str(ctx.config_path)],
        env=ctx.env,
    )
    run_command(
        [sys.executable, str(ctx.repo_dir / "src" / "plot_pca.py"), str(ctx.config_path)],
        env=ctx.env,
    )

    print("All done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
