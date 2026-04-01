#!/usr/bin/env python3
"""Shared helpers for Python pipeline scripts."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Sequence


REPO_DIR = Path(__file__).resolve().parent.parent
if str(REPO_DIR) not in sys.path:
    sys.path.insert(0, str(REPO_DIR))

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(line_buffering=True)
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(line_buffering=True)

from src.config import load_config


PLINK_PREFIX_EXTENSIONS = (".bed", ".bim", ".fam")


def _stringify(value: object) -> str:
    return str(value)


@dataclass(frozen=True)
class PipelineContext:
    """Resolved pipeline configuration shared by all step scripts."""

    repo_dir: Path
    script_dir: Path
    config_path: Path
    config: dict
    input_prefix: Path
    results_dir: Path
    log_dir: Path
    run_stamp: str
    run_dir: Path
    chr_set_args: tuple[str, ...]
    autosomes: str
    qc_geno: str
    qc_mind: str
    qc_maf: str
    qc_hwe: str
    qc_king: str
    ld_window: str
    ld_step: str
    ld_r2: str
    n_pcs: str
    fst_blocksize: str

    @property
    def env(self) -> dict[str, str]:
        env = os.environ.copy()
        env["RUN_STAMP"] = self.run_stamp
        env["ENIPRO_RUN_DIR"] = str(self.run_dir)
        return env


def load_context(run_stamp: str | None = None) -> PipelineContext:
    """Load config, resolve paths, and initialize the active run directory."""
    script_dir = Path(__file__).resolve().parent
    repo_dir = script_dir.parent
    config_path = repo_dir / "config" / "params.yaml"
    config = load_config(config_path)

    results_dir = Path(config["results_dir"])
    log_dir = Path(config["log_dir"])
    resolved_run_stamp = (
        run_stamp
        or os.environ.get("RUN_STAMP")
        or datetime.now().strftime("%Y%m%d_%H%M%S")
    )
    run_dir = results_dir / resolved_run_stamp
    run_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(config_path, run_dir / "params.yaml")
    log_dir.mkdir(parents=True, exist_ok=True)

    os.environ["RUN_STAMP"] = resolved_run_stamp
    os.environ["ENIPRO_RUN_DIR"] = str(run_dir)

    chr_set_args = [str(config["chr_set"])]
    if config.get("no_xy"):
        chr_set_args.append("no-xy")

    qc_cfg = config["qc"]
    ld_cfg = config["ld_prune"]
    pca_cfg = config["pca"]
    fst_cfg = config["fst"]

    return PipelineContext(
        repo_dir=repo_dir,
        script_dir=script_dir,
        config_path=config_path,
        config=config,
        input_prefix=Path(config["input_prefix"]),
        results_dir=results_dir,
        log_dir=log_dir,
        run_stamp=resolved_run_stamp,
        run_dir=run_dir,
        chr_set_args=tuple(chr_set_args),
        autosomes=str(config["autosomes"]),
        qc_geno=_stringify(qc_cfg["geno"]),
        qc_mind=_stringify(qc_cfg["mind"]),
        qc_maf=_stringify(qc_cfg["maf"]),
        qc_hwe=_stringify(qc_cfg["hwe"]),
        qc_king=_stringify(qc_cfg["king_cutoff"]),
        ld_window=_stringify(ld_cfg["window"]),
        ld_step=_stringify(ld_cfg["step"]),
        ld_r2=_stringify(ld_cfg["r2"]),
        n_pcs=_stringify(pca_cfg["n_pcs"]),
        fst_blocksize=_stringify(fst_cfg["blocksize"]),
    )


def ensure_dir(path: Path) -> Path:
    """Create a directory tree if needed and return the path."""
    path.mkdir(parents=True, exist_ok=True)
    return path


def run_command(
    command: Sequence[str | Path],
    *,
    env: dict[str, str] | None = None,
    cwd: Path | None = None,
) -> None:
    """Execute a subprocess and stream output directly to the terminal."""
    subprocess.run(
        [str(part) for part in command],
        check=True,
        env=env,
        cwd=str(cwd) if cwd is not None else None,
    )


def line_count(path: Path) -> int:
    """Count newline-delimited records in a text file."""
    with path.open() as handle:
        return sum(1 for _ in handle)


def count_king_removed(path: Path) -> int:
    """Count sample rows in a KING cutoff output file, excluding the header."""
    if not path.exists():
        return 0
    return max(line_count(path) - 1, 0)


def copy_plink_prefix(source_prefix: Path, destination_prefix: Path) -> None:
    """Copy a PLINK binary trio from one prefix to another."""
    destination_prefix.parent.mkdir(parents=True, exist_ok=True)
    for extension in PLINK_PREFIX_EXTENSIONS:
        shutil.copy2(
            source_prefix.with_suffix(extension),
            destination_prefix.with_suffix(extension),
        )
