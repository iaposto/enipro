"""Load and validate project configuration from params.yaml."""

import sys
from pathlib import Path

import yaml


def load_config(config_path: str | None = None) -> dict:
    """Load params.yaml and resolve paths."""
    if config_path is None:
        config_path = Path(__file__).resolve().parent.parent / "config" / "params.yaml"
    else:
        config_path = Path(config_path)

    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    # Validate required paths exist
    input_prefix = Path(cfg["input_prefix"])
    for ext in (".bed", ".bim", ".fam"):
        p = input_prefix.with_suffix(ext)
        if not p.exists():
            print(f"WARNING: input file not found: {p}", file=sys.stderr)

    return cfg


def get_breed_colors(cfg: dict) -> dict[str, str]:
    """Return {breed_code: hex_color} mapping."""
    return {code: info["color"] for code, info in cfg["breeds"].items()}


def get_breed_labels(cfg: dict) -> dict[str, str]:
    """Return {breed_code: display_label} mapping."""
    return {code: info["label"] for code, info in cfg["breeds"].items()}
