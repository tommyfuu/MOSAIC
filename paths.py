"""Paths shared by the pipeline scripts; machine-specific values live in paths.yml."""
import os

REPO_ROOT = os.path.dirname(os.path.abspath(__file__))


def _read_paths_yml():
    values = {}
    with open(os.path.join(REPO_ROOT, "paths.yml")) as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if ":" in line:
                key, value = line.split(":", 1)
                values[key.strip()] = value.strip().strip("'\"")
    return values


_values = _read_paths_yml()
SIM_OUTPUT_ROOT = os.environ.get("SIM_OUTPUT_ROOT", _values["sim_output_root"]).rstrip("/")
BENCH_DIR = os.path.join(REPO_ROOT, "benchmark")
DATA_DIR = os.path.join(REPO_ROOT, "data")
