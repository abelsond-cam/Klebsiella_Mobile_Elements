#!/usr/bin/env python3
"""Generate a limited mgefinder_dataset.txt for testing.

Thin wrapper around generate_mgefinder_dataset.py. Limiting is now native via
--limit (sublineage mode), so this just forwards --test-n as --limit. Kept for
back-compat with run_pipeline.py / Makefile callers.
"""

import argparse
import subprocess
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Generate limited MGEfinder dataset for testing")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--row", type=int, default=0)
    parser.add_argument("--test-n", type=int, required=True, help="Limit to first N samples")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--sublineage", default=None)
    args = parser.parse_args()

    print(f">>> Limited dataset: first {args.test_n} sample(s)")
    cmd = [
        "python3", "src/generate_mgefinder_dataset.py",
        "--config", str(args.config),
        "--row", str(args.row),
        "--limit", str(args.test_n),
        "--out", str(args.out),
    ]
    if args.sublineage:
        cmd += ["--sublineage", args.sublineage]
    sys.exit(subprocess.run(cmd).returncode)


if __name__ == "__main__":
    main()
