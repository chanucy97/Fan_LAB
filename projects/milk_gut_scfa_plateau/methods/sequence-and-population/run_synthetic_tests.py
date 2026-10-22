"""Run isolated synthetic suites; never load participant data or start HPC jobs."""
from pathlib import Path
import argparse
import os
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
SUITES = [
    ("comparison", "test_compare_strains.py"),
    ("sequence_integrity", "test_assembly_scan.py"),
    ("sequence_integrity", "test_profile_sequence_audit.py"),
    ("sequence_integrity/enzyme_context", "test_core.py"),
    ("sequence_integrity/readback_stage2_deploy_v1", "test_core.py"),
    ("sequence_integrity/readback_stage2_deploy_v1", "test_driver.py"),
]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", PYTHONUTF8="1")
    failures = 0
    for directory, test in SUITES:
        print(f"Running {directory}/{test}", flush=True)
        result = subprocess.run([sys.executable, "-B", "-m", "unittest", "-v", test], cwd=ROOT/directory, env=env)
        failures += result.returncode != 0
    return 1 if failures else 0

if __name__ == "__main__":
    raise SystemExit(main())
