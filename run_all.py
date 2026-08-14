"""Cross-platform equivalent of the Makefile, for machines without `make`.

    python run_all.py all           every main and supplementary panel
    python run_all.py figures       main figures 1-5
    python run_all.py supplementary supplementary S1-S6
    python run_all.py published     the submitted Fig 4B/4C (wrong control clock)
    python run_all.py notebook      rebuild the walkthrough notebook and execute it
    python run_all.py test          the test suite
    python run_all.py clean         remove caches (figures/_out/ is committed)

`make` targets of the same name do the same thing; this exists because `make`
is not usually installed on Windows.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
FIGURES = ["figure1.py", "figure2.py", "figure3.py", "figure4.py", "figure5.py"]
SUPPLEMENTARY = [
    "figS1_barcode_linear.py",
    "figS2_cluster_selection.py",
    "figS3_coclustering_support.py",
    "figS4_window_sweep.py",
    "figS5_eigenvalue_spectra.py",
    "figS6_global_coclustering.py",
]


def run(script: Path, *args: str) -> None:
    rel = script.relative_to(ROOT)
    print(f"\n=== {rel} {' '.join(args)}".rstrip(), flush=True)
    result = subprocess.run([sys.executable, str(script), *args], cwd=ROOT)
    if result.returncode != 0:
        sys.exit(f"FAILED: {rel} (exit {result.returncode})")


def do_figures() -> None:
    for name in FIGURES:
        run(ROOT / "figures" / name)


def do_supplementary() -> None:
    for name in SUPPLEMENTARY:
        run(ROOT / "figures" / "supplementary" / name)


def do_published() -> None:
    run(ROOT / "figures" / "figure4.py", "--clock", "published")


def do_test() -> None:
    result = subprocess.run([sys.executable, "-m", "pytest", "tests/", "-q"], cwd=ROOT)
    sys.exit(result.returncode)


def do_notebook() -> None:
    run(ROOT / "tools" / "build_notebook.py", "--execute")


def do_clean() -> None:
    # figures/_out/ is committed, so it is not cleaned here. `git checkout
    # figures/_out` is the way back if a run leaves it dirty.
    targets = [ROOT / ".pytest_cache"]
    targets += [p for p in ROOT.rglob("__pycache__") if p.is_dir()]
    targets += [p for p in ROOT.rglob(".ipynb_checkpoints") if p.is_dir()]
    for path in targets:
        if path.exists():
            shutil.rmtree(path)
            print(f"removed {path.relative_to(ROOT)}")


TARGETS = {
    "all": lambda: (do_figures(), do_supplementary()),
    "figures": do_figures,
    "supplementary": do_supplementary,
    "published": do_published,
    "notebook": do_notebook,
    "test": do_test,
    "clean": do_clean,
}


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("target", choices=sorted(TARGETS), nargs="?", default="all")
    TARGETS[ap.parse_args().target]()
    print("\ndone.")
