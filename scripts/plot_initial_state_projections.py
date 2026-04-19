#!/usr/bin/env python3
"""Plot initial-state projections from output/fort_118.dat.

Reads the first snapshot block (up to the first marker row with i=9999),
using columns:
  i j k e capn p w u v zmom ymom xmom tq sq amuq
and builds:
  - transverse plane: e(x,y) integrated over z
  - reaction plane:    e(x,z) integrated over y
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot initial fluidized state projections")
    p.add_argument("-i", "--input", default="output/fort_118.dat", help="Input fort_118 file")
    p.add_argument("-o", "--output", default="plots/initial_state_projections.png", help="Output image")
    p.add_argument("--dx", type=float, default=1.0, help="Cell size dx (for axis units)")
    p.add_argument("--dy", type=float, default=1.0, help="Cell size dy (for axis units)")
    p.add_argument("--dz", type=float, default=1.0, help="Cell size dz (for axis units)")
    return p.parse_args()


def load_first_snapshot(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    i_vals: list[int] = []
    j_vals: list[int] = []
    k_vals: list[int] = []
    e_vals: list[float] = []

    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            parts = s.split()
            if len(parts) < 4:
                continue
            try:
                i = int(parts[0])
                j = int(parts[1])
                k = int(parts[2])
            except ValueError:
                continue

            if i == 9999:
                break

            try:
                e = float(parts[3])
            except ValueError:
                continue

            i_vals.append(i)
            j_vals.append(j)
            k_vals.append(k)
            e_vals.append(e)

    if not i_vals:
        raise RuntimeError("No cells found in first snapshot block")

    return np.array(i_vals), np.array(j_vals), np.array(k_vals), np.array(e_vals)


def build_grid(i: np.ndarray, j: np.ndarray, k: np.ndarray, e: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, int, int, int]:
    i0 = i.min()
    j0 = j.min()
    k0 = k.min()

    ii = i - i0
    jj = j - j0
    kk = k - k0

    ni = ii.max() + 1
    nj = jj.max() + 1
    nk = kk.max() + 1

    grid = np.zeros((ni, nj, nk), dtype=np.float64)
    grid[ii, jj, kk] = e
    return grid, np.arange(ni), np.arange(nj), ni, nj, nk


def main() -> None:
    args = parse_args()
    inpath = Path(args.input)
    outpath = Path(args.output)

    i, j, k, e = load_first_snapshot(inpath)
    grid, _, _, ni, nj, nk = build_grid(i, j, k, e)

    exy = grid.sum(axis=2) * args.dz
    exz = grid.sum(axis=1) * args.dy

    x = (np.arange(ni) - 0.5 * (ni - 1)) * args.dx
    y = (np.arange(nj) - 0.5 * (nj - 1)) * args.dy
    z = (np.arange(nk) - 0.5 * (nk - 1)) * args.dz

    extent_xy = [x[0], x[-1], y[0], y[-1]]
    extent_xz = [x[0], x[-1], z[0], z[-1]]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)

    im0 = axes[0].imshow(exy.T, origin="lower", extent=extent_xy, aspect="auto", cmap="magma")
    axes[0].set_title("Initial state: transverse projection (x-y), integrated in z")
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("y")
    cb0 = fig.colorbar(im0, ax=axes[0])
    cb0.set_label(r"$\int e\,dz$")

    im1 = axes[1].imshow(exz.T, origin="lower", extent=extent_xz, aspect="auto", cmap="magma")
    axes[1].set_title("Initial state: reaction-plane projection (x-z), integrated in y")
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("z")
    cb1 = fig.colorbar(im1, ax=axes[1])
    cb1.set_label(r"$\int e\,dy$")

    fig.suptitle(f"Initial post-fluidization snapshot from {inpath.name}")

    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=180)

    print(f"Saved: {outpath}")
    print(f"Cells in first snapshot: {e.size}, grid size inferred: {ni}x{nj}x{nk}")


if __name__ == "__main__":
    main()
