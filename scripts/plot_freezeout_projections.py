#!/usr/bin/env python3
"""Plot freeze-out hypersurface projections from freezout.dat.

Reads columns:
  0: t, 1: x, 2: y, 3: z
and produces three 2D projections with geometric cuts:
  - t-z with |x|<cut and |y|<cut
  - t-x with |y|<cut and |z|<cut
  - t-y with |x|<cut and |z|<cut
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Project freezout hypersurface onto t-z, t-x, and t-y planes"
    )
    parser.add_argument(
        "-i",
        "--input",
        default="freezout.dat",
        help="Input freeze-out file (default: freezout.dat)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="freezout_projections.png",
        help="Output image path (default: freezout_projections.png)",
    )
    parser.add_argument(
        "--cut",
        type=float,
        default=0.1,
        help="Absolute cut value used in all projections (default: 0.1)",
    )
    parser.add_argument(
        "--max-points",
        type=int,
        default=250000,
        help="Maximum number of points per panel to draw (default: 250000)",
    )
    return parser.parse_args()


def load_coordinates(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    data = np.loadtxt(path, usecols=(0, 1, 2, 3), ndmin=2)
    if data.size == 0:
        raise ValueError(f"Input file is empty: {path}")

    t = data[:, 0]
    x = data[:, 1]
    y = data[:, 2]
    z = data[:, 3]
    return t, x, y, z


def maybe_downsample(a: np.ndarray, b: np.ndarray, max_points: int) -> tuple[np.ndarray, np.ndarray]:
    n = a.size
    if n <= max_points or max_points <= 0:
        return a, b
    idx = np.linspace(0, n - 1, max_points, dtype=int)
    return a[idx], b[idx]


def add_panel(
    ax: plt.Axes,
    xvals: np.ndarray,
    yvals: np.ndarray,
    title: str,
    xlabel: str,
    ylabel: str,
    max_points: int,
) -> None:
    xplot, yplot = maybe_downsample(xvals, yvals, max_points)

    if xplot.size == 0:
        ax.text(0.5, 0.5, "No points after cut", ha="center", va="center", transform=ax.transAxes)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
        return

    ax.scatter(xplot, yplot, s=2, alpha=0.5, linewidths=0)
    ax.set_title(f"{title} (N={xvals.size})")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)

    t, x, y, z = load_coordinates(input_path)
    c = abs(args.cut)

    mask_tz = (np.abs(x) < c) & (np.abs(y) < c)
    mask_tx = (np.abs(y) < c) & (np.abs(z) < c)
    mask_ty = (np.abs(x) < c) & (np.abs(z) < c)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)

    add_panel(
        axes[0],
        z[mask_tz],
        t[mask_tz],
        title=f"t-z projection |x|,|y|<{c}",
        xlabel="z",
        ylabel="t",
        max_points=args.max_points,
    )
    add_panel(
        axes[1],
        x[mask_tx],
        t[mask_tx],
        title=f"t-x projection |y|,|z|<{c}",
        xlabel="x",
        ylabel="t",
        max_points=args.max_points,
    )
    add_panel(
        axes[2],
        y[mask_ty],
        t[mask_ty],
        title=f"t-y projection |x|,|z|<{c}",
        xlabel="y",
        ylabel="t",
        max_points=args.max_points,
    )

    fig.suptitle(f"Freeze-out hypersurface projections from {input_path.name}", fontsize=12)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()
