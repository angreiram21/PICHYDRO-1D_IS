#!/usr/bin/env python3
"""Compare particle distributions before/after tau->t mapping.

This reproduces the same kinematic mapping logic used in fluidization_smash.f90,
but does not perform grid deposition/fluidization.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def momentum_rapidity(e: np.ndarray, pz: np.ndarray) -> np.ndarray:
    tiny = 1e-12
    num = e + pz
    den = e - pz
    y = np.empty_like(e)
    bad = (num <= tiny) | (den <= tiny)
    y[~bad] = 0.5 * np.log(num[~bad] / den[~bad])
    y[bad] = np.where(pz[bad] >= 0.0, 1e3, -1e3)
    return y


def parse_line(parts: list[str]):
    # OSCAR2013Extended
    if len(parts) >= 22:
        vals = list(map(float, parts[:22]))
        t, x, y, z, _m, e, px, py, pz = vals[:9]
        pdg = vals[9]
        ncoll = vals[12]
        baryon = vals[20]
        return {"kind": "oscar", "t": t, "x": x, "y": y, "z": z, "e": e, "px": px, "py": py, "pz": pz, "b": baryon, "pdg": pdg, "ncoll": ncoll}

    # SMASH_IC.dat (tau x y eta mt px py Rap ...)
    if len(parts) >= 12:
        tau, x, y, eta, mt, px, py, rap, pdg, _q, baryon, _s = map(float, parts[:12])
        t = tau * np.cosh(eta)
        z = tau * np.sinh(eta)
        e = mt * np.cosh(rap)
        pz = mt * np.sinh(rap)
        return {"kind": "milne", "t": t, "x": x, "y": y, "z": z, "e": e, "px": px, "py": py, "pz": pz, "b": baryon, "pdg": pdg, "ncoll": -1.0}

    # Legacy fallback: tau eta x y px py pz E baryon
    if len(parts) >= 9:
        tau, eta, x, y, px, py, pz, e, baryon = map(float, parts[:9])
        t = tau * np.cosh(eta)
        z = tau * np.sinh(eta)
        return {"kind": "legacy", "t": t, "x": x, "y": y, "z": z, "e": e, "px": px, "py": py, "pz": pz, "b": baryon, "pdg": 0.0, "ncoll": -1.0}

    return None


def load_particles(path: Path, apply_oscar_spectator_filter: bool, ycut: float):
    rows = []
    dropped_spectators = 0
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            row = parse_line(parts)
            if row is None:
                continue

            if apply_oscar_spectator_filter and row["kind"] == "oscar":
                pdg = abs(row["pdg"])
                if abs(row["b"]) >= 0.5 and pdg in (2212.0, 2112.0) and row["ncoll"] <= 0.5:
                    y = momentum_rapidity(np.array([row["e"]]), np.array([row["pz"]]))[0]
                    if abs(y) >= ycut:
                        dropped_spectators += 1
                        continue

            rows.append(row)

    if not rows:
        raise RuntimeError("No particles loaded after filtering")

    arr = {k: np.array([r[k] for r in rows], dtype=float) for k in ["t", "x", "y", "z", "e", "px", "py", "pz", "b"]}
    return arr, dropped_spectators


def apply_mapping(arr: dict[str, np.ndarray], target_t: float, drop_if_t_above_target: bool):
    x = arr["x"].copy()
    y = arr["y"].copy()
    z = arr["z"].copy()
    t = arr["t"].copy()
    e = arr["e"]
    px = arr["px"]
    py = arr["py"]
    pz = arr["pz"]

    valid_e = e > 1e-20
    vx = np.zeros_like(e)
    vy = np.zeros_like(e)
    vz = np.zeros_like(e)
    vx[valid_e] = px[valid_e] / e[valid_e]
    vy[valid_e] = py[valid_e] / e[valid_e]
    vz[valid_e] = pz[valid_e] / e[valid_e]

    v2 = vx * vx + vy * vy + vz * vz
    bad = v2 >= 1.0
    if np.any(bad):
        scale = 0.999999 / np.sqrt(v2[bad])
        vx[bad] *= scale
        vy[bad] *= scale
        vz[bad] *= scale

    dt = target_t - t
    keep = valid_e.copy()
    if drop_if_t_above_target:
        keep &= dt >= 0.0

    x2 = x + vx * dt
    y2 = y + vy * dt
    z2 = z + vz * dt
    t2 = np.full_like(t, target_t)

    out = {
        "x": x2[keep],
        "y": y2[keep],
        "z": z2[keep],
        "t": t2[keep],
        "e": e[keep],
        "b": arr["b"][keep],
    }
    stats = {
        "n_in": len(t),
        "n_out": int(np.sum(keep)),
        "dropped_time_or_e": int(len(t) - np.sum(keep)),
    }
    return out, stats


def weighted_hist2d(x, y, w, bins=220, xr=None, yr=None):
    h, xed, yed = np.histogram2d(x, y, bins=bins, range=[xr, yr], weights=w)
    return h.T, xed, yed


def main():
    ap = argparse.ArgumentParser(description="Compare pre/post particle mapping")
    ap.add_argument("-i", "--input", required=True, help="SMASH particle file")
    ap.add_argument("-o", "--output", default="plots/particle_mapping_comparison.png", help="Output image")
    ap.add_argument("--target-t", type=float, default=-1.0, help="Target t; <=0 uses max(t)")
    ap.add_argument("--ycut", type=float, default=1.5, help="|y| cut for OSCAR spectator filter")
    ap.add_argument("--no-spectator-filter", action="store_true", help="Disable OSCAR spectator filter")
    ap.add_argument("--bins", type=int, default=220)
    args = ap.parse_args()

    inpath = Path(args.input)
    arr, n_spec = load_particles(inpath, apply_oscar_spectator_filter=not args.no_spectator_filter, ycut=args.ycut)

    target_t = args.target_t if args.target_t > 0 else float(np.max(arr["t"]))
    mapped, st = apply_mapping(arr, target_t, drop_if_t_above_target=True)

    # Common ranges
    x_all = np.concatenate([arr["x"], mapped["x"]])
    y_all = np.concatenate([arr["y"], mapped["y"]])
    z_all = np.concatenate([arr["z"], mapped["z"]])
    xr = [np.percentile(x_all, 0.5), np.percentile(x_all, 99.5)]
    yr = [np.percentile(y_all, 0.5), np.percentile(y_all, 99.5)]
    zr = [np.percentile(z_all, 0.5), np.percentile(z_all, 99.5)]

    hxy_pre, xed, yed = weighted_hist2d(arr["x"], arr["y"], arr["e"], bins=args.bins, xr=xr, yr=yr)
    hxy_post, _, _ = weighted_hist2d(mapped["x"], mapped["y"], mapped["e"], bins=args.bins, xr=xr, yr=yr)
    hxz_pre, xed2, zed = weighted_hist2d(arr["x"], arr["z"], arr["e"], bins=args.bins, xr=xr, yr=zr)
    hxz_post, _, _ = weighted_hist2d(mapped["x"], mapped["z"], mapped["e"], bins=args.bins, xr=xr, yr=zr)

    fig, ax = plt.subplots(2, 2, figsize=(12, 10), constrained_layout=True)

    im = ax[0, 0].imshow(hxy_pre, origin="lower", extent=[xed[0], xed[-1], yed[0], yed[-1]], aspect="auto", cmap="magma")
    ax[0, 0].set_title("Pre-mapping particles: x-y (E-weighted)")
    ax[0, 0].set_xlabel("x")
    ax[0, 0].set_ylabel("y")
    fig.colorbar(im, ax=ax[0, 0], label="sum(E)")

    im = ax[0, 1].imshow(hxy_post, origin="lower", extent=[xed[0], xed[-1], yed[0], yed[-1]], aspect="auto", cmap="magma")
    ax[0, 1].set_title("Post-mapping particles: x-y (E-weighted)")
    ax[0, 1].set_xlabel("x")
    ax[0, 1].set_ylabel("y")
    fig.colorbar(im, ax=ax[0, 1], label="sum(E)")

    im = ax[1, 0].imshow(hxz_pre, origin="lower", extent=[xed2[0], xed2[-1], zed[0], zed[-1]], aspect="auto", cmap="magma")
    ax[1, 0].set_title("Pre-mapping particles: x-z (E-weighted)")
    ax[1, 0].set_xlabel("x")
    ax[1, 0].set_ylabel("z")
    fig.colorbar(im, ax=ax[1, 0], label="sum(E)")

    im = ax[1, 1].imshow(hxz_post, origin="lower", extent=[xed2[0], xed2[-1], zed[0], zed[-1]], aspect="auto", cmap="magma")
    ax[1, 1].set_title("Post-mapping particles: x-z (E-weighted)")
    ax[1, 1].set_xlabel("x")
    ax[1, 1].set_ylabel("z")
    fig.colorbar(im, ax=ax[1, 1], label="sum(E)")

    fig.suptitle(
        f"Particle mapping comparison ({inpath.name})\n"
        f"N: {st['n_in']} -> {st['n_out']} | dropped(spec): {n_spec} | dropped(time/E): {st['dropped_time_or_e']} | target_t={target_t:.3f}"
    )

    outpath = Path(args.output)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=180)

    print(f"Saved: {outpath}")
    print(f"N pre/post: {st['n_in']} -> {st['n_out']}; dropped spectators: {n_spec}; dropped time/E: {st['dropped_time_or_e']}; target_t={target_t:.6f}")


if __name__ == "__main__":
    main()
