#!/usr/bin/env python3
"""Compare hydro freeze-out fluxes with sampler event averages."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Dict, List


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Check hydro-sampler conservation of E, P and B."
    )
    parser.add_argument(
        "--fo",
        default="output/fo_flux_summary.dat",
        help="Hydro freeze-out summary file.",
    )
    parser.add_argument(
        "--oscar",
        default="output/sampler_test/particle_lists.oscar",
        help="Sampler OSCAR particle list.",
    )
    parser.add_argument(
        "--outdir",
        default="output/conservation_test",
        help="Directory for analysis outputs.",
    )
    return parser.parse_args()


def baryon_number_from_pdg(pdg: int) -> int:
    apdg = abs(pdg)
    if apdg < 1000:
        return 0
    # Ion PDG code: 10LZZZAAAI, baryon number is A.
    if apdg >= 1_000_000_000:
        a_mass = (apdg // 10) % 1000
        return int(math.copysign(a_mass, pdg))
    # Hadron baryons/antibaryons (sufficient for sampler output).
    if apdg < 10_000_000:
        return int(math.copysign(1, pdg))
    return 0


def mean_and_se(values: List[float]) -> Dict[str, float]:
    n = len(values)
    if n == 0:
        return {"mean": 0.0, "se": 0.0, "std": 0.0}
    mean = sum(values) / n
    if n == 1:
        return {"mean": mean, "se": 0.0, "std": 0.0}
    var = sum((x - mean) ** 2 for x in values) / (n - 1)
    std = math.sqrt(var)
    return {"mean": mean, "se": std / math.sqrt(n), "std": std}


def parse_fo_summary(path: Path) -> Dict[str, float]:
    rows: List[List[str]] = []
    with path.open("r", encoding="utf-8") as fobj:
        for line in fobj:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            toks = line.split()
            if len(toks) < 17:
                continue
            rows.append(toks)
    if not rows:
        raise RuntimeError(f"No data rows found in {path}")

    last = rows[-1]
    return {
        "tag": last[0],
        "ncyc": int(float(last[1])),
        "t_fm": float(last[2]),
        "cycP0": float(last[3]),
        "cycPx": float(last[4]),
        "cycPy": float(last[5]),
        "cycPz": float(last[6]),
        "cycB": float(last[7]),
        "accP0": float(last[8]),
        "accPx": float(last[9]),
        "accPy": float(last[10]),
        "accPz": float(last[11]),
        "accB": float(last[12]),
        "Eneg": float(last[13]),
        "Epos": float(last[14]),
        "Bneg": float(last[15]),
        "Bpos": float(last[16]),
    }


def parse_oscar_events(path: Path) -> List[Dict[str, float]]:
    events: List[Dict[str, float]] = []
    active = False
    current = {"E": 0.0, "Px": 0.0, "Py": 0.0, "Pz": 0.0, "B": 0.0}

    with path.open("r", encoding="utf-8") as fobj:
        for raw in fobj:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#"):
                if line.startswith("# event ") and " out " in line:
                    active = True
                    current = {"E": 0.0, "Px": 0.0, "Py": 0.0, "Pz": 0.0, "B": 0.0}
                elif line.startswith("# event ") and " end " in line and active:
                    events.append(current.copy())
                    active = False
                continue

            if not active:
                continue
            cols = line.split()
            if len(cols) < 10:
                continue
            p0 = float(cols[5])
            px = float(cols[6])
            py = float(cols[7])
            pz = float(cols[8])
            pdg = int(float(cols[9]))
            current["E"] += p0
            current["Px"] += px
            current["Py"] += py
            current["Pz"] += pz
            current["B"] += baryon_number_from_pdg(pdg)

    return events


def write_event_csv(path: Path, events: List[Dict[str, float]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fobj:
        writer = csv.writer(fobj)
        writer.writerow(["event", "E_GeV", "Px_GeV", "Py_GeV", "Pz_GeV", "B"])
        for ievt, ev in enumerate(events):
            writer.writerow([ievt, ev["E"], ev["Px"], ev["Py"], ev["Pz"], ev["B"]])


def main() -> None:
    args = parse_args()
    fo_path = Path(args.fo)
    oscar_path = Path(args.oscar)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    hydro = parse_fo_summary(fo_path)
    events = parse_oscar_events(oscar_path)
    if not events:
        raise RuntimeError(f"No events parsed from {oscar_path}")

    e_vals = [ev["E"] for ev in events]
    px_vals = [ev["Px"] for ev in events]
    py_vals = [ev["Py"] for ev in events]
    pz_vals = [ev["Pz"] for ev in events]
    b_vals = [ev["B"] for ev in events]

    e_stat = mean_and_se(e_vals)
    px_stat = mean_and_se(px_vals)
    py_stat = mean_and_se(py_vals)
    pz_stat = mean_and_se(pz_vals)
    b_stat = mean_and_se(b_vals)

    # Hydro fluxes are in MeV in this code; sampler is in GeV.
    e_fo = hydro["accP0"] / 1000.0
    px_fo = hydro["accPx"] / 1000.0
    py_fo = hydro["accPy"] / 1000.0
    pz_fo = hydro["accPz"] / 1000.0
    b_fo = hydro["accB"]

    delta_e = 100.0 * (e_stat["mean"] - e_fo) / e_fo if e_fo != 0.0 else 0.0
    delta_b = 100.0 * (b_stat["mean"] - b_fo) / b_fo if b_fo != 0.0 else 0.0
    px_norm = 100.0 * px_stat["mean"] / e_fo if e_fo != 0.0 else 0.0
    py_norm = 100.0 * py_stat["mean"] / e_fo if e_fo != 0.0 else 0.0
    pz_norm = 100.0 * pz_stat["mean"] / e_fo if e_fo != 0.0 else 0.0

    summary = {
        "inputs": {"fo": str(fo_path), "oscar": str(oscar_path)},
        "hydro_reference": {
            "ncyc": hydro["ncyc"],
            "t_fm": hydro["t_fm"],
            "E_GeV": e_fo,
            "Px_GeV": px_fo,
            "Py_GeV": py_fo,
            "Pz_GeV": pz_fo,
            "B": b_fo,
        },
        "sampler_stats": {
            "N_events": len(events),
            "E_GeV": e_stat,
            "Px_GeV": px_stat,
            "Py_GeV": py_stat,
            "Pz_GeV": pz_stat,
            "B": b_stat,
        },
        "deltas_percent": {
            "delta_E": delta_e,
            "delta_B": delta_b,
            "Px_over_Efo": px_norm,
            "Py_over_Efo": py_norm,
            "Pz_over_Efo": pz_norm,
        },
    }

    with (outdir / "summary.json").open("w", encoding="utf-8") as fobj:
        json.dump(summary, fobj, indent=2, sort_keys=True)

    write_event_csv(outdir / "event_observables.csv", events)

    with (outdir / "table_conservation.csv").open("w", encoding="utf-8", newline="") as fobj:
        writer = csv.writer(fobj)
        writer.writerow(
            [
                "N_ev",
                "E_mean_GeV",
                "E_se_GeV",
                "delta_E_percent",
                "Px_over_Efo_percent",
                "Py_over_Efo_percent",
                "Pz_over_Efo_percent",
                "B_mean",
                "B_se",
                "delta_B_percent",
            ]
        )
        writer.writerow(
            [
                len(events),
                e_stat["mean"],
                e_stat["se"],
                delta_e,
                px_norm,
                py_norm,
                pz_norm,
                b_stat["mean"],
                b_stat["se"],
                delta_b,
            ]
        )

    with (outdir / "hydro_reference.tex").open("w", encoding="utf-8") as fobj:
        fobj.write(
            "\\begin{tabular}{lccccc}\n"
            "\\toprule\n"
            "Source & $E_{\\mathrm{FO}}$ & $P_{x,\\mathrm{FO}}$ & "
            "$P_{y,\\mathrm{FO}}$ & $P_{z,\\mathrm{FO}}$ & $B_{\\mathrm{FO}}$ \\\\\n"
            "\\midrule\n"
            f"Hydro freeze-out flux & {e_fo:.8g} & {px_fo:.8g} & {py_fo:.8g} & "
            f"{pz_fo:.8g} & {b_fo:.8g} \\\\\n"
            "\\bottomrule\n"
            "\\end{tabular}\n"
        )

    with (outdir / "sampler_conservation.tex").open("w", encoding="utf-8") as fobj:
        fobj.write(
            "\\begin{tabular}{lcccccc}\n"
            "\\toprule\n"
            "$N_{\\mathrm{ev}}$ & $\\langle E\\rangle$ & $\\delta_E$ [\\%] & "
            "$100\\langle P_x\\rangle/E_{\\mathrm{FO}}$ [\\%] & "
            "$100\\langle P_y\\rangle/E_{\\mathrm{FO}}$ [\\%] & "
            "$100\\langle P_z\\rangle/E_{\\mathrm{FO}}$ [\\%] & "
            "$\\delta_B$ [\\%] \\\\\n"
            "\\midrule\n"
            f"{len(events)} & {e_stat['mean']:.8g} $\\pm$ {e_stat['se']:.3g} & "
            f"{delta_e:.6g} & {px_norm:.6g} & {py_norm:.6g} & {pz_norm:.6g} & "
            f"{delta_b:.6g} \\\\\n"
            "\\bottomrule\n"
            "\\end{tabular}\n"
        )

    print("Hydro reference (GeV):")
    print(f"  E={e_fo:.8g} Px={px_fo:.8g} Py={py_fo:.8g} Pz={pz_fo:.8g} B={b_fo:.8g}")
    print("Sampler means:")
    print(
        f"  <E>={e_stat['mean']:.8g} +- {e_stat['se']:.3g}, "
        f"<B>={b_stat['mean']:.8g} +- {b_stat['se']:.3g}"
    )
    print("Relative deltas [%]:")
    print(
        f"  delta_E={delta_e:.6g}, Px/Efo={px_norm:.6g}, "
        f"Py/Efo={py_norm:.6g}, Pz/Efo={pz_norm:.6g}, delta_B={delta_b:.6g}"
    )
    print(f"Wrote outputs in: {outdir}")


if __name__ == "__main__":
    main()
