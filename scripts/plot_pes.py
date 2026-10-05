#!/usr/bin/env python
"""plot_pes.py - render the PES scan grid written by probe_pes.f90.

Usage:  python plot_pes.py [stem]        (stem.dat -> stem.png; default pes_scan)

1D scan  -> a curve of V(s1) [eV] with a kcal/mol secondary axis and the
            minimum marked; 2D scan -> a filled-contour map with line
            contours, a colorbar [eV] and the grid minimum marked. The
            header metadata (# LABELx / TITLE / CONTAINER / FD_MAX) is
            parsed from the dat; nothing else is passed in.
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# eV -> kcal/mol (the legislated conversion literal of the library, used verbatim)
EV_KCAL = 23.0605


def read_dat(path):
    meta = {}
    rows = []
    with open(path) as fh:
        for raw in fh:
            if raw.startswith("#"):
                parts = raw[1:].split(None, 1)
                if parts:
                    meta.setdefault(parts[0], parts[1].strip() if len(parts) > 1 else "")
                continue
            if raw.strip():
                rows.append([float(x) for x in raw.split()])
    if not rows:
        raise SystemExit(f"plot_pes: fatal: no data rows in {path}")
    return meta, np.array(rows)


def main():
    stem = sys.argv[1] if len(sys.argv) > 1 else "pes_scan"
    meta, data = read_dat(stem + ".dat")
    naxes = int(meta.get("NAXES", "1"))
    lab1 = meta.get("LABEL1", "s1 [A]")
    lab2 = meta.get("LABEL2", "s2 [A]")
    title = meta.get("TITLE", stem)
    container = meta.get("CONTAINER", "")
    fd_note = meta.get("FD_MAX", "")

    if naxes == 1:
        s1, v = data[:, 0], data[:, 1]
        imin = int(np.nanargmin(v))
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(s1, v, lw=2.0, color="#1a5fb4")
        ax.plot(s1[imin], v[imin], "*", ms=16, color="#e01b24", zorder=5,
                label=f"min {v[imin]:.4f} eV at {lab1.split()[0]} = {s1[imin]:.3f}")
        ax.set_xlabel(lab1)
        ax.set_ylabel("V [eV]")
        ax.grid(alpha=0.3)
        ax.legend(loc="best", fontsize=9)
        sec = ax.secondary_yaxis(
            "right", functions=(lambda e: e * EV_KCAL, lambda k: k / EV_KCAL))
        sec.set_ylabel("[kcal/mol]")
    else:
        n1, n2 = int(meta["N1"]), int(meta["N2"])
        if data.shape[0] != n1 * n2:
            raise SystemExit(
                f"plot_pes: fatal: expected {n1*n2} rows, found {data.shape[0]}")
        s2 = data[:, 0].reshape(n2, n1)   # outer (slow) axis
        s1 = data[:, 1].reshape(n2, n1)   # inner (fast) axis
        v = data[:, 2].reshape(n2, n1)
        jmin, imin = np.unravel_index(int(np.nanargmin(v)), v.shape)
        fig, ax = plt.subplots(figsize=(8.5, 6.5))
        vmin = np.nanmin(v)
        vmax = np.nanmax(v)
        if not np.isfinite(vmin) or not np.isfinite(vmax):
            raise SystemExit("plot_pes: fatal: the grid carries no finite V values")
        # robust color scale: a steep repulsive corner must not wash out
        # the interesting window - clip at the 99.5th percentile when the
        # max towers over it, and say so in the footer
        clip = float(np.nanpercentile(v, 99.5))
        clipped = vmax > clip + 0.1 * (clip - vmin)
        vmax_plot = clip if clipped else vmax
        levels = np.linspace(vmin, vmax_plot, 41)
        cf = ax.contourf(s1, s2, v, levels=levels, cmap="viridis",
                         extend="max" if clipped else "neither")
        cs = ax.contour(s1, s2, np.clip(v, None, vmax_plot), levels=13,
                        colors="black", linewidths=0.5, alpha=0.6)
        ax.clabel(cs, fmt="%.2f", fontsize=6)
        fig.colorbar(cf, ax=ax, label="V [eV]")
        ax.plot(s1[jmin, imin], s2[jmin, imin], "*", ms=16, color="#e01b24",
                label=f"grid min {v[jmin, imin]:.4f} eV")
        ax.set_xlabel(lab1)
        ax.set_ylabel(lab2)
        ax.legend(loc="upper right", fontsize=9)

    nfin = int(np.isnan(data[:, 2]).sum()) if data.shape[1] >= 3 else 0
    ax.set_title(title, fontsize=11)
    foot = container + ("   |   CD check: " + fd_note.split(";", 1)[-1].strip()
                        if fd_note and fd_note != "none (FD_CHECK = F)" else "")
    if nfin:
        foot += f"   |   {nfin} non-finite points"
    if naxes == 2 and clipped:
        foot += f"   |   color scale clipped at {vmax_plot:.2f} eV (raw max {vmax:.1f})"
    fig.text(0.5, 0.005, foot, ha="center", fontsize=7, color="0.35")
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    out = stem + ".png"
    fig.savefig(out, dpi=160)
    print(f"plot_pes: wrote {out}")


if __name__ == "__main__":
    main()
