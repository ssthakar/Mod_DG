#!/usr/bin/env python3
"""Post-process the Mod_DG solution file (legacy ASCII VTK written by
grid::post_proc::writevtk_mesh).

Usage:
    python3 postproc.py                       # reads this.vtk, writes this.png
    python3 postproc.py this.vtk -o mach.png  # custom input/output
    python3 postproc.py this.vtk -f Mach p    # choose fields

Fields available: rho, u, v, E, Mach (from the file) and p (computed).
Needs numpy and matplotlib.
"""
import argparse

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

GAMMA = 1.4


def read_vtk(path):
    """Return (points[npoin,2], tris[nelem,3], {name: cell array})."""
    tok = open(path).read().split()
    i = 0
    pts = tris = None
    fields = {}
    ncell_data = None
    while i < len(tok):
        t = tok[i]
        if t == "POINTS":
            n = int(tok[i + 1])
            i += 3  # POINTS n type
            pts = np.array(tok[i:i + 3 * n], float).reshape(n, 3)[:, :2]
            i += 3 * n
        elif t == "CELLS":
            n = int(tok[i + 1])
            i += 3
            raw = np.array(tok[i:i + 4 * n], int).reshape(n, 4)
            tris = raw[:, 1:]
            i += 4 * n
        elif t == "CELL_DATA":
            ncell_data = int(tok[i + 1])
            i += 2
        elif t == "SCALARS":
            name = tok[i + 1]
            i += 4  # SCALARS name type ncomp
            if tok[i] == "LOOKUP_TABLE":
                i += 2
            fields[name] = np.array(tok[i:i + ncell_data], float)
            i += ncell_data
        else:
            i += 1
    if "rho" in fields and "E" in fields:
        rho, u, v, E = (fields[k] for k in ("rho", "u", "v", "E"))
        fields["p"] = (GAMMA - 1.0) * (E - 0.5 * rho * (u * u + v * v))
    return pts, tris, fields


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("vtk", nargs="?", default="this.vtk")
    ap.add_argument("-o", "--out", default=None, help="output image (default <vtk>.png)")
    ap.add_argument("-f", "--fields", nargs="+", default=["Mach", "p", "rho"])
    ap.add_argument("--mesh", action="store_true", help="overlay the mesh edges")
    args = ap.parse_args()

    pts, tris, fields = read_vtk(args.vtk)
    tri = mtri.Triangulation(pts[:, 0], pts[:, 1], tris)
    names = [f for f in args.fields if f in fields]
    missing = set(args.fields) - set(names)
    if missing:
        print("skipping unknown fields:", ", ".join(sorted(missing)))

    xr = np.ptp(pts[:, 0]); yr = np.ptp(pts[:, 1])
    fig, axes = plt.subplots(len(names), 1, figsize=(10, 10 * yr / xr * len(names) + 1),
                             squeeze=False, constrained_layout=True)
    for ax, name in zip(axes[:, 0], names):
        pc = ax.tripcolor(tri, facecolors=fields[name], cmap="viridis", shading="flat")
        if args.mesh:
            ax.triplot(tri, lw=0.2, color="k", alpha=0.3)
        ax.set_aspect("equal")
        ax.set_title(f"{name}   min {fields[name].min():.4g}   max {fields[name].max():.4g}")
        fig.colorbar(pc, ax=ax, shrink=0.9)
    out = args.out or args.vtk.rsplit(".", 1)[0] + ".png"
    fig.savefig(out, dpi=150)
    print("wrote", out)


if __name__ == "__main__":
    main()
