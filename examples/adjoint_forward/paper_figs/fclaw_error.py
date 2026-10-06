#!/usr/bin/env python3
"""
Compare a ForestClaw AMR solution against a uniformly refined reference.

Both inputs are ForestClaw fort.q ASCII files. The AMR field is prolonged
(piecewise-constant) onto the uniform fine grid, then error norms are computed.

Usage:
    python3 fclaw_error.py  fort.q_AMR  fort.q_UNIFORM
"""
import sys
import numpy as np

LDOM = 2.0   # domain is [0, LDOM]^2
# observations: location (xm, ym) and observed value dm; edit to match your setup
OBS = [
    {"xm": 0.67, "ym": 1.0, "dm": 3.0},
    {"xm": 0.80, "ym": 1.2, "dm": 3.0},
    {"xm": 0.55, "ym": 1.2, "dm": 2.0},
]


def read_fort_q(fname):
    patches = []
    with open(fname) as f:
        lines = f.readlines()
    i, n = 0, len(lines)
    while i < n:
        if 'grid_number' in lines[i]:
            i += 1
            hdr = {}
            for key in ['AMR_level', 'block_number', 'mpi_rank', 'mx', 'my',
                        'xlow', 'ylow', 'dx', 'dy']:
                hdr[key] = lines[i].split()[0]
                i += 1
            mx, my = int(hdr['mx']), int(hdr['my'])
            vals = []
            while len(vals) < mx * my and i < n:
                s = lines[i].strip()
                if s:
                    vals.append(float(s))
                i += 1
            patches.append(dict(
                level=int(hdr['AMR_level']), mx=mx, my=my,
                xlow=float(hdr['xlow']), ylow=float(hdr['ylow']),
                dx=float(hdr['dx']), dy=float(hdr['dy']),
                q=np.array(vals).reshape(my, mx)))   # q[j, i]
        else:
            i += 1
    return patches


def _leaf_at(ps, x, y):
    """Finest patch whose area contains (x, y)."""
    best = None
    for p in ps:
        if (p["xlow"] <= x < p["xlow"] + p["mx"]*p["dx"] and
                p["ylow"] <= y < p["ylow"] + p["my"]*p["dy"]):
            if best is None or p["dx"] < best["dx"]:
                best = p
    return best


def _cell_value(ps, p, i, j):
    """q in cell (i, j) of patch p (0-based).  Indices one past the patch edge
    stand in for ForestClaw's ghost cells: take the value of the leaf cell
    that contains that neighbour's centre."""
    if 0 <= i < p["mx"] and 0 <= j < p["my"]:
        return p["q"][j, i]
    x = p["xlow"] + (i + 0.5)*p["dx"]
    y = p["ylow"] + (j + 0.5)*p["dy"]
    x %= LDOM; y %= LDOM                        # periodic domain
    n = _leaf_at(ps, x, y)
    ii = int((x - n["xlow"])/n["dx"]); jj = int((y - n["ylow"])/n["dy"])
    return n["q"][jj, ii]


def forestclaw_gauge(ps, x, y):
    """This uses the same formula as ForestClaw's fclaw2d_clawpatch46_fort_gauges_update:
    find the leaf patch containing the gauge, take the cell (i, j) containing
    it, and interpolate bilinearly with cells i+1 and j+1 using the offsets of
    the gauge from the centre of (i, j).  Reproduces the gauge values printed
    by adjoint_forward (e.g. q_F = -0.0006, 0.4981, 0.9999 at t = 0.4)."""
    p = _leaf_at(ps, x, y)
    i = int((x - p["xlow"])/p["dx"]); j = int((y - p["ylow"])/p["dy"])
    xo = (x - (p["xlow"] + (i + 0.5)*p["dx"]))/p["dx"]
    yo = (y - (p["ylow"] + (j + 0.5)*p["dy"]))/p["dy"]
    q = lambda a, b: _cell_value(ps, p, a, b)
    return ((1 - xo)*(1 - yo)*q(i, j)     + xo*(1 - yo)*q(i + 1, j)
            + (1 - xo)*yo*q(i, j + 1)     + xo*yo*q(i + 1, j + 1))


def to_uniform(patches, N):
    """Prolong patches onto an N x N uniform grid (piecewise constant)."""
    h = LDOM / N
    Q = np.full((N, N), np.nan)          # Q[j, i]
    for p in patches:
        rx, ry = int(round(p['dx'] / h)), int(round(p['dy'] / h))
        i0, j0 = int(round(p['xlow'] / h)), int(round(p['ylow'] / h))
        for j in range(p['my']):
            for i in range(p['mx']):
                Q[j0 + j * ry:j0 + (j + 1) * ry,
                  i0 + i * rx:i0 + (i + 1) * rx] = p['q'][j, i]
    return Q


def dof(patches):
    return sum(p['mx'] * p['my'] for p in patches)


def main(amr_file, ref_file):
    amr = read_fort_q(amr_file)
    ref = read_fort_q(ref_file)

    dxmin = min(min(p['dx'] for p in amr), min(p['dx'] for p in ref))
    N = int(round(LDOM / dxmin))

    Qa = to_uniform(amr, N)
    Qr = to_uniform(ref, N)

    diff = Qa - Qr
    E2 = np.linalg.norm(diff) / np.linalg.norm(Qr)
    Einf = np.nanmax(np.abs(diff))

    h = LDOM / N

    # per-observation: AMR vs uniform difference, and misfit to the data dm
    rows = []
    for ob in OBS:
        # ForestClaw's gauge interpolation (bilinear, on the leaf patch), so
        # these match the gauge values printed by adjoint_forward.
        qa = forestclaw_gauge(amr, ob["xm"], ob["ym"])
        qr = forestclaw_gauge(ref, ob["xm"], ob["ym"])
        rows.append((ob["xm"], ob["ym"], ob["dm"], qa, qr,
                     abs(qa - qr), abs(qa - ob["dm"]), abs(qr - ob["dm"])))
    obs_err = max(r[5] for r in rows)          # max AMR-vs-uniform obs error
    amr_misfit = max(r[6] for r in rows)        # max |q_AMR - dm|
    ref_misfit = max(r[7] for r in rows)        # max |q_ref - dm|

    da, dr = dof(amr), dof(ref)
    print(f"uniform grid           : {N} x {N}")
    print(f"AMR active cells (DOF)  : {da}")
    print(f"uniform cells (DOF)     : {dr}")
    print(f"cell-count reduction    : {dr / da:.2f}x   (this frame only)")
    print(f"relative L2 error  E2   : {E2:.3e}")
    print(f"max error          Einf : {Einf:.3e}")
    print(f"AMR peak q              : {np.nanmax(Qa):.4f}")
    print(f"ref peak q              : {np.nanmax(Qr):.4f}")
    print()
    print("per-observation:")
    print(f"  {'(xm, ym)':<14}{'dm':>5}{'q_AMR':>9}{'q_ref':>9}"
          f"{'|A-R|':>10}{'|A-dm|':>10}{'|R-dm|':>10}")
    for (xm, ym, dm, qa, qr, dar, amf, rmf) in rows:
        print(f"  ({xm:.2f}, {ym:.2f}) {dm:>5.1f}{qa:>9.4f}{qr:>9.4f}"
              f"{dar:>10.2e}{amf:>10.4f}{rmf:>10.4f}")
    print()
    print(f"max obs-site error (AMR vs uniform) : {obs_err:.3e}")
    print(f"max data misfit |q_AMR - dm|        : {amr_misfit:.3f}")
    print(f"max data misfit |q_ref - dm|        : {ref_misfit:.3f}")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
