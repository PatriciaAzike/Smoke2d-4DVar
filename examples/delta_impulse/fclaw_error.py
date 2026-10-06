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
        ii, jj = int(ob["xm"] / h), int(ob["ym"] / h)
        qa, qr = Qa[jj, ii], Qr[jj, ii]
        rows.append((ob["xm"], ob["ym"], ob["dm"], qa, qr,
                     abs(qa - qr), abs(qa - ob["dm"]), abs(qr - ob["dm"])))
    obs_err = max(r[5] for r in rows)          # max AMR-vs-uniform obs error
    amr_misfit = max(r[6] for r in rows)        # max |q_AMR - dm|
    ref_misfit = max(r[7] for r in rows)        # max |q_ref - dm|

    da, dr = dof(amr), dof(ref)
    print(f"uniform grid           : {N} x {N}")
    print(f"AMR active cells (DOF)  : {da}")
    print(f"uniform cells (DOF)     : {dr}")
    print(f"cell-count reduction    : {dr / da:.2f}x")
    print(f"relative L2 error  E2   : {E2:.3e}")
    print(f"max error          Einf : {Einf:.3e}")
    print(f"AMR peak q              : {np.nanmax(Qa):.4f}")
    print(f"ref peak q              : {np.nanmax(Qr):.4f}")
    print()
    print("per-observation:")
    print(f"  {'(xm, ym)':<14}{'dm':>5}{'q_AMR':>9}{'q_ref':>9}"
          f"{'|A-R|':>10}{'|A-dm|':>10}{'|R-dm|':>10}")
    for (xm, ym, dm, qa, qr, dar, amf, rmf) in rows:
        print(f"  ({xm:.2f}, {ym:.2f}) {dm:>5.1f}{qa:>9.3f}{qr:>9.3f}"
              f"{dar:>10.2e}{amf:>10.3f}{rmf:>10.3f}")
    print()
    print(f"max obs-site error (AMR vs uniform) : {obs_err:.3e}")
    print(f"max data misfit |q_AMR - dm|        : {amr_misfit:.3f}")
    print(f"max data misfit |q_ref - dm|        : {ref_misfit:.3f}")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
