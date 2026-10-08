#!/usr/bin/env python3
"""
Reading and panel routines for the delta-impulse figures (alpha_m advecting
backward, r_m built from it). plot_delta_impulse.py imports this module, picks
the frames and sets the colour scale; this file is not run on its own.

Panels are 990 x 910 at 300 dpi with style_panel.m's axes and colorbar rectangles,
matching every other field panel in the paper. The AMR patch outlines are drawn on.

TIME CONVENTION -- this is the thing that is easy to get wrong.

The adjoint marches on a REVERSED clock. Its file frame N carries
tau = N*(tfinal/nout), and the physical time is t = tfinal - tau. The representer
runs on physical time directly. So the same physical instant is

    adjoint frame N        and        forward frame  nout - N

clock() and nearest() therefore work in PHYSICAL times: the nearest frame is
looked up in each directory separately, and panels are named by physical time, so

    alpha_t00.80.png   pairs with   rm_t00.80.png

Check: adjoint frame nout (physical t = 0) and forward frame 0 hold the same
field, the representer's initial condition r_m(.,0) = W_i^{-1} alpha_m(.,0)
(identical when W_i = 1); their mass and centroid should match.
"""
import glob, os, re, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from fclaw_error import read_fort_q                      # noqa: E402

LDOM = 2.0

# Computer Modern, to match the paper's body text. 
plt.rcParams.update({
    "mathtext.fontset": "cm",
    "font.family":      "serif",
    "font.serif":       ["cmr10", "DejaVu Serif"],
    "axes.unicode_minus": False,
    "axes.formatter.use_mathtext": True,   # cmr10 wants mathtext for tick labels
})
XM, YM = 1.5, 1.0                     # the impulse location

FIG_W, FIG_H, DPI = 3.30, 3.03335, 300        # style_panel.m
AX_POS = [0.165, 0.205, 0.610, 0.665]  # room for the x and y labels;
CB_POS = [0.815, 0.205, 0.038, 0.665]  # see style_panel_delta.m
FS = 16

_par = os.path.join(HERE, "parula.npy")   # MATLAB parula, sampled
CMAP = (LinearSegmentedColormap.from_list("parula", np.load(_par))
        if os.path.exists(_par) else plt.get_cmap("viridis"))


def clock(directory, tfinal, backward):
    """{physical time: frame number} for every frame present in directory."""
    out = {}
    for f in sorted(glob.glob(os.path.join(directory, "fort.t[0-9]" * 1 + "*"))):
        m = re.search(r"fort\.t(\d{4})$", f)
        if not m:
            continue
        trun = float(open(f).readline().split()[0])
        # trun is the solver clock: tau for the adjoint (backward), t otherwise.
        # t_frame is the frame's physical time: T - tau for the adjoint, t otherwise.
        t_frame = (tfinal - trun) if backward else trun
        out[round(t_frame, 6) + 0.0] = int(m.group(1))
    return out


def nearest(cl, t):
    k = min(cl, key=lambda k: abs(k - t))
    return cl[k], k


def to_uniform(patches, N):
    h = LDOM / N
    Q = np.full((N, N), np.nan)
    for p in patches:
        rx, ry = int(round(p["dx"] / h)), int(round(p["dy"] / h))
        i0, j0 = int(round(p["xlow"] / h)), int(round(p["ylow"] / h))
        for j in range(p["my"]):
            for i in range(p["mx"]):
                Q[j0 + j * ry:j0 + (j + 1) * ry,
                  i0 + i * rx:i0 + (i + 1) * rx] = p["q"][j, i]
    return Q


def panel(patches, clim, label, t, out):
    N = int(round(LDOM / min(p["dx"] for p in patches)))
    Q = to_uniform(patches, N)

    fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=DPI, facecolor="w")
    ax = fig.add_axes(AX_POS)
    im = ax.imshow(Q, origin="lower", extent=[0, LDOM, 0, LDOM], cmap=CMAP,
                   vmin=clim[0], vmax=clim[1], interpolation="nearest")
    for p in patches:
        ax.add_patch(Rectangle((p["xlow"], p["ylow"]),
                               p["mx"] * p["dx"], p["my"] * p["dy"],
                               fill=False, edgecolor="k", linewidth=0.5))
    ax.plot(XM, YM, "ko", ms=5.3, mfc="k", mec="k")   # 22 px, matching the paper panels
    ax.set_xlim(0, LDOM); ax.set_ylim(0, LDOM); ax.set_aspect("equal")
    ax.set_xticks([0, 1, 2]); ax.set_yticks([0, 1, 2])
    ax.set_xlabel("x", fontsize=FS); ax.set_ylabel("y", fontsize=FS)
    ax.set_title(r"%s at $t = %.2f$" % (label, t), fontsize=FS)
    ax.tick_params(labelsize=FS)

    cax = fig.add_axes(CB_POS)
    cb = fig.colorbar(im, cax=cax)
    cb.ax.tick_params(labelsize=FS)

    fig.savefig(os.path.join(HERE, "panels", out), dpi=DPI, facecolor="w")
    plt.close(fig)
    print("    %-20s  %dx%d  %d patches" % (out, N, N, len(patches)))
