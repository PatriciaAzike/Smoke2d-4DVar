#!/usr/bin/env python3
"""
Regenerate the paper panels:
  * pseudo-1D, W_eps = 100 and 0.01 
  * full 2D, W_eps = 100 and 0.01,

    python3 paper_figs/regen.py   # writes into paper_figs/out/

Inputs are the frames copied into this folder.

WHY: style_panel.m leaves a bottom margin of 0.132*910 = 120 px. At FontSize 16 on
a 300 dpi canvas one text line is ~67 px, and the bottom needs tick labels plus an
xlabel, about 175 px. So the xlabel fell off the canvas and the ylabel was clipped.

New rectangle, measured so nothing clips with xlabel, ylabel, title and a
two-digit colorbar label all present (clearances L 25, R 36, T 43, B 23 px):

    field panels   axes [0.165 0.205 0.610 0.665]   cb [0.815 0.205 0.038 0.665]
    line panels    axes [0.165 0.205 0.790 0.665]   (no colorbar, so wider)

Both share the left edge, the bottom and the height, so a line panel above a field
panel still aligns. The cost is the plot box: 604 x 605 px against 683 x 700, so
the data area prints ~12% smaller at the same \\includegraphics width.

clim and tick values were measured off the existing PNGs by locating the colorbar
tick labels and solving the linear map, then cross-checked against afterframe.m:

    2D / uniform / prior model   clim [-0.5, 4.5]   ticks 0:4    (afterframe.m agrees)
    pseudo-1D field              clim [-0.5, 3.06]  ticks 0:3    (measured -0.49, 3.11)
"""
import os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = HERE
sys.path.insert(0, HERE)
from fclaw_error import read_fort_q, forestclaw_gauge    # noqa: E402

LDOM = 2.0

# Computer Modern, to match the paper's body text.
# use matplotlib's bundled Computer Modern instead -- same glyphs, no
# LaTeX dependency, and much faster.
plt.rcParams.update({
    "mathtext.fontset": "cm",
    "font.family":      "serif",
    "font.serif":       ["cmr10", "DejaVu Serif"],
    "axes.unicode_minus": False,
    "axes.formatter.use_mathtext": True,   # cmr10 wants mathtext for tick labels
})
FIG_W, FIG_H, DPI = 3.30, 3.03335, 300
AX_FIELD = [0.165, 0.205, 0.610, 0.665]
CB_POS   = [0.815, 0.205, 0.038, 0.665]
FS = 16
CMAP = LinearSegmentedColormap.from_list(
    "parula", np.load(os.path.join(HERE, "parula.npy")))

G2D = [(0.67, 1.0), (0.80, 1.2), (0.55, 1.2)]      # the three gauges
G1D = [(0.55, 1.0), (0.67, 1.0), (0.80, 1.0)]      # pseudo-1D: all on y = 1

#   output name                    source                              clim        ticks     title            gauges
JOBS = [
 ("1d_amr_0004_t00.40_hw",   "fort.q0004_amr_1d_hw",   (-0.5,3.06),range(4), r"$\hat{q}$ at t = 0.40", G1D),
 ("1d_amr_0004_t00.40_lw",   "fort.q0004_amr_1d_lw",   (-0.5,3.06),range(4), r"$\hat{q}$ at t = 0.40", G1D),
 ("2d_amr_0004_t00.40_hw",   "fort.q0004_amr_hw",      (-0.5,4.5), range(5), r"$\hat{q}$ at t = 0.40", G2D),
 ("2d_amr_0010_t01.00_hw",   "fort.q0010_amr_hw",      (-0.5,4.5), range(5), r"$\hat{q}$ at t = 1.00", G2D),
 ("2d_amr_0004_t00.40_lw",   "fort.q0004_amr_lw",      (-0.5,4.5), range(5), r"$\hat{q}$ at t = 0.40", G2D),
 ("2d_amr_0010_t01.00_lw",   "fort.q0010_amr_lw",      (-0.5,4.5), range(5), r"$\hat{q}$ at t = 1.00", G2D),
 ("uniform_0004_t00.40",     "fort.q0004_uniform",     (-0.5,4.5), range(5), r"$\hat{q}$ at t = 0.40", G2D),
 ("uniform_0010_t01.00",     "fort.q0010_uniform",     (-0.5,4.5), range(5), r"$\hat{q}$ at t = 1.00", G2D),
 ("model_amr_0004_t00.40_hw","fort.q0004_model_2d_hw", (-0.5,4.5), range(5), r"$q_F$ at t = 0.40",     G2D),
 ("model_amr_0010_t01.00_hw","fort.q0010_model_2d_hw", (-0.5,4.5), range(5), r"$q_F$ at t = 1.00",     G2D),
]


def to_uniform(patches, N):
    h = LDOM / N
    Q = np.full((N, N), np.nan)
    for p in patches:
        rx, ry = int(round(p["dx"]/h)), int(round(p["dy"]/h))
        i0, j0 = int(round(p["xlow"]/h)), int(round(p["ylow"]/h))
        for j in range(p["my"]):
            for i in range(p["mx"]):
                Q[j0+j*ry:j0+(j+1)*ry, i0+i*rx:i0+(i+1)*rx] = p["q"][j, i]
    return Q


def sample_plotted_cell(Q, x, y):
    """Return the cell value displayed by the nearest-neighbour image at (x,y)."""
    ny, nx = Q.shape
    i = min(max(int(np.floor(x * nx / LDOM)), 0), nx - 1)
    j = min(max(int(np.floor(y * ny / LDOM)), 0), ny - 1)
    return Q[j, i]


def field_panel(name, src, clim, ticks, title, gauges):
    ps = read_fort_q(os.path.join(ROOT, src))
    N = int(round(LDOM / min(p["dx"] for p in ps)))
    Q = to_uniform(ps, N)

    fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=DPI, facecolor="w")
    ax = fig.add_axes(AX_FIELD)
    im = ax.imshow(Q, origin="lower", extent=[0, LDOM, 0, LDOM], cmap=CMAP,
                   vmin=clim[0], vmax=clim[1], interpolation="nearest")
    for p in ps:
        ax.add_patch(Rectangle((p["xlow"], p["ylow"]),
                               p["mx"]*p["dx"], p["my"]*p["dy"],
                               fill=False, edgecolor="k", linewidth=0.5))
    # Color each gauge by its interpolated concentration, using the same map
    # and limits as the field.  The black edge marks its location.
    norm = Normalize(vmin=clim[0], vmax=clim[1], clip=True)
    # ForestClaw's gauge interpolation, so the numbers printed below match the
    # gauge values reported by adjoint_forward.
    gauge_values = [forestclaw_gauge(ps, gx, gy) for gx, gy in gauges]
    for (gx, gy), value in zip(gauges, gauge_values):
        ax.plot(gx, gy, marker="o", ms=5.3,
                mfc=CMAP(norm(value)), mec="black", mew=1.0,
                linestyle="none", zorder=5)

    ax.set_xlim(0, LDOM); ax.set_ylim(0, LDOM); ax.set_aspect("equal")
    ax.set_xticks([0, 1, 2]); ax.set_yticks([0, 1, 2])
    ax.set_xlabel("x", fontsize=FS); ax.set_ylabel("y", fontsize=FS)
    ax.set_title(title, fontsize=FS)
    ax.tick_params(labelsize=FS)

    cax = fig.add_axes(CB_POS)
    cb = fig.colorbar(im, cax=cax, ticks=list(ticks))
    cb.ax.tick_params(labelsize=FS)

    out = os.path.join(HERE, "out", name + ".png")
    fig.savefig(out, dpi=DPI, facecolor="w"); plt.close(fig)
    print("  %-28s %4d patches  %dx%d  qmax %.4f  gauges %s"
          % (name, len(ps), N, N, max(p["q"].max() for p in ps),
             ", ".join("%.4f" % v for v in gauge_values)))


if __name__ == "__main__":
    os.makedirs(os.path.join(HERE, "out"), exist_ok=True)
    for job in JOBS:
        field_panel(*job)
