#!/usr/bin/env python3
"""Plot physical-time panels from the instantaneous delta-impulse run."""

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import plot_delta as plotting

parser = argparse.ArgumentParser()
parser.add_argument("--t", type=float, nargs="+", default=[0.0, 0.8, 1.5])
parser.add_argument("--tfinal", type=float, default=2.0)
args = parser.parse_args()

os.makedirs(os.path.join(HERE, "panels"), exist_ok=True)


def read_row(directory, backward):
    clock = plotting.clock(directory, args.tfinal, backward)
    picked = []
    for requested_time in args.t:
        frame, t_frame = plotting.nearest(clock, requested_time)
        patches = plotting.read_fort_q(
            os.path.join(directory, "fort.q%04d" % frame)
        )
        picked.append((frame, t_frame, patches))
        print(
            "  physical t=%.2f -> frame %02d (t=%.2f)"
            % (requested_time, frame, t_frame)
        )
    return picked


print("alpha_1")
alpha = read_row(os.path.join(HERE, "adjoint0"), backward=True)
print("r_1")
representer = read_row(os.path.join(HERE, "forward0"), backward=False)

# Use one range for the complete two-row figure, so equal colours denote equal
# values in alpha_m and r_m as well as across time within each row.
all_panels = alpha + representer
vmin = min(
    min(patch["q"].min() for patch in patches)
    for _, _, patches in all_panels
)
vmax = max(
    max(patch["q"].max() for patch in patches)
    for _, _, patches in all_panels
)
if vmin > -0.01*vmax:
    vmin = -0.083*vmax
clim = (vmin, vmax)
print("shared two-row scale [%.4g, %.4g]" % clim)

for _, t_frame, patches in alpha:
    plotting.panel(
        patches, clim, r"$\alpha_1$", t_frame,
        "alpha_t%05.2f.png" % t_frame,
    )
for _, t_frame, patches in representer:
    plotting.panel(
        patches, clim, r"$r_1$", t_frame,
        "rm_t%05.2f.png" % t_frame,
    )
