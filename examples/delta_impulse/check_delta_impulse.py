#!/usr/bin/env python3
"""Check the event-source frame and the adjoint-to-representer handoff."""

import glob
import configparser
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from fclaw_error import read_fort_q


def read_problem():
    config = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    config.read(os.path.join(HERE, "adjoint_options.ini"))
    user = config["adjoint-user"]
    run = config["adjoint"]

    def first(name):
        return float(user.get(name).split()[0])

    forward_config = configparser.ConfigParser(
        inline_comment_prefixes=("#", ";")
    )
    forward_config.read(os.path.join(HERE, "forward_options.ini"))

    return {
        "eps2": first("epsilon-2d"),
        "tfinal": float(run.get("tfinal").split()[0]),
        "xm": first("xm"),
        "ym": first("ym"),
        "tm": first("tm"),
        "pseudo_1d_experiment": user.get("pseudo-1d-experiment").split()[0].upper().startswith("T"),
        "wi": float(
            forward_config["forward-user"].get("W_i").split()[0]
        ),
    }


def frame_at(directory, t, tfinal, backward):
    candidates = []
    for qfile in glob.glob(os.path.join(directory, "fort.q[0-9][0-9][0-9][0-9]")):
        suffix = qfile[-4:]
        tfile = os.path.join(directory, "fort.t" + suffix)
        if not os.path.exists(tfile):
            continue
        with open(tfile) as stream:
            trun = float(stream.readline().split()[0])
        # trun is the solver clock: tau for the adjoint (backward), t otherwise.
        # t_frame is the frame's physical time: T - tau for the adjoint, t otherwise.
        t_frame = tfinal-trun if backward else trun
        candidates.append((abs(t_frame-t), t_frame, qfile))
    if not candidates:
        raise RuntimeError("no output frames in " + directory)
    return min(candidates)[1:]


def stats(qfile):
    mass = mx = my = 0.0
    qmax = -np.inf
    xmax = ymax = np.nan
    levels = set()

    for patch in read_fort_q(qfile):
        x = patch["xlow"] + (np.arange(patch["mx"]) + 0.5)*patch["dx"]
        y = patch["ylow"] + (np.arange(patch["my"]) + 0.5)*patch["dy"]
        X, Y = np.meshgrid(x, y)
        Q = patch["q"]
        area = patch["dx"]*patch["dy"]
        mass += Q.sum()*area
        mx += (Q*X).sum()*area
        my += (Q*Y).sum()*area
        levels.add(patch["level"])
        if Q.max() > qmax:
            j, i = np.unravel_index(np.argmax(Q), Q.shape)
            qmax = Q[j, i]
            xmax, ymax = X[j, i], Y[j, i]

    return {
        "mass": mass,
        "cx": mx/mass,
        "cy": my/mass,
        "qmax": qmax,
        "xmax": xmax,
        "ymax": ymax,
        "levels": sorted(levels),
    }


problem = read_problem()
time, qfile = frame_at(
    os.path.join(HERE, "adjoint0"), problem["tm"],
    problem["tfinal"], backward=True,
)
result = stats(qfile)
expected_peak = 1.0/(4*np.pi*problem["eps2"])
centroid_error = np.hypot(
    result["cx"]-problem["xm"], result["cy"]-problem["ym"]
)
peak_error = np.hypot(
    result["xmax"]-problem["xm"], result["ymax"]-problem["ym"]
)

print("instantaneous adjoint event")
print(f"  requested physical time : {problem['tm']:.6f}")
print(f"  frame physical time     : {time:.6f}")
print(f"  source location         : ({problem['xm']:.6f}, {problem['ym']:.6f})")
print(f"  centroid                : ({result['cx']:.6f}, {result['cy']:.6f})")
print(f"  centroid error          : {centroid_error:.3e}")
print(f"  maximum location        : ({result['xmax']:.6f}, {result['ymax']:.6f})")
print(f"  maximum-location error  : {peak_error:.3e}")
print(f"  mass                    : {result['mass']:.8f}  (expected 1)")
print(f"  maximum                 : {result['qmax']:.8f}")
print(f"  point-kernel maximum    : {expected_peak:.8f}")
print(f"  AMR levels present      : {result['levels']}")

if abs(result["mass"]-1.0) > 2.0e-2:
    raise SystemExit("FAIL: event mass is not approximately one")
if centroid_error > 2.0e-2:
    raise SystemExit("FAIL: event is not centered on the source")
if result["qmax"] < 0.9*expected_peak:
    raise SystemExit("FAIL: event peak is unexpectedly low")

alpha0_time, alpha0_file = frame_at(
    os.path.join(HERE, "adjoint0"), 0.0,
    problem["tfinal"], backward=True,
)
r0_time, r0_file = frame_at(
    os.path.join(HERE, "forward0"), 0.0,
    problem["tfinal"], backward=False,
)
rtm_time, rtm_file = frame_at(
    os.path.join(HERE, "forward0"), problem["tm"],
    problem["tfinal"], backward=False,
)
alpha0 = stats(alpha0_file)
r0 = stats(r0_file)
rtm = stats(rtm_file)

handoff_mass_error = abs(problem["wi"]*r0["mass"]-alpha0["mass"])
handoff_centroid_error = np.hypot(
    r0["cx"]-alpha0["cx"], r0["cy"]-alpha0["cy"]
)
refocus_error = np.hypot(
    rtm["xmax"]-problem["xm"], rtm["ymax"]-problem["ym"]
)

print("\nadjoint-to-representer handoff")
print(f"  W_i                     : {problem['wi']:.6f}")
print(f"  alpha mass at t=0       : {alpha0['mass']:.8f}")
print(f"  W_i r mass at t=0       : {problem['wi']*r0['mass']:.8f}")
print(f"  mass difference         : {handoff_mass_error:.3e}")
print(f"  centroid difference     : {handoff_centroid_error:.3e}")
print(f"  adjoint AMR levels      : {alpha0['levels']}")
print(f"  representer AMR levels  : {r0['levels']}")

print("\nforward refocusing")
print(f"  frame physical time     : {rtm_time:.6f}")
print(f"  maximum location        : ({rtm['xmax']:.6f}, {rtm['ymax']:.6f})")
print(f"  maximum-location error  : {refocus_error:.3e}")
print(f"  maximum                 : {rtm['qmax']:.8f}")

if abs(alpha0_time) > 1.0e-10 or abs(r0_time) > 1.0e-10:
    raise SystemExit("FAIL: t=0 handoff frames were not found")
if abs(rtm_time-problem["tm"]) > 1.0e-10:
    raise SystemExit("FAIL: representer output is not at the observation time")
if handoff_mass_error > 2.0e-2:
    raise SystemExit("FAIL: mass changed during the adjoint-to-forward handoff")
if handoff_centroid_error > 2.0e-2:
    raise SystemExit("FAIL: handoff displaced the field")
if refocus_error > 4.0e-2:
    raise SystemExit("FAIL: representer did not refocus at the gauge")

print("\noverall result            : PASS")
