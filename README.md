# Smoke2d-4DVar

Representer-Based Data Assimilation on Adaptively Refined Meshes for Idealized Smoke Transport

Smoke2d-4DVar is a computational model for 2d idealized wildfire smoke transport.  A key feature of the model is that it incorporates 4DVar data assimilation using data using the representer method (Bennett, 2002).   Smoke2d-4DVar is built on [ForestClaw](https://github.com/ForestClaw/forestclaw),
a parallel, patch-based adaptive mesh refinement (AMR) framework that uses
[p4est](https://www.p4est.org) for the mesh and the wave propagation algorithm finite-volume solvers (Leveque, 2002) . A smoke concentration is advected by a prescribed, velocity field,  A small number 
of point observations of the concentration are assimilated by the representer method
(Bennett, 2002), which gives the minimizer (optimal estimate) of the weak-constraint 
cost functional by combining one prior trajectory run with one adjoint and one 
representer run per observation.

Every run in the method (the prior, each adjoint, each representer) has its own 
adaptive mesh. The meshes are coupled through ForestClaw's overlap exchange mechanism, 
which interpolates a field from one mesh onto the cell centres of another.

Smoke2d-4DVar is the 2d version of a larger 3d code and is described in the paper Azike et al., (2026).

---

## Contents

1. [The smoke transport problem](#1-the-smoke-transport-problem)
2. [Repository layout](#2-repository-layout)
3. [Requirements](#3-requirements)
4. [Building Smoke2d-4DVar](#4-building-smoke2d-4dvar)
5. [Running the assimilation (`adjoint_forward`)](#5-running-the-assimilation-adjoint_forward)
6. [Configuration reference](#6-configuration-reference)
7. [Output](#7-output)
8. [Verification: the delta-impulse test (`delta_impulse`)](#8-verification-the-delta-impulse-test-delta_impulse)
9. [Plotting](#9-plotting)
10. [Uniform-mesh reference runs](#10-uniform-mesh-reference-runs)
11. [Troubleshooting](#11-troubleshooting)
12. [Known limitations](#12-known-limitations)
13. [References](#13-references)
14. [License](#14-license)

---

## 1. An idealized smoke transport model


The evolution of smoke concentration $q(x,y,t)$ in an idealized setting can be described by the transport model
$$
\frac{\partial q}{\partial t} + \nabla\cdot(\mathbf{u} q) = 0,
\qquad q(\mathbf{x},0) = q_0(\mathbf{x}),
$$

where $\mathbf u$ is a steady, velocity $\mathbf{u}(x,y) = (u(x,y),v(x,y))$. 

The code illustrates from the paper cited above. 

## 2. Repository layout

The files in this repository are laid out as follows. 

```
Smoke2d-4DVar/
├── CMakeLists.txt              top level; finds ForestClaw
├── src/                        library smoke2d-4DVar
│   ├── smoke2d_*.{c,h,f90}     options, gauge helpers, Fortran module
│   ├── 2d/                     2D transport: prescribe velocity field and store in auxiliary arrays
│   └── fortran_source/         clawpatch46_* tagging/interpolation routines used by the AMR routines
└── examples/
    ├── adjoint_forward/        the assimilation (Sections 3 and 8)
    │   ├── adjoint_forward.cpp     driver: prior, adjoints, representers, β, optimal estimate
    │   ├── user_run.c              one-step time stepper so fields can be exchanged between steps
    │   ├── af_overlap_patch.c      mesh-to-mesh interpolation callbacks
    │   ├── af_common_*.f90         shared initial condition, aux (velocity) and overlap routines
    │   ├── adjoint/                adjoint solver: options, source (heat kernels)
    │   ├── forward/                prior ("model") and representer solvers, options, sources;
    │   │                           write_gauges.py sets the observations (Section 8)
    │   ├── *_options.ini           run configuration (Section 9)
    │   └── *.m                     MATLAB plotting (Section 12)
    └── delta_impulse/          verification with an instantaneous impulse (Section 11)
```

An executable for each example (adjoint_forward and delta_impulse) is built. 

## 3. Requirements

| Requirement | Notes |
|---|---|
| C, C++17 and Fortran compilers | The same compilers used to build ForestClaw (e.g. GCC/gfortran). |
| CMake ≥ 3.19 | ForestClaw's minimum. |
| MPI | ForestClaw is built with MPI. The runs reported here used one MPI rank. |
| Python 3 with NumPy and Matplotlib | `forward/write_gauges.py` and the `delta_impulse` scripts. `write_gauges.py` also needs `fclaw_analysis.py` from ForestClaw's `python/` directory. |
| MATLAB | Plotting with the Clawpack/ForestClaw MATLAB graphics routines. |

In what follows, we describe the process that automatically downloads and builds needed libraries, including ForestClaw and p4est.  A more advanced user may wish to build these libraries separately.  

## 4. Building Smoke2d-4DVar

We recommend building Smoke2d-4DVar out of source, in a directory parallel to the repository. 

```sh
git clone https://github.com/PatriciaAzike/Smoke2d-4DVar.git
mkdir smoke2d-4DVar-build
```

Create an executable configuration file with the following commands. 

```sh
# File : config-Smoke2d-4DVar.sh

cmake \
    -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_CXX_COMPILER=g++ \
    -DCMAKE_Fortran_COMPILER=gfortran \
    -DCMAKE_C_FLAGS="-O2 -g -Wall" \
    -DCMAKE_CXX_FLAGS="-O2 -g -Wall" \
    -DCMAKE_Fortran_FLAGS="-O2 -g -Wall -cpp -Wno-unused-dummy-argument" \
    -DCMAKE_EXE_LINKER_FLAGS="-llapack" \
    -DCMake_External_build=True \
    ../Smoke2d-4DVar
```
Then, configure, build and install the executables `adjoint_forward` and `delta_impulse`. 

```sh
cd smoke2d-4DVar-build
../config-Smoke2d-4DVar
make -j4 install
```
This builds:

- `<build-directory>/examples/adjoint_forward/adjoint_forward`
- `<build-directory>/examples/delta_impulse/delta_impulse`

## 5. Running the assimilation (`adjoint_forward`)

All inputs are read from, and all outputs written below,
`examples/adjoint_forward/`.

```sh
cd Smoke2d-4DVar/examples/adjoint_forward

# 1. Link the executable here.
ln -sf ../../../<build-directory>/examples/adjoint_forward/adjoint_forward

# 2. Set the observations. The list at the top of forward/write_gauges.py
#    (xm, ym, tm, dm) is the single source: the script writes them into
#    adjoint_options.ini and writes gauges.data into forward/ and model/
#    (needs ForestClaw's python/ directory for fclaw_analysis.py).
cd forward
PYTHONPATH=/path/to/forestclaw/python python3 write_gauges.py
cd ..

# 3. Run 
./adjoint_forward
```

Re-run step 2 whenever the observations change; the driver stops with a
message if the number of gauges does not match `mdata`. The data weights
`w_eps` are not in the script's list: set them directly in
`adjoint_options.ini`.

Do not move or delete anything in `model/`, `adjoint{j}/` or `forward{j}/`
while the program runs. The final optimal estimate step reloads every checkpoint
written earlier in the run.

**What the run prints at the end** (values depend on the configuration):

```
Calling get_gauges for model
get_gauges: q_out[0] = ...           prior q_F at each observation
Calling get_gauges for representer j = 0
get_gauges: q_out[0] = ...           column j of R
...
Representer matrix: max|R_ij - R_ji| = ..., relative asymmetry = ...
frobenius norm: Representer matrix: ||R - R^T||_F / ||R||_F = ...
beta[0] = ...                        the coefficients β
...
```

## 6. Configuration reference

The driver reads three option files, one per kind of run. Each configuration file has a
sections : (`[*-user]`), a mesh section (`[*-clawpatch]`), a run section
(`[adjoint]`, `[forward]`, `[model]`) and a solver section (`[*-clawpack46]`).
The run and solver sections take the standard ForestClaw and Clawpack 4.6
options; the table lists the ones that matter here.

| File | Governs | Main options |
|---|---|---|
| `adjoint_options.ini` | the adjoint runs, **and the observations** | see below |
| `forward_options.ini` | the representer runs | `W_f`, `W_i` in `[forward-user]` |
| `model_options.ini` | the prior run | `initial-condition`, and the same `W_f`, `W_i` fields in `[model-user]` |

**Observations and experiment** (`[adjoint-user]` in `adjoint_options.ini`).
`mdata`, `xm`, `ym`, `tm` and `dm` are written by `forward/write_gauges.py`
from its observation list; edit them there, not here, or the gauge files and
the options will disagree. The script rewrites the whole file with Python's
`configparser`, which drops full-line comments.

| Option | Meaning |
|---|---|
| `pseudo-1d` | `T` = pseudo-1D experiment, `F` = full 2D (Section 1). Also selects the velocity field for all runs. |
| `mdata` | number of observations $M$ |
| `xm`, `ym`, `tm` | observation positions and times (`mdata` values each) |
| `dm` | observed values $d_m$ |
| `w_eps` | data weights $W_{\varepsilon,m}$ |
| `epsilon-1d`, `epsilon-2d` | heat-kernel width $\varepsilon$ for the pseudo-1D / 2D source |
| `x0`, `y0` | centre of the prior disk (full 2D, `initial-condition = 2`) |

**Weights** (`[forward-user]` and `[model-user]`): `W_f` (model error),
`W_i` (initial-condition error). Keep the two sections identical: the
prior and representer runs share one Fortran module
(`forward/forward_module.f90`), so whichever file was read last supplies the
values for both.

**prior initial condition** (`initial-condition` in `[model-user]`):
`0` zero, `1` Gaussian $e^{-\beta\lvert\mathbf{x}-\mathbf{x}_0\rvert^2}$
(`beta`, `x0`, `y0`), `2` strip (pseudo-1D) or disk (2D) as in Section 1.

**Mesh and time** (run and clawpatch sections of each file):

| Option | Configured value | Meaning |
|---|---|---|
| `mx`, `my` | 32 | cells per patch in each direction |
| `minlevel`, `maxlevel` | 2–4 (prior, representers), 2–5 (adjoints) | coarsest and finest refinement level; level $\ell$ has $2^\ell \times 32$ cells across the domain |
| `refine_threshold`, `coarsen_threshold` | 0.125, 0.06 | max − min of the field above which a patch is refined / below which it is coarsened (Section 3) |
| `tfinal`, `nout` | 2.0, 20 | final time and number of output intervals (output every 0.1) |
| `ax`, `bx`, `ay`, `by` | 0, 2, 0, 2 | domain |
| `periodic_x`, `periodic_y` | True | periodic boundaries |

## 7. Output

All frames are ForestClaw ASCII output (`fort.qNNNN` data, `fort.tNNNN`
header with the time and number of patches), `NNNN = 0000 … nout`.

| Location | Contents | Clock |
|---|---|---|
| `examples/adjoint_forward/fort.*` | the optimal estimate $\hat q$, on the prior mesh | physical time |
| `model/` | prior $q_F$, its gauge files, checkpoints | physical time |
| `adjoint{j}/` | adjoint $\alpha_{j+1}$, checkpoints | **reverse** time |
| `forward{j}/` | representer $r_{j+1}$, checkpoints | physical time |
| `forward/` | representer gauge files | physical time |

**Adjoint frames run on a reversed clock.** Adjoint frame $N$ is at run time
$N\,T/\texttt{nout}$, i.e. physical time $T - N\,T/\texttt{nout}$. The same
physical instant is adjoint frame $N$ and forward frame $\texttt{nout}-N$.

Checkpoints are `fort_frame_NNNN.checkpoint` and
`fort_frame_NNNN.partition`. They are needed during the run and can be
deleted afterwards.

## 8. Verification: the delta-impulse test (`delta_impulse`)

`examples/delta_impulse/` checks the adjoint–representer machinery on one
observation, with the observation applied as an **instantaneous** impulse:
the spatial heat kernel is added to the adjoint once, when the backward
sweep reaches $\tau_m = T - t_m$, instead of being spread over time by the
temporal kernel. With nothing smeared in time, the event frame can be
compared directly with the heat kernel.

It reuses the `adjoint_forward` sources, with two differences:

- `examples/delta_impulse/src/adjoint_src2_instantaneous.f90` replaces the
  adjoint source routine;
- the driver is compiled with `INSTANTANEOUS_VERIFICATION_ONLY`, so it stops
  after the representer and does not solve for $\beta$.

```sh
ln -sf <build-directory>/examples/delta_impulse/delta_impulse .
python3 write_gauges.py                         # writes gauges.data, sets up run directories
./delta_impulse                                 # runs
python3 check_delta_impulse.py                  # checks the event frame and the handoff
python3 plot_delta_impulse.py --t 0.0 0.8 1.5   # panels in panels/
```

`write_gauges.py` needs `fclaw_analysis.py` from ForestClaw's `python/` 
directory on `PYTHONPATH`.

For the configured case (one observation at $(1.5, 1.0)$, $t_m = 1.5$,
$\varepsilon = 0.01$), the event frame should have mass ≈ 1, centroid
≈ $(1.5, 1.0)$ and peak ≈ $1/(4\pi\varepsilon) = 7.96$, up to cell-centre
sampling and AMR transfer. `check_delta_impulse.py` reports these and also
checks the adjoint-to-representer handoff at $t=0$ and the location of the
representer peak at $t = t_m$. See
[`examples/delta_impulse/README.md`](examples/delta_impulse/README.md) for
details.

## 9. Plotting

**MATLAB** (`examples/adjoint_forward/*.m`). `setplot2.m` and
`afterframe.m` configure the Clawpack/ForestClaw MATLAB graphics (e.g.
`plotclaw2`, `showpatchborders`), which must be on the MATLAB path. Run
`plotclaw2` in the directory holding the frames to plot. `afterframe.m`
draws the patch borders and, for the 2D case, the exact filament of the
prior disk (`filament_soln.m`).

**Python** (`examples/delta_impulse/`). `plot_delta_impulse.py` draws the
adjoint and representer at chosen physical times with one shared colour
scale; `plot_delta.py` holds the reading and panel routines it uses.

## 10. Uniform-mesh reference runs

A uniformly refined reference is obtained by setting `minlevel = maxlevel`
in each option file (level 4 for the prior and representers, level 5
for the adjoints, i.e. $512^2$ and $1024^2$ cells). For the configuration in
this repository (full 2D, three observations, one MPI rank), the adaptive
run took $1293 \pm 15$ s and the uniform reference $23{,}667 \pm 619$ s
(mean ± sample standard deviation of three runs each), a speedup of about
18.3. Almost all of the saving is in the adjoint solves.

## 11. Troubleshooting

| Symptom | Cause and fix |
|---|---|
| `find_package(FORESTCLAW)` fails | Pass `-DFORESTCLAW_ROOT=<install prefix>`, and make sure ForestClaw was built with `-Dclawpack=ON`. |
| Undefined `dgesv_` at link time | No LAPACK was linked through ForestClaw's dependencies. Add `find_package(LAPACK REQUIRED)` and link `LAPACK::LAPACK` to `adjoint_forward` and `delta_impulse`. |
| `get_gauges: FATAL -- accumulator has N gauges but mdata = M` | `gauges.data` is out of date. Re-run `forward/write_gauges.py` (from inside `forward/`); it rewrites both `forward/gauges.data` and `model/gauges.data`. |
| The prior (`model/`) is zero everywhere and never refines | The initial condition is wrong or missing, e.g. `adjoint/adjoint_fdisc.f90` (the disk for `pseudo-1d = 2`) not compiled in. |
| `File does not exist` when "Restarting model from checkpoint file" | A checkpoint in `model/` (or `adjoint{j}/`, `forward{j}/`) was deleted or moved during the run. Re-run without touching those directories. |
| Stale results after changing options | Output directories are overwritten, not cleared; delete old `fort.*` files before a run if frame counts change. |

## 12. Known limitations

- Tested on one MPI rank. The overlap exchange is written for distributed
  meshes, but multi-rank runs have not been validated.
- The forward and adjoint discretizations are not exact discrete adjoints,
  so $R$ is symmetrized before the solve (Section 3).
- Only output style 1 (equally spaced output times) is supported by the
  one-step time stepper.
- The velocity field is analytic and steady.

## 13. References

- P. O. Azike et al., "Representer-Based Data Assimilation on Adaptively 
  Refined Meshes for Idealized Wildfire Smoke Transport", 2026.
- A. F. Bennett, *Inverse Methods in Physical Oceanography*, Cambridge
  University Press, 1992.
- A. F. Bennett, *Inverse Modeling of the Ocean and Atmosphere*, Cambridge
  University Press, 2002.
- D. Calhoun and C. Burstedde, "ForestClaw: A parallel algorithm for
  patch-based adaptive mesh refinement on a forest of quadtrees,"
  arXiv:1703.03116, 2017.
- C. Burstedde, L. C. Wilcox and O. Ghattas, "p4est: Scalable algorithms for
  parallel adaptive mesh refinement on forests of octrees," *SIAM J. Sci.
  Comput.* 33(3):1103–1133, 2011.
- R. J. LeVeque, *Finite Volume Methods for Hyperbolic Problems*, Cambridge
  University Press, 2002.

## 14. License

BSD 3-Clause; see [LICENSE](LICENSE). ForestClaw is distributed under its
own BSD-style license, and p4est under the GPL (version 2 or later).
