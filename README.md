# Smoke2d-4DVar

Representer-Based Data Assimilation on Adaptively Refined Meshes for Idealized Smoke Transport

Smoke2d_4DVar is a computational model for 2d idealized wildfire smoke transport.  A key feature of the model is that it incorporates 4DVar data assimilation using data using the representer method (Bennett, 2005).   Smoke2d-4DVar is built on [ForestClaw](https://github.com/ForestClaw/forestclaw),
a parallel, patch-based adaptive mesh refinement (AMR) framework that uses
[p4est](https://www.p4est.org) for the mesh and the wave propagation algorithm finite-volume solvers (Leveque, 2002) . A smoke concentration is advected by a prescribed, velocity field,  A small number 
of point observations of the concentration are assimilated by the representer method
(Bennett, 2005), which gives the minimizer (optimal estimate) of the weak-constraint 
cost functional by combining one prior trajectory run with one adjoint and one 
representer run per observation.

Every run in the method (the prior, each adjoint, each representer) has its own 
adaptive mesh. The meshes are coupled through ForestClaw's overlap exchange mechanism, 
which interpolates a field from one mesh onto the cell centres of another.

Smoke2d_4DVar is the 2d version of a larger 3d code and is described in the paper <fill in paper name here>

---

## Contents

1. [The smoke transport problem](#1-the-smoke-transport-problem)
2. [Weak-constraint 4D-Var and the representer method](#2-Weak-constraint-4D-Var-and-the-representer-method)
3. [How the code implements the method](#3-how-the-code-implements-the-method)
4. [Repository layout](#4-repository-layout)
5. [Requirements](#5-requirements)
6. [Building ForestClaw](#6-building-forestclaw)
7. [Building Smoke2d-4DVar](#7-building-smoke2d-4dvar)
8. [Running the assimilation (`adjoint_forward`)](#8-running-the-assimilation-adjoint_forward)
9. [Configuration reference](#9-configuration-reference)
10. [Output](#10-output)
11. [Verification: the delta-impulse test (`delta_impulse`)](#11-verification-the-delta-impulse-test-delta_impulse)
12. [Plotting](#12-plotting)
13. [Uniform-mesh reference runs](#13-uniform-mesh-reference-runs)
14. [Troubleshooting](#14-troubleshooting)
15. [Known limitations](#15-known-limitations)
16. [References](#16-references)
17. [License](#17-license)

---

## 1. An idealized smoke transport model


The evolution of smoke concentration $q(x,y,t)$ in an idealized setting can be described by the transport model
$$
\frac{\partial q}{\partial t} + \nabla\cdot(\mathbf{u} q) = 0,
\qquad q(\mathbf{x},0) = q_0(\mathbf{x}),
$$

where $\mathbf u$ is a steady, velocity $\mathbf{u}(x,y) = (u(x,y),v(x,y))$. 

The code illustrates from the paper cited above. 

## 2. How the code implements the method

The driver is `examples/adjoint_forward/adjoint_forward.cpp`. It creates one
ForestClaw "global" (a mesh plus a solver) for the prior, and one
adjoint and one representer global for each observation, and then runs:

1. **Prior.** Run $q_F$ from $t=0$ to $T$ in `model/`, writing an output
   frame and a checkpoint (restart file) at each of the `nout + 1` output
   times.
2. **For each observation $m$** (`j = m-1` in `run_program`; `m` itself 
   is the number of observations, `mdata`):
   1. **Adjoint.** Run $\alpha_m$ once through the whole backward sweep in
      `adjoint{j}/`, checkpointing at every output time.
   2. **Representer initial condition.** Interpolate $\alpha_m$ at $t=0$ onto
      the representer mesh and set $r_m(\cdot,0) = W_i^{-1}\alpha_m(\cdot,0)$.
   3. **Representer, one step at a time.** The representer needs $\alpha_m$
      at every one of its own time steps, but $\alpha_m$ was computed
      backward. After each representer step, the adjoint is restarted from
      the nearest checkpoint on the far side of the current time, advanced
      (backward) to the representer's current time, and interpolated onto the
      representer mesh, where it enters as the source term $W_f^{-1}\alpha_m$. 
      The representer is checkpointed at every output time in `forward{j}/`.
3. **Coefficients.** Read $q_F$ and each $r_j$ at the observation points from
   the ForestClaw gauges, and form the representer matrix $R$. Calculate its asymmetry ($\lVert R - R^T\rVert_F / \lVert R\rVert_F$), symmetrize it, add
   $W_\varepsilon^{-1}$ to the diagonal and solve for $\beta$ with LAPACK
   `dgesv`.
4. **Optimal estimate $\hat{q}$.** For each output time, reload the prior $q_F$ and all representers from their checkpoints, interpolate each $r_m$ onto the
prior mesh, add $\beta_m r_m$ and write $\hat q$ as an output frame in the 
example directory.

**Mesh coupling.** The meshes differ in refinement pattern and level, so
fields are transferred point by point: every cell centre of the receiving
mesh is a query point, and `fclaw_overlap_exchange` finds the leaf patch of
the sending mesh that contains it and interpolates there
(`af_overlap_patch.c`, `af_common_overlap.f90`).

**Discretization.** Every solve uses Clawpack 4.6's second-order wave
propagation method with the monotonized-centred limiter.

- The adjoint uses a quasilinear (original WPA) Riemann solver, consistent with
  its advective-form equation.
- The prior and representers use a flux-form (f-wave) Riemann solver
  with cell-centred velocities, so they conserve mass: the mass of $q_F$ is
  constant to round-off throughout a run, on the adaptive mesh as well.
- Source terms are added by a fractional-step splittling: after each transport 
  step, a forward -Euler update $q \leftarrow q + \Delta t\, S$.

The forward and adjoint equations are discretized separately (flux form versus advective form, each on its own adaptive mesh, with fields transferred between meshes by interpolation), so the discrete adjoint is not the exact transpose of the discrete forward operator. $R$ is therefore slightly asymmetric; the asymmetry is reported and $R$ is replaced by its symmetric part before solving for $\beta$.

**Adaptivity.** Each mesh refines where its own field varies. A patch is
refined when the max − min of its field (prior ($q_F$), adjoint ($\alpha_m$) or representer ($r_m$)) exceeds `refine_threshold`, and a family of
four patches is coarsened when it is below `coarsen_threshold` on all four.
These tests are the routines in `src/fortran_source/`
(`clawpatch46_fort_tag4refinement`, `clawpatch46_fort_tag4coarsening`),
which every solver installs. The mesh is regridded after every time step.

## 4. Repository layout

```
Smoke2d-4DVar/
├── CMakeLists.txt              top level; finds ForestClaw
├── src/                        library smoke2d-4DVar
│   ├── smoke2d_*.{c,h,f90}     options, gauge helpers, Fortran module
│   ├── 2d/                     2D transport: velocity from stream function, setaux, Riemann solvers
│   └── fortran_source/         clawpatch46_* tagging/interpolation routines used by the solvers
└── examples/
    ├── adjoint_forward/        the assimilation (Sections 3 and 8)
    │   ├── adjoint_forward.cpp     driver: prior, adjoints, representers, β, optimal estimate
    │   ├── user_run.c              one-step time stepper so fields can be exchanged between steps
    │   ├── af_overlap_patch.c      mesh-to-mesh interpolation callbacks
    │   ├── af_common_*.f90         shared initial condition, aux (velocity) and overlap routines
    │   ├── adjoint/                adjoint solver: options, source (heat kernels), stream function −ψ
    │   ├── forward/                prior ("model") and representer solvers, options, sources;
    │   │                           write_gauges.py sets the observations (Section 8)
    │   ├── *_options.ini           run configuration (Section 9)
    │   └── *.m                     MATLAB plotting (Section 12)
    └── delta_impulse/          verification with an instantaneous impulse (Section 11)
```

Two executables are built: `adjoint_forward` and `delta_impulse`. The
second compiles the same sources with a different adjoint source routine
and stops after the representer (Section 11).

## 5. Requirements

| Requirement | Notes |
|---|---|
| C, C++17 and Fortran compilers | The same compilers used to build ForestClaw (e.g. GCC/gfortran). |
| CMake ≥ 3.19 | ForestClaw's minimum. |
| MPI | ForestClaw is built with MPI. The runs reported here used one MPI rank. |
| ForestClaw (`develop` branch) with Clawpack 4.6 | Section 6. Provides p4est and libsc. |
| LAPACK | For `dgesv` in the coefficient solve (Section 14 if it doesn't link). |
| Python 3 with NumPy and Matplotlib | `forward/write_gauges.py` and the `delta_impulse` scripts. `write_gauges.py` also needs `fclaw_analysis.py` from ForestClaw's `python/` directory. |
| MATLAB | Plotting with the Clawpack/ForestClaw MATLAB graphics routines. |

## 6. Building ForestClaw

Smoke2d-4DVar needs ForestClaw's `develop` branch, built with CMake and the
Clawpack solvers enabled. ForestClaw's own instructions are in its
[wiki](https://github.com/ForestClaw/forestclaw/wiki).

ForestClaw needs p4est and libsc. Either build and install p4est first and
point ForestClaw at it (below), or omit `-Dsubmodules=off`, `-DP4EST_ROOT` and
`-DSC_ROOT` and let ForestClaw check out and build the bundled copies as part 
of its own build. 

Here p4est (with libsc) is built and installed first, and ForestClaw is then  
built against it. Each is configured from its own build directory next to the  
source:

```
ForestClaw/
├── p4est/              p4est-build/       (install: p4est-build/local)
└── forestclaw/         forestclaw-build/  (install: forestclaw-build/local)
```

**p4est** (libsc is checked out and installed with it):

```sh
mkdir -p $HOME/ForestClaw && cd $HOME/ForestClaw
git clone https://github.com/cburstedde/p4est.git
mkdir p4est-build && cd p4est-build
P4EST=$HOME/ForestClaw/p4est-build/local
cmake \
    -DCMAKE_INSTALL_PREFIX=${P4EST} \
    -DCMAKE_C_COMPILER=mpicc \
    -DCMAKE_C_FLAGS="-O2 -g -Wall" \
    -Dmpi=on \
    -GNinja \
    ../p4est
ninja && ninja install
```

**ForestClaw**:

```sh
cd $HOME/ForestClaw
git clone -b develop https://github.com/ForestClaw/forestclaw.git
mkdir forestclaw-build && cd forestclaw-build
P4EST=$HOME/ForestClaw/p4est-build/local
FCLAW=$HOME/ForestClaw/forestclaw-build
cmake \
    -DCMAKE_INSTALL_PREFIX=${FCLAW}/local \
    -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_CXX_COMPILER=g++ \
    -DCMAKE_Fortran_COMPILER=gfortran \
    -DCMAKE_C_FLAGS="-O2 -g -Wall" \
    -DCMAKE_CXX_FLAGS="-O2 -g -Wall" \
    -DCMAKE_Fortran_FLAGS="-O2 -g -Wall -Wno-unused-dummy-argument" \
    -DP4EST_ROOT=${P4EST} \
    -DSC_ROOT=${P4EST} \
    -Dclawpack=on \
    -Dmpi=on \
    -Dsubmodules=off \
    -GNinja \
    ../forestclaw
ninja && ninja install
```

- `-Dclawpack=ON` is required (Smoke2d links `FORESTCLAW::CLAWPACK4.6`).
- `-Dmpi=ON` builds with MPI; p4est must be built with MPI as well.
- `-Dsubmodules=off` tells ForestClaw not to build its bundled copies of
  p4est and libsc, and to use the installed ones instead; `P4EST_ROOT` and
  `SC_ROOT` say where they are.
- Optionally, `-Dapplications=OFF` skips ForestClaw's own examples, which are
  not needed here.
- On macOS, use GCC (e.g. Homebrew `gcc-14`, `g++-14`, `gfortran-14`). If the
  link step warns about compact unwind, add
  `-DCMAKE_EXE_LINKER_FLAGS="-Wl,-no_compact_unwind"` (macOS only; this flag
  breaks the link on Linux).
- For a debug build, add `-DP4EST_ENABLE_DEBUG=1 -DSC_ENABLE_DEBUG=1` to the
  p4est C flags and `-DFCLAW_ENABLE_DEBUG=1` to the ForestClaw C/C++ flags.

The ForestClaw install prefix (here `$HOME/ForestClaw/forestclaw-build/local`) is
what Smoke2d-4DVar needs as `FORESTCLAW_ROOT`.

## 7. Building Smoke2d-4DVar

Build out of source, in a directory next to the repository. The
`delta_impulse` run script expects the build directory to be called
`smoke2d-4DVar-build` and to sit next to `Smoke2d-4DVar/` (Section 11).

```sh
git clone https://github.com/PatriciaAzike/Smoke2d-4DVar.git
mkdir smoke2d-4DVar-build && cd smoke2d-4DVar-build
P4EST=$HOME/ForestClaw/p4est-build/local
FCLAW=$HOME/ForestClaw/forestclaw-build/local
cmake \
    -DFORESTCLAW_ROOT=${FCLAW} \
    -DP4EST_ROOT=${P4EST} \
    -DSC_ROOT=${P4EST} \
    -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_CXX_COMPILER=g++ \
    -DCMAKE_Fortran_COMPILER=gfortran \
    -DCMAKE_C_FLAGS="-O2 -g -Wall" \
    -DCMAKE_CXX_FLAGS="-O2 -g -Wall" \
    -DCMAKE_Fortran_FLAGS="-O2 -g -Wall -cpp -Wno-unused-dummy-argument" \
    -DCMAKE_EXE_LINKER_FLAGS="-llapack" \
    ../Smoke2d-4DVar
make -j4
```

- `-cpp` is required: some Fortran sources use `#if` blocks.
- LAPACK is required (the coefficients β are solved with `dgesv`). If it is
  not on the default library path, add `-L/path/to/lib`, e.g.
  `-L/opt/local/lib -llapack` with MacPorts. 
- `P4EST_ROOT`/`SC_ROOT` are needed when p4est was installed separately
  (Section 6).
- On macOS, add `-Wl,-no_compact_unwind` to the linker flags if the link
  warns about compact unwind.
- For a debug build, use `-O0` and add `-DFCLAW_ENABLE_DEBUG=1` to the C/C++
  flags.

This builds:

- `smoke2d-4DVar-build/examples/adjoint_forward/adjoint_forward`
- `smoke2d-4DVar-build/examples/delta_impulse/delta_impulse`

**Alternative (no separate ForestClaw install).** Configuring with
`-Dexternal_build=ON` instead of `-DFORESTCLAW_ROOT=...` makes CMake clone
ForestClaw's `develop` branch into the build directory and build it with MPI
as part of Smoke2d-4DVar. This needs internet access to GitHub.

## 8. Running the assimilation (`adjoint_forward`)

All inputs are read from, and all outputs written below,
`examples/adjoint_forward/`.

```sh
cd Smoke2d-4DVar/examples/adjoint_forward

# 1. Link the executable here.
ln -sf ../../../smoke2d-4DVar-build/examples/adjoint_forward/adjoint_forward .

# 2. Set the observations. The list at the top of forward/write_gauges.py
#    (xm, ym, tm, dm) is the single source: the script writes them into
#    adjoint_options.ini and writes gauges.data into forward/ and model/
#    (needs ForestClaw's python/ directory for fclaw_analysis.py).
cd forward
PYTHONPATH=/path/to/forestclaw/python python3 write_gauges.py
cd ..

# 3. Run (one MPI rank).
./adjoint_forward          # or: mpirun -n 1 ./adjoint_forward
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

The relative asymmetry of $R$ is about 1% in the full-2D configuration
and about 7% in the pseudo-1D one; symmetrizing $R$ changes $\beta$ by at
most 1.6% (Euclidean norm). After the coefficients, the program writes the 
optimal estimate frames and prints ForestClaw's timing summary for each run.

## 9. Configuration reference

The driver reads three option files, one per kind of run. Each has a
problem section (`[*-user]`), a mesh section (`[*-clawpatch]`), a run section
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
| `pseudo-1d` | `1` = pseudo-1D experiment, `2` = full 2D (Section 1). Also selects the velocity field for all runs. |
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

## 10. Output

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

## 11. Verification: the delta-impulse test (`delta_impulse`)

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
cd Smoke2d-4DVar/examples/delta_impulse
./run_delta_impulse.sh                             # sets up run directories, runs
python3 check_delta_impulse.py                     # checks the event frame and the handoff
python3 plot_delta_impulse.py --t 0.0 0.8 1.5      # panels in panels/
```

`run_delta_impulse.sh` looks for the executable in
`../../../smoke2d-4DVar-build/examples/delta_impulse/`; set
`DELTA_IMPULSE_EXE` to use another path. All inputs and outputs stay in
this directory.

For the configured case (one observation at $(1.5, 1.0)$, $t_m = 1.5$,
$\varepsilon = 0.01$), the event frame should have mass ≈ 1, centroid
≈ $(1.5, 1.0)$ and peak ≈ $1/(4\pi\varepsilon) = 7.96$, up to cell-centre
sampling and AMR transfer. `check_delta_impulse.py` reports these and also
checks the adjoint-to-representer handoff at $t=0$ and the location of the
representer peak at $t = t_m$. See
[`examples/delta_impulse/README.md`](examples/delta_impulse/README.md) for
details.

## 12. Plotting

**MATLAB** (`examples/adjoint_forward/*.m`). `setplot2.m` and
`afterframe.m` configure the Clawpack/ForestClaw MATLAB graphics (e.g.
`plotclaw2`, `showpatchborders`), which must be on the MATLAB path. Run
`plotclaw2` in the directory holding the frames to plot. `afterframe.m`
draws the patch borders and, for the 2D case, the exact filament of the
prior disk (`filament_soln.m`).

**Python** (`examples/delta_impulse/`). `plot_delta_impulse.py` draws the
adjoint and representer at chosen physical times with one shared colour
scale; `plot_delta.py` holds the reading and panel routines it uses.

## 13. Uniform-mesh reference runs

A uniformly refined reference is obtained by setting `minlevel = maxlevel`
in each option file (level 4 for the prior and representers, level 5
for the adjoints, i.e. $512^2$ and $1024^2$ cells). For the configuration in
this repository (full 2D, three observations, one MPI rank), the adaptive
run took $1293 \pm 15$ s and the uniform reference $23{,}667 \pm 619$ s
(mean ± sample standard deviation of three runs each), a speedup of about
18.3. Almost all of the saving is in the adjoint solves.

## 14. Troubleshooting

| Symptom | Cause and fix |
|---|---|
| `find_package(FORESTCLAW)` fails | Pass `-DFORESTCLAW_ROOT=<install prefix>`, and make sure ForestClaw was built with `-Dclawpack=ON`. |
| Undefined `dgesv_` at link time | No LAPACK was linked through ForestClaw's dependencies. Add `find_package(LAPACK REQUIRED)` and link `LAPACK::LAPACK` to `adjoint_forward` and `delta_impulse`. |
| `get_gauges: FATAL -- accumulator has N gauges but mdata = M` | `gauges.data` is out of date. Re-run `forward/write_gauges.py` (from inside `forward/`); it rewrites both `forward/gauges.data` and `model/gauges.data`. |
| The prior (`model/`) is zero everywhere and never refines | The initial condition is wrong or missing, e.g. `adjoint/adjoint_fdisc.f90` (the disk for `pseudo-1d = 2`) not compiled in. |
| `File does not exist` when "Restarting model from checkpoint file" | A checkpoint in `model/` (or `adjoint{j}/`, `forward{j}/`) was deleted or moved during the run. Re-run without touching those directories. |
| Stale results after changing options | Output directories are overwritten, not cleared; delete old `fort.*` files before a run if frame counts change. |

## 15. Known limitations

- Tested on one MPI rank. The overlap exchange is written for distributed
  meshes, but multi-rank runs have not been validated.
- The forward and adjoint discretizations are not exact discrete adjoints,
  so $R$ is symmetrized before the solve (Section 3).
- Only output style 1 (equally spaced output times) is supported by the
  one-step time stepper.
- The velocity field is analytic and steady.

## 16. References

- A. F. Bennett, *Inverse Methods in Physical Oceanography*, Cambridge
  University Press, 1992.
- A. F. Bennett, *Inverse Modeling of the Ocean and Atmosphere*, Cambridge
  University Press, 2005.
- D. Calhoun and C. Burstedde, "ForestClaw: A parallel algorithm for
  patch-based adaptive mesh refinement on a forest of quadtrees,"
  arXiv:1703.03116, 2017.
- C. Burstedde, L. C. Wilcox and O. Ghattas, "p4est: Scalable algorithms for
  parallel adaptive mesh refinement on forests of octrees," *SIAM J. Sci.
  Comput.* 33(3):1103–1133, 2011.
- R. J. LeVeque, *Finite Volume Methods for Hyperbolic Problems*, Cambridge
  University Press, 2002.

## 17. License

BSD 3-Clause; see [LICENSE](LICENSE). ForestClaw is distributed under its
own BSD-style license, and p4est under the GPL (version 2 or later).
