# delta_impulse: instantaneous delta-impulse verification

The observation is applied as an *instantaneous* event in time (a delta in
time as well as a Gaussian in space); the name is shortened to `delta_impulse`.

This example reuses the `adjoint_forward` sources and leaves the production
`adjoint_forward` executable unchanged.  It replaces only the adjoint source update: the spatial
heat kernel is added once when the reverse-time integration crosses
`tau_m = T - t_m`.  There is no temporal Gaussian and no multiplication by the
time step.  The verification executable stops after constructing the unit
representer; it does not continue into the observation-fit coefficient solve.

For the configured full-2D case, `epsilon = 0.01`, so the event frame at
physical time `t_m = 1.5` should have mass approximately one, centroid
approximately `(1.5, 1.0)`, and maximum approximately
`1/(4*pi*epsilon) = 7.9577`, subject to cell-centre sampling and AMR transfer.

The event time corresponds to reverse time `tau_m = 2.0 - 1.5 = 0.5`, which is
one of the `nout = 20` output times.  The step therefore lands on the event
before the frame is written.  This frame represents the adjoint state just after the
jump, `alpha(t_m^-)`, so it contains the full spatial Gaussian. With a temporal Gaussian, 
a frame at `t_m` would hold only about half of it, since only half the temporal kernel 
has been added by its centre.

## Build

From the Smoke2d-4DVar build directory:

```sh
cmake --build . --target delta_impulse -j2
```

## Run

From this directory:

```sh
./run_delta_impulse.sh
```

All runtime inputs and outputs stay here.  In particular, the raw files are
written to `adjoint0/`, `forward0/`, and `model/` below this directory, rather
than overwriting the corresponding parent directories.

After the run:

```sh
python3 check_delta_impulse.py
```

The check reports the mass, location, peak magnitude, and AMR levels in the
event frame.  It also checks the adjoint-to-representer handoff at `t = 0` and
the location of the refocused representer peak at `t = t_m`.

Generate the corresponding physical-time panels with:

```sh
python3 plot_delta_impulse.py --t 0.0 0.8 1.5
```

They are written to `panels/` in this directory.  One color scale is shared by
both rows and all requested times, so the same color always represents the same
field value.
