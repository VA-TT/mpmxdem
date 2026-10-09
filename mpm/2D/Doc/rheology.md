# Standalone rheology processing

The `rheology` executable processes MPM configuration files without opening the
OpenGL `see` application. It reuses `MPMbox::read()` and `MPMbox::postProcess()`
so the interpolation and smoothed stress/velocity-gradient values are the same
as the previous `see` workflow.

Build it with CMake:

```sh
cmake -S . -B BUILD
cmake --build BUILD --target rheology -j$(nproc)
```

Run it from any directory containing `confN.txt` files and the referenced DEM
configuration files:

```sh
BUILD/rheology /path/to/results --bins 80 --min-count 2
```

By default the tool uses the same smoothed stress, velocity gradient, and
out-of-plane stress as `see`. To use the stress, out-of-plane stress, and
density stored directly on each material point:

```sh
BUILD/rheology /path/to/results --mode raw
```

The configuration format does not store the raw velocity gradient. Therefore
`raw` mode still uses the velocity gradient reconstructed by
`MPMbox::postProcess()` for the inertial number `I`; the output header records
this explicitly. Use `smooth` when you want exact compatibility with the old
`see` rheology export, and `raw` when you want to inspect the effect of the
unsmoothed MP stress/density fields.

The tool writes, in the results directory:

- `rheology_confN.txt`: MP index, position, inertial number `I`,
  `sin(phi)`, rheological friction `tau/p`, and packing fraction `Phi`.
- `mean_I_time.dat`: particle-averaged `I` for every configuration.
- `binned_mu_stats.dat`: logarithmic-`I` bins for `mu(I)`.
- `binned_packing_stats.dat`: logarithmic-`I` bins for `Phi(I)`.

The calculation uses the same definitions as `See/see.cpp`: smoothed velocity
and stress, 2D or 3D pressure according to the DEM flag, representative DEM
diameter for `I`, Mohr friction, deviatoric `tau/p`, and DEM grain density for
`Phi`. `CHCL_DEM_Multi` is supported by selecting the RVE associated with each
material point.

The existing gnuplot scripts can plot the generated binned files. For example:

```sh
cd /path/to/results
gnuplot -e "datafile='rheology_conf23.txt'" mu_I_binned.plt
gnuplot -e "datafile='rheology_conf23.txt'" packing_fraction_binned.plt
gnuplot mean_I_time.plt
```
