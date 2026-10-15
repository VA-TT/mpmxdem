# Geometry-triggered RVE replacement

An optional DEM -> Mohr-Coulomb transition is enabled by:

```
model MohrCoulomb MC_RVE E nu phi c psi
CHCL.replaceRVE 1.5 MC_RVE
```

E and c are in Pa; phi and psi are in radians. The MC model must be declared.
The option is disabled by default. A zero factor disables it explicitly.
The setting is written to conf files and reset when a viewer clears a conf.

For each CHCL_DEM or CHCL_DEM_Multi MP, the code checks the present periodic
cell and the cell predicted by `F * inverse(previous F)`. The shortest distance
between parallel opposite faces is `det(h) / max(|b x c|,|c x a|,|a x b|)`,
where a,b,c are columns of h. This also handles skewed cells. Replacement
occurs before advancing DEM if either thickness is <= factor * (2 Rmax).
Invalid cell geometry is rejected. The factor 1.5 comes from the supplied
article; it is not calibrated for this sample (its Biax stress spike starts
before the 1.5 threshold).
The new model starts from the existing stress and total/plastic strains.
Mass, volume, density, position, velocity and F are preserved. The current
strain increment is applied once by MC, not by both laws. The first MC return
mapping can still project an inherited DEM stress outside the MC yield surface;
continuity at handoff does not guarantee unchanged stress after this update. DEM storage is
released and DEM tracking ends. Future steps use MC; model identifiers in
conf files retain the switch on restart. Each event prints immediately to stderr, even with the logger disabled:
`[RVE replacement] DEM -> Mohr-Coulomb: MP=... t=... s ...`.
The message includes MP number, time,
current and predicted thickness ratios and the target model. There is no
reverse switch and no reconstruction of past DEM plastic strains.

The existing MohrCoulomb is a 2D plane-strain return mapping, with sigma_zz
integrated and assumed intermediate, not a fully general 3D MC solver.
Initial sigma_xz/yz are preserved at transition but this MC model has no
subsequent evolution equation for those components. Microstructure,
coordination and actual DEM-cell porosity cease to be available after the
switch. Bulk density continues to evolve from MPM mass/volume; it is not
reset to the initial density when the model changes.

Example calibrated model (angles in radians):

```
model MohrCoulomb MC_RVE 14346192 0.4034179 0.3810024848832474 0 0.004460601499306681
CHCL.replaceRVE 1.5 MC_RVE
```

Calibration simulations and generated plots are external data, not shipped
with this solver change.

The example calibration uses rho_initial=1603.72072810721. Default residual parameters come from Biax
I=1e-4 at eps_d=60%: E=14.346192 MPa, nu=0.4034179,
phi=21.82983437 deg, psi=0.25557364 deg, c=0 assumed. E is fitted through
0.2% axial strain and nu follows the initial volumetric branch, using the
same apparent-parameter method as the previous VER8000 scripts. Residual
psi uses a total-strain fit over eps_d=57.5--62.5%; it is not a separate
measurement of plastic dilation. A single confinement does not identify c.
Triax diagnostics are also exported, with residual taken at axial strain
30%; that is not the same strain measure/state as Biax eps_d=60%.

Build and check without replacing the existing executable:

```
cmake -S . -B BUILD_RVE_REPLACEMENT
cmake --build BUILD_RVE_REPLACEMENT --target mpmbox check_rve_replacement -j 4
BUILD_RVE_REPLACEMENT/check_rve_replacement
python3 Tests/check_rve_replacement_run.py --simulation /path/to/75-MP-case
```

The C++ check is self-contained. The Python integration check additionally
requires an external 75-MP case with `input_RVE_replacement.txt` and
`VER1000.txt`; it does not depend on the plotting scripts.

Run the new executable from a simulation directory using the new input.
A full collapse run with this option has not been performed by these tests.
