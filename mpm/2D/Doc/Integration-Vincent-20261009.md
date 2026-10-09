# Current policy: prefer upstream implementations

After review on 2026-10-09, the production configuration reader and writer
were restored to Vincent's implementation. All new configurations, including
planeStrain runs, use the upstream 36-column MP layout. The custom confFormat
2/3 writer and reader are no longer used. Stress is stored as the 2x2 block plus
sigma_zz, sigma_xz and sigma_yz; deformation gradient F remains 2x2 with F_zz=1
implicit in plane strain. This is still a 2D MPM solver, not a general 3D solver.

The previous integration is retained at
backup/before-upstream-priority-20261009. The original local code and stash
also remain available. Previously generated confFormat 2/3 results require
the retained older reader; changing their header does not convert them.

Local see/rheology workflows now call Vincent's Post/ScalarTools directly for
stress invariants, principal friction, strain rate and equivalent shear strain.
The inertial number uses d*gdot*sqrt(rho_s/P), matching mpmpost. Legacy values
can change, especially for volume-changing flows, cohesion and nonzero xz/yz
shears. Local workflows are retained; the numerical definitions are upstream.

CHCL_DEM_Multi and its deterministic RVE assignment remain local extensions.
Core/MPMbox.cpp differs from upstream only in their include and registration.
HookeElasticity.cpp is identical to upstream. Build fixes and local tools remain.

The mpmxdem-vincent worktree has been removed. The active checkout is
~/mpmxdem on VA_modification; MPMrun still points to its BUILD/mpmbox.

---

The following is the historical record of the first integration, before this
upstream-priority review. Its custom-format choices are superseded above.

# Vincent upstream integration — 2026-10-09

The integration branch starts from VA_modification (179bfde48) and merges
upstream/main (148e5f796). No history is rebased and nothing is pushed.

Backups:
- backup/VA_modification-before-vincent-20261009
- backup/main-before-vincent-20261009
- stash: Before update 2026-10-08 VA_modification (kept after apply)

Retained: CHCL_DEM_Multi, plane-strain tensor output, custom see modes,
rheology runner, local spy sources and documentation. Large simulation data
and third-party dependency directories remain untracked on disk.

Integration adaptations:
- CHCL_DEM_Multi initializes sigma_xz and sigma_yz like CHCL_DEM.
- Standard conf output uses Vincent's 36-column MP layout.
- Custom plane-strain output uses confFormat 3 (60 columns): the existing
  five 3x3 tensors, followed by sigma_xz, sigma_yz and contact force. The
  reader retains confFormat 2 support within a supported file version.
- Out-of-plane plastic strain now lives in MPMbox::modelState, following
  Vincent's memory layout.
- Tests/mpmconf.py reads standard and custom layouts, including the old
  34-column layout for test-side inspection.
- CMake builds rheology against mpmbox_lib and puts dependency build files
  under BUILD to avoid the external Toofus directory target error.
- Local B-spline clamping is superseded by Vincent's locateElement and ghost
  nodes: computation stops if a particle leaves the grid; visualization clamps
  the element index. This preserves upstream validation instead of silently
  continuing an out-of-grid simulation.
- Rheology tau/P exports include the newly stored out-of-plane shears.

Old conf headers are still subject to Vincent's September 2026 version check.
Changing only the header of an old result does not convert its data layout.
The retained legacy backup is available when old results need the old reader.

The optional mpmxdem-vincent worktree is a reference checkout, not the active
working directory. MPMrun continues to point to ~/mpmxdem/mpm/2D/BUILD/mpmbox.

Validation:
- Full CMake build: mpmbox, see, mpmpost, rheology succeeded.
- Tests/runtests.py: 32 PASS, 0 FAIL after the final solver adaptation.
- Tests/benchmark.py: completed all four cases; no saved baseline is available,
  so this run does not establish a speedup or numerical equivalence to old code.
- Two-RVE smoke: standard and plane-strain output, eight MPs, two configurations
  each, finite values; C++ read and raw rheology post-processing succeeded.
- The local rheology runner initializes the profiler required by upstream's
  profiling-enabled solver library.

Validation logs for this session are in /tmp/mpmxdem-integrated-*.log and
/tmp/mpmxdem-integrated-smoke-c7duxp9z. These temporary paths are not backups.

## Upstream-priority validation

- Full build succeeded for mpmbox, see, mpmpost and rheology.
- Regression suite: 32 PASS, 0 FAIL.
- Two-RVE standard and planeStrain smoke runs wrote and read the upstream
  layout successfully; raw rheology post-processing completed.
- A separate parity check imposed an affine planar expansion on eight MPs
  with nonzero out-of-plane shears. Legacy rheology and Vincent's mpmpost
  Scalars agreed on I, sinPhi, mu and Phi within 1e-5 relative / 1e-6 absolute
  tolerance (mpmpost prints rounded values). Different RVE diameter/density
  parameters were matched to each MP's assigned RVE.
- Four benchmark cases completed. No saved old baseline is available.

Logs: /tmp/mpmxdem-priority-*.log and
/tmp/mpmxdem-upstream-compare-04w24je6 (temporary session artifacts).


## Final overlap audit

The remaining local differences were reviewed against upstream/main:
- OneStep, ShapeFunctions, BoundaryForceLaw, all upstream constitutive models,
  MaterialPoint, Post and DEM numerical sources are unchanged from upstream.
- MPMbox adds only CHCL_DEM_Multi registration. Its init now calls CHCL_DEM::init
  through a local initializer, so allocation and initial stress transfer use
  Vincent's implementation directly; stress updates are inherited unchanged.
- set_MP_grid assigns MP.nb before init and keeps IDs unique across calls. This
  is needed for the deterministic multi-RVE extension, not a constitutive change.
- Viewer pressure coloring and key 6 now use the same 3D pressure/shear-rate
  definitions as Vincent's Scalars. The sig_yy display is consistent end-to-end.
- RVE diameter comes from the selected initial RVE, avoiding stale/misindexed
  ADs data and matching the standalone rheology runner.
- Empty viewer data is guarded before minmax dereferences.
- Makefile includes ScalarTools and declares archive dependencies; CMake is the
  tested build path. The full handwritten Makefile build has not been tested.
- BeamDev and PointsInfo remain preserved local source files; neither is
  currently registered in the upstream Spy factory, so they do not alter runs.

Validation after these changes:
- Full normal CMake build succeeded.
- 32/32 regression tests passed.
- Two-RVE standard/planeStrain runs and C++ configuration reading passed.
- I, sinPhi, mu, Phi parity with mpmpost passed again.
- Optional headless test BUILD/check_viewer_formulas passed: analytical 3D stress
  with xz/yz shear, planar expansion, Hencky strain and actual pressure/rate colors.
  Rebuild with: cmake --build BUILD --target check_viewer_formulas.

No simulation datasets are staged. Existing backups and stash are retained.
This audit does not establish correctness for every possible physical setup.
