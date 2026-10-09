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
