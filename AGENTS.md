# AGENTS.md - DFTB+ Development Guide

## Main directive
Work through the `TODO.md` step by step. Do not stop before finishing all tasks, defer outputting a summary of changes until the todo is complete.

Take your time in order to do things properly. Use git — whenever an atomic task is complete, remove it from the TODO and commit it. If a task changed anything important with regard to the project, update `AGENTS.md` correspondingly.

The HSD IO migration from xmlf90 to hsd-fortran/hsd-data is **complete**. All legacy code (xmlf90, hsdcompat, re-export wrappers) has been removed. The codebase now uses hsd-fortran and hsd-data directly.

## Project Overview

DFTB+ is a density-functional tight-binding (DFTB) simulation package. It is a large Fortran project (~400+ source files) built with CMake and the Fypp preprocessor. The codebase supports optional MPI, OpenMP, GPU, and numerous external library integrations.

## Quick Reference

```bash
# Configure (minimal build)
cmake -B build -DCMAKE_BUILD_TYPE=Debug \
  -DWITH_OMP=FALSE -DWITH_MPI=FALSE -DWITH_API=TRUE \
  -DWITH_UNIT_TESTS=TRUE \
  -DCMAKE_Fortran_FLAGS="-ffree-line-length-none"
  
# Download Slater Koster data
./utils/get_opt_externals ALL

# Build
cmake --build build -j$(nproc) 2>&1 | tail -5

# Unit tests (-R "unit"), Application tests (-R "app/dftb+"):
ctest --test-dir build -j6 --stop-on-failure --output-on-failure | tee testlog.txt | tail -5
```

Application tests take a very long time to run.
Therefore prefer to run a subset whenever possible, and always employ parallelism and use --stop-on-failure --output-on-failure.
Do not run the tests multiple time to grep for things, instead direct the output to a file the first time round.
Ideally, distill iterated upon behaviour into a unit test.

## Build System

### CMake + Fypp

- **CMake** orchestrates the build. Top-level `CMakeLists.txt` has all options.
- **Fypp** (`external/fypp/bin/fypp`) is a Python-based Fortran preprocessor. Source files are `.F90` (with Fypp directives like `#:if WITH_MPI`, `@:ASSERT`). They are preprocessed to `.f90` in the build directory.
- Fypp flags are derived from build options via `dftbp_add_fypp_defines()` in `cmake/DftbPlusUtils.cmake`.
- The include file `src/dftbp/include/common.fypp` defines all conditional compilation macros.
- All source files in a subdirectory are collected into `ALL-SOURCES-FPP`, preprocessed by `dftbp_preprocess()`, then compiled into a single library `dftbplus`.

### Key Build Options

| Option | Default | Description |
|--------|---------|-------------|
| `WITH_API` | `TRUE` | Build public API (required for `isAsiCallbackEnabled` field) |
| `WITH_MPI` | `FALSE` | MPI parallelism |
| `WITH_OMP` | `TRUE` | OpenMP parallelism |
| `WITH_GPU` | `FALSE` | GPU acceleration |
| `WITH_UNIT_TESTS` | `FALSE` | Build Fortuno unit tests |
| `WITH_SOCKETS` | `FALSE` | i-PI socket communication |
| `WITH_TBLITE` | `FALSE` | xTB support |

### External Dependency Management

Dependencies use a **hybrid resolution** system via `dftbp_config_hybrid_dependency()`:

1. **Submodule** — source from `external/<pkg>/origin/` (git submodule)
2. **Find** — CMake `find_package()` for pre-installed libraries
3. **Fetch** — CMake `FetchContent` download from git

The method order is set by `HYBRID_CONFIG_METHODS` (default: `"Submodule;Find;Fetch"`).

**Always bundled:** `hsd-data` (brings hsd-fortran), `ddcosmo`, `fypp`
**Optional hybrid:** `mpifx`, `scalapackfx`, `mbd`, `libnegf`, `s-dftd3`, `dftd4`, `tblite`, `fortuno`, `chimes`, etc.
**System required:** BLAS, LAPACK

## Source Organization

All library source lives in `src/dftbp/`. Module naming convention: `dftbp_<subdir>_<module>`.

| Directory | Purpose |
|-----------|---------|
| `common/` | Accuracy, assertions, atomic data, environment, file I/O, timers, unit conversion |
| `type/` | Core data types: geometry, orbitals, Slater-Koster data, neighbor lists |
| `io/` | HSD utilities (wrappers around hsd-fortran/hsd-data), unit conversion, tagged output, line/token readers, logging |
| `math/` | Mathematical routines |
| `dftb/` | Core physics: charges, Coulomb, dispersion, DFT+U, Hamiltonian, SCC |
| `dftbplus/` | High-level: input parsing, initialization, main driver, output |
| `dftbplus/input/` | Input processing sub-modules |
| `extlibs/` | Wrapper modules for external libs (LAPACK, ARPACK, etc.) |
| `api/mm/` | Public API (C/Fortran bindings) |
| `include/` | Fypp include files (`common.fypp`, `error.fypp`) |
| Other: `derivs/`, `elecsolvers/`, `geoopt/`, `md/`, `mixer/`, `reks/`, `timedep/`, `solvation/`, `transport/`, `poisson/`, `geometry/` | Domain-specific physics modules |

## HSD Parsing Architecture (post-migration)

HSD (Human-friendly Structured Data) is the input format for DFTB+. The IO
stack has been fully migrated from the legacy xmlf90 DOM to **hsd-fortran** /
**hsd-data**.

### Layer 1: hsd-fortran + hsd-data (external libraries)

- **hsd-fortran** provides the in-memory tree (`hsd_table`, `hsd_value`),
  parser (`hsd_load`), serializer (`hsd_dump`), and typed accessors (`hsd_get`,
  `hsd_set`, `hsd_get_matrix`, etc.).
- **hsd-data** adds multi-format IO (`data_load`, `data_dump`) for XML, JSON,
  TOML, HDF5 on top of hsd-fortran.

Integration: `add_subdirectory(external/hsd-data)` →
`target_link_libraries(dftbplus PUBLIC hsd-data)`.

### Layer 2: DFTB+ IO Wrappers (src/dftbp/io/)

| Module | File | Purpose |
|--------|------|---------|
| `dftbp_io_hsdutils` | `hsdutils.F90` | DFTB+-specific convenience wrappers: `getChildValue` (22 overloads), `getChild`, `setChildValue`, `dftbp_error`/`dftbp_warning`, child-list iteration, atom/index selection, processed-flag management |
| `dftbp_io_hsderror` | `hsderror.F90` | Error bridge: `hsd_fatal_error`, `hsd_warning`, `hsd_error_from_stat` |
| `dftbp_io_unitconv` | `unitconv.F90` | Unit conversion: `convertUnitHsd` maps HSD modifier strings → DFTB+ `TUnit` system |
| `dftbp_io_unitconvfuncs` | `unitconvfuncs.F90` | Per-category converter callbacks for `hsd_get_with_unit` |
| `dftbp_io_taggedoutput` | `taggedoutput.F90` | Multi-format result output (tag/HSD/JSON/XML) from `hsd_table` trees |

### Layer 3: Application (src/dftbp/dftbplus/)

| Module | File | Purpose |
|--------|------|---------|
| `dftbp_dftbplus_hsdhelpers` | `hsdhelpers.F90` | Orchestrates input parsing: `data_load` → tree → `hsd_dump` processed |
| `dftbp_dftbplus_parser` | `parser.F90` | Maps `hsd_table` tree → `TInputData` (~3000 lines) |
| `dftbp_dftbplus_oldcompat` | `oldcompat.F90` | Backward compatibility: converts old parser versions (1–13) to current (14) |

### HSD Parsing Flow

```
dftb_in.hsd (HSD text)
    ↓  data_load() — hsd-data dispatches to hsd-fortran parser
hsd_table (hsd-fortran tree in memory)
    ↓  parser.parseHsdTree() → hsdutils.getChildValue/getChild
TInputData (Fortran data structures)
    ↓  hsd_dump()
dftb_pin.hsd (processed output)
```

### Module Dependency Graph

```
hsd-data (external/hsd-data)
  └── hsd-fortran (transitive dependency)
        │
dftbp_io_hsdutils  ← uses hsd + hsd_data types directly
├── dftbp_io_hsderror       ← error formatting via hsd_format_error
├── dftbp_io_unitconv       ← unit conversion via TUnit + hsd_set
├── dftbp_io_taggedoutput   ← result serialization
│
├── dftbp_dftbplus_parser      ← input parser
├── dftbp_dftbplus_oldcompat   ← version compat
└── dftbp_dftbplus_hsdhelpers  ← orchestration
```

### What Was Removed

- `external/xmlf90/` — legacy XML DOM library (deleted)
- `dftbp_extlibs_xmlf90` — re-export wrapper (deleted)
- `dftbp_extlibs_hsddata` — re-export wrapper (deleted)
- `hsdcompat.F90` — transitional compat shim (deleted)
- `hsdutils2.F90` / `hsdparser.F90` / `xmlutils.F90` — replaced by hsd-fortran API
- All `detailedError` / `detailedWarning` calls → `dftbp_error` / `dftbp_warning`

## Application Entry Point

`app/dftb+/dftbplus.F90`:

```
program dftbplus
  1. initGlobalEnv()         — init MPI/OpenMP/global state
  2. printDftbHeader()       — version banner
  3. parseHsdInput(input)    — parse dftb_in.hsd → TInputData
  4. TEnvironment_init(env)  — set up compute environment
  5. main%initProgramVariables(input, env) — init calculation
  6. runDftbPlus(main, env)  — execute DFTB+ calculation
  7. Cleanup
```

## Test Infrastructure

### Unit Tests (test/src/dftbp/unit/)
- Enabled with `WITH_UNIT_TESTS=TRUE` (serial builds only)
- Uses **Fortuno** framework
- 7 test suites: `common`, `dftb`, `include`, `io`, `math`, `mixer`, `type` (+ 1 integration test suite)
- Run: `ctest --test-dir build -R "unit"`

### Application Tests (test/app/dftb+/)
- ~400 regression tests across categories (scc, non-scc, spin, dispersion, md, transport, etc.)
- Each runs `dftb+` binary and compares output against reference data
- Require Slater-Koster data files in `external/slakos/`
- Run: `ctest --test-dir build -R "app/dftb+"`

