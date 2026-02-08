# AGENTS.md - DFTB+ Development Guide

## Main directive
Work through the `TODO.md` step by step. Do not stop before finishing all tasks, defer outputting a summary of changes until the todo is complete.

Take your time in order to do things properly. Use git -- Whenever an atomic task is complete, remove it from the README and commit it. If a task changed anything important with regard to the project, update `AGENTS.md` correspondingly. 
The project is currently under development / undergoing major refactoring, thus breaking changes are allowed with no further consideration required, as long as the new behaviour is documented accordingly.
If during a task it becomes apparent that some hsd-fortran behaviour is suboptimal, instead of working around it prefer to either fix it on the spot or tack it onto the `TODO.md` to resolve it later on.

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

# Unit tests
ctest --test-dir build -R "unit" 2>&1 | tail -5

# Application tests
ctest --test-dir build -R "app/dftb+" -j6 --stop-on-failure --output-on-failure | tee testlog.txt | tail -5
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

**Always bundled:** `xmlf90`, `ddcosmo`, `fypp`
**Optional hybrid:** `mpifx`, `scalapackfx`, `mbd`, `libnegf`, `s-dftd3`, `dftd4`, `tblite`, `fortuno`, `chimes`, etc.
**System required:** BLAS, LAPACK

## Source Organization

All library source lives in `src/dftbp/`. Module naming convention: `dftbp_<subdir>_<module>`.

| Directory | Purpose |
|-----------|---------|
| `common/` | Accuracy, assertions, atomic data, environment, file I/O, timers, unit conversion |
| `type/` | Core data types: geometry, orbitals, Slater-Koster data, neighbor lists |
| `io/` | **HSD parser**, HSD utilities, XML utilities, line/token readers, logging |
| `math/` | Mathematical routines |
| `dftb/` | Core physics: charges, Coulomb, dispersion, DFT+U, Hamiltonian, SCC |
| `dftbplus/` | High-level: input parsing, initialization, main driver, output |
| `dftbplus/input/` | Input processing sub-modules |
| `extlibs/` | Wrapper modules for external libs (xmlf90, LAPACK, ARPACK, etc.) |
| `api/mm/` | Public API (C/Fortran bindings) |
| `include/` | Fypp include files (`common.fypp`, `error.fypp`) |
| Other: `derivs/`, `elecsolvers/`, `geoopt/`, `md/`, `mixer/`, `reks/`, `timedep/`, `solvation/`, `transport/`, `poisson/`, `geometry/` | Domain-specific physics modules |

## HSD Parsing Architecture

HSD (Human-friendly Structured Data) is the input format for DFTB+. The parsing system is a three-layer architecture:

### Layer 1: XML DOM (external/xmlf90)

`xmlf90` provides the in-memory tree representation. Key types:
- `fnode` — universal DOM node (element, text, document)
- `fnodeList` — linked list of nodes
- `fnamedNodeMap` — node attributes (name→value)

Wrapper: `src/dftbp/extlibs/xmlf90.F90` re-exports needed symbols.

### Layer 2: HSD Parser (src/dftbp/io/)

| Module | File | Purpose |
|--------|------|---------|
| `dftbp_io_hsdparser` | `hsdparser.F90` | Core parser: HSD text → xmlf90 DOM tree; `parseHSD()`, `dumpHSD()` |
| `dftbp_io_hsdutils` | `hsdutils.F90` | Typed value extraction: `getChildValue()` (21 variants), `setChildValue()`, `getChild()`, `detailedError()` |
| `dftbp_io_hsdutils2` | `hsdutils2.F90` | Unit conversion (`convertUnitHsd`), unprocessed node tracking, node renaming |
| `dftbp_io_xmlutils` | `xmlutils.F90` | DOM tree navigation: `getFirstChildByName()`, `getChildrenByName()`, etc. |

### Layer 3: Application (src/dftbp/dftbplus/)

| Module | File | Purpose |
|--------|------|---------|
| `dftbp_dftbplus_hsdhelpers` | `hsdhelpers.F90` | Orchestrates input parsing: reads file → builds DOM → extracts data → dumps processed |
| `dftbp_dftbplus_parser` | `parser.F90` | Maps DOM tree → `TInputData` (8500+ lines, reads all DFTB+ keywords) |
| `dftbp_dftbplus_oldcompat` | `oldcompat.F90` | Backward compatibility: converts old parser versions (1–13) to current (14) |

### HSD Parsing Flow

```
dftb_in.hsd (HSD text)
    ↓  hsdparser.parseHSD() — hand-written recursive descent parser
fnode* (xmlf90 DOM tree in memory)
    ↓  parser.parseHsdTree() → hsdutils.getChildValue/getChild
TInputData (Fortran data structures)
    ↓  hsdparser.dumpHSD()
dftb_pin.hsd (processed output)
```

### HSD Node Attributes

When an HSD node is parsed into the DOM tree, metadata is stored as XML attributes:
- `"name"` — original HSD name (preserves capitalization)
- `"m"` — modifier (e.g., `"Angstrom"` from `[Angstrom]`)
- `"l"` — list flag
- `"start"` / `"end"` — source line numbers
- `"file"` — source file name
- `"proc"` — processed flag (set after parsing to detect unrecognized keywords)

### Module Dependency Graph

```
xmlf90 (external library)
  └── dftbp_extlibs_xmlf90          ← re-export wrapper
        ├── dftbp_io_xmlutils        ← XML tree utilities
        ├── dftbp_io_hsdparser       ← core HSD parser
        │     └── dftbp_io_hsdutils  ← typed value get/set
        │           └── dftbp_io_hsdutils2  ← unit conversion, tracking
        │                 ├── dftbp_dftbplus_parser      ← input parser
        │                 ├── dftbp_dftbplus_oldcompat    ← version compat
        │                 └── dftbp_dftbplus_hsdhelpers   ← orchestration
        └── dftbp_hsdapi             ← public API façade
```

### Key Integration Points for hsd-fortran

To replace the legacy HSD parser with hsd-fortran, the integration points are:

1. **`dftbp_io_hsdparser`** — Replace `parseHSD()` with hsd-fortran's `hsd_load()`
2. **`dftbp_io_hsdutils`** — Replace `getChildValue()` / `setChildValue()` with hsd-fortran's `hsd_get()` / `hsd_set()`
3. **`dftbp_io_xmlutils`** — Replace DOM navigation with hsd-fortran's `hsd_has_child()`, iteration, etc.
4. **`dftbp_io_hsdutils2`** — Replace unprocessed node tracking, unit conversion helpers
5. **`dftbp_extlibs_xmlf90`** — Can be removed once xmlf90 is fully replaced
6. **`dftbp_dftbplus_parser`** — Must be updated to use hsd-fortran API (largest change, 8500+ lines)
7. **`dftbp_dftbplus_hsdhelpers`** — Simplify orchestration to use hsd-fortran directly

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
- 7 test suites: `common`, `dftb`, `include`, `io`, `math`, `mixer`, `type`
- Run: `ctest --test-dir build -R "unit"`

### Application Tests (test/app/dftb+/)
- ~400 regression tests across categories (scc, non-scc, spin, dispersion, md, transport, etc.)
- Each runs `dftb+` binary and compares output against reference data
- Require Slater-Koster data files in `external/slakos/`
- Run: `ctest --test-dir build -R "app/dftb+"`

