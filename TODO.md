## Phase 0: Preparation & Dependency Wiring (SPECIFICATION.md §3, Phase 0)

- Add `hsd-data` as a hybrid dependency (Submodule/Find/Fetch) in the top-level
  `CMakeLists.txt`, following the existing pattern for other hybrid dependencies.
  hsd-data auto-fetches hsd-fortran, so only hsd-data needs explicit wiring.
  Place the source under `external/hsd-data/`.

- Create the wrapper module `src/dftbp/extlibs/hsddata.F90` that does
  `use hsd_data` and re-exports the public API under the DFTB+ naming convention
  (`dftbp_extlibs_hsddata`). Mirror the existing `dftbp_extlibs_xmlf90` pattern.
  Register it in `src/dftbp/extlibs/CMakeLists.txt`.

- Verify the build compiles cleanly with both old (xmlf90) and new (hsd-data)
  dependencies present: `cmake --build build -j$(nproc)`.

- Add a smoke-test under `test/src/dftbp/unit/io/` that loads a simple HSD
  string via `hsd_load_string`, retrieves a value with `hsd_get`, and checks
  the result. This proves the dependency wiring works end-to-end.

## Phase 1: Adapter Shim & Bridge Modules (SPECIFICATION.md §3, Phase 1; §4; §5)

- Create `src/dftbp/io/unitbridge.F90`: bridge module that wraps the existing
  `TUnit`-based unit conversion system (`unitconversion.F90`) behind the
  `hsd_get_with_unit` converter function interface. Provide one converter
  function per physical dimension (length, energy, time, force, etc.).

- Create `src/dftbp/io/hsderror.F90`: bridge module providing
  `hsd_fatal_error(node, msg)` and `hsd_warning(node, msg)` that extract
  `node%name` / `node%line` and delegate to `dftbp_io_message`'s `error()` /
  `warning()`.

- Create the transitional shim `src/dftbp/io/hsdcompat.F90` that exposes
  legacy call signatures (`getChildValue`, `getChild`, `setChildValue`,
  `detailedError`, `convertUnitHsd`) implemented on top of `hsd_table` /
  `hsd_value`. This enables incremental parser conversion.

## Recursive Task Directive

> Once all items above are exhausted, **replenish** a couple of new actionable
> points by considering the current project state and surveying
> `SPECIFICATION.md` and the codebase. Keep this directive at the **bottom** of
> the TODO list. Resume working through the list sequentially.
>
> **Stop condition:** no further actionable points can be generated because the
> project state fully and without question meets `SPECIFICATION.md`.
