# SPECIFICATION: Migrating DFTB+ IO from Legacy HSD/xmlf90 to hsd-fortran via hsd-data

**Version:** 1.0  
**Date:** 2026-02-07  
**Scope:** Complete replacement of the legacy HSD parser, xmlf90 DOM layer, and XML
output with hsd-fortran + hsd-data. No compatibility layer.

---

## Table of Contents

1. [Upstream Prerequisites (hsd-fortran & hsd-data)](#1-upstream-prerequisites)
2. [Architectural Overview](#2-architectural-overview)
3. [Migration Phases](#3-migration-phases)
   - Phase 0: Preparation & dependency wiring
   - Phase 1: Adapter shim layer
   - Phase 2: Core IO module replacement
   - Phase 3: Main parser conversion
   - Phase 4: Output replacement (detailed.xml → hsd-data)
   - Phase 5: Auxiliary application conversion
   - Phase 6: Public API conversion
   - Phase 7: Legacy removal & cleanup
4. [Unit Conversion Strategy](#4-unit-conversion-strategy)
5. [Error Handling Strategy](#5-error-handling-strategy)
6. [Backward Compatibility (oldcompat.F90)](#6-backward-compatibility)
7. [Testing Strategy](#7-testing-strategy)
8. [File-by-File Impact Analysis](#8-file-by-file-impact-analysis)
9. [Risk Assessment](#9-risk-assessment)

---

## 1. Upstream Prerequisites

Before any work begins in DFTB+, the following changes are required in the
upstream libraries.

### 1.1 hsd-fortran Changes Required

| # | Change | Rationale |
|---|--------|-----------|
| F1 | **`hsd_get_with_unit` — array & matrix overloads** | The legacy `convertUnitHsd` supports rank-0, rank-1, and rank-2 conversions. `hsd_get_with_unit` currently only supports scalars (`real(dp)` out). Add: `hsd_get_array_with_unit(table, path, val(:), target_unit, converter, stat)` and `hsd_get_matrix_with_unit(table, path, val(:,:), nrows, ncols, target_unit, converter, stat)`. |
| F2 | **`hsd_set_string_array`** | Already on hsd-fortran TODO. Needed for round-trip of string-list values (e.g., `TypeNames`, `SpinPerAtom` labels). |
| F3 | **`hsd_set_matrix` (integer & real)** | Already on hsd-fortran TODO. Needed for writing back matrix defaults in the processed output (`dftb_pin.hsd`). |
| F4 | **`hsd_set_attrib`** | Already on hsd-fortran TODO. The legacy parser injects/modifies modifier attributes on nodes (e.g., removing a consumed modifier). |
| F5 | **Polymorphic child retrieval** | The legacy parser heavily uses a pattern where a node has a single table child whose *name* is the selector (e.g., `Driver = ConjugateGradient { ... }`). hsd-fortran already supports this via `hsd_get_child` returning the named child. However, a convenience helper `hsd_get_choice(table, path, name, child_table, stat)` that returns both the child's name (as an output string) and a typed `hsd_table` pointer would reduce boilerplate in every `select case` dispatch throughout the 8500-line parser. |
| F6 | **Iterator with name filter** | Legacy code uses `getChildren(node, name, list)` to get all children with a given name (e.g., multiple `Atom` blocks). hsd-fortran's iterator covers all children; add a name-filtered variant: `hsd_get_children(table, path/name, children_array, stat)` returning an array of `hsd_node_ptr`. |
| F7 | **`hsd_set_or_default` pattern** | The legacy `getChildValue(node, name, val, default)` both *reads* the value and *injects the default into the tree* if absent (so the processed output contains all defaults). hsd-fortran's `hsd_get_or` returns the default to the caller but does *not* write it back to the tree. Add `hsd_get_or_set(table, path, val, default, stat)` that (a) returns the value, (b) if absent creates the node with the default value. This is critical for `dftb_pin.hsd` generation. |
| F8 | **Path normalization** | Already on hsd-fortran TODO. Trailing slashes, double slashes should be handled gracefully. |
| F9 | Fix `collect_unknown_fields` hardcoded 256-error limit | Already on hsd-fortran TODO (hsd_schema.f90). |
| F10 | Fix `schema_add_field_enum` silent truncation at 64 chars | Already on hsd-fortran TODO. |
| F11 | **Processed-node tracking / unprocessed-node detection** | The legacy system marks every accessed node with a `"proc"` attribute, then walks the tree to find unmarked nodes (→ warnings about unrecognized keywords). hsd-fortran has no equivalent. Options: (a) Add a `hsd_mark_processed(table, path)` + `hsd_get_unprocessed(table, paths)` pair, or (b) Rely on `schema_validate_strict` to detect unknown fields. Option (b) is preferred since schema validation is more principled, but requires that every section of the parser defines a schema. |
| F12 | **`hsd_rename_child(table, old_name, new_name, stat)`** | Used by backward-compatibility transformations (oldcompat.F90) and British/American spelling normalization (`localiseName`). |

### 1.2 hsd-data Changes Required

| # | Change | Rationale |
|---|--------|-----------|
| D1 | **`root_name` parameter on `data_load`** | Already on hsd-data TODO. Verify the document root tag matches expectations (e.g., `"dftbplusinput"`). |
| D2 | **Complete TOML backend** | Already on hsd-data TODO. While not strictly blocking migration, providing TOML as an alternative input format is a project goal. |
| D3 | **Ensure XML backend preserves nested attributes** | The legacy `detailed.xml` uses the xmlf90 writer which emits standard XML. The hsd-data XML writer must faithfully round-trip all attributes (including `unit`, custom attributes). Verify with a cross-format test using the exact `detailed.xml` schema. |
| D4 | **`data_dump` indentation control** | The legacy `writeDetailedXml` uses `indent=.true.` with xmlf90. hsd-data's `pretty` flag should produce equivalent formatting. Verify indent is 2-space (matching hsd-fortran convention). |
| D5 | **Remaining fixture round-trip tests** | Already on hsd-data TODO. These must all pass before DFTB+ can rely on cross-format fidelity. |

---

## 2. Architectural Overview

### 2.1 Current Architecture (Legacy)

```
dftb_in.hsd  ──[hsdparser.parseHSD()]──→  fnode* DOM tree (xmlf90)
                                                │
                          oldcompat.convertOldHSD()
                                                │
                          parser.parseHsdTree()  ─→  TInputData
                          │ uses: hsdutils.getChildValue()
                          │ uses: hsdutils2.convertUnitHsd()
                          │ uses: xmlutils (DOM navigation)
                          │
                          ├──→  dftb_pin.hsd  (via hsdparser.dumpHSD())
                          └──→  detailed.xml  (via xmlf90 wxml writer)
```

**29 source files** across the codebase directly depend on `fnode` / xmlf90.

### 2.2 Target Architecture

```
dftb_in.{hsd,json,xml,toml}  ──[data_load()]──→  hsd_table tree (hsd-fortran)
                                                       │
                              oldcompat.convertOldTree()
                                                       │
                              parser.parseTree()  ─→  TInputData
                              │ uses: hsd_get / hsd_get_or_set
                              │ uses: hsd_get_with_unit()
                              │ uses: hsd_get_choice() for dispatch
                              │ uses: hsd_schema for validation
                              │
                              ├──→  dftb_pin.hsd  (via data_dump())
                              └──→  detailed.{xml,json,hsd,hdf5}  (via data_dump())
```

**Zero** source files depend on xmlf90. The `external/xmlf90` dependency is
removed entirely.

---

## 3. Migration Phases

### Phase 0: Preparation & Dependency Wiring

**Goal:** Make hsd-data available as a build dependency; establish
infrastructure.

| Step | Description |
|------|-------------|
| 0.1 | Add `hsd-data` as a hybrid dependency (Submodule/Find/Fetch) in the top-level `CMakeLists.txt`, analogous to existing deps. Since hsd-data auto-fetches hsd-fortran, only hsd-data needs explicit wiring. |
| 0.2 | Add a new wrapper module `src/dftbp/extlibs/hsddata.F90` that does `use hsd_data` and re-exports the public API under the DFTB+ module naming convention (`dftbp_extlibs_hsddata`). This mirrors the current `dftbp_extlibs_xmlf90` pattern. |
| 0.3 | Add the new module to `src/dftbp/extlibs/CMakeLists.txt`. |
| 0.4 | Verify build: the library compiles with both old and new dependencies present. |
| 0.5 | Add a unit test under `test/src/dftbp/unit/io/` that loads a simple HSD string via `hsd_load_string` and retrieves a value — proving the wiring works end-to-end. |

### Phase 1: Adapter Shim Layer

**Goal:** Create a thin compatibility mapping between the legacy API calls and
hsd-fortran/hsd-data calls, to enable incremental conversion of the 8500-line
parser without a big-bang rewrite.

| Step | Description |
|------|-------------|
| 1.1 | Create `src/dftbp/io/hsdcompat.F90` — a temporary bridge module that provides the legacy call signatures (`getChildValue`, `getChild`, `setChildValue`, `detailedError`, `convertUnitHsd`, `getNodeName`) but implemented on top of `hsd_table`/`hsd_value`. |
| 1.2 | The shim types: `type(fnode_shim)` wrapping `type(hsd_table), pointer` or `class(hsd_node), pointer`. The shim `getChildValue(node_shim, name, val, default, modifier, child)` delegates to `hsd_get_or_set`. The shim `convertUnitHsd(modifier, units, child, val)` delegates to the DFTB+ unit converter (which remains unchanged — it operates on `TUnit` arrays and `real(dp)` values). |
| 1.3 | The shim `detailedError(node_shim, msg)` constructs an error message using `node%name` and `node%line` (from `hsd_node`) rather than the legacy xmlf90 attribute-based path/line extraction. |
| 1.4 | The shim is a *transitional* artifact. Each parser subroutine converted in Phase 3 switches from the shim to direct hsd-fortran calls, and the shim is removed when all conversions are complete (Phase 7). |

### Phase 2: Core IO Module Replacement

**Goal:** Replace the bottom layers: the HSD parser, HSD dumper, XML DOM, and
associated utilities.

| Step | Old Module | Action |
|------|-----------|--------|
| 2.1 | `dftbp_io_hsdparser` (921 lines) | **Remove entirely.** Replace all `parseHSD()` calls with `data_load()` and all `dumpHSD()` calls with `data_dump()`. `getNodeHSDName`, `getHSDPath` → derive from `hsd_node%name` / path-based navigation. `attrModifier`, `attrName`, etc. constants → no longer needed (hsd-fortran stores these as `hsd_node%attrib`, `hsd_node%name`). |
| 2.2 | `dftbp_io_hsdutils` (3651 lines) | **Remove entirely.** Every `getChildValue` call site → `hsd_get` / `hsd_get_or_set`. Every `setChildValue` call site → `hsd_set`. Every `writeChildValue` call site → `hsd_set` on an output tree then `data_dump`. `getChild` → `hsd_get_child`. `detailedError`/`detailedWarning` → new error reporting based on `hsd_error_t`. |
| 2.3 | `dftbp_io_hsdutils2` (467 lines) | **Remove entirely.** `convertUnitHsd` → `hsd_get_with_unit` (with DFTB+ converter function). `getUnprocessedNodes` / `warnUnprocessedNodes` → `schema_validate_strict` or a dedicated walk. `readHSDAsXML` → `data_load(filename, tree, fmt=DATA_FMT_XML)`. `setNodeName` / `renameChildren` / `localiseName` → `hsd_rename_child` (F12). `getDescendant` → `hsd_get_child`. |
| 2.4 | `dftbp_io_xmlutils` (265 lines) | **Remove entirely.** All DOM navigation → hsd-fortran query API. |
| 2.5 | `dftbp_extlibs_xmlf90` (24 lines) | **Remove entirely.** No more xmlf90 re-exports. |
| 2.6 | `dftbp_io_tokenreader` (615 lines) | **Evaluate.** This module is used to parse whitespace-separated tokens from text strings (e.g., reading coordinates, SK file data). It uses the xmlf90 `string` type. If the only xmlf90 dependency is the `string` type, convert to use standard Fortran `character(len=:), allocatable`. The core tokenization logic is independent of HSD and should be retained (possibly with interface changes). |
| 2.7 | `dftbp_io_indexselection` (509 lines) | **Keep but de-couple.** Currently uses `fnode` for error reporting. Change to accept a file/line context string instead. |
| 2.8 | `dftbp_io_linereader` (120 lines) | **Remove.** Only used by the legacy HSD parser. |

### Phase 3: Main Parser Conversion

**Goal:** Convert the 8500-line `parser.F90` and its sub-modules from
fnode/xmlf90 calls to hsd-fortran calls. This is the largest single task.

#### 3.1 Strategy

The parser consists of ~100 internal subroutines organized by input section.
Each subroutine follows the same pattern:

```fortran
! Legacy pattern (repeated hundreds of times):
call getChild(node, "SectionName", child, requested=.false.)
call getChildValue(node, "KeyName", value, defaultValue, modifier=modifier, child=field)
call convertUnitHsd(char(modifier), energyUnits, field, value)
```

The target pattern:

```fortran
! New pattern:
call hsd_get_or_set(table, "KeyName", value, defaultValue, stat=stat)
call hsd_get_with_unit(table, "KeyName", value, "au", dftbp_unit_converter, stat)
```

Or for polymorphic dispatch blocks:

```fortran
! Legacy:
call getChildValue(node, "", child)
call getNodeName(child, buffer)
select case (char(buffer))
  case ("conjugategradient")
    ...
end select

! New:
call hsd_get_choice(table, "", choice_name, choice_table, stat)
select case (choice_name)
  case ("conjugategradient")
    ...
end select
```

#### 3.2 Conversion Order

Convert section-by-section, each as an independently testable unit:

| Order | Subroutine(s) | Lines | Notes |
|-------|--------------|-------|-------|
| 1 | `readParserOptions` | ~60 | Simple, good warm-up |
| 2 | `handleInputVersion` | ~30 | Reads `InputVersion`, triggers oldcompat |
| 3 | `readGeometry` | ~400 | Medium complexity; tests coordinate parsing |
| 4 | `readOptions` | ~100 | Simple flags |
| 5 | `readAnalysis` | ~660 | Medium; PDOS/eigenvector options |
| 6 | `readFilling` | ~80 | Small; unit conversions |
| 7 | `readKPoints` | ~750 | Medium; matrix values, k-point grids |
| 8 | `readSolver` | ~200 | Medium; polymorphic dispatch |
| 9 | `readSKFiles` | ~500 | Complex; species-indexed maps |
| 10 | `readSpinPolarisation` | ~250 | Medium |
| 11 | `readDispersion` | ~1000 | Complex; many sub-blocks |
| 12 | `readDriver` | ~700 | Complex; many optimizer types |
| 13 | `readDFTBHam` | ~550 | Core Hamiltonian — most important |
| 14 | `readXTBHam` | ~600 | Similar to DFTB Ham |
| 15 | `readElecDynamics` | ~300 | Medium |
| 16 | `readTransportGeometry` | ~1800 | Large; NEGF-specific |
| 17 | `readReks` | ~200 | Small |
| 18 | `parseHybridBlock` | ~200 | Small |
| 19 | `parseChimes` | ~200 | Small |
| 20 | `parseHsdTree` (top-level) | ~100 | Orchestrator; convert last |

#### 3.3 Sub-Module Parsers

These files also need conversion:

| File | Lines | Notes |
|------|-------|-------|
| `src/dftbp/dftbplus/input/fileaccess.F90` | 62 | Trivial |
| `src/dftbp/dftbplus/input/geoopt.F90` | 259 | Medium; unit conversions |
| `src/dftbp/solvation/solvparser.F90` | 627 | Medium; heavy unit conversion |
| `src/dftbp/solvation/gbsafile.F90` | ~100 | Small |
| `src/dftbp/dftb/elecconstraints.F90` | ~200 | Medium |
| `src/dftbp/type/typegeometryhsd.F90` | ~200 | Medium |
| `src/dftbp/dftbplus/specieslist.F90` | 130 | Medium |

#### 3.4 inputdata.F90 (TInputData)

`TInputData` and its component types (`TControl`, `TGeometry`, `TSlater`,
`TTransPar`) are **unchanged** by this migration. They are pure Fortran data
types with no HSD/XML dependency. The migration only changes *how* they are
populated, not *what* they contain.

### Phase 4: Output Replacement

**Goal:** Replace the legacy xmlf90-based output with hsd-data.

#### 4.1 `detailed.xml` → Multi-Format Output

The current `writeDetailedXml()` in `mainio.F90` (called from `main.F90`) uses
the xmlf90 XML writer (`xmlf_t`) to emit `detailed.xml`. This is replaced with:

1. Build an `hsd_table` tree representing the detailed output data:
   ```fortran
   type(hsd_table) :: detailed
   call new_table(detailed, "detailedout")
   call hsd_set(detailed, "identity", runId)
   call hsd_set(detailed, "geometry/typenames", typeNames)
   call hsd_set(detailed, "geometry/periodic", isPeriodic)
   ! ... etc.
   ```

2. Dump using hsd-data:
   ```fortran
   call data_dump(detailed, "detailed.xml", error)
   ! Or, based on user option:
   call data_dump(detailed, "detailed.json", error)
   call data_dump(detailed, "detailed.hsd", error)
   ```

3. The output format becomes configurable via an input keyword
   (`Options/DetailedOutputFormat`), defaulting to `"xml"` for backward
   compatibility during transition, switching to `"hsd"` in a future release.

#### 4.2 `dftb_pin.hsd` (Processed Input)

Currently generated by `dumpHSD(hsdTree, "dftb_pin.hsd")` after the DOM tree
has been enriched with defaults. In the new architecture:

- The `hsd_get_or_set` calls (F7) ensure defaults are written back to the tree.
- After parsing completes: `call data_dump(root, "dftb_pin.hsd", error)`.
- Optionally also: `call data_dump(root, "dftb_pin.json", error)`.

#### 4.3 `autotest.tag` (Tagged Output)

The `TTaggedWriter` system (`taggedoutput.F90`) is independent of HSD/xmlf90
and can remain as-is. It writes a flat text format for regression testing.
Future work (not part of this migration) could replace it with HDF5 via
hsd-data.

#### 4.4 Waveplot's `detailed.xml` Reading

`app/waveplot/initwaveplot.F90` reads `detailed.xml` using
`readHSDAsXML(filename, fnode_tree)` then navigates the fnode tree. Replace
with:

```fortran
call data_load("detailed.xml", detailed, error)
call hsd_get(detailed, "geometry/typenames", typeNames, stat)
! ... etc.
```

This is cleaner and format-agnostic — waveplot would accept
`detailed.{xml,json,hsd}` automatically.

### Phase 5: Auxiliary Application Conversion

**Goal:** Convert all auxiliary apps that use the legacy HSD API.

| Application | File | Lines | Effort |
|-------------|------|-------|--------|
| modes | `app/modes/initmodes.F90` | ~800 | Medium |
| phonons | `app/phonons/initphonons.F90` | ~600 | Medium |
| transporttools | `app/transporttools/parser.F90` | ~800 | Medium |
| waveplot | `app/waveplot/initwaveplot.F90` | ~900 | Medium (+ detailed.xml reading) |
| skderivs | `app/misc/skderivs/skderivs.F90` | ~200 | Small |

Each follows the same conversion pattern as Phase 3: replace `parseHSD` with
`data_load`, replace `getChildValue`/`getChild` with `hsd_get`/`hsd_get_child`,
replace `dumpHSD` with `data_dump`, replace `destroyNode` with
`table%destroy()`.

### Phase 6: Public API Conversion

**Goal:** Convert the public DFTB+ API layer.

#### 6.1 `hsdapi.F90`

Currently re-exports `fnode`, `fnodeList`, `getChild`, `getChildValue`,
`setChildValue`, `dumpHsd`. Replace with re-exports from hsd-data:

```fortran
module dftbp_hsdapi
  use hsd_data
  implicit none
  public
end module dftbp_hsdapi
```

**This is a breaking API change.** External users who build HSD trees
programmatically must switch from `fnode` pointer manipulation to `hsd_table`
operations. Since this is a major version bump with no compatibility layer, this
is expected.

#### 6.2 `mmapi.F90`

The `TDftbPlusInput` type stores `hsdTree` as `type(fnode), pointer`. Change to
`type(hsd_table)`:

- `getInputFromFile(fileName, input)` → `call data_load(fileName, input%hsdTree, error)`
- `getEmptyInput(input)` → `call new_table(input%hsdTree, "dftbplusinput")`
- `getRootNode(root)` → return `input%hsdTree` directly
- `setupCalculator(input)` → `call parseTree(input%hsdTree, inpData, parserFlags)`

### Phase 7: Legacy Removal & Cleanup

**Goal:** Remove all legacy code and the xmlf90 dependency.

| Step | Action |
|------|--------|
| 7.1 | Delete `src/dftbp/extlibs/xmlf90.F90` |
| 7.2 | Delete `src/dftbp/io/hsdparser.F90` |
| 7.3 | Delete `src/dftbp/io/hsdutils.F90` |
| 7.4 | Delete `src/dftbp/io/hsdutils2.F90` |
| 7.5 | Delete `src/dftbp/io/xmlutils.F90` |
| 7.6 | Delete `src/dftbp/io/linereader.F90` |
| 7.7 | Delete `src/dftbp/io/hsdcompat.F90` (transitional shim from Phase 1) |
| 7.8 | Remove `external/xmlf90/` from submodules and CMake |
| 7.9 | Remove all xmlf90-related CMake wiring (`find_package`, submodule, FetchContent`) |
| 7.10 | Remove `dftbp_io_charmanip` dependencies on xmlf90 `string` type (if any remain) |
| 7.11 | Update `src/dftbp/common/filesystem.F90` and `src/dftbp/type/linkedlists0.F90` to remove xmlf90 `string` type usage (replace with `character(len=:), allocatable`) |
| 7.12 | Verify `fortitude check` / linting passes (if adopted) |
| 7.13 | Full regression test pass (`ctest --test-dir build`) |

---

## 4. Unit Conversion Strategy

### 4.1 Current Approach

The legacy system stores unit conversion tables as arrays of `TUnit` (from
`unitconversion.F90`), each containing a name and a multiplicative conversion
factor. `convertUnitHsd` looks up the modifier string in the unit table and
multiplies the value.

### 4.2 Target Approach

`hsd_get_with_unit` (from hsd-fortran `hsd_validation`) accepts a
user-provided `converter` function. DFTB+ provides a single converter
function that wraps the existing `TUnit`-based system:

```fortran
module dftbp_io_unitbridge
  use dftbp_common_unitconversion, only: TUnit, convertUnit, &
      & lengthUnits, energyUnits, timeUnits, ... etc.
  use hsd_data, only: dp
  implicit none
  private
  public :: dftbp_convert_length, dftbp_convert_energy, ... etc.

contains

  pure function dftbp_convert_length(value, from_unit, to_unit) result(res)
    real(dp), intent(in) :: value
    character(len=*), intent(in) :: from_unit, to_unit
    real(dp) :: res
    ! Convert from_unit → atomic units → to_unit using lengthUnits table
    res = value * get_factor(lengthUnits, from_unit) / get_factor(lengthUnits, to_unit)
  end function

end module
```

This preserves all existing unit definitions while adapting to hsd-fortran's
converter interface. The `unitconversion.F90` module itself is unchanged.

### 4.3 Usage

```fortran
! Legacy:
call getChildValue(node, "Temperature", temp, 0.0_dp, modifier=mod, child=ch)
call convertUnitHsd(char(mod), energyUnits, ch, temp)

! New:
call hsd_get_with_unit(table, "Temperature", temp, "au", dftbp_convert_energy, stat)
```

For values with defaults that also need unit conversion, the pattern is:

```fortran
call hsd_get_or_set(table, "Temperature", temp, 0.0_dp, stat=stat)
if (stat == HSD_STAT_OK) then
  call hsd_get_with_unit(table, "Temperature", temp, "au", dftbp_convert_energy, stat)
end if
```

---

## 5. Error Handling Strategy

### 5.1 Current Approach

The legacy parser uses `detailedError(fnode, msg)` which:
1. Extracts the HSD source path via `getHSDPath(fnode)`
2. Extracts the source line via `getAttribute(fnode, "start")`
3. Calls `error()` which prints and aborts

### 5.2 Target Approach

Replace with a centralized error reporting module:

```fortran
module dftbp_io_hsderror
  use hsd_data, only: hsd_error_t, hsd_node, hsd_table, hsd_value
  implicit none
  private
  public :: hsd_fatal_error, hsd_warning

contains

  subroutine hsd_fatal_error(node, msg)
    class(hsd_node), intent(in) :: node
    character(len=*), intent(in) :: msg
    ! Construct error with node%name, node%line
    ! Call error() from dftbp_io_message
  end subroutine

  subroutine hsd_warning(node, msg)
    class(hsd_node), intent(in) :: node
    character(len=*), intent(in) :: msg
    ! Construct warning with node%name, node%line
    ! Call warning() from dftbp_io_message
  end subroutine

end module
```

The `stat`-based error handling from hsd-fortran is used for expected
conditions (missing optional keys, type mismatches that are recoverable).
Fatal errors (invalid input that prevents continuation) still call the DFTB+
`error()` routine directly.

---

## 6. Backward Compatibility (oldcompat.F90)

### 6.1 Current Approach

`oldcompat.F90` (997 lines) contains 14 DOM tree transformation subroutines
(`convert_1_2` through `convert_13_14`). Each manipulates the `fnode` DOM tree:
renaming nodes, moving subtrees, changing attribute values, deleting/creating
nodes.

### 6.2 Target Approach

Rewrite `oldcompat.F90` to operate on `hsd_table` using hsd-fortran operations:

| Legacy Operation | hsd-fortran Equivalent |
|-----------------|----------------------|
| `setNodeName(child, newName)` | `hsd_rename_child(table, oldName, newName)` (F12) |
| `renameChildren(node, old, new)` | Loop + `hsd_rename_child` |
| `localiseName(node, local, english)` | `hsd_rename_child` |
| `getChild(node, name, child)` | `hsd_get_child(table, name, child)` |
| `appendChild(node, newChild)` | `table%add_child(newChild)` |
| `removeChild(node, child)` | `hsd_remove_child(table, name)` |
| `setAttribute(node, "name", val)` | `hsd_set_attrib(table, path, val)` (F4) |
| `setUnprocessed(node)` | Not needed with schema-based validation |

The conversion subroutines themselves are straightforward 1:1 translations.
The `versionMaps` and version-checking logic remain unchanged (pure data).

---

## 7. Testing Strategy

### 7.1 Regression Testing

All ~400 existing application regression tests (`test/app/dftb+/`) must
continue to pass. These tests run full DFTB+ calculations and compare output
against reference data. They are the primary correctness gate.

### 7.2 Processed Output Comparison

After Phase 3, `dftb_pin.hsd` output from the new parser must be byte-for-byte
identical (or semantically identical) to the legacy output. Create a comparison
script that:
1. Runs the legacy code on a test input → `dftb_pin.hsd.legacy`
2. Runs the new code on the same input → `dftb_pin.hsd.new`
3. Compares (allowing whitespace normalization)

### 7.3 Cross-Format Output Testing

After Phase 4, verify that `detailed.xml` produced by hsd-data's XML writer
is semantically equivalent to the legacy xmlf90-produced version. Parse both
into `hsd_table` trees and compare with `hsd_table_equal`.

### 7.4 New Unit Tests

Add unit tests for:

| Test | Location |
|------|----------|
| HSD load/dump round-trip | `test/src/dftbp/unit/io/` |
| Unit conversion bridge functions | `test/src/dftbp/unit/io/` |
| Error reporting (line numbers, paths) | `test/src/dftbp/unit/io/` |
| Individual parser sections (readGeometry, readHamiltonian, ...) | `test/src/dftbp/unit/dftbplus/` |
| Old compatibility transformations | `test/src/dftbp/unit/dftbplus/` |
| Multi-format input acceptance (HSD, JSON, XML) | `test/app/dftb+/` |

### 7.5 API Tests

Verify that the public API (`mmapi.F90`, `hsdapi.F90`) works with
`hsd_table` types. Update the existing API test suite in `test/src/` and
`test/tools/pythonapi/`.

---

## 8. File-by-File Impact Analysis

### Files to DELETE (Legacy IO, ~5940 lines)

| File | Lines | Reason |
|------|-------|--------|
| `src/dftbp/extlibs/xmlf90.F90` | 24 | xmlf90 wrapper |
| `src/dftbp/io/hsdparser.F90` | 921 | Legacy HSD parser |
| `src/dftbp/io/hsdutils.F90` | 3,651 | Legacy value extraction |
| `src/dftbp/io/hsdutils2.F90` | 467 | Legacy HSD utilities |
| `src/dftbp/io/xmlutils.F90` | 265 | Legacy DOM navigation |
| `src/dftbp/io/linereader.F90` | 120 | Only used by legacy parser |
| `external/xmlf90/` (all) | ~6,400 | External library |

### Files to CREATE (~300 lines)

| File | Lines (est.) | Purpose |
|------|-------------|---------|
| `src/dftbp/extlibs/hsddata.F90` | 20 | hsd-data wrapper module |
| `src/dftbp/io/unitbridge.F90` | 150 | Unit converter bridge functions |
| `src/dftbp/io/hsderror.F90` | 80 | Error reporting bridge |
| `src/dftbp/io/hsdcompat.F90` | 200 | Transitional shim (deleted in Phase 7) |

### Files to HEAVILY MODIFY (core conversion, ~21,000 lines affected)

| File | Lines | Change Description |
|------|-------|-------------------|
| `src/dftbp/dftbplus/parser.F90` | 8,528 | Every `fnode`→`hsd_table`, every `getChildValue`→`hsd_get`, etc. |
| `src/dftbp/dftbplus/mainio.F90` | 6,119 | `writeDetailedXml`→ build `hsd_table` + `data_dump`; remove xmlf90 writer |
| `src/dftbp/dftbplus/oldcompat.F90` | 997 | All DOM manipulation → hsd-fortran operations |
| `src/dftbp/dftbplus/hsdhelpers.F90` | 97 | `parseHSD`→`data_load`, `dumpHSD`→`data_dump` |
| `src/dftbp/solvation/solvparser.F90` | 627 | All fnode calls → hsd-fortran |
| `src/dftbp/io/tokenreader.F90` | 615 | Remove xmlf90 `string` type dependency |
| `src/dftbp/io/indexselection.F90` | 509 | Remove fnode dependency in error reporting |
| `src/dftbp/dftbplus/specieslist.F90` | 130 | fnode → hsd_table |
| `src/dftbp/dftbplus/input/geoopt.F90` | 259 | fnode → hsd_table |
| `src/dftbp/dftbplus/input/fileaccess.F90` | 62 | fnode → hsd_table |
| `src/dftbp/type/typegeometryhsd.F90` | ~200 | fnode → hsd_table |
| `src/dftbp/dftb/elecconstraints.F90` | ~200 | fnode → hsd_table |
| `src/dftbp/solvation/gbsafile.F90` | ~100 | fnode → hsd_table |

### Files to LIGHTLY MODIFY (~200 lines affected)

| File | Change |
|------|--------|
| `src/dftbp/common/filesystem.F90` | Remove xmlf90 `string` type (use `character`) |
| `src/dftbp/type/linkedlists0.F90` | Remove xmlf90 `string` type (use `character`) |
| `src/dftbp/api/mm/mmapi.F90` | `fnode` pointer → `hsd_table` |
| `src/dftbp/api/mm/hsdapi.F90` | Re-export hsd-data instead of xmlf90 |
| `app/dftb+/dftbplus.F90` | Minimal; `parseHsdInput` signature unchanged |

### Auxiliary Apps to CONVERT (~3300 lines)

| File | Lines |
|------|-------|
| `app/modes/initmodes.F90` | ~800 |
| `app/phonons/initphonons.F90` | ~600 |
| `app/transporttools/parser.F90` | ~800 |
| `app/waveplot/initwaveplot.F90` | ~900 |
| `app/misc/skderivs/skderivs.F90` | ~200 |

---

## 9. Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| Subtle behavioral differences between legacy parser and hsd-fortran parser (whitespace handling, quoting, escape sequences) | Medium | High | Comprehensive regression testing; add edge-case HSD inputs to test suite |
| `dftb_pin.hsd` output differs from legacy (causing user confusion) | Medium | Low | Accept cosmetic differences; verify semantic equivalence |
| Performance regression in parsing (hsd-fortran vs legacy) | Low | Low | HSD parsing is <<1% of total runtime; hash table lookup in hsd-fortran is O(1) vs legacy O(n) |
| hsd-data XML output differs from legacy xmlf90 output (breaking waveplot/external tools that parse `detailed.xml`) | Medium | Medium | Validate XML output with XML diff tools; document format changes |
| `hsd_get_with_unit` converter interface mismatches with DFTB+ unit system | Low | Medium | Thorough unit conversion testing; bridge module isolates the interface |
| oldcompat transformations break under new tree representation | Medium | High | Unit-test each conversion function individually |
| Public API break affects downstream users | High (by design) | Medium | Document API changes; provide migration guide; major version bump |
| Phase 3 regression during incremental conversion (half-migrated parser) | Medium | High | Shim layer (Phase 1) enables compiling/testing after each section conversion |

---

## Appendix A: Legacy API → hsd-fortran Cheat Sheet

| Legacy Call | New Call |
|------------|---------|
| `parseHSD(rootTag, file, doc)` | `data_load(file, table, error)` |
| `dumpHSD(doc, file)` | `data_dump(table, file, error)` |
| `getChild(node, name, child)` | `hsd_get_child(table, name, child, stat)` |
| `getChild(node, name, child, requested=.false.)` | `hsd_get_child(table, name, child, stat)` + check `stat` |
| `getChildren(node, name, list)` | `hsd_get_children(table, name, array, stat)` (F6) |
| `getChildValue(node, name, val)` | `hsd_get(table, name, val, stat)` |
| `getChildValue(node, name, val, default)` | `hsd_get_or_set(table, name, val, default, stat)` (F7) |
| `getChildValue(node, name, val, default, modifier=m, child=c)` | `hsd_get_or_set(...)` then `hsd_get_attrib(table, name, m, stat)` |
| `setChildValue(node, name, val)` | `hsd_set(table, name, val)` |
| `writeChildValue(xf, name, val)` | `hsd_set(output_table, name, val)` |
| `convertUnitHsd(mod, units, child, val)` | `hsd_get_with_unit(table, path, val, "au", converter, stat)` |
| `detailedError(node, msg)` | `hsd_fatal_error(node, msg)` |
| `detailedWarning(node, msg)` | `hsd_warning(node, msg)` |
| `warnUnprocessedNodes(node, flag)` | `schema_validate_strict(schema, table, errors)` |
| `getNodeName(node, name)` → `char(name)` | `child%name` (direct field access) |
| `getAttribute(node, "m")` | `node%attrib` / `hsd_get_attrib(table, path, attrib)` |
| `setAttribute(node, "proc", "")` | Not needed (schema validation replaces) |
| `destroyNode(doc)` | `table%destroy()` |
| `readHSDAsXML(file, doc)` | `data_load(file, table, error, fmt=DATA_FMT_XML)` |
| `setNodeName(child, name)` | `hsd_rename_child(table, oldName, newName)` (F12) |
| `localiseName(node, uk, us)` | `hsd_rename_child(table, uk, us)` (F12) |

## Appendix B: Multi-Format Input Support

After migration, DFTB+ gains automatic support for multiple input formats:

```bash
# All equivalent:
dftb+ < dftb_in.hsd
dftb+ < dftb_in.json
dftb+ < dftb_in.xml
dftb+ < dftb_in.toml
```

The `hsdhelpers.F90` orchestrator detects the format by extension and calls
`data_load` accordingly. The internal representation is always `hsd_table`.
`dftb_pin.hsd` is always generated in HSD format (user-readable). The detailed
output format is configurable.
