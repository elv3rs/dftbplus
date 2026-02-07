!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Exporting the functionality we use from the hsd-data library.
!>
!> hsd-data provides unified multi-format IO (HSD, XML, JSON, TOML, HDF5) built on top of
!> hsd-fortran's in-memory tree representation. This module re-exports the public API under the
!> DFTB+ module naming convention.
module dftbp_extlibs_hsddata
  use hsd_data, only : &
      ! Precision kinds
      & dp, sp, &
      ! Error type and status codes
      & hsd_error_t, &
      & HSD_STAT_OK, HSD_STAT_SYNTAX_ERROR, HSD_STAT_UNCLOSED_TAG, &
      & HSD_STAT_UNCLOSED_ATTRIB, HSD_STAT_UNCLOSED_QUOTE, HSD_STAT_ORPHAN_TEXT, &
      & HSD_STAT_INCLUDE_CYCLE, HSD_STAT_INCLUDE_DEPTH, HSD_STAT_FILE_NOT_FOUND, &
      & HSD_STAT_IO_ERROR, HSD_STAT_TYPE_ERROR, HSD_STAT_NOT_FOUND, &
      & HSD_STAT_SCHEMA_ERROR, &
      ! Core tree types
      & hsd_node, hsd_table, hsd_value, hsd_node_ptr, hsd_iterator, &
      & new_table, new_value, &
      ! Value type constants
      & VALUE_TYPE_NONE, VALUE_TYPE_STRING, VALUE_TYPE_INTEGER, &
      & VALUE_TYPE_REAL, VALUE_TYPE_LOGICAL, VALUE_TYPE_ARRAY, VALUE_TYPE_COMPLEX, &
      ! Visitor pattern
      & hsd_visitor_t, hsd_accept, &
      ! Accessors / mutators / query
      & hsd_get, hsd_get_or, hsd_get_matrix, &
      & hsd_set, &
      & hsd_get_child, hsd_get_table, hsd_has_child, hsd_remove_child, &
      & hsd_get_type, hsd_is_table, hsd_is_value, hsd_is_array, &
      & hsd_child_count, hsd_get_keys, hsd_get_attrib, hsd_has_attrib, &
      & hsd_merge, hsd_clone, hsd_table_equal, &
      ! Validation / schema
      & hsd_require, hsd_validate_range, hsd_validate_one_of, hsd_get_with_unit, &
      & hsd_schema_t, hsd_field_def_t, &
      & FIELD_REQUIRED, FIELD_OPTIONAL, &
      & FIELD_TYPE_ANY, FIELD_TYPE_STRING, FIELD_TYPE_INTEGER, &
      & FIELD_TYPE_REAL, FIELD_TYPE_LOGICAL, FIELD_TYPE_TABLE, &
      & FIELD_TYPE_ARRAY, FIELD_TYPE_COMPLEX, &
      & schema_init, schema_destroy, schema_add_field, schema_add_field_enum, &
      & schema_validate, schema_validate_strict, &
      ! Format constants (hsd-data)
      & DATA_FMT_AUTO, DATA_FMT_HSD, DATA_FMT_XML, DATA_FMT_JSON, &
      & DATA_FMT_TOML, DATA_FMT_HDF5, &
      ! Unified multi-format IO (hsd-data)
      & data_load, data_load_string, &
      & data_dump, data_dump_to_string, &
      & data_detect_format, data_format_available, data_convert
  implicit none
  public

end module dftbp_extlibs_hsddata
