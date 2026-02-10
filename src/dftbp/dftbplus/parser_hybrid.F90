!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

!> Parses hybrid xc-functional and ChIMES related input.
module dftbp_dftbplus_parser_hybrid
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_filesystem, only : findFile, getParamSearchPaths
  use dftbp_common_unitconversion, only : lengthUnits
  use dftbp_dftb_hybridxc, only : hybridXcAlgo, hybridXcFunc, hybridXcGammaTypes
  use dftbp_dftb_repulsive_chimesrep, only : TChimesRepInp
  use dftbp_dftbplus_inputdata, only : TControl, THybridXcInp
  use dftbp_io_charmanip, only : newline, tolower, unquote
  use hsd, only : hsd_get_or_set, hsd_get_table, hsd_get_choice, hsd_get_attrib, hsd_get, &
      & HSD_STAT_OK, hsd_schema_t, hsd_error_t, schema_init, schema_add_field, &
      & schema_validate, schema_destroy, FIELD_OPTIONAL, FIELD_TYPE_TABLE, FIELD_TYPE_STRING
  use dftbp_io_hsdutils, only : dftbp_error, dftbp_warning, getNodeName, getNodeName2
  use dftbp_io_unitconv, only : convertUnitHsd
  use dftbp_io_message, only : error
  use hsd_data, only : hsd_table
  use dftbp_dftbplus_parser_skfiles, only : TCharLcArray
  use dftbp_type_oldskdata, only : parseHybridXcTag
  use dftbp_type_typegeometry, only : TGeometry
  implicit none

  private
  public :: parseHybridBlock, parseChimes

contains


  !> Parses hybrid xc-functional input.
  subroutine parseHybridBlock(node, input, ctrl, geo, skFiles)

    !> Node to parse
    type(hsd_table), intent(in), pointer :: node

    !> Range separated data structure to fill
    type(THybridXcInp), intent(inout), allocatable :: input

    !> Geometry structure
    type(TGeometry), intent(in) :: geo

    !> General control structure
    type(TControl), intent(in) :: ctrl

    !> List of SK file names to read in for every interaction
    type(TCharLcArray), intent(inout) :: skFiles(:,:)

    !! File name of representative SK-file to read
    character(lc) :: fileName

    !! True, if hybrid xc-functional input block present
    logical :: isHybridInp

    !! True, if hybrid xc-functional extra tag found in SK-file(s)
    logical :: isHybridSk

    !! Hybrid xc-functional extra tag found in the SK-file
    integer :: hybridXcSkTag

    !! Type of hybrid xc-functional found in the SK-file
    integer :: hybridXcSkType

    !! Hybrid functional type of user input
    integer :: hybridXcInputTag

    !! Auxiliary node pointers
    type(hsd_table), pointer :: hybridChild, hybridValue, screeningChild
    type(hsd_table), pointer :: screeningValue, cmChild, cmValue, child1, child2

    !! Temporary string buffers
    character(len=:), allocatable :: buffer, modifier

    !! Temporary string buffer, that stores the gamma function type
    character(len=:), allocatable :: strBuffer
    integer :: stat

    @:ASSERT(size(skFiles, dim=1) == size(skFiles, dim=2))
    @:ASSERT((size(skFiles, dim=1) > 0))

    ! Extracting hybridXc tag from first SK-file only is a workaround and assumes that a set of
    ! given SK-files uses the same parameters (which should always be the case)!
    fileName = skFiles(1, 1)%items(1)

    ! Check if SK-files contain extra tag for hybrid xc-functionals
    call parseHybridXcTag(fileName, hybridXcTag=hybridXcSkTag, hybridXcType=hybridXcSkType)
    isHybridSk = hybridXcSkTag /= hybridXcFunc%none

    call hsd_get_table(node, "Hybrid", hybridChild, stat, auto_wrap=.true.)
    if (associated(hybridChild)) then
      call hsd_get_choice(hybridChild, "", buffer, hybridValue, stat)
      if (.not. associated(hybridValue) .and. len_trim(buffer) == 0) buffer = "none"
    else
      buffer = "none"
      hybridValue => null()
    end if

    isHybridInp = associated(hybridChild) .and. (tolower(buffer) /= "none")

    if (isHybridInp .and. .not. isHybridSk) then
      call error("Hybrid input block present, but SK-file '" // trim(fileName)&
          & // "'" // newline // "   appears to be (semi-)local.")
    elseif (isHybridSk .and. .not. isHybridInp) then
      call error("Hybrid SK-file '" // trim(fileName) // "' provided," // newline //&
          & "   but the hybrid block is missing from the HSD input file.")
    end if

    hybridXcInputTag = hybridXcFunc%none
    if (isHybridInp) then
      ! Convert hybrid functional type of user input to enumerator
      select case(tolower(buffer))
      case ("lc")
        hybridXcInputTag = hybridXcFunc%lc
      case ("cam")
        hybridXcInputTag = hybridXcFunc%cam
      case default
        call dftbp_error(hybridChild, "Unknown hybrid xc-functional type '" // buffer&
            & // "' in input.")
      end select

      ! Check if hybrid functional type is in line with SK-files
      if (hybridXcInputTag == hybridXcFunc%lc .and. hybridXcSkType /= hybridXcFunc%lc) then
        call dftbp_error(hybridChild, "Requested hybrid functional type conflict with provided&
            & SK-file(s).")
      end if

      allocate(input)
      input%hybridXcType = hybridXcSkType
      call hsd_get_table(hybridValue, "Screening", screeningChild, stat, auto_wrap=.true.)
      if (associated(screeningChild)) then
        call hsd_get_choice(screeningChild, "", buffer, screeningValue, stat)
        if (.not. associated(screeningValue) .and. len_trim(buffer) == 0) buffer = "matrixbased"
      else
        buffer = "matrixbased"
        screeningValue => null()
        screeningChild => null()
      end if

      select case(tolower(buffer))
      case ("neighbourbased")
        input%hybridXcAlg = hybridXcAlgo%neighbourBased
        call hsd_get_or_set(screeningValue, "CutoffReduction", input%cutoffRed, 0.0_dp)
        call hsd_get_attrib(screeningValue, "CutoffReduction", modifier, stat)
        if (stat /= HSD_STAT_OK) modifier = ""
        call hsd_get_table(screeningValue, "CutoffReduction", child1, stat, auto_wrap=.true.)
        if (.not. associated(child1)) child1 => screeningValue
        call convertUnitHsd(modifier, lengthUnits, child1, input%cutoffRed)
        if (geo%tPeriodic) then
          call hsd_get_or_set(screeningValue, "Threshold", input%screeningThreshold, 1e-6_dp)
        end if
      case ("thresholded")
        input%hybridXcAlg = hybridXcAlgo%thresholdBased
        call hsd_get_or_set(screeningValue, "Threshold", input%screeningThreshold, 1e-6_dp)
        call hsd_get_or_set(screeningValue, "CutoffReduction", input%cutoffRed, 0.0_dp)
        call hsd_get_attrib(screeningValue, "CutoffReduction", modifier, stat)
        if (stat /= HSD_STAT_OK) modifier = ""
        call hsd_get_table(screeningValue, "CutoffReduction", child1, stat, auto_wrap=.true.)
        if (.not. associated(child1)) child1 => screeningValue
        call convertUnitHsd(modifier, lengthUnits, child1, input%cutoffRed)
      case ("matrixbased")
        input%hybridXcAlg = hybridXcAlgo%matrixBased
        ! In this case, CutoffRedunction is not used so it should be set to zero.
        input%cutoffRed = 0.0_dp
      case default
        call dftbp_error(screeningChild, "Invalid screening method '" // buffer // "'")
      end select

      if (ctrl%tSpinOrbit) then
        call dftbp_error(hybridChild, "Spin-orbit coupling not currently supported for hybrids")
      end if
      if (ctrl%t2Component) then
        if (input%hybridXcAlg /= hybridXcAlgo%matrixBased) then
          call dftbp_error(screeningChild, "MatrixBased screening required for noncollinear spin")
        end if
      end if

      ! Additional settings for periodic sytems
      ifPeriodic: if (geo%tPeriodic) then

        ! parse gamma function type (full, truncated, mic, ...)
        call hsd_get_table(hybridValue, "CoulombMatrix", cmChild, stat, auto_wrap=.true.)
        if (associated(cmChild)) then
          call hsd_get_choice(cmChild, "", buffer, cmValue, stat)
          if (.not. associated(cmValue) .and. len_trim(buffer) == 0) buffer = "truncated"
        else
          buffer = "truncated"
          cmValue => null()
          cmChild => null()
        end if

        select case(tolower(buffer))
        case ("full")
          input%gammaType = hybridXcGammaTypes%full
        case ("minimumimage")
          input%gammaType = hybridXcGammaTypes%mic
        case ("truncated")
          input%gammaType = hybridXcGammaTypes%truncated
        case ("truncatedanddamped")
          input%gammaType = hybridXcGammaTypes%truncatedAndDamped
        case default
          call dftbp_error(cmChild, "Invalid Gamma function type '" // buffer // "'")
        end select

        ! g-Summation cutoff not needed for MIC CAM Hamiltonian
        if (input%gammaType /= hybridXcGammaTypes%mic) then
          call hsd_get_table(cmValue, "GSummationCutoff", child1, stat, auto_wrap=.true.)
          if (associated(child1)) then
            call hsd_get_attrib(cmValue, "GSummationCutoff", modifier, stat)
            if (stat /= HSD_STAT_OK) modifier = ""
            allocate(input%gSummationCutoff)
            call hsd_get(child1, "#text", input%gSummationCutoff, stat=stat)
            if (stat /= HSD_STAT_OK) call dftbp_error(child1, "Missing GSummationCutoff value")
            call convertUnitHsd(modifier, lengthUnits, child1, input%gSummationCutoff)
          end if
        end if

        if (input%gammaType == hybridXcGammaTypes%truncated&
            & .or. input%gammaType == hybridXcGammaTypes%truncatedAndDamped) then
          call hsd_get_table(cmValue, "CoulombCutoff", child1, stat, auto_wrap=.true.)
          if (associated(child1)) then
            call hsd_get_attrib(cmValue, "CoulombCutoff", modifier, stat)
            if (stat /= HSD_STAT_OK) modifier = ""
            allocate(input%gammaCutoff)
            call hsd_get(child1, "#text", input%gammaCutoff, stat=stat)
            if (stat /= HSD_STAT_OK) call dftbp_error(child1, "Missing CoulombCutoff value")
            call convertUnitHsd(modifier, lengthUnits, child1, input%gammaCutoff)
          end if
        end if

      else
        ! Always use unaltered gamma function for non-periodic systems
        input%gammaType = hybridXcGammaTypes%full
      end if ifPeriodic

      ! Number of primitive cells regarded in MIC, along each supercell folding direction
      if (input%gammaType == hybridXcGammaTypes%mic) then
        allocate(input%wignerSeitzReduction)
        call hsd_get_or_set(cmValue, "WignerSeitzReduction", input%wignerSeitzReduction, 0)
      end if

    end if

    ! -- Schema validation (proof-of-concept, warnings only) --
    block
      type(hsd_schema_t) :: schema
      type(hsd_error_t), allocatable :: schemaErrors(:)
      integer :: iErr

      if (isHybridInp .and. associated(hybridValue)) then
        call schema_init(schema, name="Hybrid")
        call schema_add_field(schema, "Screening", FIELD_OPTIONAL, FIELD_TYPE_TABLE)
        call schema_add_field(schema, "CoulombMatrix", FIELD_OPTIONAL, FIELD_TYPE_TABLE)
        call schema_validate(schema, hybridValue, schemaErrors)
        if (size(schemaErrors) > 0) then
          do iErr = 1, size(schemaErrors)
            call dftbp_warning(hybridValue, "[schema] " // schemaErrors(iErr)%message)
          end do
        end if
        call schema_destroy(schema)
      end if
    end block

  end subroutine parseHybridBlock



  !> Parses Chimes related options.
  subroutine parseChimes(root, chimesRepInput)
    type(hsd_table), pointer, intent(in) :: root
    type(TChimesRepInp), allocatable, intent(out) :: chimesRepInput

    type(hsd_table), pointer :: chimes
    integer :: stat

  #:if WITH_CHIMES
    character(len=:), allocatable :: buffer
    character(len=:), allocatable :: searchPath(:)
    character(len=:), allocatable :: chimesFile
  #:endif

    call hsd_get_table(root, "Chimes", chimes, stat, auto_wrap=.true.)
    if (.not. associated(chimes)) return
    #:if WITH_CHIMES
      allocate(chimesRepInput)
      call hsd_get_or_set(chimes, "ParameterFile", buffer, "chimes.dat")
      chimesFile = unquote(buffer)
      call getParamSearchPaths(searchPath)
      call findFile(searchPath, chimesFile, chimesRepInput%chimesFile)
      if (.not. allocated(chimesRepInput%chimesFile)) then
        call error("Could not find ChIMES parameter file '" // chimesFile // "'")
      end if
    #:else
      call dftbp_error(chimes, "ChIMES repulsive correction requested, but code was compiled&
          & without ChIMES support")
    #:endif

    ! -- Schema validation (proof-of-concept, warnings only) --
    block
      type(hsd_schema_t) :: schema
      type(hsd_error_t), allocatable :: schemaErrors(:)
      integer :: iErr

      call schema_init(schema, name="Chimes")
      call schema_add_field(schema, "ParameterFile", FIELD_OPTIONAL, FIELD_TYPE_STRING)
      call schema_validate(schema, chimes, schemaErrors)
      if (size(schemaErrors) > 0) then
        do iErr = 1, size(schemaErrors)
          call dftbp_warning(chimes, "[schema] " // schemaErrors(iErr)%message)
        end do
      end if
      call schema_destroy(schema)
    end block

  end subroutine parseChimes


end module dftbp_dftbplus_parser_hybrid
