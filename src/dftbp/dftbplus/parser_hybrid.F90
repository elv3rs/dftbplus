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
  use dftbp_io_hsdutils, only : getChild, getChildValue
  use dftbp_io_hsdutils, only : dftbp_error, getNodeName, getNodeHSDName
  use dftbp_io_unitconv, only : convertUnitHsd
  use dftbp_io_message, only : error
  use hsd_data, only : hsd_table
  use dftbp_type_linkedlist, only : get, TListCharLc
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
    type(TListCharLc), intent(inout) :: skFiles(:,:)

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

    @:ASSERT(size(skFiles, dim=1) == size(skFiles, dim=2))
    @:ASSERT((size(skFiles, dim=1) > 0))

    ! Extracting hybridXc tag from first SK-file only is a workaround and assumes that a set of
    ! given SK-files uses the same parameters (which should always be the case)!
    call get(skFiles(1, 1), fileName, 1)

    ! Check if SK-files contain extra tag for hybrid xc-functionals
    call parseHybridXcTag(fileName, hybridXcTag=hybridXcSkTag, hybridXcType=hybridXcSkType)
    isHybridSk = hybridXcSkTag /= hybridXcFunc%none

    call getChild(node, "Hybrid", child=hybridChild, requested=.false.)
    call getChildValue(node, "Hybrid", hybridValue, "None", child=hybridChild)
    call getNodeName(hybridValue, buffer)

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
      call getChildValue(hybridValue, "Screening", screeningValue, "MatrixBased",&
          & child=screeningChild)

      call getNodeName(screeningValue, buffer)
      select case(tolower(buffer))
      case ("neighbourbased")
        input%hybridXcAlg = hybridXcAlgo%neighbourBased
        call getChildValue(screeningValue, "CutoffReduction", input%cutoffRed, 0.0_dp,&
            & modifier=modifier, child=child1)
        call convertUnitHsd(modifier, lengthUnits, child1, input%cutoffRed)
        if (geo%tPeriodic) then
          call getChildValue(screeningValue, "Threshold", input%screeningThreshold, 1e-6_dp)
        end if
      case ("thresholded")
        input%hybridXcAlg = hybridXcAlgo%thresholdBased
        call getChildValue(screeningValue, "Threshold", input%screeningThreshold, 1e-6_dp)
        call getChildValue(screeningValue, "CutoffReduction", input%cutoffRed, 0.0_dp,&
            & modifier=modifier, child=child1)
        call convertUnitHsd(modifier, lengthUnits, child1, input%cutoffRed)
      case ("matrixbased")
        input%hybridXcAlg = hybridXcAlgo%matrixBased
        ! In this case, CutoffRedunction is not used so it should be set to zero.
        input%cutoffRed = 0.0_dp
      case default
        call getNodeHSdName(screeningValue, buffer)
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
        call getChildValue(hybridValue, "CoulombMatrix", cmValue, "Truncated", child=cmChild)

        call getNodeName(cmValue, buffer)
        select case(buffer)
        case ("full")
          input%gammaType = hybridXcGammaTypes%full
        case ("minimumimage")
          input%gammaType = hybridXcGammaTypes%mic
        case ("truncated")
          input%gammaType = hybridXcGammaTypes%truncated
        case ("truncatedanddamped")
          input%gammaType = hybridXcGammaTypes%truncatedAndDamped
        case default
          call getNodeHSdName(cmValue, buffer)
          call dftbp_error(cmChild, "Invalid Gamma function type '" // strBuffer // "'")
        end select

        ! g-Summation cutoff not needed for MIC CAM Hamiltonian
        if (input%gammaType /= hybridXcGammaTypes%mic) then
          call getChild(cmValue, "GSummationCutoff", child=child1, modifier=modifier,&
              & requested=.false.)
          if (associated(child1)) then
            allocate(input%gSummationCutoff)
            call getChildValue(child1, "", input%gSummationCutoff, modifier=modifier, child=child2)
            call convertUnitHsd(modifier, lengthUnits, child2, input%gSummationCutoff)
          end if
        end if

        if (input%gammaType == hybridXcGammaTypes%truncated&
            & .or. input%gammaType == hybridXcGammaTypes%truncatedAndDamped) then
          call getChild(cmValue, "CoulombCutoff", child=child1, modifier=modifier,&
              & requested=.false.)
          if (associated(child1)) then
            allocate(input%gammaCutoff)
            call getChildValue(child1, "", input%gammaCutoff, modifier=modifier, child=child2)
            call convertUnitHsd(modifier, lengthUnits, child2, input%gammaCutoff)
          end if
        end if

      else
        ! Always use unaltered gamma function for non-periodic systems
        input%gammaType = hybridXcGammaTypes%full
      end if ifPeriodic

      ! Number of primitive cells regarded in MIC, along each supercell folding direction
      if (input%gammaType == hybridXcGammaTypes%mic) then
        allocate(input%wignerSeitzReduction)
        call getChildValue(cmValue, "WignerSeitzReduction", input%wignerSeitzReduction, default=0)
      end if

    end if

  end subroutine parseHybridBlock



  !> Parses Chimes related options.
  subroutine parseChimes(root, chimesRepInput)
    type(hsd_table), pointer, intent(in) :: root
    type(TChimesRepInp), allocatable, intent(out) :: chimesRepInput

    type(hsd_table), pointer :: chimes

  #:if WITH_CHIMES
    character(len=:), allocatable :: buffer
    character(len=:), allocatable :: searchPath(:)
    character(len=:), allocatable :: chimesFile
  #:endif

    call getChild(root, "Chimes", chimes, requested=.false.)
    if (.not. associated(chimes)) return
    #:if WITH_CHIMES
      allocate(chimesRepInput)
      call getChildValue(chimes, "ParameterFile", buffer, default="chimes.dat")
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

  end subroutine parseChimes


end module dftbp_dftbplus_parser_hybrid
