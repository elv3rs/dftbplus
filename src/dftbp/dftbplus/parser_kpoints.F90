#:include 'common.fypp'
#:include 'error.fypp'

!> Parser routines for K-point sampling input blocks.
module dftbp_dftbplus_parser_kpoints
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_status, only : TStatus
  use dftbp_dftb_hybridxc, only : checkSupercellFoldingMatrix
  use dftbp_dftb_periodic, only : getSuperSampling
  use dftbp_dftbplus_inputdata, only : TControl
  use dftbp_io_charmanip, only : tolower
  use dftbp_io_hsdcompat, only : hsd_table, textNodeName, detailedError, &
      & getChild, getChildValue, getNodeName
  use dftbp_io_message, only : error, warning
  use dftbp_math_simplealgebra, only : determinant33, diagonal
  use dftbp_type_linkedlist, only : asArray, asVector, destruct, get, init, len, &
      & TListIntR1, TListRealR1
  use dftbp_type_typegeometry, only : TGeometry

  implicit none

  private
  public :: readKPoints, maxSelfConsIterations

contains

  subroutine readKPoints(node, ctrl, geo, errStatus)

    !> Relevant node in input tree
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure
    type(TGeometry), intent(in) :: geo

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    ! Assume SCC can has usual default number of steps if needed
    ctrl%poorKSampling = .false.

    ! We can omit any hybrid xc-functional related checks for helical boundary conditions, since
    ! such a calculation will nevertheless be stopped due to the incompatibility of these features
    ctrl%checkStopHybridCalc = .false.

    ! K-Points
    if (geo%tPeriodic) then
      call getEuclideanKSampling(ctrl, node, geo, errStatus)
      @:PROPAGATE_ERROR(errStatus)
    elseif (geo%tHelical) then
      call getHelicalKSampling(ctrl, node, geo)
    end if

    call maxSelfConsIterations(node, ctrl, "MaxSCCIterations", ctrl%maxSccIter)
    ! Eventually, perturbation routines should also have restart reads:
    if (ctrl%poorKSampling .and. ctrl%tSCC .and. .not.ctrl%tReadChrg) then
      call warning("It is strongly suggested you use the ReadInitialCharges option.")
    end if

    ! Check if hybrid calculation needs to be stopped due to invalid k-point sampling
    if (ctrl%checkStopHybridCalc) then
      if (ctrl%maxSccIter == 1) then
        call warning("Restarting a hybrid xc-functional run with what appears to be&
            & a poor k-point sampling that does probably" // NEW_LINE('A') // " not match the&
            & original sampling (however fine for bandstructure calculations).")
      else
        call error("Error while parsing k-point sampling for a hybrid xc-functional&
            & run." // NEW_LINE('A') // "   Only allowed for bandstructure calculations,&
            & i.e. a single SCC iteration.")
      end if
    end if

  end subroutine readKPoints


  !> Set the maximum number of SCC cycles, depending on k-point behaviour
  subroutine maxSelfConsIterations(node, ctrl, label, maxSccIter)

    !> Relevant node in input tree
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Name of the tag
    character(*), intent(in) :: label

    !> Number of self-consistent iterations
    integer, intent(out) :: maxSccIter

    ! string for error return
    character(lc) :: warningStr

    integer :: ii

    maxSccIter = 1
    if (ctrl%tSCC) then
      if (ctrl%poorKSampling) then
        ! prevent full SCC with these points
        ii = 1
      else
        ii = 100
      end if
      call getChildValue(node, trim(label), maxSccIter, ii)
    end if

    if (ctrl%poorKSampling .and. maxSccIter /= 1) then
      write(warningStr, "(A,I3)") "A self-consistent cycle with these k-points probably will&
          & not correctly calculate many properties, maximum iterations set to:", maxSccIter
      call warning(warningStr)
    end if

  end subroutine maxSelfConsIterations


  !> Tries to infer whether the k-point sampling is restricted to the Gamma-point.
  pure function isGammaOnly(nKPoint, kPoint, kWeight)

    !> Number of k-points for the calculation
    integer, intent(in) :: nKPoint

    !> The k-points for the system
    real(dp), intent(in) :: kPoint(:,:)

    !> Weights for the k-points
    real(dp), intent(in) :: kWeight(:)

    !> True, if this appears to be a Gamma-only calculation
    logical :: isGammaOnly

    if (.not. nKPoint == 1) then
      isGammaOnly = .false.
    else
      isGammaOnly = .not. ((.not. all(abs(kPoint(:, 1)) < 1.0e-08_dp))&
          & .or. (.not. abs(kWeight(1)) - 1.0_dp < 1.0e-08_dp))
    end if

  end function isGammaOnly


  !> The k-points in Euclidean space
  subroutine getEuclideanKSampling(ctrl, node, geo, errStatus)

    !> Relevant node in input tree
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure
    type(TGeometry), intent(in) :: geo

    !> Error status
    type(TStatus), intent(inout) :: errStatus

    type(hsd_table), pointer :: value1, child
    character(len=:), allocatable :: buffer, modifier
    integer :: ind, ii, jj, kk
    real(dp), target :: coeffsAndShifts(3, 4)
    real(dp) :: rTmp3(3)
    type(TListIntR1) :: li1
    type(TListRealR1) :: lr1
    integer, allocatable :: tmpI1(:)
    real(dp), allocatable :: kpts(:,:)

    !! True, if k-points should be reduced by inversion
    logical :: tReduceByInversion

    !! True, if a Gamma-only k-point sampling is requested
    logical :: tGammaOnly

    call getChildValue(node, "KPointsAndWeights", value1, child=child, modifier=modifier)
    call getNodeName(value1, buffer)

    select case(buffer)

    case ("supercellfolding")
      ctrl%poorKSampling = .false.
      if (len(modifier) > 0) then
        call detailedError(child, "No modifier is allowed, if the SupercellFolding scheme is used.")
      end if
      call getChildValue(value1, "", coeffsAndShifts)
      if (abs(determinant33(coeffsAndShifts(:,1:3))) - 1.0_dp < -1e-06_dp) then
        call detailedError(value1, "Determinant of the supercell matrix must be greater than 1")
      end if
      if (any(abs(modulo(coeffsAndShifts(:,1:3) + 0.5_dp, 1.0_dp) - 0.5_dp)&
          & > 1e-06_dp)) then
        call detailedError(value1, "The components of the supercell matrix must be integers.")
      end if
      if (allocated(ctrl%hybridXcInp)) then
        call checkSupercellFoldingMatrix(coeffsAndShifts, errStatus)
        ctrl%supercellFoldingDiag = nint(diagonal(coeffsAndShifts(:,:3)))
        @:PROPAGATE_ERROR(errStatus)
        ctrl%supercellFoldingMatrix = coeffsAndShifts
      end if
      tReduceByInversion = (.not. ctrl%tSpinOrbit)
      call getSuperSampling(coeffsAndShifts(:,1:3), modulo(coeffsAndShifts(:,4), 1.0_dp),&
          & ctrl%kPoint, ctrl%kWeight, reduceByInversion=tReduceByInversion)
      ctrl%nKPoint = size(ctrl%kPoint, dim=2)

    case ("klines")
      ! probably unable to integrate charge for SCC
      ctrl%poorKSampling = .true.
      call init(li1)
      call init(lr1)
      call getChildValue(value1, "", 1, li1, 3, lr1)
      if (len(li1) < 1) then
        call detailedError(value1, "At least one line must be specified.")
      end if
      allocate(tmpI1(len(li1)))
      allocate(kpts(3, 0:len(lr1)))
      call asVector(li1, tmpI1)
      call asArray(lr1, kpts(:,1:len(lr1)))
      kpts(:,0) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
      call destruct(li1)
      call destruct(lr1)
      if (any(tmpI1 < 0)) then
        call detailedError(value1, "Interval steps must be greater equal to &
            &zero.")
      end if
      ctrl%nKPoint = sum(tmpI1)
      if (ctrl%nKPoint < 1) then
        call detailedError(value1, "Sum of the interval steps must be greater &
            &than zero.")
      end if
      ii = 1
      do while (tmpI1(ii) == 0)
        ii = ii + 1
      end do
      allocate(ctrl%kPoint(3, ctrl%nKPoint))
      allocate(ctrl%kWeight(ctrl%nKPoint))
      ind = 1
      do jj = ii, size(tmpI1)
        if (tmpI1(jj) == 0) then
          cycle
        end if
        rTmp3(:) = (kpts(:,jj) - kpts(:,jj-1)) / real(tmpI1(jj), dp)
        do kk = 1, tmpI1(jj)
          ctrl%kPoint(:,ind) = kpts(:,jj-1) + real(kk, dp) * rTmp3
          ind = ind + 1
        end do
      end do
      ctrl%kWeight(:) = 1.0_dp
      if (len(modifier) > 0) then
        select case (tolower(modifier))
        case ("relative")
        case ("absolute")
          ctrl%kPoint(:,:) =  matmul(transpose(geo%latVecs), ctrl%kPoint)
          kpts(:,:) = matmul(transpose(geo%latVecs), kpts)
        case default
          call detailedError(child, "Invalid modifier: '" // modifier &
              &// "'")
        end select
      end if
      deallocate(tmpI1)
      deallocate(kpts)

    case (textNodeName)

      ! no idea, but assume user knows what they are doing
      ctrl%poorKSampling = .false.

      call init(lr1)
      call getChildValue(child, "", 4, lr1, modifier=modifier)
      if (len(lr1) < 1) then
        call detailedError(child, "At least one k-point must be defined.")
      end if
      ctrl%nKPoint = len(lr1)
      allocate(kpts(4, ctrl%nKPoint))
      call asArray(lr1, kpts)
      call destruct(lr1)
      if (len(modifier) > 0) then
        select case (tolower(modifier))
        case ("relative")
          continue
        case ("absolute")
          kpts(1:3,:) =  matmul(transpose(geo%latVecs), kpts(1:3,:))
        case default
          call detailedError(child, "Invalid modifier: '" // modifier &
              &// "'")
        end select
      end if
      allocate(ctrl%kPoint(3, ctrl%nKPoint))
      allocate(ctrl%kWeight(ctrl%nKPoint))
      ctrl%kPoint(:,:) = kpts(1:3, :)
      ctrl%kWeight(:) = kpts(4, :)
      deallocate(kpts)
    case default
      call detailedError(value1, "Invalid K-point scheme")
    end select

    ! Catch problematic k-point sampling in case this is a hybrid calculation
    ctrl%checkStopHybridCalc = allocated(ctrl%hybridXcInp) .and. geo%tPeriodic&
        & .and. (buffer /= "supercellfolding") .and. ctrl%tReadChrg

    ! Check for hybrid xc-functional requirements
    tGammaOnly = isGammaOnly(ctrl%nKPoint, ctrl%kPoint, ctrl%kWeight)
    if (.not. tGammaOnly) then
      if (allocated(ctrl%hybridXcInp) .and. geo%tPeriodic&
          & .and. (buffer /= "supercellfolding") .and. (.not. ctrl%tReadChrg)) then
        call detailedError(child, "Error while parsing k-point sampling for a hybrid xc-functional&
            & run. Currently only" // NEW_LINE('A') // "   the supercell folding technique (or any&
            & format specifying the Gamma-point only)" // NEW_LINE('A') // "   is supported.")
      end if
    end if

    ! Hybrid calculations expect the supercell folding coefficients/shifts to be present
    if (allocated(ctrl%hybridXcInp) .and. tGammaOnly&
        & .and. (buffer /= "supercellfolding")) then
      coeffsAndShifts(:,:) = 0.0_dp
      do ii = 1, 3
        coeffsAndShifts(ii, ii) = 1.0_dp
      end do
      ctrl%supercellFoldingMatrix = coeffsAndShifts
      ctrl%supercellFoldingDiag = nint(diagonal(coeffsAndShifts(:,:3)))
    end if

  end subroutine getEuclideanKSampling


  !> The k-points for helical boundaries
  subroutine getHelicalKSampling(ctrl, node, geo)

    !> Relevant node in input tree
    type(hsd_table), pointer :: node

    !> Control structure to be filled
    type(TControl), intent(inout) :: ctrl

    !> Geometry structure
    type(TGeometry), intent(in) :: geo

    character(len=:), allocatable :: buffer
    type(hsd_table), pointer :: value1, child
    type(TListRealR1) :: lr1
    real(dp):: rTmp3(3), rTmp22(2,2)
    integer :: iTmp, iTmp2(2), kk, ii, jj
    real(dp), allocatable :: kPts(:,:)
    character(lc) :: errorStr

    ! assume the user knows what they are doing
    ctrl%poorKSampling = .false.

    call getChildValue(node, "KPointsAndWeights", value1, child=child)
    call getNodeName(value1, buffer)
    select case(buffer)
    case ("helicaluniform")
      call getChildValue(value1, "", rTmp3(:2))
      if (abs(modulo(rTmp3(1) + 0.5_dp, 1.0_dp) - 0.5_dp) > 1e-6_dp) then
        call detailedError(value1, "The k-point grid must be integer values.")
      end if
      iTmp = nint(rTmp3(1))
      if (iTmp < 1) then
        call detailedError(node, "Number of grid points must be above 0")
      end if
      if (.not.ctrl%tSpinOrbit) then
        ctrl%nKPoint = iTmp * nint(geo%latvecs(3,1))
        allocate(ctrl%kPoint(2, ctrl%nKPoint))
        ctrl%kPoint(:,:) = 0.0_dp
        allocate(ctrl%kWeight(ctrl%nKPoint))
        ctrl%kWeight(:) = 1.0_dp / real(iTmp,dp)
        do ii = 0, iTmp-1
          ctrl%kPoint(1,ii+1) = ii * 0.5_dp*ctrl%kWeight(ii+1) + 0.5_dp*rTmp3(2)/rTmp3(1)
        end do
        ctrl%kWeight(:) = 1.0_dp / real(ctrl%nKPoint,dp)
        do ii = 2, nint(geo%latvecs(3,1))
          ctrl%kPoint(1,(ii-1)*iTmp+1:ii*iTmp) = ctrl%kPoint(1,1:iTmp)
          ctrl%kPoint(2,(ii-1)*iTmp+1:ii*iTmp) = real(ii-1,dp)/nint(geo%latvecs(3,1))
        end do
      else
        call error("Helical boundaries not yet added for spin-orbit")
      end if
    case ("helicalsampled")
      call getChildValue(value1, "", rTmp22)
      iTmp2 = nint(rTmp22(:,1))
      if (any(abs(iTmp2-rTmp22(:,1)) > 1e-6_dp)) then
        call detailedError(value1, "The k-point grid must be integers.")
      end if
      if (any(iTmp2 < 1)) then
        call detailedError(node, "Number of grid points must be above 0")
      end if
      if (iTmp2(2) > nint(geo%latvecs(3,1))) then
        write(errorStr, '("The k-point grid for the helix rotational operation (",I0,&
            & ") is larger than the rotation order (C_",I0,").")') iTmp2(2), nint(geo%latvecs(3,1))
        call detailedError(node, errorStr)
      end if
      if (mod(nint(geo%latvecs(3,1)),iTmp2(2)) /= 0) then
        write(errorStr, '("The k-point grid for the helix rotational operation (n_k=",I0,&
            & ") is not a divisor of the rotation order (C_",I0,").")') iTmp2(2),&
            & nint(geo%latvecs(3,1))
        call detailedError(node, errorStr)
      end if
      if (abs(rTmp22(2,2) * nint(geo%latvecs(3,1)) - nint(rTmp22(2,2) * nint(geo%latvecs(3,1))))&
          & > epsilon(1.0_dp)) then
        write(errorStr, '("The shift of the k-points along the rotation is incommensurate, it must&
            & be an integer multiple of 1/",I0)') nint(geo%latvecs(3,1))
        call detailedError(node, errorStr)
      end if
      if (.not.ctrl%tSpinOrbit) then
        ctrl%nKPoint = product(iTmp2)
        allocate(ctrl%kPoint(2, ctrl%nKPoint))
        ctrl%kPoint(:,:) = 0.0_dp
        allocate(ctrl%kWeight(ctrl%nKPoint))

        kk = 1
        do ii = 0, iTmp2(1)-1
          do jj = 0, iTmp2(2)-1
            ctrl%kPoint(1,kk) = ii * 0.5_dp / rTmp22(1,1) + 0.5_dp*rTmp22(1,2)/rTmp22(1,1)
            ctrl%kPoint(2,kk) = mod(jj * 1.0_dp / rTmp22(2,1) + rTmp22(2,2), 1.0_dp)
            kk = kk + 1
          end do
        end do

        ctrl%kWeight(:) = 1.0_dp / real(ctrl%nKPoint,dp)

      else
        call error("Helical boundaries not yet added for spin-orbit")
      end if

    case (textNodeName)

      call init(lr1)
      call getChildValue(child, "", 3, lr1)
      if (len(lr1) < 1) then
        call detailedError(child, "At least one k-point must be defined.")
      end if
      ctrl%nKPoint = len(lr1)
      allocate(kpts(3, ctrl%nKPoint))
      call asArray(lr1, kpts)
      call destruct(lr1)
      allocate(ctrl%kPoint(2, ctrl%nKPoint))
      allocate(ctrl%kWeight(ctrl%nKPoint))
      ! first two are point values
      ctrl%kPoint(:2,:) = kpts(:2, :)
      ! test if the second k-point is commensurate with the C_n operation
      if (any(abs(kpts(2,:)*nint(geo%latvecs(3,1)) - nint(kpts(2,:) * nint(geo%latvecs(3,1))))&
          & > epsilon(1.0_dp))) then
        call error("Specified k-value(s) incommensurate with C_n operation.")
      end if
      ! last one is the weight
      ctrl%kWeight(:) = kpts(3, :)
      deallocate(kpts)

    case default
      call detailedError(value1, "Invalid K-point scheme")
    end select

  end subroutine getHelicalKSampling


end module dftbp_dftbplus_parser_kpoints
