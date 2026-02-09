#:include 'common.fypp'
#:include 'error.fypp'

!> Parser routines for transport/NEGF/Poisson input blocks.
module dftbp_dftbplus_parser_transport
  use dftbp_common_accuracy, only : dp, lc, mc
  use dftbp_common_constants, only : Bohr__AA
  use dftbp_common_globalenv, only : stdout
  use dftbp_common_unitconversion, only : energyUnits, lengthUnits
  use dftbp_dftbplus_inputconversion, only : transformPdosRegionInfo
  use dftbp_dftbplus_inputdata, only : TInputData
#:if WITH_POISSON
  use dftbp_extlibs_poisson, only : TPoissonInfo
  use dftbp_poisson_boundaryconditions, only : poissonBCsEnum, bcPoissonNames
#:endif
  use dftbp_io_charmanip, only : i2c, tolower, unquote
  use hsd_data, only : hsd_table, new_table
  use hsd, only : hsd_get, hsd_get_attrib, hsd_get_or_set, hsd_table_ptr, hsd_get_child_tables, &
      & hsd_get_table, hsd_set, hsd_get_choice, HSD_STAT_OK
  use dftbp_io_hsdutils, only : dftbp_error, dftbp_warning, getSelectedAtomIndices,&
      & getNodeName, getNodeName2, hasInlineData
  use dftbp_io_unitconv, only : convertUnitHsd
  use dftbp_io_message, only : error, warning
  use dftbp_math_simplealgebra, only : cross3
  use dftbp_type_commontypes, only : TOrbitals

  use dftbp_type_typegeometry, only : reduce, setLattice, TGeometry
  use dftbp_type_wrappedintr, only : TWrappedInt1
#:if WITH_TRANSPORT
  use dftbp_transport_negfvars, only : ContactInfo, TElPh, TNEGFGreenDensInfo, TNEGFTunDos
#:endif
  use dftbp_transport_negfvars, only : TTransPar

  implicit none

  private
#:if WITH_TRANSPORT
  public :: readTransportGeometry, readGreensFunction, readTunAndDos, finalizeNegf
#:endif
#:if WITH_POISSON
  public :: readPoisson
#:endif

contains

#:if WITH_TRANSPORT
  !> Read geometry information for transport calculation
  subroutine readTransportGeometry(root, geom, transpar)

    !> Root node containing the current block
    type(hsd_table), pointer :: root

    !> geometry of the system, which may be modified for some types of calculation
    type(TGeometry), intent(inout) :: geom

    !> Parameters of the transport calculation
    type(TTransPar), intent(inout) :: transpar

    type(hsd_table), pointer :: pDevice, pTask, pTaskType
    character(len=:), allocatable :: buffer, modifier
    type(hsd_table), pointer :: pTmp, field
    type(hsd_table_ptr), allocatable :: pNodeList(:)
    integer :: contact
    real(dp) :: lateralContactSeparation
    logical, allocatable :: atomInRegion(:)
    integer :: ii
    character(lc) :: strTmp
    integer :: stat

    transpar%defined = .true.
    transpar%tPeriodic1D = .not. geom%tPeriodic

    !! Note: we parse first the task because we need to know it to define the
    !! mandatory contact entries. On the other hand we need to wait that
    !! contacts are parsed to resolve the name of the contact for task =
    !! contacthamiltonian
    call hsd_get_choice(root, "Task", buffer, pTaskType, stat)
    if (stat /= HSD_STAT_OK) then
      buffer = 'uploadcontacts'
      pTaskType => null()
      call hsd_set(root, "Task", 'uploadcontacts')
    end if
    call hsd_get_table(root, "Task", pTask, stat, auto_wrap=.true.)
    if (.not. associated(pTask)) pTask => root

    call hsd_get_table(root, "Device", pDevice, stat, auto_wrap=.true.)
    if (.not. associated(pDevice)) call dftbp_error(root, "Missing required block: 'Device'")
    block
      integer, allocatable :: tmpArr(:)
      call hsd_get(pDevice, "AtomRange", tmpArr, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(pDevice, "Missing required array: 'AtomRange'")
      transpar%idxdevice(1:2) = tmpArr(1:2)
    end block
    call hsd_get_table(pDevice, "FirstLayerAtoms", pTmp, stat, auto_wrap=.true.)
    call readFirstLayerAtoms(pTmp, transpar%PL, transpar%nPLs, transpar%idxdevice)
    if (.not.associated(pTmp)) then
      call hsd_set(pDevice, "FirstLayerAtoms", transpar%PL)
    end if

    call hsd_get_child_tables(root, "Contact", pNodeList)
    transpar%ncont = size(pNodeList)
    allocate(transpar%contacts(transpar%ncont))
    call readContacts(pNodeList, transpar%contacts, geom, buffer, transpar%contactLayerTol)

    ! check for atoms in multiple contact ranges/device or atoms missing from any of these regions
    allocate(atomInRegion(geom%nAtom), source=.false.)
    atomInRegion(transpar%idxdevice(1):transpar%idxdevice(2)) = .true.
    do ii = 1, transpar%nCont
      if (any(atomInRegion(transpar%contacts(ii)%idxrange(1):transpar%contacts(ii)%idxrange(2))))&
          & then
        write(strTmp, "(A,A,A,I0)")"Contact '", trim(transpar%contacts(ii)%name),&
            & "' contains an atom already in the device region or another contact: Atom nr. ",&
            & findloc(atomInRegion(transpar%contacts(ii)%idxrange(1):&
            & transpar%contacts(ii)%idxrange(2)), .true.) + transpar%contacts(ii)%idxrange(1)
        pTmp => pNodeList(ii)%ptr
        call dftbp_error(pTmp, strTmp)
      end if
      atomInRegion(transpar%contacts(ii)%idxrange(1):transpar%contacts(ii)%idxrange(2)) = .true.
    end do
    if (any(.not.atomInRegion)) then
      write(strTmp, "(A,I0,A)")"Atom ", findloc(atomInRegion, .false.),&
          & " is not in the device region or any contact"
      call dftbp_error(root, strTmp)
    end if

    transpar%taskUpload = .false.

    select case (buffer)
    case ("contacthamiltonian")

      call hsd_get(pTaskType, "ContactId", buffer, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(pTaskType, "Missing required value: 'ContactId'")
      call hsd_get_table(pTaskType, "ContactId", pTmp, stat, auto_wrap=.true.)
      if (.not. associated(pTmp)) pTmp => pTaskType
      contact = getContactByName(transpar%contacts(:)%name, tolower(trim(unquote(buffer))),&
          & pTmp)
      transpar%taskContInd = contact
      if (.not. geom%tPeriodic) then
        call hsd_get_or_set(pTaskType, "ContactSeparation", lateralContactSeparation, 1000.0_dp,&
            & child=field)
        call hsd_get_attrib(pTaskType, "ContactSeparation", modifier)
        call convertUnitHsd(modifier,lengthUnits,field,lateralContactSeparation)
      end if

      call reduceGeometry(transpar%contacts(contact)%lattice, transpar%contacts(contact)%idxrange,&
          & lateralContactSeparation, geom)

      transpar%ncont = 0

      call hsd_get_or_set(root, "writeBinaryContact", transpar%tWriteBinShift, .true.)

    case ("uploadcontacts")

      transpar%taskUpload = .true.

      call hsd_get_or_set(root, "readBinaryContact", transpar%tReadBinShift, .true.)

    case default

      call getNodeName2(pTaskType, buffer)
      call dftbp_error(pTask, "Invalid task '" // buffer // "'")

   end select

   geom%areContactsPresent = transpar%nCont > 0 .and. transpar%taskUpload

  end subroutine readTransportGeometry


  !> Reduce the geometry for the contact calculation
  subroutine reduceGeometry(contactVec, contactRange, lateralContactSeparation, geom)

    !> Vector between principle layers in the contact
    real(dp), intent(in) :: contactVec(3)

    !> Range of atoms in the contact
    integer, intent(in) :: contactRange(2)

    !> Lateral separation distance between contacts in a periodic box
    real(dp), intent(in) :: lateralContactSeparation

    !> atomic geometry
    type(TGeometry), intent(inout) :: geom

    real(dp) :: contUnitVec(3), dots(3), newLatVecs(3, 3), newOrigin(3)
    real(dp) :: minProj, maxProj
    logical :: mask(3)
    integer :: ind, indPrev, indNext, ii

    if (geom%tPeriodic) then
      contUnitVec = contactVec / sqrt(sum(contactVec**2, dim=1))
      dots = abs(matmul(contUnitVec, geom%latVecs))
      mask(:) = (abs(dots - sqrt(sum(geom%latVecs, dim=1)**2)) < 1e-8_dp)
      if (count(mask) /= 1) then
        call error("Too many lattice vectors parallel to the contact")
      end if
      ind = findloc(mask, .true., 1)
      newLatVecs = geom%latVecs
      newLatVecs(:,ind) = 2.0_dp * contactVec
      newOrigin = geom%origin
    else
      newLatVecs(:,1) = 2.0_dp * contactVec
      mask(:) = abs(contactVec) > 1e-8_dp
      ind = findloc(mask, .true., 1)
      ! Note: ind is one-based, subtract 1 before modulo and add 1 after.
      indNext = modulo(ind + 1 - 1, 3) + 1
      indPrev = modulo(ind - 1 - 1, 3) + 1
      newLatVecs(indNext, 2) = -newLatVecs(ind, 1)
      newLatVecs(ind, 2) = newLatVecs(indNext, 1)
      newLatVecs(indPrev, 2) = 0.0_dp
      newLatVecs(:,3) = cross3(newLatVecs(:,1), newLatVecs(:,2))
      newLatVecs(:,2) = newLatVecs(:,2) / sqrt(sum(newLatVecs(:,2)**2))
      newLatVecs(:,3) = newLatVecs(:,3) / sqrt(sum(newLatVecs(:,3)**2))
      newOrigin(:) = 0.0_dp
    end if
    call reduce(geom, contactRange(1), contactRange(2))
    if (.not. geom%tPeriodic) then
      do ii = 2, 3
        minProj = minval(matmul(newLatVecs(:,ii), geom%coords))
        maxProj = maxval(matmul(newLatVecs(:,ii), geom%coords))
        newLatVecs(:,ii) = ((maxProj - minProj) + lateralContactSeparation) * newLatVecs(:,ii)
      end do
    end if
    call setLattice(geom, newOrigin, newLatVecs)

  end subroutine reduceGeometry


  !> Reads settings for the first layer atoms in principal layers
  subroutine readFirstLayerAtoms(pnode, pls, npl, idxdevice, check)

    type(hsd_table), pointer, intent(in) :: pnode

    !> Start atoms in the principal layers
    integer, allocatable, intent(out) :: pls(:)

    !> Number of principal layers
    integer, intent(out) :: npl

    !> Atoms range of the device
    integer, intent(in) :: idxdevice(2)

    !> Optional setting to turn on/off check (defaults to on if absent)
    logical, optional, intent(in) :: check


    integer, allocatable :: intArr(:)
    integer :: stat
    logical :: checkidx

    checkidx = .true.
    if (present(check)) checkidx = check

    if (associated(pnode)) then
        call hsd_get(pnode, "#text", intArr, stat=stat)
        if (stat /= 0) call dftbp_error(pnode, "Error reading first layer atoms")
        npl = size(intArr)
        allocate(pls(npl))
        pls(:) = intArr
        if (checkidx) then
          if (any(pls < idxdevice(1) .or. &
                  pls > idxdevice(2))) then
             call dftbp_error(pnode, "First layer atoms must be between " &
               &// i2c(idxdevice(1)) // " and " // i2c(idxdevice(2)) // ".")
          end if
        end if
      else
         npl = 1
         allocate(pls(npl))
         pls = (/ 1 /)
      end if

  end subroutine readFirstLayerAtoms


  !> Reads Green's function settings
  subroutine readGreensFunction(pNode, greendens, transpar, tempElec)

    !> Input tree
    type(hsd_table), pointer :: pTmp

    !> Settings for Green's function solver
    type(TNEGFGreenDensInfo), intent(inout) :: greendens

    !> Transport solver settings
    type(TTransPar), intent(inout) :: transpar

    !> Electron temperature
    real(dp), intent(in) :: tempElec

    type(hsd_table), pointer :: pNode
    type(hsd_table), pointer :: field, child1, child2
    real(dp) :: Estep
    integer :: defValue, ii
    character(len=:), allocatable :: buffer, modifier
    logical :: realAxisConv, equilibrium

    real(dp), allocatable :: realArr(:)
    integer :: stat

    greendens%defined = .true.

    if (.not. transpar%defined) then
      !! Fermi level: in case of collinear spin we accept two values
      !! (up and down)
      call hsd_get(pNode, "FermiLevel", realArr, stat=stat)
      if (stat /= 0) call dftbp_error(pNode, "Error reading FermiLevel")
      call hsd_get_attrib(pNode, "FermiLevel", modifier)
      if (size(realArr) == 1) then
        greendens%oneFermi(1) = realArr(1)
        greendens%oneFermi(2) = realArr(1)
      else if (size(realArr) == 2) then
        greendens%oneFermi(1:2) = realArr(1:2)
      else
        call dftbp_error(pNode, "FermiLevel accepts 1 or 2 (for collinear spin) values")
      end if
      call convertUnitHsd(modifier, energyUnits, pNode, greendens%oneFermi)

      call hsd_get_table(pNode, "FirstLayerAtoms", pTmp, stat, auto_wrap=.true.)
      call readFirstLayerAtoms(pTmp, greendens%PL, greendens%nPLs,&
                                &transpar%idxdevice, check = .false.)
      if (.not.associated(pTmp)) then
        call hsd_set(pNode, "FirstLayerAtoms", greendens%PL)
      end if
      !call getChild(pNode, "ContactPLs", pTmp, requested=.false.)
      !if (associated(pTmp)) then
      !  call init(li)
      !  call getChildValue(pTmp, "", li)
      !  allocate(transpar%cblk(len(li)))
      !  call asArray(li,transpar%cblk)
      !  call destruct(li)
      !end if
      allocate(greendens%kbT(1))
      greendens%kbT(:) = tempElec
    else
      if (transpar%ncont > 0) then
        allocate(greendens%kbT(transpar%ncont))
        do ii = 1, transpar%ncont
          if (transpar%contacts(ii)%kbT .ge. 0.0_dp) then
            greendens%kbT(ii) = transpar%contacts(ii)%kbT
          else
            greendens%kbT(ii) = tempElec
          end if
        enddo
      end if
    end if

    call hsd_get_or_set(pNode, "LocalCurrents", greendens%doLocalCurr, .false.)
    call hsd_get_or_set(pNode, "Verbosity", greendens%verbose, 51)
    call hsd_get_or_set(pNode, "Delta", greendens%delta, 1.0e-5_dp, child=field)
    call hsd_get_attrib(pNode, "Delta", modifier)
    call convertUnitHsd(modifier, energyUnits, field, greendens%delta)
    call hsd_get_or_set(pNode, "ReadSurfaceGFs", greendens%readSGF, .false.)
    call hsd_get_or_set(pNode, "SaveSurfaceGFs", greendens%saveSGF, .not.greendens%readSGF)
    block
      integer, allocatable :: tmpArr(:)
      call hsd_get_or_set(pNode, "ContourPoints", tmpArr, [ 20, 20 ])
      greendens%nP(1:2) = tmpArr(1:min(size(tmpArr),2))
    end block
    call hsd_get_or_set(pNode, "EnclosedPoles",  greendens%nPoles, 3)
    call hsd_get_or_set(pNode, "LowestEnergy", greendens%enLow, -2.0_dp, child=field)
    call hsd_get_attrib(pNode, "LowestEnergy", modifier)
    call convertUnitHsd(modifier, energyUnits, field, greendens%enLow)
    call hsd_get_or_set(pNode, "FermiCutoff", greendens%nkT, 10)
      ! Fermi energy had not been set by other means yet

      ! Non equilibrium integration along real axis:
      ! The code will perform the integration if the number of points is larger
      ! than zero, no matter if there's bias or not.
      ! Therefore I restored the default on the energy step, as it works at zero
      ! bias and it scales flawlessy with increasing bias
      ! It is still allowed to directly set the number of points, if preferred
      ! libNEGF only wants the number of points in input
      call hsd_get_table(pNode, "RealAxisPoints", child1, stat, auto_wrap=.true.)
      call hsd_get_table(pNode, "RealAxisStep", child2, stat, auto_wrap=.true.)
      call hsd_get_attrib(pNode, "RealAxisStep", buffer)
      realAxisConv = .false.
      ! Set a bool to verify if all contacts are at the same potential (if so,
      ! no points are needed)
      equilibrium = .true.
      do ii = 2, transpar%ncont
        if (transpar%contacts(1)%potential .ne. transpar%contacts(ii)%potential &
           & .or. transpar%contacts(1)%kbT .ne. transpar%contacts(ii)%kbT ) then
           equilibrium = .false.
        end if
      end do

      ! Both Points and Step cannot be specified
      if  (associated (child1) .and. associated(child2)) then
        call dftbp_error(child1, "RealAxisPoints and RealAxisStep " &
                            &// " cannot be specified together.")
      ! If only one is specified, take it as valid value
      else if (associated(child1)) then
        call hsd_get(pNode, "RealAxisPoints", greendens%nP(3), stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(pNode, "Missing required value: 'RealAxisPoints'")
      else if (associated(child2)) then
        call hsd_get(pNode, "RealAxisStep", Estep, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(pNode, "Missing required value: 'RealAxisStep'")
        call hsd_get_attrib(pNode, "RealAxisStep", modifier)
        call convertUnitHsd(modifier, energyUnits, child2, Estep)
        realAxisConv = .true.
      ! If the system is under equilibrium we set the number of
      ! points to zero
      else if (equilibrium) then
        call hsd_get_or_set(pNode, "RealAxisPoints", greendens%nP(3), 0)
      else
        !Default is a point every 1500H
        call hsd_get_or_set(pNode, "RealAxisStep", Estep, 6.65e-4_dp, child=child2)
        call hsd_get_attrib(pNode, "RealAxisStep", modifier)
        realAxisConv = .true.
      end if
      ! RealAxisConv means that we have a step and we convert it in a number
      ! of points
      if (realAxisConv) then
        defValue = int(1.0_dp/Estep &
          & * (maxval(transpar%contacts(:)%potential) &
          & - minval(transpar%contacts(:)%potential) + &
          & 2 * greendens%nKT * maxval(greendens%kbT)))
        greendens%nP(3) = defvalue
        !call getChildValue(pNode, "RealAxisPoints", greendens%nP(3), &
        !    & defvalue, child=child1)
      end if

  end subroutine readGreensFunction
#:endif


#:if WITH_POISSON

  !> Read in Poisson related data
#:if WITH_TRANSPORT
  subroutine readPoisson(pNode, poisson, tPeriodic, transpar, latVecs, updateSccAfterDiag)
#:else
  subroutine readPoisson(pNode, poisson, tPeriodic, latVecs, updateSccAfterDiag)
#:endif

    !> Input tree
    type(hsd_table), pointer :: pNode

    !> data type for Poisson solver settings
    type(TPoissonInfo), intent(inout) :: poisson

    !> Is this a periodic calculation
    logical, intent(in) :: tPeriodic

  #:if WITH_TRANSPORT
    !> Parameters of the transport calculation
    type(TTransPar), intent(inout) :: transpar
  #:endif

    !> Lattice vectors if periodic
    real(dp), allocatable, intent(in) :: latVecs(:,:)

    !> Whether Scc should be updated with the output charges (obtained after diagonalisation)
    logical, intent(out) :: updateSccAfterDiag

    type(hsd_table), pointer :: pTmp, pTmp2, pChild, field
    character(len=:), allocatable :: buffer, modifier
    real(dp) :: denstol, gatelength_l
    logical :: needsPoissonBox
    integer :: stat

  #:if WITH_TRANSPORT
    integer :: ii
  #:endif

    poisson%defined = .true.
    needsPoissonBox = .not. tPeriodic
  #:if WITH_TRANSPORT
    needsPoissonBox = needsPoissonBox .or. transpar%tPeriodic1D .or. transpar%nCont == 1
  #:endif

    if (needsPoissonBox) then
    #:if WITH_TRANSPORT
      if (transpar%nCont == 1 .and. .not. transpar%tPeriodic1D) then
        poisson%poissBox(:) = 0.0_dp
        do ii = 1, 3
          if (ii == transpar%contacts(1)%dir) then
            call hsd_get(pNode, "PoissonThickness", poisson%poissBox(ii), stat=stat)
            if (stat /= HSD_STAT_OK) &
                & call dftbp_error(pNode, "Missing required value: 'PoissonThickness'")
            call hsd_get_attrib(pNode, "PoissonThickness", modifier)
            call hsd_get_table(pNode, "PoissonThickness", field, stat, auto_wrap=.true.)
            if (.not. associated(field)) field => pNode
            call convertUnitHsd(modifier, lengthUnits, field, poisson%poissBox)
          else
            poisson%poissBox(ii) = sqrt(sum(latVecs(:,ii)**2))
          end if
        end do
      else
        block
          real(dp), allocatable :: tmpRA(:)
          call hsd_get(pNode, "PoissonBox", tmpRA, stat=stat)
          if (stat /= HSD_STAT_OK) call dftbp_error(pNode, "Missing required array: 'PoissonBox'")
          poisson%poissBox(1:3) = tmpRA(1:3)
        end block
        call hsd_get_attrib(pNode, "PoissonBox", modifier)
        call hsd_get_table(pNode, "PoissonBox", field, stat, auto_wrap=.true.)
        if (.not. associated(field)) field => pNode
        call convertUnitHsd(modifier, lengthUnits, field, poisson%poissBox)
      end if
    #:else
      block
        real(dp), allocatable :: tmpRA(:)
        call hsd_get(pNode, "PoissonBox", tmpRA, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(pNode, "Missing required array: 'PoissonBox'")
        poisson%poissBox(1:3) = tmpRA(1:3)
      end block
      call hsd_get_attrib(pNode, "PoissonBox", modifier)
      call hsd_get_table(pNode, "PoissonBox", field, stat, auto_wrap=.true.)
      if (.not. associated(field)) field => pNode
      call convertUnitHsd(modifier, lengthUnits, field, poisson%poissBox)
    #:endif
    end if

    poisson%foundBox = needsPoissonBox
    block
      real(dp), allocatable :: tmpRA(:)
      call hsd_get_or_set(pNode, "MinimalGrid", tmpRA, [ 0.3_dp, 0.3_dp, 0.3_dp ])
      poisson%poissGrid(1:3) = tmpRA(1:min(size(tmpRA),3))
    end block
    call hsd_get_attrib(pNode, "MinimalGrid", modifier)
    call hsd_get_table(pNode, "MinimalGrid", field, stat, auto_wrap=.true.)
    if (.not. associated(field)) field => pNode
    call convertUnitHsd(modifier, lengthUnits, field, poisson%poissGrid)
    call hsd_get_or_set(pNode, "NumericalNorm", poisson%numericNorm, .false.)
    call hsd_get_table(pNode, "AtomDensityCutoff", pTmp, stat, auto_wrap=.true.)
    if (associated(pTmp)) call hsd_get_attrib(pNode, "AtomDensityCutoff", modifier)
    call hsd_get_table(pNode, "AtomDensityTolerance", pTmp2, stat, auto_wrap=.true.)
    if (associated(pTmp) .and. associated(pTmp2)) then
      call dftbp_error(pNode, "Only one of the tags AtomDensityCutoff or AtomDensityTolerance&
          & can be specified.")
    else if (associated(pTmp)) then
      call hsd_get_or_set(pTmp, "#text", poisson%maxRadAtomDens, 14.0_dp)
      call convertUnitHsd(modifier, lengthUnits, pTmp, poisson%maxRadAtomDens)
      if (poisson%maxRadAtomDens <= 0.0_dp) then
        call dftbp_error(pTmp2, "Atom density cutoff must be > 0")
      end if
    else
      call hsd_get_or_set(pNode, "AtomDensityTolerance", denstol, 1e-6_dp, child=pTmp2)
      if (denstol <= 0.0_dp) then
        call dftbp_error(pTmp2, "Atom density tolerance must be > 0")
      end if
      ! Negative value to signal automatic determination
      poisson%maxRadAtomDens = -denstol
    end if

    call hsd_get_or_set(pNode, "CutoffCheck", poisson%cutoffcheck, .true.)
    call hsd_get_or_set(pNode, "Verbosity", poisson%verbose, 51)
    call hsd_get_or_set(pNode, "SavePotential", poisson%savePotential, .false.)
    call hsd_get_or_set(pNode, "PoissonAccuracy", poisson%poissAcc, 1.0e-6_dp)
    call hsd_get_or_set(pNode, "BuildBulkPotential", poisson%bulkBC, .true.)
    call hsd_get_or_set(pNode, "ReadOldBulkPotential", poisson%readBulkPot, .false.)
    call hsd_get_or_set(pNode, "RecomputeAfterDensity", updateSccAfterDiag, .false.)
    call hsd_get_or_set(pNode, "MaxPoissonIterations", poisson%maxPoissIter, 60)

    poisson%overrideBC(:) = poissonBCsEnum%periodic
    call hsd_get_table(pNode, "OverrideDefaultBC", pTmp, stat, auto_wrap=.true.)
    if (associated(pTmp)) then
      call getPoissonBoundaryConditionOverrides(pTmp,&
          & [ poissonBCsEnum%dirichlet, poissonBCsEnum%neumann ], poisson%overrideBC)
    end if

    poisson%overrBulkBC(:) = poissonBCsEnum%unset
    call hsd_get_table(pNode, "OverrideBulkBC", pTmp, stat, auto_wrap=.true.)
    if (associated(pTmp)) then
      call getPoissonBoundaryConditionOverrides(pTmp,&
          & [ poissonBCsEnum%periodic, poissonBCsEnum%dirichlet, poissonBCsEnum%neumann ],&
          & poisson%overrBulkBC)
    end if

    call hsd_get_choice(pNode, "BoundaryRegion", buffer, pTmp, stat)
    if (stat /= HSD_STAT_OK) then
      buffer = "global"
      pTmp => null()
      call hsd_set(pNode, "BoundaryRegion", "global")
    end if
    select case(buffer)
    case ("global")
      poisson%localBCType = "G"
    case ("square")
      poisson%localBCType = "S"
      call hsd_get_or_set(pTmp, "BufferLength", poisson%bufferLocBC, 9.0_dp, child=field)
      call hsd_get_attrib(pTmp, "BufferLength", modifier)
      call convertUnitHsd(modifier, lengthUnits, field, poisson%bufferLocBC)
    case ("circle")
      poisson%localBCType = "C"
      call hsd_get_or_set(pTmp, "BufferLength", poisson%bufferLocBC, 9.0_dp, child=field)
      call hsd_get_attrib(pTmp, "BufferLength", modifier)
      call convertUnitHsd(modifier, lengthUnits, field, poisson%bufferLocBC)
    case default
      if (associated(pTmp)) call getNodeName2(pTmp, buffer)
      call dftbp_error(pTmp, "Invalid boundary region type '" // buffer // "'")
    end select

    call hsd_get_or_set(pNode, "BoxExtension", poisson%bufferBox, 0.0_dp, child=field)
    call hsd_get_attrib(pNode, "BoxExtension", modifier)
    call convertUnitHsd(modifier, lengthUnits, field, poisson%bufferBox)
    if (poisson%bufferBox.lt.0.0_dp) then
      call dftbp_error(pNode, "BoxExtension must be a positive number")
    endif

    ! PARSE GATE OPTIONS
    call hsd_get_choice(pNode, "Gate", buffer, pTmp2, stat)
    if (stat /= HSD_STAT_OK) then
      buffer = "none"
      pTmp2 => null()
      call hsd_set(pNode, "Gate", "none")
    end if
    call hsd_get_table(pNode, "Gate", pChild, stat, auto_wrap=.true.)
    if (.not. associated(pChild)) pChild => pNode

    poisson%insLength = 0.0_dp
    poisson%insRad = 0.0_dp
    select case(buffer)
    case ("none")
      poisson%gateType = "N"
    case ("planar")
      poisson%gateType = "P"
      call hsd_get_or_set(pTmp2, "GateLength", poisson%gateLength_l, 0.0_dp, child=field)
      call hsd_get_attrib(pTmp2, "GateLength", modifier)
      call convertUnitHsd(modifier, lengthUnits, field, poisson%gateLength_l)

      gatelength_l = poisson%gateLength_l !avoids a warning on intents
      call hsd_get_or_set(pTmp2, "GateLength_l", poisson%gateLength_l, gateLength_l, child=field)
      call hsd_get_attrib(pTmp2, "GateLength_l", modifier)
      call convertUnitHsd(modifier, lengthUnits, field, poisson%gateLength_l)

      call hsd_get_or_set(pTmp2, "GateLength_t", poisson%gateLength_t, poisson%gateLength_l,&
          & child=field)
      call hsd_get_attrib(pTmp2, "GateLength_t", modifier)
      call convertUnitHsd(modifier, lengthUnits, field, poisson%gateLength_t)

      call hsd_get_or_set(pTmp2, "GateDistance", poisson%gateRad, 0.0_dp, child=field)
      call hsd_get_attrib(pTmp2, "GateDistance", modifier)
      call convertUnitHsd(modifier, lengthUnits, field, poisson%gateRad)

      call hsd_get_or_set(pTmp2, "GatePotential", poisson%gatepot, 0.0_dp, child=field)
      call hsd_get_attrib(pTmp2, "GatePotential", modifier)
      call convertUnitHsd(modifier, energyUnits, field, poisson%gatepot)

      !call getChildValue(pTmp2, "GateDirection", poisson%gatedir, 2)
      poisson%gatedir = 2

    case ("cylindrical")
      poisson%gateType = "C"
      call hsd_get_or_set(pTmp2, "GateLength", poisson%gateLength_l, 0.0_dp, child=field)
      call hsd_get_attrib(pTmp2, "GateLength", modifier)
      call convertUnitHsd(modifier, lengthUnits, field, poisson%gateLength_l)

      call hsd_get_or_set(pTmp2, "GateRadius", poisson%gateRad, 0.0_dp, child=field)
      call hsd_get_attrib(pTmp2, "GateRadius", modifier)
      call convertUnitHsd(modifier, lengthUnits, field, poisson%gateRad)

      call hsd_get_or_set(pTmp2, "GatePotential", poisson%gatepot, 0.0_dp, child=field)
      call hsd_get_attrib(pTmp2, "GatePotential", modifier)
      call convertUnitHsd(modifier, energyUnits, field, poisson%gatepot)

    case default
      call getNodeName2(pTmp2, buffer)
      call dftbp_error(pTmp2, "Invalid gate type '" // buffer // "'")

    end select

    call hsd_get_or_set(pNode, "MaxParallelNodes", poisson%maxNumNodes, 1)

    poisson%scratch = "contacts"

  end subroutine readPoisson


  !> Over-rides the boundary conditions on the Poisson solver
  subroutine getPoissonBoundaryConditionOverrides(pNode, availableConditions, overrideBC)

    !> Input data tree
    type(hsd_table), pointer, intent(in) :: pNode

    !> List of conditions that can be set as choices
    integer, intent(in) :: availableConditions(:)

    !> Array of boundary condition types on the 6 faces of the box, 0 for use of default
    integer, intent(inout) :: overrideBC(:)

    integer :: bctype, iBC
    integer :: faceBC, oppositeBC
    integer :: ii, stat
    character(:), allocatable :: strArr(:)
    type(hsd_table), pointer :: pNode2, pChild
    character(lc) :: strTmp
    character(1), parameter :: sDirs(3) = ['x','y','z']

    call hsd_get_table(pNode, "none", pNode2, stat, auto_wrap=.true.)
    if (associated(pNode2)) return
    do iBC = 1, size(availableConditions)
      bctype = availableConditions(iBC)
      call hsd_get_table(pNode, trim(bcPoissonNames(bctype)), pNode2, stat, auto_wrap=.true.)
      if (associated(pNode2)) then
        call hsd_get(pNode2, "boundaries", strArr, stat=stat)
        if (stat /= 0) call dftbp_error(pNode2, "Error reading boundaries")
        call hsd_get_table(pNode2, "boundaries", pChild, stat, auto_wrap=.true.)
        if (.not. associated(pChild)) call dftbp_error(pNode2, "Missing required block: 'boundaries'")
        if (size(strArr) > 6) then
          call dftbp_error(pChild,"A maximum of 6 boundaries (or fewer) can be set")
        end if
        do ii = 1, size(strArr)
          strTmp = strArr(ii)
          select case(trim(strTmp))
          case("x")
            overrideBC(1) = bctype
            overrideBC(2) = bctype
          case("xmin")
            overrideBC(1) = bctype
          case("xmax")
            overrideBC(2) = bctype
          case("y")
            overrideBC(3) = bctype
            overrideBC(4) = bctype
          case("ymin")
            overrideBC(3) = bctype
          case("ymax")
            overrideBC(4) = bctype
          case("z")
            overrideBC(5) = bctype
            overrideBC(6) = bctype
          case("zmin")
            overrideBC(5) = bctype
          case("zmax")
            overrideBC(6) = bctype
          end select
        end do
      end if
    end do

    ! If a face is set to be periodic, the opposite one should be of the same type
    do ii = 1, 3
      faceBC = overrideBC(2 * ii)
      oppositeBC = overrideBC(2 * ii - 1)
      if (faceBC == poissonBCsEnum%periodic .neqv. oppositeBC == poissonBCsEnum%periodic) then
        call dftbp_error(pChild, "Periodic override must be set both min max along "//sDirs(ii))
      end if
    end do

  end subroutine getPoissonBoundaryConditionOverrides

#:endif


#:if WITH_TRANSPORT
  !> Correctness checking of atom ranges and returning contact vector and direction.
  subroutine getContactVector(atomrange, geom, id, name, pContact, contactLayerTol, contactVec,&
      & contactDir)

    !> Range of atoms in the contact
    integer, intent(in) :: atomrange(2)

    !> Atomic geometry, including the contact atoms
    type(TGeometry), intent(in) :: geom

    !> Index for this contact
    integer, intent(in) :: id

    !> Contact name
    character(mc), intent(in) :: name

    !> Node in the parser, needed for error handling
    type(hsd_table), pointer :: pContact

    !> Allowed discrepancy in positions of atoms between the contact's two  principle layers
    real(dp), intent(in) :: contactLayerTol

    !> Vector direction between principal layers in the contact
    real(dp), intent(out) :: contactVec(3)

    !> Which supercell vector the contact vector is parallel to
    integer, intent(out) :: contactDir

    integer :: iStart, iStart2, iEnd, ii
    logical :: mask(3)
    character(lc) :: errorStr

    !! Correctness check for the atom ranges
    iStart = atomrange(1)
    iEnd = atomrange(2)
    if (iStart < 1 .or. iEnd < 1 .or. iStart > geom%nAtom .or. iEnd > geom%nAtom) then
      call dftbp_error(pContact, "Invalid atom range '" // i2c(iStart) &
          &// " " // i2c(iEnd) // "', values should be between " // i2c(1) &
          &// " and " // i2c(geom%nAtom) // ".")
    end if
    if (iEnd < iStart) then
      call dftbp_error(pContact, "Invalid atom order in contact '" // i2c(iStart) // " " //&
          & i2c(iEnd) // "', should be asscending order.")
    end if

    if (mod(iEnd - iStart + 1, 2) /= 0) then
      call dftbp_error(pContact, "Nr. of atoms in the contact must be even")
    end if

    ! Determining intra-contact layer vector
    iStart2 = iStart + (iEnd - iStart + 1) / 2
    contactVec = geom%coords(:,iStart) - geom%coords(:,iStart2)

    if (any(sum( (geom%coords(:,iStart:iStart2-1) - geom%coords(:,iStart2:iEnd)&
        & - spread(contactVec, dim=2, ncopies=iStart2-iStart))**2, dim=1) > contactLayerTol**2))&
        & then
      write(stdout,"(1X,A,I0,A,I0)")'Contact vector defined from atoms ', iStart, ' and ',iStart2
      write(stdout,"(1X,A,I0,'-',I0)")'Contact layer 1 atoms: ',iStart, iStart2-1
      write(stdout,"(1X,A,I0,'-',I0)")'Contact layer 2 atoms: ',iStart2, iEnd
      do ii = 0, iStart2 -1 -iStart
        if (sum((geom%coords(:,ii+iStart)-geom%coords(:,ii+iStart2) - contactVec)**2)&
            & > contactLayerTol**2) then
          write(stdout,"(1X,A,I0,A,I0,A)")'Atoms ',iStart+ii, ' and ', iStart2+ii,&
              & ' inconsistent with the contact vector.'
          exit
        end if
      end do
      write(stdout,*)'Mismatches in atomic positions in the two layers:'
      write(stdout,"(3F20.12)")((geom%coords(:,iStart:iStart2-1) - geom%coords(:,iStart2:iEnd)&
          & - spread(contactVec(:), dim=2, ncopies=iStart2-iStart))) * Bohr__AA

      write (errorStr,"('Contact ',A,' (',A,') does not consist of two rigidly shifted layers')")&
          & i2c(id), trim(name)
      call error(errorStr)

    end if

    ! Determine to which axis the contact vector is parallel.
    mask(:) = (abs(abs(contactVec) - sqrt(sum(contactVec**2))) < 1.0e-8_dp)
    if (count(mask) /= 1) then
      call warning("Contact vector " // i2c(id) // " not parallel to any coordinate axis.")
      contactDir = 0
    else
      contactDir = findloc(mask, .true., 1)
    end if

  end subroutine getContactVector


  !> Read dephasing block
  subroutine readDephasing(node, orb, geom, tp, tundos)

    !> Input tree node
    type(hsd_table), pointer :: node

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Atomic geometry, including the contact atoms
    type(TGeometry), intent(in) :: geom

    !> Parameters of the transport calculation
    type(TTransPar), intent(inout) :: tp

    !> Parameters of tunneling and dos calculation
    type(TNEGFTunDos), intent(inout) :: tundos

    type(hsd_table), pointer :: value1, child
    integer :: stat
    character(len=:), allocatable :: tmpBuf

    call hsd_get_table(node, "VibronicElastic", child, stat, auto_wrap=.true.)
    if (associated(child)) then
      tp%tDephasingVE = .true.
      call readElPh(child, tundos%elph, geom, orb, tp)
    end if

    call hsd_get_table(node, "BuettikerProbes", child, stat, auto_wrap=.true.)
    if (associated(child)) then
      call hsd_get_choice(child, "", tmpBuf, value1, stat)
      if (stat /= HSD_STAT_OK) value1 => null()
    else
      value1 => null()
    end if
    if (associated(value1)) then
      tp%tDephasingBP = .true.
      call readDephasingBP(child, tundos%bp, geom, orb, tp)
    end if

    ! Lowdin transformations involve dense matrices and works only in small systems
    ! For the dftb+ official release the options are disabled
    tp%tOrthonormal = .false.
    tp%tOrthonormalDevice = .false.
    !call getChildValue(node, "Orthonormal", tp%tOrthonormal, .false.)
    !call getChildValue(node, "OrthonormalDevice", tp%tOrthonormalDevice, .false.)
    tp%tNoGeometry = .false.
    tp%NumStates = 0

  end subroutine readDephasing


  !> Read Electron-Phonon blocks (for density and/or current calculation)
  subroutine readElPh(node, elph, geom, orb, tp)

    !> Input node in the tree
    type(hsd_table), pointer :: node

    !> container for electron-phonon parameters
    type(TElPh), intent(inout) :: elph

    !> Geometry type
    type(TGeometry), intent(in) :: geom

    !> Orbitals infos
    type(TOrbitals), intent(in) :: orb

    !> Transport parameter type
    type(TTransPar), intent(in) :: tp


    logical :: block_model, semilocal_model
    integer :: stat

    elph%defined = .true.
    !! Only local el-ph model is defined (elastic for now)
    elph%model = 1

    call hsd_get_or_set(node, "MaxSCBAIterations", elph%scba_niter, 100)
    call hsd_get_or_set(node, "atomBlock", block_model, .false.)
    if (block_model) then
      elph%model = 2
    endif

    !BUG: semilocal model crashes because of access of S before its allocation
    !     this because initDephasing was moved into initprogram
    call hsd_get_or_set(node, "semiLocal", semilocal_model, .false.)
    if (semilocal_model) then
      call dftbp_error(node, "semilocal dephasing causes crash and has been "//&
           & "temporarily disabled")
      elph%model = 3
    endif

    call readCoupling(node, elph, geom, orb, tp)

  end subroutine readElPh


  !> Read Buettiker probe dephasing blocks (for density and/or current calculation)
  subroutine readDephasingBP(node, elph, geom, orb, tp)

    !> Node in input document tree
    type(hsd_table), pointer :: node

    !> container for buttiker-probes parameters
    type(TElPh), intent(inout) :: elph

    !> Geometry type
    type(TGeometry), intent(in) :: geom

    !> Orbitals infos
    type(TOrbitals), intent(in) :: orb

    !> Transport parameter type
    type(TTransPar), intent(inout) :: tp

    logical :: block_model, semilocal_model
    character(len=:), allocatable :: model
    type(hsd_table), pointer :: dephModel
    integer :: stat

    call dftbp_error(node,"Buettiker probes are still under development")

    elph%defined = .true.
    call hsd_get_choice(node, "", model, dephModel, stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(node, "Missing required block content")

    select case(model)
    case("dephasingprobes")
      !! Currently only zeroCurrent condition is implemented
      !! This corresponds to elastic dephasing probes
      tp%tZeroCurrent=.true.
      !! Only local bp model is defined (elastic for now)
    case("voltageprobes")
      call dftbp_error(dephModel,"voltageProbes have been not implemented yet")
      tp%tZeroCurrent=.false.
    case default
      call dftbp_error(dephModel,"unknown model")
    end select

    elph%model = 1

    call hsd_get_or_set(dephModel, "MaxSCBAIterations", elph%scba_niter, 100)

    call hsd_get_or_set(dephModel, "atomBlock", block_model, .false.)
    if (block_model) then
      elph%model = 2
    endif

    !BUG: semilocal model crashes because of access of S before its allocation
    !     this because initDephasing occurs in initprogram
    call hsd_get_or_set(dephModel, "semiLocal", semilocal_model, .false.)
    if (semilocal_model) then
      call dftbp_error(dephModel, "semilocal dephasing is not working yet")
      elph%model = 3
    endif

    call readCoupling(dephModel, elph, geom, orb, tp)

  end subroutine readDephasingBP


  !> Reads coupling strength and mode for dephasing
  !> 2 modes support, constant or specified per each orbital
  subroutine readCoupling(node, elph, geom, orb, tp)

    !> Node in the input tree
    type(hsd_table), pointer :: node

    !> container for buttiker-probes parameters
    type(TElPh), intent(inout) :: elph

    !> Geometry type
    type(TGeometry), intent(in) :: geom

    !> Orbitals infos
    type(TOrbitals), intent(in) :: orb

    !> Transport parameter type
    type(TTransPar), intent(in) :: tp

    character(len=:), allocatable :: buffer, method, modifier, modifier2
    type(hsd_table), pointer :: val, child, child2, child3, child4, field
    type(hsd_table_ptr), allocatable :: children(:)
    integer :: norbs, ii, jj, iAt
    integer :: atm_range(2)
    real(dp) :: rTmp
    integer, allocatable :: tmpI1(:)
    real(dp), allocatable :: atmCoupling(:)
    integer :: stat

    !! Allocate coupling array
    norbs = 0
    if (tp%defined) then
      atm_range(1) = tp%idxdevice(1)
      atm_range(2) = tp%idxdevice(2)
    else
      atm_range(1) = 1
      atm_range(2) = geom%nAtom
    endif
    do ii=atm_range(1), atm_range(2)
      norbs = norbs + orb%nOrbAtom(ii)
    enddo
    allocate(elph%coupling(norbs))
    elph%coupling(:) = 0.d0

    elph%orbsperatm = orb%nOrbAtom(atm_range(1):atm_range(2))

    call hsd_get_table(node, "Coupling", child, stat, auto_wrap=.true.)
    if (.not. associated(child)) call dftbp_error(node, "Missing required block: 'Coupling'")
    call hsd_get_attrib(node, "Coupling", modifier)
    call hsd_get_choice(child, "", method, val, stat)
    if (stat /= HSD_STAT_OK) call dftbp_error(child, "Missing coupling method")

    ! This reads also things like:  "Coupling [eV] = 0.34"
    !if (is_numeric(method)) then
    !  call getChildValue(node, "Coupling", rTmp, child=field)
    !  call convertUnitHsd(modifier, energyUnits, field, rTmp)
    !  elph%coupling = rTmp
    !  return
    !end if

    select case (method)
    case ("allorbitals")
      call hsd_get_table(child, "AllOrbitals", child2, stat, auto_wrap=.true.)
      call hsd_get(child2, "#text", elph%coupling, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(child2, "Missing required array: coupling values")
      call hsd_get_table(child2, "#text", field, stat, auto_wrap=.true.)
      if (.not. associated(field)) field => child2
      call convertUnitHsd(modifier, energyUnits, field, elph%coupling)

    case ("atomcoupling")
      call hsd_get_table(child, "AtomCoupling", child2, stat, auto_wrap=.true.)
      allocate(atmCoupling(atm_range(2)-atm_range(1)+1))
      atmCoupling = 0.d0
      call hsd_get_child_tables(child2, "AtomList", children)
      do ii = 1, size(children)
        child3 => children(ii)%ptr
        call hsd_get(child3, "Atoms", buffer, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(child3, "Missing required value: 'Atoms'")
        call hsd_get_table(child3, "Atoms", child4, stat, auto_wrap=.true.)
        if (.not. associated(child4)) child4 => child3
        call getSelectedAtomIndices(child4, buffer, geom%speciesNames, geom%species, tmpI1)
        call hsd_get(child3, "Value", rTmp, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(child3, "Missing required value: 'Value'")
        call hsd_get_attrib(child3, "Value", modifier2)
        call hsd_get_table(child3, "Value", field, stat, auto_wrap=.true.)
        if (.not. associated(field)) field => child3
        ! If not defined, use common unit modifier defined after Coupling
        if (len(modifier2)==0) then
          call convertUnitHsd(modifier, energyUnits, field, rTmp)
        else
          call convertUnitHsd(modifier2, energyUnits, field, rTmp)
        end if
        do jj=1, size(tmpI1)
          iAt = tmpI1(jj)
          if (atmCoupling(iAt) /= 0.0_dp) then
            call dftbp_warning(child3, "Previous setting of coupling &
                &for atom" // i2c(iAt) // " has been overwritten")
          end if
          atmCoupling(iAt) = rTmp
        end do
      end do

      ! Transform atom coupling in orbital coupling
      norbs = 0
      do ii=atm_range(1), atm_range(2)
        elph%coupling(norbs + 1:norbs + orb%nOrbAtom(ii)) = atmCoupling(ii)
        norbs = norbs + orb%nOrbAtom(ii)
      enddo
      deallocate(atmCoupling)

    case ("constant")
      call hsd_get(child, "Constant", rtmp, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(child, "Missing required value: 'Constant'")
      call hsd_get_table(child, "Constant", field, stat, auto_wrap=.true.)
      if (.not. associated(field)) field => child
      call convertUnitHsd(modifier, energyUnits, field, rTmp)
      elph%coupling = rTmp

    case default
      call dftbp_error(node, "Coupling definition unknown")
    end select

  end subroutine readCoupling


  !> Read Tunneling and Dos options from analysis block
  subroutine readTunAndDos(root, orb, geo, tundos, transpar, tempElec)
    type(hsd_table), pointer :: root
    type(TOrbitals), intent(in) :: orb
    type(TGeometry), intent(in) :: geo

    !> tundos is the container to be filled
    type(TNEGFTunDos), intent(inout) :: tundos
    type(TTransPar), intent(inout) :: transpar
    real(dp), intent(in) :: tempElec

    type(hsd_table), pointer :: pTmp, pNode, field
    type(hsd_table_ptr), allocatable :: pNodeList(:)
    integer :: ii, jj, ind, ncont, nKT
    real(dp) :: eRange(2), eRangeDefault(2)
    character(len=:), allocatable :: modifier
    type(TWrappedInt1), allocatable :: iAtInRegion(:)
    logical, allocatable :: tShellResInRegion(:)
    character(lc), allocatable :: regionLabelPrefixes(:)
    real(dp), allocatable :: realArr(:)
    integer :: stat

    tundos%defined = .true.

    ! ncont is needed for contact option allocation
    ncont = transpar%ncont

    call hsd_get_or_set(root, "Verbosity", tundos%verbose, 51)
    call hsd_get_or_set(root, "WriteLDOS", tundos%writeLDOS, .true.)
    call hsd_get_or_set(root, "WriteTunn", tundos%writeTunn, .true.)

    ! Read Temperature. Can override contact definition
    allocate(tundos%kbT(ncont))
    call hsd_get_table(root, "ContactTemperature", pTmp, stat, auto_wrap=.true.)
    if (associated(pTmp)) call hsd_get_attrib(root, "ContactTemperature", modifier)
    if (associated(pTmp)) then
      call hsd_get(pTmp, "#text", realArr, stat=stat)
      if (stat /= 0) call dftbp_error(root, "Error reading ContactTemperature")
      if (size(realArr) /= ncont) then
        call dftbp_error(root, "ContactTemperature does not match the number of contacts")
      end if
      tundos%kbT(:) = realArr
      call convertUnitHsd(modifier, energyUnits, pTmp, tundos%kbT)
    else
      do ii = 1, ncont
        if (transpar%contacts(ii)%kbT >= 0) then
          tundos%kbT(ii) = transpar%contacts(ii)%kbT
        else
          tundos%kbT(ii) = tempElec
        end if
      end do
    end if

    ! Parsing of energy range
    ! If the calculation is in equilibrium (all potentials to 0.0)
    ! then an energy range and step must be specified (it is assumed
    ! that the user use this filed to calculate a DOS or T(E) )
    ! If the calculation is out of equilibrium, a default similar to
    ! GreensFunction RealAxisStep is set to ensure that the current
    ! can be calculated without manually specify the energy parameters.

    if (all(transpar%contacts(:)%potential.eq.0.0)) then
      ! No default meaningful
      block
        real(dp), allocatable :: tmpRA(:)
        call hsd_get(root, "EnergyRange", tmpRA, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(root, "Missing required array: 'EnergyRange'")
        eRange(1:2) = tmpRA(1:2)
      end block
      call hsd_get_attrib(root, "EnergyRange", modifier)
      call hsd_get_table(root, "EnergyRange", field, stat, auto_wrap=.true.)
      if (.not. associated(field)) field => root
      call convertUnitHsd(modifier, energyUnits, field, eRange)
      call hsd_get(root, "EnergyStep", tundos%estep, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(root, "Missing required value: 'EnergyStep'")
      call hsd_get_attrib(root, "EnergyStep", modifier)
      call hsd_get_table(root, "EnergyStep", field, stat, auto_wrap=.true.)
      if (.not. associated(field)) field => root
      call convertUnitHsd(modifier, energyUnits, field, tundos%estep)
    else
      ! Default meaningful
      ! nKT is set to GreensFunction default, i.e. 10
      ! I avoid an explicit nKT option because I find it confusing here
      ! (it makes sense only out of equilibrium)
      ! Emin = min(-mu); Emax=max(-mu) where mu is Vi-min(Efi)
      ! Note: if Efi != min(Efi) a built in potential is added in poisson
      ! to aling the leads, we don't need to include it here
      nKT = 10
      eRangeDefault(1) = minval(-1.0*transpar%contacts(:)%potential) + &
                        & minval(1.0*transpar%contacts(:)%eFermi(1)) -   &
                        & nKT * maxval(tundos%kbT)
      eRangeDefault(2) = maxval(-1.0*transpar%contacts(:)%potential) + &
                        & minval(transpar%contacts(:)%eFermi(1)) +   &
                        & nKT * maxval(tundos%kbT)
      call hsd_get_or_set(root, "EnergyStep", tundos%estep, 6.65e-4_dp, child=field)
      call hsd_get_attrib(root, "EnergyStep", modifier)
      call convertUnitHsd(modifier, energyUnits, field, tundos%estep)
      block
        real(dp), allocatable :: tmpRA(:)
        call hsd_get_or_set(root, "EnergyRange", tmpRA, eRangeDefault)
        eRange(1:2) = tmpRA(1:min(size(tmpRA),2))
      end block
      call hsd_get_attrib(root, "EnergyRange", modifier)
      call hsd_get_table(root, "EnergyRange", field, stat, auto_wrap=.true.)
      if (.not. associated(field)) field => root
      call convertUnitHsd(modifier, energyUnits, field, eRange)
    end if

    tundos%emin = eRange(1)
    tundos%emax = eRange(2)
    ! Terminal currents
    call hsd_get_table(root, "TerminalCurrents", pTmp, stat, auto_wrap=.true.)
      if (associated(pTmp)) then
        call hsd_get_child_tables(pTmp, "EmitterCollector", pNodeList)
        allocate(tundos%ni(size(pNodeList)))
        allocate(tundos%nf(size(pNodeList)))
        do ii = 1, size(pNodeList)
          pNode => pNodeList(ii)%ptr
          call getEmitterCollectorByName(pNode, tundos%ni(ii), tundos%nf(ii),&
              & transpar%contacts(:)%name)
        end do
      else
        allocate(tundos%ni(ncont-1) )
        allocate(tundos%nf(ncont-1) )
        block
          type(hsd_table) :: newTbl
          call new_table(newTbl, name="TerminalCurrents")
          call root%add_child(newTbl)
        end block
        call hsd_get_table(root, "TerminalCurrents", pTmp, stat)
        ind = 1
        do ii = 1, 1
          do jj = ii + 1, ncont
            call hsd_set(pTmp, "EmitterCollector", &
                &(/ transpar%contacts(ii)%name, transpar%contacts(jj)%name /))
            tundos%ni(ind) = ii
            tundos%nf(ind) = jj
            ind = ind + 1
          end do
        end do
      end if
      call hsd_get_or_set(root, "Delta", tundos%delta, 1.0e-5_dp, child=field)
      call hsd_get_attrib(root, "Delta", modifier)
      call convertUnitHsd(modifier, energyUnits, field, &
          &tundos%delta)
      call hsd_get_or_set(root, "BroadeningDelta", tundos%broadeningDelta, 0.0_dp, child=field)
      call hsd_get_attrib(root, "BroadeningDelta", modifier)
      call convertUnitHsd(modifier, energyUnits, field, &
          &tundos%broadeningDelta)

      call readPDOSRegions(root, geo, transpar%idxdevice, iAtInRegion, &
          & tShellResInRegion, regionLabelPrefixes)

      if (allocated(iAtInRegion)) then
        call transformPdosRegionInfo(iAtInRegion, tShellResInRegion, &
            & regionLabelPrefixes, orb, geo%species, tundos%dosOrbitals, &
            & tundos%dosLabels)
      end if

  end subroutine readTunAndDos


  !> Read bias information, used in Analysis and Green's function eigensolver
  subroutine readContacts(pNodeList, contacts, geom, task, contactLayerTol)

    !> Node to process
    type(hsd_table_ptr), intent(in) :: pNodeList(:)

    !> Contacts
    type(ContactInfo), allocatable, dimension(:), intent(inout) :: contacts

    !> Geometry of the system
    type(TGeometry), intent(in) :: geom

    !> What type of transport-related calculation is this?
    character(*), intent(in) :: task

    !> Tolerance to distortion of contact vectors
    real(dp), intent(out) :: contactLayerTol

    integer :: ii, stat
    type(hsd_table), pointer :: field, pNode, pTmp, child1, child2
    character(len=:), allocatable :: buffer, modifier
    real(dp), allocatable :: realArr(:)

    do ii = 1, size(contacts)

      contacts(ii)%wideBand = .false.
      contacts(ii)%wideBandDos = 0.0_dp

      pNode => pNodeList(ii)%ptr
      call hsd_get(pNode, "Id", buffer, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(pNode, "Missing required value: 'Id'")
      call hsd_get_table(pNode, "Id", pTmp, stat, auto_wrap=.true.)
      if (.not. associated(pTmp)) pTmp => pNode
      buffer = tolower(trim(unquote(buffer)))
      if (len(buffer) > mc) then
        call dftbp_error(pTmp, "Contact id may not be longer than " // i2c(mc) // " characters.")
      end if
      contacts(ii)%name = buffer
      if (any(contacts(1:ii-1)%name == contacts(ii)%name)) then
        call dftbp_error(pTmp, "Contact id '" // trim(contacts(ii)%name) //  "' already in use")
      end if

      call hsd_get_or_set(pNode, "PLShiftTolerance", contactLayerTol, 1e-5_dp, child=field)
      call hsd_get_attrib(pNode, "PLShiftTolerance", modifier)
      call convertUnitHsd(modifier, lengthUnits, field, contactLayerTol)

      block
        integer, allocatable :: tmpArr(:)
        call hsd_get(pNode, "AtomRange", tmpArr, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(pNode, "Missing required array: 'AtomRange'")
        contacts(ii)%idxrange(1:2) = tmpArr(1:2)
      end block
      call hsd_get_table(pNode, "AtomRange", pTmp, stat, auto_wrap=.true.)
      if (.not. associated(pTmp)) pTmp => pNode
      call getContactVector(contacts(ii)%idxrange, geom, ii, contacts(ii)%name, pTmp,&
        & contactLayerTol, contacts(ii)%lattice, contacts(ii)%dir)
      contacts(ii)%length = sqrt(sum(contacts(ii)%lattice**2))

      ! Contact temperatures. A negative default is used so it is quite clear when the user sets a
      ! different value. In such a case this overrides values defined in the Filling block
      call hsd_get_table(pNode, "Temperature", field, stat, auto_wrap=.true.)
      if (associated(field)) then
        call hsd_get_attrib(pNode, "Temperature", modifier)
        call hsd_get_or_set(pNode, "Temperature", contacts(ii)%kbT, 0.0_dp, child=field)
        call convertUnitHsd(modifier, energyUnits, field, contacts(ii)%kbT)
      else
        contacts(ii)%kbT = -1.0_dp ! -1.0 simply means 'not defined'
      end if

      if (task .eq. "uploadcontacts") then
        call hsd_get_or_set(pNode, "Potential", contacts(ii)%potential, 0.0_dp, child=field)
        call hsd_get_attrib(pNode, "Potential", modifier)
        call convertUnitHsd(modifier, energyUnits, field, contacts(ii)%potential)

        call hsd_get_or_set(pNode, "WideBand", contacts(ii)%wideBand, .false.)

        if (contacts(ii)%wideBand) then

          ! WideBandApproximation is defined as energy spacing between levels of the contact. In the
          ! code the inverse value (Density of states) is used. Convert the negf input
          ! value. Default is 20 / e eV.
          call hsd_get_or_set(pNode, "LevelSpacing", contacts(ii)%wideBandDos, 0.735_dp,&
              & child=field)
          call hsd_get_attrib(pNode, "LevelSpacing", modifier)
          call convertUnitHsd(modifier, energyUnits, field, contacts(ii)%wideBandDos)
          contacts(ii)%wideBandDos = 1.d0 / contacts(ii)%wideBandDos

        end if


        ! Fermi level: in case of collinear spin we accept two values (up and down)
        ! call init(fermiBuffer)
        ! call getChildValue(pNode, "FermiLevel", fermiBuffer, modifier=modifier)
        ! if ( len(fermiBuffer) .eq. 1) then
        !   call asArray(fermiBuffer, contacts(ii)%eFermi)
        !   contacts(ii)%eFermi(2) = contacts(ii)%eFermi(1)
        ! else if ( len(fermiBuffer) .eq. 2) then
        !   call asArray(fermiBuffer, contacts(ii)%eFermi)
        ! else
        !   call dftbp_error(pNode, "FermiLevel accepts 1 or 2 (for collinear spin) values")
        ! end if
        ! call destruct(fermiBuffer)


        call hsd_get_table(pNode, "FermiLevel", child2, stat, auto_wrap=.true.)
        if (.not. associated(child2)) then
          block
            type(hsd_table) :: emptyTbl
            call new_table(emptyTbl, name="fermilevel")
            call pNode%add_child(emptyTbl)
          end block
          call hsd_get_table(pNode, "FermiLevel", child2, stat, auto_wrap=.true.)
          buffer = ""
          child1 => null()
        else
          call hsd_get_attrib(pNode, "FermiLevel", modifier)
          call hsd_get_choice(child2, "", buffer, child1, stat)
          if (stat /= HSD_STAT_OK) buffer = ""
        end if
        if (buffer == "" .and. .not. hasInlineData(child2)) then
          contacts(ii)%tFermiSet = .false.
          call dftbp_warning(pNode, "Missing Fermi level - required to be set in solver block or&
              & read from a contact shift file")
        else
          call hsd_get(child2, "#text", realArr, stat=stat)
          if (stat /= 0) call dftbp_error(pNode, "Error reading FermiLevel")
          select case(size(realArr))
          case (1)
            contacts(ii)%eFermi(1) = realArr(1)
            contacts(ii)%eFermi(2) = realArr(1)
          case (2)
            contacts(ii)%eFermi(1:2) = realArr(1:2)
          case default
            call dftbp_error(pNode, "FermiLevel accepts 1 or 2 (for collinear spin) values")
          end select
          call convertUnitHsd(modifier, energyUnits, child2, contacts(ii)%eFermi)

          contacts(ii)%tFermiSet = .true.

          ! NOTE: These options have been commented out: there is a problem in parallel execution
          ! since one single file is accessed by all processors causing rush conditions
          ! The options are therefore disabled for the official dftb+ release
          !call getChildValue(pNode, "WriteSelfEnergy", contacts(ii)%tWriteSelfEnergy, .false.)
          !call getChildValue(pNode, "WriteSurfaceGF", contacts(ii)%tWriteSurfaceGF, .false.)
          !call getChildValue(pNode, "ReadSelfEnergy", contacts(ii)%tReadSelfEnergy, .false.)
          !call getChildValue(pNode, "ReadSurfaceGF", contacts(ii)%tReadSurfaceGF, .false.)
          contacts(ii)%tWriteSelfEnergy = .false.
          contacts(ii)%tWriteSurfaceGF = .false.
          contacts(ii)%tReadSelfEnergy = .false.
          contacts(ii)%tReadSurfaceGF = .false.
        end if

      end if

    end do

  end subroutine readContacts


  !> Read in Fermi levels
  subroutine getFermiLevels(pNode, eFermis, nodeModifier)

    !> Document tree node to start from
    type(hsd_table), pointer :: pNode

    !> Fermi energies for contacts
    real(dp), intent(out) :: eFermis(:)

    !> Any node modifiers in action
    character(len=*), intent(in) :: nodeModifier

    real(dp) :: eFermi
    type(hsd_table), pointer :: pChild
    character(len=:), allocatable :: modifier
    integer :: stat

    call hsd_get_table(pNode, "SetForAll", pChild, stat, auto_wrap=.true.)
    if (associated(pChild)) then
      call hsd_get(pChild, "#text", eFermi, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(pChild, "Missing required value in 'SetForAll'")
      call convertUnitHsd(nodeModifier, energyUnits, pNode, eFermi)
      eFermis(:) = eFermi
    else
      block
        real(dp), allocatable :: tmpRA(:)
        call hsd_get(pNode, "#text", tmpRA, stat=stat)
        if (stat /= HSD_STAT_OK) call dftbp_error(pNode, "Missing required array: Fermi levels")
        eFermis(1:min(size(tmpRA),size(eFermis))) = tmpRA(1:min(size(tmpRA),size(eFermis)))
      end block
      call hsd_get_attrib(pNode, "#text", modifier)
      call convertUnitHsd(modifier, energyUnits, pNode, eFermis)
    end if

  end subroutine getFermiLevels


  !> Get contacts for terminal currents by name
  subroutine getEmitterCollectorByName(pNode, emitter, collector, contactNames)

    !> Node in the input tree for error reporting
    type(hsd_table), pointer :: pNode

    !> Contact number for emitting
    integer, intent(out) :: emitter

    !> Contact number for collecting
    integer, intent(out) :: collector

    !> Labels of contacts
    character(len=*), intent(in) :: contactNames(:)

    character(:), allocatable :: strArr(:)
    integer :: stat

    call hsd_get(pNode, "#text", strArr, stat=stat)
    if (stat /= 0 .or. size(strArr) /= 2) then
      call dftbp_error(pNode, "You must provide two contacts")
    end if
    emitter = getContactByName(contactNames, strArr(1), pNode)
    collector = getContactByName(contactNames, strArr(2), pNode)

  end subroutine getEmitterCollectorByName


  !> Getting the contact by name
  function getContactByName(contactNames, contName, pNode) result(contact)

    !> Node in the input tree for error reporting
    type(hsd_table), pointer :: pNode

    !> All of the contact labels
    character(len=*), intent(in) :: contactNames(:)

    !> Specific contact label to identify
    character(len=*), intent(in) :: contName

    !> Contact number
    integer :: contact

    logical :: tFound

    tFound = .false.
    do contact = 1, size(contactNames)
      tFound = (contactNames(contact) == contName)
      if (tFound) then
        exit
      end if
    end do
    if (.not. tFound) then
      call dftbp_error(pNode, "Invalid collector contact name '" // trim(contName) // "'")
    end if

  end function getContactByName


  !> Read the names of regions to calculate PDOS for
  subroutine readPDOSRegions(node, geo, idxdevice, iAtInregion, tShellResInRegion, regionLabels)

    !> Node to be parsed
    type(hsd_table), pointer, intent(in) :: node

    !> Geometry of the system
    type(TGeometry), intent(in) :: geo

    !> Is the region to be projected by shell
    integer, intent(in) :: idxdevice(2)

    !> Atoms in a given region
    type(TWrappedInt1), allocatable, intent(out) :: iAtInRegion(:)

    !> Is the region to be projected by shell
    logical, allocatable, intent(out) :: tShellResInRegion(:)

    !> Labels for the regions
    character(lc), allocatable, intent(out) :: regionLabels(:)

    integer :: nReg, iReg
    integer, allocatable :: tmpI1(:)
    type(hsd_table_ptr), allocatable :: children(:)
    type(hsd_table), pointer :: child, child2
    character(len=:), allocatable :: buffer
    character(lc) :: strTmp
    logical :: do_ldos
    integer :: stat

    call hsd_get_child_tables(node, "Region", children)
    nReg = size(children)

    if (nReg == 0) then
      call hsd_get_or_set(node, "ComputeLDOS", do_ldos, .true.)
      if (do_ldos) then
        write(strTmp,"(I0, ':', I0)") idxdevice(1), idxdevice(2)
        block
          type(hsd_table) :: newTbl
          call new_table(newTbl, name="Region")
          call node%add_child(newTbl)
        end block
        call hsd_get_table(node, "Region", child, stat)
        call hsd_set(child, "Atoms", trim(strTmp))
        call hsd_set(child, "Label", "localDOS")
        call hsd_get_child_tables(node, "Region", children)
        nReg = size(children)
      else
        return
      end if
    end if

    allocate(tShellResInRegion(nReg))
    allocate(regionLabels(nReg))
    allocate(iAtInRegion(nReg))
    do iReg = 1, nReg
      child => children(iReg)%ptr
      call hsd_get(child, "Atoms", buffer, stat=stat)
      if (stat /= HSD_STAT_OK) call dftbp_error(child, "Missing required value: 'Atoms'")
      call hsd_get_table(child, "Atoms", child2, stat, auto_wrap=.true.)
      if (.not. associated(child2)) child2 => child
      call getSelectedAtomIndices(child2, buffer, geo%speciesNames,&
          & geo%species(idxdevice(1) : idxdevice(2)), tmpI1,&
          & selectionRange=[idxdevice(1), idxdevice(2)], indexRange=[1, geo%nAtom])
      iAtInRegion(iReg)%data = tmpI1
      call hsd_get_or_set(child, "ShellResolved", tShellResInRegion(iReg), .false., child=child2)
      if (tShellResInRegion(iReg)) then
        if (.not. all(geo%species(tmpI1) == geo%species(tmpI1(1)))) then
          call dftbp_error(child2, "Shell resolved PDOS can only summed up over atoms of the same&
              & type")
        end if
      end if
      write(strTmp, "('region',I0)") iReg
      call hsd_get_or_set(child, "Label", buffer, trim(strTmp))
      regionLabels(iReg) = unquote(buffer)
    end do

  end subroutine readPDOSRegions


  !> Some assignment and consistency check in negf/poisson containers before calling initialization
  subroutine finalizeNegf(input)

    !> Input structure for DFTB+
    type(TInputData), intent(inout) :: input

    integer :: ii

    !! Check consistency between different deltas
    if (input%ginfo%tundos%defined.and.input%ginfo%greendens%defined) then
      if (input%ginfo%tundos%delta.ne.input%ginfo%greendens%delta) then
        call error("Delta parameter must be the same in GreensFunction and TunnelingAndDos")
      end if
    end if

    !! Assign spin degeneracy to every block which may use it
    if (input%ginfo%tundos%defined) then
      if (input%ctrl%tSpin) input%ginfo%tundos%gSpin = 1
      if (.not.input%ctrl%tSpin) input%ginfo%tundos%gSpin = 2
    end if
    if (input%ginfo%greendens%defined) then
      if (input%ctrl%tSpin) input%ginfo%greendens%gSpin = 1
      if (.not.input%ctrl%tSpin) input%ginfo%greendens%gSpin = 2
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !! Inheritance of first layer indexes to green solver when transport is defined
    if (input%transpar%defined .and. input%ginfo%greendens%defined) then
      input%ginfo%greendens%nPLs = input%transpar%nPLs
      input%ginfo%greendens%PL = input%transpar%PL
    end if

    #:block REQUIRES_COMPONENT('Poisson-solver', WITH_POISSON)
      !! Not orthogonal directions in transport are only allowed if no Poisson
      if (input%poisson%defined.and.input%transpar%defined) then
        do ii = 1, input%transpar%ncont
          ! If dir is  any value but x,y,z (1,2,3) it is considered oriented along
          ! a direction not parallel to any coordinate axis
          if (input%transpar%contacts(ii)%dir.lt.1 .or. &
            &input%transpar%contacts(ii)%dir.gt.3 ) then
            call error("Contact " // i2c(ii) // " not parallel to any &
              & coordinate axis and is not compatible with Poisson solver")
          end if
        end do
      end if
    #:endblock

    !! Temporarily not supporting surface green function read/load
    !! for spin polarised, because spin is handled outside of libnegf
    if (input%ginfo%greendens%defined) then
      if (input%ctrl%tSpin .and. input%ginfo%greendens%saveSGF) then
        call error("SaveSurfaceGFs must be disabled in collinear spin calculations")
      end if
      if  (input%ctrl%tSpin .and. input%ginfo%greendens%readSGF) then
        call error("ReadSurfaceGFs must be disabled in collinear spin calculations")
      end if
    end if

  end subroutine finalizeNegf
#:endif


end module dftbp_dftbplus_parser_transport
