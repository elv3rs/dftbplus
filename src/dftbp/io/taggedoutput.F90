!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

#! Template parameters: (data type, suffix name, data type name in the tagged output, format string)
#:set TEMPLATE_PARAMS = [('real(dp)', 'Real', 'real', 'formReal'),&
    & ('complex(dp)', 'Cplx', 'complex', 'formCmplx'),&
    & ('integer', 'Integer', 'integer', 'formInt'),&
    & ('logical', 'Logical', 'logical', 'formLogical')]

#! Maximal rank to include into the interface (string from 0 - scalar)
#:set MAX_RANK = 5


!> Contains routines to write out various data structures in a comprehensive tagged format.
module dftbp_io_taggedoutput
  use dftbp_common_accuracy, only : dp
  use hsd_data, only : hsd_table, hsd_value, hsd_node, &
      & VALUE_TYPE_INTEGER, VALUE_TYPE_REAL, VALUE_TYPE_LOGICAL, VALUE_TYPE_COMPLEX
  use dftbp_common_file, only : TFileDescr, openFile, closeFile
  implicit none

  private
  public :: tagLabels
  public :: TTaggedWriter, TTaggedWriter_init
  public :: writeTreeAsTag


  !> Length of permissible tag labels. Tag names should be shorter than lenLabel!
  integer, parameter :: lenLabel = 20

  !> Max length of the format strings for individual items
  integer, parameter :: lenFormStr = 20


  !> Contains a writer to write data in tagged form.
  type :: TTaggedWriter
    private
    ! Format strings
    character(len=lenFormStr) :: formReal
    character(len=lenFormStr) :: formCmplx
    character(len=lenFormStr) :: formInt
    character(len=lenFormStr) :: formLogical
    logical :: initialized = .false.
  contains
    #:for _, SUFFIX, _, _ in TEMPLATE_PARAMS
      #:for RANK in range(MAX_RANK + 1)
        procedure, private :: write${SUFFIX}$${RANK}$ => TTaggedWriter_write${SUFFIX}$${RANK}$
        generic :: write => write${SUFFIX}$${RANK}$
      #:endfor
    #:endfor
  end type TTaggedWriter


  !> Enumeration of the possible tag labels
  type :: TTagLabelsEnum

    !> unit cell volume (periodic)
    character(lenLabel) :: volume = 'cell_volume'

    !> final geometry
    character(lenLabel) :: endCoord = 'end_coords'

    !> excitation energies in Casida formalism
    character(lenLabel) :: excEgy = 'exc_energies_sqr'

    !> excited state force contributions
    character(lenLabel) :: excForce = 'exc_forces'

    !> oscillator strength for excitations
    character(lenLabel) :: excOsc = 'exc_oscillator'

    !> Transition dipole moments for excitations
    character(lenLabel) :: excDipole = 'exc_transdip'

    !> Square of transition charges for target state
    character(lenLabel) :: transQ = 'transq_sqr'

    !> nonadiabatic coupling vector, H
    character(lenLabel) :: nacH = 'coupling_vectors'

    !> nonadiabatic coupling vector TD-DFTB
    character(lenLabel) :: nacv = 'nac_vectors'

    !> ground state total forces
    character(lenLabel) :: forceTot = 'forces'

    !> forces on any external charges present
    character(lenLabel) :: chrgForces = 'forces_ext_charges'

    !> Fermi level(s)
    character(lenLabel) :: fermiLvl = 'fermi_level'

    !> number of electrons
    character(lenLabel) :: nElec = 'number_of_electrons'

    !> eigenvalues/single particle states
    character(lenLabel) :: eigvals = 'eigenvalues'

    !> filling of the eigenstates
    character(lenLabel) :: eigFill = 'filling'

    !> Gibbs free energy for finite pressure periodic systems
    character(lenLabel) :: gibbsFree = 'gibbs_energy'

    !> Gross atomic charges
    character(lenLabel) :: qOutAtGross  = 'gross_atomic_charges'

    !> Charge model 5 corrected atomic gross charges
    character(lenLabel) :: qOutAtCM5 = 'cm5_atomic_charges'

    !> Gross atomic spin polarizations
    character(lenLabel) :: spinOutAtGross  = 'gross_atomic_spins'

    !> numerically calculated second derivatives matrix
    character(lenLabel) :: hessianNum = 'hessian_numerical'

    !> numerically calculated Born charges
    character(lenLabel) :: BorndDipNum = 'born_mudrv_numerical'

    !> final energy components after real-time propagation
    character(lenLabel) :: tdenergy = 'final_energy'

    !> final dipole moment vector after real-time propagation
    character(lenLabel) :: tddipole = 'final_dipole_moment'

    !> final negative gross atomic Mulliken charges after real-time propagation
    character(lenLabel) :: tdcharges = 'final_td_charges'

    !> final forces components after real-time (Ehrenfest) propagation
    character(lenLabel) :: ehrenforces = 'final_ehrenfest_forc'

    !> final geometry after real-time (Ehrenfest) propagation
    character(lenLabel) :: ehrencoords = 'final_ehrenfest_geom'

    !> final velocities after real-time (Ehrenfest) propagation
    character(lenLabel) :: ehrenvelos = 'final_ehrenfest_velo'

    !> final molecular orbitals occupations after real-time (Ehrenfest) propagation
    character(lenLabel) :: tdprojocc = 'final_td_proj_occ'

    !> Sum of bond populaion values (should be number of electrons)
    character(lenLabel) :: sumBondPopul = 'sum_bond_pops'

    !> final atom-resolved energies
    character(lenLabel) :: atomenergies = 'atomic_energies'

    !> total energy including electron TS contribution
    character(lenLabel) :: freeEgy = 'mermin_energy'

    !> Mulliken charges
    character(lenLabel) :: qOutput = 'orbital_charges'

    !> Pipek-Mezey localisation score of single particle levels
    character(lenLabel) :: pmlocalise = 'pm_localisation'

    !> total stress tensor for periodic geometries
    character(lenLabel) :: stressTot = 'stress'

    !> total tunneling vector
    character(lenLabel) :: tunn = 'total_tunneling'

    !> total projected DOS vector
    character(lenLabel) :: ldos = 'total_localdos'

    !> total bond currents
    character(lenLabel) :: localCurrents = 'local_currents'

    !> total internal energy
    character(lenLabel) :: egyTotal   = 'total_energy'

    !> total internal energy for averaged state in REKS
    character(lenLabel) :: egyAvg   = 'averaged_energy'

    !> total internal energy extrapolated to 0 K
    character(lenLabel) :: egy0Total   = 'extrapolated0_energy'

    !> Energy, which if differentiated gives - force
    character(lenLabel) :: egyForceRelated = 'forcerelated_energy'

    !> Internal electric field
    character(lenLabel) :: internField = 'internal_efield'

    !> External electric field
    character(lenLabel) :: externField = 'external_efield'

    !> Static electric polarizability from linear response/perturbation
    character(lenLabel) :: dmudEPerturb = 'staticPolResponse'

    !> Static gross charge (Mulliken) response from linear response/perturbation
    character(lenLabel) :: dqdEPerturb = 'staticChargeReponse'

    !> Derivatives of ground state single particle eigenvalues wrt. k
    character(lenLabel) :: dEigenDE = 'dEidEfield'

    !> Number of electrons at the Fermi energy
    character(lenLabel) :: neFermi = 'neFermi'

    !> Derivative of the Fermi energy with respect to electric field
    character(lenLabel) :: dEfdE = 'dEfdE'

    !> Derivatives of ground state single particle eigenvalues wrt. onsite potentials
    character(lenLabel) :: dEigenDVons = 'dEidVons'

    !> Derivatives of ground state single particle eigenvalues wrt. potential at an atom
    character(lenLabel) :: dEigenDV = 'dEidV'

    !> Static gross charge (Mulliken) response with respect to potential at an atom
    character(lenLabel) :: dqdV = 'dqdV'

    !> Static net charge (onsite) response with respect to potential at an atom
    character(lenLabel) :: dqnetdV = 'dqnetdV'

    !> Derivatives of gross atomic charges wrt. x
    character(lenLabel) :: dqdx = 'dqdx'

    !> Born effective charges
    character(lenLabel) :: borncharges = 'borncharges'

    !> two-electron addition/removal energies in ppRPA formalism
    character(lenLabel) :: egyppRPA = '2e_add-rem_energies'

    !> atomic masses
    character(lenLabel) :: atomMass = 'atomic_masses'

    !> Total dipole moment
    character(lenLabel) :: dipoleMoment = 'dipole_moments'

    !> Rescaled dipole moment (for example if solvated)
    character(lenLabel) :: scaledDipole = 'scaled_dipole'

    !> Atomic dipole moments
    character(lenLabel) :: dipoleAtom = 'atomic_dipole_moment'

  end type TTagLabelsEnum


  !> Enum containing the tag labels used.
  type(TTagLabelsEnum), parameter :: tagLabels = TTagLabelsEnum()


contains


  !> initialise writer
  subroutine TTaggedWriter_init(this)

    !> Instance
    type(TTaggedWriter), intent(out) :: this

    integer :: nDecDigit, nExpDigit, nChar, nField

    if (this%initialized) then
      return
    end if

    !! "-3.1234567E-123 ": nDec = 7, nExpDigit = 3, nChar = 16
    nExpDigit = ceiling(log(maxexponent(1.0_dp) / log(10.0)) / log(10.0))
    nDecDigit = precision(1.0_dp)
    nChar = nDecDigit + nExpDigit + 6
    nField = 80 / nChar
    if (nField == 0) then
      nField = 1
    end if

    write (this%formReal, "('(', I2.2, 'E', I2.2, '.', I2.2, 'E', I3.3, ')')") nField, nChar,&
        & nDecDigit, nExpDigit

    write (this%formCmplx, "('(', I2.2, '(2E', I2.2, '.', I2.2, 'E', I3.3, '))')") nField / 2,&
        & nChar, nDecDigit, nExpDigit

    nChar = digits(1) + 2
    nField = 80 / nChar
    if (nField == 0) then
      nField = 1
    end if
    write (this%formInt, "('(', I2.2, 'I', I2.2, ')')") nField, nChar
    write (this%formLogical, "('(40L2)')")

    this%initialized = .true.

  end subroutine TTaggedWriter_init


#:for DATA_TYPE, SUFFIX, DATA_TYPE_TAG_NAME, FORMAT_STRING in TEMPLATE_PARAMS
  #:for RANK in range(MAX_RANK + 1)

  !> Write tagged data (data type: ${DATA_TYPE}$)
  subroutine TTaggedWriter_write${SUFFIX}$${RANK}$(this, file, tag, data, optForm)

    !> Instance
    class(TTaggedWriter), intent(inout) :: this

    !> File ID
    integer, intent(in) :: file

    !> tag label
    character(len=*), intent(in) :: tag

    !> data to print
    ${DATA_TYPE}$, intent(in) :: data${FORTRAN_ARG_DIM_SUFFIX(RANK)}$

    !> optional formatting string
    character(len=*), optional, intent(in) :: optForm

    character(len=20) :: form

    @:ASSERT(this%initialized)

    if (present(optForm)) then
      form = getLabel(optForm)
    else
      form = getLabel(this%${FORMAT_STRING}$)
    end if
    #:if RANK
      call writeTaggedHeader(file, tag, '${DATA_TYPE_TAG_NAME}$', shape(data))
    #:else
      call writeTaggedHeader(file, tag, '${DATA_TYPE_TAG_NAME}$')
    #:endif
    write(file, form) data

  end subroutine TTaggedWriter_write${SUFFIX}$${RANK}$

  #:endfor
#:endfor


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Private functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !> Writes the tagged header.
  subroutine writeTaggedHeader(file, tag, dataType, dataShape)

    !> File id to write to
    integer, intent(in) :: file

    !> Tag name
    character(*), intent(in) :: tag

    !> Form string to use
    character(*), intent(in) :: dataType

    !> Original shape of the data
    integer, intent(in), optional :: dataShape(:)

    character(100) :: buffer

    if (present(dataShape)) then
      if (size(dataShape) == 1) then
        write(buffer, "(5A,I0,4A)") '("', getLabel(tag), ":", trim(dataType), ":", size(dataShape),&
            & ":", '",', 'I0', ')'
      else
        write(buffer, "(5A,I0,3A,I0,2A)") '("', getLabel(tag), ":", trim(dataType), ":",&
            & size(dataShape), ":", '",', 'I0,', size(dataShape) - 1, '(",",I0)', ')'
      end if
      write(file, buffer) dataShape
    else
      write(file, "(4A,I0,A)") getLabel(tag), ":", trim(dataType), ":", 0, ":"
    end if

  end subroutine writeTaggedHeader


  !> Extracts the label for a tag
  function getLabel(tag)

    !> relevant tag
    character(len=*), intent(in) :: tag

    !> Label
    character(len=20) :: getLabel

    integer :: lentrim

    lentrim = len_trim(tag)
    if (lentrim >= lenLabel) then
      getLabel(:) = tag(1:lenLabel)
    else
      getLabel(1:lentrim) = tag(1:lentrim)
      getLabel(lentrim+1:lenLabel) = " "
    end if

  end function getLabel


  !> Write an hsd_table tree to a file in legacy tag format.
  !!
  !! Each hsd_value child of root is serialized as:
  !!   tagname             :datatype:rank:dim1,dim2,...
  !!     <formatted data>
  !! Table children are skipped (they are structural containers).
  subroutine writeTreeAsTag(fileName, root)

    !> Output file name (opened in append mode)
    character(len=*), intent(in) :: fileName

    !> Root table whose children are written
    type(hsd_table), intent(in) :: root

    type(TFileDescr) :: fd
    class(hsd_node), pointer :: child
    character(len=lenFormStr) :: formReal, formCmplx, formInt, formLogical
    integer :: nDecDigit, nExpDigit, nChar, nField, ii

    ! Initialise format strings (same logic as TTaggedWriter_init)
    nExpDigit = ceiling(log(maxexponent(1.0_dp) / log(10.0)) / log(10.0))
    nDecDigit = precision(1.0_dp)
    nChar = nDecDigit + nExpDigit + 6
    nField = 80 / nChar
    if (nField == 0) nField = 1
    write(formReal, "('(', I2.2, 'E', I2.2, '.', I2.2, 'E', I3.3, ')')") nField, nChar,&
        & nDecDigit, nExpDigit
    write(formCmplx, "('(', I2.2, '(2E', I2.2, '.', I2.2, 'E', I3.3, '))')") nField / 2,&
        & nChar, nDecDigit, nExpDigit
    nChar = digits(1) + 2
    nField = 80 / nChar
    if (nField == 0) nField = 1
    write(formInt, "('(', I2.2, 'I', I2.2, ')')") nField, nChar
    write(formLogical, "('(40L2)')")

    call openFile(fd, fileName, mode="a")

    do ii = 1, root%num_children
      call root%get_child(ii, child)
      if (.not. associated(child)) cycle
      select type (child)
      type is (hsd_value)
        call writeOneTagValue(fd%unit, child, formReal, formCmplx, formInt, formLogical)
      end select
    end do

    call closeFile(fd)

  end subroutine writeTreeAsTag


  !> Write a single hsd_value node in tag format to an open file unit.
  subroutine writeOneTagValue(unit, val, formReal, formCmplx, formInt, formLogical)

    !> File unit (already open)
    integer, intent(in) :: unit

    !> Value node to serialise
    type(hsd_value), intent(in) :: val

    !> Format strings for the various data types
    character(len=*), intent(in) :: formReal, formCmplx, formInt, formLogical

    integer :: rank
    integer, allocatable :: dataShape(:)

    select case (val%value_type)

    case (VALUE_TYPE_REAL)
      if (allocated(val%real_matrix)) then
        call parseShapeAttrib(val, rank, dataShape)
        if (.not. allocated(dataShape)) then
          allocate(dataShape(2))
          dataShape = shape(val%real_matrix)
        end if
        call writeTaggedHeader(unit, val%name, 'real', dataShape)
        write(unit, formReal) val%real_matrix
      else if (allocated(val%real_array)) then
        call parseShapeAttrib(val, rank, dataShape)
        if (.not. allocated(dataShape)) then
          allocate(dataShape(1))
          dataShape(1) = size(val%real_array)
        end if
        call writeTaggedHeader(unit, val%name, 'real', dataShape)
        write(unit, formReal) val%real_array
      else
        call writeTaggedHeader(unit, val%name, 'real')
        write(unit, formReal) val%real_value
      end if

    case (VALUE_TYPE_INTEGER)
      if (allocated(val%int_matrix)) then
        call parseShapeAttrib(val, rank, dataShape)
        if (.not. allocated(dataShape)) then
          allocate(dataShape(2))
          dataShape = shape(val%int_matrix)
        end if
        call writeTaggedHeader(unit, val%name, 'integer', dataShape)
        write(unit, formInt) val%int_matrix
      else if (allocated(val%int_array)) then
        call parseShapeAttrib(val, rank, dataShape)
        if (.not. allocated(dataShape)) then
          allocate(dataShape(1))
          dataShape(1) = size(val%int_array)
        end if
        call writeTaggedHeader(unit, val%name, 'integer', dataShape)
        write(unit, formInt) val%int_array
      else
        call writeTaggedHeader(unit, val%name, 'integer')
        write(unit, formInt) val%int_value
      end if

    case (VALUE_TYPE_LOGICAL)
      if (allocated(val%logical_array)) then
        call parseShapeAttrib(val, rank, dataShape)
        if (.not. allocated(dataShape)) then
          allocate(dataShape(1))
          dataShape(1) = size(val%logical_array)
        end if
        call writeTaggedHeader(unit, val%name, 'logical', dataShape)
        write(unit, formLogical) val%logical_array
      else
        call writeTaggedHeader(unit, val%name, 'logical')
        write(unit, formLogical) val%logical_value
      end if

    case (VALUE_TYPE_COMPLEX)
      if (allocated(val%complex_array)) then
        call parseShapeAttrib(val, rank, dataShape)
        if (.not. allocated(dataShape)) then
          allocate(dataShape(1))
          dataShape(1) = size(val%complex_array)
        end if
        call writeTaggedHeader(unit, val%name, 'complex', dataShape)
        write(unit, formCmplx) val%complex_array
      else
        call writeTaggedHeader(unit, val%name, 'complex')
        write(unit, formCmplx) val%complex_value
      end if

    end select

  end subroutine writeOneTagValue


  !> Parse rank and shape from the attrib field of an hsd_value.
  !!
  !! The expected format is "rank:dim1,dim2,...", e.g. "3:5,3,2".
  !! If the attrib is absent or unparseable, dataShape is left unallocated.
  subroutine parseShapeAttrib(val, rank, dataShape)

    !> Value node whose attrib field is examined
    type(hsd_value), intent(in) :: val

    !> Parsed rank (0 if not parseable)
    integer, intent(out) :: rank

    !> Parsed shape dimensions (unallocated on failure)
    integer, allocatable, intent(out) :: dataShape(:)

    integer :: colonPos, nDims, pos, commaPos, ii
    character(len=256) :: buf

    rank = 0
    if (.not. allocated(val%attrib)) return
    if (len_trim(val%attrib) == 0) return

    buf = val%attrib
    colonPos = index(buf, ':')
    if (colonPos == 0) return

    read(buf(1:colonPos-1), *, err=999, end=999) rank
    if (rank == 0) return

    ! Count dimensions
    nDims = 1
    do ii = colonPos + 1, len_trim(buf)
      if (buf(ii:ii) == ',') nDims = nDims + 1
    end do

    allocate(dataShape(nDims))
    pos = colonPos + 1
    do ii = 1, nDims
      if (ii < nDims) then
        commaPos = index(buf(pos:), ',') + pos - 1
        read(buf(pos:commaPos-1), *) dataShape(ii)
        pos = commaPos + 1
      else
        read(buf(pos:len_trim(buf)), *) dataShape(ii)
      end if
    end do
    return

    999 continue
    rank = 0

  end subroutine parseShapeAttrib


end module dftbp_io_taggedoutput
