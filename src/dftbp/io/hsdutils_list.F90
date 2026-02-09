!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

!> Linked-list getChildValue overloads (variable-length data readers).
!>
!> These were originally part of dftbp_io_hsdutils but are separated into their
!> own module so the core wrapper layer stays free of linked-list dependencies.
!>
!> Callers that need the linked-list overloads import getChildValue from this
!> module *in addition to* dftbp_io_hsdutils — Fortran merges the two generics
!> automatically.
module dftbp_io_hsdutils_list
  use dftbp_common_accuracy, only : dp
  use dftbp_io_charmanip, only : unquote
  use dftbp_io_hsdutils, only : dftbp_error, getChild, getFirstTextChild, getModifier
  use dftbp_io_tokenreader, only : getNextToken, TOKEN_EOS, TOKEN_ERROR, TOKEN_OK
  use dftbp_type_linkedlist, only : TListReal, TListRealR1, TListInt, TListIntR1, &
      & TListComplex, TListComplexR1, TListString, append, len
  use hsd_data, only : hsd_table
  implicit none

  private

  interface getChildValue
    module procedure :: getChVal_lString
    module procedure :: getChVal_lReal
    module procedure :: getChVal_lRealR1
    module procedure :: getChVal_lInt
    module procedure :: getChVal_lIntR1
    module procedure :: getChVal_lComplex
    module procedure :: getChVal_lComplexR1
    module procedure :: getChVal_lIntR1RealR1
    module procedure :: getChVal_lStringIntR1RealR1
  end interface getChildValue

  public :: getChildValue

contains


  !> Resolve the child table and its text content by name.
  !!
  !! If name is empty, operates on node directly.
  !! If name is non-empty, uses getChild (which also marks the child processed).
  subroutine resolveChild_(node, name, container, text)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    type(hsd_table), pointer, intent(out) :: container
    character(len=:), allocatable, intent(out) :: text

    if (len_trim(name) == 0) then
      container => node
      call getFirstTextChild(node, text)
    else
      call getChild(node, name, container)
      call getFirstTextChild(container, text)
    end if

  end subroutine resolveChild_


  !> Get a list of strings from a child's text content.
  subroutine getChVal_lString(node, name, variableValue, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    type(TListString), intent(inout) :: variableValue
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    character(len=:), allocatable :: text, token
    type(hsd_table), pointer :: container
    integer :: iStart, iErr

    call resolveChild_(node, name, container, text)
    iStart = 1
    call getNextToken(text, token, iStart, iErr)
    do while (iErr == TOKEN_OK)
      call append(variableValue, trim(unquote(token)))
      call getNextToken(text, token, iStart, iErr)
    end do
    if (iErr == TOKEN_ERROR) then
      call dftbp_error(container, "Invalid string value in '" // name // "'")
    end if
    if (present(modifier)) call getModifier(node, name, modifier)
    if (present(child)) child => container

  end subroutine getChVal_lString


  !> Get a list of real values from a child's text content.
  subroutine getChVal_lReal(node, name, variableValue, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    type(TListReal), intent(inout) :: variableValue
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    character(len=:), allocatable :: text
    type(hsd_table), pointer :: container
    real(dp) :: buffer
    integer :: iStart, iErr

    call resolveChild_(node, name, container, text)
    iStart = 1
    call getNextToken(text, buffer, iStart, iErr)
    do while (iErr == TOKEN_OK)
      call append(variableValue, buffer)
      call getNextToken(text, buffer, iStart, iErr)
    end do
    if (iErr == TOKEN_ERROR) then
      call dftbp_error(container, "Invalid real value in '" // name // "'")
    end if
    if (present(modifier)) call getModifier(node, name, modifier)
    if (present(child)) child => container

  end subroutine getChVal_lReal


  !> Get a list of real arrays (rank-1 of given dimension) from a child's text content.
  subroutine getChVal_lRealR1(node, name, dim, variableValue, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    integer, intent(in) :: dim
    type(TListRealR1), intent(inout) :: variableValue
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    character(len=:), allocatable :: text
    type(hsd_table), pointer :: container
    real(dp), allocatable :: buffer(:)
    integer :: iStart, iErr, nItem

    allocate(buffer(dim))
    call resolveChild_(node, name, container, text)
    iStart = 1
    call getNextToken(text, buffer, iStart, iErr, nItem)
    do while (iErr == TOKEN_OK)
      call append(variableValue, buffer)
      call getNextToken(text, buffer, iStart, iErr, nItem)
    end do
    if (iErr == TOKEN_ERROR) then
      call dftbp_error(container, "Invalid real value in '" // name // "'")
    else if (iErr == TOKEN_EOS .and. nItem /= 0) then
      call dftbp_error(container, "Unexpected end of data in '" // name // "'")
    end if
    if (present(modifier)) call getModifier(node, name, modifier)
    if (present(child)) child => container

  end subroutine getChVal_lRealR1


  !> Get a list of integer values from a child's text content.
  subroutine getChVal_lInt(node, name, variableValue, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    type(TListInt), intent(inout) :: variableValue
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    character(len=:), allocatable :: text
    type(hsd_table), pointer :: container
    integer :: buffer
    integer :: iStart, iErr

    call resolveChild_(node, name, container, text)
    iStart = 1
    call getNextToken(text, buffer, iStart, iErr)
    do while (iErr == TOKEN_OK)
      call append(variableValue, buffer)
      call getNextToken(text, buffer, iStart, iErr)
    end do
    if (iErr == TOKEN_ERROR) then
      call dftbp_error(container, "Invalid integer value in '" // name // "'")
    end if
    if (present(modifier)) call getModifier(node, name, modifier)
    if (present(child)) child => container

  end subroutine getChVal_lInt


  !> Get a list of integer arrays (rank-1 of given dimension) from a child's text content.
  subroutine getChVal_lIntR1(node, name, dim, variableValue, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    integer, intent(in) :: dim
    type(TListIntR1), intent(inout) :: variableValue
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    character(len=:), allocatable :: text
    type(hsd_table), pointer :: container
    integer, allocatable :: buffer(:)
    integer :: iStart, iErr, nItem

    allocate(buffer(dim))
    call resolveChild_(node, name, container, text)
    iStart = 1
    call getNextToken(text, buffer, iStart, iErr, nItem)
    do while (iErr == TOKEN_OK)
      call append(variableValue, buffer)
      call getNextToken(text, buffer, iStart, iErr, nItem)
    end do
    if (iErr == TOKEN_ERROR) then
      call dftbp_error(container, "Invalid integer value in '" // name // "'")
    else if (iErr == TOKEN_EOS .and. nItem /= 0) then
      call dftbp_error(container, "Unexpected end of data in '" // name // "'")
    end if
    if (present(modifier)) call getModifier(node, name, modifier)
    if (present(child)) child => container

  end subroutine getChVal_lIntR1


  !> Get a list of complex values from a child's text content.
  subroutine getChVal_lComplex(node, name, variableValue, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    type(TListComplex), intent(inout) :: variableValue
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    character(len=:), allocatable :: text
    type(hsd_table), pointer :: container
    complex(dp) :: buffer
    integer :: iStart, iErr

    call resolveChild_(node, name, container, text)
    iStart = 1
    call getNextToken(text, buffer, iStart, iErr)
    do while (iErr == TOKEN_OK)
      call append(variableValue, buffer)
      call getNextToken(text, buffer, iStart, iErr)
    end do
    if (iErr == TOKEN_ERROR) then
      call dftbp_error(container, "Invalid complex value in '" // name // "'")
    end if
    if (present(modifier)) call getModifier(node, name, modifier)
    if (present(child)) child => container

  end subroutine getChVal_lComplex


  !> Get a list of complex arrays (rank-1 of given dimension) from a child's text content.
  subroutine getChVal_lComplexR1(node, name, dim, variableValue, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    integer, intent(in) :: dim
    type(TListComplexR1), intent(inout) :: variableValue
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    character(len=:), allocatable :: text
    type(hsd_table), pointer :: container
    complex(dp), allocatable :: buffer(:)
    integer :: iStart, iErr, nItem

    allocate(buffer(dim))
    call resolveChild_(node, name, container, text)
    iStart = 1
    call getNextToken(text, buffer, iStart, iErr, nItem)
    do while (iErr == TOKEN_OK)
      call append(variableValue, buffer)
      call getNextToken(text, buffer, iStart, iErr, nItem)
    end do
    if (iErr == TOKEN_ERROR) then
      call dftbp_error(container, "Invalid complex value in '" // name // "'")
    else if (iErr == TOKEN_EOS .and. nItem /= 0) then
      call dftbp_error(container, "Unexpected end of data in '" // name // "'")
    end if
    if (present(modifier)) call getModifier(node, name, modifier)
    if (present(child)) child => container

  end subroutine getChVal_lComplexR1


  !> Get paired lists of integer and real arrays from a child's text content.
  subroutine getChVal_lIntR1RealR1(node, name, dimInt, valueInt, dimReal, valueReal, &
      & modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    integer, intent(in) :: dimInt
    type(TListIntR1), intent(inout) :: valueInt
    integer, intent(in) :: dimReal
    type(TListRealR1), intent(inout) :: valueReal
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    character(len=:), allocatable :: text
    type(hsd_table), pointer :: container
    integer, allocatable :: bufferInt(:)
    real(dp), allocatable :: bufferReal(:)
    integer :: iStart, iErr, nItem

    allocate(bufferInt(dimInt))
    allocate(bufferReal(dimReal))
    call resolveChild_(node, name, container, text)
    iStart = 1
    iErr = TOKEN_OK
    do while (iErr == TOKEN_OK)
      call getNextToken(text, bufferInt, iStart, iErr, nItem)
      if (iErr == TOKEN_ERROR) then
        call dftbp_error(container, "Invalid integer value in '" // name // "'")
      end if
      if (iErr == TOKEN_EOS .and. nItem /= 0) then
        call dftbp_error(container, "Unexpected end of data in '" // name // "'")
      end if
      if (iErr == TOKEN_OK) then
        call append(valueInt, bufferInt)
        call getNextToken(text, bufferReal, iStart, iErr, nItem)
        if (iErr == TOKEN_ERROR) then
          call dftbp_error(container, "Invalid real value in '" // name // "'")
        end if
        if (iErr == TOKEN_EOS .and. nItem /= 0) then
          call dftbp_error(container, "Unexpected end of data in '" // name // "'")
        end if
        if (iErr == TOKEN_OK) then
          call append(valueReal, bufferReal)
        end if
      end if
    end do
    if (len(valueInt) /= len(valueReal)) then
      call dftbp_error(container, "Unexpected end of data in '" // name // "'")
    end if
    if (present(modifier)) call getModifier(node, name, modifier)
    if (present(child)) child => container

  end subroutine getChVal_lIntR1RealR1


  !> Get combined string, integer-array and real-array lists from a child's text content.
  subroutine getChVal_lStringIntR1RealR1(node, name, valueStr, dimInt, valueInt, dimReal, &
      & valueReal, modifier, child)
    type(hsd_table), intent(inout), target :: node
    character(len=*), intent(in) :: name
    type(TListString), intent(inout) :: valueStr
    integer, intent(in) :: dimInt
    type(TListIntR1), intent(inout) :: valueInt
    integer, intent(in) :: dimReal
    type(TListRealR1), intent(inout) :: valueReal
    character(len=:), allocatable, intent(out), optional :: modifier
    type(hsd_table), pointer, intent(out), optional :: child

    character(len=:), allocatable :: text, bufferStr
    type(hsd_table), pointer :: container
    integer, allocatable :: bufferInt(:)
    real(dp), allocatable :: bufferReal(:)
    integer :: iStart, iErr, nItem

    allocate(bufferInt(dimInt))
    allocate(bufferReal(dimReal))
    call resolveChild_(node, name, container, text)
    iStart = 1
    iErr = TOKEN_OK
    do while (iErr == TOKEN_OK)
      call getNextToken(text, bufferStr, iStart, iErr)
      if (iErr == TOKEN_ERROR) then
        call dftbp_error(container, "Invalid string value in '" // name // "'")
      end if
      if (iErr == TOKEN_EOS) exit
      call append(valueStr, bufferStr)

      call getNextToken(text, bufferInt, iStart, iErr, nItem)
      if (iErr /= TOKEN_OK) then
        call dftbp_error(container, "Invalid integer value in '" // name // "'")
      end if
      call append(valueInt, bufferInt)

      call getNextToken(text, bufferReal, iStart, iErr, nItem)
      if (iErr /= TOKEN_OK) then
        call dftbp_error(container, "Invalid real value in '" // name // "'")
      end if
      call append(valueReal, bufferReal)
    end do
    if (len(valueStr) /= len(valueInt) .or. len(valueInt) /= len(valueReal)) then
      call dftbp_error(container, "Unexpected end of data in '" // name // "'")
    end if
    if (present(modifier)) call getModifier(node, name, modifier)
    if (present(child)) child => container

  end subroutine getChVal_lStringIntR1RealR1


end module dftbp_io_hsdutils_list
