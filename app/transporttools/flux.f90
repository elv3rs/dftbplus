!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

program flux
  implicit none

  integer, parameter :: dp = kind(1.0d0)

  integer :: nat, natoms
  integer, dimension(:,:), allocatable :: neig
  integer, dimension(:), allocatable :: nn
  real(dp), dimension(:,:), allocatable :: Inm
  real(dp), dimension(:,:), allocatable :: coord
  integer :: m, itmp, io, i, j, n, k, maxnn, maxnn_valid
  real(dp) :: Iz, Imax, frac, width, arr_len
  real(dp) :: rtmp(3)
  real(dp) :: e(3)
  real(dp) :: dr(3)
  character(64) :: arg, filename1, filename2, color
  character(4) :: id
  logical :: bondflux
  integer :: iargc
  integer :: fp

  iargc = command_argument_count()
  if (iargc < 4) then
     write(*,*) "flux supercell.xyz -a|b 'lcurrent.dat' maxneigbors [-f fraction -w width -c color]"
     write(*,*) " -a : atom currents; -b: bond currents"
     write(*,*) " maxneighbours: neighbours considered in flux calculation"
     write(*,*) " fraction [1.0]: The arrow lengths are normalized to I_max/fraction"
     write(*,*) " width [0.2]: arrows are given a width depending on length"
     write(*,*) "              arrows too long are made thicker and rescaled"
     write(*,*) " color      : colors for the arrows"
     stop
  end if

  call get_command_argument(1, filename1)

  call get_command_argument(2, arg)
  if (trim(arg) == "-b") bondflux = .true.
  if (trim(arg) == "-a") bondflux = .false.

  call get_command_argument(3, filename2)

  n = 40
  call get_command_argument(4, arg)
  read(arg,*) maxnn

  frac = 1.0
  width = 0.2
  color = "yellow"

  do m = 5, iargc, 2
    call get_command_argument(m, arg)
    if (trim(arg) == "-f") then
      call get_command_argument(m + 1, arg)
      read(arg,*) frac
    end if
    if (trim(arg) == "-w") then
      call get_command_argument(m + 1, arg)
      read(arg,*) width
    end if
    if (trim(arg) == "-c") then
      call get_command_argument(m + 1, arg)
      read(arg,*) color
    end if
  end do

  open(fp, file=trim(filename1), action="read")
  read(fp, *) nat
  allocate(neig(nat,n))
  allocate(nn(nat))
  allocate(Inm(nat,n))
  allocate(coord(3,nat))
  read(fp, *)
  do m=1, nat
    read(fp,*) id, coord(1:3,m)
  end do
  close(fp)

  ! Figure out the number of atoms
  open(fp,file=trim(filename2), action="read")
  Inm = 0.0_dp
  loop1: do m=1, nat
    read(fp,*, iostat=io) itmp, rtmp(1:3), nn(m), (neig(m,i), Inm(m,i), i=1,nn(m))
    if (io<0) then
      natoms = m-1
      exit loop1
    end if
  end do loop1
  close(fp)

  Imax=maxval(abs(Inm(1:natoms,1:maxnn)))
  print*,"# Natoms=",natoms
  print*,"# Imax=",Imax

  if(Imax==0.0_dp) Imax=1.0_dp
  k=1

  ! Bond flux
  if (bondflux) then
    do m = 1, natoms

       maxnn_valid = min(maxnn,nn(m))

       do j=1, maxnn_valid

         ! vector joining atoms
         n=neig(m,j)
         ! vector directed along the bond:
         dr(:)=coord(:,n)-coord(:,m)
         !draw the bond current only if it follows the bond
         arr_len = frac*Inm(m,j)/Imax

         if (arr_len > 0.1_dp) then
            write(*,*) "# ",n,m,Inm(m,j)/Imax
            write(id,"(i4.4)") k
            arg="draw arr"//id//" arrow width"

            e(:)=coord(:,m)+dr(:)*0.9
            !write(*,'(a24,f6.2,a2,3(f10.4),a2,3(f10.4),a14)') trim(arg),arr_len*width,' {', &
            !         & coord(1:3,m),'}{',e(1:3),'} color '//trim(color)

           if(arr_len < 0.2_dp) then
              e(:)=coord(:,m)+dr(:)*0.9
              write(*,"(a24,f6.2,a2,3(f10.4),a2,3(f10.4),a14)") trim(arg),0.25*width," {", &
                    & coord(1:3,m),"}{",e(1:3),"} color "//trim(color)
           end if
           if(arr_len >= 0.2 .and. arr_len < 0.4_dp) then
              e(:)=coord(:,m)+dr(:)*0.9
              write(*,"(a24,f6.2,a2,3(f10.4),a2,3(f10.4),a14)") trim(arg),0.5*width," {", &
                    & coord(1:3,m),"}{",e(1:3),"} color "//trim(color)
           end if
           if(arr_len >= 0.4 .and. arr_len < 1.0_dp) then
              e(:)=coord(:,m)+dr(:)*0.9
              write(*,"(a24,f6.2,a2,3(f10.4),a2,3(f10.4),a14)") trim(arg),width," {", &
                    & coord(1:3,m),"}{",e(1:3),"} color "//trim(color)
           end if
           if(arr_len >= 1.0_dp .and. arr_len < 1.5_dp) then
              e(:)=coord(:,m)+dr(:)*0.8
              write(*,"(a24,f6.2,a2,3(f10.4),a2,3(f10.4),a14)") trim(arg),1.5*width," {", &
                    & coord(1:3,m),"}{",e(1:3),"} color "//trim(color)
           end if
           if(arr_len >= 1.5_dp .and. arr_len < 2.0_dp) then
              e(:)=coord(:,m)+dr(:)*0.8
               write(*,"(a24,f6.2,a2,3(f10.4),a2,3(f10.4),a14)") trim(arg),2.0*width," {", &
                     & coord(1:3,m),"}{",e(1:3),"} color "//trim(color)
           end if
           if(arr_len >= 2.0_dp .and. arr_len < 2.5_dp) then
              e(:)=coord(:,m)+dr(:)*0.8
              write(*,"(a24,f6.2,a2,3(f10.4),a2,3(f10.4),a14)") trim(arg),2.5*width," {", &
                    & coord(1:3,m),"}{",e(1:3),"} color "//trim(color)
           end if
           if(arr_len >= 2.5_dp) then
              e(:)=coord(:,m)+dr(:)*0.8
              write(*,"(a24,f6.2,a2,3(f10.4),a2,3(f10.4),a14)") trim(arg),3.0*width," {", &
                     & coord(1:3,m),"}{",e(1:3),"} color "//trim(color)
            end if
            k=k+1
         end if

       end do
    end do

  else  ! Atomic flux

    do m = 1, natoms
       Iz =0.0_dp
       e = 0.0_dp

       do j=1, maxnn

         ! vector joining atoms
         n=neig(m,j)
         ! vector directed along the bond:
         dr(:)= coord(:,n)-coord(:,m)
         e(:) = e(:) + frac*Inm(m,j)*dr(:)/Imax

       end do

       write(id,"(i4.4)") k
       arg="draw arr"//id//" arrow width"
       write(*,"(a24,f6.2,a2,3(f10.4),a2,3(f10.4),a14)") trim(arg),width," {", &
           & coord(1:3,m),"}{",coord(:,m)+e(1:3),"} color "//trim(color)
       k=k+1

    end do

  end if

  deallocate(neig,nn,Inm,coord)

end program flux
