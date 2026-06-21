!-*- mode: F90 -*-!
!------------------------------------------------------------!
! Copyright (C) 2026 Wannier Developer Group                 !
!                                                            !
! This library is free software; you can redistribute it     !
! and/or modify it under the terms of the GNU Lesser General !
! Public License as published by the Free Software           !
! Foundation; either version 2.1 of the License, or (at your !
! option) any later version.                                 !
!                                                            !
! This library is distributed in the hope that it will be    !
! useful,but WITHOUT ANY WARRANTY; without even the implied  !
! warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR    !
! PURPOSE.  See the GNU Lesser General Public License for    !
! more details.                                              !
!                                                            !
! You should have received a copy of the GNU Lesser General  !
! Public License along with this library; if not, see        !
! <https://www.gnu.org/licenses/>.                           !
!                                                            !
! The webpage of the Wannier90 code is                       !
! <https://www.wannier.org>.                                 !
!                                                            !
! The Wannier90 code is hosted on GitHub                     !
! <https://github.com/wannier-developers/wannier90>          !
!------------------------------------------------------------!
!                                                            !
!  w90orb2orb: reformat orb files                            !
!                                                            !
!------------------------------------------------------------!

module w90orb_parameters

  use w90_types

  implicit none

  public

  integer, save :: num_kpts !BGS put in k_point_type?
  integer, save :: num_bands !! Number of bands

end module w90orb_parameters

module w90_conv_orb
  !! Module to convert orb files from formatted to unformmated
  !! and vice versa - useful for switching between computers
  use w90_constants, only: dp

  implicit none

  logical, save :: export_flag
  ! header length is identical to that of pw2wannier90.f90
  character(len=60), save :: header
  complex(kind=dp), allocatable, save :: orb_o(:, :, :, :)

contains

  subroutine io_error(error_msg, stdout)
    !================================================
    !
    !! Abort the code giving an error message
    !
    !================================================

    implicit none

    character(len=*), intent(in) :: error_msg
    integer :: stdout

    close (stdout)
    write (*, '(1x,a)') trim(error_msg)
    write (*, '(A)') "Error: examine the output/error file for details"
    stop
  end subroutine io_error

  !================================================!
  subroutine print_usage(stdout)
    !================================================!
    !
    !! Writes the usage of the program to stdout
    !
    !================================================!

    implicit none

    integer, intent(in) :: stdout

    write (stdout, '(A)') "Usage:"
    write (stdout, '(A)') "  w90orb2orb.x ACTION [SEEDNAME]"
    write (stdout, '(A)') "where ACTION can be one of the following:"
    write (stdout, '(A)') "  -export"
    write (stdout, '(A)') "  -u2f"
    write (stdout, '(A)') "      Convert from unformatted (standard) format to formatted format, to export"
    write (stdout, '(A)') "      the orb file on a different machine."
    write (stdout, '(A)') "      The seedname.orb file is read and the seedname.orb.fmt file is generated."
    write (stdout, '(A)') "  -import"
    write (stdout, '(A)') "  -f2u"
    write (stdout, '(A)') "      Convert from formatted format to unformatted (standard) format, to import"
    write (stdout, '(A)') "      the orb file seedname.orb.fmt from a different machine."
    write (stdout, '(A)') "      The seedname.orb.fmt file is read and the seedname.orb file is generated."
  end subroutine print_usage

  !================================================!
  subroutine conv_get_seedname(stdout, seedname)
    !================================================!
    !
    !! Set the seedname from the command line
    !
    !================================================!
    implicit none

    integer, intent(in) :: stdout
    character(len=50), intent(inout)  :: seedname

    integer :: num_arg
    character(len=50) :: ctemp

    num_arg = command_argument_count()
    if (num_arg == 1) then
      seedname = 'wannier'
    elseif (num_arg == 2) then
      call get_command_argument(2, seedname)
    else
      call print_usage(stdout)
      call io_error('Wrong command line arguments, see logfile for usage', stdout)
    end if

    ! If on the command line the whole seedname.win was passed, I strip the last ".win"
    if (len(trim(seedname)) .ge. 5) then
      if (seedname(len(trim(seedname)) - 4 + 1:) .eq. ".win") then
        seedname = seedname(:len(trim(seedname)) - 4)
      end if
    end if

    call get_command_argument(1, ctemp)
    if (index(ctemp, '-import') > 0) then
      export_flag = .false.
    elseif (index(ctemp, '-f2u') > 0) then
      export_flag = .false.
    elseif (index(ctemp, '-export') > 0) then
      export_flag = .true.
    elseif (index(ctemp, '-u2f') > 0) then
      export_flag = .true.
    else
      write (stdout, '(A)') 'Wrong command line action: '//trim(ctemp)
      call print_usage(stdout)
      call io_error('Wrong command line arguments, see logfile for usage', stdout)
    end if

  end subroutine conv_get_seedname

  !================================================!
  subroutine conv_read_orb(stdout, seedname)
    !================================================!
    !
    !! Read unformatted orb file
    !
    !================================================!

    use w90_constants, only: eps6, dp
    use w90orb_parameters, only: num_bands, num_kpts

    implicit none

    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname

    integer :: orb_unit, m, n, ik, ierr, s, counter
    complex(kind=dp), allocatable :: orb_temp(:, :)

    write (stdout, '(3a)') 'Reading information from unformatted file ', trim(seedname), '.orb :'

    open (newunit=orb_unit, file=trim(seedname)//'.orb', status='old', form='unformatted', err=109)

    ! Read comment line
    read (orb_unit, err=110, end=110) header
    header = ADJUSTL(header)
    write (stdout, '(1x,a)') trim(header)

    ! Consistency checks
    read (orb_unit, err=110, end=110) num_bands, num_kpts
    write (stdout, '(1x,a,i0)') "Number of bands: ", num_bands
    write (stdout, '(1x,a,i0)') "Number of k-points: ", num_kpts

    allocate (orb_o(num_bands, num_bands, num_kpts, 3), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating spm_temp in conv_read_orb', stdout)

    allocate (orb_temp(3, (num_bands*(num_bands + 1))/2), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating spm_temp in conv_read_orb', stdout)
    do ik = 1, num_kpts
      read (orb_unit) ((orb_temp(s, m), s=1, 3), m=1, (num_bands*(num_bands + 1))/2)
      counter = 0
      do m = 1, num_bands
        do n = 1, m
          counter = counter + 1
          orb_o(n, m, ik, 1) = orb_temp(1, counter)
          orb_o(n, m, ik, 2) = orb_temp(2, counter)
          orb_o(n, m, ik, 3) = orb_temp(3, counter)
          ! Although each diagonal element of orb_o should be a real number,
          ! actually it has a very small imaginary part.
          ! We skip the conjugation on the diagonal elements so that
          ! the file after formatted <==> unformatted conversions is exactly
          ! the same as the original file, otherwise the diagonal elements
          ! are the conjugations of those of the original file.
          if (m == n) cycle
          orb_o(m, n, ik, 1) = conjg(orb_temp(1, counter))
          orb_o(m, n, ik, 2) = conjg(orb_temp(2, counter))
          orb_o(m, n, ik, 3) = conjg(orb_temp(3, counter))
        end do
      end do
    end do

    close (orb_unit)
    write (stdout, '(1x,a)') "orb: read."

    deallocate (orb_temp, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating spm_temp in conv_read_orb', stdout)

    write (stdout, '(1x,a)') 'read done.'

    return

109 call io_error('Error opening '//trim(seedname)//'.orb.fmt in conv_read_orb', stdout)
110 call io_error('Error reading '//trim(seedname)//'.orb.fmt in conv_read_orb', stdout)

  end subroutine conv_read_orb

  !================================================!
  subroutine conv_read_orb_fmt(stdout, seedname)
    !================================================!
    !
    !! Read formatted orb file
    !
    !================================================!

    use w90_constants, only: eps6, dp
    use w90orb_parameters, only: num_bands, num_kpts

    implicit none

    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname

    integer :: orb_unit, m, n, ik, ierr
    real(kind=dp) :: s_real, s_img

    write (stdout, '(3a)') 'Reading information from formatted file ', trim(seedname), '.orb.fmt :'

    open (newunit=orb_unit, file=trim(seedname)//'.orb.fmt', status='old', position='rewind', form='formatted', err=109)

    ! Read comment line
    read (orb_unit, '(a)') header
    header = ADJUSTL(header)
    write (stdout, '(1x,a)') trim(header)

    ! Consistency checks
    read (orb_unit, *, err=110, end=110) num_bands, num_kpts
    write (stdout, '(1x,a,i0)') "Number of bands: ", num_bands
    write (stdout, '(1x,a,i0)') "Number of k-points: ", num_kpts

    allocate (orb_o(num_bands, num_bands, num_kpts, 3), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating orb_o in conv_read_orb_fmt', stdout)

    do ik = 1, num_kpts
      do m = 1, num_bands
        do n = 1, m
          read (orb_unit, *, err=110, end=110) s_real, s_img
          orb_o(n, m, ik, 1) = cmplx(s_real, s_img, dp)
          read (orb_unit, *, err=110, end=110) s_real, s_img
          orb_o(n, m, ik, 2) = cmplx(s_real, s_img, dp)
          read (orb_unit, *, err=110, end=110) s_real, s_img
          orb_o(n, m, ik, 3) = cmplx(s_real, s_img, dp)
          ! Although each diagonal element of orb_o should be a real number,
          ! actually it has a very small imaginary part.
          ! We skip the conjugation on the diagonal elements so that
          ! the file after formatted <==> unformatted conversions is exactly
          ! the same as the original file, otherwise the diagonal elements
          ! are the conjugations of those of the original file.
          if (m == n) cycle
          ! Read upper-triangular part, now build the rest
          orb_o(m, n, ik, 1) = conjg(orb_o(n, m, ik, 1))
          orb_o(m, n, ik, 2) = conjg(orb_o(n, m, ik, 2))
          orb_o(m, n, ik, 3) = conjg(orb_o(n, m, ik, 3))
        end do
      end do
    end do

    close (orb_unit)

    write (stdout, '(1x,a)') 'read done.'

    return

109 call io_error('Error opening '//trim(seedname)//'.orb.fmt in conv_read_orb_fmt', stdout)
110 call io_error('Error reading '//trim(seedname)//'.orb.fmt in conv_read_orb_fmt', stdout)

  end subroutine conv_read_orb_fmt

  !================================================!
  subroutine conv_write_orb(stdout, seedname)
    !================================================!
    !
    !! Write unformatted orb file
    !
    !================================================!

    use w90_io, only: io_date
    use w90orb_parameters, only: num_bands, num_kpts

    implicit none

    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname

    integer :: orb_unit, m, n, ik, counter, s, ierr
    complex(kind=dp), allocatable :: orb_temp(:, :)

    write (stdout, '(3a)') 'Writing information to unformatted file ', trim(seedname), '.orb :'

    open (newunit=orb_unit, file=trim(seedname)//'.orb', form='unformatted')

    allocate (orb_temp(3, (num_bands*(num_bands + 1))/2), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating spm_temp in conv_write_orb', stdout)

    write (orb_unit) header
    write (orb_unit) num_bands, num_kpts

    do ik = 1, num_kpts
      counter = 0
      do m = 1, num_bands
        do n = 1, m
          counter = counter + 1
          do s = 1, 3
            orb_temp(s, counter) = orb_o(n, m, ik, s)
          end do
        end do
      end do
      write (orb_unit) ((orb_temp(s, m), s=1, 3), m=1, ((num_bands*(num_bands + 1))/2))
    end do

    close (orb_unit)

    write (stdout, '(1x,a)') 'write done.'

  end subroutine conv_write_orb

  !================================================!
  subroutine conv_write_orb_fmt(stdout, seedname)
    !================================================!
    !
    !! Write formatted orb file
    !
    !================================================!

    use w90_io, only: io_date
    use w90orb_parameters, only: num_bands, num_kpts

    implicit none

    integer, intent(in) :: stdout
    character(len=50), intent(in)  :: seedname

    integer :: orb_unit, m, n, ik, s

    write (stdout, '(3a)') 'Writing information to formatted file ', trim(seedname), '.orb.fmt :'

    open (newunit=orb_unit, file=trim(seedname)//'.orb.fmt', form='formatted', status='replace', position='rewind')

    write (orb_unit, *) header
    write (orb_unit, *) num_bands, num_kpts

    do ik = 1, num_kpts
      do m = 1, num_bands
        do n = 1, m
          do s = 1, 3
            write (orb_unit, '(2es26.16)') orb_o(n, m, ik, s)
          end do
        end do
      end do
    end do

    close (orb_unit)

    write (stdout, '(1x, a)') 'write done.'

  end subroutine conv_write_orb_fmt

end module w90_conv_orb

program w90orb2orb
  !! Program to convert orb files from formatted to unformmated
  !! and vice versa - useful for switching between computers
  use w90_constants, only: dp
  use w90_conv_orb

  implicit none

  ! Export mode:
  !  TRUE:  create formatted .orb.fmt from unformatted .orb ('-export')
  !  FALSE: create unformatted .orb from formatted .orb.fmt ('-import')

  !logical :: file_found
  !integer :: file_unit
  integer :: stdout !, ierr, num_nodes
  character(len=50) :: seedname

  open (newunit=stdout, file='w90orb2orb.log')

  call conv_get_seedname(stdout, seedname)

  if (export_flag .eqv. .true.) then
    call conv_read_orb(stdout, seedname)
    write (stdout, '(a)') ''
    call conv_write_orb_fmt(stdout, seedname)
  else
    call conv_read_orb_fmt(stdout, seedname)
    write (stdout, '(a)') ''
    call conv_write_orb(stdout, seedname)
  end if

  close (unit=stdout)

end program w90orb2orb
