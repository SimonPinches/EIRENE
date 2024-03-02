      module eirmod_openfile

      private

      public :: eirene_openfile

      contains

      subroutine eirene_openfile (iunit,file,form,access,status,
     .                            position,iostat)
      implicit none
      integer, intent(inout) :: iunit
      integer, intent(out), optional :: iostat
      character(len=*) :: file
      character(len=*), optional :: form, access, status, position
      character(len=:), allocatable :: frm, acc, stat, pos
      integer :: iun, io

      if (present(form)) then
        frm = form
      else
        frm = 'UNFORMATTED'
      end if
      if (present(access)) then
        acc = access
      else
        acc = 'SEQUENTIAL'
      end if
      if (present(status)) then
        stat = status
      else
        stat = 'UNKNOWN'
      end if
      if (present(position)) then
        pos = position
      else
        pos = 'ASIS'
      end if

      io = 0
      if (iunit > 0) then
        iun = iunit
#ifndef NAGFOR
        open (iun,file=file,form=frm,access=acc,status=stat,
     .        position=pos,iostat=io)
#else
        open (iun,file=file,form=frm,access=acc,status=stat,
     .        iostat=io)
#endif

      else

#ifdef F2003
#ifndef NAGFOR
        open (newunit=iun,file=file,
     .        form=frm,access=acc,status=stat,position=pos,iostat=io)
#else
        open (newunit=iun,file=file,
     .        form=frm,access=acc,status=stat,iostat=io)
#endif
        iunit=iun
#else
        iun=eirene_newunit()
#ifndef NAGFOR
        open (iun,file=file,form=frm,access=acc,status=stat,
     .        position=pos,iostat=io)
#else
        open (iun,file=file,form=frm,access=acc,status=stat,iostat=io)
#endif
        iunit=iun
#endif

      end if

      if (present(iostat)) iostat=io

      if (io /= 0) then
        call eirene_masage ('ERROR IN OPENFILE FOR FILE')
        call eirene_masage (file)
!pb        call eirene_exit_own(1)
      end if

      return

      contains

      integer function eirene_newunit ()
      implicit none
      logical used

      eirene_newunit = 99
      used = .true.
      do while (eirene_newunit.gt.0 .and. used)
        inquire(unit=eirene_newunit,opened=used)
        if (used) eirene_newunit = eirene_newunit - 1
      end do
      if (used) then
        call eirene_masage
     .   ('No available free unit numbers for output!     ')
        call eirene_exit_own(1)
      end if
      return
      end function eirene_newunit

      end subroutine eirene_openfile

      end module eirmod_openfile
