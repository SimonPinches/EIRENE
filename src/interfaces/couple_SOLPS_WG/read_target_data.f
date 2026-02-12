
      subroutine eirene_read_target_data (iun, iunout, ntrg, nflai)

      use eirmod_precision
      use eirmod_parmmod
      use eirmod_braeir
      use eirmod_comusr, only : nmassp, ncharp, nchrgp

      implicit none

      integer, intent(in) :: iun, iunout, ntrg, nflai
!pb      integer, intent(out) :: nfc
      integer :: i, j, ipl, ntr_ion, inmass, inchar, inchrg,
     .           ntrgi, ier, ifl, ndat, itrg
      real(dp) :: pmass, pchar, pchrg
      CHARACTER(256) :: sstr, line
      external :: eirene_locstr, eirene_exit_own

!pb      nfc = 0

      OPEN (UNIT=IUN,ACCESS='SEQUENTIAL',FORM='FORMATTED',ERR=1986)
      REWIND IUN

c number of targets:
      write(sstr,'(a)')
     .          '*** TARGET DATA'
      CALL EIRENE_locstr(iun,sstr,ier)
      if(ier /=0) then
        write(iunout,*) 'READ_TARGET_DATA: ',trim(sstr),' not found'
        close(iun)
        call EIRENE_exit_own(1)
      endif
      call eirene_skip_comments(iun,line)

      read(line,*) ntrgi
      if (ntrgi /= ntrg) then
        write (iunout,*) 'READ_TARGET_DATA:', trim(sstr)
        write (iunout,*) ' wrong number of targets in plasma file'
!       write (iunout,*) ' check for correct number in file ',filename
        call EIRENE_exit_own(1)
      endif

      if (.not.allocated(trgt)) allocate(trgt(ntrg))

c read target data
      do itrg = 1, ntrg
c misc target data for target ITRG
        write (sstr,'(a,i0,a)') '*** MISC TARGET DATA #',itrg
        CALL EIRENE_locstr(iun,sstr,ier)
        if(ier /=0) then
          write(iunout,*) 'READ_TARGET_DATA: tag ',trim(sstr),
     .                    ' not found'
          close(iun)
          call EIRENE_exit_own(1)
        endif
        call eirene_skip_comments(iun,line)

        read(line,*) ndat
        call eirene_alloc_target_data (itrg, ndat, nflai)
!pb        nfc = nfc + ndat
        trgt(itrg)%ntrgdat = ndat

        do i = 1, ndat
          read (iun,*) j, trgt(itrg)%faces(i), trgt(itrg)%fcori(i),
     .                 trgt(itrg)%tet(i),
     .                 trgt(itrg)%fe(i), trgt(itrg)%fsh(i),
     .                 trgt(itrg)%flength(i)
        enddo

c read data for ion IFL and target ITRG
        do ifl = 1, nflai
          write (sstr,'(a,i0,a,i0)') '*** ION #',ifl,
     .           ' TARGET DATA #',itrg
          CALL EIRENE_locstr(iun,sstr,ier)
          if(ier /=0) then
            write(iunout,*) 'READ_TARGET_DATA: tag ',trim(sstr),
     .                      ' not found'
            close(iun)
            call EIRENE_exit_own(1)
          endif
          call eirene_skip_comments(iun,line)

          read(line,*) pmass,pchar,pchrg
          inmass = nint(pmass)
          inchar = nint(pchar)
          inchrg = nint(pchrg)

          ipl = 0
          do i=1,nflai

            if (nmassp(i) /= inmass) cycle
            if (ncharp(i) /= inchar) cycle
            if (nchrgp(i) /= inchrg) cycle

            ipl = i
            if (ipl.ne.0) exit
          end do

          read(iun,*) ntr_ion

          if (ntr_ion /= ndat) then
            write (iunout,*) 'READ_TARGET_DATA:', trim(sstr)
            write (iunout,*) ' wrong number of faces in target file'
!         write (iunout,*) ' check for correct number in file ',filename
            call EIRENE_exit_own(1)
          endif

          do i = 1, ndat
            read (31,*)
     .        j,trgt(itrg)%flux(i,ipl),trgt(itrg)%tit(i),
     .          trgt(itrg)%dnit(i,ipl),trgt(itrg)%vxt(i,ipl),
     .          trgt(itrg)%vyt(i,ipl),trgt(itrg)%vzt(i,ipl),
     .          trgt(itrg)%fi(i,ipl),trgt(itrg)%fel(i,ipl),
     .          trgt(itrg)%vpart(i,ipl),trgt(itrg)%mach(i,ipl),
     .          trgt(itrg)%usr(i,ipl),trgt(itrg)%zit(i,ipl)
          end do

        end do

      end do

      return
 1986 WRITE(IUNOUT,*)
     w      "ERROR IN READ_TARGET_DATA: CANNOT OPEN FORT.31",
     w      " (PLASMA BACKGROUND)"
      CALL EIRENE_EXIT_OWN (1)
      return

      contains

      subroutine eirene_skip_comments(iun,line)
      integer, intent(in) :: iun
      character(*), intent(out) :: line

      read(iun,'(a)') line
      do while (line(1:1) == '*')
        read(iun,'(a)') line
      enddo

      return
      end subroutine eirene_skip_comments


      end subroutine eirene_read_target_data
