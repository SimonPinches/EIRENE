!pb  21.11.06: index error corrected in defintion of ap%dte
 
      subroutine EIRENE_read_adas (ir,reac,isw,iz1)

cdr  purpose:  read a 2d table TAB2D of A&M data, and put them into REACDAT data structure
cdr            internal eirene reaction no. IR
cdr
cdr  input:
c           ir:           internal reaction number on eirene structure REACDAT
c           reac:
c           isw:   =0     data for interaction potential                 (not in use)
c                  =1     data for collision cross section               (not in use)
c                  =2-4   data for reaction rate coefficient             (only = 4  in use)
c                  =5-7   data for momentum weighted rate coefficient    (not in use)
c                  =8-10  data for energy weighted rate coefficient      (only = 10 in use)
c                  =11,12 other data, such as red. pop. coefficients     (not in use)
c           iz1:   particular charge state to be found within data file, 
c                  which containes charge states in the range  iza,....ize

c 
c  to be done: units, log-lin, scaling, asymptotics
 
      use EIRMOD_precision
      use EIRMOD_parmmod
      use EIRMOD_comxs
      use EIRMOD_comprt, only: iunout
 
      implicit none
 
      integer, intent(in) :: ir, isw, iz1
      character(len=*), intent(in) :: reac
      integer :: nz, nde, nte, iza, ize, io, lc, ind, ian, ien,
     .           ide, ite, iz
      character(132) :: zeile
      type(adas_data), pointer :: ap

cdr  error exit for unfinished options
      if (isw.ne.4 .and. isw.ne.10 .and. isw.ne.12)  goto 1000



c.............................................................
c
cdr   from here on a particular data file format is assumed.
c
c     this file format is described in ...
 
      read (29+ifoff,*,iostat=io) nz, nde, nte, iza, ize
 
      if (io .ne. 0) then
        write (iunout,*) ' ERROR READING FILE FROM TAB2D DATABASE '
        write (iunout,*) ' DIRECTORY IS ',reac
        call EIRENE_exit_own(1)
      end if
 
      if ((iz1 < iza) .or. (iz1 > ize)) then
        write (iunout,*) ' ERROR READING FILE FROM TAB2D DATABASE '
        write (iunout,*) ' REQUESTED Z1 IS NOT AVAILABLE '
        write (iunout,*) ' Z1, ZA, ZE ',IZ1, IZA, IZE
        call EIRENE_exit_own(1)
      end if

c  storage for 2d table, a rate coefficient vs. Te, ne. 
      allocate (ap)
      allocate (ap%dens(nde))
      allocate (ap%temp(nte))
      allocate (ap%dde(nde))
      allocate (ap%dte(nte))
      allocate (ap%tab2d(nte,nde))
 
      ap%ndens = nde
      ap%ntemp = nte
 
      read (29+ifoff,*)
 
      lc = len_trim(reac)
      if (reac(lc:lc) == 'r') then
        read (29+ifoff,*)
        read (29+ifoff,*)
      end if
 
! read densities
      read (29+ifoff,*) (ap%dens(ide),ide=1,nde)
 
! read temperatures
      read (29+ifoff,*) (ap%temp(ite),ite=1,nte)
 
! find appropriate Z1-block
 
      do
 
! read line between data blocks
        read (29+ifoff,'(A132)') zeile
        if (zeile(2:5) == '----') then
          ind = index(zeile,'Z1')
          if (ind == 0) cycle
          ian = ind + scan(zeile(ind+1:),'=') + 1
          ien = ian + scan(zeile(ian+1:),'/') - 1
          read (zeile(ian:ien),*) iz
          if (iz == iz1) exit
        end if
        if (zeile(2:5) == '____') then
          ind = index(zeile,'Q=')
          if (ind == 0) cycle
          ian = ind + scan(zeile(ind+1:),' ')
          ien = len_trim(zeile)
          read (zeile(ian:ien),*) iz
          if (iz+1 == iz1) exit
        end if

      end do
 
      do ite = 1, nte
        read (29+ifoff,*) (ap%tab2d(ite,ide), ide = 1,nde)
      end do
 
      close (29+ifoff)
 
! set up differenz arrays for first and 2nd independent parameter
 
      do ide=1,nde-1
        ap%dde(ide) = 1._dp / (ap%dens(ide+1) - ap%dens(ide))
      end do
 
      do ite=1,nte-1
        ap%dte(ite) = 1._dp / (ap%temp(ite+1) - ap%temp(ite))
      end do
 
 
      select case (isw)
 
      case (0)
        IF (REACDAT(IR)%LPOT) THEN
          WRITE (IUNOUT,*) ' POTENTIAL ALREADY SPECIFIED FOR REACTION',
     .                       IR
          DEALLOCATE (AP)
          WRITE (IUNOUT,*) ' CHECK SPECIFICATION OF REACTIONS'
          CALL EIRENE_EXIT_OWN(1)
        END IF
        reacdat(ir)%lpot = .true.

        call eirene_alloc_fit_form (reacdat(ir)%pot)
!  already done in eirene_alloc_fit_form
!        allocate (reacdat(ir)%pot)
!        nullify (reacdat(ir)%pot%line)
!        nullify (reacdat(ir)%pot%poly)
!        nullify (reacdat(ir)%pot%hyd)

        reacdat(ir)%pot%adas => ap

        reacdat(ir)%pot%ifit = 3
!  now initialized in eirene_alloc_fit_form
!        REACDAT(IR)%POT%RC1MIN = 0._DP
!        REACDAT(IR)%POT%RC1MAX = HUGE(1._DP)
!        REACDAT(IR)%POT%RC2MIN = 0._DP
!        REACDAT(IR)%POT%RC2MAX = HUGE(1._DP)
!        REACDAT(IR)%POT%FP1L = 0._DP
!        REACDAT(IR)%POT%FP1R = 0._DP
!        REACDAT(IR)%POT%FP2B = 0._DP
!        REACDAT(IR)%POT%FP2T = 0._DP
!        REACDAT(IR)%POT%JFEX1MN = 0
!        REACDAT(IR)%POT%JFEX1MX = 0
!        REACDAT(IR)%POT%JFEX2MN = 0
!        REACDAT(IR)%POT%JFEX2MX = 0
 
      case (1)
        IF (REACDAT(IR)%LCRS) THEN
          WRITE (IUNOUT,*) ' CROSS SECTION ALREADY SPECIFIED',
     .                     ' FOR REACTION', IR
          DEALLOCATE (AP)
          WRITE (IUNOUT,*) ' CHECK SPECIFICATION OF REACTIONS'
          CALL EIRENE_EXIT_OWN(1)
        END IF
        reacdat(ir)%lcrs = .true.

        call eirene_alloc_fit_form (reacdat(ir)%crs)
!  already done in eirene_alloc_fit_form
!        allocate (reacdat(ir)%crs)
!        nullify (reacdat(ir)%crs%line)
!        nullify (reacdat(ir)%crs%poly)
!        nullify (reacdat(ir)%crs%hyd)

        reacdat(ir)%crs%adas => ap

        reacdat(ir)%crs%ifit = 3
!  now initialized in eirene_alloc_fit_form
!        REACDAT(IR)%CRS%RC1MIN = 0._DP
!        REACDAT(IR)%CRS%RC1MAX = HUGE(1._DP)
!        REACDAT(IR)%CRS%RC2MIN = 0._DP
!        REACDAT(IR)%CRS%RC2MAX = HUGE(1._DP)
!        REACDAT(IR)%CRS%FP1L = 0._DP
!        REACDAT(IR)%CRS%FP1R = 0._DP
!        REACDAT(IR)%CRS%FP2B = 0._DP
!        REACDAT(IR)%CRS%FP2T = 0._DP
!        REACDAT(IR)%CRS%JFEX1MN = 0
!        REACDAT(IR)%CRS%JFEX1MX = 0
!        REACDAT(IR)%CRS%JFEX2MN = 0
!        REACDAT(IR)%CRS%JFEX2MX = 0
 
      case (2:4)
        IF (REACDAT(IR)%LRTC) THEN
          WRITE (IUNOUT,*) ' RATE COEFFICIENT ALREADY SPECIFIED',
     .                     ' FOR REACTION', IR
          DEALLOCATE (AP)
          WRITE (IUNOUT,*) ' CHECK SPECIFICATION OF REACTIONS'
          CALL EIRENE_EXIT_OWN(1)
        END IF
        reacdat(ir)%lrtc = .true.

        call eirene_alloc_fit_form (reacdat(ir)%rtc)
!  already done in eirene_alloc_fit_form
!        allocate (reacdat(ir)%rtc)
!        nullify (reacdat(ir)%rtc%line)
!        nullify (reacdat(ir)%rtc%poly)
!        nullify (reacdat(ir)%rtc%hyd)

        reacdat(ir)%rtc%adas => ap

        reacdat(ir)%rtc%ifit = 3
        REACDAT(IR)%RTC%RC1MIN = ap%temp(1)
        REACDAT(IR)%RTC%RC1MAX = ap%temp(nte)
        REACDAT(IR)%RTC%RC2MIN = ap%dens(1)
        REACDAT(IR)%RTC%RC2MAX = ap%dens(nde)
!  now initialized in eirene_alloc_fit_form
!        REACDAT(IR)%RTC%FP1L = 0._DP
!        REACDAT(IR)%RTC%FP1R = 0._DP
!        REACDAT(IR)%RTC%FP2B = 0._DP
!        REACDAT(IR)%RTC%FP2T = 0._DP
!        REACDAT(IR)%RTC%JFEX1MN = 0
!        REACDAT(IR)%RTC%JFEX1MX = 0
!        REACDAT(IR)%RTC%JFEX2MN = 0
!        REACDAT(IR)%RTC%JFEX2MX = 0
 
      case (5:7)
        IF (REACDAT(IR)%LRTCMW) THEN
          WRITE (IUNOUT,*) ' MOMEMTUM WEIGHTED RATE COEFFICIENT',
     .                     ' ALREADY SPECIFIED FOR REACTION', IR
          DEALLOCATE (AP)
          WRITE (IUNOUT,*) ' CHECK SPECIFICATION OF REACTIONS'
          CALL EIRENE_EXIT_OWN(1)
        END IF
        reacdat(ir)%lrtcmw = .true.

        call eirene_alloc_fit_form (reacdat(ir)%rtcmw)
!  already done in eirene_alloc_fit_form
!        allocate (reacdat(ir)%rtcmw)
!        nullify (reacdat(ir)%rtcmw%line)
!        nullify (reacdat(ir)%rtcmw%poly)
!        nullify (reacdat(ir)%rtcmw%hyd)

        reacdat(ir)%rtcmw%adas => ap

        reacdat(ir)%rtcmw%ifit = 3
        REACDAT(IR)%RTCMW%RC1MIN = ap%temp(1)
        REACDAT(IR)%RTCMW%RC1MAX = ap%temp(nte)
        REACDAT(IR)%RTCMW%RC2MIN = ap%dens(1)
        REACDAT(IR)%RTCMW%RC2MAX = ap%dens(nde)
!  now initialized in eirene_alloc_fit_form
!        REACDAT(IR)%RTCMW%FP1L = 0._DP
!        REACDAT(IR)%RTCMW%FP1R = 0._DP
!        REACDAT(IR)%RTCMW%FP2B = 0._DP
!        REACDAT(IR)%RTCMW%FP2T = 0._DP
!        REACDAT(IR)%RTCMW%JFEX1MN = 0
!        REACDAT(IR)%RTCMW%JFEX1MX = 0
!        REACDAT(IR)%RTCMW%JFEX2MN = 0
!        REACDAT(IR)%RTCMW%JFEX2MX = 0
 
      case (8:10)
        IF (REACDAT(IR)%LRTCEW) THEN
          WRITE (IUNOUT,*) ' ENERGY WEIGHTED RATE COEFFICIENT',
     .                     ' ALREADY SPECIFIED FOR REACTION', IR
          DEALLOCATE (AP)
          WRITE (IUNOUT,*) ' CHECK SPECIFICATION OF REACTIONS'
          CALL EIRENE_EXIT_OWN(1)
        END IF
        reacdat(ir)%lrtcew = .true.

        call eirene_alloc_fit_form (reacdat(ir)%rtcew)
!  already done in eirene_alloc_fit_form
!        allocate (reacdat(ir)%rtcew)
!        nullify (reacdat(ir)%rtcew%line)
!        nullify (reacdat(ir)%rtcew%poly)
!        nullify (reacdat(ir)%rtcew%hyd)

        reacdat(ir)%rtcew%adas => ap

        reacdat(ir)%rtcew%ifit = 3
        REACDAT(IR)%RTCEW%RC1MIN = ap%temp(1)
        REACDAT(IR)%RTCEW%RC1MAX = ap%temp(nte)
        REACDAT(IR)%RTCEW%RC2MIN = ap%dens(1)
        REACDAT(IR)%RTCEW%RC2MAX = ap%dens(nde)
!  now initialized in eirene_alloc_fit_form
!        REACDAT(IR)%RTCEW%FP1L = 0._DP
!        REACDAT(IR)%RTCEW%FP1R = 0._DP
!        REACDAT(IR)%RTCEW%FP2B = 0._DP
!        REACDAT(IR)%RTCEW%FP2T = 0._DP
!        REACDAT(IR)%RTCEW%JFEX1MN = 0
!        REACDAT(IR)%RTCEW%JFEX1MX = 0
!        REACDAT(IR)%RTCEW%JFEX2MN = 0
!        REACDAT(IR)%RTCEW%JFEX2MX = 0
 
      case (11:12)
        IF (REACDAT(IR)%LOTH) THEN
          WRITE (IUNOUT,*) ' OTHER 2D A&M DATA, E.G. RED. POP. COEF.',
     .                     ' ALREADY SPECIFIED FOR REACTION', IR
          DEALLOCATE (AP)
          WRITE (IUNOUT,*) ' CHECK SPECIFICATION OF REACTIONS'
          CALL EIRENE_EXIT_OWN(1)
        END IF
        reacdat(ir)%loth = .true.

        call eirene_alloc_fit_form (reacdat(ir)%oth)
!  already done in eirene_alloc_fit_form
!        allocate (reacdat(ir)%oth)
!        nullify (reacdat(ir)%oth%line)
!        nullify (reacdat(ir)%oth%poly)
!        nullify (reacdat(ir)%oth%hyd)

        reacdat(ir)%oth%adas => ap

        reacdat(ir)%oth%ifit = 3
!  now initialized in eirene_alloc_fit_form
!        REACDAT(IR)%OTH%RC1MIN = 0._DP
!        REACDAT(IR)%OTH%RC1MAX = HUGE(1._DP)
!        REACDAT(IR)%OTH%RC2MIN = 0._DP
!        REACDAT(IR)%OTH%RC2MAX = HUGE(1._DP)
!        REACDAT(IR)%OTH%FP1L = 0._DP
!        REACDAT(IR)%OTH%FP1R = 0._DP
!        REACDAT(IR)%OTH%FP2B = 0._DP
!        REACDAT(IR)%OTH%FP2T = 0._DP
!        REACDAT(IR)%OTH%JFEX1MN = 0
!        REACDAT(IR)%OTH%JFEX1MX = 0
!        REACDAT(IR)%OTH%JFEX2MN = 0
!        REACDAT(IR)%OTH%JFEX2MX = 0
 
      case default

        goto 1000
        
      end select
      return

1000  continue
      WRITE (IUNOUT,*) ' ERROR IN "READ_TAB2D" : '
      WRITE (IUNOUT,*) ' WRONG REACTION TYPE SPECIFIED FOR TAB2D OPTION'
      WRITE (IUNOUT,*) ' REACTION NO. ', IR
      WRITE (IUNOUT,*) ' REACTION TYPE H.', ISW
      CALL EIRENE_EXIT_OWN(1)
 
      return
      end subroutine EIRENE_read_adas
