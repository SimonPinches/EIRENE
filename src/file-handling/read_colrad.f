       subroutine EIRENE_read_colrad (ir,reac,isw,iz1)

cdr  purpose:  prepare usage of A&M data from an internal, built-in, 
cdr            collisional radiative code: 
cdr  1)  H_colrad, 
cdr  2)  He_colrad,
cdr  3)  H2-colrad....
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
c           iz1:   not in use 
c                  
cdr:  currently used only H.4, 2.1.5 and H.10, 2.1.5, EI,  ionisation
cdr   to be done:         H.4. 2.1.8 and H.10, 2.1.8, RC   recombination
cdr                       and  H.11, H.12: selected population
c 
c  to be done: units, log-lin, scaling, asymptotics
 
      use EIRMOD_precision
      use EIRMOD_parmmod
      use EIRMOD_comxs
      use EIRMOD_comprt, only: iunout
 
      implicit none
 
      integer, intent(in) :: ir, isw, iz1
      character(len=*), intent(in) :: reac

      close (29+ifoff)  ! nothing further to be read, currently

cdr  error exit for unfinished options
      if (isw.ne.4 .and. isw.ne.10)  goto 1000
cdr  tbd: also exit unless 2.1.5,  in particular: 
cdr       2.1.8 (recombination) is missing.
cdr  other reactions are not programmed in xsectp, rate-coeff, energy rate coef. 
      

!  ALREADY INITIALIZED IN EIRENE_INIT_CMDTA
!        REACDAT(IR)%ETH = 0._DP
!        REACDAT(IR)%RTMAX = 0._DP
!        REACDAT(IR)%ERTMAX = -HUGE(1._DP)

        SELECT CASE (ISW)
        CASE (2:4)
          IF (REACDAT(IR)%LRTC) THEN
            WRITE (IUNOUT,*) ' RATE COEFFICIENT ALREADY SPECIFIED',
     .                       ' FOR REACTION', IR
            WRITE (IUNOUT,*) ' CHECK SPECIFICATION OF REACTIONS'
            CALL EIRENE_EXIT_OWN(1)
          END IF

          CALL EIRENE_ALLOC_FIT_FORM (REACDAT(IR)%RTC)
!  ALREADY DONE IN EIRENE_ALLOC_FIT_FORM
!          ALLOCATE (REACDAT(IR)%RTC)
!          NULLIFY(REACDAT(IR)%RTC%ADAS)
!          NULLIFY(REACDAT(IR)%RTC%LINE)
!          NULLIFY(REACDAT(IR)%RTC%POLY)
!          NULLIFY(REACDAT(IR)%RTC%HYD)

          REACDAT(IR)%LRTC = .TRUE.
          REACDAT(IR)%RTC%IFIT = 5

!  ALREADY INITIALIZED IN EIRENE_ALLOC_FIT_FORM
!          REACDAT(IR)%RTC%RC1MIN = 0._DP
!          REACDAT(IR)%RTC%RC1MAX = HUGE(1._DP)
!          REACDAT(IR)%RTC%RC2MIN = 0._DP
!          REACDAT(IR)%RTC%RC2MAX = HUGE(1._DP)
!          REACDAT(IR)%RTC%FP1L = 0._DP
!          REACDAT(IR)%RTC%FP1R = 0._DP
!          REACDAT(IR)%RTC%FP2B = 0._DP
!          REACDAT(IR)%RTC%FP2T = 0._DP
!          REACDAT(IR)%RTC%JFEX1MN = 0
!          REACDAT(IR)%RTC%JFEX1MX = 0
!          REACDAT(IR)%RTC%JFEX2MN = 0
!          REACDAT(IR)%RTC%JFEX2MX = 0
          
        CASE (5:7)
          IF (REACDAT(IR)%LRTCMW) THEN
            WRITE (IUNOUT,*) ' MOMENTUM WEIGHTED RATE COEFFICIENT',
     .                       ' ALREADY SPECIFIED FOR REACTION', IR
            WRITE (IUNOUT,*) ' CHECK SPECIFICATION OF REACTIONS'
            CALL EIRENE_EXIT_OWN(1)
          END IF

          CALL EIRENE_ALLOC_FIT_FORM (REACDAT(IR)%RTCMW)
!  ALREADY DONE IN EIRENE_ALLOC_FIT_FORM
!          ALLOCATE (REACDAT(IR)%RTCMW)
!          NULLIFY(REACDAT(IR)%RTCMW%ADAS)
!          NULLIFY(REACDAT(IR)%RTCMW%LINE)
!          NULLIFY(REACDAT(IR)%RTCMW%POLY)
!          NULLIFY(REACDAT(IR)%RTCMW%HYD)

          REACDAT(IR)%LRTCMW = .TRUE.
          REACDAT(IR)%RTCMW%IFIT = 5

!  ALREADY INITIALIZED IN EIRENE_ALLOC_FIT_FORM
!          REACDAT(IR)%RTCMW%RC1MIN = 0._DP
!          REACDAT(IR)%RTCMW%RC1MAX = HUGE(1._DP)
!          REACDAT(IR)%RTCMW%RC2MIN = 0._DP
!          REACDAT(IR)%RTCMW%RC2MAX = HUGE(1._DP)
!          REACDAT(IR)%RTCMW%FP1L = 0._DP
!          REACDAT(IR)%RTCMW%FP1R = 0._DP
!          REACDAT(IR)%RTCMW%FP2B = 0._DP
!          REACDAT(IR)%RTCMW%FP2T = 0._DP
!          REACDAT(IR)%RTCMW%JFEX1MN = 0
!          REACDAT(IR)%RTCMW%JFEX1MX = 0
!          REACDAT(IR)%RTCMW%JFEX2MN = 0
!          REACDAT(IR)%RTCMW%JFEX2MX = 0
          
        CASE (8:10)
          IF (REACDAT(IR)%LRTCEW) THEN
            WRITE (IUNOUT,*) ' ENERGY WEIGHTED RATE COEFFICIENT',
     .                       ' ALREADY SPECIFIED FOR REACTION', IR
            WRITE (IUNOUT,*) ' CHECK SPECIFICATION OF REACTIONS'
            CALL EIRENE_EXIT_OWN(1)
          END IF

          CALL EIRENE_ALLOC_FIT_FORM (REACDAT(IR)%RTCEW)
!  ALREADY DONE IN EIRENE_ALLOC_FIT_FORM
!          ALLOCATE (REACDAT(IR)%RTCEW)
!          NULLIFY(REACDAT(IR)%RTCEW%ADAS)
!          NULLIFY(REACDAT(IR)%RTCEW%LINE)
!          NULLIFY(REACDAT(IR)%RTCEW%POLY)
!          NULLIFY(REACDAT(IR)%RTCEW%HYD)

          REACDAT(IR)%LRTCEW = .TRUE.
          REACDAT(IR)%RTCEW%IFIT = 5

!  ALREADY INITIALIZED IN EIRENE_ALLOC_FIT_FORM
!          REACDAT(IR)%RTCEW%RC1MIN = 0._DP
!          REACDAT(IR)%RTCEW%RC1MAX = HUGE(1._DP)
!          REACDAT(IR)%RTCEW%RC2MIN = 0._DP
!          REACDAT(IR)%RTCEW%RC2MAX = HUGE(1._DP)
!          REACDAT(IR)%RTCEW%FP1L = 0._DP
!          REACDAT(IR)%RTCEW%FP1R = 0._DP
!          REACDAT(IR)%RTCEW%FP2B = 0._DP
!          REACDAT(IR)%RTCEW%FP2T = 0._DP
!          REACDAT(IR)%RTCEW%JFEX1MN = 0
!          REACDAT(IR)%RTCEW%JFEX1MX = 0
!          REACDAT(IR)%RTCEW%JFEX2MN = 0
!          REACDAT(IR)%RTCEW%JFEX2MX = 0
          
        CASE DEFAULT
          GOTO 1000         
        END SELECT
        RETURN

1000  continue
      WRITE (IUNOUT,*) ' ERROR IN "READ_COLRAD" : '
      WRITE (IUNOUT,*) ' WRONG DATA TYPE FOR INTERNAL COLRAD OPTION'
      WRITE (IUNOUT,*) ' REACTION NO. ', IR
      WRITE (IUNOUT,'(1X,A,I0)') ' DATA TYPE H.', ISW
      CALL EIRENE_EXIT_OWN(1)
      RETURN

      end subroutine EIRENE_read_colrad
