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
cdr                       H.4. 2.1.8 and H.10, 2.1.8, RC   recombination
cdr                       and  H.11, H.12: selected population coefficients
c 
c  to be done: units, log-lin, scaling, asymptotics
 
      use EIRMOD_precision
      use EIRMOD_parmmod
      use EIRMOD_comxs
      use EIRMOD_comprt, only: iunout
 
      implicit none
 
      integer, intent(in) :: ir, isw, iz1
      character(len=*), intent(in) :: reac
      integer, save :: ifirst
      integer, save :: ihsw(21)
      integer :: ivar, i, istr
      character(8), save :: hstr(21)

      close (29+ifoff)  ! nothing further to be read, currently

      if (ifirst == 0) then
        ifirst = 1
c  identifyers for data available from intrinsic CR code
c
c  h_colrad (atomic)

c   cr rate and e_rate, coupling to ground state H(1), ionisation
        ihsw(1) = 4
        hstr(1) = '2.1.5   ' 
        ihsw(2) = 10
        hstr(2) = '2.1.5   '
c   cr rate and e_rate, coupling to continuum H+, recombination 
        ihsw(3) = 4
        hstr(3) = '2.1.8   ' 
        ihsw(4) = 10
        hstr(4) = '2.1.8   ' 
c   cr rate and e_rate, coupling to radiation field, photo-excitation 
        ihsw(5) = 4
        hstr(5) = '2.1.5PH ' 
        ihsw(6) = 10
        hstr(6) = '2.1.5PH '
c  reduced population coefficient, H(n=3,2,4,5,6) states, coupling to ground state H(1), ionisation
        ihsw(7) = 12
        hstr(7) = '2.1.5a  '     ! n=3
        ihsw(8) = 12
        hstr(8) = '2.1.5b  '     ! n=2
        ihsw(9) = 12
        hstr(9) = '2.1.5c  '     ! n=4
        ihsw(10) = 12
        hstr(10) = '2.1.5d  '    ! n=5
        ihsw(11) = 12
        hstr(11) = '2.1.5e  '    ! n=6
c  reduced population coefficient, H(n=3,2,4,5,6) states, coupling to H+, recombination 
        ihsw(12) = 12
        hstr(12) = '2.1.8a  ' 
        ihsw(13) = 12
        hstr(13) = '2.1.8b  ' 
        ihsw(14) = 12
        hstr(14) = '2.1.8c  ' 
        ihsw(15) = 12
        hstr(15) = '2.1.8d  ' 
        ihsw(16) = 12
        hstr(16) = '2.1.8e  ' 
c  reduced population coefficient, H(n=3,2,4,5,6) states, coupling to radiation field 
        ihsw(17) = 12
        hstr(17) = '2.1.5PHa' 
        ihsw(18) = 12
        hstr(18) = '2.1.5PHb' 
        ihsw(19) = 12
        hstr(19) = '2.1.5PHc' 
        ihsw(20) = 12
        hstr(20) = '2.1.5PHd' 
        ihsw(21) = 12
        hstr(21) = '2.1.5PHe' 
      end if

cdr  error exit for unfinished options
      if (isw.ne.4 .and. isw.ne.10 .and. isw.ne.12)  goto 1000
cdr  tbd: also exit unless 2.1.5,  in particular: 
cdr       2.1.8 (recombination) is missing.
cdr  other reactions are not programmed in xsectp, rate-coeff, energy rate coef. 
      
cdr  IDENTIFY THE NUMBER IVAR (between 1:21) OF THE VARIABLE HSTR(IVAR) 
cdr  TO BE STORED ON M_HCOL(1:NHCOL_STORE).
        IVAR = 0
        DO I = 1, 21
           IF (ISW /= IHSW(I)) CYCLE
           IF (REAC(1:8) == HSTR(I)(1:8)) THEN
             IVAR = I
             EXIT
           END IF
        END DO
        IF (IVAR == 0) GOTO 1000

!  CHECK IF VARIABLE HAS ALREADY BEEN MARKED FOR STORING EARLIER
        ISTR = NHCOL_STORE + 1
        DO I = 1, NHCOL_STORE
          IF (IVAR == M_HCOL(I)) THEN
            ISTR = I
            EXIT
          END IF
        END DO
!  VARIABLE NOT YET MARKED FOR STORING --> MARK
        IF (ISTR > NHCOL_STORE) THEN
          NHCOL_STORE = ISTR
          M_HCOL(NHCOL_STORE) = IVAR
        END IF

!  ALREADY INITIALIZED IN EIRENE_INIT_CMDTA
!       REACDAT(IR)%ETH = 0._DP
!       REACDAT(IR)%RTMAX = 0._DP
!       REACDAT(IR)%ERTMAX = -HUGE(1._DP)

        SELECT CASE (ISW)
        CASE (2:4)
          IF (REACDAT(IR)%LRTC) THEN
            WRITE (IUNOUT,*) ' RATE COEFFICIENT ALREADY SPECIFIED',
     .                       ' FOR REACTION', IR
            WRITE (IUNOUT,*) ' CHECK SPECIFICATION OF REACTIONS'
            CALL EIRENE_EXIT_OWN(1)
          END IF

          CALL EIRENE_ALLOC_FIT_FORM (REACDAT(IR)%RTC)

          REACDAT(IR)%LRTC = .TRUE.
          REACDAT(IR)%RTC%IFIT = 5

          ALLOCATE (REACDAT(IR)%RTC%CRM)
          REACDAT(IR)%RTC%CRM%IFLAV = 1
          REACDAT(IR)%RTC%CRM%IVARST = ISTR
          
        CASE (5:7)
          IF (REACDAT(IR)%LRTCMW) THEN
            WRITE (IUNOUT,*) ' MOMENTUM WEIGHTED RATE COEFFICIENT',
     .                       ' ALREADY SPECIFIED FOR REACTION', IR
            WRITE (IUNOUT,*) ' CHECK SPECIFICATION OF REACTIONS'
            CALL EIRENE_EXIT_OWN(1)
          END IF

          CALL EIRENE_ALLOC_FIT_FORM (REACDAT(IR)%RTCMW)

          REACDAT(IR)%LRTCMW = .TRUE.
          REACDAT(IR)%RTCMW%IFIT = 5


          ALLOCATE (REACDAT(IR)%RTCMW%CRM)
          REACDAT(IR)%RTCMW%CRM%IFLAV = 1
          REACDAT(IR)%RTCMW%CRM%IVARST = ISTR  
        
        CASE (8:10)
          IF (REACDAT(IR)%LRTCEW) THEN
            WRITE (IUNOUT,*) ' ENERGY WEIGHTED RATE COEFFICIENT',
     .                       ' ALREADY SPECIFIED FOR REACTION', IR
            WRITE (IUNOUT,*) ' CHECK SPECIFICATION OF REACTIONS'
            CALL EIRENE_EXIT_OWN(1)
          END IF

          CALL EIRENE_ALLOC_FIT_FORM (REACDAT(IR)%RTCEW)

          REACDAT(IR)%LRTCEW = .TRUE.
          REACDAT(IR)%RTCEW%IFIT = 5


          ALLOCATE (REACDAT(IR)%RTCEW%CRM)
          REACDAT(IR)%RTCEW%CRM%IFLAV = 1
          REACDAT(IR)%RTCEW%CRM%IVARST = ISTR

        CASE (11:12)
          IF (REACDAT(IR)%LOTH) THEN
            WRITE (IUNOUT,*) ' OTHER RATE COEFFICIENT',
     .                       ' ALREADY SPECIFIED FOR REACTION', IR
            WRITE (IUNOUT,*) ' CHECK SPECIFICATION OF REACTIONS'
            CALL EIRENE_EXIT_OWN(1)
          END IF

          CALL EIRENE_ALLOC_FIT_FORM (REACDAT(IR)%OTH)

          REACDAT(IR)%LOTH = .TRUE.
          REACDAT(IR)%OTH%IFIT = 5

          ALLOCATE (REACDAT(IR)%OTH%CRM)
          REACDAT(IR)%OTH%CRM%IFLAV = 1
          REACDAT(IR)%OTH%CRM%IVARST = ISTR
          
        CASE DEFAULT
          GOTO 1000         
        END SELECT
        RETURN

1000  continue
      CALL EIRENE_LEER(1)
      WRITE (IUNOUT,*) ' ERROR IN "READ_COLRAD" : '
      WRITE (IUNOUT,*) ' WRONG DATA TYPE FOR INTERNAL COLRAD OPTION'
      WRITE (IUNOUT,*) ' REACTION NO. ', IR
      WRITE (IUNOUT,'(1X,A,I0)') ' DATA TYPE H.', ISW
      WRITE (IUNOUT,'(1X,A,A8)') ' DATA NR.    ', REAC(1:8)
      CALL EIRENE_EXIT_OWN(1)
      RETURN

      end subroutine EIRENE_read_colrad
