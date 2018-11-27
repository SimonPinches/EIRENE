!csw
!csw ADAS EIRENE MODULE
!csw (parts taken from NIMBUS)
!csw [19sep05] creation
!csw [09jan06] _eirene suffixes to all subroutines added
!csw
!csw s.wiesen
!csw
module mod_adas
  implicit none

  public :: adas_eirene

  character(len=80), public, save :: adas_eirene_userid

  private

  ! common cadas1
  INTEGER*4 :: IYEAR , IDYEAR
  CHARACTER(len=80) :: USERID

  ! common cadas2
  integer, parameter ::  ITDIM=55 , IDDIM=30 , IZDIM=28 , ISDIM=5 , IPDIM=7
  real*4 :: ZRAL  (ITDIM,IDDIM,IZDIM,ISDIM)
  real*4 ::  ZCAL  (ITDIM,IDDIM,IZDIM,ISDIM)
  real*4 ::  ZSAL  (ITDIM,IDDIM,IZDIM,ISDIM)
  real*4 ::  ZPLTL (ITDIM,IDDIM,IZDIM,ISDIM)
  real*4 ::  ZPLSL (ITDIM,IDDIM,IZDIM,ISDIM)
  real*4 :: ZPRBL (ITDIM,IDDIM,IZDIM,ISDIM)
  real*4 :: ZPRCL (ITDIM,IDDIM,IZDIM,ISDIM)
  real*4 :: ZSAL0 (ITDIM,IDDIM,ISDIM)
  real*4 :: ZPLTL0(ITDIM,IDDIM,ISDIM)
  real*4 :: ZPLSL0(ITDIM,IDDIM,ISDIM)
  real*4 :: ZTEL(ITDIM,IPDIM,ISDIM)
  real*4 :: ZNEL(IDDIM,IPDIM,ISDIM)

  INTEGER :: ITMAX(IPDIM,ISDIM),  IDMAX(IPDIM,ISDIM),  IZ0A(ISDIM),  ISMAX
  CHARACTER(len=8) :: ZLINFO(IZDIM,IPDIM,ISDIM)

contains

  SUBROUTINE ADAS_eirene( IZ0    , IYR  , IDYR  , USER , IZMAX, &
       IFORCE , ECUT , ICHAN , IOUT , IADAS, &
       TE     , TI, &
       DE     , NH   , DH    , HMASS, &
       SA0    , SA   , RTA   , PTA0 , PTA, &
       IER    )
    IMPLICIT NONE
    !C
    !C.......................................................................
    !C
    !C ROUTINE : ADAS
    !C
    !C VERSION : V1.R2.M0
    !C
    !C PURPOSE : TO ACCESS THE ADAS ATOMIC DATABASE.
    !C
    !C INPUT   : (I*4) IZ0          = NUCLEAR CHARGE OF ELEMENT.
    !C           (I*4) IYR          = YEAR OF ADAS FILES BEING ACCESSED.
    !C           (I*4) IDYR         = DEFAULT YEAR TO BE USED IF A REQUESTED
    !C                                ADAS 'IYR' FILE DOES NOT EXIST.
    !C           (C**) USER         = USERID WHERE ADAS FILES ARE STORED.
    !C           (I*4) IZMAX        = MAX. NUMBER OF CHARGED STATES ALLOWED
    !C           (I*4) IFORCE       = 0 --- DISALLOW 'ADASRE' ROUTINE TO READ
    !C                                      THE SAME 'IZ0' DATA TWICE.
    !C                              = 1 --- ALLOW 'ADASRE' ROUTINE TO READ
    !C                                      THE SAME 'IZ0' DATA AGAIN AND SO
    !C                                      OVERWRITE THE PREVIOUS READING.
    !C           (R*4) ECUT         = ENERGY CUT-OFF (EV) (USED TO ACCESS
    !C                                ENERGY FILTERED RADIATED POWER DATA
    !C                                SETS - SEE IEVCUT BELOW)
    !C           (I*4) ICHAN        = INPUT CHANNEL FOR USE IN READ ROUTINE
    !C           (I*4) IOUT         = OUTPUT CHANNEL FOR MESSAGES/ERRORS
    !C           (I*4) IADAS        = 0 --- USE 'ADASITX'
    !C                              = 1 --- USE 'ADASIT'
    !C
    !C           (R*4) TE           = ELECTRON TEMPERATURE (EV) AT WHICH THE
    !C                                ADAS COEFFICIENTS ARE REQUIRED
    !C           (R*4) TI           = ION TEMPERATURE (EV) AT WHICH THE
    !C                                ADAS COEFFICIENTS ARE REQUIRED
    !C           (R*4) DE           = ELECTRON DENSITY (CM-3) AT WHICH THE
    !C                                ADAS COEFFICIENTS ARE REQUIRED
    !C           (I*4) NH           = NO. OF NEUTRAL H ISOTOPES
    !C           (R*4) DH(NH)       = HYDROGEN DENSITY (CM-3) AT WHICH THE
    !C                                ADAS COEFFICIENTS ARE REQUIRED
    !C           (R*4) HMASS(NH)    = MASS OF HYDROGEN ISOTOPE (AMU)
    !C
    !C OUTPUT  : (R*4) SA0          = NEUTRAL IONIS. RATE COEFFT.(CM3 S-1)
    !C           (R*4) SA()         = ION IONIS. RATE COEFFT.    (CM3 S-1)
    !C           (R*4) RTA()        = ION RECOMB. RATE COEFFT.   (CM3 S-1)
    !C           (R*4) PTA0         = NEUTRAL RADIATED POWER  (WATTS CM3 S-1)
    !C           (R*4) PTA()        = ION RADIATED POWER      (WATTS CM3 S-1)
    !C
    !C           (I*4) IER          =     0 --- ROUTINE SUCCESSFUL
    !C                              = 10000 --- IZMAX .GT. IZDIM
    !C                              = 10001 --- IADAS ROUTINE DOES NOT EXIST
    !C                    (ADASRE)  = 10010 --- IZ0(INPUT) <> IZ(ADAS)
    !C                    (   "  )  = 10020 --- ITDIM SET TOO LOW
    !C                    (   "  )  = 10030 --- IDDIM SET TOO LOW
    !C                    (   "  )  = 10040 --- ISDIM SET TOO LOW
    !C                    (   "  )  = 10050 --- REQUESTED ELEMENT
    !C                                          NOT SUPPORTED
    !C                    (   "  )  = 10060 --- IYEAR AND IDYEAR FILES
    !C                                          NOT FOUND
    !C                    (ADASITX) = 10110 --- REQUIRED ELEMENT (IZ0) HAS
    !C                                          NOT BEEN EXTRACTED FROM ADAS
    !C                                          'ISONUCLEAR MASTER FILES'
    !C                    (   "   ) = 10120 --- FRACTH(1)+...+FRACTH(3)
    !C                                          DOES NOT SUM TO 1.0
    !C
    !C PROGRAM : (I*4) IZDIM        <= IZMAX
    !C           (R*4) HFRAC()      = HFRAC(1) : HYDROGEN  FUEL FRACTION
    !C                                HFRAC(2) : DEUTERIUM FUEL FRACTION
    !C                                HFRAC(3) : TRITIUM   FUEL FRACTION
    !C
    !C ROUTINE : ADASRE, ADASITX, ADASIT
    !C
    !C AUTHOR  : J.SPENCE  (K1/0/80)  EXT.4866
    !C           JET
    !C
    !C DATE    : V1.R1.M0 --- 07/04/94 --- CREATION
    !C           V1.R2.M0 --- 22/05/96 --- MORE IFAIL TRAPPING FROM 'ADASRE'
    !C
    !C.......................................................................
    !C
    INTEGER*4, parameter  ::  IZDIM=28
    !C
    character(len=80), intent(in) :: user
    INTEGER*4, intent(in) :: IZ0    , IYR        , IDYR       , IZMAX, &
         IFORCE , ICHAN      , IOUT       , IADAS, nh
    integer*4 :: I, IH
    !C
    real*4, intent(in) :: ecut,te,ti,de,dh(nh),hmass(nh)
    real*4, intent(out) :: sa0,sa(IZMAX),rta(IZMAX),pta0,pta(IZMAX)
    integer*4, intent(out) :: ier
    REAL*4 ::     HFRAC(3), dhtot, tran
    !C
    !C----------------------------- CHECKS & INITIALISE ---------------------
    !C
    IF( IZMAX.GT.IZDIM ) THEN
       WRITE(IOUT,1000) IZMAX , IZDIM
       IER = 10000
       return
    END IF
    !C
    IYEAR     = IYR
    IDYEAR    = IDYR
    USERID(1:80)    = USER(1:80)
    !C
    HFRAC(1) = 0.0E0
    HFRAC(2) = 0.0E0
    HFRAC(3) = 0.0E0
    DHTOT    = 0.0E0
    DO IH=1,NH
       IF( HMASS(IH).GT.0.9 .AND. HMASS(IH).LT.1.1 ) THEN
          !C         H
          DHTOT = DHTOT + DH(IH)
          HFRAC(1) = HFRAC(1)+DH(IH)
       ELSE IF( HMASS(IH).GT.1.9 .AND. HMASS(IH).LT.2.1 ) THEN
          !C         D
          DHTOT = DHTOT + DH(IH)
          HFRAC(2) = HFRAC(2)+DH(IH)
       ELSE IF( HMASS(IH).GT.2.4 .AND. HMASS(IH).LT.2.6 ) THEN
          !C         DT
          DHTOT = DHTOT + DH(IH)
          HFRAC(2) = HFRAC(2)+DH(IH)*0.5
          HFRAC(3) = HFRAC(3)+DH(IH)*0.5
       ELSE IF( HMASS(IH).GT.2.9 .AND. HMASS(IH).LT.3.1 ) THEN
          !C         T
          DHTOT = DHTOT + DH(IH)
          HFRAC(3) = HFRAC(3)+DH(IH)
       ENDIF
    ENDDO
    IF(NH.LE.0 .OR. DHTOT.LE.0.0) THEN
       !C       DUMMY
       HFRAC(1)=1.0
    ELSE
       TRAN = HFRAC(1)+HFRAC(2)+HFRAC(3)
       DO I=1,3
          HFRAC(I)=HFRAC(I)/TRAN
       ENDDO
    ENDIF
    !C
    !C----------------------------- READ ADAS -------------------------------
    !C
    CALL ADASRE_eirene( IZ0 , ECUT , ICHAN , IOUT , IFORCE , IER )
    !C
    IF( IER.NE.0 ) THEN
       IF( IER.EQ.1 ) THEN
          WRITE(IOUT,1010) IZ0
       ELSE IF( IER.EQ.2 ) THEN
          WRITE(IOUT,1020) 'ITDIM'
       ELSE IF( IER.EQ.3 ) THEN
          WRITE(IOUT,1020) 'IDDIM'
       ELSE IF( IER.EQ.4 ) THEN
          WRITE(IOUT,1020) 'ISDIM'
       ELSE IF( IER.EQ.5 ) THEN
          WRITE(IOUT,1030) IZ0
       ELSE IF( IER.EQ.6 ) THEN
          WRITE(IOUT,1040) IYR , IDYR
       ELSE IF( IER.EQ.7 ) THEN
          WRITE(IOUT,1080)
       ELSE IF( IER.EQ.8 ) THEN
          WRITE(IOUT,1090)
       ELSE
          WRITE(IOUT,1099) 'ADASRE' , IER
       END IF
       IER = 10000 + 10*IER
       return
    END IF
    !C
    !C----------------------------- USE  ADAS -------------------------------
    !C
    IF( IADAS.EQ.0 ) THEN
       CALL ADASITX_eirene( IZ0 , IOUT, &
            TE  , TI    , DE     , DHTOT, HFRAC(1), &
            SA0 , SA(1) , RTA(1) , PTA0 , PTA(1), &
            IER  )
    ELSE IF( IADAS.EQ.1 ) THEN
       CALL ADASIT_eirene ( IZ0 , IOUT, &
            TE  , TI    , DE     , DHTOT, HFRAC(1),&
            SA0 , SA(1) , RTA(1) , PTA0 , PTA(1),&
            IER  )
    ELSE
       WRITE(IOUT,1050)
       IER = 10001
       return
    END IF
    !C
    IF( IER.NE.0 ) THEN
       IF( IER.EQ.1 ) THEN
          WRITE(IOUT,1060) IZ0
       ELSE IF( IER.EQ.2 ) THEN
          WRITE(IOUT,1070) HFRAC(1)+HFRAC(2)+HFRAC(3)
       ELSE
          WRITE(IOUT,1099) 'ADASITX' , IER
       END IF
       IER = 10100 + 10*IER
       return
    END IF
    !C
    !C----------------------------- FORMAT ----------------------------------
    !C
1000 FORMAT( / ' *** ERROR *** ADAS : IZMAX (',I3,') .GT. IZDIM (',I3 , ')' / '                      POSSIBLE OVERWRITE PROBLEM' / )
1010 FORMAT( / ' *** ERROR *** ADASRE : IZ0 (',I3,') <> IZ (ADAS)' / )
1020 FORMAT( / ' *** ERROR *** ADASRE : ',A,' SET TOO LOW.' / '                        CONSULT ADAS SUPERVISOR.' / )
1030 FORMAT( / ' *** ERROR *** ADASRE : IZ0 (',I3,') NOT SUPPORTED' / )
1040 FORMAT( / ' *** ERROR *** ADASRE : IYEAR (',I3,') & IDYEAR (',I3, ') FILES NOT FOUND' / )
1050 FORMAT( / ' *** ERROR *** ADAS : IADAS (',I3,') ROUTINE IS NOT SUPPORTED.' / )
1060 FORMAT( / ' *** ERROR *** ADASITX : IZ0 (',I3,') HAS NOT BEEN EXTRACTED FROM ADAS FILES' / )
1070 FORMAT( / ' *** ERROR *** ADASITX : HYDROGEN ISOTOPIC FRACTIONS (', 1P , E12.4 , ') DO NOT SUM TO 1.0' / )
1080 FORMAT( / ' *** ERROR *** ADASRE : UNABLE TO OPEN ADAS DATAFILE'/)
1090 FORMAT( / ' *** ERROR *** ADASRE : NO. OF SPECIES > MAX ALLOWED'/)
1099 FORMAT( / ' *** ERROR *** ' , A , ' : IFAIL = ' , I6, ' CONSULT ADAS SUPERVISOR.' / )
    !C
    !C----------------------------- FORMAT ----------------------------------
    !C
    IER   = 0
    !C
    RETURN
  END SUBROUTINE ADAS_EIRENE


  SUBROUTINE ADASRE_eirene( IZ0 , ECUT , ICHAN , IOUT , IFORCE , IER )
    IMPLICIT NONE
    !C
    !C+ .....................................................................
    !C
    !C ROUTINE : ADAS READING OF DATAFILES
    !C           ---  --
    !C VERSION : V2.R1.M0
    !C
    !C PURPOSE : TO READ A COMPLETE SET OF 'ADAS' MASTER FILES FOR AN ELEMENT
    !C
    !C INPUT   : (I*4) IZ0          = NUCLEAR CHARGE OF ELEMENT TO BE READ
    !C           (R*4) ECUT         = ENERGY CUT-OFF (EV) (USED TO ACCESS
    !C                                ENERGY FILTERED RADIATED POWER DATA
    !C                                SETS - SEE IEVCUT BELOW)
    !C           (I*4) ICHAN        = N/A
    !C           (I*4) IOUT         = OUTPUT CHANNEL FOR MESSAGES/ERRORS
    !C           (I*4) IFORCE       = 0 --- ALLOW ROUTINE TO DETERMINE
    !C                                      WHETHER TO READ DATA OR NOT.
    !C                                      ( E.G DISALLOWING THE READING
    !C                                            OF THE SAME SPECIES TWICE)
    !C                              = 1 --- FORCE ROUTINE TO READ DATA AND
    !C                                      OVERWRITE THE OLD SPECIES DATA
    !C                                      WITH THE NEW ONE.
    !C                                      ( E.G READ THE SAME SPECIES AGAIN
    !C                                            BUT FROM A DIFFERENT YEAR )
    !C
    !C OUTPUT  : (I*4) IER          = 0 --- ROUTINE SUCCESSFUL
    !C                              = 2 --- NO. ADAS TEMP'S > MAX ALLOWED
    !C                              = 3 --- NO. ADAS DENS'S > MAX ALLOWED
    !C                              = 4 --- NO. ADAS STATES > MAX ALLOWED
    !C                              = 7 --- UNABLE TO OPEN A DATAFILE
    !C                              = 8 --- SPECIES EXCEEDS MAXIMUM ALLOWED
    !C
    !C PROGRAM : (I*4)  ICLMAX      = MAX. ICLASS ALLOWED
    !C           (C*2)  YEAR        = IYEAR
    !C           (C*2)  YEARDF      = IDYEAR
    !C           (I*4)  ISOLD       = ORIGINAL INDEX OF SPECIES WHICH WILL
    !C                                BE OVERWRITTEN BY NEW DATA (IFORCE=1)
    !C           (I*4)  IEVCUT      = ABSOLUTE INTEGER PART OF ENERGY CUTOF
    !C           (I*4)  IFAIL       = D2LINK ERROR CODES
    !C           (C*100)MESS        = MESSAGE FROM D2DATA
    !C           (LOG)  LINIT       = .TRUE. --- INITIALISE /CADAS2/
    !C
    !C NOTES   : (1) FROM COMMON BLOCK 'CADAS1' :-
    !C               (I*4)  IYEAR  = 'ADAS' YEAR IDENTIFIER REQUESTED
    !C               (I*4)  IDYEAR = 'ADAS' YEAR IDENTIFIER DEFAULT
    !C               (C*6)  USERID = ON WHO'S FILE SPACE ARE 'ADAS' FILES
    !C
    !C           (2) IF IYEAR = 0 THEN IDYEAR IS USED
    !C
    !C           (3) IF THERE IS A MISSING 'ADAS' FILE FOR IYEAR THEN
    !C               THE 'ADAS' FILE FROM IDYEAR WILL BE USED INSTEAD.
    !C
    !C           (4) A WARNING MESSAGE WILL BE GIVEN IF A IDYEAR IS
    !C               USED AS A REPLACEMENT FOR IYEAR WHEN IYEAR <> 0.
    !C
    !C           (5) ADAS DATAFILE NAMES ;
    !C    ICLASS=1 ... ACD    = COLLISIONAL DIELECTRONIC
    !C                          RECOMBINATION RATE COEFFICIENTS (CM**3 SEC-1)
    !C    ICLASS=2 ... SCD    = COLLISIONAL DIELECTRONIC
    !C                          IONISATION RATE COEFFICIENTS (CM**3 SEC-1)
    !C    ICLASS=3 ... CCD    = COLLISIONAL RADIATIVE CHARGE EXCHANGE
    !C                          RECOMBINATION RATE COEFFICIENTS CM**3 SEC-1)
    !C    ICLASS=4 ... PLT    = TOTAL LOW LEVEL LINE POWER
    !C                          LOSS RATE COEFFICIENTS (ERG CM**3 SEC-1)
    !C    ICLASS=5 ... PRB    = RECOMBINATION/CASCADE + BREMSSTRAHLUNG POWER
    !C                          LOSS RATE COEFFICIENTS (ERG CM**3 SEC-1)
    !C    ICLASS=6 ... PRC    = CHARGE EXCHANGE POWER
    !C                          LOSS RATE COEFFICIENTS (ERG CM**3 SEC-1)
    !C    ICLASS=7 ... PLS    = SPECIFIC LINE POWER
    !C                          LOSS RATE COEFFICIENTS (ERG CM**3 SEC-1)
    !C
    !C ROUTINE : D2LINK       = CALL D2DATA AND THEN FILL /CADAS2/
    !C           I4UNIT       = SET MESSAGE OUTPUT CHANNEL FOR D2DATA
    !C           XXUID        = SET ADAS DATA SOURCE USER ID
    !C
    !C AUTHOR  : JAMES SPENCE  (K1/0/80)  EXT. 4865
    !C           JET
    !C
    !C HISTORY : V2.R1.M0 --- 16/05/96 --- CREATION
    !C
    !C- .....................................................................
    !C
    !C..INPUT
    INTEGER*4 :: IZ0 , ICHAN , IOUT , IFORCE
    REAL*4 :: ECUT
    !C
    !C..OUTPUT
    INTEGER*4 :: IER
    !C..PROGRAM
    INTEGER*4, parameter :: ICLMAX=7
    INTEGER*4 :: IS     , ISOLD    , IEVCUT , IFAIL , ICLASS
    CHARACTER(len=2) :: YEAR , YEARDF
    character(len=100) :: MESS
    LOGICAL, save :: LINIT=.true.
    !C
    !C-------------------------------- INITIALISE ---------------------------
    !C
    IF( LINIT ) THEN
       ISMAX = 0
       LINIT = .FALSE.
    END IF
    !C
    IER      = 0
    ISOLD    = 0
    !C
    CALL I4UNIT(IOUT)
    !C
    CALL XXUID(USERID)
    !C
    !C--------------------- CONVERT YEAR INTEGER TO STRING ------------------
    !C
    WRITE(YEARDF,'(I2)') IDYEAR
    !C
    IF( IYEAR.GT.0 ) THEN
       WRITE(YEAR  ,'(I2)') IYEAR
    ELSE
       YEAR = YEARDF
    END IF
    !C
    !C---------------------- CONVERT ECUT REAL TO INTEGER -------------------
    !C
    IEVCUT = ABS(ECUT)
    !C
    !C----------------------------- SPECIES LOADING -------------------------
    !C
    !C.. SKIP ROUTINE IF THE SPECIES HAS ALREADY BEEN READ OR FORCE RE-READ.
    DO  IS = 1 , ISMAX
       IF( IZ0A(IS).EQ.IZ0 ) THEN
          IF( IFORCE.GE.1 ) THEN
             ISOLD = IS
          ELSE
             return
          END IF
       END IF
    enddo
    !C
    !C.. IS THERE STORAGE AVAILABLE FOR THE SPECIES ?
    IF( ISOLD.EQ.0 ) THEN
       IS   = ISMAX + 1
       IF( IS.GT.ISDIM ) THEN
          write(iout,'(a,i3,a)') ' *** ERROR *** SPECIES DIMENSION (ISDIM=',isdim,') : TOO LOW'
          write(iout,'(a)')      '               THE DIMENSION MUST BE, AT LEAST, EQUAL TO '
          write(iout,'(a)')      '               THE NUMBER OF SPECIES YOU WANT TO ACCESS.'
          IER = 7
          return
       END IF
    ELSE
       IS   = ISOLD
       write(iout,'(a,i3,a,i2,a)') ' *** WARNING *** SPECIES INDEX =',ISOLD,' (IZ0 = ',IZ0,')'
       write(iout,'(a)')           '                 OVERWRITTEN WITH NEW ADAS DATA'
    END IF
    !C
    !C----------------------------- READ ADAS FILES -------------------------
    !C
    DO ICLASS = 1 , ICLMAX
       CALL D2LINK_eirene( ICLASS , YEAR  , YEARDF , IZ0 , IS , IEVCUT , IOUT , MESS   , IFAIL )
       IF( IFAIL.NE.0 ) THEN
          WRITE(IOUT,*) MESS
          return
       END IF
    enddo
    !C
    !C----------------------------- SPECIES STORAGE -------------------------
    !C
    IF( ISOLD.EQ.0 ) THEN
       ISMAX       = IS
       IZ0A(ISMAX) = IZ0
    END IF
    RETURN
  END SUBROUTINE ADASRE_EIRENE


  !C=======================================================================
  SUBROUTINE D2LINK_eirene( ICLASS , YEAR  , YEARDF , IZ0    , IS    , IEVCUT, IOUT , MESS   , IFAIL )
    IMPLICIT NONE
    !C
    !C+ .....................................................................
    !C
    !C ROUTINE : ADAS READING OF DATAFILES
    !C           ---  --
    !C VERSION : V1.R1.M0
    !C
    !C PURPOSE : TO READ A COMPLETE SET OF 'ADAS' MASTER FILES FOR AN ELEMENT
    !C           USING D2DATA AND STORE IN /CADAS2/
    !C
    !C INPUT   : (I*4) ICLASS       = 1 --- READ 'ACD' DATA
    !C                              = 2 --- READ 'SCD' DATA
    !C                              = 3 --- READ 'CCD' DATA
    !C                              = 4 --- READ 'PRB' DATA
    !C                              = 5 --- READ 'PLT' DATA
    !C                              = 6 --- READ 'PRC' DATA
    !C                              = 7 --- READ 'PLS' DATA
    !C           (C*2) YEAR         = 'ADAS' YEAR REQUESTED (YY)
    !C           (C*2) YEARDF       = 'ADAS' YEAR DEFAULT (YY)
    !C           (I*4) IZ0          = NUCLEAR CHARGE OF ELEMENT TO BE READ
    !C           (I*4) IS           = SPECIES INDEX
    !C           (I*4) IEVCUT       = ENERGY CUT-OFF (EV) (USED TO ACCESS
    !C                                ENERGY FILTERED RADIATED POWER DATA
    !C                                SETS (FOR ICLASS=4,5,6)
    !C           (I*4) IOUT         = OUTPUT CHANNEL FOR MESSAGES/ERRORS
    !C
    !C OUTPUT  : (C**) MESS         = MESSAGE FROM D2DATA
    !C           (I*4) IFAIL        = 0 --- ROUTINE SUCCESSFUL
    !C                              = 2 --- NO. ADAS TEMP'S > MAX ALLOWED
    !C                              = 3 --- NO. ADAS DENS'S > MAX ALLOWED
    !C                              = 4 --- NO. ADAS STATES > MAX ALLOWED
    !C                              = 6 --- UNABLE TO OPEN A DATAFILE
    !C                              = 8 --- ICLASS NOT SUPPORTED
    !C                              = 9 --- NOT ALL IZ0 CHARGED STATES READ
    !C
    !C PARAM   : (I*4) ITDIM        = MAX. NUMBER OF TEMPERATURES    ALLOWED
    !C           (I*4) IDDIM        = MAX. NUMBER OF DENSITIES       ALLOWED
    !C           (I*4) IZDIM        = MAX. NUMBER OF CHARGED STATES  ALLOWED
    !C           (I*4) IPDIM        = MAX. NUMBER OF ADAS FILE TYPES ALLOWED
    !C
    !C /CADAS2/: (R*4) ZRAL(I,J,Z,S) =RADIATIVE + DIELECTRONIC RECOMBINATION
    !C           (R*4) ZCAL(I,J,Z,S) =CHARGE EXCHANGE RECOMBINATION
    !C           (R*4) ZSAL0(I,J,S)  =IONISATION (NEUTRAL)
    !C           (R*4) ZSAL(I,J,Z,S) =IONISATION (NON-NEUTRAL)
    !C           (R*4) ZPLTL0(I,J,S) =TOTAL LINE RADIATED POWER (NEUTRAL)
    !C           (R*4) ZPLTL(I,J,Z,S)=TOTAL LINE RADIATED POWER (NON-NEUTRAL
    !C           (R*4) ZPRBL(I,J,Z,S)=RAD. + DIEL RECOM. + BREMS. POWER
    !C           (R*4) ZPRCL(I,J,Z,S)=CX. RECOM. POWER
    !C           (R*4) ZPLSL0(I,J,S) =SPECIFIC LINE POWER (NEUTRAL)
    !C           (R*4) ZPLSL(I,J,Z,S)=SPECIFIC LINE POWER (NON-NEUTRAL)
    !C           (I*4) ITMAX(P,S)    =NUMBER OF TEMPERATURES (EV)
    !C           (I*4) IDMAX(P,S)    =NUMBER OF TEMPERATURES (CM-3)
    !C           (I*4) IZ0A(S)       =NUCLEAR CHARGE
    !C           (I*4) ISMAX         =NUMBER OF SPECIES READ FROM ADAS
    !C           (C*8) ZLINFO(Z,P,S) =WAVELENGTH INFORMATION LINE
    !C
    !C                 PARAMETERS:
    !C                 -----------
    !C           (I*4) ITDIM        = MAX. NUMBER OF TEMPERATURES    ALLOWED
    !C           (I*4) IDDIM        = MAX. NUMBER OF DENSITIES       ALLOWED
    !C           (I*4) IZDIM        = MAX. NUMBER OF CHARGED STATES  ALLOWED
    !C           (I*4) IPDIM        = MAX. NUMBER OF ADAS FILE TYPES ALLOWED
    !C
    !C                 INDEXING:
    !C                 ---------
    !C                 I            =  TE        (1 <= I >= ITMAX(P,S))
    !C                 J            =  DENS      (1 <= J >= IDMAX(P,S))
    !C                 Z            =  STAGE     (1 <= Z >= IZ0A(S))
    !C                 S            =  SPECIES   (1 <= S >= ISMAX)
    !C                 P            =  ADAS FILE (1 <= K >= IPDIM)
    !C
    !C                 *    ALL TABLES ARE LOG10  IN CGS UNITS
    !C
    !C                 *    ENERGY LOWER LIMIT  EGR  OPERATES
    !C
    !C PROGRAM : (R*8)  DTEVD()     = ADAS DATAFILE TEMPERATURES (EV)
    !C           (R*8)  DDENSD()    = ADAS DATAFILE DENSITIES (CM-3)
    !C       (R*8) DCOEFD(IZ,IT,ID) = ADAS DATAFILE COEFFICIENTS
    !C           (R*8)  ZDATA()     = ADAS DATAFILE CHARGED STATE
    !C
    !C           (I*4)  INTP        = NO. OF (TE,DE) INTERPOLATION POINTS
    !C                                FOR S.R. D2DATA
    !C           (R*8)  TEINT()     = TEMPERATURE INTERPOLATION POINTS (EV)
    !C           (R*8)  DEINT()     = DENSITY INTERPOLATION POINTS (CM-3)
    !C           (R*8)  COFINT()    = INTERPOLATED COEFFICIENT FOR (TE,DE)
    !C
    !C           (I*4)  IT          = LOOP ARRAY FOR TEMPERATURES
    !C           (I*4)  ID          = LOOP ARRAY FOR DENSITIES
    !C           (I*4)  IZ          = LOOP ARRAY FOR CHARGED STATES
    !C           (I*4)  IP          = LOOP ARRAY FOR POINTS
    !C           (I*4)  ICTOIP()    = MAP ICLASS TO IPDIM ELEMENT IN /CADAS2/
    !C
    !C NOTES   : (1) A WARNING MESSAGE WILL BE GIVEN IF A IDYEAR IS
    !C               USED AS A REPLACEMENT FOR IYEAR WHEN IYEAR <> 0.
    !C
    !C           (2) ADAS DATAFILE NAMES ;
    !C    ICLASS=1 ... ACD    = COLLISIONAL DIELECTRONIC
    !C                          RECOMBINATION RATE COEFFICIENTS (CM**3 SEC-1)
    !C    ICLASS=2 ... SCD    = COLLISIONAL DIELECTRONIC
    !C                          IONISATION RATE COEFFICIENTS (CM**3 SEC-1)
    !C    ICLASS=3 ... CCD    = COLLISIONAL RADIATIVE CHARGE EXCHANGE
    !C                          RECOMBINATION RATE COEFFICIENTS CM**3 SEC-1)
    !C    ICLASS=4 ... PLT    = TOTAL LOW LEVEL LINE POWER
    !C                          LOSS RATE COEFFICIENTS (ERG CM**3 SEC-1)
    !C    ICLASS=5 ... PRB    = RECOMBINATION/CASCADE + BREMSSTRAHLUNG POWER
    !C                          LOSS RATE COEFFICIENTS (ERG CM**3 SEC-1)
    !C    ICLASS=6 ... PRC    = CHARGE EXCHANGE POWER
    !C                          LOSS RATE COEFFICIENTS (ERG CM**3 SEC-1)
    !C    ICLASS=7 ... PLS    = SPECIFIC LINE POWER
    !C                          LOSS RATE COEFFICIENTS (ERG CM**3 SEC-1)
    !C
    !C ROUTINE : D2DATA       = READ MASTER FILE
    !C           I4UNIT       = SET MESSAGE OUTPUT CHANNEL FOR D2DATA
    !C
    !C AUTHOR  : JAMES SPENCE  (K1/0/80)  EXT. 4865
    !C           JET
    !C
    !C HISTORY : V1.R1.M0 --- 17/05/96 --- CREATION
    !C
    !C- .....................................................................
    !C
    !C..INPUT
    INTEGER*4 ::  ICLASS   , IZ0 , IS   , IEVCUT , IOUT
    CHARACTER(len=2) ::  YEAR, YEARDF
    !C
    !C..OUTPUT
    INTEGER*4 ::   IFAIL
    CHARACTER( len=*) ::  MESS
    !C..PROGRAM
    INTEGER*4, parameter ::  INTP =1
    REAL*8 ::     DTEVD(ITDIM)                , DDENSD(ITDIM)
    real*8 ::     DCOEFD(ITDIM,ITDIM,ITDIM)   , ZDATA(ITDIM)
    real*8 ::     COFINT(INTP)
    INTEGER*4 ::  IT     , ID       , IZ      , IP ,  ITMAXD , IDMAXD   , IZMAXD
    !C
    !C--- FIXED INTERPOLATION VALUES AS D2DATA INTERPOLATION IS NOT WANTED --
    !C
    real*8 ::     TEINT(INTP)=1.0e00                 , DEINT(INTP)=1.e13
    integer*4 ::  ICTOIP(7)=(/ 1 , 3 , 2 , 5 , 4 , 6 , 7 /)
    !C
    !C-------------------------------- INITIALISE ---------------------------
    !C
    IFAIL  = 0
    !C
    ITMAXD = 0
    IDMAXD = 0
    IZMAXD = 0
    !C
    !C-------------------------------- READ ADAS ----------------------------
    !C
    CALL D2DATA_eirene( YEAR          , YEARDF    , MESS   , IFAIL &
         , IZ0           , 0         , ICLASS , INTP   , IEVCUT &
         , ITDIM         , ITMAXD    , IDMAXD , IZMAXD &
         , TEINT(1)      , DEINT(1) &
         , DTEVD(1)      , DDENSD(1) &
         , DCOEFD(1,1,1) , ZDATA(1) &
         , COFINT(1)     )
    !C
    !C-------------------------------- CHECK DATA ---------------------------
    !C
    IF( IFAIL.LT.0 ) IFAIL = 0
    !C
    IF( IFAIL.EQ.1 ) IFAIL = 6
    !C
    IF( ITMAXD.GT.ITDIM) THEN
       IFAIL = 2
       WRITE(MESS,9000) 'ITMAXD' , ITMAXD , 'ITDIM' , ITDIM
    ELSE IF( IDMAXD.GT.IDDIM ) THEN
       IFAIL = 3
       WRITE(MESS,9000) 'IDMAXD' , IDMAXD , 'IDDIM' , IDDIM
    ELSE IF( IZMAXD.GT.IZDIM ) THEN
       IFAIL = 4
       WRITE(MESS,9000) 'IZMAXD' , IZMAXD , 'IZDIM' , IZDIM
    END IF
    !C
    IF( ZDATA(1).NE.1 .OR. ZDATA(IZMAXD).NE.IZ0 ) THEN
       IFAIL = 9
       WRITE(MESS,9020)
    END IF
    !C
    IF( IFAIL.GT.0 ) return
    !C
    !C------------------------------ FILL /CADAS2/ --------------------------
    !C
    IP           = ICTOIP(ICLASS)
    !C
    ITMAX(IP,IS) = ITMAXD
    IDMAX(IP,IS) = IDMAXD
    !C
    DO IT = 1 , ITMAXD
       ZTEL(IT,IP,IS) = DTEVD(IT)
    enddo
    !C
    DO ID = 1 , IDMAXD
       ZNEL(ID,IP,IS) = DDENSD(ID)
    enddo
    !C
    DO IZ = 1 , IZMAXD
       READ(MESS,'(44X,A)') ZLINFO(IZ,IP,IS)
       !C.. <<PATCH FOR LACK OF HYDROGEN ISOTOPE MASS IN SOME CX MASTER FILES>>
       IF( ICLASS.EQ.3 .AND. ZLINFO(IZ,IP,IS)(1:4).NE.' MH=' ) ZLINFO(IZ,IP,IS) = ' MH=2.00'
       !C..
    enddo
    !C
    DO IT       = 1 , ITMAXD
       DO ID    = 1 , IDMAXD
          DO IZ = 1 , IZMAXD
             !C
             !C.. COLLISIONAL DILELECTRONIC - RADIATIVE RECOMBINATION RATE
             IF( ICLASS.EQ.1 ) THEN
                ZRAL(IT,ID,IZ,IS)    = DCOEFD(IZ,IT,ID)
                !C
                !C.. IONISATION RATE COEFFICIENT
             ELSE IF( ICLASS.EQ.2 .AND. IZ.EQ.1 ) THEN
                ZSAL0(IT,ID,IS)      = DCOEFD(1,IT,ID)
             ELSE IF( ICLASS.EQ.2 .AND. IZ.GT.1 ) THEN
                ZSAL(IT,ID,IZ-1,IS)  = DCOEFD(IZ,IT,ID)
                !C
                !C.. CHARGE EXCHANGE RECOMBINATION
             ELSE IF( ICLASS.EQ.3 ) THEN
                ZCAL(IT,ID,IZ,IS)    = DCOEFD(IZ,IT,ID)
                !C
                !C.. RAD. + DIEL RECOMB. + BREMS. POWER
             ELSE IF( ICLASS.EQ.4 ) THEN
                ZPRBL(IT,ID,IZ,IS)   = DCOEFD(IZ,IT,ID)
                !C
                !C.. TOTAL LINE RADIATIVE POWER
             ELSE IF( ICLASS.EQ.5 .AND. IZ.EQ.1 ) THEN
                ZPLTL0(IT,ID,IS)     = DCOEFD(1,IT,ID)
             ELSE IF( ICLASS.EQ.5 .AND. IZ.GT.1 ) THEN
                ZPLTL(IT,ID,IZ-1,IS) = DCOEFD(IZ,IT,ID)
                !C
                !C.. CHARGE EXCHANGE RECOMBINATION POWER
             ELSE IF( ICLASS.EQ.6 ) THEN
                ZPRCL(IT,ID,IZ,IS)   = DCOEFD(IZ,IT,ID)
                !C
                !C.. SPECIFIC LINE POWER
             ELSE IF( ICLASS.EQ.7 .AND. IZ.EQ.1 ) THEN
                ZPLSL0(IT,ID,IS)     = DCOEFD(1,IT,ID)
             ELSE IF( ICLASS.EQ.7 .AND. IZ.GT.1 ) THEN
                ZPLSL(IT,ID,IZ-1,IS) = DCOEFD(IZ,IT,ID)
                !C
                !C.. ICLASS VALUE UNSUPPORTED !
             ELSE
                IFAIL = 8
                WRITE(MESS,9010) ICLASS
                return
             END IF
          enddo
       enddo
    enddo
    !C
    !C-----------------------------------------------------------------------
    !C
9000 FORMAT( ' *** ERROR *** D2DATA : ' , A , ' (' , I5, ') > ' , A , ' (' , I5 ,')' )
9010 FORMAT( ' *** ERROR *** D2LINK : ICLASS (' , I5, ') NOT SUPPORTED' )
9020 FORMAT( ' *** ERROR *** D2LINK : INCOMPLETE CHARGE STATE RANGE', ' READ FROM ADAS' )
    !C
    !C-----------------------------------------------------------------------
    !C
    RETURN
  END SUBROUTINE D2LINK_EIRENE

!C=======================================================================
  SUBROUTINE ADASIT_eirene( IZ0     , IOUT , TE      , TI     , DENS   , DENSH , FRACTH, SA0     , SA     , RTA    , PTA0  , PTA, IFAIL)
    IMPLICIT REAL*4(A-H,O-Z)
    !C
    !C-----------------------------------------------------------------------
    !C
    !C
    !C PURPOSE : TO INTERPOLATE 'ADAS' ISONUCLEAR MASTER FILES AT THE
    !C           REQUIRED ELECTRON TEMPERATURE, ION TEMPERATURE & ELECTRON
    !C           DENSITY.  THE REQUIRED ELEMENT MUST ALREADY HAVE BEEN
    !C           EXTRACTED FROM THE 'ADAS' ISONUCLEAR MASTER FILES.
    !C
    !C
    !C TYPE    : SINGLE PRECISION
    !C
    !C
    !C INPUT   : (I*4) IZ0          = NUCLEAR CHARGE OF REQUIRED ELEMENT
    !C           (I*4) IOUT         = OUTPUT CHANNEL FOR MESSAGES/ERRORS
    !C           (R*4) TE           = ELECTRON TEMPERATURE (EV) AT WHICH THE
    !C                                ADAS COEFFICIENTS ARE REQUIRED
    !C           (R*4) TI           = ION TEMPERATURE (EV) AT WHICH THE
    !C                                ADAS COEFFICIENTS ARE REQUIRED
    !C           (R*4) DENS         = ELECTRON DENSITY (CM-3) AT WHICH THE
    !C                                ADAS COEFFICIENTS ARE REQUIRED
    !C           (R*4) DENSH        = HYDROGEN DENSITY (CM-3) AT WHICH THE
    !C                                ADAS COEFFICIENTS ARE REQUIRED
    !C           (R*4) FRACTH()     = FRACTH(1) : HYDROGEN  FUEL FRACTION
    !C                                FRACTH(2) : DEUTERIUM FUEL FRACTION
    !C                                FRACTH(3) : TRITIUM   FUEL FRACTION
    !C
    !C
    !C OUTPUT  : (R*4) SA0          = NEUTRAL IONIS. RATE COEFFT.(CM3 S-1)
    !C           (R*4) SA()         = ION IONIS. RATE COEFFT.    (CM3 S-1)
    !C           (R*4) RTA()        = ION RECOMB. RATE COEFFT.   (CM3 S-1)
    !C           (R*4) PTA0         = NEUTRAL RADIATED POWER  (WATTS CM3 S-1)
    !C           (R*4) PTA()        = ION RADIATED POWER      (WATTS CM3 S-1)
    !C           (I*4) IFAIL        = 0 --- ROUTINE SUCCESSFUL
    !C                              = 1 --- THE REQUIRED ELEMENT (IZ0) HAS
    !C                                      NOT BEEN EXTRACTED FROM THE ADAS
    !C                                      'ISONUCLEAR MASTER FILES'
    !C                              = 2 --- FRACTH(1)+...+FRACTH(MHYD) IS
    !C                                      NOT A REAL NUMBER.
    !C
    !C
    !C COMMON  : (I*4) ITDIM        = MAX. NUMBER OF TEMPERATURES   ALLOWED
    !C           (I*4) IDDIM        = MAX. NUMBER OF DENSITIES      ALLOWED
    !C           (I*4) IZDIM        = MAX. NUMBER OF CHARGED STATES ALLOWED
    !C           (I*4) ISDIM        = MAX. NUMBER OF SPECIES        ALLOWED
    !C           (I*4) IPDIM        = MAX. NUMBER OF ADAS FILES     ALLOWED
    !C
    !C
    !C PROGRAM : (I*4) MHYD         = NUMBER OF HYDROGEN ISOTOPES
    !C           (I*4) IS           = INDEX OF ELEMENT STORED IN /CADAS2/
    !C           (R*4) FSUM         = SUM OF HYDROGEN ISOTOPIC FRACTIONS
    !C           (R*4) HMADAS       = HYDROGEN ISOTOPE USED IN MASTER CX FILE
    !C                                (1=PROTIUM, 2=DEUTERIUM, 3=TRITIUM)
    !C           (C*8) STRING       = WORK STRING
    !C
    !C
    !C NOTES   : (1) THE HYDROGEN ISOTOPIC FRACTIONS ARE NORMALISED FOR
    !C               THE CODE IN CASE THEY ARE NOT ON INPUT.
    !C           (2) THE ARRAY STACKING IS UNUSUAL!
    !C               FOR RA, CA, PRB, PRC DATA IS STORED AS FOLLOWS:
    !C                RECOMBINING ION CHARGE   RECOMBINED ION CHARGE  INDEX
    !C                ----------------------   ---------------------  -----
    !C                     1                       0 (NEUTRAL)          1
    !C                     IZ                      IZ-1                 IZ
    !C                     IZ0 (BARE NUCLEUS)      IZ0-1                IZ0
    !C
    !C               FOR SA, PLT, PLS DATA IS STORED AS FOLLOWS:
    !C                IONISING  ION  CHARGE    IONISED  ION  CHARGE   INDEX
    !C                ----------------------   ---------------------  -----
    !C                     0 (NEUTRAL)             1                 SEPARATE
    !C                     IZ                      IZ+1                 IZ
    !C                     IZ0-1                   IZ0 (BARE NUCLEUS)   IZ0-1
    !C
    !C
    !C               THIS IS DIFFERENT FROM THE CONVENTIONS ACCORDING TO
    !C               ADF12 IN THE ADAS MASTER FILES.
    !C
    !C AUTHOR  : JAMES SPENCE  (K1/0/80)  EXT. 4866
    !C           JET
    !C
    !C
    !C DATE    : 25/07/91
    !C
    !C
    !C CHANGES : 31/07/91, H.P.SUMMERS - PROVIDE HYDROGEN ISOTOPE MASS USED
    !C                                   IN ADAS MASTER CX FILE PRODUCTION
    !C                                   IN INTEGER VARIABLE 'MH'.
    !C                                   CHANGE CX INTERPOLATION TO OBTAIN
    !C                                   COLL-RAD DEPENDENCE AT TE/NE AND
    !C                                   SCALE TO TI AT NE=0
    !C
    !C
    !C-----------------------------------------------------------------------
    !C
    !C
    !C
    integer, parameter :: mhyd=3
    integer*4, intent(in) :: iz0,iout
    real*4, intent(in) :: te,ti,dens,densh,fracth(mhyd)
    real*4, intent(out) :: sa0,sa(izdim),rta(izdim),pta0,pta(izdim)
    integer*4, intent(out) :: ifail
    !C
    !C-----------------------------------------------------------------------
    !C
    integer :: ITE(IPDIM) , INE(IPDIM)
    real*4 :: ZINT(IPDIM), ZINT1(IPDIM)
    !C
    real*4 :: CXC(IZDIM), PROCXC(IZDIM)
    real*4 :: CA(IZDIM) , PLT(IZDIM), PRB(IZDIM), PLS(IZDIM)
    real*4 :: PRC(IZDIM), RA(IZDIM)
    character(len=8) :: string
    integer :: is,i,jh,iz,iz01,k,iti2,iti6
    !C
    !C-----------------------------------------------------------------------
    !C
    real*4 :: HMASAU(mhyd) = (/ 1.007825D+00 , 2.0140D+00 , 3.01605D+00 /)
    !C
    !C-----CHECK THAT REQUIRED ELEMENT HAS BEEN READ-------------------------
    !C
    IS        = 0
    DO I    = 1 , ISMAX
       IF( IZ0A(I).EQ.IZ0 ) IS = I
    enddo
    !C
    IF( IS.LE.0 ) THEN
       WRITE(IOUT,9000) IZ0
       IFAIL = 1
       return
    END IF
    !C
    !C-----IDENTIFY HYDROGEN ISOTOPE USED IN ADAS CX MASTER FILE PRODUCTION--
    !C
    STRING   = ZLINFO(1,2,IS)
    READ(STRING,'(4X,F4.2)')HMADAS
    !C
    !C-----MEAN HYDROGEN ISOTOPE MASS FOR CHARGE EXCHANGE--------------------
    !C
    FSUM      = 0.0D+00
    DO JH   = 1 , MHYD
       FSUM    = FSUM + FRACTH(JH)
    enddo
    !C
    IF( FSUM.LE.0.0D+00 ) THEN
       WRITE(IOUT,9010) FSUM
       IFAIL = 2
       return
    END IF
    !C
    ZAVMAS    = 0.0D+00
    DO JH   = 1 , MHYD
       ZAVMAS = ZAVMAS + FRACTH(JH)*HMASAU(JH)/FSUM
    enddo
    !C
    !C=======================================================================
    !C
    !C  TABLE INTERPOLATION SECTION
    !C
    !C       1. RA(IZ)         ---  RADIATIVE + DIELECTRONIC RECOMBINATION
    !C       2. CA(IZ)         ---  CHARGE EXCHANGE RECOMBINATION
    !C
    !C       3. SA0            ---  IONISATION (NEUTRAL)
    !C          SA(IZ)         ---  IONISATION (NON-NEUTRAL)
    !C
    !C       4. PLT0           ---  TOTAL LINE RADIATED POWER (NEUTRAL)
    !C          PLT(IZ)        ---  TOTAL LINE RADIATED POWER (NON-NEUTRAL)
    !C
    !C       5. PRB(IZ)        ---  RAD. + DIEL RECOM. + BREMS. POWER
    !C       6. PRC(IZ)        ---  CX. RECOM. POWER
    !C
    !C       7. PLS0           ---  SPECIFIC LINE POWER     (NEUTRAL)
    !C          PLS(IZ)        ---  SPECIFIC LINE POWER     (NON-NEUTRAL)
    !C
    !C
    !C       INDEXING:
    !C                  IZ =  STAGE    (1 .LE. IZ .LE. IZ0 )
    !C=======================================================================
    !C-----------------------------------------------------------------------
    !C     ZERO THE RATES
    !C-----------------------------------------------------------------------
    SA0   = 0.0E0
    PLT0  = 0.0E0
    PLS0  = 0.0E0
    DO IZ=1,IZ0
       RA(IZ) = 0.0E0
       CA(IZ) = 0.0E0
       SA(IZ) = 0.0E0
       PLT(IZ)= 0.0E0
       PRB(IZ)= 0.0E0
       PRC(IZ)= 0.0E0
       PLS(IZ)= 0.0E0
    enddo
    !C
    IZ01 = IZ0 - 1
    !C-----------------------------------------------------------------------
    !C     FIND LOOKUP TABLE POSITIONS FOR ELECTRON TEMP AND DENSITY
    !C-----------------------------------------------------------------------
    ZTELG = ALOG10(TE)
    CALL SRCHA_eirene(ZTELG,ZTEL(1,1,IS),ITDIM,ITMAX(1,IS),ITE)
    !C
    ZNELG = ALOG10(DENS)
    CALL SRCHA_eirene(ZNELG,ZNEL(1,1,IS),IDDIM,IDMAX(1,IS),INE)
    !C
    DO K=1,IPDIM
       ZINT(K) = (ZTEL(ITE(K),K,IS)-ZTELG) / (ZTEL(ITE(K),K,IS)-ZTEL(ITE(K)-1,K,IS))
       ZINT1(K) = (ZNEL(INE(K),K,IS)-ZNELG) / (ZNEL(INE(K),K,IS)-ZNEL(INE(K)-1,K,IS))
    enddo
    !C-------------------------------------
    !C  INTERPOLATE NEUTRAL STAGE TABLES
    !C-------------------------------------
    !C
    !C----------
    !C  2*LINEAR
    !C----------
    ZSL1 = ZSAL0(ITE(3)-1,INE(3)-1,IS)*ZINT(3) + ZSAL0(ITE(3),INE(3)-1,IS)*(1.0-ZINT(3))
    ZSL2 = ZSAL0(ITE(3)-1,INE(3),IS)*ZINT(3) + ZSAL0(ITE(3),INE(3),IS)*(1.0-ZINT(3))
    ZSLG = ZSL1 * ZINT1(3) + ZSL2 * (1.0 - ZINT1(3))
    SA0  = 10.0**ZSLG
    !C
    ZPLT1 = ZPLTL0(ITE(4)-1,INE(4)-1,IS)*ZINT(4) + ZPLTL0(ITE(4),INE(4)-1,IS)*(1.0-ZINT(4))
    ZPLT2 = ZPLTL0(ITE(4)-1,INE(4),IS)*ZINT(4) + ZPLTL0(ITE(4),INE(4),IS)*(1.0-ZINT(4))
    ZPLTLG = ZPLT1 * ZINT1(4) + ZPLT2 * (1.0 - ZINT1(4))
    PLT0   = 10.0**ZPLTLG
    !C
    ZPLS1 = ZPLSL0(ITE(7)-1,INE(7)-1,IS)*ZINT(7) + ZPLSL0(ITE(7),INE(7)-1,IS)*(1.0-ZINT(7))
    ZPLS2 = ZPLSL0(ITE(7)-1,INE(7),IS)*ZINT(7) + ZPLSL0(ITE(7),INE(7),IS)*(1.0-ZINT(7))
    ZPLSLG = ZPLS1 * ZINT1(7) + ZPLS2 * (1.0 - ZINT1(7))
    PLS0   = 10.0**ZPLSLG
    !C-------------------------------------
    !C  INTERPOLATE ION  STAGE TABLES
    !C-------------------------------------
    !C----------
    !C  2*LINEAR
    !C----------
    DO IZ=1,IZ01
       ZSL1 = ZSAL(ITE(3)-1,INE(3)-1,IZ,IS)*ZINT(3) + ZSAL(ITE(3),INE(3)-1,IZ,IS)*(1.0-ZINT(3))
       ZSL2 = ZSAL(ITE(3)-1,INE(3),IZ,IS)*ZINT(3) + ZSAL(ITE(3),INE(3),IZ,IS)*(1.0-ZINT(3))
       ZSLG = ZSL1 * ZINT1(3) + ZSL2 * (1.0 - ZINT1(3))
       SA(IZ) = 10.0**ZSLG
       !C
       ZPLT1 = ZPLTL(ITE(4)-1,INE(4)-1,IZ,IS)*ZINT(4) + ZPLTL(ITE(4),INE(4)-1,IZ,IS)*(1.0-ZINT(4))
       ZPLT2 = ZPLTL(ITE(4)-1,INE(4),IZ,IS)*ZINT(4) + ZPLTL(ITE(4),INE(4),IZ,IS)*(1.0-ZINT(4))
       ZPLTLG = ZPLT1 * ZINT1(4) + ZPLT2 * (1.0 - ZINT1(4))
       PLT(IZ) = 10.0**ZPLTLG
       !C
       ZPLS1 = ZPLSL(ITE(7)-1,INE(7)-1,IZ,IS)*ZINT(7) + ZPLSL(ITE(7),INE(7)-1,IZ,IS)*(1.0-ZINT(7))
       ZPLS2 = ZPLSL(ITE(7)-1,INE(7),IZ,IS)*ZINT(7) + ZPLSL(ITE(7),INE(7),IZ,IS)*(1.0-ZINT(7))
       ZPLSLG = ZPLS1 * ZINT1(7) + ZPLS2 * (1.0 - ZINT1(7))
       PLS(IZ) = 10.0**ZPLSLG
!C
       ZRA1 = ZRAL(ITE(1)-1,INE(1)-1,IZ,IS)*ZINT(1) + ZRAL(ITE(1),INE(1)-1,IZ,IS)*(1.0-ZINT(1))
       ZRA2 = ZRAL(ITE(1)-1,INE(1),IZ,IS)*ZINT(1) + ZRAL(ITE(1),INE(1),IZ,IS)*(1.0-ZINT(1))
       ZRLG = ZRA1 * ZINT1(1) + ZRA2 * (1.0 - ZINT1(1))
       RA(IZ) = 10.0**ZRLG
!C
       ZPA1 = ZPRBL(ITE(5)-1,INE(5)-1,IZ,IS)*ZINT(5) + ZPRBL(ITE(5),INE(5)-1,IZ,IS)*(1.0-ZINT(5))
       ZPA2 = ZPRBL(ITE(5)-1,INE(5),IZ,IS)*ZINT(5) + ZPRBL(ITE(5),INE(5),IZ,IS)*(1.0-ZINT(5))
       ZPRBLG = ZPA1 * ZINT1(5) + ZPA2 * (1.0 - ZINT1(5))
       PRB(IZ) = 10.0**ZPRBLG
       !C--------------------------------------------------------------
    enddo
    !C-----------------------------------------------------------------------
    !C  INTERPOLATE LAST STAGE TABLES
    !C-----------------------------------------------------------------------
    SA(IZ0) = 0.0E0
    PLT(IZ0) = 0.0E0
    PLS(IZ0) = 0.0E0
    !C
    !C----------
    !C  2*LINEAR
    !C----------
    ZRA1 = ZRAL(ITE(1)-1, INE(1)-1, IZ0, IS) * ZINT(1) + ZRAL(ITE(1), INE(1)-1, IZ0, IS) * (1.0-ZINT(1))
    ZRA2 = ZRAL(ITE(1)-1, INE(1), IZ0, IS) * ZINT(1) + ZRAL(ITE(1), INE(1), IZ0, IS) * (1.0-ZINT(1))
    ZRLG = ZRA1 * ZINT1(1) + ZRA2 * (1.0 - ZINT1(1))
    RA(IZ0) = 10.0**ZRLG
    !C
    ZPA1 = ZPRBL(ITE(5)-1, INE(5)-1, IZ0, IS) * ZINT(5) + ZPRBL(ITE(5), INE(5)-1, IZ0, IS) * (1.0-ZINT(5))
    ZPA2 = ZPRBL(ITE(5)-1, INE(5), IZ0, IS) * ZINT(5) + ZPRBL(ITE(5), INE(5), IZ0, IS) * (1.0-ZINT(5))
    ZPRBLG = ZPA1 * ZINT1(5) + ZPA2 * (1.0 - ZINT1(5))
    PRB(IZ0) = 10.0**ZPRBLG

    !C-----------------------------------------------------------------------
    !C  INTERPOLATE CHARGE EXCHANGE TABLE
    !C-----------------------------------------------------------------------
    IF ( DENSH .GT. 0.0000) THEN
       !C--------------------------------------------------
       !C  SHIFT TI ACCORDING TO MEAN HYDROGEN ISOTOPE MASS
       !C---------------------------------------------------
       TISH  = TI*HMADAS/ZAVMAS
       !C---------------------------------------------------------
       !C  FIND LOOKUP TABLE POSITION FOR ION/NEUTRAL TEMPERATURES
       !C---------------------------------------------------------
       ZTILG = ALOG10(TISH)
       CALL SRCHB_eirene(ZTILG,ZTEL(1,2,IS),ITMAX(2,IS),ITI2)
       CALL SRCHB_eirene(ZTILG,ZTEL(1,6,IS),ITMAX(6,IS),ITI6)
       !C
       ZINTI2=(ZTEL(ITI2,2,IS)-ZTILG) / (ZTEL(ITI2,2,IS)-ZTEL(ITI2-1,2,IS))
       ZINTI6=(ZTEL(ITI6,6,IS)-ZTILG) / (ZTEL(ITI6,6,IS)-ZTEL(ITI6-1,6,IS))
       !C----------------------------
       !C  1*LINEAR AT TI AND TE AT LOWEST DENSITY
       !C----------------------------
       DO IZ = 1, IZ0
          ZCAI0= ZCAL(ITI2-1,1,IZ,IS)*ZINTI2 + ZCAL(ITI2,1,IZ,IS)*(1.0-ZINTI2)
          ZCAE0= ZCAL(ITE(2)-1,1,IZ,IS)*ZINT(2) + ZCAL(ITE(2),1,IZ,IS)*(1.0-ZINT(2))
          !C
          ZPAI0= ZPRCL(ITI6-1,1,IZ,IS)*ZINTI6 + ZPRCL(ITI6,1,IZ,IS)*(1.0-ZINTI6)
          ZPAE0= ZPRCL(ITE(6)-1,1,IZ,IS)*ZINT(2) + ZPRCL(ITE(6),1,IZ,IS)*(1.0-ZINT(6))
          !C---------------------
          !C  2*LINEAR AT TE/DENS
          !C---------------------
          ZCA1 = ZCAL(ITE(2)-1,INE(2)-1,IZ,IS)*ZINT(2) + ZCAL(ITE(2),INE(2)-1,IZ,IS)*(1.0-ZINT(2))
          ZCA2 = ZCAL(ITE(2)-1,INE(2),IZ,IS)*ZINT(2) + ZCAL(ITE(2),INE(2),IZ,IS)*(1.0-ZINT(2))
          ZCLG = ZCA1 * ZINT1(2) + ZCA2 * (1.0 - ZINT1(2))
          !C---------------------------------------------------------------------
          !C  FORM SCALED RESULT & MULTIPLY BY DENSH/DENS FOR COMBINATION WITH RA
          !C---------------------------------------------------------------------
          ZCLG=10.0**(ZCLG+ZCAI0-ZCAE0)
          CA(IZ)=ZCLG*DENSH/DENS
          !C
          ZPA1 = ZPRCL(ITE(6)-1,INE(6)-1,IZ,IS)*ZINT(6) + ZPRCL(ITE(6),INE(6)-1,IZ,IS)*(1.0-ZINT(6))
          ZPA2 = ZPRCL(ITE(6)-1,INE(6),IZ,IS)*ZINT(6) + ZPRCL(ITE(6),INE(6),IZ,IS)*(1.0-ZINT(6))
          ZPRCLG = ZPA1 * ZINT1(6) + ZPA2 * (1.0 - ZINT1(6))
          !C---------------------------------------------------------------------
          !C  FORM SCALED RESULT & MULTIPLY BY DENSH/DENS FOR COMBINATION WITH PRB
          !C---------------------------------------------------------------------
          ZPRCLG=10.0**(ZPRCLG+ZPAI0-ZPAE0)
          PRC(IZ)=ZPRCLG*DENSH/DENS
       enddo
    ENDIF
    !C------------------------------------
    !C    COMBINE RA & CA , PLT, PRB & PRC
    !C------------------------------------
    PTA0=PLT0
    IZ=0
    DO IZ=1,IZ0
       RTA(IZ)=RA(IZ)+CA(IZ)
       PTA(IZ)=PLT(IZ)+PRB(IZ)+PRC(IZ)
    enddo
    !C
    !C-----------------------------------------------------------------------
    !C
    RETURN
    !C
    !C-----------------------------------------------------------------------
    !C
9000 FORMAT( ' *** ERROR *** ELEMENT IZ0=',I2,' HAS NOT BEEN EXTRACTED' / '               FROM THE ADAS ISONUCLEAR MASTER FILES.'  )
9010 FORMAT( ' *** ERROR *** HYDROGEN ISOTOPIC FRACTIONS SUM TO ', 1PE12.4 )
    !C
    !C-----------------------------------------------------------------------
    !C
  END SUBROUTINE ADASIT_EIRENE

  !C=======================================================================
  SUBROUTINE ADASITX_eirene( IZ0   , IOUT , TE    , TI   , DENS , DENSH , FRACTH , SA0   , SA   , RTA  , PTA0  , PTA  , IFAIL )
    IMPLICIT REAL*4(A-H,O-Z)
    !C
    !C-----------------------------------------------------------------------
    !C
    !C
    !C PURPOSE : TO INTERPOLATE 'ADAS' ISONUCLEAR MASTER FILES AT THE
    !C           REQUIRED ELECTRON TEMPERATURE, ION TEMPERATURE & ELECTRON
    !C           DENSITY.  THE REQUIRED ELEMENT MUST ALREADY HAVE BEEN
    !C           EXTRACTED FROM THE 'ADAS' ISONUCLEAR MASTER FILES.
    !C
    !C
    !C TYPE    : SINGLE PRECISION
    !C
    !C
    !C INPUT   : (I*4) IZ0          = NUCLEAR CHARGE OF REQUIRED ELEMENT
    !C           (I*4) IOUT         = OUTPUT CHANNEL FOR MESSAGES/ERRORS
    !C           (R*4) TE           = ELECTRON TEMPERATURE (EV) AT WHICH THE
    !C                                ADAS COEFFICIENTS ARE REQUIRED
    !C           (R*4) TI           = ION TEMPERATURE (EV) AT WHICH THE
    !C                                ADAS COEFFICIENTS ARE REQUIRED
    !C           (R*4) DENS         = ELECTRON DENSITY (CM-3) AT WHICH THE
    !C                                ADAS COEFFICIENTS ARE REQUIRED
    !C           (R*4) DENSH        = HYDROGEN DENSITY (CM-3) AT WHICH THE
    !C                                ADAS COEFFICIENTS ARE REQUIRED
    !C           (R*4) FRACTH()     = FRACTH(1) : HYDROGEN  FUEL FRACTION
    !C                                FRACTH(2) : DEUTERIUM FUEL FRACTION
    !C                                FRACTH(3) : TRITIUM   FUEL FRACTION
    !C
    !C
    !C OUTPUT  : (R*4) SA0          = NEUTRAL IONIS. RATE COEFFT.(CM3 S-1)
    !C           (R*4) SA()         = ION IONIS. RATE COEFFT.    (CM3 S-1)
    !C           (R*4) RTA()        = ION RECOMB. RATE COEFFT.   (CM3 S-1)
    !C           (R*4) PTA0         = NEUTRAL RADIATED POWER  (WATTS CM3 S-1)
    !C           (R*4) PTA()        = ION RADIATED POWER      (WATTS CM3 S-1)
    !C           (I*4) IFAIL        = 0 --- ROUTINE SUCCESSFUL
    !C                              = 1 --- THE REQUIRED ELEMENT (IZ0) HAS
    !C                                      NOT BEEN EXTRACTED FROM THE ADAS
    !C                                      'ISONUCLEAR MASTER FILES'
    !C                              = 2 --- FRACTH(1)+...+FRACTH(MHYD) IS
    !C                                      NOT A REAL NUMBER.
    !C
    !C
    !C COMMON  : (I*4) ITDIM        = MAX. NUMBER OF TEMPERATURES   ALLOWED
    !C           (I*4) IDDIM        = MAX. NUMBER OF DENSITIES      ALLOWED
    !C           (I*4) IZDIM        = MAX. NUMBER OF CHARGED STATES ALLOWED
    !C           (I*4) ISDIM        = MAX. NUMBER OF SPECIES        ALLOWED
    !C           (I*4) IPDIM        = MAX. NUMBER OF ADAS FILES     ALLOWED
    !C
    !C
    !C PROGRAM : (I*4) MHYD         = NUMBER OF HYDROGEN ISOTOPES
    !C           (I*4) IS           = INDEX OF ELEMENT STORED IN /CADAS2/
    !C           (R*4) FSUM         = SUM OF HYDROGEN ISOTOPIC FRACTIONS
    !C           (R*4) HMADAS       = HYDROGEN ISOTOPE USED IN MASTER CX FILE
    !C                                (1=PROTIUM, 2=DEUTERIUM, 3=TRITIUM)
    !C           (C*8) STRING       = WORK STRING
    !C
    !C
    !C NOTES   : (1) THE HYDROGEN ISOTOPIC FRACTIONS ARE NORMALISED FOR
    !C               THE CODE IN CASE THEY ARE NOT ON INPUT.
    !C           (2) THE ARRAY STACKING IS UNUSUAL!
    !C               FOR RA, CA, PRB, PRC DATA IS STORED AS FOLLOWS:
    !C                RECOMBINING ION CHARGE   RECOMBINED ION CHARGE  INDEX
    !C                ----------------------   ---------------------  -----
    !C                     1                       0 (NEUTRAL)          1
    !C                     IZ                      IZ-1                 IZ
    !C                     IZ0 (BARE NUCLEUS)      IZ0-1                IZ0
    !C
    !C               FOR SA, PLT, PLS DATA IS STORED AS FOLLOWS:
    !C                IONISING  ION  CHARGE    IONISED  ION  CHARGE   INDEX
    !C                ----------------------   ---------------------  -----
    !C                     0 (NEUTRAL)             1                 SEPARATE
    !C                     IZ                      IZ+1                 IZ
    !C                     IZ0-1                   IZ0 (BARE NUCLEUS)   IZ0-1
    !C
    !C
    !C               THIS IS DIFFERENT FROM THE CONVENTIONS ACCORDING TO
    !C               ADF12 IN THE ADAS MASTER FILES.
    !C
    !C AUTHOR  : JAMES SPENCE  (K1/0/80)  EXT. 4866
    !C           JET
    !C
    !C
    !C DATE    : 25/07/91
    !C
    !C
    !C CHANGES : 31/07/91, H.P.SUMMERS - PROVIDE HYDROGEN ISOTOPE MASS USED
    !C                                   IN ADAS MASTER CX FILE PRODUCTION
    !C                                   IN INTEGER VARIABLE 'MH'.
    !C                                   CHANGE CX INTERPOLATION TO OBTAIN
    !C                                   COLL-RAD DEPENDENCE AT TE/NE AND
    !C                                   SCALE TO TI AT NE=0
    !C            5/08/91, H.P.SUMMERS - ADD IONISATION  POTENTIAL ENERGY
    !C                                   PART TO THE POWER TO MATCH SIMONINI
    !C                                   USE - REMOVE LATER AS NOT PHYSICAL !
    !C
    !C
    !C-----------------------------------------------------------------------
    !C
    !C
    integer, parameter :: mhyd=3
    real*4, intent(in) :: FRACTH(MHYD),te,ti,dens,densh
    integer*4, intent(in) :: iz0,iout
    real*4, intent(out) ::  SA(IZDIM)    , RTA(IZDIM)   , PTA(IZDIM), sa0,pta0
    integer, intent(out) :: ifail
    !C
    !C-----------------------------------------------------------------------
    !C
    integer :: ITE(IPDIM) , INE(IPDIM)
    real*4 ::  ZINT(IPDIM), ZINT1(IPDIM)
    !C
    real*4 ::  CXC(IZDIM), PROCXC(IZDIM)
    real*4 ::  CA(IZDIM) , PLT(IZDIM), PRB(IZDIM), PLS(IZDIM)
    real*4 ::  PRC(IZDIM), RA(IZDIM) , PIA(IZDIM)
    character(len=8) :: string
    integer :: is,i,jh,iz,iz01,k,iti2,iti6,ixp
    !C
    !C-----------------------------------------------------------------------
    !C
    real*4 :: HMASAU(MHYD) =(/ 1.007825D+00 , 2.0140D+00 , 3.01605D+00 /)
    !C
    !C-----------------------------------------------------------------------
    !C
    real*4 :: BWNOA(55) = (/   109679.0 , &
         198311.0 ,   438909.0 , &
         43487.0 ,   610079.0 ,   987661.0 , &
         75192.0 ,   146883.0 ,  1241242.0 ,  1756019.0 , &
         66928.0 ,   202887.0 ,   305931.0 ,  2092001.0 , &
         2744108.0 , &
         90820.0 ,   196645.0 ,   386241.0 ,   520178.0 ,&
         3162395.0 ,  3952062.0 , &
         117225.0 ,   238751.0 ,   382704.0 ,   624866.0 ,&
         789537.0 ,  4452745.0 ,  5380089.0 , &
         109837.0 ,   283240.0 ,   443086.0 ,   624384.0 ,&
         918657.0 ,  1114008.0 ,  5963114.0 ,  7026394.0 ,&
         140525.0 ,   282059.0 ,   505777.0 ,   702830.0 ,&
         921430.0 ,  1267622.0 ,  1493629.0 ,  7693780.0 ,&
         8897240.0 ,&
         173932.0 ,   330391.0 ,   511800.0 ,   783300.0 ,&
         1018000.0 ,  1273800.0 ,  1671792.0 ,  1928462.0 , &
         9644960.0 , 10986873.0 /)
    !C
    !C-----CHECK THAT REQUIRED ELEMENT HAS BEEN READ-------------------------
    !C
    IS        = 0
    DO I    = 1 , ISMAX
       IF( IZ0A(I).EQ.IZ0 ) IS = I
    enddo
    !C
    IF( IS.LE.0 ) THEN
       WRITE(IOUT,9000) IZ0
       IFAIL = 1
       return
    END IF
    !C
    !C-----IDENTIFY HYDROGEN ISOTOPE USED IN ADAS CX MASTER FILE PRODUCTION--
    !C
    STRING   = ZLINFO(1,2,IS)
    READ(STRING,'(4X,F4.2)')HMADAS
    !C
    !C-----MEAN HYDROGEN ISOTOPE MASS FOR CHARGE EXCHANGE--------------------
    !C
    FSUM      = 0.0D+00
    DO JH   = 1 , MHYD
       FSUM    = FSUM + FRACTH(JH)
    enddo
    !C
    IF( FSUM.LE.0.0D+00 ) THEN
       WRITE(IOUT,9010) FSUM
       IFAIL = 2
       return
    END IF
    !C
    ZAVMAS    = 0.0D+00
    DO JH   = 1 , MHYD
       ZAVMAS = ZAVMAS + FRACTH(JH)*HMASAU(JH)/FSUM
    enddo
    !C
    !C=======================================================================
    !C
    !C  TABLE INTERPOLATION SECTION
    !C
    !C       1. RA(IZ)         ---  RADIATIVE + DIELECTRONIC RECOMBINATION
    !C       2. CA(IZ)         ---  CHARGE EXCHANGE RECOMBINATION
    !C
    !C       3. SA0            ---  IONISATION (NEUTRAL)
    !C          SA(IZ)         ---  IONISATION (NON-NEUTRAL)
    !C
    !C       4. PLT0           ---  TOTAL LINE RADIATED POWER (NEUTRAL)
    !C          PLT(IZ)        ---  TOTAL LINE RADIATED POWER (NON-NEUTRAL)
    !C
    !C       5. PRB(IZ)        ---  RAD. + DIEL RECOM. + BREMS. POWER
    !C       6. PRC(IZ)        ---  CX. RECOM. POWER
    !C
    !C       7. PLS0           ---  SPECIFIC LINE POWER     (NEUTRAL)
    !C          PLS(IZ)        ---  SPECIFIC LINE POWER     (NON-NEUTRAL)
    !C
    !C
    !C       INDEXING:
    !C                  IZ =  STAGE    (1 .LE. IZ .LE. IZ0 )
    !C=======================================================================
    !C-----------------------------------------------------------------------
    !C     ZERO THE RATES
    !C-----------------------------------------------------------------------
    SA0   = 0.0E0
    PLT0  = 0.0E0
    PLS0  = 0.0E0
    DO IZ=1,IZ0
       RA(IZ) = 0.0E0
       CA(IZ) = 0.0E0
       SA(IZ) = 0.0E0
       PLT(IZ)= 0.0E0
       PRB(IZ)= 0.0E0
       PRC(IZ)= 0.0E0
       PLS(IZ)= 0.0E0
    enddo
    !C
    IZ01 = IZ0 - 1
    !C-----------------------------------------------------------------------
    !C     FIND LOOKUP TABLE POSITIONS FOR ELECTRON TEMP AND DENSITY
    !C-----------------------------------------------------------------------
    ZTELG = ALOG10(TE)
    CALL SRCHA_eirene(ZTELG,ZTEL(1,1,IS),ITDIM,ITMAX(1,IS),ITE)
    !C
    ZNELG = ALOG10(DENS)
    CALL SRCHA_eirene(ZNELG,ZNEL(1,1,IS),IDDIM,IDMAX(1,IS),INE)
    !C
    DO K=1,IPDIM
       ZINT(K) = (ZTEL(ITE(K),K,IS)-ZTELG) / (ZTEL(ITE(K),K,IS)-ZTEL(ITE(K)-1,K,IS))
       ZINT1(K) = (ZNEL(INE(K),K,IS)-ZNELG) / (ZNEL(INE(K),K,IS)-ZNEL(INE(K)-1,K,IS))
    enddo
    !C-------------------------------------
    !C  INTERPOLATE NEUTRAL STAGE TABLES
    !C-------------------------------------
    !C
    !C----------
    !C  2*LINEAR
    !C----------
    ZSL1 = ZSAL0(ITE(3)-1,INE(3)-1,IS)*ZINT(3) + ZSAL0(ITE(3),INE(3)-1,IS)*(1.0-ZINT(3))
    ZSL2 = ZSAL0(ITE(3)-1,INE(3),IS)*ZINT(3) + ZSAL0(ITE(3),INE(3),IS)*(1.0-ZINT(3))
    ZSLG = ZSL1 * ZINT1(3) + ZSL2 * (1.0 - ZINT1(3))
    SA0  = 10.0**ZSLG
    !C
    ZPLT1 = ZPLTL0(ITE(4)-1,INE(4)-1,IS)*ZINT(4) + ZPLTL0(ITE(4),INE(4)-1,IS)*(1.0-ZINT(4))
    ZPLT2 = ZPLTL0(ITE(4)-1,INE(4),IS)*ZINT(4) + ZPLTL0(ITE(4),INE(4),IS)*(1.0-ZINT(4))
    ZPLTLG = ZPLT1 * ZINT1(4) + ZPLT2 * (1.0 - ZINT1(4))
    PLT0   = 10.0**ZPLTLG
    !C
    ZPLS1 = ZPLSL0(ITE(7)-1,INE(7)-1,IS)*ZINT(7) + ZPLSL0(ITE(7),INE(7)-1,IS)*(1.0-ZINT(7))
    ZPLS2 = ZPLSL0(ITE(7)-1,INE(7),IS)*ZINT(7) + ZPLSL0(ITE(7),INE(7),IS)*(1.0-ZINT(7))
    ZPLSLG = ZPLS1 * ZINT1(7) + ZPLS2 * (1.0 - ZINT1(7))
    PLS0   = 10.0**ZPLSLG
    !C-------------------------------------
    !C  INTERPOLATE ION  STAGE TABLES
    !C-------------------------------------
    !C----------
    !C  2*LINEAR
    !C----------
    DO IZ=1,IZ01
       ZSL1 = ZSAL(ITE(3)-1,INE(3)-1,IZ,IS)*ZINT(3) + ZSAL(ITE(3),INE(3)-1,IZ,IS)*(1.0-ZINT(3))
       ZSL2 = ZSAL(ITE(3)-1,INE(3),IZ,IS)*ZINT(3) + ZSAL(ITE(3),INE(3),IZ,IS)*(1.0-ZINT(3))
       ZSLG = ZSL1 * ZINT1(3) + ZSL2 * (1.0 - ZINT1(3))
       SA(IZ) = 10.0**ZSLG
       !C
       ZPLT1 = ZPLTL(ITE(4)-1,INE(4)-1,IZ,IS)*ZINT(4) + ZPLTL(ITE(4),INE(4)-1,IZ,IS)*(1.0-ZINT(4))
       ZPLT2 = ZPLTL(ITE(4)-1,INE(4),IZ,IS)*ZINT(4) + ZPLTL(ITE(4),INE(4),IZ,IS)*(1.0-ZINT(4))
       ZPLTLG = ZPLT1 * ZINT1(4) + ZPLT2 * (1.0 - ZINT1(4))
       PLT(IZ) = 10.0**ZPLTLG
       !C
       ZPLS1 = ZPLSL(ITE(7)-1,INE(7)-1,IZ,IS)*ZINT(7) + ZPLSL(ITE(7),INE(7)-1,IZ,IS)*(1.0-ZINT(7))
       ZPLS2 = ZPLSL(ITE(7)-1,INE(7),IZ,IS)*ZINT(7) + ZPLSL(ITE(7),INE(7),IZ,IS)*(1.0-ZINT(7))
       ZPLSLG = ZPLS1 * ZINT1(7) + ZPLS2 * (1.0 - ZINT1(7))
       PLS(IZ) = 10.0**ZPLSLG
       !C
       ZRA1 = ZRAL(ITE(1)-1,INE(1)-1,IZ,IS)*ZINT(1) + ZRAL(ITE(1),INE(1)-1,IZ,IS)*(1.0-ZINT(1))
       ZRA2 = ZRAL(ITE(1)-1,INE(1),IZ,IS)*ZINT(1) + ZRAL(ITE(1),INE(1),IZ,IS)*(1.0-ZINT(1))
       ZRLG = ZRA1 * ZINT1(1) + ZRA2 * (1.0 - ZINT1(1))
       RA(IZ) = 10.0**ZRLG
       !C
       ZPA1 = ZPRBL(ITE(5)-1,INE(5)-1,IZ,IS)*ZINT(5) + ZPRBL(ITE(5),INE(5)-1,IZ,IS)*(1.0-ZINT(5))
       ZPA2 = ZPRBL(ITE(5)-1,INE(5),IZ,IS)*ZINT(5) + ZPRBL(ITE(5),INE(5),IZ,IS)*(1.0-ZINT(5))
       ZPRBLG = ZPA1 * ZINT1(5) + ZPA2 * (1.0 - ZINT1(5))
       PRB(IZ) = 10.0**ZPRBLG
       !C--------------------------------------------------------------
    enddo
    !C-----------------------------------------------------------------------
    !C  INTERPOLATE LAST STAGE TABLES
    !C-----------------------------------------------------------------------
    SA(IZ0) = 0.0E0
    PLT(IZ0) = 0.0E0
    PLS(IZ0) = 0.0E0
    !C
    !C----------
    !C  2*LINEAR
    !C----------
    ZRA1 = ZRAL(ITE(1)-1, INE(1)-1, IZ0, IS) * ZINT(1) + ZRAL(ITE(1), INE(1)-1, IZ0, IS) * (1.0-ZINT(1))
    ZRA2 = ZRAL(ITE(1)-1, INE(1), IZ0, IS) * ZINT(1) + ZRAL(ITE(1), INE(1), IZ0, IS) * (1.0-ZINT(1))
    ZRLG = ZRA1 * ZINT1(1) + ZRA2 * (1.0 - ZINT1(1))
    RA(IZ0) = 10.0**ZRLG
    !C
    ZPA1 = ZPRBL(ITE(5)-1, INE(5)-1, IZ0, IS) * ZINT(5) + ZPRBL(ITE(5), INE(5)-1, IZ0, IS) * (1.0-ZINT(5))
    ZPA2 = ZPRBL(ITE(5)-1, INE(5), IZ0, IS) * ZINT(5) + ZPRBL(ITE(5), INE(5), IZ0, IS) * (1.0-ZINT(5))
    ZPRBLG = ZPA1 * ZINT1(5) + ZPA2 * (1.0 - ZINT1(5))
    PRB(IZ0) = 10.0**ZPRBLG
    !C-----------------------------------------------------------------------
    !C  INTERPOLATE CHARGE EXCHANGE TABLE
    !C-----------------------------------------------------------------------
    IF ( DENSH .GT. 0.0000) THEN
       !C--------------------------------------------------
       !C  SHIFT TI ACCORDING TO MEAN HYDROGEN ISOTOPE MASS
       !C---------------------------------------------------
       TISH  = TI*HMADAS/ZAVMAS
       !C---------------------------------------------------------
       !C  FIND LOOKUP TABLE POSITION FOR ION/NEUTRAL TEMPERATURES
       !C---------------------------------------------------------
       ZTILG = ALOG10(TISH)
       CALL SRCHB_eirene(ZTILG,ZTEL(1,2,IS),ITMAX(2,IS),ITI2)
       CALL SRCHB_eirene(ZTILG,ZTEL(1,6,IS),ITMAX(6,IS),ITI6)
       !C
       ZINTI2=(ZTEL(ITI2,2,IS)-ZTILG) / (ZTEL(ITI2,2,IS)-ZTEL(ITI2-1,2,IS))
       ZINTI6=(ZTEL(ITI6,6,IS)-ZTILG) / (ZTEL(ITI6,6,IS)-ZTEL(ITI6-1,6,IS))
       !C----------------------------
       !C  1*LINEAR AT TI AND TE AT LOWEST DENSITY
       !C----------------------------
       DO IZ = 1, IZ0
          ZCAI0= ZCAL(ITI2-1,1,IZ,IS)*ZINTI2 + ZCAL(ITI2,1,IZ,IS)*(1.0-ZINTI2)
          ZCAE0= ZCAL(ITE(2)-1,1,IZ,IS)*ZINT(2) + ZCAL(ITE(2),1,IZ,IS)*(1.0-ZINT(2))
          !C
          ZPAI0= ZPRCL(ITI6-1,1,IZ,IS)*ZINTI6 + ZPRCL(ITI6,1,IZ,IS)*(1.0-ZINTI6)
          ZPAE0= ZPRCL(ITE(6)-1,1,IZ,IS)*ZINT(2) + ZPRCL(ITE(6),1,IZ,IS)*(1.0-ZINT(6))
          !C---------------------
          !C  2*LINEAR AT TE/DENS
          !C---------------------
          ZCA1 = ZCAL(ITE(2)-1,INE(2)-1,IZ,IS)*ZINT(2) + ZCAL(ITE(2),INE(2)-1,IZ,IS)*(1.0-ZINT(2))
          ZCA2 = ZCAL(ITE(2)-1,INE(2),IZ,IS)*ZINT(2) + ZCAL(ITE(2),INE(2),IZ,IS)*(1.0-ZINT(2))
          ZCLG = ZCA1 * ZINT1(2) + ZCA2 * (1.0 - ZINT1(2))
          !C---------------------------------------------------------------------
          !C  FORM SCALED RESULT & MULTIPLY BY DENSH/DENS FOR COMBINATION WITH RA
          !C---------------------------------------------------------------------
          ZCLG=10.0**(ZCLG+ZCAI0-ZCAE0)
          CA(IZ)=ZCLG*DENSH/DENS
          !C
          ZPA1 = ZPRCL(ITE(6)-1,INE(6)-1,IZ,IS)*ZINT(6) + ZPRCL(ITE(6),INE(6)-1,IZ,IS)*(1.0-ZINT(6))
          ZPA2 = ZPRCL(ITE(6)-1,INE(6),IZ,IS)*ZINT(6) + ZPRCL(ITE(6),INE(6),IZ,IS)*(1.0-ZINT(6))
          ZPRCLG = ZPA1 * ZINT1(6) + ZPA2 * (1.0 - ZINT1(6))
          !C---------------------------------------------------------------------
          !C  FORM SCALED RESULT & MULTIPLY BY DENSH/DENS FOR COMBINATION WITH PRB
          !C---------------------------------------------------------------------
          ZPRCLG=10.0**(ZPRCLG+ZPAI0-ZPAE0)
          PRC(IZ)=ZPRCLG*DENSH/DENS
       enddo
    ENDIF
    !C------------------------------------
    !C    COMBINE RA & CA , PLT, PRB & PRC
    !C------------------------------------
    IXP=(IZ0*(IZ0-1))/2+1
    !C
    IF(IZ0.LE.10) THEN
       PIA0=BWNOA(IXP)*1.98618E-23*SA0
    ELSE
       PIA0=0.0
    ENDIF
    !C
    PTA0=PLT0+PIA0
    IZ=0
    DO IZ=1,IZ0
       IF(IZ.LT.IZ0.AND.IZ0.LE.10) THEN
          PIA(IZ)=BWNOA(IXP+IZ)*1.98618E-23*SA(IZ)
       ELSE
          PIA(IZ)=0.0
       ENDIF
       !C
       RTA(IZ)=RA(IZ)+CA(IZ)
       PTA(IZ)=PLT(IZ)+PRB(IZ)+PRC(IZ)+PIA(IZ)
    enddo
    !C
    !C-----------------------------------------------------------------------
    !C
    RETURN
    !C
    !C-----------------------------------------------------------------------
    !C
9000 FORMAT( ' *** ERROR *** ELEMENT IZ0=',I2,' HAS NOT BEEN EXTRACTED' / '               FROM THE ADAS ISONUCLEAR MASTER FILES.'  )
9010 FORMAT( ' *** ERROR *** HYDROGEN ISOTOPIC FRACTIONS SUM TO ' , 1PE12.4 )
    !C
    !C-----------------------------------------------------------------------
    !C
  END SUBROUTINE ADASITX_EIRENE

!!$  !C=======================================================================
!!$  SUBROUTINE ADASR1(  ICHAN , IPOINT , IZ0    , IS , ITDIM , IDDIM  , IZDIM  , ISDIM  , IPDIM , ITMAX , IDMAX  , IZ0MAX , TEMPE , DENSE  , ZDATA  , ZLINFO , IFAIL )
!!$    IMPLICIT REAL*4(A-H,O-Z)
!!$    !C
!!$    !C-----------------------------------------------------------------------
!!$    !C
!!$    !C PURPOSE : TO READ 'ADAS' ELEMENT MASTER FILES FORMATTED ACCORDING TO
!!$    !C           ADAS DATA FORMAT 'ADF01'.  NO SEPARATION OF THE NEUTRAL
!!$    !C           IONISATION STAGE DATA INTO SPECIAL ARRAYS IS DONE.
!!$    !C
!!$    !C           THE DATA RETURNED IS IN THE FORM :-
!!$    !C
!!$    !C                        ZDATA( IT , ID , IZ , IS )
!!$    !C           WITH
!!$    !C
!!$    !C                        TEMPE(IT,IPOINT,IS)
!!$    !C                        DENSE(IT,IPOINT,IS)
!!$    !C           WHERE,
!!$    !C                  IT       =  TEMPERATURE INDEX ( 1 - ITMAX  )
!!$    !C                  ID       =  DENSITY     INDEX ( 1 - IDMAX  )
!!$    !C                  IZ       =  ION. STAGE  INDEX ( 1 - IZ0MAX )
!!$    !C                  IS       =  SPECIES     INDEX ( 1 - ISDIM  )
!!$    !C                  IPOINT   =  DATA TYPE   INDEX ( 1 - IPDIM  )
!!$    !C                  TEMPE(,,)=  LOG10(ELEC. TEMPERATURE)(EV)
!!$    !C                  DENSE(,) =  LOG10(ELEC. DENSITY)    (CM**3)
!!$    !C
!!$    !C
!!$    !C TYPE    : SINGLE PRECISION
!!$    !C
!!$    !C INPUT   : (I*4)  ICHAN    =  STREAM NUMBER (PREVIOUSLY ALLOCATED)
!!$    !C           (I*4)  IPOINT   =  DATA TYPE POINTER FOR ARRAY STACKING
!!$    !C           (I*4)  ITDIM    =  ARRAY DIMENSION LIMIT FOR ITMAX
!!$    !C           (I*4)  IDDIM    =  ARRAY DIMENSION LIMIT FOR IDMAX
!!$    !C           (I*4)  IZDIM    =  ARRAY DIMENSION LIMIT FOR IZ0MAX
!!$    !C           (I*4)  ISDIM    =  ARRAY DIMENSION LIMIT FOR ISMAX
!!$    !C           (I*4)  IPDIM    =  ARRAY DIMENSION LIMIT FOR IPOINT
!!$    !C           (I*4)  IZ0      =  NUCLEAR CHARGE OF SPECIES
!!$    !C
!!$    !C OUTPUT  : (I*4)  ITMAX(,) =  NUMBER OF TEMPERATURES FOR IPOINT
!!$    !C           (I*4)  IDMAX(,) =  NUMBER OF DENSITIES FOR IPOINT
!!$    !C           (I*4)  IZ0MAX   =  NUMBER OF CHARGED STATES (MUST EQUAL IZ0)
!!$    !C           (R*4)  TEMPE(,,)=  LOG10(ELECTRON TEMPERATURES)
!!$    !C           (R*4)  DENSE(,,)=  LOG10(ELECTRON DENSITIES)
!!$    !C           (R*4)  ZDATA(,,,)= RETURNED COEFFICIENT ARRAY (SEE ABOVE)
!!$    !C                              (LOG10 VALUES)
!!$    !C           (C*8)  ZLINFO(,,)= TEXT INFORMATION FOR EACH ION AND DATA
!!$    !C                              TYPE.  USED FOR WAVELENGTH OF SPECIFIC
!!$    !C                              LINE AND HYDROGEN MASS FOR CX DATA
!!$    !C                              ONLY AT THIS STAGE
!!$    !C           (I*4)  IFAIL    =  0 --- IF ROUTINE SUCCESSFUL
!!$    !C                           =  1 --- MISMATCH OF NUCLEAR CHARGES
!!$    !C                           =  2 --- ITDIM SET TOO LOW
!!$    !C                           =  3 --- IDDIM SET TOO LOW
!!$    !C
!!$    !C NOTE    : THE ROUTINE ABORTS IF IZ0 DOES NOT MATCH THE NUCLEAR CHARGE,
!!$    !C           IZ0MAX, AND NUMBER OF STAGES,IZMAX, READ FROM THE ADAS FILE.
!!$    !C
!!$    !C AUTHOR  : JAMES SPENCE  (K1/0/80)  EXT. 4866
!!$    !C           JET
!!$    !C
!!$    !C DATE    : 23/03/90
!!$    !C
!!$    !C
!!$    !C CHANGES : 12/12/90, H.P.SUMMERS - PUT ADAS DATA FORM CODE IN POSITIONS
!!$    !C                                   76-80 OF FIRST LINE OF OUTPUT FILES.
!!$    !C                                   CODE IS 'ADF01'.
!!$    !C
!!$    !C           20/ 2/91, H.P.SUMMERS - ADD ZLINFO
!!$    !C           20/ 2/91, H.P.SUMMERS - ADD IPOINT
!!$    !C           23/07/91, J.SPENCE    - REMOVE WRITE OPTION & 'SANC0' STUFF.
!!$    !C           25/07/91, J.SPENCE    - ALLOW MULTIPLE SPECIES
!!$    !C
!!$    !C-----------------------------------------------------------------------
!!$    !C
!!$    DIMENSION ZDATA( ITDIM , IDDIM , IZDIM , ISDIM ) &
!!$         , TEMPE( ITDIM , IPDIM , ISDIM ) &
!!$         , DENSE( IDDIM , IPDIM , ISDIM ) &
!!$         , ITMAX( IPDIM , ISDIM ) &
!!$         , IDMAX( IPDIM , ISDIM ) &
!!$         , ZLINFO( IZDIM , IPDIM , ISDIM )
!!$    !C
!!$    !C-----------------------------------------------------------------------
!!$    !C
!!$    CHARACTER* :: ZLINFO
!!$    CHARACTER*80 ::  STRING
!!$    !C
!!$    !C-----------------------------------------------------------------------
!!$    !C
!!$    !C
!!$    IFAIL = 0
!!$    !C
!!$    !C-----------------------------------------------------------------------
!!$    !C
!!$    READ(ICHAN,1600) IZADAS , IDMAX(IPOINT,IS) , ITMAX(IPOINT,IS) , IDUMMY , IZ0MAX
!!$    !C
!!$    IF( ITDIM.LT.ITMAX(IPOINT,IS) ) THEN
!!$       IFAIL=2
!!$       RETURN
!!$    ENDIF
!!$    !C
!!$    IF( IDDIM.LT.IDMAX(IPOINT,IS) ) THEN
!!$       IFAIL=3
!!$       RETURN
!!$    ENDIF
!!$    !C
!!$    IF( IZ0.EQ.IZADAS.AND.IZ0.EQ.IZ0MAX ) THEN
!!$       !C
!!$       READ(ICHAN,1000) STRING
!!$       !C
!!$       READ(ICHAN,1200) ( DENSE(ID,IPOINT,IS),ID=1,IDMAX(IPOINT,IS))
!!$       !C
!!$       READ(ICHAN,1200) ( TEMPE(IT,IPOINT,IS),IT=1,ITMAX(IPOINT,IS))
!!$       !C
!!$       DO IZ = 1 , IZ0MAX
!!$          READ(ICHAN,1800)  ZLINFO(IZ,IPOINT,IS)
!!$          DO 300 IT = 1 , ITMAX(IPOINT,IS)
!!$             READ(ICHAN,1400)(ZDATA(IT,ID,IZ,IS),ID=1,IDMAX(IPOINT,IS))
!!$          enddo
!!$       enddo
!!$       !C
!!$    ELSE
!!$       !C
!!$       IFAIL = 1
!!$       !C
!!$    END IF
!!$    !C
!!$    RETURN
!!$    !C
!!$    !C-----------------------------------------------------------------------
!!$    !C
!!$1000 FORMAT( A )
!!$1200 FORMAT( 8(1X , F9.5) )
!!$1400 FORMAT( 8F10.5 )
!!$1600 FORMAT( 5I5 )
!!$1800 FORMAT( 44X , A )
!!$    !C
!!$    !C-----------------------------------------------------------------------
!!$    !C
!!$  END SUBROUTINE ADASR1
!!$
!!$  !C=======================================================================
!!$  SUBROUTINE ADASR2(  ICHAN , IPOINT , IZ0    , IS , ITDIM , IDDIM  , IZDIM  , ISDIM , IPDIM , ITMAX , IDMAX  , IZ0MAX , TEMPE , DENSE  , ZDATA0 , ZDATA , ZLINFO , IFAIL )
!!$    IMPLICIT REAL*4(A-H,O-Z)
!!$    !C
!!$    !C-----------------------------------------------------------------------
!!$    !C
!!$    !C PURPOSE : TO READ 'ADAS' ELEMENT MASTER FILES FORMATTED ACCORDING TO
!!$    !C           ADAS DATA FORMAT 'ADF01'.  SEPARATION OF THE NEUTRAL
!!$    !C           IONISATION STAGE DATA INTO SPECIAL ARRAYS IS DONE
!!$    !C
!!$    !C           THE DATA RETURNED IS IN THE FORM :-
!!$    !C
!!$    !C                        ZDATA0( IT , ID ,IS )
!!$    !C           AND
!!$    !C                        ZDATA( IT , ID , IZ-1 , IS )
!!$    !C           WITH
!!$    !C
!!$    !C                        TEMPE(IT,IPOINT,IS)
!!$    !C                        DENSE(IT,IPOINT,IS)
!!$    !C           WHERE,
!!$    !C                  IT       =  TEMPERATURE INDEX ( 1 - ITMAX  )
!!$    !C                  ID       =  DENSITY     INDEX ( 1 - IDMAX  )
!!$    !C                  IZ       =  ION. STAGE  INDEX ( 2 - IZ0MAX )
!!$    !C                  IS       =  SPECIES     INDEX ( 1 - ISDIM  )
!!$    !C                  IPOINT   =  DATA TYPE   INDEX ( 1 - IPDIM  )
!!$    !C                  TEMPE(,,)=  LOG10(ELEC. TEMPERATURE)(EV)
!!$    !C                  DENSE(,,)=  LOG10(ELEC. DENSITY)    (CM**3)
!!$    !C
!!$    !C
!!$    !C TYPE    : SINGLE PRECISION
!!$    !C
!!$    !C INPUT   : (I*4)  ICHAN    =  STREAM NUMBER (PREVIOUSLY ALLOCATED)
!!$    !C           (I*4)  IPOINT   =  DATA TYPE POINTER FOR ARRAY STACKING
!!$    !C           (I*4)  ITDIM    =  ARRAY DIMENSION LIMIT FOR ITMAX
!!$    !C           (I*4)  IDDIM    =  ARRAY DIMENSION LIMIT FOR IDMAX
!!$    !C           (I*4)  IZDIM    =  ARRAY DIMENSION LIMIT FOR IZ0MAX
!!$    !C           (I*4)  ISDIM    =  ARRAY DIMENSION LIMIT FOR ISMAX
!!$    !C           (I*4)  IPDIM    =  ARRAY DIMENSION LIMIT FOR IPOINT
!!$    !C           (I*4)  IZ0      =  NUCLEAR CHARGE OF SPECIES
!!$    !C
!!$    !C OUTPUT  : (I*4)  ITMAX(,) =  NUMBER OF TEMPERATURES FOR IPOINT
!!$    !C           (I*4)  IDMAX(,) =  NUMBER OF DENSITIES FOR IPOINT
!!$    !C           (I*4)  IZ9MAX   =  NUMBER OF CHARGE STATES (MUST EQUAL IZ0)
!!$    !C           (R*4)  TEMPE(,,)=  ELECTRON TEMPERATURES
!!$    !C           (R*4)  DENSE(,,)=  ELECTRON DENSITIES
!!$    !C           (R*4)  ZDATA0(,,)= RETURNED NEUTRAL COEFFICIENT ARRAY
!!$    !C           (R*4)  ZDATA(,,,)= RETURNED COEFFICIENT ARRAY (SEE ABOVE)
!!$    !C                              (LOG10 VALUES)
!!$    !C           (C*8)  ZLINFO(,,)= TEXT INFORMATION FOR EACH ION AND DATA
!!$    !C                              TYPE.  USED FOR WAVELENGTH OF SPECIFIC
!!$    !C                              LINE AND HYDROGEN MASS FOR CX DATA
!!$    !C                              ONLY AT THIS STAGE
!!$    !C           (I*4)  IFAIL    =  0 --- IF ROUTINE SUCCESSFUL
!!$    !C                           =  1 --- MISMATCH OF NUCLEAR CHARGES
!!$    !C                           =  2 --- ITDIM SET TOO LOW
!!$    !C                           =  3 --- IDDIM SET TOO LOW
!!$    !C
!!$    !C NOTE    : THE ROUTINE ABORTS IF IZ0 DOES NOT MATCH THE NUCLEAR CHARGE,
!!$    !C           IZ0MAX, AND NUMBER OF STAGES,IZMAX, READ FROM THE ADAS FILE.
!!$    !C
!!$    !C AUTHOR  : JAMES SPENCE  (K1/0/80)  EXT. 4866
!!$    !C           JET
!!$    !C
!!$    !C DATE    : 23/03/90
!!$    !C
!!$    !C
!!$    !C CHANGES : 12/12/90, H.P.SUMMERS - PUT ADAS DATA FROM CODE IN POSITIONS
!!$    !C                                   76-80 OF FIRST LINE OF OUTPUT FILES.
!!$    !C                                   CODE IS 'ADF01'.
!!$    !C           20/ 2/91, H.P.SUMMERS - ADD ZLINDO
!!$    !C           20/ 2/91, H.P.SUMMERS - ADD IPOINT
!!$    !C           23/07/91, J.SPENCE    - REMOVE WRITE OPTION & 'SANC0' STUFF.
!!$    !C           25/07/91, J.SPENCE    - ALLOW MULTIPLE SPECIES
!!$    !C
!!$    !C-----------------------------------------------------------------------
!!$    !C
!!$    DIMENSION ZDATA0( ITDIM , IDDIM , ISDIM )&
!!$         , ZDATA( ITDIM , IDDIM , IZDIM , ISDIM )&
!!$         , TEMPE( ITDIM , IPDIM , ISDIM )&
!!$         , DENSE( IDDIM , IPDIM , ISDIM )&
!!$         , ITMAX( IPDIM , ISDIM ) , IDMAX( IPDIM , ISDIM )&
!!$         , ZLINFO( IZDIM , IPDIM , ISDIM )
!!$    !C
!!$    !C-----------------------------------------------------------------------
!!$    !C
!!$    CHARACTER* ::  ZLINFO
!!$    CHARACTER*80 ::  STRING
!!$    !C
!!$    !C-----------------------------------------------------------------------
!!$    !C
!!$    !C
!!$    IFAIL = 0
!!$    !C
!!$    !C-----------------------------------------------------------------------
!!$    !C
!!$    READ(ICHAN,1600) IZADAS , IDMAX(IPOINT,IS) , ITMAX(IPOINT,IS) , IDUMMY , IZ0MAX
!!$    !C
!!$    IF( ITDIM.LT.ITMAX(IPOINT,IS) )THEN
!!$       IFAIL=2
!!$       RETURN
!!$    ENDIF
!!$    !C
!!$    IF( IDDIM.LT.IDMAX(IPOINT,IS) )THEN
!!$       IFAIL=3
!!$       RETURN
!!$    ENDIF
!!$    !C
!!$    IF( IZ0.EQ.IZADAS.AND.IZ0.EQ.IZ0MAX ) THEN
!!$       !C
!!$       READ(ICHAN,1000) STRING
!!$       !C
!!$       READ(ICHAN,1200) ( DENSE(ID,IPOINT,IS),ID=1,IDMAX(IPOINT,IS))
!!$       !C
!!$       READ(ICHAN,1200) ( TEMPE(IT,IPOINT,IS),IT=1,ITMAX(IPOINT,IS))
!!$       !C
!!$       DO IZ1 = 1 , IZ0MAX
!!$          IZ = IZ1 - 1
!!$          READ(ICHAN,1800)  ZLINFO(IZ1,IPOINT,IS)
!!$          DO IT = 1 , ITMAX(IPOINT,IS)
!!$             IF( IZ1.EQ.1 ) THEN
!!$                READ(ICHAN,1400)(ZDATA0(IT,ID,IS) , ID = 1 , IDMAX(IPOINT,IS))
!!$             ELSE
!!$                READ(ICHAN,1400)(ZDATA(IT,ID,IZ,IS) , ID = 1 , IDMAX(IPOINT,IS))
!!$             END IF
!!$          enddo
!!$       enddo
!!$       !C
!!$    ELSE
!!$       !C
!!$       IFAIL = 1
!!$       !C
!!$    END IF
!!$    !C
!!$    RETURN
!!$    !C
!!$    !C-----------------------------------------------------------------------
!!$    !C
!!$1000 FORMAT( A )
!!$1200 FORMAT( 8(1X,F9.5) )
!!$1400 FORMAT( 8F10.5 )
!!$1600 FORMAT( 5I5 )
!!$1800 FORMAT( 44X, A )
!!$    !C
!!$    !C-----------------------------------------------------------------------
!!$    !C
!!$  END SUBROUTINE ADASR2

  !C=======================================================================
  SUBROUTINE SRCHA_eirene(X,XA,NXDIM,NX,IX)
    implicit none
    !C-----------------------------------------------------------------------
    !C  ARRAY VERSION OF CODE SRCH
    !C
    !C  SEARCHES FOR THE POSITION OF A VALUE IN AN MONOTONIC INCREASING
    !C  VECTOR OF VALUES.
    !C  THE INDEX RETURNED IS THAT OF THE VECTOR VALUE AT THE TOP OF THE
    !C  RANGE WITHIN WHICH THE SEARCH VALUE LIES.
    !C
    !C  A BINARY SEARCH IS USED
    !C
    !C  INPUT
    !C      X     =  SEARCH VALUE TO BE POSITIONED
    !C      XA(,) =  MONOTONIC INCREASING VECTOR  OF VALUES TO BE SCANNED
    !C      NX() = NUMBER OF VALUES IN THE VECTOR
    !C
    !C  OUTPUT
    !C      IX() = POSITION MARKER OF SEARCH VALUE
    !C
    !C
    !C  ******  H.P.SUMMERS, JET           20 FEB 1991      *********
    !C-----------------------------------------------------------------------
    real*4, intent(in) :: x,xa(nxdim,7)
    integer, intent(in) :: nxdim,nx(7)
    integer, intent(out) :: ix(7)
    integer :: k,ir,iix
    DO  K=1,7
       IX(K)=(NX(K)+1)/2
       IR=(IX(K)+1)/2
       l1: do
          if(ix(k) .le. 0) then
             !write(*,*) ' mod_adas: warning, ix(k) <= 0 (1)'
             ix(k) = 0
             exit l1
          endif

          IF(X.EQ.XA(IX(K),K)) THEN
             exit l1
          ELSEIF (X.LT.XA(IX(K),K))THEN
             IX(K)=IX(K)-IR
          ELSE
             IX(K)=IX(K)+IR
          ENDIF

          IF(IR.GT.1)THEN
             IR=(IR+1)/2
             cycle l1
          endif

          if(ix(k) .le. 0)  then
             !write(*,*) ' mod_adas: warning, ix(k) <= 0 (2)'
             ix(k) = 0
             exit l1
          endif

          if(ix(k) .gt. nxdim)  then
             !write(*,*) ' mod_adas: warning, ix(k) > nxdim (3)'
             ix(k) = nxdim
             exit l1
          endif

          IF (X.LT.XA(IX(K),K))THEN
             IX(K)=IX(K)-IR
             exit l1
          ENDIF
       enddo l1

       IX(K)=IX(K)+1
       IX(K)=MAX0(2,IX(K))
       IX(K)=MIN0(IX(K),NX(K))
    enddo
    RETURN
  END SUBROUTINE SRCHA_EIRENE

  !C=======================================================================
  SUBROUTINE SRCHB_eirene(X,XA,NX,IX)
    implicit none
    !C-----------------------------------------------------------------------
    !C  SEARCHES FOR THE POSITION OF A VALUE IN AN MONOTONIC INCREASING
    !C  VECTOR OF VALUES.
    !C  THE INDEX RETURNED IS THAT OF THE VECTOR VALUE AT THE TOP OF THE
    !C  RANGE WITHIN WHICH THE SEARCH VALUE LIES.
    !C
    !C  A BINARY SEARCH IS USED
    !C
    !C  INPUT
    !C      X     =  SEARCH VALUE TO BE POSITIONED
    !C      XA()  =  MONOTONIC INCREASING VECTOR  OF VALUES TO BE SCANNED
    !C      NX    = NUMBER OF VALUES IN THE VECTOR
    !C
    !C  OUTPUT
    !C      IX    = POSITION MARKER OF SEARCH VALUE
    !C
    !C
    !C
    !C  AUTHOR:   H.P.SUMMERS, JET
    !C  DATE:     20 FEB 1991
    !C
    !C  UPDATE:   HPS 8 JULY 1993  - CHANGED SUBROUTINE NAME FROM SRCH TO
    !C                               SRCHB AT REQUEST OF JAMES SPENCE
    !C                               CALLS ALSO ALTERED IN ADASIT AND ADASITX
    !C-----------------------------------------------------------------------
    real*4, intent(in) :: x,xa(nx)
    integer, intent(in) :: nx
    integer, intent(out) :: ix
    integer :: ir
    IX=(NX+1)/2
    IR=(IX+1)/2
    do
       IF(X.EQ.XA(IX)) THEN
          exit
       ELSEIF (X.LT.XA(IX))THEN
          IX=IX-IR
       ELSE
          IX=IX+IR
       ENDIF

       IF(IR.GT.1)THEN
          IR=(IR+1)/2
          cycle
       ELSEIF (X.LT.XA(IX))THEN
          IX=IX-IR
       ENDIF
       exit
    enddo

    IX=IX+1
    IX=MAX0(2,IX)
    IX=MIN0(IX,NX)
    RETURN
  END SUBROUTINE SRCHB_EIRENE

  !CX  Port of JET3090 version to UNIX by L. Horton 3/8/95
  !CX
  SUBROUTINE D2DATA_eirene( YEAR   , YEARDF , TITLF  , IFAIL, &
       IZ0    , IZ1    , ICLASS , ITMAX  , IEVCUT,&
       ITDIMD , ITMAXD , IDMAXD , IZMAXD ,&
       DTEV   , DDENS ,&
       DTEVD  , DDENSD , DRCOFD , ZDATA ,&
       DRCOFI )
    IMPLICIT none
    !C
    !C-----------------------------------------------------------------------
    !C
    !C PURPOSE : TO EXTRACT 'SANC0' COLLISIONAL DIELECTRONIC DATA
    !C
    !C NOTE    : THE SOURCE DATA IS STORED AS FOLLOWS:
    !C
    !C   (1) $ADASUSER/<DEFADF>/acd<YR>/acd<YR>_<ELEMENT SYMBOL>.dat
    !C   (2) $ADASUSER/<DEFADF>/scd<YR>/scd<YR>_<ELEMENT SYMBOL>.dat
    !C   (3) $ADASUSER/<DEFADF>/ccd<YR>/ccd<YR>_<ELEMENT SYMBOL>.dat
    !C   (4) $ADASUSER/<DEFADF>/prb<YR>/prb<YR>_<ELEMENT SYMBOL>_ev<CUT>.dat
    !C   (5) $ADASUSER/<DEFADF>/plt<YR>/plt<YR>_<ELEMENT SYMBOL>_ev<CUT>.dat
    !C   (6) $ADASUSER/<DEFADF>/prc<YR>/prc<YR>_<ELEMENT SYMBOL>_ev<CUT>.dat
    !C   (7) $ADASUSER/<DEFADF>/pls<YR>/pls<YR>_<ELEMENT SYMBOL>.dat
    !C
    !C           IF <CUT> = 0 THEN _ev<CUT> IS DELETED FROM ABOVE FILES.
    !C
    !C INPUT  : (C*2)  YEAR      = YEAR OF DATA
    !C          (C*2)  YEARDF    = DEFAULT YEAR OF DATA IF REQUESTED YEAR
    !C                             DOES NOT EXIST.
    !C          (I*4)  IZ0       = NUCLEAR CHARGE
    !C          (I*4)  IZ1       = MINIMUM ION CHARGE + 1
    !C          (I*4)  ICLASS    = CLASS OF DATA (1 - 7)
    !C          (I*4)  ITMAX     = NUMBER OF ( DTEV() , DDENS() ) PAIRS
    !C          (I*4)  IEVCUT    = ENERGY CUT-OFF (EV)
    !C          (I*4)  ITDIMD    = MAXIMUM NUMBER OF DATA TEMP & DENS
    !C          (R*8)  DTEV()    = DLOG10(ELECTRON TEMPERATURES (EV))
    !C          (R*8)  DDENS()   = DLOG10(ELECTRON DENSITIES (CM-3))
    !C
    !C OUTPUT : (C**)  TITLF     = INFORMATION STRING
    !C          (I*4)  ITMAXD    = NUMBER OF DATA DTEVD()
    !C          (I*4)  IDMAXD    = NUMBER OF DATA DDENS()
    !C          (I*4)  IZMAXD    = NUMBER OF DATA ZDATA()
    !C          (I*4)  ZDATA()   = Z1 CHARGES IN DATASET
    !C          (I*4)  IFAIL     = -1   IF ROUTINE SUCCESSFUL BUT THE DEFAULT
    !C                                  YEAR FOR THE DATA WAS USED.
    !C                           = 0    IF ROUTINE SUCCESSFUL - DATA FOR THE
    !C                                  REQUESTED YEAR USED.
    !C                           = 1    IF ROUTINE OPEN STATEMENT FAILED
    !C          (R*8)  DTEVD()   = DLOG10(DATA ELECTRON TEMPERATURES (EV))
    !C          (R*8)  DDENSD()  = DLOG10(DATA ELECTRON DENSITIES (CM-3))
    !C          (R*8)  DRCOFD()  = DLOG10(DATA RATE COEFFICIENTS (CM-3/S))
    !C          (R*8)  DRCOFI()  = INTERPOLATION OF DRCOFD(,,) FOR
    !C                             DTEV() & DDENS()
    !C
    !C PROGRAM: (C*2)  XFESYM    = FUNCTION - SEE ROUTINES SECTION BELOW
    !C          (C*2)  ESYM      = ELEMENT SYMBOL FOR NUCLEAR CHARGE IZ0
    !C          (C*60) USERID    = USER ID UNDER WHICH ADAS DATA IS STORED
    !C          (C*60) DSNAME    = FILE NAME ( SEE ABOVE TYPES )
    !C          (C*80) STRING    = GENERAL VARIABLE
    !C          (C*80) BLANK     = BLANK STRING
    !C          (C*2)  YEARSV    = LAST YEAR USED IN THIS ROUTINE
    !C          (I*4)  IREAD     = INPUT STREAM FOR OPEN STATEMENT
    !C          (I*4)  IZ0SV     = LAST IZ0 USED IN THIS ROUTINE
    !C          (I*4)  ICLSV     = LAST ICLASS USED IN THIS ROUTINE
    !C          (I*4)  INDXZ1    = LOCATION OF IZ1 IN ZDATA()
    !C          (I*4)  LCK       = MUST BE GREATER THAN 'ITMAXD' & 'IDMAXD'
    !C                             & 'ITMAX' - ARRAY SIZE FOR SPLINE CALCS.
    !C          (R*8)  A()       = GENERAL ARRAY
    !C          (R*8)  DRCOF0(,) = INTERPOLATION OF DRCOFD(,,) W.R.T DTEV()
    !C          (L*8)  LEXIST    = TRUE --- FILE TO OPEN EXISTS ELSE NOT
    !C
    !C PE BRIDEN = ADDED VARIABLES (14/01/91)
    !C
    !C          (I*4)  L1      = PARAMETER = 1
    !C          (I*4)  IOPT    = DEFINES THE BOUNDARY DERIVATIVES FOR THE
    !C                             SPLINE ROUTINE 'XXSPLE', SEE 'XXSPLE'.
    !C
    !C          (L*4)  LSETX   = .TRUE.  => SET UP SPLINE PARAMETERS RELATING
    !C                                      TO X-AXIS.
    !C                           .FALSE. => DO NOT SET UP SPLINE PARAMETERS
    !C                                      RELATING TO X-AXIS.
    !C                                      (I.E. THEY WERE SET IN A PREVIOUS
    !C                                            CALL )
    !C                           (VALUE SET TO .FALSE. BY 'XXSPLE')
    !C
    !C
    !C          (R*8)  DY()    = SPLINE INTERPOLATED DERIVATIVES
    !C
    !C          (R*8 ADAS FUNCTION - 'R8FUN1' ( X -> X) )
    !C
    !C PE BRIDEN = ADDED VARIABLES (23/04/93)
    !C
    !C          (I*4 ADAS FUNCTION - 'I4UNIT' (OUTPUT STREAM))
    !C
    !C AUTHOR : JAMES SPENCE (TESSELLA SUPPORT SERVICES PLC)
    !C          K1/0/80
    !C          JET  EXT. 4866
    !C
    !C DATE   : 22/02/90
    !C
    !C DATE   : 21/08/90 PE BRIDEN - REVISION: SEQUA(43) CHANGED ('TE'->'TC')
    !C
    !C DATE   : 08/10/90 PE BRIDEN - REVISION: RENAMED SUBROUTINE
    !C
    !C DATE   : 12/11/90 PE BRIDEN - CORRECTION: MOVE THE SETTING OF 'INDXZ1'
    !C                                           TO AFTER THE  '20 CONTINUE'
    !C                                           STATEMENT.   ALSO SAVE  THE
    !C                                           VALUE OF 'IZ1MIN'.
    !C
    !C DATE   : 14/01/91 PE BRIDEN - ADAS91:     CALLS TO NAG SPLINE ROUTINES
    !C                                           'E01BAF' & 'E02BBF' REPLACED
    !C                                           BY  CALLS   TO  ADAS  SPLINE
    !C                                           ROUTINE 'XXSPLN'.
    !C
    !C DATE   : 25/06/91 PE BRIDEN - CORRECTION: CHANGED FOLLOWING DIMENSION:
    !C                                            'DIMENSION DRCOFI(ITDIMD)'
    !C                                           TO
    !C                                            'DIMENSION DRCOFI(ITMAX)'
    !C
    !C DATE   : 07/08/91 PE BRIDEN - ADDED ERROR HANDLING IF THE OPEN STATE-
    !C                               MENT FAILS. (IFAIL=1 RETURNED)
    !C
    !C DATE   : 27/04/92 PE BRIDEN - ADDED DEFAULT YEAR FOR DATA IF REQUESTED
    !C                               YEAR DOES NOT EXIST. (ADDED 'YEARDF')
    !C                               INTRODUCED IFAIL = -1 IF DEFAULT YEAR
    !C                               WAS USED AND NOT THE REQUESTED YEAR.
    !C
    !C DATE   : 10/03/93 PE BRIDEN - ALLOWED INPUT DATA SETS TO BE ACCESSED
    !C                               FROM ANY USERID (DEFAULT = JETSHP)
    !C                               - INTRODUCED USERID VARIABLE AND CALL
    !C                                 TO XXUID.
    !C
    !C DATE   : 23/04/93 PE BRIDEN - ADDED I4UNIT FUNCTION TO WRITE
    !C                               STATEMENTS FOR SCREEN MESSAGES
    !C
    !C UPDATE:  24/05/93 - PE BRIDEN - ADAS91: CHANGED I4UNIT(0)-> I4UNIT(-1)
    !C
    !C UPDATE:  14/09/94 - PE BRIDEN - ADAS91: ADDED CHECK TO MAKE SURE THAT
    !C                                         ITMAX, ITMAXD AND IDMAXD ARE
    !C                                         IN RANGE (I.E. <= LCK).
    !C
    !C DATE   : 17/03/95 LD HORTON - MODIFIED FOR UNIX.  CLEANED UP FILE
    !C                               HANDLING
    !C
    !C DATE   :  3/08/95 LD HORTON - REPLACED XXSPLN WITH XXSPLE
    !C
    !C DATE   :  6/12/95 LD HORTON - MOVED LINTRP TO BE LOCAL TO MAINTAIN
    !C                               COMPATIBILITY WITH MAINFRAME
    !C
    !C-----------------------------------------------------------------------
    !C
    INTEGER, parameter ::  L1=1, MCLASS=7,iread=4999,lck=100
    !C
    !C
    integer*4, intent(in) :: iz0,iz1,iclass,itmax,ievcut,itdimd
    CHARACTER(len=2), intent(in) :: YEAR, YEARDF
    real*8, intent(in) :: DTEV(ITMAX)   , DDENS(ITMAX)
    character(len=*), intent(out) ::  TITLF
    integer*4, intent(out) :: itmaxd,idmaxd,izmaxd,ifail
    real*8, intent(out) :: ZDATA(ITDIMD)
    real*8, intent(out) ::  DTEVD(ITDIMD) , DDENSD(ITDIMD)
    real*8, intent(out) ::  DRCOFD(ITDIMD,ITDIMD,ITDIMD), DRCOFI(ITMAX)

    INTEGER ::  I4UNIT,izmax,iz1max,id,it,iz,indxz1
    INTEGER ::  IOPT
    INTEGER ::   LENF1, LENF2, LENF3, LENF4, LENF5, LENF6
    real*8 ::  A(LCK)
    REAL*8 ::  DY(LCK)
    real*8 ::  DRCOF0(LCK,LCK)
    LOGICAL ::  LINTRP(LCK)
    !C
    character(len=2) :: YEARSV='  ',ESYM,XFESYM*2
    character(len=6) :: EVCUT
    character(len=5) :: defadf
    character(len=80) :: USERID, DSNAME, STRING, BLANKS
    CHARACTER(len=3) :: CLASS(MCLASS) = (/'acd', 'scd', 'ccd', 'prb', 'plt', 'prc', 'pls'/)

    LOGICAL ::  LEXIST  , LSETX
    !C
    EXTERNAL  ::  R8FUN1
    !C
    integer, save :: IZ1MIN
    !C
    !C------SET DEFAULT DIRECTORY--------------------------------------------------------
    !C
    PARAMETER (DEFADF='adf11')
    !C
    !C-----------------------------------------------------------------------
    !C
    integer, save ::  IZ0SV=0,ICLSV=0,IEVSV=0

    !C
    !C------DIMENSION CHECK--------------------------------------------------
    !C
    IF (LCK.LT.ITMAX) STOP ' D2DATA ERROR: ITMAX > 100 (LCK): DECREASE ITMAX'
    !C
    !C-----------------------------------------------------------------------
    !C
    IF(.not.(YEAR.EQ.YEARSV.AND.IZ0.EQ.IZ0SV.AND.ICLASS.EQ.ICLSV.AND.IEVCUT.EQ.IEVSV)) then
       YEARSV = YEAR
       IZ0SV  = IZ0
       ICLSV  = ICLASS
       IEVSV  = IEVCUT
       IFAIL  = 0
       !C
       !C------ENERGY CUTOFF----------------------------------------------------
       !C
       IF( IEVCUT.GT.0 ) THEN
          WRITE(EVCUT,1050) IEVCUT
          CALL XXSLEN(EVCUT,LENF5,LENF6)
       END IF
       !C
       !C------GET ADAS DATA SOURCE USERID--------------------------------------
       !C
       USERID = '?'
       CALL XXUID(USERID)
       CALL XXSLEN(USERID,LENF1,LENF2)
       !C
       !C------ELEMENT NAME-----------------------------------------------------
       !C
       ESYM = XFESYM(IZ0)
       CALL XXSLEN(ESYM,LENF3,LENF4)
       !C
       !C------FILE NAME--------------------------------------------------------
       !C
       do
          DSNAME=USERID(LENF1:LENF2)//'/'//DEFADF//'/'//CLASS(ICLASS)// &
               YEARSV//'/'//CLASS(ICLASS)//YEARSV//'_'// &
               ESYM(LENF3:LENF4)//'.dat'
          IF ((ICLASS.GE.4.AND.ICLASS.LE.6) .AND. IEVCUT.NE.0) THEN
             DSNAME=USERID(LENF1:LENF2)//'/'//DEFADF//'/'//CLASS(ICLASS)// &
                  YEARSV//'/'//CLASS(ICLASS)//YEARSV//'_'// &
                  ESYM(LENF3:LENF4)//'_ev'//EVCUT(LENF5:LENF6)//'.dat'
          ENDIF
          !C
          !C------PE BRIDEN - MODIFICATION 27/04/92 - INCLUSION OF DEFAULT YEAR -
          !C
          !C------DOES FILE TO BE OPEN EXIST OR NOT--------------------------------
          !C
          INQUIRE(FILE=DSNAME,EXIST=LEXIST)
          !C
          IF ( (.NOT.LEXIST) .AND. (YEARSV.NE.YEARDF) ) THEN
             WRITE(I4UNIT(-1),1060) DSNAME , YEARDF
             IFAIL  = -1
             YEARSV = YEARDF
             cycle
          ENDIF
          exit
       enddo
       !C
       IF( .NOT.LEXIST ) GOTO 9999
       !C
       TITLF=BLANKS
       WRITE(TITLF,1000) DSNAME
       !C
       !C------PE BRIDEN - END OF MODIFICATION 27/04/92
       !C
       !C------READ FILE # IREAD------------------------------------------------
       !C
       !C       OPEN(UNIT=IREAD,FILE=DSNAME,ACTION='READ',ERR=9999)
       OPEN(UNIT=IREAD,FILE=DSNAME,ERR=9999)
       !C
       READ(IREAD,1010) IZMAX , IDMAXD , ITMAXD , IZ1MIN , IZ1MAX
       !C
       READ(IREAD,1020) STRING
       READ(IREAD,1040) ( DDENSD(ID) , ID = 1 , IDMAXD )
       READ(IREAD,1040) ( DTEVD(IT)  , IT = 1 , ITMAXD )
       !C
       IZMAXD = 0
       DO IZ = IZ1MIN , IZ1MAX
          IZMAXD = IZMAXD + 1
          ZDATA(IZMAXD) = IZ
          !C         IF( IZ .EQ. IZ1 ) INDXZ1 = IZ
          READ(IREAD,1020)STRING
          DO IT = 1 , ITMAXD
             READ(IREAD,1040) ( DRCOFD(IZMAXD,IT,ID) , ID = 1 , IDMAXD )
          enddo
       enddo
       !C
       CLOSE(IREAD)
       !C
       TITLF = STRING
       !C
       !C------INTERPOLATE USING SPLINES (NAG ALGORITHM)------------------------
       !C
       IF ( (LCK.LT.ITMAXD) .OR. (LCK.LT.IDMAXD) ) STOP ' D2DATA ERROR: ITMAXD AND/OR IDMAXD > 100 (LCK): INCREASE LCK'
       !C

    endif
    !C
    !C------PE BRIDEN - CORRECTION 12/11/90 - SET INDXZ1 AFTER '20 CONTINUE'-
    !C
    INDXZ1 = IZ1 - IZ1MIN + 1
    !C
    !C-----------------------------------------------------------------------
    !C
    !C
    !C>>>>>>INTERPOLATE DRCOFD(,,,) W.R.T TEMPERATURE
    !C
    LSETX = .TRUE.
    IOPT  = -1
    !C
    DO ID = 1 , IDMAXD
       !csw
       if(indxz1 <= 0) cycle
       !csw
       !C
       DO IT = 1 , ITMAXD
         A(IT) = DRCOFD(INDXZ1,IT,ID)
       enddo
       !C
       CALL XXSPLE( LSETX  , IOPT  , R8FUN1       , &
            ITMAXD , DTEVD , A            ,&
            ITMAX  , DTEV  , DRCOF0(1,ID) ,&
            DY     , LINTRP )
       !C
    enddo
    !C
    !C>>>>>>INTERPOLATE ABOVE RESULT W.R.T DENSITY
    !C
    LSETX = .TRUE.
    IOPT  = -1
    !C
    DO IT = 1 , ITMAX
       !C
       DO ID = 1 , IDMAXD
          A(ID) = DRCOF0(IT,ID)
       enddo
       !C
       CALL XXSPLE( LSETX  , IOPT      , R8FUN1     , &
            IDMAXD , DDENSD    , A          ,&
            L1     , DDENS(IT) , DRCOFI(IT) ,&
            DY     , LINTRP )
       !C
    enddo
    !C
    RETURN
    !C
    !C-----------------------------------------------------------------------
    !C DATA SET OPENING/EXISTENCE ERROR HANDLING
    !C-----------------------------------------------------------------------
    !C
9999 IFAIL  = 1
    YEARSV = '  '
    IZ0SV  = 0
    ICLSV  = 0
    RETURN
    !C
    !C-----------------------------------------------------------------------
    !C
1000 FORMAT('FILE = ',1A60)
1010 FORMAT(5I5)
1020 FORMAT(1A80)
1040 FORMAT(8F10.5)
1050 FORMAT(I6)
1060 FORMAT(1X,'NOTE: REQUESTED DATASET - ',A50,' DOES NOT EXIST.'/ 7 X,      'USING DEFAULT YEAR (',A2,') DATASET INSTEAD'/)
    !C
    !C-----------------------------------------------------------------------
    !C
  END SUBROUTINE D2DATA_EIRENE



end module mod_adas
