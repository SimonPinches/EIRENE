!  03.08.06:  data structure for reaction data redefined
!  25.04.07:  reading of rate coefficients from HYDKIN database added
c  changed in 2011:  new atomic/molecular data structure introduced,
c                       REACDAT(IR)% ...
C
c  at the end of this routine, for each reaction card, call: SET_REACTION_DATA.F
cdr  jan.14: started to comment, cleanup
cdr  april 2015: further commenting cleanup, nov. 15: continued
cdr  jan 16: started to document options for asymptotics

cdr:  possible conflict with file fort.29, which is also used in coupling to B2
cdr:  subr. infcop.f, there to provide extra information regarding grid distortion
C
C
      SUBROUTINE EIRENE_SLREAC (IR,FILNAM,H123,REAC,CRC,
     .                   RCMIN, RCMAX, FP, JFEXMN, JFEXMX, ELNAME, IZ1)
c
c  open data stream 29 and read atomic data set no. IR
c          (note: general input-stream/output-stream no. offset ifoff
c           may have been set (for entire eirene run),
c           then stream is "29+ifoff".  default: ifoff=0)
c  and at the end: call to SET_REACTION_DATA.F (in module COMXS) 
c                  to fill REACDAT data structure
c
c
C  input
c    IR    : store data on eirene array CREAC(...,...,IR)

c
c
c
C    FILNAM: read a&m data from file filnam,
c            FILNAM=AMJUEL, HYDHEL, METHAN, H2VIBR, CONST
CC           FILNAM=ADAS:  special treatment, see below.
C            FILNAM=H-COL: nothing to be done here, use internal CR code h-colrad.f
C            FILNAM=HYDRTC: nothing to be done here  ??
C
c    H123  : identifyer for data type in filnam, e.g. H.1, H.2, H.3, ...


c    REAC  : in case FILNAM.ne.CONST:
c               number of reaction in data file "filnam", e.g. 2.2.5
c               and parameter fit-flag is found from the datafile (if available)
c    REAC  : in case FILNAM.eq.CONST:
c               reac IS MIS-USED AS  fit-flag: iftflg.
C               not NICE, VERY CONFUSING.
C               BETTER MAKE AN OWN INPUT PARAMETER IFTFLG IN CASE OPTION FILNAM= "CONST"

c  what does that mean for H-COL? ADAS ?  what about "spectral database"?  where described, where read ?

c  is iftflg not known in case of AMJUEL?

C            in case FILNAM=ADAS:  the file name DSN = REAC_ELNAME.dat is opened (stream 29+ifoff)
C                                  and then subroutine read_adas.f is called.
c    CRC   : type of process, e.g. EI, CX, OT, etc.
c
c  parameters for extrapolation beyond specified range [RMN,RMX] of data (asymptotics),
c  these asymptotics parameters may already have been read from input file, block 4., subroutine input.f
c  or, if not, they will be searched for on the A&M data files read here.
c
c    RCMIN: LOG(RMN), RMN: lower boundary for indep. dependent variable (energy, temperature, density)
c           Default: RMX = exp(-20.)
c    RCMAX: LOG(RMX), RMX: upper boundary for indep. dependent variable (energy, temperature, density)
c           Default: RMX = exp(20.)
c    FP     Fitting coefficients for extrapolation (three for MIN and three for MAX
c
c    JFEXMN Flag for selecting extrapolation expression, left end (minimum)
c           =0  :  no data yet, try to read extrapolation from atomic data file here
c           else:  extrapolation is set explicitly in input file, block 4a
c                  skip reading extrapolation data from data file, even if they are available
c    JFEXMX Flag for selecting extrapolation expression, right end (maximum)
c           =0  :  no data yet, try to read extrapolation from atomic data file here
c           else:  extrapolation is set explicitly in input file, block 4a
c                  skip reading extrapolation data from data file, even if they are available

c  specific input, only available in case FILNAM=ADAS
c    ELNAME:  only in case FILNAM=ADAS: the new file name REAC_ELNAME is construced
C    IZ1   :  only in case FILNAM=ADAS:    ???????????????  ion charge in filename ?????

C  internal
C    ISW   <-- H123
C    IO    derived from ISW, initial value of 2nd index in CREAC-array

C  output
c    ISWR  : eirene flag for type of process  (1,2,...7), coding EI,CX,EL,...

c    CREAC : (old version) eirene storage array for a&m data CREAC(9,-1:9,IR)
c    REACDAT(IR)%.... (new version) eirene atomic data structure.

c    MODCLF: see below: further information on input a&m data structure
c    DELPOT: ionisation potential (for H.10 data),
c            currently handeled in input.f. not nice! also missing still for: H.8, H.9
c
C    IFTFLG=IFTFLG(IR,IFLG)
C    IFLG  derived from ISW
C          0 for potential, 1 for cross section, 2 for rate-coeff,
C          3 for mom-weighted rate coeff. 4 for energy weighted rate coeff.
c    IFTFLG: eirene flag for type of fitting expression ("fit-flag=...")
c
c   iftflg only kown for filname=const ?????
c
c
C            DEFAULTS: =2 IFLG=0
C                         FOR POTENTIAL (GEN. MORSE)
C                         IFLG > 0:
C                         NOT IN USE
C                      =0 FOR ALL OTHERS (POLYNOM, DOUBLE POLYNOM)
C                      =3, IFLG=1:
C                          ionisation/excitation cross section formula (METHANE,...)
C                          IFLG > 1, IFLG=0:
C                          NOT IN USE
C
C                      =L10 (L=0,1): ONLY ONE CONSTANT RATE OR RATE-COEFF.
C                      =LMN L=0: rate coefficient.
c                                multiply with density
C                      =LMN L=1: rate, not rate coefficient.
c                                no multiplication of density
c
C  READ A&M DATA FROM THE FILES INTO EIRENE ARRAY CREAC
C
C
C  OUTPUT (IN COMMON COMXS):
C    READ DATA FROM "FILNAM" INTO ARRAY "CREAC"
C    DEFINE PARAMETER MODCLF(IR) (5 DIGITS NMLKJ)
C    FIRST DEZIMAL  J           =1  POTENTIAL AVAILABLE
C                                   (ON CREAC(..,-1,IR))
C                   J           =0  ELSE
C    SECOND DEZIMAL K           =1  CROSS SECTION AVAILABLE
C                                   (ON CREAC(..,0,IR))
C                   K           =0  ELSE
C    THIRD  DEZIMAL L           =1  <SIGMA V> FOR ONE
C                                   PARAMETER E (E.G.
C                                   PROJECTILE ENERGY OR ELECTRON
C                                   DENSITY) AVAILABLE
C                                   (ON CREAC(..,1,IR))
C                               =2  <SIGMA V> FOR
C                                   9 PROJECTILE ENERGIES AVAILABLE
C                                   (ON CREAC(..,J,IR),J=1,9)
C                               =3  <SIGMA V> FOR
C                                   9 ELECTRON DENSITIES  AVAILABLE
C                                   (ON CREAC(..,J,IR),J=1,9)
C                   L           =0  ELSE
C    FOURTH DEZIMAL M               DATA FOR MOMENTUM EXCHANGE
C                                   TO BE WRITTEN
C    FIFTH  DEZIMAL N           =1  DELTA E FOR ONE PARAMETER E (E.G.
C                                   PROJECTILE ENERGY OR ELECTRON
C                                   DENSITY) AVAILABLE
C                                   (ON CREAC(..,1,IR))
C                               =2  DELTA E FOR
C                                   9 PROJECTILE ENERGIES AVAILABLE
C                                   (ON CREAC(..,J,IR),J=1,9)
C                               =3  DELTA E FOR
C                                   9 ELECTRON DENSITIES  AVAILABLE
C                                   (ON CREAC(..,J,IR),J=1,9)
C                   N           =0  ELSE
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMPRT
      USE EIRMOD_COMXS
      USE EIRMOD_CINIT
      USE EIRMOD_PHOTON

      IMPLICIT NONE

      INTEGER,      INTENT(IN) :: IR, IZ1
      INTEGER,      INTENT(IN OUT) :: JFEXMN, JFEXMX
      CHARACTER(8), INTENT(IN) :: FILNAM
      CHARACTER(4), INTENT(IN) :: H123
      CHARACTER(LEN=*), INTENT(IN) :: REAC, ELNAME
      CHARACTER(3), INTENT(IN) :: CRC
      REAL(DP), INTENT(IN OUT) :: RCMIN, RCMAX, FP(6)
      CHARACTER(50) :: REACSTR
      REAL(DP) :: CONST, E_EL, E_K
      REAL(DP) :: CREACD(9,9)  ! INTERMEDIATE STORAGE FOR FIT PARAMETERS
      INTEGER :: I, IND, J, K, IH, I0P1, I0, IC, IREAC, ISW, INDFF,
     .           IFLG, INC, IANF, IFILE, IL
      CHARACTER(80) :: ZEILE
      CHARACTER(2) :: CHR
      CHARACTER(3) :: CHRL, CHRR
      CHARACTER(200) :: DSN, DIR
      CHARACTER(1) :: CUT
      CHARACTER(4) :: CH123
      CHARACTER(3) :: CCRC
      LOGICAL :: LCONST,LGEMIN,LGEMAX
C
!   set some defaults
      LGEMIN=.FALSE.
      LGEMAX=.FALSE.
      ISWR(IR)=0
      CONST=0.
      CHR='l0'
      CHRL='ll0'
      CHRR='lr0'
      I0=0
      CREACD = 0._DP
C
      IF (INDEX(CRC,'EI').NE.0.OR.
     .    INDEX(CRC,'DS').NE.0) ISWR(IR)=1
      IF (INDEX(CRC,'CX').NE.0) ISWR(IR)=3
      IF (INDEX(CRC,'II').NE.0.OR.
     .    INDEX(CRC,'PI').NE.0) ISWR(IR)=4
      IF (INDEX(CRC,'EL').NE.0) ISWR(IR)=5
      IF (INDEX(CRC,'RC').NE.0) ISWR(IR)=6
      IF (INDEX(CRC,'OT').NE.0) ISWR(IR)=7
C
      IF (INDEX(FILNAM,'CONST').NE.0) THEN
        LCONST=.TRUE.
!  nothing to be done
      ELSEIF (INDEX(FILNAM,'H-COL').NE.0) THEN
        LCONST=.FALSE.
!  nothing to be done
      ELSE
!  open data file, stream 29+ifoff.
        DO IFILE=1,NDBNAMES
          IF (INDEX(FILNAM,DBHANDLE(IFILE)).NE.0) EXIT
        END DO
        IF (IFILE <= NDBNAMES) THEN
C  proper filnam found
          LCONST=.FALSE.
          IF (INDEX(FILNAM,'ADAS') == 0) THEN
! FILNAM=AMJUEL, HYDHEL, METHAN, H2VIBR, PHOTON....: open data file
            OPEN (UNIT=29+ifoff,FILE=DBFNAME(IFILE))

          ELSEIF (INDEX(FILNAM,'ADAS').NE.0) THEN
! FILNAM=ADAS: open data file
! FIND NAME OF SPECIFIC ADAS-FILE TO BE READ,  DSN=abc.dat
!           reconstruct 'DSN' from:  reac, elname
            DIR = ' '
            IL = 0
            IF (VERIFY(DBFNAME(IFILE),' ') .NE. 0) THEN
              IC = SCAN(DBFNAME(IFILE),'/\\')
              CUT = DBFNAME(IFILE)(IC:IC)
              DIR=TRIM(DBFNAME(IFILE)) // CUT //
     .            ADJUSTL(TRIM(REAC)) // CUT
              IL = INDEX(DIR,CUT,.TRUE.)
            END IF
            IF (IL == 0) THEN
              DSN = ADJUSTL(TRIM(REAC)) // '_' // TRIM(ELNAME) // '.dat'
            ELSE
              DSN = DIR(1:IL) // ADJUSTL(TRIM(REAC)) // '_' //
     .              TRIM(ELNAME) // '.dat'
            END IF
            OPEN (UNIT=29+ifoff,FILE=DSN)
          END IF

C  THE A&M DATA FILE FILNAM IS NOW OPENDED, ON STREAM 29 (+ifoff)

        ELSE
          WRITE (iunout,*)
     .      ' NO VALID FILENAME IN REACTION CARD'
          WRITE (iunout,*) ' CHOOSE EITHER '
          WRITE (iunout,*) ' AMJUEL, METHAN, HYDHEL, H2VIBR, PHOTON '
          WRITE (iunout,*) ' OR '
          WRITE (iunout,*) ' ADAS '
          WRITE (iunout,*) ' OR '
          WRITE (iunout,*) ' H-COL'
          WRITE (iunout,*) ' OR '
          WRITE (iunout,*) ' CONST FOR ENTERING REACTION DATA VIA '
          WRITE (iunout,*) ' EIRENE INPUT-FILE '
          CALL EIRENE_EXIT_OWN(1)
        END IF
      ENDIF
C
      IF (H123(4:4).EQ.' ') THEN
        READ (H123(3:3),'(I1)') ISW
      ELSE
        READ (H123(3:4),'(I2)') ISW
      ENDIF

C
      IF (INDEX(FILNAM,'PHOTON').NE.0) THEN
        CALL EIRENE_READ_PHOTDBK (IR,REAC,ISW)
        RETURN
      END IF

      REACSTR=REPEAT(' ',11)
      IANF=VERIFY(REAC,' ')
      IF (IANF > 0) THEN
        IREAC=INDEX(REAC(IANF:),' ')-1
        IF (IREAC.LT.0) IREAC=LEN(REAC(IANF:))
        REACSTR(2:IREAC+1)=REAC(IANF:IREAC+IANF-1)
C  ADD ONE MORE BLANK, IF POSSIBLE
        IREAC=IREAC+2
      ELSE
        IF (.NOT.LCONST) THEN
          WRITE (iunout,*) ' NO REACTION SPECIFIED IN REACTION CARD ',IR
          CALL EIRENE_EXIT_OWN (1)
        END IF
      END IF

      REAC_NAME(IR) = REACSTR(2:)

C  H.0
      IF (ISW.EQ.0) THEN
        CHR='p0'
        CHRL='pl0'
        CHRR='pr0'
        I0=-1
        MODCLF(IR)=MODCLF(IR)+1
        IFLG=0
C  DEFAULT POTENTIAL: GENERALISED MORSE
        IFTFLG(IR,IFLG)=2
C  H.1
      ELSEIF (ISW.EQ.1) THEN
        CHR='a0'
        CHRL='al0'
        CHRR='ar0'
        I0=0
        MODCLF(IR)=MODCLF(IR)+10
        IFLG=1
C  DEFAULT CROSS SECTION: 8TH ORDER POLYNOM OF LN(SIGMA)
        IFTFLG(IR,IFLG)=0
C  H.2
      ELSEIF (ISW.EQ.2) THEN
        CHR='b0'
        CHRL='bl0'
        CHRR='br0'
        I0=1
        MODCLF(IR)=MODCLF(IR)+100
        IFLG=2
C  DEFAULT RATE COEFFICIENT: 8TH ORDER POLYNOM OF LN(<SIGMA V>) FOR E0=0.
        IFTFLG(IR,IFLG)=0
C  H.3
      ELSEIF (ISW.EQ.3) THEN
        MODCLF(IR)=MODCLF(IR)+200
        I0=1
        IFLG=2
        IFTFLG(IR,IFLG)=0
C  H.4
      ELSEIF (ISW.EQ.4) THEN
        MODCLF(IR)=MODCLF(IR)+300
        I0=1
        IFLG=2
        IFTFLG(IR,IFLG)=0
C  H.5
      ELSEIF (ISW.EQ.5) THEN
        CHR='e0'
        CHRL='el0'
        CHRR='er0'
        I0=1
        MODCLF(IR)=MODCLF(IR)+1000
        IFLG=3
        IFTFLG(IR,IFLG)=0
C  H.6
      ELSEIF (ISW.EQ.6) THEN
        MODCLF(IR)=MODCLF(IR)+2000
        I0=1
        IFLG=3
        IFTFLG(IR,IFLG)=0
C  H.7
      ELSEIF (ISW.EQ.7) THEN
        MODCLF(IR)=MODCLF(IR)+3000
        I0=1
        IFLG=3
        IFTFLG(IR,IFLG)=0
C  H.8
      ELSEIF (ISW.EQ.8) THEN
        CHR='h0'
        CHRL='hl0'
        CHRR='hr0'
        I0=1
        MODCLF(IR)=MODCLF(IR)+10000
        IFLG=4
        IFTFLG(IR,IFLG)=0
C  H.9
      ELSEIF (ISW.EQ.9) THEN
        MODCLF(IR)=MODCLF(IR)+20000
        I0=1
        IFLG=4
        IFTFLG(IR,IFLG)=0
C  H.10
      ELSEIF (ISW.EQ.10) THEN
        MODCLF(IR)=MODCLF(IR)+30000
        I0=1
        IFLG=4
        IFTFLG(IR,IFLG)=0
C  H.11
      ELSEIF (ISW.EQ.11) THEN
        CHR='k0'
        CHRL='kl0'
        CHRR='kr0'
        I0=1
        IFLG=5
        IFTFLG(IR,IFLG)=0
C  H.12
      ELSEIF (ISW.EQ.12) THEN
        I0=1
        IFLG=5
        IFTFLG(IR,IFLG)=0
      ENDIF

      IF (INDEX(FILNAM,'H-COL').NE.0) THEN
        SELECT CASE (ISW)
        CASE (2:4)
          IF (REACDAT(IR)%LRTC) THEN
            WRITE (IUNOUT,*) ' RATE COEFFICIENT ALREADY SPECIFIED',
     .                       ' FOR REACTION', IR
            WRITE (IUNOUT,*) ' CHECK SPECIFICATION OF REACTIONS'
            CALL EIRENE_EXIT_OWN(1)
          END IF
          ALLOCATE (REACDAT(IR)%RTC)
          NULLIFY(REACDAT(IR)%RTC%ADAS)
          NULLIFY(REACDAT(IR)%RTC%LINE)
          NULLIFY(REACDAT(IR)%RTC%POLY)
          NULLIFY(REACDAT(IR)%RTC%HYD)
          REACDAT(IR)%LRTC = .TRUE.
          REACDAT(IR)%RTC%IFIT = 5
        CASE (8:10)
          IF (REACDAT(IR)%LRTCEW) THEN
            WRITE (IUNOUT,*) ' ENERGY WEIGHTED RATE COEFFICIENT',
     .                       ' ALREADY SPECIFIED FOR REACTION', IR
            WRITE (IUNOUT,*) ' CHECK SPECIFICATION OF REACTIONS'
            CALL EIRENE_EXIT_OWN(1)
          END IF
          ALLOCATE (REACDAT(IR)%RTCEW)
          NULLIFY(REACDAT(IR)%RTCEW%ADAS)
          NULLIFY(REACDAT(IR)%RTCEW%LINE)
          NULLIFY(REACDAT(IR)%RTCEW%POLY)
          NULLIFY(REACDAT(IR)%RTCEW%HYD)
          REACDAT(IR)%LRTCEW = .TRUE.
          REACDAT(IR)%RTCEW%IFIT = 5
        CASE DEFAULT
          WRITE (IUNOUT,*) ' WRONG REACTION TYPE SPECIFIED '
          WRITE (IUNOUT,*) ' REACTION NO. ', IR
          WRITE (IUNOUT,*) ' REACTION TYPE H.', ISW
          CALL EIRENE_EXIT_OWN(1)
        END SELECT
        RETURN
      END IF

      IF (INDEX(FILNAM,'ADAS').NE.0) THEN
        CALL EIRENE_READ_ADAS (IR,REAC,ISW,IZ1)
c  close unit=29+ifoff:   done in READ_TABLE2.f
        RETURN
      END IF

      IF (INDEX(FILNAM,'HYDRTC').NE.0) THEN
        CLOSE (UNIT=29+ifoff)
        CH123 = H123
        CCRC = CRC
        CALL EIRENE_READ_HYDKIN
     .      (IR,DBFNAME(IFILE),CH123,REAC,CCRC,RCMIN,RCMAX,
     .       E_EL,E_K,.FALSE.)
        RETURN
      END IF
C
      IF (LCONST) THEN
        IND=INDEX(REACSTR,'FT')
        IF (IND /= 0) THEN
          READ (REACSTR(IND+2:),*) IFTFLG(IR,IFLG)
        END IF

        IF (MOD(IFTFLG(IR,IFLG),100) == 10) THEN
C
C  READ ONLY ONE FIT COEFFICIENT FROM INPUT FILE 'iunin'
          READ (IUNIN,6664) CREACD(1,1)
          REACLINES(IRLINES)%NCONST = 1
          REACLINES(IRLINES)%CONST(1) = CREACD(1,1)
        ELSE
C
C  READ 9 FIT COEFFICIENTS FROM INPUT FILE 'iunin'
          READ (IUNIN,6664) (CREACD(IC,1),IC=1,9)
          REACLINES(IRLINES)%NCONST = 9
          REACLINES(IRLINES)%CONST(1:9) = CREACD(1:9,1)

        END IF
        CALL EIRENE_SET_REACTION_DATA(IR,ISW,IFTFLG(IR,IFLG),CREACD,
     .                         IUNOUT,.FALSE.)
        RETURN
      ENDIF
C
C  READ FROM DATA FILE, stream 29
C
C  already ruled out here (done at this point):
C  FILNAM= "H-COL", "CONST", "ADAS", "HYDRTC", "PHOTON"
C  in these cases: already returned to calling program
C
!dr   ELSEIF (.NOT.LCONST) THEN
C  AT THIS POINT: FILNAM= AMJUEL, HYDHEL, H2VIBR, METHAN, i.e. single or double polynomial fits

CC  now identify proper dataset within file FILNAM

100     READ (29+ifoff,'(A80)',END=990) ZEILE
        IF (INDEX(ZEILE,'##BEGIN DATA HERE##').EQ.0) GOTO 100

1       READ (29+ifoff,'(A80)',END=990) ZEILE
        IF (INDEX(ZEILE,H123).EQ.0) GOTO 1     !  infinite loop passible !
C
2       READ (29+ifoff,'(A80)',END=990) ZEILE
        IF (INDEX(ZEILE,'H.').NE.0) GOTO 990
        IF (INDEX(ZEILE,'Reaction ').EQ.0.or.
     .      INDEX(ZEILE,REACSTR(1:ireac)).EQ.0) GOTO 2  ! infinite loop possible  !
!dr   ENDIF
C
C  SINGLE PARAM. FIT, ISW=0,1,2,5,8,11
C
      IF (ISW.EQ.0.OR.ISW.EQ.1.OR.ISW.EQ.2.OR.ISW.EQ.5.OR.ISW.EQ.8.OR.
     .    ISW.EQ.11) THEN
C       IF (.NOT.LCONST) THEN  ! this is redundant, LCONST = TRUE is completely done already,
c                                already returned to calling program
3       READ (29+ifoff,'(A80)',END=990) ZEILE
        INDFF=INDEX(ZEILE,'fit-flag')
        IF (INDEX(ZEILE,CHR)+INDFF.EQ.0) GOTO 3
c  input line found which either contains fit-flag, or the reaction identifier a0,b0,...k0
        IF (INDFF > 0) THEN  !OTHERWISE: use DEFAULT FOR FIT-FLAG: iftflg = 0
c  read parameter for type of fitting expression from data file
          READ (ZEILE((INDFF+8):80),*) IFTFLG(IR,IFLG)
          GOTO 3
        ENDIF
C  read only one constant:  (FIT-FLAG = 10, 110, ....)
        IF (MOD(IFTFLG(IR,IFLG),100) == 10) THEN
          IND=INDEX(ZEILE,CHR(1:1))
          READ (ZEILE((IND+2):80),'(E20.12)') CREACD(1,1)
        ELSE
C  READ 9 FIT COEFFICIENTS, SEPARATED BY 'CHR'  FIXED FORMAT E20.12
C  THREE LINES WITH THREE DATA PER LINE
          DO 9 J=0,2
            IND=0
            DO 4 I=1,3
              IND=IND+INDEX(ZEILE((IND+1):80),CHR(1:1))
              READ (ZEILE((IND+2):80),'(E20.12)') CREACD(J*3+I,1)
4           CONTINUE
            READ (29+ifoff,'(A80)',END=990) ZEILE
9         CONTINUE
        END IF
C
C  READ ASYMPTOTICS FROM DATA FILE, IF AVAILABLE
C  I0P1=1 FOR CROSS SECTION
C  I0P1=2 FOR (WEIGHTED) RATE COEFFICIENT
        I0P1=I0+1
        IF (ISW.EQ.0) GOTO 12 ! NO ASYMPTOTICS FOR POTENTIALS

c  CHRL is label of left (low E,T) extrapolation fit parameters, a0l, b0l,....
        IF (INDEX(ZEILE,CHRL).NE.0.AND.JFEXMN.EQ.0) THEN
c  at this point: left extrapolation fit found in dataset fort.29, and
c                 left extrapolation was not overruled explicitly in input file 'fort.iunin'.

c  read three parameters FP(i), i=1,3 for 'left' extrapolation
          IND=0
          DO 5 I=1,3
            INC=INDEX(ZEILE((IND+1):80),CHR(1:1))
            IF (INC.GT.0) THEN
              IND=IND+INDEX(ZEILE((IND+1):80),CHR(1:1))
              READ (ZEILE((IND+3):80),'(E20.12)') FP(I)
            ENDIF
5         CONTINUE
          LGEMIN=.true.
          READ (29+ifoff,'(A80)',END=990) ZEILE
        ENDIF

c  same as above. for right (high E,T) extraploation fit
        IF (INDEX(ZEILE,CHRR).NE.0.AND.JFEXMX.EQ.0) THEN
c  read three parameters FP(i), i=4,6 for 'right' extrapolation
          IND=0
          DO 7 I=4,6
            INC=INDEX(ZEILE((IND+1):80),CHR(1:1))
            IF (INC.GT.0) THEN
              IND=IND+INDEX(ZEILE((IND+1):80),CHR(1:1))
              READ (ZEILE((IND+3):80),'(E20.12)') FP(I)
            ENDIF
7         CONTINUE
          LGEMAX=.true.
          READ (29+ifoff,'(A80)',END=990) ZEILE
        ENDIF
c


        if (lgemin.and.jfexmn.eq.0) then
c  at this point:  low end extrapolation parameter FP(1:3) have been read from atomic data file.
c  read value of lower validity bound, e.g. ELABMIN,...
c  and return as RCMIN
          IND=INDEX(ZEILE,'=')
          READ (ZEILE((IND+2):80),'(E12.5)') rcmin
          rcmin=log(rcmin)
          jfexmn=5
          READ (29+ifoff,'(A80)',END=990) ZEILE
        endif
        if (lgemax.and.jfexmx.eq.0) then
c  at this point:  high end extrapolation parameter FP(4:6) have been read from atomic data file.
c  read value of upper validity bound, e.g. ELABMAX,...
c  and return as RCMAX
          IND=INDEX(ZEILE,'=')
          READ (ZEILE((IND+2):80),'(E12.5)') rcmax
          rcmax=log(rcmax)
          jfexmx=5
          READ (29+ifoff,'(A80)',END=990) ZEILE
        endif
C
CDR  this part needs to be re-written and/or documented

C  ANY OTHER ASYMPTOTICS INFO ON FILE?  SEARCH FOR Tmin, or Emin
        IF ((INDEX(ZEILE,'Tmin').NE.0.and.I0P1==2).or.
     .      (INDEX(ZEILE,'Emin').NE.0.and.I0P1==1)) then
          IND=INDEX(ZEILE,'n')
          READ (ZEILE((IND+2):80),'(E9.2)') rcmin
          rcmin=log(rcmin)
C  extrapolation from subr. CROSS
          if (I0P1.eq.1.and.iswr(ir).eq.1) jfexmn=1
          if (I0P1.eq.1.and.iswr(ir).eq.3) jfexmn=-1
          if (I0P1.eq.1.and.iswr(ir).eq.5) jfexmn=-1
C  extrapolation from subr. CDEF
C   ??    if (I0PT.eq.2) jfexmn=-1
          READ (29+ifoff,'(A80)',END=990) ZEILE
        ENDIF
12      CONTINUE
C       ELSEIF (LCONST) THEN
C  NOTHING TO BE DONE
C       ENDIF

C   AT THIS POINT WE HAVE STORED FOR REACTION ir, DATA TYPE iflg:
C   IFTFLG(IR,iflg)   (DEFAUT:   =0)
C   9 FIT COEFFICIENTS ON INTERMEDIATE ARRAY CREACD(1...9,1)
C   AND POSSIBLY (SOME OF) THE EXTRAPOLATION PARAMETERS RCMIN,RCMAX, FP(1:6)
C
C  TWO PARAM. FIT, ISW=3,4,6,7,9,10,12
      ELSEIF (ISW.EQ.3.OR.ISW.EQ.4.OR.ISW.EQ.6.OR.ISW.EQ.7.OR.
     .        ISW.EQ.9.OR.ISW.EQ.10.OR.ISW.EQ.12) THEN
        DO 11 J=0,2
16        READ (29+ifoff,'(A80)',END=990) ZEILE
C  SEARCH FOR STRING 'fit-flag'  or 'Index'
          INDFF=INDEX(ZEILE,'fit-flag')
          IF (INDEX(ZEILE,'Index')+INDFF.EQ.0) GOTO 16
          IF (INDFF > 0) THEN
            READ (ZEILE((INDFF+8):80),*) IFTFLG(IR,IFLG)
            GOTO 16
          ENDIF
          READ (29+ifoff,'(1X)')
          IF (MOD(IFTFLG(IR,IFLG),100) == 10) THEN
C  IFTFLG = 10, 110,  210,....ETC:  READ ONLY ONE CONSTANT PARAMETER
            READ (29+ifoff,*) IH,CREACD(1,1)
            EXIT
          ELSE
            DO 17 I=1,9
C   READ BLOCK OF 9 LINES, THREE DATA EACH, UNFORMATTED
              READ (29+ifoff,*) IH,(CREACD(I,K),K=J*3+1,J*3+3)
17          CONTINUE
          END IF
11      CONTINUE
C   AT THIS POINT WE HAVE STORED FOR REACTION ir:
C   IFTFLG(IR)   (DEFAUT:   =0)
C   81 FIT COEFFICIENTS ON INTERMEDIATE ARRAY CREACD(1...9,1...9)

C   NO ASYMPTOTICS AVAILABLE YET FOR 2 PARAMETER FIT
C
      ENDIF

      CALL
     .  EIRENE_SET_REACTION_DATA(IR,ISW,IFTFLG(IR,IFLG),CREACD,IUNOUT,
     .                       .TRUE.,RCMIN,RCMAX,FP,JFEXMN,JFEXMX)
C
      CLOSE (UNIT=29+ifoff)
C
      RETURN
C
990   WRITE (iunout,*) ' NO DATA FOUND FOR REACTION ',H123,' ',REAC,
     .                 ' IN DATA SET ',FILNAM
      WRITE (iunout,*) ' IR,MODCLF(IR) ',IR,MODCLF(IR)
      CLOSE (UNIT=29+ifoff)
      CALL EIRENE_EXIT_OWN(1)
991   WRITE (iunout,*) ' INVALID CONSTANT IN SLREAC. CONST= ',CONST
      WRITE (iunout,*) ' CHECK "REACTION CARDS" FOR REACTION NO. ',IR
      CLOSE (UNIT=29+ifoff)
      CALL EIRENE_EXIT_OWN(1)
6664  FORMAT (6E12.4)
      END
