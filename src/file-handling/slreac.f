!  03.08.06:  data structure for reaction data redefined
!  25.04.07:  reading of rate coefficients from HYDKIN database added
c  changed in 2011:  new atomic/molecular data structure introduced,
c                    REACDAT(IR)% ..., replaces array CREAC(...)
C
c    at the end of this routine, for each reaction card, call: SET_REACTION_DATA(IR,..)
cdr  jan.14: started to comment, cleanup
cdr  april 2015: further commenting, cleanup, nov. 15: continued
cdr  jan 16: started to document options for asymptotics
!pb  apr 16: extensions to allow more precise comments in AMJUEL, HYDHEL, METHAN and H2VIBR
c            data files,
cdr          such as character strings H.xxx
cdr          taken over from ITER-IO branch
!pb  may 16: bug fix to the extensions (resolving problem reading HYDHEL H.3)

cdr:  possible conflict with file fort.29, which is also used in coupling to B2
cdr:  subr. infcop.f, there to provide extra information regarding grid distortion
cdr:  june 16:  added H.5 - H.7 options for H_COL case.
cdr            started to clarify extrapolation options for polynomial fits. Not ready
cdr            some comments corrected
cdr   Aug. 16: reading Tmin, Emin from hydhel disabled.
cdr            May have corrupted extrapolation in some cases
c     Sept.16: two new internal subroutines,
c              a) READ_RANGE:  to read validity range information,
c              b) READ_COEFFS: three parameters for each validity boundary, for extrapolation options
c    June  17: read_colrad (for old H-COL option (now CRM)) moved to separate routine.
cdr  Jan   19:  filnam=CRM --> CR... to prepare merge with branch ...emis....,
cdr             H, He internal CR codes, formulation I, II (MS resolved or not)
cdr  Feb   19:  remove obsolete (and unfinished) option HYDRTC 
C
C
      SUBROUTINE EIRENE_SLREAC (IR,FILNAM,H123,REAC,CRC,
     .                          RC1MIN, RC1MAX, FP1, JFEX1MN, JFEX1MX,
     .                          RC2MIN, RC2MAX, FP2, JFEX2MN, JFEX2MX,
     .                          ELNAME, IZ1,
     .                          IROW_ESC, ICOL_ESC, POP_ESC)
c
c  open data stream 29 and read atomic dataset no. IR
c          (note: general input-stream/output-stream no. offset ifoff
c           may have been set (for entire eirene run),
c           then stream is "29+ifoff".  default: ifoff=0)
c  and at the end: call to SET_REACTION_DATA.F (in module COMXS)
c                  to fill REACDAT data structure
c
c
C  input
c    IR    : store data on eirene data structure REACDAT(...,...,IR)

c
c
c
C    FILNAM: read A&M data from file filnam,
c            FILNAM=AMJUEL, HYDHEL, METHAN, H2VIBR, CONST: polynomial fits
CC           FILNAM=TAB2D, ADAS:  special treatment, see below.
C            FILNAM=CR...: nothing to be done here, use internal CR code xx_colrad.f
c                          currently available: h_colrad.f
C
c    H123  : identifier for data type in filnam, e.g. H.1, H.2, H.3, ...


c    REAC  : in case FILNAM = AMJUEL, HYDHEL, METHAN, H2VIBR:
c               number of reaction in data file "filnam", e.g. 2.2.5
c               and parameter fit-flag is found from the datafile (if available)
c    REAC  : in case FILNAM.eq.CONST:
c               reac IS MISUSED AS fit-flag: iftflg.
C               not NICE, VERY CONFUSING.
C               BETTER MAKE AN OWN INPUT PARAMETER IFTFLG IN CASE OPTION FILNAM= "CONST"

cdr what does that mean for CRM? ADAS ?  what about "spectral database"?
cdr where described, where read ?

C            in case FILNAM = TAB2D, ADAS:  the file name DSN = REAC_ELNAME.dat is opened
C                                           (stream 29+ifoff)
C                                  and then subroutine read_tab2d.f is called.

c    CRC   : type of process, e.g. EI, CX, EL, PI, RC, OT, etc.
c
c  parameters for extrapolation beyond specified range [RiMN,RiMX, i=1,2] of data (asymptotics),
c  these asymptotics parameters may already have been read from input file, block 4., subroutine input.f
c  or, if not, they will be searched for on the A&M data files read here.
c
c    for i=1,2:
c    RCiMIN: LOG(RiMN), RiMN: lower boundary for indep. dependent variable (energy, temperature, density)
c           Default: RiMN = exp(-20.)
c    RCiMAX: LOG(RiMX), RiMX: upper boundary for indep. dependent variable (energy, temperature, density)
c           Default: RiMX = exp(20.)
c    FPi    Fitting coefficients for extrapolation (three for MIN and three for MAX, each)
c
c    JFEXiMN Flag for selecting extrapolation expression, left end (minimum)
c            =0  :  no data yet, try to read extrapolation from atomic data file here
c            else:  extrapolation is set explicitly in input file, block 4a
c                  skip reading extrapolation data from data file, even if they are available
c    JFEXiMX Flag for selecting extrapolation expression, right end (maximum)
c            =0  :  no data yet, try to read extrapolation from atomic data file here
c            else:  extrapolation is set explicitly in input file, block 4a
c                  skip reading extrapolation data from data file, even if they are available

c  specific input, only available in case FILNAM=ADAS
c    ELNAME:  only in case FILNAM=ADAS: the new file name REAC_ELNAME is construced
C    IZ1   :  only in case FILNAM=ADAS: ionization stage number within ADAS file

C  internal
C    ISW   <-- H123:   H.0: ISW=0, H.1: ISW=1, H.2: ISW=2, ...,H.12: ISW=12
C    I0    derived from ISW, initial value of 2nd index in old CREAC arrays

C  output
c    ISWR  : eirene flag for type of process (1,2,...7), coding EI,CX,EL,PI,...

c    CREAC :          (old version) eirene storage array for a&m data CREAC(9,-1:9,IR)
c    REACDAT(IR)%.... (new version) eirene atomic data structure.

c    MODCLF: see below: further information on input a&m data structure
c    DELPOT: ionisation potential (for H.10 data),
c            currently handled in input.f. not nice! also missing still for: H.8, H.9
c
C    IFTFLG=IFTFLG(IR,IH): flag for type of fitting expression ("fit-flag=...")
C    IH  internally derived from ISW, for different types of data:
C          0 for interaction potential, or differential cross-sections (ISW=0)
C          1 for cross-section, (ISW=1)
C          2 for rate coeff, (ISW=2,3,4)
C          3 for mom-weighted rate coeff. (ISW=5,6,7)
C          4 for energy-weighted rate coeff. (ISW=8,9,10)
C          5 for other quantities, population densities, etc.. (ISW=11,12)
c
c
C       IFTFLG(IR,IH) DEFAULTS:

C           CASE IH=0  (H.0):
C       IFTFLG(IR,0)   =2,  FOR INTERACTION POTENTIAL (GEN. MORSE)

C           CASE  IH=1  (H.1):
C       IFTFLG(IR,1)   =0,  CROSS-SECTION (9-POLYNOMIAL)
C                      =3,  cross-section (ionisation/excitation cross-section
C                           formula (METHANE,...)
C           CASE  IH=2,3....,10 (H.2, H.3,....H.10)
C       IFTFLG(IR,...  =0,  FOR RATE COEFFICIENTS (9-POLYNOMIAL, 9X9-DOUBLE POLYNOMIAL)
C                      =10, FOR RATE COEFFICIENTS (CONSTANT)
C                      =100 FOR RATE, not rate coefficient,
c                      =110 FOR RATE, not rate coefficient, (CONSTANT)
c
C  READ A&M DATA FROM THE FILES INTO EIRENE ARRAY CREAC
C
C
C  OUTPUT (IN COMMON COMXS):
C    READ DATA FROM "FILNAM" INTO ARRAY "CREAC"
C    DEFINE PARAMETER MODCLF(IR) (5 DIGITS NMLKJ)
C    FIRST DECIMAL  J           =1  POTENTIAL AVAILABLE
C                                   (ON CREAC(..,-1,IR))
C                   J           =0  ELSE
C    SECOND DECIMAL K           =1  CROSS-SECTION AVAILABLE
C                                   (ON CREAC(..,0,IR))
C                   K           =0  ELSE
C    THIRD  DECIMAL L           =1  <SIGMA V> FOR ONE
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
C    FOURTH DECIMAL M               DATA FOR MOMENTUM EXCHANGE
C                                   TO BE WRITTEN
C    FIFTH  DECIMAL N           =1  DELTA E FOR ONE PARAMETER E (E.G.
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
      INTEGER,      INTENT(IN), OPTIONAL :: IROW_ESC, ICOL_ESC
      REAL(DP),     INTENT(IN), OPTIONAL :: POP_ESC

      CHARACTER(8), INTENT(IN) :: FILNAM
      CHARACTER(4), INTENT(IN) :: H123
      CHARACTER(LEN=*), INTENT(IN) :: REAC, ELNAME
      CHARACTER(3), INTENT(IN) :: CRC
cdr  asymptotics parameters already read from input block 4?
cdr  if not:  try to read from external A&M data file
cdr  in either case: store these on data structure REACDAT, in call to: set_reaction_data(IR,...)
      INTEGER,  INTENT(IN OUT) :: JFEX1MN, JFEX1MX, JFEX2MN, JFEX2MX
      REAL(DP), INTENT(IN OUT) :: RC1MIN, RC1MAX, FP1(6),
     .                            RC2MIN, RC2MAX, FP2(6)
cdr
      REAL(DP) :: RTMAX, ERTMAX, ETH
      CHARACTER(50) :: REACSTR
      REAL(DP) :: CONST, E_EL, E_K
      LOGICAL :: LCONST
      REAL(DP) :: CREACD(9,9)  ! INTERMEDIATE STORAGE FOR FIT PARAMETERS

cdr  for reading asymptotics parameters from data files
      INTEGER ::  IF1MN, IF1MX, IF2MN, IF2MX
      REAL(DP) :: FP1L(3), FP1R(3), FP2L(3), FP2R(3)
      REAL(DP) :: R1MN, R1MX, R2MN, R2MX

      INTEGER :: I, IND, J, K, IH, I0, IC, IREAC, ISW, INDFF,
     .           IFLG, IANF, IFILE, IL, INDG

      CHARACTER(80) :: ZEILE, LAST_TEX, ULINE
      CHARACTER(4) :: CHR, CETH, BEND
      CHARACTER(3) :: CH1L, CH1R, CH2L, CH2R
      CHARACTER(200) :: DSN, DIR
      CHARACTER(1) :: CUT, BACK
      CHARACTER(4) :: CH123
      CHARACTER(3) :: CCRC
      CHARACTER(8) :: SECTION, FITFLAG
      CHARACTER(7) :: C1L, C1R, C2L, C2R, CMR, CEMR
      LOGICAL :: LGC1MIN,LGC1MAX,LGC2MIN,LGC2MAX,
     .           LGR1MIN,LGR1MAX,LGR2MIN,LGR2MAX
C
! defining backslash character
      BACK="\\"
      SECTION=BACK // 'section'
      BEND=BACK // 'end'
C
!   set some defaults

      ISWR(IR)=0
      CONST=0.
      CHR=' l0 '

      I0=0
      CREACD = 0._DP
      FITFLAG ='fit-flag'

c  some additional  (optional) reaction data:  threshold energy,
c                                              max rate coeff sigma*v_rel,
c                                              at E_rel=ERTMAX
      CMR  = 'MAXRATE'
      CEMR = 'ELAB'
      CETH = 'ETH'
      RTMAX = 0._DP
      ERTMAX = -HUGE(1._DP)
      ETH = 0._DP

C     Defaults: no asymptotics
      CH1L='ll0'
      CH1R='lr0'
      CH2L='ll0'
      CH2R='lr0'
      C1L = 'XXMIN'
      C1R = 'XXMAX'
      C2L = 'YYMIN'
      C2R = 'YYMAX'
C
c  type (class) of reaction process

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
      ELSEIF (INDEX(FILNAM,'CR').NE.0) THEN
        LCONST=.FALSE.
!  nothing to be done
      ELSE   ! in all other cases: open data file, stream 29+ifoff
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
! FIND NAME OF SPECIFIC TAB2D FILE TO BE READ, DSN=abc.dat
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
          WRITE (iunout,*) ' AMJUEL, METHAN, HYDHEL, H2VIBR '
          WRITE (iunout,*) ' OR '
          WRITE (iunout,*) ' TAB1D, TAB2D '
          WRITE (iunout,*) ' OR '
          WRITE (iunout,*) ' CR'
          WRITE (iunout,*) ' OR '
          WRITE (iunout,*) ' CONST '
          WRITE (iunout,*) ' OR '
          WRITE (iunout,*) ' PHOTON'
          WRITE (iunout,*) ' FOR ENTERING REACTION DATA VIA '
          WRITE (iunout,*) ' EIRENE INPUT FILE '
          WRITE (iunout,*) ' FILNAM WAS : ', FILNAM
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
C
C Set character string identifiers CHR to search coefficients in data files.
C  H.0
      IF (ISW.EQ.0) THEN
        CHR=' p0 '
        I0=-1
        MODCLF(IR)=MODCLF(IR)+1
        IFLG=0
C  DEFAULT POTENTIAL: GENERALISED MORSE
        IFTFLG(IR,0)=2
C  no asymptotics yet for interaction potentials

C  H.1
      ELSEIF (ISW.EQ.1) THEN
        CHR=' a0 '
        I0=0
        MODCLF(IR)=MODCLF(IR)+10
        IFLG=1
C  DEFAULT CROSS-SECTION: 8TH-ORDER POLYNOMIAL OF LN(SIGMA) VS LN(E)
        IFTFLG(IR,1)=0
c  (laboratory)  energy range, asymptotics
        CH1L='al0'
        CH1R='ar0'
        C1L = 'ELABMIN'
        C1R = 'ELABMAX'

C  H.2
      ELSEIF (ISW.EQ.2) THEN
        CHR=' b0 '
        CH1L='bl0'
        CH1R='br0'
        I0=1
        MODCLF(IR)=MODCLF(IR)+100
        IFLG=2
C  DEFAULT RATE COEFFICIENT: 8TH-ORDER POLYNOMIAL OF LN(<SIGMA V>) VS LN(T), FOR E0=0.
        IFTFLG(IR,2)=0
c  temperature range, asymptotics
        C1L = 'T1MIN'
        C1R = 'T1MAX'

C  H.3
      ELSEIF (ISW.EQ.3) THEN
        CHR=' c0 '
        CH1L='cl0'
        CH1R='cr0'
        CH2L='cb0'
        CH2R='ct0'
        MODCLF(IR)=MODCLF(IR)+200
        I0=1
        IFLG=2
C  DEFAULT RATE COEFFICIENT: DOUBLE POLYNOMIAL OF LN(<SIGMA V>) VS LN(T) AND LN(E0)
        IFTFLG(IR,2)=0
c  temperature range, asymptotics
c  beam energy range, asymptotics
        C1L = 'T1MIN'
        C1R = 'T1MAX'
        C2L = 'E2MIN'
        C2R = 'E2MAX'
C  H.4
      ELSEIF (ISW.EQ.4) THEN
        CHR=' d0 '
        CH1L='dl0'
        CH1R='dr0'
        CH2L='db0'
        CH2R='dt0'
        MODCLF(IR)=MODCLF(IR)+300
        I0=1
        IFLG=2
C  DEFAULT RATE COEFFICIENT: DOUBLE POLYNOMIAL OF LN(<SIGMA V>) VS LN(T) AND LN(NE)
        IFTFLG(IR,2)=0
c  temperature range, asymptotics
c  beam energy range, asymptotics
        C1L = 'T1MIN'
        C1R = 'T1MAX'
        C2L = 'N2MIN'
        C2R = 'N2MAX'
C  H.5
      ELSEIF (ISW.EQ.5) THEN
        CHR=' e0 '
        CH1L='el0'
        CH1R='er0'
        I0=1
        MODCLF(IR)=MODCLF(IR)+1000
        IFLG=3
C  MOMENTUM-WEIGHTED RATE COEFFICIENT
        IFTFLG(IR,3)=0
        C1L = 'T1MIN'
        C1R = 'T1MAX'

C  H.6
      ELSEIF (ISW.EQ.6) THEN
        CHR=' f0 '
        CH1L='fl0'
        CH1R='fr0'
        CH2L='fb0'
        CH2R='ft0'
        MODCLF(IR)=MODCLF(IR)+2000
        I0=1
        IFLG=3
C  MOMENTUM-WEIGHTED RATE COEFFICIENT
        IFTFLG(IR,3)=0
        C1L = 'T1MIN'
        C1R = 'T1MAX'
        C2L = 'E2MIN'
        C2R = 'E2MAX'
C  H.7
      ELSEIF (ISW.EQ.7) THEN
        CHR=' g0 '
        CH1L='gl0'
        CH1R='gr0'
        CH2L='gb0'
        CH2R='gt0'
        MODCLF(IR)=MODCLF(IR)+3000
        I0=1
        IFLG=3
C  MOMENTUM-WEIGHTED RATE COEFFICIENT
        IFTFLG(IR,3)=0
        C1L = 'T1MIN'
        C1R = 'T1MAX'
        C2L = 'N2MIN'
        C2R = 'N2MAX'
C  H.8
      ELSEIF (ISW.EQ.8) THEN
        CHR=' h0 '
        CH1L='hl0'
        CH1R='hr0'
        I0=1
        MODCLF(IR)=MODCLF(IR)+10000
        IFLG=4
C  ENERGY-WEIGHTED RATE COEFFICIENT
        IFTFLG(IR,4)=0
        C1L = 'T1MIN'
        C1R = 'T1MAX'

C  H.9
      ELSEIF (ISW.EQ.9) THEN
        CHR=' i0 '
        CH1L='il0'
        CH1R='ir0'
        CH2L='ib0'
        CH2R='it0'
        MODCLF(IR)=MODCLF(IR)+20000
        I0=1
        IFLG=4
C  ENERGY-WEIGHTED RATE COEFFICIENT
        IFTFLG(IR,4)=0
        C1L = 'T1MIN'
        C1R = 'T1MAX'
        C2L = 'E2MIN'
        C2R = 'E2MAX'
C  H.10
      ELSEIF (ISW.EQ.10) THEN
        CHR=' j0 '
        CH1L='jl0'
        CH1R='jr0'
        CH2L='jb0'
        CH2R='jt0'
        MODCLF(IR)=MODCLF(IR)+30000
        I0=1
        IFLG=4
C  ENERGY-WEIGHTED RATE COEFFICIENT
        IFTFLG(IR,4)=0
        C1L = 'T1MIN'
        C1R = 'T1MAX'
        C2L = 'N2MIN'
        C2R = 'N2MAX'
C  H.11
      ELSEIF (ISW.EQ.11) THEN
        CHR=' k0 '
        CH1L='kl0'
        CH1R='kr0'
        I0=1
        IFLG=5
        IFTFLG(IR,IFLG)=0
        C1L = 'P1MIN'
        C1R = 'P1MAX'
C  H.12
      ELSEIF (ISW.EQ.12) THEN
        CHR=' l0 '
        CH1L='ll0'
        CH1R='lr0'
        CH2L='lb0'
        CH2R='lt0'
        I0=1
        IFLG=5
        IFTFLG(IR,IFLG)=0
        C1L = 'P1MIN'
        C1R = 'P1MAX'
        C2L = 'P2MIN'
        C2R = 'P2MAX'
      ENDIF


      IF (INDEX(FILNAM,'CR').NE.0) THEN
        CALL EIRENE_READ_COLRAD (IR,REAC,ISW,IZ1,
     .                           IROW_ESC,ICOL_ESC,POP_ESC)
c  close unit=29+ifoff:   done in READ_COLRAD.f
        RETURN
      END IF

      IF (INDEX(FILNAM,'TAB2D').NE.0 .OR.
     .    INDEX(FILNAM,'ADAS') .NE.0) THEN
        CALL EIRENE_READ_TAB2D (IR,REAC,ISW,IZ1)
c  close unit=29+ifoff:   done in READ_TAB2D.f
        RETURN
      END IF

      IF (INDEX(FILNAM,'CONST').NE.0) THEN
        IND=INDEX(REACSTR,'FT')
        IF (IND /= 0) THEN
c  read parameter for type of fitting expression from data file
          READ (REACSTR(IND+2:),*) IFTFLG(IR,IFLG)
        END IF

        IF (MOD(IFTFLG(IR,IFLG),100) == 10) THEN
C
C  SET A REACTION CONSTANT
C  READ ONLY ONE FIT COEFFICIENT FROM INPUT FILE 'iunin'
          READ (IUNIN,6664) CREACD(1,1)
        ELSE
C
C  READ 9 FIT COEFFICIENTS (2 CARDS) FROM INPUT FILE 'iunin'
          READ (IUNIN,6664) (CREACD(IC,1),IC=1,9)

        END IF
        CALL EIRENE_SET_REACTION_DATA    ! this routine sets only "POLY" data
     .          (IR,ISW,IFTFLG(IR,IFLG),CREACD,IUNOUT,.FALSE.)
c  no optional extrapolation flags here
        RETURN
      ENDIF
C
C  READ FROM DATA FILE, stream 29
C
C  already ruled out here (done at this point):
C  FILNAM= "CR...", "CONST", "ADAS",  "PHOTON"
C  in all these cases: already returned to calling program
C
C......................................................................
C  AT THIS POINT: FILNAM= AMJUEL, HYDHEL, H2VIBR, METHAN, i.e. single or double polynomial fits

CC  now identify proper dataset within file FILNAM

  100 READ (29+ifoff,'(A80)',END=990) ZEILE
      IF (INDEX(ZEILE,'##BEGIN DATA HERE##').EQ.0) GOTO 100

      LAST_TEX=REPEAT(' ',80)
    1 READ (29+ifoff,'(A80)',END=990) ZEILE
!ITER IF (INDEX(ZEILE,H123).EQ.0) GOTO 1
!PB      IF (INDEX(ZEILE,H123).EQ.0 .or.
!PB  .      INDEX(ZEILE,'section').EQ.0) GOTO 1   !  infinite loop possible !
      IF (INDEX(ZEILE,H123).EQ.0) THEN
        IF (INDEX(ZEILE,BACK).NE.0) LAST_TEX=ZEILE
        GOTO 1
      ELSE
        IF ((INDEX(LAST_TEX,SECTION) .EQ. 0) .AND.
     .      (INDEX(ZEILE,SECTION) .EQ. 0)) GOTO 1
      END IF
C
    2 READ (29+ifoff,'(A80)',END=990) ZEILE
!ITER IF (INDEX(ZEILE,'H.').NE.0) GOTO 990
      IF (INDEX(ZEILE,'H.').NE.0 .and.
     .    INDEX(ZEILE,'section').NE.0) GOTO 990
      IF (INDEX(ZEILE,'Reaction ').EQ.0.or.
     .    INDEX(ZEILE,REACSTR(1:ireac)).EQ.0) GOTO 2  ! infinite loop possible  !

C
C  SINGLE PARAM. FIT, ISW=0,1,2,5,8,11
C
      IF (ISW.EQ.0.OR.ISW.EQ.1.OR.ISW.EQ.2.OR.ISW.EQ.5.OR.ISW.EQ.8.OR.
     .    ISW.EQ.11) THEN

    3   READ (29+ifoff,'(A80)',END=990) ZEILE
        INDFF=INDEX(ZEILE,FITFLAG)
        IF (INDEX(ZEILE,CHR)+INDFF.EQ.0) GOTO 3
c  input line found which either contains fit-flag, or the reaction identifier a0,b0,...k0
        IF (INDFF > 0) THEN  ! OTHERWISE: use DEFAULT FOR FIT-FLAG: iftflg = 0
c  read parameter for type of fitting expression from data file
          READ (ZEILE((INDFF+8):80),*) IFTFLG(IR,IFLG)
          GOTO 3
        ENDIF
C  read only one constant:  (FIT-FLAG = 10, 110, ....)
        IF (MOD(IFTFLG(IR,IFLG),100) == 10) THEN
          IND=INDEX(ZEILE,CHR(2:2))
          READ (ZEILE((IND+2):80),'(E20.12)') CREACD(1,1)
        ELSE
C  READ 9 FIT COEFFICIENTS, SEPARATED BY 'CHR'  FIXED FORMAT E20.12
C  THREE LINES WITH THREE DATA PER LINE
          DO 9 J=0,2
            IND=0
            DO 4 I=1,3
              IND=IND+INDEX(ZEILE((IND+1):80),CHR(2:2))
              READ (ZEILE((IND+2):80),'(E20.12)') CREACD(J*3+I,1)
    4       CONTINUE
            READ (29+ifoff,'(A80)',END=990) ZEILE
    9     CONTINUE
        END IF
C

C  SINGLE PARAMETER POLYNOMIAL FITS: DONE

C   AT THIS POINT WE HAVE STORED FOR REACTION ir, DATA TYPE iflg: 0,...,5
C   IFTFLG(IR,iflg)   (DEFAUT:   =0)
C   9 FIT COEFFICIENTS ON INTERMEDIATE ARRAY CREACD(1...9,1)
C   AND POSSIBLY (SOME OF) THE EXTRAPOLATION PARAMETERS RCMIN,RCMAX, FP(1:6)
C


        GOTO 1000

C  TWO PARAM. FIT, ISW=3,4,6,7,9,10,12
      ELSEIF (ISW.EQ.3.OR.ISW.EQ.4.OR.ISW.EQ.6.OR.ISW.EQ.7.OR.
     .        ISW.EQ.9.OR.ISW.EQ.10.OR.ISW.EQ.12) THEN
C READ 3 BLOCKS "J" OF DATA. each block contains 9 LINES, 3 numbers per line, i.e. 3 sub blocks
        DO 11 J=0,2
   16     READ (29+ifoff,'(A80)',END=990) ZEILE
C  SEARCH FOR STRING 'fit-flag'  or 'Index'
          INDFF=INDEX(ZEILE,'fit-flag')
          IF (INDEX(ZEILE,'Index')+INDFF.EQ.0) GOTO 16
          IF (INDFF > 0) THEN
c  read parameter for type of fitting expression from data file
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
C   READ 9 LINES, THREE DATA EACH LINE, UNFORMATTED I.E. READ 3 SUB-BLOCKS K,K+1,K+2
              READ (29+ifoff,*) IH,(CREACD(I,K),K=J*3+1,J*3+3)
   17       CONTINUE
c    first  index I: I-th block, vertical, Temp. dependence
c    second index K:  from sub block to sub-block (horizontal), ne, eb dependence.
c  d.h. erster sub block entspricht ln(ne/1e8))=0, oder ne=1e8, corona rate vs. T
          END IF
   11   CONTINUE
        READ (29+ifoff,'(A80)',END=990) ZEILE
C   AT THIS POINT WE HAVE STORED FOR REACTION ir: , DATA TYPE iflg: 0,...,5
C   IFTFLG(IR,iflg)   (DEFAUT:   =0)
C   81 FIT COEFFICIENTS ON INTERMEDIATE ARRAY CREACD(1...9,1...9)

C
C  DOUBLE PARAMETER POLYNOMIAL FITS: DONE

        GOTO 1000
C
      ENDIF

 1000 CONTINUE

C  NEXT: READ ASYMPTOTICS INFORMATION FROM ATOMIC DATA FILE
C        HYDHEL, AMJUEL, H2VIBR, METHANE.

C FOR 1D OR 2D DATASETS. 4 BOUNDARIES, LEFT1, RIGHT1, LEFT2, RIGHT2.
C FOR 1D: ONLY "LEFT1" AND "RIGHT1" ARE USED

      IF (ISW.EQ.0) GOTO 2000    ! NO ASYMPTOTICS FOR POTENTIALS

c  INDICATE, IF EXTRAPOLATION INFORMATION (r1mn,.., if1mn,...) IS FOUND ON DATA FILE
c  DEFAULT:  NO DATA FOUND
      LGR1MIN=.FALSE.
      LGR1MAX=.FALSE.
      LGR2MIN=.FALSE.
      LGR2MAX=.FALSE.

c  INDICATE, IF EXTRAPOLATION COEFFICIENTS (fp1l,fp1r,fp2l,fp2r) ARE FOUND ON DATA FILE
c  DEFAULT:  NO DATA FOUND
      LGC1MIN=.FALSE.
      LGC1MAX=.FALSE.
      LGC2MIN=.FALSE.
      LGC2MAX=.FALSE.

C  FLAG FOR CHOICE OF EXTRAPOLATION OPTION:
C  DEFAULT ASYMPTOTIC EXPRESSION  (...=0): NO ASYMPTOTICS
C  DEFAULT ASYMPTOTIC EXPRESSION  (...=1): SET TO ZERO BEYOND LAST VALID POINT
C  DEFAULT ASYMPTOTIC EXPRESSION  (...=4): TAKE LAST VALID POINT AT r1mn,r1mx,....
C  DEFAULT ASYMPTOTIC EXPRESSION  (...=5): 2ND-ORDER POLYNOMIAL BEYOND LAST VALID POINT
C                                          (OLD DEFAULT FOR CROSS-SECTIONS)
c  AND EXTRAPOLATE CONSTANT FROM THERE
      IF1MN = 0
      IF1MX = 0
      IF2MN = 0
      IF2MX = 0

C  FURTHER PARAMETERS, NOT RELATED TO ASYMPTOTICS
      RTMAX = 0._DP
      ERTMAX = -HUGE(1._DP)

cdr Try to read asymptotics, unless already read: jfeximn,jfeximx /= 0
cdr Even if jfeximn,jfeximx /= 0, read anyway, but later: do not use

c  data for one "reaction IR" are located between "\begin" and "\end"
!     BEND="\end" stops looking for asymptotics
      DO WHILE(INDEX(ZEILE,BEND) == 0)

        ULINE = ZEILE
        CALL EIRENE_UPPERCASE(ULINE)
c
c  at this point we have found a card ZEILE which contains
c  one of the extrapolation parameter identifiers al0,ar0,....,k0l,k0r
c  that corresponds to the H.1, ....H.12 type of data IR.
c  next: read up to three fit coefficients FP.L OR FP.R.
c  CHR(2:2) is set to either character a,b,c,....,or k
c
        IF (INDEX(ZEILE,CH1L) /= 0) THEN
          CALL EIRENE_READ_COEFFS (ZEILE,CHR(2:2),FP1L)
          LGC1MIN = .TRUE.
        ENDIF
        IF (INDEX(ZEILE,CH1R) /= 0) THEN
          CALL EIRENE_READ_COEFFS (ZEILE,CHR(2:2),FP1R)
          LGC1MAX = .TRUE.
        END IF

        IF (INDEX(ZEILE,CH2L) /= 0) THEN
          CALL EIRENE_READ_COEFFS (ZEILE,CHR(2:2),FP2L)
          LGC2MIN = .TRUE.
        END IF
        IF (INDEX(ZEILE,CH2R) /= 0) THEN
          CALL EIRENE_READ_COEFFS (ZEILE,CHR(2:2),FP2R)
          LGC2MAX = .TRUE.
        END IF
c
c  currently foreseen asymptotic data identifiers in data files:
c  c1l,c2l,c1r,c2r:  ELABMIN, ELABMAX,
c                    T1MIN,T1MAX,E2MIN,E2MAX, N2MIN, N2MAX,
c                    P1MIN,P1MAX,P2MIN,P2MAX
        IF (INDEX(ULINE,TRIM(C1L)) /= 0) THEN
          CALL EIRENE_READ_RANGE (ULINE,C1L,'EXT-FLG',R1MN,IF1MN)
          LGR1MIN = .TRUE.
c  old default: 2nd-order polynomial beyond valid range, with coefs. FPL1
          IF (IF1MN == 0 .AND. LGC1MIN) IF1MN = 5
c  default extrapolation from r1mn (by constant continuation) will be: jfex1mn=4
          IF (IF1MN == 0) IF1MN = 4
        END IF
        IF (INDEX(ULINE,TRIM(C1R)) /= 0) THEN
          CALL EIRENE_READ_RANGE (ULINE,C1R,'EXT-FLG',R1MX,IF1MX)
          LGR1MAX = .TRUE.
c  old default: 2nd-order polynomial beyond valid range, with coefs. FPR1
          IF (IF1MX == 0 .AND. LGC1MAX) IF1MX = 5
c  default extrapolation from r1mx (by constant continuation) will be: jfex1mx=4
          IF (IF1MX == 0) IF1MX = 4
        END IF
        IF (INDEX(ULINE,TRIM(C2L)) /= 0) THEN
          CALL EIRENE_READ_RANGE (ULINE,C2L,'EXT-FLG',R2MN,IF2MN)
          LGR2MIN = .TRUE.
c  old default: 2nd-order polynomial beyond valid range, with coefs. FPL2
          IF (IF2MN == 0 .AND. LGC2MIN) IF2MN = 5
c  default extrapolation from r2mn (by constant continuation) will be: jfex2mn=4
          IF (IF2MN == 0) IF2MN = 4
        END IF
        IF (INDEX(ULINE,TRIM(C2R)) /= 0) THEN
          CALL EIRENE_READ_RANGE (ULINE,C2R,'EXT-FLG',R2MX,IF2MX)
          LGR2MAX = .TRUE.
c  old default: 2nd-order polynomial beyond valid range, with coefs. FPr2
          IF (IF2MX == 0 .AND. LGC2MAX) IF2MX = 5
c  default extrapolation from r2mx (by constant continuation) will be: jfex2mx=4
          IF (IF2MX == 0) IF2MX = 4
        END IF
C
C  ...AND FURTHER REACTION PARAMETERS, NOT RELATED TO ASYMPTOTICS
C     ETH
C     RTMAX
C     ERTMAX
C
        IND = INDEX(ULINE,TRIM(CETH))
        IF (IND /= 0) THEN
          READ (ULINE(IND+3:80),*) ETH
        END IF
        IND = INDEX(ULINE,TRIM(CMR))
        IF (IND /= 0) THEN
          INDG = INDEX(ULINE,'=')
          READ (ULINE(INDG+1:80),*) RTMAX
          IND = INDEX(ULINE,TRIM(CEMR))
          IF (IND /= 0) THEN
            INDG = IND + INDEX(ULINE(IND:80),'=')
            READ (ULINE(INDG+1:),*) ERTMAX
          END IF

        END IF

        READ (29+ifoff,'(A80)',END=990) ZEILE
      END DO

c  unless asymptotics are already explicitly defined in input block 4a
c  put asymptotics information into proper (intermediate) data structure:
c  flags:      jfex1mn,jfex1mx,jfex2mn,jfex2mx
c  boundaries: rc1min,rc1max,rc2min,rc2max
c  parameters: fp1(1:3),fp1(4:6),fp2(1:3),fp2(4:6)

      IF (JFEX1MN == 0) THEN
        IF (LGR1MIN .AND. .NOT. LGC1MIN.and.if1mn.ge.3.) THEN
          WRITE (IUNOUT,*) ' WARNING FROM SLREAC'
          WRITE (IUNOUT,*) ' REACTION ',IR, 'TYPE ',H123
          WRITE (IUNOUT,*) ' LOWER RANGE FOR 1ST PARAMETER OF FIT',
     .          ' SPECIFIED BUT',
     .          ' NO COEFFICIENTS FOR EXTRAPOLATION PROVIDED'
          CALL EIRENE_MASJ1R('IF1MN,R1MN      ',if1mn,r1mn)
          IF (IF1MN.EQ.4)
     .      WRITE (IUNOUT,*) 'I.E.: CONTINUATION AS CONSTANT'
          CALL EIRENE_LEER(1)
        ELSEIF (LGR1MIN) THEN
          WRITE (IUNOUT,*) 'ASYMPTOTICS FROM SLREAC'
          WRITE (IUNOUT,*) 'REACTION ',IR, 'TYPE ',H123
          WRITE (IUNOUT,*) 'LOWER RANGE FOR 1ST PARAMETER OF FIT'
          CALL EIRENE_MASJ1R('IF1MN,R1MN      ',if1mn,r1mn)
          if (if1mn.ge.3)
     .      CALL EIRENE_MASRR1('PARAMETERS ',fp1l,3,3)
          CALL EIRENE_LEER(1)
        END IF

        IF (LGR1MIN) RC1MIN = LOG(R1MN)
        IF (LGC1MIN) FP1(1:3) = FP1L
        JFEX1MN = IF1MN
        IF (LGR1MIN .AND. LGC1MIN .AND. (JFEX1MN == 0))
! DEFAULT EXTRAPOLATION=EXP(FP(1)+FP(2)*PARM+FP(3)*PARM**2), 2ND ORDER ON LOG SCALE
     .          JFEX1MN = 5
      END IF

      IF (JFEX1MX == 0) THEN
        IF (LGR1MAX .AND. .NOT. LGC1MAX.and.if1mx.ge.3.) THEN
          WRITE (IUNOUT,*) ' WARNING FROM SLREAC'
          WRITE (IUNOUT,*) ' REACTION ',IR, ' TYPE ',H123
          WRITE (IUNOUT,*) ' UPPER RANGE FOR 1ST PARAMETER OF FIT',
     .          ' SPECIFIED BUT',
     .          ' NO COEFFICIENTS FOR EXTRAPOLATION PROVIDED'
          CALL EIRENE_MASJ1R('IF1MX,R1MX      ',if1mx,r1mx)
          IF (IF1MX.EQ.4)
     .      WRITE (IUNOUT,*) 'I.E.: CONTINUATION AS CONSTANT'
          CALL EIRENE_LEER(1)
        ELSEIF (LGR1MAX) THEN
          WRITE (IUNOUT,*) 'ASYMPTOTICS FROM SLREAC'
          WRITE (IUNOUT,*) 'REACTION ',IR, ' TYPE ',H123
          WRITE (IUNOUT,*) 'UPPER RANGE FOR 1ST PARAMETER OF FIT'
          CALL EIRENE_MASJ1R('IF1MX,R1MX      ',if1mx,r1mx)
          if (if1mx.ge.3)
     .      CALL EIRENE_MASRR1('PARAMETERS ',fp1r,3,3)
          CALL EIRENE_LEER(1)
        END IF

        IF (LGR1MAX) RC1MAX = LOG(R1MX)
        IF (LGC1MAX) FP1(4:6) = FP1R
        JFEX1MX = IF1MX
        IF (LGR1MAX .AND. LGC1MAX .AND. (JFEX1MX == 0))
! DEFAULT EXTRAPOLATION=EXP(FP(1)+FP(2)*PARM+FP(3)*PARM**2), 2ND ORDER ON LOG SCALE
     .          JFEX1MX = 5
      END IF

      IF (JFEX2MN == 0) THEN
        IF (LGR2MIN .AND. .NOT. LGC2MIN.and.if2mn.ge.3.) THEN
          WRITE (IUNOUT,*) ' WARNING FROM SLREAC'
          WRITE (IUNOUT,*) ' REACTION ',IR, ' TYPE ',H123
          WRITE (IUNOUT,*) ' LOWER RANGE FOR 2ND PARAMETER OF FIT',
     .          ' SPECIFIED BUT',
     .          ' NO COEFFICIENTS FOR EXTRAPOLATION PROVIDED'
          CALL EIRENE_MASJ1R('IF2MN,R2MN      ',if2mn,r2mn)
          IF (IF2MN.EQ.4)
     .      WRITE (IUNOUT,*) 'I.E.: CONTINUATION AS CONSTANT'
          CALL EIRENE_LEER(1)
        ELSEIF (LGR2MIN) THEN
          WRITE (IUNOUT,*) 'ASYMPTOTICS FROM SLREAC'
          WRITE (IUNOUT,*) 'REACTION ',IR, ' TYPE ',H123
          WRITE (IUNOUT,*) 'LOWER RANGE FOR 2ND PARAMETER OF FIT'
          CALL EIRENE_MASJ1R('IF2MN,R2MN      ',if2mn,r2mn)
          if (if2mn.ge.3)
     .      CALL EIRENE_MASRR1('PARAMETERS ',fp2l,3,3)
          CALL EIRENE_LEER(1)
        END IF
        IF (LGR2MIN) RC2MIN = LOG(R2MN)
        IF (LGC2MIN) FP2(1:3) = FP2L
        JFEX2MN = IF2MN
        IF (LGR2MIN .AND. LGC2MIN .AND. (JFEX2MN == 0))
! DEFAULT EXTRAPOLATION=EXP(FP(1)+FP(2)*PARM+FP(3)*PARM**2), 2ND ORDER ON LOG SCALE
     .         JFEX2MN = 5
      END IF

      IF (JFEX2MX == 0) THEN
        IF (LGR2MAX .AND. .NOT. LGC2MAX.and.if2mx.ge.3.) THEN
          WRITE (IUNOUT,*) ' WARNING FROM SLREAC'
          WRITE (IUNOUT,*) ' REACTION ',IR, ' TYPE ',H123
          WRITE (IUNOUT,*) ' UPPER RANGE FOR 2ND PARAMETER OF FIT',
     .          ' SPECIFIED BUT',
     .          ' NO COEFFICIENTS FOR EXTRAPOLATION PROVIDED'
          CALL EIRENE_MASJ1R('IF2MX,R2MX      ',if2mx,r2mx)
          IF (IF2MX.EQ.4)
     .      WRITE (IUNOUT,*) 'I.E.: CONTINUATION AS CONSTANT'
          CALL EIRENE_LEER(1)
        ELSEIF (LGR2MAX) THEN
          WRITE (IUNOUT,*) 'ASYMPTOTICS FROM SLREAC'
          WRITE (IUNOUT,*) 'REACTION ',IR, ' TYPE ',H123
          WRITE (IUNOUT,*) 'UPPER RANGE FOR 2ND PARAMETER OF FIT'
          CALL EIRENE_MASJ1R('IF2MX,R2MX      ',if2mx,r2mx)
          if (if2mx.ge.3)
     .      CALL EIRENE_MASRR1('PARAMETERS ',fp2r,3,3)
          CALL EIRENE_LEER(1)
        END IF
        IF (LGR2MAX) RC2MAX = LOG(R2MX)
        IF (LGC2MAX) FP2(4:6) = FP2R
        JFEX2MX = IF2MX
        IF (LGR2MAX .AND. LGC2MAX .AND. (JFEX2MX == 0))
! DEFAULT EXTRAPOLATION=EXP(FP(1)+FP(2)*PARM+FP(3)*PARM**2), 2ND ORDER ON LOG SCALE
     .        JFEX2MX = 5
      END IF
C
 2000 CONTINUE

      CALL
     .  EIRENE_SET_REACTION_DATA   ! this routine sets only "POLY" data
     .            (IR,ISW,IFTFLG(IR,IFLG),CREACD,IUNOUT,.TRUE.,
c  from here on: optional input to SET_REACTION_DATA
     .                       RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .                       RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .                       RTMAX,ERTMAX,ETH)
C
      CLOSE (UNIT=29+ifoff)
C
      RETURN
C
  990 WRITE (iunout,*) ' NO DATA FOUND FOR REACTION ',H123,' ',REAC,
     .                 ' IN DATA SET ',FILNAM
      WRITE (iunout,*) ' IR,MODCLF(IR) ',IR,MODCLF(IR)
      CLOSE (UNIT=29+ifoff)
      CALL EIRENE_EXIT_OWN(1)
  991 WRITE (iunout,*) ' INVALID CONSTANT IN SLREAC. CONST= ',CONST
      WRITE (iunout,*) ' CHECK "REACTION CARDS" FOR REACTION NO. ',IR
      CLOSE (UNIT=29+ifoff)
      CALL EIRENE_EXIT_OWN(1)
 6664 FORMAT (6E12.4)

      CONTAINS

      SUBROUTINE EIRENE_READ_COEFFS (ZEILE,CH,FP)
c  ZEILE is in upper case, and CH is found.
c  CH is: a,b,c,.....k,l, (depending on H.1, H.2, ...H.12)
c  read up to three parameters FP(i), i=1,3 formatted: '(E20.12)'
c  to be used for extrapolation
      CHARACTER(80), INTENT(IN) :: ZEILE
      CHARACTER(1), INTENT(IN) :: CH
      REAL(DP), INTENT(OUT) :: FP(3)
      INTEGER :: I, IND, INC

      IND=0
      DO I=1,3
        INC=INDEX(ZEILE((IND+1):80),CH)
        IF (INC.GT.0) THEN
           IND=IND+INDEX(ZEILE((IND+1):80),CH)
           READ (ZEILE((IND+3):80),'(E20.12)') FP(I)
        ENDIF
      END DO

      END SUBROUTINE EIRENE_READ_COEFFS


      SUBROUTINE EIRENE_READ_RANGE (ZEILE,KEY1,KEY2,RNG,IFX)
c  called from slreac, after the original fit coefficients for reaction IR
c  are read.
c  At this point an extrapolation card ZEILE belonging to this reaction IR
c  has already been found.
c
c  This routine:
c  Reads validity range from atomic data file:
c  Search in ZEILE for key1, key2 and return: RNG, IFX

c  key1:  'ELABMIN', 'ELABMAX',   'T1MIN','T1MAX','E2MIN','E2MAX',
C         'N2MIN','N2MAX','P1MIN','P1MAX','P1MIN','P1MAX'=,
C                     read RNG (unformatted, real)
c  key2:  'EXT-FLG'= ,read IFX (unformatted, integer)

      CHARACTER(80), INTENT(IN) :: ZEILE
      CHARACTER(7), INTENT(IN) :: KEY1, KEY2
      REAL(DP), INTENT(OUT) :: RNG
      INTEGER, INTENT(OUT) :: IFX
      INTEGER :: IND1, IND2, INDE, INDP, INDX, INDA, INDG
      CHARACTER(20) :: FORM

      IND1 = INDEX(ZEILE,TRIM(KEY1))
      IND2 = INDEX(ZEILE,TRIM(KEY2))

      RNG = 0._DP
      IFX = 0

      IF (IND1 > 0) THEN
        INDG = INDEX(ZEILE,'=')
        INDE = INDG + VERIFY(ZEILE(INDG+1:),'+-0123456789DEed. ') - 1
        INDP = SCAN(ZEILE(INDG+1:),'.')
        INDX = SCAN(ZEILE(INDG+1:),'EDed')
        INDA = SCAN(ZEILE(INDG+1:),'+-0123456789.')
        FORM=REPEAT(' ',20)
        WRITE (FORM,'(A2,I0,A1,I0,A1)')
     .         '(E',INDX+3-INDA+1,'.',INDX-INDP-1,')'
        READ (ZEILE(INDG+INDA:INDE),FORM) RNG
      END IF

      IF (IND2 > 0) READ (ZEILE(IND2+7:),*) IFX

      END SUBROUTINE EIRENE_READ_RANGE

      END
