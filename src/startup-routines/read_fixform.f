
      SUBROUTINE EIRENE_READ_FIXFORM (NL, SAREA_SAVE, IERROR)
C
C   READ INPUT DATA AND SET DEFAULT VALUES
c   IN CASE IITER.GT.1 OR ITIMV.GT.1 : SKIP READING NEW INPUT FROM IUNIN.
C                                      ONLY INPUT DATA PROCESSING (STATEMENT 4000 FF)
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CADGEO
      USE EIRMOD_CCONA
      USE EIRMOD_CGRPTL
      USE EIRMOD_CLOGAU
      USE EIRMOD_CPL3D
      USE EIRMOD_CPLOT
      USE EIRMOD_CINIT
      USE EIRMOD_COMSIG
      USE EIRMOD_CREF
      USE EIRMOD_CPOLYG
      USE EIRMOD_CGRID
      USE EIRMOD_CZT1
      USE EIRMOD_CTRCEI
      USE EIRMOD_CCOUPL
      USE EIRMOD_CGEOM
      USE EIRMOD_CSDVI
      USE EIRMOD_COMPRT
      USE EIRMOD_CPES, ONLY: NPRS,NLIDENT
      USE EIRMOD_COMNNL
      USE EIRMOD_COMSOU
      USE EIRMOD_CSTEP
      USE EIRMOD_COMSPL
      USE EIRMOD_CTEXT
      USE EIRMOD_CLGIN
      USE EIRMOD_COMXS
      USE EIRMOD_CSPEI
      USE EIRMOD_CTRIG
      USE EIRMOD_CTETRA
      USE EIRMOD_CESTIM
      USE EIRMOD_CUPD
      USE EIRMOD_PHOTON
      USE EIRMOD_MPI
      USE EIRMOD_SECOND_OWN, ONLY: EIRENE_SECOND_OWN
      USE EIRMOD_TIMEA, ONLY: EIRENE_TIMEA0
      USE EIRMOD_PROFILES, ONLY: EIRENE_PROFR
      USE EIRMOD_JSON

      IMPLICIT NONE

C
      TYPE(TEMPERATURE),POINTER :: TEMPCUR
      TYPE(DENSITY),POINTER :: DENCUR
      TYPE(VELOCITY),POINTER :: VELCUR
      TYPE(VOLUMEP),POINTER :: VOLCUR
C
      TYPE(TSURFACE), POINTER :: SURFCUR, SURFCUR2
      TYPE(REFMODEL), POINTER :: REFCUR, REFLAST
      TYPE(SPEC_REF), POINTER :: SPR

      TYPE REFFILE
        CHARACTER(500) :: RFILE
        TYPE(REFFILE), POINTER :: NEXT
      END TYPE REFFILE

      TYPE(REFFILE), POINTER :: REFFILES, CURFILE

      TYPE(EIRENE_SPECTRUM), POINTER :: ESPEC, SSPEC
      TYPE(TCONTRIB) :: CNT
cdr  optional arguments for FILNAM=CONST options
      INTERFACE
        SUBROUTINE EIRENE_SLREAC (IR,FILNAM,H123,REAC,CRC,
     .             RC1MIN, RC1MAX, FP1, JFEX1MN, JFEX1MX,
     .             RC2MIN, RC2MAX, FP2, JFEX2MN, JFEX2MX,
     .             ELNAME, IZ1, 
     .             IROW_ESC, ICOL_ESC, POP_ESC,  ! for internal CR models, line emission etc..
     .             IFTFL, NCOEF, COEF)  ! for filnam=const
        USE EIRMOD_PRECISION
        INTEGER,      INTENT(IN) :: IR, IZ1
        INTEGER,      INTENT(IN), OPTIONAL :: IROW_ESC, ICOL_ESC,
     .                                        IFTFL, NCOEF
        REAL(DP),     INTENT(IN), OPTIONAL :: POP_ESC
        REAL(DP),     INTENT(IN), OPTIONAL :: COEF(9)      
        CHARACTER(8), INTENT(IN) :: FILNAM
        CHARACTER(4), INTENT(IN) :: H123
        CHARACTER(LEN=*), INTENT(IN) :: REAC
        CHARACTER(2), INTENT(IN) :: ELNAME
        CHARACTER(3), INTENT(IN) :: CRC
        INTEGER,  INTENT(IN OUT) :: JFEX1MN, JFEX1MX,JFEX2MN, JFEX2MX
        REAL(DP), INTENT(IN OUT) :: RC1MIN, RC1MAX, FP1(6),
     .                              RC2MIN, RC2MAX, FP2(6)
        END SUBROUTINE EIRENE_SLREAC
      END INTERFACE


      TYPE(TRANSFORM), POINTER :: TRF
C
      INTEGER, INTENT(IN) :: NL
      REAL(DP), INTENT(INOUT) :: SAREA_SAVE(NL)
      INTEGER, INTENT(INOUT) :: IERROR
      REAL(DP) :: AFF(3,3), AFFI(3,3), FP1(6), FP2(6), COEF(9)
      REAL(DP) :: RP1, SA, SI, THMAX, SM, SPP, DTIMVO, SAVE, VOLTOT_TAL,
     .          SPH, ALR, ROTNRM, RPSDL, XSH, YSH, ZSH, REFNRM,
     .          DPP, R1MN, R1MX, R2MN, R2MX,
     .          SPCMN, SPCMX, SPC_SHIFT,
     .          SPCPLT_X, SPCPLT_Y, SPCPLT_SAME, SPCVX, SPCVY, SPCVZ,
     .          VNORM, ESCD2A, ESCD2M, ESCD2I, ESCD2PH, ESCD2P,
     .          RC1MIN, RC1MAX, RC2MIN, RC2MAX, POP_ESC

C  RUN TIME STATISTICS IN INITIALIZATION PHASE, WITHIN INPUT.F
cdr   REAL(DP) :: tpb1, tpb2, EIRENE_SECOND_OWN, timea
c
C  MULTIPLIER FOR BOTH CPU TIME NTCPU AND MAX NUMBER OF MC HISTORIES NPTS, ....
      REAL(DP) :: AMPTS

      INTEGER :: IHELP(NLIMPS), IPRSF(12), NUMTAL(12), IADTYP(0:4)
      INTEGER :: JP, KT, IM, II, IP, IA, NFLGS, NSPZS1, NSPZS2, IPH,
     .           NTLVF, NSRF, NTLS, I1000, NSP, NTL, ICHORI, IRAD,
     .           ILIMPS, ISS, ILA, IB, INC, NSOPT, IIN, IEN, NSMSTRA1,
     .           EIRENE_ILLZ, IZ, ISTREAM, ITALI, IBEND,
     .           JATM, JMOL, JION, JPHOT, JPLS, JTRJ, JSPZ,
     .           NO, IGO, IRPTA3, IRPTE2, IRPTE3, IH, IDIMP,
     .           JDUMMY, IRPTA1, IRPTA2, IRPTE1, NLJ, I1, I2, I3,
     .           NTIME0, IREAD, I, ISTS, J, IST,
     .           IPOS2, IPOS0, NM, K, JJ, IPOS1, NRGEN, INUM, NTLSF,
     .           L, INILGJ, INI, ICO, IS, NTLV, ID, IRE, NSC,
     .           IN, INELGJ, NPRCSF, MXL, NSPZV1, NSPZV2, NFLGV,
     .           IPRCSF, IR, MT, MP, NDUMM, NUMSEC, NDUMM1, NDUMM2,
     .           NRTAL1, NCOPIE, NFR, IPLN,
     .           NDUMM3, NDUMM4, NRE, IFLR,
     .           ISPSRF, ISPTYP, NSPS, NSPSA, IPTYP, IPSPZ,
     .           IANF, IEND, IDEFLT_SPUT, IDEFLT_SPEZ,
     .           IPLSTI, IPLSV, IFILE, ISRFCLL,
     .           IDIREC, ISTCHR,  ITOK, IER, ILOGS, IO,
     .           IUNIN_SAVE, NLOGIN, IFLG, IDUM,
     .           JFEX1MN, JFEX1MX, JFEX2MN, JFEX2MX,
     .           NB, NS, NA, ISTR, IRC, IOPT, ITAL,
     .           NRC, IADV, NUM_COMPO, NUM_CONTRIB, ICNT, IDMDL, IND,
     .           ILINE, JCOMP, KCONTR, IROW_ESC, ICOL_ESC,
     .           IFTFL, NCOEF, IRET, ITLV, ITLS, JL, IUSR
      INTEGER, SAVE :: NITER0, IUSROUT=0
      INTEGER, EXTERNAL :: EIRENE_IDEZ
      INTEGER, DIMENSION(1) :: ISTR_A
      LOGICAL :: LHELP(NLIMPS), NLSRON_SAVE(NSTRA)
      LOGICAL :: LRPS3D, LRPSCN, LHYDDEF, LINCL45, LMULTI, LRDMLTI
      LOGICAL, ALLOCATABLE :: LOGRDH(:)
      CHARACTER(10) :: CDATE, CTIME
      CHARACTER(12) :: CHR, HYDKIN_DEFAULT, CADAPT
      CHARACTER(420) :: ZEILE, ULINE
      CHARACTER(500) :: FILE, FILE45
      CHARACTER(8) :: FILNAM, varname, spcname
      CHARACTER(4) :: H123, CLAB
      CHARACTER(50) :: REAC
      CHARACTER(3) :: CRC
      CHARACTER(400) :: PATH, RFILNM
      CHARACTER(6) :: HANDLE
      CHARACTER(2) :: ELNAME
C
C  INITIALISE SOME DATA AND SET DEFAULTS
C
      IREAD=0
      IERROR=0

C  UNIT NUMBER FOR INPUT FILE: MUST BE DIFFERENT FROM: 5,8,10,11,12
C  13,14, AND 15
!pb IUNIN set in COMPRT
!pb   IUNIN=1+ifoff
      IUNIN_SAVE = IUNIN
      IUSROUT = 0
C
C  UNIT NUMBER FOR OUTPUT FILE: MUST BE DIFFERENT FROM: 5,8,10,11,12
C  13,14, AND 15 AND IUNIN
!pb IUNOUT has already been set in subroutine EIRENE
!pb   IUNOUT=6+ifoff
C

      NULLIFY(SURFLIST)
      NULLIFY(REFLIST)
C
C  SET DEFAULT DATA FOR BLOCK 13
C
      DTIMV=1.D30
      TIME0=0.
      NSNVI=0
      NTMSTP=1

      LMULTI = .FALSE.
C
C  READ TEXT DESCRIBING THE RUN, 100--199
C
      READ (IUNIN,'(A72)') TXTRUN
      WRITE (iunout,'(1X,A)') trim(TXTRUN)
      CALL EIRENE_LEER(1)
  109 READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1).EQ.'*') THEN
        IF (ZEILE(1:3).NE.'***') THEN
          WRITE (iunout,'(1X,A)') trim(ZEILE)
        CALL EIRENE_LEER(1)
!  STORE COMMENT LINES FOR OUTPUTING TO JSON FILE 
          call eirene_push_string_stack(cm_stack,trim(zeile)) 
        END IF  
        GOTO 109
      ELSE
        READ (ZEILE,6666) NPRLL,NMODE,NTCPU,NFILE,NITER0,NITER,
     .                    NTIME0,NTIME
      ENDIF
      CALL EIRENE_LEER(1)
      IITER=MAX0(1,NITER0)
      ITIMV=MAX0(1,NTIME0)

      READ (IUNIN,'(A72)') ZEILE
C......................................................................
C  IS THERE AN ADDITIONAL INPUT LINE FOR STORAGE OPTIMIZATION?

C  THESE DEFAULTS HAVE ALREADY BEEN SET IN FIND_PARAM
!     NOPTIM = ...
!     NOPTM1 = 1
!     NGEOM_USR = 0
!     NCOUP_INPUT = 1
!     NSMSTRA = 1    ADDITIONAL STORAGE FOR SUM OVER STRATA IS MADE AVAILABLE.
!     NSTORAM = 9    FULL STORAGE FOR ATOMIC DATA ARRAYS.
!                    SOME OF THAT CAN BE ELIMINATED BY "ON THE FLY A&M COMPUTATIONS"
!     NGSTAL = 0
!
!     NRTAL1 = 0     ONLY FOR FIND_PARAM, NOT USED ANY FURTHER

      IF ((INDEX(ZEILE,'F') + INDEX(ZEILE,'f') + INDEX(ZEILE,'T') +
     .     INDEX(ZEILE,'t')) == 0) THEN
C  YES
C THIS LINE:  NOPTIM, NSTORAM,.... is read in FIND_PARAM and modified in SET_PARMMOD(1)
!pb do not overwrite here
        READ (ZEILE,6666) IDUM,NOPTM1,NGEOM_USR,NCOUP_INPUT,
     .                    NSMSTRA1,NSTORAM,NGSTAL,NRTAL1
        READ (IUNIN,'(A72)') ZEILE
      END IF

cdr  repeat the same code as in find_param....?
C  NSTORAM IS REDEFINED, FINALLY EITHER =0  (A&M STORAGE SAVE MODE)
C                                    OR =9  (FULL A&M STORAGE MODE, =DEFAULT)
      NSTORAM = MIN(NSTORAM,9)
      IF (NSTORAM < 9) NSTORAM = 0
      NOPTM1 = MAX(NOPTM1,1)
C...........................................................................
c   done with this optional "storage save mode card"
C
* For gfortran: it does not accept empty field for logicals
      call fix_logical_input(zeile,16)
      READ (ZEILE,6665) NLSCL,NLTEST,NLANA,NLDRFT,NLCRR,
     .                  NLERG,NLIDENT,NLONE,NLMOVIE,NLDFST,
     .                  NLOLDRAN,NLCASCAD,NLOCTREE,NLWRMSH,NEXVS

C  OPTIONAL INPUT CARDS, FOR PATHWAYS AND NAME DEFINITIONS
C                        FOR EXTERNAL DATABASES: AMJUEL, HYDHEL,.....

C  IF DATA FILES ARE NOT ENTERED HERE:
c              THEN BY DEFAULT: THE EXTERNAL DATA FILES ARE EXPECTED IN THE SAME
C              DIRECTORY AS THE ONE FROM WHICH THE RUN IS STARTED.
C
cdr scan for optional CFILE lines: path to external database files
cdr Parse input card for PATH = DBFNAME

cdr  probably identical code as just done in find_param.f, repeated here.
cdr  begin

      READ (IUNIN,'(A420)') ZEILE
      IREAD=1
      I1 = INDEX(ZEILE,'CFILE')
      DO WHILE (I1 /= 0)
        I2 = VERIFY(ZEILE(I1+5:),' ') + I1 + 4
        I3 = SCAN(ZEILE(I2+1:),' ')
        HANDLE=REPEAT(' ',6)
        HANDLE(1:I3) = ZEILE(I2:I2+I3-1)
c   cfile card found. Is this one of the permitted external files?
        DO IFILE = 1,NDBNAMES
          IF (INDEX(DBHANDLE(IFILE),HANDLE) /= 0) EXIT
        END DO
        IF (IFILE <= NDBNAMES) THEN
c   yes, file type no 'ifile' as stored on dbhandle, in eirmod_cinit.
c   currently: 16 types of files are recognized
          IANF = I2+I3+VERIFY(ZEILE(I2+I3:),' ')-1
          IEND = IANF+SCAN(ZEILE(IANF+1:),' ')-1
          DBFNAME(IFILE)(1:IEND-IANF+1) = ZEILE(IANF:IEND)
          LDBREAD(IFILE) = .TRUE.
          WRITE (IUNOUT,*) 'PATH SET FOR FILE ',HANDLE
          WRITE (IUNOUT,*) 'PATH = ',ZEILE(IANF:IEND)
        ELSE
          WRITE (IUNOUT,*) ' WRONG NAME FOR DATABASE ENTERED '
          WRITE (IUNOUT,*) ' DATABASE DEFINITION FOR ',HANDLE,
     .                     ' IGNORED '
        END IF
        READ (IUNIN,'(A420)') ZEILE
        I1 = INDEX(ZEILE,'CFILE')
      END DO

C  READING OF OPTIONAL DATABASE NAME AND PATH-CARDS DONE.

C  CONTINUE WITH MANDATORY INPUT.  IREAD=1 AT THIS POINT

C  READING OF INPUT BLOCK 1 DONE

      NFILEN=EIRENE_IDEZ(NFILE,1,5)
      NFILEM=EIRENE_IDEZ(NFILE,2,5)
      NFILEL=EIRENE_IDEZ(NFILE,3,5)
      NFILEK=EIRENE_IDEZ(NFILE,4,5)
      NFILEJ=EIRENE_IDEZ(NFILE,5,5)
      CALL EIRENE_LEER(2)
      CALL EIRENE_MASAGE('*** 1. DATA FOR OPERATING MODE')
      CALL EIRENE_LEER(1)
      CALL EIRENE_MASAGE('       PARALLELISATION MODE:')
      SELECT CASE( NPRLL )
        CASE( -1 )
          CALL EIRENE_MASAGE('       MPI USER-DEFINED')
C       CASE( 0 )
C         Reserved for default, see below
        CASE( 1 )
          CALL EIRENE_MASAGE('       MPI PROPORTIONAL ALLOCATION')
        CASE DEFAULT
          CALL EIRENE_MASAGE('       MPI "EMBARRASSINGLY PARALLEL"')
          NPRLL = 0
      END SELECT
      WRITE (IUNOUT,*) '       NUMBER OF PROCESSES NPRS= ',NPRS
      CALL EIRENE_LEER(1)
      IF (NMODE.NE.0) THEN
        CALL EIRENE_MASAGE
     .  ('       EIRENE READS ADDITIONAL DATA FROM FILE')
        CALL EIRENE_MASAGE
     .  ('       INTERFACING ROUTINE INFCOP IS CALLED')
        CALL EIRENE_MASAGE
     .  ('       AT ENTRIES IF0COP (GEOMETRY)')
        CALL EIRENE_MASAGE
     .  ('                  IF1COP (BACKGROUND MEDIUM)')
        CALL EIRENE_MASAGE
     .  ('                  IF2COP (BOUNDARY CONDITIONS)')
        IF (NMODE.GT.0) THEN
          CALL EIRENE_MASAGE('       RETURN DATA TO EXTERNAL CODE:')
          CALL EIRENE_MASAGE('       ENTRIES IF3COP AND IF4COP ARE')
          CALL EIRENE_MASAGE('       CALLED AT THE END OF EACH STRATUM')
          CALL EIRENE_MASAGE('       AND AT THE END OF THE RUN, RESP.')
        ELSE
          CALL EIRENE_MASAGE
     .  ('       NO RETURN OF DATA TO EXTERNAL CODE:')
          CALL EIRENE_MASAGE('       ENTRIES IF3COP AND IF4COP')
          CALL EIRENE_MASAGE('       ARE NOT CALLED')
        ENDIF
        WRITE (iunout,*) '       NMODE= ',NMODE
      ELSE
        CALL EIRENE_MASAGE
     .  ('       EIRENE RUN AS STAND ALONE CODE')
        CALL EIRENE_MASAGE
     .  ('       INTERFACING ROUTINE INFCOP IS NOT CALLED')
      ENDIF
      CALL EIRENE_LEER(1)
      WRITE (iunout,*) '       EIRENE ASSUMES A TOTAL CPUTIME'
      WRITE (iunout,*) '       OF ',NTCPU,' SECONDS'
      CALL EIRENE_LEER(1)

      IF (NFILEN.EQ.1) THEN
        WRITE (iunout,*) '       EIRENE SAVES OUTPUT DATA'
        WRITE (iunout,*) '       ON FILES FT10 AND FT11 AFTER'
        WRITE (iunout,*)
     .    '       HAVING COMPUTED THE PARTICLE HISTORIES'
      ELSEIF (NFILEN.EQ.2) THEN
        WRITE (iunout,*) '       EIRENE READS OUTPUT DATA FROM'
        WRITE (iunout,*)
     .    '       AN EARLIER RUN FROM FILES FT10 AND FT11'
        WRITE (iunout,*) '       NO NEW HISTORIES ARE COMPUTED'
      ELSEIF (NFILEN.EQ.6) THEN
        WRITE (iunout,*) '       EIRENE SAVES OUTPUT DATA'
        WRITE (iunout,*) '       ON FILES FT10 AND FT11 AFTER'
        WRITE (iunout,*)
     .    '       HAVING COMPUTED THE PARTICLE HISTORIES'
        WRITE (iunout,*) '       FOR THE SUM OVER STRATA TALLIES ONLY'
      ELSEIF (NFILEN.EQ.7) THEN
        WRITE (iunout,*) '       EIRENE READS OUTPUT DATA FROM'
        WRITE (iunout,*)
     .    '       AN EARLIER RUN FROM FILES FT10 AND FT11'
        WRITE (iunout,*) '       FOR THE SUM OVER STRATA TALLIES ONLY'
        WRITE (iunout,*) '       NO NEW HISTORIES ARE COMPUTED'
      ENDIF
      IF (NFILEN.NE.0) CALL EIRENE_LEER(1)

      IF (NFILEM.EQ.1) THEN
        WRITE (iunout,*) '       EIRENE SAVES GEOMETRICAL DATA'
        WRITE (iunout,*) '       ON FILE FT12'
      ELSEIF (NFILEM.EQ.2) THEN
        WRITE (iunout,*) '       EIRENE READS GEOMETRICAL DATA FROM'
        WRITE (iunout,*) '       AN EARLIER RUN FROM FILE FT12'
      ENDIF
      IF (NFILEM.NE.0) CALL EIRENE_LEER(1)

      IF (NFILEL.GE.6.AND.NFILEL.LE.9) THEN
        NFILEL = NFILEL - 5
        NLSHRT13 = .TRUE.
      END IF
      SELECT CASE (NFILEL)
        CASE (1)
          IF (NLSHRT13) THEN
            WRITE (iunout,*) '       EIRENE SAVES ONLY PLASMA DATA,' 
            WRITE (iunout,*) '       STARTING FROM SPECIES NUMBER'
            WRITE (iunout,*) '       NFLA+1 (COUPLING VARIABLE, BLOCK'
            WRITE (iunout,*) '       14),'
          ELSE
            WRITE (iunout,*) '       EIRENE SAVES PLASMA DATA, A&M DATA'
            WRITE (iunout,*) '       AND SOURCE DISTRIBUTION DATA'
          END IF
          WRITE (iunout,*) '       ON FILE FT13 AT END OF RUN, I.E.,'
          WRITE (iunout,*) '       AFTER LAST TIMESTEP OR ITERATION'
        CASE (2)
          IF (NLSHRT13) THEN
            WRITE (iunout,*) '       EIRENE READS (AND EXPECTS) ONLY'
            WRITE (iunout,*) '       PLASMA DATA, STARTING FROM SPECIES'
            WRITE (iunout,*) '       NUMBER NFLA+1 (COUPLING VARIABLE,'
            WRITE (iunout,*) '       BLOCK 14),'
          ELSE
            WRITE (iunout,*) '       EIRENE READS PLASMA, A&M DATA'
            WRITE (iunout,*) '       AND SOURCE DISTRIBUTION DATA'
          ENDIF
          WRITE (iunout,*) '       FROM FILE FT13 '
        CASE (3)
          IF (NLSHRT13) THEN
            WRITE (iunout,*) '       EIRENE READS (AND EXPECTS) ONLY'
            WRITE (iunout,*) '       PLASMA DATA, STARTING FROM SPECIES'
            WRITE (iunout,*) '       NUMBER NFLA+1 (COUPLING VARIABLE,'
            WRITE (iunout,*) '       BLOCK 14),'
          ELSE
            WRITE (iunout,*) '       EIRENE READS PLASMA, A&M DATA'
            WRITE (iunout,*) '       AND SOURCE DISTRIBUTION DATA'
          ENDIF
          WRITE (iunout,*) '       FROM FILE FT13 AND'
          IF (NLSHRT13) THEN
            WRITE (iunout,*) '       SAVES ONLY PLASMA DATA, STARTING'
            WRITE (iunout,*) '       FROM SPECIES NUMBER NFLA+1,'
          ELSE
            WRITE (iunout,*) '       SAVES PLASMA DATA, A&M DATA'
            WRITE (iunout,*) '       AND SOURCE DISTRIBUTION DATA'
          END IF
          WRITE (iunout,*) '       ON FILE FT13 AT END OF RUN, I.E.,'
          WRITE (iunout,*) '       AFTER LAST TIMESTEP OR ITERATION'
        CASE (4)
          IF (NLSHRT13) THEN
            WRITE (iunout,*) '       EIRENE READS (AND EXPECTS) ONLY'
            WRITE (iunout,*) '       PLASMA DATA, STARTING FROM SPECIES'
            WRITE (iunout,*) '       NUMBER NFLA+1 (COUPLING VARIABLE,'
            WRITE (iunout,*) '       BLOCK 14),'
          ELSE
            WRITE (iunout,*) '       EIRENE READS PLASMA AND A&M DATA'
          END IF
          WRITE (iunout,*) '       FROM FILE FT13 AND'
          IF (NLSHRT13) THEN
            WRITE (iunout,*) '       SAVES ONLY PLASMA DATA, STARTING'
            WRITE (iunout,*) '       FROM SPECIES NUMBER NFLA+1,'
          ELSE
            WRITE (iunout,*) '       SAVES PLASMA DATA, A&M DATA'
            WRITE (iunout,*) '       AND SOURCE DISTRIBUTION DATA'
          END IF
          WRITE (iunout,*) '       ON FILE FT13 AT END OF RUN, I.E.,'
          WRITE (iunout,*) '       AFTER LAST TIMESTEP OR ITERATION'
          WRITE (iunout,*) '       SOURCE DISTRIBUTION IS NEWLY'
          WRITE (iunout,*) '       DETERMINED'
      END SELECT
      IF (NFILEL.NE.0) CALL EIRENE_LEER(1)

      SELECT CASE (NFILEK)
        CASE (1)
          WRITE (iunout,*)
     .      '       EIRENE SAVES DATA FOR RECOMMENDED INPUT'
          WRITE (iunout,*) 
     .      '       MODIFICATIONS AND RUN TIME PER STRATUM ON FILE FT14'
        CASE (2)
          WRITE (iunout,*)
     .      '       EIRENE READS DATA FOR RECOMMENDED INPUT'
          WRITE (iunout,*)
     .      '       MODIFICATIONS AND RUN TIME PER STARTUM FROM FILE'
          WRITE (iunout,*) 
     .      '       FT14, AND CARRIES OUT INPUT MODIFICATIONS THIS RUN'
        CASE (3)
          WRITE (iunout,*)
     .      '       EIRENE READS OLD DATA FOR RECOMMENDED INPUT'
          WRITE (iunout,*)
     .      '       MODIFICATIONS AND RUN TIME PER STARTUM FROM FILE'
          WRITE (iunout,*) 
     .      '       FT14, AND CARRIES OUT INPUT MODIFICATIONS THIS RUN'
          WRITE (iunout,*)
     .      '       EIRENE SAVES NEW DATA FOR RECOMMENDED INPUT'
          WRITE (iunout,*)
     .      '       MODIFICATIONS AND RUN TIME PER STRATUM ON FILE FT14'
          WRITE (iunout,*) '       FOR NEXT RUN'
        CASE (4)
          WRITE (iunout,*)
     .      '       EIRENE READS OLD RUN TIME PER STRATUM FROM FILE'
          WRITE (iunout,*) '       FT14'
        CASE (5)
          WRITE (iunout,*)
     .      '       EIRENE READS OLD RUN TIME PER STARTUM FROM FILE'
          WRITE (iunout,*) '       FT14'
          WRITE (iunout,*)
     .      '       EIRENE SAVES NEW RUN TIME PER STRATUM ON FILE FT14'
          WRITE (iunout,*) '       FOR NEXT RUN'
      END SELECT
      IF (NFILEK.NE.0) CALL EIRENE_LEER(1)

      IF (NFILEJ.EQ.1.AND.NTIME.GT.0) THEN
        WRITE (iunout,*) '       EIRENE SAVES SNAPSHOT POPULATION AT'
        WRITE (iunout,*) '       END OF LAST TIMESTEP ON FILE FT15'
      ELSEIF (NFILEJ.EQ.2) THEN
        WRITE (iunout,*) '       EIRENE READS SNAPSHOT POPULATION FROM'
        WRITE (iunout,*) '       FILE FT15'
      ELSEIF (NFILEJ.EQ.3.AND.NTIME.GT.0) THEN
        WRITE (iunout,*) '       EIRENE READS SNAPSHOT POPULATION FOR'
        WRITE (iunout,*) '       STRATUM NSTRAI+1 FOR FIRST TIMESTEP'
        WRITE (iunout,*) '       FROM FILE FT15'
        WRITE (iunout,*) '       EIRENE SAVES NEW SNAPSHOT POPULATION'
        WRITE (iunout,*) '       AT END OF LAST TIMESTEP ON FILE FT15'
      ELSEIF (NFILEJ.EQ.3.AND.NTIME.EQ.0) THEN
        WRITE (iunout,*) '       EIRENE READS SNAPSHOT POPULATION FOR'
        WRITE (iunout,*) '       FIRST TIMESTEP FROM FILE FT15'
        WRITE (iunout,*) '       NO FURTHER SNAPSHOP PRODUCED'
        WRITE (iunout,*) '       DUE TO NTIME=0'
      ENDIF
      IF (NFILEJ.NE.0) CALL EIRENE_LEER(1)

      IF (NITER.GE.1) THEN
        WRITE (iunout,*) '       EIRENE RUN IN ITERATIVE MODE.'
        WRITE (iunout,*) '       ITERATIONS: ',IITER,' TO ',NITER
        WRITE (iunout,*)
     .    '       SUBROUTINE "MODUSR" IS CALLED AFTER EACH'
        WRITE (iunout,*) '       ITERATION'
      ELSE
        WRITE (iunout,*) '       EIRENE RUN IN NON-ITERATIVE MODE'
      ENDIF
      CALL EIRENE_LEER(1)
      IF (NTIME.GE.1) THEN
        WRITE (iunout,*) '       EIRENE RUN IN TIME DEP. MODE.'
        WRITE (iunout,*) '       TIME-CYCLES: ',ITIMV,' TO ',NTIME
        WRITE (iunout,*)
     .    '       SUBROUTINE "TMSUSR" IS CALLED AFTER EACH'
        WRITE (iunout,*) '       TIME-CYCLE'
      ELSE
        WRITE (iunout,*) '       EIRENE RUN IN STATIONARY MODE'
      ENDIF
      CALL EIRENE_LEER(1)
C
C  to be done:
c  printout:  which default storage settings have we overruled (consequences thereof?)
C
C
C  READ DATA FOR STANDARD MESH, 200---299
C
      IF (IREAD == 0) READ (IUNIN,*)
      CALL EIRENE_MASAGE('*** 2. DATA FOR STANDARD MESH')
      CALL EIRENE_LEER(1)
      IREAD=0
C
  201 READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1) .EQ. '*') THEN
        WRITE (iunout,'(1X,A)') trim(ZEILE)
        CALL EIRENE_LEER(1)
        GOTO 201
      ENDIF
      READ (ZEILE,6666) (INDGRD(J),J=1,3)
C
C INPUT SUB-BLOCK 2A
C
C  RADIAL MESH
  210 READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1) .EQ. '*') THEN
        GOTO 210
      ENDIF
      READ (ZEILE,6665) NLRAD
      IF (NLRAD) THEN
C
* For gfortran: it does not accept empty field for logicals
        READ (IUNIN,'(A72)') ZEILE
        call fix_logical_input(zeile,8)
        READ (ZEILE,6665) NLSLB,NLCRC,NLELL,NLTRI,NLPLG,NLFEM,NLTET,
     .                    NLGEN
        NLPLG_IN = NLPLG
        NLFEM_IN = NLFEM

! CHECK GEOMETRY SWITCHES
       ILOGS = 0
       IF (NLSLB) ILOGS = ILOGS + 1
       IF (NLCRC) ILOGS = ILOGS + 1
       IF (NLELL) ILOGS = ILOGS + 1
       IF (NLTRI) ILOGS = ILOGS + 1
       IF (NLPLG) ILOGS = ILOGS + 1
       IF (NLFEM) ILOGS = ILOGS + 1
       IF (NLTET) ILOGS = ILOGS + 1
       IF (NLGEN) ILOGS = ILOGS + 1
       IF ((ILOGS == 0) .OR. (ILOGS > 1)) THEN
         WRITE (IUNOUT,*) ' ERROR IN GEOMETRY SPECIFICATION'
         WRITE (IUNOUT,*) ' ONE AND ONLY ONE OF THE FOLLOWING',
     .                    ' FLAGS MAY BE .TRUE. '
         WRITE (IUNOUT,*) 'NLSLB = ',NLSLB
         WRITE (IUNOUT,*) 'NLCRC = ',NLCRC
         WRITE (IUNOUT,*) 'NLELL = ',NLELL
         WRITE (IUNOUT,*) 'NLTRI = ',NLTRI
         WRITE (IUNOUT,*) 'NLPLG = ',NLPLG
         WRITE (IUNOUT,*) 'NLFEM = ',NLFEM
         WRITE (IUNOUT,*) 'NLTET = ',NLTET
         WRITE (IUNOUT,*) 'NLGEN = ',NLGEN
         CALL EIRENE_EXIT_OWN(1)
       END IF

       READ (IUNIN,6666) NR1ST,NRSEP,NRPLG,NPPLG,NRKNOT,NCOOR
       NR1ST_IN = NR1ST
       IF (NR1ST < 0) THEN
         NR1ST = N1ST
       ENDIF
        IF (INDGRD(1).LE.5) THEN
          IF (NLSLB.OR.NLCRC.OR.NLELL.OR.NLTRI) THEN
            READ (IUNIN,6664) RIA,RGA,RAA,RRA
            IF (NLELL.OR.NLTRI) THEN
              READ (IUNIN,6664) EP1IN,EP1OT,EP1CH,EXEP1
              READ (IUNIN,6664) ELLIN,ELLOT,ELLCH,EXELL
              IF (NLTRI) THEN
                READ (IUNIN,6664) TRIIN,TRIOT,TRICH,EXTRI
              ENDIF
            ENDIF
          ENDIF
          IF (NLPLG) THEN
            READ(IUNIN,6664) XPCOR,YPCOR,ZPCOR,PLREFL
            READ (IUNIN,6666) (NPOINT(1,K),NPOINT(2,K),K=1,NPPLG)
            DO 212 I=1,NR1ST
              READ (IUNIN,6664) (XPOL(I,J),YPOL(I,J),J=1,NRPLG)
  212       CONTINUE
            IF (PLREFL.GT.0.D0) NR1ST=NR1ST+1
          ENDIF
          IF (NLFEM .OR. NLTET) THEN
            READ (IUNIN,'(A72)') ZEILE
            IREAD = 1
            CLAB = ZEILE(1:4)
            CALL EIRENE_UPPERCASE(CLAB)
            IF (INDEX(ZEILE,'CASE')==0) THEN
              WRITE (IUNOUT,*) ' ERROR IN GEOMETRY SPECIFICATION'
              WRITE (IUNOUT,*)
     .          ' TRIANGLE OR TETRAHEDRON GRID SWITCHED ON'
              WRITE (IUNOUT,*) ' BUT NO CASENAME SPECIFIED'
              CALL EIRENE_EXIT_OWN(1)
            END IF

            READ (ZEILE(6:),'(A66)') CASENAME
            IREAD = 0
            CASENAME=ADJUSTL(CASENAME)
            I2=INDEX(CASENAME,' ')

            IF (NLFEM) CALL EIRENE_READ_TRIANG (CASENAME(1:I2))
            IF (NLTET) CALL EIRENE_READ_TETRA (CASENAME(1:I2))

            READ (IUNIN,'(A72)') ZEILE
            IREAD = 1
            IF ( (ZEILE(1:1) .NE. '*') .AND.
C  DOES THE NEXT LINE START WITH A 'LOGICAL' INPUT
     .           (INDEX('FTft',ZEILE(1:1)) == 0) ) THEN
C  NO: READ COORDINATE TRANSLATION CARD
              READ (ZEILE,6664) XPCOR,YPCOR,ZPCOR
              IREAD = 0
            END IF
          END IF
          IF (NLGEN) THEN
            NRGEN=NR1ST
          ENDIF
        ELSEIF (INDGRD(1).EQ.6) THEN
C  IS THERE ONE MORE LINE, OR IS NLPOL THE NEXT VARIABLE
          READ (IUNIN,'(A72)') ZEILE

CDR  SKIP READING COMMENT INPUT CARDS STARTING WITH *
          do while (zeile(1:1) == '*')
            WRITE (iunout,'(1X,A)') trim(ZEILE)
            READ (IUNIN,'(A72)') ZEILE
          end do

          IREAD = 1
          IPOS1=INDEX(ZEILE,'T')
          IPOS2=INDEX(ZEILE,'F')
          IF (IPOS1.GT.0.OR.IPOS2.GT.0) THEN
            WRITE (iunout,*) 'ONE INPUT LINE MISSING IN BLOCK 2A'
            WRITE (iunout,*) 'AUTOMATIC CORRECTION PERFORMED'
            READ (ZEILE,'(L1)') NLPOL
            IREAD = 0
            GOTO 222
          ENDIF
          IF (NLSLB.OR.NLCRC.OR.NLELL.OR.NLTRI) THEN
            READ (ZEILE,6664) RIA,RGA,RAA
            IREAD = 0
          ELSEIF (NLPLG.OR.NLFEM.OR.NLTET) THEN
            READ (ZEILE,6664) XPCOR,YPCOR,ZPCOR
            IREAD = 0
          ENDIF
        ENDIF
      ENDIF
C
C  POLOIDAL MESH
C
C INPUT SUB-BLOCK 2B

!pb   do not change previously defined IREAD
!pb   IREAD=0
      CALL EIRENE_SKIP_READ_COMMENT(IREAD,IUNIN,ZEILE)

      READ (ZEILE,6665) NLPOL
      IREAD = 0
C
 222  CONTINUE
      NLPOL_IN = NLPOL
      READ (IUNIN,'(A72)') ZEILE
* For gfortran: it does not accept empty field for logicals
      call fix_logical_input(zeile,3)
      READ (ZEILE,6665) NLPLY,NLPLA,NLPLP
      READ (IUNIN,6666) NP2ND,NPSEP,NPPLA,NPPER
      NP2ND_IN = NP2ND
      IF (INDGRD(2).LE.5) THEN
        READ (IUNIN,6664) YIA,YGA,YAA,YYA
      ELSEIF (INDGRD(2).EQ.6) THEN
      ENDIF
C
C  TOROIDAL MESH
C
C INPUT SUB-BLOCK 2C
C
      IREAD=0
      CALL EIRENE_SKIP_READ_COMMENT(IREAD,IUNIN,ZEILE)

      READ (ZEILE,6665) NLTOR
      IREAD=0
C
      READ (IUNIN,'(A80)') ZEILE
* For gfortran: it does not accept empty field for logicals
      call fix_logical_input(zeile,3)
      READ (ZEILE,6665) NLTRZ,NLTRA,NLTRT
      READ (IUNIN,6666) NT3RD,NTSEP,NTTRA,NTPER
      NT3RD_IN = NT3RD
      IF (INDGRD(3).LE.5) THEN
        READ (IUNIN,6664) ZIA,ZGA,ZAA,ZZA,ROA
      ELSEIF (INDGRD(3).EQ.6) THEN
      ENDIF
C
C  MESH MULTIPLICATION
C
C INPUT SUB-BLOCK 2D
C
      IREAD=0
      CALL EIRENE_SKIP_READ_COMMENT(IREAD,IUNIN,ZEILE)
      READ (ZEILE,6665) NLMLT
      IREAD=0
C
      IF (NLMLT) THEN
        READ (IUNIN,6666) NBMLT
        READ (IUNIN,6664) (VOLCOR(NM),NM=1,NBMLT)
        READ (IUNIN,'(A72)') ZEILE
      ELSE
        NBMLT=1
        VOLCOR(1)=1.D0
C  FIND START OF NEXT INPUT BLOCK: 2E. SEARCH FOR *, T OR F
  241   READ (IUNIN,'(A72)') ZEILE
        CALL EIRENE_UPPERCASE(ZEILE)
        IPOS0=INDEX(ZEILE,'*')
        IPOS1=INDEX(ZEILE,'T')
        IPOS2=INDEX(ZEILE,'F')
        IF (IPOS0.EQ.0.AND.IPOS1.EQ.0.AND.IPOS2.EQ.0) GOTO 241
      ENDIF
      IREAD=1
C
C  ADDITIONAL CELLS OUTSIDE STANDARD MESH
C
  250 IF (IREAD.EQ.0) READ (IUNIN,'(A72)') ZEILE
C
C INPUT SUB-BLOCK 2E
C
      IF (ZEILE(1:1) .EQ. '*') THEN
        IREAD=0
        GOTO 250
      ENDIF

      READ (ZEILE,6665) NLADD
      IREAD=0
C
      IF (NLADD) THEN
        READ (IUNIN,6666) NRADD
        READ (IUNIN,6664) (VOLADD(NM),NM=1,NRADD)
      ELSE
C  FIND START OF NEXT INPUT BLOCK: 3A
  252   READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:3) .NE. '***') GOTO 252
        IREAD=1
      ENDIF
C
C
C  READING FOR INPUT BLOCK 2 DONE
C
      NBMLT=MAX0(1,NBMLT)
      NR1ST=MAX0(1,NR1ST)
      NP2ND=MAX0(1,NP2ND)
      IF (.NOT.NLPOL) NP2ND=1
      NT3RD=MAX0(1,NT3RD)
      IF (.NOT.NLTOR) NT3RD=1
C
      IF (NLPLG.AND.INDGRD(1).LE.4) THEN
        IF (NLPOL.AND.NRPLG.NE.NP2ND) GOTO 994
      ENDIF
C
C  * 3A: READ DATA FOR NON-DEFAULT SURFACE MODELS ON STANDARD SURFACES
C  300--349
C
      IF (IREAD.EQ.0) READ (IUNIN,*)
      CALL EIRENE_MASAGE
     .  ('*** 3A. DATA FOR NON-DEFAULT STANDARD SURFACES')
  310 READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1) .EQ. '*') GOTO 310
      IREAD=1
      READ(ZEILE,6666) NSTSI
      IREAD=0
      WRITE (iunout,*) '        NSTSI= ',NSTSI
      CALL EIRENE_LEER(1)
      DO 311 ISTS=1,NSTSI
        NLJ=NLIM+ISTS
        IF (IREAD.EQ.0) THEN
          READ (IUNIN,'(A72)') TXTSFL(NLJ)
        ELSE
          READ (ZEILE,'(A72)') TXTSFL(NLJ)
          IREAD=0
        ENDIF
        READ (IUNIN,6666) JDUMMY,IDIMP,INUMP(ISTS,IDIMP),IRPTA1,
     .                    IRPTE1,IRPTA2,IRPTE2,IRPTA3,IRPTE3

        IF (.NOT.(NLFEM.OR.NLTET.OR.NLGEN)) THEN

          IF ((IDIMP == 1) .AND. (INUMP(ISTS,IDIMP) > N1ST)) THEN
            WRITE (iunout,*) ' ERROR IN SPECIFICATION OF NON-DEFAULT'
            WRITE (iunout,*) ' SURFACE ',ISTS
            WRITE (iunout,*) ' NUMBER OF RADIAL SURFACE > N1ST'
            WRITE (iunout,*) ' CHECK INPUT FILE'
            CALL EIRENE_EXIT_OWN(1)
          ELSEIF ((IDIMP == 2) .AND. (INUMP(ISTS,IDIMP) > N2ND)) THEN
            WRITE (iunout,*) ' ERROR IN SPECIFICATION OF NON-DEFAULT'
            WRITE (iunout,*) ' SURFACE ',ISTS
            WRITE (iunout,*) ' NUMBER OF POLOIDAL SURFACE > N2ND'
            WRITE (iunout,*) ' CHECK INPUT FILE'
            CALL EIRENE_EXIT_OWN(1)
          ELSEIF ((IDIMP == 3) .AND.
     .            ((NLTOR.AND.(INUMP(ISTS,IDIMP) > N3RD)) .OR.
     .             (NLTRA.AND.(INUMP(ISTS,IDIMP) > NTTRA)))) THEN
            WRITE (iunout,*) ' ERROR IN SPECIFICATION OF NON-DEFAULT'
            WRITE (iunout,*) ' SURFACE ',ISTS
            WRITE (iunout,*) ' NUMBER OF TOROIDAL SURFACE > N3RD'
            WRITE (iunout,*) ' CHECK INPUT FILE'
            CALL EIRENE_EXIT_OWN(1)
          END IF

        END IF
C
C  OLD INPUT VERSION BEGIN
        IF (IDIMP.EQ.1.AND.IRPTA1.NE.IRPTE1) THEN
          WRITE (iunout,*) 'WARNING FROM INPUT BLOCK 3A, ISTS= ',ISTS
          WRITE (iunout,*) 'NEW INPUT FOR IRPTA,IRPTE....'
          WRITE (iunout,*) 'AUTOMATIC CORRECTION CARRIED OUT'
          IRPTA2=IRPTA1
          IRPTE2=IRPTE1
          IRPTA1=INUMP(ISTS,1)
          IRPTE1=INUMP(ISTS,1)
        ELSEIF (IDIMP.EQ.1.AND.IRPTA1.EQ.IRPTE1.AND.
     .          IRPTA2+IRPTE2+IRPTA3+IRPTE3.EQ.0) THEN
          IRPTA2=IRPTA1
          IRPTE2=IRPTE1
          IRPTA1=INUMP(ISTS,1)
          IRPTE1=INUMP(ISTS,1)
        ELSEIF (IDIMP.EQ.1.AND.IRPTA1.EQ.0.AND.IRPTE1.EQ.0) THEN
          IRPTA1=INUMP(ISTS,1)
          IRPTE1=INUMP(ISTS,1)
        ELSEIF (IDIMP.EQ.2.AND.IRPTA2.EQ.0.AND.IRPTE2.EQ.0) THEN
          IRPTA2=INUMP(ISTS,2)
          IRPTE2=INUMP(ISTS,2)
        ELSEIF (IDIMP.EQ.3.AND.IRPTA3.EQ.0.AND.IRPTE3.EQ.0) THEN
          IRPTA3=INUMP(ISTS,3)
          IRPTE3=INUMP(ISTS,3)
        ENDIF
C  OLD INPUT VERSION DONE
C
C  OVERWRITE DEFAULTS FOR IRPTA, IRPTE ARRAYS
        IF (IDIMP.NE.1) THEN
          IF (IRPTA1.GT.1) IRPTA(ISTS,1)=IRPTA1
          IF (IRPTE1.LE.1.OR.IRPTE1.GT.NR1ST) THEN
            IRPTE(ISTS,1)=MAX(2,NR1ST)
          ELSE
            IRPTE(ISTS,1)=IRPTE1
          ENDIF
          IF ((IDIMP.EQ.2.AND.IRPTE2.GT.NP2ND).OR.
     .        (IDIMP.EQ.3.AND.IRPTE3.GT.NT3RD)) THEN
            WRITE (iunout,*) 'WARNING:'
            WRITE (iunout,*) 'NON-DEFAULT STANDARD SURFACE NO. ',ISTS
            WRITE (iunout,*) 'OUT OF RANGE. INVISIBLE!'
          ENDIF
        ELSE ! IDIMP .EQ. 1
          IRPTA(ISTS,1) = IRPTA1
          IRPTE(ISTS,1) = IRPTE1
        ENDIF
        IF (IDIMP.NE.2) THEN
          IF (IRPTA2.GT.1) IRPTA(ISTS,2)=IRPTA2
          IF (IRPTE2.LE.1.OR.IRPTE2.GT.MAX(NP2ND,NRPLG)) THEN
            IRPTE(ISTS,2)=MAX(2,NP2ND,NRPLG)
          ELSE
            IRPTE(ISTS,2)=IRPTE2
          ENDIF
          IF ((IDIMP.EQ.1.AND.IRPTE1.GT.NR1ST).OR.
     .        (IDIMP.EQ.3.AND.IRPTE3.GT.NT3RD)) THEN
            WRITE (iunout,*) 'WARNING:'
            WRITE (iunout,*) 'NON-DEFAULT STANDARD SURFACE NO. ',ISTS
            WRITE (iunout,*) 'OUT OF RANGE. INVISIBLE!'
          ENDIF
        ELSE ! IDIMP .EQ. 2
          IRPTA(ISTS,2) = IRPTA2
          IRPTE(ISTS,2) = IRPTE2
        ENDIF
        IF (IDIMP.NE.3) THEN
          IF (IRPTA3.GT.1) IRPTA(ISTS,3)=IRPTA3
          IF (IRPTE3.LE.1.OR.IRPTE3.GT.NT3RD) THEN
            IRPTE(ISTS,3)=MAX(2,NT3RD)
          ELSE
            IRPTE(ISTS,3)=IRPTE3
          ENDIF

          IF ((IDIMP.EQ.1.AND.IRPTE1.GT.NR1ST).OR.
     .        (IDIMP.EQ.2.AND.IRPTE2.GT.MAX(NP2ND,NRPLG))) THEN
            WRITE (iunout,*) 'WARNING:'
            WRITE (iunout,*) 'NON-DEFAULT STANDARD SURFACE NO. ',ISTS
            WRITE (iunout,*) 'OUT OF RANGE. INVISIBLE!'
          ENDIF
        ELSE ! IDIMP .EQ. 3
          IRPTA(ISTS,3) = IRPTA3
          IRPTE(ISTS,3) = IRPTE3
        ENDIF

        READ (IUNIN,6666) ILIIN(NLJ),ILSIDE(NLJ),ILSWCH(NLJ),
     .                    ILEQUI(NLJ),ILTOR(NLJ),ILCOL(NLJ),
     .                    ILFIT(NLJ),ILCELL(NLJ),ILBOX(NLJ),
     .                    ILPLG(NLJ)
        ILCOL_IN(NLJ) = ILCOL(NLJ)
        IF (ABS(ILCOL(NLJ)).EQ.7) THEN
          WRITE (iunout,*) 'COLOUR FLAG ILCOL CHANGED FOR SURFACE NO. ',
     .                      NLJ
          WRITE (iunout,*) 'COLOUR NO. 7 IS RESERVED FOR "NON-ANALOG'
          WRITE (iunout,*)
     .      'SURFACES" (SPLITTING, R.R., WEIGHT WINDOWS,..)'
          ILCOL(NLJ)=ILCOL(NLJ)-2
        ENDIF
        READ (IUNIN,'(A72)') ZEILE
        IREAD=1
        IF (ZEILE(1:1).EQ.'*') THEN
C  NO LOCAL SURFACE INTERACTION MODEL FOUND, USE: DEFAULT
          GOTO 314
!pb     ELSEIF (ILIIN(NLJ).LE.0) THEN
!pb       GOTO 312
C  ASSIGN ONE OF THE SURFACE MODELS (BLOCK 6) TO THIS SURFACE
        ELSEIF (ZEILE(1:8).EQ.'SURFMOD_') THEN
          ALLOCATE(SURFCUR)
          SURFCUR%MODNAME = TRIM(ADJUSTL(ZEILE(9:)))
          SURFCUR%NOSURF = NLJ
          SURFCUR%NEXT => SURFLIST
          SURFLIST => SURFCUR
          IREAD=0
          SMOD_NAME(NLJ)=SURFCUR%MODNAME
          GOTO 314
C  READ LOCAL (FOR THIS SURFACE) SURFACE INTERACTION MODEL, NEXT 3 OR 4 INPUT CARDS
        ELSE
          READ (ZEILE,6666) ILREF(NLJ),ILSPT(NLJ),ISRS(1,NLJ),
     .                      ISRC(1,NLJ)
          IREAD=0
          READ (IUNIN,6664) ZNML(NLJ),EWALL(NLJ),EWBIN(NLJ),
     .                      TRANSP(1,1,NLJ),TRANSP(1,2,NLJ),FSHEAT(NLJ)
          READ (IUNIN,6664) RECYCF(1,NLJ),RECYCT(1,NLJ),RECPRM(1,NLJ),
     .                      EXPPL(1,NLJ),EXPEL(1,NLJ),EXPIL(1,NLJ)
          READ (IUNIN,'(A72)') ZEILE
          IREAD=1
C  READ ONE MORE LINE FOR NON-DEFAULT SPUTTER MODEL
          IF (ZEILE(1:1).NE.'*') THEN
            READ (ZEILE,6664) RECYCS(1,NLJ),RECYCC(1,NLJ),SPTPRM(1,NLJ),
     .                        ESPUTS(1,NLJ),ESPUTC(1,NLJ)
            IREAD=0
          ELSEIF (ILSPT(NLJ).NE.0) THEN
            WRITE (iunout,*) 'WARNING: SPUTTERING AT NON-DEF. SURFACE ',
     .                        ISTS
            WRITE (iunout,*)
     .        'BUT NO PARAMETERS RECYCS, RECYCC ARE READ'
            WRITE (iunout,*)
     .        'DEFAULT MODEL: "NO SPUTTERING" IS USED.'
            WRITE (iunout,*) 'DO YOU REALLY WANT THIS?'
            ILSPT(NLJ)=0
          ENDIF
        ENDIF
C
  314   CONTINUE
C
  311 CONTINUE
C
C  * 3B: READ DATA FOR ADDITIONAL SURFACES 350--399
C
      IF (IREAD.EQ.0) READ (IUNIN,*)
      CALL EIRENE_MASAGE('*** 3B. DATA FOR ADDITIONAL SURFACES')
  350 READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1).EQ.'*') THEN
        WRITE (iunout,'(1X,A)') trim(ZEILE)
        GOTO 350
      ENDIF
      IREAD=0
      READ (ZEILE,6666) NLIMI
      WRITE (iunout,*) '        NLIMI= ',NLIMI
      CALL EIRENE_LEER(1)
C
      IF (NLIMI.GT.0) THEN
        DO 353 I=1,NLIMI
          IHELP(I)=IGJUM0(I)
  353   CONTINUE
  351   READ (IUNIN,'(A72)') ZEILE
        IREAD=1
        IF (ZEILE(1:3).EQ.'CH0') THEN
          IREAD=0
          
!  STORE CH0 LINES FOR OUTPUTING TO JSON FILE 
          call eirene_push_string_stack(ch0_stack,zeile(1:72))

          CALL EIRENE_DEKEY (ZEILE(4:72),IHELP,1,1,1,NLIMPS)
          GOTO 351
        ENDIF
        DO 352 I=1,NLIMI
          IGJUM0(I)=IHELP(I)
  352   CONTINUE
      ENDIF
C
      DO 360 I=1,NLIMI
        IF (IREAD.NE.0) THEN
          READ (ZEILE,'(A72)') TXTSFL(I)
        ELSE
          READ (IUNIN,'(A72)') TXTSFL(I)
        ENDIF
        WRITE (iunout,'(1X,A)') trim(TXTSFL(I))
        IREAD=0
!pb read data for switched off surfaces in order to be able to write
!   data onto json files
!        IF (IGJUM0(I).NE.0) THEN
!361       READ (IUNIN,'(A72)') ZEILE
!          IREAD=1
!          IF (ZEILE(1:1).EQ.'*') GOTO 369
!          IF (ZEILE(1:9).EQ.'TRANSFORM') GOTO 368
!          GOTO 361
!        ENDIF
C
C   GENERAL SURFACE DATA
  362   READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:3).EQ.'CH1') THEN
          
!  STORE CH1 LINES FOR OUTPUTING TO JSON FILE 
           call eirene_push_string_stack(ch1_stack(i),zeile(1:72))

          IF (NLIMPB >= NLIMPS) THEN
            CALL EIRENE_DEKEY (ZEILE(4:72),IGJUM1,0,NLIMPS,I,NLIMPS)
          ELSE
            CALL EIRENE_DEKEYB
     .  (ZEILE(4:72),IGJUM1,0,NLIMPS,I,NLIMPB,NBITS)
          END IF
          GOTO 362

        ELSEIF (ZEILE(1:3).EQ.'CH2') THEN
        
!  STORE CH2 LINES FOR OUTPUTING TO JSON FILE 
          call eirene_push_string_stack(ch2_stack(i),zeile(1:72))
          IF (NLIMPB >= NLIMPS) THEN
            CALL EIRENE_DEKEY (ZEILE(4:72),IGJUM2,0,NLIMPS,I,NLIMPS)
          ELSE
            CALL EIRENE_DEKEYB
     .  (ZEILE(4:72),IGJUM2,0,NLIMPS,I,NLIMPB,NBITS)
          END IF
          GOTO 362
        ELSE
          READ (ZEILE,6664) RLB(I),SAREA_SAVE(I),RLWMN(I),RLWMX(I)
          IF (SAREA_SAVE(I).LE.0.D0) SAREA_SAVE(I)=666.
          READ (IUNIN,6666) ILIIN(I),ILSIDE(I),ILSWCH(I),
     .                      ILEQUI(I),ILTOR(I),ILCOL(I),
     .                      ILFIT(I),ILCELL(I),ILBOX(I),
     .                      ILPLG(I)
          IF (ABS(ILCOL(I)).EQ.7) THEN
            WRITE (iunout,*)
     .        'COLOUR FLAG ILCOL CHANGED FOR SURFACE NO. ',I
            WRITE (iunout,*) 'COLOUR NO. 7 IS RESERVED FOR "NON-ANALOG'
            WRITE (iunout,*)
     .        'SURFACES" (SPLITTING, R.R., WEIGHT WINDOWS,..)'
            ILCOL(I)=ILCOL(I)-2
          ENDIF
C  READ SURFACE COEFFICIENTS
          IF (RLB(I).LT.2.) THEN
            READ (IUNIN,6664) A0LM(I),A1LM(I),A2LM(I),A3LM(I),A4LM(I),
     .                        A5LM(I)
            READ (IUNIN,6664) A6LM(I),A7LM(I),A8LM(I),A9LM(I)
          ELSEIF (RLB(I).LT.3.) THEN
            READ (IUNIN,6664) P1(1,I),P1(2,I),P1(3,I),P2(1,I),P2(2,I),
     .                        P2(3,I)
          ELSEIF (RLB(I).LT.5.) THEN
            READ (IUNIN,6664) P1(1,I),P1(2,I),P1(3,I),P2(1,I),P2(2,I),
     .                        P2(3,I)
            READ (IUNIN,6664) P3(1,I),P3(2,I),P3(3,I),P4(1,I),P4(2,I),
     .                        P4(3,I)
          ELSEIF (RLB(I).LT.7.) THEN
            READ (IUNIN,6664) P1(1,I),P1(2,I),P1(3,I),P2(1,I),P2(2,I),
     .                        P2(3,I)
            READ (IUNIN,6664) P3(1,I),P3(2,I),P3(3,I),P4(1,I),P4(2,I),
     .                        P4(3,I)
            READ (IUNIN,6664) P5(1,I),P5(2,I),P5(3,I),P6(1,I),P6(2,I),
     .                        P6(3,I)
          ENDIF
        ENDIF
C  READ BOUNDARY DATA
        IF (RLB(I).GT.0..AND.RLB(I).LT.2.) THEN
          READ (IUNIN,6664) XLIMS1(1,I),YLIMS1(1,I),ZLIMS1(1,I),
     .                      XLIMS2(1,I),YLIMS2(1,I),ZLIMS2(1,I)
        ELSEIF (RLB(I).LE.0.D0) THEN
          IH=INT(-RLB(I))
          ILIN(I)=EIRENE_IDEZ(IH,1,2)
          ISCN(I)=EIRENE_IDEZ(IH,2,2)
          DO 363 J=1,ILIN(I)
            READ (IUNIN,6664) ALIMS(J,I),XLIMS(J,I),YLIMS(J,I),
     .                        ZLIMS(J,I)
  363     CONTINUE
          DO 364 J=1,ISCN(I)
            READ (IUNIN,6664) ALIMS0(J,I),XLIMS1(J,I),YLIMS1(J,I),
     .                        ZLIMS1(J,I),XLIMS2(J,I),YLIMS2(J,I)
            READ (IUNIN,6664) ZLIMS2(J,I),XLIMS3(J,I),YLIMS3(J,I),
     .                        ZLIMS3(J,I)
  364     CONTINUE
        ENDIF


!  STORE VERTICES OF ADDITIONAL SURFACES FOR OUTPUTING THEM ONTO JSON FILE
!  ORIGINAL VALUES MIGHT BE CHANGED LATER ON
        CALL EIRENE_COPY_ADDSRF(I)
C
C READ LOCAL SURFACE INTERACTION MODEL FOR SURFACE NO. I
C
  367   READ (IUNIN,'(A72)') ZEILE
        IREAD=1
        IF (ZEILE(1:1).EQ.'*') THEN
C  SKIP READING LOCAL SURFACE INTERACTION MODEL, USE: DEFAULT
          GOTO 369
        ELSEIF (ZEILE(1:9).EQ.'TRANSFORM') THEN
C  TRANSFORM BLOCK
          GOTO 368
        ELSEIF (ILIIN(I).LE.0) THEN
C  TRANSPARENT: SKIP READING SURFACE INTERACTION MODEL
          GOTO 367
        ELSEIF (ZEILE(1:8).EQ.'SURFMOD_') THEN
          ALLOCATE(SURFCUR)
          SURFCUR%MODNAME = TRIM(ADJUSTL(ZEILE(9:)))
          SURFCUR%NOSURF = I
          SURFCUR%NEXT => SURFLIST
          SURFLIST => SURFCUR
          IREAD=0
          SMOD_NAME(I)=SURFCUR%MODNAME
          GOTO 368
        ELSE
C  READ LOCAL REFLECTION MODEL
          READ (ZEILE,6666) ILREF(I),ILSPT(I),ISRS(1,I),ISRC(1,I)
          IREAD=0
          READ (IUNIN,6664) ZNML(I),EWALL(I),EWBIN(I),
     .                      TRANSP(1,1,I),TRANSP(1,2,I),FSHEAT(I)
          READ (IUNIN,6664) RECYCF(1,I),RECYCT(1,I),RECPRM(1,I),
     .                      EXPPL(1,I),EXPEL(1,I),EXPIL(1,I)
          READ (IUNIN,'(A72)') ZEILE
          IREAD=1
C  READ ONE MORE LINE FOR NON-DEFAULT SPUTTER MODEL
          IF (ZEILE(1:1).NE.'*'.AND.ZEILE(1:9).NE.'TRANSFORM') THEN
            READ (ZEILE,6664) RECYCS(1,I),RECYCC(1,I),SPTPRM(1,I),
     .                        ESPUTS(1,I),ESPUTC(1,I)
            IREAD=0
          ELSEIF (ILSPT(I).NE.0) THEN
            WRITE (iunout,*) 'WARNING: SPUTTERING FOR ADD. SURFACE ',I
            WRITE (iunout,*)
     .        'BUT NO PARAMETERS RECYCS, RECYCC ARE READ'
            WRITE (iunout,*) 'DEFAULT MODEL: "NO SPUTTERING" IS USED.'
            WRITE (iunout,*) 'DO YOU REALLY WANT THIS?'
            ILSPT(I)=0
          ENDIF
        ENDIF
C
  368   IF (IREAD.EQ.0) READ (IUNIN,'(A72)') ZEILE
        IREAD=1
        IF (ZEILE(1:9).NE.'TRANSFORM') GOTO 369
C
C  TRANSFORM FROM CONVENIENT COORDINATE SYSTEM TO EIRENE COORDINATE
C  SYSTEM. THIS IS POSSIBLE FOR ALL SURFACES READ BY NOW, I.E. FROM
C  ITINI=1 TO ITEND=I
        IREAD=0
        READ (IUNIN,6666) ITINI(I),ITEND(I)
        IF (ITINI(I).LT.1) ITINI(I)=I
        IF (ITEND(I).GT.I) ITEND(I)=I
        READ (IUNIN,6664) XLCOR(I),YLCOR(I),ZLCOR(I)
        READ (IUNIN,6664) XLREF(I),YLREF(I),ZLREF(I)
        READ (IUNIN,6664) XLROT(I),YLROT(I),ZLROT(I),ALROT(I)
C
c  STORE TRANSFORMATIONS FOR OUTPUT TO JSON FILE
        ALLOCATE(TRF)
        TRF%ITINI = ITINI(I)
        TRF%ITEND = ITEND(I)
        TRF%XLCOR = XLCOR(I)
        TRF%YLCOR = YLCOR(I)
        TRF%ZLCOR = ZLCOR(I)
        TRF%XLREF = XLREF(I)
        TRF%YLREF = YLREF(I)
        TRF%ZLREF = ZLREF(I)
        TRF%XLROT = XLROT(I)
        TRF%YLROT = YLROT(I)
        TRF%ZLROT = ZLROT(I)
        TRF%ALROT = ALROT(I)
        NULLIFY(TRF%NEXT)
        IF (.NOT.ASSOCIATED(TRNSFRM(I)%HEAD)) THEN
          TRNSFRM(I)%HEAD => TRF
          TRNSFRM(I)%LAST => TRF
        ELSE
          TRNSFRM(I)%LAST%NEXT => TRF
          TRNSFRM(I)%LAST => TRF
        END IF
c  SHIFT
        XSH=XLCOR(I)
        IF (XSH.NE.0.D0) CALL EIRENE_XSHADD(XSH,ITINI(I),ITEND(I))
        YSH=YLCOR(I)
        IF (YSH.NE.0.D0) CALL EIRENE_YSHADD(YSH,ITINI(I),ITEND(I))
        ZSH=ZLCOR(I)
        IF (ZSH.NE.0.D0) CALL EIRENE_ZSHADD(ZSH,ITINI(I),ITEND(I))
C  REFLECTION
        REFNRM=XLREF(I)*XLREF(I)+YLREF(I)*YLREF(I)+ZLREF(I)*ZLREF(I)
        IF (REFNRM.GT.EPS10) THEN
          CALL EIRENE_SETREF(AFF,AFFI,1,XLREF(I),YLREF(I),ZLREF(I))
          CALL EIRENE_ROTADD(AFF,AFFI,ITINI(I),ITEND(I))
        ENDIF
C  ROTATION
        ROTNRM=XLROT(I)*XLROT(I)+YLROT(I)*YLROT(I)+ZLROT(I)*ZLROT(I)
        ALR=ALROT(I)
        IF (ALR.NE.0..AND.ROTNRM.GT.EPS10) THEN
          CALL EIRENE_SETROT(AFF,AFFI,1,XLROT(I),YLROT(I),ZLROT(I),
     .                                  ALROT(I))
          CALL EIRENE_ROTADD(AFF,AFFI,ITINI(I),ITEND(I))
        ENDIF
        GOTO 368
C
  369   CONTINUE
C
  360 CONTINUE

      CALL EIRENE_LEER(1)
C
C  READ DATA FOR SPECIES SPECIFICATION AND ATOMIC PHYSICS MODULE
C  400--499
C
C  AT THIS POINT THE INPUT LINE *** 4.  .... IS EXPECTED
      IF (IREAD.EQ.0) READ (IUNIN,'(A72)') ZEILE
      IREAD=0
      IF ((ZEILE(1:3) .NE. '***') .OR. (INDEX(ZEILE,'4.') == 0)) THEN
        WRITE (IUNOUT,*) 'INPUT ERROR, BLOCK *** 4. NOT FOUND'
        CALL EIRENE_EXIT_OWN(1)
      ENDIF

      CALL EIRENE_MASAGE('*** 4. DATA FOR SPECIES SPECIFICATION AND')
      CALL EIRENE_MASAGE('       ATOMIC PHYSICS MODULE')
      CALL EIRENE_LEER(1)

! CHECK FOR INCLUDE LINE
      READ (IUNIN,'(A420)') ZEILE
      IREAD=1  !  next input line is already read from iunin, now on 'ZEILE'
      ULINE = ZEILE
      CALL EIRENE_UPPERCASE(ULINE)
      I1 = INDEX(ULINE,'INCLUDE')
      LINCL45 = .FALSE.
      IUNIN_SAVE = IUNIN

      IF (I1 > 0) THEN

C  "Include" found. Skip all the rest of block 4 and 5 of input file (fort.iunin),
C   and read this information only from the "include" file instead
C   stream: 2+ifoff
C   Zeile  = INCLUDE 'FILE45'
C
        LINCL45 = .TRUE.
C

        CALL EIRENE_READ_TOKEN(ZEILE(I1+7:),' ',FILE45,ITOK,IER,.FALSE.)
        IREAD = 0

        IUNIN_SAVE = IUNIN
        IUNIN = 2+ifoff !  FILE45 is the include file. read block 4 and 5 from
!                          FILE45 (stream fort.2) rather than from stream fort.iunin
!                          close fort.2 at end of block 5.
        WRITE (IUNOUT,*) 'EXTERNAL A&M INPUT BLOCK 4 AND 5 FOUND'
        WRITE (IUNOUT,*) 'FILE45 = ',TRIM(FILE45)
        CALL EIRENE_LEER(1)
        OPEN (IUNIN,FILE=FILE45,FORM='FORMATTED',ACCESS='SEQUENTIAL')
c  read comment lines on external A&M data file FILE45, stream fort.2
  401   READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 401
        IREAD=1
        GOTO 402
      END IF

C   no "Include" found. Read block 4 and 5 from input
C   stream: IUNIN
C
      IF (IREAD == 0) READ (IUNIN,*)

      WRITE (iunout,*)
     .  '       ATOMIC REACTION CARDS, NREACI DATA FIELDS'

      READ (IUNIN,'(A72)') ZEILE  ! THIS IS THE FIRST NON-COMMENT LINE
  402 CALL EIRENE_UPPERCASE(ZEILE)

chk.................................................................
C  special only in case of HYDKIN INTERFACE: find string "DEFAULT"
      IEND=INDEX(ZEILE,'DEFAULT')
      LHYDDEF =.FALSE.
      IF (IEND > 0) THEN
        LHYDDEF =.TRUE.
CDR  lhyddef is currently only supported in proprietary versions of eirene.
        WRITE (IUNOUT,*) 'LHYDDEF: OPTION IS NOT READY, EXIT CALLED'
        CALL EIRENE_EXIT_OWN(1)
CDR
        CALL
     .  EIRENE_READ_TOKEN(ZEILE(IEND+7:),' ',HYDKIN_DEFAULT,ITOK,IER,
     .                  .FALSE.)

        IEND = IEND + 7 + ITOK
        CALL EIRENE_READ_TOKEN(ZEILE(IEND+1:),' ',CADAPT,ITOK,IER,
     .                  .FALSE.)

        READ (IUNIN,'(A72)') ZEILE

      END IF
chk..................................................................

C  Normal start of reading database A&M processes

      READ (ZEILE,*) NREACI
      WRITE (iunout,*) '       NREACI= ',NREACI
      CALL EIRENE_LEER(1)

      IF (NPHOTI > 0) CALL EIRENE_PH_INIT(0)
!pbcrm
      IREAD = 0
C
  411 IF (IREAD == 0) READ (IUNIN,'(A420)') ZEILE
      IF (ZEILE(1:1).NE.'*') THEN

!  STORE REACTION LINES FOR OUTPUTING TO JSON FILE 
        call eirene_push_string_stack(crs_stack,trim(zeile)) 
C
C  READ ONE REACTION FROM FILE "FILNAM" AT A TIME. Input card is on "ZEILE"
C

        CALL EIRENE_READ_REACLINES (ZEILE, IUNIN, IUNOUT, IREAD,
     .       IR, FILNAM, H123, REAC, CRC, MP, MT, DPP,
     .       RC1MIN, RC1MAX, JFEX1MN, JFEX1MX, FP1,
     .       RC2MIN, RC2MAX, JFEX2MN, JFEX2MX, FP2,
     .       IZ, ELNAME, IROW_ESC, ICOL_ESC, POP_ESC, 
     .       IFTFL, NCOEF, COEF)

C  SAVE SOME OF THE INPUT FLAGS FOR LATER
C  PROCESSING (MASS SCALING, POTENTIAL ENERGY INCREMENT) IN XSTCX,XSTEI,...

        MASSP(IR)=MP
        MASST(IR)=MT
        DELPOT(IR)=DPP
C
        CALL EIRENE_SLREAC (IR,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX, ! additional input card: asymptotics P1
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX, ! additional input card: asymptotics P2
     .               ELNAME,IZ,                 ! additional input card read for TAB2D/ADAS format
     .               IROW_ESC,ICOL_ESC,POP_ESC, ! (optional) additional input card read CR  internal models
     .               IFTFL, NCOEF, COEF)        ! (optional) for "CONST models"
        GOTO 411
      ENDIF
      IF (TRCAMD) CALL EIRENE_LEER(1)
C
      WRITE (iunout,*)
     . '*** 4A. NEUTRAL ATOMS SPECIES CARDS, NATMI SPECIES'
      READ (IUNIN,*) NATMI
      WRITE (iunout,*) '       NATMI= ',NATMI
      CALL EIRENE_LEER(1)
      NATMI_IN = NATMI
C
      NSPH=NPHOTI
      DO 421 JATM=1,NATMI
        ISPZ=NSPH+JATM
        READ (IUNIN,66666) I,TEXTS(ISPZ),NMASSA(JATM),NCHARA(JATM),
     .                       NDUMM1,NDUMM2,  !NPRT=1, NCHRGA=0, DEFAULT
     .                       ISRF(ISPZ,1),ISRT(ISPZ,1),NUMSEC,
     .                       NRCA(JATM),NFOLA(JATM),NGENA(JATM),
     .                       NHSTS(ISPZ)
C  DEFAULTS FOR ATOMIC SPECIES:
        NPRT(ISPZ)=1
        DO 422 K=1,NRCA(JATM)
          IF (NUMSEC .LT. 3) THEN
            READ (IUNIN,6666) IREACA(JATM,K),IBULKA(JATM,K),
     .                        ISCD1A(JATM,K),ISCD2A(JATM,K),
     .                        ISCDEA(JATM,K),IESTMA(JATM,K),
     .                        IBGKA(JATM,K)
          ELSE IF (NUMSEC .EQ. 3) THEN
            READ (IUNIN,6666) IREACA(JATM,K),IBULKA(JATM,K),
     .                        ISCD1A(JATM,K),ISCD2A(JATM,K),
     .                        ISCD3A(JATM,K),
     .                        ISCDEA(JATM,K),IESTMA(JATM,K),
     .                        IBGKA(JATM,K)
            WRITE (iunout,*) ' WARNING !!! '
            WRITE (iunout,*)
     .        ' THREE SECONDARY GROUPS USED FOR SPECIES ',
     .          TEXTS(ISPZ)
            WRITE (iunout,*) ' ISCD1A = ',ISCD1A(JATM,K)
            WRITE (iunout,*) ' ISCD2A = ',ISCD2A(JATM,K)
            WRITE (iunout,*) ' ISCD3A = ',ISCD3A(JATM,K)
            WRITE (iunout,*) ' ISCDEA = ',ISCDEA(JATM,K)
            WRITE (iunout,*) ' IESTMA = ',IESTMA(JATM,K)
            WRITE (iunout,*) ' IBGKA  = ',IBGKA(JATM,K)
          ELSE IF (NUMSEC .EQ. 4) THEN
            READ (IUNIN,6666) IREACA(JATM,K),IBULKA(JATM,K),
     .                        ISCD1A(JATM,K),ISCD2A(JATM,K),
     .                        ISCD3A(JATM,K),ISCD4A(JATM,K),
     .                        ISCDEA(JATM,K),IESTMA(JATM,K),
     .                        IBGKA(JATM,K)
            WRITE (iunout,*) ' WARNING !!! '
            WRITE (iunout,*)
     .        ' FOUR SECONDARY GROUPS USED FOR SPECIES ',
     .          TEXTS(ISPZ)
            WRITE (iunout,*) ' ISCD1A = ',ISCD1A(JATM,K)
            WRITE (iunout,*) ' ISCD2A = ',ISCD2A(JATM,K)
            WRITE (iunout,*) ' ISCD3A = ',ISCD3A(JATM,K)
            WRITE (iunout,*) ' ISCD4A = ',ISCD4A(JATM,K)
            WRITE (iunout,*) ' ISCDEA = ',ISCDEA(JATM,K)
            WRITE (iunout,*) ' IESTMA = ',IESTMA(JATM,K)
            WRITE (iunout,*) ' IBGKA  = ',IBGKA(JATM,K)
          END IF
          LMULTI = LMULTI .OR. (IBGKA(JATM,K) /= 0)
          READ (IUNIN,6664) EELECA(JATM,K),EBULKA(JATM,K),
     .                      ESCD1A(JATM,K),ESCD2A,
     .                      FREACA(JATM,K),EDPOTA(JATM,K)
cdr  .                      FLDLMA(JATM,K)  removed, now controlled by negative ngena
          ESCD1A(JATM,K) = ESCD1A(JATM,K)+ESCD2A
          NSC = 0
          IF (ISCD3A(JATM,K) > 0) NSC = 3
          IF (ISCD4A(JATM,K) > 0) NSC = 4
          IF (NSC > 0) THEN
            IF ((REACDAT(IREACA(JATM,K))%NOSEC > 0) .AND.
     .          (REACDAT(IREACA(JATM,K))%NOSEC /= NSC)) THEN
              WRITE (IUNOUT,*) ' INCONSISTENCY FOUND CONCERNING',
     .              ' REACTION ',IREACA(JATM,K)
              WRITE (IUNOUT,*) ' NUMBER OF SECONDARIES FOUND',
     .              ' PREVIOUSLY WAS ',REACDAT(IREACA(JATM,K))%NOSEC
              WRITE (IUNOUT,*) ' NUMBER OF SECONDARIES FOUND',
     .              ' NOW IS         ',NSC
              WRITE (IUNOUT,*) ' USE',
     .              MAX(REACDAT(IREACA(JATM,K))%NOSEC,NSC),
     .              ' SECONDARIES'
              REACDAT(IREACA(JATM,K))%NOSEC =
     .              MAX(REACDAT(IREACA(JATM,K))%NOSEC,NSC)
            ELSE
              REACDAT(IREACA(JATM,K))%NOSEC = NSC
            END IF
          END IF
  422   CONTINUE
  421 CONTINUE
C
C  READ NEUTRAL MOLECULES SPECIES CARDS
C
      READ (IUNIN,*)
      WRITE (iunout,*)
     . '*** 4B. NEUTRAL MOLECULE SPECIES CARDS, NMOLI SPECIES'
      READ (IUNIN,*) NMOLI
      WRITE (iunout,*) '       NMOLI= ',NMOLI
      CALL EIRENE_LEER(1)
      NMOLI_IN=NMOLI
      NSPA=NSPH+NATMI
      DO 431 JMOL=1,NMOLI
        ISPZ=NSPA+JMOL
        READ (IUNIN,66666) I,TEXTS(ISPZ),NMASSM(JMOL),NCHARM(JMOL),
     .                       NPRT(ISPZ),NDUMM,   !NCHRGM=0, DEFAULT
     .                       ISRF(ISPZ,1),ISRT(ISPZ,1),NUMSEC,
     .                       NRCM(JMOL),NFOLM(JMOL),NGENM(JMOL),
     .                       NHSTS(ISPZ)
        IF (ISRT(ISPZ,1).LT.0) THEN
          WRITE (iunout,*) 'INPUT ERROR IN BLOCK 4B '
          WRITE (iunout,*) 'MOLECULAR SPECIES ',JMOL,':'
          WRITE (iunout,*) 'ISRT LT 0 OPTION IS NOT AVAILABLE ANYMORE'
          WRITE (iunout,*) 'PROBABLY YOU MEAN: ISRT= ',NMOLI+1
          IERROR=IERROR+1
        ENDIF
        DO 432 K=1,NRCM(JMOL)
          IF (NUMSEC .LT. 3) THEN
            READ (IUNIN,6666) IREACM(JMOL,K),IBULKM(JMOL,K),
     .                        ISCD1M(JMOL,K),ISCD2M(JMOL,K),
     .                        ISCDEM(JMOL,K),IESTMM(JMOL,K),
     .                        IBGKM(JMOL,K)
          ELSE IF (NUMSEC .EQ. 3) THEN
            READ (IUNIN,6666) IREACM(JMOL,K),IBULKM(JMOL,K),
     .                        ISCD1M(JMOL,K),ISCD2M(JMOL,K),
     .                        ISCD3M(JMOL,K),
     .                        ISCDEM(JMOL,K),IESTMM(JMOL,K),
     .                        IBGKM(JMOL,K)
            WRITE (iunout,*) ' WARNING !!! '
            WRITE (iunout,*)
     .        ' THREE SECONDARY GROUPS USED FOR SPECIES ',
     .          TEXTS(ISPZ)
            WRITE (iunout,*) ' ISCD1M = ',ISCD1M(JMOL,K)
            WRITE (iunout,*) ' ISCD2M = ',ISCD2M(JMOL,K)
            WRITE (iunout,*) ' ISCD3M = ',ISCD3M(JMOL,K)
            WRITE (iunout,*) ' ISCDEM = ',ISCDEM(JMOL,K)
            WRITE (iunout,*) ' IESTMM = ',IESTMM(JMOL,K)
            WRITE (iunout,*) ' IBGKM  = ',IBGKM(JMOL,K)
          ELSE IF (NUMSEC .EQ. 4) THEN
            READ (IUNIN,6666) IREACM(JMOL,K),IBULKM(JMOL,K),
     .                        ISCD1M(JMOL,K),ISCD2M(JMOL,K),
     .                        ISCD3M(JMOL,K),ISCD4M(JMOL,K),
     .                        ISCDEM(JMOL,K),IESTMM(JMOL,K),
     .                        IBGKM(JMOL,K)
            WRITE (iunout,*) ' WARNING !!! '
            WRITE (iunout,*)
     .        ' FOUR SECONDARY GROUPS USED FOR SPECIES ',
     .          TEXTS(ISPZ)
            WRITE (iunout,*) ' ISCD1M = ',ISCD1M(JMOL,K)
            WRITE (iunout,*) ' ISCD2M = ',ISCD2M(JMOL,K)
            WRITE (iunout,*) ' ISCD3M = ',ISCD3M(JMOL,K)
            WRITE (iunout,*) ' ISCD4M = ',ISCD4M(JMOL,K)
            WRITE (iunout,*) ' ISCDEM = ',ISCDEM(JMOL,K)
            WRITE (iunout,*) ' IESTMM = ',IESTMM(JMOL,K)
            WRITE (iunout,*) ' IBGKM  = ',IBGKM(JMOL,K)
          END IF
          LMULTI = LMULTI .OR. (IBGKM(JMOL,K) /= 0)
          READ (IUNIN,6664) EELECM(JMOL,K),EBULKM(JMOL,K),
     .                      ESCD1M(JMOL,K),ESCD2M,
     .                      FREACM(JMOL,K),EDPOTM(JMOL,K)
cdr  for backward compatibility:  formerly: two KER values, now one total is used.
          ESCD1M(JMOL,K) = ESCD1M(JMOL,K)+ESCD2M
          NSC = 0
          IF (ISCD3M(JMOL,K) > 0) NSC = 3
          IF (ISCD4M(JMOL,K) > 0) NSC = 4
          IF (NSC > 0) THEN
            IF ((REACDAT(IREACM(JMOL,K))%NOSEC > 0) .AND.
     .          (REACDAT(IREACM(JMOL,K))%NOSEC /= NSC)) THEN
              WRITE (IUNOUT,*) ' INCONSISTENCY FOUND CONCERNING',
     .              ' REACTION ',IREACM(JMOL,K)
              WRITE (IUNOUT,*) ' NUMBER OF SECONDARIES FOUND',
     .              ' PREVIOUSLY WAS ',REACDAT(IREACM(JMOL,K))%NOSEC
              WRITE (IUNOUT,*) ' NUMBER OF SECONDARIES FOUND',
     .              ' NOW IS         ',NSC
              WRITE (IUNOUT,*) ' USE ',
     .              MAX(REACDAT(IREACM(JMOL,K))%NOSEC,NSC),
     .              ' SECONDARIES'
              REACDAT(IREACM(JMOL,K))%NOSEC =
     .              MAX(REACDAT(IREACM(JMOL,K))%NOSEC,NSC)
            ELSE
              REACDAT(IREACM(JMOL,K))%NOSEC = NSC
            END IF
          END IF
  432   CONTINUE
  431 CONTINUE
C
C  READ TEST PARTICLE IONS SPECIES CARDS
C
      READ (IUNIN,*)
      WRITE (iunout,*) '*** 4C. TEST IONS SPECIES CARDS, NIONI SPECIES'
      READ (IUNIN,*) NIONI
      WRITE (iunout,*) '       NIONI= ',NIONI
      CALL EIRENE_LEER(1)
      NIONI_IN=NIONI
      NSPAM=NSPH+NATMI+NMOLI
      DO 441 JION=1,NIONI
        ISPZ=NSPAM+JION
        READ (IUNIN,66666) I,TEXTS(ISPZ),NMASSI(JION),NCHARI(JION),
     .                       NPRT(ISPZ),NCHRGI(JION),
     .                       ISRF(ISPZ,1),ISRT(ISPZ,1),NUMSEC,
     .                       NRCI(JION),NFOLI(JION),NGENI(JION),
     .                       NHSTS(ISPZ)
        DO 442 K=1,NRCI(JION)
          IF (NUMSEC .LT. 3) THEN
            READ (IUNIN,6666) IREACI(JION,K),IBULKI(JION,K),
     .                        ISCD1I(JION,K),ISCD2I(JION,K),
     .                        ISCDEI(JION,K),IESTMI(JION,K),
     .                        IBGKI(JION,K)
          ELSE IF (NUMSEC .EQ. 3) THEN
            READ (IUNIN,6666) IREACI(JION,K),IBULKI(JION,K),
     .                        ISCD1I(JION,K),ISCD2I(JION,K),
     .                        ISCD3I(JION,K),
     .                        ISCDEI(JION,K),IESTMI(JION,K),
     .                        IBGKI(JION,K)
            WRITE (iunout,*) ' WARNING !!! '
            WRITE (iunout,*) ' THREE SECONDARY GROUPS FOR SPECIES ',
     .                    TEXTS(ISPZ)
            WRITE (iunout,*) ' ISCD1I = ',ISCD1I(JION,K)
            WRITE (iunout,*) ' ISCD2I = ',ISCD2I(JION,K)
            WRITE (iunout,*) ' ISCD3I = ',ISCD3I(JION,K)
            WRITE (iunout,*) ' ISCDEI = ',ISCDEI(JION,K)
            WRITE (iunout,*) ' IESTMI = ',IESTMI(JION,K)
            WRITE (iunout,*) ' IBGKI  = ',IBGKI(JION,K)
          ELSE IF (NUMSEC .EQ. 4) THEN
            READ (IUNIN,6666) IREACI(JION,K),IBULKI(JION,K),
     .                        ISCD1I(JION,K),ISCD2I(JION,K),
     .                        ISCD3I(JION,K),ISCD4I(JION,K),
     .                        ISCDEI(JION,K),IESTMI(JION,K),
     .                        IBGKI(JION,K)
            WRITE (iunout,*) ' WARNING !!! '
            WRITE (iunout,*) ' FOUR SECONDARY GROUPS FOR SPECIES ',
     .                    TEXTS(ISPZ)
            WRITE (iunout,*) ' ISCD1I = ',ISCD1I(JION,K)
            WRITE (iunout,*) ' ISCD2I = ',ISCD2I(JION,K)
            WRITE (iunout,*) ' ISCD3I = ',ISCD3I(JION,K)
            WRITE (iunout,*) ' ISCD4I = ',ISCD4I(JION,K)
            WRITE (iunout,*) ' ISCDEI = ',ISCDEI(JION,K)
            WRITE (iunout,*) ' IESTMI = ',IESTMI(JION,K)
            WRITE (iunout,*) ' IBGKI  = ',IBGKI(JION,K)
          END IF
          LMULTI = LMULTI .OR. (IBGKI(JION,K) /= 0)
          READ (IUNIN,6664) EELECI(JION,K),EBULKI(JION,K),
     .                      ESCD1I(JION,K),ESCD2I,
     .                      FREACI(JION,K),EDPOTI(JION,K)
          ESCD1I(JION,K) = ESCD1I(JION,K)+ESCD2I
          NSC = 0
          IF (ISCD3I(JION,K) > 0) NSC = 3
          IF (ISCD4I(JION,K) > 0) NSC = 4
          IF (NSC > 0) THEN
            IF ((REACDAT(IREACI(JION,K))%NOSEC > 0) .AND.
     .          (REACDAT(IREACI(JION,K))%NOSEC /= NSC)) THEN
              WRITE (IUNOUT,*) ' INCONSISTENCY FOUND CONCERNING',
     .              ' REACTION ',IREACI(JION,K)
              WRITE (IUNOUT,*) ' NUMBER OF SECONDARIES FOUND',
     .              ' PREVIOUSLY WAS ',REACDAT(IREACI(JION,K))%NOSEC
              WRITE (IUNOUT,*) ' NUMBER OF SECONDARIES FOUND',
     .              ' NOW IS         ',NSC
              WRITE (IUNOUT,*) ' USE ',
     .              MAX(REACDAT(IREACI(JION,K))%NOSEC,NSC),
     .              ' SECONDARIES'
              REACDAT(IREACI(JION,K))%NOSEC =
     .              MAX(REACDAT(IREACI(JION,K))%NOSEC,NSC)
            ELSE
              REACDAT(IREACI(JION,K))%NOSEC = NSC
            END IF
          END IF
  442   CONTINUE
  441 CONTINUE
C
      READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:3) == '***') THEN
        IREAD=1
        NPHOTI=0
        NPHOTI_IN=NPHOTI
        GOTO 500
      END IF
      IREAD=0
      WRITE (iunout,*)
     . '*** 4D. NEUTRAL PHOTONS SPECIES CARDS, NPHOTI SPECIES'
      READ (IUNIN,*) NPHOTI
      WRITE (iunout,*) '       NPHOTI= ',NPHOTI
      CALL EIRENE_LEER(1)
      NPHOTI_IN=NPHOTI
C
      DO 451 JPHOT=1,NPHOTI
        ISPZ=JPHOT
        READ (IUNIN,66666) I,TEXTS(ISPZ),NDUMM1,NDUMM2,
     .                       NDUMM3,NDUMM4,
     .                       ISRF(ISPZ,1),ISRT(ISPZ,1),NUMSEC,
     .                       NRCPH(JPHOT),NFOLPH(JPHOT),NGENPH(JPHOT),
     .                       NHSTS(ISPZ)
C  DEFAULTS FOR PHOTONIC SPECIES:
        NPRT(ISPZ)=1
        DO 452 K=1,NRCPH(JPHOT)
          IF (NUMSEC .LT. 3) THEN
            READ (IUNIN,6666) IREACPH(JPHOT,K),IBULKPH(JPHOT,K),
     .                        ISCD1PH(JPHOT,K),ISCD2PH(JPHOT,K),
     .                        ISCDEPH(JPHOT,K),IESTMPH(JPHOT,K),
     .                        IBGKPH(JPHOT,K)
          ELSE IF (NUMSEC .EQ. 3) THEN
            READ (IUNIN,6666) IREACPH(JPHOT,K),IBULKPH(JPHOT,K),
     .                        ISCD1PH(JPHOT,K),ISCD2PH(JPHOT,K),
     .                        ISCD3PH(JPHOT,K),
     .                        ISCDEPH(JPHOT,K),IESTMPH(JPHOT,K),
     .                        IBGKPH(JPHOT,K)
            WRITE (iunout,*) ' WARNING !!! '
            WRITE (iunout,*)
     .        ' THREE SECONDARY GROUPS USED FOR SPECIES ',
     .          TEXTS(ISPZ)
            WRITE (iunout,*) ' ISCD1PH = ',ISCD1PH(JPHOT,K)
            WRITE (iunout,*) ' ISCD2PH = ',ISCD2PH(JPHOT,K)
            WRITE (iunout,*) ' ISCD3PH = ',ISCD3PH(JPHOT,K)
            WRITE (iunout,*) ' ISCDEPH = ',ISCDEPH(JPHOT,K)
            WRITE (iunout,*) ' IESTMPH = ',IESTMPH(JPHOT,K)
            WRITE (iunout,*) ' IBGKPH  = ',IBGKPH(JPHOT,K)
          ELSE IF (NUMSEC .EQ. 4) THEN
            READ (IUNIN,6666) IREACPH(JPHOT,K),IBULKPH(JPHOT,K),
     .                        ISCD1PH(JPHOT,K),ISCD2PH(JPHOT,K),
     .                        ISCD3PH(JPHOT,K),ISCD4PH(JPHOT,K),
     .                        ISCDEPH(JPHOT,K),IESTMPH(JPHOT,K),
     .                        IBGKPH(JPHOT,K)
            WRITE (iunout,*) ' WARNING !!! '
            WRITE (iunout,*)
     .        ' FOUR SECONDARY GROUPS USED FOR SPECIES ',
     .          TEXTS(ISPZ)
            WRITE (iunout,*) ' ISCD1PH = ',ISCD1PH(JPHOT,K)
            WRITE (iunout,*) ' ISCD2PH = ',ISCD2PH(JPHOT,K)
            WRITE (iunout,*) ' ISCD3PH = ',ISCD3PH(JPHOT,K)
            WRITE (iunout,*) ' ISCD4PH = ',ISCD4PH(JPHOT,K)
            WRITE (iunout,*) ' ISCDEPH = ',ISCDEPH(JPHOT,K)
            WRITE (iunout,*) ' IESTMPH = ',IESTMPH(JPHOT,K)
            WRITE (iunout,*) ' IBGKPH  = ',IBGKPH(JPHOT,K)
          END IF
          LMULTI = LMULTI .OR. (IBGKPH(JPHOT,K) /= 0)
          READ (IUNIN,6664) EELECPH(JPHOT,K),EBULKPH(JPHOT,K),
     .                      ESCD1PH(JPHOT,K),ESCD2PH,
     .                      FREACPH(JPHOT,K),EDPOTPH(JPHOT,K)
cdr  .                      FLDLMPH(JPHOT,K)  removed. now controlled by negative ngenph
          ESCD1PH(JPHOT,K) = ESCD1PH(JPHOT,K)+ESCD2PH
          NSC = 0
          IF (ISCD3PH(JPHOT,K) > 0) NSC = 3
          IF (ISCD4PH(JPHOT,K) > 0) NSC = 4
          IF (NSC > 0) THEN
            IF ((REACDAT(IREACPH(JPHOT,K))%NOSEC > 0) .AND.
     .          (REACDAT(IREACPH(JPHOT,K))%NOSEC /= NSC)) THEN
              WRITE (IUNOUT,*) ' INCONSISTENCY FOUND CONCERNING',
     .              ' REACTION ',IREACPH(JPHOT,K)
              WRITE (IUNOUT,*) ' NUMBER OF SECONDARIES FOUND',
     .              ' PREVIOUSLY WAS ',REACDAT(IREACPH(JPHOT,K))%NOSEC
              WRITE (IUNOUT,*) ' NUMBER OF SECONDARIES FOUND',
     .              ' NOW IS         ',NSC
              WRITE (IUNOUT,*) ' USE ',
     .              MAX(REACDAT(IREACPH(JPHOT,K))%NOSEC,NSC),
     .              ' SECONDARIES'
              REACDAT(IREACPH(JPHOT,K))%NOSEC =
     .              MAX(REACDAT(IREACPH(JPHOT,K))%NOSEC,NSC)
            ELSE
              REACDAT(IREACPH(JPHOT,K))%NOSEC = NSC
            END IF
          END IF
  452   CONTINUE
  451 CONTINUE
C
C  READ DATA FOR PLASMA BACKGROUND, 500--599
C
  500 CONTINUE
C
      IF (IREAD == 0) READ (IUNIN,'(A72)') ZEILE
C     WRITE (iunout,'(1X,A)') trim(ZEILE)
      CALL EIRENE_MASAGE('*** 5. DATA FOR PLASMA BACKGROUND')
      CALL EIRENE_LEER(1)
C
C  READ BULK IONS SPECIES CARDS
C
      READ (IUNIN,'(A72)') ZEILE
C     WRITE (iunout,'(1X,A)') trim(ZEILE)
      WRITE (iunout,*) '*** 5A. BULK ION SPECIES CARDS, NPLSI SPECIES'
  510 READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1) .EQ. '*') GOTO 510
      READ(ZEILE,6666) NPLSI
      WRITE (iunout,*) '       NPLSI= ',NPLSI
      CALL EIRENE_LEER(1)
      NPLSI_IN=NPLSI
      NSPAMI=NSPAM+NIONI
      NSPTOT=NSPAMI+NPLSI
! counter for additional reaction cards possibly needed for "density models".
      IDMDL = 0
c  loop over background species (except: electrons)
      DO 511 JPLS=1,NPLSI
        ISPZ=NSPAMI+JPLS
        READ (IUNIN,66666) I,TEXTS(ISPZ),NMASSP(JPLS),NCHARP(JPLS),
     .                       NPRT(ISPZ),NCHRGP(JPLS),
     .                       ISRF(ISPZ,1),ISRT(ISPZ,1),NUMSEC,
     .                       NRCP(JPLS),NDUMM1,NDUMM2,
     .                       NHSTS(ISPZ),NDUMM4,
     .                       CDENMODEL(JPLS),NRE
        CALL EIRENE_UPPERCASE (CDENMODEL(JPLS))
        DO 512 K=1,NRCP(JPLS)
          IF (NUMSEC .LT. 3) THEN
            READ (IUNIN,6666) IREACP(JPLS,K),IBULKP(JPLS,K),
     .                        ISCD1P(JPLS,K),ISCD2P(JPLS,K),
     .                        ISCDEP(JPLS,K)
          ELSE IF (NUMSEC .EQ. 3) THEN
            READ (IUNIN,6666) IREACP(JPLS,K),IBULKP(JPLS,K),
     .                        ISCD1P(JPLS,K),ISCD2P(JPLS,K),
     .                        ISCD3P(JPLS,K),ISCDEP(JPLS,K)
            WRITE (iunout,*) ' WARNING !!! '
            WRITE (iunout,*)
     .        ' THREE SECONDARY GROUPS USED FOR SPECIES ',
     .          TEXTS(ISPZ)
            WRITE (iunout,*) ' ISCD1P = ',ISCD1P(JPLS,K)
            WRITE (iunout,*) ' ISCD2P = ',ISCD2P(JPLS,K)
            WRITE (iunout,*) ' ISCD3P = ',ISCD3P(JPLS,K)
            WRITE (iunout,*) ' ISCDEP = ',ISCDEP(JPLS,K)
          ELSE IF (NUMSEC .EQ. 4) THEN
            READ (IUNIN,6666) IREACP(JPLS,K),IBULKP(JPLS,K),
     .                        ISCD1P(JPLS,K),ISCD2P(JPLS,K),
     .                        ISCD3P(JPLS,K),ISCD4P(JPLS,K),
     .                        ISCDEP(JPLS,K)
            WRITE (iunout,*) ' WARNING !!! '
            WRITE (iunout,*)
     .        ' FOUR SECONDARY GROUPS USED FOR SPECIES ',
     .          TEXTS(ISPZ)
            WRITE (iunout,*) ' ISCD1P = ',ISCD1P(JPLS,K)
            WRITE (iunout,*) ' ISCD2P = ',ISCD2P(JPLS,K)
            WRITE (iunout,*) ' ISCD3P = ',ISCD3P(JPLS,K)
            WRITE (iunout,*) ' ISCD4P = ',ISCD4P(JPLS,K)
            WRITE (iunout,*) ' ISCDEP = ',ISCDEP(JPLS,K)
          END IF
          READ (IUNIN,6664) EELECP(JPLS,K),EBULKP(JPLS,K),
     .                      ESCD1P(JPLS,K),ESCD2P,
     .                      FREACP(JPLS,K),EDPOTP(JPLS,K)
          ESCD1P(JPLS,K) = ESCD1P(JPLS,K)+ESCD2P
c
cdr  deal with non-default number of secondaries, NSC > 2
          NSC = 0
          IF (ISCD3P(JPLS,K) > 0) NSC = 3
          IF (ISCD4P(JPLS,K) > 0) NSC = 4
          IF (NSC > 0) THEN
            IF ((REACDAT(IREACP(JPLS,K))%NOSEC > 0) .AND.
     .          (REACDAT(IREACP(JPLS,K))%NOSEC /= NSC)) THEN
              WRITE (IUNOUT,*) ' INCONSISTENCY FOUND CONCERNING',
     .              ' REACTION ',IREACP(JPLS,K)
              WRITE (IUNOUT,*) ' NUMBER OF SECONDARIES FOUND',
     .              ' PREVIOUSLY WAS ',REACDAT(IREACP(JPLS,K))%NOSEC
              WRITE (IUNOUT,*) ' NUMBER OF SECONDARIES FOUND',
     .              ' NOW IS         ',NSC
              WRITE (IUNOUT,*) ' USE ',
     .              MAX(REACDAT(IREACP(JPLS,K))%NOSEC,NSC),
     .              ' SECONDARIES'
              REACDAT(IREACP(JPLS,K))%NOSEC =
     .              MAX(REACDAT(IREACP(JPLS,K))%NOSEC,NSC)
            ELSE

              REACDAT(IREACP(JPLS,K))%NOSEC = NSC
            END IF
cdr  NSC: number of secondaries. But here: 0,1,2 secondaries all have NSC=0 ??
          END IF
c
  512   CONTINUE
        IF (LEN_TRIM(CDENMODEL(JPLS)) > 0) THEN
          ALLOCATE (TDMPAR(JPLS)%TDM)
          TDMPAR(JPLS)%TDM%NRE=MAX(NRE,1)
          ALLOCATE(TDMPAR(JPLS)%TDM%ISP(TDMPAR(JPLS)%TDM%NRE),SOURCE=0)
          ALLOCATE(TDMPAR(JPLS)%TDM%ITP(TDMPAR(JPLS)%TDM%NRE),SOURCE=0)
          ALLOCATE(TDMPAR(JPLS)%TDM%ISTR(TDMPAR(JPLS)%TDM%NRE),SOURCE=0)
          ALLOCATE(TDMPAR(JPLS)%TDM%IRC(TDMPAR(JPLS)%TDM%NRE),SOURCE=0)

C needs trim, check if possible...
          SELECT CASE (CDENMODEL(JPLS))
          CASE (FORT//'13')
            IDMDL = IDMDL + 1
            READ (IUNIN,6666) TDMPAR(JPLS)%TDM%ISP(1)
c  default: only for bulk ions
                              TDMPAR(JPLS)%TDM%ITP(1)=4
          CASE (FORT//'10')
            IDMDL = IDMDL + 1
            READ (IUNIN,6666) TDMPAR(JPLS)%TDM%ISP(1),
     .                        TDMPAR(JPLS)%TDM%ITP(1),
     .                        TDMPAR(JPLS)%TDM%ISTR(1)
          CASE ('CONSTANT')
            IDMDL = IDMDL + 1
            READ (IUNIN,6664) TDMPAR(JPLS)%TDM%TVAL,
     .                        TDMPAR(JPLS)%TDM%DVAL,
     .                        TDMPAR(JPLS)%TDM%VXVAL,
     .                        TDMPAR(JPLS)%TDM%VYVAL,
     .                        TDMPAR(JPLS)%TDM%VZVAL
          CASE ('MULTIPLY')
            IDMDL = IDMDL + 1
            READ (IUNIN,'(3I6,6x,3E12.4)')
     .           TDMPAR(JPLS)%TDM%ISP(1),
     .           TDMPAR(JPLS)%TDM%ITP(1),
     .           TDMPAR(JPLS)%TDM%ISTR(1),
     .           TDMPAR(JPLS)%TDM%DFACTOR,
     .           TDMPAR(JPLS)%TDM%TFACTOR,
     .           TDMPAR(JPLS)%TDM%VFACTOR
                 TDMPAR(JPLS)%TDM%ITP(1)=4
          CASE ('SAHA')
            IDMDL = IDMDL + 1
!PB   TO BE WRITTEN
          CASE ('BOLTZMANN')
            IDMDL = IDMDL + 1
            READ (IUNIN,'(3I6,6x,2E12.4)')
     .           TDMPAR(JPLS)%TDM%ISP(1),
     .           TDMPAR(JPLS)%TDM%ITP(1),
     .           TDMPAR(JPLS)%TDM%ISTR(1),
     .           TDMPAR(JPLS)%TDM%G_BOLTZ,
     .           TDMPAR(JPLS)%TDM%DELTAE
          CASE ('CORONA')
            IDMDL = IDMDL + 1  !  ONE MORE H.2 REACTION data set
            READ (IUNIN,'(A72)') ZEILE
            IF (VERIFY(ZEILE,'1234567890eEdD+- ') > 0) THEN
              WRITE (IUNOUT,*) ' OLD VERSION OF INPUT FOR DENSITY',
     .                         ' MODEL CORONA USED'
              WRITE (IUNOUT,*) ' NEW INPUT FORMAT:',
     .                         ' REACTIONS SPECIFIED IN BLOCK 4 ONLY '
              WRITE (IUNOUT,*) ' CHECK IN THE EIRENE MANUAL '
              CALL EIRENE_EXIT_OWN(1)
            END IF
            READ (ZEILE,'(4I6,E12.4)')
     .           TDMPAR(JPLS)%TDM%ISP(1),
     .           TDMPAR(JPLS)%TDM%ITP(1),
     .           TDMPAR(JPLS)%TDM%ISTR(1),
     .           TDMPAR(JPLS)%TDM%IRC(1),
     .           TDMPAR(JPLS)%TDM%A_CORONA
            IR = TDMPAR(JPLS)%TDM%IRC(1)
            IF (.NOT.REACDAT(IR)%LRTC) THEN
              WRITE (iunout,*)
     .          ' WRONG REACTION SPECIFIED FOR CORONA MODEL'
              WRITE (iunout,*) ' ONLY H.2 REACTION DATA ARE PERMITTED'
              WRITE (iunout,*) ' IPLS = ',JPLS
              CALL EIRENE_EXIT_OWN(1)
            END IF
          CASE ('COLRAD')
            IDMDL = IDMDL + 1  !  ONE MORE H.11 or H.12 REACTION data set
            DO I=1, TDMPAR(JPLS)%TDM%NRE
              READ (IUNIN,'(A72)') ZEILE
              IF (VERIFY(ZEILE,'1234567890eEdD+- ') > 0) THEN
                WRITE (IUNOUT,*) ' OLD VERSION OF INPUT FOR DENSITY',
     .                           ' MODEL COLRAD USED'
                WRITE (IUNOUT,*) ' NEW INPUT FORMAT:',
     .                           ' REACTIONS SPECIFIED IN BLOCK 4 ONLY '
                WRITE (IUNOUT,*) ' CHECK IN THE EIRENE MANUAL '
                CALL EIRENE_EXIT_OWN(1)
              END IF
              READ (ZEILE,6666)
     .             TDMPAR(JPLS)%TDM%ISP(I),
     .             TDMPAR(JPLS)%TDM%ITP(I),
     .             TDMPAR(JPLS)%TDM%ISTR(I),
     .             TDMPAR(JPLS)%TDM%IRC(I)

              IR = TDMPAR(JPLS)%TDM%IRC(I)
              IF (.NOT.REACDAT(IR)%LOTH) THEN
                WRITE (iunout,*)
     .            ' WRONG REACTION SPECIFIED FOR COLRAD MODEL'
                WRITE (iunout,*) ' ONLY H.11 OR H.12 (OT) REACTION DATA'
                WRITE (IUNOUT,*) ' ARE PERMITTED'
                WRITE (iunout,*) ' IPLS = ',JPLS
                CALL EIRENE_EXIT_OWN(1)
              END IF
            END DO
          CASE DEFAULT
!PB  NOTHING TO BE DONE
          END SELECT
        END IF
  511 CONTINUE   ! nplsi
cdr  additional reaction cards due to density models: idmdl
cdr   nreaci=nreaci=idmdl  ??
      CALL EIRENE_LEER(1)
C
      DO I=1,NSPZ
        CALL EIRENE_UPPERCASE (TEXTS(I))
      END DO

C
  520 READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1) .EQ. '*') THEN
C       WRITE (iunout,'(1x,a)') trim(ZEILE)
        GOTO 520
      ENDIF
      CALL EIRENE_MASAGE('*** 5B. PLASMA BACKGROUND DATA')
      CALL EIRENE_LEER(1)
      READ (ZEILE,6666) (INDPRO(J),J=1,12)

!     SAVE VALUES READ FROM INPUT FILE
      INDPRO_IN = INDPRO

      DO J=1,12
C  Indicate "smoothed" input tallies (interpolation into cells)
C  amongst the 12 input tallies read here.
C  formerly: LSMOPRO(J) flag.
C  Smoothing is currently only for tallies 1 to 7.
        IF (ABS(INDPRO(J)) > 100) THEN
          INDPRO(J) = MOD(INDPRO(J),100)
          SELECT CASE(J)
          CASE (1)
            LTESMO = .TRUE.
          CASE (2)
            LTISMO = .TRUE.
          CASE (3)
            LDESMO = .TRUE.
            LDISMO = .TRUE.
          CASE (4)
            LVXSMO = .TRUE.
            LVYSMO = .TRUE.
            LVZSMO = .TRUE.
          CASE (5)
            LBXSMO = .TRUE.
            LBYSMO = .TRUE.
            LBZSMO = .TRUE.
            LBFSMO = .TRUE.
          CASE (6)
            LADSMO = .TRUE.
          CASE (7)
            LEXSMO = .TRUE.
            LEYSMO = .TRUE.
            LEZSMO = .TRUE.
            LEFSMO = .TRUE.
            LPOTSMO = .TRUE.     ! smoothed electric potential
          END SELECT
        END IF
      ENDDO

C  Te profile
      IF (INDPRO(1).LE.5.AND.NPLSI.GT.0)
     .  READ (IUNIN,6664) TE0,TE1,TE2,TE3,TE4,TE5

C  Ti profile(s)
!pb   NPLSTI = 1
!pb   IF (LMULTI .OR. (INDPRO(2) < 0) .OR. (INDPRO(2) > 9)) NPLSTI=NPLS

      NPLSTI = NPLS
      IF (INDPRO(2) < 0) NPLSTI=1

      IF ((NPLS > 1) .AND. (NPLSTI == 1)) THEN
        WRITE (IUNOUT,*) 'WARNING !'
        WRITE (IUNOUT,*) 'TIIN PROVIDED FOR ONE SPECIES ONLY ',
     .                   'DUE TO INDPRO(2) < 0'
        IF (LMULTI) THEN
          WRITE (IUNOUT,*) 'DIMENSION OF TIIN OVERWRITTEN ',
     .                     'BECAUSE BGK REACTIONS PRESENT'
          NPLSTI = NPLS
        END IF
        WRITE (IUNOUT,*) ' NPLSTI = ',NPLSTI
      END IF

!pb   ion temperature for each bulk ion
      NLMLTI = (NPLSTI > 1)
      MPLSTI=1
      IF (NLMLTI)     MPLSTI = (/ (I,I=1,NPLS) /)
!pb      IF (NPLSTI > 1) MPLSTI = (/ (I,I=1,NPLS) /)

      INDPRO(2)=IABS(INDPRO(2))
!pb      NLMLTI= INDPRO(2) > 9
!pb   read parameters TI0...TI5 for each temperature
      LRDMLTI= INDPRO(2) > 9
      IF (INDPRO(2) > 9) INDPRO(2) = MOD(INDPRO(2),10)

      IF (INDPRO(2).LE.5.AND.NPLSI.GT.0) THEN
!pb        IF (NLMLTI) THEN
        IF (LRDMLTI) THEN
          READ (IUNIN,6664) (TI0(I),TI1(I),TI2(I),TI3(I),TI4(I),TI5(I),
     .                       I=1,NPLSI)
        ELSE
C  ONLY ONE COMMON ION TEMPERATURE FOR ALL SPECIES
          READ (IUNIN,6664)  TI0(1),TI1(1),TI2(1),TI3(1),TI4(1),TI5(1)
          DO 530 JPLS=2,NPLSTI
            TI0(JPLS)=TI0(1)
            TI1(JPLS)=TI1(1)
            TI2(JPLS)=TI2(1)
            TI3(JPLS)=TI3(1)
            TI4(JPLS)=TI4(1)
            TI5(JPLS)=TI5(1)
  530     CONTINUE
        ENDIF
      ENDIF

c  di profiles
      IF (INDPRO(3).LE.5)
     .  READ (IUNIN,6664) (DI0(I),DI1(I),DI2(I),DI3(I),DI4(I),DI5(I),
     .                     I=1,NPLSI)
c  v_in profile(s)
cdr  default is: cm/s units for flow field(s)
      NLMACH=INDPRO(4).LT.0  ! Mach number units instead, rather than cm/s
      INDPRO(4)=IABS(INDPRO(4))

cdr  default is: |indpro| < 10:  npls flow fields, one per background species
cdr              |indpro| > 10:  only one common flow field
      NPLSV = NPLS
      IF (INDPRO(4) > 9) NPLSV = 1

      IF (INDPRO(4) > 9) INDPRO(4) = MOD(INDPRO(4),10)
      NLMLV = (NPLSV > 1)
      IF (.NOT.NLMLV) THEN
        MPLSV = 1                    ! find all flow fields on ipls=1 storage
      ELSE
        MPLSV = (/ (I,I=1,NPLS) /)   ! find each individual flow field on its ipls storage
      ENDIF
      IF (INDPRO(4).LE.5.AND.NPLSI.GT.0) THEN
        IF (.NOT. NLMLV) THEN
C  ONLY ONE COMMON FLOW VELOCITY FOR ALL NPLS SPECIES
          READ (IUNIN,6664) VX0(1),VX1(1),VX2(1),VX3(1),VX4(1),VX5(1)
          READ (IUNIN,6664) VY0(1),VY1(1),VY2(1),VY3(1),VY4(1),VY5(1)
          READ (IUNIN,6664) VZ0(1),VZ1(1),VZ2(1),VZ3(1),VZ4(1),VZ5(1)
        ELSE
C  READ NPLSI SETS OF INPUT PARAMETERS, ONE FOR EACH SPECIES IPLS
          READ (IUNIN,6664) (VX0(I),VX1(I),VX2(I),VX3(I),VX4(I),VX5(I),
     .                       I=1,NPLSI)
          READ (IUNIN,6664) (VY0(I),VY1(I),VY2(I),VY3(I),VY4(I),VY5(I),
     .                       I=1,NPLSI)
          READ (IUNIN,6664) (VZ0(I),VZ1(I),VZ2(I),VZ3(I),VZ4(I),VZ5(I),
     .                       I=1,NPLSI)
        END IF
      ENDIF

c  pitch or B-field profile
c                              !  default:                B-FIELD WITH BX=0
      NLPITCH=INDPRO(5).LT.0   !  for 1D parallel B runs: B-FIELD WITH BY=0

      INDPRO(5)=IABS(INDPRO(5))
      IF (INDPRO(5).LE.5)
     .  READ (IUNIN,6664) B0,B1,B2,B3,B4,B5

c  cell volume -profile card, OPTIONAL

      IF (INDPRO(12).LE.5) THEN
        READ (IUNIN,'(A72)',IOSTAT=IO) ZEILE
        CALL EIRENE_UPPERCASE(ZEILE)
        IREAD=1
        IF ((IO /= 0) .OR. (ZEILE(1:3) .EQ. '***')
     .                .OR. (ZEILE(1:3) .EQ. 'OPT')) THEN
          WRITE (iunout,*) 'ONE INPUT LINE (VOL) MISSING IN BLOCK 5'
          WRITE (iunout,*) 'AUTOMATIC CORRECTION PERFORMED'
          VL0=0
        ELSE
          READ (ZEILE,6664) VL0,VL1,VL2,VL3,VL4,VL5
          IREAD=0
        ENDIF
      ENDIF

c  check for "OPTIONAL" input cards.
c  Here: explicitly switch off input tallies,
c        or activate gradient tallies: INTLOPTS(ITAL)
c  ITAL:  -ntali<ital<0, then: ital=iabs(ital).
c  IOPT:   turn off or add gradient tally

      IO = 0
      IF (IREAD == 0) READ (IUNIN,'(A72)',IOSTAT=IO) ZEILE
      IREAD = 1
      IF ((IO == 0) .AND. (ZEILE(1:3) .NE. '***')) THEN
        CALL EIRENE_UPPERCASE(ZEILE)
        IF (INDEX(ZEILE,'OPTIONAL') > 0) THEN
          DO 
            READ (IUNIN,'(A72)',IOSTAT=IO) ZEILE
            IF ((IO /= 0) .OR. (ZEILE(1:3) .EQ. '***')) EXIT
            READ (ZEILE,6666) ITAL, IOPT
            IF ((ITAL >= 0) .OR. (ITAL < -NTALI)) CYCLE         
            INTLOPTS(IABS(ITAL)) = IOPT
          END DO
        END IF
      END IF

      IF (LINCL45) THEN
        CLOSE (IUNIN)
        IUNIN = IUNIN_SAVE
        LINCL45 =.FALSE.

        DO
          READ (IUNIN,'(A72)') ZEILE
          IF ((ZEILE(1:3) == '***') .AND.
     ,        (INDEX(ZEILE,'6.') > 0)) EXIT
        END DO
        IREAD = 1
      END IF

C  TRY TO GENERATE INPUT FOR EIRENE CORRESPONDING TO A HYDKIN RUN.
C  NOT READY
      IF (LHYDDEF) THEN
        WRITE (IUNOUT,*)
     .    'LHYDDEF IS TRUE: OPTION NOT READY, EXIT CALLED'
        CALL EIRENE_EXIT_OWN(1)
cdr  subr. EIRENE_SETUP_HYDKIN_REACTIONS moved to folder: unfinished_business
cdr     CALL EIRENE_SETUP_HYDKIN_REACTIONS(HYDKIN_DEFAULT,CADAPT)
      ENDIF
C
C  READ DATA FOR REFLECTION MODEL  600--699
C
      IF (IREAD.EQ.0) READ (IUNIN,*)
      IREAD=0
      NULLIFY(REFFILES)
      CALL EIRENE_MASAGE('*** 6. GENERAL DATA FOR REFLECTION MODEL')
      CALL EIRENE_LEER(1)
C
  610 READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1) .EQ. '*') GOTO 610
      READ (ZEILE,6665) NLTRIM
      IREAD=0
      NFR = 0
      IF (NLTRIM) THEN

c  read TRIM reflection datasets A_on_B
        READ (IUNIN,'(A420)') ZEILE
        IREAD=1
        IF (INDEX(ZEILE,'PATH')+INDEX(ZEILE,'path').EQ.0) THEN
C  NO PATH SPECIFIED FOR REFLECTION DATABASE
          WRITE (iunout,*)
     .      ' NO PATH SPECIFIED FOR TRIM REFLECTION DATABASE'
          WRITE (iunout,*) ' OLD TRIM DATABASE VERSION USED'
          CALL EIRENE_LEER(1)
          LTRIM_OLD=.TRUE.
C  TAKE OLD "TRIM.DAT" FILE WITH NHD6=12 TARGET-PROJECTILE COMBINATIONS
C           "TRIM.DAT" IS EXPECTED ON INPUT STREAM IUN=(21+IFOFF) (SUBR. REFDAT.F)
C  SKIP LINES CONTAINING SPECIFICATIONS FOR DATABASES AND
C  CONTINUE READING with DATD
  615     IF (INDEX(ZEILE,'ON')+INDEX(ZEILE,'on').NE.0) THEN
            READ (IUNIN,'(A72)') ZEILE
            GOTO 615
          ELSE
C  NEXT VALID INPUT CARD FOUND
            IREAD=1
            GOTO 620
          ENDIF
        ELSE
C  PATH SPECIFICATION FOR DATABASE FOUND
          LTRIM_OLD=.FALSE.
          READ (ZEILE(7:LEN(ZEILE)),'(A400)') PATH
          PATH=ADJUSTL(PATH)
          IF (LEN_TRIM(PATH) == 0) PATH='./ '
          I2=INDEX(PATH,' ')
          IF ((I2 == 0) .AND. (ZEILE(408:408) /= ' ')) THEN
            WRITE (iunout,*) ' PATH FOR TRIM REFLECTION DATABASE IS',
     .                  ' TOO LONG !'
            CALL EIRENE_EXIT_OWN(1)
          END IF
          FILE=REPEAT(' ',420)
          FILE(1:I2-1)=PATH(1:I2-1)
          WRITE (iunout,*) ' PATH = ',FILE(1:I2-1)
C  PATH FOUND. NEXT: READ ONE OR MORE CARDS FILNAM A_ON_B
C         NFR=0  !dr now already set above
          READ (IUNIN,'(A72)') ZEILE
  625     IF (INDEX(ZEILE,'ON')+INDEX(ZEILE,'on').NE.0) THEN
            NFR=NFR+1
            READ (ZEILE,'(A72)') RFILNM
            FILE(I2:)=RFILNM
            WRITE (iunout,*) ' NFR =',NFR,' FILNAM = ',FILE(I2:I2+20)
C           WRITE (iunout,'(A,A)') ' FILE = ',FILE
            ALLOCATE (CURFILE)
            CURFILE%RFILE = FILE
            CURFILE%NEXT => REFFILES
            REFFILES => CURFILE
            READ (IUNIN,'(A72)') ZEILE
            GOTO 625
          ELSE
            IREAD=1
            GOTO 620
          ENDIF
        ENDIF
      ENDIF

  620 CONTINUE  !  READING OF REFLECTION DATASETS 'A_ON_B' COMPLETED

      IF (ASSOCIATED(REFFILES)) THEN
c  TRIM files for NFR target-projectile combinations are requested.
c  read them one by one in subr. RDTRIM
        NHD6 = NFR
      ELSE
c  old default: all TRIM data in one single file.
c  read this (single) file "TRIM.DAT" in subr. REFDAT
        NHD6 = 12
      END IF

      NH0=NHD1*NHD2*NHD6
      NH1=NH0*NHD3
      NH2=NH1*NHD4
      NH3=NH2*NHD5
      CALL EIRENE_ALLOC_CREF

      IFLR = 0
      DO WHILE (ASSOCIATED(REFFILES))
        IFLR = IFLR + 1
        REFFIL(IFLR) = REFFILES%RFILE
        CURFILE => REFFILES
        REFFILES => REFFILES%NEXT
        DEALLOCATE(CURFILE)
      END DO
cdr AT THIS POINT: IFLR=NFR OR IFLR=0 ??
      IF (IFLR.NE.NFR.AND.IFLR.NE.0) GOTO 993

c  next: read species index sampling distributions datm, dmol, dion, dpls, and in case nphot > 0, also dphot

      IF (IREAD.EQ.0) READ (IUNIN,'(A72)') ZEILE
      READ (ZEILE,6664) (DATD(JATM),JATM=1,MIN(NATMI_IN,6))
      IREAD=0
      IF (NATMI_IN.GT.6) THEN
        READ (IUNIN,6664) (DATD(JATM),JATM=7,NATMI_IN)
      ENDIF
      READ (IUNIN,6664) (DMLD(JMOL),JMOL=1,NMOLI_IN)
      READ (IUNIN,6664) (DIOD(JION),JION=1,NIONI_IN)
      READ (IUNIN,6664) (DPLD(JPLS),JPLS=1,NPLSI_IN)
      IF (NPHOTI > 0)   !dr try to make this more logic: always read dphd.
cdr                     !dr backward compatible ?
     .  READ (IUNIN,6664) (DPHD(JPHOT),JPHOT=1,NPHOTI_IN)
      
      DATD_IN = DATD
      DMLD_IN = DMLD
      DIOD_IN = DIOD
      DPLD_IN = DPLD
      DPHD_IN = DPHD
      
c  next: read universal surface reflection model flags
      READ (IUNIN,6664) ERMIN,ERCUT,RPROB0,RINTEG(1),EINTEG(1),AINTEG(1)

c  read surface models identified by character string 'SURFMOD_...'

      DO
        IF (IREAD.EQ.0) READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:3) .EQ. '***') THEN
          IREAD = 1
          EXIT
        END IF
        IF (ZEILE(1:8) /= 'SURFMOD_') CYCLE
        ALLOCATE (REFCUR)
        ALLOCATE (REFCUR%JSRS(nspz))
        ALLOCATE (REFCUR%JSRC(nspz))
        ALLOCATE (REFCUR%TRANSPR(nspz,2))
        ALLOCATE (REFCUR%RCYCFR(nspz))
        ALLOCATE (REFCUR%RCYCTR(nspz))
        ALLOCATE (REFCUR%RCPRMR(nspz))
        ALLOCATE (REFCUR%EXPPLR(nspz))
        ALLOCATE (REFCUR%EXPELR(nspz))
        ALLOCATE (REFCUR%EXPILR(nspz))
        ALLOCATE (REFCUR%RCYCSR(nspz))
        ALLOCATE (REFCUR%RCYCCR(nspz))
        ALLOCATE (REFCUR%STPRMR(nspz))
        ALLOCATE (REFCUR%ESPTSR(nspz))
        ALLOCATE (REFCUR%ESPTCR(nspz))
        NULLIFY (REFCUR%SPEC_LINES)

        REFCUR%JSRS = 0
        REFCUR%JSRC = 0
        REFCUR%TRANSPR = 0._DP
        REFCUR%RCYCFR  = 0._DP
        REFCUR%RCYCTR  = 0._DP
        REFCUR%RCPRMR  = 0._DP
        REFCUR%EXPPLR  = 0._DP
        REFCUR%EXPELR  = 0._DP
        REFCUR%EXPILR  = 0._DP
        REFCUR%RCYCSR  = 0._DP
        REFCUR%RCYCCR  = 0._DP
        REFCUR%STPRMR  = 0._DP
        REFCUR%ESPTSR  = 0._DP
        REFCUR%ESPTCR  = 0._DP

        REFCUR%REFNAME = TRIM(ADJUSTL(ZEILE(9:)))
        IREAD=0
        READ (IUNIN,6666) REFCUR%JLREF,REFCUR%JLSPT,
     .                    REFCUR%JSRS(1),REFCUR%JSRC(1)
        READ (IUNIN,6664) REFCUR%ZNMLR,REFCUR%EWALLR,REFCUR%EWBINR,
     .                    REFCUR%TRANSPR(1,1),
     .                    REFCUR%TRANSPR(1,2),REFCUR%FSHEATR
        READ (IUNIN,6664) REFCUR%RCYCFR(1),REFCUR%RCYCTR(1),
     .                    REFCUR%RCPRMR(1),REFCUR%EXPPLR(1),
     .                    REFCUR%EXPELR(1),REFCUR%EXPILR(1)
        DO I=2,NSPZ
          REFCUR%JSRS(I) = REFCUR%JSRS(1)
          REFCUR%JSRC(I) = REFCUR%JSRC(1)
          REFCUR%TRANSPR(I,1:2)=REFCUR%TRANSPR(1,1:2)
          REFCUR%RCYCFR(I) = REFCUR%RCYCFR(1)
          REFCUR%RCYCTR(I) = REFCUR%RCYCTR(1)
          REFCUR%RCPRMR(I) = REFCUR%RCPRMR(1)
          REFCUR%EXPPLR(I) = REFCUR%EXPPLR(1)
          REFCUR%EXPELR(I) = REFCUR%EXPELR(1)
          REFCUR%EXPILR(I) = REFCUR%EXPILR(1)
        END DO
C
C  DEFAULT
        REFCUR%RCYCSR=RECYCS(1,0)
        REFCUR%RCYCCR=RECYCC(1,0)
        REFCUR%STPRMR=SPTPRM(1,0)
        REFCUR%ESPTSR=ESPUTS(1,0)
        REFCUR%ESPTCR=ESPUTC(1,0)
C
        READ (IUNIN,'(A72)') ZEILE
        IREAD=1
        ideflt_sput=-1
        ideflt_spez=-1
C  READ ONE MORE LINE FOR NON-DEFAULT SPUTTER MODEL
        do while (ZEILE(1:1).NE.'*'.AND.ZEILE(1:8).NE.'SURFMOD_')
          varname = zeile(1:8)
          call EIRENE_uppercase (varname)
          spcname = zeile(10:17)
          call EIRENE_uppercase (spcname)
          ispz=-1
          ico=0
          do is=1,nspz
            if (texts(is) == spcname) then
              if (ico == 0) then
                ispz = is
                ico = 1
              else
                call EIRENE_leer(2)
                write (iunout,*) ' WARNING !!'
                write (iunout,*)
     .           ' ambiguous species names found in surface model "',
     .                       REFCUR%REFNAME,'"'
                write (iunout,*)
     .           varname,' is modified for species ',ispz
              end if
            end if
          end do     ! end of "is"-loop
          ideflt_spez=1   !tentatively assume: species card
          if (ispz > 0) then
            allocate(spr)
            spr%varname = varname
            spr%spcname = spcname
            if ((varname == 'ISRS') .or. (varname == 'ISRC')) then           
              read (zeile(18:),'(I6)') spr%ival
            else
              read (zeile(18:),'(E12.4)') spr%rval
            end if
            spr%next => refcur%spec_lines
            refcur%spec_lines => spr
          end if
          select case (varname)
          case ('ISRS')
            if (ispz > 0)
     .          read (zeile(18:),'(I6)') REFCUR%JSRS(ispz)
          case ('ISRC')
            if (ispz > 0)
     .          read (zeile(18:),'(I6)') REFCUR%JSRC(ispz)
          case ('TRANSP1')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%TRANSPR(ispz,1)
          case ('TRANSP2')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%TRANSPR(ispz,2)
          case ('RECYCF')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%RCYCFR(ispz)
          case ('RECYCT')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%RCYCTR(ispz)
          case ('RECPRM')
            if (ispz >0)
     .          read (zeile(18:),'(E12.4)') REFCUR%RCPRMR(ispz)
          case ('EXPPL')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%EXPPLR(ispz)
          case ('EXPEL')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%EXPELR(ispz)
          case ('EXPIL')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%EXPILR(ispz)
          case ('RECYCS')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%RCYCSR(ispz)
          case ('RECYCC')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%RCYCCR(ispz)
          case ('SPTPRM')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%STPRMR(ispz)
          case ('ESPUTS')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%ESPTSR(ispz)
          case ('ESPUTC')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%ESPTCR(ispz)

          case default
c  not a species card, hence: a sputer model card
            if (ispz < 0) then
              READ (ZEILE,6664) REFCUR%RCYCSR(1),REFCUR%RCYCCR(1),
     .                          REFCUR%STPRMR(1),REFCUR%ESPTSR(1),
     .                          REFCUR%ESPTCR(1)
              DO I=2,NSPZ
                REFCUR%RCYCSR(I) = REFCUR%RCYCSR(1)
                REFCUR%RCYCCR(I) = REFCUR%RCYCCR(1)
                REFCUR%STPRMR(I) = REFCUR%STPRMR(1)
                REFCUR%ESPTSR(I) = REFCUR%ESPTSR(1)
                REFCUR%ESPTCR(I) = REFCUR%ESPTCR(1)
              ENDDO
              ideflt_sput=1
              ideflt_spez=-1
            end if
          end select
          if ((ideflt_spez > 0) .and. (ispz < 0)) then
            write (iunout,*)
     .        ' wrong card in species dep. reflection model'
            write (iunout,'(1x,a)') trim(zeile)
          end if
          READ (IUNIN,'(A72)') ZEILE
          IREAD=1
C
        END do    ! end of do-while loop (search for * or for SURFMOD_)
C all lines of this particular SURFMOD_... are read now
C
        IF (ideflt_sput.le.0.and.REFCUR%JLSPT.NE.0) THEN
          WRITE (iunout,*) 'WARNING: SPUTTERING FOR MODEL ',
     .                      REFCUR%REFNAME
          WRITE (iunout,*) 'BUT NO PARAMETERS RECYCS, RECYCC ARE READ'
          WRITE (iunout,*) 'DEFAULT MODEL: "NO SPUTTERING" IS USED.'
          WRITE (iunout,*) 'DO YOU REALLY WANT THIS?'
          REFCUR%JLSPT=0
        ENDIF
        L=LEN_TRIM(REFCUR%REFNAME)
        WRITE (iunout,*) 'SURFACE MODEL "',REFCUR%REFNAME(1:L),
     .                   '" DEFINED'
!pb     REFCUR%NEXT => REFLIST
!pb     REFLIST => REFCUR
!  add new reflection model at the end of the list
          nullify(refcur%next)
          if (.not.associated(reflist)) then
            reflist => refcur
            reflast => refcur
          else
            reflast%next => refcur
            reflast => reflast%next
          end if
      END DO
      CALL EIRENE_LEER(1)
C
C  READ DATA FOR PRIMARY SOURCE  700--799
C
      IF (IREAD.EQ.0) READ (IUNIN,*)
      CALL EIRENE_MASAGE
     .  ('*** 7. DATA FOR PRIMARY SOURCES, NSTRAI STRATA')
C
  710 READ (IUNIN,'(A72)') ZEILE
      IREAD=1
      IF (ZEILE(1:1) .EQ. '*') GOTO 710
      READ (ZEILE,6666) NSTRAI
      IREAD=0
      WRITE (iunout,*) '       NSTRAI= ',NSTRAI
      CALL EIRENE_LEER(1)
      READ (IUNIN,6666) (INDSRC(IST),IST=1,NSTRAI)
      READ (IUNIN,6664) ALLOC, AMPTS

C  AMPTS: (option added 2014) common multiplier for max. allowed cpu time NTCPU,
C          and for number of histories NPTS (see below)
      IF (AMPTS.EQ.0.0_DP) AMPTS=1.0_DP
      IF(AMPTS.NE.1.0_DP) THEN
        MPTS_COMSOU=AMPTS
        NTCPU=INT(REAL(NTCPU)*AMPTS)
        WRITE (iunout,*) '       NTCPU, NPTS ENHANCED BY FACTOR AMPTS ',
     .                     AMPTS
      ELSE
        MPTS_COMSOU=1.0
      END IF
C
      DO 712 ISTRA=1,NSTRAI
        IF (INDSRC(ISTRA).EQ.6) GOTO 712  ! SKIP READING INPUT FOR THIS STRATUM.

        I=ISTRA
        IF (IREAD.EQ.0) THEN
          READ (IUNIN,'(A72)') TXTSOU(ISTRA)
        ELSEIF (IREAD.EQ.1) THEN
          READ (ZEILE,'(A72)') TXTSOU(ISTRA)
          IREAD=0
        ENDIF
  713   READ (IUNIN,'(A72)') ZEILE
        IREAD=1
        IF (ZEILE(1:1) .EQ. '*') GOTO 713
        call fix_logical_input(zeile,4)
        READ (ZEILE,6665) NLAVRP(I),NLAVRT(I),NLSYMP(I),NLSYMT(I)
C    .                   ,NLRAY(I)
                          NLRAY(I)=.FALSE.  ! CDR: UNFINISHED PROPRIETARY OPTION
        IREAD=0
        READ (IUNIN,6666) NPTS(ISTRA),NINITL(I),NEMODS(I),NAMODS(I),

c.............................................
c  further proprietary input flags. Not ready for external use.
c  NMINPTS:   ENFORCE MINIMUM NUMBER OF HISTORIES FOR STRATUM.
cdr           currently set in various places, but not used anywhere.
     .                    NMINPTS(ISTRA),
c  NPTSDEL: =0   (VALUES NPTSDEL .GT.0 ONLY FOR INTERNAL TESTTING PROCEDURES
C                 OF MPI parallelization. Requires 2 runs with special input settings
cdr See comments in MCARLO.F
     .                    NPTSDEL(ISTRA)
C    .                   ,NRAYEN(ISTRA)   !CDR  NOT IN USE
c...................................................................

C APPLY MULTIPLIER AMPTS TO SELECTED MC PARTICLE NUMBER RANGE
        IF(AMPTS.GT.0) THEN
         NPTS(ISTRA)=INT(NPTS(ISTRA)*AMPTS)
         NMINPTS(ISTRA)=INT(NMINPTS(ISTRA)*AMPTS)
        END IF
C

csw if npts < 0 --> npts = infinity
        if(npts(istra).lt.0) npts(istra) = huge(1)-1

        READ (IUNIN,66662) FLUX(ISTRA),SCALV(ISTRA),IVLSF(ISTRA),
     .                    ISCLS(ISTRA),ISCLT(ISTRA),ISCL1(ISTRA),
     .                    ISCL2(ISTRA),ISCL3(ISTRA),ISCLB(ISTRA),
     .                    ISCLA(ISTRA)
        READ (IUNIN,'(A72)') ZEILE
        call fix_logical_input(zeile,5)
        READ (ZEILE,6665) NLATM(ISTRA),NLMOL(ISTRA),NLION(I),NLPLS(I),
     .                    NLPHOT(ISTRA)
        READ (IUNIN,6666) NSPEZ(ISTRA)
        READ (IUNIN,'(A72)') ZEILE
        call fix_logical_input(zeile,5)
        READ (ZEILE,6665) NLPNT(ISTRA),NLLNE(ISTRA),
     .                    NLSRF(ISTRA),NLVOL(ISTRA),NLCNS(ISTRA)
C  SAME FOR POINT, LINE, SURFACE AND VOLUME SOURCES
        READ (IUNIN,6666) NSRFSI(ISTRA)
        IF (NSRFSI(ISTRA).GT.NSRFS)
     .    CALL EIRENE_MASPRM(
     .         'NSRFS',5,NSRFS,'NSRFSI(I)',9,NSRFSI(I),IERROR)
        DO 715 J=1,NSRFSI(ISTRA)
          READ (IUNIN,6666) INUM,INDIM(J,I),INSOR(J,I),
     .                      INGRDA(J,I,1),INGRDE(J,I,1),
     .                      INGRDA(J,I,2),INGRDE(J,I,2),
     .                      INGRDA(J,I,3),INGRDE(J,I,3)
          READ (IUNIN,6664) SORWGT(J,I),SORLIM(J,I),
     .                      SORIND(J,I),SOREXP(J,I),SORIFL(J,I)
          READ (IUNIN,6666) NRSOR(J,I),NPSOR(J,I),NTSOR(J,I),
     .                      NBSOR(J,I),NASOR(J,I),NISOR(J,I),ISTOR(J,I)
          READ (IUNIN,6664) SORAD1(J,I),SORAD2(J,I),SORAD3(J,I),
     .                      SORAD4(J,I),SORAD5(J,I),SORAD6(J,I)
  715   CONTINUE
C  VELOCITY SPACE DISTRIBUTION
        READ (IUNIN,6664) SORENI(ISTRA),SORENE(I),SORVDX(I),SORVDY(I),
     .                    SORVDZ(ISTRA)
        READ (IUNIN,6664) SORCOS(I),SORMAX(I),SORCTX(I),SORCTY(I),
     .                    SORCTZ(ISTRA),RAYFRAC(ISTRA)
!  STORE VALUES FOR JSON OUTPUT
        SORCOS_IN(I) = SORCOS(I) 
        SORMAX_IN(I) = SORMAX(I) 
C
  712 CONTINUE
C
      READ (IUNIN,*)
C     READ ADDITIONAL DATA FOR SOME SPECIFIC ZONES
C
      CALL EIRENE_MASAGE ('*** 8. ADDITIONAL DATA FOR SPECIFIC ZONES')
  810 READ (IUNIN,'(A72)') ZEILE
      IREAD=1
      IF (ZEILE(1:1) .EQ. '*') GOTO 810
      READ (ZEILE,6666) NZADD
      IREAD=0
      WRITE (iunout,*) '       NZADD= ',NZADD
      CALL EIRENE_LEER(1)
      NULLIFY(TEMPLIST)
      NULLIFY(DENLIST)
      NULLIFY(VELLIST)
      NULLIFY(VOLLIST)
      ALLOCATE (INI_ZONE(NZADD))
      ALLOCATE (INE_ZONE(NZADD))
      ALLOCATE (CH3_STACK(NZADD))
      INI_ZONE = 0
      INE_ZONE = 0
      DO 811 I=1,NZADD
        NULLIFY(CH3_STACK(I)%HEAD)
        NULLIFY(CH3_STACK(I)%LAST)
CDR  skip reading optional comment lines
  814   READ (IUNIN,'(A72)') ZEILE
        IREAD=1
        IF (ZEILE(1:1) .EQ. '*') GOTO 814
CDR
        READ (ZEILE,6666) INI,INE
        INI_ZONE(I) = INI
        INE_ZONE(I) = INE
        IREAD=0
        IF (INI.GT.NRAD.OR.INI.LE.0) GOTO 998
        IF (INE.GT.NRAD) GOTO 998
        IF (INE.LE.0) INE=INI
  812   READ (IUNIN,'(A72)') ZEILE
C  IGJUM3 FLAG
        IF (ZEILE(1:3).EQ.'CH3') THEN
          INILGJ=MAX(1,MIN(NOPTIM,INI))
          if (inilgj.ne.ini) write (iunout,*)
     .      'CH3 option failed. ini, noptim =',ini,NOPTIM
          INELGJ=MIN(NOPTIM,INE)
          if (inelgj.ne.ine) write (iunout,*)
     .      'CH3 option failed. ine, noptim =',ine,NOPTIM
          DO 821 IN=INILGJ,INELGJ
            IF (NLIMPB >= NLIMPS) THEN
              CALL EIRENE_DEKEY (ZEILE(4:72),IGJUM3,0,NOPTIM,IN,NLIMPS)
            ELSE
              CALL EIRENE_DEKEYB
     .                    (ZEILE(4:72),IGJUM3,0,NOPTIM,IN,NLIMPB,NBITS)
            END IF
  821     CONTINUE
          call eirene_push_string_stack(ch3_stack(i),zeile(1:72))
          GOTO 812
C  TEMPERATURE
        ELSEIF (ZEILE(1:1).EQ.'T') THEN
          DO 822 IN=INI,INE
            ALLOCATE(TEMPCUR)
            READ (ZEILE(2:72),66664) TEMPCUR%IDION,TEMPCUR%TE,
     .                               TEMPCUR%TI
            TEMPCUR%IN = IN
            TEMPCUR%NEXT => TEMPLIST
            TEMPLIST => TEMPCUR
  822     CONTINUE
          GOTO 812
C  DENSITY
        ELSEIF (ZEILE(1:1).EQ.'D') THEN
          DO 823 IN=INI,INE
            ALLOCATE(DENCUR)
            READ (ZEILE(2:72),66664) DENCUR%IDION,DENCUR%DI
            DENCUR%IN = IN
            DENCUR%NEXT => DENLIST
            DENLIST => DENCUR
  823     CONTINUE
          GOTO 812
C  VELOCITY (CM/SEC OR MACH)
        ELSEIF ((ZEILE(1:1).EQ.'V'.AND.ZEILE(2:2).NE.'L')
     .       .OR.ZEILE(1:1).EQ.'M') THEN
          DO 824 IN=INI,INE
            ALLOCATE(VELCUR)
            READ (ZEILE(2:72),66664) VELCUR%IDION,VELCUR%VX,VELCUR%VY,
     .                               VELCUR%VZ
            IF (ZEILE(1:1).EQ.'M') THEN
              VELCUR%IZ=1
            ELSEIF (ZEILE(1:1).EQ.'V') THEN
              VELCUR%IZ=-1
            ENDIF
            VELCUR%IN = IN
            VELCUR%NEXT => VELLIST
            VELLIST => VELCUR
  824     CONTINUE
          GOTO 812
C  VOLUME
        ELSEIF (ZEILE(1:2).EQ.'VL') THEN
          DO 825 IN=INI,INE
            ALLOCATE(VOLCUR)
            READ (ZEILE(3:72),6664) VOLCUR%VOL
            VOLCUR%IN = IN
            VOLCUR%NEXT => VOLLIST
            VOLLIST => VOLCUR
  825     CONTINUE
          GOTO 812
        ELSEIF (ZEILE(1:1).EQ.'*') THEN
          IREAD=1
          GOTO 811
        ELSE
          GOTO 998
        ENDIF
  811 CONTINUE

C
C  READ DATA FOR STATISTICS AND NON-ANALOG MODEL, 900--999
C
      IF (IREAD.EQ.0) READ (IUNIN,*)
      CALL EIRENE_MASAGE
     .  ('*** 9. DATA FOR STATISTIC AND NON-ANALOG MODEL')
C
  910 IF (IREAD.EQ.0) READ (IUNIN,'(A72)') ZEILE
      IREAD = 0
      IF (ZEILE(1:1) .EQ. '*') GOTO 910
C  DATA FOR CONDITIONAL EXPECTATION ESTIMATOR
      NLOGIN = MAX(1,NATMI_IN + NMOLI_IN + NIONI_IN + NPHOTI_IN)
C  read 60 or more logical flags
      ALLOCATE (LOGRDH(NLOGIN))
      LOGRDH = .FALSE.
      DO J=1, NLOGIN, 60
        call fix_logical_input(zeile,60)
        READ (ZEILE,6665) LOGRDH(J:MIN(J+59,NLOGIN))
        READ (IUNIN,'(A72)') ZEILE
      END DO
c
      IF (NATMI_IN > 0) NLPRCA(1:NATMI_IN) = LOGRDH(1:NATMI_IN)
      IF (NMOLI_IN > 0)
     .    NLPRCM(1:NMOLI_IN) = LOGRDH(NATMI_IN+1 : NATMI_IN+NMOLI_IN)
      IF (NIONI_IN > 0) NLPRCI(1:NIONI_IN) =
     .       LOGRDH(NATMI_IN+NMOLI_IN+1 : NATMI_IN+NMOLI_IN+NIONI_IN)
      IF (NPHOTI_IN > 0) NLPRCPH(1:NPHOTI_IN) =
     .       LOGRDH(NATMI_IN+NMOLI_IN+NIONI_IN+1 : NLOGIN)
      DEALLOCATE (LOGRDH)
C  READ NON-DEFAULT OR ADDITIONAL SURFACES, TOWARDS WHICH COND. EXP. EST. IS USED
C  DEFAULT: ???
      READ (ZEILE,6666) NPRCSF

      NPRCSF=MIN0(NLIMPS,NPRCSF)
c  read 12 or more integer flags
      IPRCSF=1
  911 CONTINUE
      IF (IPRCSF.LE.NPRCSF) THEN
        READ (IUNIN,6666) (IPRSF(J),J=1,12)
        DO J=1,12
        IF (IPRSF(J).GT.0.AND.IPRSF(J).LE.NLIMPS)
     .     NLPRCS(IPRSF(J))=.TRUE.
        ENDDO
        IPRCSF=IPRCSF+12
        GOTO 911
      ENDIF

C  DATA FOR SPLITTING AND RUSSIAN ROULETTE
      READ (IUNIN,6666) MAXLEV,MAXRAD,MAXPOL,MAXTOR,MAXADD
      MXL=15
      IF (MAXLEV.GT.MXL)
     .    CALL EIRENE_MASPRM('MXL',3,MXL,'MAXLEV',6,MAXLEV,IERROR)
      IF (ABS(MAXRAD).GT.N1ST)
     .    CALL EIRENE_MASPRM('N1ST',4,N1ST,'MAXRAD',6,MAXRAD,IERROR)
      IF (MAXPOL.GT.N2ND)
     .    CALL EIRENE_MASPRM('N2ND',4,N2ND,'MAXPOL',6,MAXPOL,IERROR)
      IF (MAXTOR.GT.N3RD)
     .    CALL EIRENE_MASPRM('N3RD',4,N3RD,'MAXTOR',6,MAXTOR,IERROR)
      IF (MAXADD.GT.NLIM)
     .    CALL EIRENE_MASPRM('NLIM',4,NLIM,'MAXADD',6,MAXADD,IERROR)
      DO IN=1,MAXRAD
        READ (IUNIN,66665) ID,NSSPL(IN),PRMSPL(IN)
      ENDDO
      DO IN=1,MAXPOL
        READ (IUNIN,66665) ID,NSSPL(N1ST+IN),PRMSPL(N1ST+IN)
      ENDDO
      DO IN=1,MAXTOR
        READ (IUNIN,66665) ID,NSSPL(N1ST+N2ND+IN),PRMSPL(N1ST+N2ND+IN)
      ENDDO
      DO IN=1,MAXADD
        READ (IUNIN,66665) ID,NSSPL(N1ST+N2ND+N3RD+IN),
     .                        PRMSPL(N1ST+N2ND+N3RD+IN)
      ENDDO
C  DATA FOR BIAS SAMPLING
      READ (IUNIN,6664) WMINV,WMINS,WMINC,WMINL
      WMINV=MAX(WMINV,EPS60)
      WMINS=MAX(WMINS,EPS60)
      WMINC=MAX(WMINC,EPS60)
      READ (IUNIN,6664) SPLPAR
C  DATA FOR STANDARD DEVIATION
      READ (IUNIN,*)
      CALL EIRENE_LEER(1)
      WRITE (iunout,*) '       CARDS FOR STANDARD DEVIATION'
      READ (IUNIN,6666) NSIGVI,NSIGSI,NSIGCI,NSIGI_SPC
      WRITE (iunout,*) '       NSIGVI,NSIGSI,NSIGCI = ',
     .                         NSIGVI,NSIGSI,NSIGCI
      WRITE (iunout,*) '       NSIGI_SPC            = ',NSIGI_SPC
      CALL EIRENE_LEER(1)
      DO 913 J=1,NSIGVI
        READ (IUNIN,6666) IGH(J),IIH(J)
  913 CONTINUE
      DO 914 J=1,NSIGSI
        READ (IUNIN,6666) IGHW(J),IIHW(J)
  914 CONTINUE
      DO 915 J=1,NSIGCI
        READ (IUNIN,6666) IGHC(1,J),IIHC(1,J),IGHC(2,J),IIHC(2,J)
  915 CONTINUE
C
C   READ DATA FOR ADDITIONAL AND SURFACE-AVERAGED TALLIES
      READ (IUNIN,*)
      CALL EIRENE_MASAGE
     .  ('*** 10. DATA FOR ADDITIONAL TALLIES, COLLISION')
      CALL EIRENE_MASAGE
     .  ('        ESTIMATORS AND ALGEBRAIC EXPRESSIONS')
C
 1010 READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1) .EQ. '*') GOTO 1010
      READ (ZEILE,6666) NADVI,NCLVI,NALVI,NADSI,NALSI,NADSPC
      CALL EIRENE_LEER(1)
      WRITE (iunout,*) '       NADVI,NCLVI,NALVI    = ',
     .                         NADVI,NCLVI,NALVI
      WRITE (iunout,*) '       NADSI,NALSI,NADSPC   = ',
     .                         NADSI,NALSI,NADSPC
      CALL EIRENE_LEER(1)
C
C
      READ (IUNIN,*)
      CALL EIRENE_MASAGE('*** 10A. DATA FOR ADDITIONAL TALLIES')

      DO 1020 J=1,NADVI
 1021   READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 1021
        READ (ZEILE,6666) IADVE(J),IADVS(J),IADVT(J),IADRC(J)
        READ (IUNIN,'(A72)') TXTTLA(J)
        READ (IUNIN,'(2A24)') TXTSCA(J),TXTUTA(J)
 1020 CONTINUE

      READ (IUNIN,*)
      CALL EIRENE_MASAGE
     .     ('*** 10B. DATA FOR COLLISION ESTIMATORS')
      DO 1030 J=1,NCLVI
 1031   READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 1031
        READ (ZEILE,6666) ICLVE(J),ICLVS(J),ICLVT(J),ICLRC(J)
        READ (IUNIN,'(A72)') TXTTLC(J)
        READ (IUNIN,'(2A24)') TXTSCC(J),TXTUTC(J)
 1030 CONTINUE
      READ (IUNIN,*)
      CALL EIRENE_MASAGE('*** 10C. DATA FOR ALGEBRAIC EXPRESSIONS')
      DO 1040 J=1,NALVI
 1041   READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 1041
        READ (ZEILE,'(A72)') CHRTAL(J)
        READ (IUNIN,'(A72)') TXTTLR(J)
        READ (IUNIN,'(2A24)') TXTSCR(J),TXTUTR(J)
 1040 CONTINUE
C
      READ (IUNIN,*)
      CALL EIRENE_MASAGE
     .  ('*** 10D. DATA FOR ADDITIONAL SURFACE TALLIES')
      DO 1050 J=1,NADSI
 1051   READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 1051
        READ (ZEILE,6666) IADSE(J),IADSS(J),IADST(J),IADSC(J)
        READ (IUNIN,'(A72)') TXTTLW(J,NTLSA)
        READ (IUNIN,'(2A24)') TXTSPW(J,NTLSA),TXTUNW(J,NTLSA)
 1050 CONTINUE
C
      READ (IUNIN,*)
      CALL EIRENE_MASAGE('*** 10E. DATA FOR ALGEBRAIC SURFACE TALLIES')
      DO 1060 J=1,NALSI
 1061   READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 1061
        READ (ZEILE,'(A72)') CHRTLS(J)
        READ (IUNIN,'(A72)') TXTTLW(J,NTLSR)
        READ (IUNIN,'(2A24)') TXTSPW(J,NTLSR),TXTUNW(J,NTLSR)
 1060 CONTINUE
C
      READ (IUNIN,'(A72)') ZEILE
      IREAD=1
      IF (ZEILE(1:3) == '** ') THEN
        CALL EIRENE_MASAGE('*** 10F. DATA FOR SPECTRA')
        IF (NADSPC > 0) THEN
          ALLOCATE(ESTIML(NADSPC))
          IF (NSMSTRA > 0) ALLOCATE(SMESTL(NADSPC))
        ELSE
          ALLOCATE(ESTIML(1))
        END IF
C
        DO J=1,NADSPC
          DO
            READ (IUNIN,'(A72)') ZEILE
            IF (ZEILE(1:1) .NE. '*') EXIT
          END DO
cdr  this input and this code segment should be modified:
c  first read isrfcll: "type of spectrum", then ispsrf "where is this spectrum scored"
c  then further details (species: iptyp,ipspz, and then: direction,....)
c

          READ (ZEILE,'(12I6)') ISPSRF, IPTYP, IPSPZ, ISPTYP, NSPS,
     .                          ISRFCLL, IDIREC
          READ (IUNIN,'(6E12.4)') SPCMN, SPCMX, SPC_SHIFT,
     .                            SPCPLT_X, SPCPLT_Y, SPCPLT_SAME
          SPCVX = 0._DP
          SPCVY = 0._DP
          SPCVZ = 0._DP
          IF (IDIREC /= 0) THEN
            READ (IUNIN,'(6E12.4)') SPCVX, SPCVY, SPCVZ
            VNORM = SQRT(SPCVX**2+SPCVY**2+SPCVZ**2)+EPS60
            SPCVX = SPCVX / VNORM
            SPCVY = SPCVY / VNORM
            SPCVZ = SPCVZ / VNORM
            IF (ISRFCLL == 0) THEN
              WRITE (IUNOUT,*) ' SPECTRUM NUMBER ',J
              WRITE (IUNOUT,*) ' DEFINITION OF LINE OF SIGHT FOR',
     .                         ' SURFACE SPECTRUM IS NOT FORESEEN'
              WRITE (IUNOUT,*) ' SPECIFICATION OF DIRECTION IS',
     .                         ' IGNORED'
              IDIREC = 0
            END IF
          END IF

cdr  better: first discriminate by isrfcll,  then, for each value of isrfcll: do the rest
c   isrcfll=0 : surface-averaged tally
c   isrfcll=1 : volume-averaged tally, integrated over all directions
c   isrfcll=2 : volume-averaged tally, along a specific direction
c  it seems: surface-averaged directional tallies: not yet forseen
c  one may try to score coefficients of orthogonal angular expansion, per energy bin.

          IF (ISPSRF > 0) THEN
            IF ((ISRFCLL == 0) .AND. (ISPSRF > NLIMI)) THEN
C  SPECTRUM AT AN ADDITIONAL SURFACE
              WRITE (iunout,*)
     .          ' SURFACE INDEX FOR SPECTRUM OUT OF BOUNDS'
              WRITE (iunout,*) ' SPECTRUM NUMBER = ',J
              WRITE (iunout,*) ' SURFACE NUMBER  = ',ISPSRF
              IERROR = IERROR + 1
            END IF
C  SPECTRUM IN CELL, POSSIBLY ALONG A CERTAIN DIRECTION
cdr:  next 2 lines: why different condition for cell based spectra cell numbers ??
            IF (((ISRFCLL == 1) .AND. (ISPSRF > NRTAL)) .OR.   ! cell based spectrum in scoring cell ISPSRF
     .          ((ISRFCLL == 2) .AND. (ISPSRF > NRAD))) THEN   ! directional cell based spectrum in scoring cell ISPSRF
              WRITE (iunout,*)
     .          ' CELL INDEX FOR SPECTRUM OUT OF BOUNDS'
              WRITE (iunout,*) ' SPECTRUM NUMBER = ',J
              WRITE (iunout,*) ' CELL NUMBER     = ',ISPSRF
              IERROR = IERROR + 1
            END IF
          ELSEIF (ISPSRF < 0) THEN
            ISPSRF = ABS(ISPSRF)
            IF ((ISRFCLL == 0) .AND. (ISPSRF > NSTSI)) THEN
C  SPECTRUM AT A NON-DEFAULT STANDARD SURFACE
              WRITE (iunout,*)
     .          ' SURFACE INDEX FOR SPECTRUM OUT OF BOUNDS'
              WRITE (iunout,*) ' SPECTRUM NUMBER = ',J
              WRITE (iunout,*) ' SURFACE NUMBER  = ',ISPSRF
              IERROR = IERROR + 1
            ELSEIF (ISRFCLL == 0) THEN
              ISPSRF = NLIM+ISPSRF
            ELSE
C           IF (((ISRFCLL == 1) .AND. (ISPSRF < 0)) .OR.
C    .          ((ISRFCLL == 2) .AND. (ISPSRF < 0))) THEN
C    .      ... WRONG INPUT !
            END IF
          ELSEIF (ISPSRF == 0) THEN
            WRITE (iunout,*) ' SURFACE OR CELL INDEX = 0: NOT FORESEEN'
            WRITE (iunout,*) ' SPECTRUM NUMBER = ',J
            IERROR = IERROR + 1
          END IF

          IF ((IPTYP < 0) .OR. (IPTYP > 4)) THEN
            WRITE (iunout,*) ' PARTICLE TYPE ',IPTYP,' NOT FORESEEN'
            WRITE (iunout,*) ' SPECTRUM NUMBER = ',J
            IERROR = IERROR + 1
          ELSE
            IF (((IPTYP == 0)
     .        .AND.((IPSPZ < 0).OR.(IPSPZ > NPHOTI))) .OR.
     .          ((IPTYP == 1)
     .        .AND.((IPSPZ < 0).OR.(IPSPZ > NATMI))) .OR.
     .          ((IPTYP == 2)
     .        .AND.((IPSPZ < 0).OR.(IPSPZ > NMOLI))) .OR.
     .          ((IPTYP == 3)
     .        .AND.((IPSPZ < 0).OR.(IPSPZ > NIONI))) .OR.
     .          ((IPTYP == 4)
     .        .AND.((IPSPZ < 0).OR.(IPSPZ > NPLSI)))) THEN
              WRITE (iunout,*) ' PARTICLE SPECIES INDEX OUT OF BOUNDS'
              WRITE (iunout,*) ' SPECTRUM NUMBER = ',J
              WRITE (iunout,*) ' SPECIES NUMBER  = ',IPSPZ
              IERROR = IERROR + 1
            END IF
          END IF

          IF ((ISRFCLL == 0) .AND. (ISPTYP > 2)) THEN
            WRITE (IUNOUT,*) ' WRONG TYPE OF SPECTRUM SPECIFIED'
            WRITE (IUNOUT,*) ' SPECTRUM NUMBER = ',J
            WRITE (IUNOUT,*) ' SPECTRUM TYPE   = ',ISPTYP
          END IF

          ESPEC => ESTIML(J)

          ESPEC%ISPCSRF = ISPSRF
          ESPEC%IPRTYP = IPTYP
          ESPEC%IPRSP = IPSPZ
          ESPEC%ISPCTYP = ISPTYP
          ESPEC%NSPC = ABS(NSPS)
          ESPEC%ISRFCLL = ISRFCLL
          ESPEC%IDIREC = IDIREC
cdr
cdr       ESPEC%LOG = .FALSE. ! this was too restrictive !
cdr  Option     LOG = .TRUE. WAS ALREADY AVAILABLE IN SCORING/UPDATE_SPECTRUM

cdr   X.B. correction Sept 17, from SOLPS-ITER branch,
          IF (NSPS > 0) THEN
            NSPSA=NSPS
            ESPEC%LOG = .FALSE.
            ESPEC%SPCMIN=SPCMN
            ESPEC%SPCMAX=SPCMX
            ESPEC%SPCDEL=(SPCMX-SPCMN)/REAL(NSPSA,DP)
          ELSEIF (NSPS < 0) THEN
            NSPSA=-NSPS
            ESPEC%LOG = .TRUE.
cdr  missing here: check for SPCMN, SPCMX > 0.
cdr  In particular in case IDIREC > 0 scores may fall into negative bins.
cdr  tbd: Rule out combination of IDIREC > 0 and NSPS < 0
            ESPEC%SPCMIN = log10(SPCMN)
            ESPEC%SPCMAX = log10(SPCMX)
            ESPEC%SPCDEL = log10(SPCMX/SPCMN)/REAL(NSPSA,DP)
          ELSE
            WRITE (iunout,*) ' SPECTRUM TALLY NUMBER = ',J
            WRITE (iunout,*) ' ZERO NO. OF SPECTRAL BINS "NSPS". EXIT '
            CALL EIRENE_EXIT_OWN(1)
          ENDIF
          ESPEC%ESP_00=SPC_SHIFT
          ESPEC%SPC_XPLT=SPCPLT_X
          ESPEC%SPC_YPLT=SPCPLT_Y
          ESPEC%SPC_SAME=SPCPLT_SAME
          ESPEC%SPCVX=SPCVX
          ESPEC%SPCVY=SPCVY
          ESPEC%SPCVZ=SPCVZ
C  MIN AND MAX ENERGY SCORE (EV) ON THIS TALLY
          ESPEC%ESP_MIN=1.E30
          ESPEC%ESP_MAX=-1.E30

          ESPEC%SPCDELI=1._DP/(ESPEC%SPCDEL+EPS60)
          ALLOCATE(ESPEC%SPC(0:NSPSA+1))
c  standard deviation of spectrally resolved tallies
          IF (NSIGI_SPC > 0) THEN
            ALLOCATE(ESPEC%SDV(0:NSPSA+1))
            ALLOCATE(ESPEC%SGM(0:NSPSA+1))
            ALLOCATE(ESPEC%STV(0:NSPSA+1))
            ALLOCATE(ESPEC%GG(0:NSPSA+1))
          ELSE
            NULLIFY(ESPEC%SDV, ESPEC%SGM, ESPEC%STV, ESPEC%GG)
          END IF
          ESPEC%SPC(0:NSPSA+1) = 0._DP
          ESPEC%IMETSP = 0

          IF (NSMSTRA > 0) THEN
c  sum over strata
            SSPEC => SMESTL(J)
            ALLOCATE(SSPEC%SPC(0:NSPSA+1))
c  standard deviation of spectra tallies, sum over strata intermediate storage
            IF (NSIGI_SPC > 0) THEN
              ALLOCATE(SSPEC%SDV(0:NSPSA+1))
              ALLOCATE(SSPEC%SGM(0:NSPSA+1))
              ALLOCATE(SSPEC%STV(0:NSPSA+1))
              ALLOCATE(SSPEC%GG(0:NSPSA+1))
            ELSE
              NULLIFY(SSPEC%SDV, SSPEC%SGM, SSPEC%STV, SSPEC%GG)
            END IF
            SMESTL(J) = ESTIML(J)
          END IF
C
        END DO
        IREAD=0
      ELSE
Cdr  No block 10F found, i.e. no spectrally resolved tallies at all.
cdr  Why do we allocate estiml in input.f and not in eirmod_cestim ?
        IF (.NOT.ALLOCATED(ESTIML)) THEN
          ALLOCATE(ESTIML(1))
        END IF
      END IF
C
C   READ DATA FOR NUMERICAL AND GRAPHICAL OUTPUT 1100--1199
C
      IF (IREAD == 0) READ (IUNIN,*)
      CALL EIRENE_LEER(1)
      CALL EIRENE_MASAGE
     .  ('*** 11. DATA FOR NUMERICAL AND GRAPHICAL OUTPUT')
c
c  search for input block 11a
 1110 READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1) .EQ. '*') GOTO 1110
* For gfortran: it does not accept empty field for logical
      call fix_logical_input(zeile,35)
      READ (ZEILE,6665) TRCPLT,TRCHST,TRCNAL,TRCMOD,TRCSIG,
     .                  TRCGRD,TRCSUR,TRCREF,TRCFLE,TRCAMD,
     .                  TRCINT,TRCLST,TRCSOU,TRCREC,TRCTIM,
     .                  TRCBLA,TRCBLM,TRCBLI,TRCBLP,TRCBLE,
     .                  TRCBLPH,TRCTAL,TRCOCT,TRCCEN,TRCRNF,
CVK TRACING FOR DEBUGGING, V.Kotov:  not in use in present EIRENE version
cdr  .                  TRCDBG2,TRCDBGE,TRCDBGM,TRCDBGF,TRCDBGL,
cdr  .                  TRCDBGS,TRCDBGG,TRCDBGMPI,TRCDBGC,
CPB  ACTIVATE SPECIES-RESOLVED CPU CONSUMPTION OPTION
     .                  TRCHKTIM
      do I = 0, NSTRA, 60
        READ (IUNIN,'(A72)') ZEILE
        call fix_logical_input(zeile,60)
        READ (ZEILE,6665) (TRCSRC(J),J=I,MIN(I+59,NSTRA))
      end do
C
      READ (IUNIN,6666) NVOLPR, NSPCPR
      CALL EIRENE_LEER(1)
      IF(NVOLPR > NVLPR) THEN
         WRITE (iunout,*) 'NVOLPR > NVLPR, EXIT.'
         WRITE (iunout,*) ' NVOLPR, NVLPR ', NVOLPR, NVLPR
         CALL EIRENE_EXIT_OWN(1)
      ENDIF
      WRITE (iunout,*) '        NVOLPR= ',NVOLPR
      DO 1120 J=1,NVOLPR
        READ (IUNIN,6666) NTLV,NFLGV,NSPZV1,NSPZV2,NTLVF
        IF (NTLV.LT.-NTALI.OR.NTLV.GT.NTALV) GOTO 990
cdr  output stream for particular tallies: allow only ntlv>=70
        IF ((NTLVF.NE.0).AND.(NTLVF.LT.70)) GOTO 995
        NPRTLV(J)  =NTLV
        NFLAGV(J)  =NFLGV
        NSPEZV(J,1)=NSPZV1
        NSPEZV(J,2)=NSPZV2
        NTLVFL(J)  =NTLVF
 1120 CONTINUE
C
      READ (IUNIN,6666) NSURPR
      IF(NSURPR > NSRPR) THEN
         WRITE (iunout,*) 'NSURPR > NSRPR, EXIT.'
         WRITE (iunout,*) ' NSURPR, NSRPR ', NSURPR, NSRPR
         CALL EIRENE_EXIT_OWN(1)
      ENDIF
C
      WRITE (iunout,*) '        NSURPR= ',NSURPR
      DO 1130 J=1,NSURPR
        READ (IUNIN,6666) NSRF,NTLS,NFLGS,NSPZS1,NSPZS2,NTLSF
c  negative surface numbers may be used for non-default standard surfaces
        IF (NSRF.LT.0) NSRF=NLIM+IABS(NSRF)
c  zero as surface number for time horizon
        IF (NSRF.EQ.0.AND.NLIM+NSTSI.LT.NLIMPS) NSRF=NLIM+NSTSI+1

        IF (NSRF.LE.0.OR.NSRF.GT.NLIMPS) GOTO 991
        NPRSRF(J)=NSRF
        NPRTLS(J)=NTLS
        NFLAGS(J)=NFLGS
        NSPEZS(J,1)=NSPZS1
        NSPEZS(J,2)=NSPZS2
        NTLSFL(J)=NTLSF
 1130 CONTINUE

c  overrule default switching on/off of volume-averaged tallies

 1131 READ (IUNIN,'(A72)') ZEILE
      CALL EIRENE_UPPERCASE(ZEILE)
      IREAD=1
      IF ((ZEILE(1:1) .NE.'*') .AND. (SCAN(ZEILE,'FT') == 0)) THEN

c  there are ntlvout volume tallies to be dealt with
c  and ntlsout surface tallies

        READ (ZEILE,6666) NVTLOUT
        IREAD=0
        ITLV=0
        IF (NVTLOUT > 0) THEN
          ALLOCATE (IVTLOUT(NVTLOUT))
          IVTLOUT = 0
        END IF
        DO WHILE ((NVTLOUT > 0) .AND. (ITLV < NVTLOUT))
          READ (IUNIN,6666) (NUMTAL(J),J=1,12)
          DO J=1,12
            IF ((NUMTAL(J) /= 0) .AND. (ABS(NUMTAL(J)) <= NTALV)) THEN
              IF (NUMTAL(J) < 0) THEN
                LMISTALV(ABS(NUMTAL(J))) = .TRUE.  ! SWITCHED OFF
              ELSE
                LMISTALV(NUMTAL(J)) = .FALSE.      ! SWITCHED ON
              END IF
              IVTLOUT(ITLV+J) = NUMTAL(J)
            END IF
          END DO
          ITLV=ITLV+12
        END DO

c  overrule default switching on/off of surface-averaged tallies

        READ (IUNIN,6666) NSTLOUT
        ITLS=0
        IF (NSTLOUT > 0) THEN
          ALLOCATE (ISTLOUT(NSTLOUT))
          ISTLOUT = 0
        END IF
        DO WHILE ((NSTLOUT > 0) .AND. (ITLS < NSTLOUT))
          READ (IUNIN,6666) (NUMTAL(J),J=1,12)
          DO J=1,12
            IF ((NUMTAL(J) /= 0) .AND. (ABS(NUMTAL(J)) <= NTALS)) THEN
              IF (NUMTAL(J) < 0) THEN
                LMISTALS(ABS(NUMTAL(J))) = .TRUE.  ! SWITCHED OFF
              ELSE
                LMISTALS(NUMTAL(J)) = .FALSE.      ! SWITCHED ON
              END IF
              ISTLOUT(ITLS+J) = NUMTAL(J)
            END IF
          END DO
          ITLS=ITLS+12
        END DO

      ELSE
C EITHER A COMMENT LINE, OR THE NEXT INPUT CARD (LOGICALS FOR PLOTTING) IS ON 'zeile'
        IF (ZEILE(1:1) .eq.'*')  goto 1131
c  reading 'switch tallies off' done
      END IF

C  search for input block 11b

 1132 IF (IREAD == 0) READ (IUNIN,'(A72)') ZEILE
      IREAD=0
      IF (ZEILE(1:1) .EQ. '*') GOTO 1132
C  2D GEOMETRY PLOT
      call fix_logical_input(zeile,16)
      READ (ZEILE,6665) PL1ST,PL2ND,PL3RD,PLADD,PLHST,
     .                  PLCUT(1),PLCUT(2),PLCUT(3),PLBOX,PLSTOR,
     .                  PLNUMV,PLNUMS,PLARR,LRPSCUT,PLIDL,
     .                  PLVTK
      READ (IUNIN,6666) NPLINR,NPLOTR,NPLDLR,NPLINP,NPLOTP,NPLDLP,
     .                  NPLINT,NPLOTT,NPLDLT
C  3D GEOMETRY PLOT
      DO 1140 J=1,5
        READ (IUNIN,6662) PL3A(J),TEXTLA(J),IPLTA(J),
     .                (IPLAA(J,I),IPLEA(J,I),I=1,IPLTA(J))
 1140 CONTINUE
      DO 1141 J=1,3
        READ (IUNIN,6662) PL3S(J),TEXTLS(J),IPLTS(J),
     .                (IPLAS(J,I),IPLES(J,I),I=1,IPLTS(J))
 1141 CONTINUE
C
      READ (IUNIN,6664) CH2MX,CH2MY,      CH2X0,CH2Y0,CH2Z0
      READ (IUNIN,6664) CH3MX,CH3MY,CH3MZ,CH3X0,CH3Y0,CH3Z0
      READ (IUNIN,6664) ANGLE1,ANGLE2,ANGLE3
C
C  PARTICLE HISTORY PLOTS IN 2D OR 3D GEOMETRY PLOTS
      READ (IUNIN,6666) I1TRC,I2TRC,(ISYPLT(J),J=1,8),ILINIE
C

c  search for input block 11c

 1151 READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1) .EQ. '*') GOTO 1151
C  DATA FOR PLOTS OF VOLUME-AVERAGED TALLIES
      READ (ZEILE,6666) NVOLPL
      WRITE (iunout,*) '        NVOLPL= ',NVOLPL
      CALL EIRENE_LEER(1)
C
      IF (NVOLPL.GT.NPTAL)
     .  CALL EIRENE_MASPRM('NPTAL',5,NPTAL,'NVOLPT',6,NVOLPL,IERROR)
C
      IF (NVOLPL.LE.0) GOTO 1175
      DO I = 0, NSTRA, 60
        READ (IUNIN,'(A72)') ZEILE
        call fix_logical_input(zeile,60)
        READ (ZEILE,6665) (PLTSRC(J),J=I,MIN(I+59,NSTRA))
      END DO
C
      IF (LRPSCUT) THEN
        READ (IUNIN,6664) CUTPLANE(1:4)
!  TEST IF CUTPLANE IS AXIS-PARALLEL
        IF ((SUM(ABS(CUTPLANE(2:4)))-1._DP > EPS10) .OR.
     .      (CUTPLANE(2)*CUTPLANE(3) > 0._DP) .OR.
     .      (CUTPLANE(2)*CUTPLANE(4) > 0._DP) .OR.
     .      (CUTPLANE(3)*CUTPLANE(4) > 0._DP)) THEN
          WRITE (iunout,*) ' PLANE IS NOT PARALLEL TO ONE AXIS'
          WRITE (iunout,*) ' CUTTING OF TETRAHEDRA ABANDONED'
          WRITE (iunout,*) ' LRPSCUT SET TO .FALSE.'
          LRPSCUT = .FALSE.
        END IF
      END IF
C
      IPLANE=1
      LRAPS3D=.FALSE.
      LR3DCON=.FALSE.
      RAPSDEL=-HUGE(1._DP)
C
      DO 1150 J=1,NVOLPL
 1152   READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 1152
        READ (ZEILE,6666) NSP
        IF (NSP.GT.NPLT)
     .    CALL EIRENE_MASPRM('NPLT',4,NPLT,'NSP',3,NSP,IERROR)
        NSPTAL(J)=NSP
        READ (IUNIN,'(A72)') ZEILE
        call fix_logical_input(zeile,4)
        READ (ZEILE,6665) PLTL2D(J),PLTL3D(J),PLTLLG(J),PLTLER(J)
        READ (IUNIN,6664) TALZMI(J),TALZMA(J),TALXMI(J),TALXMA(J),
     .                    TALYMI(J),TALYMA(J)
        IF (PLTL2D(J)) THEN
          READ (IUNIN,'(A72)') ZEILE
          call fix_logical_input(zeile,2)
          READ (ZEILE,6665) LHIST2(J),LSMOT2(J)
          DO 1160 I=1,NSPTAL(J)
            READ (IUNIN,6666) ISPTAL(J,I),NTL,
     .                        NPLIN2(J,I),NPLOT2(J,I),NPLDL2(J,I)
            IF (NTL.LT.-NTALI.OR.NTL.GT.NTALV.OR.NTL.EQ.0) GOTO 990
            NPTALI(J,I)=NTL
            NPLIN2(J,I) = MAX(NPLIN2(J,I),1)
            NPLDL2(J,I) = MAX(NPLDL2(J,I),1)
 1160     CONTINUE
        ENDIF
        IF (PLTL3D(J)) THEN
          READ (IUNIN,'(A72)') ZEILE
          call fix_logical_input(zeile,8)
          READ (ZEILE,6665) LHIST3(J),LCNTR3(J),LSMOT3(J),
     .                      LRAPS3(J),LVECT3(J),LRPVC3(J),
     .                      LRPS3D,LRPSCN
          READ (IUNIN,'(A72)') ZEILE
          call fix_logical_input(zeile,3)
          READ (ZEILE,6665) LPRAD3(J),LPPOL3(J),LPTOR3(J)
          DO 1161 I=1,NSPTAL(J)
            READ (IUNIN,6666) ISPTAL(J,I),NTL,IPROJ3(J,I),
     .                        NPLI13(J,I),NPLO13(J,I),
     .                        NPLI23(J,I),NPLO23(J,I),IPLN
            IF (NTL.LT.-NTALI.OR.NTL.GT.NTALV.OR.NTL.EQ.0) GOTO 990
            NPTALI(J,I)=NTL
 1161     CONTINUE
          READ (IUNIN,6664) TALW1(J),TALW2(J),FCABS1(J),FCABS2(J),
     .                      RPSDL
          IF (FCABS1(J).LE.0.0D0) FCABS1(J)=1.0D0
          IF (FCABS2(J).LE.0.0D0) FCABS2(J)=1.0D0
          LRAPS3D=LRAPS3D.OR.LRPS3D
          LR3DCON=LR3DCON.OR.LRPSCN
          IPLANE=MAX(IPLANE,IPLN)
          RAPSDEL=MAX(RAPSDEL,RPSDL)
        ENDIF
 1150 CONTINUE
      IF (NLTRA) RAPSDEL=RAPSDEL*DEGRAD
C
C  SKIP INPUT LINES, UNTIL INPUT BLOCK 12 STARTS
 1175 READ (IUNIN,'(A72)') ZEILE
      IREAD=1
      IF (ZEILE(1:3) .EQ. '***') GOTO 1205
      GOTO 1175
C
C  READ DATA FOR DIAGNOSTIC MODULE  1200--1299
C
 1205 CALL EIRENE_MASAGE('*** 12. DATA FOR DIAGNOSTIC MODULE')
 1210 READ (IUNIN,'(A72)') ZEILE
      IREAD=1
      IF (ZEILE(1:1) .EQ. '*') GOTO 1210

c June 18: new: more general option for definition of emission lines (in NCHTAL=2 option)

cdr  read further atomic/molecular data: population coefficients, QSS ratios, etc
cdr       needed for setting up volumetric line emissivity profiles
cdr       as further add. tally ADDV (additional to those already defined in block 10a)
cdr       The ADDV tallies may then be used for line of sight integration, along
cdr       the CHORDS defined further below.
cdr  Distinct from the other addv tallies from block 10a these
cdr  further addv tallies (beyond NADVI) are currently apparently neither
cdr  scaled, averaged, integrated, nor do they have text (units, species) assigned.

cdr  check if such emissivity tallies are defined: search for 'DEFINE LINES' in next card
      ULINE=ZEILE
      CALL EIRENE_UPPERCASE(ULINE)
cdr  careful: slightly different from find_param
      NLEMIS = NCHOR > 0  !dr  and: NCHTAL=2 ??
      LDEF_LINES = .FALSE.
      
      IF (INDEX(ULINE,'DEFINE_LINES') > 0) THEN
cdr read volumetric emission profile data
        LDEF_LINES = .TRUE.       
        IADV = NADVI
        NLEMIS = .TRUE.
        READ (IUNIN,6666) NUM_LINES, MOD_ADDV
        IF (NUM_LINES > 0) THEN
          ALLOCATE (EMIS_LINES(NUM_LINES))
          EMIS_LINES%LINE_NAME = REPEAT(' ',80)
          EMIS_LINES%NUM_COMPO = 0
          WRITE (IUNOUT,*) 'FURTHER EMISSIVITY PROFILES SET:'
        END IF
        DO ILINE=1, NUM_LINES
          READ (IUNIN,'(A80)') ZEILE
          DO WHILE (ZEILE(1:1) == '*')
            READ (IUNIN,'(A80)') ZEILE
            write (IUNOUT,'(A80)') ZEILE
          END DO
cdr  further below, a particular emission profile is identified
cdr  (e.g. for a chord ichori) either by its name  (and ch_line%...)
cdr  or, if that fails, by its energy (and EMIN1 flag)
          READ (ZEILE,'(A80)') EMIS_LINES(ILINE)%LINE_NAME
          IREAD=0
          READ (IUNIN,6666) NUM_COMPO
          READ (IUNIN,6664) EMIS_LINES(ILINE)%EINSTEIN,
     .                      EMIS_LINES(ILINE)%TRANS_EN
C
          EMIS_LINES(ILINE)%NUM_COMPO = NUM_COMPO

c  Deal with ADDV storage for sum over components.
c  The storage on ADDV for individual components is done below (JCOMP loop).
cdr       IADV = IADV + 1
c  if mod_addv=0: reset storage needs for each line
c                 back to nadvi+1..nadvi+num_compo+1,
c                 i.e. addv tallies are only saved for one line at a time.
cdr       IF (MOD_ADDV == 0) IADV = NADVI + 1
          IF (MOD_ADDV == 0) IADV = NADVI

          EMIS_LINES(ILINE)%IADV_TOTAL = IADV + NUM_COMPO+1

          IF (NUM_COMPO > 0) THEN
            ALLOCATE (EMIS_LINES(ILINE)%COMPO(NUM_COMPO))

            DO JCOMP=1,NUM_COMPO

              READ (IUNIN,'(A72)')
     .              EMIS_LINES(ILINE)%COMPO(JCOMP)%COMPO_NAME
              READ (IUNIN,6666) NUM_CONTRIB, IRC
              EMIS_LINES(ILINE)%COMPO(JCOMP)%NUM_CONTRIB = NUM_CONTRIB
              EMIS_LINES(ILINE)%COMPO(JCOMP)%IRC = IRC
              ALLOCATE
     .         (EMIS_LINES(ILINE)%COMPO(JCOMP)%CONTRIB(NUM_CONTRIB))
cdr  for each new component: store emission profile on addv(iadv)
              IADV = IADV + 1
              EMIS_LINES(ILINE)%COMPO(JCOMP)%IADV = IADV

cdr  now we dwell on the contributions:
cdr  e.g. different isotopes,... but same rates, same population factors in each component
              DO KCONTR = 1, NUM_CONTRIB
                CNT%ISP = -1
                CNT%ITP = -1
                CNT%IRATIO = 0
                CNT%IRC_RAT = 0
                CNT%ISP_RAT = -1
                CNT%ITP_RAT = -1
                READ (IUNIN,6666) CNT%ISP, CNT%ITP, CNT%IRATIO

cdr  read QSS ratio between two densities, e.g.:  H2+/H2, if nfoli(H2+)=-1.
cdr  If density n_B of parent state (e.g. H2+) for upper level is not amongst the densities
cdr  known in this run,
cdr  but is in an (QSS) equilibrium (nfol<0) with such a density n_A (e.g. H2) instead.
                IF (CNT%IRATIO > 0) THEN
cdr  read QSS density ratio n_B/n_A(Te,ne). Then the
cdr  upper state population is n_B *pop_B(upper)= n_A * ratio * pop_B(upper)
                  READ (IUNIN,6666) CNT%IRC_RAT(1)
cdr   Read a second density ratio.
cdr   It may turn out that, after reading the first density ratio,
cdr   that now the new parent density n_A is still not amongst the density
cdr   known to eirene in this run. Then read a second density ratio.
cdr   Example:  n_B=n_H3+, ratio1=nH3+/nH2+. I.e. n_A =nH2+.
cdr   See routine emissivity.f for further explanations.
                  IF (CNT%IRATIO == 2) THEN
                    READ (IUNIN,6666) CNT%IRC_RAT(2),
     .                                CNT%ISP_RAT(1),CNT%ITP_RAT(1),
     .                                CNT%ISP_RAT(2),CNT%ITP_RAT(2)
                    END IF
                  END IF

                EMIS_LINES(ILINE)%COMPO(JCOMP)%CONTRIB(KCONTR) = CNT
              END DO  !  kcontr

            END DO    ! jcomp
c  increase counter iadv, for next line, because sum over comp.
c                         is stored on num_compo+1
            IADV=IADV+1
          END IF      ! if jcomp.gt.0
        END DO        ! iline
        READ (IUNIN,'(A80)') ZEILE
        DO WHILE (ZEILE(1:1) == '*')
          READ (IUNIN,'(A80)') ZEILE
        END DO
        IREAD=1
      ELSEIF (INDEX(ULINE,'DEFAULT_LINES') > 0) THEN
        IREAD=0
        NLEMIS=.TRUE.
        READ (IUNIN,'(A80)') ZEILE
        DO WHILE (ZEILE(1:1) == '*')
          READ (IUNIN,'(A80)') ZEILE
        END DO
        IREAD=1
      ENDIF

! no definition of emissivity lines was read in
! define OLD default emissivity model for chords, for backward compatibility.
cdr  allocate storage and fill structure EMIS_LINES
cdr  such that old default options are recovered.
cdr This is exclusive: as soon as at least one line emission profile is
cdr read from block "12.0", no default emissivities are set at all.
      IF (NLEMIS.AND..NOT.ALLOCATED(EMIS_LINES))
     .   CALL EIRENE_SETUP_DEFAULT_EMISSIVITY

c  June 18: new:  end of new code, further modifications below, for NCHTAL=2 option
      CALL EIRENE_MASAGE
     .  ('*** 12A. DATA FOR LINES OF SIGHT              ')
      IF (IREAD == 0) READ (IUNIN,'(A80)') ZEILE
      READ (ZEILE,6666) NCHORI,NCHENI
      WRITE (iunout,*) '       NCHORI= ',NCHORI
      IREAD=0
      NCHOR = NCHORI
      NCHEN = NCHENI
      CALL EIRENE_LEER(1)
      WRITE (iunout,'(1x,a,2i5)') 'NCHORI,NCHENI= ',NCHORI,NCHENI
      CALL EIRENE_LEER(1)
      IF (IABS(NCHENI).GT.NCHEN)
     .    CALL EIRENE_MASPRM
     .         ('NCHEN',5,NCHEN,'IABS(NCHENI)',12,IABS(NCHENI),
     .                 IERROR)
      IF (NCHORI.LE.0) GOTO 1230
      CALL EIRENE_ALLOC_COMSIG
      CALL EIRENE_INIT_COMSIG
      DO 1220 ICHORI=1,NCHORI
        READ (IUNIN,'(A72)') TXTSIG(ICHORI)
        READ (IUNIN,6666) NCHTAL(ICHORI),NSPSCL(ICHORI),NSPNEW(ICHORI),
     .                    ISTCHR
cdr  new input option: generalized side on line emissivities
cdr                    for nchtal=2 option.
        READ (IUNIN,'(A400)') ZEILE
        IREAD=1

        IF (NCHTAL(ICHORI) == 2) THEN
cdr  old default: neither "define_lines" nor "default_lines" found.
cdr  still:  chords with NCHTAL=2 present. Set default now, latest. 
          IF (.NOT.ALLOCATED(EMIS_LINES))
     .      CALL EIRENE_SETUP_DEFAULT_EMISSIVITY
cdr  for nchtal=2: one optional extra input card may be read:  search for 'USE_LINE'
cdr  and fill CH_LINE_NAME(ICHORI) with the name of that line.
cdr  Alternatively the energy parameters EMIN1 may be used.
          ULINE = ZEILE
          CALL EIRENE_UPPERCASE(ULINE)
          IND = INDEX(ULINE,'USE_LINE')
          IF (IND > 0) THEN
            READ (ZEILE(IND+8:),'(A80)') CH_LINE_NAME(ICHORI)
            READ (IUNIN,'(A400)') ZEILE
           IREAD=1
          END IF
        END IF

        READ (ZEILE,6666) NSPSTR(ICHORI),NSPSPZ(ICHORI),  ! here should come: NSPTP(..), TYPE
     .                    NSPINI(ICHORI),NSPEND(ICHORI),
     .                    NSPBLC(ICHORI),NSPADD(ICHORI)
        IREAD=0
        READ (IUNIN,6664) EMIN1(ICHORI),EMAX1(ICHORI),ESHIFT(ICHORI)
        READ (IUNIN,66664) IPIVOT(ICHORI),
     .                     XPIVOT(ICHORI),YPIVOT(ICHORI),ZPIVOT(ICHORI)
        READ (IUNIN,66664) ICHORD(ICHORI),
     .                     XCHORD(ICHORI),YCHORD(ICHORI),ZCHORD(ICHORI)
        NLSTCHR(ICHORI) = ISTCHR > 0  ! automatically add directional cell-based spectra, along line of sight
 1220 CONTINUE   ! ichori

      READ (IUNIN,'(A72)') ZEILE
      call fix_logical_input(zeile,5)
      READ (ZEILE,6665) PLCHOR,PLSPEC,PRSPEC,PLARGL,PRARGL

 1230 CONTINUE  ! nchori > 0

C  SKIP READING REST OF THIS BLOCK
      READ (IUNIN,'(A72)') ZEILE
      IREAD=0
      IF (ZEILE(1:3).NE.'***') GOTO 1230
      IREAD=1
C
C  READ DATA FOR TIME-DEPENDENT AND NONLINEAR MODE  1300--1399
C
      IF (IREAD.EQ.0) READ (IUNIN,*)
      IREAD=0
      CALL EIRENE_MASAGE
     .  ('*** 13. DATA FOR ITERATIVE AND TIME DEP. OPTION')
C
      NINITL_READ = 0
      READ (IUNIN,6666) NPRNLI_IN, NINITL_READ, NPRMUL
      NPRNLI = NPRNLI_IN
      IF (NPRMUL > 1) NPRNLI=NPRNLI * NPRMUL

      CALL EIRENE_LEER(1)

      IF (NLERG.AND.NPRNLI.LE.0) THEN
C  NO TIME HORIZON DEFINED, DESPITE NLERG=.TRUE.
C  THEREFORE: SET A DEFAULT TIME HORIZON HERE
        NPRNLI=100
        IF (NTIME.EQ.0) NTIME=1
        WRITE (iunout,*) '        NPRNLI= ',NPRNLI,
     .                   ' (MODIFIED DUE TO NLERG)'
      ELSE
        WRITE (iunout,*) '        NPRNLI= ',NPRNLI
      ENDIF
      CALL EIRENE_LEER(1)

      IF (NPRNLI.LE.0.OR.NTIME.EQ.0) THEN
C  TURN OFF TIME DEP MODE IF EITHER NTIME=0 OR NPRNLI=0
        IF (NPRNLI.GT.0) THEN
          WRITE (IUNOUT,*) 'TIME DEP. MODE TURNED OFF, BECAUSE NTIME=0'
          NPRNLI=0
        ENDIF
        IF (NTIME.GT.0) THEN
          WRITE (IUNOUT,*) 'TIME DEP. MODE TURNED OFF, BECAUSE NPRNLI=0'
          NTIME=0
        ENDIF
        GOTO 1350
      ENDIF

      LDEF_TIME_HORIZON = .FALSE.
      READ (IUNIN,'(A72)') ZEILE
      IREAD=1
      IF (ZEILE(1:1).EQ.'*') THEN
C  DATA FOR DEFAULT TIME HORIZON
        NTIME0=0
        NTIME=1
        DTIMV=1.D0
        TIME0=0.D0
        NPTST=0
        NTMSTP=1
        NSNVI=0
        IREAD=0
        LDEF_TIME_HORIZON = .TRUE.
        GOTO 1350
      ENDIF
cdr  skip reading further comment lines...tbd
c
      READ (ZEILE,6666) NPTST,NTMSTP
      IREAD=0

C   ENFORCE ONE-BY-ONE RELAUNCH FROM CENSUS, IN CASE NLMOVIE
      IF (NLMOVIE.AND.NPTST.GE.0) THEN
        WRITE (IUNOUT,*) 'NPTST RESET TO -1, BECAUSE OF NLMOVIE OPTION'
        NPTST=-1
      ENDIF

      READ (IUNIN,6664) DTIMV,TIME0
      IREAD=0
c
C   READ DATA FOR SNAPSHOT TALLIES
      READ (IUNIN,*)
      CALL EIRENE_MASAGE('*** 13A. DATA FOR SNAPSHOT TALLIES')
      READ (IUNIN,6666) NSNVI
      IREAD=0
      WRITE (iunout,*) '        NSNVI= ',NSNVI
      CALL EIRENE_LEER(1)
      DO 1320 J=1,NSNVI
 1321   READ (IUNIN,'(A72)') ZEILE
        IREAD=1
        IF (ZEILE(1:1) .EQ. '*') GOTO 1321
        READ (ZEILE,6666) ISNVE(J),ISNVS(J),ISNVT(J),ISNRC(J)
        IREAD=0
!       READ (IUNIN,'(A72)') TXTTAL(J,NTALT)
!       READ (IUNIN,'(2A24)') TXTSPC(J,NTALT),TXTUNT(J,NTALT)
        READ (IUNIN,'(A72)') TXTTLT(J)
        READ (IUNIN,'(2A24)') TXTSCT(J),TXTUTT(J)
 1320 CONTINUE
C
      IF (NTIME.LE.0) THEN
        WRITE (iunout,*) 'ERROR IN INPUT: TIME DEP. MODE BUT NTIME.LE.0'
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
 1350 CONTINUE
C  SKIP READING REST OF THIS BLOCK
      IF (IREAD.EQ.0) READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:3).NE.'***') GOTO 1350
      IREAD=1

      IF (IREAD.EQ.0) READ (IUNIN,*)
C
C
 6662 FORMAT (L1,1X,A24,1X,I1,1X,4(2I3,1X))
 6664 FORMAT (6E12.4)
 6665 FORMAT (12(5L1,1X))
 6666 FORMAT (12I6)
66661 FORMAT (I3,1X,A6,1X,A4,A9,A3,2I3,3E12.4)
66662 FORMAT (2E12.4,10I6)
66664 FORMAT (I6,6X,5E12.4)
66665 FORMAT (2I6,4E12.4)
66666 FORMAT (I2,1X,A8,12(I3),1X,A10,1X,I2)

      RETURN
C
C  ERROR EXITS
C
  990 CONTINUE
      WRITE (iunout,*) 'TALLY NUMBER FOR PRINTOUT OR PLOT OF'
      WRITE (iunout,*) 'VOLUME-AVERAGED TALLIES OUT OF RANGE'
      CALL EIRENE_EXIT_OWN(1)
  991 CONTINUE
      WRITE (iunout,*) 'TALLY NUMBER FOR PRINTOUT OF SURFACE-AVERAGED'
      WRITE (iunout,*) 'TALLIES OUT OF RANGE'
      CALL EIRENE_EXIT_OWN(1)
  993 CONTINUE
      WRITE (iunout,*)
     .  'INCONSISTENCY WRT. SURFACE REFLECTION DATABASE, BLOCK 6 '
      CALL EIRENE_EXIT_OWN(1)
  994 CONTINUE
      WRITE (iunout,*) 'ERROR IN INPUT: NRPLG.NE.NP2ND, BUT NLPOL=TRUE'
      WRITE (iunout,*) 'NRPLG,NP2ND ',NRPLG,NP2ND
      CALL EIRENE_EXIT_OWN(1)
  995 CONTINUE
      WRITE (iunout,*) 'ERROR IN INPUT: ',
     .            'WRONG UNIT NUMBER FOR OUTPUT OF TALLY SPECIFIED'
      WRITE (iunout,*) 'NTLV,NTLVF',NTLV,NTLVF
      CALL EIRENE_EXIT_OWN(1)
  998 CONTINUE
      WRITE (iunout,*) 'ERROR IN INPUT BLOCK FOR ADDITIONAL DATA FOR'
      WRITE (iunout,*) 'SPECIFIC ZONES FOUND AT ZONE NO. ',I
      CALL EIRENE_EXIT_OWN(1)
      
      END SUBROUTINE EIRENE_READ_FIXFORM
      

      SUBROUTINE EIRENE_READ_BLK14_FIXFORM(IERROR)
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT
      USE EIRMOD_CTEXT
      USE EIRMOD_CSPEI
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: IERROR 
      CHARACTER(420) :: ZEILE, ULINE
      INTEGER :: IDUM, NCOPIE, J
C
C
C  ADDITIONAL INPUT FOR THIS RUN COMES FROM EITHER
C  ANOTHER CODE (DATA FILE) OR FROM AN EARLIER RUN OF EIRENE
C
C  INPUT BLOCK 14 BEGIN
C
C  READ DATA IN INTERFACING SUBROUTINE INFCOP  1400 -- 1499
C
      CALL EIRENE_LEER(1)
      CALL EIRENE_MASAGE
     .  ('*** 14. DATA FOR INTERFACING ROUTINE "INFCOP"')
      IF (NMODE.EQ.0) THEN
C  STAND ALONE RUN, READ BLOCK *** 14 HERE
        WRITE (iunout,*) '        SUBR. INFCOP NOT CALLED.'
!pb NCOPII is completely redundant        
!       READ (IUNIN,6666) NAINI,NCOPII,NCOPIE
        READ (IUNIN,6666) NAINI,IDUM,NCOPIE
        NCPVI=NCOPIE
        WRITE (iunout,*) '        NAINI, NCPVI = ',NAINI,NCPVI
        IF (NAINI.GT.NAIN) THEN
          CALL EIRENE_MASPRM('NAIN',4,NAIN,'NAINI',5,NAINI,IERROR)
!PB          GOTO 1500
          WRITE (iunout,*) IERROR,' INPUT OR PARAMETER ERRORS DETECTED'
          WRITE (iunout,*)
     .      ' SEE THE ERROR MESSAGES LISTED ABOVE AND CORRECT'
          WRITE (iunout,*) ' THE ERRORS BEFORE RE-EXECUTION'
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
        
        IF (NCPVI.GT.NCPV) THEN
          CALL EIRENE_MASPRM('NCPV',4,NCPV,'NCPVI',5,NCPVI,IERROR)
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
!pb     CALL EIRENE_ALLOC_CCOUPL(2)
        DO 3020 J=1,NAINI
 3021     READ (IUNIN,'(A72)') ZEILE
          IF (ZEILE(1:1) .EQ. '*') GOTO 3021
          READ (ZEILE,6666) NAINS(J),NAINT(J)
          READ (IUNIN,'(A72)') TXTPLS(J,NTALN)
          READ (IUNIN,'(2A24)') TXTPSP(J,NTALN),TXTPUN(J,NTALN)
 3020   CONTINUE

      ELSEIF (NMODE.NE.0) THEN
C  COUPLED RUN, READ BLOCK *** 14 IN INTERFACING ROUTINE INFCOP (ENTRY IF0COP)
        NCPVI=0
        NAINI=0
C  READ BLOCK 14 AND GEOMETRY FROM EXTERNAL DATABASE (FT30), also set NAINI, NCOPII, NCOPIE there
        CALL EIRENE_IF0COP(.TRUE.)
      ENDIF

      RETURN
 6666 FORMAT (12I6)
      END SUBROUTINE EIRENE_READ_BLK14_FIXFORM
