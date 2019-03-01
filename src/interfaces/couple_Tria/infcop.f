!PB   17.11.05  NEMODS=-2  =>  NEMODS=8
!PB             NEMODS=-3  =>  NEMODS=9
!PB   14.07.11  OVERWRITE NEMODS SPECIFIED BY THE EIRENE INPUT BY
!PB             SETTINGS CONSISTENT WITH SHEATH SETTINGS IN B2
!PB             IF SHEATH PARAMETERS ARE GIVEN BY B2 NEMODS=9
!PB             IF NO SHEATH PARAMETERS WERE TRANSFERRED BY B2 NEMODS=8
CDR  09. 2014 GENERAL RELATIONS BETWEEN FINE (UNSTRUCTURED)  AND COARSE (STRUCTURED)
C             GRID MADE MORE EXPLICIT.  NEW NSBOX FOR FINE GRID SET.
C             NP2NDQ REMOVED (REDUNDANT: =NP2TAL)
C             NR1STQ REMOVED (REDUNDANT: =NR1TAL),
C             AND ERROR: NR1STQ WAS USED BEFORE DEFINITION --> PROBLEMS WITH NSTGRD ARRAYS
C             VIA FILES FROM FORT.29?

c             plus minor notational cleanup, comments added

cdr  23.04.2015
cdr  issue: underlying coarse grid format for printout?
cdr  see previous version.
cdr  now:
cdr  ncelln test for coarse grid cell number replaced, by ixtri(itri)=0 test
cdr
cdr  distinct from solps-iter version:
cdr  nr1tal, np2tal still set such that underlying structured grid is enabled for printout
cdr  with PRTTAL, and SYMET option.  To be checked: is this consistent with NCLTAL coarse-graining
cdr  to be done:
cdr  can one always find a "ncltal map" such that underlying structured grid is preserved
cdr  for printing and plotting?
c

cpb  15.09.15:  added: default bfield =1 (tesla), if bfield=0, cell-wise.


cdr  start to use species-resolved energy tallies.
cdr  nov.15:  eapl,empl,eipl: now species-resolved
cdr           --> eppl_cop, eploda get an additional species index ipls
C
cdr  dec. 15: not ready.  started to comment, and to extend copv tallies
cdr           for total and internal energy with species index.  Not ready...

cdr  tbd:  check storage on copv, whenever used.

cdr  jan.17
cdr  epel(nrad)      --> epel_cop(nrad)
cdr  cppv(ncpv,nrad) --> mppl_cop(npls,nrad), for notational consistency
c    also done        :  cpvoda(ncpv,...) --> cpvoda(npls,...). perhaps 'ncpv' can go out??
C    to be done          CPPVS            --> mppls_cop
c    new: write file of triang. data for idl plotting tool.

cdr  a number of changes involving internal energy source options have been changed.
cdr  to be checked.
cdr  distinct from solps_iter version:  all pppl_cop, mppl_cop, eppl_cop, epel_cop
cdr  are all set.
cdr  in solps_iter: only pppl_cop is set. eppl_cop and epel_cop contributions are directly
cdr                 removed from eapl, eael, and not set. --> short cycle not possibel there.
cdr                 mppl_cop is set, but not on mapl.  what does that mean in wrneutrals ??

cdr   cheim needs species index?
cdr   a number of 'save' missing, compared to solps_iter
cdr   llcut (FZJ versions) rather than lcut (solps_iter), lcut is broadcasted and on common
cdr   eirmod_cpolyg. Also set again in timer. same meaning?
c
cdr      additional input read from block 14:
c        LCHKQUD   uncommented read variable, taken out from read ***14,
c                  now set to .false. Fills array IREVERS. Still needed ?
cdr
c             IMF,ITCO
c  unused ?

cdr feb 17:  nradd_tal corrected: NRADD was added twice
c            set vacuum parameters also in additional cells from block 2e,
c            not just on those from tria-grid generation.
c  made allocatable (as in solps_iter): chpm, cheem, cheim,chmom
c  default setting of nemod1=9 and 8 at targets removed (if(false)....)
c  bug fix re lgjum3 implemented from couple_solps

cdr Feb 17:  remove duplicated code re call to geousr.
c
c            double grid structure: triangles and polygons, for tallies.
c            parameters nr1st,  np2nd,  ..... for geometry
c            parameters nr1tal, np2tal, ..... for scoring
c            parameters nr1tal_save, ....     for interfacing tallies between b2 and eirene
c                                             is always the b2 (structured) coarse grid

cdr June 17: NFL dependence added:  SFNIT(0:NSTEP,NFL) in global particle balance: from couple_SOLPS
c                            sfni..(..,nfl), ssns,ssni,balann,totn,rn.
cdr Nov 17 : precision checks with triangles vs. quadrangles: use relative units.
cdr          only partially done

cdr Jan 18 : bug fix re vol.rec., only one ipls per stratum is supported
c            code was correct in solps4.3, and garching versions of couple_b2/b2.5

cdr March 18: new variable LCOARSE: maintain underlying coarse structured grid, scoring
cdr           on coarse grid (NCLTAL array). Otherwise: only fine (triangular) grid structure
cdr           remove unused array: scpveii
cdr Oct.  18: bug fix: bvin(iplsv) rather than bvin(ipls) in one place
c......................................................................................


C   EIRENE CODE SEGMENT COUPLE_$, $ MAY CURRENTLY STAND FOR B2,
C                                                           B2.5,
C                                                           DIVIMP,
C                                                           TRIA,
C                                                           TETRA,
C                                                           TRANSP,
C                                                           DUMMY
C
C   THIS VERSION: $=Tria,  NOV. 2002
c                 proprietary version of FZJ, for local B2 code versions
C
c  geometry data not any longer via work array into eirene
c                due to module structure
c  eliminate cut cells from balances ( plus: rename lcut to llcut(..)... why?)
c  new input: ncopib, ncopeb
C             fniprt

C   UPDATES:
C   OPTION TO EVALUATE B-FIELD VECTORS FROM GRIDADAP FILE FT29
C   FOR NON-ORTHOGONAL GRIDS
C
C   THIS CODE SEGMENT CONTAINS VARIOUS SUBROUTINES NEEDED FOR
C   INTERFACING THE EIRENE CODE TO PLASMA FLUID CODES.
C   IT READS GEOMETRICAL DATA (MESHES) FROM FILE FT30
C   AND PRODUCES THE EIRENE INPUT DATA (BLOCK 2).
C   IT READS PLASMA BACKGROUND DATA FROM FILE (FT31) OR COMMON BLOCKS,
C   IT MAY (OPTIONAL) ALSO READ PLASMA DATA FROM FILE FT13
C   WRITTEN IN A PREVIOUS EIRENE RUN (E.G. IN ORDER TO ITERATE
C   IN SOME BACKGROUND SPECIES)
C   IT THEN PRODUCES INPUT DATA FOR EIRENE
C   INPUT BLOCK 5 (PLASMA DATA) AND BLOCK 7 (SURFACE RECYCLING SOURCES)
C
C
C   THIS PARTICULAR VERSION LINKS EIRENE TO THE B2_Tria   2D MULTIFLUID EDGE
C   PLASMA TRANSPORT CODE.
C
C   IT WAS WRITTEN BY D.REITER AND P.BOERNER, FZ-JUELICH
C   E-MAIL: D.REITER @ FZ-JUELICH.DE, AND: www.eirene.de
C
C
C
C   MOST OF THE FORTRAN IN THIS CODE SEGMENT HAS BEEN DEVELOPED
C   UNDER KFA-NET CONTRACT NO. 428/90-8/FU-D
C
C   FINAL REPORT BY: D.REITER(1), P.BOERNER(1), B.KUEPPERS(1),
C                    M.BAELMANS(2) AND G.P.MADDISON(3)
C                    (1992)
C   1): KFA-JUELICH GMBH
C   2): UNIV. LEUVEN, ERM, KFA-JUELICH
C   3): AEA TECHNOLOGY, FUSION, CULHAM, UKAEA FUSION ASSOCIATION
C
*DK COUPLE
C
      SUBROUTINE EIRENE_INFCOP
C
C     THIS SUBROUTINE DEFINES THE PLASMA MODEL IN CASE OF A COUPLED
C     NEUTRAL-PLASMA CALCULATION
C
C     THE ENTRY "IF0COP" RECEIVES GEOMETRICAL INPUT DATA FROM AN
C     EXTERNAL FILE (E.G. OTHER PLASMA CODES)
C     AND PREPARES THEM FOR AN EIRENE RUN
C
C     THE ENTRY "IF1COP" RECEIVES PLASMA INPUT DATA FROM AN
C     EXTERNAL FILE (E.G. OTHER PLASMA CODES)
C     AND PREPARES THEM FOR AN EIRENE RUN
C
C     THE ENTRY "IF2COP" PREPARES THE SOURCE SAMPLING DISTRIBUTION
C     FROM THE EXTERNAL DATA, AND MAY OVERWRITE OTHER INPUT
C     DATA FROM BLOCKS 1 TO 13 AS WELL
C
C     THE ENTRIES "IF3COP, IF4COP" RETURN  RESULTS TO AN EXTERNAL CODE
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_BRASPOI
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_CADGEO
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CPLOT
      USE EIRMOD_CINIT
      USE EIRMOD_CPOLYG
      USE EIRMOD_CGRID
      USE EIRMOD_CZT1
      USE EIRMOD_CTRCEI
      USE EIRMOD_CCOUPL
      USE EIRMOD_CGEOM
      USE EIRMOD_CSDVI
      USE EIRMOD_CSDVI_BGK
      USE EIRMOD_CSDVI_COP
      USE EIRMOD_COMPRT
      USE EIRMOD_COMNNL
      USE EIRMOD_COMSOU
      USE EIRMOD_CSTEP
      USE EIRMOD_CTEXT
      USE EIRMOD_CLGIN
      USE EIRMOD_COUTAU
      USE EIRMOD_COMXS
      USE EIRMOD_CSPEI
      USE EIRMOD_CTRIG
C  PLASMA DATA: NI,TE,TI,VV,UU,PR,UP,RR,FNIX,FNIY.. (BRAAMS ---> EIRENE)
      USE EIRMOD_BRAEIR
C  NEUTRAL SOURCE TERMS: SNI,SMO,SEE,SEI (EIRENE ---> BRAAMS)
      USE EIRMOD_EIRBRA
      USE EIRMOD_BRASCL

      IMPLICIT NONE
C
C
C  GEOMETRICAL DATA FROM GRIDADAP
      REAL(DP), ALLOCATABLE, SAVE ::
     R  ALPHXB(:,:), ALPHYB(:,:), XAISO(:,:)

      REAL(DP), ALLOCATABLE, SAVE ::
     R  PUX(:),      PUY(:),      PVX(:),      PVY(:),
     R  PUXE(:), PUYE(:), PUXN(:), PUYN(:),
c  the 4 pv... arrays. Needed?
     R  PVXE(:), PVYE(:), PVXN(:), PVYN(:)

      INTEGER, ALLOCATABLE, SAVE ::
     I  IAISO(:,:)
C
      TYPE(CELLSIM), POINTER :: CPSIM
      TYPE(CELLMUL), POINTER :: CPMUL
C
      REAL(DP) :: SEES0(NSTRA), SEIS0(NSTRA)
      REAL(DP), ALLOCATABLE, SAVE ::
     .            CHPM(:,:), CHEEM(:), CHEIM(:),
     .            CHMOM(:,:)
      REAL(DP) :: DI(NPLS), VP(NPLS)

cdr for species-dependent global particle balance
      REAL(DP) :: SFNISY(NFL),SFNINY(NFL),SFNIWX(NFL),SFNIEX(NFL)
      REAL(DP) :: SSN(NFL),SSNI(NFL),BALANN(NFL),TOTN(NFL),RN(NFL)

C pppl_cop, mppl_cop, eppl_cop and epel_cop are the exact
c volumetric source tallies,
c while default tallies pppl, mppl, eppl and epel are
c the corresponding tallies scored from random sampling in eirene
      REAL(DP), ALLOCATABLE, SAVE ::
     .            PPPL_COP(:,:), MPPL_COP(:,:),
     .            EPPL_COP(:,:), EPEL_COP(:)
C
c  for short cycle correction terms, in vol. rec. strata.
      REAL(DP), ALLOCATABLE, SAVE ::
     .            PPLODA(:,:), CPVODA(:,:),
     .            EPLODA(:,:), EPEODA(:)
C
      REAL(DP) :: EFLX(NSTRA),
     R          DUMMY(0:NDXP,0:NDYP),
     R          SFNIT(0:NSTEP,NFL), SFEIT(0:NSTEP),
     R          SFEET(0:NSTEP), SHEAE(0:NSTEP), SHEAI(0:NSTEP)

      INTEGER :: NRWL(NSTRA)

      REAL(DP), SAVE :: SCALM, SCALE, SCALI, SEES, SEIS,
     .          SFEISY, SFEESY,
     .          VPARA, RECADD, RECTOT,
     .          EEADD, EIADD, PIADD,
     .          SMOCL, CHEES, CHEIS, SNICL,
     .          SIGNUM, FNIYB0,
     .          SSE, BALANI, BALANE, SSEE, SSI, RE, RI, RNT, TOT,
     .          TOTI, TOTE, SFEENY, SFEIWX, BALAN,
     .          TIFLX, PIFLX,
     .          SSEI, SFEIEX, SFEEEX, SFEEWX, SFEINY,
     .          XCOOR, YCOOR, ZCOOR, VSX, VSY, VS, VTX, VL, BVAC, V, T,
     .          BN, BX, BY, BZ,
     .          DELTE_PARA, DELTE_PERP, DELTI_PERP, DELTI_PARA,
     .          TX, TY, VPRO, VTY, VT, TEST, XMUE, PX,
     .          PY,
     .          ALX, ALE, ALW, ALS, ALN, AL, UUBC, VTEST, VTEST2,
     .          VR, CS,
     .          PERW, PARW, PARWI, PERWI, DRR, EADD, EMAXW, ESHEATH,
     .          TE, CUR,
     .          VPZ, PM1, VPY, PN1, VPX, GAMMA, ESUM,
     .          CHI, CHP, CHE, SUMEI, SUMEE, SUMM, SUMN,
     .          ETOT, FLXI, FLX, OR, DELY, XANF, YANF, UDBC,
     .          RBC, UPBC, VVBC, DELX, PUYS, RRBS,
     .          EEMAX, EESHT, RP1, THMAX, TIS, TES,
     .          PIPV, PUXS, PNORM, PVXS, PVYS, PUPV,
     .          FLX_EIR,
     .          SUMN_OLD,
     .          SNIRES, SMORES, SEERES, SEIRES, UU, PITB,
     .          DXPOL, DYPOL, PAR,
     .          fniprt, fltt, e0b2, frac, celdel, dd, cfac

      INTEGER, SAVE :: J, IRC, JC, INC,
     .           K, IADD, NAS, IPUNKT, NSSIR,
     .           NUMSI, NBAR, ISNR, ISC, IS, NASMOD, NRS, NADMOD,
     .           NBARSI, IP1, IFL, IS1, IR1, IAIN, IAOT, IREAD,
     .           NTGPRI, IPRT, IO29, NEND, NINI, NSSIP,
     .           LTARG, I, IPL, IERROR, IMODE, NPLP,
     .           NRED, IDUMMY,  ISTS, ITRI, ISTR,
     .           NREC11, NEM, MINSPEZ, MAXSPEZ,
     .           IR, IP, IT, IA, IB, JUN,
     .           IN, IX, IY, NDX2, IIRC, NDXY,
     .           IFIRST, ISTRAI, NPES, MTRI, NPEC, NPBS,
     .           NPBC, IACT, IANF, IO, IIPLS, IEPLS, IG, ITARG, IGITT,
     .           I34, IRRC, ISC1, ISC2, ISCS, ICOU, ISP,
     .           IXI, IXE, NCOPIB, NCOPEB, IPLSTI, IPLSV, IPLV,
     .           IST_RATE, MSHFRM, IMF, ITCO,
!  ADDITIONAL STORAGE FOR LIN. COMB. OF TALLIES (E.G. INTERNAL ENERGY SOURCES)
     .           ICPV, icp1, icp2, icp3,
     .           icoadd, icoscr, icog,
     .           istat_cop, IXM1, IYM1,
     .           ntrfrm,
     .           nr1tal_save,np2tal_save,nt3tal_save,nsbox_tal_save,
     .           nsurf_tal_save,nradd_tal_save

      INTEGER, INTENT(IN) :: ISTRAA, ISTRAE, NEW_ITER, IFRST, ITRG
      REAL(DP) :: EIRENE_STEP, EIRENE_FTABRC1, EIRENE_FEELRC1,
     .            EIRENE_SHEATH, EIRENE_EMAXW
      INTEGER, EXTERNAL :: EIRENE_IDEZ
C
      LOGICAL, INTENT(INOUT) :: LSTP
      LOGICAL, SAVE :: LSHORT, LSTOP, LTEST, LSTP3,
     .                 LNONREC_SY,LNONREC_NY,LNONREC_WX,LNONREC_EX
     .                ,IFBOUND, lchkqud,LCOARSE

      LOGICAL, ALLOCATABLE, SAVE :: LLCUT(:)

      logical :: l1, l2, lxsrf
CTRIG A
      TYPE :: CELL
        INTEGER :: TRIANGLE
        TYPE(CELL),POINTER :: NEXT
      END TYPE CELL

      TYPE :: POIFELD
        TYPE (CELL),POINTER :: P
      END TYPE POIFELD

      TYPE (POIFELD), ALLOCATABLE, SAVE :: HEADS(:,:)
      TYPE (CELL),POINTER :: CURPOI
CTRIG E
C

      REAL(DP), ALLOCATABLE, SAVE ::
     . CHPS(:),    SNIS(:),    CHMOS(:),  SMOS(:),  SCALN(:),
     . SNIS0(:,:), SMOS0(:,:),
c
     . RESSNI(:,:),  RESSMO(:,:),
     . RESSEE(:), RESSEI(:)
     .,FLXEIR(:)

      REAL(DP), ALLOCATABLE, SAVE ::
     . TORL(:,:), ESHT(:,:), ELTEST(:,:), ORI(:,:)

      real(dp), allocatable :: helpw(:)

      REAL(DP), ALLOCATABLE, save :: uuba(:,:,:), upba(:,:,:)
      REAL(DP), ALLOCATABLE, save :: uubh(:,:,:), upbh(:,:,:)

      logical :: lhit(nrad)


      INTEGER, ALLOCATABLE, SAVE :: IHELP(:)
C
      CHARACTER(10) :: CHR
      CHARACTER(6)  :: CITARG
      CHARACTER(72) :: ZEILE
      CHARACTER(20) :: FORM
      CHARACTER(1) :: NSEW

      TYPE(RATE_STORE), POINTER :: RTIS
C
      DATA LTARG/0/
C
CTRIG A
C
      INTEGER, ALLOCATABLE :: NUMTRI(:), NUMSID(:)
CTRIG E
C
C
      ENTRY EIRENE_IF0COP
C
      LSHORT=.FALSE.
C
      GOTO 99990
C
C  TO INITIALIZE THE SHORT CYCLING, THE GEOMETRY HAS TO BE
C  DEFINED ONCE (ENTRY: INTER0)
C
      ENTRY EIRENE_INTER0
      LSHORT=.TRUE.
99990 CONTINUE

      call eirene_leer(2)
      write (iunout,*) 'Subr. INFCOP called '
      write (iunout,*) 'This is a proprietary FZJ version of an '
      write (iunout,*) 'interfacing code to B2, B2.5 plasma solvers.'
      write (iunout,*) 'NOT ready for 3rd parties'
      call eirene_leer(2)
C
      IERROR=0
C
      IMODE=IABS(NMODE)
C
      IF (.NOT.ALLOCATED(CHPM)) THEN
        ALLOCATE (CHPM(NPLS,NRAD))
        ALLOCATE (CHEEM(NRAD))
        ALLOCATE (CHEIM(NRAD))
        ALLOCATE (CHMOM(NPLS,NRAD))

        ALLOCATE (PPPL_COP(NPLS,NRAD))
        ALLOCATE (MPPL_COP(NPLS,NRAD))
        ALLOCATE (EPPL_COP(NPLS,NRAD))
        ALLOCATE (EPEL_COP(NRAD))

        ALLOCATE (PPLODA(NPLS,NRAD))
        ALLOCATE (CPVODA(NPLS,NRAD))
        ALLOCATE (EPLODA(NPLS,NRAD))
        ALLOCATE (EPEODA(NRAD))
      END IF
cdr March 18: removed from input block 14. Unclear meaning.
      lchkqud = .false.
cdr
      mshfrm = 0   !  optional flag for geometry file format: linda, carre, sonnet
      ntrfrm = 0
C
      IF (.NOT.LSHORT.AND.ITIMV.LE.1) THEN
        WRITE (iunout,*) '        SUBROUTINE INFCOP IS CALLED  '
C  READ INPUT DATA OF BLOCK 14
C  SAVE INPUT DATA OF BLOCK 14 FOR SHORT CYCLE ON COMMON CCOUPL
        CALL EIRENE_LEER(1)
        CALL EIRENE_ALLOC_CCOUPL(1)
        READ (IUNIN,'(5L1)') LSYMET,LBALAN,LCOARSE
        IF (TRCINT)
     .  WRITE (iunout,*) ' LSYMET,LBALAN,LCOARSE = ',
     .                     LSYMET,LBALAN,LCOARSE

        READ (IUNIN,'(5I6)') NFLA,NCUTB,NCUTL,IMF,
cdr  imf  flag for different formats of geometry file: linda, sonnet, carree. What is What?
     .                       ntrfrm
        IF (IMF /= 0) MSHFRM=IMF
        NCUTB_SAVE=NCUTB
        IF (TRCINT) THEN
          WRITE (iunout,*) ' NFLA,NCUTB,NCUTL,IMF,NTRFRM = ',
     .                       NFLA,NCUTB,NCUTL,IMF,NTRFRM
          WRITE (iunout,*) ' IPLS,IFLB(IPLS),FCTE(IPLS),BMASS(IPLS)'
        ENDIF
        DO 20 IPL=1,NPLSI
          READ (IUNIN,'(2I6,2E12.4)') I,IFLB(IPL),FCTE(IPL),BMASS(IPL)
          IF (TRCINT) WRITE (iunout,'(2I6,2E12.4)')
     .                              IPL,IFLB(IPL),FCTE(IPL),BMASS(IPL)
   20   CONTINUE
        READ (IUNIN,'(2I6)') NDXA,NDYA
        IF (TRCINT) WRITE (iunout,*) 'NDXA,NDYA= ',NDXA,NDYA
C  NUMBER OF TARGET SOURCES ON B2 SURFACES: NTARGI
        READ (IUNIN,'(I6)') NTARGI
        IF (TRCINT) WRITE (iunout,*) 'NTARGI=    ',NTARGI
        CALL EIRENE_LEER(1)
        IF (NTARGI.GT.NSTEP) THEN
          CALL EIRENE_MASPRM ('NSTEP',5,NSTEP,'NTARGI',6,NTARGI,IERROR)
          WRITE (iunout,*) 'EXIT CALLED FROM SUBR. INFCOP '
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
C  NUMBER OF PARTS PER TARGET SOURCE
        IF (NTARGI.GT.0) READ (IUNIN,'(12I6)') (NTGPRT(IT),IT=1,NTARGI)
        DO 22 IT=1,NTARGI
          IF (NTGPRT(IT).GT.NGITT) THEN
            NTGPRI=NTGPRT(IT)
            CALL EIRENE_MASPRM('NGITT',5,NGITT,'NTGPRT',6,NTGPRI,IERROR)
            WRITE (iunout,*) 'EXIT CALLED FROM SUBR. INFCOP '
            CALL EIRENE_EXIT_OWN(1)
          ENDIF
   22   CONTINUE
        IREAD=0
C  ALL INDICES: AFTER INDEX MAPPING
C  NDT: INDEX OF X-CELL (EAST OR NORTH SURFACE OF BRAAMS CELL) OF TARGET
C  NINCT: DIRECTION OF OUTER TARGET NORMAL WITH RESPECT TO POSITIVE DIR.
C  NIXY: SOURCE ON Y SURFACE: NIXY=1; SOURCE ON X SURFACE: NIXY=2
C  NTIN,NTEN: SOURCE RANGE FROM GRIDPOINT NTIN TO GRIDPOINT NTEN
        IF (TRCINT)
     .  WRITE (iunout,*) '    IT,  NDT,NINCT, NIXY, NTIN, NTEN',
     .              ',NIFLG, NPTC, NPTCM,NSPZI,NSPZE,NEMOD'
        DO 30 IT=1,NTARGI
          DO 33 IPRT=1,NTGPRT(IT)
            CALL EIRENE_SKIP_READ_COMMENT(IREAD,IUNIN,ZEILE)
            READ (ZEILE,'(12I6)') I,NDT(IT,IPRT),NINCT(IT,IPRT),
     .                              NIXY(IT,IPRT),NTIN(IT,IPRT),
     .                              NTEN(IT,IPRT),NIFLG(IT,IPRT),
     .                              NPTC(IT,IPRT),NPTCM(IT,IPRT),
     .                              NSPZI(IT,IPRT),NSPZE(IT,IPRT),
     .                              NEMOD(IT,IPRT)
            IREAD=0
            NSPZI(IT,IPRT)=MAX0(1,NSPZI(IT,IPRT))
            NSPZE(IT,IPRT)=MIN0(NFLA,NSPZE(IT,IPRT))
            IF (NSPZE(IT,IPRT).LT.NSPZI(IT,IPRT)) THEN
              WRITE (iunout,*) 'WARNING FROM INFCOP: '
              WRITE (iunout,*) 'ITARG,IPRT : ',IT,IPRT
              WRITE (iunout,*) 'NSPZI,NSPZE MODIFIED TO 1,NFLA, RESP.'
              NSPZI(IT,IPRT)=1
              NSPZE(IT,IPRT)=NFLA
            ENDIF
            IF (TRCINT)
     .      WRITE (iunout,'(1X,7I6,2I7,3I6)')
     .                               IT,NDT(IT,IPRT),NINCT(IT,IPRT),
     .                               NIXY(IT,IPRT),NTIN(IT,IPRT),
     .                               NTEN(IT,IPRT),NIFLG(IT,IPRT),
     .                               NPTC(IT,IPRT),NPTCM(IT,IPRT),
     .                               NSPZI(IT,IPRT),NSPZE(IT,IPRT),
     .                               NEMOD(IT,IPRT)
            IF (NIXY(IT,IPRT).EQ.1) THEN
              IF (NTIN(IT,IPRT).LE.0.OR.NTIN(IT,IPRT).GE.NR1ST.OR.
     .            NTEN(IT,IPRT).GT.NR1ST) THEN
                WRITE (iunout,*) 'ERROR IN INPUT BLOCK 14, NTIN, NTEN '
                CALL EIRENE_EXIT_OWN(1)
              ENDIF
            ELSEIF (NIXY(IT,IPRT).EQ.2) THEN
              IF (NTIN(IT,IPRT).LE.0.OR.NTIN(IT,IPRT).GE.NP2ND.OR.
     .            NTEN(IT,IPRT).GT.NP2ND) THEN
                WRITE (iunout,*) 'ERROR IN INPUT BLOCK 14, NTIN, NTEN '
                CALL EIRENE_EXIT_OWN(1)
              ENDIF
            ENDIF
   33     CONTINUE
          IF (TRCINT) CALL EIRENE_LEER(1)
   30   CONTINUE
        READ (IUNIN,'(6E12.4)')  CHGP,CHGEE,CHGEI,CHGMOM
        IF (TRCINT) CALL EIRENE_MASR4
     .                         ('CHGP,CHGEE,CHGEI,CHGMOM         ',
     .                           CHGP,CHGEE,CHGEI,CHGMOM)
C  READ ADDITIONAL DATA TO BE TRANSFERRED FROM B2_TRIA INTO EIRENE
C  HERE: B2_TRIA VOLUME TALLIES
        READ (IUNIN,'(I6)') NAINB
C  ADDITIONAL INPUT TALLY ADIN:  ITAL=12
        NAIN = MAX(NAIN,NAINB)
        CALL EIRENE_ALLOC_CCOUPL(2)
        WRITE (iunout,*) '        NAINI = ',NAINB
        IF (NAINB.GT.NAIN) THEN
          CALL EIRENE_MASPRM ('NAIN',4,NAIN,'NAINB',5,NAINB,IERROR)
          WRITE (iunout,*) 'EXIT CALLED FROM SUBR. INFCOP '
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
        IF (TRCINT.AND.NAINB.GT.0)
     .      WRITE (iunout,*) 'I,NAINS(IAIN),NAINT(IAIN)'
        DO 40 IAIN=1,NAINB
          READ (IUNIN,'(6I6)') I,NAINS(IAIN),NAINT(IAIN)
          READ (IUNIN,'(A72)') TXTPLS(IAIN,12)
          READ (IUNIN,'(2A24)') TXTPSP(IAIN,12),TXTPUN(IAIN,12)
          IF (TRCINT) THEN
            WRITE (iunout,'(6I6)') I,NAINS(IAIN),NAINT(IAIN)
            WRITE (iunout,'(1X,A72)') TXTPLS(IAIN,12)
            WRITE (iunout,'(1X,2A24)') TXTPSP(IAIN,12),TXTPUN(IAIN,12)
          ENDIF
   40   CONTINUE
C  READ ADDITIONAL DATA TO BE TRANSFERRED FROM EIRENE INTO B2
C  HERE: EIRENE SURFACE TALLIES
        READ (IUNIN,'(I6)') NAOTB
        WRITE (iunout,*) '        NAOTI = ',NAOTB
        IF (NAOTB.GT.NLIMPS) THEN
          CALL EIRENE_MASPRM ('NLIMPS',6,NLIMPS,'NAOTB',5,NAOTB,IERROR)
          WRITE (iunout,*) 'EXIT CALLED FROM SUBR. INFCOP '
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
        IF (TRCINT.AND.NAOTB.GT.0)
     .      WRITE (iunout,*) 'I,NAOTS(IAOT),NAOTT(IAOT)'
        DO 50 IAOT=1,NAOTB
          READ (IUNIN,'(6I6)') I,NAOTS(IAOT),NAOTT(IAOT)
          IF (TRCINT) THEN
            WRITE (iunout,'(6I6)') I,NAOTS(IAOT),NAOTT(IAOT)
          ENDIF
   50   CONTINUE
      ENDIF
C
C READING BLOCK 14 FROM FORMATTED INPUT FILE (IUNIN) FINISHED
C
C
C  DEFINE ADDITIONAL TALLIES FOR COUPLING (UPDATED IN SUBR. UPTCOP
C                                              AND IN SUBR. COLLIDE)

cdr  already done in if0prm. Hidden link, must be removed....
      NCPVI = NPLSI
      NCPV  = MAX(NCPV,NCPVI)
C
C SAVE SOME MORE INPUT DATA FOR SHORT CYCLE ON COMMON CCOUPL
      LNLPLG=NLPLG
      LNLDRF=NLDRFT
      LTRCFL=TRCFLE
      NSTRI=NSTRAI
      DO 60 ISTR=1,NSTRAI
        LNLVOL(ISTR)=NLVOL(ISTR)
   60 CONTINUE
      NMODEI=NMODE
      NFILNN=NFILEN
C
C  DEFINE ADDITIONAL TALLIES FOR COUPLING (UPDATED IN SUBR. UPTCOP
C                                              AND IN SUBR. COLLIDE)
      IF (NCPVI.EQ.0) GOTO 70

CDR  SET THE NCPVI= NPLSI COUPLE TALLIES

      DO IPLS=1,NPLSI
        ICPVE(IPLS)=1
        ICPRC(IPLS)=1
        TXTTAL(IPLS,NTALM)=
     .  'ENERGY-WEIGHTED CX RATE OF ATOMS WITH IPLS                  '
        TXTSPC(IPLS,NTALM)=TEXTS(NSPAMI+IPLS)
        TXTUNT(IPLS,NTALM)='AMP                     '
      ENDDO
C
   70 CONTINUE
C
C
C  TRANSFER GEOMETRY
C
      IF (.NOT.(INDGRD(1).EQ.6.OR.INDGRD(2).EQ.6.OR.INDGRD(3).EQ.6))
     .RETURN
C
      OPEN (UNIT=29,ACCESS='SEQUENTIAL',FORM='FORMATTED')
      REWIND 29
C
      OPEN (UNIT=30,ACCESS='SEQUENTIAL',FORM='FORMATTED')
      REWIND 30
C
C  READ IN DATA TO SET UP GEOMETRY FOR NEUTRAL GAS TRANSPORT CODE
C  STATEMENT NUMBER 1000 ---> 1999
C
C  AT PRESENT THE DATA COME FROM THE FILE FT30
C  THIS PART WILL HAVE TO BE MODIFIED AS SOON AS BRAAMS PROVIDES
C  CELL VERTICES AND CUT DESCRIPTION
C
 1000 CONTINUE
C
C  ACTUAL MESH SIZE USED IN THIS RUN: FIRST CARD OF GEOMETRY DATA FILE
C
C
      IF (NDYA.NE.NR1ST-1.OR.NDYA.GT.NDY) THEN
        WRITE (iunout,*) ' PARAMETER ERROR DETECTED IN INTFCE '
        WRITE (iunout,*) ' NDYA MUST BE = NR1ST-1 AND <= NDY'
        WRITE (iunout,*) ' NDYA,NR1ST-1,NDY = ',NDYA,NR1ST-1,NDY
        CALL EIRENE_EXIT_OWN(1)
      ELSEIF (NDXA.NE.NP2ND-1.OR.NDXA.GT.NDX) THEN
        WRITE (iunout,*) ' PARAMETER ERROR DETECTED IN INTFCE '
        WRITE (iunout,*) ' NDXA MUST BE = NP2ND-1 AND <= NDX'
        WRITE (iunout,*) ' NDXA,NP2ND-1,NDX = ',NDXA,NP2ND-1,NDX
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
C
C  EACH FLUXSURFACE IS GIVEN BY A POLYGON OF LENGTH NDXA+1, I.E.
C  WITH NDXA SEGMENTS. THERE ARE NDYA+1 POLYGONS
C  READ IN POLYGONS CELL BY CELL. IX IS INDEX ALONG ONE POLYGON
C                                 IY IS INDEX PERP. TO THE POLYGONS
C
C
C         DIRECTION  OF INCREASING IY ("RADIAL")
C PERP.      ^                 ^ PERP.
C POLYG.NO.IX|      (VV,SY,    | POLYG.NO.IX+1
C            |       FNIY)     |
C            |       ^         |
C            |       |         |
C         X3,Y3      |        X4,Y4
C            ________X__________   -----> ALONG POLYGON NO. (IY+1)
C            |  CELL NO.(IX,IY)|
C            |                 |
C            |       X         X-----> (UU,UP,FNIX,SX)
C            |                 |
C            | (TE,TI,NI,PR,RR,|
C            |  VOL,GX,GY)     |
C            -------------------   -----> ALONG POLYGON NO. (IY)
C         X1,Y1               X2,Y2
C                                         DIRECTION OF INCREASING IX ("POLOIDAL")
C
C
      NPOINT = 0
      IF (.NOT.ALLOCATED(PUX)) THEN
c  unit vector parallel to B field, in poloidal section
        ALLOCATE (PUX(NRAD))
        ALLOCATE (PUY(NRAD))
c  only for inclined target option:
        ALLOCATE (PUXE(NRAD))
        ALLOCATE (PUYE(NRAD))
        ALLOCATE (PUXN(NRAD))
        ALLOCATE (PUYN(NRAD))
c  unit vector perp. to B field, in poloidal section
        ALLOCATE (PVX(NRAD))
        ALLOCATE (PVY(NRAD))
c  only for inclined target option:
        ALLOCATE (PVXE(NRAD))
        ALLOCATE (PVYE(NRAD))
        ALLOCATE (PVXN(NRAD))
        ALLOCATE (PVYN(NRAD))
        PUXE = 0._DP
        PUYE = 0._DP
        PUXN = 0._DP
        PUYN = 0._DP
        PVXE = 0._DP
        PVYE = 0._DP
        PVXN = 0._DP
        PVYN = 0._DP
      END IF
C
      CALL EIRENE_GEOMD (NDXA,NDYA,NPLP,NR1ST,
     .                   PUX,PUY,PVX,PVY,MSHFRM)
C
      IF (NDXA+1.NE.NRPLG) THEN
        WRITE (iunout,*) 'ERROR IN INFCOP: NRPLG.NE.NDXA+1'
        WRITE (iunout,*) 'NDXA+1,NRPLG ',NDXA+1,NRPLG
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
      IF (NPLP.NE.NPPLG) THEN
        WRITE (iunout,*) 'ERROR IN INFCOP: NPPLG.NE.NPLP'
        WRITE (iunout,*) 'NPLP,NPPLG ',NPLP,NPPLG
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
C
      CALL EIRENE_LEER(1)
      ALLOCATE(LLCUT(0:NDXP))
      DO IX=0,NDXP
        LLCUT(IX)=.FALSE.
        DO IPRT=1,NPLP-1
          IXI=NPOINT(2,IPRT)
          IXE=NPOINT(1,IPRT+1)
          IF (IX.GE.IXI.AND.IX.LT.IXE) THEN
            LLCUT(IX)=.TRUE.
            WRITE (iunout,*) 'POLOIDAL CUT CELL INTRODUCED AT IP= ',IX
          ENDIF
        ENDDO
      ENDDO
      CALL EIRENE_LEER(1)
CDR UP TO NOW: UNDERLYING STRUCTURED (COARSE) GRID OF POLYGONS
      levgeo = 3
      nr1stm = nr1st-1
      call eirene_sneigh
C
C  READ ADDITIONAL GEOMETRICAL DATA (MESH DISTORTION, DEAD CELLS)
C  SAME FORMAT AS FORT.31, I.E., INDEX MAPPING MAY BE NECESSARY
      READ (29,'(A)',IOSTAT=IO29) CHR
      IF (IO29.EQ.0) THEN
        REWIND 29
        NRED=(NPPLG-1)*(NCUTL-NCUTB)
        NDX2=NDXA-NRED

        ALLOCATE (ALPHXB(0:NDXP,0:NDYP))
        ALLOCATE (ALPHYB(0:NDXP,0:NDYP))
        ALLOCATE (XAISO(0:NDXP,0:NDYP))
        ALLOCATE (IAISO(0:NDXP,0:NDYP))

C  DEFAULT: ALL CELLS ARE VALID
        IAISO = 1

        CALL EIRENE_PLASM (29,NDX2,NDYA,1,NDX,NDY,1,ALPHXB)
        CALL EIRENE_PLASM (29,NDX2,NDYA,1,NDX,NDY,1,ALPHYB)

        FORM=REPEAT(' ',20)
        FORM(1:10) = '(      I1)'
        WRITE (FORM(2:7),'(I6)') NDX2+2
        DO IY=NDYA+1,0,-1
          READ (29,FORM) (IAISO(IX,IY),IX=0,NDX2+1)
        ENDDO
C
        IF (NCUTL.EQ.NCUTB) GOTO 1020
C
C CONVERT IAISO TO REAL, FOR INDMAP
        XAISO=IAISO
        CALL EIRENE_INDMAP (XAISO,DUMMY,NDX,NDY,1,NDXA,NDYA,1,NCUTB,
     .               NCUTL,NPOINT,NPLP)
C IAISO BACK TO INTEGER
        IAISO=XAISO
C
        CALL EIRENE_INDMAP (ALPHXB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,NCUTB,
     .               NCUTL,NPOINT,NPLP)
        CALL EIRENE_INDMAP (ALPHYB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,NCUTB,
     .               NCUTL,NPOINT,NPLP)
 1020   CONTINUE
C
!  ALPHXB, ALPHYB GIVE THE DIRECTION OF THE B-FIELD IN THE
!  CARTESIAN PLANE
        write (iunout,*) 'testoutput from fort.29 in infcop'
        write (iunout,*) 'irad,ipol, angles.....'
        DO IY=1,NDYA
          DO IX =1,NDXA
            IN=IY+(IX-1)*NR1ST
            ALE=ALPHXB(IX,IY)
            ALW=ALPHXB(IX-1,IY)
            IF (MAX(ALE,ALW)-MIN(ALE,ALW) > PIA) THEN
              write (iunout,*) 'modulus 2PI used', ale,alw
              AL=MIN(ALE,ALW)
              ALW=MAX(ALE,ALW)
              ALE=AL+PI2A
              write (iunout,*) 'new values ale,alw ',ale,alw
            END IF
            ALN=ALPHYB(IX,IY)
            ALS=ALPHYB(IX,IY-1)
! cell-centered angle of B_pol (psi-contour line) against eirene x-coordinate
            ALX=0.25D0*(ALE+ALW+ALN+ALS)
c           write (iunout,'(1x,2i3,1P,5e12.3)') iy,ix,
c    .                                          ale,alw,aln,als,alx
! cell-centered unit vector along poloidal direction
            PUX(IN)=COS(ALX)
            PUY(IN)=SIN(ALX)
! cell-centered unit vector along "radial" (grad psi) direction,
!                    strictly orthonormal to  PU (poloidal) direction
            PVX(IN)=-PUY(IN)
            PVY(IN)=PUX(IN)
! surface-centered: nothing to be done, the values on fort.29 are already surface-centered
!                   on east and north sides of a cell, respectively
            PUXE(IN)=COS(ALE)
            PUYE(IN)=SIN(ALE)
            PUXN(IN)=COS(ALN)
            PUYN(IN)=SIN(ALN)
! again: radial (grad PSI) unit vector, strictly orthonormal to PU, by construction
            PVXE(IN)=-SIN(ALE)
            PVYE(IN)=COS(ALE)
            PVXN(IN)=-SIN(ALN)
            PVYN(IN)=COS(ALN)
          END DO
        END DO
C
        DEALLOCATE (ALPHXB)
        DEALLOCATE (ALPHYB)
        DEALLOCATE (XAISO)
        DEALLOCATE (IAISO)
C
      ELSE
        CALL EIRENE_LEER(1)
        WRITE (iunout,*)
     .    ' NO FILE FORT.29 WITH MODIFIED GRID INFO. FOUND '
        WRITE (iunout,*) ' OLD VERSION CALCULATION MAGN. FIELD FROM ',
     .                   ' GRID IS USED '
        WRITE (iunout,*) ' GRID IS ASSUMED TO BE ORTHOGONAL '
        WRITE (iunout,*) ' NO INFO RE. ISOLATED CELLS FROM THIS FILE '
        CALL EIRENE_LEER(1)
      END IF
C
C  TRANSFER FLAGS
C
      NAINI=NAINB
C
CTRIG A
C  READ DATA FOR TRIANGULAR MESH
C
      OPEN (UNIT=33,ACCESS='SEQUENTIAL',FORM='FORMATTED')
      OPEN (UNIT=34,ACCESS='SEQUENTIAL',FORM='FORMATTED')
      OPEN (UNIT=35,ACCESS='SEQUENTIAL',FORM='FORMATTED')
C
      READ(33,*) NRKNOT
      WRITE(iunout,*) 'NRKNOT = ',NRKNOT


C
C     READ IN THE NUMBER OF TRIANGLES AND ATTRIBUTES OF THE TRIANGLES
      READ(34,*) NTRII
      WRITE(iunout,*) 'NTRII = ',NTRII

C
C  EACH ELEMENT (TRIANGLE) IS GIVEN BY 3 POINTS
C
C                   3
C                 /   \
C            3  /       \  2
C             /           \
C           /               \
C          1.................2
C                  1
C
C
      IF (NTRFRM == 0) THEN
        READ(33,*) (XTRIAN(I),I=1,NRKNOT)
        READ(33,*) (YTRIAN(I),I=1,NRKNOT)

         if (plidl) then
c write file 'triang_new.npco_char' for triang-grid, for idl tool, in appropriate format.
c first: fetch a free file unit number
           open(newunit=jun,file='triang_new.npco_char',
     .          access='SEQUENTIAL',form='FORMATTED')
c next: write file 'triang_new.npco_char'
           write (jun,'(i9)') NRKNOT
           DO I=1,NRKNOT
             WRITE(jun,'(i9,2es24.16)') I,XTRIAN(I),YTRIAN(I)
           ENDDO
           close (unit=jun)
         end if
      ELSE
        DO I=1,NRKNOT
          READ(33,*) J,XTRIAN(I),YTRIAN(I)
        ENDDO
      END IF
C
      IF (NTRII.GT.NRAD.OR.NTRII.GT.NTRI) THEN
        WRITE (iunout,*) ' PARAMETER ERROR DETECTED IN INFUSR '
        WRITE (iunout,*) ' NTRII MUST BE < NRAD AND <= NTRI'
        WRITE (iunout,*) ' NTRII,NRAD,NTRI = ',NTRII,NRAD,NTRI
        CALL EIRENE_EXIT_OWN(1)
      ENDIF

      DO I=1,NTRII
        READ(34,*) J,NECKE(1,I),NECKE(2,I),NECKE(3,I)
      ENDDO

      READ (35,*) IDUMMY
      IF (IDUMMY /= NTRII) THEN
        WRITE (IUNOUT,*) ' NUMBERS OF TRIANGLES DO NOT MATCH '
        WRITE (IUNOUT,*) ' IN FILES ...ELEMENTE AND ...NEIGHBORS'
        WRITE (IUNOUT,*) ' CHECK THE GEOMETRY, EXIT CALLED '
        CALL EIRENE_EXIT_OWN(1)
      END IF

      DO I=1,NTRII
cdr  june 17:
cdr: careful: I ne J possible. Unless triangles are sorted as J= 1,2,3... on fort.35
cdr           This is implicitly assumed here ??
          READ(35,*) J,NCHBAR(1,I),NSEITE(1,I),IDUMMY,
     >                 NCHBAR(2,I),NSEITE(2,I),IDUMMY,
     >                 NCHBAR(3,I),NSEITE(3,I),IDUMMY,
     >                 IXTRI(I),IYTRI(I)

C       WRITE (iunout,*) J,NECKE(1,J),NECKE(2,J),NECKE(3,J),
C    >                   NCHBAR(1,J),NSEITE(1,J),
C    >                   NCHBAR(2,J),NSEITE(2,J),NCHBAR(3,J),NSEITE(3,J)

C THE SPECIAL SURFACE PROPERTY (IF ANY) IS ON INMTI ARRAY, AND TRANSFERED INTO
C EIRENE VIA COMMON.
      ENDDO

C  FOR ALL QUADRANGLES BUILD LIST OF TRIANGLES BELONGING
C  TO THE QUADRANGLE
      IF(.NOT.ALLOCATED(HEADS)) ALLOCATE (HEADS(N1ST,N2ND))
      DO IR=1,NR1ST
        DO IP=1,NP2ND
          NULLIFY(HEADS(IR,IP)%P)
        ENDDO
      ENDDO

      DO ITRI=1,NTRII
        IF (IXTRI(ITRI).GT.0) THEN
          IR=IYTRI(ITRI)
          IP=IXTRI(ITRI)
          ALLOCATE(CURPOI)
          CURPOI%TRIANGLE = ITRI
          CURPOI%NEXT => HEADS(IR,IP)%P
          HEADS(IR,IP)%P => CURPOI
        ENDIF
      ENDDO

C  BUILD NSTGRD ARRAY OF "BLOCKED" TRIANGLES FROM XAISO ARRAY FROM FORT.29
      IF (IO29.EQ.0) THEN
        DO ITRI=1,NTRII
          IY=IYTRI(ITRI)
          IX=IXTRI(ITRI)
          IF (IX .GT. 0) THEN
            IN=IY+(IX-1)*NR1ST
            NSTGRD(ITRI)=ABS(XAISO(IX,IY)-1.)
          ENDIF
        ENDDO
        DEALLOCATE (XAISO)
        DEALLOCATE (IAISO)
      ENDIF
C
C
C  DETERMINE THE ARRAY INMTI FOR ALL NON-DEFAULT STD. SURFACES
C  ISTS=INMTI(ISIDE,NRCELL), ISIDE=1, 2, OR 3
C
      ICOG = 0
      DO ISTS=1,NSTSI

C  FIRST: RADIAL SURFACES

        DO IR=1,NR1ST
          IF (IR.EQ.INUMP(ISTS,1)) THEN
            IR1=IR+1
            IF (IR1.GT.NR1ST) IR1=IR-1
            DO IP=IRPTA(ISTS,2),IRPTE(ISTS,2)-1
              IF (IR.LT.NR1ST) THEN
                CURPOI => HEADS(IR,IP)%P
              ELSE
                CURPOI => HEADS(IR1,IP)%P
              ENDIF
              DO WHILE (ASSOCIATED(CURPOI))
                ITRI=CURPOI%TRIANGLE
CVKG TO FIX A BUG WITH GEOMETRY
                ISC1=0
                ISC2=0
                DO IS=1,3
                  IF(EIRENE_POINT_ON_INTERVAL(XTRIAN(NECKE(IS,ITRI)),
     f                             YTRIAN(NECKE(IS,ITRI)),
     f                             XPOL(IR,IP),YPOL(IR,IP),
     f                             XPOL(IR,IP+1),YPOL(IR,IP+1))) THEN
                  IF(ISC1.GT.0) THEN
                    ISC2=IS
                  ELSE
                    ISC1=IS
                  END IF
                END IF
              ENDDO

C  NODES ISC1 AND ISC2 OF TRIANGLE ITRI ARE LOCATED ON RADIAL SURFACE IR
C  THAT MEANS  SIDE "NUMSI" OF TRIANGLE "ITRI" BELONGS TO NDS
              IF (ISC1.GT.0.AND.ISC2.GT.0) THEN
                NUMSI=MIN(ISC1,ISC2)
                IF (NUMSI.EQ.1.AND.MAX(ISC1,ISC2).EQ.3) NUMSI=3
C
                  ICOG=ICOG+1
                  INSPAT(NUMSI,ITRI)=ICOG
                  INMTI(NUMSI,ITRI)=NLIM+ISTS
                  NBAR=NCHBAR(NUMSI,ITRI)
                  IF (NBAR.GT.0) THEN
                    NBARSI=NSEITE(NUMSI,ITRI)
                    ICOG=ICOG+1
                    INSPAT(NBARSI,NBAR)=ICOG
                    INMTI(NBARSI,NBAR)=NLIM+ISTS
                    LXSRF=.FALSE.
                    CALL CORRECTNSS !VK
                  ENDIF
                ENDIF
                CURPOI => CURPOI%NEXT
              ENDDO
            ENDDO
          ENDIF
        ENDDO

C  NEXT: POLOIDAL SURFACES

        DO IP=1,NP2ND
          IF (IP.EQ.INUMP(ISTS,2)) THEN
            IP1=IP+1
            IF (IP1.GT.NP2ND) IP1=IP-1
            DO IR=IRPTA(ISTS,1),IRPTE(ISTS,1)-1
              IF (IP.LT.NP2ND) THEN
                CURPOI => HEADS(IR,IP)%P
              ELSE
                CURPOI => HEADS(IR,IP1)%P
              ENDIF
              DO WHILE (ASSOCIATED(CURPOI))
                ITRI=CURPOI%TRIANGLE
CVKG TO FIX GEOMETRY BUG
                ISC1=0
                ISC2=0
                DO IS=1,3
                  IF(EIRENE_POINT_ON_INTERVAL(XTRIAN(NECKE(IS,ITRI)),
     f                             YTRIAN(NECKE(IS,ITRI)),
     f                             XPOL(IR,IP),YPOL(IR,IP),
     f                             XPOL(IR+1,IP),YPOL(IR+1,IP))) THEN
                  IF(ISC1.GT.0) THEN
                    ISC2=IS
                  ELSE
                    ISC1=IS
                  END IF
                END IF
              ENDDO

C  NODES ISC1 AND ISC2 OF TRIANGLE ITRI ARE LOCATED ON POLOIDAL SURFACE IP
C  THAT MEANS  SIDE "NUMSI" OF TRIANGLE "ITRI" BELONGS TO NDS
              IF (ISC1.GT.0.AND.ISC2.GT.0) THEN
                  NUMSI=MIN(ISC1,ISC2)
                  IF (NUMSI.EQ.1.AND.MAX(ISC1,ISC2).EQ.3) NUMSI=3
C
                  ICOG=ICOG+1
                  INSPAT(NUMSI,ITRI)=ICOG
                  INMTI(NUMSI,ITRI)=NLIM+ISTS
                  NBAR=NCHBAR(NUMSI,ITRI)
                  IF (NBAR.GT.0) THEN
                    NBARSI=NSEITE(NUMSI,ITRI)
                    ICOG=ICOG+1
                    INSPAT(NBARSI,NBAR)=ICOG
                    INMTI(NBARSI,NBAR)=NLIM+ISTS
                    LXSRF=.TRUE.
                    CALL CORRECTNSS !VK
                  ENDIF
                ENDIF
                CURPOI => CURPOI%NEXT
              ENDDO
            ENDDO
          ENDIF
        ENDDO

C  NEXT TOROIDAL SURFACES

C  ADDED HERE ONLY FOR OTHER INFCOP, IN CASE OF 3D GRID RESOLUTION. IRRELEVANT FOR B2_TRIA

        DO IT = 1, NT3RD
          IF (IT.EQ.INUMP(ISTS,3)) THEN
            DO ITRI = 1, NTRII
              IR = IXTRI(ITRI)
              IP = IXTRI(ITRI)
              IF ((IR >= IRPTA(ISTS,1)) .AND.
     .            (IR <= IRPTE(ISTS,1)) .AND.
     .            (IP >= IRPTA(ISTS,2)) .AND.
     .            (IP <= IRPTE(ISTS,2))) THEN
                INMTI3(ITRI,IT) = ISTS
              END IF
            END DO
          END IF
        END DO

      ENDDO
C
C  NOW THE ADJUSTMENTS, WHICH ARE AUTOMATICALLY DONE IN GEOUSR
C  (INPUT BLOCK 15, B2-CODE SPECIFIC)
C

      CALL EIRENE_GEOUSR

C
C  DETERMINE THE ARRAY INMTI FOR ALL ADDITIONAL SURFACES
C  ISTS=INMTI(ISIDE,NRCELL), ISIDE=1, 2, OR 3
C
      DO I=1,NLIMI
        IF (IGJUM0(I)==0) THEN
C  SURFACE I IS ACTIVE
          IF (ILPLG(I).NE.0) THEN
C  SURFACE I IS PART IF A CONTOUR USED FOR THE MESHGENERATOR
            VSX=P2(1,I)-P1(1,I)
            VSY=P2(2,I)-P1(2,I)
            VS=SQRT(VSX**2+VSY**2)+EPS60
C
            DO 1111 ITRI=1,NTRII
              DO IS=1,3
                IF (IS.EQ.1) THEN
                  VTX=XTRIAN(NECKE(2,ITRI))-XTRIAN(NECKE(1,ITRI))
                  VTY=YTRIAN(NECKE(2,ITRI))-YTRIAN(NECKE(1,ITRI))
                  ISCS=1
                  ISC1=1
                  ISC2=2
                ELSEIF (IS.EQ.2) THEN
                  VTX=XTRIAN(NECKE(3,ITRI))-XTRIAN(NECKE(2,ITRI))
                  VTY=YTRIAN(NECKE(3,ITRI))-YTRIAN(NECKE(2,ITRI))
                  ISCS=2
                  ISC1=2
                  ISC2=3
                ELSEIF (IS.EQ.3) THEN
                  VTX=XTRIAN(NECKE(1,ITRI))-XTRIAN(NECKE(3,ITRI))
                  VTY=YTRIAN(NECKE(1,ITRI))-YTRIAN(NECKE(3,ITRI))
                  ISCS=3
                  ISC1=3
                  ISC2=1
                ENDIF
                VT=SQRT(VTX**2+VTY**2)+EPS60
C
                VPRO=(VSX*VTY-VTX*VSY)/(VT*VS)
                IF (ABS(VPRO).LT.1.E-2) THEN
C  SURFACES ARE PARALLEL,
C  TEST IF ONE POINT OF THE APPROPRIATE TRIANGLE SIDE BELONGS TO
C  THE SURFACE
                  PX=P1(1,I)
                  PY=P1(2,I)
                  ICOU=1
                  ISC=ISC1
 1112             TX=XTRIAN(NECKE(ISC,ITRI))
                  TY=YTRIAN(NECKE(ISC,ITRI))

                  if(.true.) then
                    l1=eirene_point_on_interval(tx,ty,p1(1,i),p1(2,i),
     .                                                p2(1,i),p2(2,i))
                    TX=XTRIAN(NECKE(ISC2,ITRI))
                    TY=YTRIAN(NECKE(ISC2,ITRI))
                    l2=eirene_point_on_interval(tx,ty,p1(1,i),p1(2,i),
     .                                                p2(1,i),p2(2,i))
                    if(l1 .and. l2) then
                      IGJUM0(I)=1
                      ICOG=ICOG+1
                      INSPAT(ISCS,ITRI)=ICOG
                      INMTI(ISCS,ITRI)=I
                      IF (LCHKQUD)
     .                IREVERS(ISCS,ITRI)=
     .                        INT(SIGN(1._DP,PX*VTRIX(ISCS,ITRI)+
     .                                       PY*VTRIY(ISCS,ITRI)))
                    endif
                  else

                  IF (ABS(VSX).GT.ABS(VSY)) THEN
                    XMUE=(TX-PX)/VSX
                    IF (XMUE.GE.-1.D-5 .AND. XMUE.LE.1.+1.D-5) THEN
                      TEST=(TY-PY-XMUE*VSY)/VS
!pb                      IF (ABS(TEST).LT.1.D-5) THEN
                      IF (ABS(TEST).LT.1.D-4) THEN
                        IF (ICOU.EQ.2) THEN
C  TAKE CORRESPONDING ADDITIONAL SURFACE "I" OUT
C  AND REPLACE IT BY NON-DEFAULT STD. SURFACE
                          IGJUM0(I)=1
                          ICOG=ICOG+1
                          INSPAT(ISCS,ITRI)=ICOG
                          INMTI(ISCS,ITRI)=I
                          IF (LCHKQUD)
     .                    IREVERS(ISCS,ITRI)=
     .                        INT(SIGN(1._DP,PX*VTRIX(ISCS,ITRI)+
     .                                       PY*VTRIY(ISCS,ITRI)))
                          GOTO 1111
                        ELSE
                          ICOU=2
                          ISC=ISC2
                          GOTO 1112
                        ENDIF
                      ENDIF
                    ENDIF
                  ELSE
                    XMUE=(TY-PY)/VSY
                    IF (XMUE.GE.-1.D-5 .AND. XMUE.LE.1.+1.D-5) THEN
                      TEST=(TX-PX-XMUE*VSX)/VS
!pb                      IF (ABS(TEST).LT.1.D-5) THEN
                      IF (ABS(TEST).LT.1.D-4) THEN
                        IF (ICOU.EQ.2) THEN
C  TAKE CORRESPONDING ADDITIONAL SURFACE "I" OUT
C  AND REPLACE IT BY NON-DEFAULT STD. SURFACE
                          IGJUM0(I)=1
                          ICOG=ICOG+1
                          INSPAT(ISCS,ITRI)=ICOG
                          INMTI(ISCS,ITRI)=I
                          IF (LCHKQUD)
     .                    IREVERS(ISCS,ITRI)=
     .                        INT(SIGN(1._DP,PX*VTRIX(ISCS,ITRI)+
     .                                       PY*VTRIY(ISCS,ITRI)))
                          GOTO 1111
                        ELSE
                          ICOU=2
                          ISC=ISC2
                          GOTO 1112
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDIF

                  endif

                ENDIF
              ENDDO
 1111       CONTINUE  ! END OF ITRI LOOP
          ENDIF
        ENDIF
      ENDDO

C
      NLPLG=.FALSE.
      NLPOL=.FALSE.
      NLFEM=.TRUE.
      LEVGEO=4

C  SAVE OLD, STRUCTURED, 2D GRID INFO FOR TALLY OUTPUT, ETC...
      NR1TAL=NR1ST
      NP2TAL=NP2ND
      NT3TAL=NT3RD
      NSURF_TAL=NR1TAL*NP2TAL*NT3TAL*NBMLT

C     NSBOX_TAL, NRADD_TAL: SEE BELOW

C  SET NEW, FINER GRID STRUCTURE
      NR1ST=NTRII+1  !  add one cell, for volume average in TRIA region
      NP2ND=1
      NT3RD=1
      NSURF=NR1ST*NP2ND*NT3RD*NBMLT
      NSBOX=NSURF+NRADD


C  counter for additional cells outside original COARSE standard grid
      icoadd = NSURF_TAL

C  NCLTAL(ITRI):  TRIANGLE ITRI IS PART OF ORIGINAL STRUCTURED GRID CELL IX,IY,
C                 WITH IX,IY, CODED IN THE 1D ARRAY FORM (NCELL) OF EIRENE STANDARD GRIDS
C                 NCELL=NCLTAL(ITRI)
      IF (LCOARSE) THEN
      	write (iunout,*) 'SCORING OF VOLUME-AVERAGED TALLIES  '
        write (iunout,*) 'IS ON COARSE GRID CELLS NCELL ONLY. '
      DO ITRI=1,NTRII
        IY=IYTRI(ITRI)
        IX=IXTRI(ITRI)
        IF (IX .GT. 0) THEN
           IN=IY+(IX-1)*NR1TAL
           NCLTAL(ITRI)=IN
        ELSE   ! ADDITIONAL TRIA CELLS, OUTSIDE OLD STRUCTURED GRID
           icoadd = icoadd+1
           NCLTAL(ITRI)=icoadd
        ENDIF
      ENDDO

      ELSEIF (.NOT.LCOARSE) THEN
        write (iunout,*) 'SCORING OF VOLUME-AVERAGED TALLIES  '
        write (iunout,*) 'IS ON FINE (TRIA) GRID ONLY. '
C                 NCLTAL(ITRI)=ITRI
        DO ITRI=1,NTRII
          IY=IYTRI(ITRI)
          IX=IXTRI(ITRI)
          IF (IX .GT. 0) THEN
            IN=IY+(IX-1)*NR1TAL
            NCLTAL(ITRI)=ITRI
          ELSE   ! ADDITIONAL TRIA CELLS, OUTSIDE OLD STRUCTURED GRID
            icoadd = icoadd+1
            NCLTAL(ITRI)=itri
          ENDIF
        ENDDO
      ENDIF

c   ntrii+1 is the storage for summed/integrated tallies, over the standard grid
c           i.e. not including the additional cell region (if any)
      NCLTAL(NTRII+1) = 0

      ICOSCR= icoadd   ! NUMBER OF SCORING CELLS, SO FAR.
!
      if (nsbox > nsurf) then
c  There are NRADD further additional cells, from input block 2e.
c  These are not part of triangular grid. Add them now to grid
        IF (LCOARSE) THEN
          do itri = nsurf+1, nsbox
            icoscr = icoscr + 1
            ncltal(itri) = icoscr
          end do
        ELSEIF (.NOT.LCOARSE) THEN
          do itri = nsurf+1, nsbox
            icoscr = icoscr + 1
            ncltal(itri) = itri
          end do
        ENDIF
      end if

C  ADDITIONAL CELLS: THOSE CELLS THAT ARE ADDED FROM OUTSIDE ORIGINAL GRID,
C  PLUS THOSE FROM INPUT FILE BLOCK 2E (NRADD)
      NRADD_TAL = (ICOSCR - NSURF_TAL)
C  TOTAL NUMBER OF CELLS OF COARSE GRID: STRUCTURED GRID PLUS ALL ADDITIONAL CELLS
      NSBOX_TAL=NSURF_TAL+NRADD_TAL


C  save old COURSE grid structure

      nr1tal_save=nr1tal
      np2tal_save=np2tal
      nt3tal_save=nt3tal
      nsurf_tal_save=nsurf_tal
      nsbox_tal_save=nsbox_tal
      nradd_tal_save=nradd_tal

      IF (.NOT.LCOARSE) THEN
c  if no coarser grid has been set for scoring:
c  reset scoring grid parameters back to fine tally grid

        nr1tal=nr1st
        np2tal=np2nd
        nt3tal=nt3rd
        nsurf_tal=nsurf
        nsbox_tal=nsbox
        nradd_tal=nradd
      ENDIF

      NGITT = COUNT(INMTI(1:3,1:NTRII) .NE. 0)

C  CARRY OUT SOME CONSISTENCY CHECKS ON NEW TRIAGULAR GRID
      DO ITRI=1,NTRII
        DO IS=1,3
          IF (NCHBAR(IS,ITRI).EQ.0.AND.INMTI(IS,ITRI).EQ.0) THEN
            WRITE (iunout,*) ' ERROR IN INFCOP '
            WRITE (iunout,*) ' OPEN SIDE OF TRIANGLE ',ITRI,' SIDE ',IS

            if (ixtri(itri) == 0) then
              write (iunout,*)
     .               'triangle outside original structured grid'
              call eirene_masj1('itri  ',itri)
            else
              write (iunout,*)
     .               'triangle inside original structured grid'
              call eirene_masj2('nr,np          ',
     .                           iytri(itri),ixtri(itri))
            endif

            write (iunout,*) ' necke ',necke(1:3,itri)
            write (iunout,*) ' xtrian,ytrian(1) ',xtrian(necke(1,itri)),
     .                                       ytrian(necke(1,itri))
            write (iunout,*) ' xtrian,ytrian(2) ',xtrian(necke(2,itri)),
     .                                       ytrian(necke(2,itri))
            write (iunout,*) ' xtrian,ytrian(3) ',xtrian(necke(3,itri)),
     .                                       ytrian(necke(3,itri))
            IS1=IS+1
            IF (IS.EQ.3) IS1=1
            WRITE (iunout,*) ' XTRIAN,YTRIAN ',XTRIAN(NECKE(IS,ITRI)),
     .                                    YTRIAN(NECKE(IS,ITRI))
            WRITE (iunout,*) ' XTRIAN,YTRIAN ',XTRIAN(NECKE(IS1,ITRI)),
     .                                    YTRIAN(NECKE(IS1,ITRI))
          ENDIF
        ENDDO
      ENDDO
C

csw 09jan2012 missing IGJUM3 correction!
CVK CORRECT  IGJUM1 AND IGJUM3 FOR TRIANGLES OUTSIDE STANDARD (B2) MESH


      WRITE(iunout,*) 'IF0COP: CHECKING IGJUM3 FROM SETTINGS IN GEOUSR'
      DO ITRI=1,NTRII

C SET FLAG FOR BOUNDARY CELLS
c ifbound=t means: TRIANGULAR CELL NO. ITRI has AT LEAST ONE neighbor IT, which
c                  is not part of the STANDARD (polygonal) grid

        IFBOUND=.FALSE.
        DO IS=1,3
          IT=NCHBAR(IS,ITRI)
          IF(IT.GT.0) THEN
            IF(IYTRI(IT).LE.0.OR.IXTRI(IT).LE.0) THEN
              IFBOUND=.TRUE.
              EXIT
            END IF
          END IF
        END DO


        DO I=1,NLIM
C SWITCHING ON ADDITIONAL SURFACE CHECKING FOR TRIANGULAR CELL ITRI (OUTSIDE B2 MESH)
          IF(IXTRI(ITRI).LE.0.AND.IYTRI(ITRI).LE.0
     .                      .AND.IGJUM0(I).EQ.0) THEN
            if (IGJUM3(ITRI,I).ne.0) then
              if (trcint.and.trcsur)
     .        write (iunout,*) 'activated itri,i, 1', itri,i
              IGJUM3(ITRI,I)=0
            endif
          END IF
C SWITCHING ON ADDITIONAL SURFACE CHECKING FOR BOUNDARY CELL ITRI OF B2 GRID
          IF(IFBOUND.AND.IGJUM0(I).EQ.0) THEN
            if (IGJUM3(ITRI,I).ne.0) then
              if (trcint.and.trcsur)
     .         write (iunout,*) 'activated itri,i, 2', itri,i
              IGJUM3(ITRI,I)=0
            endif
          END IF
        END DO
      END DO
CVK END
csw



      CALL EIRENE_LEER(2)
      CALL EIRENE_HEADNG(' CASE REDEFINED IN COUPLE_TRIA:  ',33)
      WRITE (iunout,*) 'NLPLG,NLFEM ',NLPLG,NLFEM
      WRITE (iunout,*) 'NLPOL       ',NLPOL

      WRITE (IUNOUT,*) 'NEW (FINE, UN-STRUCTURED) GRID '
      CALL EIRENE_MASJ5('NR1ST,   NP2ND,  NSURF,  NRADD,  NSBOX  ',
     .                   NR1ST,   NP2ND,  NSURF,  NRADD,  NSBOX)

      WRITE (IUNOUT,*) 'OLD (COARSE, STRUCTURED) GRID '
      CALL EIRENE_MASJ5('NR1_TAL,NP2_TAL,NSF_TAL, NRA_TAL,NBX_TAL',
     .                   NR1TAL_SAVE, NP2TAL_SAVE, NSURF_TAL_SAVE,
     .                   NRADD_TAL_SAVE, NSBOX_TAL_SAVE)
      CALL EIRENE_LEER(2)
CTRIG E
      RETURN
C
C   GEOMETRY DEFINITION PART FINISHED
C
      ENTRY EIRENE_IF1COP
C
C   NOW READ THE PLASMA STATE GIVEN BY BRAAMS
C   AT PRESENT THE DATA COME FROM THE FILE FT31
C   FURTHERMORE: SCALING TO EIRENE UNITS AND INDEX MAPPING
C   STATEMENT NO. 2000 ---> 2999
C
C  IN CASE OF "SHORT CYCLE" THE PLASMA STATE IS TRANSFERRED VIA COMMON
C
      LSHORT=.FALSE.
      CALL EIRENE_LEER(1)
      WRITE (iunout,*) 'IF1COP CALLED '
      IF (NLPLAS) THEN
        WRITE (IUNOUT,*) 'PLASMA DATA EXPECTED ON BRAEIR'
      ELSE
        WRITE (IUNOUT,*) 'PLASMA DATA EXPECTED ON FORT.31'
      ENDIF
C  SKIP READING PLASMA, IF NLPLAS
      IF (NLPLAS) GOTO 2100
C
      GOTO 99991
C
C  IN CASE OF "SHORT CYCLE" OR TIME DEP. MODE
C  THE PLASMA STATE IS TRANSFERRED VIA COMMON
C  ONLY SCALING TO EIRENE UNITS AND INDEX MAPPING NEEDS TO BE DONE HERE
C
      ENTRY EIRENE_INTER1
      LSHORT=.TRUE.
      GOTO 2100
C
99991 CONTINUE
C
C  TRANSFER PROFILES
C
      IF (.NOT.(INDPRO(1).EQ.6.OR.INDPRO(2).EQ.6.OR.INDPRO(3).EQ.6.OR.
     .          INDPRO(4).EQ.6)) RETURN
C
C
      OPEN (UNIT=31,ACCESS='SEQUENTIAL',FORM='FORMATTED')
      REWIND 31
C
      IF (NFLA.GT.NFL) THEN
        WRITE (iunout,*) ' PARAMETER ERROR DETECTED IN INFCOP '
        WRITE (iunout,*) ' NFLA MUST BE <= NFL'
        WRITE (iunout,*) ' NFLA,NFL = ',NFLA,NFL
        CALL EIRENE_EXIT_OWN(1)
      ENDIF

      CALL EIRENE_ALLOC_BRAEIR(NDX,NDY,NFL,IFOFF)
C
C  B2-BRAAMS CODE SPECIFIC BEGIN
      NRED=(NPPLG-1)*(NCUTL-NCUTB)
      NDX2=NDXA-NRED
      write(iunout,*) 'NDX2, NDYA, NFLA = ',ndx2,ndya,nfla
      CALL EIRENE_PLASM (31,NDX2,NDYA,NFLA,NDX,NDY,NFL,DNIB)
      CALL EIRENE_PLASM (31,NDX2,NDYA,NFLA,NDX,NDY,NFL,UUB)
      CALL EIRENE_PLASM (31,NDX2,NDYA,NFLA,NDX,NDY,NFL,VVB)
      CALL EIRENE_PLASM (31,NDX2,NDYA,1,NDX,NDY,1,TEB)
      CALL EIRENE_PLASM (31,NDX2,NDYA,1,NDX,NDY,1,TIB)
      CALL EIRENE_PLASM (31,NDX2,NDYA,1,NDX,NDY,1,PRB)
      CALL EIRENE_PLASM (31,NDX2,NDYA,NFLA,NDX,NDY,NFL,UPB)
      CALL EIRENE_PLASM (31,NDX2,NDYA,1,NDX,NDY,1,RRB)
      CALL EIRENE_PLASM (31,NDX2,NDYA,NFLA,NDX,NDY,NFL,FNIXB)
      CALL EIRENE_PLASM (31,NDX2,NDYA,NFLA,NDX,NDY,NFL,FNIYB)
      CALL EIRENE_PLASM (31,NDX2,NDYA,1,NDX,NDY,1,FEIXB)
      CALL EIRENE_PLASM (31,NDX2,NDYA,1,NDX,NDY,1,FEIYB)
      CALL EIRENE_PLASM (31,NDX2,NDYA,1,NDX,NDY,1,FEEXB)
      CALL EIRENE_PLASM (31,NDX2,NDYA,1,NDX,NDY,1,FEEYB)
C
C  OPTIONAL ARRAYS: VOLB, BFELDB
      VOLB = 0.D0
      BFELDB = 0.D0

C  CELL VOLUMES AS USED IN B2
      CALL EIRENE_PLASM (31,NDX2,NDYA,1,NDX,NDY,1,VOLB)
C  MAGNETIC FIELD STRENGTH (TESLA)
      CALL EIRENE_PLASM (31,NDX2,NDYA,1,NDX,NDY,1,BFELDB)

      CALL EIRENE_PLASM (31,NDX2,NDYA,NFLA,NDX,NDY,NFL,VPARXB)
      CALL EIRENE_PLASM (31,NDX2,NDYA,NFLA,NDX,NDY,NFL,VPARYB)
      CALL EIRENE_PLASM (31,NDX2,NDYA,NFLA,NDX,NDY,NFL,VRADXB)
      CALL EIRENE_PLASM (31,NDX2,NDYA,NFLA,NDX,NDY,NFL,VRADYB)
      CALL EIRENE_PLASM (31,NDX2,NDYA,1,NDX,NDY,1,DELTAE_PARXB)
      CALL EIRENE_PLASM (31,NDX2,NDYA,1,NDX,NDY,1,DELTAE_PARYB)
      CALL EIRENE_PLASM (31,NDX2,NDYA,1,NDX,NDY,1,DELTAE_RADXB)
      CALL EIRENE_PLASM (31,NDX2,NDYA,1,NDX,NDY,1,DELTAE_RADYB)
      CALL EIRENE_PLASM (31,NDX2,NDYA,1,NDX,NDY,1,DELTAI_PARXB)
      CALL EIRENE_PLASM (31,NDX2,NDYA,1,NDX,NDY,1,DELTAI_PARYB)
      CALL EIRENE_PLASM (31,NDX2,NDYA,1,NDX,NDY,1,DELTAI_RADXB)
      CALL EIRENE_PLASM (31,NDX2,NDYA,1,NDX,NDY,1,DELTAI_RADYB)
      CALL EIRENE_PLASM (31,NDX2,NDYA,1,NDX,NDY,1,DELTA_SHEATHXB)
      CALL EIRENE_PLASM (31,NDX2,NDYA,1,NDX,NDY,1,DELTA_SHEATHYB)
      CALL EIRENE_PLASM (31,NDX2,NDYA,1,NDX,NDY,1,AISOB)

!pb  for the time being set default values
!csw 26sep2011
      if ((maxval(abs(delta_sheathxb)) < eps30) .and.
     .    (maxval(abs(delta_sheathyb)) < eps30)) then
        delta_sheathxb=3.1
        delta_sheathyb=3.1
      end if
!csw
      if (maxval(aisob) < eps30) then
        aisob = 1.d0
      end if

C
 2100 CONTINUE
C
C  NO INDEX MAPPING REQUIRED, IF NEW TIMESTEP ON SAME PLASMA
C  BRAEIR NOT MODIFIED SINCE LAST CALL TO IF1COP
      IF (NCUTB_SAVE.EQ.NCUTL) THEN
        WRITE (iunout,*) 'NO INDEX MAPPING DONE'
      ELSE
        WRITE (iunout,*) 'DO INDEX MAPPING   ', NCUTB_SAVE,NCUTB,NCUTL
      ENDIF
      CALL EIRENE_LEER(1)
C
C  INDEX MAPPING : NDY DIRECTION
C
C  INDEX MAPPING : NDX DIRECTION
C  SET THE NUMBER OF COLUMNS PER CUT FROM NCUTB (BRAAMS IMPLEMENTATION)
C  TO WHAT IS FOUND FROM THE EIRENE GEOMETRY FILE (NCUTL)
C
      IF (NCUTL.EQ.NCUTB_SAVE) GOTO 2101
C  FIRST THE ZONE-CENTERED DATA
      CALL EIRENE_INDMAP (DNIB,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (TEB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,NCUTB,NCUTL,
     .             NPOINT,NPLP)
      CALL EIRENE_INDMAP (TIB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,NCUTB,NCUTL,
     .             NPOINT,NPLP)
      CALL EIRENE_INDMAP (RRB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,NCUTB,NCUTL,
     .             NPOINT,NPLP)
      CALL EIRENE_INDMAP (PRB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,NCUTB,NCUTL,
     .             NPOINT,NPLP)
C  NOW THE SURFACE-CENTERED DATA
      CALL EIRENE_INDMAP (FNIXB,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (FNIYB,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .             NCUTB,NCUTL,NPOINT,NPLP)
C  distinct from B2.5: these velocities are surface-centered in b2
      CALL EIRENE_INDMAP (UUB,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (VVB,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (UPB,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .             NCUTB,NCUTL,NPOINT,NPLP)
C  same in B2 and in B2.5: these ENERGY fluxes are surface-centered
      CALL EIRENE_INDMAP (FEIXB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (FEIYB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (FEEXB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (FEEYB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,
     .             NCUTB,NCUTL,NPOINT,NPLP)
C  UNUSED INPUT TALLIES
      CALL EIRENE_INDMAP (VOLB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (BFELDB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,
     .             NCUTB,NCUTL,NPOINT,NPLP)

C  ADDITIONAL INPUT TALLIES FOR BOUNDARY CONDITIONS
      CALL EIRENE_INDMAP (VPARXB,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (VPARYB,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (VRADXB,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (VRADYB,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (DELTAE_PARXB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (DELTAE_PARYB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (DELTAE_RADXB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (DELTAE_RADYB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (DELTAI_PARXB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (DELTAI_PARYB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (DELTAI_RADXB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (DELTAI_RADYB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (DELTA_SHEATHXB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (DELTA_SHEATHYB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,
     .             NCUTB,NCUTL,NPOINT,NPLP)
C  FLAGS FOR ISOLATED CELL REGIONS
      CALL EIRENE_INDMAP (AISOB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,
     .             NCUTB,NCUTL,NPOINT,NPLP)
c
      WRITE (iunout,*) 'INDEX MAPPING DONE '
C
 2101 CONTINUE
C
C  INDICATE, THAT NOW BRAEIR CONTAINS DATA AFTER INDEX-MAPPING
      NCUTB_SAVE=NCUTL
C
C  RESET 2D ARRAYS ONTO 1D EIRENE ARRAYS, RESCALE TO EIRENE UNITS
C  AND CONVERT BRAAMS VECTORS INTO CARTESIAN EIRENE VECTORS
C
C  UNITS CONVERSION FACTORS
      T=1./ELCHA
      V=1.D2                             !pb 1.e2 -> 1.d2
      VL=1.D6                            !pb 1.e6 -> 1.d6
CTRIG A
C  VACUUM DATA NEEDED FOR REGION OUTSIDE B2-MESH
      TVAC=0.02_DP                       !pb 0.02 -> 0.02_dp
      DVAC=1.D2
      VVAC=0.
      BVAC=1.
CTRIG E
      DO 2105 IPLS=1,NPLSI
        D(IPLS)=1.D-6*FCTE(IPLS)         !pb 1.e-6 -> 1.d-6
        FL(IPLS)=ELCHA*FCTE(IPLS)        !  SCALING FOR FLUX: 1/S --> AMP
 2105 CONTINUE
C
C  SET PLASMA BACKGROUND ON TRIANGULAR GRID
C  A TRIANGULAR CELL ITRI RECEIVES THE PLASMA DATA FROM ITS LARGER HOST CELL IX,IY,
C  WITHOUT ANY WEIGHTING/INTERPOLATION ETC...
C

      BZINTF = 1._DP
      BFINTF = 1._DP
      DO ITRI=1,NTRII
C  INSIDE ORIGINAL 2D GRID, IX,IY
        IY=IYTRI(ITRI)
        IX=IXTRI(ITRI)
        IF (IX .GT. 0) THEN
          IN=IY+(IX-1)*NR1TAL_SAVE
          TEINTF(ITRI)=TEB(IX,IY)*T
C
C  ONLY ONE ION TEMPERATURE AVAILABLE FROM PLASMA FLUID CODE,
C  SEE LOOP 2150 BELOW
          TIINTF(1,ITRI)=TIB(IX,IY)*T
C
C  poloidal field
          BX=PUX(IN)*RRB(IX,IY)   ! +PVX(IN)*0., but radial field is zero
          BY=PUY(IN)*RRB(IX,IY)   ! +PVY(IN)*0.
c  toroidal field
          BZ=SQRT(1.-RRB(IX,IY)**2)
c  normalize B-field vector to length 1 (one)
c
          BN=SQRT(BX*BX+BY*BY+BZ*BZ)
          BXINTF(ITRI)=BX/BN
          BYINTF(ITRI)=BY/BN
          BZINTF(ITRI)=BZ/BN
          BFINTF(ITRI)=BN
c
          VLINTF(ITRI)=VOLB(IX,IY)*VL
        ELSE
c  outside the original b2 grid:
C    set default vacuum temperatures TVAC
C    set default b-field: (0,0,1)
          TEINTF(ITRI)=TVAC
          TIINTF(1,ITRI)=TVAC
          BXINTF(ITRI)=0.
          BYINTF(ITRI)=0.
          BZINTF(ITRI)=1.
          BFINTF(ITRI)=1.
c
C         VLINTF(ITRI)=1.
        ENDIF
      ENDDO  !ITRI LOOP

C  ntrii+1:  for average over ntrii grid

C  additional cells from input block 2e
      DO itri=ntrii+2, nsbox
        TEINTF(ITRI)=TVAC
        TIINTF(1,ITRI)=TVAC
        BXINTF(ITRI)=0.
        BYINTF(ITRI)=0.
        BZINTF(ITRI)=1.
        BFINTF(ITRI)=1.
c
C       VLINTF(ITRI)=1.
      ENDDO
C
C  SET SAME ION TEMPERATURE FOR ALL EIRENE BACKGROUND SPECIES
C
      DO 2150 IPLSTI=1,NPLSTI
      DO 2150 ITRI=1,NSBOX
        TIINTF(IPLSTI,ITRI)=TIINTF(1,ITRI)
 2150 CONTINUE
C
CDR  set density from B2 array DNIB, for each fluid
CDR  set plasma flow velocity field from B2 arrays UPB (parallel velocity)
c  without drifts:
c  upb * pitch:  poloidal velocity (i.e. cartesian x,y direction).
c  poloidal field direction is given by that of the poloidal cell face PU..(in),
C  i.e. along a flux surface. (PU(...) is cell-centered)
c  and upb*(1-pitch^2): toroidal velocity  (i.e. cartesian z direction (nltrz) or
c                                                toroidal phi direction (nltra)
c  sign of flowfield follows the sign of poloidal grid in B2.
c
c  with drifts:
c  uudia and vvdia are additional flow velocities in b2.5 only.
c
c
c
c
      IREAD=0
      DO 2200 IPLS=1,NPLSI
        IF (IFLB(IPLS).GT.0) THEN  ! DEAL WITH B2 ION SPECIES ONLY, EXCLUDE VIRTUAL EIRENE BACKGROUND
          IPLSV=MPLSV(IPLS)
          DO 2201 IFL=1,NFLA
            IF (IFLB(IPLS).NE.IFL) GOTO 2201
            DO ITRI=1,NTRII
              IY=IYTRI(ITRI)
              IYM1=IY-1
              IX=IXTRI(ITRI)
              IF (IX .GT. 0) THEN
                IXM1 = IX-1
                IF (LLCUT(IXM1)) THEN
                  IXM1 = NGHPOL(4,IY,IX)
                END IF

C  SET PLASMA DENSITY
                IN=IY+(IX-1)*NR1TAL_SAVE
                DIINTF(IPLS,ITRI)=DNIB(IX,IY,IFL)*D(IPLS)
C  INTERPOLATE PLASMA VELOCITIES TO CELL CENTERS.
                UUBC=0.5*(UUB(IXM1,IY,IFL)+UUB(IX,IY,IFL))
                UPBC=0.5*(UPB(IXM1,IY,IFL)+UPB(IX,IY,IFL))
cdr  radial velocity
                VVBC=0.5*(VVB(IX,IYM1,IFL)+VVB(IX,IY,IFL))

cdr  rrb: pitch  B_pol/B_total
C  UDBC: DIAMAGNETIC VELOCITY, SHOULD BE ZERO IN B2, AND NONZERO IN EB2
C        NOTE: IN LINEAR DEVICES: RRB=1
                IF (RRB(IX,IY).LT.1.D0) THEN
                  RBC=SQRT(1.-RRB(IX,IY)**2)
                  UDBC=1./RBC*(UUBC-RRB(IX,IY)*UPBC)
                ELSE
                  RBC=0.
                  UDBC=0.
                ENDIF
cdr  now set cartesian flow velocity components
                VXINTF(IPLSV,ITRI)=(PUX(IN)*UUBC+PVX(IN)*VVBC)*V
                VYINTF(IPLSV,ITRI)=(PUY(IN)*UUBC+PVY(IN)*VVBC)*V
                VZINTF(IPLSV,ITRI)=(RBC*UPBC-RRB(IX,IY)*UDBC)*V
              ELSE
c  region outside  B2 grid, but inside triangular grid
                DIINTF(IPLS,ITRI)=DVAC
                VXINTF(IPLSV,ITRI)=VVAC
                VYINTF(IPLSV,ITRI)=VVAC
                VZINTF(IPLSV,ITRI)=VVAC
              ENDIF
            ENDDO  !itri loop
c  further additional cells from input block 2e ?
            DO itri=ntrii+2, nsbox
              DIINTF(IPLS,ITRI)=DVAC
              VXINTF(IPLSV,ITRI)=VVAC
              VYINTF(IPLSV,ITRI)=VVAC
              VZINTF(IPLSV,ITRI)=VVAC
            ENDDO

C  EIRENE BACKGROUND SPECIES "IPLS" IS NOW FILLED WITH B2 DATA "IFL"
 2201     CONTINUE   !IFL loop

C  NO DATA FOR "IPLS" IN B2 FILES
        ELSEIF (IFLB(IPLS).EQ.-13) THEN
C  READ DATA FOR "IPLS" FROM EIRENE DUMP FILE FT13

          IF (IREAD.EQ.0) THEN
            OPEN (UNIT=13,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
            REWIND 13
            READ (13,IOSTAT=IO) TEIN,TIIN,DEIN,DIIN,VXIN,VYIN,VZIN
            IREAD=1
            IF (TRCFLE) WRITE (iunout,*) 'READ 13: RCMUSR, IO= ',IO
            CLOSE (UNIT=13)
          ENDIF
          IF (IO.EQ.0) THEN
            IF (TRCINT) THEN
               WRITE(IUNOUT,*) 'PLASMA DATA FOR IPLS READ FROM FORT.13'
               WRITE(IUNOUT,*) 'IPLS, IPLSV,IPLSTI ',IPLS,IPLSV,IPLSTI
            ENDIF

            IPLSTI = MPLSTI(IPLS)
            IPLSV = MPLSV(IPLS)
            DIINTF(IPLS,:)=DIIN(IPLS,:)
            VXINTF(IPLSV,:)=VXIN(IPLSV,:)
            VYINTF(IPLSV,:)=VYIN(IPLSV,:)
            VZINTF(IPLSV,:)=VZIN(IPLSV,:)
            TIINTF(IPLSTI,:)=TIIN(IPLSTI,:)
          ENDIF
        ELSE
C  SET PARAMETERS FOR SPECIES IPLS TO ZERO
C  NOTHING TO BE DONE HERE
        ENDIF
 2200 CONTINUE
C  B2-BRAAMS CODE SPECIFIC END
C
C
C  READ OTHER B2_TRIA ARRAYS INTO EIRENE, FOR PRINTOUT AND PLOTTING
C
c  density, species index as in B2 code, cell-centered
      DO 2300 IAIN=1,NAINB
        IF (NAINT(IAIN).EQ.1.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO 2321 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=DNIB(IXTRI(IN),IYTRI(IN),NAINS(IAIN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
 2321     CONTINUE
c  poloidal (projection) flow velocity, species index as in B2 code, north surface-centered
        ELSEIF (NAINT(IAIN).EQ.2.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO 2322 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=UUB(IXTRI(IN),IYTRI(IN),NAINS(IAIN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
 2322     CONTINUE
c  radial drift velocity, species index as in B2 code, east surface-centered
        ELSEIF (NAINT(IAIN).EQ.3.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO 2323 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=VVB(IXTRI(IN),IYTRI(IN),NAINS(IAIN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
 2323     CONTINUE
c  plasma pressure, cell-centered, no species index
        ELSEIF (NAINT(IAIN).EQ.6) THEN
          DO 2326 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=PRB(IXTRI(IN),IYTRI(IN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
 2326     CONTINUE
c  parallel velocity, species index as in B2 code, north surface-centered
        ELSEIF (NAINT(IAIN).EQ.7.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO 2327 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=UPB(IXTRI(IN),IYTRI(IN),NAINS(IAIN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
 2327     CONTINUE
c  pitch angle, no species index
        ELSEIF (NAINT(IAIN).EQ.8) THEN
          DO 2328 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=RRB(IXTRI(IN),IYTRI(IN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
 2328     CONTINUE
        ELSEIF (NAINT(IAIN).EQ.9.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO 2329 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=FNIXB(IXTRI(IN),IYTRI(IN),NAINS(IAIN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
 2329     CONTINUE
        ELSEIF (NAINT(IAIN).EQ.10.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO 2330 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=FNIYB(IXTRI(IN),IYTRI(IN),NAINS(IAIN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
 2330     CONTINUE
        ELSEIF (NAINT(IAIN).EQ.11) THEN
          DO 2331 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=FEIXB(IXTRI(IN),IYTRI(IN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
 2331     CONTINUE
        ELSEIF (NAINT(IAIN).EQ.12) THEN
          DO 2332 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=FEIYB(IXTRI(IN),IYTRI(IN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
 2332     CONTINUE
        ELSEIF (NAINT(IAIN).EQ.13) THEN
          DO 2333 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=FEEXB(IXTRI(IN),IYTRI(IN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
 2333     CONTINUE
        ELSEIF (NAINT(IAIN).EQ.14) THEN
          DO 2334 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=FEEYB(IXTRI(IN),IYTRI(IN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
 2334     CONTINUE
C   NAINT=15,16:  USED ONLY IN B2.5 COUPLING: UUDIAG, VVDIAG
c   cell volume as in b2 code, no species index (cell-centered)
        ELSEIF (NAINT(IAIN).EQ.17) THEN
          DO 2335 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=VOLB(IXTRI(IN),IYTRI(IN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
 2335     CONTINUE
cdr  magnetic field strength, Tesla
        ELSEIF (NAINT(IAIN).EQ.18) THEN
          DO 2336 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=BFELDB(IXTRI(IN),IYTRI(IN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
 2336     CONTINUE


cdr  free: NAINT=20 --29:  reserved for AMDIAG:  scaled atomic/molecular rate coefficients
cdr                        evaluated on computational grid. See Manual.

        ENDIF
 2300 CONTINUE
C
      RETURN
C
 2999 CONTINUE
C
C  PLASMA PROFILES ARE NOW READ IN
C
      ENTRY EIRENE_IF2COP(ITRG)
      ITARG=ITRG
      IF (ITARG.GT.NTARGI) THEN
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) 'SOURCE DATA FOR STRATUM ISTRA= ',ITARG
        WRITE (iunout,*)
     .    'CANNOT BE DEFINED IN IF2COP. CHANGE INDSRC(ISTRA)'
        CALL EIRENE_LEER(1)
        RETURN
      ENDIF
C
C  NEXT DEFINE FLUXES, TEMPERATURES AND VELOCITIES AT THE TARGETS
C  (FLUXES IN AMP/(CM ALONG TARGET), TEMPERATURES IN EV, VELOCITIES IN CM/SEC)
C   FNIXB*FL (FNIYB*FL) ARE GIVEN IN AMP
C  STATEMENT NO 3000 ---> 3999
C
 3000 CONTINUE
C
      IF (TRCINT.AND.LTARG.EQ.0) THEN
        LTARG=1
        WRITE (iunout,*) 'ITARG: TARGET NUMBER '
        WRITE (iunout,*) 'IPRT : SUBSECTION OF TARGET '
        WRITE (iunout,*)
     .        'NPBS : BRAAMS (SURFACE) X-CELL INDEX OF TARGET '
        WRITE (iunout,*) 'NPBC : BRAAMS (ZONE) P-CELL INDEX OF TARGET '
        WRITE (iunout,*) 'NPES : POLOIDAL SURFACE INDEX OF TARGET'
        WRITE (iunout,*) '       IN EIRENE MESH'
        WRITE (iunout,*) 'NPEC : 1ST POLOIDAL CELL INDEX OF EIRENE MESH'
        WRITE (iunout,*) '       SEEN BY MONTE CARLO HISTORIES'
      ENDIF
C
      DO 3005 IPLS=1,NPLSI
      DO 3005 IGITT=1,NGITT
        FLSTEP(IPLS,ITARG,IGITT)=0.
 3005 CONTINUE
C
      ALLOCATE (NUMSID(NGITT))
      ALLOCATE (NUMTRI(NGITT))
      ALLOCATE (TORL(NSTRA,NGITT))
      ALLOCATE (ESHT(NSTEP,NGITT))
      ALLOCATE (ELTEST(NSTEP,NGITT))
      ALLOCATE (ORI(NSTEP,NGITT))

      RRSTEP(ITARG,1)=0.
      IG=0
      IIPLS=NPLSI
      IEPLS=1
      DO 3040 IPRT=1,NTGPRT(ITARG)
C  NINCT= 1: PLASMA FLUX IN SAME   DIRECTION AS B2 COORDINATE
C  NINCT=-1: PLASMA FLUX IN OPPOS. DIRECTION AS B2 COORDINATE
C  BRAAMS X-CELL CONTAINING THE TARGET DATA (BOUNDARY CONDITIONS)
C  (SURFACE-CENTERED, EAST OR NORTH) (AFTER INDEX MAPPING)
C  (E.G. SURFACE NO.0 AND SURFACE NO. NX) AT TARGETS.
        NPBS=NDT(ITARG,IPRT)
C  BRAAMS P-CELL CONTAINING THE TARGET DATA (BOUNDARY CONDITIONS)
C  (ZONE-CENTERED) (AFTER INDEX MAPPING)
C
C  THIS LINE, IF B2-BOUNDARY CONDITIONS ARE COMPUTED FROM GUARD CELLS
C  (E.G. CELL NO.0 AND CELL NO. NX+1) AT TARGETS.
        NPBC=NPBS+MAX0(0,NINCT(ITARG,IPRT))
C
C  1ST EIRENE CELL ALONG TARGET
        NPEC=NPBS-MIN0(0,NINCT(ITARG,IPRT))
C  EIRENE SURFACE NUMBER AT TARGET
        NPES=NPBS+1
        IF (TRCINT) THEN
          WRITE (iunout,'(a,6(i4))') 'ITARG,IPRT,NPBS,NPBC,NPES,NPEC ',
     .                                ITARG,IPRT,NPBS,NPBC,NPES,NPEC
        ENDIF
C
C  FIRST: SOURCES AT POLOIDAL (Y) SURFACES (EAST OR WEST CELL FACES)
        IF (NIXY(ITARG,IPRT).EQ.2) GOTO 3020
C
        ITRI=0
        DO IY=NTIN(ITARG,IPRT),NTEN(ITARG,IPRT)-1
          ICOU = 0
          CURPOI => HEADS(IY,NPEC)%P
          DO WHILE (ASSOCIATED(CURPOI))
            IT=CURPOI%TRIANGLE
C  TEST WHETHER TRIANGLE BELONGS TO QUADRANGULAR CELLS ALONG THE TARGET
            IF (IXTRI(IT).EQ.NPEC .AND.
     .         (IYTRI(IT).GE.NTIN(ITARG,IPRT) .AND.
     .          IYTRI(IT).LT.NTEN(ITARG,IPRT))) THEN
              dxpol = xpol(iy+1,npes) - xpol(iy,npes)
              dypol = ypol(iy+1,npes) - ypol(iy,npes)
              dd = sqrt(dxpol*dxpol + dypol*dypol)

              do is = 1, 3
!  test if side IS of triangle IT is parallel to B2 cell face
                par = vtrix(is,it)*dypol-vtriy(is,it)*dxpol
cdr  scale to relative units of triangle coordinates
                par=par/dd

                if (abs(par) < 5*eps5) then
                  is1 = is + 1
                  if (is1 > 3) is1 = 1
!  test if the vertices of the triangle side lie on B2 cell face
                  L1 = EIRENE_POINT_ON_INTERVAL(XTRIAN(NECKE(IS,IT)),
     f                             YTRIAN(NECKE(IS,IT)),
     f                             XPOL(IY,NPES),YPOL(IY,NPES),
     f                             XPOL(IY+1,NPES),YPOL(IY+1,NPES))
                  L2 = EIRENE_POINT_ON_INTERVAL(XTRIAN(NECKE(IS1,IT)),
     f                             YTRIAN(NECKE(IS1,IT)),
     f                             XPOL(IY,NPES),YPOL(IY,NPES),
     f                             XPOL(IY+1,NPES),YPOL(IY+1,NPES))
                  IF ( L1 .AND. L2 ) THEN
                    ITRI=ITRI+1
                    IF (ITRI.GT.NGITT) THEN
                      WRITE (iunout,*)
     .                        ' NOT ENOUGH GRIDPOINTS FOR DEFINING',
     .                        ' STEP-FUNCTION '
                      WRITE (iunout,*) ' INCREASE PARAMETER NGITT '
                      WRITE (iunout,*) ' NGITT = ',NGITT
                      CALL EIRENE_EXIT_OWN(1)
                    ENDIF
!  triangle found
                    NUMTRI(ITRI)=IT
                    NUMSID(ITRI)=IS
                    IF (LCHKQUD)
     .              IREVERS(IS,IT) = INT(SIGN(1._DP,DXPOL*VTRIX(IS,IT) +
     .                                              DYPOL*VTRIY(IS,IT)))
                    ICOU = ICOU + 1
                    EXIT
                  END IF
                end if
              end do
            ENDIF
            CURPOI => CURPOI%NEXT
          ENDDO
!  there has to be at least one triangle per B2 cell
          IF (ICOU == 0) THEN
            WRITE (iunout,*) ' NO TRIANGLE FOUND FOR B2 CELL ',NPBC, IY
            CALL EIRENE_EXIT_OWN(1)
          END IF
        ENDDO
        MTRI=ITRI
C SORT TRIANGLES ALONG TARGET
         XANF=XPOL(NTIN(ITARG,IPRT),NPES)
         YANF=YPOL(NTIN(ITARG,IPRT),NPES)
         IANF=1
         IACT=1
         DO WHILE (IANF .LT. MTRI)
           DO IT=IANF,MTRI
             ITRI=NUMTRI(IT)
             IS=NUMSID(IT)
             IS1=IS+1
             IF (IS1.GT.3) IS1=1
             IF (((XANF-XTRIAN(NECKE(IS,ITRI)))**2+
     .           (YANF-YTRIAN(NECKE(IS,ITRI)))**2) .LT. 5*EPS5) THEN
              NUMTRI(IT)=NUMTRI(IACT)
              NUMSID(IT)=NUMSID(IACT)
              NUMTRI(IACT)=ITRI
              NUMSID(IACT)=IS
              IACT=IACT+1
              XANF=XTRIAN(NECKE(IS1,ITRI))
              YANF=YTRIAN(NECKE(IS1,ITRI))
            ELSEIF (((XANF-XTRIAN(NECKE(IS1,ITRI)))**2+
     .               (YANF-YTRIAN(NECKE(IS1,ITRI)))**2).LT.5*EPS5) THEN
              NUMTRI(IT)=NUMTRI(IACT)
              NUMSID(IT)=NUMSID(IACT)
              NUMTRI(IACT)=ITRI
              NUMSID(IACT)=IS
              IACT=IACT+1
              XANF=XTRIAN(NECKE(IS,ITRI))
              YANF=YTRIAN(NECKE(IS,ITRI))
            ENDIF
          ENDDO
          IANF=IACT
        ENDDO
C
        DO IT=1,MTRI
          ITRI=NUMTRI(IT)
          IS=NUMSID(IT)
          IS1=IS+1
          IF (IS1.GT.3) IS1=1
          IX=IXTRI(ITRI)
          IY=IYTRI(ITRI)
          IF ((IX < 0) .OR. (IY < 0)) CYCLE
          IG=IG+1
          IF (IG.GT.NGITT) GOTO 999
C  TESTEP, TISTEP: ZONE-CENTERED TEMPERATURE IN BOUNDARY ZONE (EV)
          ORI(ITARG,IG) = NINCT(ITARG,IPRT)
          TESTEP(ITARG,IG) = TEB(NPBC,IY)*T
C  RRSTEP,IRSTEP,IPSTEP: GEOMETRICAL INFORMATION ALONG TARGET
C  RRSTEP IS THE ARC LENGTH ALONG THE TARGET (CM)
          RRSTEP(ITARG,IG+1)=RRSTEP(ITARG,IG) +
     .        SQRT((XTRIAN(NECKE(IS,ITRI))-XTRIAN(NECKE(IS1,ITRI)))**2+
     .             (YTRIAN(NECKE(IS,ITRI))-YTRIAN(NECKE(IS1,ITRI)))**2)
C  EIRENE CELL NUMBER INFORMATION ALONG TARGET
          IRSTEP(ITARG,IG)=ITRI
          IPSTEP(ITARG,IG)=IS
          ITSTEP(ITARG,IG)=1
          IASTEP(ITARG,IG)=0
          IBSTEP(ITARG,IG)=1
          IGSTEP(ITARG,IG)=200000+NPES
          IF (INMTI(IS,ITRI).EQ.0) THEN
            WRITE (iunout,*) 'ERROR IN INFCOP '
            WRITE (iunout,*) 'SOURCE NOT ON A KNOWN SURFACE'
            WRITE (iunout,*) 'ITARG,IG,IPRT ',ITARG,IG,IPRT
          ENDIF
C FACTOR FOR SHEATH POTENTIAL: Vsh=Te*SHSTEP(ITARG,IG)
          SHSTEP(ITARG,IG)=DELTA_SHEATHXB(NPBS,IY)
C  TORL: TOROIDAL LENGTH (CM) AT TARGET SEGMENT IY: CENTER OF GRAVITY
          TORL(ITARG,IG)=2.*PIA*0.5*(XTRIAN(NECKE(IS,ITRI))+
     .                               XTRIAN(NECKE(IS1,ITRI)))
          DO 3013 IPLS=1,NPLSI
            IPLSTI=MPLSTI(IPLS)
            IPLSV=MPLSV(IPLS)
            ELSTEP(IPLS,ITARG,IG)=0.
            TISTEP(IPLSTI,ITARG,IG) = TIB(NPBC,IY)*T
C  DISTEP: ZONE-CENTERED DENSITY IN BOUNDARY ZONE
            IFL=IFLB(IPLS)
            IF (IFL.LE.0.OR.IFL.GT.NFLA) GOTO 3013
            DISTEP(IPLS,ITARG,IG)=DNIB(NPBC,IY,IFL)*D(IPLS)
C  FLSTEP: SURFACE-CENTERED FLUX (AMP/CM ALONG TARGET)
            IF (NSPZI(ITARG,IPRT).LE.IFL.AND.
     .                               IFL.LE.NSPZE(ITARG,IPRT)) THEN
              IIPLS=MIN0(IIPLS,IPLS)
              IEPLS=MAX0(IEPLS,IPLS)
!  LENGTH OF CELL FACE OF B2 CELL
              CELDEL=SQRT((XPOL(IY+1,NPES)-XPOL(IY,NPES))**2+
     .                    (YPOL(IY+1,NPES)-YPOL(IY,NPES))**2)
!  LENGTH OF TRIANGLE SIDE ALONG B2 CELL FACE
              DELY =
     .        SQRT((XTRIAN(NECKE(IS,ITRI))-XTRIAN(NECKE(IS1,ITRI)))**2+
     .             (YTRIAN(NECKE(IS,ITRI))-YTRIAN(NECKE(IS1,ITRI)))**2)
              FLSTEP(IPLS,ITARG,IG)=0.
              IF (DELY.GT.0.) THEN
                FRAC = DELY / CELDEL
                FLSTEP(IPLS,ITARG,IG)=MAX(0._DP,ORI(ITARG,IG)*
     .                          FNIXB(NPBS,IY,IFL))*FL(IPLS)/DELY*FRAC

C  SET DEFAULT ION ENERGY FLUXES FROM B2_TRIA BOUNDARY CONDITIONS
                delti_para=3
                delte_para=0.5
                delti_perp=2
                delte_perp=0
!  only one of the next two is different from 0
!pb                delti_para=deltai_parxb(npbs,iy)
!pb                delti_perp=deltai_radxb(npbs,iy)
!  only one of the next two is different from 0
!pb                delte_para=deltae_parxb(npbs,iy)
!pb                delte_perp=deltae_radxb(npbs,iy)

                   tis=TISTEP(IPLSTI,ITARG,IG)
                   tes=TESTEP(ITARG,IG)
!pb                ELSTEP(IPLS,ITARG,IG) = ELSTEP(IPLS,ITARG,IG)+
!pb     .                                  FL(IPLS)/DELY*
!pb     .            (TIS*delti_perp*ABS(Fnix_yb(npbs,iy,ifl))+
!pb     .             TIS*delti_para*ABS(fnixb  (npbs,iy,ifl))+
!pb     .             TES*delte_para*ABS(fnixb  (npbs,iy,ifl)))

                   ELSTEP(IPLS,ITARG,IG) = ELSTEP(IPLS,ITARG,IG)+
     .                                     FL(IPLS)/DELY*
     .            (TIS*delti_para+TES*delte_para)
     .             *ABS(fnixb(npbs,iy,ifl))
              ENDIF
            ENDIF

C  VXSTEP,VYSTEP,VZSTEP: SURFACE-CENTERED FLOW VELOCITY (CM/S)
C  NOTE: PV VECTOR IS CELL-CENTERED, BUT EXACT VECTOR CAN BE FOUND FROM
C        DATA FOR POLOIDAL POLYGON NPES
            IN=IY+(NPEC-1)*NR1TAL_SAVE

            PVXS=PVXE(IN)
            PVYS=PVYE(IN)
            PUXS=PUXE(IN)
            PUYS=PUYE(IN)
            PITB=RRB(NPBS,IY)
            UU=PITB*VPARXB(NPBS,IY,IFL)
            VXSTEP(IPLSV,ITARG,IG)=
     .            (PUXS*UU+PVXS*VRADXB(NPBS,IY,IFL))*V
            VYSTEP(IPLSV,ITARG,IG)=
     .            (PUYS*UU+PVYS*VRADXB(NPBS,IY,IFL))*V
            VZSTEP(IPLSV,ITARG,IG)=
     .            (SQRT(1.-PITB**2)*VPARXB(NPBS,IY,IFL))*V
            RRBS=UUB(NPBS,IY,IFL)/(UPB(NPBS,IY,IFL)+EPS60)
            IF(ABS(RRBS).GT.1) THEN                            !VK
              WRITE(IUNOUT,*) "IF2COP WARNING,UUB>UPB!"        !VK
              VZSTEP(IPLS,ITARG,IG)=0                          !VK
            ELSE                                               !VK
              VZSTEP(IPLSV,ITARG,IG)=
     .              (SQRT(1.-RRBS**2)*UPB(NPBS,IY,IFL))*V
            END IF                                             !VK
 3013     CONTINUE

          IF (.FALSE.) THEN
!  intermediate model for testing:
!  within a given cell all incident ions (of all species) have same
!  energy E0B2 (delta function approximation) (NEMOD1=8 or NEMOD1=9)

          DELY=RRSTEP(ITARG,IG+1)-RRSTEP(ITARG,IG)
          IF (DELY.GT.0.) THEN
            fltt = sum(flstep(1:nplsi,itarg,ig))
            e0b2 = ABS(feixb(npbs,iy))/dely/(fltt+eps60)*frac
            do ipls=1,nplsi
              IFL=IFLB(IPLS)
              IF (IFL.LE.0.OR.IFL.GT.NFLA) CYCLE
              IF (NSPZI(ITARG,IPRT).LE.IFL.AND.
     .                               IFL.LE.NSPZE(ITARG,IPRT)) THEN
                ELSTEP(IPLS,ITARG,IG) = ELSTEP(IPLS,ITARG,IG)+
     .                                  FLSTEP(IPLS,ITARG,IG)*E0B2
              END IF
            end do
          END IF
          END IF
        ENDDO
C
        GOTO 3030
C
 3020   CONTINUE
C
C  SECOND: SOURCES AT RADIAL (X) SURFACES
C
        ITRI=0
        DO IX=NTIN(ITARG,IPRT),NTEN(ITARG,IPRT)-1
          IF (LLCUT(IX)) CYCLE
          ICOU = 0
          CURPOI => HEADS(NPEC,IX)%P
          DO WHILE (ASSOCIATED(CURPOI))
            IT=CURPOI%TRIANGLE
C  TEST WHETHER TRIANGLE BELONGS TO QUADRANGULAR CELLS ALONG THE TARGET
            IF (IYTRI(IT).EQ.NPEC .AND.
     .         (IXTRI(IT).GE.NTIN(ITARG,IPRT) .AND.
     .          IXTRI(IT).LT.NTEN(ITARG,IPRT))) THEN
              ISC=0
              dxpol = xpol(npes,ix+1) - xpol(npes,ix)
              dypol = ypol(npes,ix+1) - ypol(npes,ix)
              dd = sqrt(dxpol*dxpol + dypol*dypol)
c
              do is = 1, 3
!  test if side IS of triangle IT is parallel to B2 cell face
                par = vtrix(is,it)*dypol-vtriy(is,it)*dxpol
cdr  scale to relative units of triangle coordinates
                par=par/dd


                if (abs(par) < 5.E-3_dp) then
                  is1 = is + 1
                  if (is1 > 3) is1 = 1
!  test if the vertices of the triangle side lie on B2 cell face
                  L1 = EIRENE_POINT_ON_INTERVAL(XTRIAN(NECKE(IS,IT)),
     f                             YTRIAN(NECKE(IS,IT)),
     f                             XPOL(NPES,IX),YPOL(NPES,IX),
     f                             XPOL(NPES,IX+1),YPOL(NPES,IX+1))
                  L2 = EIRENE_POINT_ON_INTERVAL(XTRIAN(NECKE(IS1,IT)),
     f                             YTRIAN(NECKE(IS1,IT)),
     f                             XPOL(NPES,IX),YPOL(NPES,IX),
     f                             XPOL(NPES,IX+1),YPOL(NPES,IX+1))
                  IF ( L1 .AND. L2 ) THEN
                    ITRI=ITRI+1
                    IF (ITRI.GT.NGITT) THEN
                      WRITE (iunout,*)
     .                        ' NOT ENOUGH GRIDPOINTS FOR DEFINING',
     .                        ' STEP-FUNCTION '
                      WRITE (iunout,*) ' INCREASE PARAMETER NGITT '
                      WRITE (iunout,*) ' NGITT = ',NGITT
                      CALL EIRENE_EXIT_OWN(1)
                    ENDIF
!  triangle found
                    NUMTRI(ITRI)=IT
                    NUMSID(ITRI)=IS
                    ICOU = ICOU + 1
                    IF (LCHKQUD)
     .              IREVERS(IS,IT) = INT(SIGN(1._DP,DXPOL*VTRIX(IS,IT) +
     .                                              DYPOL*VTRIY(IS,IT)))
                    EXIT
                  END IF
                end if
              end do
            ENDIF
            CURPOI => CURPOI%NEXT
          ENDDO
!  there has to be at least one triangle per B2 cell
          IF (ICOU == 0) THEN
            WRITE (iunout,*) ' NO TRIANGLE FOUND FOR B2 CELL ',IX, NPBC
            CALL EIRENE_EXIT_OWN(1)
          END IF
        ENDDO
        MTRI=ITRI
C SORT TRIANGLES ALONG TARGET
        XANF=XPOL(NPES,NTIN(ITARG,IPRT))
        YANF=YPOL(NPES,NTIN(ITARG,IPRT))
        IANF=1
        IACT=1
        DO WHILE (IANF .LT. MTRI)
          DO IT=IANF,MTRI
            ITRI=NUMTRI(IT)
            IS=NUMSID(IT)
            IS1=IS+1
            IF (IS1.GT.3) IS1=1
            IF (((XANF-XTRIAN(NECKE(IS,ITRI)))**2+
     .           (YANF-YTRIAN(NECKE(IS,ITRI)))**2) .LT. 5*EPS5) THEN
              NUMTRI(IT)=NUMTRI(IACT)
              NUMSID(IT)=NUMSID(IACT)
              NUMTRI(IACT)=ITRI
              NUMSID(IACT)=IS
              IACT=IACT+1
              XANF=XTRIAN(NECKE(IS1,ITRI))
              YANF=YTRIAN(NECKE(IS1,ITRI))
            ELSEIF (((XANF-XTRIAN(NECKE(IS1,ITRI)))**2+
     .               (YANF-YTRIAN(NECKE(IS1,ITRI)))**2).LT.5*EPS5) THEN
              NUMTRI(IT)=NUMTRI(IACT)
              NUMSID(IT)=NUMSID(IACT)
              NUMTRI(IACT)=ITRI
              NUMSID(IACT)=IS
              IACT=IACT+1
              XANF=XTRIAN(NECKE(IS,ITRI))
              YANF=YTRIAN(NECKE(IS,ITRI))
            ENDIF
          ENDDO
          IANF=IACT
        ENDDO
C
        DO IT=1,MTRI
          ITRI=NUMTRI(IT)
          IS=NUMSID(IT)
          IS1=IS+1
          IF (IS1.GT.3) IS1=1
          IX=IXTRI(ITRI)
          IY=IYTRI(ITRI)
          IF ( (IX < 0) .OR. (IY < 0) ) CYCLE
          IG=IG+1
          IF (IG.GT.NGITT) GOTO 999
C  TESTEP, TISTEP: ZONE-CENTERED TEMPERATURE IN BOUNDARY ZONE (EV)
          ORI(ITARG,IG) = NINCT(ITARG,IPRT)
          TESTEP(ITARG,IG) = TEB(IX,NPBC)*T
C  RRSTEP,IRSTEP,IPSTEP: GEOMETRICAL INFORMATION ALONG TARGET
C  EIRENE CELL NUMBER INFORMATION ALONG TARGET
          RRSTEP(ITARG,IG+1)=RRSTEP(ITARG,IG) +
     .       SQRT((XTRIAN(NECKE(IS,ITRI))-XTRIAN(NECKE(IS1,ITRI)))**2+
     .            (YTRIAN(NECKE(IS,ITRI))-YTRIAN(NECKE(IS1,ITRI)))**2)
          IRSTEP(ITARG,IG)=ITRI
          IPSTEP(ITARG,IG)=IS
          ITSTEP(ITARG,IG)=1
          IASTEP(ITARG,IG)=0
          IBSTEP(ITARG,IG)=1
          IGSTEP(ITARG,IG)=100000+NPES
          IF (INMTI(IS,ITRI).EQ.0) THEN
            WRITE (iunout,*) 'ERROR IN INFCOP '
            WRITE (iunout,*) 'SOURCE NOT ON A KNOWN SURFACE'
            WRITE (iunout,*) 'ITARG,IG,IPRT ',ITARG,IG,IPRT
          ENDIF
C FACTOR FOR SHEATH POTENTIAL: Vsh=Te*SHSTEP(ITARG,IG)
          SHSTEP(ITARG,IG)=DELTA_SHEATHYB(IX,NPBS)
C  TORL: TOROIDAL LENGTH (CM) AT TARGET SEGMENT IY: CENTER OF GRAVITY
          TORL(ITARG,IG)=2.*PIA*0.5*(XTRIAN(NECKE(IS,ITRI))+
     .                               XTRIAN(NECKE(IS1,ITRI)))
          DO 3023 IPLS=1,NPLSI
            IPLSTI=MPLSTI(IPLS)
            IPLSV=MPLSV(IPLS)
            ELSTEP(IPLS,ITARG,IG)=0.
            TISTEP(IPLSTI,ITARG,IG) = TIB(IX,NPBC)*T
C  DISTEP: ZONE-CENTERED DENSITY IN BOUNDARY ZONE (EV)
            IFL=IFLB(IPLS)
            IF (IFL.LE.0.OR.IFL.GT.NFLA) GOTO 3023
            DISTEP(IPLS,ITARG,IG)=DNIB(IX,NPBC,IFL)*D(IPLS)
C  FLSTEP: SURFACE-CENTERED FLUX (AMP/CM ALONG TARGET)
            IF (NSPZI(ITARG,IPRT).LE.IFL.AND.
     .                               IFL.LE.NSPZE(ITARG,IPRT)) THEN
              IIPLS=MIN0(IIPLS,IPLS)
              IEPLS=MAX0(IEPLS,IPLS)
!  LENGTH OF CELL FACE OF B2 CELL
              CELDEL=SQRT((XPOL(NPES,IX+1)-XPOL(NPES,IX))**2+
     .                  (YPOL(NPES,IX+1)-YPOL(NPES,IX))**2)
!  LENGTH OF TRIANGLE SIDE ALONG B2 CELL FACE
              DELX =
     .        SQRT((XTRIAN(NECKE(IS,ITRI))-XTRIAN(NECKE(IS1,ITRI)))**2+
     .             (YTRIAN(NECKE(IS,ITRI))-YTRIAN(NECKE(IS1,ITRI)))**2)
              FLSTEP(IPLS,ITARG,IG)=0.
              IF (DELX.GT.0.) THEN
cdr why should celdel and dely be different ??
                FRAC = DELX / CELDEL
                FLSTEP(IPLS,ITARG,IG)=MAX(0._DP,ORI(ITARG,IG)*
     .                            FNIYB(IX,NPBS,IFL))*FL(IPLS)/DELX*FRAC

C  SET DEFAULT ION ENERGY FLUXES FROM B2 BOUNDARY CONDITIONS
                delti_para=3
                delte_para=0.5
                delti_perp=2
                delte_perp=0
!  only one of the next two is different from 0
!pb                delti_para=deltai_paryb(ix,npbs)
!pb                delti_perp=deltai_radyb(ix,npbs)
!  only one of the next two is different from 0
!pb                delte_para=deltae_paryb(ix,npbs)
!pb                delte_perp=deltae_radyb(ix,npbs)
                tis=TISTEP(IPLSTI,ITARG,IG)
                tes=TESTEP(ITARG,IG)
!pb                ELSTEP(IPLS,ITARG,IG) = ELSTEP(IPLS,ITARG,IG) +
!pb     .                                  FL(IPLS)/DELX*
!pb     .            (TIS*delti_perp*ABS(Fniyb  (ix,npbs,ifl))+
!pb     .             TIS*delti_para*ABS(fniy_xb(ix,npbs,ifl))+
!pb     .             TES*delte_para*ABS(fniy_xb(ix,npbs,ifl)))

                ELSTEP(IPLS,ITARG,IG) = ELSTEP(IPLS,ITARG,IG) +
     .                                  FL(IPLS)/DELX*
     .            (TIS*delti_perp*ABS(Fniyb(ix,npbs,ifl)))


!dr  try ion energy fluxes from B2 directly. But then velocs-sampling inconsistency
!                ELSTEP(IPLS,ITARG,IG) = ELSTEP(IPLS,ITARG,IG) +
!     .                                  ABS(Feiyb(ix,npbs))/DELX
              ENDIF
            ENDIF
C
C  VXSTEP,VYSTEP,VZSTEP: SURFACE-CENTERED FLOW VELOCITY (CM/S)
C  NOTE: PU VECTOR IS CELL-CENTERED, BUT EXACT VECTOR CAN BE FOUND FROM
C        RADIAL POLYGON NPES DATA
            IN=NPEC+(IX-1)*NR1TAL_SAVE
            PVXS=PVXN(IN)
            PVYS=PVYN(IN)
            PUXS=PUXN(IN)
            PUYS=PUYN(IN)
!pb
            IXM1 = IX-1
            IF (LLCUT(IXM1)) THEN
               IXM1 = NGHPOL(4,NPEC,IX)
            END IF
!pb
            UUBC=0.5*(UUB(IXM1,NPBC,IFL)+UUB(IX,NPBC,IFL))
            UPBC=0.5*(UPB(IXM1,NPBC,IFL)+UPB(IX,NPBC,IFL))
            PITB=RRB(NPBS,IY)
            UU=PITB*VPARYB(IX,NPBS,IFL)
            VXSTEP(IPLSV,ITARG,IG)=
     .            (PUXS*UU+PVXS*VRADYB(IX,NPBS,IFL))*V
            VYSTEP(IPLSV,ITARG,IG)=
     .            (PUYS*UU+PVYS*VRADYB(IX,NPBS,IFL))*V
            VZSTEP(IPLSV,ITARG,IG)=
     .            (SQRT(1.-PITB**2)*VPARYB(IX,NPBS,IFL))*V

!pb         RRBS=UUB(IX,NPBS,IFL)/(UPB(IX,NPBS,IFL)+EPS60)
            RRBS=RRB(IX,NPBC)
            IF(ABS(RRBS).GT.1) THEN                            !VK
              WRITE(IUNOUT,*) "IF2COP WARNING,UUB>UPB!"        !VK
              VZSTEP(IPLS,ITARG,IG)=0                          !VK
            ELSE                                               !VK
!pb 17.09.2012
             VZSTEP(IPLS,ITARG,IG)=
     .             (SQRT(1.-RRBS**2)*UPBC)*V
!pb              VZSTEP(IPLSV,ITARG,IG)=
!pb     .              (SQRT(1.-RRBS**2)*UPB(IX,NPBS,IFL))*V
            END IF                                             !VK
 3023     CONTINUE

          IF (.FALSE.) THEN
!  intermediate model for testing:
!  within a given cell all incident ions (of all species) have same
!  energy E0B2 (delta function approximation) (NEMOD1=8 or NEMOD1=9)

          DELX=RRSTEP(ITARG,IG+1)-RRSTEP(ITARG,IG)
          IF (DELX.GT.0.) THEN
            fltt = sum(flstep(1:nplsi,itarg,ig))
            e0b2 = ABS(feiyb(ix,npbs))/delx/(fltt+eps60)*frac
            do ipls=1,nplsi
              IFL=IFLB(IPLS)
              IF (IFL.LE.0.OR.IFL.GT.NFLA) CYCLE
              IF (NSPZI(ITARG,IPRT).LE.IFL.AND.
     .                               IFL.LE.NSPZE(ITARG,IPRT)) THEN
                ELSTEP(IPLS,ITARG,IG) = ELSTEP(IPLS,ITARG,IG)+
     .                                  FLSTEP(IPLS,ITARG,IG)*E0B2
              END IF
            end do
          END IF
          END IF
        ENDDO
 3030   CONTINUE
C
 3040 CONTINUE
      NRWL(ITARG)=IG+1

C
      IF (TRCINT) CALL EIRENE_LEER(2)
C
C  INITIALIZE FUNCTION STEP (FOR RANDOM SAMPLING ALONG TARGET)
C  SET SOME SOURCE PARAMETERS EXPLICITLY TO ENFORCE INPUT CONSISTENCY
C  also: sum over species: flstep(0,...), elstep(0,...) will be set.
C
      FLUX(ITARG)=EIRENE_STEP(IIPLS,IEPLS,NRWL(ITARG),ITARG)
C
      NLPLS(ITARG)=.TRUE.
      NLATM(ITARG)=.FALSE.
      NLMOL(ITARG)=.FALSE.
      NLION(ITARG)=.FALSE.
C
      NLSRF(ITARG)=.TRUE.
      NLPNT(ITARG)=.FALSE.
      NLLNE(ITARG)=.FALSE.
      NLVOL(ITARG)=.FALSE.
      NLCNS(ITARG)=.FALSE.
C
      NSRFSI(ITARG)=1
CTRIG A
C  IN TRIA OPTION INDIM=4 (for LEVGEO=3) CORRESPONDS TO INDIM=1 (for LEVGEO=4)
      INDIM(1,ITARG)=1
CTRIG E
      IF (INDSRC(ITARG).NE.6) THEN
        I34=EIRENE_IDEZ(INT(SORLIM(1,ITARG)),3,3)
        SORLIM(1,ITARG)=I34*100+40   !sample with step fct.
      ELSEIF (INDSRC(ITARG).EQ.6) THEN
C  SORLIM DEFAULT WAS 0.D0
        SORLIM(1,ITARG)=0240
      ENDIF
      SORIND(1,ITARG)=ITARG

!  intermediate model for testing:
      IF (.FALSE.) THEN
!  SELECT THE PREPROGRAMMED SOURCE ENERGY CONDITIONAL DISTRIBUTION CONSISTENT
!  TO B2 SETTINGS:  USE ELSTEP DIRECTLY, NO SAMPLING
C  USE ENERGY FLUXES ELSTEP SPECIFIED ABOVE, IE., SORENE, SORENI ARE REDUNDANT
      IF (ABS(SUM(SHSTEP(ITARG,1:NRWL(ITARG)))) > EPS10) THEN
!  USE SHEATH PARAMETER GIVEN BY B2
        NEMODS(ITARG)=9
      ELSE
!  NO SHEATH IS CALCULATED
        NEMODS(ITARG)=8
      END IF
      ENDIF

C  IN CASE INDIM=4: INSOR,INDGRD,... ARE REDUNDANT
      NRSOR(1,ITARG)=-1
      NPSOR(1,ITARG)=-1
      IF (INDSRC(ITARG).LT.6) THEN
        WRITE (iunout,*) 'MESSAGE FROM IF2COP: '
        WRITE (iunout,*) 'SOURCE STRENGTH AND SPATIAL DISTRIBUTION FOR '
        WRITE (iunout,*) 'STRATUM ',ISTRA,' MODIFIED.'
        CALL EIRENE_MASR1('FLUX=   ',FLUX(ISTRA))
        WRITE (iunout,*) 'USE STEP FUNCTION ISTEP= ',ITARG,
     .                   ' FROM BLOCK 14'
        WRITE (iunout,*) ' FROM BLOCK 14: SORLIM= ',SORLIM(1,ITARG)
        CALL EIRENE_LEER(1)
      ENDIF
C
      IF (INDSRC(ITARG).EQ.6) THEN
C  DEFINE SOURCE FOR TARGET RECYCLING STRATUM ITARG
C  ASSUME NOW: ITARG=ISTRA
C  DEFAULTS ARE ALREADY SET IN SUBR. INPUT.
C
        CALL EIRENE_FTCRI(ITARG,CITARG)
        TXTSOU(ITARG)= 'SURFACE RECYCLING SOURCE NO.'//CITARG
        NPTS(ITARG)=NPTC(ITARG,1)*MPTS_COMSOU
        IF (NPTS(ITARG) < 0) NPTS(ITARG) = HUGE(1)
        NMINPTS(ITARG)=NPTCM(ITARG,1)*MPTS_COMSOU !VK MINIMUM NUMBER OF HYSTORIES FOR THE STRATUM
        NINITL(ITARG)=ITARG*1001
        NSPEZ(ITARG)=0
        SORIFL(1,ITARG)=NIFLG(ITARG,1)
        SORWGT(1,ITARG)=1.
C  USE ENERGY FLUXES SPECIFIED HERE, IE., SORENE, SORENI ARE REDUNDANT
        NEMODS(ITARG)=9
        NAMODS(ITARG)=1
C
C  USE POLYGON MESH, IE., SORAD1,...,SORAD4 ARE REDUNDANT.
        SORAD5(1,ITARG)=ZIA
        SORAD6(1,ITARG)=ZAA
C
C  VELOCITY SPACE DISTRIBUTION
        SORCOS(ITARG)=1.
        SORMAX(ITARG)=0.
C
C
C  DO 2028 LOOP FROM SUBR. INPUT
        THMAX=MAX(0._DP,MIN(PIHA,SORMAX(ITARG)*DEGRAD))
        IF (NAMODS(ITARG).EQ.1) THEN
          RP1=SORCOS(ITARG)+1.
          SORCOS(ITARG)=1./RP1
          IF (ABS(COS(THMAX)).LE.EPS10) THEN
            SORMAX(ITARG)=1.
          ELSE
            SORMAX(ITARG)=1.-COS(THMAX)**RP1
          ENDIF
        ELSEIF (NAMODS(ITARG).EQ.2) THEN
          SORCOS(ITARG)=SORCOS(ITARG)*DEGRAD
          SORMAX(ITARG)=THMAX
        ENDIF
        NLSYMT(0)=NLSYMT(0).AND.NLSYMT(ITARG)
        NLSYMP(0)=NLSYMP(0).AND.NLSYMP(ITARG)
C
      ENDIF
C
C  SOURCE DEFINITION FOR TARGET RECYCLING STRATUM ITARG COMPLETED
C
 3999 CONTINUE
C
C  TARGET DATA ITARG ARE DEFINED NOW
C
C
C  COMPUTE MACH NUMBERS, SHEATH POTENTIAL, AND EXACT SURFACE ENERGY FLUXES
C  FOR COMPARISON WITH SAMPLED ENERGY FLUXES
C  E-FLUX "ETOTP". THIS IS ONLY FOR DIAGNOSTICS PURPOSES
C  E.G. TO CHECK CONSISTENCY OF BOUNDARY CONDITIONS
C  STATEMENT NO. 6000 ---> 6500
C
      IF (.NOT.TRCINT) GOTO 6500
C
      EEMAX=0.
      EESHT=0.
C
      NEM=NEMODS(ITARG)
      DO 6011 IG=1,NRWL(ITARG)-1
        OR=ORI(ITARG,IG)
C
C  COMPUTE SHEATH POTENTIAL ESHT(ITARG,IG)
C  USE ALL NPLSI SPECIES, NOT JUST IFL=NSPZI,NSPZE
C
        ESHT(ITARG,IG)=0.D0
        IF (NEM.EQ.3.OR.NEM.EQ.5.OR.NEM.EQ.7) THEN
          DO 6005 IPL=1,NPLSI
            IPLV=MPLSV(IPL)
            VPX=VXSTEP(IPLV,ITARG,IG)
            VPY=VYSTEP(IPLV,ITARG,IG)
            VPZ=VZSTEP(IPLV,ITARG,IG)
            VP(IPL)=SQRT(VPX**2+VPY**2+VPZ**2)
            DI(IPL)=DISTEP(IPL,ITARG,IG)
 6005     CONTINUE
          TE=TESTEP(ITARG,IG)
          CUR=0.
          GAMMA=0.
          ESHT(ITARG,IG)=EIRENE_SHEATH(TE,DI,VP,NCHRGP,GAMMA,CUR,
     .                         NPLSI,-ITARG)
        ELSE IF (NEM == 9) THEN
          IF (IGSTEP(ITARG,IG).GT.200000) THEN
            IY=IYTRI(IRSTEP(ITARG,IG))
            IF (IY < 0) THEN
              ESHT(ITARG,IG)=0._DP
            ELSE
              NPES=IGSTEP(ITARG,IG)-200000
              NPBS=NPES-1
              TE=TESTEP(ITARG,IG)
              ESHT(ITARG,IG)=DELTA_SHEATHXB(NPBS,IY)*TE
            END IF
          ELSEIF (IGSTEP(ITARG,IG).LT.200000) THEN
            IX=IXTRI(IPSTEP(ITARG,IG))
            IF (IX < 0) THEN
              ESHT(ITARG,IG)=0._DP
            ELSE
              NPES=IGSTEP(ITARG,IG)-100000
              NPBS=NPES-1
              TE=TESTEP(ITARG,IG)
              ESHT(ITARG,IG)=DELTA_SHEATHYB(IX,NPBS)*TE
            END IF
          ENDIF
        ENDIF
C
        ELTEST(ITARG,IG)=0.  ! ELSTEP MAY ALREADY HAVE BEEN SET IN CALL TO FCT. STEP
        DO 6009 IPLS=1,NPLSI
          IF (FLSTEP(IPLS,ITARG,IG).EQ.0.D0) GOTO 6009
C
          IPLSTI=MPLSTI(IPLS)
          IPLSV=MPLSV(IPLS)
          IF (IGSTEP(ITARG,IG).GT.200000) THEN
C  CHECK BOHM CRITERION AT "POLOIDAL" TARGET SURFACE COMPONENTS
            ITRI=IRSTEP(ITARG,IG)
!pb            NPES=IGSTEP(ITARG,IG)-200000
            NPES=IPSTEP(ITARG,IG)
          ELSEIF (IGSTEP(ITARG,IG).LT.200000) THEN
C  CHECK BOHM CRITERION AT "POLOIDAL" TARGET SURFACE COMPONENTS
            ITRI=IRSTEP(ITARG,IG)
!pb            NPES=IGSTEP(ITARG,IG)-100000
            NPES=IPSTEP(ITARG,IG)
          END IF
          VT=SQRT(2.*TISTEP(IPLSTI,ITARG,IG)/BMASS(IPLS))*CVEL2A
C  VELOCITY COMPONENT NORMAL TO POLOIDAL TARGET SURFACE
C  I.E., POLOIDAL COMPONENT V-POL
C  ASSUMING ORTHOGONAL TARGET
          PM1=(PTRIX(NPES,ITRI)*VXSTEP(IPLSV,ITARG,IG)+
     .         PTRIY(NPES,ITRI)*VYSTEP(IPLSV,ITARG,IG))
C  VELOCITY COMPONENT PARALLEL TO POLOIDAL TARGET SURFACE
C  I.E., RADIAL PLUS TOROIDAL COMPONENT, V-RAD + V-TOR
C  AGAIN: ASSUMING ORTHOGONAL TARGET
          VPX=VXSTEP(IPLSV,ITARG,IG)-PM1*PTRIX(NPES,ITRI)
          VPY=VYSTEP(IPLSV,ITARG,IG)-PM1*PTRIY(NPES,ITRI)
          VPZ=VZSTEP(IPLSV,ITARG,IG)-0.
          PN1=SQRT(VPX**2+VPY**2+VPZ**2)
          PERW=0.
          PARW=0.
          IF (VT.GT.0.) THEN
            PERW=PM1/VT
            PARW=PN1/VT
          ENDIF
C
          CS=SQRT((1.*TISTEP(IPLSTI,ITARG,IG)+
     .                TESTEP(ITARG,IG))/BMASS(IPLS))*CVEL2A
C THE MACH NUMBER BOUNDARY CONDITION ONLY AFFECTS THE PARALLEL TO B
C MOMENTUM, I.E., NOT THE RADIAL VELOCITY
          VTEST=SQRT(PM1**2+VPZ**2)
          VTEST=VTEST/(CS+EPS60)
          VR=SQRT(VPX**2+VPY**2)
          VTEST2=VPZ/(CS+EPS60)
          IF (TRCINT) THEN
            WRITE (iunout,*) 'IPL,ITG,IG,MACH_PAR',
     .                        IPLS,ITARG,IG,VTEST
C           WRITE (iunout,*) 'POL., TOR., RAD. (CM/S) ',PM1,VPZ,VR
C           CALL EIRENE_LEER(1)
C
          ENDIF
C
C  BOHM CRITERION CHECK DONE
C
C  ELTEST: TOTAL ION ENERGY FLUX ONTO TARGET:EMAXW + ESHET
C
C  NEXT: TARGET MAXW. ENERGY FLUXES
C  EADD=  IN EV, SUCH THAT EADD*PARTICLE FLUX = ENERGY FLUX
          DRR=RRSTEP(ITARG,IG+1)-RRSTEP(ITARG,IG)
C  ENERGY FLUX DEFINED WITH PARAMETERS IN INPUT BLOCK 7
          IF (NEM.EQ.1) THEN
            EADD=SORENI(ITARG)
          ELSEIF (NEM.EQ.2.OR.NEM.EQ.3) THEN
            EADD=SORENI(ITARG)*TISTEP(IPLSTI,ITARG,IG)+SORENE(ITARG)*
     .           TESTEP(ITARG,IG)
          ELSEIF (NEM.GE.4 .AND. NEM.LE.7) THEN
            PERWI=PERW/SQRT(BMASS(IPLS)/RMASSP(IPLS))
            PARWI=PARW/SQRT(BMASS(IPLS)/RMASSP(IPLS))
            EADD=EIRENE_EMAXW(TISTEP(IPLSTI,ITARG,IG),PERWI,PARWI)
          ELSEIF (NEM.EQ.8 .OR. NEM.EQ.9) THEN
C  ENERGY FLUX ELSTEP IS ALREADY DEFINED BY B2-BOUNDARY CONDITIONS (SUM: IPLS=0?)
            EADD=ELSTEP(IPLS,ITARG,IG)/FLSTEP(IPLS,ITARG,IG)
          ENDIF
          EMAXW=EADD
          ESUM=EMAXW*FLSTEP(IPLS,ITARG,IG)
          ELTEST(ITARG,IG)=ELTEST(ITARG,IG)+ESUM
          EEMAX=EEMAX+ESUM*DRR

C  ADD ENERGY GAIN BY SHEATH ACCELERATION TO TOTAL
          IF (NEM.EQ.3.OR.NEM.EQ.5.OR.NEM.EQ.7.OR.NEM.EQ.9) THEN
            ESHEATH=NCHRGP(IPLS)*ESHT(ITARG,IG)
            ESUM=ESHEATH*FLSTEP(IPLS,ITARG,IG)
            ELTEST(ITARG,IG)=ELTEST(ITARG,IG)+ESUM
            EESHT=EESHT+ESUM*DRR
          ENDIF

 6009   CONTINUE  ! IPLS loop

 6011 CONTINUE    ! IG,  CELL ALONG TARGET
C
      CALL EIRENE_LEER(1)
      WRITE (iunout,*) 'TARGET DATA: TARGET NO. ITARG=ISTRA= ',ITARG
      WRITE (iunout,*) TXTSOU(ISTRA)
      WRITE (iunout,'(1X,A3,9A11,A3)')
     .'IG','ARC','P-FLUX','E-FLUX','TE','TI','SHEATH/TE',
     . 'VXSTEP','VYSTEP','VZSTEP'
      DO 6100 IG=1,NRWL(ITARG)-1

        IF (IGSTEP(ITARG,IG).GT.200000) THEN
          IF (ORI(ITARG,IG).LT.0) NSEW='W'
          IF (ORI(ITARG,IG).GT.0) NSEW='E'
        ENDIF
        IF (IGSTEP(ITARG,IG).LT.200000) THEN
          IF (ORI(ITARG,IG).LT.0) NSEW='S'
          IF (ORI(ITARG,IG).GT.0) NSEW='N'
        ENDIF
        WRITE (iunout,'(1X,I3,1P,9E11.3,3X,A1)')
     .             IG,RRSTEP(ITARG,IG),FLSTEP(0,ITARG,IG),
     .             ELTEST(ITARG,IG),
     .             TESTEP(ITARG,IG),TISTEP(1,ITARG,IG),
     .             ESHT(ITARG,IG)/(TESTEP(ITARG,IG)+EPS60),
     .             VXSTEP(1,ITARG,IG),VYSTEP(1,ITARG,IG),
     .             VZSTEP(1,ITARG,IG),NSEW
 6100 CONTINUE
      WRITE (iunout,'(1X,I3,1P,1E11.3)') NRWL(ITARG),
     .                                 RRSTEP(ITARG,NRWL(ITARG))
      CALL EIRENE_MASR1 ('EEMAX   ',EEMAX)
      CALL EIRENE_MASR1 ('EESHT   ',EESHT)
C
      ETOT=EEMAX+EESHT
      EFLX(ITARG)=EEMAX+EESHT
      WRITE (iunout,*) 'PARTICLE FLUX(IPLS), IPLS=1,NPLSI '
      WRITE (iunout,'(1X,1P,6E12.4)') (FLTOT(ISP,ITARG),ISP=1,NPLSI)
      CALL EIRENE_LEER(1)
      WRITE (iunout,*) 'ENERGY FLUX '
      WRITE (iunout,'(1X,1P,1E12.4)') EFLX(ITARG)
      CALL EIRENE_LEER(2)
C
 6300 CONTINUE
C
C  SET SOME OTHER DATA SPECIFIC FOR EIRENE CODE REQUIREMENTS
C  STATEMENT NO. 6500 ---> 6999
C
 6500 CONTINUE

      DEALLOCATE (NUMSID)
      DEALLOCATE (NUMTRI)
      DEALLOCATE (TORL)
      DEALLOCATE (ESHT)
      DEALLOCATE (ELTEST)
      DEALLOCATE (ORI)
C
C
      RETURN
  999 CONTINUE
      WRITE (iunout,*) 'ERROR IN IF2COP: NGITT TOO SMALL '
      CALL EIRENE_EXIT_OWN(1)
      RETURN
C
C
      ENTRY EIRENE_IF3COP(ISTRAA,ISTRAE,NEW_ITER)
C
C
      WRITE (iunout,*) ' IF3COP IS CALLED, ISTRAA,ISTRAE '
      WRITE (iunout,*) ISTRAA,ISTRAE
      LSHORT=.FALSE.
      LSTP3=.TRUE.
      LSTOP=LSTP3
      IFIRST=0
      NDXY=(NDXA-1)*NR1TAL_SAVE+NDYA
      GOTO 99992
C
      ENTRY EIRENE_INTER3(LSTP,IFRST,ISTRAA,ISTRAE,NEW_ITER)
C
C  ENTRY FOR SHORT CYCLE FROM SUBR. EIRSRT
C
C  IFIRST=0: RESTORE DATA FROM A PREVIOUS EIRENE RUN, SET REFERENCE
C            DATA FOR "STOP-CRITERION" SNIS,SEES,SEIS
C  IFIRST>0: MODIFY SOURCE TERMS ACCORDING TO NEW PLASMA CONDITIONS,
C            COMPARE INTEGRALS WITH SNIS,...., AND DECIDE TO STOP OR
C            CONTINUE SHORT CYCLE (LSTOP)
C
      LSHORT=.TRUE.
      LSTP3=LSTP
      LSTOP=LSTP
      IFIRST=IFRST
      NDXY=(NDXA-1)*NR1TAL_SAVE+NDYA
C
99992 CONTINUE

      IF (.NOT.ALLOCATED(CHPS)) THEN
        ALLOCATE (CHPS(NFL))
        ALLOCATE (SNIS(0:NFL))
        ALLOCATE (CHMOS(NFL))
        ALLOCATE (SMOS(0:NFL))
        ALLOCATE (SCALN(0:NFL))
        ALLOCATE (SNIS0(NSTRA,0:NFL))
        ALLOCATE (SMOS0(NSTRA,0:NFL))

        ALLOCATE (RESSNI(0:NSTRA,NFL))
        ALLOCATE (RESSMO(0:NSTRA,NFL))
        ALLOCATE (RESSEE(0:NSTRA))
        ALLOCATE (RESSEI(0:NSTRA))

        ALLOCATE (FLXEIR(NSTRA))

        CALL EIRENE_ALLOC_BRASPOI
        CALL EIRENE_ALLOC_EIRBRA(NDX, NDY, NFL, NSTRA,IFOFF)
C
        RESSNI = 0._DP
        RESSMO = 0._DP
        RESSEE = 0._DP
        RESSEI = 0._DP
      ENDIF
C
      IF (.NOT.LSHORT) THEN
        RESSNI(ISTRAA:ISTRAE,:) = 0._DP
        RESSMO(ISTRAA:ISTRAE,:) = 0._DP
        RESSEE(ISTRAA:ISTRAE) = 0._DP
        RESSEI(ISTRAA:ISTRAE) = 0._DP
      ENDIF

      DO 10000 ISTRAI=ISTRAA,ISTRAE
C
C  FIRSTLY INITIALIZE SOURCE TERM ARRAYS
C
        sni(:,:,:,istrai) = 0.d0
        smo(:,:,:,istrai) = 0.d0
        see(:,:,istrai) = 0.d0
        sei(:,:,istrai) = 0.d0
C
        IF (XMCP(ISTRAI).LE.1.) GOTO 10000
C
        IF (LSHORT) GOTO 7000
C
C  READ DATA FROM STRATUM NO. ISTRAI BACK INTO WORKING SPACE
C  IF REQUIRED
C
        IF (ISTRAI.EQ.IESTR) THEN
C  NOTHING TO BE DONE
        ELSEIF (NFILEN.EQ.1.OR.NFILEN.EQ.2) THEN
          IESTR=ISTRAI
          CALL EIRENE_RSTRT(ISTRAI,NSTRAI,NESTM1,NESTM2,NADSPC,
     .               ESTIMV,ESTIMS,ESTIML,
     .               NSDVI1,SDVI1,NSDVI2,SDVI2,
     .               NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .               NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .               NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .               NSIGI_SPC,TRCFLE)
        ELSE
          WRITE (iunout,*) 'ERROR IN INFCOP: STRATUM ISTRAI= ',ISTRAI
          WRITE (iunout,*) 'IS NOT AVAILABLE. EXIT CALLED'
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
C
C  DATA TRANSFER BACK FROM EIRENE TO EXTERNAL CODE
C  STATEMENT NO 7000 ---> 7999
C
 7000   CONTINUE
C
        IF (.NOT.LSHORT) CALL EIRENE_SAVE_TALLIES(ISTRAI)
C
C  SCALE SURFACE SOURCES PER UNIT FLUX, FOR OTHER SOURCES USE
C  EIRENE SCALINGS
        IF (ISTRAI.LE.NTARGI.AND.WTOTP(0,ISTRAI).NE.0.) THEN
C  FLUX FROM EIRENE TO PLASMA CODE: NEGATIVE
          FLX=-WTOTP(0,ISTRAI)
          FLXI=1./FLX
          FLXEIR(ISTRAI)=1._DP
        ELSEIF (ISTRAI.LE.NTARGI.AND.WTOTP(0,ISTRAI).EQ.0.) THEN
          WRITE (iunout,*) 'NO PLASMA FLUX FROM STRATUM NO. ISTRAI= ',
     .                      ISTRAI
          WRITE (iunout,*)
     .       'NO DATA RETURNED TO PLASMA CODE FOR THIS STRATUM'
          GOTO 7999
        ELSEIF (ISTRAI.GT.NTARGI) THEN
          FLXI=1.

C  IF THE SOURCE STRENGTH IS TO BE CHANGED DURING THE SHORT CYCLE (E.G.: VOL-REC)
C  THEN FLXEIR HAS TO BE RESET TO SCALE TO NEW SOURCE STRENGTH DURING SHORT CYCLE
          FLXEIR(ISTRAI)=1._DP
        ENDIF

C
        CHPM  = 0._DP
        CHMOM = 0._DP
        CHEEM = 0._DP
        CHEIM = 0._DP
C
        IF (.NOT.LSHORT) GOTO 7400

        IST_RATE = ITS(ISTRAI)
        RTIS => RTS(IST_RATE)%RTA

        COPV=0.D0
        CPMUL => COPVS(ISTRAI)%PMUL
        DO WHILE (ASSOCIATED(CPMUL))
          ICPV=CPMUL%IART
          IN=CPMUL%ICM
          COPV(ICPV,IN)=CPMUL%VALUEM
          CPMUL => CPMUL%NXTMUL
        ENDDO

        MAPL=0.D0
        CPMUL => MAPLS(ISTRAI)%PMUL
        DO WHILE (ASSOCIATED(CPMUL))
          IPLS=CPMUL%IART
          IN=CPMUL%ICM
          MAPL(IPLS,IN)=CPMUL%VALUEM
          CPMUL => CPMUL%NXTMUL
        ENDDO

        MMPL=0.D0
        CPMUL => MMPLS(ISTRAI)%PMUL
        DO WHILE (ASSOCIATED(CPMUL))
          IPLS=CPMUL%IART
          IN=CPMUL%ICM
          MMPL(IPLS,IN)=CPMUL%VALUEM
          CPMUL => CPMUL%NXTMUL
        ENDDO

        MIPL=0.D0
        CPMUL => MIPLS(ISTRAI)%PMUL
        DO WHILE (ASSOCIATED(CPMUL))
          IPLS=CPMUL%IART
          IN=CPMUL%ICM
          MIPL(IPLS,IN)=CPMUL%VALUEM
          CPMUL => CPMUL%NXTMUL
        ENDDO

        MPHPL=0.D0
        CPMUL => MPHPLS(ISTRAI)%PMUL
        DO WHILE (ASSOCIATED(CPMUL))
          IPLS=CPMUL%IART
          IN=CPMUL%ICM
          MPHPL(IPLS,IN)=CPMUL%VALUEM
          CPMUL => CPMUL%NXTMUL
        ENDDO
C
C  SHORT LOOP CORRECTION FOR ELECTRON IMPACT IONISATION OF ATOMS
C                        AND BULK ION CHARGE EXCHANGE WITH ATOMS:
C                        PARTICLE AND ENERGY SOURCES
C
C  PARTICLE SOURCE: SPLIT FOR MULTIPLE IPLS SPECIES
        PAPL=0.D0
        CPMUL => PAPLS(ISTRAI)%PMUL
        DO WHILE (ASSOCIATED(CPMUL))
          IPLS=CPMUL%IART
          IN=CPMUL%ICM
          PAPL(IPLS,IN)=CPMUL%VALUEM
          CPMUL => CPMUL%NXTMUL
        ENDDO
C  ELECTRON ENERGY: SINGLE (ELECTRON) SPECIES ARRAY
        EAEL=0.D0
        CPSIM => EAELS(ISTRAI)%PSIM
        DO WHILE (ASSOCIATED(CPSIM))
          IN=CPSIM%ICS
          EAEL(IN)=CPSIM%VALUES
          CPSIM => CPSIM%NXTSIM
        END DO

C  ION ENERGY: SPLIT FOR MULTIPLE IPLS SPECIES
        EAPL=0.D0
        CPMUL => EAPLS(ISTRAI)%PMUL
        DO WHILE (ASSOCIATED(CPMUL))
          IPLS=         CPMUL%IART
          IN=           CPMUL%ICM
          EAPL(IPLS,IN)=CPMUL%VALUEM
          CPMUL =>      CPMUL%NXTMUL
        END DO

        IF (IFIRST.EQ.0) GOTO 7310

        CPMUL => PDENAS(ISTRAI)%PMUL
        DO WHILE (ASSOCIATED(CPMUL))
          IATM=CPMUL%IART
          IN=CPMUL%ICM
          DO IPLS=1,NPLSI
            CHP=CPMUL%VALUEM*
     .          (SPLNWA(IN,IATM,IPLS)-RTIS%SPLODA(IN,IATM,IPLS))*ELCHA
            PAPL(IPLS,IN)=PAPL(IPLS,IN)+CHP
            CHPM(IPLS,IN)=CHPM(IPLS,IN)+CHP
          ENDDO
          CHE=CPMUL%VALUEM*
     .        (SEENWA(IN,IATM)-RTIS%SEEODA(IN,IATM))*ELCHA
          EAEL(IN)=EAEL(IN)+CHE
          CHEEM(IN)=CHEEM(IN)+CHE
          CPMUL => CPMUL%NXTMUL
        ENDDO

        CPMUL => COPVS(ISTRAI)%PMUL
        DO WHILE (ASSOCIATED(CPMUL))
          ICPV=CPMUL%IART
          IF (ICPV.LE.NPLSI) THEN
            IPLS=ICPV
            IN=CPMUL%ICM
            CHI=CPMUL%VALUEM*
     .          (SEINWA(IN,IPLS)-RTIS%SEIODA(IN,IPLS))*ELCHA
            EAPL(IPLS,IN)=EAPL(IPLS,IN)+CHI
            CHEIM(IN)=CHEIM(IN)+CHI
          ENDIF
          CPMUL => CPMUL%NXTMUL
        ENDDO

 7310   CONTINUE

C
C  CORRECTION FOR ELECTRON IMPACT IONISATION AND CX OF ATOMS FINISHED
C
C
C  SHORT LOOP CORRECTION FOR ELECTRON IMPACT DISSOCIATION OF TEST IONS
C

        PIPL=0.D0
        CPMUL => PIPLS(ISTRAI)%PMUL
        DO WHILE (ASSOCIATED(CPMUL))
          IPLS=CPMUL%IART
          IN=CPMUL%ICM
          PIPL(IPLS,IN)=CPMUL%VALUEM
          CPMUL => CPMUL%NXTMUL
        ENDDO

        EIEL=0.D0
        CPSIM => EIELS(ISTRAI)%PSIM
        DO WHILE (ASSOCIATED(CPSIM))
          IN=CPSIM%ICS
          EIEL(IN)=CPSIM%VALUES
          CPSIM => CPSIM%NXTSIM
        END DO

        EIPL=0.D0
        CPMUL => EIPLS(ISTRAI)%PMUL
        DO WHILE (ASSOCIATED(CPMUL))
          IPLS=CPMUL%IART
          IN=CPMUL%ICM
          EIPL(IPLS,IN)=CPMUL%VALUEM
          CPMUL => CPMUL%NXTMUL
        END DO

        IF (IFIRST.EQ.0) GOTO 7330

        CPMUL => PDENIS(ISTRAI)%PMUL
        DO WHILE (ASSOCIATED(CPMUL))
          IION=CPMUL%IART
          IN=CPMUL%ICM
          DO IPLS=1,NPLSI
            CHP=CPMUL%VALUEM *
     .          (SPLNWI(IN,IION,IPLS)-RTIS%SPLODI(IN,IION,IPLS))*ELCHA
            PIPL(IPLS,IN)=PIPL(IPLS,IN)+CHP
            CHPM(IPLS,IN)=CHPM(IPLS,IN)+CHP
            CHI=CPMUL%VALUEM *
     .          (SEINWI(IN,IION)-RTIS%SEIODI(IN,IION))*ELCHA
            EIPL(IPLS,IN)=EIPL(IPLS,IN)+CHI
            CHEIM(IN)=CHEIM(IN)+CHI
          END DO

          CHE=CPMUL%VALUEM *
     .        (SEENWI(IN,IION)-RTIS%SEEODI(IN,IION))*ELCHA
          EIEL(IN)=EIEL(IN)+CHE
          CHEEM(IN)=CHEEM(IN)+CHE
          CPMUL => CPMUL%NXTMUL
        ENDDO

 7330   CONTINUE
C
C
C  CORRECTION FOR TEST IONS FINISHED
C
C
C  SHORT LOOP CORRECTION FOR ELECTRON IMPACT COLLISIONS
C             OF MOLECULES
C

        PMPL=0.D0
        CPMUL => PMPLS(ISTRAI)%PMUL
        DO WHILE (ASSOCIATED(CPMUL))
          IPLS=CPMUL%IART
          IN=CPMUL%ICM
          PMPL(IPLS,IN)=CPMUL%VALUEM
          CPMUL => CPMUL%NXTMUL
        ENDDO

        EMEL=0.D0
        CPSIM => EMELS(ISTRAI)%PSIM
        DO WHILE (ASSOCIATED(CPSIM))
          IN=CPSIM%ICS
          EMEL(IN)=CPSIM%VALUES
          CPSIM => CPSIM%NXTSIM
        END DO

        EMPL=0.D0
        CPMUL => EMPLS(ISTRAI)%PMUL
        DO WHILE (ASSOCIATED(CPMUL))
          IPLS=CPMUL%IART
          IN=CPMUL%ICM
          EMPL(IPLS,IN)=CPMUL%VALUEM
          CPMUL => CPMUL%NXTMUL
        END DO

        IF (IFIRST.EQ.0) GOTO 7350

        CPMUL => PDENMS(ISTRAI)%PMUL
        DO WHILE (ASSOCIATED(CPMUL))
          IMOL=CPMUL%IART
          IN=CPMUL%ICM
          DO IPLS=1,NPLSI
            CHP=CPMUL%VALUEM*
     .          (SPLNWM(IN,IMOL,IPLS)-RTIS%SPLODM(IN,IMOL,IPLS))*ELCHA
            PMPL(IPLS,IN)=PMPL(IPLS,IN)+CHP
            CHPM(IPLS,IN)=CHPM(IPLS,IN)+CHP
          END DO
          CHE=CPMUL%VALUEM*
     .        (SEENWM(IN,IMOL)-RTIS%SEEODM(IN,IMOL))*ELCHA
          EMEL(IN)=EMEL(IN)+CHE
          CHEEM(IN)=CHEEM(IN)+CHE
          CPMUL => CPMUL%NXTMUL
        ENDDO

 7350   CONTINUE
C  CORRECTION FOR ELECTRON IMPACT DISSOCIATION OF MOLECULES FINISHED

C
C  SHORT LOOP CORRECTION FOR VOLUME RECOMBINATION PROCESSES  (UNFINISHED)
C

        PPLODA=0.D0
        CPMUL => PPPL_COPS(ISTRAI)%PMUL
        DO WHILE (ASSOCIATED(CPMUL))
          IPLS=CPMUL%IART
          IN=CPMUL%ICM
          PPLODA(IPLS,IN)=CPMUL%VALUEM
          CPMUL => CPMUL%NXTMUL
        ENDDO

        CPVODA=0.D0
        CPMUL => CPPVS(ISTRAI)%PMUL
        DO WHILE (ASSOCIATED(CPMUL))
          IPLS=CPMUL%IART
          IN=CPMUL%ICM
          CPVODA(IPLS,IN)=CPMUL%VALUEM
          CPMUL => CPMUL%NXTMUL
        ENDDO

        EPLODA=0.D0
        CPMUL => EPPL_COPS(ISTRAI)%PMUL
        DO WHILE (ASSOCIATED(CPMUL))
          IPLS=CPMUL%IART
          IN=CPMUL%ICM
          EPLODA(IPLS,IN)=CPMUL%VALUEM
          CPMUL => CPMUL%NXTMUL
        END DO

        EPEODA=0.D0
        CPSIM => EPELS(ISTRAI)%PSIM
        DO WHILE (ASSOCIATED(CPSIM))
          IN=CPSIM%ICS
          EPEODA(IN)=CPSIM%VALUES
          CPSIM => CPSIM%NXTSIM
        END DO
C
C
C  SHORT LOOP CORRECTION FINISHED
C
 7400   CONTINUE
C
C
C  ADD CONTRIBUTIONS TO SOURCE RATES, FROM PRIMARY VOLUME RECOMBINATION SOURCE
C
        PPPL_COP = 0.D0
        MPPL_COP = 0.D0
        EPPL_COP = 0.D0
        EPEL_COP = 0.D0

        IF (NLVOL(ISTRAI).AND.NLPLS(ISTRAI)) THEN
C
          RECTOT = 0._DP
cdr  only one bulk ion species per volume source stratum supported
          IPLS=NSPEZ(ISTRAI)  ! RANGE CHECK FOR IPLS ALREADY DONE IN SAMVOL

          CNDYNP=AMUA*RMASSP(IPLS)
          IPLSTI = MPLSTI(IPLS)
          DO 7472 IIRC=1,NPRCI(IPLS)
              IRRC=LGPRC(IPLS,IIRC)
              SUMN=0.0
              SUMM=0.0
              SUMEI=0.0
              SUMEE=0.0
              lhit = .false.
              DO 7471 IN=1,NTRII
                INC=NCLTAL(IN)
                if ((nstgrd(in) /= 0) .or. lgvac(in,ipls)) cycle
                IF (NSTORDR >= NRAD) THEN
                  RECADD=-TABRC1(IRRC,IN)*DIIN(IPLS,IN)*ELCHA
                  EEADD=  EELRC1(IRRC,IN)*DIIN(IPLS,IN)*ELCHA
                ELSE
                  RECADD=-EIRENE_FTABRC1(IRRC,IN)*DIIN(IPLS,IN)*ELCHA
                  EEADD=  EIRENE_FEELRC1(IRRC,IN)*DIIN(IPLS,IN)*ELCHA
                END IF
                PIADD=PARMOM(IPLS,IN)*RECADD
                EIADD=(1.5*TIIN(IPLSTI,IN)+EDRIFT(IPLS,IN))*RECADD

!pb 21012013
!  if lcoarse: add contribution to tallies only once per fine grid cell
!  tallies are to be scaled with voltal
                if (.not.lhit(inc)) then
                  PPPL_COP(IPLS,INC)=PPPL_COP(IPLS,INC)+RECADD
                  MPPL_COP(IPLS,INC)=MPPL_COP(IPLS,INC)+PIADD
                  EPPL_COP(IPLS,INC)=EPPL_COP(IPLS,INC)+EIADD
                  EPEL_COP(INC)=EPEL_COP(INC)+EEADD
                  lhit(inc) = .true.
                end if

! ALL INPUT TALLIES, TAB.., FTAB.., VOL, ARE ON UNDERLYING FINE GRID.
! and hence: RECADD,PIADD,EIADD,EEADD also on fine grid.
!  here we scale with volume of triangle cell
!  this needs to be done per triangle
                SUMN=SUMN+RECADD*VOL(IN)
                SUMM=SUMM+PIADD*VOL(IN)
                SUMEI=SUMEI+EIADD*VOL(IN)
                SUMEE=SUMEE+EEADD*VOL(IN)

 7471         CONTINUE  ! loop over grid

              RECTOT = RECTOT + SUMN
              WRITE (iunout,*) 'IPLS,IRRC ',IPLS,IRRC
              CALL EIRENE_MASR4('SUMN, SUMM, SUMEI, SUMEE        ',
     .                     SUMN,SUMM,SUMEI,SUMEE)
 7472     CONTINUE

cdr
CC SUMN_OLD=WTOTP ???
          IF (.NOT.LSHORT) SUMN_OLD=RECTOT

C  RESCALE PPPL_COP,.... TO SOURCE STRENGTH FROM LAST FULL EIRENE RUN
C  BECAUSE ALSO PAPL,.... ARE SCALED LIKE THIS
          PPPL_COP=PPPL_COP*SUMN_OLD/RECTOT
          MPPL_COP=MPPL_COP*SUMN_OLD/RECTOT
          EPPL_COP=EPPL_COP*SUMN_OLD/RECTOT
          EPEL_COP=EPEL_COP*SUMN_OLD/RECTOT

C  RESCALE SOURCES FOR B2 (PPPL_COP,PAPL,....) TO NEW SOURCE STRENGTH
          FLXEIR(ISTRAI)=RECTOT/SUMN_OLD
cdr
          IF (LSHORT) THEN
            DO IPLS=1,NPLSI
              CHPM(IPLS,1:NSBOX_TAL) = CHPM(IPLS,1:NSBOX_TAL) +
     .             PPPL_COP(IPLS,1:NSBOX_TAL) - PPLODA(IPLS,1:NSBOX_TAL)
              CHEIM(1:NSBOX_TAL) = CHEIM(1:NSBOX_TAL) +
     .             EPPL_COP(IPLS,1:NSBOX_TAL) - EPLODA(IPLS,1:NSBOX_TAL)
            END DO

            CHEEM(1:NSBOX_TAL) = CHEEM(1:NSBOX_TAL) +
     .            EPEL_COP(1:NSBOX_TAL) - EPEODA (1:NSBOX_TAL)
          END IF

C
C  SAVE RECOMBINATION SOURCES FOR USE BY SHORT CYCLE
C
          IF (.NOT.LSHORT) THEN
            DO IPLS=1,NPLSI
              DO IN=1,NSBOX_TAL
                IF (PPPL_COP(IPLS,IN) .NE. 0.D0) THEN
!pb               ALLOCATE(CPMUL)
                  CPMUL => EIRENE_NEW_MULARR()
                  CPMUL%IART = IPLS
                  CPMUL%ICM = IN
                  CPMUL%VALUEM = PPPL_COP(IPLS,IN)
                  CPMUL%NXTMUL => PPPL_COPS(ISTRAI)%PMUL
                  PPPL_COPS(ISTRAI)%PMUL => CPMUL
                ENDIF
                IF (MPPL_COP(IPLS,IN) .NE. 0.D0) THEN
!pb               ALLOCATE(CPMUL)
                  CPMUL => EIRENE_NEW_MULARR()
                  CPMUL%IART = IPLS
                  CPMUL%ICM = IN
                  CPMUL%VALUEM = MPPL_COP(IPLS,IN)
                  CPMUL%NXTMUL => CPPVS(ISTRAI)%PMUL
                  CPPVS(ISTRAI)%PMUL => CPMUL
                ENDIF
                IF (EPPL_COP(IPLS,IN) .NE. 0.D0) THEN
!pb               ALLOCATE(CPMUL)
                  CPMUL => EIRENE_NEW_MULARR()
                  CPMUL%IART = IPLS
                  CPMUL%ICM = IN
                  CPMUL%VALUEM = EPPL_COP(IPLS,IN)
                  CPMUL%NXTMUL => EPPL_COPS(ISTRAI)%PMUL
                  EPPL_COPS(ISTRAI)%PMUL => CPMUL
                ENDIF
              ENDDO
            ENDDO
            DO IN=1,NSBOX_TAL
              IF (EPEL_COP(IN) .NE. 0.D0) THEN
!pb             ALLOCATE(CPSIM)
                CPSIM => EIRENE_NEW_SIMARR()
                CPSIM%ICS = IN
                CPSIM%VALUES = EPEL_COP(IN)
                CPSIM%NXTSIM => EPELS(ISTRAI)%PSIM
                EPELS(ISTRAI)%PSIM => CPSIM
              ENDIF
            ENDDO
          END IF

        ENDIF
C
        IF (.NOT.LSYMET) GOTO 7500
C
C  SECONDLY SYMMETRISE EIRENE ARRAYS IN CASE OF UP-DOWN SYMMETRY IN MODEL
C
C
C   THIRDLY WRITE EIRENE ARRAYS (1D) ONTO BRAAMS ARRAYS (2D)
C   AND RESCALE TO PROPER UNITS: #/CELL/STRATUM FLUX
C   # STANDS FOR PARTICLES (SNI), MOMENTUM (SMO)
C   AND ENERGY (SEE,SEI)
C
 7500   CONTINUE

        IF (ISTRAI <= NTARGI) THEN
          FLX_EIR = 1._DP
        ELSE
          FLX_EIR = FLXEIR(ISTRAI)
        END IF

!pb 22012013 copv
cdr  copv(icp1+1,...) to copv(icp5+1,...)  linear algebraic tallies
cdr                                       scored in upfcop, so also their variances
cdr                                       can be set without need for covariance matrices
        icp1=nplsi
        icp2=2*nplsi
        icp3=3*nplsi

cdr  fill bulk particle source rate sni(...ifl) from all contributing
cdr  test particle sources papl,pmpl,pipl,pppl (...,ipls)
cdr  test particle may have scored on a finer mesh (lcoarse=.false)

        DO 7510 IFL=1,NFLA
          CHPS(IFL)=0.
          SNIS(IFL)=0.
          CHMOS(IFL)=0.
          SMOS(IFL)=0.

cdr  fill bulk particle source rate sni(...ifl) from all contributing
cdr  test particle sources papl,pmpl,pipl,pppl (...,ipls)
cdr  test particle may have scored on a finer mesh (lcoarse=.false)

          DO 7510 IPLS=1,NPLSI
            IF (IFLB(IPLS).NE.IFL) GOTO 7510
            IPLSV=MPLSV(IPLS)
c  ipls contributes to plasma code species ifl
            DO IX=1,NDXA
            IF (LLCUT(IX)) CYCLE
              DO IY=1,NDYA
C  multiple grid options
C  either (not lcoarse) LOOP OVER ALL TRIANGLES CONTRIBUTING TO ix,iy CELL
C  or     (    lcoarse) already scored on B2.5 grid cell INC=IY+(IX-1)*NR1TAL
                IF (.NOT.LCOARSE) THEN

                CURPOI => HEADS(IY,IX)%P
                DO WHILE (ASSOCIATED(CURPOI))
                  IT=CURPOI%TRIANGLE
                  INC=NCLTAL(IT)  !here: inc=it: all tally data are on fine grid
                  CURPOI=>CURPOI%NEXT

                  SNICL=(PAPL(IPLS,INC)+PMPL(IPLS,INC)+PIPL(IPLS,INC)+
     .                   PPPL_COP(IPLS,INC))*VOLTAL(INC)*FLX_EIR
                  SNI(IX,IY,IFL,ISTRAI)=SNI(IX,IY,IFL,ISTRAI)+SNICL
                  SNIS(IFL)=SNIS(IFL)+ SNICL
                  CHPS(IFL)=CHPS(IFL)+CHPM(IPLS,INC)*VOLTAL(INC)
                ENDDO
                ELSEIF (LCOARSE) THEN
                  INC=IY+(IX-1)*NR1TAL_SAVE

                  SNICL=(PAPL(IPLS,INC)+PMPL(IPLS,INC)+PIPL(IPLS,INC)+
     .                   PPPL_COP(IPLS,INC))*VOLTAL(INC)*FLX_EIR
                  SNI(IX,IY,IFL,ISTRAI)=SNI(IX,IY,IFL,ISTRAI)+SNICL
                  SNIS(IFL)=SNIS(IFL)+ SNICL
                  CHPS(IFL)=CHPS(IFL)+CHPM(IPLS,INC)*VOLTAL(INC)
                ENDIF  ! LCOARSE OPTION
              ENDDO  !IY LOOP
            ENDDO  !IX LOOP

cdr  build alternative source rates, from corresponding copv tallies
cdr  scored in UPDLIN.
cdr  check storage on copv tallies, ncpv ??
            if (NCPV.GE.ICP3+NPLS) THEN
c  skip working on internal lin. comb. of tallies unless sufficient storage


            lhit = .false.
            DO ITRI=1,NTRII
              IY=IYTRI(ITRI)
              IX=IXTRI(ITRI)
              if ((ix<=0).or.(iy <= 0)) cycle
              IN=IY+(IX-1)*NR1TAL_SAVE
              if (lhit(in)) cycle
cdr  scale the particle sources from copv(icp1+ipls)
cdr  then add pppl_cop: vol. recombination contribution to part. sources
              cfac = (PAPL(IPLS,IN)+PMPL(IPLS,IN)+PIPL(IPLS,IN)) /
     .               (copv(icp1+ipls,in) + eps60)
              copv(icp1+ipls,in)=(copv(icp1+ipls,in)*cfac +
     .             PPPL_COP(IPLS,IN)) * VOLTAL(IN)*FLX_EIR
cdr  is this now any different from sni set above?


cdr  add pppl_cop contribution to internal energy sources rate
!pb              copv(icp3+3,in)=copv(icp3+3,in) +
!pb     .            0.5_dp*rmassp(ipls)*bvin(ipls,in)**2*PPPL_COP(IPLS,IN)
              copv(icp3+3,in)=copv(icp3+3,in) +
     .            cvrssp(ipls)*bvin(iplsv,itri)**2*PPPL_COP(IPLS,IN)
              lhit(in) = .true.
            end do  ! ITRI LOOP
            copv(icp1+ipls,:) = copv(icp1+ipls,:) * flxi

            ENDIF  ! STORAGE ON COPV
!pb

            IF (.NOT.LSHORT) THEN

!pb replace sigma_cop
              istat_cop = 0
              do i = 1, nsigvi
                if ((iih(i) == ntalm).and.(igh(i) == NPLSI+IPLS)) then
                  istat_cop = i
                  exit
                end if
              end do

              if (istat_cop > 0) then
                DO IX=1,NDXA
                  DO IY=1,NDYA
cdr  scaling was already on coarse grained grid for B2: in=ncltal(it)

                    INC=IY+(IX-1)*NR1TAL_SAVE

                    SNIRES=(PAPL(IPLS,INC)+PMPL(IPLS,INC)+
     .                      PIPL(IPLS,INC))*VOLTAL(INC)*FLX_EIR
                    RESSNI(ISTRAI,IFL)=RESSNI(ISTRAI,IFL)+
     .                                 ABS(SIGMA(ISTAT_COP,INC)*
     .                                 SNIRES/100.D0)
                  END DO
                END DO
              end if

            END IF

cdr   particle sources done.

cdr   next: dwell on parallel momentum sources. still inside ifl and ipls loop
cdr   ipls contributes to plasma code species ifl

            DO IX=1,NDXA
              IF (LLCUT(IX)) CYCLE
              DO IY=1,NDYA
                IF (.NOT.LCOARSE) THEN
                CURPOI => HEADS(IY,IX)%P
                DO WHILE (ASSOCIATED(CURPOI))
                  IT=CURPOI%TRIANGLE
                  INC=NCLTAL(IT)
                  CURPOI=>CURPOI%NEXT

                  SIGNUM=SIGN(1._DP,BVIN(IPLSV,IT))
                  SMOCL=(MAPL(IPLS,INC)+MMPL(IPLS,INC)+MIPL(IPLS,INC)+
     .                   MPPL_COP(IPLS,INC))*
     .                   VOLTAL(INC)*1.D-5*SIGNUM*FLX_EIR
                  SMO(IX,IY,IFL,ISTRAI)=SMO(IX,IY,IFL,ISTRAI)+SMOCL
                  SMOS(IFL)=SMOS(IFL)+SMOCL
                  CHMOS(IFL)=CHMOS(IFL)+CHMOM(IPLS,INC)*VOLTAL(INC)
                ENDDO
              ELSEIF (LCOARSE) THEN
! use BVIN from first triangle belonging the quadrangular cell
cdr  associated(curpoi) is taken for granted here !?

                CURPOI => HEADS(IY,IX)%P
                IT=CURPOI%TRIANGLE
                INC=IY+(IX-1)*NR1TAL_SAVE

                SIGNUM=SIGN(1._DP,BVIN(IPLSV,IT))
                SMOCL=(MAPL(IPLS,INC)+MMPL(IPLS,INC)+MIPL(IPLS,INC)+
     .                 MPPL_COP(IPLS,INC))*
     .                 VOLTAL(INC)*1.D-5*SIGNUM*FLX_EIR
                SMO(IX,IY,IFL,ISTRAI)=SMO(IX,IY,IFL,ISTRAI)+SMOCL
                SMOS(IFL)=SMOS(IFL)+SMOCL
                CHMOS(IFL)=CHMOS(IFL)+CHMOM(IPLS,INC)*VOLTAL(INC)
              ENDIF  ! LCOARSE OPTION

              ENDDO  !IY LOOP
            ENDDO  !IX LOOP


cdr  build alternative source rates, from corresponding copv tallies
cdr  scored in updlin.
cdr  check storage on copv tallies, ncpv ??
            if (NCPV.GE.ICP3+NPLS) THEN
c  skip working on internal lin. comb. of tallies, unless sufficient storage

            lhit = .false.
            DO ITRI=1,NTRII
              IY=IYTRI(ITRI)
              IX=IXTRI(ITRI)
              if ((ix<=0).or.(iy <= 0)) cycle
              IN=IY+(IX-1)*NR1TAL_SAVE
              if (lhit(in)) cycle
              SIGNUM=SIGN(1._DP,BVIN(IPLSV,ITRI))
              cfac = (MAPL(IPLS,IN)+MMPL(IPLS,IN)+MIPL(IPLS,IN))/
     .               (copv(icp2+ipls,in) + eps60)
              copv(icp2+ipls,in)=(copv(icp2+ipls,in)*cfac +
     .                MPPL_COP(IPLS,IN))*
     .                VOLTAL(IN)*1.D-5*SIGNUM*FLX_EIR
!pb 30012013 sei internal
              copv(icp3+3,in)=copv(icp3+3,in) -
     .            bvin(iplsv,itri)*MPPL_COP(IPLS,IN)*SIGNUM*
     .            cveli2/amua*2._DP
              lhit(in) = .true.
            end do  !itri loop
            copv(icp2+ipls,:) = copv(icp2+ipls,:) * flxi

            endif  ! STORAGE ON COPV

            IF (.NOT.LSHORT) THEN

!pb replace sigma_cop
              istat_cop = 0
              do i = 1, nsigvi
                if ((iih(i) == ntalm).and.(igh(i) == 2*NPLSI+IPLS)) then
                  istat_cop = i
                  exit
                end if
              end do

              if (istat_cop > 0) then
                DO IX=1,NDXA
                  DO IY=1,NDYA
!pb 21012013 ncltal
!                  CURPOI => HEADS(IY,IX)%P
!                  DO WHILE (ASSOCIATED(CURPOI))
!                    IT=CURPOI%TRIANGLE
!                    INC=NCLTAL(IT)
                     INC=IY+(IX-1)*NR1TAL_SAVE

                     SMORES=(MAPL(IPLS,INC)+MMPL(IPLS,INC)+
     .                       MIPL(IPLS,INC))*
     .                       VOLTAL(INC)*1.D-5*SIGNUM*FLX_EIR
                     RESSMO(ISTRAI,IFL)=RESSMO(ISTRAI,IFL)+
     .                                  ABS(SIGMA(ISTAT_COP,INC)*
     .                                  SMORES/100.D0*1.D5)
                  END DO  !IY LOOP
                END DO !IX LOOP
              end if

            END IF
 7510   CONTINUE  ! close loops nfla and npls

C
        CHEES=0.
        SEES=0.
        DO IX=1,NDXA
          IF (LLCUT(IX)) CYCLE
          DO IY=1,NDYA
C  multiple grid options
C  either (not lcoarse) LOOP OVER ALL TRIANGLES CONTRIBUTING TO ix,iy CELL
C  or     (    lcoarse) already scored on B2.5 grid cell INC=IY+(IX-1)*NR1TAL
            IF (.NOT.LCOARSE) THEN
            CURPOI => HEADS(IY,IX)%P
            DO WHILE (ASSOCIATED(CURPOI))
              IT=CURPOI%TRIANGLE
              INC=NCLTAL(IT) !here: inc=it: all tally data are on fine grid
              CURPOI=>CURPOI%NEXT

              SEE(IX,IY,ISTRAI)=SEE(IX,IY,ISTRAI)+
     .           (EAEL(INC)+EMEL(INC)+EIEL(INC)+
     .            EPEL_COP(INC))*VOLTAL(INC)*ELCHA
              CHEES=CHEES+CHEEM(INC)*VOLTAL(INC)
              SEES=SEES+(EAEL(INC)+EMEL(INC)+EIEL(INC)+
     .                   EPEL_COP(INC))*VOLTAL(INC)
            ENDDO
            ELSEIF (LCOARSE) THEN
              INC=IY+(IX-1)*NR1TAL_SAVE
              SEE(IX,IY,ISTRAI)=SEE(IX,IY,ISTRAI)+
     .           (EAEL(INC)+EMEL(INC)+EIEL(INC)+
     .            EPEL_COP(INC))*VOLTAL(INC)*ELCHA
              CHEES=CHEES+CHEEM(INC)*VOLTAL(INC)
              SEES=SEES+(EAEL(INC)+EMEL(INC)+EIEL(INC)+
     .                   EPEL_COP(INC))*VOLTAL(INC)
            ENDIF  ! LCOARSE OPTION

          ENDDO  !IY LOOP
        ENDDO  !IX LOOP


cdr  build alternative source rates, from corresponding copv tallies
cdr  scored in updlin.
cdr  check storage on copv tallies, ncpv ??
        if (NCPV.GE.ICP3+NPLS) THEN
        lhit = .false.
        DO ITRI=1,NTRII
          IY=IYTRI(ITRI)
          IX=IXTRI(ITRI)
              if ((ix<=0).or.(iy <= 0)) cycle
              IN=IY+(IX-1)*NR1TAL_SAVE
              if (lhit(in)) cycle
!pb 22012013 copv
cdr  electron energy source rate, now on copv icp3+1
              cfac = (EAEL(IN)+EMEL(IN)+EIEL(IN)) /
     .               (copv(icp3+1,in) + eps60)
              copv(icp3+1,in)=(copv(icp3+1,in)*cfac +
     .             EPEL_COP(IN)) * VOLTAL(IN)*ELCHA
cdr  this is now identical to see above ?

          lhit(in) = .true.
        end do  ! itri loop

        copv(icp3+1,:) = copv(icp3+1,:) * flxi

        endif  ! STORAGE ON COPV

        IF (.NOT.LSHORT) THEN

!pb replace sigma_cop
          istat_cop = 0
          do i = 1, nsigvi
            if ((iih(i) == ntalm).and.(igh(i) == ICP3+1)) then
              istat_cop = i
              exit
            end if
          end do

          if (istat_cop > 0) then
            DO IX=1,NDXA
              DO IY=1,NDYA
!pb 21012013 ncltal
!             CURPOI => HEADS(IY,IX)%P
!             DO WHILE (ASSOCIATED(CURPOI))
!               IT=CURPOI%TRIANGLE
!               IN=NCLTAL(IT)

                IN=IY+(IX-1)*NR1TAL_SAVE

                SEERES=(EAEL(IN)+EMEL(IN)+EIEL(IN))*VOLTAL(IN)*FLX_EIR
!pb              RESSEE(ISTRAI)=RESSEE(ISTRAI)+
!pb     .                       ABS(SIGMA_COP(2*NPLSI+1,IN)*
!pb     .                       SEERES/100.D0)
                RESSEE(ISTRAI)=RESSEE(ISTRAI)+
     .                         ABS(SIGMA(ISTAT_COP,IN)*
     .                         SEERES/100.D0)
              END DO
            END DO
          end if
        END IF

C
        CHEIS=0.
        SEIS=0.
        DO IFL=1,NFLA
          DO  IPLS=1,NPLSI
            IF (IFLB(IPLS).NE.IFL) CYCLE

            DO IX=1,NDXA
              IF (LCUT(IX)) CYCLE
              DO IY=1,NDYA
C  multiple grid options
C  either (not lcoarse) LOOP OVER ALL TRIANGLES CONTRIBUTING TO ix,iy CELL
C  or     (    lcoarse) already scored on B2.5 grid cell INC=IY+(IX-1)*NR1TAL
              IF (.NOT.LCOARSE) THEN
                CURPOI => HEADS(IY,IX)%P
                DO WHILE (ASSOCIATED(CURPOI))
                  IT=CURPOI%TRIANGLE
                  INC=NCLTAL(IT) !here: inc=it: all tally data are on fine grid
                  CURPOI=>CURPOI%NEXT

                  SEI(IX,IY,ISTRAI)=SEI(IX,IY,ISTRAI)+
     .                 (EAPL(IPLS,INC)+EMPL(IPLS,INC)+
     .                  EIPL(IPLS,INC)+
     .                  EPPL_COP(IPLS,INC))*VOLTAL(INC)*ELCHA
                  CHEIS=CHEIS+CHEIM(INC)*VOLTAL(IN)
                  SEIS=SEIS+(EAPL(IPLS,INC)+EMPL(IPLS,INC)+
     .                       EIPL(IPLS,INC)+
     .                       EPPL_COP(IPLS,INC))*VOLTAL(INC)
                ENDDO
              ELSEIF (LCOARSE) THEN
                INC=IY+(IX-1)*NR1TAL_SAVE
                SEI(IX,IY,ISTRAI)=SEI(IX,IY,ISTRAI)+
     .               (EAPL(IPLS,INC)+EMPL(IPLS,INC)+
     .                EIPL(IPLS,INC)+
     .                EPPL_COP(IPLS,INC))*VOLTAL(INC)*ELCHA
                CHEIS=CHEIS+CHEIM(INC)*VOLTAL(INC)
                SEIS=SEIS+(EAPL(IPLS,INC)+EMPL(IPLS,INC)+
     .                     EIPL(IPLS,INC)+
     .                     EPPL_COP(IPLS,INC))*VOLTAL(INC)
              ENDIF  ! LCOARSE OPTION

              ENDDO ! IY LOOP
            ENDDO ! IX LOOP

          ENDDO  ! IPLS LOOP
        ENDDO ! IFL LOOP

        if (NCPV.GE.ICP3+NPLS) THEN
c  skip working on internal lin. comb. of tallies, unless sufficient storage

           lhit = .false.
           DO ITRI=1,NTRII
             IY=IYTRI(ITRI)
             IX=IXTRI(ITRI)
             if ((ix<=0).or.(iy <= 0)) cycle
             IN=IY+(IX-1)*NR1TAL_SAVE
             if (lhit(in)) cycle

             cfac = (EAPL(IPLS,IN)+EMPL(IPLS,IN)+EIPL(IPLS,IN)) /
     .              (copv(icp3+2,in) + eps60)
             copv(icp3+2,in)=(copv(icp3+2,in)*cfac +
     .                  EPPL_COP(IPLS,IN)) * VOLTAL(IN)*ELCHA
!pb 30012013 sei internal
             copv(icp3+3,in)=(copv(icp3+3,in) +
     .                  EPPL_COP(IPLS,IN)) * VOLTAL(IN)*ELCHA
             lhit(in) = .true.
           end do

           copv(icp3+2,:) = copv(icp3+2,:) * flxi
           copv(icp3+3,:) = copv(icp3+3,:) / elcha

         endif  ! STORAGE ON COPV


      IF (.NOT.LSHORT) THEN

!pb replace sigma_cop
          istat_cop = 0
          do i = 1, nsigvi
            if ((iih(i) == ntalm).and.(igh(i) == ICP3+2)) then
              istat_cop = i
              exit
            end if
          end do

          if (istat_cop > 0) then
            DO IX=1,NDXA
              DO IY=1,NDYA
!pb 21012013 ncltal
!            CURPOI => HEADS(IY,IX)%P
!            DO WHILE (ASSOCIATED(CURPOI))
!              IT=CURPOI%TRIANGLE
!              IN=NCLTAL(IT)
               IN=IY+(IX-1)*NR1TAL_SAVE
               SEIRES=(EAPL(IPLS,IN)+EMPL(IPLS,IN)+
     .                 EIPL(IPLS,IN))*VOLTAL(IN)*FLX_EIR
!pb                RESSEI(ISTRAI)=RESSEI(ISTRAI)+
!pb     .                       ABS(SIGMA_COP(2*NPLSI+2,IN)*
!pb     .                       SEIRES/100.D0)
                  RESSEI(ISTRAI)=RESSEI(ISTRAI)+
     .                           ABS(SIGMA(ISTAT_COP,IN)*
     .                           SEIRES/100.D0)
!pb 21012013 ncltal
!              CURPOI=>CURPOI%NEXT
!            END DO ! CURPOI
              END DO   !NDYA
            END DO    !NDXA
          end if   !ISTAT_COP
      END IF  !lshort

      WRITE (iunout,*) 'RECYCLING SOURCE FROM IF3COP ',ISTRAI
      WRITE (iunout,8888) sum(SNIS(1:nfla)), seis, sees
C
C   NEXT:
C   IF LSHORT OR IFIRST.GT.0     : CRITERION TO STOP SHORT CYCLE,
C   IF NOT LSHORT OR IFIRST.EQ.0 : ONLY RESCALE SURFACE SOURCE STRATA
C                  UNITS: # PER UNIT TARGET PLATE FLUX
C
        IF (IFIRST.EQ.0) THEN
C
          SNIS0(ISTRAI,0)=0.
          SMOS0(ISTRAI,0)=0.
          DO 7550 IFL=1,NFLA
            SNIS0(ISTRAI,0)=SNIS0(ISTRAI,0)+SNIS(IFL)*FLXI
            SNIS0(ISTRAI,IFL)=SNIS(IFL)*FLXI
            SMOS0(ISTRAI,0)=SMOS0(ISTRAI,0)+SMOS(IFL)*FLXI
            SMOS0(ISTRAI,IFL)=SMOS(IFL)*FLXI
 7550     CONTINUE
          SEES0(ISTRAI)=SEES*FLXI
          SEIS0(ISTRAI)=SEIS*FLXI
C
        ELSEIF (LSHORT.AND.IFIRST.GT.0) THEN
C
          SNIS(0)=0.
          DO 7551 IFL=1,NFLA
            SNIS(0)=SNIS(0)+SNIS(IFL)
            SCALN(IFL)=SNIS0(ISTRAI,IFL)/(SNIS(IFL)+EPS60)
 7551     CONTINUE

          SCALN(0)=SNIS0(ISTRAI,0)/(SNIS(0)+EPS60)
C
C         SCALM=SMOS0(ISTRAI,0)/(SMOS(0)+EPS60)  ???
C         SCALE=SEES0(ISTRAI)/(SEES+EPS60)       ???
C         SCALI=SEIS0(ISTRAI)/(SEIS+EPS60)       ???
C
          SCALM=1.
          SCALE=1.
          SCALI=1.
          DO 7555 IX=0,NDXA+1
            DO 7552 IY=0,NDYA+1
              SEE(IX,IY,ISTRAI)=SEE(IX,IY,ISTRAI)*SCALN(0)
              SEI(IX,IY,ISTRAI)=SEI(IX,IY,ISTRAI)*SCALN(0)
 7552       CONTINUE
 7555     CONTINUE
          DO 7556 IFL=1,NFLA
            DO 7553 IX=0,NDXA+1
              DO 7554 IY=0,NDYA+1
                SNI(IX,IY,IFL,ISTRAI)=SNI(IX,IY,IFL,ISTRAI)*SCALN(IFL)
                SMO(IX,IY,IFL,ISTRAI)=SMO(IX,IY,IFL,ISTRAI)*SCALN(IFL)
 7554         CONTINUE
 7553       CONTINUE
 7556     CONTINUE
C
          LTEST=.TRUE.  ! tactically assume: this stratum will continue in short cycle mode

          IF (LSTOP)
     .      WRITE (iunout,*) 'STOP SHORT CYCLE: ALL B2 TIMESTEPS DONE '

!  DO AT LEAST ONE MORE TIME STEP
CDR  DECIDE FOR THIS CURRENT STRATUM ISTRAI:
CDR  SHORT CYCLE (IMPLICIT CORRECTION) ONLY, OR FULL MONTE CARLO


          DO 7558 IFL=1,NFLA
            TEST=CHPS(IFL)/(SNIS(IFL)+1.D-60)*100.
            write (iunout,*) ' global change in sni,ifl ',test,ifl
            IF (ABS(TEST).GT.CHGP) THEN
              LSTP3=.TRUE.
              LTEST=.FALSE.  !  stop short cycle mode. Full new set of trajectories.
              WRITE (iunout,*) 'STOP SHORT CYCLE: PART. SOURCES: ',
     .                     SNIS(IFL),CHPS(IFL),TEST
              WRITE (iunout,*) 'STRATUM ISTRAI, SPECIES IFL ',
     .                          ISTRAI,IFL
            ENDIF
            TEST=CHMOS(IFL)/(SMOS(IFL)+1.D-60)*100.
            write (iunout,*) ' global change in smo,ifl ',test,ifl
            IF (ABS(TEST).GT.CHGMOM) THEN
              LSTP3=.TRUE.
              LTEST=.FALSE. !  stop short cycle mode. Full new set of  trajectories.
              WRITE (iunout,*) 'STOP SHORT CYCLE: MOMENTUM SOURCE: ',
     .                     SMOS(IFL),CHMOS(IFL),TEST
              WRITE (iunout,*) 'STRATUM ISTRAI, SPECIES IFL ',
     .                          ISTRAI,IFL
            ENDIF
 7558     CONTINUE

          TEST=CHEES/(SEES+1.D-60)*100.
          write (iunout,*) ' global change in see,ifl ',test,ifl
          IF (ABS(TEST).GT.CHGEE) THEN
            LSTP3=.TRUE.
            LTEST=.FALSE. !  stop short cycle mode. Full new set of  trajectories.
            WRITE (iunout,*) 'STOP SHORT CYCLE: EL EN. SOURCE: ',SEES,
     .                   CHEES,TEST
            WRITE (iunout,*) 'STRATUM ISTRAI ',ISTRAI
          ENDIF
          TEST=CHEIS/(SEIS+1.D-60)*100.
          write (iunout,*) ' global change in sei,ifl ',test,ifl
          IF (ABS(TEST).GT.CHGEI) THEN
            LSTP3=.TRUE.
            LTEST=.FALSE. !  stop short cycle mode. Full new set of  trajectories.
            WRITE (iunout,*) 'STOP SHORT CYCLE: ION EN. SOURCE: ',
     .                        SEIS,CHEIS,TEST
            WRITE (iunout,*) 'STRATUM ISTRAI ',ISTRAI
          ENDIF
          IF (LSHORT) LSTP=LSTP3

        ELSE
cdr  here: .not.lshort .and. ifirst.gt.0
        ENDIF

        NLSRON(ISTRAI) = .NOT.LTEST
C
        IF (.NOT.LSHORT) THEN
C
          DO 7560 IX=0,NDXA+1
            DO 7565 IY=0,NDYA+1
              SEE(IX,IY,ISTRAI)=SEE(IX,IY,ISTRAI)*FLXI
              SEI(IX,IY,ISTRAI)=SEI(IX,IY,ISTRAI)*FLXI
 7565       CONTINUE
 7560     CONTINUE
          DO 7570 IFL=1,NFLA
            DO 7580 IX=0,NDXA+1
              DO 7590 IY=0,NDYA+1
                SNI(IX,IY,IFL,ISTRAI)=SNI(IX,IY,IFL,ISTRAI)*FLXI
                SMO(IX,IY,IFL,ISTRAI)=SMO(IX,IY,IFL,ISTRAI)*FLXI
 7590         CONTINUE
 7580       CONTINUE
 7570     CONTINUE
C
        ENDIF
C
C   THIRDLY:
C   INDEX MAPPING BACK TO BRAAMS IMPLEMENTATION OF LINDA GEOMETRY
C
        IF (NCUTL.EQ.NCUTB) GOTO 7600
C
        CALL EIRENE_INDMPI (SNI,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .               NCUTB,NCUTL,NPOINT,NPPLG,NSTRA,ISTRAI)
        CALL EIRENE_INDMPI (SMO,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .               NCUTB,NCUTL,NPOINT,NPPLG,NSTRA,ISTRAI)
        CALL EIRENE_INDMPI (SEE,DUMMY,NDX,NDY,1  ,NDXA,NDYA,1   ,
     .               NCUTB,NCUTL,NPOINT,NPPLG,NSTRA,ISTRAI)
        CALL EIRENE_INDMPI (SEI,DUMMY,NDX,NDY,1  ,NDXA,NDYA,1   ,
     .               NCUTB,NCUTL,NPOINT,NPPLG,NSTRA,ISTRAI)
C
 7600   CONTINUE
C
 7700   CONTINUE
C
 7999 CONTINUE
C
C  DATA TRANSFER BACK TO PLASMA CODE FINISHED FOR STRATUM NO. ISTRAI
C
10000 CONTINUE
C
      RETURN
C
      ENTRY EIRENE_IF4COP
C
      NREC11=NOUTAU
      OPEN (UNIT=11,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=8*NREC11)
      IRC=3
C  WRITE RCCPL
      WRITE (11,REC=IRC) RCCPL
      IF (TRCINT.OR.TRCFLE)
     .    WRITE (iunout,*) 'WRITE 11  RCCPL,   IRC= ',IRC
C     IRC=3   STILL
C  WRITE ICCPL1
      ALLOCATE (IHELP(NOUTAU))
      JC=0
      DO K=1,NPTRGT
        DO J=1,10*NSTEP
          JC=JC+1
          IHELP(JC)=ICCPL1(J,K)
          IF (JC == NOUTAU) THEN
            IRC=IRC+1
            WRITE (11,REC=IRC) IHELP
            IF (TRCINT.OR.TRCFLE)
     .          WRITE (iunout,*) 'WRITE 11  ICCPL1,  IRC= ',IRC
            JC=0
          END IF
        END DO
      END DO
c  write last (incomplete) record of ICCPL1
      IF (JC > 0) THEN
        IHELP(JC+1:NOUTAU) = 0   ! fill up last record, up to full length NOUTAU
        IRC=IRC+1
        WRITE (11,REC=IRC) IHELP
        IF (TRCINT.OR.TRCFLE)
     .      WRITE (iunout,*) 'WRITE 11  ICCPL1,  IRC= ',IRC
      END IF
      DEALLOCATE (IHELP)
C  WRITE ICCPL2
      IRC=IRC+1
      WRITE (11,REC=IRC) ICCPL2
      IF (TRCINT.OR.TRCFLE)
     .    WRITE (iunout,*) 'WRITE 11  ICCPL2,  IRC= ',IRC
      IRC=IRC+1
      WRITE (11,REC=IRC) LCCPL
      IF (TRCINT.OR.TRCFLE)
     .    WRITE (iunout,*) 'WRITE 11  LCCPL,   IRC= ',IRC
C
!pb  LSTP is dummy argument to entry IF3COP, thus not available here
!pb  LSTP3 is stored in IF3COP
!pb   IF (LSHORT) LSTOP=LSTP
      IF (LSHORT) LSTOP=LSTP3
C
      IF (.NOT.LSTOP) RETURN
C
      IF (.NOT.(LBALAN)) GOTO 11000
C
C  BALANCES, SHOULD BE DONE ONLY AT THE END OF B2 RUN
C  AT THE END OF AN EIRENE RUN THE BALANCES MAY BE OFF AT LEAST AT
C  THE BEGINNING OF THE CYCLING PROCEDURE, BECAUSE THE PLASMA STILL
C  HAS TO ADJUST TO THE NEW SOURCES
C
C
C...............................................................................
C
C  CHECK FLUXES AT THE GRID BOUNDARY, SIGN, RECYCLING OR NON-RECYCLING BOUNDARY
C  COUNT FLUXES FROM OUTSIDE INTO GRID AS POSITIVE
C
C...............................................................................
C
C  FIRST: SOUTH EDGE: IY=0
C
C NON RECYCLING FLUXES AT SOUTH EDGE: SFEISY,SFEESY,SFNISY
      SFEISY=0.
      SFEESY=0.
      SFNISY=0.
      NPBS=0
      NPES=NPBS+1
C
      DO 10113 IX=1,NDXA
C
C IS (IX,0) A RECYCLING SOURCE? IF YES, DO NOT COUNT HERE
C
        DO ITARG=1,NTARGI
          DO IPRT=1,NTGPRT(ITARG)
            IF (NIXY(ITARG,IPRT).EQ.2.AND.NDT(ITARG,IPRT).EQ.0) THEN
              IF (IX.GE.NTIN(ITARG,IPRT).AND.
     .            IX.LT.NTEN(ITARG,IPRT)) GOTO 10110
            ENDIF
          ENDDO
        ENDDO

        ITARG=0
        IF (LLCUT(IX)) GOTO 10113
C
C SURFACE NORMAL IS INWARD. HENCE: TAKE ALL FLUXES F...YB POSITIVE
C SIGN OF ADDITIONAL COMPONENT DUE TO INCLINED GRID AS SIGN OF F...YB
C
        SFEISY=SFEISY+FEIYB(IX,0)
        SFEESY=SFEESY+FEEYB(IX,0)
        DO 10111 IFL=1,NFLA
          SFNISY(IFL)=SFNISY(IFL)+FNIYB(IX,0,IFL)
10111   CONTINUE
        GOTO 10113
10110   CONTINUE

C  DO NOT RECYCLE TARGET FLUXES WITH FALSE ORIENTATION
        IF (ITARG.GT.0) THEN
C  ITARG, IPRT KNOWN FROM ABOVE
          FLX=0.
          TIFLX=0.
          PIFLX=0.
          DO 10112 IFL=NSPZI(ITARG,IPRT),NSPZE(ITARG,IPRT)
            FLX=FLX+FNIYB(IX,0,IFL)
            IF (NINCT(ITARG,IPRT)*FNIYB(IX,0,IFL).LT.0) THEN
              WRITE (iunout,*)
     .          'RECYCLING TARGET, BUT WRONG FLOW DIRECTION: '
              WRITE (iunout,*)
     .          'SOUTH,IX,ITARG,IPRT,IF ',IX,ITARG,IPRT,IFL
              SFNISY(IFL)=SFNISY(IFL)+FNIYB(IX,0,IFL)
              PIFLX      =PIFLX      +FNIYB(IX,0,IFL)
              TIFLX=TIFLX+FEIYB(IX,0)*FNIYB(IX,0,IFL)
            ENDIF
10112     CONTINUE

          TIFLX=TIFLX/(FLX+EPS60)
          SFNISY=SFNISY+PIFLX
          SFEISY=SFEISY+TIFLX
        ENDIF

10113 CONTINUE
C
      SFNISY=SFNISY*ELCHA
C
      LNONREC_SY=ANY(SFNISY(1:nfla).NE.0.0).OR.SFEISY.NE.0.0.OR.
     .                                         SFEESY.NE.0.0
      WRITE (37,*) 'NON-RECYCLING FLUXES FROM SOUTH EDGE '
      WRITE (37,8888) sum(SFNISY(1:nfla)),SFEISY,SFEESY
C
C
C  SECOND: NORTH EDGE: IY=NDYA
C
C NON RECYCLING FLUXES AT NORTH EDGE: SFEINY,SFEENY,SFNINY
      SFEINY=0.
      SFEENY=0.
      SFNINY=0.
      NPBS=NDYA
      NPES=NPBS+1
      DO 10118 IX=1,NDXA
C
C IS (IX,NDYA) A RECYCLING SOURCE? IF YES, DO NOT COUNT HERE
C
        DO ITARG=1,NTARGI
          DO IPRT=1,NTGPRT(ITARG)
            IF (NIXY(ITARG,IPRT).EQ.2.AND.NDT(ITARG,IPRT).EQ.NDYA) THEN
              IF (IX.GE.NTIN(ITARG,IPRT).AND.
     .            IX.LT.NTEN(ITARG,IPRT)) GOTO 10115
            ENDIF
          ENDDO
        ENDDO
        ITARG=0
        IF (LLCUT(IX)) GOTO 10118
C
C SURFACE NORMAL IS OUTWARD. HENCE: TAKE ALL FLUXES F...YB NEGATIVE
C SIGN OF ADDITIONAL COMPONENT DUE TO INCLINED GRID AS SIGN OF F...YB
C
        SFEINY=SFEINY-FEIYB(IX,NDYA)
        SFEENY=SFEENY-FEEYB(IX,NDYA)
        DO 10116 IFL=1,NFLA
          SFNINY(IFL)=SFNINY(IFL)-FNIYB(IX,NDYA,IFL)
10116   CONTINUE
        GOTO 10118
10115   CONTINUE

C  DO NOT RECYCLE TARGET FLUXES WITH FALSE ORIENTATION
        IF (ITARG.GT.0) THEN
C  ITARG, IPRT KNOWN FROM ABOVE
          FLX=0.
          TIFLX=0.
          PIFLX=0.
          DO 10117 IFL=NSPZI(ITARG,IPRT),NSPZE(ITARG,IPRT)
            FLX=FLX+FNIYB(IX,NDYA,IFL)
            IF (NINCT(ITARG,IPRT)*FNIYB(IX,NDYA,IFL).LT.0) THEN
              WRITE (iunout,*)
     .          'RECYCLING TARGET, BUT WRONG FLOW DIRECTION: '
              WRITE (iunout,*)
     .          'NORTH,IX,ITARG,IPRT,IFL ',IX,ITARG,IPRT,IFL
              SFNINY(IFL)=SFNINY(IFL)-FNIYB(IX,NDYA,IFL)
              PIFLX      =PIFLX      +FNIYB(IX,NDYA,IFL)
              TIFLX=TIFLX+FEIYB(IX,NDYA)*(-FNIYB(IX,NDYA,IFL))
            ENDIF
10117     CONTINUE
          TIFLX=TIFLX/(FLX+EPS60)
          SFNINY=SFNINY-PIFLX
          SFEINY=SFEINY-TIFLX
        ENDIF

10118 CONTINUE
C
      SFNINY=SFNINY*ELCHA
C
      LNONREC_NY=ANY(SFNINY(1:nfla).NE.0.0).OR.SFEINY.NE.0.0.OR.
     .                                         SFEENY.NE.0.0
      WRITE (37,*) 'NON-RECYCLING FLUXES TO NORTH EDGE '
      WRITE (37,8888) sum(SFNINY(1:nfla)),SFEINY,SFEENY
C
C
C  THIRD: WEST EDGE: IX=0
C
C NON RECYCLING FLUXES AT WEST EDGE: SFEIWX,SFEEWX,SFNIWX
      SFEIWX=0.
      SFEEWX=0.
      SFNIWX=0.
      DO 10123 IY=1,NDYA
C
C IS (0,IY) A RECYCLING SOURCE? IF YES, DO NOT COUNT HERE
C
        DO ITARG=1,NTARGI
          DO IPRT=1,NTGPRT(ITARG)
             IF (NIXY(ITARG,IPRT).EQ.1.AND.NDT(ITARG,IPRT).EQ.0) THEN
               IF (IY.GE.NTIN(ITARG,IPRT).AND.
     .             IY.LT.NTEN(ITARG,IPRT)) GOTO 10120
            ENDIF
          ENDDO
        ENDDO
        ITARG=0
C
        SFEIWX=SFEIWX+FEIXB(0,IY)
        SFEEWX=SFEEWX+FEEXB(0,IY)
        DO 10121 IFL=1,NFLA
          SFNIWX(IFL)=SFNIWX(IFL)+FNIXB(0,IY,IFL)
10121   CONTINUE

cdr init  ... something seems fundamentally wrong here ...
10120   CONTINUE

C  DO NOT RECYCLE TARGET FLUXES WITH FALSE ORIENTATION
        IF (ITARG.GT.0) THEN
C  ITARG, IPRT KNOWN FROM ABOVE
          FLX=0.
          TIFLX=0.
          PIFLX=0.
          DO 10122 IFL=NSPZI(ITARG,IPRT),NSPZE(ITARG,IPRT)
            FLX=FLX+FNIXB(0,IY,IFL)
            IF (NINCT(ITARG,IPRT)*FNIXB(0,IY,IFL).LT.0) THEN
              WRITE (iunout,*)
     .          'RECYCLING TARGET, BUT WRONG FLOW DIRECTION: '
              WRITE (iunout,*) 'WEST,IY,ITARG,IPRT,IFL ',
     .                               IY,ITARG,IPRT,IFL
              SFNIWX(IFL)=SFNIWX(IFL)+FNIXB(0,IY,IFL)
              PIFLX      =PIFLX      +FNIXB(0,IY,IFL)
              TIFLX=TIFLX+FEIXB(0,IY)*FNIXB(0,IY,IFL)
            ENDIF
10122     CONTINUE
          TIFLX=TIFLX/(FLX+EPS60)
          SFNIWX=SFNIWX+PIFLX
          SFEIWX=SFEIWX+TIFLX
        ENDIF
cdr end

10123 CONTINUE
C
      SFNIWX=SFNIWX*ELCHA
C
      LNONREC_WX=ANY(SFNIWX(1:nfla).NE.0.0).OR.SFEIWX.NE.0.0.OR.
     .                                         SFEEWX.NE.0.0
      WRITE (37,*) 'NON-RECYCLING FLUXES FROM WEST EDGE '
      WRITE (37,8888) sum(SFNIWX(1:nfla)),SFEIWX,SFEEWX
C
C
C  FOURTH: EAST EDGE: IX=NDXA
C
C NON-RECYCLING FLUXES AT EAST EDGE: SFEIEX,SFEEEX,SFNIEX
      SFEIEX=0.
      SFEEEX=0.
      SFNIEX=0.
      DO 10128 IY=1,NDYA
C
C IS (NDXA,IY) A RECYCLING SOURCE? IF YES, DO NOT COUNT HERE
C
        DO ITARG=1,NTARGI
          DO IPRT=1,NTGPRT(ITARG)
            IF (NIXY(ITARG,IPRT).EQ.1.AND.NDT(ITARG,IPRT).EQ.NDXA) THEN
              IF (IY.GE.NTIN(ITARG,IPRT).AND.
     .            IY.LT.NTEN(ITARG,IPRT)) GOTO 10125
            ENDIF
          ENDDO
        ENDDO
        ITARG=0
C
        WRITE (iunout,*) 'EAST,IY ',IY
        SFEIEX=SFEIEX-FEIXB(NDXA,IY)
        SFEEEX=SFEEEX-FEEXB(NDXA,IY)
        DO 10126 IFL=1,NFLA
          SFNIEX(IFL)=SFNIEX(IFL)-FNIXB(NDXA,IY,IFL)
10126   CONTINUE

10125   CONTINUE

C  DO NOT RECYCLE TARGET FLUXES WITH FALSE ORIENTATION
        IF (ITARG.GT.0) THEN
C  ITARG, IPRT KNOWN FROM ABOVE
          FLX=0.
          TIFLX=0.
          PIFLX=0.
          DO 10127 IFL=NSPZI(ITARG,IPRT),NSPZE(ITARG,IPRT)
            FLX=FLX+FNIXB(NDXA,IY,IFL)
            IF (NINCT(ITARG,IPRT)*FNIXB(NDXA,IY,IFL).LT.0) THEN
              WRITE (iunout,*)
     .          'RECYCLING TARGET, BUT WRONG FLOW DIRECTION: '
              WRITE (iunout,*) 'EAST,IY,ITARG,IPRT,IFL ',
     .                               IY,ITARG,IPRT,IFL
              SFNIEX(IFL)=SFNIEX(IFL)-FNIXB(NDXA,IY,IFL)
              PIFLX      =PIFLX      +FNIXB(NDXA,IY,IFL)
              TIFLX=TIFLX+FEIXB(NDXA,IY)*(-FNIXB(NDXA,IY,IFL))
            ENDIF
10127     CONTINUE
          TIFLX=TIFLX/(FLX+EPS60)
          SFNIEX=SFNIEX-PIFLX
          SFEIEX=SFEIEX-TIFLX
        ENDIF

10128 CONTINUE
C
      SFNIEX=SFNIEX*ELCHA
C
      LNONREC_EX=ANY(SFNIEX(1:nfla).NE.0.0).OR.SFEIEX.NE.0.0.OR.
     .                                         SFEEEX.NE.0.0
      WRITE (37,*) 'NON-RECYCLING FLUXES TO EAST EDGE '
      WRITE (37,8888) sum(SFNIEX(1:nfla)),SFEIEX,SFEEEX
C
C  NEXT: FLUXES TO THOSE SURFACES, AT WHICH RECYCLING BOUNDARY
C        CONDITIONS ARE SPECIFIED
      CALL EIRENE_LEER(2)
C
10130 CONTINUE
C
      SFEIT(0)=0.
      SFEET(0)=0.
      SFNIT(0,:)=0.
      SHEAE(0)=0.
      SHEAI(0)=0.
      DO 10139 I=1,NTARGI
        SFEIT(I)=0.
        SFEET(I)=0.
        SFNIT(I,:)=0.
        SHEAE(I)=0.
        SHEAI(I)=0.
        DO IPRT=1,NTGPRT(I)
          fniprt = 0.
          IF (NIXY(I,IPRT).EQ.1) THEN
C  BALANCE CONTRIB. X-GRID REC. SOURCE
            NPBS=NDT(I,IPRT)
            NPBC=NPBS+MAX0(0,NINCT(I,IPRT))
            IF (NINCT(I,IPRT) > 0) THEN
! EAST TARGET IN B2
              IF (NPBC > 0) THEN
              IF ((NGHPOL(4,NTIN(I,IPRT),NPBC) > 0) .AND.
     .            (ABS(NPBC-NGHPOL(4,NTIN(I,IPRT),NPBC)) > 1)) THEN
                NPBS=NGHPOL(4,NTIN(I,IPRT),NPBC)
              END IF
              END IF
            ELSE IF  (NINCT(I,IPRT) < 0) THEN
! WEST TARGET IN B2
              IF (NPBC > 0) THEN
              IF ((NGHPOL(2,NTIN(I,IPRT),NPBC) > 0) .AND.
     .            (ABS(NPBC-NGHPOL(2,NTIN(I,IPRT),NPBC)) > 1)) THEN
                NPBS=NGHPOL(2,NTIN(I,IPRT),NPBC)
              END IF
              END IF
            END IF
            DO 10132 IY=NTIN(I,IPRT),NTEN(I,IPRT)-1
              DO 10131 IFL=NSPZI(I,IPRT),NSPZE(I,IPRT)
                IF (NINCT(I,IPRT)*FNIXB(NPBS,IY,IFL).GT.0) THEN
                  SFNIT(I,IFL)=SFNIT(I,IFL)-
     .                   NINCT(I,IPRT)*FNIXB(NPBS,IY,IFL)

                  write (37,*) iprt, npbs, iy,
     .                         FNIXB(NPBS,IY,IFL)
                  fniprt = fniprt + FNIXB(NPBS,IY,IFL)

cdr sheath contributions: count negative for electrons, positive for ions
cdr unfinished:  need to account for charge state of ion species IFL
                  SHEAE(I)=SHEAE(I)+TEB(NPBC,IY)*
     .             NINCT(I,IPRT)*FNIXB(NPBS,IY,IFL)*
     .             (-DELTA_SHEATHXB(NPBS,IY))
                  SHEAI(I)=SHEAI(I)+TEB(NPBC,IY)*
     .             NINCT(I,IPRT)*FNIXB(NPBS,IY,IFL)*
     .             DELTA_SHEATHXB(NPBS,IY)
cdr  sheath done
                ELSE
                  WRITE (iunout,*)
     .              'WRONG ORIENTATION OF W/E-TARGET RECYCLING FLUX '
                  WRITE (iunout,*) 'ITARG, IPRT, IPLS, NDT, IY ',
     .                         I    , IPRT, IFL,   NDT(I,IPRT), IY
                  WRITE (iunout,*) 'FNIX(NDT,IY) ',
     .                              FNIXB(NPBS,IY,IFL)
                ENDIF
10131         CONTINUE  ! IFL

CDR  feixb, feiyb are only defined for sum over species.
cdr  try a proportional allocation of energy fluxes to species.

              SFEIT(I)=SFEIT(I)-NINCT(I,IPRT)*FEIXB(NPBS,IY)
              SFEET(I)=SFEET(I)-NINCT(I,IPRT)*FEEXB(NPBS,IY)
10132       CONTINUE  !  IY=1,ntin.nten

C  BALANCE CONTRIB. FROM Y-GRID RECYCLING SOURCE
          ELSEIF (NIXY(I,IPRT).EQ.2) THEN
            DO 10135 IX=NTIN(I,IPRT),NTEN(I,IPRT)-1
              IF (LLCUT(IX)) GOTO 10135
              DO 10136 IFL=NSPZI(I,IPRT),NSPZE(I,IPRT)
                IF (NINCT(I,IPRT)*FNIYB(IX,NDT(I,IPRT),IFL).GT.0.) THEN
                  SFNIT(I,IFL)=SFNIT(I,IFL)-
     .                   NINCT(I,IPRT)*FNIYB(IX,NDT(I,IPRT),IFL)

                  write (37,*) iprt, ix,NDT(I,IPRT),
     .                         FNIYB(IX,NDT(I,IPRT),IFL)
                  fniprt = fniprt + FNIYB(IX,NDT(I,IPRT),IFL)

cdr sheath contributions: count negative for electrons, positive for ions
cdr unfinished:  need to account for charge state of ion species IFL
                SHEAE(I)=SHEAE(I)+TEB(IX,NDT(I,IPRT))*
     .           NINCT(I,IPRT)*FNIYB(IX,NDT(I,IPRT),IFL)*
     .           (-DELTA_SHEATHYB(IX,NDT(I,IPRT)))
                SHEAI(I)=SHEAI(I)+TEB(IX,NDT(I,IPRT))*
     .           NINCT(I,IPRT)*FNIYB(IX,NDT(I,IPRT),IFL)*
     .           DELTA_SHEATHYB(IX,NDT(I,IPRT))
cdr  sheath done
                ELSE
                  WRITE (iunout,*)
     .              'WRONG ORIENTATION OF S/N-TARGET RECYCLING FLUX '
                  WRITE (iunout,*) 'ITARG, IPRT, IPLS, IX, NDT ',
     .                         I    , IPRT, IFL,   IX, NDT(I,IPRT)
                  WRITE (iunout,*) 'FNIY(IX,NDT) ',
     .                              FNIYB(IX,NDT(I,IPRT),IFL)
                ENDIF
10136         CONTINUE
              SFEIT(I)=SFEIT(I)-NINCT(I,IPRT)*FEIYB(IX,NDT(I,IPRT))
              SFEET(I)=SFEET(I)-NINCT(I,IPRT)*FEEYB(IX,NDT(I,IPRT))
10135       CONTINUE
          ENDIF
        ENDDO
C
        SFNIT(I,:)=SFNIT(I,:)*ELCHA
C
C
        WRITE (37,*) 'FLUXES TO TARGET NO. ',I
        WRITE (37,8888) (SFNIT(I,IFL),IFL=1,NFL),SFEIT(I),SFEET(I)
C
        SFEIT(0)=SFEIT(0)+SFEIT(I)
        SFEET(0)=SFEET(0)+SFEET(I)
        SFNIT(0,:)=SFNIT(0,:)+SFNIT(I,:)
        SHEAE(0)=SHEAE(0)+SHEAE(I)
        SHEAI(0)=SHEAI(0)+SHEAI(I)
10139 CONTINUE
C
      SSNI=0.
      SSEI=0.
      SSEE=0.

      DO 10150 ISTR=1,NSTRAI
        ISTRA = ISTR
        IF (XMCP(ISTRA).LE.1) GOTO 10150
        FLX=0.
        IF (ISTRA.LE.NTARGI) THEN
          FLX=SUM(ABS(SFNIT(ISTRA,1:NFLA)))
        ELSE
          FLX=1.
        ENDIF

        SSN=0.
        SSI=0.
        SSE=0.
        DO 10140 IX=1,NDXA
           DO 10140 IY=1,NDYA
             DO 10141 IFL=1,NFLA
               SSN(IFL)=SSN(IFL)+SNI(IX,IY,IFL,ISTRA)
10141        CONTINUE
             SSI=SSI+SEI(IX,IY,ISTRA)
             SSE=SSE+SEE(IX,IY,ISTRA)
10140   CONTINUE
C
      WRITE (37,*) 'RECYCLING SOURCE RATES, POTENTIAL+RAD. EN. ',ISTRA
      WRITE (37,8888) SSN*FLX,SSI*FLX/ELCHA,SSE*FLX/ELCHA
C
C  TRENNEN VON RAD. UND POTENTIELLER ENERGY IM ELECTRONENKANAL.
C  DAZU ABER TEILCHENQUELLE SPEZIESAUFGELOEST NOETIG.
C
C
C
C     WRITE (37,*) 'RADIATION LOSSES VIA NEUTRAL CHANNEL ',ISTRA
C     WRITE (37,8888) 0.,0.,0.
C

        SSNI(1:NFLA)=SSNI(1:NFLA)+SSN(1:NFLA)*FLX
        SSEI=SSEI+SSI*FLX/ELCHA
        SSEE=SSEE+SSE*FLX/ELCHA
10150 CONTINUE
C
      WRITE (37,*) 'EQUILIBRATION '
      WRITE (37,8888) 0.,B2QIE,-B2QIE
C
C
      WRITE (37,*) 'BREMSSTRAHLUNG '
      WRITE (37,8888) 0.,0.,B2BREM
C
      WRITE (37,*) 'CHARGED IMPURITY RAD.,IONIZ. AND RECOMB. '
      WRITE (37,8888) 0.,0.,B2RAD
C
      WRITE (37,*) 'ELECTRIC FIELD TERMS (PRESSURE GRADIENTS)'
      WRITE (37,8888) 0.,B2VDP,-B2VDP
C
      BALANI=SFEISY+SFEINY+SFEIT(0)+SSEI+B2QIE+B2VDP+
     .       SFEIWX+SFEIEX
      BALANE=SFEESY+SFEENY+SFEET(0)+SSEE+B2BREM+B2RAD-B2QIE+
     .       SFEEWX+SFEEEX-B2VDP
      BALANN(1:NFLA)=SFNISY(1:NFLA)+SFNINY(1:NFLA)+SFNIWX(1:NFLA)+
     .               SFNIEX(1:NFLA)+SFNIT(0,1:NFLA)+SSNI(1:NFLA)
C
      TOTI=ABS(SFEISY+SFEINY)+ABS(SFEIT(0))+
     .     ABS(SHEAI(0))+ABS(SSEI)
      TOTE=ABS(SFEESY+SFEENY)+ABS(SFEET(0))+
     .     ABS(SHEAE(0))+ABS(SSEE)
      TOTN(1:NFLA)=ABS(SFNISY(1:NFLA))+SFNINY(1:NFLA)+
     .             ABS(SFNIT(0,1:NFLA))+ABS(SSNI(1:NFLA))
      RE=BALANE/(TOTE+EPS60)*100.
      RI=BALANI/(TOTI+EPS60)*100.
      RN(1:NFLA)=BALANN(1:NFLA)/(TOTN(1:NFLA)+EPS60)*100.

c  residuals, contributions from noise in source terms.
c  sum over strata
      DO IFL=1,NFLA
        RESSNI(0,IFL) = SUM(RESSNI(1:NSTRAI,IFL))
        RESSMO(0,IFL) = SUM(RESSMO(1:NSTRAI,IFL))
      END DO
      RESSEE(0) = SUM(RESSEE(1:NSTRAI))
      RESSEI(0) = SUM(RESSEI(1:NSTRAI))
C
      CALL EIRENE_LEER (1)
      IF (LBALAN) THEN
        WRITE (iunout,*) 'B2-EIRENE GLOBAL BALANCES '
        WRITE (iunout,*) 'PARTICLE FLUXES (SFNI..) IN AMP'
        WRITE (iunout,*) 'ENERGY FLUXES (SFEI..,SFEE..,) IN WATT'
        CALL EIRENE_LEER(1)
        IF (LNONREC_SY) THEN
          WRITE (iunout,*) ' NON-RECYCLING FLUXES AT SOUTH EDGE '
          CALL EIRENE_MASR2(' SFEISY,SFEESY  ',SFEISY,SFEESY)
          DO IFL=1,NFLA
            WRITE(iunout,'(A,I0,A,ES12.4)') 'SFNISY(IFL)=',IFL,') ',
     .                                       SFNISY(IFL)
          ENDDO
        ENDIF
        IF (LNONREC_NY) THEN
          WRITE (iunout,*) ' NON-RECYCLING FLUXES AT NORTH EDGE'
          CALL EIRENE_MASR2(' SFEINY,SFEENY  ',SFEINY,SFEENY)
          DO IFL=1,NFLA
            WRITE(iunout,'(A,I0,A,ES12.4)') 'SFNINY(IFL)=',IFL,') ',
     .                                       SFNINY(IFL)
          ENDDO
        ENDIF
        IF (LNONREC_WX) THEN
          WRITE (iunout,*) ' NON-RECYCLING FLUXES AT WEST EDGE '
          CALL EIRENE_MASR2(' SFEIWX,SFEEWX  ',SFEIWX,SFEEWX)
          DO IFL=1,NFLA
            WRITE(iunout,'(A,I0,A,ES12.4)') 'SFNIWX(IFL)=',IFL,') ',
     .                                       SFNIWX(IFL)
          ENDDO
        ENDIF
        IF (LNONREC_EX) THEN
          WRITE (iunout,*) ' NON-RECYCLING FLUXES AT EAST EDGE '
          CALL EIRENE_MASR2(' SFEIEX,SFEEEX  ',SFEIEX,SFEEEX)
          DO IFL=1,NFLA
            WRITE(iunout,'(A,I0,A,ES12.4)') 'SFNIEX(IFL)=',IFL,') ',
     .                                       SFNIEX(IFL)
          ENDDO
        ENDIF
        CALL EIRENE_MASRR1 (' TARGETS,EI',SFEIT(1),NTARGI,5)
        CALL EIRENE_MASRR1 (' TARGETS,EE',SFEET(1),NTARGI,5)
        DO ITARG=1,NTARGI
          IF (ANY(SFNIT(ITARG,1:NFLA).NE.0.0)) THEN
          DO IFL=1,NFLA
             WRITE(iunout,'(A,I0,A,I0,A,ES12.4)') 'TARGET ', ITARG,
     .                      ', NI(IFL =',IFL,') ',SFNIT(ITARG,IFL)
          ENDDO
          ENDIF
          CALL EIRENE_LEER(1)
        ENDDO
        CALL EIRENE_LEER(1)
        CALL EIRENE_MASR2(' TOTALS, EI,EE  ',SFEIT(0),SFEET(0))
        DO IFL=1,NFLA
           WRITE(iunout,'(A,I0,A,ES12.4)') 'TOTALS, NI(IFL)=',IFL,') ',
     .                                      SFNIT(0,IFL)
        ENDDO

        CALL EIRENE_LEER(2)

        WRITE (iunout,*) ' NEUTRAL PLASMA INTERACTION: '
        CALL EIRENE_MASR2(' SSEI,SSEE      ',SSEI,SSEE)
        DO IFL=1,NFLA
           WRITE(iunout,'(A,I0,A,ES12.4)') 'SSNI(IFL=',IFL,') ',
     .                                      SSNI(IFL)
        ENDDO
        CALL EIRENE_LEER(2)

        WRITE (iunout,*)
     .    ' VOLUMETRIC ENERGY SINKS FOR ELECTRONS, FROM B2 '
        CALL EIRENE_MASR4(' B2BREM,B2RAD,-B2QIE,-B2VDP     ',
     .               B2BREM,B2RAD,-B2QIE,-B2VDP)
        WRITE (iunout,*)
     .    ' TARGET SHEATH CONTRIBUTIONS,ELECTRONS AND IONS '
        CALL EIRENE_MASRR1 (' TARGETS,EI',SHEAI(1),NTARGI,5)
        CALL EIRENE_MASRR1 (' TARGETS,EE',SHEAE(1),NTARGI,5)
        CALL EIRENE_MASR2(' TOTALS,EI,EE    ',SHEAI(0),SHEAE(0))
        CALL EIRENE_LEER(2)

        CALL EIRENE_MASR2(' BALANI,BALANE  ',BALANI,BALANE)
        DO IFL=1,NFLA
           WRITE(iunout,'(A,I0,A,ES12.4)') 'BALANN(IFL)=',IFL,') ',
     .                                      BALANN(IFL)
        ENDDO
        CALL EIRENE_LEER(1)

        CALL EIRENE_MASR2('REL.ERR.(%)RI,RE',RI,RE)
        DO IFL=1,NFLA
           WRITE(iunout,'(A,I0,A,ES12.4)') 'RN(IFL)=',IFL,') ',RN(IFL)
        ENDDO
        CALL EIRENE_LEER(1)
        MINSPEZ=99
        MAXSPEZ=-1
        DO ITARG=1,NTARGI
          DO IPRT=1,NTGPRT(ITARG)
             MINSPEZ=MIN(MINSPEZ,NSPZI(ITARG,IPRT))
             MAXSPEZ=MAX(MAXSPEZ,NSPZE(ITARG,IPRT))
          ENDDO
        ENDDO
        BALAN=0.
        TOT=0.
        DO IFL=MINSPEZ,MAXSPEZ
          BALAN=BALAN+SFNISY(IFL)+SFNINY(IFL)+SFNIWX(IFL)
     .               +SFNIEX(IFL)+SFNIT(0,IFL)+SSNI(IFL)
          TOT=TOT+ABS(SFNISY(IFL)+SFNINY(IFL))+ABS(SFNIT(0,IFL))+
     .            ABS(SSNI(IFL))
        ENDDO
        RNT=BALAN/(TOT+EPS60)*100.
        CALL EIRENE_MASJ2('SUMMED OVER     ',MINSPEZ,MAXSPEZ)
        CALL EIRENE_MASR3('BALAN,TOT,RNT           ',BALAN,TOT,RNT)

        CALL EIRENE_LEER(1)
        WRITE (iunout,*) ' NOISE FROM SOURCE TERMS '

        RESSNI(0,1:NFLA) = RESSNI(0,1:NFLA)/ELCHA
        RESSMO(0,1:NFLA) = RESSMO(0,1:NFLA)/ELCHA
        CALL EIRENE_MASR4(' RESSEE,RESSEI,RESSNI,RESSMO    ',
     .        RESSEE(0),RESSEI(0),SUM(RESSNI(0,1:NFLA)),
     .        SUM(RESSMO(0,1:NFLA)))
        CALL EIRENE_LEER(1)

        WRITE (iunout,*) ' RESSNI-CONTRIBUTIONS BY DIFFERENT SPECIES '
cdr  wrong format in call to masrr1
cdr     CALL EIRENE_MASRR1 (' RESSNI    ',RESSNI(0,1:NFLA),NFLA,5)
        if (.not.allocated(helpw)) allocate (helpw(nfla))
        helpw(1:nfla) = RESSNI(0,1:NFLA)
        CALL EIRENE_MASRR1 (' RESSNI    ',HELPW,NFLA,5)

        WRITE (iunout,*) ' RESSMO-CONTRIBUTIONS BY DIFFERENT SPECIES '
cdr  wrong format in call to masrr1
cdr     CALL EIRENE_MASRR1 (' RESSMO    ',RESSMO(0,1:NFLA),NFLA,5)
        helpw(1:nfla) = RESSMO(0,1:NFLA)
        CALL EIRENE_MASRR1 (' RESSMO    ',HELPW,NFLA,5)
        if (allocated(helpw)) deallocate (helpw)

      ENDIF  !LBALAN
C
      CALL EIRENE_LEER (1)
C
11000 CONTINUE
C
      RETURN
C
 8888 FORMAT (3E14.6)
C
C
      CONTAINS

       FUNCTION EIRENE_POINT_ON_INTERVAL(X,Y,X1,Y1,X2,Y2)
C
C FINDS IF THE POINT X;Y BELONGS TO INTERVAL (X1,Y1)..(X2,Y2)
C
       LOGICAL :: EIRENE_POINT_ON_INTERVAL
       REAL(DP),INTENT(IN) :: X,Y,X1,Y1,X2,Y2
       REAL(DP) :: A,B,C,D,L, EPS4

       EIRENE_POINT_ON_INTERVAL=.FALSE.
       EPS4 = 10._DP * EPS5

       IF(X1.LT.X2.AND.
     f    (X.GT.X2+EPS4.OR.X.LT.X1-EPS4)) RETURN
       IF(X1.GT.X2.AND.
     f    (X.GT.X1+EPS4.OR.X.LT.X2-EPS4)) RETURN
       IF(Y1.LT.Y2.AND.
     f    (Y.GT.Y2+EPS4.OR.Y.LT.Y1-EPS4)) RETURN
       IF(Y1.GT.Y2.AND.
     f    (Y.GT.Y1+EPS4.OR.Y.LT.Y2-EPS4)) RETURN

       A=Y2-Y1
       B=X1-X2
       C=Y1*X2-Y2*X1
       D=A*X+B*Y+C
       D=D*D
       L=A*A+B*B
       IF(L.LT.EPS5) THEN
        WRITE(iunout,*) "WARNING FROM  POINT_ON_INTERVAL"
        WRITE(iunout,*) "THE LENGTH OF THE INTERVAL IS TOO SMALL"
        WRITE(iunout,*) "L,X1,Y1,X2,Y2 ",L,X1,Y1,X2,Y2
        RETURN
       END IF
       IF(D.LT.1.D-8*L) EIRENE_POINT_ON_INTERVAL=.TRUE.

       RETURN

       END FUNCTION EIRENE_POINT_ON_INTERVAL



C DEFINE NORMAL DIRECTION FOR SURFACE-AVERAGED TALLIES (SEE FOLNEUT.F)
       SUBROUTINE CORRECTNSS

         IF(ITRI.GT.NTRIS.OR.NBAR.GT.NTRIS.OR.
     .     NUMSI.GT.3.OR.NBARSI.GT.3)
     .     WRITE(iunout,*) "ERROR IN CORRECTNSS",
     .                  "ITRI,NTRIS,NBAR,NTRIS,NUMSI,NBARSI",
     .                   ITRI,NTRIS,NBAR,NTRIS,NUMSI,NBARSI

         INMTINSS(NUMSI,ITRI)=1
         INMTINSS(NBARSI,NBAR)=1

         IF ((IXTRI(ITRI) > 0) .AND. (IYTRI(ITRI) > 0) .AND.
     .       (IXTRI(NBAR) > 0) .AND. (IYTRI(NBAR) > 0)) THEN
!  both triangles inside mesh
            IF (IXTRI(ITRI) == IXTRI(NBAR)) THEN
              IF (IYTRI(ITRI) > IYTRI(NBAR)) THEN
                INMTINSS(NUMSI,ITRI) = -1
              ELSE
                INMTINSS(NBARSI,NBAR) = -1
              END IF

            ELSE IF (IYTRI(ITRI) == IYTRI(NBAR)) THEN
              IF (IXTRI(ITRI) > IXTRI(NBAR)) THEN
                INMTINSS(NUMSI,ITRI) = -1
              ELSE
                INMTINSS(NBARSI,NBAR) = -1
              END IF
            END IF


!  triangle IT inside mesh, triangle NBAR outside mesh
          ELSE IF ((IXTRI(ITRI) > 0) .AND. (IYTRI(ITRI) > 0)) THEN

!  poloidal surface
            IF (LXSRF) THEN
              IF (IXTRI(ITRI).EQ.1) THEN        ! 'W SURFACE'
                INMTINSS(NUMSI,ITRI ) = -1
              ELSE                            ! 'E SURFACE'
                INMTINSS(NBARSI,NBAR) = -1
              END IF

!  radial surface
            ELSE
              IF (IYTRI(ITRI).EQ.1) THEN        ! 'S SURFACE'
                INMTINSS(NUMSI,ITRI ) = -1
              ELSE                            ! 'N SURFACE'
                INMTINSS(NBARSI,NBAR) = -1
              END IF

            END IF

!  triangle NBAR inside mesh, triangle IT outside mesh
          ELSEIF((IXTRI(NBAR) > 0) .AND. (IYTRI(NBAR) > 0)) THEN

!  poloidal surface
            IF (LXSRF) THEN
              IF (IXTRI(NBAR).EQ.1) THEN      ! 'W SURFACE'
                INMTINSS(NBARSI,NBAR ) = -1
              ELSE                            ! 'E SURFACE'
                INMTINSS(NUMSI,ITRI) = -1
              END IF

!  radial surface
            ELSE
              IF (IYTRI(NBAR).EQ.1) THEN      ! 'S SURFACE'
                INMTINSS(NBARSI,NBAR ) = -1
              ELSE                            ! 'N SURFACE'
                INMTINSS(NUMSI,ITRI) = -1
              END IF

            END IF

!  both triangles outside mesh --> do not know - do nothing
          END IF

       END  SUBROUTINE CORRECTNSS
CVK END

      END

C> \brief Any property requiring hand-over in parallel part.
C>
C> This interfacing routine is called in the parallel part of EIRENE
C> after the broadcast of any other quantity and before MCARLO.
      SUBROUTINE EIRENE_INFCOP_PRE_MCARLO
      RETURN
      END SUBROUTINE EIRENE_INFCOP_PRE_MCARLO

C> \brief Return data at the end of an EIRENE stratum
C>
C> At this entry data are transferred back from EIRENE to the
C> interfacing module (and from there, after possibly further
C> processing, to the external code). This entry is called from the
C> "strata-loop" in subr. MCARLO, after all trajectories for a
C> particular stratum ISTRA have been sampled and after all volume and
C> surface tallies have been scaled and processed to their final form.
C>
C> The call to IF3COP is controlled by the flag NMODE (input block 1).
C> At call data are expected for stratum ISTRA. There they may be
C> further prepared (e.g. normalized, or scaled to other units) for
C> transfer to the external code
      SUBROUTINE EIRENE_INFCOP_POST_STRATUM(ISTRA)
      integer, intent(in) :: istra
      RETURN
      END SUBROUTINE EIRENE_INFCOP_POST_STRATUM
