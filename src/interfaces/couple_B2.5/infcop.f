cdr Nov 2017: clean up, comments. Sync. between couple_B2.5, vs. 2008, 
cdr           with couple_B2, Nov. 17 git master.
c             1) remove FNIX_YB, FNIY_XB 
c             2) add input: mshfrm
c             3) remove printout of sputter fluxes, after end of global balances
c             4) additional species index (ipls) in eapl,empl,eipl
c             5) LLCUT included
c             6) BFINTF set  (as in couple_B2, still not from B2.5 plasma files)
c             7) if not lshort: call eirene_save_tallies
c             8) sigma_cop removed. Ressni, ressee,..., now: per stratum. Allocatable
c             9) clarify: dimensions (allocations) in braeir: b2 vs. b2.5 ??
c            10) add. elstep, eemax,esheath in step fct. NEMODS=2,3, RATHER 8,9
c            11) EPEL --> EPEL_COP  (also in couple_B2)
c             CPPV --> MPPL_COP
c             ELTEST, EMAXW,... for a target energ y flux as interpreted from B2 output. 
c            12) bug fix re vol.rec., only one ipls per stratum is supported
c                code was correct in solps4.3, and garching versions of couple_b2/b2.5
c            
c

cdr Nov. 17: removed dead option LINDIM: here and in couple_b2_parallel
cdr Dec. 17:
c    trcsou --> trcint: consistency checks for target recycling step functions
c               trcsou: print step functions from samsrf.f, as finally used in eirene.
C
C   EIRENE CODE SEGMENT COUPLE_$, $ MAY CURRENTLY STAND FOR B2,
C                                                           B2.5,
C                                                           DIVIMP,
C                                                           TRIA,
C                                                           TETRA,
C                                                           TRANSP,
C                                                           DUMMY
C
C   THIS VERSION: $COUPLE_B2.5  JAN. 2018
c                 proprietary version of FZJ, for local B2.5 code versions   
C
C   UPDATES:
C   OPTION TO EVALUATE B-FIELD VECTORS FROM GRIDADAP FILE FT29
C   FOR NON-ORTHOGONAL GRIDS
C
C   THIS CODE SEGMENT CONTAINES VARIOUS SUBROUTINES NEEDED FOR
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
C   THIS PARTICULAR VERSION LINKS EIRENE TO THE B2.5 2D MULTIFLUID EDGE
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

cdr for species dependent global particle balance
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
     .          SFEISY, SFEESY, RECADD, RECTOT,
     .          EEADD, EIADD, PIADD, 
     .          SMOCL, CHEES, CHEIS, SNICL,
     .          SIGNUM, 
     .          SSE, BALANI, BALANE, SSEE, SSI, RE, RI, RNT, TOT,
     .          TOTI, TOTE, BALAN, RRBC,
     .          SSEI, SFEIEX, SFEEEX, SFEENY, SFEIWX, VVBC,
     .          UUBC, UPBC, RBC, UDBC, VL, V, T, 
     .          BX, BY, BZ, BN, TEST,
     .          DELTE_PARA, DELTI_PARA, DELTE_PERP, DELTI_PERP,
     .          TES,TIS,
     .          DELY, ALX, ALE, ALW, ALS, ALN, AL, ETOT,
     .          FLX, ESUM, DR, VR, VTEST, VTEST2, EADD, 
     .          EMAXW, ESHEATH, SI,
     .          PARWI, PERWI, SUMM, SUMN, SUMEI, SUMEE, FLXI,
     .          CHI,  CHP, CHE, CS, THMAX, EESHT, EEMAX,
     .          RP1, DELX, PNORM, PVYS, PVXS, PUPV, RRBS, PUYS, PUXS,
     .          VPX, VPY, VT, PARW, PERW, PN1, OR, VPZ, GAMMA, CUR, TE,
     .          SFEEWX, SFEINY, PM1, DRR, VDBC     

      INTEGER, SAVE :: J, IRC, JC, INC, IADD, 
     .           IP, ITARG, IO, IFL, NPES,
     .           IIPLS, IG, IGITT, IEPLS, NPEC, NPBC, NPBS, NTGPRI,
     .           IT, I, IPRT, IAOT, IAIN, IREAD, IPL, INN,
     .           IMODE, IERROR, LTARG, IN, IX, IY,
     .           NPLP, NDX2, NRED, IO29, NDXY, IFIRST,
     .           ISTRAI, IRRC, K, IR, IIRC, ICPV, I34,
     .           NREC11, NEM, ISTR,
     .           IXI, IXE, IPLSTI, IPLSV, IPLV, ISP,
     .           mshfrm, imf, ixm1,iym1,
     .           MINSPEZ, MAXSPEZ

      INTEGER, INTENT(IN) :: ISTRAA, ISTRAE, NEW_ITER, IFRST, ITRG
      REAL(DP) :: EIRENE_STEP, EIRENE_FTABRC1, EIRENE_FEELRC1, 
     .            EIRENE_SHEATH, EIRENE_EMAXW
      INTEGER, EXTERNAL :: EIRENE_IDEZ
C
      LOGICAL, INTENT(INOUT) :: LSTP
      LOGICAL, SAVE :: LSHORT, LSTOP, LTEST, LSTP3,
     .                 LNONREC_SY,LNONREC_NY,LNONREC_WX,LNONREC_EX
     .                ,LCOARSE

      LOGICAL, ALLOCATABLE, SAVE :: LLCUT(:)

      REAL(DP), ALLOCATABLE, SAVE ::
     . CHPS(:),    SNIS(:),    CHMOS(:),  SMOS(:),  SCALN(:),
     . SNIS0(:,:), SMOS0(:,:),
c
     . RESSNI(:,:),  RESSMO(:,:), 
     . RESSEE(:), RESSEI(:)
     .,FLXEIR(:)

      REAL(DP), ALLOCATABLE, SAVE ::
     . TORL(:,:), ESHT(:,:), ELTEST(:,:), ORI(:,:)

      real(dp),allocatable :: helpw(:)

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
C
C
      ENTRY EIRENE_IF0COP
C
      LSHORT=.FALSE.
C
      GOTO 99990
C
C  TO INITIALISE THE SHORT CYCLING, THE GEOMETRY HAS TO BE
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
cdr
cdr   lchkqud = .false.  !  only needed for triangular grid options
cdr
      mshfrm = 0   !  optional flag for geometry file format: linda, carree, sonnet
      NLSHRT13 = .TRUE.  !  only short version of fort13 is used: calls WRPLAM_SHRT, RPLAM_SHRT
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

        READ (IUNIN,'(5I6)') NFLA,NCUTB,NCUTL,IMF
cdr  imf  flag for different formats of geometry file: linda, sonnet, carree. What is What?
        IF (IMF /= 0) MSHFRM=IMF
        NCUTB_SAVE=NCUTB
        IF (TRCINT) THEN
          WRITE (iunout,*) ' NFLA,NCUTB,NCUTL,IMF = ',
     .                       NFLA,NCUTB,NCUTL,IMF
          WRITE (iunout,*) ' IPLS,IFLB(IPLS),FCTE(IPLS),BMASS(IPLS)'
        ENDIF
        DO 20 IPL=1,NPLSI
          READ (IUNIN,'(2I6,2E12.4)') I,IFLB(IPL),FCTE(IPL),BMASS(IPL)
          IF (TRCINT)
     .    WRITE (iunout,*)          IPL,IFLB(IPL),FCTE(IPL),BMASS(IPL)
20      CONTINUE
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
22      CONTINUE
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
33        CONTINUE
          IF (TRCINT) CALL EIRENE_LEER(1)
30      CONTINUE
        READ (IUNIN,'(6E12.4)')  CHGP,CHGEE,CHGEI,CHGMOM
        IF (TRCINT) CALL EIRENE_MASR4
     .                         ('CHGP,CHGEE,CHGEI,CHGMOM         ',
     .                           CHGP,CHGEE,CHGEI,CHGMOM)
C  READ ADDITIONAL DATA TO BE TRANSFERRED FROM B2.5 INTO EIRENE
C  HERE: B2.5 VOLUME TALLIES
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
40      CONTINUE
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
50      CONTINUE
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
      NDX = NDXA
      NDY = NDYA
      NFL = NFLA
      NDXP = NDX+1
      NDYP = NDY+1
C
      LNLPLG=NLPLG
      LNLDRF=NLDRFT
      LTRCFL=TRCFLE
      NSTRI=NSTRAI
      DO 60 ISTR=1,NSTRAI
        LNLVOL(ISTR)=NLVOL(ISTR)
60    CONTINUE
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
     .  'ENERGY WEIGHTED CX RATE OF ATOMS WITH IPLS                  '
        TXTSPC(IPLS,NTALM)=TEXTS(NSPAMI+IPLS)
        TXTUNT(IPLS,NTALM)='AMP                     '
      ENDDO
C
70    CONTINUE
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
1000  CONTINUE
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
1020    CONTINUE
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
! cell centered angle of B_pol (psi-contour line) against eirene x-coordinate
            ALX=0.25D0*(ALE+ALW+ALN+ALS)
            PUXE(IN)=COS(ALE)
            PUYE(IN)=SIN(ALE)
            PUXN(IN)=COS(ALN)
            PUYN(IN)=SIN(ALN)
! cell centered unit vector along poloidal direcion
            PUX(IN)=COS(ALX)
            PUY(IN)=SIN(ALX)
! cell centered unit vector along "radial" (grad psi) direcion,
!                    strictly orthonormal to  PU (poloidal) direction
            PVX(IN)=-PUY(IN)
            PVY(IN)=PUX(IN)
          END DO
        END DO
C
        DO IY=1,NDYA
          DO IX =1,NDXA
            IN=IY+(IX-1)*NR1ST
            NSTGRD(IN)=ABS(IAISO(IX,IY)-1)
          END DO
        END DO

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
cdr
cdr  check with if0prm, values of ndxp,.....
cdr  for b2:       0 ...NDXP ?
c     CALL EIRENE_ALLOC_BRAEIR(NDXP,NDYP,NFL,IFOFF)
cdr  for b2.5:    -1 ...NDX  ?
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
      CALL EIRENE_PLASM (31,NDX2,NDYA,NFLA,NDX,NDY,NFL,UUDIAB)
      CALL EIRENE_PLASM (31,NDX2,NDYA,NFLA,NDX,NDY,NFL,VVDIAB)
      CALL EIRENE_PLASM (31,NDX2,NDYA,1,NDX,NDY,1,POB)
C
C  OPTIONAL ARRAYS: VOLB, BFELDB
      VOLB = 0.D0
      BFELDB = 0.D0

C  CELL VOLUMES AS USED IN B2
      CALL EIRENE_PLASM (31,NDX2,NDYA,1,NDX,NDY,1,VOLB)
C  MAGNETIC FIELD STRENGTH (TESLA)
      CALL EIRENE_PLASM (31,NDX2,NDYA,1,NDX,NDY,1,BFELDB)
C
2100  CONTINUE
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
C  FIRST THE ZONE CENTERED DATA
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
C  NOW THE SURFACE CENTERED DATA
      CALL EIRENE_INDMAP (FNIXB,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (FNIYB,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .             NCUTB,NCUTL,NPOINT,NPLP)
C  distinct from B2: these velocities are cell centered in b2.5
      CALL EIRENE_INDMAP (UUB,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (VVB,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (UPB,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .             NCUTB,NCUTL,NPOINT,NPLP)
C   same in B2 and in B2.5:  these ENERGY fluxes are surface centered
      CALL EIRENE_INDMAP (FEIXB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (FEIYB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (FEEXB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (FEEYB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,
     .             NCUTB,NCUTL,NPOINT,NPLP)
c  b2.5 only: additional velocities from plasma drifts, cell centered
      CALL EIRENE_INDMAP (UUDIAB,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .             NCUTB,NCUTL,NPOInt,NPLP)
      CALL EIRENE_INDMAP (VVDIAB,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .             NCUTB,NCUTL,NPOINT,NPLP)
C  UNUSED INPUT TALLIES
c  presumably:  POB  (electric potential ??) and BFELDB  are cell centered
      CALL EIRENE_INDMAP (POB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (VOLB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (BFELDB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,
     .             NCUTB,NCUTL,NPOINT,NPLP)
C
2101  CONTINUE
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
      DO 2105 IPLS=1,NPLSI
        D(IPLS)=1.D-6*FCTE(IPLS)         !pb 1.e-6 -> 1.d-6
        FL(IPLS)=ELCHA*FCTE(IPLS)        !  SCALING FOR FLUX: 1/S --> AMP
2105  CONTINUE
C
      BZINTF = 1._DP
      BFINTF = 1._DP
      DO 2110 IY=1,NDYA
        DO 2120 IX =1,NDXA
          IN=IY+(IX-1)*NR1ST
          TEINTF(IN)=TEB(IX,IY)*T
C
C  ONLY ONE ION TEMPERATURE AVAILABLE FROM PLASMA FLUID CODE,
C  SEE LOOP 2150 BELOW
          TIINTF(1,IN)=TIB(IX,IY)*T
C
C  polodial field
          BX=PUX(IN)*RRB(IX,IY)   ! +PVX(IN)*0., but radial field is zero
          BY=PUY(IN)*RRB(IX,IY)   ! +PVY(IN)*0.
c  toroidal field
          BZ=SQRT(1.-RRB(IX,IY)**2)
c  normalize B-field vector to length 1 (one)
c
          BN=SQRT(BX*BX+BY*BY+BZ*BZ)
          BXINTF(IN)=BX/BN
          BYINTF(IN)=BY/BN
          BZINTF(IN)=BZ/BN
          BFINTF(IN)=BN
          VLINTF(IN)=VOLB(IX,IY)*VL
2120    CONTINUE
2110  CONTINUE
C
C  SET SAME ION TEMPERATURE FOR ALL EIRENE BACKGROUND SPECIES
C
      DO 2150 IPLS=1,NPLSI
      IPLSTI=MPLSTI(IPLS)
      DO 2150 IY=1,NDYA
        DO 2150 IX =1,NDXA
          IN=IY+(IX-1)*NR1ST
          TIINTF(IPLSTI,IN)=TIINTF(1,IN)
2150  CONTINUE
C
CDR  set density from B2 array DNIB, for each fluid
CDR  set plasma flow velocity field from B2 arrays UPB (parallel velocity)
c  without drifts:
c  upb * pitch:  poloidal velocity (i.e. cartesian x,y direction).
c  poloidal field direction is given by that of the poloidal cell face PU..(in),
C  i.e. along a flux surface. (PU(...) is cell centered)
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
            DO 2210 IY=1,NDYA
              IYM1=IY-1
              DO 2220 IX = 1,NDXA
                IN=IY+(IX-1)*NR1ST
                DIINTF(IPLS,IN)=DNIB(IX,IY,IFL)*D(IPLS)
c
                UPBC=UPB(IX,IY,IFL)
                UDBC=UUDIAB(IX,IY,IFL)
cdr  radial velocity
                VVBC=0.5*(VVB(IX,IYM1,IFL)+VVB(IX,IY,IFL))

cdr  rrb: pitch  B_pol/B_total
C  UDBC: DIAMAGNETIC VELOCITY, SHOULD BE ZERO IN B2, AND NONZERO IN EB2
C        NOTE: IN LINEAR DEVICES: RRB=1
                VDBC=VVDIAB(IX,IY,IFL)
                RRBC=RRB(IX,IY)
cdr  now set cartesian flow velocity components
                VXINTF(IPLSV,IN)=
     &           (PUX(IN)*(UPBC*RRBC-UDBC*SQRT(1.-RRBC**2))+
     &            PVX(IN)*VVBC)*V
                VYINTF(IPLSV,IN)=
     &           (PUY(IN)*(UPBC*RRBC-UDBC*SQRT(1.-RRBC**2))+
     &            PVY(IN)*VVBC)*V
                VZINTF(IPLSV,IN)=(UPBC*SQRT(1.-RRBC**2)+UDBC*RRBC)*V
2220          CONTINUE
2210        CONTINUE
C  EIRENE BACKGROUND SPECIES "IPLS" IS NOW FILLED WITH B2 DATA "IFL"
2201      CONTINUE   !IFL loop

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
            IPLSV  = MPLSV(IPLS)
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
2200  CONTINUE
C  B2.5-BRAAMS CODE SPECIFIC END
C
C
C  READ OTHER B2.5 ARRAYS INTO EIRENE, FOR PRINTOUT AND PLOTTING
C
c  density, species index as in B2 code, cell centered
      DO 2300 IAIN=1,NAINB
        IF (NAINT(IAIN).EQ.1.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO 2321 IY=1,NDYA
          DO 2321 IX=1,NDXA
            IN=IY+(IX-1)*NR1ST
            ADINTF(IAIN,IN)=DNIB(IX,IY,NAINS(IAIN))
2321      CONTINUE
c  poloidal (projection) flow velocity, species index as in B2 code, north surface centered
        ELSEIF (NAINT(IAIN).EQ.2.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO 2322 IY=1,NDYA
          DO 2322 IX=1,NDXA
            IN=IY+(IX-1)*NR1ST
            ADINTF(IAIN,IN)=UUB(IX,IY,NAINS(IAIN))
2322      CONTINUE
c  radial drift velocity, species index as in B2 code, east surface centered
        ELSEIF (NAINT(IAIN).EQ.3.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO 2323 IY=1,NDYA
          DO 2323 IX=1,NDXA
            IN=IY+(IX-1)*NR1ST
            ADINTF(IAIN,IN)=VVB(IX,IY,NAINS(IAIN))
2323      CONTINUE
c  plasma pressure, cell centered, no species index
        ELSEIF (NAINT(IAIN).EQ.6) THEN
          DO 2326 IY=1,NDYA
          DO 2326 IX=1,NDXA
            IN=IY+(IX-1)*NR1ST
            ADINTF(IAIN,IN)=PRB(IX,IY)
2326      CONTINUE
c  parallel velocity, species index as in B2 code, north surface centered
        ELSEIF (NAINT(IAIN).EQ.7.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO 2327 IY=1,NDYA
          DO 2327 IX=1,NDXA
            IN=IY+(IX-1)*NR1ST
            ADINTF(IAIN,IN)=UPB(IX,IY,NAINS(IAIN))
2327      CONTINUE
c  pitch angle, no species index
        ELSEIF (NAINT(IAIN).EQ.8) THEN
          DO 2328 IY=1,NDYA
          DO 2328 IX=1,NDXA
            IN=IY+(IX-1)*NR1ST
            ADINTF(IAIN,IN)=RRB(IX,IY)
2328      CONTINUE
        ELSEIF (NAINT(IAIN).EQ.9.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO 2329 IY=1,NDYA
          DO 2329 IX=1,NDXA
            IN=IY+(IX-1)*NR1ST
            ADINTF(IAIN,IN)=FNIXB(IX,IY,NAINS(IAIN))
2329      CONTINUE
        ELSEIF (NAINT(IAIN).EQ.10.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO 2330 IY=1,NDYA
          DO 2330 IX=1,NDXA
            IN=IY+(IX-1)*NR1ST
            ADINTF(IAIN,IN)=FNIYB(IX,IY,NAINS(IAIN))
2330      CONTINUE
        ELSEIF (NAINT(IAIN).EQ.11) THEN
          DO 2331 IY=1,NDYA
          DO 2331 IX=1,NDXA
            IN=IY+(IX-1)*NR1ST
            ADINTF(IAIN,IN)=FEIXB(IX,IY)
2331      CONTINUE
        ELSEIF (NAINT(IAIN).EQ.12) THEN
          DO 2332 IY=1,NDYA
          DO 2332 IX=1,NDXA
            IN=IY+(IX-1)*NR1ST
            ADINTF(IAIN,IN)=FEIYB(IX,IY)
2332      CONTINUE
        ELSEIF (NAINT(IAIN).EQ.13) THEN
          DO 2333 IY=1,NDYA
          DO 2333 IX=1,NDXA
            IN=IY+(IX-1)*NR1ST
            ADINTF(IAIN,IN)=FEEXB(IX,IY)
2333      CONTINUE
        ELSEIF (NAINT(IAIN).EQ.14) THEN
          DO 2334 IY=1,NDYA
          DO 2334 IX=1,NDXA
            IN=IY+(IX-1)*NR1ST
            ADINTF(IAIN,IN)=FEEYB(IX,IY)
2334      CONTINUE
        ELSEIF (NAINT(IAIN).EQ.15.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO IY=1,NDYA
          DO IX=1,NDXA
            IN=IY+(IX-1)*NR1ST
            ADINTF(IAIN,IN)=UUDIAB(IX,IY,NAINS(IAIN))
          ENDDO
          ENDDO
        ELSEIF (NAINT(IAIN).EQ.16.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO IY=1,NDYA
          DO IX=1,NDXA
            IN=IY+(IX-1)*NR1ST
            ADINTF(IAIN,IN)=VVDIAB(IX,IY,NAINS(IAIN))
          ENDDO
          ENDDO
c   cell volume as in b2 code, no species index (cell centered)
        ELSEIF (NAINT(IAIN).EQ.17) THEN
          DO 2335 IY=1,NDYA
          DO 2335 IX=1,NDXA
            IN=IY+(IX-1)*NR1ST
            ADINTF(IAIN,IN)=VOLB(IX,IY)
2335      CONTINUE
cdr  magnetic field strength, Tesla
        ELSEIF (NAINT(IAIN).EQ.18) THEN
          DO 2336 IY=1,NDYA
          DO 2336 IX=1,NDXA
            IN=IY+(IX-1)*NR1ST
            ADINTF(IAIN,IN)=BFELDB(IX,IY)
2336      CONTINUE


cdr  free: NAINT=20 --29:  reserved for AMDIAG:  scaled atomic/molecular rate coefficients
cdr                        evaluated on computational grid. See Manual.

        ENDIF
2300  CONTINUE
C
      RETURN
C
2999  CONTINUE
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
3000  CONTINUE
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
3005  CONTINUE
C
      ALLOCATE (TORL(NSTRA,NGITT))
      ALLOCATE (ESHT(NSTEP,NGITT))
      ALLOCATE (ELTEST(NSTEP,NGITT))
      ALLOCATE (ORI(NSTEP,NGITT))

      RRSTEP(ITARG,1)=0.
      IG=0
      IIPLS=NPLSI
      IEPLS=1
      DO 3040 IPRT=1,NTGPRT(ITARG)
C  NINCT= 1: PLASMA FLUX IN SAME   DIRECTION AS B2 CO-ORDINATE
C  NINCT=-1: PLASMA FLUX IN OPPOS. DIRECTION AS B2 CO-ORDINATE
C  BRAAMS X-CELL CONTAINING THE TARGET DATA (BOUNDARY CONDITIONS)
C  (SURFACE CENTERED, EAST OR NORTH) (AFTER INDEX MAPPING)
C  (E.G. SURFACE NO.0 AND SURFACE NO. NX) AT TARGETS.
        NPBS=NDT(ITARG,IPRT)
C  BRAAMS P-CELL CONTAINING THE TARGET DATA (BOUNDARY CONDITIONS)
C  (ZONE CENTERED) (AFTER INDEX MAPPING)
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
        DO IY=NTIN(ITARG,IPRT),NTEN(ITARG,IPRT)-1
          IG=IG+1
          IF (IG.GT.NGITT) GOTO 999
C  TESTEP, TISTEP: ZONE CENTERED TEMPERATURE IN BOUNDARY ZONE (EV)
          ORI(ITARG,IG) = NINCT(ITARG,IPRT)
          TESTEP(ITARG,IG) = TEB(NPBC,IY)*T
C  RRSTEP,IRSTEP,IPSTEP: GEOMETRICAL INFORMATION ALONG TARGET
C  RRSTEP IS THE ARC LENGTH ALONG THE TARGET (CM)
          RRSTEP(ITARG,IG+1)=RRSTEP(ITARG,IG) +
     .                       SQRT((XPOL(IY+1,NPES)-XPOL(IY,NPES))**2+
     .                            (YPOL(IY+1,NPES)-YPOL(IY,NPES))**2)
C  EIRENE CELL NUMBER INFORMATION ALONG TARGET
          IRSTEP(ITARG,IG)=IY
          IPSTEP(ITARG,IG)=NPEC
          ITSTEP(ITARG,IG)=1
          IASTEP(ITARG,IG)=0
          IBSTEP(ITARG,IG)=1
          IGSTEP(ITARG,IG)=200000+NPES
C  TORL: TOROIDAL LENGTH (CM) AT TARGET SEGMENT IY: CENTER OF GRAVITY
          TORL(ITARG,IG)=2.*PIA*0.5*(XPOL(IY+1,NPES)+
     .                               XPOL(IY,NPES))
          DO 3013 IPLS=1,NPLSI
            IPLSTI=MPLSTI(IPLS)
            IPLSV=MPLSV(IPLS)
            ELSTEP(IPLS,ITARG,IG)=0.
            TISTEP(IPLSTI,ITARG,IG) = TIB(NPBC,IY)*T
C  DISTEP: ZONE CENTERED DENSITY IN BOUNDARY ZONE
            IFL=IFLB(IPLS)
            IF (IFL.LE.0.OR.IFL.GT.NFLA) GOTO 3013
            DISTEP(IPLS,ITARG,IG)=DNIB(NPBC,IY,IFL)*D(IPLS)
C  FLSTEP: SURFACE CENTERED FLUX (AMP/CM ALONG TARGET)
            IF (NSPZI(ITARG,IPRT).LE.IFL.AND.
     .                               IFL.LE.NSPZE(ITARG,IPRT)) THEN
              IIPLS=MIN0(IIPLS,IPLS)
              IEPLS=MAX0(IEPLS,IPLS)
!  LENGTH OF CELL FACE OF B2.5 CELL
              DELY=RRSTEP(ITARG,IG+1)-RRSTEP(ITARG,IG)
              FLSTEP(IPLS,ITARG,IG)=0.
              IF (DELY.GT.0.) THEN
                FLSTEP(IPLS,ITARG,IG)=MAX(0._DP,ORI(ITARG,IG)*
     .                                FNIXB(NPBS,IY,IFL))*FL(IPLS)/DELY

C  SET DEFAULT ION ENERGY FLUXES FROM B2.5 BOUNDARY CONDITIONS
                delti_para=3
                delte_para=0.5
                delti_perp=2
                delte_perp=0
                tis=TISTEP(IPLSTI,ITARG,IG)
                tes=TESTEP(ITARG,IG)
                ELSTEP(IPLS,ITARG,IG) = ELSTEP(IPLS,ITARG,IG) +
     .                                  FL(IPLS)/DELY*
     .            (TIS*delti_para+TES*delte_para)
     .             *ABS(fnixb(npbs,iy,ifl))
              ENDIF
            ENDIF

C  VXSTEP,VYSTEP,VZSTEP: SURFACE CENTERED FLOW VELOCITY (CM/S)
C  NOTE: PV VECTOR IS CELL CENTERED, BUT EXACT VECTOR CAN BE FOUND FROM
C        DATA FOR POLOIDAL POLYGON NPES
            IN=IY+(NPEC-1)*NR1ST
            PVXS=XPOL(IY+1,NPES)-XPOL(IY,NPES)
            PVYS=YPOL(IY+1,NPES)-YPOL(IY,NPES)
            PNORM=SQRT(PVXS**2+PVYS**2)
            PVXS=PVXS/(PNORM+EPS60)
            PVYS=PVYS/(PNORM+EPS60)
C  ORTHONORMALISE PU VECTOR WITH RESPECT TO PV
            PUPV=PUX(IN)*PVXS+PUY(IN)*PVYS
            PUXS=PUX(IN)-PUPV*PVXS
            PUYS=PUY(IN)-PUPV*PVYS
            PNORM=SQRT(PUXS**2+PUYS**2)
            PUXS=PUXS/(PNORM+EPS60)
            PUYS=PUYS/(PNORM+EPS60)
! Detlev und Xavier
            VXSTEP(IPLSV,ITARG,IG)=
     .            (PUXS*UUB(NPBS,IY,IFL)+PVXS*VVB(NPBS,IY,IFL))*V
            VYSTEP(IPLSV,ITARG,IG)=
     .            (PUYS*UUB(NPBS,IY,IFL)+PVYS*VVB(NPBS,IY,IFL))*V
            RRBS=0.5*(RRB(NPBC,IY)+
     .                RRB(NPBC-NINCT(ITARG,IPRT),IY))
            VZSTEP(IPLSV,ITARG,IG)=
     .            (SQRT(1.-RRBS**2)*UPB(NPBC,IY,IFL)+
     .             RRBS*UUDIAB(NPBC,IY,IFL))*V
!
3013      CONTINUE
        ENDDO
C
        GOTO 3030
C
3020    CONTINUE
C
C  SECOND: SOURCES AT RADIAL (X) SURFACES
C
        DO IX=NTIN(ITARG,IPRT),NTEN(ITARG,IPRT)-1
          IG=IG+1
          IF (IG.GT.NGITT) GOTO 999
C  TESTEP, TISTEP: ZONE CENTERED TEMPERATURE IN BOUNDARY ZONE (EV)
          ORI(ITARG,IG) = NINCT(ITARG,IPRT)
          TESTEP(ITARG,IG) = TEB(IX,NPBC)*T
C  RRSTEP,IRSTEP,IPSTEP: GEOMETRICAL INFORMATION ALONG TARGET
C  EIRENE CELL NUMBER INFORMATION ALONG TARGET
          RRSTEP(ITARG,IG+1)=RRSTEP(ITARG,IG) +
     .                     SQRT((XPOL(NPES,IX+1)-XPOL(NPES,IX))**2+
     .                          (YPOL(NPES,IX+1)-YPOL(NPES,IX))**2)
          IRSTEP(ITARG,IG)=NPEC
          IPSTEP(ITARG,IG)=IX
          ITSTEP(ITARG,IG)=1
          IASTEP(ITARG,IG)=0
          IBSTEP(ITARG,IG)=1
          IGSTEP(ITARG,IG)=100000+NPES
C  TORL: TOROIDAL LENGTH (CM) AT TARGET SEGMENT IY: CENTER OF GRAVITY
          TORL(ITARG,IG)=2.*PIA*0.5*(XPOL(NPES,IX+1)+
     .                               XPOL(NPES,IX))
          DO 3023 IPLS=1,NPLSI
            IPLSTI=MPLSTI(IPLS)
            IPLSV=MPLSV(IPLS)
            ELSTEP(IPLS,ITARG,IG)=0.
            TISTEP(IPLSTI,ITARG,IG) = TIB(IX,NPBC)*T
C  DISTEP: ZONE CENTERED DENSITY IN BOUNDARY ZONE (EV)
            IFL=IFLB(IPLS)
            IF (IFL.LE.0.OR.IFL.GT.NFLA) GOTO 3023
            DISTEP(IPLS,ITARG,IG)=DNIB(IX,NPBC,IFL)*D(IPLS)
C  FLSTEP: SURFACE CENTERED FLUX (AMP/CM ALONG TARGET)
            IF (NSPZI(ITARG,IPRT).LE.IFL.AND.
     .                               IFL.LE.NSPZE(ITARG,IPRT)) THEN
              IIPLS=MIN0(IIPLS,IPLS)
              IEPLS=MAX0(IEPLS,IPLS)
!  LENGTH OF CELL FACE OF B2 CELL
              DELX=RRSTEP(ITARG,IG+1)-RRSTEP(ITARG,IG)
              FLSTEP(IPLS,ITARG,IG)=0.
              IF (DELX.GT.0.) THEN
                FLSTEP(IPLS,ITARG,IG)=MAX(0._DP,ORI(ITARG,IG)*
     .                                FNIYB(IX,NPBS,IFL))*FL(IPLS)/DELX
C  SET DEFAULT ION ENERGY FLUXES FROM B2.5 BOUNDARY CONDITIONS
                delti_para=3
                delte_para=0.5
                delti_perp=2
                delte_perp=0
                tis=TISTEP(IPLSTI,ITARG,IG)
                tes=TESTEP(ITARG,IG)
                ELSTEP(IPLS,ITARG,IG) = ELSTEP(IPLS,ITARG,IG) +
     .                                  FL(IPLS)/DELX*
     .            (TIS*delti_perp*ABS(Fniyb(ix,npbs,ifl)))
              ENDIF
            ENDIF
C
C  VXSTEP,VYSTEP,VZSTEP: SURFACE CENTERED FLOW VELOCITY (CM/S)
C  NOTE: PU VECTOR IS CELL CENTERED, BUT EXACT VECTOR CAN BE FOUND FROM
C        RADIAL POLYGON NPES DATA
            IN=NPEC+(IX-1)*NR1ST
            PUXS=XPOL(NPES,IX+1)-XPOL(NPES,IX)
            PUYS=YPOL(NPES,IX+1)-YPOL(NPES,IX)
            PNORM=SQRT(PUXS**2+PUYS**2)
            PUXS=PUXS/(PNORM+EPS60)
            PUYS=PUYS/(PNORM+EPS60)
C  ORTHONORMALISE PV VECTOR WITH RESPECT TO PU
            PUPV=PUXS*PVX(IN)+PUYS*PVY(IN)
            PVXS=PVX(IN)-PUPV*PUXS
            PVYS=PVY(IN)-PUPV*PUYS
            PNORM=SQRT(PVXS**2+PVYS**2)
            PVXS=PVXS/(PNORM+EPS60)
            PVYS=PVYS/(PNORM+EPS60)
! Detlev und Xavier
            VXSTEP(IPLSV,ITARG,IG)=
     .            (PUXS*UUB(IX,NPBS,IFL)+PVXS*VVB(IX,NPBS,IFL))*V
            VYSTEP(IPLSV,ITARG,IG)=
     .            (PUYS*UUB(IX,NPBS,IFL)+PVYS*VVB(IX,NPBS,IFL))*V
            RRBS=0.5*(RRB(IX,NPBC)+
     .                RRB(IX,NPBC-NINCT(ITARG,IPRT)))
cxpb            RRBS=UUB(IX,NPBS,IFL)/(UPB(IX,NPBS,IFL)+EPS60)
            VZSTEP(IPLSV,ITARG,IG)=
     .            (SQRT(1.-RRBS**2)*UPB(IX,NPBC,IFL)+
     .             RRBS*UUDIAB(IX,NPBC,IFL))*V
3023      CONTINUE
        ENDDO
3030    CONTINUE
C
3040  CONTINUE
      NRWL(ITARG)=IG+1

C
      IF (TRCINT) CALL EIRENE_LEER(2)
C
C  INITIALISE FUNCTION STEP (FOR RANDOM SAMPLING ALONG TARGET)
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
      INDIM(1,ITARG)=4
      IF (INDSRC(ITARG).NE.6) THEN
        I34=EIRENE_IDEZ(INT(SORLIM(1,ITARG)),3,3)
        SORLIM(1,ITARG)=I34*100+04   !sample with step fct. 
      ELSEIF (INDSRC(ITARG).EQ.6) THEN
C  SORLIM DEFAULT WAS 0.D0
        SORLIM(1,ITARG)=0204
      ENDIF
      SORIND(1,ITARG)=ITARG

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


        NINITL(ITARG)=ITARG*1001
        NSPEZ(ITARG)=-1
        SORIFL(1,ITARG)=NIFLG(ITARG,1)
        SORWGT(1,ITARG)=1.
Cdr  USE ENERGY FLUXES SPECIFIED HERE, IE., SORENE, SORENI ARE REDUNDANT
cdr     NEMODS(ITARG)=9

        IF (NIXY(ITARG,1).EQ.1) THEN
C TARGET RECYCLING SOURCE AT POLOIDAL SURFACE NPES
          NEMODS(ITARG)=3
          NAMODS(ITARG)=1
          SORENI(ITARG)=3.
          SORENE(ITARG)=0.5
        ELSEIF (NIXY(ITARG,1).EQ.2) THEN
C WALL RECYCLING SOURCE AT RADIAL SURFACE NPES
          NEMODS(ITARG)=2
          NAMODS(ITARG)=1
          SORENI(ITARG)=2.
          SORENE(ITARG)=0.
        ENDIF
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
3999  CONTINUE
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
        IF (IGSTEP(ITARG,IG).GT.200000) THEN
          IY=IRSTEP(ITARG,IG)
          NPES=IGSTEP(ITARG,IG)-200000
          IF (NEM.EQ.3.OR.NEM.EQ.5.OR.NEM.EQ.7) THEN
            DO 6005 IPL=1,NPLSI
              IPLV=MPLSV(IPL)
              PM1=(PPLNX(IY,NPES)*VXSTEP(IPLV,ITARG,IG)+
     .             PPLNY(IY,NPES)*VYSTEP(IPLV,ITARG,IG))*OR
              VPZ=VZSTEP(IPLV,ITARG,IG)
              VP(IPL)=SQRT(PM1**2+VPZ**2)
              DI(IPL)=DISTEP(IPL,ITARG,IG)
6005        CONTINUE
            TE=TESTEP(ITARG,IG)
            CUR=0.
            GAMMA=0.
            ESHT(ITARG,IG)=EIRENE_SHEATH(TE,DI,VP,NCHRGP,GAMMA,CUR,
     .                          NPLSI,-ITARG)
          ELSE IF (NEM == 9) THEN
            TE=TESTEP(ITARG,IG)
            ESHT(ITARG,IG)=DELTA_SHEATHXB(NPBS,IY)*TE
          END IF
        ELSEIF (IGSTEP(ITARG,IG).LT.200000) THEN
          IX=IPSTEP(ITARG,IG)
          NPES=IGSTEP(ITARG,IG)-100000
          IF (NEM.EQ.3.OR.NEM.EQ.5.OR.NEM.EQ.7) THEN
C  in order to find sheath potential, we need ALL plasma particle flux components
            DO 6006 IPL=1,NPLSI
              IPLV=MPLSV(IPL)
              PM1=(PPLNX(NPES,IX)*VXSTEP(IPLV,ITARG,IG)+
     .             PPLNY(NPES,IX)*VYSTEP(IPLV,ITARG,IG))*OR
              VPZ=VZSTEP(IPLV,ITARG,IG)
              VP(IPL)=SQRT(PM1**2+VPZ**2)
              DI(IPL)=DISTEP(IPL,ITARG,IG)
6006        CONTINUE
            TE=TESTEP(ITARG,IG)
            CUR=0.
            GAMMA=0.
            ESHT(ITARG,IG)=EIRENE_SHEATH(TE,DI,VP,NCHRGP,GAMMA,CUR,
     .                           NPLSI,-ITARG)
          ELSE IF (NEM == 9) THEN
            TE=TESTEP(ITARG,IG)
            ESHT(ITARG,IG)=DELTA_SHEATHYB(IX,NPBS)*TE
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
            IY=IRSTEP(ITARG,IG)
            NPES=IGSTEP(ITARG,IG)-200000
            VT=SQRT(2.*TISTEP(IPLSTI,ITARG,IG)/BMASS(IPLS))*CVEL2A
C  VELOCITY COMPONENT NORMAL TO POLOIDAL TARGET SURFACE
C  I.E., POLOIDAL COMPONENT V-POL
C  ASSUMING ORTHOGONAL TARGET
            PM1=(PPLNX(IY,NPES)*VXSTEP(IPLSV,ITARG,IG)+
     .           PPLNY(IY,NPES)*VYSTEP(IPLSV,ITARG,IG))*OR
C  VELOCITY COMPONENT PARALLEL TO POLOIDAL TARGET SURFACE
C  I.E., RADIAL PLUS TOROIDAL COMPONENT, V-RAD + V-TOR
C  AGAIN: ASSUMING ORTHOGONAL TARGET
            VPX=VXSTEP(IPLSV,ITARG,IG)-PM1*PPLNX(IY,NPES)*OR
            VPY=VYSTEP(IPLSV,ITARG,IG)-PM1*PPLNY(IY,NPES)*OR
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
     .                  TESTEP(ITARG,IG))/BMASS(IPLS))*CVEL2A
C THE MACH NUMBER BOUNDARY CONDITION ONLY AFFECTS THE PARALLEL TO B
C MOMENTUM, I.E., NOT THE RADIAL VELOCITY
            VTEST=SQRT(PM1**2+VPZ**2)
            VTEST=VTEST/(CS+EPS60)
            VR=SQRT(VPX**2+VPY**2)
            VTEST2=VPZ/(CS+EPS60)
            IF (TRCINT) THEN
              WRITE (iunout,*) 'IPL,ITG,IG,MACH_PAR',
     .                          IPLS,ITARG,IG,VTEST
C             WRITE (iunout,*) 'POL., TOR., RAD. (CM/S)',PM1,VPZ,VR
C             CALL EIRENE_LEER(1)
            ENDIF
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
            IF (IGSTEP(ITARG,IG).LT.200000) THEN
              EADD=0.
              WRITE (iunout,*) 'INVALID OPTION FOUND IN IF2COP '
              WRITE (iunout,*) 'POSSIBLE ERROR IN TARGET ENERGY FLUX '
              WRITE (iunout,*) 'ITARG,IG ',ITARG,IG
            ELSE
              PERWI=PERW/SQRT(BMASS(IPLS)/RMASSP(IPLS))
              PARWI=PARW/SQRT(BMASS(IPLS)/RMASSP(IPLS))
              EADD=EIRENE_EMAXW(TISTEP(IPLSTI,ITARG,IG),PERWI,PARWI)
            ENDIF
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

6009    CONTINUE  ! IPLS loop
        GOTO 6011
6010    CONTINUE
C  TO BE WRITTEN
6011  CONTINUE    ! IG,  CELL ALONG TARGET
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
6100  CONTINUE
      WRITE (iunout,'(1X,I3,1P,1E11.3)') NRWL(ITARG),
     .                                 RRSTEP(ITARG,NRWL(ITARG))
      CALL EIRENE_MASR1 ('EEMAX    ',EEMAX)
      CALL EIRENE_MASR1 ('EESHT    ',EESHT)
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
6300  CONTINUE
C
C  SET SOME OTHER DATA SPECIFIC FOR EIRENE CODE REQUIREMENTS
C  STATEMENT NO. 6500 ---> 6999
C
6500  CONTINUE

      DEALLOCATE (TORL)
      DEALLOCATE (ESHT)
      DEALLOCATE (ELTEST)
      DEALLOCATE (ORI)
C
C
      RETURN
999   CONTINUE
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
      NDXY=(NDXA-1)*NR1ST+NDYA
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
      NDXY=(NDXA-1)*NR1ST+NDYA
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
        CALL EIRENE_ALLOC_EIRBRA(NDXP, NDYP, NFL, NSTRA,IFOFF)
      END IF
C
      IF (NEW_ITER == 0) THEN
        RESSNI(ISTRAA:ISTRAE,:) = 0._DP
        RESSMO(ISTRAA:ISTRAE,:) = 0._DP
        RESSEE(ISTRAA:ISTRAE) = 0._DP
        RESSEI(ISTRAA:ISTRAE) = 0._DP
      ENDIF

      volSUMN(ISTRAA:ISTRAE)=0.0           !dpc
      volSUMM(ISTRAA:ISTRAE)=0.0           !dpc
      volSUMEI(ISTRAA:ISTRAE)=0.0          !dpc
      volSUMEE(ISTRAA:ISTRAE)=0.0          !dpc

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
7000    CONTINUE
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

7310    CONTINUE

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

7330    CONTINUE
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

7350    CONTINUE
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
7400    CONTINUE
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
              DO 7471 IR=1,NR1ST-1
              DO 7471 K=1,NPPLG
              DO 7471 IP=NPOINT(1,K),NPOINT(2,K)-1
c DPC 1997.05.20 added check for valid point in specified vol. rec. domain
                if(ir.ge.INGRDA(1,ISTRAI,1).and.
     1             ir.lt.INGRDE(1,ISTRAI,1).and.
     2             ip.ge.INGRDA(1,ISTRAI,2).and.
     3             ip.lt.INGRDE(1,ISTRAI,2)) then
c dpc
                IN=(IP-1)*NR1ST+IR
                INC=NCLTAL(IN)
                IF (NSTORDR >= NRAD) THEN
                  RECADD=-TABRC1(IRRC,IN)*DIIN(IPLS,IN)*ELCHA
                  EEADD=  EELRC1(IRRC,IN)*DIIN(IPLS,IN)*ELCHA
                ELSE
                  RECADD=-EIRENE_FTABRC1(IRRC,IN)*DIIN(IPLS,IN)*ELCHA
                  EEADD=  EIRENE_FEELRC1(IRRC,IN)*DIIN(IPLS,IN)*ELCHA
                END IF
                PPPL_COP(IPLS,INC)=PPPL_COP(IPLS,INC)+RECADD
                SUMN=SUMN+RECADD*VOL(IN)
                PIADD=PARMOM(IPLS,IN)*RECADD
                MPPL_COP(IPLS,INC)=MPPL_COP(IPLS,INC)+PIADD
                SUMM=SUMM+PIADD*VOL(IN)
                EIADD=(1.5*TIIN(IPLSTI,IN)+EDRIFT(IPLS,IN))*RECADD
                EPPL_COP(IPLS,INC)=EPPL_COP(IPLS,INC)+EIADD
                SUMEI=SUMEI+EIADD*VOL(IN)
                EPEL_COP(INC)=EPEL_COP(INC)+EEADD
                SUMEE=SUMEE+EEADD*VOL(IN)
                END IF

7471          CONTINUE  ! loop over grid


              WRITE (iunout,*) 'IPLS,IRRC ',IPLS,IRRC
              CALL EIRENE_MASR4('SUMN, SUMM, SUMEI, SUMEE        ',
     .                     SUMN,SUMM,SUMEI,SUMEE)
c  now sum over IRRC rec processes for bulk ion IPLS
              volSUMN(ISTRAI)=volSUMN(ISTRAI)+SUMN             ! dpc
              volSUMM(ISTRAI)=volSUMM(ISTRAI)+SUMM             ! dpc
              volSUMEI(ISTRAI)=volSUMEI(ISTRAI)+SUMEI          ! dpc
              volSUMEE(ISTRAI)=volSUMEE(ISTRAI)+SUMEE          ! dpc
7472        CONTINUE
7473      CONTINUE
        ENDIF
C
        IF (.NOT.LSYMET) GOTO 7500
C
C  SECONDLY SYMMETRISE EIRENE ARRAYS ACCORDING TO SYMMETRY IN MODEL
C
C
C   THIRDLY WRITE EIRENE ARRAYS (1D) ONTO BRAAMS ARRAYS (2D)
C   AND RESCALE TO PROPER UNITS: #/CELL/STRATUM FLUX
C   # STANDS FOR PARTICLES (SNI), MOMENTUM (SMO)
C   AND ENERGY (SEE,SEI)
C
7500    CONTINUE
        DO 7510 IFL=1,NFLA
          CHPS(IFL)=0.
          SNIS(IFL)=0.
          CHMOS(IFL)=0.
          SMOS(IFL)=0.

          DO 7510 IPLS=1,NPLSI
            IF (IFLB(IPLS).NE.IFL) GOTO 7510
            IPLSV=MPLSV(IPLS)
c  ipls contributes to plasma code species ifl
            DO 7520 IX=1,NDXA
              DO 7530 IY=1,NDYA
                IN=IY+(IX-1)*NR1ST
                INC=NCLTAL(IN)
                SNICL=(PAPL(IPLS,INC)+PMPL(IPLS,INC)+PIPL(IPLS,INC)+
     .                 PPPL_COP(IPLS,INC))*VOLTAL(INC)
                SNI(IX,IY,IFL,ISTRAI)=SNI(IX,IY,IFL,ISTRAI)+SNICL
                SNIS(IFL)=SNIS(IFL)+SNICL
                CHPS(IFL)=CHPS(IFL)+CHPM(IPLS,INC)*VOLTAL(INC)
7530          CONTINUE
7520        CONTINUE
            DO 7536 IX=1,NDXA
              IF (LLCUT(IX)) CYCLE
              DO 7533 IY=1,NDYA
                IN=IY+(IX-1)*NR1ST
                INC=NCLTAL(IN)
                SIGNUM=SIGN(1._DP,BVIN(IPLSV,IN))
                SMOCL=(MAPL(IPLS,INC)+MMPL(IPLS,INC)+MIPL(IPLS,INC)+
     .                 MPPL_COP(IPLS,INC))*
     .                 VOLTAL(INC)*1.D-5*SIGNUM
                SMO(IX,IY,IFL,ISTRAI)=SMO(IX,IY,IFL,ISTRAI)+SMOCL
                SMOS(IFL)=SMOS(IFL)+SMOCL
                CHMOS(IFL)=CHMOS(IFL)+CHMOM(IPLS,INC)*VOLTAL(INC)
7533          CONTINUE
7536        CONTINUE
!pb 7539        CONTINUE
7510    CONTINUE

C
        CHEES=0.
        SEES=0.
        DO 7540 IX=1,NDXA
          DO 7545 IY=1,NDYA
            INN=IY+(IX-1)*NR1ST
            IN=NCLTAL(INN)
            SEE(IX,IY,ISTRAI)=(EAEL(IN)+EMEL(IN)+
     .                         EIEL(IN)+EPEL_COP(IN))*
     .                         VOLTAL(IN)
            CHEES=CHEES+CHEEM(IN)*VOLTAL(IN)
            SEES=SEES+SEE(IX,IY,ISTRAI)
            SEE(IX,IY,ISTRAI)=SEE(IX,IY,ISTRAI)*ELCHA
7545      CONTINUE
7540    CONTINUE
C
        CHEIS=0.
        SEIS=0.
        DO 7544 IFL=1,NFLA
          DO  7543 IPLS=1,NPLSI
            IF (IFLB(IPLS).NE.IFL) GOTO 7543
            DO 7542 IX=1,NDXA
              IF (LLCUT(IX)) CYCLE
              DO 7541 IY=1,NDYA
                INN=IY+(IX-1)*NR1ST
                IN=NCLTAL(INN)
                SEI(IX,IY,ISTRAI)=SEI(IX,IY,ISTRAI) + 
     .                           (EAPL(IPLS,IN)+EMPL(IPLS,IN)+
     .                            EIPL(IPLS,IN)+EPPL_COP(IPLS,IN))*
     .                            VOLTAL(IN)
 7541         CONTINUE
 7542       CONTINUE

 7543     CONTINUE
 7544   CONTINUE

        DO IX = 1, NDXA
          IF (LLCUT(IX)) CYCLE
          DO IY=1,NDYA
            INN=IY+(IX-1)*NR1ST
            IN=NCLTAL(INN)
            CHEIS=CHEIS+CHEIM(IN)*VOLTAL(IN)
            SEIS=SEIS+SEI(IX,IY,ISTRAI)
            SEI(IX,IY,ISTRAI)=SEI(IX,IY,ISTRAI)*ELCHA
          END DO
        END DO
C
C   NEXT:
C   IF LSHORT: CRITERION TO STOP SHORT CYCLE,
C   IF NOT LSHORT: RESCALE SURFACE SOURCE STRATA
C                  UNITS: # PER UNIT TARGET PLATE FLUX
C
        IF (LSHORT.AND.IFIRST.EQ.0) THEN
C
          SNIS0(ISTRAI,0)=0.
          SMOS0(ISTRAI,0)=0.
          DO 7550 IFL=1,NFLA
            SNIS0(ISTRAI,0)=SNIS0(ISTRAI,0)+SNIS(IFL)
            SNIS0(ISTRAI,IFL)=SNIS(IFL)
            SMOS0(ISTRAI,0)=SMOS0(ISTRAI,0)+SMOS(IFL)
            SMOS0(ISTRAI,IFL)=SMOS(IFL)
7550      CONTINUE
          SEES0(ISTRAI)=SEES
          SEIS0(ISTRAI)=SEIS
C
        ELSEIF (LSHORT.AND.IFIRST.GT.0) THEN
C
          SNIS(0)=0.
          DO 7551 IFL=1,NFLA
            SNIS(0)=SNIS(0)+SNIS(IFL)
            SCALN(IFL)=SNIS0(ISTRAI,IFL)/(SNIS(IFL)+EPS60)
7551      CONTINUE

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
7552        CONTINUE
7555      CONTINUE
          DO 7556 IFL=1,NFLA
            DO 7553 IX=0,NDXA+1
              DO 7554 IY=0,NDYA+1
                SNI(IX,IY,IFL,ISTRAI)=SNI(IX,IY,IFL,ISTRAI)*SCALN(IFL)
                SMO(IX,IY,IFL,ISTRAI)=SMO(IX,IY,IFL,ISTRAI)*SCALN(IFL)
7554          CONTINUE
7553        CONTINUE
7556      CONTINUE
C
          LTEST=.TRUE.  ! tactically assume: this stratum will continue in short cycle mode

          IF (LSTOP) THEN
            WRITE (iunout,*) 'STOP SHORT CYCLE: ALL B2 TIMESTEPS DONE '

          ELSE   ! DO AT LEAST ONE MORE TIME STEP
CDR  DECIDE FOR THIS CURRENT STRATUM ISTRAI: 
CDR     SHORT CYCLE (IMPLICIT CORRECTION) ONLY, OR FULL MONTE CARLO


            DO 7558 IFL=1,NFLA
              TEST=CHPS(IFL)/(SNIS(IFL)+1.D-60)*100.
              write (iunout,*) ' global change in sni,ifl ',test,ifl
              IF (ABS(TEST).GT.CHGP) THEN
                LSTP3=.TRUE.
                LTEST=.FALSE.  !  stop short cycle mode. Full new set of  trajectories.
                WRITE (iunout,*) 'STOP SHORT CYCLE: PART. SOURCES: ',
     .                            SNIS(IFL),CHPS(IFL),TEST
                WRITE (iunout,*) 'STRATUM ISTRAI, SPECIES IFL ',
     .                            ISTRAI,IFL
              ENDIF
              TEST=CHMOS(IFL)/(SMOS(IFL)+1.D-60)*100.
              write (iunout,*) ' global change in smo,ifl ',test,ifl
              IF (ABS(TEST).GT.CHGMOM) THEN
                LSTP3=.TRUE.
                LTEST=.FALSE. !  stop short cycle mode. Full new set of  trajectories.
                WRITE (iunout,*) 'STOP SHORT CYCLE: MOMENTUM SOURCE: ',
     .                            SMOS(IFL),CHMOS(IFL),TEST
                WRITE (iunout,*) 'STRATUM ISTRAI, SPECIES IFL ',
     .                            ISTRAI,IFL
              ENDIF
7558        CONTINUE

            TEST=CHEES/(SEES+1.D-60)*100.
            write (iunout,*) ' global change in see,ifl ',test,ifl
            IF (ABS(TEST).GT.CHGEE) THEN
              LSTP3=.TRUE.
              LTEST=.FALSE. !  stop short cycle mode. Full new set of  trajectories.
              WRITE (iunout,*) 'STOP SHORT CYCLE: EL EN. SOURCE: ',SEES,
     .                          CHEES,TEST
              WRITE (iunout,*) 'STRATUM ISTRAI ',ISTRAI
            ENDIF
            TEST=CHEIS/(SEIS+1.D-60)*100.
            write (iunout,*) ' global change in sei,ifl ',test,ifl
            IF (ABS(TEST).GT.CHGEI) THEN
              LSTP3=.TRUE.
              LTEST=.FALSE. !  stop short cycle mode. Full new set of  trajectories.
              WRITE (iunout,*) 'STOP SHORT CYCLE: ION EN. SOURCE: ',
     .                          SEIS,CHEIS,TEST
              WRITE (iunout,*) 'STRATUM ISTRAI ',ISTRAI
            ENDIF
            IF (LSHORT) LSTP=LSTP3
          ENDIF
C
        ELSEIF (.NOT.LSHORT) THEN
C
          DO 7560 IX=0,NDXA+1
            DO 7565 IY=0,NDYA+1
              SEE(IX,IY,ISTRAI)=SEE(IX,IY,ISTRAI)*FLXI
              SEI(IX,IY,ISTRAI)=SEI(IX,IY,ISTRAI)*FLXI
7565        CONTINUE
7560      CONTINUE
          DO 7570 IFL=1,NFLA
            DO 7580 IX=0,NDXA+1
              DO 7590 IY=0,NDYA+1
                SNI(IX,IY,IFL,ISTRAI)=SNI(IX,IY,IFL,ISTRAI)*FLXI
                SMO(IX,IY,IFL,ISTRAI)=SMO(IX,IY,IFL,ISTRAI)*FLXI
7590          CONTINUE
7580        CONTINUE
7570      CONTINUE
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
7600    CONTINUE
C
7700    CONTINUE
C
7999  CONTINUE
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
C  BALANCES, SHOULD BE DONE ONLY AT THE END OF B2.5 RUN
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
C
      DO 10113 IX=1,NDXA
C
C IS (IX,0) A RECYCLING SOURCE? IF YES, DO NOT COUNT HERE
C
        DO ITARG=1,NTARGI
          DO IPRT=1,NTGPRT(ITARG)
            IF (NIXY(ITARG,IPRT).EQ.2.AND.NDT(ITARG,IPRT).EQ.0) THEN
              IF (IX.GE.NTIN(ITARG,IPRT).AND.
     .            IX.LT.NTEN(ITARG,IPRT)) GOTO 10113
            ENDIF
          ENDDO
        ENDDO

        ITARG=0
C
C SURFACE NORMAL IS INWARD. HENCE: TAKE ALL FLUXES F...YB POSITIVE
C SIGN OF ADDITIONAL COMPONENT DUE TO INCLINED GRID AS SIGN OF F...YB
C
        SFEISY=SFEISY+FEIYB(IX,0)
        SFEESY=SFEESY+FEEYB(IX,0)
        DO 10111 IFL=1,NFLA
          SFNISY(IFL)=SFNISY(IFL)+FNIYB(IX,0,IFL)
10111   CONTINUE
10113 CONTINUE
C
      SFNISY=SFNISY*ELCHA
C
      LNONREC_SY=ANY(SFNISY(1:nfla).NE.0.0).OR.SFEISY.NE.0.0.OR.
     .                                         SFEESY.NE.0.0
      WRITE (37,*) 'NON-RECYCLING FLUXES FROM SOUTH EDGE '
      WRITE (37,8888) SFNISY,SFEISY,SFEESY
C
C
C  SECOND: NORTH EDGE: IY=NDYA
C
C NON RECYCLING FLUXES AT NORTH EDGE: SFEINY,SFEENY,SFNINY
      SFEINY=0.
      SFEENY=0.
      SFNINY=0.
      DO 10118 IX=1,NDXA
C
C IS (IX,NDYA) A RECYCLING SOURCE? IF YES, DO NOT COUNT HERE
C
        DO ITARG=1,NTARGI
          DO IPRT=1,NTGPRT(ITARG)
            IF (NIXY(ITARG,IPRT).EQ.2.AND.NDT(ITARG,IPRT).EQ.NDYA) THEN
              IF (IX.GE.NTIN(ITARG,IPRT).AND.
     .            IX.LT.NTEN(ITARG,IPRT)) GOTO 10118
            ENDIF
          ENDDO
        ENDDO
        ITARG=0
C
C SURFACE NORMAL IS OUTWARD. HENCE: TAKE ALL FLUXES F...YB NEGATIVE
C SIGN OF ADDITIONAL COMPONENT DUE TO INCLINED GRID AS SIGN OF F...YB
C
        SFEINY=SFEINY-FEIYB(IX,NDYA)
        SFEENY=SFEENY-FEEYB(IX,NDYA)
        DO 10116 IFL=1,NFLA
          SFNINY(IFL)=SFNINY(IFL)-FNIYB(IX,NDYA,IFL)
10116   CONTINUE
10118 CONTINUE
C
      SFNINY=SFNINY*ELCHA
C
      LNONREC_NY=ANY(SFNINY(1:nfla).NE.0.0).OR.SFEINY.NE.0.0.OR.
     .                                         SFEENY.NE.0.0
      WRITE (37,*) 'NON-RECYCLING FLUXES TO NORTH EDGE '
      WRITE (37,8888) SFNINY,SFEINY,SFEENY
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
     .             IY.LT.NTEN(ITARG,IPRT)) GOTO 10123
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
10123 CONTINUE
C
      SFNIWX=SFNIWX*ELCHA
C
      LNONREC_WX=ANY(SFNIWX(1:nfla).NE.0.0).OR.SFEIWX.NE.0.0.OR.
     .                                         SFEEWX.NE.0.0
      WRITE (37,*) 'NON-RECYCLING FLUXES FROM WEST EDGE '
      WRITE (37,8888) SFNIWX,SFEIWX,SFEEWX
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
     .            IY.LT.NTEN(ITARG,IPRT)) GOTO 10128
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
10128 CONTINUE
C
      SFNIEX=SFNIEX*ELCHA
C
      LNONREC_EX=ANY(SFNIEX(1:nfla).NE.0.0).OR.SFEIEX.NE.0.0.OR.
     .                                         SFEEEX.NE.0.0
      WRITE (37,*) 'NON-RECYCLING FLUXES TO EAST EDGE '
      WRITE (37,8888) SFNIEX,SFEIEX,SFEEEX
C
C  NEXT: FLUXES TO THOSE SURFACES, AT WHICH RECYCLING BOUNDARY
C        CONDITIONS ARE SPECIFIED
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
          IF (NIXY(I,IPRT).EQ.1) THEN
C  BALANCE CONTRIB. X-GRID REC. SOURCE
            DO 10132 IY=NTIN(I,IPRT),NTEN(I,IPRT)-1
              SFEIT(I)=SFEIT(I)-NINCT(I,IPRT)*FEIXB(NDT(I,IPRT),IY)
              SFEET(I)=SFEET(I)-NINCT(I,IPRT)*FEEXB(NDT(I,IPRT),IY)
              DO 10131 IFL=NSPZI(I,IPRT),NSPZE(I,IPRT)
                IF (NINCT(I,IPRT)*FNIXB(NDT(I,IPRT),IY,IFL).GT.0) THEN
                SFNIT(I,IFL)=SFNIT(I,IFL)-
     .                   NINCT(I,IPRT)*FNIXB(NDT(I,IPRT),IY,IFL)

cdr sheath contributions: count negative for electrons, positive for ions 
cdr unfinished:  need to account for charge state of ion species IFL
                SHEAE(I)=SHEAE(I)+TEB(NDT(I,IPRT),IY)*
     .           NINCT(I,IPRT)*FNIXB(NDT(I,IPRT),IY,IFL)*
     .          (-2.8)
                SHEAI(I)=SHEAI(I)+TEB(NDT(I,IPRT),IY)*
     .           NINCT(I,IPRT)*FNIXB(NDT(I,IPRT),IY,IFL)*
     .           2.8
cdr  sheath done
                ELSE
                  WRITE (iunout,*)
     .              'WRONG ORIENTATION OF W/E-TARGET RECYCLING FLUX '
                  WRITE (iunout,*) 'ITARG, IPRT, IPLS, NDT, IY ',
     .                         I    , IPRT, IFL,   NDT(I,IPRT), IY
                  WRITE (iunout,*) 'FNIX(NDT,IY) ',
     .                              FNIXB(NDT(I,IPRT),IY,IFL)
                ENDIF
10131         CONTINUE
10132       CONTINUE

C  BALANCE CONTRIB. FROM Y-GRID RECYCLING SOURCE
          ELSEIF (NIXY(I,IPRT).EQ.2) THEN
            DO 10135 IX=NTIN(I,IPRT),NTEN(I,IPRT)-1
              SFEIT(I)=SFEIT(I)-NINCT(I,IPRT)*FEIYB(IX,NDT(I,IPRT))
              SFEET(I)=SFEET(I)-NINCT(I,IPRT)*FEEYB(IX,NDT(I,IPRT))
              DO 10136 IFL=NSPZI(I,IPRT),NSPZE(I,IPRT)
                IF (NINCT(I,IPRT)*FNIYB(IX,NDT(I,IPRT),IFL).GT.0.) THEN
                SFNIT(I,IFL)=SFNIT(I,IFL)-
     .                   NINCT(I,IPRT)*FNIYB(IX,NDT(I,IPRT),IFL)

cdr sheath contributions: count negative for electrons, positive for ions 
cdr unfinished:  need to account for charge state of ion species IFL
                SHEAE(I)=SHEAE(I)+TEB(IX,NDT(I,IPRT))*
     .           NINCT(I,IPRT)*FNIYB(IX,NDT(I,IPRT),IFL)*
     .           (-2.8)
                SHEAI(I)=SHEAI(I)+TEB(IX,NDT(I,IPRT))*
     .           NINCT(I,IPRT)*FNIYB(IX,NDT(I,IPRT),IFL)*
     .           (2.8)
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
      WRITE (37,*) 'CHARGED IMPURITY RAD.,IONIS. AND RECOMB. '
      WRITE (37,8888) 0.,0.,B2RAD
C
      WRITE (37,*) 'ELECTRIC FIELD TERMS (PRESSURE GRADIENTS)'
      WRITE (37,8888) 0.,B2VDP,-B2VDP
C
      BALANI=SFEISY+SFEINY+SFEIT(0)+SHEAI(0)+SSEI+B2QIE+B2VDP+
     .       SFEIWX+SFEIEX
      BALANE=SFEESY+SFEENY+SFEET(0)+SHEAE(0)+SSEE+B2BREM+B2RAD-B2QIE+
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
8888  FORMAT (3E14.6)
      END
