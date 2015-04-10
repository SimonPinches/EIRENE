cdr  150407:  orientation of B field made optional, additional input flags ibrad,ibpol,ibtor in block 14.
cdr           magnitude of bfield (T) tranfered.
cdr           do be checked: orientation of uudiag for reconstruction of carthesian flow velocity components.
cdr           apparently not used: vvdiag.



C   EIRENE CODE SEGMENT COUPLE_$, $ MAY CURRENTLY STAND FOR B2,
C                                                           B2.5,
C                                                           DIVIMP,
C                                                           TRIA,
C                                                           TETRA,
C                                                           TRANSP,
C                                                           DUMMY
C
C   THIS VERSION: $B2.5/$TRIA combined by s.wiesen@fz-juelich.de, 2011
C
c  geometry data not any longer via work array into eirene
c                due to module structure
c  eliminate cut cells from balances (lcut(..))
c  new input: ncopib, ncopeb
c
c  additionally modified for MPI use, s.wiesen@fz-juelich.de, apr2011
c
C
C   UPDATES:
C   OPTION TO EVALUATE B-FIELD VECTORS FROM GRIDADAP FILE FT29
C   FOR NON-ORTHOGONAL GRIDS
C
C   THIS CODE SEGMENT CONTAINES VARIOUS SUBROUTINES NEEDED FOR
C   INTERFACING THE EIRENE CODE TO PLASMA FLUID CODES.
C   IT READS GEOMETRICAL DATA (MESHES) FROM FILE FT30
C   AND PRODUCES THE EIRENE INPUT DATA (BLOCK 2).
C   IT READS PLASMA BACKGROUND DATA FROM FILE FT31 OR COMMON BLOCKS,
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
C   IT WAS WRITTEN BY D.REITER AND P.BOERNER, KFA-JUELICH
C   E-MAIL: D.REITER @ EIRENE.DE
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
      USE EIRMOD_BRAEIR
      USE EIRMOD_EIRBRA
      USE EIRMOD_BRASCL
      USE EIRMOD_CSPEZ
csw mpi
      use eirmod_cpes
csw
csw 26jan2011 extra B25
      use eirmod_extraB25
csw
      IMPLICIT NONE
csw mpi
      include 'mpif.h'
      integer :: ier,istrx,irank,istrr, irnk
      real*8, allocatable :: dumvec(:)
      real*8, allocatable :: save_estimv(:,:),save_estims(:,:),
     .                       save_sigma_cop(:,:),wtotp_dum(:,:)
csw
C
C  GEOMETRICAL DATA FROM GRIDADAP
      REAL(DP), ALLOCATABLE ::
     R  ALPHXB(:,:), ALPHYB(:,:), XAISO(:,:)

      REAL(DP), ALLOCATABLE, SAVE ::
     R  PUX(:),      PUY(:),      PVX(:),      PVY(:),
     R  PUXE(:), PUYE(:), PUXN(:), PUYN(:),
     R  PVXE(:), PVYE(:), PVXN(:), PVYN(:)

      INTEGER, ALLOCATABLE ::
     I  IAISO(:,:)
C
      TYPE(CELLSIM), POINTER :: CPSIM
      TYPE(CELLMUL), POINTER :: CPMUL
C
      REAL(DP) :: SEES0(NSTRA), SEIS0(NSTRA)
      REAL(DP) :: CHPM(NPLS,NRAD), CHEEM(NRAD), CHEIM(NRAD),
     .            CHMOM(NPLS,NRAD)
      REAL(DP) :: DI(NPLS), VP(NPLS)
      REAL(DP) :: SFNISY(NFL),SFNINY(NFL),SFNIWX(NFL),SFNIEX(NFL)
      REAL(DP) :: SSN(NFL),SSNI(NFL),BALANN(NFL),TOTN(NFL),RN(NFL)
      REAL(DP) :: PPPL_COP(NPLS,NRAD), CPPV(NCPV,NRAD),
     .            EPPL_COP(NRAD), EPEL(NRAD)
      REAL(DP) :: PPLODA(NPLS,NRAD), CPVODA(NCPV,NRAD), 
     .            EPLODA(NRAD), EPEODA(NRAD)
C
      REAL(DP) :: EFLX(NSTRA),
     R          DUMMY(0:NDXP,0:NDYP),
     R          SFNIT(0:NSTEP,NFL), SFEIT(0:NSTEP),
     R          SFEET(0:NSTEP), SHEAE(0:NSTEP), SHEAI(0:NSTEP)

      INTEGER :: NRWL(NSTRA)

      REAL(DP), SAVE :: SCALM, SCALE, SCALI, CHEIS, SEES, SEIS, TEST,
     .          SFEISY, SFEESY, RECADD, RECTOT,
     .          EEADD, PIADD, SIGNUM, SMOCL, CHEES, EIADD, SNICL,
     .          SSE, BALANI, BALANE, SSEE, SSI, RE, RI, RNT, TOT,
     .          TOTI, TOTE, SFEENY, SFEIWX, BALAN, RRBC,
     .          SSEI, SFEIEX, SFEEEX, VVBC,
     .          UUBC, UPBC, RBC, UDBC, VL, V, T, BX, BY, BZ, BN,
     .          DELTE_PARA, DELTI_PARA, DELY, DELTI_PERP, TES, TIS,
     .          DELTE_PERP, ALX, ALE, ALW, ALS, ALN, AL, ETOT,
     .          FLX, ESUM, DR, VR, VTEST, EADD, SI,VTEST2,
     .          PARWI, PERWI, SUMM, SUMN, SUMEI, SUMEE, FLXI, CHP,
     .          CNDYNP, CHI, CHE, CS, THMAX, EESHT, EEMAX,
     .          RP1, DELX, PNORM, PVYS, PVXS, PUPV, RRBS, PUYS, PUXS,
     .          VPX, VPY, VT, PARW, PERW, PN1, OR, VPZ, GAMMA, CUR, TE,
     .          SFEEWX, SFEINY, PM1, DRR, VDBC,
     .          TIFLX,XCOOR,YCOOR,ZCOOR,VSX,VSY,VS,VTX,
     .          BVAC,TX,TY,VPRO,VTY,XMUE,PX,PY,
     .          XANF,YANF,PIPV,FLX_EIR,
     .          SUMN_OLD,SNIRES,SMORES,SEERES,SEIRES,UU,PITB,
     .          DXPOL,DYPOL,PAR,dx,dy,brad,bpol,btor

       INTEGER, SAVE :: J, IRC, JC, INC, IADD, IP, ITARG, IO, IFL, NPES,
     .           IIPLS, IG, IGITT, IEPLS, NPEC, NPBC, NPBS, NTGPRI,
     .           IT, I, IPRT, IAOT, IAIN, IREAD, IPL, INN,
     .           IMODE, IERROR, LTARG, IN, IX, IY,
     .           NCOPI, NPLP, NDX2, NRED, IO29, NDXY, IFIRST,
     .           ISTRAI, IRRC, K, IR, IIRC, ICPV, IF, I34,
     .           NREC11, NEM, MINSPEZ, MAXSPEZ, ISP, IPLSTI, IPLSV,
     .           IPLV,l,
     .           NAS,IPUNKT,NSSIR,NUMSI,NBAR,ISNR,ISC,IS,NASMOD,
     .           NRS,NADMOD,NBARSI,IP1,NP2NDQ,IS1,IR1,
     .           NEND,NINI,NSSIP,MTRI,
     .           IDUMMY,NR1STQ,ISTS,ITRI,IACT,IANF,ICOG,
     .           ISC1,ISC2,ISCS,ICOU,IXI,IXE,NCOPIB,NCOPEB,
     .           IST_RATE, MSHFRM, IMF, istat_cop,ibrad,ibpol,ibtor,
     .           ntrfrm
      INTEGER, INTENT(IN) :: ISTRAA, ISTRAE, NEW_ITER, IFRST, ITRG
      REAL(DP) :: EIRENE_STEP, EIRENE_FTABRC1, EIRENE_FEELRC1, 
     .            EIRENE_SHEATH, EIRENE_EMAXW  
      INTEGER, EXTERNAL :: EIRENE_IDEZ
C
      LOGICAL, INTENT(INOUT) :: LSTP
      LOGICAL, SAVE :: LSHORT, LSTOP, LTEST, LSTP3, IFBOUND, lchkqud
csw 14apr2011, LCUT now in EIRMOD_CPOLYG (broadcasted)
csw      LOGICAL, ALLOCATABLE, SAVE :: LCUT(:)
      logical :: l1, l2, lxsrf
!pb qq not needed/used
!pb      real(DP) :: DUMVAL,ud,vv,up,qq
      real(DP) :: DUMVAL,ud,vv,up
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
C  READ PLASMA PARAMETERS, RESCALE THEM IF NECESSARY
C  AND TRANSFER THEM TO EIRENE VIA
C  EIRENE FUNCTION "PROFR"  THROUGH EQUIVALENCE ON ARRAY SMESTV
C
      REAL(DP), ALLOCATABLE, SAVE ::
     . CHPS(:),    SNIS(:),    CHMOS(:),  SMOS(:),  SCALN(:),
     . SNIS0(:,:), SMOS0(:,:),
     . RESSNI(:,:),  RESSMO(:,:), RESSEE(:), RESSEI(:), FLXEIR(:)
      REAL(DP) :: SPAT(0:NATM,0:NSTRA), SPML(0:NMOL,0:NSTRA),
     .            SPIO(0:NION,0:NSTRA), SPPL(0:NPLS,0:NSTRA) 

      REAL(DP), ALLOCATABLE ::
     . TORL(:,:), ESHT(:,:), ORI(:,:)

      REAL(DP) :: OUTHELP(NFL)


      INTEGER, ALLOCATABLE :: IHELP(:)
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
!pb      SAVE
C
csw 11apr2011
      real(dp) :: dd
csw 21feb2012
      integer :: iat
csw
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
C
      IERROR=0
C
      IMODE=IABS(NMODE)
!pb 
      lchkqud = .false.
      mshfrm = 0
      NLSHRT13 = .TRUE.
      ntrfrm = 0
C
      IF (.NOT.LSHORT.AND.ITIMV.LE.1) THEN
        WRITE (iunout,*) '        SUBROUTINE INFCOP IS CALLED  '
C  READ INPUT DATA OF BLOCK 14
C  SAVE INPUT DATA OF BLOCK 14 FOR SHORT CYCLE ON COMMON CCOUPL
        CALL EIRENE_LEER(1)
        CALL EIRENE_ALLOC_CCOUPL(1)
        READ (IUNIN,'(5L1)') LSYMET,LBALAN,LCHKQUD
        IF (TRCINT)
     .  WRITE (iunout,*) ' LSYMET,LBALAN = ',LSYMET,LBALAN
        READ (IUNIN,'(9I6)') NFLA,NCUTB,NCUTL,IMF,
     .                       ntrfrm, nfull,ibrad,ibpol,ibtor
        IF (IMF /= 0) MSHFRM=IMF

cdr added in april 2015:
c  flags for orientation of radial (not in use), poloidal and toroidal magnetic field components
        brad = 1._dp
        if (ibrad < 0) brad = -brad
        bpol = 1._dp
        if (ibpol < 0) bpol = -bpol
        btor = 1._dp
        if (ibtor < 0) btor = -btor


        NCUTB_SAVE=NCUTB
        IF (TRCINT) THEN
          WRITE (iunout,*) ' NFLA,NCUTB,NCUTL = ',NFLA,NCUTB,NCUTL
          WRITE (iunout,*) ' IPLS,IFLB(IPLS),FCTE(IPLS),BMASS(IPLS)'
        ENDIF
        DO 20 IPL=1,NPLSI
          READ (IUNIN,'(2I6,2E12.4)') I,IFLB(IPL),FCTE(IPL),BMASS(IPL)
          IF (TRCINT)
     .    WRITE (iunout,*) IPL,IFLB(IPL),FCTE(IPL),BMASS(IPL)
20      CONTINUE
        READ (IUNIN,'(2I6)') NDXA,NDYA
        IF (TRCINT) WRITE (iunout,*) 'NDXA,NDYA ',NDXA,NDYA
C  NUMBER OF TARGET SOURCES ON B2 SURFACES: NTARGI
        READ (IUNIN,'(I6)') NTARGI
        WRITE (iunout,*) '        NTARGI= ',NTARGI
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
331         READ (IUNIN,'(A72)') ZEILE
            IREAD=1
            IF (ZEILE(1:1).EQ.'*') THEN
C             WRITE (iunout,......)
              GOTO 331
            ENDIF
            READ (ZEILE,'(12I6)') I,NDT(IT,IPRT),NINCT(IT,IPRT),
     .                              NIXY(IT,IPRT),NTIN(IT,IPRT),
     .                              NTEN(IT,IPRT),NIFLG(IT,IPRT),
     .                              NPTC(IT,IPRT),NPTCM(IT,IPRT),
     .                              NSPZI(IT,IPRT),NSPZE(IT,IPRT),
     .                              NEMOD(IT,IPRT) !VK NPTCM
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
     .      WRITE (iunout,'(1X,12I6)') IT,NDT(IT,IPRT),NINCT(IT,IPRT),
     .                               NIXY(IT,IPRT),NTIN(IT,IPRT),
     .                               NTEN(IT,IPRT),NIFLG(IT,IPRT),
     .                               NPTC(IT,IPRT),NPTCM(IT,IPRT),
     .                               NSPZI(IT,IPRT),NSPZE(IT,IPRT),
     .                               NEMOD(IT,IPRT) !VK NPTCM
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
        IF (TRCINT) WRITE (iunout,*) 'CHGP,CHGEE,CHGEI,CHGMOM ',
     .                           CHGP,CHGEE,CHGEI,CHGMOM
C  READ ADDITIONAL DATA TO BE TRANSFERRED FROM B2 INTO EIRENE
C  HERE: B2 VOLUME TALLIES
        READ (IUNIN,'(3I6)') NAINB,NCOPIB,NCOPEB
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
csw 26jan2011 extra B25, not used anymore (called in eirene_mc)
c      if(my_pe == 0) then
c        call eirene_extrab25_eirpbls_init(2,0,0,0,0,0,nfla)
c      endif
csw
C
C  DEFINE ADDITIONAL TALLIES FOR COUPLING (UPDATED IN SUBR. UPTCOP
C                                              AND IN SUBR. COLLIDE)
      NCOPI=0
      IF (NMODE.GT.0) NCOPI=4
      IF (NCOPEB.NE.0) NCOPI=MAX(0,NCOPEB)
      NCPVI=NCOPI*NPLSI
      NCOP = NCOPI
C
C SAVE SOME MORE INPUT DATA FOR SHORT CYCLE ON COMMON CCOUPL
      NDX = NDXA
      NDY = NDYA
      NFL = NFLA
      NDXP = NDX+1
      NDYP = NDY+1
      LNLPLG=NLPLG
      LNLDRF=NLDRFT
      LTRCFL=TRCFLE
      NSTRI=NSTRAI
      DO 60 ISTRA=1,NSTRAI
        LNLVOL(ISTRA)=NLVOL(ISTRA)
60    CONTINUE
      NMODEI=NMODE
      NFILNN=NFILEN
C
      IF (NCPVI.EQ.0) GOTO 70
      DO IPLS=1,NPLSI
        ICPVE(IPLS)=1
        ICPRC(IPLS)=1
        TXTTAL(IPLS,NTALM)=
     .  'ENERGY WEIGHTED CX RATE OF ATOMS WITH IPLS                  '
        TXTSPC(IPLS,NTALM)=TEXTS(NSPAMI+IPLS)
        TXTUNT(IPLS,NTALM)='AMP                       '
C
        ICPVE(NPLSI+IPLS)=3
        ICPRC(NPLSI+IPLS)=1
        TXTTAL(NPLSI+IPLS,NTALM)=
     .  'PAR. MOM. SOURCE, FROM ATOMS, FOR IPLS             '
        TXTSPC(NPLSI+IPLS,NTALM)=TEXTS(NSPAMI+IPLS)
        TXTUNT(NPLSI+IPLS,NTALM)='G*CM/S* AMP * CM**-3       '
C
        ICPVE(2*NPLSI+IPLS)=3
        ICPRC(2*NPLSI+IPLS)=2
        TXTTAL(2*NPLSI+IPLS,NTALM)=
     .  'PAR. MOM. SOURCE, FROM MOLECULES, FOR IPLS         '
        TXTSPC(2*NPLSI+IPLS,NTALM)=TEXTS(NSPAMI+IPLS)
        TXTUNT(2*NPLSI+IPLS,NTALM)='G*CM/S* AMP * CM**-3       '
C
        ICPVE(3*NPLSI+IPLS)=3
        ICPRC(3*NPLSI+IPLS)=3
        TXTTAL(3*NPLSI+IPLS,NTALM)=
     .  'PAR. MOM. SOURCE, FROM TEST IONS, FOR IPLS         '
        TXTSPC(3*NPLSI+IPLS,NTALM)=TEXTS(NSPAMI+IPLS)
        TXTUNT(3*NPLSI+IPLS,NTALM)='G*CM/S* AMP * CM**-3       '
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
        ALLOCATE (PUX(NRAD))
        ALLOCATE (PUY(NRAD))
        ALLOCATE (PUXE(NRAD))
        ALLOCATE (PUYE(NRAD))
        ALLOCATE (PUXN(NRAD))
        ALLOCATE (PUYN(NRAD))
        ALLOCATE (PVX(NRAD))
        ALLOCATE (PVY(NRAD))
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
     .            PUX,PUY,PVX,PVY,MSHFRM)
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
csw 14apr2011      ALLOCATE(LCUT(0:NDXP))
!pb      DO IX=0,NDXP
      DO IX=0,NDXA+1
        LCUT(IX)=.FALSE.
        DO IPRT=1,NPLP-1
          IXI=NPOINT(2,IPRT)
          IXE=NPOINT(1,IPRT+1)
          IF (IX.GE.IXI.AND.IX.LT.IXE) THEN
            LCUT(IX)=.TRUE.
            WRITE (iunout,*) 'POLOIDAL CUT CELL INTRODUCED AT IP= ',IX
          ENDIF
        ENDDO
      ENDDO
      CALL EIRENE_LEER(1)

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
        IAISO = 0

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
!  CARTHESIAN PLANE
        DO IY=1,NDYA
          DO IX =1,NDXA
            IN=IY+(IX-1)*NR1ST
            ALE=ALPHXB(IX,IY)
            ALW=ALPHXB(IX-1,IY)
            IF (MAX(ALE,ALW)-MIN(ALE,ALW) > PIA) THEN
              AL=MIN(ALE,ALW)
              ALW=MAX(ALE,ALW)
              ALE=AL+PI2A
            END IF
            ALN=ALPHYB(IX,IY)
            ALS=ALPHYB(IX,IY-1)
            ALX=0.25D0*(ALE+ALW+ALN+ALS)
! cell centered
            PUX(IN)=COS(ALX)
            PUY(IN)=SIN(ALX)
            PVX(IN)=-PUY(IN)
            PVY(IN)=PUX(IN)
! surface centered
            PUXE(IN)=COS(ALE)
            PUYE(IN)=SIN(ALE)
            PUXN(IN)=COS(ALN)
            PUYN(IN)=SIN(ALN)
            PVXE(IN)=-SIN(ALE)
            PVYE(IN)=COS(ALE)
            PVXN(IN)=-SIN(ALN)
            PVYN(IN)=COS(ALN)
          END DO
        END DO
C
        DEALLOCATE (ALPHXB)
        DEALLOCATE (ALPHYB)
C
      ELSE
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) 
     .    ' NO FILE FORT.29 WITH MODIFIED GRID INFO. FOUND '
        WRITE (iunout,*) ' OLD VERSION CALCULATION MAGN. FIELD FROM ',
     .               'GRID IS USED '
        WRITE (iunout,*) ' GRID IS ASSUMED TO BE ORTHOGONAL '
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
      if (ntrfrm == 0) then
         READ(33,*) (XTRIAN(I),I=1,NRKNOT)
         READ(33,*) (YTRIAN(I),I=1,NRKNOT)
      else
         DO I=1,NRKNOT
           READ(33,*) J,XTRIAN(I),YTRIAN(I)
         ENDDO
      end if
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
        WRITE (IUNOUT,*) ' NUMBER OF TRIANGLES DO NOT MATCH '
        WRITE (IUNOUT,*) ' IN ELEMENTE AND NEIGHBOR FILES'
        WRITE (IUNOUT,*) ' PLEASE CHECK THE GEOMETRY '
        CALL EIRENE_EXIT_OWN(1)
      END IF

      DO I=1,NTRII
        READ(35,*) J,NCHBAR(1,I),NSEITE(1,I),IDUMMY,
     >               NCHBAR(2,I),NSEITE(2,I),IDUMMY,
     >               NCHBAR(3,I),NSEITE(3,I),IDUMMY,
C
     >               IXTRI(I),IYTRI(I)
C       WRITE(iunout,*) J,NECKE(1,J),NECKE(2,J),NECKE(3,J),
C    >               NCHBAR(1,J),NSEITE(1,J),
C    >               NCHBAR(2,J),NSEITE(2,J),NCHBAR(3,J),NSEITE(3,J)

C THE SPECIAL SURFACE PROPERTY (IF ANY) IS ON INMTI ARRAY, AND TRANSFERED INTO
C EIRENE VIA COMMON.
      ENDDO

      ALLOCATE (HEADS(N1ST,N2ND))
      DO IR=1,NR1ST
        DO IP=1,NP2ND
          NULLIFY(HEADS(IR,IP)%P)
        ENDDO
      ENDDO

C  FOR ALL QUADRANGLES BUILD LIST OF TRIANGLES BELONGING
C  TO THE QUADRANGLE
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

      IF (IO29.EQ.0) THEN
        DO ITRI=1,NTRII
          IY=IYTRI(ITRI)
          IX=IXTRI(ITRI)
          IF (IX .GT. 0) THEN
            IN=IY+(IX-1)*NR1STQ
            NSTGRD(ITRI)=ABS(XAISO(IX,IY)-1.)
          ENDIF
        ENDDO
        DEALLOCATE (XAISO)
        DEALLOCATE (IAISO)
      ENDIF
C
C
C  DETERMINE THE ARRAY INMTI FOR ALL NON DEFAULT STD. SURFACES
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
                IT=CURPOI%TRIANGLE
CVKG TO FIX A BUG WITH GEOMETRY
                ISC1=0
                ISC2=0
                DO IS=1,3
                  IF(EIRENE_POINT_ON_INTERVAL(XTRIAN(NECKE(IS,IT)),
     f                             YTRIAN(NECKE(IS,IT)),
     f                             XPOL(IR,IP),YPOL(IR,IP),
     f                             XPOL(IR,IP+1),YPOL(IR,IP+1))) THEN
                  IF(ISC1.GT.0) THEN
                    ISC2=IS 
                  ELSE
                    ISC1=IS
                  END IF
                END IF 
              ENDDO

C  NODES ISC1 AND ISC2 OF TRIANGLE IT ARE LOCATED ON RADIAL SURFACE IR
C  THAT MEANS  SIDE "NUMSI" OF TRIANGLE "IT" BELONGS TO NDS
              IF (ISC1.GT.0.AND.ISC2.GT.0) THEN
                NUMSI=MIN(ISC1,ISC2) 
                IF (NUMSI.EQ.1.AND.MAX(ISC1,ISC2).EQ.3) NUMSI=3
C
                  ICOG=ICOG+1
                  INSPAT(NUMSI,IT)=ICOG
                  INMTI(NUMSI,IT)=NLIM+ISTS
                  NBAR=NCHBAR(NUMSI,IT)
                  IF (NBAR.GT.0) THEN
                    NBARSI=NSEITE(NUMSI,IT)
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
                IT=CURPOI%TRIANGLE
CVKG TO FIX GEOMETRY BUG
                ISC1=0
                ISC2=0
                DO IS=1,3
                  IF(EIRENE_POINT_ON_INTERVAL(XTRIAN(NECKE(IS,IT)),
     f                             YTRIAN(NECKE(IS,IT)),
     f                             XPOL(IR,IP),YPOL(IR,IP),
     f                             XPOL(IR+1,IP),YPOL(IR+1,IP))) THEN
                  IF(ISC1.GT.0) THEN
                    ISC2=IS 
                  ELSE
                    ISC1=IS
                  END IF
                END IF 
              ENDDO

C  NODES ISC1 AND ISC2 OF TRIANGLE IT ARE LOCATED ON POLOIDAL SURFACE IP
C  THAT MEANS  SIDE "NUMSI" OF TRIANGLE "IT" BELONGS TO NDS
              IF (ISC1.GT.0.AND.ISC2.GT.0) THEN
                  NUMSI=MIN(ISC1,ISC2) 
                  IF (NUMSI.EQ.1.AND.MAX(ISC1,ISC2).EQ.3) NUMSI=3
C
                  ICOG=ICOG+1
                  INSPAT(NUMSI,IT)=ICOG
                  INMTI(NUMSI,IT)=NLIM+ISTS
                  NBAR=NCHBAR(NUMSI,IT)
                  IF (NBAR.GT.0) THEN
                    NBARSI=NSEITE(NUMSI,IT)
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
      ENDDO
C
C  NOW THE ADJUSTMENTS, WHICH ARE AUTOMATICALLY DONE IN GEOUSR OTHERWISE 
C  (INPUT BLOCK 15, B2-CODE SPECIFIC)
C
csw
      if(.true.) then
        call eirene_geousr
      else
      READ (IUNIN,*)
      READ (IUNIN,'(2I6)') NADMOD,NASMOD

      DO I=1,NADMOD
        READ (IUNIN,'(2I6,3E12.4)') NRS,IPUNKT,XCOOR,YCOOR,ZCOOR

        GOTO (1,2,3,4,5,6),IPUNKT
        WRITE (iunout,*) 'WRONG POINTNUMBER IN INFCOP '
        WRITE (iunout,*) 'INPUT LINE READING'
        WRITE (iunout,'(2I6,1P,3E12.4)') NRS,IPUNKT,XCOOR,YCOOR,ZCOOR
        WRITE (iunout,*) ' IS IGNORED '
        GOTO 10

    1   CONTINUE
        P1(1,NRS)=XCOOR
        P1(2,NRS)=YCOOR
        P1(3,NRS)=ZCOOR
        GOTO 10

    2   CONTINUE
        P2(1,NRS)=XCOOR
        P2(2,NRS)=YCOOR
        P2(3,NRS)=ZCOOR
        GOTO 10

    3   CONTINUE
        P3(1,NRS)=XCOOR
        P3(2,NRS)=YCOOR
        P3(3,NRS)=ZCOOR
        GOTO 10

    4   CONTINUE
        P4(1,NRS)=XCOOR
        P4(2,NRS)=YCOOR
        P4(3,NRS)=ZCOOR
        GOTO 10

    5   CONTINUE
        P5(1,NRS)=XCOOR
        P5(2,NRS)=YCOOR
        P5(3,NRS)=ZCOOR
        GOTO 10

    6   CONTINUE
        P6(1,NRS)=XCOOR
        P6(2,NRS)=YCOOR
        P6(3,NRS)=ZCOOR

   10   CONTINUE
      ENDDO

      DO I=1,NASMOD
        READ (IUNIN,'(5I6)') NAS,IPUNKT,NSSIR,NSSIP
        IF (IPUNKT.EQ.1) THEN
          P1(1,NAS)=XPOL(NSSIR,NSSIP)
          P1(2,NAS)=YPOL(NSSIR,NSSIP)
        ELSEIF (IPUNKT.EQ.2) THEN
          P2(1,NAS)=XPOL(NSSIR,NSSIP)
          P2(2,NAS)=YPOL(NSSIR,NSSIP)
        ELSE
          WRITE (iunout,*) 'WRONG POINTNUMBER IN GEOUSR '
          WRITE (iunout,*) 'INPUT LINE READING'
          WRITE (iunout,'(5I6)') NAS,IPUNKT,NSSIR,NSSIP
          WRITE (iunout,*) ' IS IGNORED '
        ENDIF
      ENDDO
csw
      endif
C
C  DETERMINE THE ARRAY INMTI FOR ALL ADDITIONAL SURFACES
C  ISTS=INMTI(ISIDE,NRCELL), ISIDE=1, 2, OR 3
C
      DO I=1,NLIMI
        IF (IGJUM0(I)==0) THEN
C  SURFACE I IS ACTIV
          IF (ILPLG(I).NE.0) THEN
C  SURFACE I IS PART IF A CONTOUR USED FOR THE MESHGENERATOR
            VSX=P2(1,I)-P1(1,I)
            VSY=P2(2,I)-P1(2,I)
            VS=SQRT(VSX**2+VSY**2)+EPS60
C
            DO 1111 IT=1,NTRII
              DO IS=1,3
                IF (IS.EQ.1) THEN
                  VTX=XTRIAN(NECKE(2,IT))-XTRIAN(NECKE(1,IT))
                  VTY=YTRIAN(NECKE(2,IT))-YTRIAN(NECKE(1,IT))
                  ISCS=1
                  ISC1=1
                  ISC2=2
                ELSEIF (IS.EQ.2) THEN
                  VTX=XTRIAN(NECKE(3,IT))-XTRIAN(NECKE(2,IT))
                  VTY=YTRIAN(NECKE(3,IT))-YTRIAN(NECKE(2,IT))
                  ISCS=2
                  ISC1=2
                  ISC2=3
                ELSEIF (IS.EQ.3) THEN
                  VTX=XTRIAN(NECKE(1,IT))-XTRIAN(NECKE(3,IT))
                  VTY=YTRIAN(NECKE(1,IT))-YTRIAN(NECKE(3,IT))
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
1112              TX=XTRIAN(NECKE(ISC,IT))
                  TY=YTRIAN(NECKE(ISC,IT))
csw
                  if(.true.) then
                    l1=eirene_point_on_interval(tx,ty,p1(1,i),p1(2,i),
     .                                                p2(1,i),p2(2,i))
                    TX=XTRIAN(NECKE(ISC2,IT))
                    TY=YTRIAN(NECKE(ISC2,IT))
                    l2=eirene_point_on_interval(tx,ty,p1(1,i),p1(2,i),
     .                                                p2(1,i),p2(2,i))
                    if(l1 .and. l2) then
                      IGJUM0(I)=1
                      ICOG=ICOG+1
                      INSPAT(ISCS,IT)=ICOG
                      INMTI(ISCS,IT)=I
                      IF (LCHKQUD)
     .                IREVERS(ISCS,IT)=INT(SIGN(1._DP,PX*VTRIX(ISCS,IT)+ 
     .                                          PY*VTRIY(ISCS,IT)))
                    endif
                  else
csw
                  IF (ABS(VSX).GT.ABS(VSY)) THEN
                    XMUE=(TX-PX)/VSX
                    IF (XMUE.GE.-1.D-5 .AND. XMUE.LE.1.+1.D-5) THEN
                      TEST=(TY-PY-XMUE*VSY)/VS
!pb                      IF (ABS(TEST).LT.1.D-5) THEN
                      IF (ABS(TEST).LT.1.D-4) THEN
                        IF (ICOU.EQ.2) THEN
C  TAKE CORRESPONDING ADDITIONAL SURFACE "I" OUT
C  AND REPLACE IT BY NON DEFAULT STD. SURFACE
                          IGJUM0(I)=1
                          ICOG=ICOG+1
                          INSPAT(ISCS,IT)=ICOG
                          INMTI(ISCS,IT)=I
                           IF (LCHKQUD)
     .                     IREVERS(ISCS,IT)=
     .                        INT(SIGN(1._DP,PX*VTRIX(ISCS,IT)+ 
     .                                       PY*VTRIY(ISCS,IT)))
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
C  AND REPLACE IT BY NON DEFAULT STD. SURFACE
                          IGJUM0(I)=1
                          ICOG=ICOG+1
                          INSPAT(ISCS,IT)=ICOG
                          INMTI(ISCS,IT)=I
                          IF (LCHKQUD) 
     .                    IREVERS(ISCS,IT)=
     .                        INT(SIGN(1._DP,PX*VTRIX(ISCS,IT)+ 
     .                                       PY*VTRIY(ISCS,IT)))
                           GOTO 1111
                        ELSE
                          ICOU=2
                          ISC=ISC2
                          GOTO 1112
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDIF
csw
                  endif
csw
                ENDIF
              ENDDO
1111        CONTINUE
          ENDIF
        ENDIF
      ENDDO

      DO IT=1,NTRII
        DO IS=1,3
          IF (NCHBAR(IS,IT).EQ.0.AND.INMTI(IS,IT).EQ.0) THEN
            WRITE (iunout,*) ' ERROR IN INFCOP '
            WRITE (iunout,*) ' OPEN SIDE OF TRIANGLE ',IT,' SIDE ',IS
            write (iunout,*) ' necke ',necke(1:3,it)
            write (iunout,*) ' xtrian,ytrian(1) ',xtrian(necke(1,it)),
     .                                       ytrian(necke(1,it))
            write (iunout,*) ' xtrian,ytrian(2) ',xtrian(necke(2,it)),
     .                                       ytrian(necke(2,it))
            write (iunout,*) ' xtrian,ytrian(3) ',xtrian(necke(3,it)),
     .                                       ytrian(necke(3,it))
            IS1=IS+1
            IF (IS.EQ.3) IS1=1
            WRITE (iunout,*) ' XTRIAN,YTRIAN ',XTRIAN(NECKE(IS,IT)),
     .                                    YTRIAN(NECKE(IS,IT))
            WRITE (iunout,*) ' XTRIAN,YTRIAN ',XTRIAN(NECKE(IS1,IT)),
     .                                    YTRIAN(NECKE(IS1,IT))
          ENDIF
        ENDDO
      ENDDO
C

csw 09jan2012 missing IGJUM3 correction!
CVK CORRECT  IGJUM1 AND IGJUM3 FOR TRIANGLES OUTSIDE STANDART (B2) MESH
CVK
      WRITE(iunout,*) 'IF0COP: CORRECTING IGJUM3'
      DO J=1,NTRII
C SET FLAG FOR BOUNDARY CELLS
        IFBOUND=.FALSE.
        DO IS=1,3
         IT=NCHBAR(IS,J)
         IF(IT.GT.0) THEN
          IF(IYTRI(IT).LE.0.OR.IXTRI(IT).LE.0) THEN
           IFBOUND=.TRUE.
           EXIT
          END IF
         END IF
        END DO
        DO I=1,NLIM
C SWITCHING ON ADDITIONAL SURFACES CHECKING FOR TRIANGULAR GRID (OUTSIDE B2 MESH
         IF(IXTRI(J).LE.0.AND.IYTRI(J).LE.0.AND.IGJUM0(I).EQ.0) THEN
           IGJUM3(J,I)=0
          END IF
C SWITCHING ON ADDITIONAL SURFACES CHECKING FOR BOUNDARY CELLS OF B2 GRID
          IF(IFBOUND.AND.IGJUM0(I).EQ.0) THEN
           IGJUM3(J,I)=0
          END IF
        END DO
      END DO
CVK END
csw

      NLPLG=.FALSE.
      NLFEM=.TRUE.
      LEVGEO=4
      NR1STQ=NR1ST
      NP2NDQ=NP2ND
      NR1ST=NTRII+1
      NLPOL=.FALSE.
      NP2ND=1
      NR1TAL=NR1ST
      NP2TAL=NP2ND
      NT3TAL=NT3RD
      NSBOX_TAL=NR1TAL*NP2TAL*NT3TAL*NBMLT+NRADD
      NGITT = COUNT(INMTI(1:3,1:NTRII) .NE. 0)

      CALL EIRENE_LEER(2)
      CALL EIRENE_HEADNG(' CASE REDEFINED IN COUPLE_TRIA: ',32)
      WRITE (iunout,*) 'NLPLG,NLFEM ',NLPLG,NLFEM
      WRITE (iunout,*) 'NLPOL       ',NLPOL
      WRITE (iunout,*) 'NR1ST,NP2ND ',NR1ST,NP2ND
      WRITE (iunout,*) 'NR1TAL,NP2TAL,NSBOX_TAL ',
     .                  NR1TAL,NP2TAL,NSBOX_TAL
      CALL EIRENE_LEER(2)
CTRIG E

csw 08mar2013, switch over to save out tallies per stratum on file-system
csw necessary when collecting data from MPI nodes in if3cop in case of nprs < nstrai
csw
      if(nprs > 1 .and. nprs < nstrai) then
        nfilen=1
      endif
csw
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
      LSHORT=.FALSE.
      CALL EIRENE_LEER(1)
      WRITE (iunout,*) 'IF1COP CALLED '
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
      IF (NLPLAS) WRITE (iunout,*) 'PLASMA DATA EXPECTED ON BRAEIR'
      IF (.NOT.NLPLAS)
     .   WRITE (iunout,*) 'PLASMA DATA EXPECTED ON FORT.31'
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
      CALL EIRENE_PLASM (31,NDX2,NDYA,NFLA,NDX,NDY,NFL,UUDIAB)
      CALL EIRENE_PLASM (31,NDX2,NDYA,NFLA,NDX,NDY,NFL,VVDIAB)
      CALL EIRENE_PLASM (31,NDX2,NDYA,1,NDX,NDY,1,POB)
C
C  OPTIONAL ARRAYS: VOLB, BFELDB,FNIX_YB, FNIY_XB
      VOLB = 0.D0
      BFELDB = 0.D0
!pb  FNIX_YB, FNIY_XB were never really used  
!pb      FNIX_YB = 0.D0
!pb      FNIY_XB = 0.D0
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
 
c  now removed again: fluxes: y-fluxes across x-cell-faces, and vice versa
C  X-SURFACE MAY BE INCLINED, HENCE: IT MAY RECEIVE A Y-FLUX TOO
!      CALL PLASM (31,NDX2,NDYA,NFLA,NDX,NDY,NFL,FNIX_YB)
C  Y-SURFACE MAY BE INCLINED, HENCE: IT MAY RECEIVE A X-FLUX TOO
!      CALL PLASM (31,NDX2,NDYA,NFLA,NDX,NDY,NFL,FNIY_XB)

C
2100  CONTINUE

!pb  for the time being set default values
csw 26sep2011
c      delta_sheathxb=3.1
c      delta_sheathyb=3.1
csw 06jan2012 NO! not for direct comparison SOLSP4.3!!!
      delta_sheathxb=0.0
      delta_sheathyb=0.0
csw
C
C  NO INDEX MAPPING REQUIRED, IF NEW TIMESTEP ON SAME PLASMA
C  BRAEIR NOT MODIFIED SINCE LAST CALL TO IF1COP
      IF (NCUTB_SAVE.EQ.NCUTL) THEN
        WRITE (iunout,*) 'NO INDEX MAPPING DONE'
      ELSE
        WRITE (iunout,*) 'INDEX MAPPING DONE ', NCUTB_SAVE,NCUTB,NCUTL
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
C  NOW THE SURFACE CENTERED PARTICLE FLUXES
      CALL EIRENE_INDMAP (FNIXB,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (FNIYB,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .             NCUTB,NCUTL,NPOINT,NPLP)
C  distinct from B2: these velocities are now cell centered in b2.5
      CALL EIRENE_INDMAP (UUB,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (VVB,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (UPB,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .             NCUTB,NCUTL,NPOINT,NPLP)
C   same as in B2:  these ENERGY fluxes are surface centered
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

c  presumably:  POB  (electric potential ??) and BFELDB  are cell centered
      CALL EIRENE_INDMAP (POB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (VOLB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,
     .             NCUTB,NCUTL,NPOINT,NPLP)
      CALL EIRENE_INDMAP (BFELDB,DUMMY,NDX,NDY,1,NDXA,NDYA,1,
     .             NCUTB,NCUTL,NPOINT,NPLP)

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
c
!pb      CALL EIRENE_INDMAP (FNIX_YB,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
!pb     .             NCUTB,NCUTL,NPOINT,NPLP)
!pb      CALL EIRENE_INDMAP (FNIY_XB,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
!pb     .             NCUTB,NCUTL,NPOINT,NPLP)
C
2101  CONTINUE
C
C  INDICATE, THAT NOW BRAEIR CONTAINS DATA AFTER INDEX-MAPPING
      NCUTB_SAVE=NCUTL
C
C  RESET 2D ARRAYS ONTO 1D EIRENE ARRAYS, RESCALE TO EIRENE UNITS
C  AND CONVERT BRAAMS VECTORS INTO CARTHESIAN EIRENE VECTORS
C
C  UNITS CONVERSION FACTORS
      T=1./ELCHA
      V=1.E2
      VL=1.E6
CTRIG A
C  VACCUM DATA NEEDED FOR REGION OUTSIDE B2-MESH
      TVAC=0.02
      DVAC=1.D2
      VVAC=0.
      BVAC=1.
CTRIG E
      DO 2105 IPLS=1,NPLSI
        D(IPLS)=1.E-6*FCTE(IPLS)
        FL(IPLS)=ELCHA*FCTE(IPLS)
2105  CONTINUE
C
      BZINTF = 1._DP
      BFINTF = 1._DP
      DO ITRI=1,NTRII
        IY=IYTRI(ITRI)
        IX=IXTRI(ITRI)
        IF (IX .GT. 0) THEN
          IN=IY+(IX-1)*NR1STQ
          TEINTF(ITRI)=TEB(IX,IY)*T
C
C  ONLY ONE ION TEMPERATURE AVAILABLE FROM PLASMA FLUID CODE,
C  SEE LOOP 2150 BELOW
          TIINTF(1,ITRI)=TIB(IX,IY)*T
C
CDR Construct magnetic field from:
C    a) poloidal field direction is given by that of the poloidal cell face PU..(in),
C      (PU(...) is cell centered)
C       and in the direction of increasing poloidal B2 cell index
c    b) toroidal field is in eirene positive z-direction (periodic cylinder, nltrz-option)
c                                (or positive 3rd coodinate "phi", in case nltra-option)
c    modulus of the ratio poloidal to poloidal field is given by the B2-array pitch RRB
C    magnitude of B-field is given by B2-array BFELDB

C  polodial field
          BX=PUX(IN)*RRB(IX,IY)   ! +PVX(IN)*0., but radial field is zero
          BY=PUY(IN)*RRB(IX,IY)   ! +PVY(IN)*0.
c  toroidal field
          BZ=SQRT(1.-RRB(IX,IY)**2)
c  normalize B-field vector to length 1 (one)
c  and apply input flagts for b-field orientation
          BN=SQRT(BX*BX+BY*BY+BZ*BZ)
          BXINTF(ITRI)=BX/BN*bpol
          BYINTF(ITRI)=BY/BN*bpol
          BZINTF(ITRI)=BZ/BN*btor
          BFINTF(ITRI)=BFELDB(IX,IY)
c
          VLINTF(ITRI)=VOLB(IX,IY)*VL
        ELSE
c  outside the b2 grid: set default b-field: (0,0,1)
          TEINTF(ITRI)=TVAC
          TIINTF(1,ITRI)=TVAC
          BXINTF(ITRI)=0.
          BYINTF(ITRI)=0.
          BZINTF(ITRI)=1.
          BFINTF(ITRI)=1.
c
C         VLINTF(ITRI)=1.
        ENDIF
      ENDDO
C
C  SET SAME ION TEMPERATURE FOR ALL EIRENE BACKGROUND SPECIES
C

      DO 2150 IPLS=1,NPLSTI
      DO 2150 ITRI=1,NTRII
        TIINTF(IPLS,ITRI)=TIINTF(1,ITRI)
2150  CONTINUE
C
CDR  set density from B2 array DNIB, for each fluid
CDR  set plasma flow velocity field from B2 arrays UPB (parallel velocity)
c  without drifts:
c  upb * pitch:  poloidal velocity (i.e. carthesian x,y direction).
c  poloidal field direction is given by that of the poloidal cell face PU..(in),
C  i.e. along a flux surface. (PU(...) is cell centered)
c  and upb*(1-pitch^2): toroidal velocity  (i.e. carthesian z direction (nltrz) or
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
        IF (IFLB(IPLS).GT.0) THEN
          IPLSV=MPLSV(IPLS)
          DO 2201 IFL=1,NFLA
            IF (IFLB(IPLS).NE.IFL) GOTO 2201
            DO ITRI=1,NTRII
              IY=IYTRI(ITRI)
              IX=IXTRI(ITRI)
              IF (IX .GT. 0) THEN
                IN=IY+(IX-1)*NR1STQ
                DIINTF(IPLS,ITRI)=DNIB(IX,IY,IFL)*D(IPLS)
cdr  parallel velocity, cell centered
                UPBC=UPB(IX,IY,IFL)
cdr  diamagnetic velocity, i.e. in (B x grad-PSI) direction.
cdr  take grad PSI ("radial") to be in direction of B2 iy grid.
cdr  unclear: orientation
                UDBC=UUDIAB(IX,IY,IFL)
cdr  radial velocity
                VVBC=0.5*(VVB(IX,IY-1,IFL)+VVB(IX,IY,IFL))
cdr  ???  perhaps a radial component of drift velocities?  not used any further currently
                VDBC=VVDIAB(IX,IY,IFL)
cdr  pitch  B_pol/B_total
                RRBC=RRB(IX,IY)
cdr  now set carthesina flow velocity components
                VXINTF(IPLSV,ITRI)=
     &           (PUX(IN)*(UPBC*RRBC-UDBC*SQRT(1.-RRBC**2))+
     &            PVX(IN)*VVBC)*V
                VYINTF(IPLSV,ITRI)=
     &           (PUY(IN)*(UPBC*RRBC-UDBC*SQRT(1.-RRBC**2))+
     &            PVY(IN)*VVBC)*V
                VZINTF(IPLSV,ITRI)=(UPBC*SQRT(1.-RRBC**2)+UDBC*RRBC)*V
              ELSE
c  region outside  B2 grid
                DIINTF(IPLS,ITRI)=DVAC
                VXINTF(IPLSV,ITRI)=VVAC
                VYINTF(IPLSV,ITRI)=VVAC
                VZINTF(IPLSV,ITRI)=VVAC
              ENDIF
            ENDDO
C  EIRENE BACKGROUND SPECIES "IPLS" IS NOW FILLED WITH B2 DATA "IFL"
2201      CONTINUE

C  NO DATA FOR "IPLS" IN B2 FILES
        ELSEIF (IFLB(IPLS).EQ.-13) THEN
C  READ DATA FOR "IPLS" FROM EIRENE DUMP FILE FT13
csw 09jan2012 NO!!!
csw          IF (IREAD.EQ.0) THEN
csw            OPEN (UNIT=13,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
csw            REWIND 13
csw            READ (13,IOSTAT=IO) TEIN,TIIN,DEIN,DIIN,VXIN,VYIN,VZIN
csw            IREAD=1
csw            IF (TRCFLE) WRITE (iunout,*) 'READ 13: RCMUSR, IO= ',IO
csw            CLOSE (UNIT=13)
csw          ENDIF
csw          IF (IO.EQ.0) THEN

            IPLSTI = MPLSTI(IPLS)
            IPLSV = MPLSV(IPLS)
            DO ITRI=1,NTRII
              DIINTF(IPLS,ITRI)=DIIN(IPLS,ITRI)
              VXINTF(IPLSV,ITRI)=VXIN(IPLSV,ITRI)
              VYINTF(IPLSV,ITRI)=VYIN(IPLSV,ITRI)
              VZINTF(IPLSV,ITRI)=VZIN(IPLSV,ITRI)
              TIINTF(IPLSTI,ITRI)=TIIN(IPLSTI,ITRI)
            ENDDO
csw          ENDIF
        ELSE
C  SET PARAMETERS FOR SPECIES IPLS TO ZERO
C  NOTHING TO BE DONE HERE
        ENDIF
2200  CONTINUE
C  B2-BRAAMS CODE SPECIFIC END
C
C
C  READ OTHER B2 ARRAYS INTO EIRENE, FOR PRINTOUT AND PLOTTING
C
      DO 2300 IAIN=1,NAINB
        IF (NAINT(IAIN).EQ.1.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO 2321 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=DNIB(IXTRI(IN),IYTRI(IN),NAINS(IAIN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
2321      CONTINUE
        ELSEIF (NAINT(IAIN).EQ.2.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO 2322 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=UUB(IXTRI(IN),IYTRI(IN),NAINS(IAIN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
2322      CONTINUE
        ELSEIF (NAINT(IAIN).EQ.3.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO 2323 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=VVB(IXTRI(IN),IYTRI(IN),NAINS(IAIN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
2323      CONTINUE
        ELSEIF (NAINT(IAIN).EQ.6) THEN
          DO 2326 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=PRB(IXTRI(IN),IYTRI(IN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
2326      CONTINUE
        ELSEIF (NAINT(IAIN).EQ.7.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO 2327 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=UPB(IXTRI(IN),IYTRI(IN),NAINS(IAIN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
2327      CONTINUE
        ELSEIF (NAINT(IAIN).EQ.8) THEN
          DO 2328 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=RRB(IXTRI(IN),IYTRI(IN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
2328      CONTINUE
        ELSEIF (NAINT(IAIN).EQ.9.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO 2329 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
!pb              ADINTF(IAIN,IN)=FNIXB(IXTRI(IN),IYTRI(IN),NAINS(IAIN))+
!pb     .                      FNIX_YB(IXTRI(IN),IYTRI(IN),NAINS(IAIN))
              ADINTF(IAIN,IN)=FNIXB(IXTRI(IN),IYTRI(IN),NAINS(IAIN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
2329      CONTINUE
        ELSEIF (NAINT(IAIN).EQ.10.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO 2330 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
!pb              ADINTF(IAIN,IN)=FNIYB(IXTRI(IN),IYTRI(IN),NAINS(IAIN))+
!pb     .                      FNIY_XB(IXTRI(IN),IYTRI(IN),NAINS(IAIN))
              ADINTF(IAIN,IN)=FNIYB(IXTRI(IN),IYTRI(IN),NAINS(IAIN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
2330      CONTINUE
        ELSEIF (NAINT(IAIN).EQ.11) THEN
          DO 2331 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=FEIXB(IXTRI(IN),IYTRI(IN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
2331      CONTINUE
        ELSEIF (NAINT(IAIN).EQ.12) THEN
          DO 2332 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=FEIYB(IXTRI(IN),IYTRI(IN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
2332      CONTINUE
        ELSEIF (NAINT(IAIN).EQ.13) THEN
          DO 2333 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=FEEXB(IXTRI(IN),IYTRI(IN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
2333      CONTINUE
        ELSEIF (NAINT(IAIN).EQ.14) THEN
          DO 2334 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=FEEYB(IXTRI(IN),IYTRI(IN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
2334      CONTINUE
        ELSEIF (NAINT(IAIN).EQ.15.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=UUDIAB(IXTRI(IN),IYTRI(IN),NAINS(IAIN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
          ENDDO
        ELSEIF (NAINT(IAIN).EQ.16.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=VVDIAB(IXTRI(IN),IYTRI(IN),NAINS(IAIN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
          ENDDO
        ELSEIF (NAINT(IAIN).EQ.17) THEN
          DO 2335 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=VOLB(IXTRI(IN),IYTRI(IN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
2335      CONTINUE
        ELSEIF (NAINT(IAIN).EQ.18) THEN
          DO 2336 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=BFELDB(IXTRI(IN),IYTRI(IN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
2336      CONTINUE
        ENDIF
2300  CONTINUE
C

csw 04dec2014 collecting normals for B2.5/B2 cells per triangle
      write(iunout,*) 'IF0COP: collecting normals'
      if(allocated(plnxtri)) then
        deallocate(plnxtri, plnytri, pplnxtri, pplnytri)
      endif
      allocate(plnxtri(ntrii))
      allocate(plnytri(ntrii))
      allocate(pplnxtri(ntrii))
      allocate(pplnytri(ntrii))
      do it=1,ntrii
        ix=ixtri(it)
        iy=iytri(it)
        ir = iy
        ip = ix
        if (ix<=0 .or. iy<=0) cycle

!  set radial unit vector  e_r:
!  set poloidal unit vector from underlying polygon grid, then take normal
!  to that vector to be the radial unit vector
        dx=xpol(ir,ip+1) - xpol(ir,ip)
        dy=ypol(ir,ip+1) - ypol(ir,ip)
        dd=sqrt(dx**2 + dy**2)+1.d-60
c  normal to that vector dx,dy. But: orientation ???
        plnxtri(it) = dy/dd
        plnytri(it) =-dx/dd

!  set poloidal unit vector  e_p:
        dx=xpol(ir+1,ip) - xpol(ir,ip)
        dy=ypol(ir+1,ip) - ypol(ir,ip)
        dd=sqrt(dx**2 + dy**2)+1.d-60
        if(dy .lt. 1.d-12) then
          pplnxtri(it) = 0.d0
          pplnytri(it) = 1.d0
        else
          pplnxtri(it) = 1.d0
          pplnytri(it) =-dx/dy
        endif
        dd=(pplnxtri(it)**2 + pplnytri(it)**2)+1.d-60
        pplnxtri(it)=pplnxtri(it)/dd
        pplnytri(it)=pplnytri(it)/dd
      enddo !it
C
!pb 22.01.2014 change of source strength for gas puff as required
!pb            from SOLPS moved here from subroutine EIRENE_MCARLO
csw
csw 24oct2011
      if(my_pe==0) then
        if (.not.allocated(flux_save)) then
           allocate(flux_save(nstra))
           flux_save = 0._dp
        end if
        DO ISTRAI=NTARGI+1,NSTRAI
          ISTRA = ISTRAI
csw 20mar2013
csw          if(.not.nlvol(istra)) then
          if(.not.nlvol(istra) .and. .not. nlcns(istra)) then
            IF (FLUX_save(ISTRA).NE.0.d0 ) THEN
              FLUX(ISTRA)=FLUX_save(ISTRA)*ELCHA
!pb 22.01.2014 Do NOT overwrite sources which are not mentioned in 
!pb            SOLPS input

!pb            ELSE
!pb              FLUX(ISTRA)=1.
            ENDIF
          endif
        ENDDO
      endif
csw 


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
      ALLOCATE (NUMSID(NGITT))
      ALLOCATE (NUMTRI(NGITT))
      ALLOCATE (TORL(NSTRA,NGITT))
      ALLOCATE (ESHT(NSTEP,NGITT))
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
          WRITE (iunout,*) 'ITARG,IPRT,NPBS,NPBC,NPES,NPEC ',
     .                 ITARG,IPRT,NPBS,NPBC,NPES,NPEC
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
              do is = 1, 3
!  test if side IS of triangle IT is parallel to B2 cell face
                par = vtrix(is,it)*dypol-vtriy(is,it)*dxpol
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
     .           (YANF-YTRIAN(NECKE(IS,ITRI)))**2). LT. 5*EPS5) THEN
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
          IG=IG+1
          IF (IG.GT.NGITT) GOTO 999
C  TESTEP, TISTEP: ZONE CENTERED TEMPERATURE IN BOUNDARY ZONE (EV)
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
!pb FACTOR FOR SHEATH POTENTIAL
          SHSTEP(ITARG,IG)=DELTA_SHEATHXB(NPBS,IY)
C  TORL: TOROIDAL LENGTH (CM) AT TARGET SEGMENT IY: CENTER OF GRAVITY
          TORL(ITARG,IG)=2.*PIA*0.5*(XTRIAN(NECKE(IS,ITRI))+
     .                               XTRIAN(NECKE(IS1,ITRI)))
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
              DELY=SQRT((XPOL(IY+1,NPES)-XPOL(IY,NPES))**2+
     .                  (YPOL(IY+1,NPES)-YPOL(IY,NPES))**2)
              FLSTEP(IPLS,ITARG,IG)=0.
              IF (DELY.GT.0.) THEN
                FLSTEP(IPLS,ITARG,IG)=MAX(0._DP,ORI(ITARG,IG)*
     .                                FNIXB(NPBS,IY,IFL))*FL(IPLS)/DELY
C  CORRECT FOR INCLINED TARGETS: ADD FLUXES FROM SECOND DIRECTION
C  USE SIGN FROM "MAIN" CONTRIBUTION TO DECIDE ORIENTATION OF SEC. CONTR.
!pb                IF (FLSTEP(IPLS,ITARG,IG).GT.0.) THEN
!pb                  FLSTEP(IPLS,ITARG,IG)=
!pb     .            FLSTEP(IPLS,ITARG,IG)+ABS(FNIX_YB(NPBS,IY,IFL))*
!pb     .                   FL(IPLS)/DELY
!pb                ENDIF
C  SET ION ENERGY FLUXES FROM B2-BOUNDARY CONDITIONS
                delti_para=3
                delte_para=0.5
                delti_perp=2
                delte_perp=0
!pb 20.09.2011 moved up
!pb csw 26sep2011
!pb                delta_sheathxb=3.1
!pb                delta_sheathyb=3.1
!pb csw
!  only one of the next two is different from 0
!                delti_para=deltai_parxb(npbs,iy)
!                delti_perp=deltai_radxb(npbs,iy)
!  only one of the next two is different from 0
!                delte_para=deltae_parxb(npbs,iy)
!                delte_perp=deltae_radxb(npbs,iy)
                tis=TISTEP(IPLSTI,ITARG,IG)
                tes=TESTEP(ITARG,IG)
!pb                ELSTEP(IPLS,ITARG,IG) = ELSTEP(IPLS,ITARG,IG)+
!pb     .                                  FL(IPLS)/DELY*
!pb     .            (TIS*delti_perp*ABS(Fnix_yb(npbs,iy,ifl))+
!pb     .             TIS*delti_para*ABS(fnixb  (npbs,iy,ifl))+
!pb     .             TES*delte_para*ABS(fnixb  (npbs,iy,ifl)))
!pb                ELSTEP(IPLS,ITARG,IG) = ELSTEP(IPLS,ITARG,IG)+
!pb     .                                  FL(IPLS)/DELY*
!pb     .            (TIS*delti_perp*ABS(Fnix_yb(npbs,iy,ifl))+
!pb     .             TIS*delti_para*ABS(fnixb  (npbs,iy,ifl))+
!pb     .             TES*delte_para*ABS(fnixb  (npbs,iy,ifl)))
                ELSTEP(IPLS,ITARG,IG) = ELSTEP(IPLS,ITARG,IG)+
     .                                  FL(IPLS)/DELY*
     .            (TIS*delti_para*ABS(fnixb  (npbs,iy,ifl))+
     .             TES*delte_para*ABS(fnixb  (npbs,iy,ifl)))
              ENDIF
            ENDIF
C  VXSTEP,VYSTEP,VZSTEP: SURFACE CENTERED FLOW VELOCITY (CM/S)
C  NOTE: PV VECTOR IS CELL CENTERED, BUT EXACT VECTOR CAN BE FOUND FROM
C        DATA FOR POLOIDAL POLYGON NPES
            IN=IY+(NPEC-1)*NR1STQ
            PVXS=XPOL(IY+1,NPES)-XPOL(IY,NPES)
            PVYS=YPOL(IY+1,NPES)-YPOL(IY,NPES)
            PNORM=SQRT(PVXS**2+PVYS**2)
            PVXS=PVXS/(PNORM+EPS60)
            PVYS=PVYS/(PNORM+EPS60)
C  ORTHONORMALIZE PU VECTOR WITH RESPECT TO PV
            PUPV=PUX(IN)*PVXS+PUY(IN)*PVYS
            PUXS=PUX(IN)-PUPV*PVXS
            PUYS=PUY(IN)-PUPV*PVYS
            PNORM=SQRT(PUXS**2+PUYS**2)
            PUXS=PUXS/(PNORM+EPS60)
            PUYS=PUYS/(PNORM+EPS60)
! Detlev und Xavier und Sven
!pb qq, qqc not needed/used
!pb            qq=qqc(npbs,iy)
!pb            uu=uub(npbs,iy,ifl)
!pb            uu=uub(npbc,iy,ifl)
!pb 13.9.2012  uub is a cell centered quantity like upb
!            RRBS=0.5*(RRB(NPBC,IY)+
!     .                RRB(NPBC-NINCT(ITARG,IPRT),IY))
            RRBS=RRB(NPBC-NINCT(ITARG,IPRT),IY)
!pb 21.9.2012            uu=uub(npbc,iy,ifl)
            uu=upb(npbc,iy,ifl)*rrbs
            up=upb(npbc,iy,ifl)
            ud=uudiab(npbc,iy,ifl)
            vv=vvb(npbs,iy,ifl)
            VXSTEP(IPLSV,ITARG,IG)=
     .            (PUXS*UU+PVXS*VV)*V
            VYSTEP(IPLSV,ITARG,IG)=
     .            (PUYS*UU+PVYS*VV)*V
            VZSTEP(IPLSV,ITARG,IG)=
     .            (SQRT(1.-RRBS**2)*UP+
     .             RRBS*UD)*V
!
! Detlev und Xavier
!            VXSTEP(IPLSV,ITARG,IG)=
!     .            (PUXS*UUB(NPBS,IY,IFL)+PVXS*VVB(NPBS,IY,IFL))*V
!            VYSTEP(IPLSV,ITARG,IG)=
!     .            (PUYS*UUB(NPBS,IY,IFL)+PVYS*VVB(NPBS,IY,IFL))*V
!            RRBS=0.5*(RRB(NPBC,IY)+
!     .                RRB(NPBC-NINCT(ITARG,IPRT),IY))
!cxpb            RRBS=UUB(NPBS,IY,IFL)/(UPB(NPBS,IY,IFL)+EPS60)
!            VZSTEP(IPLSV,ITARG,IG)=
!     .            (SQRT(1.-RRBS**2)*UPB(NPBC,IY,IFL)+
!     .             RRBS*UUDIAB(NPBC,IY,IFL))*V

! Detlevs Versuch?
!            UU=FNIXB(NPBS,IY,IFL)/DNIB(NPBC,IY,IFL)*V
!C           UP=
!C           VV=
!            VXSTEP(IPLSV,ITARG,IG)=
!     .            (PUXS*UU              +PVXS*VVB(NPBS,IY,IFL))*V
!            VYSTEP(IPLSV,ITARG,IG)=
!     .            (PUYS*UU              +PVYS*VVB(NPBS,IY,IFL))*V
!            RRBS=RRB(NPBC,IY)
!            UP=UU/(RRBS+eps60)
!            VZSTEP(IPLSV,ITARG,IG)=
!     .            (SQRT(1.-RRBS**2)*UP)*V
C Original
C
C            VXSTEP(IPLSV,ITARG,IG)=
C     .            (PUXS*UUB(NPBS,IY,IFL)+PVXS*VVB(NPBS,IY,IFL))*V
C            VYSTEP(IPLSV,ITARG,IG)=
C     .            (PUYS*UUB(NPBS,IY,IFL)+PVYS*VVB(NPBS,IY,IFL))*V
C            RRBS=UUB(NPBS,IY,IFL)/(UPB(NPBS,IY,IFL)+EPS60)
C            VZSTEP(IPLSV,ITARG,IG)=
C     .            (SQRT(1.-RRBS**2)*UPB(NPBS,IY,IFL))*V
3013      CONTINUE
        ENDDO
C
        GOTO 3030
C
3020    CONTINUE
C
C  SECOND: SOURCES AT RADIAL (X) SURFACES
C
        ITRI=0
        DO IX=NTIN(ITARG,IPRT),NTEN(ITARG,IPRT)-1
csw 11apr2011
          IF(LCUT(IX)) CYCLE
csw
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
csw 11apr2011
              dd = sqrt(dxpol*dxpol + dypol*dypol)
csw
              do is = 1, 3
!  test if side IS of triangle IT is parallel to B2 cell face
                par = vtrix(is,it)*dypol-vtriy(is,it)*dxpol
csw 11apr2011
                par=par/dd
csw
csw                if (abs(par) < 5*eps5) then
                if (abs(par) < 5.d-4) then
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
     .           (YANF-YTRIAN(NECKE(IS,ITRI)))**2). LT. 5*EPS5) THEN
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
          IG=IG+1
          IF (IG.GT.NGITT) GOTO 999
C  TESTEP, TISTEP: ZONE CENTERED TEMPERATURE IN BOUNDARY ZONE (EV)
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
!pb FACTOR FOR SHEATH POTENTIAL
          SHSTEP(ITARG,IG)=DELTA_SHEATHYB(IX,NPBS)
C  TORL: TOROIDAL LENGTH (CM) AT TARGET SEGMENT IY: CENTER OF GRAVITY
          TORL(ITARG,IG)=2.*PIA*0.5*(XTRIAN(NECKE(IS,ITRI))+
     .                               XTRIAN(NECKE(IS1,ITRI)))
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
              DELX=SQRT((XPOL(NPES,IX+1)-XPOL(NPES,IX))**2+
     .                  (YPOL(NPES,IX+1)-YPOL(NPES,IX))**2)
              FLSTEP(IPLS,ITARG,IG)=0.
              IF (DELX.GT.0.) THEN
                FLSTEP(IPLS,ITARG,IG)=MAX(0._DP,ORI(ITARG,IG)*
     .                            FNIYB(IX,NPBS,IFL))*FL(IPLS)/DELX
C  CORRECT FOR INCLINED TARGETS: ADD FLUXES FROM SECOND DIRECTION
!pb                IF (FLSTEP(IPLS,ITARG,IG).GT.0.) THEN
!pb                  FLSTEP(IPLS,ITARG,IG)=
!pb     .            FLSTEP(IPLS,ITARG,IG)+ABS(FNIY_XB(IX,NPBS,IFL))*
!pb     .                   FL(IPLS)/DELX
!pb                ENDIF
C  SET ION ENERGY FLUXES FROM B2 BOUNDARY CONDITIONS
                delti_para=3
                delte_para=0.5
                delti_perp=2
                delte_perp=0
!  only one of the next two is different from 0
!                delti_para=deltai_paryb(ix,npbs)
!                delti_perp=deltai_radyb(ix,npbs)
!  only one of the next two is different from 0
!                delte_para=deltae_paryb(ix,npbs)
!                delte_perp=deltae_radyb(ix,npbs)
                tis=TISTEP(IPLSTI,ITARG,IG)
                tes=TESTEP(ITARG,IG)
!pb                ELSTEP(IPLS,ITARG,IG) = ELSTEP(IPLS,ITARG,IG) +
!pb     .                                  FL(IPLS)/DELX*
!pb     .            (TIS*delti_perp*ABS(Fniyb  (ix,npbs,ifl))+
!pb     .             TIS*delti_para*ABS(fniy_xb(ix,npbs,ifl))+
!pb     .             TES*delte_para*ABS(fniy_xb(ix,npbs,ifl)))
!pb                ELSTEP(IPLS,ITARG,IG) = ELSTEP(IPLS,ITARG,IG) +
!pb     .                                  FL(IPLS)/DELX*
!pb     .            (TIS*delti_perp*ABS(Fniyb  (ix,npbs,ifl))+
!pb     .             TIS*delti_para*ABS(fniy_xb(ix,npbs,ifl))+
!pb     .             TES*delte_para*ABS(fniy_xb(ix,npbs,ifl)))
                ELSTEP(IPLS,ITARG,IG) = ELSTEP(IPLS,ITARG,IG) +
     .                                  FL(IPLS)/DELX*
     .              TIS*delti_perp*ABS(Fniyb(ix,npbs,ifl))
              ENDIF
            ENDIF
C  VXSTEP,VYSTEP,VZSTEP: SURFACE CENTERED FLOW VELOCITY (CM/S)
C  NOTE: PU VECTOR IS CELL CENTERED, BUT EXACT VECTOR CAN BE FOUND FROM
C        RADIAL POLYGON NPES DATA
            IN=NPEC+(IX-1)*NR1STQ
            PUXS=XPOL(NPES,IX+1)-XPOL(NPES,IX)
            PUYS=YPOL(NPES,IX+1)-YPOL(NPES,IX)
            PNORM=SQRT(PUXS**2+PUYS**2)
            PUXS=PUXS/(PNORM+EPS60)
            PUYS=PUYS/(PNORM+EPS60)
C  ORTHONORMALIZE PV VECTOR WITH RESPECT TO PU
            PUPV=PUXS*PVX(IN)+PUYS*PVY(IN)
            PVXS=PVX(IN)-PUPV*PUXS
            PVYS=PVY(IN)-PUPV*PUYS
            PNORM=SQRT(PVXS**2+PVYS**2)
            PVXS=PVXS/(PNORM+EPS60)
            PVYS=PVYS/(PNORM+EPS60)
! Detlev und Xavier
!22102012            VXSTEP(IPLSV,ITARG,IG)=
!22102012     .            (PUXS*UUB(IX,NPBS,IFL)+PVXS*VVB(IX,NPBS,IFL))*V
!22102012            VYSTEP(IPLSV,ITARG,IG)=
!22102012     .            (PUYS*UUB(IX,NPBS,IFL)+PVYS*VVB(IX,NPBS,IFL))*V
            VXSTEP(IPLSV,ITARG,IG)=
     .            (PUXS*UUB(IX,NPBC,IFL)+PVXS*VVB(IX,NPBS,IFL))*V
            VYSTEP(IPLSV,ITARG,IG)=
     .            (PUYS*UUB(IX,NPBC,IFL)+PVYS*VVB(IX,NPBS,IFL))*V
            RRBS=RRB(IX,NPBC)
!22102012            RRBS=0.5*(RRB(IX,NPBC)+
!22102012     .                RRB(IX,NPBC-NINCT(ITARG,IPRT)))
cxpb            RRBS=UUB(IX,NPBS,IFL)/(UPB(IX,NPBS,IFL)+EPS60)
            VZSTEP(IPLSV,ITARG,IG)=
     .            (SQRT(1.-RRBS**2)*UPB(IX,NPBC,IFL)+
     .             RRBS*UUDIAB(IX,NPBC,IFL))*V
C  SURFACE CENTERED VELOCITIES
!            VXSTEP(IPLSV,ITARG,IG)=
!     .            (PUXS*UUB(IX,NPBS,IFL)+PVXS*VVB(IX,NPBS,IFL))*V
!            VYSTEP(IPLSV,ITARG,IG)=
!     .            (PUYS*UUB(IX,NPBS,IFL)+PVYS*VVB(IX,NPBS,IFL))*V
!            RRBS=UUB(IX,NPBS,IFL)/(UPB(IX,NPBS,IFL)+EPS60)
!            VZSTEP(IPLSV,ITARG,IG)=
!     .            (SQRT(1.-RRBS**2)*UPB(IX,NPBS,IFL))*V
3023      CONTINUE
        ENDDO
3030    CONTINUE
C
3040  CONTINUE
      NRWL(ITARG)=IG+1

C
      IF (TRCSOU) CALL EIRENE_LEER(2)
C
C  INITIALIZE FUNCTION STEP (FOR RANDOM SAMPLING ALONG TARGET)
C  SET SOME SOURCE PARAMETERS EXPLICITLY TO ENFORCE INPUT CONSISTENCY
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
C  IN TRIA OPTION INDIM=4 (LEVGEO=3) CORRESPONDS TO INDIM=1 (LEVGEO=4)
      INDIM(1,ITARG)=1
CTRIG E
      IF (INDSRC(ITARG).NE.6) THEN
        I34=EIRENE_IDEZ(INT(SORLIM(1,ITARG)),3,3)
        SORLIM(1,ITARG)=I34*100+40
      ELSEIF (INDSRC(ITARG).EQ.6) THEN
C  SORLIM DEFAULT WAS 0.D0
        SORLIM(1,ITARG)=0240
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
        NMINPTS(ITARG)=NPTCM(ITARG,1)*MPTS_COMSOU !VK MINIMUM NUMBER OF HYSTORIES FOR THE STRATUM
        NINITL(ITARG)=ITARG*1001
!pb     NSPEZ(ITARG)=-1
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
3999  CONTINUE
C
C  TARGET DATA ARE DEFINED NOW
C
C
C  COMPUTE EXACT SURFACE ENERGY FLUXES FOR COMPARISON WITH SAMPLED
C  E-FLUX "ETOTP". THIS IS ONLY FOR DIAGNOSTICS PURPOSES
C  E.G. TO CHECK CONSISTENCY OF BOUNDARY CONDITIONS
C  STATEMENT NO. 6000 ---> 6500
C
      IF (.NOT.TRCSOU) GOTO 6500
C
      EEMAX=0.
      EESHT=0.
C
      DO 6011 IG=1,NRWL(ITARG)-1
        OR=ORI(ITARG,IG)
C
C  COMPUTE SHEATH POTENTIAL ESHT(ITARG,IG)
C  USE ALL NPLSI SPECIES, NOT JUST IFL=NSPZI,NSPZE
C
        ESHT(ITARG,IG)=0.D0
        NEM=IABS(NEMODS(ITARG))
        IF (NEM.EQ.3.OR.NEM.EQ.5.OR.NEM.EQ.7) THEN
          DO 6005 IPL=1,NPLSI
            IPLV=MPLSV(IPL)
            VPX=VXSTEP(IPLV,ITARG,IG)
            VPY=VYSTEP(IPLV,ITARG,IG)
            VPZ=VZSTEP(IPLV,ITARG,IG)
            VP(IPL)=SQRT(VPX**2+VPY**2+VPZ**2)
            DI(IPL)=DISTEP(IPL,ITARG,IG)
6005      CONTINUE
          TE=TESTEP(ITARG,IG)
          CUR=0.
          GAMMA=0.
          ESHT(ITARG,IG)=EIRENE_SHEATH(TE,DI,VP,NCHRGP,GAMMA,CUR,
     .                         NPLSI,-ITARG)
        ELSE IF (NEM == 9) THEN
          IF (IGSTEP(ITARG,IG).GT.200000) THEN
            NPES=IGSTEP(ITARG,IG)-200000
            NPBS=NPES-1
            TE=TESTEP(ITARG,IG)
            ESHT(ITARG,IG)=DELTA_SHEATHXB(NPBS,IY)*TE
          ELSEIF (IGSTEP(ITARG,IG).LT.200000) THEN
            NPES=IGSTEP(ITARG,IG)-100000
            NPBS=NPES-1
            TE=TESTEP(ITARG,IG)
            ESHT(ITARG,IG)=DELTA_SHEATHYB(IX,NPBS)*TE
          ENDIF
        ENDIF
C

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
            WRITE (iunout,*) 'IPLS,ITARG,IG,MACH ',IPLS,ITARG,IG,VTEST,
     .                                             VTEST2
C           WRITE (iunout,*) 'POL., TOR., RAD. ',PM1,VPZ,VR
            CALL EIRENE_LEER(1)
          END IF
C
C  BOHM CRITERION CHECK DONE
C
C  NEXT: TARGET ENERGY FLUXES
          DRR=RRSTEP(ITARG,IG+1)-RRSTEP(ITARG,IG)
C  ENERGY FLUX DEFINED IN INPUT BLOCK 7
          IF (NEMODS(ITARG).EQ.1) THEN
            EADD=SORENI(ITARG)
            ESUM=EADD*FLSTEP(IPLS,ITARG,IG)
            EEMAX=EEMAX+ESUM*DRR
          ELSEIF (NEMODS(ITARG).EQ.2.OR.NEMODS(ITARG).EQ.3) THEN
            EADD=SORENI(ITARG)*TISTEP(IPLSTI,ITARG,IG)+SORENE(ITARG)*
     .           TESTEP(ITARG,IG)
            ESUM=EADD*FLSTEP(IPLS,ITARG,IG)
            EEMAX=EEMAX+ESUM*DRR
          ELSEIF (NEMODS(ITARG).GE.4 .AND. NEMODS(ITARG).LE.7) THEN
            PERWI=PERW/SQRT(BMASS(IPLS)/RMASSP(IPLS))
            PARWI=PARW/SQRT(BMASS(IPLS)/RMASSP(IPLS))
            EADD=EIRENE_EMAXW(TISTEP(IPLSTI,ITARG,IG),PERWI,PARWI)
            ESUM=EADD*FLSTEP(IPLS,ITARG,IG)
            EEMAX=EEMAX+ESUM*DRR
          ELSEIF (NEMODS(ITARG).EQ.8 .OR. NEMODS(ITARG).EQ.9) THEN
C  ENERGY FLUX DEFINED BY B2-BOUNDARY CONDITIONS
            EEMAX=EEMAX+ELSTEP(IPLS,ITARG,IG)*DRR
          ENDIF

C  ENERGY GAIN BY SHEATH ACCELERATION
          EADD=NCHRGP(IPLS)*ESHT(ITARG,IG)
          ESUM=EADD*FLSTEP(IPLS,ITARG,IG)
          EESHT=EESHT+ESUM*DRR
6009    CONTINUE
6011  CONTINUE
C
      CALL EIRENE_LEER(1)
      WRITE (iunout,*) 'TARGET DATA: TARGET NO. ITARG=ISTRA= ',ITARG
      WRITE (iunout,*) TXTSOU(ISTRA)
      WRITE (iunout,'(1X,A3,9A11,A3)') 
     .'IG','ARC','P-FLUX','E-FLUX','TE','TI','SHEATH/TE',
     . 'VXSTEP','VYSTEP','VZSTEP'
      DO 6100 IG=1,NRWL(ITARG)-1

        IF (IGSTEP(ITARG,IG).GT.200000) THEN
          IF (ORI(ITARG,IG).LT.0) NSEW='E'
          IF (ORI(ITARG,IG).GT.0) NSEW='W'
        ENDIF
        IF (IGSTEP(ITARG,IG).LT.200000) THEN
          IF (ORI(ITARG,IG).LT.0) NSEW='S'
          IF (ORI(ITARG,IG).GT.0) NSEW='N'
        ENDIF
        WRITE (iunout,'(1X,I3,1P,9E11.3,A3)')
     .             IG,RRSTEP(ITARG,IG),FLSTEP(0,ITARG,IG),
     .             ELSTEP(0,ITARG,IG),
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

      DEALLOCATE (NUMSID)
      DEALLOCATE (NUMTRI)
      DEALLOCATE (TORL)
      DEALLOCATE (ESHT)
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
      NDXY=(NDXA-1)*NR1STQ+NDYA
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
      NDXY=(NDXA-1)*NR1STQ+NDYA
C
99992 CONTINUE
csw 14jul2011
      if(my_pe == 0)then
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
        CALL EIRENE_ALLOC_EIRBRA(NDX,NDY,NFL,NSTRA,IFOFF)
        RESSNI = 0._DP
        RESSMO = 0._DP
        RESSEE = 0._DP
        RESSEI = 0._DP
      END IF
C
      IF (NEW_ITER == 0) THEN
        RESSNI = 0.D0
        RESSMO = 0.D0
        RESSEE = 0.D0
        RESSEI = 0.D0
      ENDIF
C
      IF (.NOT.LSHORT) THEN
        RESSNI(ISTRAA:ISTRAE,:) = 0._DP
        RESSMO(ISTRAA:ISTRAE,:) = 0._DP
        RESSEE(ISTRAA:ISTRAE) = 0._DP
        RESSEI(ISTRAA:ISTRAE) = 0._DP
      ENDIF

      volSUMN(ISTRAA:ISTRAE)=0.0           !dpc
      volSUMM(ISTRAA:ISTRAE)=0.0           !dpc
      volSUMEI(ISTRAA:ISTRAE)=0.0          !dpc
      volSUMEE(ISTRAA:ISTRAE)=0.0          !dpc
csw
      endif


csw 08mar2013, make sure that all strata have written out onto fort.10, ie.
csw            all cpu's have come to this point
      if(nprs > 1 .and. nprs < nstrai .and. .not.lshort) then
        call mpi_barrier(MPI_COMM_WORLD,ier)
      endif
csw
csw 26jan2011 extra B25
      if(my_pe == 0) then
        call eirene_extraB25_wneuinit
csw 20jun2013 bugfix, check for nprs>1
        if(nprs > 1) then
          call eirene_extraB25_wneuclean
        endif
      endif
csw
      DO 10000 ISTRAI=ISTRAA,ISTRAE
C
csw 20mar2013
        IF (XMCP(ISTRAI).LE.1.) then
         sni(:,:,:,istrai) = 0.d0
         smo(:,:,:,istrai) = 0.d0
         see(:,:,istrai) = 0.d0
         sei(:,:,istrai) = 0.d0
         GOTO 10000
        endif
csw

csw 14jul2011
csw 08mar2013 check for nprs >= nstrai too
        IF (nprs > 1 .and. nprs >= nstrai .and. .not.lshort) then
          irnk = npesta(istrai)
          if(irnk >= nprs) then
            goto 10000

          elseif(my_pe ==0 .and.my_pe == irnk) then
            if(xmcp(istrai).le.1) then
              goto 10000
            else
              goto 7000
            endif

          elseif(my_pe == 0.and.my_pe /= irnk) then
            allocate(save_estimv(nvoltl,nrtal))
            allocate(save_estims(nsrftl,nlmpgs))
!pb 08.01.2014            allocate(save_sigma_cop(ncpv_stat,nrtals))
            save_estimv=estimv
            save_estims=estims
!pb 08.01.2014            save_sigma_cop=sigma_cop
            call mpi_recv(estimv,nrtal*nvoltl,
     .               MPI_DOUBLE_PRECISION,
     .               irnk,istrai, MPI_COMM_WORLD,
     .               MPI_STATUS_IGNORE,ier)
            call mpi_recv(estims,nsrftl*nlmpgs,
     .               MPI_DOUBLE_PRECISION,
     .               irnk,istrai, MPI_COMM_WORLD,
     .               MPI_STATUS_IGNORE,ier)
!pb 08.01.2014            call mpi_recv(sigma_cop,ncpv_stat*nrtals,
!pb 08.01.2014     .               MPI_DOUBLE_PRECISION,
!pb 08.01.2014     .               irnk,istrai, MPI_COMM_WORLD,
!pb 08.01.2014     .               MPI_STATUS_IGNORE,ier)
            allocate(wtotp_dum(0:npls,0:nstra))
            call mpi_recv(wtotp_dum(0,0),(npls+1)*(nstra+1),
     .               MPI_DOUBLE_PRECISION,
     .               irnk,istrai, MPI_COMM_WORLD,
     .               MPI_STATUS_IGNORE,ier)
            wtotp(0:npls,istrai) = wtotp_dum(0:npls,istrai)
            deallocate(wtotp_dum)

          else 
            if(my_pe/=0.and.my_pe == irnk) then
              call mpi_send(estimv,nrtal*nvoltl,
     .                 MPI_DOUBLE_PRECISION,0,istrai,
     .                 MPI_COMM_WORLD,ier)
              call mpi_send(estims,nsrftl*nlmpgs,
     .                 MPI_DOUBLE_PRECISION,0,istrai,
     .                 MPI_COMM_WORLD,ier)
!pb 08.01.2014              call mpi_send(sigma_cop,ncpv_stat*nrtals,
!pb 08.01.2014     .                 MPI_DOUBLE_PRECISION,0,istrai,
!pb 08.01.2014     .                 MPI_COMM_WORLD,ier)
              call mpi_send(wtotp(0,0),(npls+1)*(nstra+1),
     .                 MPI_DOUBLE_PRECISION,0,istrai,
     .                 MPI_COMM_WORLD,ier)
            endif
            goto 10000
          endif
        endif

C
        IF (LSHORT) GOTO 7000
C
C  READ DATA FROM STRATUM NO. ISTRAI BACK INTO WORKING SPACE
C  IF REQUIRED
C
csw 15jul2011 mpi
csw 08mar2013 check for nprs < nstrai too
        if(nprs == 1 .or. nprs < nstrai) then
csw 08mar2013 do this only for my_pe=0
          if(my_pe /=0) goto 10000
csw
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
        endif
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
     .                     ISTRAI
          WRITE (iunout,*) 
     .       'NO DATA RETURNED TO PLASMA CODE FOR THIS STRATUM'
          GOTO 7999
        ELSEIF (ISTRAI.GT.NTARGI) THEN
          FLXI=1.
C  FLXEIR HAS TO BE RESET TO SCALE TO NEW SOURCE STRENGTH DURING SHORT CYCLE
C  IF THE SOURCE STRENGTH IS TO CHANGE DURING THE SHORT CYCLE (E.G.: VOL-REC)
          FLXEIR(ISTRAI)=1._DP
        ENDIF
C
C  FIRSTLY INITIALIZE SOURCE TERM ARRAYS
C
        DO 7100 IX=0,NDXA+1
          DO 7150 IY=0,NDYA+1
            SEE(IX,IY,ISTRAI)=0.
            SEI(IX,IY,ISTRAI)=0.
7150      CONTINUE
7100    CONTINUE
        DO 7210 IF=1,NFLA
          DO 7220 IX=0,NDXA+1
            DO 7230 IY=0,NDYA+1
              SNI(IX,IY,IF,ISTRAI)=0.
              SMO(IX,IY,IF,ISTRAI)=0.
7230        CONTINUE
7220      CONTINUE
7210    CONTINUE
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
C                        AND BULK ION CHARGE EXCHANGE WITH ATOMS
C

        PAPL=0.D0
        CPMUL => PAPLS(ISTRAI)%PMUL
        DO WHILE (ASSOCIATED(CPMUL))
          IPLS=CPMUL%IART
          IN=CPMUL%ICM
          PAPL(IPLS,IN)=CPMUL%VALUEM
          CPMUL => CPMUL%NXTMUL
        ENDDO

        EAEL=0.D0
        CPSIM => EAELS(ISTRAI)%PSIM
        DO WHILE (ASSOCIATED(CPSIM))
          IN=CPSIM%ICS
          EAEL(IN)=CPSIM%VALUES
          CPSIM => CPSIM%NXTSIM
        END DO

        EAPL=0.D0
        CPSIM => EAPLS(ISTRAI)%PSIM
        DO WHILE (ASSOCIATED(CPSIM))
          IN=CPSIM%ICS
          EAPL(IN)=CPSIM%VALUES
          CPSIM => CPSIM%NXTSIM
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
            EAPL(IN)=EAPL(IN)+CHI
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
        CPSIM => EIPLS(ISTRAI)%PSIM
        DO WHILE (ASSOCIATED(CPSIM))
          IN=CPSIM%ICS
          EIPL(IN)=CPSIM%VALUES
          CPSIM => CPSIM%NXTSIM
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
          END DO
          CHE=CPMUL%VALUEM *
     .        (SEENWI(IN,IION)-RTIS%SEEODI(IN,IION))*ELCHA
          EIEL(IN)=EIEL(IN)+CHE
          CHEEM(IN)=CHEEM(IN)+CHE
          CHI=CPMUL%VALUEM *
     .        (SEINWI(IN,IION)-RTIS%SEIODI(IN,IION))*ELCHA
          EIPL(IN)=EIPL(IN)+CHI
          CHEIM(IN)=CHEIM(IN)+CHI
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
        CPSIM => EMPLS(ISTRAI)%PSIM
        DO WHILE (ASSOCIATED(CPSIM))
          IN=CPSIM%ICS
          EMPL(IN)=CPSIM%VALUES
          CPSIM => CPSIM%NXTSIM
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
        CPSIM => EPPL_COPS(ISTRAI)%PSIM
        DO WHILE (ASSOCIATED(CPSIM))
          IN=CPSIM%ICS
          EPLODA(IN)=CPSIM%VALUES
          CPSIM => CPSIM%NXTSIM
        END DO

        EPEODA=0.D0
        CPSIM => EPELS(ISTRAI)%PSIM
        DO WHILE (ASSOCIATED(CPSIM))
          IN=CPSIM%ICS
          EPEODA(IN)=CPSIM%VALUES
          CPSIM => CPSIM%NXTSIM
        END DO
C
C  CORRECTION FOR ELECTRON IMPACT DISSOCIATION OF TEST IONS FINISHED
C
C  SHORT LOOP CORRECTION FINISHED
C
7400    CONTINUE
C
C
C  ADD CONTRIBUTIONS FROM VOLUME RECOMBINATION SOURCE
C
        PPPL_COP =0.D0
        CPPV = 0.D0
        EPPL_COP = 0.D0
        EPEL = 0.D0

        IF (NLVOL(ISTRAI)) THEN
C
          RECTOT = 0._DP
          DO 7473 IPLS=1,NPLSI
            CNDYNP=AMUA*RMASSP(IPLS)
            IPLSTI = MPLSTI(IPLS)
            DO 7472 IIRC=1,NPRCI(IPLS)
              IRRC=LGPRC(IPLS,IIRC)
              SUMN=0.0
              SUMM=0.0
              SUMEI=0.0
              SUMEE=0.0
              DO 7471 IN=1,NTRII
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
csw 21feb2012, taken from SOLPS4.3
C CORRECT PAAT FOR CALCULATING RADIATION  FOR ATOMS IN WNEUTRALS
                IAT=NATPRC(IRRC)
                IF(IAT.GT.0) PAAT(IAT,INC)=PAAT(IAT,INC)-RECADD 
csw
                PIADD=PARMOM(IPLS,IN)*RECADD
                CPPV(IPLS,INC)=CPPV(IPLS,INC)+PIADD
                SUMM=SUMM+PIADD*VOL(IN)
                EIADD=(1.5*TIIN(IPLSTI,IN)+EDRIFT(IPLS,IN))*RECADD
csw 08feb2013
csw                EPPL_COP(INC)=EPPL_COP(INC)+EIADD
                EAPL(INC)=EAPL(INC)+EIADD 
csw
                SUMEI=SUMEI+EIADD*VOL(IN)
csw 08feb2013
csw                EPEL(INC)=EPEL(INC)+EEADD
                EAEL(INC)=EAEL(INC)+EEADD 
csw
                SUMEE=SUMEE+EEADD*VOL(IN)
7471          CONTINUE
              RECTOT = RECTOT + SUMN
              WRITE (iunout,*) 'IPLS,IRRC ',IPLS,IRRC
              CALL EIRENE_MASR4('SUMN, SUMM, SUMEI, SUMEE        ',
     .                     SUMN,SUMM,SUMEI,SUMEE)
            volSUMN(ISTRAI)=volSUMN(ISTRAI)+SUMN             ! dpc
            volSUMM(ISTRAI)=volSUMM(ISTRAI)+SUMM             ! dpc
            volSUMEI(ISTRAI)=volSUMEI(ISTRAI)+SUMEI          ! dpc
            volSUMEE(ISTRAI)=volSUMEE(ISTRAI)+SUMEE          ! dpc
7472        CONTINUE
7473      CONTINUE
cdr
CC SUMN_OLD=WTOTP ???
          IF (.NOT.LSHORT) SUMN_OLD=RECTOT
csw 08jun2010
          sumn_old=rectot
csw
C  RESCALE PPPL_COP,.... TO SOURCE STRENGTH FROM LAST FULL EIRENE RUN
C  BECAUSE ALSO PAPL,.... ARE SCALED LIKE THIS
          PPPL_COP=PPPL_COP*SUMN_OLD/RECTOT
          CPPV=CPPV*SUMN_OLD/RECTOT
          EPPL_COP=EPPL_COP*SUMN_OLD/RECTOT
          EPEL=EPEL*SUMN_OLD/RECTOT

C  RESCALE SOURCES FOR B2 (PPPL_COP,PAPL,....) TO NEW SOURCE STRENGTH
          FLXEIR(ISTRAI)=RECTOT/SUMN_OLD
cdr
          IF (LSHORT) THEN
            DO IPLS=1,NPLSI
              CHPM(IPLS,1:NSBOX_TAL) = CHPM(IPLS,1:NSBOX_TAL) +
     .             PPPL_COP(IPLS,1:NSBOX_TAL) - PPLODA(IPLS,1:NSBOX_TAL) 
!pb            ICPV=NPLSI+IPLS
!pb            CHMOM(IPLS,1:NSBOX_TAL) = CHMOM(IPLS,1:NSBOX_TAL) +
!pb     .           CPPV(ICPV,1:NSBOX_TAL) - CPVODA(ICPV,1:NSBOX_TAL) 
            END DO

            CHEEM(1:NSBOX_TAL) = CHEEM(1:NSBOX_TAL) +
     .            EPEL(1:NSBOX_TAL) - EPEODA (1:NSBOX_TAL) 
            CHEIM(1:NSBOX_TAL) = CHEIM(1:NSBOX_TAL) +
     .            EPPL_COP(1:NSBOX_TAL) - EPLODA (1:NSBOX_TAL)     
          END IF
          
C
C  SAVE RECOMBINATION SOURCES FOR USE BY SHORT CYCLE
C
          IF (.NOT.LSHORT) THEN
            DO IPLS=1,NPLSI
              DO IN=1,NSBOX_TAL
                IF (PPPL_COP(IPLS,IN) .NE. 0.D0) THEN
!pb                  ALLOCATE(CPMUL)
                  CPMUL => EIRENE_NEW_MULARR()
                  CPMUL%IART = IPLS
                  CPMUL%ICM = IN
                  CPMUL%VALUEM = PPPL_COP(IPLS,IN)
                  CPMUL%NXTMUL => PPPL_COPS(ISTRAI)%PMUL
                  PPPL_COPS(ISTRAI)%PMUL => CPMUL
                ENDIF
                IF (CPPV(IPLS,IN) .NE. 0.D0) THEN
!pb                  ALLOCATE(CPMUL)
                  CPMUL => EIRENE_NEW_MULARR()
                  CPMUL%IART = IPLS
                  CPMUL%ICM = IN
                  CPMUL%VALUEM = CPPV(IPLS,IN)
                  CPMUL%NXTMUL => CPPVS(ISTRAI)%PMUL
                  CPPVS(ISTRAI)%PMUL => CPMUL
                ENDIF
              ENDDO
            ENDDO
            DO IN=1,NSBOX_TAL
              IF (EPPL_COP(IN) .NE. 0.D0) THEN
!pb                ALLOCATE(CPSIM)
                CPSIM => EIRENE_NEW_SIMARR()
                CPSIM%ICS = IN
                CPSIM%VALUES = EPPL_COP(IN)
                CPSIM%NXTSIM => EPPL_COPS(ISTRAI)%PSIM
                EPPL_COPS(ISTRAI)%PSIM => CPSIM
              ENDIF
              IF (EPEL(IN) .NE. 0.D0) THEN
!pb                ALLOCATE(CPSIM)
                CPSIM => EIRENE_NEW_SIMARR()
                CPSIM%ICS = IN
                CPSIM%VALUES = EPEL(IN)
                CPSIM%NXTSIM => EPELS(ISTRAI)%PSIM
                EPELS(ISTRAI)%PSIM => CPSIM
              ENDIF
            ENDDO
          END IF

        ENDIF
C
csw 26jan2011 extra B25
csw 08mar2013 moved to here, i.e. after correction from volume recombination
        if(my_pe == 0 .and..not. lshort) then
          call eirene_extraB25_wneufill(istrai)
        endif
csw
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

        IF (ISTRAI <= NTARGI) THEN
          FLX_EIR = 1._DP
        ELSE
          FLX_EIR = FLXEIR(ISTRAI)
        END IF

        DO 7510 IFL=1,NFLA
          CHPS(IFL)=0.
          SNIS(IFL)=0.
          CHMOS(IFL)=0.
          SMOS(IFL)=0.
          DO 7510 IPLS=1,NPLSI
            IF (IFLB(IPLS).NE.IFL) GOTO 7510
            IPLSV=MPLSV(IPLS)
            DO 7520 IX=1,NDXA
              DO 7530 IY=1,NDYA
                CURPOI => HEADS(IY,IX)%P
                DO WHILE (ASSOCIATED(CURPOI))
                  IT=CURPOI%TRIANGLE
                  IN=NCLTAL(IT)
                  SNICL=(PAPL(IPLS,IN)+PMPL(IPLS,IN)+PIPL(IPLS,IN)+
     .                  PPPL_COP(IPLS,IN))*VOLTAL(IN)*FLX_EIR
                  SNI(IX,IY,IFL,ISTRAI)=SNI(IX,IY,IFL,ISTRAI)+SNICL
                  SNIS(IFL)=SNIS(IFL)+ SNICL
                  CHPS(IFL)=CHPS(IFL)+CHPM(IPLS,IN)*VOLTAL(IN)
                  CURPOI=>CURPOI%NEXT
                ENDDO
7530          CONTINUE
7520        CONTINUE

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
                    CURPOI => HEADS(IY,IX)%P
                    DO WHILE (ASSOCIATED(CURPOI))
                      IT=CURPOI%TRIANGLE
                      IN=NCLTAL(IT)
                      SNIRES=(PAPL(IPLS,IN)+PMPL(IPLS,IN)+
     .                        PIPL(IPLS,IN))*VOLTAL(IN)*FLX_EIR
!pb                    RESSNI(ISTRAI,IFL)=RESSNI(ISTRAI,IFL)+
!pb     .                                 ABS(SIGMA_COP(NCPVI+IPLS,IN)*
!pb     .                                 SNIRES/100.D0)
                      RESSNI(ISTRAI,IFL)=RESSNI(ISTRAI,IFL)+
     .                                   ABS(SIGMA(ISTAT_COP,IN)*
     .                                   SNIRES/100.D0)
                      CURPOI=>CURPOI%NEXT
                    END DO
                  END DO
                END DO
              end if
            END IF

!pb            DO 7539 IADD=NPLSI,3*NPLSI,NPLSI
!pb            IF (NCPVI.LT.IADD+NPLSI) GOTO 7539
            DO 7536 IX=1,NDXA
              DO 7533 IY=1,NDYA
                CURPOI => HEADS(IY,IX)%P
                DO WHILE (ASSOCIATED(CURPOI))
                  IT=CURPOI%TRIANGLE
                  INC=NCLTAL(IT)
                  SIGNUM=SIGN(1._DP,BVIN(IPLSV,IT))
!pb                  SMOCL=(COPV(IADD+IPLS,INC)+CPPV(IADD+IPLS,INC))*
                  SMOCL=(MAPL(IPLS,INC)+MMPL(IPLS,INC)+MIPL(IPLS,INC)+
     .                   CPPV(IPLS,INC))*
     .                   VOLTAL(INC)*1.D-5*SIGNUM*FLX_EIR
                  SMO(IX,IY,IFL,ISTRAI)=SMO(IX,IY,IFL,ISTRAI)+SMOCL
                  SMOS(IFL)=SMOS(IFL)+SMOCL
                  CHMOS(IFL)=CHMOS(IFL)+CHMOM(IPLS,INC)*VOLTAL(INC)
                  CURPOI=>CURPOI%NEXT
                ENDDO
7533          CONTINUE
7536        CONTINUE

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
                     CURPOI => HEADS(IY,IX)%P
                     DO WHILE (ASSOCIATED(CURPOI))
                       IT=CURPOI%TRIANGLE
                       IN=NCLTAL(IT)
                       SIGNUM=SIGN(1._DP,BVIN(IPLSV,IT))
                       SMORES=(MAPL(IPLS,IN)+MMPL(IPLS,IN)+
     .                         MIPL(IPLS,IN))*
     .                        VOLTAL(IN)*1.D-5*SIGNUM*FLX_EIR
!pb                    RESSMO(ISTRAI,IFL)=RESSMO(ISTRAI,IFL)+
!pb     .                                 ABS(SIGMA_COP(IPLS,IN)*
!pb     .                                 SMORES/100.D0*1.D5)
                       RESSMO(ISTRAI,IFL)=RESSMO(ISTRAI,IFL)+
     .                                    ABS(SIGMA(ISTAT_COP,IN)*
     .                                    SMORES/100.D0*1.D5)
                       CURPOI=>CURPOI%NEXT
                     END DO
                   END DO
                 END DO
              end if
            END IF
!pb 7539        CONTINUE
7510    CONTINUE
C
        CHEES=0.
        CHEIS=0.
        SEES=0.
        SEIS=0.
        DO 7540 IX=1,NDXA
          DO 7545 IY=1,NDYA
            CURPOI => HEADS(IY,IX)%P
            DO WHILE (ASSOCIATED(CURPOI))
              IT=CURPOI%TRIANGLE
              IN=NCLTAL(IT)
              SEE(IX,IY,ISTRAI)=SEE(IX,IY,ISTRAI)+
     .           (EAEL(IN)+EMEL(IN)+EIEL(IN)+EPEL(IN))*VOLTAL(IN)*ELCHA
              CHEES=CHEES+CHEEM(IN)*VOLTAL(IN)
              SEES=SEES+(EAEL(IN)+EMEL(IN)+EIEL(IN)+EPEL(IN))*VOLTAL(IN)
C
              SEI(IX,IY,ISTRAI)=SEI(IX,IY,ISTRAI)+
     .           (EAPL(IN)+EMPL(IN)+EIPL(IN)+EPPL_COP(IN))*
     .           VOLTAL(IN)*ELCHA
              CHEIS=CHEIS+CHEIM(IN)*VOLTAL(IN)
              SEIS=SEIS+(EAPL(IN)+EMPL(IN)+EIPL(IN)+EPPL_COP(IN))*
     .                  VOLTAL(IN)
              CURPOI=>CURPOI%NEXT
            ENDDO
7545      CONTINUE
7540    CONTINUE

        IF (.NOT.LSHORT) THEN

!pb replace sigma_cop 
          istat_cop = 0
          do i = 1, nsigvi
            if ((iih(i) == ntalm).and.(igh(i) == 3*NPLSI+1)) then
              istat_cop = i
              exit
            end if
          end do
           
          if (istat_cop > 0) then 
            DO IX=1,NDXA
              DO IY=1,NDYA
                CURPOI => HEADS(IY,IX)%P
                DO WHILE (ASSOCIATED(CURPOI))
                  IT=CURPOI%TRIANGLE
                  IN=NCLTAL(IT)
                  SEERES=(EAEL(IN)+EMEL(IN)+EIEL(IN))*VOLTAL(IN)*FLX_EIR
!pb              RESSEE(ISTRAI)=RESSEE(ISTRAI)+
!pb     .                       ABS(SIGMA_COP(2*NPLSI+1,IN)*
!pb     .                       SEERES/100.D0)
                  RESSEE(ISTRAI)=RESSEE(ISTRAI)+
     .                           ABS(SIGMA(ISTAT_COP,IN)*
     .                           SEERES/100.D0)
                  CURPOI=>CURPOI%NEXT
                END DO
              END DO
            END DO
          end if

!pb replace sigma_cop 
          istat_cop = 0
          do i = 1, nsigvi
            if ((iih(i) == ntalm).and.(igh(i) == 3*NPLSI+2)) then
              istat_cop = i
              exit
            end if
          end do
           
          if (istat_cop > 0) then 
            DO IX=1,NDXA
              DO IY=1,NDYA
                CURPOI => HEADS(IY,IX)%P
                DO WHILE (ASSOCIATED(CURPOI))
                  IT=CURPOI%TRIANGLE
                  IN=NCLTAL(IT)
                  SEIRES=(EAPL(IN)+EMPL(IN)+EIPL(IN))*VOLTAL(IN)*FLX_EIR
!     pb              RESSEI(ISTRAI)=RESSEI(ISTRAI)+
!     pb     .                       ABS(SIGMA_COP(2*NPLSI+2,IN)*
!     pb     .                       SEIRES/100.D0)
                  RESSEI(ISTRAI)=RESSEI(ISTRAI)+
     .                           ABS(SIGMA(ISTAT_COP,IN)*
     .                           SEIRES/100.D0)
                  CURPOI=>CURPOI%NEXT
                END DO
              END DO
            END DO
          end if
        END IF
C
C   NEXT:
C   IF LSHORT: CRITERION TO STOP SHORT CYCLE,
C   IF NOT LSHORT: RESCALE SURFACE SOURCE STRATA
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
7550      CONTINUE
          SEES0(ISTRAI)=SEES*FLXI
          SEIS0(ISTRAI)=SEIS*FLXI
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
          LTEST=.TRUE.
          IF (LSTOP)
     .      WRITE (iunout,*) 'STOP SHORT CYCLE: ALL B2 TIMESTEPS DONE '
          DO 7558 IFL=1,NFLA
            TEST=CHPS(IFL)/(SNIS(IFL)+1.D-60)*100.
            IF (ABS(TEST).GT.CHGP) THEN
              LSTP3=.TRUE.
              LTEST=.FALSE.
              WRITE (iunout,*) 'STOP SHORT CYCLE: PART. SOURCES: ',
     .                     SNIS(IFL),CHPS(IFL),TEST
              WRITE (iunout,*) 'STRATUM ISTRAI, SPECIES IFL ',ISTRAI,IFL
            ENDIF
            TEST=CHMOS(IFL)/(SMOS(IFL)+1.D-60)*100.
            IF (ABS(TEST).GT.CHGMOM) THEN
              LSTP3=.TRUE.
              LTEST=.FALSE.
              WRITE (iunout,*) 'STOP SHORT CYCLE: MOMENTUM SOURCE: ',
     .                     SMOS(IFL),CHMOS(IFL),TEST
              WRITE (iunout,*) 'STRATUM ISTRAI, SPECIES IFL ',ISTRAI,IFL
            ENDIF
7558      CONTINUE
          TEST=CHEES/(SEES+1.D-60)*100.
          IF (ABS(TEST).GT.CHGEE) THEN
            LSTP3=.TRUE.
            LTEST=.FALSE.
            WRITE (iunout,*) 'STOP SHORT CYCLE: EL EN. SOURCE: ',SEES,
     .                   CHEES,TEST
            WRITE (iunout,*) 'STRATUM ISTRAI ',ISTRAI
          ENDIF
          TEST=CHEIS/(SEIS+1.D-60)*100.
          IF (ABS(TEST).GT.CHGEI) THEN
            LSTP3=.TRUE.
            LTEST=.FALSE.
            WRITE (iunout,*) 'STOP SHORT CYCLE: ION EN. SOURCE: ',SEIS,
     .                   CHEIS,TEST
            WRITE (iunout,*) 'STRATUM ISTRAI ',ISTRAI
          ENDIF
          if(lshort) LSTP=LSTP3
        ENDIF          
csw ignore short cycle attempt
        LSTP3=.true.
        LSTOP=LSTP3
        LTEST=.false.
        if(lshort) then
          LSTP=LSTP3
        end if
csw
        NLSRON(ISTRAI) = .NOT.LTEST
C
        IF (.NOT.LSHORT) THEN
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
csw 26jan2011 extra B25, correct particle sources for recycling strata if necessary
        if(my_pe == 0) then
!pb          call eirene_extrab25_eirpbls(istrai)
          srcstrn(istrai)=flux(istrai)
        endif
csw
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
csw 14jul2011
        if(allocated(save_estimv)) then
          estimv = save_estimv
          deallocate(save_estimv)
        endif
        if(allocated(save_estims)) then
          estims = save_estims
          deallocate(save_estims)
        endif
!pb 08.01.2014        if(allocated(save_sigma_cop)) then
!pb 08.01.2014          sigma_cop = save_sigma_cop
!pb 08.01.2014          deallocate(save_sigma_cop)
!pb 08.01.2014         endif
csw
10000 CONTINUE
C
      RETURN
C
      ENTRY EIRENE_IF4COP
C
      NREC11=NOUTAU
      OPEN (UNIT=11,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=8*NREC11)
      IRC=3
      WRITE (11,REC=IRC) RCCPL
      IF (TRCINT.OR.TRCFLE)   WRITE (iunout,*) 'WRITE 11  IRC= ',IRC
      IRC=3
      ALLOCATE (IHELP(NOUTAU))
      JC=0
      DO K=1,NPTRGT
        DO J=1,10*NSTEP
          JC=JC+1
          IHELP(JC)=ICCPL1(J,K)
          IF (JC == NOUTAU) THEN
            IRC=IRC+1
            WRITE (11,REC=IRC) IHELP
            IF (TRCINT.OR.TRCFLE) WRITE (iunout,*) 'WRITE 11  IRC= ',IRC
            JC=0
          END IF
        END DO
      END DO
      IF (JC > 0) THEN
        IHELP(JC+1:NOUTAU) = 0
        IRC=IRC+1
        WRITE (11,REC=IRC) IHELP
        IF (TRCINT.OR.TRCFLE)   WRITE (iunout,*) 'WRITE 11  IRC= ',IRC
      END IF
      DEALLOCATE (IHELP)
      IRC=IRC+1
      WRITE (11,REC=IRC) ICCPL2
      IRC=IRC+1
      WRITE (11,REC=IRC) LCCPL
      IF (TRCINT.OR.TRCFLE)   WRITE (iunout,*) 'WRITE 11  IRC= ',IRC
C
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
        IF (LCUT(IX)) GOTO 10113
C
C SURFACE NORMAL IS INWARD. HENCE: TAKE ALL FLUXES F...YB POSITIVE
C SIGN OF ADDITIONAL COMPONENT DUE TO INCLINED GRID AS SIGN OF F...YB
        SFEISY=SFEISY+FEIYB(IX,0)
        SFEESY=SFEESY+FEEYB(IX,0)
        DO 10111 IF=1,NFLA
!pb          SI=SIGN(1._DP,FNIYB(IX,0,IF))
!pb          SFNISY(IF)=SFNISY(IF)+FNIYB(IX,0,IF)+SI*ABS(FNIY_XB(IX,0,IF))
          SFNISY(IF)=SFNISY(IF)+FNIYB(IX,0,IF)
10111   CONTINUE
        GOTO 10113
10110   CONTINUE
C  DO NOT RECYCLE TARGET FLUXES WITH FALSE ORIENTATION
        IF (ITARG.GT.0) THEN
C  ITARG, IPRT KNOWN FROM ABOVE
          FLX=0.
          TIFLX=0.
          DO 10112 IF=NSPZI(ITARG,IPRT),NSPZE(ITARG,IPRT)
!pb            SI=SIGN(1._DP,FNIYB(IX,0,IF))
!pb            DUMVAL = FNIYB(IX,0,IF)+SI*ABS(FNIY_XB(IX,0,IF))
            DUMVAL = FNIYB(IX,0,IF)
            FLX=FLX+DUMVAL
            IF (NINCT(ITARG,IPRT)*DUMVAL.LT.0) THEN
              WRITE (iunout,*) 
     .          'RECYCLING TARGET, BUT WRONG FLOW DIRECTION: '
              WRITE (iunout,*) 
     .          'SOUTH,IX,ITARG,IPRT,IF ',IX,ITARG,IPRT,IF
              SFNISY(IF)=SFNISY(IF)+DUMVAL
              TIFLX=TIFLX+FEIYB(IX,0)*DUMVAL
            ENDIF
10112     CONTINUE
          TIFLX=TIFLX/(FLX+EPS60)
          SFEISY=SFEISY+TIFLX
        ENDIF
10113 CONTINUE

      SFNISY=SFNISY*ELCHA
C
      WRITE (37,*) 'NON RECYCLING FLUXES FROM SOUTH EDGE '
      WRITE (37,8888) SFNISY,SFEISY,SFEESY
8888  FORMAT (3E14.6)
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
        IF (LCUT(IX)) GOTO 10118
C
C SURFACE NORMAL IS OUTWARD. HENCE: TAKE ALL FLUXES F...YB NEGATIVE
C SIGN OF ADDITIONAL COMPONENT DUE TO INCLINED GRID AS SIGN OF F...YB

        SFEINY=SFEINY-FEIYB(IX,NDYA)
        SFEENY=SFEENY-FEEYB(IX,NDYA)
        DO 10116 IF=1,NFLA
          SI=SIGN(1._DP,FNIYB(IX,NDYA,IF))
!pb          SFNINY(IF)=SFNINY(IF)-FNIYB(IX,NDYA,IF)-
!pb     .               SI*ABS(FNIY_XB(IX,NDYA,IF))
          SFNINY(IF)=SFNINY(IF)-FNIYB(IX,NDYA,IF)
10116   CONTINUE
        GOTO 10118
10115   CONTINUE
        IF (ITARG.GT.0) THEN
C  ITARG, IPRT KNOWN FROM ABOVE
          FLX=0.
          TIFLX=0.
          DO 10117 IF=NSPZI(ITARG,IPRT),NSPZE(ITARG,IPRT)
!pb            SI=SIGN(1._DP,FNIYB(IX,NDYA,IF))
!pb            DUMVAL=FNIYB(IX,NDYA,IF)+SI*ABS(FNIY_XB(IX,NDYA,IF))
            DUMVAL=FNIYB(IX,NDYA,IF)
            FLX=FLX+DUMVAL
            IF (NINCT(ITARG,IPRT)*DUMVAL.LT.0) THEN
              WRITE (iunout,*) 
     .          'RECYCLING TARGET, BUT WRONG FLOW DIRECTION: '
              WRITE (iunout,*) 
     .          'NORTH,IX,ITARG,IPRT,IF ',IX,ITARG,IPRT,IF
              SFNINY(IF)=SFNINY(IF)-DUMVAL
              TIFLX=TIFLX+FEIYB(IX,NDYA)*(-DUMVAL)
            ENDIF
10117     CONTINUE
          TIFLX=TIFLX/(FLX+EPS60)
          SFEINY=SFEINY-TIFLX
        ENDIF
10118 CONTINUE
C
      SFNINY=SFNINY*ELCHA
C
      WRITE (37,*) 'NON RECYCLING FLUXES TO NORTH EDGE '
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
     .             IY.LT.NTEN(ITARG,IPRT)) GOTO 10120
            ENDIF
          ENDDO
        ENDDO
        ITARG=0
C
        SFEIWX=SFEIWX+FEIXB(0,IY)
        SFEEWX=SFEEWX+FEEXB(0,IY)
        DO 10121 IF=1,NFLA
!pb          SI=SIGN(1._DP,FNIXB(0,IY,IF))
!pb          SFNIWX(IF)=SFNIWX(IF)+FNIXB(0,IY,IF)+SI*ABS(FNIX_YB(0,IY,IF))
          SFNIWX(IF)=SFNIWX(IF)+FNIXB(0,IY,IF)
10121   CONTINUE
10120   CONTINUE
        IF (ITARG.GT.0) THEN
C  ITARG, IPRT KNOWN FROM ABOVE
          FLX=0.
          TIFLX=0.
          DO 10122 IF=NSPZI(ITARG,IPRT),NSPZE(ITARG,IPRT)
!pb          SI=SIGN(1._DP,FNIXB(0,IY,IF))
!pb          DUMVAL=FNIXB(0,IY,IF)+SI*ABS(FNIX_YB(0,IY,IF))
          DUMVAL=FNIXB(0,IY,IF)
            FLX=FLX+DUMVAL
            IF (NINCT(ITARG,IPRT)*DUMVAL.LT.0) THEN
              WRITE (iunout,*) 
     .          'RECYCLING TARGET, BUT WRONG FLOW DIRECTION: '
              WRITE (iunout,*) 'WEST,IY,ITARG,IPRT,IF ',IY,ITARG,IPRT,IF
              SFNIWX(IF)=SFNIWX(IF)+DUMVAL
              TIFLX=TIFLX+FEIXB(0,IY)*DUMVAL
            ENDIF
10122     CONTINUE
          TIFLX=TIFLX/(FLX+EPS60)
          SFEIWX=SFEIWX+TIFLX
        ENDIF
10123 CONTINUE
C
      SFNIWX=SFNIWX*ELCHA
C
      WRITE (37,*) 'NON RECYCLING FLUXES FROM WEST EDGE '
      WRITE (37,8888) SFNIWX,SFEIWX,SFEEWX
C
C
C  FOURTH: EAST EDGE: IX=NDXA
C
C NON RECYCLING FLUXES AT EAST EDGE: SFEIEX,SFEEEX,SFNIEX
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
        SFEIEX=SFEIEX-FEIXB(NDXA,IY)
        SFEEEX=SFEEEX-FEEXB(NDXA,IY)
        DO 10126 IF=1,NFLA
!pb          SI=SIGN(1._DP,FNIXB(NDXA,IY,IF))
!pb          SFNIEX(IF)=SFNIEX(IF)-FNIXB(NDXA,IY,IF)-
!pb     .               SI*ABS(FNIX_YB(NDXA,IY,IF))
          SFNIEX(IF)=SFNIEX(IF)-FNIXB(NDXA,IY,IF)
10126   CONTINUE
10125   CONTINUE
        IF (ITARG.GT.0) THEN
C  ITARG, IPRT KNOWN FROM ABOVE
          FLX=0.
          TIFLX=0.
          DO 10127 IF=NSPZI(ITARG,IPRT),NSPZE(ITARG,IPRT)
!pb            SI=SIGN(1._DP,FNIXB(NDXA,IY,IF))
!pb            DUMVAL=FNIXB(NDXA,IY,IF)+SI*ABS(FNIX_YB(NDXA,IY,IF))
            DUMVAL=FNIXB(NDXA,IY,IF)
            FLX=FLX+DUMVAL
            IF (NINCT(ITARG,IPRT)*DUMVAL.LT.0) THEN
              WRITE (iunout,*) 
     .          'RECYCLING TARGET, BUT WRONG FLOW DIRECTION: '
              WRITE (iunout,*) 'EAST,IY,ITARG,IPRT,IF ',IY,ITARG,IPRT,IF
              SFNIEX(IF)=SFNIEX(IF)-DUMVAL
              TIFLX=TIFLX+FEIXB(NDXA,IY)*(-DUMVAL)
            ENDIF
10127     CONTINUE
          TIFLX=TIFLX/(FLX+EPS60)
          SFEIEX=SFEIEX-TIFLX
        ENDIF
10128 CONTINUE
C
      SFNIEX=SFNIEX*ELCHA
C
      WRITE (37,*) 'NON RECYCLING FLUXES TO EAST EDGE '
      WRITE (37,8888) SFNIEX,SFEIEX,SFEEEX
C
C  NEXT: FLUXES TO THOSE SURFACES, AT WHICH RECYCLING BOUNDARY 
C        CONDITIONS ARE SPECIFIED
C
10130 CONTINUE
C
      SFEIT(0)=0.
      SFEET(0)=0.
      SFNIT=0.
      SHEAE(0)=0.
      SHEAI(0)=0.
      DO 10139 I=1,NTARGI
        SFEIT(I)=0.
        SFEET(I)=0.
        SHEAE(I)=0.
        SHEAI(I)=0.
        DO IPRT=1,NTGPRT(I)
          IF (NIXY(I,IPRT).EQ.1) THEN
C  BALANCE CONTRIB. X-GRID REC. SOURCE
            DO 10132 IY=NTIN(I,IPRT),NTEN(I,IPRT)-1
              DO 10131 IF=NSPZI(I,IPRT),NSPZE(I,IPRT)
                IF (NINCT(I,IPRT)*FNIXB(NDT(I,IPRT),IY,IF).GT.0) THEN
!pb                SFNIT(I,IF)=SFNIT(I,IF)-
!pb     .                   NINCT(I,IPRT)*FNIXB(NDT(I,IPRT),IY,IF)-
!pb     .                   ABS(FNIX_YB(NDT(I,IPRT),IY,IF))
                SFNIT(I,IF)=SFNIT(I,IF)-
     .                   NINCT(I,IPRT)*FNIXB(NDT(I,IPRT),IY,IF)
!pb                SHEAE(I)=SHEAE(I)+TEB(NDT(I,IPRT),IY)*
!pb     .           (NINCT(I,IPRT)*FNIXB(NDT(I,IPRT),IY,IF)+
!pb     .            ABS(FNIX_YB(NDT(I,IPRT),IY,IF)))*
!pb     .           (-DELTA_SHEATHXB(NDT(I,IPRT),IY))
                SHEAE(I)=SHEAE(I)+TEB(NDT(I,IPRT),IY)*
     .            NINCT(I,IPRT)*FNIXB(NDT(I,IPRT),IY,IF)*
     .           (-DELTA_SHEATHXB(NDT(I,IPRT),IY))
!pb                SHEAI(I)=SHEAI(I)+TEB(NDT(I,IPRT),IY)*
!pb     .           (NINCT(I,IPRT)*FNIXB(NDT(I,IPRT),IY,IF)+
!pb     .            ABS(FNIX_YB(NDT(I,IPRT),IY,IF)))*
!pb     .           DELTA_SHEATHXB(NDT(I,IPRT),IY)
                SHEAI(I)=SHEAI(I)+TEB(NDT(I,IPRT),IY)*
     .            NINCT(I,IPRT)*FNIXB(NDT(I,IPRT),IY,IF)*
     .            DELTA_SHEATHXB(NDT(I,IPRT),IY)
                ELSE
                  WRITE (iunout,*)
     .              'WRONG ORIENTATION OF W/E-TARGET RECYCLING FLUX '
                  WRITE (iunout,*) 'ITARG, IPRT, IPLS, NDT, IY ',
     .                         I    , IPRT, IF,   NDT(I,IPRT), IY
                  WRITE (iunout,*) 'FNIX(NDT,IY) ',
     .                              FNIXB(NDT(I,IPRT),IY,IF)
                ENDIF
10131         CONTINUE
              SFEIT(I)=SFEIT(I)-NINCT(I,IPRT)*FEIXB(NDT(I,IPRT),IY)
              SFEET(I)=SFEET(I)-NINCT(I,IPRT)*FEEXB(NDT(I,IPRT),IY)
10132       CONTINUE
C  BALANCE CONTRIB. FROM Y-GRID REC. SOURCE
          ELSEIF (NIXY(I,IPRT).EQ.2) THEN
            DO 10135 IX=NTIN(I,IPRT),NTEN(I,IPRT)-1
             IF (LCUT(IX)) GOTO 10135
             DO 10136 IF=NSPZI(I,IPRT),NSPZE(I,IPRT)
                IF (NINCT(I,IPRT)*FNIYB(IX,NDT(I,IPRT),IF).GT.0.) THEN
!pb                SFNIT(I,IF)=SFNIT(I,IF)-
!pb     .                   NINCT(I,IPRT)*FNIYB(IX,NDT(I,IPRT),IF)-
!pb     .                   ABS(FNIY_XB(IX,NDT(I,IPRT),IF))
                SFNIT(I,IF)=SFNIT(I,IF)-
     .                   NINCT(I,IPRT)*FNIYB(IX,NDT(I,IPRT),IF)
!pb                SHEAE(I)=SHEAE(I)+TEB(IX,NDT(I,IPRT))*
!pb     .           (NINCT(I,IPRT)*FNIYB(IX,NDT(I,IPRT),IF)+
!pb     .            ABS(FNIY_XB(IX,NDT(I,IPRT),IF)))*
!pb     .           (-DELTA_SHEATHYB(IX,NDT(I,IPRT)))
                SHEAE(I)=SHEAE(I)+TEB(IX,NDT(I,IPRT))*
     .            NINCT(I,IPRT)*FNIYB(IX,NDT(I,IPRT),IF)*
     .           (-DELTA_SHEATHYB(IX,NDT(I,IPRT)))
!pb                SHEAI(I)=SHEAI(I)+TEB(IX,NDT(I,IPRT))*
!pb     .           (NINCT(I,IPRT)*FNIYB(IX,NDT(I,IPRT),IF)+
!pb     .            ABS(FNIY_XB(IX,NDT(I,IPRT),IF)))*
!pb     .           DELTA_SHEATHYB(IX,NDT(I,IPRT))
                SHEAI(I)=SHEAI(I)+TEB(IX,NDT(I,IPRT))*
     .            NINCT(I,IPRT)*FNIYB(IX,NDT(I,IPRT),IF)*
     .           DELTA_SHEATHYB(IX,NDT(I,IPRT))
                ELSE
                  WRITE (iunout,*)
     .              'WRONG ORIENTATION OF S/N-TARGET RECYCLING FLUX '
                  WRITE (iunout,*) 'ITARG, IPRT, IPLS, IX, NDT ',
     .                         I    , IPRT, IF,   IX, NDT(I,IPRT)
                  WRITE (iunout,*) 'FNIY(IX,NDT) ',
     .                              FNIYB(IX,NDT(I,IPRT),IF)
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
        WRITE (37,8888) (SFNIT(I,IF),IF=1,NFL),SFEIT(I),SFEET(I)
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
      DO 10150 ISTRA=1,NSTRAI
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
             DO 10141 IF=1,NFLA
               SSN(IF)=SSN(IF)+SNI(IX,IY,IF,ISTRA)
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
!pb      BALANI=SFEISY+SFEINY+SFEIT(0)+SHEAI(0)+SSEI+B2QIE+B2VDP+
!pb     .       SFEIWX+SFEIEX
!pb      BALANE=SFEESY+SFEENY+SFEET(0)+SHEAE(0)+SSEE+B2BREM+B2RAD-B2QIE+
!pb     .       SFEEWX+SFEEEX-B2VDP
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
        WRITE (iunout,*) ' NON RECYCLING FLUXES AT SOUTH EDGE '
        CALL EIRENE_MASR2(' SFEISY,SFEESY  ',SFEISY,SFEESY)
        DO IF=1,NFLA
           WRITE(iunout,*) 'SFNISY(IF =',IF,') ',SFNISY(IF)
        ENDDO
        WRITE (iunout,*) ' NON RECYCLING FLUXES AT NORTH EDGE'
        CALL EIRENE_MASR2(' SFEINY,SFEENY  ',SFEINY,SFEENY)
        DO IF=1,NFLA
           WRITE(iunout,*) 'SFNINY(IF =',IF,') ',SFNINY(IF)
        ENDDO
        WRITE (iunout,*) ' NON RECYCLING FLUXES AT WEST EDGE '
        CALL EIRENE_MASR2(' SFEIWX,SFEEWX  ',SFEIWX,SFEEWX)
        DO IF=1,NFLA
           WRITE(iunout,*) 'SFNIWX(IF =',IF,') ',SFNIWX(IF)
        ENDDO
        WRITE (iunout,*) ' NON RECYCLING FLUXES AT EAST EDGE '
        CALL EIRENE_MASR2(' SFEIEX,SFEEEX  ',SFEIEX,SFEEEX)
        DO IF=1,NFLA
           WRITE(iunout,*) 'SFNIEX(IF =',IF,') ',SFNIEX(IF)
        ENDDO
        CALL EIRENE_MASRR1 (' TARGETS,EI',SFEIT(1),NTARGI,5)
        CALL EIRENE_MASRR1 (' TARGETS,EE',SFEET(1),NTARGI,5)
        DO ITARG=1,NTARGI
          DO IF=1,NFLA
             WRITE(iunout,*) 'TARGETS, NI(IF =',IF,') ',
     .             SFNIT(ITARG,IF),ITARG
          ENDDO
        ENDDO
        CALL EIRENE_MASR2(' TOTALS, EI,EE  ',SFEIT(0),SFEET(0))
        DO IF=1,NFLA
           WRITE(iunout,*) 'TOTALS, NI(IF =',IF,') ',SFNIT(0,IF)
        ENDDO
        WRITE (iunout,*) ' NEUTRAL PLASMA INTERACTION: '
        CALL EIRENE_MASR2(' SSEI,SSEE      ',SSEI,SSEE)
        DO IF=1,NFLA
           WRITE(iunout,*) 'SSNI(IF =',IF,') ',SSNI(IF)
        ENDDO
        WRITE (iunout,*) 
     .    ' VOLUMETRIC ENERGY SINKS FOR ELECTRONS, FROM B2 '
        CALL EIRENE_MASR4(' B2BREM,B2RAD,-B2QIE,-B2VDP     ',
     .               B2BREM,B2RAD,-B2QIE,-B2VDP)
        WRITE (iunout,*) 
     .    ' TARGET SHEATH CONTRIBUTIONS,ELECTRONS AND IONS '
        CALL EIRENE_MASRR1 (' TARGETS,EI',SHEAI(1),NTARGI,5)
        CALL EIRENE_MASRR1 (' TARGETS,EE',SHEAE(1),NTARGI,5)
        CALL EIRENE_MASR2(' TOTALS,EI,EE    ',SHEAI(0),SHEAE(0))
        CALL EIRENE_LEER(1)
        CALL EIRENE_MASR2(' BALANI,BALANE  ',BALANI,BALANE)
        DO IF=1,NFLA
           WRITE(iunout,*) 'BALANN(IF =',IF,') ',BALANN(IF)
        ENDDO
        CALL EIRENE_MASR2('REL.ERR.(%)RI,RE',RI,RE)
        DO IF=1,NFLA
           WRITE(iunout,*) 'RN(IF =',IF,') ',RN(IF)
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
        DO IF=MINSPEZ,MAXSPEZ
          BALAN=BALAN+SFNISY(IF)+SFNINY(IF)+SFNIWX(IF)
     .               +SFNIEX(IF)+SFNIT(0,IF)+SSNI(IF)
          TOT=TOT+ABS(SFNISY(IF)+SFNINY(IF))+ABS(SFNIT(0,IF))+
     .            ABS(SSNI(IF))
        ENDDO
        RNT=BALAN/(TOT+EPS60)*100.
        CALL EIRENE_MASJ2('SUMMED OVER     ',MINSPEZ,MAXSPEZ)
        CALL EIRENE_MASR3('BALAN,TOT,RNT           ',BALAN,TOT,RNT)
        WRITE (iunout,*) ' NOISE FROM SOURCE TERMS '
        RESSNI(0,1:NFLA) = RESSNI(0,1:NFLA)/ELCHA
        CALL EIRENE_MASR4(' RESSEE,RESSEI,RESSNI,RESSMO    ',
     .        RESSEE(0),RESSEI(0),SUM(RESSNI(0,1:NFLA)),
     .        SUM(RESSMO(0,1:NFLA)))
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) ' RESSNI-CONTRIBUTIONS BY DIFFERENT SPECIES '
!pb  copy RESSNI to OUTHELP to avoid warnings from Intel compiler
!PB        CALL EIRENE_MASRR1 (' RESSNI    ',RESSNI(0,1:NFLA),NFLA,5)
        OUTHELP(1:NFLA) = RESSNI(0,1:NFLA)
        CALL EIRENE_MASRR1 (' RESSNI    ',OUTHELP,NFLA,5)
        WRITE (iunout,*) ' RESSMO-CONTRIBUTIONS BY DIFFERENT SPECIES '
!pb  copy RESSMO to OUTHELP to avoid warnings from Intel compiler
!PB        CALL EIRENE_MASRR1 (' RESSMO    ',RESSMO(0,1:NFLA),NFLA,5)
        OUTHELP(1:NFLA) = RESSMO(0,1:NFLA)
        CALL EIRENE_MASRR1 (' RESSMO    ',OUTHELP,NFLA,5)
      ENDIF
C
c sputtering
!      WRITE(iunout,*) 'Sputtering: Total'
!      CALL EIRENE_MASYR1('ATOMS    ',SPTATI,LOGATM,0,0,NATM,0,NSTRA,
!     .   TEXTS(NSPH+1))
!      CALL EIRENE_MASYR1('MOLECULES',SPTMLI,LOGMOL,0,0,NMOL,0,NSTRA,
!     .   TEXTS(NSPA+1))
!      CALL EIRENE_MASYR1('TEST IONS',SPTIOI,LOGION,0,0,NION,0,NSTRA,
!     .   TEXTS(NSPAM+1))
!      CALL EIRENE_MASYR1('BULK IONS',SPTPLI,LOGPLS,0,0,NPLS,0,NSTRA,
!     .   TEXTS(NSPAMI+1))
!      CALL EIRENE_MASR1('TOT. FLX',
!     .   SPTATI(0,0)+SPTMLI(0,0)+SPTIOI(0,0)+SPTPLI(0,0))

      spat = 0._dp
      spml = 0._dp
      spio = 0._dp
      sppl = 0._dp
      spat(0:natmi,0:nstrai) = sptaati(0:natmi,0:nstrai) +
     .                         sptmati(0:natmi,0:nstrai) +
     .                         sptiati(0:natmi,0:nstrai) +
     .                         sptphati(0:natmi,0:nstrai) +
     .                         sptpati(0:natmi,0:nstrai) 
      spml(0:nmoli,0:nstrai) = sptamli(0:nmoli,0:nstrai) +
     .                         sptmmli(0:nmoli,0:nstrai) +
     .                         sptimli(0:nmoli,0:nstrai) +
     .                         sptphmli(0:nmoli,0:nstrai) +
     .                         sptpmli(0:nmoli,0:nstrai) 
      spio(0:nioni,0:nstrai) = sptaioi(0:nioni,0:nstrai) +
     .                         sptmioi(0:nioni,0:nstrai) +
     .                         sptiioi(0:nioni,0:nstrai) +
     .                         sptphioi(0:nioni,0:nstrai) +
     .                         sptpioi(0:nioni,0:nstrai) 
      sppl(0:nplsi,0:nstrai) = sptapli(0:nplsi,0:nstrai) +
     .                         sptmpli(0:nplsi,0:nstrai) +
     .                         sptipli(0:nplsi,0:nstrai) +
     .                         sptphpli(0:nplsi,0:nstrai) +
     .                         sptppli(0:nplsi,0:nstrai) 
      WRITE(iunout,*) 'Sputtering: Total'
      CALL EIRENE_MASYR1('ATOMS    ',SPAT,LOGATM,0,0,NATM,0,NSTRA,
     .   TEXTS(NSPH+1))
      CALL EIRENE_MASYR1('MOLECULES',SPML,LOGMOL,0,0,NMOL,0,NSTRA,
     .   TEXTS(NSPA+1))
        CALL EIRENE_MASYR1('TEST IONS',SPIO,LOGION,0,0,NION,0,NSTRA,
     .   TEXTS(NSPAM+1))
        CALL EIRENE_MASYR1('BULK IONS',SPPL,LOGPLS,0,0,NPLS,0,NSTRA,
     .   TEXTS(NSPAMI+1))
      CALL EIRENE_MASR1('TOT. FLX',
     .   SPAT(0,0)+SPML(0,0)+SPIO(0,0)+SPPL(0,0))
C
      CALL EIRENE_LEER (1)
C
11000 CONTINUE
C
!pb      IF (.not.LSHORT) call EIRENE_wneutrals        ! DPC-ADD
csw 07feb2011 extra B2.5
        if(my_pe == 0) then
          call eirene_extraB25_wneusave
        endif
csw
      RETURN
C
csw mpi 07apr2010
      ENTRY EIRENE_BROADCAST_INFCOP
        if(my_pe > 0) then
          if(.not.allocated(lcut)) allocate(lcut(0:NDXP))

          if(.not.allocated(pux)) then
            ALLOCATE (PUX(NRAD))
            ALLOCATE (PUY(NRAD))
            ALLOCATE (PUXE(NRAD))
            ALLOCATE (PUYE(NRAD))
            ALLOCATE (PUXN(NRAD))
            ALLOCATE (PUYN(NRAD))
            ALLOCATE (PVX(NRAD))
            ALLOCATE (PVY(NRAD))
            ALLOCATE (PVXE(NRAD))
            ALLOCATE (PVYE(NRAD))
            ALLOCATE (PVXN(NRAD))
            ALLOCATE (PVYN(NRAD))
          endif
        endif
        call mpi_bcast(lcut(0),ndxp+1,MPI_LOGICAL,
     .                 0,MPI_COMM_WORLD,ier)

        call mpi_bcast(pux,nrad,MPI_DOUBLE_PRECISION,
     .                 0,MPI_COMM_WORLD,ier)
        call mpi_bcast(puy,nrad,MPI_DOUBLE_PRECISION,
     .                 0,MPI_COMM_WORLD,ier)
        call mpi_bcast(puxe,nrad,MPI_DOUBLE_PRECISION,
     .                 0,MPI_COMM_WORLD,ier)
        call mpi_bcast(puye,nrad,MPI_DOUBLE_PRECISION,
     .                 0,MPI_COMM_WORLD,ier)
        call mpi_bcast(puxn,nrad,MPI_DOUBLE_PRECISION,
     .                 0,MPI_COMM_WORLD,ier)
        call mpi_bcast(puyn,nrad,MPI_DOUBLE_PRECISION,
     .                 0,MPI_COMM_WORLD,ier)
        call mpi_bcast(pvx,nrad,MPI_DOUBLE_PRECISION,
     .                 0,MPI_COMM_WORLD,ier)
        call mpi_bcast(pvy,nrad,MPI_DOUBLE_PRECISION,
     .                 0,MPI_COMM_WORLD,ier)
        call mpi_bcast(pvxe,nrad,MPI_DOUBLE_PRECISION,
     .                 0,MPI_COMM_WORLD,ier)
        call mpi_bcast(pvye,nrad,MPI_DOUBLE_PRECISION,
     .                 0,MPI_COMM_WORLD,ier)
        call mpi_bcast(pvxn,nrad,MPI_DOUBLE_PRECISION,
     .                 0,MPI_COMM_WORLD,ier)
        call mpi_bcast(pvyn,nrad,MPI_DOUBLE_PRECISION,
     .                 0,MPI_COMM_WORLD,ier)
      return
csw
csw mpi 09jun2010
      entry eirene_broadcast_if3cop

        if(.not.allocated(snis0)) then
          ALLOCATE (SNIS0(NSTRAi,0:NFL))
          ALLOCATE (SMOS0(NSTRAi,0:NFL))
          ALLOCATE (RESSNI(0:NSTRAi,NFL))
          ALLOCATE (RESSMO(0:NSTRAi,NFL))
          ALLOCATE (RESSEE(0:NSTRAi))
          ALLOCATE (RESSEI(0:NSTRAi))
          ALLOCATE (FLXEIR(NSTRAi))
        endif

        do istrx=1,nstrai
          if(procforstra(istrx,my_pe) .and. my_pe /= 0) then
            SNIS0(istrx,0) = sum(snis0(istrx,1:nfl))
            SMOS0(istrx,0) = sum(smos0(istrx,1:nfl))

            DO IFL=1,NFL
              RESSNI(0,IFL) = ressni(0,ifl) + RESSNI(istrx,IFL)
              RESSMO(0,IFL) = ressni(0,ifl) + RESSMO(istrx,IFL)
            END DO
            RESSEE(0) = ressee(0) + RESSEE(istrx)
            RESSEI(0) = ressei(0) + RESSEI(istrx)
          endif
        enddo

        allocate(dumvec(0:nfl))
        do istrx=1,nstrai

          dumvec(0:nfl) = snis0(istrx,0:nfl)
          call mpi_bcast(dumvec, nfl+1,MPI_DOUBLE_PRECISION,
     .                 0, MPI_COMM_WORLD,ier)
          snis0(istrx,0:nfl) = dumvec(0:nfl)

          dumvec(0:nfl) = smos0(istrx,0:nfl)
          call mpi_bcast(dumvec, nfl+1,MPI_DOUBLE_PRECISION,
     .                 0, MPI_COMM_WORLD,ier)
          smos0(istrx,0:nfl) = dumvec(0:nfl)

          call mpi_bcast(sees0(istrx),    1,MPI_DOUBLE_PRECISION,
     .                 0, MPI_COMM_WORLD,ier)

          call mpi_bcast(seis0(istrx),    1,MPI_DOUBLE_PRECISION,
     .                 0, MPI_COMM_WORLD,ier)

          call mpi_bcast(flxeir(istrx),   1,MPI_DOUBLE_PRECISION,
     .                 0, MPI_COMM_WORLD,ier)

          dumvec(1:nfl) = ressni(istrx,1:nfl)
          call mpi_bcast(dumvec, nfl+1,MPI_DOUBLE_PRECISION,
     .                 0, MPI_COMM_WORLD,ier)
          ressni(istrx,1:nfl) = dumvec(1:nfl)

          dumvec(1:nfl) = ressmo(istrx,1:nfl)
          call mpi_bcast(dumvec, nfl+1,MPI_DOUBLE_PRECISION,
     .                 0, MPI_COMM_WORLD,ier)
          ressmo(istrx,1:nfl) = dumvec(1:nfl)

          call mpi_bcast(ressee(istrx),   1,MPI_DOUBLE_PRECISION,
     .                 0, MPI_COMM_WORLD,ier)

          call mpi_bcast(ressei(istrx),   1,MPI_DOUBLE_PRECISION,
     .                 0, MPI_COMM_WORLD,ier)

        enddo
        deallocate(dumvec)
        return

      entry eirene_mpirecv_if3cop(istrr,irank)
        allocate(dumvec(0:nfl))

        call mpi_recv(dumvec, (nfl+1), MPI_DOUBLE_PRECISION,
     .                irank, istrr, MPI_COMM_WORLD, MPI_STATUS_IGNORE,
     .                ier)
        snis0(istrr,0:nfl) = dumvec(0:nfl)
        
        call mpi_recv(dumvec, (nfl+1), MPI_DOUBLE_PRECISION,
     .                irank, istrr, MPI_COMM_WORLD, MPI_STATUS_IGNORE,
     .                ier)
        smos0(istrr,0:nfl) = dumvec(0:nfl)

        call mpi_recv(sees0(istrr),    1, MPI_DOUBLE_PRECISION,
     .                irank, istrr, MPI_COMM_WORLD, MPI_STATUS_IGNORE,
     .                ier)

        call mpi_recv(seis0(istrr),    1, MPI_DOUBLE_PRECISION,
     .                irank, istrr, MPI_COMM_WORLD, MPI_STATUS_IGNORE,
     .                ier)

        call mpi_recv(flxeir(istrr),   1, MPI_DOUBLE_PRECISION,
     .                irank, istrr, MPI_COMM_WORLD, MPI_STATUS_IGNORE,
     .                ier)

        call mpi_recv(dumvec, nfl+1, MPI_DOUBLE_PRECISION,
     .                irank, istrr, MPI_COMM_WORLD, MPI_STATUS_IGNORE,
     .                ier)
        ressni(istrr,1:nfl) = dumvec(1:nfl)

        call mpi_recv(dumvec, nfl+1, MPI_DOUBLE_PRECISION,
     .                irank, istrr, MPI_COMM_WORLD, MPI_STATUS_IGNORE,
     .                ier)
        ressmo(istrr,1:nfl) = dumvec(1:nfl)

        call mpi_recv(ressee(istrr),   1, MPI_DOUBLE_PRECISION,
     .                irank, istrr, MPI_COMM_WORLD, MPI_STATUS_IGNORE,
     .                ier)

        call mpi_recv(ressei(istrr),   1, MPI_DOUBLE_PRECISION,
     .                irank, istrr, MPI_COMM_WORLD, MPI_STATUS_IGNORE,
     .                ier)

        deallocate(dumvec)
        return

      entry eirene_mpisend_if3cop(istrr,irank)

        allocate(dumvec(0:nfl))

        dumvec(0:nfl) = snis0(istrr,0:nfl)
        call mpi_send(dumvec, (nfl+1), MPI_DOUBLE_PRECISION,
     .                irank,istrr, MPI_COMM_WORLD, ier)

        dumvec(0:nfl) = smos0(istrr,0:nfl)
        call mpi_send(dumvec, (nfl+1), MPI_DOUBLE_PRECISION,
     .                irank,istrr, MPI_COMM_WORLD, ier)

        call mpi_send(sees0(istrr),    1, MPI_DOUBLE_PRECISION,
     .                irank,istrr, MPI_COMM_WORLD, ier)

        call mpi_send(seis0(istrr),    1, MPI_DOUBLE_PRECISION,
     .                irank,istrr, MPI_COMM_WORLD, ier)

        call mpi_send(flxeir(istrr),   1, MPI_DOUBLE_PRECISION,
     .                irank,istrr, MPI_COMM_WORLD, ier)

        dumvec(1:nfl) = ressni(istrr,1:nfl)
        call mpi_send(dumvec, nfl+1, MPI_DOUBLE_PRECISION,
     .                irank,istrr, MPI_COMM_WORLD, ier)

        dumvec(1:nfl) = ressmo(istrr,1:nfl)
        call mpi_send(dumvec, nfl+1, MPI_DOUBLE_PRECISION,
     .                irank,istrr, MPI_COMM_WORLD, ier)

        call mpi_send(ressee(istrr),   1, MPI_DOUBLE_PRECISION,
     .                irank,istrr, MPI_COMM_WORLD, ier)

        call mpi_send(ressei(istrr),   1, MPI_DOUBLE_PRECISION,
     .                irank,istrr, MPI_COMM_WORLD, ier)

        deallocate(dumvec)
        return
csw
       contains
C
C FINDS IF THE POINT X;Y BELONGS TO INTERVAL (X1,Y1)..(X2,Y2)
C
       FUNCTION EIRENE_POINT_ON_INTERVAL(X,Y,X1,Y1,X2,Y2)
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



C DEFINE  NORMAL DIRECTION FOR SURFACE AVERAGED TALLIES (SEE FOLNEUT.F)
       SUBROUTINE CORRECTNSS

         IF(IT.GT.NTRIS.OR.NBAR.GT.NTRIS.OR.
     .     NUMSI.GT.3.OR.NBARSI.GT.3)
     .     WRITE(iunout,*) "ERROR IN CORRECTNSS",
     .                  "IT,NTRIS,NBAR,NTRIS,NUMSI,NBARSI",
     .                   IT,NTRIS,NBAR,NTRIS,NUMSI,NBARSI

         INMTINSS(NUMSI,IT)=1
         INMTINSS(NBARSI,NBAR)=1

         IF ((IXTRI(IT) > 0) .AND. (IYTRI(IT) > 0) .AND.
     .       (IXTRI(NBAR) > 0) .AND. (IYTRI(NBAR) > 0)) THEN
!  both triangles inside mesh 
            IF (IXTRI(IT) == IXTRI(NBAR)) THEN
              IF (IYTRI(IT) > IYTRI(NBAR)) THEN
                INMTINSS(NUMSI,IT) = -1
              ELSE
                INMTINSS(NBARSI,NBAR) = -1
              END IF

            ELSE IF (IYTRI(IT) == IYTRI(NBAR)) THEN
              IF (IXTRI(IT) > IXTRI(NBAR)) THEN
                INMTINSS(NUMSI,IT) = -1
              ELSE
                INMTINSS(NBARSI,NBAR) = -1
              END IF
            END IF


!  triangle IT inside mesh, triangle NBAR outside mesh  
          ELSE IF ((IXTRI(IT) > 0) .AND. (IYTRI(IT) > 0)) THEN

!  poloidal surface
            IF (LXSRF) THEN
              IF (IXTRI(IT).EQ.1) THEN        ! 'W SURFACE'
                INMTINSS(NUMSI,IT ) = -1
              ELSE                            ! 'E SURFACE'
                INMTINSS(NBARSI,NBAR) = -1
              END IF

!  radial surface
            ELSE
              IF (IYTRI(IT).EQ.1) THEN        ! 'S SURFACE'
                INMTINSS(NUMSI,IT ) = -1
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
                INMTINSS(NUMSI,IT) = -1
              END IF

!  radial surface
            ELSE
              IF (IYTRI(NBAR).EQ.1) THEN      ! 'S SURFACE'
                INMTINSS(NBARSI,NBAR ) = -1
              ELSE                            ! 'N SURFACE'
                INMTINSS(NUMSI,IT) = -1
              END IF
            
            END IF

!  both triangles outside mesh --> don't know - do nothing 
          END IF

       END  SUBROUTINE CORRECTNSS
CVK END

      END
