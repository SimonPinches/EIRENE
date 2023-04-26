cdr Jan 2016: clean-up, comments.
cdr           This is master version for all other versions of infcop.f

cdr Nov. 17: removed dead option LINDIM: here and in couple_b2_parallel
cdr Dec. 17:
c    trcsou --> trcint:  consistency checks for target recycling step functions
c               trcsou: print step functions from samsrf.f, as finally used in eirene.
c   Jan. 18:  SSNI --> SSNI(ifl), rn, balann, etc... species index in global part. bal.
c             transfered to here from couple_b2.5
c             EPEL --> EPEL_COP  (also in couple_B2.5)
c             CPPV --> MPPL_COP
c             ELTEST, EMAXW,... for a target energy flux as interpreted from B2 output.
cdr Oct. 18:  bug fix: bvin(iplsv) instead bvin(ipls) in 2 places

C
C   EIRENE CODE SEGMENT COUPLE_$, $ MAY CURRENTLY STAND FOR B2,
C                                                           B2.5,
C                                                           DIVIMP,
C                                                           TRIA,
C                                                           TETRA,
C                                                           TRANSP,
C                                                           DUMMY
C
C   THIS VERSION: $=B2,  NOV. 2002
c                 proprietary version of FZJ, for local B2 code versions
C
c  geometry data not any longer via work array into eirene
c                due to module structure
c  eliminate cut cells from balances (lcut(..))
c  removed: ncopib, ncopeb

C   UPDATES:
C   OPTION TO EVALUATE B FIELD VECTORS FROM GRIDADAP FILE FT29
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
C   THIS PARTICULAR VERSION LINKS EIRENE TO THE B2   2D MULTIFLUID EDGE
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
      MODULE EIRMOD_INFCOP

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
      USE EIRMOD_SHEATH
      USE EIRMOD_JSON
      USE EIRMOD_OPENFILE, ONLY: EIRENE_OPENFILE
      
      use json_module       

      IMPLICIT NONE

      PRIVATE

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
      REAL(DP), ALLOCATABLE, SAVE ::
     .            CHPM(:,:), CHEEM(:), CHEIM(:),
     .            CHMOM(:,:)


C pppl_cop, mppl_cop, eppl_cop and epel_cop are the exact
c volumetric source tallies,
c while default tallies pppl, mppl, eppl and epel are
c the corresponding tallies scored from random sampling in eirene
      REAL(DP), ALLOCATABLE, SAVE ::
     .            PPPL_COP(:,:), MPPL_COP(:,:),
     .            EPPL_COP(:,:), EPEL_COP(:)
C
     .           ,CPV_CMP(:,:,:)
      REAL(DP), ALLOCATABLE, SAVE ::
     .            PPLODA(:,:), CPVODA(:,:),
     .            EPLODA(:,:), EPEODA(:)
C

      REAL(DP), SAVE :: SCALM, SCALE, SCALI, CHEIS, SEES, SEIS, TEST,
     .          SFEISY, SFEESY, RECADD, RECTOT,
     .          EEADD, PIADD, SIGNUM, SMOCL, CHEES, EIADD, SNICL,
     .          SSE, BALANI, BALANE, SSEE, SSI, RE, RI, RNT, TOT,
     .          TOTI, TOTE, BALAN,
     .          SSEI, SFEIEX, SFEEEX, SFEENY, SFEIWX, VVBC,
     .          UUBC, UPBC, RBC, UDBC, VL, V, T, BX, BY, BZ, BN,
c    .          DELTE_PARA, DELTI_PARA, DELTE_PERP, DELTI_PERP, TES,TIS,
     .          DELY, ALX, ALE, ALW, ALS, ALN, AL, ETOT,
     .          FLX, ESUM, VR, VTEST, EADD, EMAXW, ESHEATH,
     .          PARWI, PERWI, SUMM, SUMN, SUMEI, SUMEE, FLXI, CHP,
     .          CHI, CHE, CS, THMAX, EESHT, EEMAX,
     .          RP1, DELX, PVYS, PVXS, PUPV, RRBS, PUYS, PUXS,
     .          VPX, VPY, VT, PARW, PERW, PN1, OR, VPZ, GAMMA, CUR, TE,
     .          SFEEWX, SFEINY, PM1, DRR, UU, PITB,
     .          FLX_EIR, SUMN_OLD, SNIRES, SMORES, SEERES, SEIRES,
     .          fltt, e0b2, dmaxiso, dminiso, CFAC, 
     .          eamisum, eplsum, bv

      INTEGER, SAVE :: J, IRC, JC, INC,
     .           IADD, IP, ITARG, IO, IFL, NPES,
     .           IIPLS, IG, IGITT, IEPLS, NPEC, NPBC, NPBS, NTGPRI,
     .           IT, I, IPRT, IAOT, IAIN, IREAD, IPL, INN,
     .           IERROR, LTARG, IN, IX, IY,
     .           NPLP, NDX2, NRED, IO29, NDXY, IFIRST,
     .           ISTRAI, IRRC, K, IR, IIRC, ICPV, IF, I34,
     .           NREC11, NEM, MINSPEZ, MAXSPEZ,
     .           JPLS, ISR, ISTEP, IST_RATE, ISTR,
     .           IXI, IXE, IPLSTI, IPLSV, IPLV, ISP,
     .           IR1, IR2,  imf,
     .           icp, icp2, icp3,
     .           istat_cop, ixm1, iym1, icp4, icp5, 
     .           js, iunin_save, iusrout

      INTEGER, ALLOCATABLE, SAVE :: IZDEN(:)
      REAL(DP), EXTERNAL :: EIRENE_STEP, EIRENE_FTABRC1, EIRENE_FEELRC1,
     .            EIRENE_EMAXW
      INTEGER, EXTERNAL :: EIRENE_IDEZ
C
      LOGICAL, SAVE :: LSHORT, LSTOP, LTEST, LSTP3
      LOGICAL, ALLOCATABLE, SAVE :: LLCUT(:)
      LOGICAL, ALLOCATABLE, SAVE :: LZDENA(:,:), LZDENM(:,:),
     .                              LZDENI(:,:)
C

      REAL(DP), ALLOCATABLE, SAVE ::
     . CHPS(:),    SNIS(:),    CHMOS(:),  SMOS(:),  SCALN(:),
     . SNIS0(:,:), SMOS0(:,:),
c
     . RESSNI(:,:),  RESSMO(:,:),
     . RESSEE(:), RESSEI(:)
     .,FLXEIR(:)

c  ion energy sources for internal (rather than total) ion energy balance
      REAL(DP), ALLOCATABLE, SAVE :: SCPVEII(:)

      REAL(DP), ALLOCATABLE, SAVE ::
     . TORL(:,:), ESHT(:,:), ELTEST(:,:), ORI(:,:)

      real(dp), allocatable :: helpw(:)

      INTEGER, ALLOCATABLE, SAVE :: IHELP(:)
C
      CHARACTER(10) :: CHR
      CHARACTER(6)  :: CITARG
      CHARACTER(72) :: ZEILE
      CHARACTER(20) :: FORM
      CHARACTER(1) :: NSEW

      TYPE(RATE_STORE), POINTER :: RTIS

      DATA LTARG/0/

      PUBLIC :: EIRENE_INFCOP, EIRENE_IF0COP, EIRENE_IF1COP,
     .          EIRENE_IF2COP, EIRENE_IF3COP, EIRENE_IF4COP,
     .          EIRENE_INTER0, EIRENE_INTER3, EIRENE_INFCOP_PRE_MCARLO,
     .          EIRENE_INFCOP_PRE_STRATA, EIRENE_INFCOP_POST_STRATUM,
     .          EIRENE_IF3COP_SUM
      
      
      CONTAINS

      SUBROUTINE EIRENE_INFCOP
C     Previously ENTRY statements were used in this function and LSHORT
C     would be set TRUE when calling the subroutine itself. This may or
C     may not be the intended behaviour
      LOGICAL :: LFIXED,LSHRT
      CALL EIRENE_INFCOP_GENERIC(LFIXED,LSHRT)
      END SUBROUTINE
           
      SUBROUTINE EIRENE_IF0COP(LFIXED,LSHRT)
      LOGICAL, INTENT(IN) :: LFIXED,LSHRT
      LSHORT=.FALSE.
      CALL EIRENE_INFCOP_GENERIC(LFIXED,LSHRT)
      END SUBROUTINE

      SUBROUTINE EIRENE_INTER0
      LOGICAL :: LFIXED,LSHRT
      LSHORT=.TRUE.      
      CALL EIRENE_INFCOP_GENERIC(LFIXED,LSHRT)
      END SUBROUTINE

      SUBROUTINE EIRENE_INTER1
      LSHORT=.TRUE.
      CALL EIRENE_INTER1_GENERIC()
      END SUBROUTINE
      
      SUBROUTINE EIRENE_IF3COP(IENTRY,LSTP,
     .     IFRST,ISTRAA,ISTRAE,NEW_ITER)
      INTEGER, INTENT(IN) :: IENTRY
      LOGICAL, INTENT(INOUT) :: LSTP
      INTEGER, INTENT(IN) :: ISTRAA, ISTRAE, NEW_ITER, IFRST
      WRITE (iunout,*) ' IF3COP IS CALLED, ISTRAA,ISTRAE '
      WRITE (iunout,*) ISTRAA,ISTRAE
      LSHORT=.FALSE.
      LSTP3=.TRUE.
      LSTOP=.TRUE.
      IFIRST=0
      NDXY=(NDXA-1)*NR1ST+NDYA
      CALL EIRENE_INTF_GENERIC(LSTP,IFRST,ISTRAA,ISTRAE,NEW_ITER)
      END SUBROUTINE

      SUBROUTINE EIRENE_INTER3(LSTP,IFRST,ISTRAA,ISTRAE,NEW_ITER)
      LOGICAL, INTENT(INOUT) :: LSTP
      INTEGER, INTENT(IN) :: ISTRAA, ISTRAE, NEW_ITER, IFRST
      LSHORT=.TRUE.
      LSTP3=LSTP
      LSTOP=LSTP
      IFIRST=IFRST
      NDXY=(NDXA-1)*NR1ST+NDYA
      CALL EIRENE_INTF_GENERIC(LSTP,IFRST,ISTRAA,ISTRAE,NEW_ITER)
      END SUBROUTINE     
      
      SUBROUTINE EIRENE_INFCOP_GENERIC(LFIXED,LSHRT)
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


      LOGICAL, INTENT(IN) :: LFIXED,LSHRT
      REAL(DP) :: DUMMY(0:NDXP,0:NDYP)

!HJL These entry points replaced with subroutines above
!      ENTRY EIRENE_IF0COP(LFIXED,LSHRT)
!      ENTRY EIRENE_INTER0

      call eirene_leer(2)
      write (iunout,*) 'Subr. INFCOP called '
      write (iunout,*) 'This is a proprietary FZJ version of an '
      write (iunout,*) 'interfacing code to B2, B2.5 plasma solvers.'
      write (iunout,*) 'NOT ready for 3rd parties'
      call eirene_leer(2)
C
      IERROR=0
      IUSROUT = 0
      IUNIN_SAVE = IUNIN
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

      mshfrm = 0   !  optional flag for geometry file format: linda, carre, sonnet
      NLSHRT13 = .TRUE.  !  only short version of fort13 is used: calls WRPLAM_SHRT, RPLAM_SHRT
C
      IF (.NOT.LSHORT.AND.ITIMV.LE.1) THEN
        IF (LFIXED) THEN
          CALL EIRENE_READ14_FIXED
        ELSE
          js = itree_num(14)
          CALL EIRENE_READ14_JSON(jtrees(js),blks(14)%p)
        ENDIF
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

      IF (LZDEN) THEN
        ALLOCATE (LZDENA(NRTAL,NSTRA))
        ALLOCATE (LZDENM(NRTAL,NSTRA))
        ALLOCATE (LZDENI(NRTAL,NSTRA))
        ALLOCATE (IZDEN(NSTRA))
        LZDENA = .FALSE.
        LZDENM = .FALSE.
        LZDENI = .FALSE.
        IZDEN = 0
      END IF
C
C
C  TRANSFER GEOMETRY
C
      IF (.NOT.(INDGRD(1).EQ.6.OR.INDGRD(2).EQ.6.OR.INDGRD(3).EQ.6))THEN
        IUNIN = IUNIN_SAVE
        IF (IUSROUT /= 0) CLOSE(IUSROUT)
        RETURN
      END IF
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
!  ALPHXB, ALPHYB GIVE THE DIRECTION OF THE B FIELD IN THE
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
     .  ' NO FILE '//FORT//'29 WITH MODIFIED GRID INFO. FOUND '
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

      IUNIN = IUNIN_SAVE
      IF (IUSROUT /= 0) CLOSE(IUSROUT)
C
      RETURN
      END SUBROUTINE
C
C   GEOMETRY DEFINITION PART FINISHED
C
!HJL Replaced ENTRY point      
!     ENTRY EIRENE_IF1COP(IENTRY)
      SUBROUTINE EIRENE_IF1COP(IENTRY)

      INTEGER, INTENT(IN) :: IENTRY
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
!     IF (NLPLAS) GOTO 2100
!HJL  Replaced gotos with subroutine calls. 1st == 2100, 2nd == 99991
!     This does not seem to match the logic in the original ENTRY calls
!     However it works and what seems to be the correct logic does not
      IF (NLPLAS) CALL EIRENE_INTER1()
      CALL EIRENE_IF1COP_EXT()
      
      RETURN
      END SUBROUTINE
C
C  IN CASE OF "SHORT CYCLE" OR TIME DEP. MODE
C  THE PLASMA STATE IS TRANSFERRED VIA COMMON
C  ONLY SCALING TO EIRENE UNITS AND INDEX MAPPING NEEDS TO BE DONE HERE
C
!HJL Replaced entry with subroutine
!     ENTRY EIRENE_INTER1

!      LSHORT=.TRUE.
!      GOTO 2100
C
!99991 CONTINUE
      SUBROUTINE EIRENE_IF1COP_EXT
      REAL(DP) :: DUMMY(0:NDXP,0:NDYP)
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
c     CALL EIRENE_ALLOC_BRAEIR(NDXP,NDYP,NFL)
cdr  for b2.5:    -1 ...NDX  ?
      CALL EIRENE_ALLOC_BRAEIR(NDX,NDY,NFL)
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

      END SUBROUTINE
C
!     2100 CONTINUE
      SUBROUTINE EIRENE_INTER1_GENERIC()
      REAL(DP) :: DUMMY(0:NDXP,0:NDYP)

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
!     IF (NCUTL.EQ.NCUTB_SAVE) GOTO 2101
!HJL GOTO replaced with do while      
      IF (NCUTL.NE.NCUTB_SAVE) THEN
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
C  distinct from B2: these velocities are cell-centered in b2.5
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
      END IF
!HJL GOTO replaced with IF
!     2101 CONTINUE
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
 2105 CONTINUE
C

!pb find density in isolated cells

      dmaxiso = -huge(1._dp)
      dminiso = huge(1._dp)

      do iy = 1, ndya
        do ix = 1, ndxa
          IN=IY+(IX-1)*NR1ST
          if (nstgrd(in) == 0) cycle
          dmaxiso = max(dmaxiso,dnib(ix,iy,1))
          dminiso = min(dminiso,dnib(ix,iy,1))
        end do
      end do
      write (iunout,*) 'dminiso = ',dminiso
      write (iunout,*) 'dmaxiso = ',dmaxiso
      do iy = 1, ndya
        do ix = 1, ndxa
          if (dnib(ix,iy,1) > dminiso) aisob(ix,iy) = 1._dp
        end do
      end do


      BZINTF = 1._DP
      BFINTF = 1._DP
      DO 2110 IY=1,NDYA
        DO 2120 IX =1,NDXA
          IN=IY+(IX-1)*NR1ST
          TEINTF(IN)=TEB(IX,IY)*T*AISOB(IX,IY)
C
C  ONLY ONE ION TEMPERATURE AVAILABLE FROM PLASMA FLUID CODE,
C  SEE LOOP 2150 BELOW
          TIINTF(1,IN)=TIB(IX,IY)*T*AISOB(IX,IY)
C
          BX=PUX(IN)*RRB(IX,IY)+PVX(IN)*0.
          BY=PUY(IN)*RRB(IX,IY)+PVY(IN)*0.
          BZ=SQRT(1.-RRB(IX,IY)**2)
c  normalize B field vector to length 1 (one)
c
          BN=SQRT(BX*BX+BY*BY+BZ*BZ)
          BXINTF(IN)=BX/BN
          BYINTF(IN)=BY/BN
          BZINTF(IN)=BZ/BN
          BFINTF(IN)=BN
          VLINTF(IN)=VOLB(IX,IY)*VL
 2120   CONTINUE
 2110 CONTINUE
C
C  SET SAME ION TEMPERATURE FOR ALL EIRENE BACKGROUND SPECIES
C
      DO IPLS=1,NPLSTI
        DO IY=1,NDYA
          DO IX =1,NDXA
            IN=IY+(IX-1)*NR1ST
            TIINTF(IPLS,IN)=TIINTF(1,IN)
          END DO
        END DO
      END DO
C
CDR  set density from B2 array DNIB, for each fluid
CDR  set plasma flow velocity field from B2 arrays UPB (parallel velocity)
c  without drifts:
c  upb * pitch: poloidal velocity (i.e. cartesian x,y direction).
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
            DO 2210 IY=1,NDYA
              IYM1=IY-1
              DO 2220 IX = 1,NDXA
                IXM1 = IX-1
                IF (LLCUT(IXM1)) THEN
                  IXM1 = NGHPOL(4,IY,IX)
                END IF

C  SET PLASMA DENSITY
                IN=IY+(IX-1)*NR1ST
                DIINTF(IPLS,IN)=DNIB(IX,IY,IFL)*D(IPLS)*AISOB(IX,IY)
C  INTERPOLATE PARALLEL PLASMA VELOCITIES TO CELL CENTERS.
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
                VXINTF(IPLSV,IN)=(PUX(IN)*UUBC+PVX(IN)*VVBC)*V
                VYINTF(IPLSV,IN)=(PUY(IN)*UUBC+PVY(IN)*VVBC)*V
                VZINTF(IPLSV,IN)=(RBC*UPBC-RRB(IX,IY)*UDBC)*V
 2220         CONTINUE
 2210       CONTINUE
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

            IPLSTI= MPLSTI(IPLS)
            IPLSV= MPLSV(IPLS)
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
C  READ OTHER B2 ARRAYS INTO EIRENE, FOR PRINTOUT AND PLOTTING
C
c  density, species index as in B2 code, cell-centered
      DO 2300 IAIN=1,NAINB
        IF (NAINT(IAIN).EQ.1.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO 2321 IY=1,NDYA
            DO IX=1,NDXA
              IN=IY+(IX-1)*NR1ST
              ADINTF(IAIN,IN)=DNIB(IX,IY,NAINS(IAIN))
            END DO
 2321     CONTINUE
c  poloidal (projection) flow velocity, species index as in B2 code, north surface-centered
        ELSEIF (NAINT(IAIN).EQ.2.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO 2322 IY=1,NDYA
            DO IX=1,NDXA
              IN=IY+(IX-1)*NR1ST
              ADINTF(IAIN,IN)=UUB(IX,IY,NAINS(IAIN))
            END DO
 2322     CONTINUE
c  radial drift velocity, species index as in B2 code, east surface-centered
        ELSEIF (NAINT(IAIN).EQ.3.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO 2323 IY=1,NDYA
            DO IX=1,NDXA
              IN=IY+(IX-1)*NR1ST
              ADINTF(IAIN,IN)=VVB(IX,IY,NAINS(IAIN))
            END DO
 2323     CONTINUE
c  plasma pressure, cell-centered, no species index
        ELSEIF (NAINT(IAIN).EQ.6) THEN
          DO 2326 IY=1,NDYA
            DO IX=1,NDXA
              IN=IY+(IX-1)*NR1ST
              ADINTF(IAIN,IN)=PRB(IX,IY)
            END DO
 2326     CONTINUE
c  parallel velocity, species index as in B2 code, north surface-centered
        ELSEIF (NAINT(IAIN).EQ.7.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO 2327 IY=1,NDYA
            DO IX=1,NDXA
              IN=IY+(IX-1)*NR1ST
              ADINTF(IAIN,IN)=UPB(IX,IY,NAINS(IAIN))
            END DO
 2327     CONTINUE
c  pitch angle, no species index
        ELSEIF (NAINT(IAIN).EQ.8) THEN
          DO 2328 IY=1,NDYA
            DO IX=1,NDXA
              IN=IY+(IX-1)*NR1ST
              ADINTF(IAIN,IN)=RRB(IX,IY)
            END DO
 2328     CONTINUE
        ELSEIF (NAINT(IAIN).EQ.9.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO 2329 IY=1,NDYA
            DO IX=1,NDXA
              IN=IY+(IX-1)*NR1ST
              ADINTF(IAIN,IN)=FNIXB(IX,IY,NAINS(IAIN))
            END DO
 2329     CONTINUE
        ELSEIF (NAINT(IAIN).EQ.10.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO 2330 IY=1,NDYA
            DO IX=1,NDXA
              IN=IY+(IX-1)*NR1ST
              ADINTF(IAIN,IN)=FNIYB(IX,IY,NAINS(IAIN))
            END DO
 2330     CONTINUE
        ELSEIF (NAINT(IAIN).EQ.11) THEN
          DO 2331 IY=1,NDYA
            DO IX=1,NDXA
              IN=IY+(IX-1)*NR1ST
              ADINTF(IAIN,IN)=FEIXB(IX,IY)
             END DO
 2331     CONTINUE
        ELSEIF (NAINT(IAIN).EQ.12) THEN
          DO 2332 IY=1,NDYA
            DO IX=1,NDXA
              IN=IY+(IX-1)*NR1ST
              ADINTF(IAIN,IN)=FEIYB(IX,IY)
            END DO
 2332     CONTINUE
        ELSEIF (NAINT(IAIN).EQ.13) THEN
          DO 2333 IY=1,NDYA
            DO IX=1,NDXA
              IN=IY+(IX-1)*NR1ST
              ADINTF(IAIN,IN)=FEEXB(IX,IY)
            END DO
 2333     CONTINUE
        ELSEIF (NAINT(IAIN).EQ.14) THEN
          DO 2334 IY=1,NDYA
            DO IX=1,NDXA
              IN=IY+(IX-1)*NR1ST
              ADINTF(IAIN,IN)=FEEYB(IX,IY)
            END DO
 2334     CONTINUE
C   NAINT=15,16:  USED ONLY IN B2.5 COUPLING: UUDIAG, VVDIAG
c   cell volume as in b2 code, no species index (cell-centered)
        ELSEIF (NAINT(IAIN).EQ.17) THEN
          DO 2335 IY=1,NDYA
            DO IX=1,NDXA
              IN=IY+(IX-1)*NR1ST
              ADINTF(IAIN,IN)=VOLB(IX,IY)
            END DO
 2335     CONTINUE
cdr  magnetic field strength, Tesla
        ELSEIF (NAINT(IAIN).EQ.18) THEN
          DO 2336 IY=1,NDYA
            DO IX=1,NDXA
              IN=IY+(IX-1)*NR1ST
              ADINTF(IAIN,IN)=BFELDB(IX,IY)
            END DO
 2336     CONTINUE


cdr  free: NAINT=20 --29:  reserved for AMDIAG:  scaled atomic/molecular rate coefficients
cdr                        evaluated on computational grid. See Manual.

        ENDIF
 2300 CONTINUE
C
      RETURN

      END SUBROUTINE
C
! 2999 CONTINUE
C
C  PLASMA PROFILES ARE NOW READ IN
C
!      ENTRY EIRENE_IF2COP(ITRG)
      SUBROUTINE EIRENE_IF2COP(ITRG)
      INTEGER, INTENT(IN) :: ITRG
      REAL(DP) :: DI(NPLS), VP(NPLS), ZI(NPLS)
      REAL(DP) :: EFLX(NSTRA),
     R          DUMMY(0:NDXP,0:NDYP),
     R          SFNIT(0:NSTEP,NFL), SFEIT(0:NSTEP),
     R          SFEET(0:NSTEP), SHEAE(0:NSTEP), SHEAI(0:NSTEP)
      INTEGER :: NRWL(NSTRA)

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
! 3000 CONTINUE
C
      IF (TRCINT.AND.LTARG.EQ.0) THEN
        LTARG=1
        WRITE (iunout,*) 'ITARG: TARGET NUMBER '
        WRITE (iunout,*) 'IPRT : SUBSECTION OF TARGET '
        WRITE (iunout,*)
     .     'NPBS : BRAAMS (SURFACE) X-CELL INDEX OF TARGET '
        WRITE (iunout,*) 'NPBC : BRAAMS (ZONE) P-CELL INDEX OF TARGET '
        WRITE (iunout,*) 'NPES : POLOIDAL SURFACE INDEX OF TARGET'
        WRITE (iunout,*) '       IN EIRENE MESH'
        WRITE (iunout,*) 'NPEC : 1ST POLOIDAL CELL INDEX OF EIRENE MESH'
        WRITE (iunout,*) '       SEEN BY MONTE CARLO HISTORIES'
      ENDIF
C
      DO 3005 IPLS=1,NPLSI
        DO IGITT=1,NGITT
          FLSTEP(IPLS,ITARG,IGITT)=0.
        END DO
 3005 CONTINUE
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
!        IF (NIXY(ITARG,IPRT).EQ.2) GOTO 3020
        IF (NIXY(ITARG,IPRT).NE.2) THEN

        IF (NINCT(ITARG,IPRT) > 0) THEN
! EAST TARGET IN B2
          IF (NPBC > 0) THEN
          IF ((NGHPOL(4,NTIN(ITARG,IPRT),NPBC) > 0) .AND.
     .        (ABS(NPBC-NGHPOL(4,NTIN(ITARG,IPRT),NPBC)) > 1)) THEN
            NPBS=NGHPOL(4,NTIN(ITARG,IPRT),NPBC)
            NPEC=NGHPOL(4,NTIN(ITARG,IPRT),NPBC)
            NPES=NPEC+1
          END IF
          END IF
        ELSE IF  (NINCT(ITARG,IPRT) < 0) THEN
! WEST TARGET IN B2
          IF (NPBC > 0) THEN
          IF ((NGHPOL(2,NTIN(ITARG,IPRT),NPBC) > 0) .AND.
     .        (ABS(NPBC-NGHPOL(2,NTIN(ITARG,IPRT),NPBC)) > 1)) THEN
            NPBS=NGHPOL(2,NTIN(ITARG,IPRT),NPBC)
            NPEC=NGHPOL(2,NTIN(ITARG,IPRT),NPBC)
            NPES=NPEC
          END IF
          END IF
        END IF
C
        DO IY=NTIN(ITARG,IPRT),NTEN(ITARG,IPRT)-1
          IG=IG+1
!          IF (IG.GT.NGITT) GOTO 999
          IF (IG.GT.NGITT) THEN
             WRITE (iunout,*) 'ERROR IN IF2COP: NGITT TOO SMALL '
             CALL EIRENE_EXIT_OWN(1)
             RETURN
          ENDIF

C  TESTEP, TISTEP: ZONE-CENTERED TEMPERATURE IN BOUNDARY ZONE (EV)
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
C FACTOR FOR SHEATH POTENTIAL: Vsh=Te*SHSTEP(ITARG,IG)
          SHSTEP(ITARG,IG)=DELTA_SHEATHXB(NPBS,IY)
C  TORL: TOROIDAL LENGTH (CM) AT TARGET SEGMENT IY: CENTER OF GRAVITY
          TORL(ITARG,IG)=2.*PIA*0.5*(XPOL(IY+1,NPES)+
     .                               XPOL(IY,NPES))
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
              DELY=RRSTEP(ITARG,IG+1)-RRSTEP(ITARG,IG)
              FLSTEP(IPLS,ITARG,IG)=0.
              IF (DELY.GT.0.) THEN
                FLSTEP(IPLS,ITARG,IG)=MAX(0._DP,ORI(ITARG,IG)*
     .                                FNIXB(NPBS,IY,IFL))*FL(IPLS)/DELY
C.......................................................................
C  CORRECT FOR INCLINED TARGETS: ADD FLUXES FROM SECOND DIRECTION
C  USE SIGN FROM "MAIN" CONTRIBUTION TO DECIDE ORIENTATION OF SEC. CONTR.
cdr             IF (FLSTEP(IPLS,ITARG,IG).GT.0.) THEN
cdr               FLSTEP(IPLS,ITARG,IG)=
cdr  .            FLSTEP(IPLS,ITARG,IG)+ABS(FNIX_YB(NPBS,IY,IFL))*
cdr  .                   FL(IPLS)/DELY
cdr             ENDIF
C  SET ION ENERGY FLUXES FROM B2-BOUNDARY CONDITIONS
!               delti_para=3
!               delte_para=0.5
!               delti_perp=2
!               delte_perp=0
!pb                tis=TISTEP(IPLSTI,ITARG,IG)
!pb                tes=TESTEP(ITARG,IG)
!pb                ELSTEP(IPLS,ITARG,IG) = ELSTEP(IPLS,ITARG,IG)+
!pb     .                                  FL(IPLS)/DELY*
!pb     .            (TIS*delti_perp*ABS(Fnix_yb(npbs,iy,ifl))+
!pb     .             TIS*delti_para*ABS(fnixb  (npbs,iy,ifl))+
!pb     .             TES*delte_para*ABS(fnixb  (npbs,iy,ifl)))

! 2008
!  only one of the next two is different from 0
!pb                delti_para=deltai_parxb(npbs,iy)
!pb                delti_perp=deltai_radxb(npbs,iy)
!  only one of the next two is different from 0
!pb                delte_para=deltae_parxb(npbs,iy)
!pb                delte_perp=deltae_radxb(npbs,iy)
!pb                tis=TISTEP(IPLSTI,ITARG,IG)
!pb                tes=TESTEP(ITARG,IG)
!pb                ELSTEP(IPLS,ITARG,IG) = ELSTEP(IPLS,ITARG,IG)+
!pb     .                                  FL(IPLS)/DELY*
!pb     .            (TIS*(delti_para+delti_perp)*ABS(fnixb(npbs,iy,ifl))+
!pb     .             TES*(delte_para+delte_perp)*ABS(fnixb(npbs,iy,ifl)))
C........................................................................

! 15.09.2010
!               ELSTEP(IPLS,ITARG,IG) = ELSTEP(IPLS,ITARG,IG)+
!    .                                  ABS(feixb(npbs,iy))/DELY
              ENDIF
            ENDIF

C  VXSTEP,VYSTEP,VZSTEP: SURFACE-CENTERED FLOW VELOCITY (CM/S)
C  NOTE: PV VECTOR IS CELL-CENTERED, BUT EXACT VECTOR CAN BE FOUND FROM
C        DATA FOR POLOIDAL POLYGON NPES
            IN=IY+(NPEC-1)*NR1ST
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

 3013     CONTINUE

!  present model:
!  within a given cell all incident ions (of all species) have same
!  energy E0B2 (delta function approximation) (NEMOD1=8 or NEMOD1=9)

          DELY=RRSTEP(ITARG,IG+1)-RRSTEP(ITARG,IG)
          IF (DELY.GT.0.) THEN
            fltt = sum(flstep(1:nplsi,itarg,ig))
            e0b2 = ABS(feixb(npbs,iy))/dely/(fltt+eps60)
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
        ENDDO
C
!     GOTO 3030
        ELSE
C
! 3020   CONTINUE
C
C  SECOND: SOURCES AT RADIAL (X) SURFACES
C
        DO IX=NTIN(ITARG,IPRT),NTEN(ITARG,IPRT)-1
          IG=IG+1
!          IF (IG.GT.NGITT) GOTO 999
          IF (IG.GT.NGITT) THEN
             WRITE (iunout,*) 'ERROR IN IF2COP: NGITT TOO SMALL '
             CALL EIRENE_EXIT_OWN(1)
             RETURN
          ENDIF
C  TESTEP, TISTEP: ZONE-CENTERED TEMPERATURE IN BOUNDARY ZONE (EV)
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
C FACTOR FOR SHEATH POTENTIAL: Vsh=Te*SHSTEP(ITARG,IG)
          SHSTEP(ITARG,IG)=DELTA_SHEATHYB(IX,NPBS)
C  TORL: TOROIDAL LENGTH (CM) AT TARGET SEGMENT IY: CENTER OF GRAVITY
          TORL(ITARG,IG)=2.*PIA*0.5*(XPOL(NPES,IX+1)+
     .                               XPOL(NPES,IX))
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
              DELX=RRSTEP(ITARG,IG+1)-RRSTEP(ITARG,IG)
              FLSTEP(IPLS,ITARG,IG)=0.
              IF (DELX.GT.0.) THEN
                FLSTEP(IPLS,ITARG,IG)=MAX(0._DP,ORI(ITARG,IG)*
     .                                FNIYB(IX,NPBS,IFL))*FL(IPLS)/DELX
C  SET DEFAULT ION ENERGY FLUXES FROM B2 BOUNDARY CONDITIONS
!                delti_para=3
!                delte_para=0.5
!                delti_perp=2
!                delte_perp=0
!  only one of the next two is different from 0
!pb                delti_para=deltai_paryb(ix,npbs)
!pb                delti_perp=deltai_radyb(ix,npbs)
!  only one of the next two is different from 0
!pb                delte_para=deltae_paryb(ix,npbs)
!pb                delte_perp=deltae_radyb(ix,npbs)
!pb                tis=TISTEP(IPLSTI,ITARG,IG)
!pb                tes=TESTEP(ITARG,IG)
!pb                ELSTEP(IPLS,ITARG,IG) = ELSTEP(IPLS,ITARG,IG) +
!pb     .                                  FL(IPLS)/DELX*
!pb     .            (TIS*delti_perp*ABS(Fniyb  (ix,npbs,ifl))+
!pb     .             TIS*delti_para*ABS(fniy_xb(ix,npbs,ifl))+
!pb     .             TES*delte_para*ABS(fniy_xb(ix,npbs,ifl)))

!pb                ELSTEP(IPLS,ITARG,IG) = ELSTEP(IPLS,ITARG,IG) +
!pb     .                                  FL(IPLS)/DELX*
!pb     .           (TIS*(delti_para+delti_perp)*ABS(Fniyb(ix,npbs,ifl))+
!pb     .            TES*(delte_para+delte_perp)*ABS(Fniyb(ix,npbs,ifl)))
! 15.09.2010
!                ELSTEP(IPLS,ITARG,IG) = ELSTEP(IPLS,ITARG,IG) +
!     .                                  ABS(Feiyb(ix,npbs))/DELX
              ENDIF
            ENDIF
C
C  VXSTEP,VYSTEP,VZSTEP: SURFACE-CENTERED FLOW VELOCITY (CM/S)
C  NOTE: PU VECTOR IS CELL-CENTERED, BUT EXACT VECTOR CAN BE FOUND FROM
C        RADIAL POLYGON NPES DATA
            IN=NPEC+(IX-1)*NR1ST
            PVXS=PVXN(IN)
            PVYS=PVYN(IN)
            PUXS=PUXN(IN)
            PUYS=PUYN(IN)
            PITB=RRB(NPBS,IY)
            UU=PITB*VPARYB(IX,NPBS,IFL)
            VXSTEP(IPLSV,ITARG,IG)=
     .            (PUXS*UU+PVXS*VRADYB(IX,NPBS,IFL))*V
            VYSTEP(IPLSV,ITARG,IG)=
     .            (PUYS*UU+PVYS*VRADYB(IX,NPBS,IFL))*V
            VZSTEP(IPLSV,ITARG,IG)=
     .            (SQRT(1.-PITB**2)*VPARYB(IX,NPBS,IFL))*V

 3023     CONTINUE

!  present model:
!  within a given cell all incident ions (of all species) have same
!  energy E0B2 (delta function approximation) (NEMOD1=8 or NEMOD1=9)

          DELX=RRSTEP(ITARG,IG+1)-RRSTEP(ITARG,IG)
          IF (DELX.GT.0.) THEN
            fltt = sum(flstep(1:nplsi,itarg,ig))
            e0b2 = ABS(feiyb(ix,npbs))/delx/(fltt+eps60)
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
        ENDDO
!     3030   CONTINUE
        ENDIF
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
      INDIM(1,ITARG)=4
      IF (INDSRC(ITARG).NE.6) THEN
        I34=EIRENE_IDEZ(INT(SORLIM(1,ITARG)),3,3)
        SORLIM(1,ITARG)=I34*100+04   !sample with step fct.
      ELSEIF (INDSRC(ITARG).EQ.6) THEN
C  SORLIM DEFAULT WAS 0.D0
        SORLIM(1,ITARG)=0204
      ENDIF
      SORIND(1,ITARG)=ITARG
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
        NSPEZ(ITARG)=-1
        SORIFL(1,ITARG)=NIFLG(ITARG,1)
        SORWGT(1,ITARG)=1.
C  USE ENERGY FLUXES SPECIFIED HERE, IE., SORENE, SORENI ARE REDUNDANT
!pb 14072011        NEMODS(ITARG)=9
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
! 3999 CONTINUE
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
          NPBS=NPES-1
          IF (NEM.EQ.3.OR.NEM.EQ.5.OR.NEM.EQ.7) THEN
            DO 6005 IPL=1,NPLSI
              IFL=IFLB(IPL)
              VPX=VXSTEP(IPL,ITARG,IG)
              VPY=VYSTEP(IPL,ITARG,IG)
              VPZ=VZSTEP(IPL,ITARG,IG)
              VP(IPL)=SQRT(VPX**2+VPY**2+VPZ**2)
              DI(IPL)=DISTEP(IPL,ITARG,IG)
              ZI(IPL)=ZISTEP(IPL,ITARG,IG)
 6005       CONTINUE
            TE=TESTEP(ITARG,IG)
            CUR=0.
            GAMMA=0.
            ESHT(ITARG,IG)=EIRENE_SHEATH(TE,DI,VP,ZI,GAMMA,CUR,
     .                         NPLSI,-ITARG)
          ELSE IF (NEM == 9) THEN
            TE=TESTEP(ITARG,IG)
            ESHT(ITARG,IG)=DELTA_SHEATHXB(NPBS,IY)*TE
          END IF
        ELSEIF (IGSTEP(ITARG,IG).LT.200000) THEN
          IX=IPSTEP(ITARG,IG)
          NPES=IGSTEP(ITARG,IG)-100000
          NPBS=NPES-1
          NEM=IABS(NEMODS(ITARG))
          IF (NEM.EQ.3.OR.NEM.EQ.5.OR.NEM.EQ.7) THEN
C  in order to find sheath potential, we need ALL plasma particle flux components
            DO 6006 IPL=1,NPLSI
              IFL=IFLB(IPL)
              VPX=VXSTEP(IPL,ITARG,IG)
              VPY=VYSTEP(IPL,ITARG,IG)
              VPZ=VZSTEP(IPL,ITARG,IG)
              VP(IPL)=SQRT(VPX**2+VPY**2+VPZ**2)
              DI(IPL)=DISTEP(IPL,ITARG,IG)
              ZI(IPL)=ZISTEP(IPL,ITARG,IG)             
 6006       CONTINUE
            TE=TESTEP(ITARG,IG)
            CUR=0.
            GAMMA=0.
            ESHT(ITARG,IG)=EIRENE_SHEATH(TE,DI,VP,ZI,GAMMA,CUR,
     .                           NPLSI,-ITARG)
          ELSE IF (NEM == 9) THEN
            TE=TESTEP(ITARG,IG)
            ESHT(ITARG,IG)=DELTA_SHEATHYB(IX,NPBS)*TE
          END IF
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
            WRITE (iunout,*) 'IPLS,ITARG,IG,MACH '
     .                       ,IPLS,ITARG,IG,VTEST
C           WRITE (iunout,*) 'POL., TOR., RAD. (CM/S) ',PM1,VPZ,VR
            CALL EIRENE_LEER(1)
C
          ELSEIF (IGSTEP(ITARG,IG).LT.200000) THEN
C  CHECK BOHM CRITERION AT "POLOIDAL" TARGET SURFACE COMPONENTS
            IX=IPSTEP(ITARG,IG)
            NPES=IGSTEP(ITARG,IG)-100000
            VT=SQRT(2.*TISTEP(IPLSTI,ITARG,IG)/BMASS(IPLS))*CVEL2A
C  VELOCITY COMPONENT NORMAL TO POLOIDAL TARGET SURFACE
C  I.E., POLOIDAL COMPONENT V-POL
            PM1=(PLNX(NPES,IX)*VXSTEP(IPLSV,ITARG,IG)+
     .           PLNY(NPES,IX)*VYSTEP(IPLSV,ITARG,IG))*OR
C  VELOCITY COMPONENT PARALLEL TO POLOIDAL TARGET SURFACE
C  I.E., RADIAL PLUS TOROIDAL COMPONENT, V-RAD + V-TOR
            VPX=VXSTEP(IPLSV,ITARG,IG)-PM1*PLNX(NPES,IX)*OR
            VPY=VYSTEP(IPLSV,ITARG,IG)-PM1*PLNY(NPES,IX)*OR
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
            WRITE (iunout,*) 'IPLS,ITARG,IG,MACH ',IPLS,ITARG,IG,VTEST
C           WRITE (iunout,*) 'POL., TOR., RAD. ',PM1,VPZ,VR
C           CALL EIRENE_LEER(1)
C
          ENDIF
C
C  BOHM CRITERION CHECK DONE
C
C  ELTEST: TOTAL ION ENERGY FLUX ONTO TARGET:EMAXW + ESHET

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
!        GOTO 6011
! 6010   CONTINUE
C  TO BE WRITTEN
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
! 6300 CONTINUE
C
C  SET SOME OTHER DATA SPECIFIC FOR EIRENE CODE REQUIREMENTS
C  STATEMENT NO. 6500 ---> 6999
C
 6500 CONTINUE

      DEALLOCATE (TORL)
      DEALLOCATE (ESHT)
      DEALLOCATE (ELTEST)
      DEALLOCATE (ORI)
C
C
      RETURN

      END SUBROUTINE
!  999 CONTINUE
!      WRITE (iunout,*) 'ERROR IN IF2COP: NGITT TOO SMALL '
!      CALL EIRENE_EXIT_OWN(1)
!      RETURN
C
C
!      SUBROUTINE EIRENE_IF3COP(IENTRY,LSTP,
!     .                    IFRST,ISTRAA,ISTRAE,NEW_ITER)
C
C
!      WRITE (iunout,*) ' IF3COP IS CALLED, ISTRAA,ISTRAE '
!      WRITE (iunout,*) ISTRAA,ISTRAE
!      LSHORT=.FALSE.
!      LSTP3=.TRUE.
!      LSTOP=.TRUE.
!      IFIRST=0
!      NDXY=(NDXA-1)*NR1ST+NDYA
!      GOTO 99992
C
!      ENTRY EIRENE_INTER3(LSTP,IFRST,ISTRAA,ISTRAE,NEW_ITER)
C
C  ENTRY FOR SHORT CYCLE FROM SUBR. EIRSRT
C
C  IFIRST=0: RESTORE DATA FROM A PREVIOUS EIRENE RUN, SET REFERENCE
C            DATA FOR "STOP-CRITERION" SNIS,SEES,SEIS
C  IFIRST>0: MODIFY SOURCE TERMS ACCORDING TO NEW PLASMA CONDITIONS,
C            COMPARE INTEGRALS WITH SNIS,...., AND DECIDE TO STOP OR
C            CONTINUE SHORT CYCLE (LSTOP)
C
!      LSHORT=.TRUE.
!      LSTP3=LSTP
!      LSTOP=LSTP
!      IFIRST=IFRST
!      NDXY=(NDXA-1)*NR1ST+NDYA
C
!99992 CONTINUE

      SUBROUTINE EIRENE_INTF_GENERIC(LSTP,IFRST,ISTRAA,ISTRAE,NEW_ITER)
      LOGICAL, INTENT(INOUT) :: LSTP
      INTEGER, INTENT(IN) :: ISTRAA, ISTRAE, NEW_ITER, IFRST
      REAL(DP) :: SEES0(NSTRA), SEIS0(NSTRA)
      REAL(DP) :: DUMMY(0:NDXP,0:NDYP)
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
        ALLOCATE (scpveii(NSTRA))

        CALL EIRENE_ALLOC_BRASPOI
        CALL EIRENE_ALLOC_EIRBRA(NDX,NDY,NFL,NSTRA)
C
        RESSNI = 0._DP
        RESSMO = 0._DP
        RESSEE = 0._DP
        RESSEI = 0._DP
c
        scpveii = 0._dp
      ENDIF
C
      IF (.NOT.LSHORT) THEN
        RESSNI(ISTRAA:ISTRAE,:) = 0._DP
        RESSMO(ISTRAA:ISTRAE,:) = 0._DP
        RESSEE(ISTRAA:ISTRAE) = 0._DP
        RESSEI(ISTRAA:ISTRAE) = 0._DP
      ENDIF

      if (.not.allocated(cpv_cmp)) then
        allocate (cpv_cmp(ncpv,nrtal,nstra))
        cpv_cmp = 0._dp
      end if

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
     .      'NO DATA RETURNED TO PLASMA CODE FOR THIS STRATUM'
          GOTO 7999
        ELSEIF (ISTRAI.GT.NTARGI) THEN
          FLXI=1.

C  IF THE SOURCE STRENGTH IS TO BE CHANGED DURING THE SHORT CYCLE (E.G.: VOL-REC)
C  THEN FLXEIR HAS TO BE RESET TO SCALE TO NEW SOURCE STRENGTH DURING SHORT CYCLE
          FLXEIR(ISTRAI)=1._DP
        ENDIF
C
C  FIRSTLY INITIALIZE SOURCE TERM ARRAYS

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
     .          (SEINW(IN,IPLS)-RTIS%SEIOD(IN,IPLS))*ELCHA
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

cdr  we are still in stratum istra.
cdr  what is this next  'lzden' option doing ? Why here, after short loop corrections

        IF (LZDEN) THEN
          IF (IZDEN(ISTRAI) == 0) THEN
            DO IN=1,NSBOX_TAL
              LZDENA(IN,ISTRAI) = SUM(ABS(PDENA(1:NATMI,IN))) < EPS10
              LZDENM(IN,ISTRAI) = SUM(ABS(PDENM(1:NMOLI,IN))) < EPS10
              LZDENI(IN,ISTRAI) = SUM(ABS(PDENI(1:NIONI,IN))) < EPS10
            END DO
            IZDEN(ISTRAI) = 1
          ELSE
            DO IN=1,NSBOX_TAL
              IF (LZDENA(IN,ISTRAI)) THEN
C  NO ATOMS IN THIS STRATUM
                PAPL(1:NPLSI,IN) = 0._DP
                MAPL(1:NPLSI,IN) = 0._DP
                EAPL(1:NPLSI,IN) = 0._DP
                EAEL(IN) = 0._DP
              END IF
              IF (LZDENM(IN,ISTRAI)) THEN
C  NO MOLECULES IN THIS STRATUM
                PMPL(1:NPLSI,IN) = 0._DP
                MMPL(1:NPLSI,IN) = 0._DP
                EMPL(1:NPLSI,IN) = 0._DP
                EMEL(IN) = 0._DP
              END IF
              IF (LZDENI(IN,ISTRAI)) THEN
C  NO TEST IONS IN THIS STRATUM
                PIPL(1:NPLSI,IN) = 0._DP
                MIPL(1:NPLSI,IN) = 0._DP
                EIPL(1:NPLSI,IN) = 0._DP
                EIEL(IN) = 0._DP
              END IF
            END DO
          END IF   ! IZDEN(ISTRAI) == 0
        END IF   ! LZDEN
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
          DO 7473 IPLS=1,NPLSI
            JPLS = NSPEZ(ISTRAI)
            IF ((JPLS > 0) .AND. (JPLS <= NPLSI) .AND.
     .          (IPLS /= JPLS)) CYCLE
            IPLSTI = MPLSTI(IPLS)
            DO ISR=1, NSRFSI(ISTRAI)
              ISTEP = SORIND(ISR,ISTRAI)
              DO 7472 IIRC=1,NPRCI(IPLS)
                IRRC=LGPRC(IPLS,IIRC)
                IF ((ISTEP > 0) .AND. (ISTEP /= IRRC)) CYCLE
                SUMN=0.0
                SUMM=0.0
                SUMEI=0.0
                SUMEE=0.0
                IR1 = MAX(1,INGRDA(ISR,ISTRAI,1))
                IR2 = MIN(NR1ST,INGRDE(ISR,ISTRAI,1))
!pb                DO 7471 IR=1,NR1ST-1
                DO 7471 IR = IR1, IR2-1
                DO K=1,NPPLG
                DO IP=NPOINT(1,K),NPOINT(2,K)-1
                  IF ((IP < INGRDA(ISR,ISTRAI,2)) .OR.
     .                (IP >= INGRDE(ISR,ISTRAI,2))) CYCLE
                  IN=(IP-1)*NR1ST+IR
                  INC=NCLTAL(IN)
C  EXCLUDE DEAD CELLS (GRID CUTS, ISOLATED CELLS FROM COUPLE_.., ETC)
C  EXCLUDE IPLS-VACUUM CELLS
                  IF (NSTGRD(IN).EQ.0.AND..NOT.LGVAC(IN,IPLS)) THEN
                    IF (NSTORDR >= NRAD) THEN
                     RECADD=-TABRC1(IRRC,IN)*DIIN(IPLS,IN)*ELCHA
                     EEADD=  EELRC1(IRRC,IN)*DIIN(IPLS,IN)*ELCHA
                    ELSE
                     RECADD=-EIRENE_FTABRC1(IRRC,IN)*DIIN(IPLS,IN)*ELCHA
                     EEADD=  EIRENE_FEELRC1(IRRC,IN)*DIIN(IPLS,IN)*ELCHA
                    END IF
                    PPPL_COP(IPLS,INC)=PPPL_COP(IPLS,INC)+RECADD
                    SUMN=SUMN+RECADD*VOL(IN)
                    PIADD=0._DP
                    IF (LPARMOM) PIADD=PARMOM(IPLS,IN)*RECADD
                    MPPL_COP(IPLS,INC)=MPPL_COP(IPLS,INC)+PIADD
                    SUMM=SUMM+PIADD*VOL(IN)
                    EIADD=1.5*TIIN(IPLSTI,IN)*RECADD
                    IF (LEDRIFT) EIADD=EIADD+EDRIFT(IPLS,IN)*RECADD
                    EPPL_COP(IPLS,INC)=EPPL_COP(IPLS,INC)+EIADD
                    SUMEI=SUMEI+EIADD*VOL(IN)
                    EPEL_COP(INC)=EPEL_COP(INC)+EEADD
                    SUMEE=SUMEE+EEADD*VOL(IN)
                  END IF
                END DO
                END DO
 7471           CONTINUE
                RECTOT = RECTOT + SUMN
                WRITE (iunout,*) 'IPLS,IRRC ',IPLS,IRRC
                CALL EIRENE_MASR4('SUMN, SUMM, SUMEI, SUMEE        ',
     .                      SUMN,SUMM,SUMEI,SUMEE)
 7472         CONTINUE
            END DO
 7473     CONTINUE
cdr
CC SUMN_OLD=WTOTP ???
          IF (.NOT.LSHORT) SUMN_OLD=RECTOT
csw 08jun2010
          sumn_old=rectot

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
     .             EPEL_COP(1:NSBOX_TAL) - EPEODA (1:NSBOX_TAL)
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
cdr  copv(icp+1,...) to copv(icp5+1,...)  linear algebraic tallies
cdr                                       scored in upfcop, so also their variances
cdr                                       can be set without need for covariance matrices
        icp=nplsi
        icp2=2*nplsi
        icp3=3*nplsi
cdr
        icp4=4*nplsi
        icp5=5*nplsi


        scpveii(istrai) = 0._dp

        cpv_cmp(:,:,istrai) = copv (:,:)

cdr  fill bulk particle source rate sni(...ifl) from all contributing
cdr  test particle sources papl,pmpl,pipl,pppl (...,ipls)
        DO 7510 IFL=1,NFLA
          CHPS(IFL)=0.
          SNIS(IFL)=0.
          CHMOS(IFL)=0.
          SMOS(IFL)=0.

          DO 7515 IPLS=1,NPLSI
            IF (IFLB(IPLS).NE.IFL) GOTO 7515
            IPLSV=MPLSV(IPLS)
c  ipls contributes to plasma code species ifl
            DO 7520 IX=1,NDXA
              IF (LLCUT(IX)) CYCLE
              DO 7530 IY=1,NDYA
                IN=IY+(IX-1)*NR1ST
                INC=NCLTAL(IN)
                SNICL=(PAPL(IPLS,INC)+PMPL(IPLS,INC)+PIPL(IPLS,INC)+
     .                 PPPL_COP(IPLS,INC))*VOLTAL(INC)*FLX_EIR
                SNI(IX,IY,IFL,ISTRAI)=SNI(IX,IY,IFL,ISTRAI)+SNICL
                SNIS(IFL)=SNIS(IFL)+SNICL
                CHPS(IFL)=CHPS(IFL)+CHPM(IPLS,INC)*VOLTAL(INC)
 7530         CONTINUE
 7520       CONTINUE

!pb 22012013 copv
cdr  build alternative source rates, from corresponding copv tallies
cdr  scored in upfcop.
cdr  tbd:  check storage on copv tallies, ncpv ??
            DO IX=1,NDXA
              IF (LLCUT(IX)) CYCLE
              DO IY=1,NDYA
                INN=IY+(IX-1)*NR1ST
                IN=NCLTAL(INN)
cdr  scale the particle sources from copv(icp+ipls)
cdr  then add pppl
                cfac = (PAPL(IPLS,IN)+PMPL(IPLS,IN)+PIPL(IPLS,IN)) /
     .                 (cpv_cmp(icp+ipls,in,istrai) + eps60)
                cpv_cmp(icp+ipls,in,istrai)=
     .                  (cpv_cmp(icp+ipls,in,istrai)*cfac +
     .                  PPPL_COP(IPLS,IN)) * VOLTAL(IN)*FLX_EIR
cdr  is this now any different from sni set above?


cdr  add pppl contribution to internal energy sources rate
                bv = 0._dp
                if (lbvin) bv = bvin(iplsv,inn)
                cpv_cmp(icp4+ipls,in,istrai)=
     .                  cpv_cmp(icp4+ipls,in,istrai) + 
     .                  cvrssp(ipls)*bv**2*PPPL_COP(IPLS,IN)
              end do  ! iy
            end do    ! ix



            IF (.NOT.LSHORT) THEN

cdr   try to find stat. variance for particle balance sources
cdr   ntalm is a copv tally.
              istat_cop = 0
              do i = 1, nsigvi
                if ((iih(i) == ntalm).and.(igh(i) == NPLSI+IPLS)) then
                  istat_cop = i
                  exit
                end if
              end do

              if (istat_cop > 0) then
                DO IX=1,NDXA
                  IF (LLCUT(IX)) CYCLE
                  DO IY=1,NDYA
                    IN=IY+(IX-1)*NR1ST
                    INC=NCLTAL(IN)
                    SNIRES=(PAPL(IPLS,INC)+PMPL(IPLS,INC)+
     .                      PIPL(IPLS,INC))*VOLTAL(INC)*FLX_EIR
                    RESSNI(ISTRAI,IFL)=RESSNI(ISTRAI,IFL)+
     .                                 ABS(SIGMA(ISTAT_COP,IN)*
     .                                 SNIRES/100.D0)
                  END DO
                END DO
              end if

            END IF

cdr   particle sources done.

cdr   next: dwell on parallel momentum sources. still inside ifl and ipls loop
cdr   ipls contributes to plasma code species ifl

            DO 7536 IX=1,NDXA
              IF (LLCUT(IX)) CYCLE
              DO 7533 IY=1,NDYA
                IN=IY+(IX-1)*NR1ST
                INC=NCLTAL(IN)
                SIGNUM=1._DP
                IF (LBVIN) SIGNUM=SIGN(1._DP,BVIN(IPLSV,IN))
                SMOCL=(MAPL(IPLS,INC)+MMPL(IPLS,INC)+MIPL(IPLS,INC)+
     .                 MPPL_COP(IPLS,INC))*
     .                 VOLTAL(INC)*1.D-5*SIGNUM*FLX_EIR
                SMO(IX,IY,IFL,ISTRAI)=SMO(IX,IY,IFL,ISTRAI)+SMOCL
                SMOS(IFL)=SMOS(IFL)+SMOCL
                CHMOS(IFL)=CHMOS(IFL)+CHMOM(IPLS,INC)*VOLTAL(INC)
 7533         CONTINUE
 7536       CONTINUE

!pb 22012013 copv
cdr  build alternative source rates, from corresponding copv tallies
cdr  scored in updlin.
cdr  tbd:  check storage on copv tallies, ncpv ??
            DO IX=1,NDXA
              IF (LLCUT(IX)) CYCLE
              DO IY=1,NDYA
                IN=IY+(IX-1)*NR1ST
                INC=NCLTAL(IN)
                SIGNUM=1._DP
                IF (LBVIN) SIGNUM=SIGN(1._DP,BVIN(IPLSV,IN))
                cfac = (MAPL(IPLS,INC)+MMPL(IPLS,INC)+MIPL(IPLS,INC))/
     .                 (cpv_cmp(icp2+ipls,inc,istrai) + eps60)
                cpv_cmp(icp2+ipls,inc,istrai)=
     .                 (cpv_cmp(icp2+ipls,inc,istrai)*cfac +
     .                  MPPL_COP(IPLS,INC))*
     .                  VOLTAL(INC)*1.D-5*SIGNUM*FLX_EIR
!pb 30012013 sei internal
                bv = 0._dp
                if (lbvin) bv = bvin(iplsv,in)
                cpv_cmp(icp3+3,inc,istrai)=cpv_cmp(icp3+3,inc,istrai) -
     .                  bv*MPPL_COP(IPLS,INC)*SIGNUM*
     .                  cveli2/amua*2._DP
              end do
            end do

            IF (.NOT.LSHORT) THEN

cdr   try to find stat. variance for momentum balance sources
cdr   ntalm is a copv tally.
              istat_cop = 0
              do i = 1, nsigvi
                if ((iih(i) == ntalm).and.(igh(i) == 2*NPLSI+IPLS)) then
                  istat_cop = i
                  exit
                end if
              end do

              if (istat_cop > 0) then
                DO IX=1,NDXA
                  IF (LLCUT(IX)) CYCLE
                  DO IY=1,NDYA
                    IN=IY+(IX-1)*NR1ST
                    INC=NCLTAL(IN)
                    SIGNUM=1._DP
                    IF (LBVIN) SIGNUM=SIGN(1._DP,BVIN(IPLSV,IN))
                    SMORES=(MAPL(IPLS,INC)+MMPL(IPLS,INC)+
     .                 MIPL(IPLS,INC))*VOLTAL(INC)*1.D-5*SIGNUM*FLX_EIR
                    RESSMO(ISTRAI,IFL)=RESSMO(ISTRAI,IFL)+
     .                                 ABS(SIGMA(ISTAT_COP,INC)*
     .                                 SMORES/100.D0*1.D5)
                  END DO
                END DO
              end if

            END IF
!pb 7539        CONTINUE
 7515   CONTINUE
 7510   CONTINUE

C
        CHEES=0.
        SEES=0.
        DO 7540 IX=1,NDXA
          IF (LLCUT(IX)) CYCLE
          DO 7545 IY=1,NDYA
            INN=IY+(IX-1)*NR1ST
            IN=NCLTAL(INN)
            SEE(IX,IY,ISTRAI)=(EAEL(IN)+EMEL(IN)+
     .                         EIEL(IN)+EPEL_COP(IN))*
     .                         VOLTAL(IN)*FLX_EIR
            CHEES=CHEES+CHEEM(IN)*VOLTAL(IN)
            SEES=SEES+SEE(IX,IY,ISTRAI)
            SEE(IX,IY,ISTRAI)=SEE(IX,IY,ISTRAI)*ELCHA
 7545     CONTINUE
 7540   CONTINUE
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
     .                            VOLTAL(IN)*FLX_EIR
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

!pb 22012013 copv
        DO IX=1,NDXA
          IF (LLCUT(IX)) CYCLE
          DO IY=1,NDYA
            INN=IY+(IX-1)*NR1ST
            IN=NCLTAL(INN)
!pb 22012013 copv
cdr  electron energy source rate, now on copv icp5+1
            cfac = (EAEL(IN)+EMEL(IN)+EIEL(IN)) /
     .             (cpv_cmp(icp3+1,in,istrai) + eps60)
            cpv_cmp(icp5+1,in,istrai)=(cpv_cmp(icp5+1,in,istrai)*cfac +
     .              EPEL_COP(IN)) * VOLTAL(IN)*ELCHA
cdr  this is now identical to see above ?

cdr  ion energy source rate,
cdr  to be redone:  eamisum and eplsum to be made species-dependent again????
cdr  i.e. eppl_cop(ipls,in),....

            eamisum = 0._dp
            eplsum = 0._dp
            DO IFL=1,NFLA
              DO IPLS=1,NPLSI
                IF (IFLB(IPLS).EQ.IFL) THEN
                  eamisum = eamisum +
     .                      EAPL(IPLS,IN)+EMPL(IPLS,IN)+EIPL(IPLS,IN)
                  eplsum = eplsum + EPPL_COP(IPLS,IN)
                END IF
              END DO
            END DO
cdr  total ion energy source rate, now on copv icp3+ipls
            cfac = EAMISUM / (cpv_cmp(icp3+ipls,in,istrai) + eps60)
            cpv_cmp(icp3+ipls,in,istrai)=(cpv_cmp(icp3+ipls,in,istrai)*
     .               cfac + EPLSUM) * VOLTAL(IN)*ELCHA
cdr  this is now identical to sei above ?

!pb 30012013 sei internal energy
            cpv_cmp(icp4+ipls,in,istrai)=(cpv_cmp(icp4+ipls,in,istrai) +
     .              EPLSUM) * VOLTAL(IN)*ELCHA
            scpveii(istrai)=scpveii(istrai)+cpv_cmp(icp4+ipls,in,istrai)
          end do
        end do

!pb        cpv_cmp(icp3+1,:,istrai) = cpv_cmp(icp3+1,:,istrai) * flxi
!pb        cpv_cmp(icp3+2,:,istrai) = cpv_cmp(icp3+2,:,istrai) * flxi
!pb        cpv_cmp(icp3+3,:,istrai) = cpv_cmp(icp3+3,:,istrai) / elcha
!pb

        IF (.NOT.LSHORT) THEN

cdr   try to find stat. variance for electron energy balance sources
cdr   ntalm is a copv tally.
          istat_cop = 0
          do i = 1, nsigvi
            if ((iih(i) == ntalm).and.(igh(i) == 3*NPLSI+1)) then
              istat_cop = i
              exit
            end if
          end do

          if (istat_cop > 0) then
            DO IX=1,NDXA
              IF (LLCUT(IX)) CYCLE
              DO IY=1,NDYA
                INN=IY+(IX-1)*NR1ST
                IN=NCLTAL(INN)
                SEERES=(EAEL(IN)+EMEL(IN)+EIEL(IN))*VOLTAL(IN)*FLX_EIR
                RESSEE(ISTRAI)=RESSEE(ISTRAI)+
     .                         ABS(SIGMA(ISTAT_COP,IN)*
     .                         SEERES/100.D0)
              END DO
            END DO
          end if

cdr   try to find stat. variance for ion energy balance sources
cdr   ntalm is a copv tally.
          istat_cop = 0
          do i = 1, nsigvi
            if ((iih(i) == ntalm).and.(igh(i) == 3*NPLSI+2)) then
              istat_cop = i
              exit
            end if
          end do

          if (istat_cop > 0) then
            DO IX=1,NDXA
              IF (LLCUT(IX)) CYCLE
              DO IY=1,NDYA
                INN=IY+(IX-1)*NR1ST
                IN=NCLTAL(INN)
                eamisum = 0._dp
                DO IFL=1,NFLA
                  DO IPLS=1,NPLSI
                    IF (IFLB(IPLS).EQ.IFL) THEN
                       eamisum = eamisum +
     .                      EAPL(IPLS,IN)+EMPL(IPLS,IN)+EIPL(IPLS,IN)
                    END IF
                  END DO
                END DO
                SEIRES=EAMISUM*VOLTAL(IN)*FLX_EIR
                RESSEI(ISTRAI)=RESSEI(ISTRAI)+
     .                         ABS(SIGMA(ISTAT_COP,IN)*
     .                         SEIRES/100.D0)
              END DO
            END DO
          end if
        END IF

        WRITE (iunout,*) 'RECYCLING SOURCE FROM IF3COP ',ISTRAI
        WRITE (iunout,8888) sum(SNIS(1:nfla)), seis, sees
 8888 FORMAT (3E14.6)
C
C   NEXT:
C   IF LSHORT: CRITERION TO STOP SHORT CYCLE,
C   IF NOT LSHORT: ONLY RESCALE SURFACE SOURCE STRATA
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
        ENDIF

        NLSRON(ISTRAI) = (.NOT.LTEST) .OR. (XMCP(ISTRAI) <= 2) .OR.
     .                   NLVOL(ISTRAI)
C
        ELSEIF (.NOT.LSHORT) THEN
C
          DO 7560 IX=0,NDXA+1
            DO 7565 IY=0,NDYA+1
              SEE(IX,IY,ISTRAI)=SEE(IX,IY,ISTRAI)*FLXI
              SEI(IX,IY,ISTRAI)=SEI(IX,IY,ISTRAI)*FLXI

CDR  NO IDEA WHAT THIS IS.....
              if (ix*iy > 0) then
                INN=IY+(IX-1)*NR1ST
                IN=NCLTAL(INN)
                cpv_cmp(icp3+1,in,istrai) =
     .                  cpv_cmp(icp3+1,in,istrai) * flxi
                cpv_cmp(icp3+2,in,istrai) =
     .                  cpv_cmp(icp3+2,in,istrai) * flxi
                cpv_cmp(icp3+3,in,istrai) =
     .                  cpv_cmp(icp3+3,in,istrai) / elcha
              end if
CDR ???
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

CDR  NO IDEA WHAT THIS IS.....
          DO IPLS=1,NPLSI
            DO IX=0,NDXA+1
              DO IY=0,NDYA+1
                if (ix*iy > 0) then
    1             INN=IY+(IX-1)*NR1ST
                  IN=NCLTAL(INN)
                  cpv_cmp(icp+ipls,in,istrai)=
     .                    cpv_cmp(icp+ipls,in,istrai)*flxi
                  cpv_cmp(icp2+ipls,in,istrai)=
     .                    cpv_cmp(icp2+ipls,in,istrai)*flxi
                end if
              end do
            end do
          end do
CDR ???
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
      END SUBROUTINE
C
!     ENTRY EIRENE_IF4COP
!HJL REplaces entry with subroutine      
      SUBROUTINE EIRENE_IF4COP
C
cdr for species-dependent global particle balance
      REAL(DP) :: SFNISY(NFL),SFNINY(NFL),SFNIWX(NFL),SFNIEX(NFL)
      REAL(DP) :: SSN(NFL),SSNI(NFL),BALANN(NFL),TOTN(NFL),RN(NFL)
      REAL(DP) :: SFNIT(0:NSTEP,NFL), SFEIT(0:NSTEP),
     R            SFEET(0:NSTEP), SHEAE(0:NSTEP), SHEAI(0:NSTEP)

      NREC11=NOUTAU
      OPEN (UNIT=11,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=8*NREC11)
      IRC=3
C  WRITE RCCPL
      WRITE (11,REC=IRC) RCCPL
      IF (TRCINT.OR.TRCFLE)
     .    WRITE (iunout,*) 'WRITE 11  RCCPL,   IRC= ',IRC
#ifdef CHECKBIN
      write (111+ifoff,*) ' RCCPL '
      write (111+ifoff,*) RCCPL
#endif
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
#ifdef CHECKBIN
      write (111+ifoff,*) ' ICCPL1 '
      write (111+ifoff,*) IHELP
#endif
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
#ifdef CHECKBIN
      write (111+ifoff,*) IHELP
#endif
        IF (TRCINT.OR.TRCFLE)
     .      WRITE (iunout,*) 'WRITE 11  ICCPL1,  IRC= ',IRC
      END IF
      DEALLOCATE (IHELP)
C  WRITE ICCPL2
      IRC=IRC+1
      WRITE (11,REC=IRC) ICCPL2
#ifdef CHECKBIN
      write (111+ifoff,*) ' ICCPL2 '
      write (111+ifoff,*) ICCPL2
#endif
      IF (TRCINT.OR.TRCFLE)
     .    WRITE (iunout,*) 'WRITE 11  ICCPL2,  IRC= ',IRC
      IRC=IRC+1
      WRITE (11,REC=IRC) LCCPL
#ifdef CHECKBIN
      write (111+ifoff,*) ' LCCPL '
      write (111+ifoff,*) LCCPL
#endif
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
10113 CONTINUE
C
      SFNISY=SFNISY*ELCHA
C
      WRITE (37,*) 'NON-RECYCLING FLUXES FROM SOUTH EDGE '
      WRITE (37,8888) SFNISY,SFEISY,SFEESY
 8888 FORMAT (3E14.6)
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
10118 CONTINUE
C
      SFNINY=SFNINY*ELCHA
C
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
      WRITE (37,*) 'NON-RECYCLING FLUXES TO EAST EDGE '
      WRITE (37,8888) SFNIEX,SFEIEX,SFEEEX
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

        write (37,*)
        write (37,*) 'target no. ',i
        write (37,*) 'iprt, ix, iy, fni '

        DO IPRT=1,NTGPRT(I)
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
              DO 10131 IF=NSPZI(I,IPRT),NSPZE(I,IPRT)
                IF (NINCT(I,IPRT)*FNIXB(NPBS,IY,IF).GT.0) THEN
                SFNIT(I,IF)=SFNIT(I,IF)-
     .                   NINCT(I,IPRT)*FNIXB(NPBS,IY,IF)

cdr sheath contributions: count negative for electrons, positive for ions
cdr unfinished:  need to account for charge state of ion species IFL
                SHEAE(I)=SHEAE(I)+TEB(NPBC,IY)*
     .           NINCT(I,IPRT)*FNIXB(NPBS,IY,IF)*
     .           (-DELTA_SHEATHXB(NPBS,IY))
                SHEAI(I)=SHEAI(I)+TEB(NPBC,IY)*
     .           NINCT(I,IPRT)*FNIXB(NPBS,IY,IF)*
     .           DELTA_SHEATHXB(NPBS,IY)
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
              DO 10136 IF=NSPZI(I,IPRT),NSPZE(I,IPRT)
                IF (NINCT(I,IPRT)*FNIYB(IX,NDT(I,IPRT),IF).GT.0.) THEN
                SFNIT(I,IF)=SFNIT(I,IF)-
     .                   NINCT(I,IPRT)*FNIYB(IX,NDT(I,IPRT),IF)

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
        DO 10142 IX=1,NDXA
           DO 10140 IY=1,NDYA
             DO 10141 IFL=1,NFLA
               SSN(IFL)=SSN(IFL)+SNI(IX,IY,IFL,ISTRA)
10141        CONTINUE
             SSI=SSI+SEI(IX,IY,ISTRA)
             SSE=SSE+SEE(IX,IY,ISTRA)
10140     CONTINUE
10142   CONTINUE
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
        WRITE (iunout,*) ' NON-RECYCLING FLUXES AT SOUTH EDGE '
        CALL EIRENE_MASR2(' SFEISY,SFEESY  ',SFEISY,SFEESY)
          DO IFL=1,NFLA
            WRITE(iunout,'(A,I0,A,ES12.4)') 'SFNISY(IFL=',IFL,') ',
     .                                       SFNISY(IFL)
        ENDDO
        WRITE (iunout,*) ' NON-RECYCLING FLUXES AT NORTH EDGE'
        CALL EIRENE_MASR2(' SFEINY,SFEENY  ',SFEINY,SFEENY)
          DO IFL=1,NFLA
            WRITE(iunout,'(A,I0,A,ES12.4)') 'SFNINY(IFL=',IFL,') ',
     .                                       SFNINY(IFL)
        ENDDO
        WRITE (iunout,*) ' NON-RECYCLING FLUXES AT WEST EDGE '
        CALL EIRENE_MASR2(' SFEIWX,SFEEWX  ',SFEIWX,SFEEWX)
          DO IFL=1,NFLA
            WRITE(iunout,'(A,I0,A,ES12.4)') 'SFNIWX(IFL=',IFL,') ',
     .                                       SFNIWX(IFL)
        ENDDO
        WRITE (iunout,*) ' NON-RECYCLING FLUXES AT EAST EDGE '
        CALL EIRENE_MASR2(' SFEIEX,SFEEEX  ',SFEIEX,SFEEEX)
          DO IFL=1,NFLA
            WRITE(iunout,'(A,I0,A,ES12.4)') 'SFNIEX(IFL=',IFL,') ',
     .                                       SFNIEX(IFL)
        ENDDO
        CALL EIRENE_MASRR1 (' TARGETS,EI',SFEIT(1),NTARGI,5)
        CALL EIRENE_MASRR1 (' TARGETS,EE',SFEET(1),NTARGI,5)
        DO ITARG=1,NTARGI
          DO IFL=1,NFLA
            WRITE(iunout,'(A,I0,A,ES12.4)') 
     .            'TARGETS, NI(IFL =',IFL,') ',
     .             SFNIT(ITARG,IFL)
          ENDDO
          CALL EIRENE_LEER(1)
        ENDDO
        CALL EIRENE_LEER(1)
        CALL EIRENE_MASR2(' TOTALS, EI,EE  ',SFEIT(0),SFEET(0))
        DO IFL=1,NFLA
           WRITE(iunout,'(A,I0,A,ES12.4)')
     .          'TOTALS, NI(IFL =',IFL,') ',
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

c  new: ion energy source terms for internal rather than total energy balance
        SCPVEII = SCPVEII / elcha
        CALL EIRENE_MASRR1 (' SSEI_INTERNAL',SCPVEII(1),NSTRAI,5)
        write (iunout,*) ' TOTAL ',sum(SCPVEII,NSTRAI)

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
           WRITE(iunout,'(A,I0,A,ES12.4)') 'BALANN(IFL=',IFL,') ',
     .                                      BALANN(IFL)
        ENDDO
        CALL EIRENE_LEER(1)

        CALL EIRENE_MASR2('REL.ERR.(%)RI,RE',RI,RE)
        DO IFL=1,NFLA
           WRITE(iunout,'(A,I0,A,ES12.4)') 'RN(IFL=',IFL,') ',RN(IFL)
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

      END SUBROUTINE

!      CONTAINS

      SUBROUTINE EIRENE_READ14_FIXED
      INTEGER :: JL

      WRITE (iunout,*) '        SUBROUTINE INFCOP IS CALLED  '
C  READ INPUT DATA OF BLOCK 14
C  SAVE INPUT DATA OF BLOCK 14 FOR SHORT CYCLE ON COMMON CCOUPL
      CALL EIRENE_LEER(1)
      CALL EIRENE_ALLOC_CCOUPL(1)
      READ (IUNIN,'(5L1)') LSYMET,LBALAN,LZDEN
      IF (TRCINT)
     .  WRITE (iunout,*) ' LSYMET,LBALAN = ',LSYMET,LBALAN
      READ (IUNIN,'(5I6)') NFLA,NCUTB,NCUTL,IMF,nfull
      NPLS_FIX = NFLA
cdr  imf  flag for different formats of geometry file: linda, sonnet, carre. What is What?
      if (imf /= 0) mshfrm = imf
      NCUTB_SAVE=NCUTB
      IF (TRCINT) THEN
        WRITE (iunout,*) ' NFLA,NCUTB,NCUTL = ',
     .                     NFLA,NCUTB,NCUTL
        WRITE (iunout,*) ' IPLS,IFLB(IPLS),FCTE(IPLS),BMASS(IPLS)'
      ENDIF
      DO IPL=1,NPLSI
        READ (IUNIN,'(2I6,2E12.4)') I,IFLB(IPL),FCTE(IPL),BMASS(IPL)
        IF (TRCINT)
     .    WRITE (iunout,*)          IPL,IFLB(IPL),FCTE(IPL),BMASS(IPL)
      END DO  ! IPL
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
      DO IT=1,NTARGI
        IF (NTGPRT(IT).GT.NGITT) THEN
          NTGPRI=NTGPRT(IT)
          CALL EIRENE_MASPRM('NGITT',5,NGITT,'NTGPRT',6,NTGPRI,IERROR)
          WRITE (iunout,*) 'EXIT CALLED FROM SUBR. INFCOP '
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
      END DO  ! IT
      IREAD=0
C  ALL INDICES: AFTER INDEX MAPPING
C  NDT: INDEX OF X-CELL (EAST OR NORTH SURFACE OF BRAAMS CELL) OF TARGET
C  NINCT: DIRECTION OF OUTER TARGET NORMAL WITH RESPECT TO POSITIVE DIR.
C  NIXY: SOURCE ON Y SURFACE: NIXY=1; SOURCE ON X SURFACE: NIXY=2
C  NTIN,NTEN: SOURCE RANGE FROM GRIDPOINT NTIN TO GRIDPOINT NTEN
      IF (TRCINT)
     .  WRITE (iunout,*) '    IT,  NDT,NINCT, NIXY, NTIN, NTEN',
     .              ',NIFLG, NPTC, NPTCM,NSPZI,NSPZE,NEMOD'
      DO IT=1,NTARGI
        DO IPRT=1,NTGPRT(IT)
          CALL EIRENE_SKIP_READ_COMMENT(IREAD,IUNIN,ZEILE)
          READ (ZEILE,'(12I6)') I,NDT(IT,IPRT),NINCT(IT,IPRT),
     .                            NIXY(IT,IPRT),NTIN(IT,IPRT),
     .                            NTEN(IT,IPRT),NIFLG(IT,IPRT),
     .                            NPTC(IT,IPRT),NPTCM(IT,IPRT),
     .                            NSPZI(IT,IPRT),NSPZE(IT,IPRT),
     .                            NEMOD(IT,IPRT)
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
     .          NTEN(IT,IPRT).GT.NR1ST) THEN
              WRITE (iunout,*) 'ERROR IN INPUT BLOCK 14, NTIN, NTEN '
              CALL EIRENE_EXIT_OWN(1)
            ENDIF
          ELSEIF (NIXY(IT,IPRT).EQ.2) THEN
            IF (NTIN(IT,IPRT).LE.0.OR.NTIN(IT,IPRT).GE.NP2ND.OR.
     .          NTEN(IT,IPRT).GT.NP2ND) THEN
              WRITE (iunout,*) 'ERROR IN INPUT BLOCK 14, NTIN, NTEN '
              CALL EIRENE_EXIT_OWN(1)
            ENDIF
          ENDIF
        END DO  ! IPRT
        IF (TRCINT) CALL EIRENE_LEER(1)
      END DO  ! IT
      READ (IUNIN,'(6E12.4)')  CHGP,CHGEE,CHGEI,CHGMOM
      IF (TRCINT) CALL EIRENE_MASR4
     .                         ('CHGP,CHGEE,CHGEI,CHGMOM         ',
     .                           CHGP,CHGEE,CHGEI,CHGMOM)
C  READ ADDITIONAL DATA TO BE TRANSFERRED FROM B2 INTO EIRENE
C  HERE: B2 VOLUME TALLIES
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
      DO IAIN=1,NAINB
        READ (IUNIN,'(6I6)') I,NAINS(IAIN),NAINT(IAIN)
        READ (IUNIN,'(A72)') TXTPLS(IAIN,12)
        READ (IUNIN,'(2A24)') TXTPSP(IAIN,12),TXTPUN(IAIN,12)
        IF (TRCINT) THEN
          WRITE (iunout,'(6I6)') I,NAINS(IAIN),NAINT(IAIN)
          WRITE (iunout,'(1X,A72)') TXTPLS(IAIN,12)
          WRITE (iunout,'(1X,2A24)') TXTPSP(IAIN,12),TXTPUN(IAIN,12)
        ENDIF
      END DO  ! IAIN
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
      DO IAOT=1,NAOTB
        READ (IUNIN,'(6I6)') I,NAOTS(IAOT),NAOTT(IAOT)
        IF (TRCINT) THEN
          WRITE (iunout,'(6I6)') I,NAOTS(IAOT),NAOTT(IAOT)
        ENDIF
      END DO  ! IAOT

C  COPY USER SPECIFIC DATA TO FILE user_data.input
      JL = 0
      IO = 0
      IUNIN_SAVE = IUNIN
      DO WHILE (IO == 0)
        READ (IUNIN,'(A72)',IOSTAT=IO) ZEILE
        IF (IO == 0) THEN
          JL = JL + 1
!pb       IF (JL == 1) OPEN(NEWUNIT=IUSROUT,FILE='user_data.input')
          IF (JL == 1) THEN
            IUSROUT = -9999
            CALL EIRENE_OPENFILE(IUSROUT,FILE='user_data.input',
     .                           FORM='FORMATTED',ACCESS='SEQUENTIAL')
          END IF
          WRITE (IUSROUT,'(A)') TRIM(ZEILE)
        END IF
      END DO
      IF (JL > 0) THEN
        REWIND IUSROUT
        IUNIN = IUSROUT
      ELSE
        IUSROUT = 0
      END IF
      
      RETURN
      END SUBROUTINE EIRENE_READ14_FIXED


      SUBROUTINE EIRENE_READ14_JSON(json,me)
      USE EIRMOD_JSON     
      use json_module           !IGNORE
     .    , lk => json_lk, rk => json_rk, ik => json_ik, ck => json_ck

      IMPLICIT NONE

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: me
      type(json_value), pointer :: pflds, pfld, ptrgs, ptrg,
     .                             prts, prt, padds, padd
      character(kind=CK,len=:),allocatable :: txt                       
      integer :: j, npl, ntrg
      integer, allocatable :: ihelp(:)
      logical :: found, foundi, foundo
      
      WRITE (iunout,*) '        SUBROUTINE INFCOP IS CALLED  '
C  READ INPUT DATA OF BLOCK 14
C  SAVE INPUT DATA OF BLOCK 14 FOR SHORT CYCLE ON COMMON CCOUPL
      CALL EIRENE_LEER(1)
      CALL EIRENE_ALLOC_CCOUPL(1)

      call json%get(me,'LSYMET',lsymet,found)
      call json%get(me,'LBALAN',lbalan,found)
      call json%get(me,'LZDEN',lzden,found)
!     call json%get(me,'LCHKQUD',lchkqud,found)
      IF (TRCINT)
     .  WRITE (iunout,*) ' LSYMET,LBALAN = ',LSYMET,LBALAN

      call json%get(me,'NFLA',nfla,found)
      call json%get(me,'NCUTB',ncutb,found)
      call json%get(me,'NCUTL',ncutl,found)
      call json%get(me,'NCUTL',ncutl,found)
      call json%get(me,'IMF',imf,found)
!     call json%get(me,'NTRFRM',ntrfrm,found)
      call json%get(me,'NFULL',nfull,found)
!     call json%get(me,'IBRAD',ibrad,found)
!     call json%get(me,'IBPOL',ibpol,found)
!     call json%get(me,'IBTOR',ibtor,found)
      NPLS_FIX = NFLA
cdr  imf  flag for different formats of geometry file: linda, sonnet, carre. What is What?
      if (imf /= 0) mshfrm = imf
      NCUTB_SAVE=NCUTB
      
      IF (TRCINT) THEN
        WRITE (iunout,*) ' NFLA,NCUTB,NCUTL = ',
     .                       NFLA,NCUTB,NCUTL
        WRITE (iunout,*) ' IPLS,IFLB(IPLS),FCTE(IPLS),BMASS(IPLS)'
      ENDIF
      
      call json%get_child(me,'B2FLUIDS',pflds)
      call json%info(pflds,n_children=npl)
      if (npl /= NPLSI) then
        write (iunout,*) ' NUMBER OF BULK IONS DOES',
     .          ' NOT MATCH NUMBER OF FLUIDS FOUND IN FILE '
        write (iunout,*) 'NPLSI = ',nplsi
        write (iunout,*) 'NPL =   ',npl
        call eirene_exit_own(1)
      end if
        
      DO IPL=1,NPLSI
        call json%get_child(pflds,ipl,pfld)

        call json%get(pfld,'IPL',i,found)
        call json%get(pfld,'IFLB',iflb(ipl),found)
        call json%get(pfld,'FCTE',fcte(ipl),found)
        call json%get(pfld,'BMASS',bmass(ipl),found)
        nullify(pfld)
        IF (TRCINT)
     .    WRITE (iunout,*)          IPL,IFLB(IPL),FCTE(IPL),BMASS(IPL)
      END DO
      nullify(pflds)
      
      call json%get(me,'NDXA',ndxa,found)
      call json%get(me,'NDYA',ndya,found)
      IF (TRCINT) WRITE (iunout,*) 'NDXA,NDYA= ',NDXA,NDYA

C     NUMBER OF TARGET SOURCES ON B2 SURFACES: NTARGI
      call json%get(me,'NTARGI',ntargi,found)
      IF (TRCINT) WRITE (iunout,*) 'NTARGI=    ',NTARGI
      CALL EIRENE_LEER(1)
      IF (NTARGI.GT.NSTEP) THEN
        CALL EIRENE_MASPRM ('NSTEP',5,NSTEP,'NTARGI',6,NTARGI,IERROR)
        WRITE (iunout,*) 'EXIT CALLED FROM SUBR. INFCOP '
        CALL EIRENE_EXIT_OWN(1)
      ENDIF

C     NUMBER OF PARTS PER TARGET SOURCE
      IF (NTARGI.GT.0) THEN
        call json%get(me,'NTGPRT',ihelp,found)
        NTGPRT(1:NTARGI) = ihelp(1:NTARGI)
        deallocate(ihelp)

        DO IT=1,NTARGI
          IF (NTGPRT(IT).GT.NGITT) THEN
            NTGPRI=NTGPRT(IT)
            CALL EIRENE_MASPRM('NGITT',5,NGITT,'NTGPRT',6,NTGPRI,IERROR)
            WRITE (iunout,*) 'EXIT CALLED FROM SUBR. INFCOP '
            CALL EIRENE_EXIT_OWN(1)
          ENDIF
        END DO
C  ALL INDICES: AFTER INDEX MAPPING
C  NDT: INDEX OF X-CELL (EAST OR NORTH SURFACE OF BRAAMS CELL) OF TARGET
C  NINCT: DIRECTION OF OUTER TARGET NORMAL WITH RESPECT TO POSITIVE DIR.
C  NIXY: SOURCE ON Y SURFACE: NIXY=1; SOURCE ON X SURFACE: NIXY=2
C  NTIN,NTEN: SOURCE RANGE FROM GRIDPOINT NTIN TO GRIDPOINT NTEN
        IF (TRCINT)
     .  WRITE (iunout,*) '    IT,  NDT,NINCT, NIXY, NTIN, NTEN',
     .              ',NIFLG, NPTC, NPTCM,NSPZI,NSPZE,NEMOD'

        call json%get_child(me,'TARGETS',ptrgs,found)
        call json%info(ptrgs,n_children=ntrg)
        if (ntrg /= NTARGI) then
          write (iunout,*) ' NUMBER OF TARGETS DOES',
     .          ' NOT MATCH NUMBER OF TARGETS FOUND IN FILE '
          write (iunout,*) 'NTARGI = ',ntargi
          write (iunout,*) 'NTRG =   ',ntrg
          call eirene_exit_own(1)
        end if
        DO IT=1,NTARGI
          call json%get_child(ptrgs,it,ptrg,found)
          call json%get_child(ptrg,'PARTS',prts,found)
          DO IPRT=1,NTGPRT(IT)
            call json%get_child(prts,iprt,prt)
            call json%get(prt,'NDT',ndt(it,iprt),found)
            call json%get(prt,'NINCT',ninct(it,iprt),found)
            call json%get(prt,'NIXY',nixy(it,iprt),found)
            call json%get(prt,'NTIN',ntin(it,iprt),found)
            call json%get(prt,'NTEN',nten(it,iprt),found)
            call json%get(prt,'NIFLG',niflg(it,iprt),found)
            call json%get(prt,'NPTC',nptc(it,iprt),found)
            call json%get(prt,'NPTCM',nptcm(it,iprt),found)
            call json%get(prt,'NSPZI',nspzi(it,iprt),found)
            call json%get(prt,'NSPZE',nspze(it,iprt),found)
            call json%get(prt,'NEMOD',nemod(it,iprt),found)
            nullify(prt)
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
     .        WRITE (iunout,'(1X,7I6,2I7,3I6)')
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
          END DO
          nullify(prts)
          nullify(ptrg)
          IF (TRCINT) CALL EIRENE_LEER(1)
        END DO
        nullify(ptrgs)
      END IF

      call json%get(me,'CHGP',chgp,found)
      call json%get(me,'CHGEE',chgee,found)
      call json%get(me,'CHGEI',chgei,found)
      call json%get(me,'CHGMOM',chgmom,found)
      IF (TRCINT) CALL EIRENE_MASR4
     .                         ('CHGP,CHGEE,CHGEI,CHGMOM         ',
     .                           CHGP,CHGEE,CHGEI,CHGMOM)

C  READ ADDITIONAL DATA TO BE TRANSFERRED FROM B2 INTO EIRENE
C  HERE: B2 VOLUME TALLIES
      call json%get(me,'NAINB',nainb,found)   
C  ADDITIONAL INPUT TALLY ADIN:  ITAL=12
      NAIN = MAX(NAIN,NAINB)
      CALL EIRENE_ALLOC_CCOUPL(2)
      WRITE (iunout,*) '        NAINI = ',NAINB
      IF (NAINB.GT.NAIN) THEN
        CALL EIRENE_MASPRM ('NAIN',4,NAIN,'NAINB',5,NAINB,IERROR)
        WRITE (iunout,*) 'EXIT CALLED FROM SUBR. INFCOP '
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
      IF (NAINB > 0) THEN
        IF (TRCINT)
     .      WRITE (iunout,*) 'I,NAINS(IAIN),NAINT(IAIN)'
        call json%get(me,'ADD_IN_TAL',padds,foundi)
        IF (FOUNDI) THEN
          DO IAIN=1,NAINB
            call json%get_child(padds,iain,padd,found)
            call json%get(padd,'NAINS',nains(iain),found)
            call json%get(padd,'NAINT',naint(iain),found)
            call json%get(padd,'TXTPLS',txt,found)
            txtpls(iain,12) = txt
            deallocate(txt)
            call json%get(padd,'TXTPSP',txt,found)
            txtpsp(iain,12) = txt
            deallocate(txt)
            call json%get(padd,'TXTPUN',txt,found)
            txtpun(iain,12) = txt
            deallocate(txt)
            nullify(padd)
            IF (TRCINT) THEN
              WRITE (iunout,'(6I6)') IAIN,NAINS(IAIN),NAINT(IAIN)
              WRITE (iunout,'(1X,A72)') TXTPLS(IAIN,12)
              WRITE (iunout,'(1X,2A24)') TXTPSP(IAIN,12),TXTPUN(IAIN,12)
            ENDIF
          END DO
          nullify(padds)
        END IF
      END IF
      
C  READ ADDITIONAL DATA TO BE TRANSFERRED FROM EIRENE INTO B2
C  HERE: EIRENE SURFACE TALLIES
      call json%get(me,'NAOTB',naotb,found)   
      WRITE (iunout,*) '        NAOTI = ',NAOTB
      IF (NAOTB.GT.NLIMPS) THEN
        CALL EIRENE_MASPRM ('NLIMPS',6,NLIMPS,'NAOTB',5,NAOTB,IERROR)
        WRITE (iunout,*) 'EXIT CALLED FROM SUBR. INFCOP '
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
      IF (NAOTB > 0) THEN
        IF (TRCINT)
     .      WRITE (iunout,*) 'I,NAOTS(IAOT),NAOTT(IAOT)'
        call json%get(me,'ADD_OUT_TAL',padds,foundo)
        IF (FOUNDO) THEN
          DO IAOT=1,NAOTB
            call json%get_child(padds,iaot,padd,found)
            call json%get(padd,'NAOTS',naots(iaot),found)
            call json%get(padd,'NAOTT',naott(iaot),found)
            nullify(padd)
            IF (TRCINT) THEN
              WRITE (iunout,'(6I6)') I,NAOTS(IAOT),NAOTT(IAOT)
            ENDIF
          END DO
          nullify(padds)
        END IF
      END IF
C
C  INPUT BLOCK 14 DONE
C
!pb   OPEN(NEWUNIT=IUSROUT,FILE='user_data.input',STATUS='OLD',
!pb  .     IOSTAT=IO)
      IUSROUT = -9999
      CALL EIRENE_OPENFILE(IUSROUT,FILE='user_data.input',STATUS='OLD',
     .     FORM='FORMATTED',ACCESS='SEQUENTIAL',IOSTAT=IO)
      IUNIN_SAVE = IUNIN
      IF (IO == 0) THEN
        IUNIN = IUSROUT
        CALL EIRENE_LEER(1)
        WRITE (IUNOUT,*) 'USER SPECIFIC INPUT READ FROM ',
     .       'user_data.input'
      ELSE
        CALL EIRENE_LEER(1)
        WRITE (IUNOUT,*) 'NO FILE FOR USR SPECIFIC INPUT FOUND'
      END IF
        
      RETURN
!      END SUBROUTINE EIRENE_READ14_JSON
C
!      END SUBROUTINE EIRENE_INFCOP
      END SUBROUTINE 

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

C> \brief Prepare some data prior to calculation of strata but after
C> the distribution of processors has been updated 
C>
      SUBROUTINE EIRENE_INFCOP_PRE_STRATA

      RETURN
      END SUBROUTINE EIRENE_INFCOP_PRE_STRATA


      SUBROUTINE EIRENE_IF3COP_SUM

      RETURN
      END SUBROUTINE EIRENE_IF3COP_SUM
 
      END MODULE EIRMOD_INFCOP
