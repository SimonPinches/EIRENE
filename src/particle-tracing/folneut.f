cdr Nov. 17   unification of update, fpath.
cdr           tbd: photon routines, static loop, logatm,mol.ion in static loop.
cdr Oct. 17   minor sync with folion
cdr           started: implementation of QSS branch: folstat_neut.f  not ready

cdr Sept.17   conditional exp. estim: external function funexp, rather than inline.
cdr           PR = prob to reach the next cell boundary.
cdr           In case of geometrical multi-steps within one macro step
cdr           (NCOU.GT.1) use PR rather than AX(2)=1, when leaving the NCOU loop
c
cdr Sept.15   Bug fix: generation limit, xgener moved in front of 100 continue
Cdr Nov.14    evaluation of NUPC(1) in static loop corrected (for 1D applications)
Cdr Oct 14 TO BE DONE: clarify role of iflag. now also used for calc-spectrum?
cdr Oct.14             spectra scoring only called if cell-based spectra are defined
c
c   ???     LDAMCEL(icell) introduced ??  "damaged cell" ??, comes via eirmod_cgeom.
c   ???     checking for 3rd grid intersection in case levgeo=4
c           (triangles plus resolution in z-direction)
C
!PB 30.01.08: optimization of calculation of intersection with additional surfaces
!             corrected
!PB 22.03.07: LEVGEO=6 --> LEVGEO=10
!PB 12.01.06: calls to update_spectrum introduced for cell-based spectra
!dr 2016:     now conditional on NADSPC_CD >= 1
!PB 02.03.06: Store trajectory from birth place to first collision with the
!             wall. It is assumed that conditional epectation estimator is
!             switched on.
!dr 2017:     conditional on ...  NLTRJ ?

!PB 18.04.06: xstorv=0 in "vacuum region" added
!PB 26.09.06: sg corrected for levgeo=4 and levgeo=5

C  MAY05: CALL UPDATE FROM STATIC LOOP WITH IFLAG=4 (RATHER =1)
C         WG. COLL EST. ON 1ST FLIGHT AFTER BIRTH.
C
      SUBROUTINE EIRENE_FOLNEUT
C
C     NEUTRAL PARTICLE, LAUNCHED AT X0,Y0,Z0, IN CELL NRCELL, IPOLG,
C     IPERID, NPCELL, NTCELL, NACELL, NBLOCK, WITH VELOCITY VELX,VELY,VELX
C     IS FOLLOWED.
C     (MODULE: COMPRT.F)
c
c
c
c
c
c
c
C
C  ON INPUT:
C     ITYP=0 OR ITYP=1 OR ITYP=2
C     IC_ION = 0  NEWBORN NEUTRAL PARTICLE, OR CONTINUATION FROM TEST ION FULL TRACK
C     IC_ION > 0  CONTINUATION FROM TEST ION WHICH WAS IN STATIC LOOP
C     IC_ION < 0  CONTINUATION FROM TEST PARTICLE IN DIFFUSION MODE (TO BE WRITTEN)
C
C  ON OUTPUT:
C
C     LGPART=TRUE
C           ITYP=3      A NEXT GENERATION TEST ION IION IS GENERATED

C     LGPART=FALSE
C           ITYP=4  NO NEXT GENERATION TEST PARTICLE IS GENERATED
C                   (PARTICLE ABSORBED IN BULK ION SPECIES)
c
c  at 100 :   start a new neutral particle, velocity is given as full cartesian vector, lcart=true
c  at 1004:   reduced (guiding centre) velocities and B-field are now set for particle. lcart=false.
C  at 1001:   particle enters static loop
C  at 1002:   particle leaves static loop
c  at 101 :   full new trajectory starts here.
c  at 104 :   an earlier track continues here.
c             initial position of track and cumulated integral for mfp sampling is not refreshed.
c             meant for continuing a track across a transparent surface
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_CADGEO
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CRAND
      USE EIRMOD_CUPD
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      USE EIRMOD_CSPEZ
      USE EIRMOD_CZT1
      USE EIRMOD_CTETRA
      USE EIRMOD_COMPRT
      USE EIRMOD_CPES
      USE EIRMOD_COMNNL
      USE EIRMOD_COMSPL
      USE EIRMOD_CLGIN
      USE EIRMOD_COUTAU
      USE EIRMOD_COMXS
      USE EIRMOD_CTRIG
      USE EIRMOD_CTRCEI
      USE EIRMOD_TIMEA, ONLY: EIRENE_TIMEA1
      USE EIRMOD_RANF, ONLY: RANF_EIRENE
      USE EIRMOD_ADDCOL, ONLY: EIRENE_ADDCOL
      USE EIRMOD_STDCOL, ONLY: EIRENE_STDCOL
      USE EIRMOD_SWITCH_PARTINFO, ONLY: EIRENE_SWITCH_PARTINFO
      USE EIRMOD_PLT2D, ONLY: EIRENE_CHCTRC
      use eirmod_timer
      use eirmod_timep
      use eirmod_colatm  
      use eirmod_colmol

      use OMP_LIB

      IMPLICIT NONE

      REAL(DP) :: CFLAG(7,MSTOR0)
      REAL(DP) :: AX(2)
      REAL(DP) :: XSTOR2(MSTOR1,MSTOR2,N2ND+N3RD),
     .            XSTORV2(NSTORV,N2ND+N3RD)
      REAL(DP) :: XSTORC(MSTOR1,MSTOR2), XSTORVC(NSTORV)
      REAL(DP) :: VELXC, TIMEC, WS, COLTYP, X0C, Y0C, Z0C, ZDT1C,
     .          X0ERR, Y0ERR, Z0ERR, VELC, E0C, VELYC, VELZC, SG,
     .          GENRC, PHIC, WEIGHC, ZLI, XLI, YLI, T, ZTS,
     .          ZMFP, ZEP1, ZLOG, ZTST, ZINT1, ZINT2, Z0S, TIMES,
     .          X0S, Y0S, PHIS, DIST, ZTC, PSAVE, TSAVE,
     .          EX, EXPM, FF, WMINC_LOCAL, PR, PPR,   ! cond exp. est
     .          EIRENE_FPATH,
     .          SCOS_NEW
ctk      REAL(DP), EXTERNAL :: RANF_EIRENE, EIRENE_FUNEXP
      REAL(DP), EXTERNAL :: EIRENE_FUNEXP
      INTEGER :: NBLCKC, NCELLC, NRCLLC, NACLLC, ITIMEC, IPERIDC,
     .           IFPTHC, IUPDTC, NPCLLC, NTCLLC, NTSAVE, NPSAVE,
     .           EIRENE_LEARC2, J, NCOUS, NLE, NRC, JCOL, NLI, ISTS,
     .           NPCOLC,
     .           JJ, NPCELC, NTCELC, NTCOLC, IFLAG, I, IM,
     .           NCLLN, IRET, IRT_STAT
      LOGICAL :: NLPR, LCNDEXP
      TYPE(CELL_INFO), POINTER :: NEW_CELL


c  all cell indices must be known at this point
c  tentatively assume: a next generation particle will be born

c  IC_NEUT, IC_ION: counter for generations within static loop
      IC_NEUT=IC_ION
C  XGENER: COUNTER FOR GENERATION LIMIT
      XGENER=0.D0

  100 LGPART=.TRUE.
      IC_NEUT=IC_NEUT+1

C  INITIALIZE COND. EXP. ESTIMATOR
      NLPR=.FALSE.
      AX(1)=1.
      AX(2)=1.
      PR=AX(2)
      WMINC_LOCAL=WMINC

C  CHECK FOR VALID SPECIES INDEX
      IF (ITYP.EQ.0) THEN
        IF (IPHOT.LE.0.OR.IPHOT.GT.NPHOTI) GOTO 998
      ELSEIF (ITYP.EQ.1) THEN
        IF (IATM.LE.0.OR.IATM.GT.NATMI) GOTO 998
      ELSEIF (ITYP.EQ.2) THEN
        IF (IMOL.LE.0.OR.IMOL.GT.NMOLI) GOTO 998
      ENDIF

      IF (ITYP.EQ.1) THEN
        LOGATM(IATM,ISTRA)=.TRUE.
        NLPR=NLPRCA(IATM)
        NRC=NRCA(IATM)
      ELSEIF (ITYP.EQ.2) THEN
        LOGMOL(IMOL,ISTRA)=.TRUE.
        NLPR=NLPRCM(IMOL)
        NRC=NRCM(IMOL)
      ELSEIF (ITYP.EQ.0) then
        LOGPHOT(IPHOT,ISTRA)=.TRUE.
        NLPR=NLPRCPH(IPHOT)
        NRC=NRCPH(IPHOT)
      ENDIF
C
C  THE CELL NUMBER NRCELL, IPOLG, IPERID, NPCELL, NTCELL, NACELL, NBLOCK
C  WAS ALREADY SET IN CALLING SUBROUTINE MCARLO
C
C  IF NLSRFX, SURFACE INDEX MRSURF MUST BE DEFINED AT THIS POINT
C  IF NLSRFY, SURFACE INDEX MPSURF MUST BE DEFINED AT THIS POINT
C  IF NLSRFZ, SURFACE INDEX MTSURF MUST BE DEFINED AT THIS POINT
C  IF NLSRFA, SURFACE INDEX MASURF MUST BE DEFINED AT THIS POINT
C
C  FOLLOW MOTION OF NEUTRAL PARTICLE OR "STATIC APPROXIMATION"?
      IF (IFPATH.NE.1.OR.NRC.LT.0) GOTO 1002
C  STOP STATIC LOOP AFTER 100 GENERATIONS LATEST TO AVOID ACCIDENTAL INFINITE LOOPS
      IF (IC_NEUT.GT.100) GOTO 1002

      IF (ITYP.EQ.1) THEN
        IF (NFOLA(IATM).NE.-1) GOTO 1002 ! go to static loop
      ELSEIF (ITYP.EQ.2) THEN
        IF (NFOLM(IMOL).NE.-1) GOTO 1002
      ELSEIF (ITYP.EQ.0) THEN
        IF (NFOLPH(IPHOT).NE.-1) GOTO 1002
      ENDIF
C
C***********************************************************************
C  STATIC APPROXIMATION
C  SIMULATE NEXT COLLISION INSTANTANEOUSLY
C***********************************************************************

      CALL EIRENE_FOLSTAT_NEUT(IC_NEUT,VELX,VELY,VELZ,CFLAG,
     .                         IRT_STAT)
      IF (IRT_STAT == 1) GOTO 101
      IF (IRT_STAT == 2) GOTO 230
      IF (IRT_STAT == 3) GOTO 380
C.............................................................

C
 1002 CONTINUE

C  AT THIS POINT: EITHER CONTINUE A FULL TRAJECTORY
C                 OR PARTICLE WAS IN STATIC APPROXIMATION,
C                 BUT NOW IT RETURNS TO FULL MOTION
C
      IF (IC_NEUT.GT.1.AND.NLTRC.AND.TRCHST)
     .  WRITE (iunout,*) 'TRAJECTORY LEAVES STATIC LOOP, ITYP=',ITYP

C  IN CASE THAT THE PARTICLE WAS IN STATIC LOOP AND ON A SURFACE,
C  SOME MORE WORK NEEDS TO BE DONE, TO REVIVE IT TO FULL KINETIC ORBIT MODE.
      IF (IC_NEUT.GT.1.AND.
     .   (NLSRFX.OR.NLSRFY.OR.NLSRFZ.OR.NLSRFA)) THEN

C  PARTICLE CONTINUES FROM SURFACE AND FROM PREVIOUS "STATIC LOOP"
C  PREPARE CELL NUMBERS FOR FIRST FLIGHT
        IC_ION=0
        IC_NEUT=0
        SCOS_NEW = SIGN(1.D0,VELX*CRTXG+VELY*CRTYG+VELZ*CRTZG)
        IF (SCOS_SAVE.NE.SCOS_NEW) THEN
          SCOS=SCOS_NEW
          ZT=0.D0
          TL=0.D0
          IPOLGN=IPOLG
          IF (NLSRFA) THEN
            CALL EIRENE_ADDCOL (X0,Y0,Z0,SCOS,IRET)
            if (IRET .EQ. 1) GOTO 101
            if (IRET .EQ. 2) GOTO 380
          ELSEIF (NLSRFX) THEN
            select case (LEVGEO)
            case (:3)
              ISTS=INMP1I(MRSURF,IPCELL,ITCELL)
              MSURFG=NPCELL+(NTCELL-1)*NP2T3
              IF (ILIIN(NLIM+ISTS) .NE. 0) THEN
                 CALL EIRENE_STDCOL (ISTS,1,SCOS,IRET)
                 if (IRET .EQ. 1) GOTO 101
                 if (IRET .EQ. 2) GOTO 380
              ENDIF
            case (4)
              ISTS=ABS(INMTI(IPOLGN,MRSURF))  !dr NLIM already added in ISTS ?
              MSURFG=INSPAT(IPOLGN,MRSURF)
              IF (ILIIN(ISTS) .NE. 0) THEN
                 CALL EIRENE_STDCOL (ISTS,1,SCOS,IRET)
                 if (IRET .EQ. 1) GOTO 101
                 if (IRET .EQ. 2) GOTO 380
              ENDIF
            case (5)
              ISTS=ABS(INMTIT(IPOLGN,MRSURF)) !dr NLIM already added in ISTS ?
C             MSURFG= ??
              IF (ILIIN(ISTS) .NE. 0) THEN
                 CALL EIRENE_STDCOL (ISTS,1,SCOS,IRET)
                 if (IRET .EQ. 1) GOTO 101
                 if (IRET .EQ. 2) GOTO 380
              ENDIF
            case (10)
              ISTS=INMP1I(MRSURF,IPCELL,ITCELL)
C             MSURFG= ??
              IF (ILIIN(NLIM+ISTS) .NE. 0) THEN
                 CALL EIRENE_STDCOL (ISTS,1,SCOS,IRET)
                 if (IRET .EQ. 1) GOTO 101
                 if (IRET .EQ. 2) GOTO 380
              ENDIF
            end select
          ELSEIF (NLSRFY) THEN
            ISTS=INMP2I(IRCELL,MPSURF,ITCELL)
            MSURFG=NRCELL+(NTCELL-1)*NR1P2
            IF (ILIIN(NLIM+ISTS) .NE. 0) THEN
               CALL EIRENE_STDCOL (ISTS,2,SCOS,IRET)
               if (IRET .EQ. 1) GOTO 101
               if (IRET .EQ. 2) GOTO 380
            ENDIF
          ELSEIF (NLSRFZ) THEN
            ISTS=INMP3I(IRCELL,IPCELL,MTSURF)
            MSURFG=NRCELL+(NPCELL-1)*NR1P2
            IF (ILIIN(NLIM+ISTS) .NE. 0) THEN
               CALL EIRENE_STDCOL (ISTS,3,SG,IRET)
               if (IRET .EQ. 1) GOTO 101
               if (IRET .EQ. 2) GOTO 380
            ENDIF
          ENDIF
C         WRITE (IUNOUT,*) 'FOLNEUT: I SHOULD NOT BE HERE'
        ENDIF
      ENDIF

C**********************************************************************
C   STATIC LOOP FINISHED. REGULAR PARTICLE TRACKING CONTINUES
C**********************************************************************

      IC_ION=0
      IC_NEUT=0
C
C  PARTICLE IN VOLUME OR ON SURFACE BUT NOT FROM "STATIC LOOP"
C
C  EACH NEUTRAL PARTICLE TRACK STARTS AT THIS POINT, IC_NEUT=0 HERE
C
  101 CONTINUE
      IF (ITYP.EQ.1) THEN
        LOGATM(IATM,ISTRA)=.TRUE.
        NLPR=NLPRCA(IATM)
        NRC=NRCA(IATM)
      ELSEIF (ITYP.EQ.2) THEN
        LOGMOL(IMOL,ISTRA)=.TRUE.
        NLPR=NLPRCM(IMOL)
        NRC=NRCM(IMOL)
      ELSEIF (ITYP.EQ.0) then
        LOGPHOT(IPHOT,ISTRA)=.TRUE.
        NLPR=NLPRCPH(IPHOT)
        NRC=NRCPH(IPHOT)
      ENDIF
      IF (NLTRJ .AND. .NOT.NLPR) THEN
        WRITE (IUNOUT,*) ' STORING OF TRAJECTORIES SWITCHED OFF',
     .           ' BECAUSE NO CONDITIONAL EXPECTATION ESTIMATOR'
        NLTRJ = .FALSE.
      ENDIF
C  WEIGHT TOO SMALL? STOP HISTORY
      IF (WEIGHT.LT.EPS30) THEN
        LGPART=.FALSE.
        RETURN
      ENDIF
      ICOL=0
      JCOL=0
C
      ZEP1=RANF_EIRENE( )
      ZLOG=-LOG(ZEP1)
      ZINT1=0.0
      ZINT2=ZINT1

C  COORDINATES FOR SUB-STEPS IN 2ND AND/OR 3RD GRID, (ONLY IF NEEDED: NL2ND; NL3RD, NLTRA)
      IF (NLTRA) X01=X0+RMTOR
      X00=X0
      Y00=Y0
      Z00=Z0
      Z01=Z0
C
C  CLEAR WORK VARIABLES AND: CONTINUE FLIGHTS ACROSS TRANSPARENT
C                            SURFACES FROM THIS POINT: NEW POSITION; OLD VELOCITY
C                            REFRESH MFP SAMPLING
  104 CONTINUE
      NJUMP=0
      DO I=1,NIMINT
        IM=IIMINT(I)
        TIMINT(IM)=0._DP
        IIMINT(I)=0
      END DO
      NIMINT = 0
      PHIS=0
      TT=1.D30
      TL=1.D30
      TS=1.D30
      ZTST=1.D30
      ZT=0.0
C
  110 CONTINUE

      NCOU=1
      NUPC(1)=0
      NCOUNT(1)=1
      NCOUNP(1)=1
      ISRFCL=-1

      NCELL=NRCELL+((NPCELL-1)+(NTCELL-1)*NP2T3)*NR1P2+NBLCKA
      IF (LDAMCEL(NCELL)) GOTO 9912
C
C TL: DISTANCE TO NEXT ADDITIONAL SURFACE
c   nli,nle: index range of additional surfaces, visible from cell no. ncell
      IF (NCELL.GT.0.AND.NCELL.LE.NOPTIM) THEN
        NLI=NLIMII(NCELL)
        NLE=NLIMIE(NCELL)
      ELSEIF (NCELL.GT.0) THEN
        NLI=1
        NLE=NLIMI
      ELSE
C  NEGATIVE CELL INDEX. STOP THIS PARTICLE
        GOTO 990
      ENDIF
      IF (NLI.LE.NLE) THEN
c  check all additional surfaces in index range nli, nle
        CALL EIRENE_TIMEA1
     .  (MSURF,NCELL,NLI,NLE,NTCELL,IPERID,X0,Y0,Z0,TIME,
     .               VELX,VELY,VELZ,VEL,
     .               MASURF,XLI,YLI,ZLI,SG,TL,NLTRC,LCNDEXP)
c  activate conditional expectation estimator,
c  if history points towards one of selected add. surfaces
!pb     NLPR=NLPRCS(MASURF).OR.NLPR
        NLPR=LCNDEXP.OR.NLPR
c
        ZDT1=TL
        ZTST=TL
        CLPD(1)=ZDT1
        IF (MASURF.NE.0) ISRFCL=1
      ENDIF
C
C TT: DISTANCE UNTIL NEXT TIMESTEP LIMIT IS REACHED
      IF (LGTIME) THEN
        TT=(DTIMVI-TIME)*VEL
        IF (TT.LT.TL) THEN
          ZDT1=TT
          ZTST=TT
          CLPD(1)=ZDT1
          ISRFCL=2
        ENDIF
      ENDIF
C
C  SCAN OVER SEGMENT
C
  210 CONTINUE
C
C  TS:   DISTANCE TO NEXT SURFACE OF STANDARD MESH
C  ZDT1: DISTANCE TRAVELLED IN CURRENT CELL
C  ZT: DISTANCE ALREADY TRAVELLED IN PREVIOUS PARTS OF THIS TRACK
C
      IF (ITIME.EQ.1) THEN
        IF (NLRAD) THEN
          CALL EIRENE_TIMER(TS)
          IF (.NOT.LGPART) GOTO 9911
C
          T=TS/TL-1.0_DP
          IF (ABS(T).LE.EPS10.AND.TL.NE.1.E30_DP) GOTO 992
          IF (TL.LT.TS.OR.TT.LT.TS) THEN
            MRSURF=0
            IPOLGN=0
C  CHECK FOR INTERSECTION WITH ADDITIONAL SURFACE
            IF (TL.LE.TT) THEN
              ZDT1=TL-ZT
              TL=ZT+ZDT1
              ZTST=TL
              ISRFCL=1
C  INTERSECTION WITH TIME SURFACE. TIME LIMIT REACHED ?
            ELSEIF (TT.LT.TL) THEN
              ZDT1=TT-ZT
              TT=ZT+ZDT1
              ZTST=TT
              ISRFCL=2
            ENDIF
          ELSE
C  INTERSECTION A  WITH 1-ST (RADIAL) GRID SURFACE
            ZDT1=TS-ZT
            ZTST=TS
            ISRFCL=0
          ENDIF
        ENDIF
C
        NCOU=1
        NUPC(1)=0
        CLPD(1)=ZDT1
        NCOUNT(1)=1
        NCOUNP(1)=1

C  CHECK SUB-GRIDS.  FOR OPTIONAL Y,Z RESOLUTION, ON 1D BACKGROUND MEDIUM
c  (can be switched off: nlpol, nltor, nltra)
C  sub-cells have the same background parameters as the parent (x-or-radial) cell,
C  i.e. same collision rates, same mean free path.

C  3RD Z (OR TOROIDAL) SUB-GRID, ALSO:  TOROIDAL PERIODICITY SURFACES
C  SUBDIVIDE GIVEN TRACK INTO Z (OR TOROIDAL) SMALLER SEGMENTS
        IF (NLTOR.OR.NLTRA) THEN
          CALL EIRENE_TIMET (ZDT1)
          TS=ZT+ZDT1
          ZTST=TS
        ENDIF
C  2ND (OR POLOIDAL) SUB-GRID
        IF (NLPOL) THEN
          CALL EIRENE_TIMEP(ZDT1)
          TS=ZT+ZDT1
          ZTST=TS
        ENDIF
C
        IF (ZDT1.LE.0.D0) GOTO 990
C
      ELSEIF (ITIME.NE.1) THEN
C
        IF (NLTOR.OR.NLTRA) THEN
          CALL EIRENE_TIMET (ZDT1)
          TS=ZT+ZDT1
          ZTST=TS
        ENDIF

      ENDIF
C
      IF (ZTST.GE.1.D30) GOTO 990
C
C  LOCAL MEAN FREE PATH
C
C  ONE "RADIAL" (FIRST GRID) AND IN TOTAL
C  NCOU CELLS ARE CROSSED BY THE CURRENT TRACK.
C  EVALUATE REACTION RATES, MEAN FREE PATH, ETC. IN THESE CELLS
C
      IFLAG=3

      IF (NLTRJ) THEN
C  STORE THIS TRAJECTORY, FOR LATER USE IN CORRELATED SAMPLING
        TRAJ(ITRJ)%TRJ%NCOU_CELL = TRAJ(ITRJ)%TRJ%NCOU_CELL + NCOU
        DO J=1,NCOU
          NCELL=NRCELL+NUPC(J)*NR1P2+NBLCKA
          ALLOCATE(NEW_CELL)
          NEW_CELL%NO_CELL = NCELL
          NEW_CELL%FLIGHT = CLPD(J)
          CALL EIRENE_CELL_INSERT(ITRJ,NEW_CELL)
        END DO
      END IF

      IF (IFPATH.NE.1.OR.NRC.LT.0) THEN
C  USE VACUUM VALUES FOR REACTION RATES, MFP, ETC..
        XSTORV(:) =0.D0
        IF (NCOU.GT.1) THEN
          XSTOR2(:,:,1:NCOU)=0.D0
          XSTORV2(:,1:NCOU) =0.D0
        ENDIF
        ZMFP=1.D10
        NCELL=NRCELL+((NPCELL-1)+(NTCELL-1)*NP2T3)*NR1P2+NBLCKA
        IF (LDAMCEL(NCELL)) GOTO 9912
      ELSE
        NCOUS=NCOU
        ZTS=ZT
        DO 212 J=1,NCOU
          JJ=J
          NCELL=NRCELL+NUPC(J)*NR1P2+NBLCKA
          IF (LDAMCEL(NCELL)) GOTO 9912
          ZMFP=EIRENE_FPATH(NCELL,CFLAG,J,NCOU)

c  so far for photons only: local (WMINL) criterion for cond. exp.est.
c  if mfp smaller than geometrical step size times WMINL, turn off
c  cond. exp. est.
          IF ((ITYP.EQ.0).AND.(ZMFP < WMINL*CLPD(J))) WMINC_LOCAL=1._DP

          IF (NCOU.GT.1) THEN
            XSTOR2(:,:,J) = XSTOR(:,:)
            XSTORV2(:,J) = XSTORV(:)
          ENDIF

C  UPDATE INTEGRAL
          ZINT1=ZINT1+CLPD(J)*ZMFPI
          IF (.NOT.NLPR) THEN
            IF (ZINT1.GE.ZLOG) THEN
C  COLLISION IN SECTION J OF CURRENT TRACK
              IF (NLPOL) NPCELL=NCOUNP(J)
              IF (NLTOR) NTCELL=NCOUNT(J)
              GO TO 213
            ENDIF
            ZINT2=ZINT1
            ZT=ZT+CLPD(J)
          ELSEIF (NLPR) THEN
            IF (JCOL.EQ.0) THEN
C   NO COLLISION YET ON THIS SEGMENT
              IF (ZINT1.GE.ZLOG) THEN
C   NOW FIRST COLLISION FOUND; IN SUB-SECTION NO. J
                JCOL=J
                IF (NLPOL) NPCOLC=NCOUNP(J)
                IF (NLTOR) NTCOLC=NCOUNT(J)
              ELSE
c   STILL UNCOLLIDED FLUX
                ZINT2=ZINT1
                ZT=ZT+CLPD(J)
              ENDIF
            ENDIF
c
c  conditional expexctation estimator for flight segment J
            AX(1)=AX(2)
            EX=CLPD(J)*ZMFPI
c
c  find new AX(1)= AX(1) * (1-exp(-ex))/ex)
c           AX(1)= AX(1) * funexp(ex)
c  funexp(x), with 0 <= x <= infinity
            FF=EIRENE_FUNEXP(EX,EXPM)
            AX(1)=AX(1)*FF

cdr old: inline-version
c           IF (EX.LE.1.D-10) THEN
c             EXPM=1.
Cc            AX(1)=AX(1)
c           ELSEIF (EX.GT.1.D2) THEN
c             EXPM=0.D0
c             AX(1)=AX(1)/EX
c           ELSE
c             EXPM=EXP(-EX)
c             AX(1)=AX(1)*(1.-EXPM)/EX
c           ENDIF
c  done: AX(1) is adapted, and EXPM is set.
c
            ZTS=ZTS+CLPD(J)
            IF (NLPOL) NPCELC=NCOUNP(J)
            IF (NLTOR) NTCELC=NCOUNT(J)
            CLPD(J)=CLPD(J)*AX(1)
C  PROB. FOR REACHING NEXT CELL BOUNDARY
            AX(2)=AX(2)*EXPM
            PR=AX(2)
C  COND. EXP.EST: STOP BECAUSE OF WMINC CRITERION
            IF (.NOT.NLTRJ.AND.(AX(2).LE.WMINC_LOCAL)) THEN
C    RESTORE POINT OF COLLISION ?
              IF (JCOL.NE.0) GOTO 213
C  JCOL=0: NO COLLISION YET; CONTINUE LOOP 212
C          REFRESH COND. EXP. EST. FOR NEXT SECTION JJ=J+1.
C          but what if this was the last section of this flight, J=NCOU
C          we might then need old PR for surface tallies.
              AX(1)=1.
              AX(2)=1.
            ENDIF
          ENDIF
  212   CONTINUE   ! NCOU LOOP OVER SUB-STEPS ICOU IN BIG RADIAL STEP: DONE
C
  213   CONTINUE   ! EXIT FROM NCOU LOOP DUE TO COLLISION AT JCOL=JJ
        NCOU=JJ
      ENDIF
C
C  CHECK FOR EVENT
C
      IF (NLPR) THEN
C  CHECK FOR 1.ST COLLISION ALONG TRACK
        IF (ICOL.EQ.0.AND.ZINT1.GE.ZLOG) GOTO 505
C  STOP CONDITIONAL TRACK AT EARLIER SECTION, BECAUSE JCOL AND WMINC CRITERION?
        IF (NCOU.LT.NCOUS) THEN
          IFLAG=2
          IF (IUPDTE.GE.1) THEN
            CALL EIRENE_UPDATE(XSTOR2,XSTORV2,IFLAG)
            IF (NADSPC_CD >= 1)
     .        CALL EIRENE_UPDATE_SPECTRUM (WEIGHT,IFLAG,1)
          ENDIF
          ZT=ZTS
          GOTO 216
        ENDIF
C  STOP TRACK BECAUSE OF COLLISION
      ELSEIF (ZINT1.GE.ZLOG) THEN
        GO TO 220
      ENDIF
C
  215 CONTINUE
C
      ZINT2=ZINT1

C  SET NEW ACCUMULATED FLIGHT LENGTH, TENTATIVE
      ZT=ZTST
C
C  UPDATE CONTRIBUTION TO VOLUME-AVERAGED ESTIMATORS
C
      IF (IUPDTE.GE.1) THEN
cdr     IFLAG= ?
        CALL EIRENE_UPDATE(XSTOR2,XSTORV2,IFLAG)
        IF (NADSPC_CD >= 1) CALL EIRENE_UPDATE_SPECTRUM (WEIGHT,IFLAG,1)
      ENDIF
C
C  STOP TRACK ?
C
      IF (ISRFCL.EQ.1) THEN
        CALL EIRENE_ADDCOL (XLI,YLI,ZLI,SG,IRET)
        IF (IRET .EQ. 1) GOTO 104
        IF (IRET .EQ. 2) GOTO 380
      ENDIF
!pb   IF (ISRFCL.EQ.2) CALL EIRENE_TIMCOL (PR,            *104,*800)
      IF (ISRFCL.EQ.2)THEN
        CALL EIRENE_TIMCOL (PR,IRET)
        IF (IRET .EQ. 1) GOTO 104
        IF (IRET .EQ. 2) GOTO 800
      ENDIF
!pb   IF (ISRFCL.EQ.3) CALL EIRENE_TORCOL (               *104)
      IF (ISRFCL.EQ.3) THEN
        CALL EIRENE_TORCOL (IRET)
        IF (IRET .EQ. 1) GOTO 104
      ENDIF
C
C  NO, CONTINUE TRACK
C
C  NEXT CELL - CHECK FOR ESCAPE OR NON-DEFAULT ACTING STANDARD SURFACE
C
      select case (LEVGEO)
      case (:3)
C  ESCAPE AT 1ST GRID SURFACE (X OR RADIAL) MRSURF
        ISTS=INMP1I(MRSURF,IPCELL,ITCELL)
        IF (NLRAD.AND.ISTS.NE.0) THEN
          SG=ISIGN(1,NINCX)
          NLSRFX=.TRUE.
          MSURFG=NPCELL+(NTCELL-1)*NP2T3
          IF (ILIIN(NLIM+ISTS) .NE. 0) THEN
             CALL EIRENE_STDCOL (ISTS,1,SG,IRET)
             if (IRET .EQ. 1) GOTO 104
             if (IRET .EQ. 2) GOTO 380
          ENDIF
        ENDIF

C  ESCAPE AT 2ND GRID SURFACE (Y OR POLOIDAL) NO. MPSURF
        ISTS=INMP2I(IRCELL,MPSURF,ITCELL)
        IF (NLPOL.AND.ISTS.NE.0) THEN
          SG=ISIGN(1,NINCY)
          NLSRFY=.TRUE.
          MSURFG=NRCELL+(NTCELL-1)*NR1P2
          IF (ILIIN(NLIM+ISTS) .NE. 0) THEN
             CALL EIRENE_STDCOL (ISTS,2,SG,IRET)
             if (IRET .EQ. 1) GOTO 104
             if (IRET .EQ. 2) GOTO 380
          ENDIF
        ENDIF

C  ESCAPE AT 3RD GRID SURFACE (Z OR TOROIDAL) MTSURF
        ISTS=INMP3I(IRCELL,IPCELL,MTSURF)
        IF (NLTOR.AND.ISTS.NE.0) THEN
          SG=ISIGN(1,NINCZ)
          NLSRFZ=.TRUE.
          MSURFG=NRCELL+(NPCELL-1)*NR1P2
          IF (ILIIN(NLIM+ISTS) .NE. 0) THEN
             CALL EIRENE_STDCOL (ISTS,3,SG,IRET)
             if (IRET .EQ. 1) GOTO 104
             if (IRET .EQ. 2) GOTO 380
          ENDIF
        ENDIF
C
C  ESCAPE AT GRID SURFACE BUILT FROM TRIANGLE SIDES IN X-Y PLANE: MRSURF
      case (4)
        ISTS=ABS(INMTI(IPOLGN,MRSURF)) !dr NLIM already added in ISTS ?
        IF (NLRAD.AND.ISTS.NE.0) THEN
!pb          SG=ISIGN(1,NINCX)
          SG=SIGN(1._DP,VELX*PTRIX(IPOLGN,MRSURF)+
     .                  VELY*PTRIY(IPOLGN,MRSURF))
          NLSRFX=.TRUE.
          MSURFG=INSPAT(IPOLGN,MRSURF)
          IF (ILIIN(ISTS) .NE. 0) THEN
             CALL EIRENE_STDCOL (ISTS,1,SG,IRET)
             if (IRET .EQ. 1) GOTO 104
             if (IRET .EQ. 2) GOTO 381
          ENDIF
  381     CONTINUE
          SG=INMTINSS(IPOLGN,MRSURF)                            !VK
          GOTO 380                                              !VK
        ENDIF

C  ESCAPE AT 3RD (Z OR TOROIDAL) GRID SURFACE FOR TRIANGULAR X-Y GRID OPTION: MTSURF
        IF (MTSURF > 0) THEN
          ISTS=INMTI3(IRCELL,MTSURF)
          IF (NLTOR.AND.ISTS.NE.0) THEN
            SG=ISIGN(1,NINCZ)
            NLSRFZ=.TRUE.
            MSURFG=NRCELL+(NPCELL-1)*NR1P2
            IF (ILIIN(NLIM+ISTS) .NE. 0) THEN
               CALL EIRENE_STDCOL(ISTS,3,SG,IRET)
             if (IRET .EQ. 1) GOTO 104
             if (IRET .EQ. 2) GOTO 380
            ENDIF
          ENDIF
        END IF
C
C  ESCAPE AT GRID SURFACE BUILD FROM TETRAHEDRA SIDES: MRSURF
      case (5)
        ISTS=ABS(INMTIT(IPOLGN,MRSURF))  !dr NLIM already added in ISTS ?
        IF (NLRAD.AND.ISTS.NE.0) THEN
          SG=SIGN(1._DP,VELX*PTETX(IPOLGN,MRSURF)+
     .                  VELY*PTETY(IPOLGN,MRSURF)+
     .                  VELZ*PTETZ(IPOLGN,MRSURF))
          NLSRFX=.TRUE.
C         MSURFG= ??
          IF (ILIIN(ISTS) .NE. 0) THEN
             CALL EIRENE_STDCOL (ISTS,1,SG,IRET)
             if (IRET .EQ. 1) GOTO 104
             if (IRET .EQ. 2) GOTO 380
          ENDIF
       ENDIF

C  ESCAPE TO GRID SURFACE ON USER-DEFINED GEOMETRY BLOCK: MRSURF
      case (10)
        ISTS=INMP1I(MRSURF,IPCELL,ITCELL)
        IF (NLRAD.AND.ISTS.NE.0) THEN
          SG=ISIGN(1,NINCX)
          NLSRFX=.TRUE.
          IF (ILIIN(NLIM+ISTS) .NE. 0) THEN
             CALL EIRENE_STDCOL (ISTS,1,SG,IRET)
             if (IRET .EQ. 1) GOTO 104
             if (IRET .EQ. 2) GOTO 380
          ENDIF
        ENDIF
      end select
C
      NRCELL=NRCELL+NINCX
      IF (NRCELL.GT.NR1STM) GOTO 990
      IF (NACELL.LT.1.AND.NRCELL.LT.1) GOTO 990
!pb
      if (.not.lgpart) goto 9912
C
C  PARTICLE ON SURFACE MRSURF BELONGING TO 1ST (RADIAL OR X-) GRID
C
C  IF NOT, THEN IT MUST, FOR SOME REASON,
C  HAVE BEEN STOPPED IN THE MIDDLE OF A TRACK.
      IF (MRSURF.EQ.0) THEN
C  ADVANCE IN SAME CELL, AND CONTINUE TRACK
        X0=X0+VELX*ZT
        Y0=Y0+VELY*ZT
        Z0=Z0+VELZ*ZT
        TIME=TIME+ZT/VEL
        IPOLG=IPOLGN
        MASURF=0
        MSURF=0
        IF (NLTRA) THEN
          PHI=MOD(PHI-ATAN2(Z01,X01)+ATAN2(Z0,(RMTOR+X0)),PI2A)
          X01=X0+RMTOR
        ENDIF
        X00=X0
        Y00=Y0
        Z00=Z0
        Z01=Z0
        GOTO 104
      ENDIF
C
C  CHECK IF WE HAVE ENCOUNTERED A SPLITTING ZONE
C  SPLITTING AND RR NOT READY FOR LEVGEO.GE.4
      IF (LEVGEO.LE.3) THEN
        IF (NLSPLT(MRSURF).AND.NLEVEL.LT.MAXLEV.AND.ICOL.EQ.0) THEN
!PB       CALL EIRENE_SPLTRR(1,MRSURF,NINCX,*210,*700)
          CALL EIRENE_SPLTRR(1,MRSURF,NINCX,IRET)
          IF (IRET == 1) GOTO 210
          IF (IRET == 2) GOTO 700
        ENDIF
      ENDIF
C
  216 CONTINUE
C
C  SWITCH OFF CONDITIONAL EXP. ESTIMATOR
C  AT AN INTERNAL SURFACE ?
      IF (NLPR.AND..NOT.NLTRJ.AND.(AX(2).LT.WMINC_LOCAL)) THEN
        IF (NLTRC) THEN
C  TEMPORARILY
          X0S=X0+VELX*ZT
          Y0S=Y0+VELY*ZT
          Z0S=Z0+VELZ*ZT
          TIMES=TIME+ZT/VEL
          IF (NLTRA)
     .      PHIS=MOD(PHI-ATAN2(Z01,X01)+ATAN2(Z0S,(X0S+RMTOR)),PI2A)
          PSAVE=PHI
          PHI=PHIS
          TSAVE=TIME
          TIME=TIMES
cym
!$OMP CRITICAL
          CALL EIRENE_CHCTRC(X0S,Y0S,Z0S,16,19)
!$OMP END CRITICAL
          PHI=PSAVE
          TIME=TSAVE
        ENDIF
        IF (ICOL.EQ.1) GOTO 512
C  ICOL=0: NO COLLISION YET; RESTART AGAIN WITH FRESH COND. EXP. ESTIMATOR
C                            IN NEXT CELL
c       IF (NLTRC)
c    .     write (iunout,*) 'continue without restart, npanu ',
c    .                       npanu,ICOL
        AX(1)=1.
        AX(2)=1.
        JCOL=0
      ENDIF
C  EITHER: GOTO 101, NEW RANDOM NUMBER, ZINT1=0, X=XS
C  OR    : GOTO 210, CONTINUE TRACK,
C  THIS IS THE SAME, BECAUSE OF EXPONENTIAL DISTRIBUTION OF PATH LENGTHS
      IF (NCELL.LE.NOPTIM) THEN
C  FROM NEW CELL NCLLN MORE ADDITIONAL SURFACES MIGHT BE VISIBLE BY THE PARTICLES
        NCLLN=NRCELL+((NPCELL-1)+(NTCELL-1)*NP2T3)*NR1P2+NBLCKA
        IF (NCLLN <= NOPTIM) THEN
          IF ((NLIMII(NCLLN) < NLIMII(NCELL)) .OR.
     .        (NLIMIE(NCLLN) > NLIMIE(NCELL))) GOTO 110
        END IF
      END IF
      GOTO 210
C
C  POINT OF COLLISION  220 -- 240
C
  220 CONTINUE
C
      DIST=CLPD(NCOU)
      CLPD(NCOU)=(ZLOG-ZINT2)*ZMFP
      ZTC=ZT+CLPD(NCOU)
      IFLAG=4
      IF (IUPDTE.GE.1) THEN
        CALL EIRENE_UPDATE(XSTOR2,XSTORV2,IFLAG)
        IF (NADSPC_CD >= 1) CALL EIRENE_UPDATE_SPECTRUM (WEIGHT,IFLAG,1)
      ENDIF

C  PUSH PARTICLE TO POINT OF COLLISION, EITHER DELTA OR REAL


      X0=X0+VELX*ZTC
      Y0=Y0+VELY*ZTC
      Z0=Z0+VELZ*ZTC
      TIME=TIME+ZTC/VEL
      IF (LEVGEO.LE.3.AND.NLPOL) THEN
        IPOLG=NPCELL
      ELSEIF (NLPLG) THEN
        IPOLG=EIRENE_LEARC2(X0,Y0,NRCELL,NPANU,'FOLNEUT 2    ')
      ELSEIF (NLFEM) THEN
        IPOLG=0
      ELSEIF (NLTET) THEN
        IPOLG=0
      ENDIF
      NLSRFX=.FALSE.
      NLSRFY=.FALSE.
      NLSRFZ=.FALSE.
      NLSRFA=.FALSE.
      MRSURF=0
      MPSURF=0
      MTSURF=0
      MASURF=0
      MSURF=0
      IF (NLTRA) PHI=MOD(PHI-ATAN2(Z01,X01)+ATAN2(Z0,(RMTOR+X0)),PI2A)
C
  230 CONTINUE
C
C  PRE-COLLISION ESTIMATOR
C
      IF (NCLVI.GT.0) THEN
        WS=WEIGHT/SIGTOT
        CALL EIRENE_UPCUSR(WS,1)
      ENDIF
C
C
C  TEST FOR CORRECT CELL NUMBER AT COLLISION POINT
C  KILL PARTICLE, IF TOO LARGE ROUND-OFF ERRORS DURING
C  PARTICLE TRACING
C
      IF (NLTEST) THEN
        CALL EIRENE_CLLTST(IRET)
        IF (IRET .EQ. 1) GOTO 997
      ENDIF
C
C  SAMPLE FROM COLLISION KERNEL FOR NEUTRAL PARTICLES
C  AT PRESENT: NO SUPPRESSION OF ABSORPTION AT IONISATION
C  FIND NEW WEIGHT, SPECIES INDEX, VELOCITY AND RETURN
C
      IF (ITYP.EQ.1) THEN
        CALL EIRENE_COLATM(CFLAG,COLTYP,DIST)
      ELSEIF (ITYP.EQ.2) THEN
        CALL EIRENE_COLMOL(CFLAG,COLTYP)
      ELSEIF (ITYP.EQ.0) THEN
        CALL EIRENE_COLPHOT(CFLAG,COLTYP)
      ENDIF
      ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)

!  PARTICLE TYPE AND SPECIES MIGHT HAVE CHANGED
!  PREPARE POINTER FOR UNIFIED SUBROUTINES
      CALL EIRENE_SWITCH_PARTINFO
C
C  POST-COLLISION ESTIMATOR
C
      IF (LGPART.AND.(NCLVI.GT.0)) THEN
        WS=WEIGHT/SIGTOT
        CALL EIRENE_UPCUSR(WS,2)
        IF (NADSPC_CD >= 1) CALL EIRENE_UPDATE_SPECTRUM (WS,2,1)
      ENDIF
C
      IF (COLTYP.EQ.2.) GOTO 700
C
      GOTO 100
C
C  SIMULATION OF COLLISION EVENT FINISHED
C
C
C  ..............................................................
C  .
C  .  INCIDENT ONTO SURFACE
C  ..............................................................
C
  380 CONTINUE
C
C
c     PR= cond. exp probability to reach this surface. PR=1. by default
      PPR=PR
      IF (ILIIN(MSURF).LE.-2) PPR=PR*SG
C
C  UPDATE EFFLUXES ONTO SURFACE AND REFLECT PARTICLE
C
!pb   CALL EIRENE_ESCAPE(PPR,SG,*100,*104,*512)
      CALL EIRENE_ESCAPE(PPR,SG,IRET)
      IF (IRET == 1) GOTO 100
      IF (IRET == 2) GOTO 104
      IF (IRET == 3) GOTO 512
C
C   GOTO 100: START NEW TRACK OF NEUTRAL PARTICLE
C   GOTO 104: CONTINUE THIS TRACK, TRANSPARENT SURFACE IS CROSSED
C   GOTO 512: RESTORE PREVIOUS COLLISION DATA,
C             CONDITIONAL EXPECTATION ESTIMATOR WAS USED
C
      NLTRJ = .FALSE.
      TRAJ(ITRJ)%TRJ%NO_SURF = MSURF
      RETURN
C
C
C  ...................................................
C  .                                                 .
C  .  CONDITIONAL EXPECTATION ESTIMATOR  500 -- 599  .
C  ...................................................
C
C
C
C
C   SAVE PRE-COLLISION DATA OF FIRST COLLISION ALONG CONDITIONAL TRACK
C   AT THIS POINT:  JCOL= NO. OF TRACK SEGMENT IN RANGE J= 1:NCOU,
C   IN WHICH 1ST COLLISION FOUND.
  505 CONTINUE
      IF (NCOU.GT.1) THEN
        ZMFP=1./XSTORV2(NSTORV,JCOL)
      ELSE
C  IN CASE NCOU.EQ.1: XSTORV HAS NOT BEEN STORED ONTO XSTORV2
        ZMFP=1./XSTORV(NSTORV)
      ENDIF

      NPCLLC=1
      NTCLLC=1
      IF (NLPOL) NPCLLC=NPCOLC
      IF (NLTOR) NTCLLC=NTCOLC
      ZDT1C=(ZLOG-ZINT2)*ZMFP
      ZTC=ZT+ZDT1C
      X0C=X0+VELX*ZTC
      Y0C=Y0+VELY*ZTC
      Z0C=Z0+VELZ*ZTC
      TIMEC=TIME+ZTC/VEL
      NRCLLC=NRCELL
      NACLLC=NACELL
      NBLCKC=NBLOCK
      NCELLC=NCELL
      ITIMEC=ITIME
      IFPTHC=IFPATH
      IUPDTC=IUPDTE
      IPERIDC=IPERID
      VELXC=VELX
      VELYC=VELY
      VELZC=VELZ
      VELC=VEL
      E0C=E0
      GENRC=XGENER
      WEIGHC=WEIGHT
      IF (NLTRA)
     .  PHIC=MOD(PHI-ATAN2(Z01,X01)+ATAN2(Z0C,(X0C+RMTOR)),PI2A)
      IF (NCOU.GT.1) THEN
        XSTORC(:,:) = XSTOR2(:,:,JCOL)
        XSTORVC(:)  = XSTORV2(:,JCOL)
      ELSE
C  IN CASE NCOU.EQ.1: XSTORV HAS NOT BEEN STORED ONTO XSTORV2
        XSTORC(:,:) = XSTOR(:,:)
        XSTORVC(:)  = XSTORV(:)
      ENDIF
      IF (NLTRC) THEN
        PSAVE=PHI
        TSAVE=TIME
        NPSAVE=NPCELL
        NTSAVE=NTCELL
        IF (NLTRA) PHI=PHIC
        TIME=TIMEC
        NPCELL=NPCLLC
        NTCELL=NTCLLC
cym
!$OMP CRITICAL
        CALL EIRENE_CHCTRC(X0C,Y0C,Z0C,16,13)
!$OMP END CRITICAL
        PHI=PSAVE
        TIME=TSAVE
        NPCELL=NPSAVE
        NTCELL=NTSAVE
      ENDIF
      ICOL=1
      IFLAG=5
C  TRACK COMPLETED ?
      IF (NCOU.GE.NCOUS) GOTO 215
C
C  TRACK NOT COMPLETED BECAUSE OF WMINC CRITERION
C  UPDATE CONTRIBUTION TO VOLUME-AVERAGED ESTIMATORS
C
      IF (IUPDTE.GE.1) THEN
        CALL EIRENE_UPDATE(XSTOR2,XSTORV2,IFLAG)
        IF (NADSPC_CD >= 1) CALL EIRENE_UPDATE_SPECTRUM (WEIGHT,IFLAG,1)
      ENDIF
      GOTO 216
C
C   RESTORE PRE-COLLISION DATA AND SAMPLE FROM COLLISION KERNEL
  512 X0=X0C
      Y0=Y0C
      Z0=Z0C
      TIME=TIMEC
      NLSRFX=.FALSE.
      NLSRFY=.FALSE.
      NLSRFZ=.FALSE.
      NLSRFA=.FALSE.
      MSURF=0
      MRSURF=0
      MPSURF=0
      MTSURF=0
      MASURF=0
      NRCELL=NRCLLC
      NPCELL=NPCLLC
      NTCELL=NTCLLC
      NACELL=NACLLC
      NBLOCK=NBLCKC
      NBLCKA=NSTRD*(NBLOCK-1)+NACELL
      NCELL=NCELLC
      ITIME=ITIMEC
      IFPATH=IFPTHC
      IUPDTE=IUPDTC
      IPERID=IPERIDC
      VELX=VELXC
      VELY=VELYC
      VELZ=VELZC
      VEL=VELC
      E0=E0C
      XGENER=GENRC
      WEIGHT=WEIGHC
      IF (LEVGEO.LE.3.AND.NLPOL) THEN
        IPOLG=NPCELL
      ELSEIF (NLPLG) THEN
        IPOLG=EIRENE_LEARC2(X0,Y0,NRCELL,NPANU,'FOLNEUT 3    ')
      ELSEIF (NLFEM) THEN
        IPOLG=0
      ELSEIF (NLTET) THEN
        IPOLG=0
      ENDIF
      IF (NLTRA) PHI=PHIC
      XSTOR(:,:) = XSTORC(:,:)
      XSTORV(:)  = XSTORVC(:)
cym
      IF (NLTRC) THEN
!$OMP CRITICAL
        CALL EIRENE_CHCTRC(X0,Y0,Z0,0,14)
!$OMP END CRITICAL
      ENDIF
      ICOL=0
      LGPART=.TRUE.
      NLTRJ = .FALSE.
      GOTO 230
C
  700 CONTINUE
C  REGULAR STOP IN SUBR. FOLNEUT, CONTINUE IN SUBR. MCARLO
      RETURN
C
  800 CONTINUE
C  REGULAR STOP IN SUBR. FOLNEUT, STOP HISTORY, CENSUS ARRAY FULL
      IF (ICOL.EQ.1.AND..NOT.LGLAST) GOTO 512
      LGPART=.FALSE.
      WEIGHT=0.
      RETURN
C
  990 CONTINUE
cym
!$OMP CRITICAL
      CALL EIRENE_LEER(1)
      CALL EIRENE_MASAGE('ERROR IN FOLNEUT, ZDT1 OR NCELL OUT OF RANGE')
      CALL EIRENE_MASAGE('PARTICLE IS KILLED')
     
      write(iunout,*) 'ERROR for NPANU,thread =',NPANU,
     .                omp_get_thread_num()
      WRITE (iunout,*) 'ERROR NPANU,NCELL,ZDT1,ZTST,TL,TS '
      WRITE (iunout,'(I8,1X,I6,1P,4(1X,1E14.7))')
     .                  NPANU,NCELL,ZDT1,ZTST,TL,TS
      CALL EIRENE_MASJ4('NRCELL,NPCELL,NTCELL,NACELL     ',
     .                   NRCELL,NPCELL,NTCELL,NACELL)
!$OMP END CRITICAL
cym
      GOTO 995
C
 9911 CONTINUE
!$OMP CRITICAL
      CALL EIRENE_LEER(1)
      CALL EIRENE_MASAGE('ERROR IN FOLNEUT, NO INTERSECTION FOUND')
      CALL EIRENE_MASAGE('PARTICLE IS KILLED')
      WRITE (iunout,*) 'NPANU,NCELL,NRCELL,NPCELL,NTCELL '
      WRITE (iunout,*)  NPANU,NCELL,NRCELL,NPCELL,NTCELL
!$OMP END CRITICAL
      GOTO 995
C
 9912 CONTINUE
!$OMP CRITICAL
      CALL EIRENE_LEER(1)
      CALL EIRENE_MASAGE('ERROR IN FOLNEUT, DAMAGED CELL HIT')
      CALL EIRENE_MASAGE('PARTICLE IS KILLED')
      WRITE (iunout,*) 'NPANU,NCELL,NRCELL,NPCELL,NTCELL '
      WRITE (iunout,*)  NPANU,NCELL,NRCELL,NPCELL,NTCELL
!$OMP END CRITICAL
      GOTO 995
C
  992 CONTINUE
!$OMP CRITICAL
      CALL EIRENE_LEER(1)
      CALL EIRENE_MASAGE('ERROR IN FOLNEUT, SURFACE CONFLICT')
      CALL EIRENE_MASR2('TL,TS           ',TL,TS)
      WRITE (iunout,*) 'NPANU ',NPANU
!$OMP END CRITICAL
      ZT=TL
      GOTO 995
C
  995 WRITE (iunout,*) 'MRSURF,MPSURF,MTSURF,MASURF ',
     .             MRSURF,MPSURF,MTSURF,MASURF
      X0ERR=X0+ZT*VELX
      Y0ERR=Y0+ZT*VELY
      Z0ERR=Z0+ZT*VELZ
      IF (NLTRC) THEN
cym
!$OMP CRITICAL
        CALL EIRENE_CHCTRC(X0ERR,Y0ERR,Z0ERR,16,18)
!$OMP END CRITICAL
      ELSE
!$OMP CRITICAL
        WRITE (iunout,'(A,1P,4(1X,1E14.7))') 'X0,Y0,Z0,ZT ',X0,Y0,Z0,ZT
        WRITE (iunout,'(A,1P,3(1X,1E14.7))') 'VELX,VELY,VELZ ',
     .                                        VELX,VELY,VELZ
        WRITE (iunout,'(A,1P,3(1X,1E14.7))') 'X0ERR,Y0ERR,Z0ERR ',
     .                                        X0ERR,Y0ERR,Z0ERR
!$OMP END CRITICAL
      ENDIF
      GOTO 999
  997 CALL EIRENE_MASAGE('ERROR IN FOLNEUT, DETECTED IN SUBR. CLLTST')
      CALL EIRENE_MASAGE('PARTICLE IS KILLED')
C   DETAILED PRINTOUT ALREADY DONE FROM SUBR. CLLTST
cym
      IF (NLTRC) THEN
!$OMP CRITICAL
        CALL EIRENE_CHCTRC(X0,Y0,Z0,16,18)
!$OMP END CRITICAL
      ENDIF
      GOTO 999
C
  998 CALL EIRENE_MASAGE('ERROR IN FOLNEUT, SPECIES INDEX OUT OF RANGE')
      CALL EIRENE_MASJ5
     . ('NPANU, ITYP, IATM, IMOL, IPHOT          ',
     .   NPANU, ITYP, IATM, IMOL, IPHOT)
      GOTO 999
C
  999 CONTINUE
      PTRASH(ISTRA)=PTRASH(ISTRA)-WEIGHT
      ETRASH(ISTRA)=ETRASH(ISTRA)-WEIGHT*E0
      LGPART=.FALSE.
      CALL EIRENE_LEER(1)
      RETURN
      END
