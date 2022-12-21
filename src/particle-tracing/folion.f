cdr aug. 18   bug fix: remove virtual neutral background species
cdr           from coulomb collision frequency evaluation
cdr aprl.18   bug fix re. parallel distance (zt,ztc,mfp,...) and
cdr           scoring distance clpd (full gyro motion distance)
cdr           clpd  is switched back and forth. Needs clean-up.
cdr Oct. 17   minor sync with folneut
cdr           started: implementation of QSS branch: folstat_ion.f  not ready

cdr Nov. 15:  check again bgk solution for energy relaxation: mass factor, exponent ??
cdr           also: manual. to be done: remove static loop from folneut and folion.
c  nov. 2015:  fnui collision frequency: retain individual frequencies, for
c              all background species: fnuiar(ipls)
c  April 2015:  call to escape at periodicidy surfaces:  with reduced velocity, lcart=f
c               no gyro phase sampling then.
c               also for proper printout from chctrc for trace ions.

c  njump=3, for internal grid surface und timusr. reset time=0
c  error exit from fpkcol: goto 9991, da alles bereits in fpkcol erledigt (ptrash....)
C  OCT 14.:  cell-based spectra scoring called only if cell-based spectra are defined
!DR  eps12 --> eps6 for testing cosine of angle of incidence.
!DR  levgeo=4:  if nlsrfx: correction of nrcell for SG gt.0 SG lt.eps6
c
c....
c
C  MAY05: CALL UPDATE FROM STATIC LOOP WITH IFLAG=4 (RATHER =1)
C         WG. COLL EST. ON 1ST FLIGHT AFTER BIRTH.
C  Sept 05: also vel=velpar before call  to ...col  routines.
!PB 12.01.06: calls to UPDATE_SPECTRUM introduced for cell-based spectra
!PB 18.04.06: xstorv=0 in "vacuum region" added
!DR  4.08.06: check v_par=0, otherwise stop trajectory (lable 992)
!DR 10.08.06: cut-off Ti with T_vac for collision frequency, for
!             ion tracing in vacuum region
!PB 28.09.06: sg corrected for levgeo=4 and levgeo=5
!DR 09.02.07: not only the direction, but also the magnitute of velocity
!             is reset to full cartesian velocity in subr. NEWFIELD
!PB 22.03.07: LEVGEO=6 --> LEVGEO=10

!DR: introduce LCART=TRUE:
!                              velx,vely,velx,vel: "true particle velocities"
!                              in this case the reduced "guiding centre" velocity vector
!                              is stored in: velxgs, velygs, velzgs, velgs, velg(3)
!              LCART=FALSE:
!                              velx,vely,velx,vel: reduced "guiding centre velocities"
!                              i.e. excluding the gyromotion.
!                              velperp and vrelpar are parameters to solve
!                              (e.g. numerically) for
!                              guiding center equation.
!                              in this latter case the last "true" velocity vector
!                              is stored in: velxts, velyts, velzts, velts, velt(3)
!  TRUE (full) VELOCITIES ARE NEEDED IN FPATHI ROUTINES, AS WELL AS AT SOLID BOUNDARIES.
!  ONLY REDUCED VELOCITIES (GUIDING CENTRE) AT ALL TRANSPARENT BOUNDARIES AND TO
!  PUSH PARTICLES
!DR  eps12 --> eps6 for testing cosine of angle of incidence.
!DR  levgeo=4:  if nlsrfx: correction of nrcell for SG gt.0 SG lt.eps6

C  .......................................................................................
C  DIFFERENCES FROM SUBR. FOLNEUT:
C    0) INTRODUCE PARAMETERS VELPAR, VELPER:
C       VELOCITY PARALLEL AND PERP TO B FIELD, RESP.
C    1) REDUCED EQ. OF MOTION: A) MOTION ALONG B FIELD: VEL= VELPAR
C                              B) GUIDING CENTRE, INCL DRIFTS (EXPL. EULER: JOSEF)
C                              C) FULL GYRO MOTION (CORRECTIONS) NEAR TARGETS (TO BE DONE)
C    2) ADDITIONALLY: "FOKKER-PLANCK COLLISIONS", ISRFCL=4
C                              A) LANGER MODEL NF, ANALYTICAL
C                              B) TRUBNIKOV REFINED, SEMI-ANALYTICAL
C                              C) BINARY: TAKIZUKA  (BENJAMIN)
C                              D) HYBRID: PARTICLE-FLUID-FOKKER-PLANCK (JOSEF)
C  .......................................................................................

C
      SUBROUTINE EIRENE_FOLION
C
C     CHARGED PARTICLE, LAUNCHED AT X0,Y0,Z0 IN CELL NRCELL, IPOLG,
C     IPERID, NPCELL, NTCELL, NACELL, NBLOCK, WITH VELOCITY VELX,VELY,VELX
C     IS FOLLOWED.
C     (MODULE: COMPRT.F)
c
c
c
C  ON INPUT:
C     ITYP=3
C     IC_NEUT = 0  NEWBORN CHARGED PARTICLE, OR CONTINUATION FROM NEUTRAL PARTICLE FULL TRACK
C     IC_NEUT > 0  CONTINUATION FROM NEUTRAL PARTICLE WHICH WAS IN STATIC LOOP
C     IC_NEUT < 0  CONTINUATION FROM TEST PARTICLE IN DIFFUSION MODE (TO BE WRITTEN)
C
C  ON OUTPUT:
C
C     LGPART=TRUE
C           ITYP=0  NEXT GENERATION PHOTON IPHOT IS GENERATED
C           ITYP=1  NEXT GENERATION ATOM IATM IS GENERATED
C           ITYP=2  NEXT GENERATION MOLECULE IMOL IS GENERATED
C     LGPART=FALSE
C           ITYP=4  NO NEXT GENERATION TEST PARTICLE IS GENERATED
C                   (PARTICLE ABSORBED IN BULK ION SPECIES)
c
c  at 100 :   start a new trace ion, velocity is given as full cartesian vector, lcart=true
c  at 1004:   reduced (guiding centre) velocities and B field are now set for particle. lcart=false.
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
      USE EIRMOD_CFPLK
      USE EIRMOD_CLOGAU
      USE EIRMOD_CRAND
      USE EIRMOD_CINIT
      USE EIRMOD_CUPD
      USE EIRMOD_CPOLYG
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      USE EIRMOD_CSPEZ
      USE EIRMOD_CZT1
      USE EIRMOD_CTETRA
      USE EIRMOD_COMPRT
      USE EIRMOD_COMNNL
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
      use eirmod_colion

      IMPLICIT NONE

C     REAL(DP) :: a,aa,aaa
c     REAL(DP) :: fnueqi,fnueqi_1,fnueqi_2
      REAL(DP) :: CFLAG(7,MSTOR0)
      REAL(DP) :: AX(2)
      REAL(DP) :: XSTOR2(MSTOR1,MSTOR2,N2ND+N3RD),
     .            XSTORV2(NSTORV,N2ND+N3RD)
      REAL(DP) :: COSIN, XLI, YLI, ZLI,
     .          PR, WS, X0ERR, Y0ERR, Z0ERR,
     .          FNUI,
     .          VELXS, VELYS, VELZS, VELS,
     .          SG,
     .          VCOS,
     .          ZLOG, ZINT1, ZEP1, ZTST, ZINT2,
     .          ZMFP, EIRENE_FPATH, ZTC,
     .          DELFAC, TIFAC,
     .          SCOS_NEW, EPSLIM
      INTEGER :: ISTS, NCOUS, ICOU, J, JJ, IPL,
     .           ICO, NLI, NLE, JCOL, NRC, 
     .           NRCOLD, IPLTI, I, IM, ICOUN,
     .           EIRENE_LEARC2, COLTYP,
     .           indf, NJUMP_EMC3 = 0, IRET, IRT_STAT, KK
      LOGICAL :: LCNDEXP


c  no conditional expectation estimators for test ions

c  all cell indices must be known at this point
c  tentatively assume: a next generation particle will be born

c  IC_NEUT, IC_ION: counter for generations within static loop
      IC_ION=IC_NEUT
C  XGENER: COUNTER FOR GENERATION LIMIT
      LCART=.TRUE.

  100 LGPART=.TRUE.
      IC_ION=IC_ION+1
c  full cartesian velocity vector VEL,VELX,VELY,VELZ at this point
c  either a new particle, or back to here from collide, escape, fpkcol,
c  with a new full (cartesian) velocity vector.
c
      IF (.NOT.LCART) GOTO 9921

      XGENER=0.D0

C  CHECK FOR VALID SPECIES INDEX
      IF (ITYP.EQ.3.AND.(IION.LE.0.OR.IION.GT.NIONI)) GOTO 998
C
C  THE CELL NUMBER NRCELL, IPOLG, IPERID, NPCELL, NTCELL, NACELL, NBLOCK
C  WAS ALREADY SET IN CALLING SUBROUTINE MCARLO
C
C  IF NLSRFX, SURFACE INDEX MRSURF MUST BE DEFINED AT THIS POINT
C  IF NLSRFY, SURFACE INDEX MPSURF MUST BE DEFINED AT THIS POINT
C  IF NLSRFZ, SURFACE INDEX MTSURF MUST BE DEFINED AT THIS POINT
C  IF NLSRFA, SURFACE INDEX MASURF MUST BE DEFINED AT THIS POINT
C
      ICO=0  ! counter for particles sitting on surface: allow for two attemps in srfchk.
      EPSLIM=EPS6
 1005 NUPC(1)=NPCELL-1+(NTCELL-1)*NP2T3
      NCELL=NRCELL+NUPC(1)*NR1P2+NBLCKA
      IF (LDAMCEL(NCELL)) GOTO 9912
      IF (NCELL.GT.NSBOX.OR.NCELL.LT.1) GOTO 991

c  find direction parallel and perpendicular to B field, and velocity components
c  i.e. convert cartesian velocity unit vector VELX,VELY,VELX into
c  parallel and perpendicular unit velocity components VELPAR
c  find B field in cell NCELL
      CALL EIRENE_NEWFIELD(X0,Y0,Z0,VELS,0)

      VELXS=VELX
      VELYS=VELY
      VELZS=VELZ
      VELS=VEL

C  SIGPAR: SIGN OF PARALLEL VELOCITY WITH RESPECT TO B
c  calculating the angle between full velocity and B field
c  BBX, BBY, BBZ are normalized!
      VCOS = VELX*BBX + VELY*BBY + VELZ*BBZ
      IF (ABS(VCOS).LT.EPS30) GOTO 992
      SIGPAR=SIGN(1._DP,VCOS)
      VELPAR=ABS(VEL*VCOS)
      VELPER=SQRT(MAX(0._DP,VEL**2 - VELPAR**2))
c  VELOCITY WITH RESPECT TO B FIELD IS NOW DEFINED:
c  VELPAR: full parallel velocity, absolute value
c  VELPER: full perpendicular velocity, always non-negative
c  SIGPAR: sign of parallel velocity with respect to B

C  NOW REDUCED VELOCITY: GUIDING CENTRE APPROXIMATION

c  APPROXIMATION A)
c  use B field line as trajectory
c  VLXPAR,VLYPAR,VLZPAR gives the direction of the full parallel velocity
c  in Cartesian coordinates - absolute value is not correct!!!
      VLXPAR=SIGPAR*BBX
      VLYPAR=SIGPAR*BBY
      VLZPAR=SIGPAR*BBZ
c  VL_PAR: parallel unit speed vector, VL_PAR = SIG*B
c  VL_PAR = (/ VLXPAR, VLYPAR, VLZPAR /)

c  set ion energy = parallel energy of the ionized test particle
c
      E0PAR=CVRSSI(IION)*VELPAR*VELPAR
C
 1004 CONTINUE

C  FOLLOW MOTION OF TEST ION OR "STATIC APPROXIMATION"?
C     IF (ITYP.EQ.3) THEN
        LOGION(IION,ISTRA)=.TRUE.
C       NLPR=   : NOT AVAILABLE
        NRC=NRCI(IION)
C     ENDIF
        
      IF (IFPATH.NE.1.OR.NRC.LT.0) GOTO 1002
C  STOP STATIC LOOP AFTER 100 GENERATIONS LATEST, TO AVOID ACCIDENTAL INFINITE LOOPS
      IF (IC_ION.GT.100) GOTO 1002

      IF (NFOLI(IION).EQ.-1) GOTO 1001 ! go to static loop
C
C  the particle may be sitting exactly on a surface (nlsrf...=.true.).
C
C  this part is special for ions: due to projection of velocity
C  onto Guiding Center Motion (or even onto B field) the correct
C  angle relative to surface may be lost (e.g. cosin lt 0 may result).
c  Also NINC may be different, depending on whether computed with full
c  or with reduced (guiding centre) velocity
C
C  Check if the correct new cell number NCELL is used, in case of nlsrf?
C  Fiddle around a bit with cell number and flight direction in this case.
c  using the reduced (guiding centre) velocity to find orientation
C  relative to surface, and possibly correct side of surface, i.e. cell
c  number. In that case: goto 1005 and try again with new cell number.
c  Else: continue at 1002

      CALL EIRENE_SRFCHK(VLXPAR,VLYPAR,VLZPAR,SG,ICO,EPSLIM,IRET)
      IF (IRET == 1) GOTO 1005
      IF (ic_ion.gt.2) THEN
        write (iunout,*) 'error re static loop, ic_ion=', ic_ion
        call eirene_exit_own(1)
      ENDIF
      GOTO 1002

c***********************************************************************
c  CORRECTIONS FOR PARTICLES SITTING EXACTLY ON SURFACES DONE.
c***********************************************************************

c  at this point: V_PARALLEL, V_PERP known,
c                 gyrophase: to be sampled, if needed
C
 1001 CONTINUE

C***********************************************************************
C  STATIC APPROXIMATION
C  SIMULATE NEXT COLLISION INSTANTANEOUSLY
C***********************************************************************

      CALL EIRENE_FOLSTAT_ION(IC_ION,VLXPAR,VLYPAR,VLZPAR,CFLAG,
     .                        IRT_STAT)
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
      IF (IC_ION.GT.1.AND.NLTRC.AND.TRCHST)
     .  WRITE (iunout,*) 'TRAJECTORY LEAVES STATIC LOOP, ITYP=',ITYP

C  IN CASE THAT THE PARTICLE WAS IN STATIC LOOP AND ON A SURFACE,
C  SOME MORE WORK NEEDS TO BE DONE, TO REVIVE IT TO FULL KINETIC ORBIT MODE.
      IF (IC_ION.GT.1.AND.
     .   (NLSRFX.OR.NLSRFY.OR.NLSRFZ.OR.NLSRFA)) THEN

C  PARTICLE CONTINUES FROM SURFACE AND FROM PREVIOUS "STATIC LOOP"
C  PREPARE CELL NUMBERS FOR FIRST FLIGHT
        IC_ION=0
        IC_NEUT=0
        SCOS_NEW = SIGN(1.D0,VLXPAR*CRTXG+VLYPAR*CRTYG+VLZPAR*CRTZG)
        IF (SCOS_SAVE.NE.SCOS_NEW) THEN
          SCOS=SCOS_NEW
          ZT=0.D0
          TL=0.D0
          IPOLGN=IPOLG
C PUSH PARTICLE TO SURFACE, USE REDUCED (GC) VELOCITY
          IF (LCART) THEN
            VELXS=VELX
            VELYS=VELY
            VELZS=VELZ
            VELS=VEL

            VELX=VLXPAR
            VELY=VLYPAR
            VELZ=VLZPAR
            VEL =VELPAR
            LCART=.FALSE.
          ENDIF
          IF (NLSRFA) THEN
            CALL EIRENE_ADDCOL (X0,Y0,Z0,SCOS,IRET)
            IF (IRET .EQ. 1) GOTO 101
            IF (IRET .EQ. 2) GOTO 380
          ELSEIF (NLSRFX) THEN
            select case (LEVGEO)
            case (:3)
              ISTS=INMP1I(MRSURF,IPCELL,ITCELL)
              MSURFG=NPCELL+(NTCELL-1)*NP2T3
              IF (ILIIN(NLIM+ISTS) .NE. 0) THEN
                 CALL EIRENE_STDCOL (ISTS,1,SCOS,IRET)
                 IF (IRET .EQ. 1) GOTO 101
                 IF (IRET .EQ. 2) GOTO 380
              ENDIF
            case (4)
              ISTS=ABS(INMTI(IPOLGN,MRSURF))  !dr NLIM already added in ISTS ?
              MSURFG=INSPAT(IPOLGN,MRSURF)
              IF (ILIIN(ISTS) .NE. 0) THEN
                 CALL EIRENE_STDCOL (ISTS,1,SCOS,IRET)
                 IF (IRET .EQ. 1) GOTO 101
                 IF (IRET .EQ. 2) GOTO 380
              ENDIF
             case (5)
              ISTS=ABS(INMTIT(IPOLGN,MRSURF)) !dr NLIM already added in ISTS ?
C             MSURFG= ??
              IF (ILIIN(ISTS) .NE. 0) THEN
                 CALL EIRENE_STDCOL (ISTS,1,SCOS,IRET)
                 IF (IRET .EQ. 1) GOTO 101
                 IF (IRET .EQ. 2) GOTO 380
              ENDIF
            case (10)
              ISTS=INMP1I(MRSURF,IPCELL,ITCELL)
C             MSURFG= ??
              IF (ILIIN(NLIM+ISTS) .NE. 0) THEN
                 CALL EIRENE_STDCOL (ISTS,1,SCOS,IRET)
                 IF (IRET .EQ. 1) GOTO 101
                 IF (IRET .EQ. 2) GOTO 380
              ENDIF
            end select
          ELSEIF (NLSRFY) THEN
            ISTS=INMP2I(IRCELL,MPSURF,ITCELL)
            MSURFG=NRCELL+(NTCELL-1)*NR1P2
            IF (ILIIN(NLIM+ISTS) .NE. 0) THEN
               CALL EIRENE_STDCOL (ISTS,2,SCOS,IRET)
               IF (IRET .EQ. 1) GOTO 101
               IF (IRET .EQ. 2) GOTO 380
            ENDIF
          ELSEIF (NLSRFZ) THEN
            ISTS=INMP3I(IRCELL,IPCELL,MTSURF)
            MSURFG=NRCELL+(NPCELL-1)*NR1P2
            IF (ILIIN(NLIM+ISTS) .NE. 0) THEN
               CALL EIRENE_STDCOL (ISTS,3,SG,IRET)
               IF (IRET .EQ. 1) GOTO 101
               IF (IRET .EQ. 2) GOTO 380
            ENDIF
          ENDIF
C         WRITE (IUNOUT,*) 'FOLION: I SHOULD NOT BE HERE'
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
C  EACH TEST ION TRACK STARTS AT THIS POINT, IC_ION=0 HERE

  101 CONTINUE
C     IF (ITYP.EQ.3) THEN
        LOGION(IION,ISTRA)=.TRUE.
C       NLPR=   : NOT AVAILABLE
        NRC=NRCI(IION)
C     ENDIF
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
      AX(1)=1.
      AX(2)=1.

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
      NCELL=NRCELL+((NPCELL-1)+(NTCELL-1)*NP2T3)*NR1P2+NBLCKA
      IF (LDAMCEL(NCELL)) GOTO 9912
      CALL EIRENE_NEWFIELD(X0,Y0,Z0,VELS,1)
C  AT THIS POINT: LCART=F

      NJUMP=0
      IF (NJUMP_EMC3 == 3) THEN
        NJUMP = 3
        NJUMP_EMC3 = 0
      ENDIF
      DO I=1,NIMINT
        IM=IIMINT(I)
        TIMINT(IM)=0._DP
        IIMINT(I)=0
      END DO
      NIMINT = 0
      TT=1.D30
      TL=1.D30
      TS=1.D30
      ZTST=1.D30
      ZT=0.0
C
      NCOU=1
      NUPC(1)=0
      NCOUNT(1)=1
      NCOUNP(1)=1
      ISRFCL=-1
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
     .               VLXPAR,VLYPAR,VLZPAR,VELPAR,
     .               MASURF,XLI,YLI,ZLI,SG,TL,NLTRC,LCNDEXP)
C       NLPR: NOT AVAILABLE FOR TEST IONS
c
        ZDT1=TL
        ZTST=TL
        CLPD(1)=ZDT1
        IF (MASURF.NE.0) ISRFCL=1
      ENDIF
C
C TT: DISTANCE UNTIL NEXT TIMESTEP LIMIT IS REACHED
C     USE VEL_GC INSTEAD OF VEL, BECAUSE ORBIT IS COMPUTED WITH REDUCED (GC) VELOCITY
C     LATER: VELPAR --> VEL_GC
      IF (LGTIME) THEN
        TT=(DTIMVI-TIME)*VELPAR
        IF (TT.LT.ZTST) THEN
          ZDT1=TT
          ZTST=TT
          CLPD(1)=ZDT1
          ISRFCL=2
        ENDIF
      ENDIF
C
C FNUI: collision frequency with background ions.

      FNUI   = 1.D-30
      IF (NRC.GE.0) THEN
        DO IPL=1,NPLSI
          IPLTI=MPLSTI(IPL)
ctest     ti=200
ctest     ni=1e14
ctest     ea=0.1
ctest     iion=1
ctest     ipls=1
ctest     nmassi(1)=16.
ctest     nmassp(1)=1.
ctest     a=fnueqi(1.d14,200.d0)
ctest     a=a*(1.+1./16.)**0.5-a*1.5*200./0.1
ctest     aa=fnueqi_1(0.1d0,1.d14,200.d0,1,1)
ctest     aaa=fnueqi_2(0.1d0,1.d14,200.d0,1,1)
ctest     write (iunout,*) 'a,aa,aaa', a,aa,aaa
ctest     stop

C  default Coulomb collision model (simple energy relaxation, e.g. also: NRC=0)
C  Set Coulomb collisions (energy relaxation) frequencies.
C  Exclude vacuum region and virtual neutral background species
          IF (.NOT.LGVAC(NCELL,IPL) .AND. (NCHRGP(IPL) > 0) ) THEN
            FNUIAR(IPL) = FNUEQI(DIIN(IPL,NCELL),TIIN(IPLTI,NCELL))
            FNUI=FNUI+FNUIAR(IPL)
          END IF

        ENDDO
      ENDIF
C TAUE: RELAXATION TIME
      TAUE=1./FNUI
C STEPSIZE=0.1*VEL_PARALLEL*TAUE, I.E. 10 COULOMB COLLISIONS PER RELAX.TIME
C TF: DISTANCE UNTIL NEXT COULOMB COLLISION
C DELFAC: INCREASE STEPSIZE AS E0 APPROACHES 1.5 * TI
      TIFAC  = MAX(TVAC,TIIN(1,NCELL))
      DELFAC =1.5_DP*TIFAC/ABS(E0-1.5_DP*TIFAC+EPS60)
C  DELTA_T = TAUE*0.1*DELFAC
C  DELTA_S = DELTA_T * VELPAR  ! = TF
C     USE VELGS INSTEAD OF VEL, BECAUSE ORBIT IS COMPUTED WITH REDUCED (GC) VELOCITY
C     LATER: VELPAR --> VEL_GC
      TF=TAUE*VELPAR*0.1*DELFAC
      if (nldfst) tf=1.E-5_DP*vel

      IF (TF.LT.ZTST) THEN
        ZTST=TF
        ZDT1=TF
        CLPD(1)=ZDT1
        ISRFCL=4
      ENDIF
C
C  SCAN OVER SEGMENT

C  BEFORE THIS SCAN: ZTST, ZDT1, CLPD(1):  MAX. POSSIBLE DISTANCE, DUE TO TIME STEP, FP_COL OR ADD. SURF.
C
CCC  210 CONTINUE
C
C
C  TS:   DISTANCE TO NEXT RADIAL SURFACE OF STANDARD MESH
C  ZDT1: DISTANCE TRAVELLED IN CURRENT RADIAL CELL
C  ZT:   ACCUMULATED DISTANCE, UNTIL THIS SEGMENT
C
C  USE PARALLEL VELOCITY, I.E., COMPUTE PARALLEL DISTANCES IN GRID
C  THUS ZT,TS,ZTST,ZDT1,CLPD ETC. ARE PARALLEL DISTANCES
C  I.E., LCART=F AT THIS POINT
C
      IF (ITIME.EQ.1) THEN
c  switch to gc velocity
        IF (LCART) THEN
          VELXS=VELX
          VELYS=VELY
          VELZS=VELZ
          VELS =VEL
          VELX=VLXPAR
          VELY=VLYPAR
          VELZ=VLZPAR
          VEL =VELPAR
          LCART=.FALSE.
        ENDIF

        IF (NLRAD) THEN
          CALL EIRENE_TIMER(TS)
          IF (.NOT.LGPART) GOTO 9911
C
          IF (TL.LT.TS.OR.TT.LT.TS.OR.TF.LT.TS) THEN
            MRSURF=0
            IPOLGN=0
C  CHECK FOR INTERSECTION WITH ADDITIONAL SURFACE
            IF (TL.LE.TT.AND.TL.LE.TF) THEN
              ZDT1=TL-ZT
              TL=ZT+ZDT1
              ZTST=TL
              ISRFCL=1
C  INTERSECTION WITH TIME SURFACE. TIME LIMIT REACHED ?
            ELSEIF (TT.LT.TL.AND.TL.LE.TF) THEN
              ZDT1=TT-ZT
              TT=ZT+ZDT1
              ZTST=TT
              ISRFCL=2
C  Fokker-Planck collision, DIFFUSIVE STEP
            ELSEIF (TF.LT.TL.AND.TF.LE.TT) THEN
              ZDT1=TF-ZT
              TF=ZT+ZDT1
              ZTST=TF
              ISRFCL=4
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
C
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

c  switch to full velocity but gc velocity is not saved
        IF (.NOT.LCART) THEN
          VELX=VELXS
          VELY=VELYS
          VELZ=VELZS
          VEL =VELS
          LCART=.TRUE.
        ENDIF

      ENDIF
C
      IF (ZTST.GE.1.D30) GOTO 990
C
C  LOCAL MEAN FREE PATH
C
c  use parallel velocity, i.e., compute parallel mean free path
c  because clpd is the parallel distance in each cell (exclud. gyro)
c  etc. ... e.g LAMBDA(PARALLEL) = VEL(PARALLEL)/SIGV.
c  the collision frequency SIGV, however, must be computed using the
c  full test ion velocity vector, because it may depend upon the relativ
c  interaction energy: to be written
c  for interactions with electrons this is usually irrelevant

      IF (IFPATH.NE.1.OR.NRC.LT.0) THEN
        XSTORV(:)=0.D0
        DO J=1,NCOU
          JJ=J
cdr  next 2 lines added, Aug. 18. Strictly not necessary, but safer
cdr  (allows using NCELL later also in this case).
          NCELL=NRCELL+NUPC(J)*NR1P2+NBLCKA
          IF (LDAMCEL(NCELL)) GOTO 9912
          XSTOR2(:,:,J)=0.D0
          XSTORV2(:,J)=0.D0
          ZMFP=1.D10
          IF (NLPOL) NPCELL=NCOUNP(J)
          IF (NLTOR) NTCELL=NCOUNT(J)
C         VEL=VELS
          IF (.TRUE.) GOTO 213
        END DO
      ELSE
c switch to parallel gc velocity
        IF (LCART) THEN
          VELXS=VELX
          VELYS=VELY
          VELZS=VELZ
          VELS =VEL
          VELX=VLXPAR
          VELY=VLYPAR
          VELZ=VLZPAR
          VEL =VELPAR
          LCART=.FALSE.
        ENDIF
        DO J=1,NCOU
          JJ=J
          NCELL=NRCELL+NUPC(J)*NR1P2+NBLCKA
          IF (LDAMCEL(NCELL)) GOTO 9912
          ZMFP=EIRENE_FPATH(NCELL,CFLAG,J,NCOU)
          IF (NCOU.GT.1) THEN
            XSTOR2(:,:,J)=XSTOR(:,:)
            XSTORV2(:,J)=XSTORV(:)
          ENDIF

C  UPDATE INTEGRAL
          ZINT1=ZINT1+CLPD(J)*ZMFPI
C         IF (.NOT.NLPR) THEN
CCC         IF (ZINT1.GE.ZLOG) THEN
C  COLLISION IN SECTION J OF CURRENT TRACK
              IF (NLPOL) NPCELL=NCOUNP(J)
              IF (NLTOR) NTCELL=NCOUNT(J)

              VELX=VELXS
              VELY=VELYS
              VELZ=VELZS
              VEL =VELS
              LCART=.TRUE.
              IF (.TRUE.) GOTO 213
CCC         ENDIF
C  THESE NEXT TWO LINES CAN NEVER BE REACHED, BECAUSE ONLY ONE
C  CELL FOR EACH TRACK OF IONS (DISTINCT FROM FOLNEUT).
C  THEN (AT THE LATEST): ROTATION OF VELOCITY DUE TO NEW B FIELD
C           ZINT2=ZINT1
C           ZT=ZT+CLPD(J)
C         ELSEIF (JCOL.EQ.0) THEN
C   CONDITIONAL EXPECTATION ESTIMATOR FOR TEST IONS: TO BE WRITTEN
C         ENDIF
C
        END DO
        VELX=VELXS
        VELY=VELYS
        VELZ=VELZS
        VEL =VELS
        LCART=.TRUE.
      ENDIF
C
  213 CONTINUE
      NCOUS=NCOU
      NCOU=JJ

CCC  IF NO COLLISION, THEN: ENFORCE ONLY ONE STEP AT A TIME
      IF (ZINT1.LT.ZLOG.AND.NCOUS.GT.1) THEN
        MRSURF=0
        MPSURF=0
        MTSURF=0
        MASURF=0
        ISRFCL=0
        NINCX=0
        NINCY=0
        NINCZ=0
      ENDIF
CCC
C
C  CHECK FOR EVENT
C
C     IF (NLPR)    ......
      IF (ZINT1.GE.ZLOG) GO TO 220
C
      ZINT2=ZINT1

C  SET NEW ACCUMULATED FLIGHT LENGTH, TENTATIVE
      ZT=ZTST
C
C  RESET CLPD TO REAL PATH LENGTH OF FULL GYRO MOTION FOR SCORING
C  vel is the full velocity, velpar is the parallel velocity only

      DO 217 ICOU=1,NCOU
        CLPD(ICOU)=CLPD(ICOU)*VEL/VELPAR
  217 CONTINUE
C
C  UPDATE CONTRIBUTION TO VOLUME-AVERAGED ESTIMATORS
C
      IF (IUPDTE.GE.1) THEN
cdr     IFLAG= 3 ?
        CALL EIRENE_UPDATE(XSTOR2,XSTORV2,3)
        IF (NADSPC_CD >= 1) CALL EIRENE_UPDATE_SPECTRUM (WEIGHT,3,1)
      ENDIF
C
C  STOP TRACK ?
C
CDR: ALLE DISTANZEN IN ...COL routines sind parallele distanzen
CDR: Daher auch wg. x = x + dist/vel  parallele geschwindigkeiten.

c  switch to parallel gc velocity
      IF (LCART) THEN
        VELXS=VELX
        VELYS=VELY
        VELZS=VELZ
        VELS =VEL
        VELX=VLXPAR
        VELY=VLYPAR
        VELZ=VLZPAR
        VEL =VELPAR
        LCART=.FALSE.
      ENDIF
      IF (ISRFCL.EQ.1) THEN
c  will fpkcol change the collision with additional surface?
        CALL EIRENE_ADDCOL(XLI,YLI,ZLI,SG,IRET)
        IF (IRET .EQ. 1) GOTO 104
        IF (IRET .EQ. 2) GOTO 380
      ELSEIF (ISRFCL.EQ.2) THEN
!pb     CALL EIRENE_FPKCOL(               *104,*2215,*9991,3)
        CALL EIRENE_FPKCOL(IRET,3)
        IF (IRET .EQ. 1) GOTO 104
        IF (IRET .EQ. 2) GOTO 2215
        IF (IRET .EQ. 3) GOTO 9991
!PB 2215 CALL EIRENE_TIMCOL(AX(2),         *104,*800)
 2215   CALL EIRENE_TIMCOL(AX(2),IRET)
        IF (IRET .EQ. 1) GOTO 104
        IF (IRET .EQ. 2) GOTO 800
      ELSEIF (ISRFCL.EQ.3) THEN
!pb     CALL EIRENE_FPKCOL(               *104,*2216,*9991,3)
        CALL EIRENE_FPKCOL(IRET,3)
        IF (IRET .EQ. 1) GOTO 104
        IF (IRET .EQ. 2) GOTO 2216
        IF (IRET .EQ. 3) GOTO 9991
!PB 2216 CALL EIRENE_TORCOL(               *104)
 2216   CALL EIRENE_TORCOL(IRET)
        IF (IRET .EQ. 1) GOTO 104
      ELSEIF (ISRFCL.EQ.4) THEN
!pb     CALL EIRENE_FPKCOL(               *104,*100,*9991,0)
        CALL EIRENE_FPKCOL(IRET,0)
        IF (IRET .EQ. 1) GOTO 104
        IF (IRET .EQ. 2) GOTO 100
        IF (IRET .EQ. 3) GOTO 9991
      ENDIF

C changing back to full cartesian velocities
      VELX=VELXS
      VELY=VELYS
      VELZ=VELZS
      VEL =VELS
      LCART=.TRUE.
C
C  NO, CONTINUE TRACK
C
C  NEXT CELL - CHECK FOR ESCAPE OR NON-DEFAULT ACTING STANDARD SURFACE

c  DO THIS WITH REDUCED VELOCITY:
      IF (LCART) THEN
            VELXS=VELX
            VELYS=VELY
            VELZS=VELZ
            VELS =VEL
            VELX=VLXPAR
            VELY=VLYPAR
            VELZ=VLZPAR
            VEL =VELPAR
            LCART=.FALSE.
      ENDIF
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
             CALL EIRENE_STDCOL(ISTS,1,SG,IRET)
             IF (IRET .EQ. 1) GOTO 104
             IF (IRET .EQ. 2) GOTO 380
          ENDIF
        ENDIF

C  ESCAPE AT 2ND GRID SURFACE (Y OR POLOIDAL) NO. MPSURF
        ISTS=INMP2I(IRCELL,MPSURF,ITCELL)
        IF (NLPOL.AND.ISTS.NE.0) THEN
          SG=ISIGN(1,NINCY)
          NLSRFY=.TRUE.
          MSURFG=NRCELL+(NTCELL-1)*NR1P2
          IF (ILIIN(NLIM+ISTS) .NE. 0) THEN
             CALL EIRENE_STDCOL(ISTS,2,SG,IRET)
             IF (IRET .EQ. 1) GOTO 104
             IF (IRET .EQ. 2) GOTO 380
          ENDIF
        ENDIF

C  ESCAPE AT 3RD GRID SURFACE (Z OR TOROIDAL) MTSURF
        ISTS=INMP3I(IRCELL,IPCELL,MTSURF)
        IF (NLTOR.AND.ISTS.NE.0) THEN
          SG=ISIGN(1,NINCZ)
          NLSRFZ=.TRUE.
          MSURFG=NRCELL+(NPCELL-1)*NR1P2
          IF (ILIIN(NLIM+ISTS) .NE. 0) THEN
             CALL EIRENE_STDCOL(ISTS,3,SG,IRET)
             IF (IRET .EQ. 1) GOTO 104
             IF (IRET .EQ. 2) GOTO 380
          ENDIF
        ENDIF
C
C  ESCAPE AT GRID SURFACE BUILT FROM TRIANGLE SIDES IN X-Y PLANE: MRSURF
      case (4)
        IF (MRSURF > 0) THEN
          ISTS=ABS(INMTI(IPOLGN,MRSURF))
          IF (NLRAD.AND.ISTS.NE.0) THEN
            NLSRFX=.TRUE.
            MSURFG=INSPAT(IPOLGN,MRSURF)
            SG=SIGN(1._DP,VELX*PTRIX(IPOLGN,MRSURF)+
     .                    VELY*PTRIY(IPOLGN,MRSURF))
            IF (ILIIN(ISTS) .NE. 0) THEN
               CALL EIRENE_STDCOL(ISTS,1,SG,IRET)
               IF (IRET .EQ. 1) GOTO 104
               IF (IRET .EQ. 2) GOTO 380
            ENDIF
          ENDIF
        END IF

C  ESCAPE AT 3RD (Z OR TOROIDAL) GRID SURFACE FOR TRIANGULAR X-Y GRID OPTION: MTSURF
        IF (MTSURF > 0) THEN
          ISTS=INMTI3(IRCELL,MTSURF)
          IF (NLTOR.AND.ISTS.NE.0) THEN
            SG=ISIGN(1,NINCZ)
            NLSRFZ=.TRUE.
            MSURFG=NRCELL+(NPCELL-1)*NR1P2
            IF (ILIIN(NLIM+ISTS) .NE. 0) THEN
               CALL EIRENE_STDCOL(ISTS,3,SG,IRET)
               IF (IRET .EQ. 1) GOTO 104
               IF (IRET .EQ. 2) GOTO 380
            ENDIF
          ENDIF
        END IF
C
C  ESCAPE AT GRID SURFACE BUILT FROM TETRAHEDRA SIDES: MRSURF
      case (5)
        ISTS=ABS(INMTIT(IPOLGN,MRSURF))  !dr NLIM already added in ISTS ?
        IF (NLRAD.AND.ISTS.NE.0) THEN
          SG=SIGN(1._DP,VELX*PTETX(IPOLGN,MRSURF)+
     .                  VELY*PTETY(IPOLGN,MRSURF)+
     .                  VELZ*PTETZ(IPOLGN,MRSURF))
          NLSRFX=.TRUE.
C         MSURFG= ??
          IF (ILIIN(ISTS) .NE. 0) THEN
             CALL EIRENE_STDCOL(ISTS,1,SG,IRET)
             IF (IRET .EQ. 1) GOTO 104
             IF (IRET .EQ. 2) GOTO 380
          ENDIF
        ENDIF

C  ESCAPE TO GRID SURFACE ON USER-DEFINED GEOMETRY BLOCK: MRSURF
      case (10)
        ISTS=INMP1I(MRSURF,IPCELL,ITCELL)
        IF (NLRAD.AND.ISTS.NE.0) THEN
          SG=ISIGN(1,NINCX)
          NLSRFX=.TRUE.
          IF (ILIIN(NLIM+ISTS) .NE. 0) THEN
             CALL EIRENE_STDCOL(ISTS,1,SG,IRET)
             IF (IRET .EQ. 1) GOTO 104
             IF (IRET .EQ. 2) GOTO 380
          ENDIF
        ENDIF
      end select
C
      NRCELL=NRCELL+NINCX
      IF (NRCELL.GT.NR1STM.OR.NRCELL.LT.1) GOTO 990
C
CDR: SPLITTING AND COND.EXP.EST. NOT AVAILABLE FOR TEST IONS
C
C  CHECK IF WE HAVE ENCOUNTERED A SPLITTING ZONE
C     IF (NLSPLT(MRSURF).AND.NLEVEL.LT.MAXLEV.AND.ICOL.EQ.0) GOTO 330
C
C  SWITCH OFF CONDITIONAL EXP. ESTIMATOR ?
C     IF (AX(2).LT.WMINC) THEN
C       IF (ICOL.EQ.1) GOTO 512
C  NO COLLISION YET; RESTART AGAIN WITH COND. EXP. ESTIMATOR
C                    IN NEW CELL
C       AX(1)=1.
C       AX(2)=1.
C       JCOL=0
C     ENDIF
CCC

C  EARLIER CLPD WAS DISTANCE ALONG FULL GYROMOTION, FOR SCORING.
C  NOW WE NEED AGAIN THE PARALLEL DISTANCE, FOR TRACKING TO
C  POINT OF COLLISION OR TO SURFACE EVENT. (I.E. LCART=F)
      IF (.NOT.LCART) THEN
c        WRITE (IUNOUT,*) 'FROM FOLION: LCART ?'
c        WRITE (IUNOUT,*) 'WARNING: VEL IS ALREADY = VELPAR HERE'
c        WRITE (IUNOUT,*) VEL,VELPAR,VELS
         CLPD(1)=CLPD(1)*VELPAR/VELS
      ENDIF
      ZTC=CLPD(1)*VELPAR/VEL   !   this now does nothing: Velpar=vel here
      IF (LCART) THEN
        VELXS=VELX
        VELYS=VELY
        VELZS=VELZ
        VELS =VEL
        VELX=VLXPAR
        VELY=VLYPAR
        VELZ=VLZPAR
        VEL =VELPAR
        LCART=.FALSE.
      ENDIF
      GOTO 2211
CCC
CCC   GOTO 210
C
C  POINT OF COLLISION  220 -- 240
C
  220 CONTINUE
C
      CLPD(NCOU)=(ZLOG-ZINT2)*ZMFP
      ZTC=ZT+CLPD(NCOU)
C  RESET CLPD TO REAL (FULL) PATH LENGTH OF FULL GYRO MOTION FOR SCORING
cdr I do not understand: for scoring clpd should be full (gyro) distance.
cdr but if I rescale clpd with vels/velpar, then trace ion balances become
cdr much worse.
cdr   if (.not.lcart) then
        DO 221 ICOU=1,NCOU
cdr       CLPD(ICOU)=CLPD(ICOU)*VELS/VELPAR
          CLPD(ICOU)=CLPD(ICOU)*VEL/VELPAR
  221   CONTINUE
cdr   endif

      IF (IUPDTE.GE.1) THEN
        CALL EIRENE_UPDATE (XSTOR2,XSTORV2,4)
        IF (NADSPC_CD >= 1) CALL EIRENE_UPDATE_SPECTRUM (WEIGHT,4,1)
      ENDIF

C  PUSH PARTICLE TO POINT OF COLLISION, EITHER DELTA OR REAL. ZTC: PARALLEL (gc) DISTANCE)

 2211 CONTINUE
      X0=X0+VLXPAR*ZTC
      Y0=Y0+VLYPAR*ZTC
      Z0=Z0+VLZPAR*ZTC
      TIME=TIME+ZTC/VELPAR
      IF (LEVGEO.LE.3.AND.NLPOL) THEN
        IPOLG=NPCELL
      ELSEIF (NLPLG) THEN
        IPOLG=EIRENE_LEARC2(X0,Y0,NRCELL,NPANU,'FOLION 2     ')
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
CCC

C  DELTA EVENT AT CELL BOUNDARY: STOP TEST ION, AND RESTART WITH REFRESHED E AND B FIELDS

      IF (ZINT1.LT.ZLOG) THEN
C  CELL SURFACE HAS BEEN REACHED BEFORE COLLISION EVENT

        IF (NINCX.NE.0) THEN
C  IT WAS A "RADIAL" (1 ST) GRID SURFACE
          NLSRFX=.TRUE.
C  AT THIS POINT: NRCELL IS THE NEW CELL TO BE ENTERED
C                 FIND MRSURF: SURFACE OF CELL BOUNDARY
C                 BETWEEN OLD AND NEW CELL.
          select case (LEVGEO)
          case (:3)
            MRSURF=NRCELL
            IF (NINCX.EQ.-1) MRSURF=NRCELL+1
          case (4)
            NRCOLD=NRCELL-NINCX
            MRSURF=NCHBAR(IPOLGN,NRCOLD)
            IPOLG=NSEITE(IPOLGN,NRCOLD)
          case (5)
            NRCOLD=NRCELL-NINCX
            MRSURF=NTBAR(IPOLGN,NRCOLD)
            IPOLG=NTSEITE(IPOLGN,NRCOLD)
          case (10)
!PB         EXPLICITLY ALLOW FOR LEVGEO=10
!PB         NOTHING DONE FOR DELTA EVENT AT CELL BOUNDARY
            NJUMP_EMC3 = 3
          case default
            WRITE (iunout,*) 'DELTA EVENT AT CELL BOUNDARY '
            WRITE (iunout,*) 'FOR INVALID LEVGEO IN SUBR. FOLION. '
            CALL EIRENE_EXIT_OWN(1)
          end select

        ELSEIF (NINCZ.NE.0) THEN
C  IT WAS A "TOROIDAL" (3 RD) GRID SURFACE
          NLSRFZ=.TRUE.
          NTCELL=KUPC(1)+NINCZ
          IF (NINCZ == 1) THEN
            MTSURF=NTCELL
          ELSEIF (NINCZ.EQ.-1) THEN
            MTSURF=NTCELL+1
          ENDIF

C  IT WAS A "POLOIDAL" (2 ND) GRID SURFACE
        ELSEIF (NINCY.NE.0) THEN
          NLSRFY=.TRUE.
          select case (LEVGEO)
          case (1)
            NPCELL=JUPC(1)+NINCY
            IF (NINCY == 1) THEN
              MPSURF=NPCELL
            ELSEIF (NINCY.EQ.-1) THEN
              MPSURF=NPCELL+1
            ENDIF
          case (2:3)
            MPSURF=LUPC(1)
            IF (MUPC(1).EQ.1) NPCELL=NGHPLS(2,NRCELL,MPSURF)
            IF (MUPC(1).NE.1) NPCELL=NGHPLS(4,NRCELL,MPSURF)
C  PERIODICITY FOR LEVGEO=2 (TO BE WRITTEN IN MORE GENERAL TERMS)
            IF (NPCELL.EQ.0.AND.LEVGEO.EQ.2) THEN
              WRITE (iunout,*) 'should not be here '
              MPSURF=NP2ND
              NPCELL=NP2NDM
            ELSEIF (NPCELL.EQ.NP2ND.AND.LEVGEO.EQ.2) THEN
              WRITE (iunout,*) 'should not be here '
              MPSURF=1
              NPCELL=1
            ENDIF
            IF (LEVGEO.LE.3.AND.NLPOL) THEN
              IPOLG=NPCELL
            ELSEIF (LEVGEO.EQ.3.AND..NOT.NLPOL) THEN
              IPOLG=EIRENE_LEARC2(X0,Y0,NRCELL,NPANU,'FOLION neu   ')
            ENDIF
          end select

        ELSE   !NONE OF THE ninc_x,y,z flags are set,
cdr  all the nincx,...y,...z=0. This can happen only in levgeo=10,
cdr  for an internal surface which is only known to external geometry block but not to eirene
cdr  try to tell external code: particle on surface, but it is an old particle, which continues.
          NLSRFX=.TRUE.
          IF (LEVGEO .NE. 10) GOTO 994
          NJUMP_EMC3 = 3
        ENDIF

        IF (NLTRC) THEN
!$OMP CRITICAL
          CALL EIRENE_CHCTRC(X0,Y0,Z0,16,19)
!$OMP END CRITICAL
        ENDIF
        NUPC(1)=NPCELL-1+(NTCELL-1)*NP2T3
        NCELL=NRCELL+NUPC(1)*NR1P2+NBLCKA
        IF (LDAMCEL(NCELL)) GOTO 9912
C  DELTA COLLISION AT SURFACE DONE, NEW CELL FOUND (except in case levgeo 10 ?)

!pb     CALL EIRENE_FPKCOL(*104,*229,*9991,3)
        CALL EIRENE_FPKCOL(IRET,3)
        IF (IRET == 1) GOTO 104
        IF (IRET == 2) GOTO 229
        IF (IRET == 3) GOTO 9991

C  FIND NEW B FIELD, NEW REDUCED (GC) VELOCITY
  229   CONTINUE
C STORE NEW FULL VELOCITY
        VELS = VEL
        CALL EIRENE_NEWFIELD(X0,Y0,Z0,VELS,1)
cdr Warning: this call is probably incorrect in case of levgeo=10.
!   Jump to external (e.g. emc3) routine)
!   but there the cell number may be set only later.
!   In fpkcol a new B field may already have been set.
!   Futhermore: a new vel vector from fpkcol may get corrupted here.

        ICO = 0
        GOTO 1004
      ENDIF
CCC
C
  230 CONTINUE
C
C  PRE-COLLISION ESTIMATOR
C
      IF (NCLVI.GT.0) THEN
        WS=WEIGHT/SIGTOT
        CALL EIRENE_UPCUSR(WS,1,KK)
        IF (NADSPC_CD >= 1) CALL EIRENE_UPDATE_SPECTRUM (WS,1,1)
      ENDIF
C
C
C  TEST FOR CORRECT CELL NUMBER AT COLLISION POINT
C  KILL PARTICLE, IF TOO LARGE ROUND-OFF ERRORS DURING
C  PARTICLE TRACING
C
      IF (NLTEST) THEN
        CALL EIRENE_CLLTST(IRET)
        IF (IRET.EQ.1) GOTO 997
      ENDIF
C
C  SAMPLE FROM COLLISION KERNEL FOR TEST IONS
C  AT PRESENT: NO SUPPRESSION OF ABSORPTION AT IONIZATION
C  FIND NEW WEIGHT, SPECIES INDEX, VELOCITY AND RETURN
C
      CALL EIRENE_COLION(CFLAG,COLTYP,KK)
      ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)

!  PARTICLE TYPE AND SPECIES MIGHT HAVE CHANGED
!  PREPARE POINTER FOR UNIFIED SUBROUTINES
      CALL EIRENE_SWITCH_PARTINFO
C
C  POST-COLLISION ESTIMATOR
C
      IF (LGPART.AND.(NCLVI.GT.0)) THEN
        WS=WEIGHT/SIGTOT
        CALL EIRENE_UPCUSR(WS,2,KK)
        IF (NADSPC_CD >= 1) CALL EIRENE_UPDATE_SPECTRUM (WS,2,1)
      ENDIF
C
      IF (COLTYP.EQ.2) GOTO 700
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
C   NEXT: REFLECTION FROM  SURFACE
C   USE FULL VELOCITY, NOT ONLY THE REDUCED PARALLEL VELOCITY.
C   THIS IS DONE BY SAMPLING THE GYRO-PHASE IN SUBR. NEWFIELD
C   REJECT THOSE GYROPHASES WHICH WOULD LEAD TO NEGATIVE ANGLE OF INCIDENCE
C
C   EXCEPTION: PERIODICITY SURFACE. THEN: NO NEED TO CONVERT TO
C              FULL CARTESIAN VELOCITY COMPONENTS
      IF (ILIIN(MSURF).GE.4) THEN
        PR=1.0
        ICO=0
        IF (.NOT.LGPART) THEN
          WRITE (IUNOUT,*) 'ERROR AT PERIODICITY SURFACE, LGPART=FALSE'
          RETURN
        ENDIF

        IF (NLTRC) THEN
!$OMP CRITICAL
          CALL EIRENE_CHCTRC(X0,Y0,Z0,0,11)
!$OMP END CRITICAL
        ENDIF
        GOTO 1004
      ENDIF
C
      IF (.NOT.LCART) THEN
        NUPC(1)=NPCELL-1+(NTCELL-1)*NP2T3
        NCELL=NRCELL+NUPC(1)*NR1P2+NBLCKA
C  ???
        IF (LDAMCEL(NCELL)) GOTO 9912  ! damaged cell, stop particle

cdr:  try to distinguish: transparent or not. Use arrays "transp(ispz...)"
cdr:  indf=1: transparent, indf=2: non-transparent

        ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
cdr  for solid surface: produce a full cartesian velocity vector, lcart=.true.
        indf=2
cdr  for transparent surface: stick to reduced (GC) velocity, lcart=false
cdr: check here: are any of "transp" flags ne. zero ???
        if (abs(transp(ispz,1,msurf))+abs(transp(ispz,2,msurf)) > 0)
cdr  what about other transparency options: iliin < 0 here ?
cdr  perhaps for those code segment 380 ...ff and call to escape is not reached?
     .     indf = 1
c
c  add gyro velocity (with random phase) to GC velocity:
        ICOUN=0
        DO
!pb       CALL EIRENE_NEWFIELD(X0,Y0,Z0,VELS,2)
          CALL EIRENE_NEWFIELD(X0,Y0,Z0,VELS,indf)
          COSIN=VELX*CRTX+VELY*CRTY+VELZ*CRTZ
C  DOES THE PARTICLE SPEED UNIT VECTOR NOW POINT TOWARDS THE SURFACE ?
          IF (.NOT.LGPART) EXIT  ! DO NOT CARE ABOUT GYRO MOTION, ABSORBED PARTICLE ANYWAY
          IF (ILIIN(MSURF) < 0) EXIT ! DO NOT CARE ABOUT GYRO MOTION, TRANSPARENT SURFACE
          IF (COSIN.GT.0.) EXIT
C  NO, TRY ANOTHER GYRO PHASE
          ICOUN=ICOUN+1
          IF (ICOUN.EQ.100) THEN
            WRITE (IUNOUT,*) 'PARTICLE KILLED AT SURFACE IN FOLION'
            WRITE (IUNOUT,*) 'NO PROPER GYRO ANGLE FOUND'
            WRITE (IUNOUT,*) 'NPANU, MSURF ',NPANU, MSURF
            WRITE (IUNOUT,*) 'VELPER,VELPAR ',VELPER,VELPAR
            LGPART=.FALSE.
            WEIGHT=0.
            ZT=0.0
            GOTO 9951
          ENDIF

        ENDDO
C  NOW A PARTICLE WITH FULL CARTESIAN VELOCITY VECTOR (LCART=T)
C  IS SET. ITS SPEED VECTOR POINTS TOWARDS THE SURFACE (COSIN.GT.0)
      ENDIF
C
C  UPDATE EFFLUXES ONTO SURFACE AND REFLECT PARTICLE
      PR=1.
      IF (ILIIN(MSURF).LE.-2) PR=SG
C
C  FOR NONTRANSPARENT SURFACES:
C  ACCELERATION IN SHEATH IS DONE IN SUBR. ESCAPE
C
!pb   CALL EIRENE_ESCAPE(PR,SG,*100,*104,*996)
      CALL EIRENE_ESCAPE(PR,SG,IRET)
      IF (IRET == 1) GOTO 100
      IF (IRET == 2) GOTO 104
      IF (IRET == 3) GOTO 996
      RETURN
C
C   100: START NEW ION TRACK AFTER SURFACE EVENT
C   104: CONTINUE THIS TRACK, TRANSPARENT SURFACE IS CROSSED
C
C
  700 CONTINUE
C  REGULAR STOP IN SUBR. FOLION, CONTINUE IN SUBR. MCARLO
      RETURN
C
  800 CONTINUE
C  REGULAR STOP IN SUBR. FOLION, STOP HISTORY, CENSUS ARRAY FULL
C     IF (ICOL.EQ.1.AND..NOT.LGLAST) GOTO 512
      LGPART=.FALSE.
      WEIGHT=0.
      RETURN
C
  990 CONTINUE
      CALL EIRENE_LEER(1)
      CALL EIRENE_MASAGE('ERROR IN FOLION, ZDT1 OR NRCELL OUT OF RANGE')
      CALL EIRENE_MASAGE('PARTICLE IS KILLED')
      WRITE (iunout,*) 'NPANU,NRCELL,ZDT1,ZTST ',NPANU,NRCELL,ZDT1,ZTST
      WRITE (iunout,*) 'TL,TS,ZINT1,ZLOG ',TL,TS,ZINT1,ZLOG
      GOTO 995
  991 CONTINUE
      CALL EIRENE_LEER(1)
      CALL EIRENE_MASAGE('ERROR IN FOLION, NCELL OUT OF RANGE')
      CALL EIRENE_MASAGE('PARTICLE IS KILLED')
      WRITE (iunout,*) 'NPANU,NCELL,NRCELL,NPCELL,NTCELL '
      WRITE (iunout,*)  NPANU,NCELL,NRCELL,NPCELL,NTCELL
      GOTO 995
 9911 CONTINUE
      CALL EIRENE_LEER(1)
      CALL EIRENE_MASAGE('ERROR IN FOLION, NO INTERSECTION FOUND')
      CALL EIRENE_MASAGE('PARTICLE IS KILLED')
      WRITE (iunout,*) 'NPANU,NCELL,NRCELL,NPCELL,NTCELL '
      WRITE (iunout,*)  NPANU,NCELL,NRCELL,NPCELL,NTCELL
      GOTO 995
C
 9912 CONTINUE
      CALL EIRENE_LEER(1)
      CALL EIRENE_MASAGE('ERROR IN FOLION, DAMAGED CELL HIT')
      CALL EIRENE_MASAGE('PARTICLE IS KILLED')
      WRITE (iunout,*) 'NPANU,NCELL,NRCELL,NPCELL,NTCELL '
      WRITE (iunout,*)  NPANU,NCELL,NRCELL,NPCELL,NTCELL
      GOTO 995
C
  992 CONTINUE
      CALL EIRENE_LEER(1)
      CALL EIRENE_MASAGE('ERROR IN FOLION, PROJECTION TO V_PAR, V_PERP')
      CALL EIRENE_MASAGE
     .  ('PROBABLY ILL-DEFINED B FIELD WRT. PARTICLE SPEED')
      WRITE (iunout,*) 'BBX,BBY,BBZ ',BBX,BBY,BBZ
      ZT=0.
      GOTO 9951

 9921 CONTINUE
      CALL EIRENE_LEER(1)
      CALL EIRENE_MASAGE('ERROR IN FOLION, LCART HAS WRONG VALUE')
      WRITE (IUNOUT,*) 'NPANU,LCART ',NPANU,LCART
      IF (NLTRC) THEN
!$OMP CRITICAL
        CALL EIRENE_CHCTRC(X0,Y0,Z0,16,18)
!$OMP END CRITICAL
      ENDIF
      GOTO 999
C
  994 CALL EIRENE_MASAGE('ERROR IN FOLION, AT SURFACE DELTA EVENT')
      WRITE (iunout,*) 'IION,NPANU ',IION,NPANU
      GOTO 999
C
  995 WRITE (iunout,*) 'MRSURF,MPSURF,MTSURF,MASURF ',
     .                  MRSURF,MPSURF,MTSURF,MASURF
 9951 X0ERR=X0+ZT*VELX
      Y0ERR=Y0+ZT*VELY
      Z0ERR=Z0+ZT*VELZ
      IF (NLTRC) THEN
!$OMP CRITICAL
        CALL EIRENE_CHCTRC(X0ERR,Y0ERR,Z0ERR,16,18)
!$OMP END CRITICAL
      ELSE
!$OMP CRITICAL
        WRITE (iunout,'(A,1P,4(1X,1E14.7))') 'X0,Y0,Z0,ZT ',X0,Y0,Z0,ZT
        WRITE (iunout,'(A,1P,4(1X,1E14.7))') 'VELX,VELY,VELZ,VEL ',
     .                                        VELX,VELY,VELZ,VEL
        WRITE (iunout,'(A,1P,3(1X,1E14.7))') 'X0ERR,Y0ERR,Z0ERR ',
     .                                        X0ERR,Y0ERR,Z0ERR
!$OMP END CRITICAL
      ENDIF
      GOTO 999
  996 CALL EIRENE_MASAGE
     .  ('ERROR IN FOLION, COND. EXP. ESTIM. NOT IN USE')
      GOTO 999
  997 CALL EIRENE_MASAGE('ERROR IN FOLION, DETECTED IN SUBR. CLLTST')
      CALL EIRENE_MASAGE('PARTICLE IS KILLED')
C   DETAILED PRINTOUT ALREADY DONE FROM SUBR. CLLTST

      IF (NLTRC) THEN
!$OMP CRITICAL
        CALL EIRENE_CHCTRC(X0,Y0,Z0,16,18)
!$OMP END CRITICAL
      ENDIF
      GOTO 999
C
  998 CALL EIRENE_MASAGE('ERROR IN FOLION, SPECIES INDEX OUT OF RANGE')
      WRITE (iunout,*) ' NPANU,IION ',NPANU,IION
      GOTO 999
C
  999 CONTINUE
!$OMP ATOMIC
      PTRASH(ISTRA)=PTRASH(ISTRA)-WEIGHT
!$OMP ATOMIC
      ETRASH(ISTRA)=ETRASH(ISTRA)-WEIGHT*E0
      LGPART=.FALSE.
      WEIGHT=0.
 9991 CALL EIRENE_LEER(1)
      RETURN

      CONTAINS
C  ION-ION ENERGY LOSS FREQUENCY (LANGER APPROXIMATION) (1/SEC)
C  NUCL.FUS. 22, NO. 6, (1986) P754, FOR CH4+ (mA=16) ON H+ (mB=1)
      FUNCTION FNUEQI(XNI,TI)
      REAL(DP) ::  FNUEQI,XNI,TI
c     FNUEQI=8.8E-8*XNI*TI**(-1.5)
c  This is not exactly the relaxation time, but instead a time
c  which appears in the analytical (BGK-like) solution EA(t).
c  to obtain an effective  nu(Ti) that can be compared with a
c  "relaxation time"
c  in the dgl dEA/dt=-nu(Ti,EA,...) times EA
c  In the present limit: this must be multiplied by a factor(EA,Ti)
      FNUEQI=8.5E-8*XNI*TI**(-1.5)
c  in calling program: FNUEQI = FNUEQI*(1.+mB/mA)**0.5-1.5*Ti/EA
c  but this is already implicitly contained in the analytic BGK solution
c  written for fnueqi without that factor.
      RETURN
      END FUNCTION FNUEQI

C  ION-ION ENERGY LOSS FREQUENCY (LOW ENERGY LIMIT, NRL) (1/SEC)
C  GENERALIZATION OF LANGER EXPRESSION TO ARBITRARY IONS (MASS, CHARGE)
C  note: for an intermediate period (1995 --2013) the mass factor
c  (1+mB/mA) had an incorrect exponent -1/2, in the NRL formularies.
c  2016: back to the correct formula (as in eighties) without that exponent

      FUNCTION FNUEQI_1(EA,XNI,TI,ION,IPL)
      REAL(DP) ::  FNUEQI_1,EA,XNI,TI
      INTEGER ::  ION,IPL
      REAL(DP) ::  Coullog,fact,za,zb,XMUA,XMUB
      Coullog=10.
      ZA=NCHRGI(ION)
      ZB=NCHRGP(IPL)
      XMUA=nMASSI(ION)
      XMUB=nMASSP(IPL)
      FACT=XNI*ZA**2*ZB**2*COULLOG*6.8E-8*XMUB**0.5/XMUA/TI**0.5
      FNUEQI_1=FACT*(2./TI*(1.+XMUB/XMUA)-2/EA-1/EA)
      RETURN
      END FUNCTION FNUEQI_1

C  ION-ION ENERGY LOSS FREQUENCY (FULL EXPRESSION, NRL) (1/SEC)
C  INVOLVING THE CHANDRASEKHAR FUNCTIONS

      FUNCTION FNUEQI_2(EA,XNI,TI,ION,IPL)
      REAL(DP) ::  FNUEQI_2,EA,XNI,TI
      INTEGER ::  ION,IPL
      REAL(DP) ::  Coullog,XNUE0,za,zb,XMUA,XMUB,XAB,
     .             vela,xma,xmb,eza,ezb
      COULLOG=10.
      ZA=NCHRGI(ION)
      ZB=NCHRGP(IPL)
      eZA=ZA*4.8032d-10  ! charge in statcoul (cgs)
      eZB=ZB*4.8032d-10
      XMUA=nMASSI(ION)
      XMUB=nMASSP(IPL)
      XMA=xmua*amua   ! in g
      XMB=xmub*amua   ! in g
      VELA=CVELAA*SQRT(EA/XMUA) ! cm/s
      XNUE0=XNI*eZA**2*eZB**2*COULLOG*4.*PIA/XMA**2/VELA**3
      XAB=XMUB/(2.*TI)*EA/XMUA*2  ! DIMENSIONLESS
      XAB=XMB/(2.*TI*1.6e-12)*vela*vela  ! DIMENSIONLESS
      FNUEQI_2=2.*XNUE0*(XMUA/XMUB*PSI_CHAND(XAB)-DPSI_CHAND(XAB))
      RETURN
      END FUNCTION FNUEQI_2

      FUNCTION PSI_CHAND(X)
      REAL(DP) :: PSI_CHAND,X
#ifdef PGF90
      REAL(DP) :: DERF
      PSI_CHAND=-DERF(SQRT(X))+2./SQRT(PIA)*EXP(-X)*SQRT(X)
#else
      PSI_CHAND=-ERF(SQRT(X))+2./SQRT(PIA)*EXP(-X)*SQRT(X)
#endif
      RETURN
      END FUNCTION PSI_CHAND

      FUNCTION DPSI_CHAND(X)
      REAL(DP) :: DPSI_CHAND,X
      DPSI_CHAND=2./SQRT(PIA)*EXP(-X)*SQRT(X)
      RETURN
      END FUNCTION DPSI_CHAND

      END SUBROUTINE EIRENE_FOLION

      SUBROUTINE EIRENE_NEWFIELD(X,Y,Z,VELS,IND)
C  FIND NEW MAGNETIC FIELD AT NEW POINT X,Y,Z IN CELL NCELL
C  IF (IND.EQ.0) RETURN WITH NEW LOCAL B FIELD BVEC
C
C  IF (IND.GE.1) ADDITIONALLY ALSO PROVIDE REDUCED (GC) VELOCITY VECTOR (SPEED UNIT VECTOR)
C    BUT RETAIN PREVIOUS MODULI: V_PARALLEL, V_PERP.
C    NEW REDUCED SPEED VECTOR:  LCART=FALSE AND VELX,VELY,VELY, SPEED: VEL (=VELPAR),
C    CHECKS DONE THAT VELPER AND VERPAR ARE PRESERVED, CHECKS REMOVED.

C  IF (IND.GE.2) ADDITIONALLY ALSO PROVIDE NEW CARTESIAN VELOCITY
C  BY SAMPLING THE GYRO PHASE, AND A COORDINATE TRANSFORMATION IN
C  VEL-SPACE.
C  NEW CARTESIAN VELOCITY VECTOR: LCART=.TRUE., VELX,VELZ,VELZ, VEL
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CCONA
      USE EIRMOD_CFPLK
      USE EIRMOD_COMPRT
      USE EIRMOD_CRAND
      USE EIRMOD_CINIT
      USE EIRMOD_RANF, ONLY: RANF_EIRENE
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: X,Y,Z,VELS
      REAL(DP) :: BVEC_1(3), VVEC(3), GYRO, BBF
      INTEGER :: IND

      CALL EIRENE_BFIELD (NCELL, X, Y, Z, BBX, BBY, BBZ, BBF,.TRUE.)
      BVEC = (/ BBX, BBY, BBZ /)

      IF (IND.LT.1) RETURN

C  FIND NEW REDUCED (GUIDING CENTRE) VELOCITY, LCART=F
C  RETAIN MODULI VEL, V_PARALLEL, V_PERP, SIGPAR,
C  ONLY THE NEW DIRECTION (REDUCED SPEED UNIT VECTORS) ARE EVALUATED
      VLXPAR=SIGPAR*BBX
      VLYPAR=SIGPAR*BBY
      VLZPAR=SIGPAR*BBZ
      VELX = VLXPAR
      VELY = VLYPAR
      VELZ = VLZPAR
      VEL  = VELPAR
      LCART=.FALSE.

      IF (IND.LT.2) RETURN

C  FIND NEW CARTESIAN VELX,VELY,VELZ (SAME VEL=VELS), LCART=T
C  NEW GYRO PHASE
      GYRO=RANF_EIRENE()*PI2A
C  BACK TO CARTESIAN COORDINATES
      CALL EIRENE_B_PROJI (BVEC,BVEC_1,VVEC,SIGPAR*VELPAR,VELPER,GYRO)
      VELX = VVEC(1)
      VELY = VVEC(2)
      VELZ = VVEC(3)
      VEL  = VELS
      LCART=.TRUE.
      RETURN
      END SUBROUTINE EIRENE_NEWFIELD
