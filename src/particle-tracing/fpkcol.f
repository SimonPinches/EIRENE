C  sept. 05: use only parallel velocity, vel=velpar,....
C              this is now made consistent in calling program folion.f
C  aug 06:   move call to chctrc: pre-collision status at point of collision
C            and: return new velocity vector in full cartesian coord.
C               fetch new BVEC at point of collision
C  mai 10:   flag ind:  if ind=2, only FP collision, but no push to
C                               new position.
c  july 15:  set E0PAR,  (was missing).
cdr nov. 15:  multiple bulk ion species, new array fnuiar(ipl)
cdr           to be done:  proper definition of eipl, and e0new, in cases
cdr                        of multiple background ion species 
C
      SUBROUTINE EIRENE_FPKCOL(*,*,*,IND)
C
C  IF IND=0 (DEFAULT)
C  1.) ADVANCE PARTICLE BY TIMESTEP DUR, AND THEN
C  2.) CARRY OUT FOKKER PLANCK ELASTIC COLLISION
C  3.) TURN VELOCITY INTO FULL CARTESIAN PARTICLE VELOCITY

C  IF IND=1, SKIP STEP 1, ONLY FOKKER PLANCK COLLISION AND FULL
C                         VELOCITY, 
C                         LCART=T
C  IF IND=2, SKIP STEP 3, ONLY PUSH AND FOKKER PLANCK COLLISION,
C                         BUT STAY IN REDUCED (GUIDING CENTRE) VELOCITY
C                         LCART=F
C  IF IND=3, SKIP STEPS 1,3, ONLY FOKKER PLANCK COLLISION,
C                        BUT STAY IN REDUCED (GUIDING CENTRE) VELOCITY
C                        LCART=F

C  ON INPUT: LCART=FALSE, ENTER WITH REDUCED (GUIDING CENTRE) VELOCITY
C            ZT = DISTANCE ALONG REDUCED (GUIDING CENTRE) PATH
C            VELX,VELY,VELZ: REDUCED (GUIDING CENTRE) SPEED UNIT VECTOR
C
C  RETURN  :  NOT IN USE
C  RETURN 1:  NOT IN USE
C  RETURN 2:  START COMPLETELY NEW TEST ION TRACK, SAME SPECIES
C             LCART=TRUE
C  RETURN 3:  ERROR, STOP TRACK IN CALLING PROGRAM. PTRASH AND ETRASH ALREADY DONE HERE
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CINIT
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_CCONA
      USE EIRMOD_CFPLK
      USE EIRMOD_CLOGAU
      USE EIRMOD_CUPD
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      USE EIRMOD_CZT1
      USE EIRMOD_COMPRT
      USE EIRMOD_CLGIN
      USE EIRMOD_COUTAU
      USE EIRMOD_COMXS
      USE EIRMOD_CVARUSR
 
      IMPLICIT NONE
 
      REAL(DP) :: DUR, E0OLD, E0NEW, VNEW, WS, FAC, GYRO,
     .            BVEC_1(3), VVEC(3), VELS, FNUI, EWG,
cdr  old code was:
c    .            D_VELPAR(1:NPLSI),
c    .            D_VELPER(1:NPLSI),
c    .            D_E0NEW_tmp(1:NPLSI)
cdr: I beliefe this is now wrong  (old code correct ?) 
cdr  argument IPL in these arrays is in the range 1:NPLS, 
cdr  because LGIEL(..,..1) is.
     .            D_VELPAR(1:NIELI(IION)),
     .            D_VELPER(1:NIELI(IION)),
     .            D_E0NEW_tmp(1:NIELI(IION))

      INTEGER :: IOLD, EIRENE_LEARC2, NCELLT, IND, IPL, IPLTI, IDSC
      REAL(DP), EXTERNAL :: RANF_EIRENE

C  SAVE INCIDENT SPECIES: IOLD
      IOLD=IION
      E0OLD=E0
      NCELLT=NCLTAL(NCELL)
      DUR=ZT/VEL   ! TIMESTEP UNTIL THIS FP COLLISION, [S]

C
      IF (LCART) GOTO 991

C  SKIP PUSH ?
      IF (IND.EQ.1.OR.IND.EQ.3) GOTO 200

C  1.) PUSH TO NEW POSITION ALONG REDUCED (GUIDING CENTRE) TRACK
C     write (6,*) 'fpkcol push, zt used', zt
      X0=X0+VELX*ZT
      Y0=Y0+VELY*ZT
      Z0=Z0+VELZ*ZT
      TIME=TIME+DUR
      IF (LEVGEO.LE.3.AND.NLPOL) THEN
        IPOLG=NPCELL
      ELSEIF (NLPLG) THEN
        IPOLG=EIRENE_LEARC2(X0,Y0,NRCELL,NPANU,'FOLION 2     ')
      ELSEIF (NLFEM) THEN
        IPOLG=0
      ENDIF
      NLSRFX=.FALSE.
      NLSRFY=.FALSE.
      NLSRFZ=.FALSE.
      MRSURF=0
      MPSURF=0
      MTSURF=0
      MASURF=0
      MSURF=0
      IF (NLTRA) PHI=MOD(PHI-ATAN2(Z01,X01)+ATAN2(Z0,(RMTOR+X0)),PI2A)
      

C  TEST FOR CORRECT CELL NUMBER AT COLLISION POINT
C  KILL PARTICLE, IF TOO LARGE ROUND OFF ERRORS DURING
C  PARTICLE TRACING
C
      IF (NLTEST) CALL EIRENE_CLLTST(*997)

C
C  FROM THIS POINT: FOKKER PLANCK COLLISION MODEL

200   CONTINUE

C
C
C  PRE COLLISION ESTIMATOR
C
      IF (NCLVI.GT.0) THEN
        WS=WEIGHT/SIGTOT
        CALL EIRENE_UPCUSR(WS,1)
      ENDIF
C
      IF (NLTRC) CALL EIRENE_CHCTRC(X0,Y0,Z0,16,7)
C
      IF (DUR.GT.0.D0) THEN

c  second minimal collision model: change velocities according to
c  the expectation values of the change
c  SIGPAR is the sign of the parallel velocity with respect to the
c  magnetic field

        D_VELPAR    = 0.0
        D_VELPER    = 0.0
        D_E0NEW_tmp = 0.0

        DO IDSC = 1, NIELI(IION) ! loop over number of elastic collisions
          IPL = LGIEL(IION,IDSC,1)
          D_VELPAR(IPL) = dVelPrl_dt(IPL)*DUR   ! in cm/s
          D_VELPER(IPL) = dVelPerp_dt(IPL)*DUR  ! in cm/s

          D_E0NEW_tmp(IPL)= CVELI2*RMASSI(IION)*
     >                      ((SIGPAR*VELPAR + D_VELPAR(IPL))**2
     >                    + (VELPER + D_VELPER(IPL))**2)
     >                    - E0
        END DO

        VELPAR = SIGPAR*VELPAR + SUM(D_VELPAR)
        SIGPAR = SIGN(1.0,VELPAR)
        VELPAR = ABS(VELPAR)
        VELPER = VELPER + SUM(D_VELPER)
        veltotal = SQRT(VELPAR**2 + VELPER**2)

        E0NEW = CVELI2*RMASSI(IION)*veltotal**2
        VNEW = RSQDVI(IOLD)*SQRT(E0NEW)
C
C  UPDATE ESTIMATORS EIIO,EIPL
        EWG = WEIGHT*(E0NEW-E0OLD)
        EIIO(NCELLT) = EIIO(NCELLT)+EWG
        DO IDSC = 1, NIELI(IION) ! loop over number of elastic collisions
          IPL = LGIEL(IION,IDSC,1)
          IF (D_E0NEW_tmp(IPL).NE.0.0) THEN
            EIPL(IPL,NCELLT) = EIPL(IPL,NCELLT)
     >                   - EWG*D_E0NEW_tmp(IPL)/(SUM(D_E0NEW_tmp)+EPS60)
          END IF
        END DO

        FAC    = SQRT(E0NEW/E0OLD)
        E0PAR  = E0PAR*FAC*FAC ! ratio of the particle velocity before and after the collision

      ENDIF
C  FP COLLISION DONE, LCART=F STILL, I.E. VEL = V_GC
c  gets new B-field

!pb VELS is not used in NEWFIELD with option 1
!pb but for the sake of decent programming set VELS
      VELS = VEL
      CALL EIRENE_NEWFIELD(X0,Y0,Z0,VELS,1)

C  SKIP TRANSFORM TO FULL VELOCITY AND RETURN WITH LCART=F  ?
c  FHa IND.EQ.2 or IND.EQ.3 means that the velocity is updated in
c  EIRENE_NEWFIELD
      IF (IND.EQ.2.OR.IND.EQ.3) GOTO 300

C  RETURN WITH FULL CARTESIAN VELOCITY VECTOR V = V_FULL
 
C  NEW B-FIELD
!pb VELS is not used in NEWFIELD with option 0
!pb but for the sake of decent programming set VELS
      VELS = VEL
      CALL EIRENE_NEWFIELD(X0,Y0,Z0,VELS,0)

C  NEW GYRO PHASE
      GYRO=RANF_EIRENE()*PI2A
C  BACK TO CARTESIAN COORDIANTES
      CALL EIRENE_B_PROJI
     .                 (BVEC,BVEC_1,VVEC,SIGPAR*VELPAR,VELPER,GYRO)
      VELX = VVEC(1)
      VELY = VVEC(2)
      VELZ = VVEC(3)
      LCART=.TRUE.

300   CONTINUE

      VEL=VNEW
      E0=E0NEW
      RETURN 2
C
C  POST COLLISION ESTIMATOR
C
C     IF (NCLVI.GT.0) THEN
C       WS=WEIGHT/SIGTOT
C       CALL UPCUSR(WS,2)
C     ENDIF
      RETURN 2
C
991   CALL EIRENE_MASAGE
     .  ('ERROR IN FPKCOL,  CALLED WITH LCART=TRUE  ')
      CALL EIRENE_MASAGE
     .  ('PARTICLE IS KILLED                            ')
      GOTO 999
997   CALL EIRENE_MASAGE
     .  ('ERROR IN FPKCOL,  DETECTED IN SUBR. CLLTST    ')
      CALL EIRENE_MASAGE
     .  ('PARTICLE IS KILLED                            ')
C   DETAILED PRINTOUT ALREADY DONE FROM SUBR. CLLTST
      IF (NLTRC) CALL EIRENE_CHCTRC(X0,Y0,Z0,16,18)
      GOTO 999
C
999   PTRASH(ISTRA)=PTRASH(ISTRA)-WEIGHT
      ETRASH(ISTRA)=ETRASH(ISTRA)-WEIGHT*E0
      LGPART=.FALSE.
      WEIGHT=0.
      CALL EIRENE_LEER(1)
      RETURN 3
      END
