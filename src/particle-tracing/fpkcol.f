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
!pb   SUBROUTINE EIRENE_FPKCOL(*,*,*,IND)
      SUBROUTINE EIRENE_FPKCOL(IRET,IND)
C
C  IF IND=0 (DEFAULT)
C  1.) ADVANCE PARTICLE BY TIMESTEP DUR, AND THEN
C  2.) CARRY OUT FOKKER-PLANCK ELASTIC COLLISION
C  3.) TURN VELOCITY INTO FULL CARTESIAN PARTICLE VELOCITY

C  IF IND=1, SKIP STEP 1, ONLY FOKKER-PLANCK COLLISION AND FULL
C                         VELOCITY,
C                         LCART=T
C  IF IND=2, SKIP STEP 3, ONLY PUSH AND FOKKER-PLANCK COLLISION,
C                         BUT STAY IN REDUCED (GUIDING CENTRE) VELOCITY
C                         LCART=F
C  IF IND=3, SKIP STEPS 1,3, ONLY FOKKER-PLANCK COLLISION,
C                        BUT STAY IN REDUCED (GUIDING CENTRE) VELOCITY
C                        LCART=F

C  ON INPUT: LCART=FALSE, ENTER WITH REDUCED (GUIDING CENTRE) VELOCITY
C            ZT = DISTANCE ALONG REDUCED (GUIDING CENTRE) PATH
C            VELX,VELY,VELZ: REDUCED (GUIDING CENTRE) SPEED UNIT VECTOR
C
C  RETURN  :  NOT IN USE
C  RETURN IRET=1:  NOT IN USE
C  RETURN IRET=2:  START COMPLETELY NEW TEST ION TRACK, SAME SPECIES
C             LCART=TRUE
C  RETURN IRET=3:  ERROR, STOP TRACK IN CALLING PROGRAM. PTRASH AND ETRASH ALREADY DONE HERE
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
      USE EIRMOD_RANF, ONLY: RANF_EIRENE
      USE EIRMOD_PLT2D, ONLY: EIRENE_CHCTRC
      USE EIRMOD_CSDVI, ONLY: LMETSP

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: IND
      INTEGER, INTENT(OUT) :: IRET
      
      REAL(DP) :: DUR, E0OLD, E0NEW, VNEW, WS, FAC, GYRO,
     .            BVEC_1(3), VVEC(3), VELS, FNUI, EWG,
     .            BX, BY, BZ, BF, DIRPROJ, VDEL, VPLASP, SIG,
     .            VX, VY, VZ
      INTEGER :: IOLD, EIRENE_LEARC2, NCELLT, IPL, IPLSV, IPLSLOC
ctk      REAL(DP), EXTERNAL :: RANF_EIRENE
C     SAVE INCIDENT SPECIES: IOLD
      
      IRET = 0
      IOLD=IION
      E0OLD=E0
      NCELLT=NCLTAL(NCELL)
      DUR=ZT/VEL   ! TIMESTEP UNTIL THIS FP COLLISION, [S]

C
      IF (LCART) GOTO 991

C  SKIP PUSH ?
      IF ((IND.EQ.1).OR.(IND.EQ.3)) GOTO 200

C  1.) PUSH TO NEW POSITION ALONG REDUCED (GUIDING CENTRE) TRACK
C     write (iunout,*) 'fpkcol push, zt used', zt
      X0=X0+VELX*ZT
      Y0=Y0+VELY*ZT
      Z0=Z0+VELZ*ZT
      TIME=TIME+DUR
      IF ((LEVGEO.LE.3).AND.NLPOL) THEN
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

      IF (NLTRC) THEN
!$OMP CRITICAL
        CALL EIRENE_CHCTRC(X0,Y0,Z0,16,7)
!$OMP END CRITICAL
      ENDIF

C  TEST FOR CORRECT CELL NUMBER AT COLLISION POINT
C  KILL PARTICLE, IF TOO LARGE ROUND-OFF ERRORS DURING
C  PARTICLE TRACING
C
      IF (NLTEST) THEN
        CALL EIRENE_CLLTST(IRET)
        IF (IRET.EQ.1) GOTO 997
      ENDIF

C
C  FROM THIS POINT: FOKKER-PLANCK COLLISION MODEL

  200 CONTINUE

C
C
C  PRE-COLLISION ESTIMATOR
C
      IF (NCLVI.GT.0) THEN
        IF (SIGTOT.GT.0) THEN ! SIGTOT can be zero here (probably empty cell, so no reactions, but in that case we should not even be here, Fokker-Planck should be deactivated as well.)
        WS=WEIGHT/SIGTOT
        ELSE
          WS=0._DP
        ENDIF
CNR     NREACI+1 : REACTION INDEX FOR FOKKER-PLANCK (NOT IN REACTION LIST IN INPUT FILE)
        CALL EIRENE_UPCUSR(WS,1,NREACI+1)
      ENDIF
C
C
      IF (DUR.GT.0.D0) THEN
C
C  FLIGHT WITH PARALLEL VELOCITY VEL=VELPAR (CM/SEC)
C  PARALLEL DISTANCE ZT (CM)
C  ENERGY RELAXATION CONSTANT TAUE
C
cdr to be done: proper new energy, according to weighting by fnuiar(ipl)
cdr relaxation towards a weighted mean background energy
cdr currently: arbitrary 1.5*Tiin(1,...)
        E0NEW=E0OLD*EXP(-DUR/TAUE)+1.5*TIIN(1,NCELL)*(1.-EXP(-DUR/TAUE))
        VNEW=RSQDVI(IOLD)*SQRT(E0NEW)
C
C  UPDATE ESTIMATORS EIIO,EIPL

!$OMP ATOMIC
        EIIO(NCELLT)=EIIO(NCELLT)+WEIGHT*(E0NEW-E0OLD)

cdr  for the time being: distribute bulk ion energy loss proportional to collision frequency
cdr  strictly bulk ipls1 and ipls2 can have different gains/losses, depending on their
cdr  temprature(ipls), even different sign.
cdr
        EWG = WEIGHT*(E0NEW-E0OLD)
        FNUI = SUM(FNUIAR(1:NPLSI))  ! CDR THIS SUM SHOULD BE KNOWN FROM CALLING ROUTINE
        DO IPL = 1, NPLSI
!$OMP ATOMIC
          EIPL(IPL,NCELLT)=EIPL(IPL,NCELLT)-EWG*FNUIAR(IPL)/FNUI
        END DO
C

        FAC=SQRT(E0NEW/E0OLD)
        IF (NLSOLEDGE) THEN
        
CNR Now update the parallel momentum tallies (*** WARNING: Only for main background ion ***)
        IF (LMIPL) THEN
C  SET THE POST-COLLISION TEST PARTICLE PARALLEL VELOCITY
          IPLSLOC = 1 ! Hardcoded only for main background ion
          CALL EIRENE_BFIELD (NCELL,X0,Y0,Z0,BX,BY,BZ,BF,.TRUE.)
          DIRPROJ = VELX*BX+VELY*BY+VELZ*BZ ! The test ion does not change (parallel) direction during the collision
          VDEL=VELPAR*(1._DP-FAC)*DIRPROJ*AMUA*RMASSI(IION)*WEIGHT ! VELPAR_old-VNEW_parallel
C
          IF (INDPRO(4) == 8) THEN
            CALL EIRENE_VECUSR(2,NCELL,X0,Y0,Z0,VX,VY,VZ,IPLSLOC,
     .                         .TRUE.)
            VPLASP=VX*BX+VY*BY+VZ*BZ
            SIG=SIGN(1._DP,VPLASP)
          ELSE
            IPLSV=MPLSV(IPLSLOC)
            SIG=1._DP
            IF (LBVIN) SIG=SIGN(1._DP,BVIN(IPLSV,NCELL))
          ENDIF
          MIPL(IPLSLOC,NCELL)=MIPL(IPLSLOC,NCELL)+VDEL*SIG
          LMETSP(NSPAMI+IPLSLOC)=.TRUE.
        ENDIF
        END IF
C in this particlar case: retain old pitch: velpar/velper. No pitch angle scattering so far.
        VELPAR=VELPAR*FAC
        VELPER=VELPER*FAC
        E0PAR=E0PAR*FAC*FAC

      ELSE  ! flight time DUR IS NOT GT 0.0
        GOTO 998
      ENDIF
C  FP COLLISION DONE, LCART=F STILL, I.E. VEL = V_GC
c  gets new B field

!pb VELS is not used in NEWFIELD with option 1
!pb but for the sake of decent programming set VELS
      VELS = VEL
CNR   ONLY UPDATE VELOCITY VECTOR IF THE PARTICLE HAS BEEN PUSHED, NO NEED
CNR   OTHERWISE, AND BREAKS TRAJECTORY IN CASE OF INT. GRID SURFACE: THE
CNR   VELOCITY NEEDS TO BE KEPT FROM THE PREVIOUS CELL, OR THE FACE/VELOCITY
CNR   INTERSECTION IN TIMER MAY NOT BE FOUND.
      IF (IND.EQ.0.OR.IND.EQ.2) THEN
      CALL EIRENE_NEWFIELD(X0,Y0,Z0,VELS,1)
      ENDIF

C  SKIP TRANSFORM TO FULL VELOCITY AND RETURN WITH LCART=F  ?

      IF (IND.EQ.2.OR.IND.EQ.3) GOTO 300

C  RETURN WITH FULL CARTESIAN VELOCITY VECTOR V = V_FULL

C  NEW B FIELD

!pb VELS is not used in NEWFIELD with option 0
!pb but for the sake of decent programming set VELS
      VELS = VEL
      CALL EIRENE_NEWFIELD(X0,Y0,Z0,VELS,0)

C  NEW GYRO PHASE
      GYRO=RANF_EIRENE()*PI2A
C  BACK TO CARTESIAN COORDINATES
      CALL EIRENE_B_PROJI
     .                 (BVEC,BVEC_1,VVEC,SIGPAR*VELPAR,VELPER,GYRO)
      VELX = VVEC(1)
      VELY = VVEC(2)
      VELZ = VVEC(3)
      LCART=.TRUE.
c strictly: e0new, vnew should be modified, due to new gyro phase.

  300 CONTINUE

      VEL=VNEW
      E0=E0NEW
      IRET = 2
CNR   RETURN
C
C  POST-COLLISION ESTIMATOR
C
      IF (NCLVI.GT.0) THEN
        IF (SIGTOT.GT.0) THEN ! SIGTOT can be zero here (probably empty cell, so no reactions, but in that case we should not even be here, Fokker-Planck should be deactivated as well.)
          WS=WEIGHT/SIGTOT
        ELSE
          WS=0._DP
        ENDIF
CNR     NREACI+1 : REACTION INDEX FOR FOKKER-PLANCK (NOT IN REACTION LIST IN INPUT FILE)
        CALL EIRENE_UPCUSR(WS,2,NREACI+1)
      ENDIF
      IRET = 2
      RETURN
C
  991 CALL EIRENE_MASAGE
     .  ('ERROR IN FPKCOL, CALLED WITH LCART=TRUE')
      CALL EIRENE_MASAGE('PARTICLE IS KILLED')
      GOTO 999
  997 CALL EIRENE_MASAGE
     .  ('ERROR IN FPKCOL, DETECTED IN SUBR. CLLTST')
      CALL EIRENE_MASAGE('PARTICLE IS KILLED')
C   DETAILED PRINTOUT ALREADY DONE FROM SUBR. CLLTST
      IF (NLTRC) THEN
!$OMP CRITICAL
        CALL EIRENE_CHCTRC(X0,Y0,Z0,16,18)
!$OMP END CRITICAL
      ENDIF
      GOTO 999
998   CALL EIRENE_MASAGE
     .  ('ERROR IN FPKCOL, NEG. FLIGHT TIME ENCOUNTERED')
      CALL EIRENE_MASAGE('PARTICLE IS KILLED')
C   DETAILED PRINTOUT ALREADY DONE FROM SUBR. CLLTST
      WRITE (IUNOUT,'(A,I6)') ' NPANU = ',NPANU
      WRITE (IUNOUT,'(A,2ES12.4)') ' ZT, VEL = ', ZT, VEL
      IF (NLTRC) THEN
!$OMP CRITICAL
        CALL EIRENE_CHCTRC(X0,Y0,Z0,16,18)
!$OMP END CRITICAL
      ENDIF 
      GOTO 999
C
  999 LGPART=.FALSE. 
!$OMP ATOMIC 
      PTRASH(ISTRA)=PTRASH(ISTRA)-WEIGHT
!$OMP ATOMIC
      ETRASH(ISTRA)=ETRASH(ISTRA)-WEIGHT*E0
      LGPART=.FALSE.
      WEIGHT=0.
      CALL EIRENE_LEER(1)
      IRET = 3
      RETURN
      END SUBROUTINE EIRENE_FPKCOL
