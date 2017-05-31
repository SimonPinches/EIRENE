cdr:  njump=3 option: continue new flight from previously found
cdr                   intersection with cell boundary. used from levgeo=10
C icts introduced, lcut --> llcut  (sometime after 2014)
!pb  22.03.07:  LEVGEO=6 --> LEVGEO=10
cdr   may 17:  commenting. tbd: commenting in particular of tetraheda-grid part.
C
C  FULL EIRENE GEOMETRY BLOCK  (GEO3D)
C
C
      SUBROUTINE EIRENE_TIMER (PT)
C
C  THIS SUBROUTINE CALCULATES DISTANCES TO INTERSECTION POINTS IN THE 1ST STANDARD
C  MESH "RSURF" (X- OR RADIAL DIRECTION)
C
C  INPUT:
C       NRCELL = CELL NUMBER FOR WHICH NEXT INTERSECTION
C                TIME IS TO BE CALCULATED (I.E. NOT NECESSARLY THE
C                CELL WHICH CONTAINS THE STARTING POINT X0,Y0,Z0)
C       NRCELL=0 IF PARTICLE OUTSIDE STANDARD MESH
C
C       NJUMP = 0 MEANS: THIS IS THE FIRST CALL OF TIMER FOR THIS TRACK
C                        IN THIS CASE , NLSRFX MUST BE KNOWN
C           NLSRFX = .TRUE. :PARTICLE ON A SURFACE, IN THIS CASE
C                            THE SUBROUTINE NEEDS 'MRSURF'
C                            MRSURF = NUMBER OF THIS SURFACE
C
C           NLSRFX = .FALSE.:PARTICLE NOT ON A SURFACE
C
C           X0,Y0,Z0 = STARTING POINT OF THIS TRACK
C           VELX,VELY,VELZ = SPEED UNIT VECTOR OF PARTICLE
C       NJUMP = 1  X0,Y0,Z0,VELX,VELY,VELZ ARE THE SAME
C                  AS IN THE PREVIOUS CALL
C       NJUMP = 2  ONLY VELX,VELY,VELZ ARE THE SAME
C                  AS IN THE PREVIOUS CALL, I.E. PARTICLE HAS
C                  BEEN MOVED BUT VELOCITY HAS NOT BEEN CHANGED
C                  (TO BE WRITTEN)
C       NJUMP = 3  PARTICLE HAS
C                  BEEN MOVED, VELOCITY HAS CHANGED, BUT PARTICLE CONTINUES FLIGHT
C                  FROM PREVIOUSLY FOUND CELL BOUNDARY INTERSECTION POINT
C  IF (LEVGEO=1 OR LEVGEO=2):
C       TIMINT NE 0. MEANS: INTERSECTION TIME IS KNOWN FROM AN EARLIER
C                CALL AND ITS VALUE IS TIMINT.
C                OTHERWISE (TIMINT=0) IT HAS TO BE CALCULATED IN THIS
C                CALL
C  IF (LEVGEO=3):
C       TIMINT NE 0. MEANS: INTERSECTION TIME IS KNOWN FROM AN EARLIER
C                CALL AND IS TO BE FOUND IN THE ARRAYS TIMPOL,IIMPOL.
C                OTHERWISE (TIMINT=0) IT HAS TO BE CALCULATED IN THIS
C                CALL
C  IF (LEVGEO=4):
C       TIMINT NE 0. MEANS: NOTHING
C  IF (LEVGEO=5):
C       TIMINT NE 0. MEANS: NOTHING
C  IF (LEVGEO=10):
C       TIMINT NE 0. MEANS: NOTHING
C
C  OUTPUT :
C       NJUMP = 1
C       NLSRFX=.FALSE.
C       MRSURF = INDEX OF NEXT SURFACE ALONG TRACK
C              = 0 IF NO NEXT SURFACE IS FOUND
C       PT = TIME (DISTANCE) TO REACH THIS SURFACE
C          = 1.D30 IF NO NEXT SURFACE IS FOUND
C       TIMINT(MRSURF) NE 0 INDICATES THAT FURTHER INTERSECTION WITH SAME SURFACE MRSURF WAS FOUND
C                      ALREADY IN THIS CALL, WHICH MAY BE USED IN A LATER CALL
C       NINCX = INDEX FOR DIRECTION IN GRID "RSURF": +1 OR -1
C             = 0 IF NO NEXT SURFACE FOUND
C       NRCELL:  NUMBER OF FINAL RADIAL CELL (IE. NOT MODIFIED)
C       IRCELL: NOT NEEDED ON RADIAL SURFACE. SET IRCELL=NRCELL
C
C  ADDITIONALLY IF NLPLG:
C  INPUT  :
C       IPOLG  = INDEX ON POLYGON OF INITIAL POINT X=(X0,Y0,Z0)
C  OUTPUT :
C       IPOLGN = INDEX ON POLYGON OF THE POINT X+PT*VEL
C         ( IF NO VALID INTERSECTION FOUND IN THIS CALL, THEN
C           IPOLGN=IPOLGO IS RETURNED, THE INDEX OF THE LAST
C           VALID POINT OF INTERSECTION FOUND IN EARLIER CALLS
C           (WHICH MAY BE THE INPUT VALUE IPOLG ITSELF) )
C
C
C  ADDITIONALLY, IF NLFEM:
C
C    TO BE WRITTEN
C
C  ADDITIONALLY, IF NLTET:
C
C    TO BE WRITTEN
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CTSURF
      USE EIRMOD_CADGEO
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CUPD
      USE EIRMOD_CPOLYG
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      USE EIRMOD_CTETRA
      USE EIRMOD_COMPRT
      USE EIRMOD_CLGIN
      USE EIRMOD_CTRIG
 
      IMPLICIT NONE
 
      REAL(DP), INTENT(OUT) :: PT
      REAL(DP) :: A(3,3), B(3), AB(3,3), PTS(4)
      REAL(DP) :: SIG1, SIG2, DET, XMU, XETA, EIRENE_SARRUS, AX, AY, V,
     .          TIMTET,
     .          ZZ, RICHTY, HELP, TM, RICHTX, YY, XX, S3, TIMT, TIM,
     .          S2, A1, A2, A3, B1, B2, B3, S1, V2, V3, ZR, ZEP1,
     .          ZSQRT, PS, VELYQ, XX0, VVELX, VELXQ, ZT1, ZT2, YVY,
     .          ZC1, DXA, TST, XA, ZB2, ZAB, ZAB2, ZB, Z0TEST,
     .          Y0TEST, ZA, T, PT1, PT2, PT3, PT4, V1, ESURF,
     .          XTEST, YTEST, X0SURF, DSRF, Y0Q, T1, T2, T3, T4, PNORMI
      INTEGER :: IRICH(2,4), ITSIDE(3,4)
      INTEGER :: IZELLO, NTIMT, IPOLGOO, IOB, I1, I2, ISW, KAN, KEN,
     .           EIRENE_ILLZ, IHELP, IZELL, NTMS, MXSF, NTMZ, NRMSRF,
     .           ICOS, IERR, IRS, NEWCEL, ITET, IT, IL, IS, NRI, MS, IR,
     .           ICALL, ITFRST, ISTS, MMSURF, ICOUP, J, K, I, JPOL,
     .           MPOL, IPOLGO, EIRENE_LEARC2, IP, ICELLR, MSAVE, ITRI,
     .           ISTS_CELL, ISD, nclpb, icts
!pb      INTEGER, ALLOCATABLE, SAVE :: ITRINO(:), ISIDNO(:)
!pb      INTEGER, SAVE :: NSTS_CELL
      INTEGER, ALLOCATABLE :: ITRINO(:), ISIDNO(:)
      INTEGER :: NSTS_CELL,nrtest,ist
      LOGICAL :: LLCUT(N2NDPLGS), LCTS(4)
      LOGICAL :: LNGB1, LNGB2, LNGB3, LNGB4,
     .           LCT1, LCT2, LCT3, LCT4, EIRENE_BITGET
 
      EQUIVALENCE (LCTS(1),LCT1), (LCTS(2),LCT2),
     .            (LCTS(3),LCT3), (LCTS(4),LCT4)
      EQUIVALENCE (PTS(1),PT1),   (PTS(2),PT2),
     .            (PTS(3),PT3),   (PTS(4),PT4)
 
      DATA ITSIDE /1,2,3,
     .             1,4,2,
     .             2,4,3,
     .             3,4,1/
      DATA IRICH / 1, -3,
     .             4,  1,
     .             5,  2,
     .             6,  3 /
      DATA ICALL /0/
      DATA ITFRST /0/

      SAVE
C
#ifdef TRACE
      IF (NLTRC) THEN
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) 'TIMER, INIT.: NRCELL,NLSRFX,MRSURF,TT,TL'
        WRITE (iunout,*)                NRCELL,NLSRFX,MRSURF,TT,TL
        WRITE (iunout,*) '              NPCELL,NLSRFY,MPSURF'
        WRITE (iunout,*)                NPCELL,NLSRFY,MPSURF
      ENDIF
#endif
C
      IRCELL=NRCELL
      PT=1.D30
      IF ((LEVGEO <= 4) .AND. (ABS(VELZ).EQ.1.D0)) THEN
c  zero velocity in x-y plane, i.e. flight in z-direction
        IPOLGN=IPOLG
        NCOUP=1
        ALPD(NCOUP)=PT
        LUPC(NCOUP)=0
        MUPC(NCOUP)=0
        JUPC(NCOUP)=IPOLGN
        RETURN
      ENDIF
C
C-------------------------------------------------------------------
      IF (LEVGEO.GT.1) GOTO 100
C
C****SLAB-MODEL IN X DIRECTION: FOR TESTING, ANALYTICAL REFERENCE CASES, ETC..
C
      IF (ABS(VELY).EQ.1.D0) THEN  
!dr  RULE OUT VELX=0.0, IN WHICH CASE THERE IS NO INTERSECTION WITH X-GRID
!dr  SAME TEST ALSO FOR Z DIRECTION IS ALREADY DONE ABOVE
        IPOLGN=IPOLG
        RETURN
      ENDIF
C
      IF (NJUMP.EQ.0) THEN
C  SAVE STARTING POINT, FLIGHT DIRECTION, FOR LATER CALLS ALONG SAME TRACK.
        XA=X0
        IF (NLSRFX) XA=RSURF(MRSURF)
        NINCX=1
        IF (VELX.LT.0.) NINCX=-1
C
C   RESET BRANCH MARKER
        NJUMP=1
        NLSRFX=.FALSE.
      ENDIF
C
      MRSURF=0
      IF (NRCELL.GT.0) THEN
        IR=NRCELL
        IF (NINCX.EQ.1) IR=IR+1
        DXA=RSURF(IR)-XA
C  TENTATIVE FLIGHT TIME (DISTANCE) TO NEXT X-GRID BOUNDARY
        TST=DXA/(VELX+EPS60)
C  IN CASE OF 2D OR 3D: CHECK IF WE HAVE LEFT THE GRID ACROSS ANOTHER BOUNDARY
C                       BEFORE REACHING THIS NEXT X-GRID BOUNDARY
        IF (NLTOR.AND.NLTRZ) THEN
          Z0TEST=Z0+TST*VELZ
          IF (ZSURF(1).GT.Z0TEST.OR.ZSURF(NT3RD).LT.Z0TEST) GOTO 5
        ENDIF
        IF (NLPOL) THEN
          Y0TEST=Y0+TST*VELY
          IF (PSURF(1).GT.Y0TEST.OR.PSURF(NP2ND).LT.Y0TEST) GOTO 5
        ENDIF
        PT=TST
        MRSURF=IR
5       CONTINUE
      ELSE
C  TRY TO FIND REENTRY SURFACE. CHECK ONLY NON DEFAULT RADIAL SURFACES
        DO 10 IR=1,NR1ST
          IF (INMP1I(IR,0,0).NE.0) THEN
            DXA=RSURF(IR)-XA
            TST=DXA/(VELX+EPS60)
            IF (TST.LE.0..OR.TST.GT.PT) GOTO 10
            IF (NLTOR.AND.NLTRZ) THEN
              Z0TEST=Z0+TST*VELZ
              IF (ZSURF(1).GT.Z0TEST.OR.ZSURF(NT3RD).LT.Z0TEST) GOTO 10
            ENDIF
            IF (NLPOL) THEN
              Y0TEST=Y0+TST*VELY
              IF (PSURF(1).GT.Y0TEST.OR.PSURF(NP2ND).LT.Y0TEST) GOTO 10
            ENDIF
            PT=TST
            MRSURF=IR
          ENDIF
10      CONTINUE
        IF (MRSURF.EQ.0) NINCX=0
      ENDIF
#ifdef TRACE
      IF (NLTRC) THEN
        WRITE (iunout,*) 'TIMER, OUT: PT,MRSURF,NINCX'
        WRITE (iunout,*)              PT,MRSURF,NINCX
      ENDIF
#endif
      RETURN
C
C
100   CONTINUE
C
C---------------------------------------------------------------------
      IF (LEVGEO.GT.2) GOTO 6000
C
      IF (NLELL) GOTO 1000
C
C****CIRCULAR MESH   (CONCENTRIC)
C
#ifdef TRACE
      IF (NLTRC) THEN
        WRITE (iunout,*) 'NJUMP,NINCX,NRCELL '
        WRITE (iunout,*)  NJUMP,NINCX,NRCELL
      ENDIF
#endif
      IF (NJUMP.EQ.0) THEN
C  SAVE STARTING POINT, FLIGHT DIRECTION, FOR LATER CALLS ALONG SAME TRACK.
        ZA=VELX*VELX+VELY*VELY
        ZB=X0*VELX+Y0*VELY
        ZB2=ZB*ZB
        ZAB=-ZB/(ZA+EPS60)
        ZAB2=ZA/(ZB2+EPS60)
C
C  TEST FOR DIRECTION
        NINCX=1
        IF (ZB.LT.0) NINCX=-NINCX
C
        IF (NLSRFX) THEN
          ZC1=RQ(MRSURF)
C  IN THIS CASE: DON'T TRUST THE VALUE OF NRCELL. RECOMPUTE FROM MRSURF
          IF (NINCX.EQ.1) NRCELL=MRSURF
          IF (NINCX.EQ.-1) NRCELL=MRSURF-1
        ELSE
          ZC1=X0*X0+Y0*Y0
        ENDIF
      ENDIF
C
      IF (NRCELL.GT.0) THEN
C  PARTICLE INSIDE STANDARD MESH
C  FIND NEXT SURFACE MRSURF
        NRI=0
        MRSURF=NRCELL
        IF (NINCX.EQ.1) MRSURF=MRSURF+1
        GOTO 204
      ENDIF
C
C  PARTICLE OUTSIDE STANDARD MESH
C  FIND NEXT SURFACE MRSURF
      NRI=1
      PS=PT
      MS=0
      IS=0
200   DO 205 IR=NRI,NR1ST
        IF (INMP1I(IR,0,0).NE.0) THEN
          MRSURF=IR
          GOTO 204
        ENDIF
205   CONTINUE
      PT=PS
      MRSURF=MS
      RETURN
C
C  CHECK SURFACE: MRSURF
C
204   IF (MRSURF.EQ.1) GO TO 201
      IF (TIMINT(MRSURF).GT.0.0) GO TO 303
203   ZR=RQ(MRSURF)
C  CHECK FOR ROOT
      ZEP1=ZAB2*(ZC1-ZR)
      IF (ZEP1.LT.1.0) GO TO 202
C  NO ROOT - PATH IN OTHER DIRECTION. THIS MUST BE THE
C  OUTWARD DIRECTION NOW, DUE TO CONVEXITY OF THE MESH
      IF (NRI.GT.0) GOTO 310
201   CONTINUE
      NINCX=1
      MRSURF=MRSURF+1
      IF (NJUMP.EQ.0) GO TO 203
      GO TO 303
202   CONTINUE
C
C  RESET BRANCH MARKER
      NJUMP=1
C
C  INTERSECTION TIMES
      ZSQRT=1.0+SQRT(1.0-ZEP1)
      ZT1=ZAB*ZEP1/ZSQRT
      ZT2=ZAB*ZSQRT
      NLSRFX=.FALSE.
C
CL              3.         RETURN RESULT
C
300   CONTINUE
C
C  CHECK SIGN OF RESULT
      IF(ZEP1.GT.0.0) GO TO 301
C  ONE ROOT NEGATIVE OR ZERO
      PT=MAX(ZT1,ZT2)
      GOTO 310
C
301   CONTINUE
      PT=ZT1
      TIMINT(MRSURF)=ZT2

c  ONE FURTHER VALID INTERSECTION FOUND: SAVE THIS FOR NEXT TIME:
      NIMINT = NIMINT+1
      IIMINT(NIMINT) = MRSURF
      GOTO 310
C
C  ROOT ALREADY KNOWN
303   PT=TIMINT(MRSURF)
      GOTO 310
C
310   CONTINUE
      IF (NRI.EQ.0) THEN
        RETURN
      ELSE
        IF (PT.GT.0.AND.PT.LT.PS) THEN
          IS=MRSURF-MAX(0,NINCX)
          MS=MRSURF
          PS=PT
        ENDIF
        NRI=IR+1
        GOTO 200
      ENDIF
C
1000  CONTINUE
C
C****ELLIPTICAL MESH  (NOT NECESSARILY CONCENTRIC NOR CONFOCAL)
C
      IF (NRCELL.LE.0) THEN
        WRITE (iunout,*) 'NRCELL.LE.0 IN TIMER: NOT READY '
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
 
C  PARTICLE OUTSIDE GRID OPTION NOT READY, SEE NLCRC
      NRI=0
 
C  COEFFICIENTS OF QUADRATIC
      IF (NJUMP.EQ.0) THEN
C  SAVE STARTING POINT, FLIGHT DIRECTION, FOR LATER CALLS ALONG SAME TRACK.
        YVY=Y0*VELY
        VELXQ=VELX*VELX
        VELYQ=VELY*VELY
        XX0=X0
        VVELX=VELX
        IF (NLSRFX) THEN
          X0SURF=X0-EP1(MRSURF)
          DSRF=ELLQ(MRSURF)
          ZR=RQ(MRSURF)
          Y0Q=(ZR-X0SURF*X0SURF)*DSRF
          MSAVE=MRSURF
C  IN THIS CASE: DON'T TRUST THE VALUE OF NRCELL. RECOMPUTE FROM MRSURF
          ZB=X0SURF*VELX+Y0*VELY/ELLQ(MRSURF)
C  TEST FOR DIRECTION. CASE: NLSRFX, MRSURF KNOWN
          NINCX=1
          IF (ZB.LT.0)     NINCX =-NINCX
          IF (NINCX.EQ.1)  NRCELL=MRSURF
          IF (NINCX.EQ.-1) NRCELL=MRSURF-1
C  NEXT SURFACE
          MRSURF=NRCELL
          IF (NINCX.EQ.1) MRSURF=MRSURF+1
          GOTO 1100
        ELSE
          MSAVE=0
          Y0Q=Y0*Y0
C   TEST FOR DIRECTION. CASE:.NOT.NLSRFX, NRCELL KNOWN
          NINCX=-1
          MRSURF=NRCELL
          ZB=(XX0-EP1(MRSURF))*VVELX+YVY/ELLQ(MRSURF)
          IF(ZB.LE.0.) GOTO 1200
          NINCX=1
          MRSURF=MRSURF+1
          GOTO 1100
        ENDIF
      ENDIF
C
2000  CONTINUE
      MRSURF=NRCELL
      IF (NINCX.EQ.1) MRSURF=MRSURF+1
1100  ZB=(XX0-EP1(MRSURF))*VVELX+YVY/ELLQ(MRSURF)
C
1200  IF(MRSURF.EQ.1) GO TO 2010
      IF(TIMINT(MRSURF).GT.0.0) GO TO 3030
      ESURF=XX0-EP1(MRSURF)
      DSRF=ELLQ(MRSURF)
      ZA=VELXQ+VELYQ/DSRF
      ZB2=ZB*ZB
      ZAB=-ZB/(ZA+EPS60)
      ZAB2=ZA/(ZB2+EPS60)
      ZC1=ESURF*ESURF+Y0Q/DSRF
      ZR=RQ(MRSURF)
      ZEP1=ZAB2*(ZC1-ZR)
C  CHECK FOR ROOT
      IF(ZEP1.LT.1.0) GO TO 2020
C  NO ROOT - PATH IN OTHER DIRECTION
C  MUST BE OUTWARD, BECAUSE OF CONVEXITY OF MESH
2010  CONTINUE
      NINCX=1
      MRSURF=MRSURF+1
      IF (NJUMP.EQ.0.) GO TO 1100
      GO TO 3030
2020  CONTINUE
C
C  RESET BRANCH MARKER
      NJUMP=1
C
C  INTERSECTION TIMES
      ZSQRT=1.0+SQRT(1.0-ZEP1)
      IF (MRSURF.EQ.MSAVE) THEN
        ZT1=0.
        ZEP1=0.
      ELSE
        ZT1=ZAB*ZEP1/ZSQRT
      ENDIF
      ZT2=ZAB*ZSQRT
      NLSRFX=.FALSE.
C
C---------------------------------------------------------------------
CL              3.         RETURN RESULT
C
3000  CONTINUE
C
C  CHECK SIGN OF RESULT
      IF(ZEP1.GT.0.0) GO TO 3010
C  ONE ROOT NEGATIVE OR ZERO
      PT=MAX(ZT1,ZT2)
      GOTO 3100
C
3010  CONTINUE
C  BOTH ROOTS POSITIVE - RETURN ONE AND SAVE OTHER
      PT=ZT1
      TIMINT(MRSURF)=ZT2
      NIMINT = NIMINT+1
      IIMINT(NIMINT) = MRSURF
      GOTO 3100
C
C  ROOT ALREADY KNOWN
3030  PT=TIMINT(MRSURF)
      GOTO 3100
 
C
3100  CONTINUE
      XTEST=X0+PT*VELX
      YTEST=Y0+PT*VELY
      IF (NLPOL) IPOLGN=EIRENE_LEARC2(XTEST,YTEST,NRCELL,NPANU,'TIMER')
#ifdef TRACE
      IF (NLTRC)
     .  WRITE (iunout,*) 'IPOLGN,NRCELL FROM TIMER ',IPOLGN,NRCELL
#endif
      IF (NRI.EQ.0) THEN
        RETURN
      ELSE
        WRITE (iunout,*) ' OPTION NRI > 0 NOT READY in TIMER '
        CALL EIRENE_EXIT_OWN(1)
CPB        IF (PT.GT.0.AND.PT.LT.PS) THEN
CPB          IS=MRSURF-MAX(0,NINCX)
CPB          MS=MRSURF
CPB          PS=PT
CPB        ENDIF
CPB        NRI=IR+1
CPB        GOTO 2000
      ENDIF
C
C     POLYGONS
C
6000  CONTINUE
C
      IF (LEVGEO.GT.3) GOTO 8000
C
C  POLYGON MESH
C
      IF (NRCELL.NE.0) THEN
#ifdef TRACE
        IF (NLTRC) WRITE (iunout,*) ' TIMER: IN NEIGHBOR PART'
#endif
        LNGB1=.TRUE.
        LNGB2=.TRUE.
        LNGB3=.TRUE.
        LNGB4=.TRUE.
C  SET INDICES OF STARTING CELL
        IR=NRCELL
        IP=NPCELL
        IF (.NOT.NLPOL) IP=IPOLG

        IF (NJUMP.EQ.1) THEN
C  CONTINUATION OF A TRAJECTORY WITH SAME STARTING POINT AND SAME DIRECTION 
          IR=ICELLR
          IP=IPOLGN
          IF (NINCX.EQ.1) LNGB1=.FALSE.
          IF (NINCX.EQ.-1) LNGB3=.FALSE.
        ENDIF

        IF (NLSRFX) THEN
!pb05032015          NLSRFX=.FALSE.
          IF (NRCELL.EQ.MRSURF) THEN
            LNGB1=.FALSE.
          ELSE
            LNGB3=.FALSE.
          ENDIF
        ENDIF
        IF (NLSRFY) THEN
          IF (IPOLG.EQ.MPSURF) THEN
            IF (NPCELL.EQ.MPSURF) THEN
              LNGB4=.FALSE.
              IP=NGHPLS(2,IR,MPSURF)
            ELSE
              LNGB2=.FALSE.
              IP=NGHPLS(4,IR,MPSURF)
            ENDIF
          ELSE
            LNGB2=.FALSE.
          ENDIF
        ENDIF
C
        NCOUP=0
C  CALCULATE INTERSECTIONS OF FLIGHT WITH CELL BOUNDARIES
6001    CONTINUE
#ifdef TRACE
        IF (NLTRC) WRITE (iunout,*) ' IR,IP,LNGB1,LNGB2,LNGB3,LNGB4',
     .                           IR,IP,LNGB1,LNGB2,LNGB3,LNGB4
        IF (NLTRC) WRITE (iunout,*) ' velx, vely, velz ',
     .                           velx, vely, velz
#endif
        T1=-1.D30
        T2=-1.D30
        T3=-1.D30
        T4=-1.D30
        PT1=-1.D30
        PT2=-1.D30
        PT3=-1.D30
        PT4=-1.D30
        IF (LNGB1)
     .  T1=((XPOL(IR,IP)-X0)*VELY-(YPOL(IR,IP)-Y0)*VELX)/
     .     (VELX*VPLY(IR,IP)-VELY*VPLX(IR,IP)+EPS60)
        IF (LNGB2)
     .  T2=((XPOL(IR,IP+1)-X0)*VELY-(YPOL(IR,IP+1)-Y0)*VELX)/
     .      (VELX*VVTY(IR,IP+1)-VELY*VVTX(IR,IP+1)+EPS60)
        IF (LNGB3)
     .  T3=((XPOL(IR+1,IP)-X0)*VELY-(YPOL(IR+1,IP)-Y0)*VELX)/
     .     (VELX*VPLY(IR+1,IP)-VELY*VPLX(IR+1,IP)+EPS60)
        IF (LNGB4)
     .  T4=((XPOL(IR,IP)-X0)*VELY-(YPOL(IR,IP)-Y0)*VELX)/
     .      (VELX*VVTY(IR,IP)-VELY*VVTX(IR,IP)+EPS60)
#ifdef TRACE
        IF (NLTRC) WRITE (iunout,*) ' T1,T2,T3,T4 ',T1,T2,T3,T4
#endif
        LCT1 = T1.GE.0.D0 .AND. T1.LE.1.D0
        LCT2 = T2.GE.0.D0 .AND. T2.LE.1.D0
        LCT3 = T3.GE.0.D0 .AND. T3.LE.1.D0
        LCT4 = T4.GE.0.D0 .AND. T4.LE.1.D0
C  CALCULATE TIME OF FLIGHT FROM STARTING POINT X=(X0,Y0,Z0)
C  TO THE BOUNDARY OF THE ACTUELL CELL
        IF (ABS(VELX).GT.ABS(VELY)) THEN
          IF (LCT1) PT1=(XPOL(IR,IP)-X0+VPLX(IR,IP)*T1)/VELX
          IF (LCT2) PT2=(XPOL(IR,IP+1)-X0+VVTX(IR,IP+1)*T2)/VELX
          IF (LCT3) PT3=(XPOL(IR+1,IP)-X0+VPLX(IR+1,IP)*T3)/VELX
          IF (LCT4) PT4=(XPOL(IR,IP)-X0+VVTX(IR,IP)*T4)/VELX
        ELSE
          IF (LCT1) PT1=(YPOL(IR,IP)-Y0+VPLY(IR,IP)*T1)/VELY
          IF (LCT2) PT2=(YPOL(IR,IP+1)-Y0+VVTY(IR,IP+1)*T2)/VELY
          IF (LCT3) PT3=(YPOL(IR+1,IP)-Y0+VPLY(IR+1,IP)*T3)/VELY
          IF (LCT4) PT4=(YPOL(IR,IP)-Y0+VVTY(IR,IP)*T4)/VELY
        ENDIF
        LCT1 = LCT1 .AND. PT1.GE.0.D0
        LCT2 = LCT2 .AND. PT2.GE.0.D0
        LCT3 = LCT3 .AND. PT3.GE.0.D0
        LCT4 = LCT4 .AND. PT4.GE.0.D0
#ifdef TRACE
        IF (NLTRC) WRITE (iunout,*) ' LCT1,LCT2,LCT3,LCT4 ',
     .                           LCT1,LCT2,LCT3,LCT4
        IF (NLTRC) WRITE (iunout,*) ' PT1,PT2,PT3,PT4 ',PT1,PT2,PT3,PT4
#endif

        ICTS=COUNT(LCTS)
C  MORE THAN ONE TENTATIVE INTERECTION FOUND
        IF (ICTS > 1) THEN
          WHERE (PTS < 0._DP)
            PTS = 1.E30_DP
          END WHERE
          ISD = MINLOC(PTS,DIM=1)
          LCTS = .FALSE.
          LCTS(ISD) = .TRUE.
          T = PTS(ISD)
C  EXACTLY ONE TENTATIVE INTERESECTION FOUND
        ELSE IF (ICTS == 1) THEN
          T=MAX(PT1,PT2,PT3,PT4)
C  NO TENTATIVE INTERSECTION FOUND, T IS NOT SET.
C       ELSE
C  THIS MAY BE A RESULT OF NUMERICAL ROUNDING ERRORS,
C  SO SOME FURTHER ATTEMPTS WILL BE MADE BELOW.....
        END IF
#ifdef TRACE
        IF (NLTRC) WRITE (iunout,*) ' T = ',T
#endif

C  IF INTERSECTION WITH POLOIDAL BOUNDARY CONTINUE WITH NEIGHBORING CELL
        IF (LCT2.OR.LCT4) THEN
          LNGB1=.TRUE.
          LNGB2=.TRUE.
          LNGB3=.TRUE.
          LNGB4=.TRUE.
          NCOUP=NCOUP+1
          ALPD(NCOUP)=T
          IF (LCT2) THEN
            LUPC(NCOUP)=IP+1
            MUPC(NCOUP)=1
            JUPC(NCOUP)=IP
            IP=NGHPOL(2,IR,IP)
            LNGB4=.FALSE.
          ENDIF
          IF (LCT4) THEN
            LUPC(NCOUP)=IP
            MUPC(NCOUP)=-1
            JUPC(NCOUP)=IP
            IP=NGHPOL(4,IR,IP)
            LNGB2=.FALSE.
          ENDIF
          ISTS=INMP2I(IR,LUPC(NCOUP),0)
!pb          IF ((.not.NLPOL.or.ISTS.eq.0).and.ip.ne.0) goto 6001
          if (ip > 0) then
            nclpb = IR+((IP-1)+(NTCELL-1)*NP2T3)*NR1P2+NBLCKA
          else
            nclpb = ncell
          end if
          IF ((ityp.ne.3.and.(.not.NLPOL.or.ISTS.eq.0).and.ip.ne.0)
     .       .and. (.not.ldamcel(nclpb)))
     .       goto 6001
C  NO NEIGHBORING CELL: PARTICLE HAS HIT A POLOIDAL BOUNDARY OF THE MESH
          MRSURF=0
          PT=1.D30
          NINCX=0
        ELSEIF (LCT1.OR.LCT3) THEN
C  INTERSECTION WITH RADIAL CELL BOUNDARY FOUND
          PT=T
          IPOLGN=IP
          LNGB1=.TRUE.
          LNGB2=.TRUE.
          LNGB3=.TRUE.
          LNGB4=.TRUE.
          NCOUP=NCOUP+1
          ALPD(NCOUP)=T
          LUPC(NCOUP)=0
          MUPC(NCOUP)=0
          JUPC(NCOUP)=IP
          IF (LCT1) THEN
            MRSURF=IR
            NINCX=-1
            ICELLR=IR-1
          ELSE
            MRSURF=IR+1
            ICELLR=IR+1
            NINCX=1
          ENDIF
        ELSE
C..................................................................................
C  NO INTERSECTION FOUND
          IF (NLSRFX.AND.NJUMP.EQ.0) THEN
C  PLAY SAVE: TRY ONCE AGAIN, IF PARTICLE ON RAD. SURFACE
#ifdef TRACE
            IF (NLTRC)
     .        WRITE (iunout,*)
     .        ' NO INTERSECTION IN TIMER. TRY ONCE AGAIN '
#endif
            MMSURF=MSURF
            IF (MSURF.GT.NLIM) MMSURF=-MSURF+NLIM
#ifdef TRACE
            IF (NLTRC) WRITE (iunout,*) 'NPANU, MSURF ',NPANU,MMSURF
#endif
            IF (.NOT.LNGB1) THEN
              ICELLR=IR-1
              IRCELL=IR-1
              IR=IR-1
              LNGB1=.TRUE.
              LNGB3=.FALSE.
            ELSEIF (.NOT.LNGB3) THEN
              ICELLR=IR+1
              IRCELL=IR+1
              IR=IR+1
              LNGB1=.FALSE.
              LNGB3=.TRUE.
            ENDIF
            LNGB2=.TRUE.
            LNGB4=.TRUE.
            NJUMP=1
            NRCELL=IR
            GOTO 6001

          ELSEIF (NLSRFY.AND.NJUMP.EQ.0) THEN
C  PLAY SAVE: TRY ONCE AGAIN, IF PARTICLE ON POL. SURFACE
#ifdef TRACE
            IF (NLTRC)
     .        WRITE (iunout,*)
     .        ' NO INTERSECTION IN TIMER. TRY ONCE AGAIN '
#endif
            MMSURF=MSURF
            IF (MSURF.GT.NLIM) MMSURF=-MSURF+NLIM
#ifdef TRACE
            IF (NLTRC) WRITE (iunout,*) 'NPANU, MSURF ',NPANU,MMSURF
#endif
            IF (.NOT.LNGB4) THEN
              IP=NGHPLS(4,IR,MPSURF)
              LNGB4=.TRUE.
              LNGB2=.FALSE.
            ELSEIF (.NOT.LNGB2) THEN
              IP=NGHPLS(2,IR,MPSURF)
              LNGB2=.TRUE.
              LNGB4=.FALSE.
            ENDIF
            LNGB1=.TRUE.
            LNGB3=.TRUE.
            NJUMP=1
            GOTO 6001
          ENDIF
          WRITE (iunout,*) ' ERROR: NO INTERSECTION FOUND IN TIMER '
          WRITE (iunout,*) ' NPANU: ',NPANU
          MRSURF=0
          PT=1.D30
          NINCX=0
          ICELLR=0
          IPOLGN=IP
        ENDIF
        NJUMP=1
#ifdef TRACE
        IF (NLTRC) THEN
          WRITE (iunout,*) ' PT,MRSURF,NINCX,IPOLGN ',
     .                  PT,MRSURF,NINCX,IPOLGN
          WRITE (iunout,*) 'NCOUP ',NCOUP
          DO ICOUP=1,NCOUP
            WRITE (iunout,*) 'ICOUP,ALPD(ICOUP) ',ICOUP,ALPD(ICOUP)
          ENDDO
        ENDIF
#endif
        RETURN
      ENDIF
C
C  PARTICLE OUTSIDE STANDARD MESH, IN ADDITIONAL CELL NACELL
C  NRCELL=0
C
      IF (NLSRFX) THEN
        NLSRFX=.FALSE.
        JPOL=IPOLG
        MPOL=MRSURF
      ELSE
        JPOL=0
        MPOL=0
      ENDIF
C
      IF (NJUMP.EQ.0) IPOLGO=IPOLG
C
      PT=1.D30
      MRSURF=0
      IPOLGN=IPOLGO
C
      DO 6100 I=1,NR1ST
        ISTS=INMP1I(I,0,0)
C  SURFACE I IS NOT A NON DEFAULT RADIAL STANDARD SURFACE
        IF (ISTS.EQ.0) GOTO 6100

        IF (TIMINT(I).EQ.0) THEN
          TIMINT(I)=1
          NIMINT = NIMINT+1
          IIMINT(NIMINT) = I
          NTIM(I)=0
C
C   SEARCH FOR ALL POSSIBLE INTERSECTIONS WITHIN THE CELL
C
          DO 6011 J=1,NRPLG
6011        LLCUT(J)=.FALSE.
C
          DO 6012 J=1,NPPLG
            DO 6012 K=NPOINT(1,J),NPOINT(2,J)-1
              V1=(YPOL(I,K+1)-Y0)*VELX-(XPOL(I,K+1)-X0)*VELY
              V2=(YPOL(I,K)-Y0)*VELX-(XPOL(I,K)-X0)*VELY
!pb allow only intersection with nondefault standard surfaces
!pb no intersection with transparent parts
!pb              LLCUT(K)=V1*V2.LE.0.
              LLCUT(K)=(V1*V2.LE.0.) .AND. (INMP1I(I,K,0) /= 0)
6012      CONTINUE
C
          IF (I.EQ.MPOL) LLCUT(JPOL)=.FALSE.
          KAN=EIRENE_ILLZ(NRPLG,LLCUT,1)+1
          KEN=NRPLG-EIRENE_ILLZ(NRPLG,LLCUT,-1)
C
C   COMPUTE THE FLIGHT TIMES TO THE INTERSECTION POINTS
C
          DO 6013 K=KAN,KEN
            IF (LLCUT(K)) THEN
              T1=((XPOL(I,K)-X0)*VPLY(I,K)-(YPOL(I,K)-Y0)*VPLX(I,K))
     .           /(VELX*VPLY(I,K)-VELY*VPLX(I,K)+EPS60)
              IF (T1.GT.0.) THEN
                NTIM(I)=NTIM(I)+1
                TIMPOL(I,NTIM(I))=T1
                IIMPOL(I,NTIM(I))=K
              ENDIF
            ENDIF
6013      CONTINUE
C
6015      ISW=0
          DO 6020 J=1,NTIM(I)-1
            IF (TIMPOL(I,J).LT.TIMPOL(I,J+1)) THEN
              ISW=ISW+1
              HELP=TIMPOL(I,J)
              TIMPOL(I,J)=TIMPOL(I,J+1)
              TIMPOL(I,J+1)=HELP
              IHELP=IIMPOL(I,J)
              IIMPOL(I,J)=IIMPOL(I,J+1)
              IIMPOL(I,J+1)=IHELP
            ENDIF
6020      CONTINUE
          IF (ISW.GT.0.AND.NTIM(I).GT.2) GOTO 6015
#ifdef TRACE
          IF (NLTRC) THEN
            WRITE (iunout,*) ' SURFACE NO. = ',I
            WRITE (iunout,*) ' TIMPOL ',(TIMPOL(I,J),J=1,NTIM(I))
            WRITE (iunout,*) ' IIMPOL ',(IIMPOL(I,J),J=1,NTIM(I))
          ENDIF
#endif
C
        ENDIF
C
C  FIND PT, MRSURF, IPOLGN
        IF (NTIM(I).GT.0) THEN
          IF (TIMPOL(I,NTIM(I)).LT.PT) THEN
            MRSURF=I
            PT=TIMPOL(I,NTIM(I))
            IPOLGN=IIMPOL(I,NTIM(I))
          ENDIF
        ENDIF
6100  CONTINUE
C
C  FIND NINCX
      NINCX=0
      IF (MRSURF.NE.0) THEN
        NTIM(MRSURF)=NTIM(MRSURF)-1
        NINCX=SIGN(1._DP,VELX*PLNX(MRSURF,IPOLGN)+
     .        VELY*PLNY(MRSURF,IPOLGN))
      ENDIF
C
      NJUMP=1
      IPOLGO=IPOLGN
C
#ifdef TRACE
      IF (NLTRC) WRITE (iunout,*) 'PT,MRSURF,IPOLGN,NINCX ',
     .                        PT,MRSURF,IPOLGN,NINCX
#endif
      RETURN
C
C
8000  CONTINUE
C
      IF (LEVGEO.GT.4) GOTO 10000
C
C   FINITE ELEMENT DISCRETISATION (MESH TRIANGLES IN X-Y- PLANE)
C
#ifdef TRACE
      IF (NLTRC) THEN
        WRITE (iunout,*) ' TIMER,NJUMP,NRCELL ',NJUMP,NRCELL
        WRITE (iunout,*) '       NLSRFX,IPOLG ',NLSRFX,IPOLG
      ENDIF
#endif
C
      IF (NJUMP.EQ.0) THEN
8001    CONTINUE
        IZELL = NRCELL
        IPOLGO=0
        IF (NLSRFX) THEN
          IPOLGO=IPOLG
        ENDIF
C     ELSEIF (NJUMP.EQ.1) THEN
      ENDIF
 
      IF (NRCELL > 0) THEN
#ifdef TRACE
        IF (NLTRC) THEN
          WRITE (iunout,*) ' IZELL,IPOLGO ',IZELL,IPOLGO
          WRITE (iunout,*) ' X0,Y0 ',X0,Y0
          WRITE (iunout,*) ' VELX,VELY ',VELX,VELY
          CALL EIRENE_LEER(1)
        ENDIF
#endif
        XX = X0
        YY = Y0
        TM=0.
        IF (IZELL.EQ.0) GOTO 9999
C
8020    CONTINUE
        DO 8010,J=1,3
          IF (J.NE.IPOLGO) THEN
            RICHTX = VTRIX(J,IZELL)
            RICHTY = VTRIY(J,IZELL)
            AX = XTRIAN(NECKE(J,IZELL))-XX
            AY = YTRIAN(NECKE(J,IZELL))-YY
            V  = (AX*VELY-AY*VELX)/(RICHTY*VELX-RICHTX*VELY+EPS60)
            IF (V.GE.0..AND.V.LE.1.) THEN
              IF (ABS(VELX).GT.ABS(VELY)) THEN
                T=(AX+V*RICHTX)/VELX
              ELSE
                T=(AY+V*RICHTY)/VELY
              ENDIF
              IF (T .GT. 0.) THEN
                XX = XX + T * VELX
                YY = YY + T * VELY
                TM = TM + T
                TIMINT(IZELL) = TM
                NIMINT = NIMINT+1
                IIMINT(NIMINT) = IZELL
C  NUMMER DER GESCHNITTENEN SEITE DES ALTEN DREIECKS
                NTIM(IZELL) = J
#ifdef TRACE
                IF (NLTRC) THEN
                  WRITE (iunout,*) 'IZELL,XX,YY,J ',IZELL,XX,YY,J
                ENDIF
#endif
C  ZELLENNUMMER DES NEUEN DREIECKS
                IZELLO=IZELL
                IZELL = NCHBAR(J,IZELLO)
                IF (IZELL .EQ. 0) GOTO 8050
C  SEITENNUMMER DES NEUEN DREIECKS
                IPOLGO = NSEITE(J,IZELLO)
#ifdef TRACE
                IF (NLTRC) THEN
                  WRITE(iunout,*) 'GEHE IN ZELLE ',IZELL,' TM= ',TM
                  WRITE(iunout,*) 'DORT AUF SEITE IPOLGO ',IPOLGO
                ENDIF
#endif
                GOTO 8050
              ENDIF
            ENDIF
          ENDIF
8010    CONTINUE
C
#ifdef TRACE
        IF (NLTRC)
     .    WRITE (iunout,*) ' NO INTERSECTION FOUND IN TRIANGLE '
#endif
        PT=1.D30
        NINCX=0
C  DURCHSUCHE NACHBARDREIECK
        IF (NLSRFX) THEN
          IZELL=NCHBAR(IPOLG,NRCELL)
          IPOLGO=NSEITE(IPOLG,NRCELL)
          NRCELL=IZELL
          NLSRFX=.FALSE.
          GOTO 8020
        ENDIF

C  KEIN SCHNITTPUNKT GEFUNDEN UND NLSRFX=.FALSE.
        write (iunout,*) ' NO INTERSECTION FOUND IN TRIANGLE '
        write (iunout,*) ' PARTICLE TRAJECTORY STOPPED'
        write (iunout,*) ' NPANU ',NPANU
        LGPART = .FALSE.

        RETURN
C
C   IZELL and IPOLGO of next triangle are now given
C
C   THIS IS DONE FOR NJUMP=0 AND NJUMP=1
C
8050  CONTINUE
      NLSRFX=.FALSE.
      NJUMP=1
      MRSURF=NRCELL
      PT=TIMINT(MRSURF)
      TIMINT(MRSURF)=0.
      NIMINT = NIMINT+1
      IIMINT(NIMINT) = MRSURF
      IPOLGN=NTIM(NRCELL)
      ISTS=ABS(INMTI(IPOLGN,NRCELL))
      IF (ISTS.EQ.0) THEN
c  on internal grid surface:  ISTS=0
        NINCX=NCHBAR(IPOLGN,NRCELL)-NRCELL
      ELSEIF (ISTS.GT.0.AND.
     .        ISTS.LE.NLIM+NSTSI) THEN
C  on non-default surface (add. or. std.): ISTS=INMTI(IPOLGN,NRCELL)
!pb        IF (ILIIN(ISTS) == 0) THEN
        IF (ILIIN(ISTS) <= 0) THEN
          NINCX=NCHBAR(IPOLGN,NRCELL)-NRCELL
        ELSE
          NINCX=SIGN(1,INMTI(IPOLGN,NRCELL))
        END IF
      ELSE
        GOTO 9999
      ENDIF
#ifdef TRACE
      IF (NLTRC) WRITE (iunout,*) ' NRCELL,MRSURF,PT,NINCX,IPOLGN ',
     .                         NRCELL,MRSURF,PT,NINCX,IPOLGN
#endif
 
      RETURN
      END IF
 
C PARTICLE OUTSIDE STANDARD MESH IN CELL NACELL, NRCELL=0
 
      IF (.NOT.ALLOCATED(ITRINO)) THEN
        NSTS_CELL = COUNT(INMTI(1:3,1:NTRII) .NE. 0)
        ALLOCATE (ITRINO(NSTS_CELL))
        ALLOCATE (ISIDNO(NSTS_CELL))
        ISTS_CELL = 0
        DO ITRI=1,NTRII
          DO J=1,3
            IF (INMTI(J,ITRI) .NE. 0) THEN
              ISTS_CELL = ISTS_CELL + 1
              ITRINO(ISTS_CELL) = ITRI
              ISIDNO(ISTS_CELL) = J
            END IF
          END DO
        END DO
      END IF
 
      XX = X0
      YY = Y0
      PT=1.E30_DP
      ISTS_CELL = 0
      IF (NLSRFX) THEN
        NLSRFX=.FALSE.
        MPOL=MRSURF
        JPOL=IPOLG
      ELSE
        MPOL=0
        JPOL=0
      END IF
      MRSURF=0
 
      DO I=1, NSTS_CELL
        IZELL=ITRINO(I)
        J=ISIDNO(I)
        IF ((IZELL == MPOL) .AND. (J == JPOL)) CYCLE
        RICHTX = VTRIX(J,IZELL)
        RICHTY = VTRIY(J,IZELL)
        AX = XTRIAN(NECKE(J,IZELL))-XX
        AY = YTRIAN(NECKE(J,IZELL))-YY
        V  = (AX*VELY-AY*VELX)/(RICHTY*VELX-RICHTX*VELY+EPS60)
        IF (V.GE.0..AND.V.LE.1.) THEN
          IF (ABS(VELX).GT.ABS(VELY)) THEN
            T=(AX+V*RICHTX)/VELX
          ELSE
            T=(AY+V*RICHTY)/VELY
          ENDIF
#ifdef TRACE
          IF (NLTRC)
     .      WRITE (IUNOUT,*) ' INTERSECTION WITH TRIANGLE ',IZELL,
     .                       ' SIDE ',J,' FOUND, PT,T = ',PT,T
#endif
          IF ((T .GT. 0.) .AND. ( T < PT)) THEN
            PT = T
            ISTS_CELL = I
          ENDIF
        ENDIF
      END DO
 
      IF (ISTS_CELL > 0) THEN
!  REENTRY FOUND
        NLSRFX=.TRUE.
        NJUMP=1
        MRSURF = ITRINO(ISTS_CELL)
        IPOLGN = ISIDNO(ISTS_CELL)
        NINCX = ITRINO(ISTS_CELL)
#ifdef TRACE
        IF (NLTRC) WRITE (iunout,*) ' NRCELL,MRSURF,PT,NINCX,IPOLGN ',
     .                                NRCELL,MRSURF,PT,NINCX,IPOLGN
#endif
      END IF
 
      RETURN
C
10000 CONTINUE
 
      IF (LEVGEO > 5) GOTO 12000
C
C  MESH OF TETRAHEDRA, FULL 3D DISCRETIZATION OF COMPUTATIONAL VOLUME
C
      TIMTET = 1.E30
      NTIMT = 0
      IF (NRCELL .NE. 0) THEN
        IF (NJUMP.EQ.0) THEN
          IZELL = NRCELL
          IPOLGO=0
          IF (NLSRFX) THEN
            IPOLGO=IPOLG
          ENDIF
C       ELSEIF (NJUMP.EQ.1) THEN
        ENDIF
#ifdef TRACE
        IF (NLTRC) THEN
          WRITE (iunout,*) ' IZELL,IPOLGO ',IZELL,IPOLGO
          WRITE (iunout,*) ' X0,Y0,Z0 ',X0,Y0,Z0
          WRITE (iunout,*) ' VELX,VELY,VELZ ',VELX,VELY,VELZ
        ENDIF
#endif
 
        XX = X0
        YY = Y0
        ZZ = Z0
        TM=0.
        IF (IZELL.EQ.0) THEN
          WRITE (iunout,*) ' IZELL = 0 IN TIMER '
          CALL EIRENE_EXIT_OWN(1)
        END IF
C
C  CHECK TETRAHEDRON IZELL FOR INTERSECTION WITH TRAJECTORY
C
11020   CONTINUE
        DO J=1,4
          IF (J.NE.IPOLGO) THEN
            SIG1=SIGN(1,IRICH(1,J))
            I1=ABS(IRICH(1,J))
            SIG2=SIGN(1,IRICH(2,J))
            I2=ABS(IRICH(2,J))
            PNORMI=RINCRC(J,IZELL)
            A(1:3,1) = (/ VTETX(I1,IZELL), VTETY(I1,IZELL),
     .                    VTETZ(I1,IZELL)/) * SIG1 * PNORMI
            A(1:3,2) = (/ VTETX(I2,IZELL), VTETY(I2,IZELL),
     .                    VTETZ(I2,IZELL)/) * SIG2 * PNORMI
            A(1:3,3) = (/ -VELX, -VELY, -VELZ /) * PNORMI
            B(1:3) = (/ XX-XTETRA(NTECK(ITSIDE(1,J),IZELL)),
     .                  YY-YTETRA(NTECK(ITSIDE(1,J),IZELL)),
     .                  ZZ-ZTETRA(NTECK(ITSIDE(1,J),IZELL)) /) * PNORMI
            DET = A(1,1) * A(2,2) * A(3,3)
     .          + A(1,2) * A(2,3) * A(3,1)
     .          + A(1,3) * A(2,1) * A(3,2)
     .          - A(3,1) * A(2,2) * A(1,3)
     .          - A(3,2) * A(2,3) * A(1,1)
     .          - A(3,3) * A(2,1) * A(1,2)
            IF (ABS(DET) < 1.D-5) CYCLE
            AB = A
            AB(:,1) = B
            XMU = ( AB(1,1) * AB(2,2) * AB(3,3)
     .            + AB(1,2) * AB(2,3) * AB(3,1)
     .            + AB(1,3) * AB(2,1) * AB(3,2)
     .            - AB(3,1) * AB(2,2) * AB(1,3)
     .            - AB(3,2) * AB(2,3) * AB(1,1)
     .            - AB(3,3) * AB(2,1) * AB(1,2) ) / DET
            IF (XMU .GE.0.D0 .AND. XMU .LE.1.D0) THEN
              AB = A
              AB(:,2) = B
              XETA = ( AB(1,1) * AB(2,2) * AB(3,3)
     .               + AB(1,2) * AB(2,3) * AB(3,1)
     .               + AB(1,3) * AB(2,1) * AB(3,2)
     .               - AB(3,1) * AB(2,2) * AB(1,3)
     .               - AB(3,2) * AB(2,3) * AB(1,1)
     .               - AB(3,3) * AB(2,1) * AB(1,2) ) / DET
              IF ((XETA.GE.0.D0 .AND. XETA.LE.1.D0) .AND.
     .            (XMU+XETA <= 1.D0)) THEN
                AB = A
                AB(:,3) = B
                T = ( AB(1,1) * AB(2,2) * AB(3,3)
     .              + AB(1,2) * AB(2,3) * AB(3,1)
     .              + AB(1,3) * AB(2,1) * AB(3,2)
     .              - AB(3,1) * AB(2,2) * AB(1,3)
     .              - AB(3,2) * AB(2,3) * AB(1,1)
     .              - AB(3,3) * AB(2,1) * AB(1,2) ) / DET
                IF (T .GT. 0.) THEN
!            IF ((XMU .GE.0.D0 .AND. XMU .LE.1.D0) .AND.
!     .          (XETA.GE.0.D0 .AND. XETA.LE.1.D0) .AND.
!     .          (XMU+XETA <= 1.D0) .AND. (T .GT. 0.)) THEN
                  XX = XX + T * VELX
                  YY = YY + T * VELY
                  ZZ = ZZ + T * VELZ
#ifdef TRACE
                  if (nltrc) then
                     write(iunout,*) ' i1, i2 ',i1,i2
                     write(iunout,'(1x,a,3es14.7)') ' velxyz ',
     .                     velx,vely,velz
                     write(iunout,'(1x,a,3es14.7)') ' a11, a12, a13 ',
     .                     a(1,:)
                     write(iunout,'(1x,a,3es14.7)') ' a21, a22, a23 ',
     .                     a(2,:)
                     write(iunout,'(1x,a,3es14.7)') ' a31, a32, a33 ',
     .                     a(3,:)
                     write(iunout,'(1x,a,3es14.7)') ' b1, b2, b3 ',b(:)
                     write (iunout,'(1x,a,3es14.7)') ' det, xmu ',
     .                     det,xmu
                  end IF
#endif
                  TM = TM + T
                  TIMTET = TM
C  NUMMER DER GESCHNITTENEN SEITE DES ALTEN DREIECKS
                  NTIMT = J
#ifdef TRACE
                  IF (NLTRC) THEN
                    WRITE (iunout,*) 'IZELL,XX,YY,ZZ,J ',
     .                                IZELL,XX,YY,ZZ,J
                  ENDIF
#endif
C  ZELLENNUMMER DES NEUEN DREIECKS
                  IZELLO=IZELL
                  IPOLGOO=IPOLGO
                  IZELL = NTBAR(J,IZELLO)
                  PT = TM
                  IF (IZELL .EQ. 0) GOTO 11050
C  SEITENNUMMER DES NEUEN DREIECKS
                  IPOLGO = NTSEITE(J,IZELLO)
#ifdef TRACE
                  IF (NLTRC) THEN
                    WRITE(iunout,*) 'GEHE IN ZELLE ',IZELL,' TM= ',TM
                    WRITE(iunout,*) 'DORT AUF SEITE IPOLGO ',IPOLGO
                  ENDIF
#endif
                  GOTO 11050
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        END DO
C
#ifdef TRACE
        IF (NLTRC)
     .    WRITE (iunout,*) ' NO INTERSECTION FOUND IN TETRAHEDRON '
#endif
        PT=1.D30
        NINCX=0
C  DURCHSUCHE NACHBARDREIECK
        IF (NLSRFX) THEN
          IZELL=NTBAR(IPOLG,NRCELL)
          IPOLGO=NTSEITE(IPOLG,NRCELL)
          NRCELL=IZELL
          NLSRFX=.FALSE.
          GOTO 11020
        ENDIF
C  KEIN SCHNITTPUNKT GEFUNDEN UND NLSRFX=.FALSE.
        write (iunout,*) ' NO INTERSECTION FOUND IN TETRAHEDRON '
        write (iunout,*) ' PARTICLE TRAJECTORY STOPPED'
        write (iunout,*) ' NPANU ',NPANU
        LGPART = .FALSE.
        RETURN
C
C
C   THIS IS DONE FOR NJUMP=0 AND NJUMP=1
C
11050   CONTINUE
        NLSRFX=.FALSE.
        NJUMP=1
        MRSURF=NRCELL
        PT=TIMTET
        IPOLGN=NTIMT
        ISTS=ABS(INMTIT(IPOLGN,NRCELL))
        IF (ISTS.EQ.0) THEN
          NINCX=NTBAR(IPOLGN,NRCELL)-NRCELL
        ELSEIF (ISTS.GT.0.AND.
     .          ISTS.LE.NLIM+NSTSI) THEN
C  ON NON DEFAULT SURFACE (ADD. OR STD.) ISTS=INMTI(IPOLGN,NRCELL)
          IF (ILIIN(ISTS) == 0) THEN
            NINCX=NTBAR(IPOLGN,NRCELL)-NRCELL
          ELSE
            NINCX=SIGN(1,INMTIT(IPOLGN,NRCELL))
          END IF
        ELSE
          GOTO 9999
        ENDIF
#ifdef TRACE
        IF (NLTRC) WRITE (iunout,*) ' NRCELL,MRSURF,PT,NINCX,IPOLGN ',
     .                                NRCELL,MRSURF,PT,NINCX,IPOLGN
#endif
 
        RETURN
 
      ELSE
C
C  PARTICLE IS OUTSIDE STANDARD MESH, IN ADDITIONAL CELL NACELL
C  NRCELL=0
C
        IF (ITFRST == 0) THEN
CDR  PREPARE SOME GRID DATA FOR FINDING RE-ENTRY INTO GRID OF TETRAHEDA
CDR  ONLY ONCE FOR ALL FUTURE TRAJECTORIES.

C    SET THE FOLLOWING VARIABLES: .....  (COMMENTING TO BE DONE).   
          ITFRST = 1
          NTETSUR = COUNT(NTBAR(1:4,1:NTET) == 0)
          ALLOCATE (IDROB(NTETSUR,2))
          IOB = 0
          DO ITET = 1,NTET
            DO IS = 1,4
              IF (NTBAR(IS,ITET) == 0) THEN
                IOB = IOB + 1
                IDROB(IOB,:) = (/ ITET,IS /)
              END IF
            END DO
          END DO
 
          ALLOCATE (ISEEOB(0:NLIMI,NTETSUR/NBITS+1))
          ISEEOB = 0
          ISEEOB(0,:) = -1
          DO I=1,NTETSUR
            IT=IDROB(I,1)
            IS=IDROB(I,2)
            SIG1=SIGN(1,IRICH(1,IS))
            I1=ABS(IRICH(1,IS))
            SIG2=SIGN(1,IRICH(2,IS))
            I2=ABS(IRICH(2,IS))
            A1 = VTETX(I1,IT) * SIG1
            A2 = VTETY(I1,IT) * SIG1
            A3 = VTETZ(I1,IT) * SIG1
            B1 = VTETX(I2,IT) * SIG2
            B2 = VTETY(I2,IT) * SIG2
            B3 = VTETZ(I2,IT) * SIG2
            V1 = A2*B3 - A3*B2
            V2 = A3*B1 - A1*B3
            V3 = A1*B2 - A2*B1
CDR  LOOP OVER ADDITIONAL SURFACES
            DO IL=1,NLIMI
              IF (RLB(IL) == 3.) THEN
                S1=(P1(1,IL)-XTETRA(NTECK(ITSIDE(1,IS),IT)))*V1 +
     .             (P1(2,IL)-YTETRA(NTECK(ITSIDE(1,IS),IT)))*V2 +
     .             (P1(3,IL)-ZTETRA(NTECK(ITSIDE(1,IS),IT)))*V3
                S2=(P2(1,IL)-XTETRA(NTECK(ITSIDE(1,IS),IT)))*V1 +
     .             (P2(2,IL)-YTETRA(NTECK(ITSIDE(1,IS),IT)))*V2 +
     .             (P2(3,IL)-ZTETRA(NTECK(ITSIDE(1,IS),IT)))*V3
                S3=(P3(1,IL)-XTETRA(NTECK(ITSIDE(1,IS),IT)))*V1 +
     .             (P3(2,IL)-YTETRA(NTECK(ITSIDE(1,IS),IT)))*V2 +
     .             (P3(3,IL)-ZTETRA(NTECK(ITSIDE(1,IS),IT)))*V3
                IF (ANY( (/ S1,S2,S3 /) > 0.D0)) THEN
                  CALL EIRENE_BITSET (ISEEOB,0,NLIMI,IL,I,1,NBITS)
                ELSE
                  IF (NOPTIM >= NTET) THEN
                    IF (NLIMPB >= NLIMPS) THEN
                      IGJUM3(ITET,IL) = 0
                    ELSE
                      CALL EIRENE_BITSET
     .                    (IGJUM3,0,NOPTIM,ITET,IL,0,NBITS)
                    END IF
                  END IF
                END IF
              ELSE
                CALL EIRENE_BITSET (ISEEOB,0,NLIMI,IL,I,1,NBITS)
              END IF
            END DO
          END DO
        END IF   ! ITFRST == 0 PART. NOW: ITFRST = 1 
 
        PT=1.D30
        MRSURF=0
 
        TIMT=1.D30
        NTMZ=0
        NTMS=0
        DO IT=1,NTETSUR
          IF ((MRSURF == 0) .AND.
     .        .NOT.EIRENE_BITGET(ISEEOB,0,NLIMI,MSURF,IT,NBITS)) CYCLE
          I = IDROB(IT,1)
          J = IDROB(IT,2)
          IF (I == IZELLO)  CYCLE
          SIG1=SIGN(1,IRICH(1,J))
          I1=ABS(IRICH(1,J))
          SIG2=SIGN(1,IRICH(2,J))
          I2=ABS(IRICH(2,J))
          A(1:3,1) = (/ VTETX(I1,I), VTETY(I1,I),
     .                  VTETZ(I1,I)/) * SIG1
          A(1:3,2) = (/ VTETX(I2,I), VTETY(I2,I),
     .                  VTETZ(I2,I)/) * SIG2
          A(1:3,3) = (/ -VELX, -VELY, -VELZ /)
          B(1:3) = (/ X0-XTETRA(NTECK(ITSIDE(1,J),I)),
     .                Y0-YTETRA(NTECK(ITSIDE(1,J),I)),
     .                Z0-ZTETRA(NTECK(ITSIDE(1,J),I)) /)
          DET = EIRENE_SARRUS(A)
          IF (ABS(DET) < 1.D-5) CYCLE
          AB = A
          AB(:,1) = B
          XMU = EIRENE_SARRUS(AB)/DET
          AB = A
          AB(:,2) = B
          XETA = EIRENE_SARRUS(AB)/DET
          AB = A
          AB(:,3) = B
          T = EIRENE_SARRUS(AB)/DET
          IF ((XMU .GE.0.D0 .AND. XMU .LE.1.D0) .AND.
     .        (XETA.GE.0.D0 .AND. XETA.LE.1.D0) .AND.
     .        (XMU+XETA <= 1.D0) .AND. (T .GT. 0.)) THEN
            IF (T < TIMT) THEN
              TIMT = T
              NTMZ = I
              NTMS = J
            END IF
          END IF
        END DO
        PT = TIMT
        IZELLO = 0
        IPOLGOO = 0
        IF (NTMZ .NE. 0) THEN
C  SEITENNUMMER DES NEUEN DREIECKS
          IPOLGO = NTMS
#ifdef TRACE
          IF (NLTRC) THEN
            WRITE(iunout,*) 'OUTSIDE; GEHE IN ZELLE ',NTMZ,' TM= ',TM
            WRITE(iunout,*) 'DORT AUF SEITE IPOLGO ',IPOLGO
          ENDIF
#endif
        END IF
 
        NLSRFX=.FALSE.
        NJUMP=1
        MRSURF=NTMZ
        IPOLGN=NTMS
        NINCX=NTMZ-NRCELL
#ifdef TRACE
        IF (NLTRC) WRITE (iunout,*) ' NRCELL,MRSURF,PT,NINCX,IPOLGN ',
     .                                NRCELL,MRSURF,PT,NINCX,IPOLGN
#endif
      END IF

      RETURN
 
12000 CONTINUE
C
C  GENERAL GEOMETRY OPTION: PROVIDE FLIGHT TIME IN CURRENT CELL
C
      IF (ICALL == 0) THEN
        ICALL = 1
c  set nrmsrf:  = a "free virtual radial" surface index, between 1 and nsurf, such that none
c              of the block 3a defined nstsi "real radial" non-default standard surfaces
c              coincides with this "virtual radial" surface.

        MXSF = MAXVAL(INUMP(1:NSTSI,1))
        IF (MXSF < NSURF) THEN
          NRMSRF = MXSF+1
        ELSE
          DO IRS=1,NSURF
            IF (ALL(INUMP(1:NSTSI,1).NE.IRS)) EXIT
          END DO
          NRMSRF = IRS
          IF (NRMSRF > NSURF) THEN
            WRITE (iunout,*) ' ERROR IN TIMER, NRMSRF WRONG '
            CALL EIRENE_EXIT_OWN(1)
          END IF
        END IF
      END IF
C

      
C  NJUMP=3:  INTERNAL GRID SURFACE, STOP AND GO.
      IF (NJUMP.EQ.3.or.njump.eq.0) then
        TIM=0.
      endif

#ifdef TRACE
      IF (NLTRC) THEN
        WRITE (iunout,*) 'TIMER, LEVGEO=10, IN: NJUMP,TIM,ZT'
        WRITE (iunout,*)                        NJUMP,TIM,ZT
      ENDIF
#endif
      
c     if (nltrc.and.tim.ne.zt) write (6,*) 'error, npanu ',npanu,tim,zt

      CALL EIRENE_TIMUSR(NRCELL,X0,Y0,Z0,VELX,VELY,VELZ,NJUMP,
     .                   NEWCEL,TIM,ICOS,IERR,NPANU,NLSRFX)

      IF (IERR.NE.0) GOTO 9999
      PT=TIM
      IF (NEWCEL.GT.0) THEN
cdr  an intersection with the surface of cell nrcell has been found.
cdr  but this is not one of the "radial" non-default standard surfaces defined in 3a.
        NINCX=NEWCEL-NRCELL
cdr  make sure that mrsurf is not pointing to any of the defined non-default standard surfaces
cdr  but why not just mrsurf=0 ???
cdr  in cases levgeo ne 10, mrsurf is the next grid surface label, 
c    no matter if non-default (3a) or not.
        MRSURF=NRMSRF
c
      ELSEIF (NEWCEL.LT.0) THEN
cdr  one of the non-default surfaces has been hit. set this surface index to mrsurf.
        NINCX=ICOS
        MRSURF=INUMP(-NEWCEL,1)
cdr
      ELSEIF (NEWCEL.EQ.0) THEN
        WRITE (iunout,*) 'NEWCEL=0, EXIT FROM MESH '
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
      NLSRFX=.FALSE.

#ifdef TRACE
      IF (NLTRC) THEN
        WRITE (iunout,*) 'TIMER, OUT: PT,MRSURF,NINCX,NRCELL'
        WRITE (iunout,*)              PT,MRSURF,NINCX,NRCELL
      ENDIF
#endif

      RETURN
 
C
9999  CONTINUE
      WRITE (iunout,*) 'ERROR IN TIMER, EXIT CALLED AT 9999'
      WRITE (iunout,*) 'NPANU ', NPANU
      CALL EIRENE_EXIT_OWN(1)
      END
