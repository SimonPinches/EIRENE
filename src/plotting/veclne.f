C   may 15:  plot all triangle sides, not only those with no neighbor.
c            changed by adding .true.or. ...  in "itr" loop for levgeo=4
C
      SUBROUTINE EIRENE_VECLNE (AORIG,BORIG,IBLD,ICURV,
     .                   IXXI,IXXE,IYYI,IYYE,
     .                   TEXT1,TEXT2,TEXT3,
     .                   LOGL,ZMA,ZMI,
     .                   HEAD,RUNID,TXHEAD,TRC)
C
C  THIS SUBROUTINE PRODUCES A VECTOR FIELD PLOT
C  MODIFIED 160600: LEVGEO=3 OPTION ADDED
C                   SCALING OF ARROWS AUTOMATICALLY
C  ARGUMENTS: XX,YY OUT, IXX,IYY --> IXXI,IXXE,IYYI,IYYE
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_CLOGAU
      USE EIRMOD_CPLOT
      USE EIRMOD_CPOLYG
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      USE EIRMOD_CTRIG
      IMPLICIT NONE
C
C
      REAL(DP), INTENT(IN) :: AORIG(*), BORIG(*)
      REAL(DP), INTENT(IN) :: ZMA, ZMI
      INTEGER, INTENT(IN) :: IBLD, ICURV
      INTEGER, INTENT(INOUT) :: IXXI, IXXE, IYYI, IYYE
      LOGICAL, INTENT(IN) :: LOGL, TRC
      CHARACTER(72), INTENT(IN) :: HEAD, RUNID, TXHEAD
      CHARACTER(72), INTENT(OUT) :: TEXT1
      CHARACTER(24), INTENT(IN) :: TEXT2, TEXT3

      REAL(DP) :: SCLFCX, SCLFCY, VMIN, VMAX, DX, DY, CM, FAK, BRFL,
     .          D, DIAMETER, XM, YM, VX, VY, VABS, XMIN, XMAX, YMIN,
     .          YMAX, PLFL
      REAL(SP) YH
      INTEGER :: ITR, IRAD, I, IR, IPART, NP1, NP2, IP
      CHARACTER*17 CH
C
C  SEARCH FOR XMIN,XMAX,YMIN,YMAX
C
      XMIN=1.D60
      YMIN=1.D60
      XMAX=-1.D60
      YMAX=-1.D60
C
      IF ((LEVGEO.EQ.1.OR.LEVGEO.EQ.2).AND.NLTOR) THEN
        IXXI=MAX(1,IXXI)
        IXXE=MIN(IXXE,NR1ST)
        IF (IXXE.LE.0) IXXE=NR1ST
        IYYI=MAX(1,IYYI)
        IYYE=MIN(IYYE,NT3RD)
        IF (IYYE.LE.0) IYYE=NT3RD
        XMIN = RHOSRF(IXXI)
        XMAX = RHOSRF(IXXE)
        YMIN = ZSURF(IYYI)
        YMAX = ZSURF(IYYE)
      ELSEIF (LEVGEO.EQ.1.AND.NLPOL) THEN
        IXXI=MAX(1,IXXI)
        IXXE=MIN(IXXE,NR1ST)
        IF (IXXE.LE.0) IXXE=NR1ST
        IYYI=MAX(1,IYYI)
        IYYE=MIN(IYYE,NP2ND)
        IF (IYYE.LE.0) IYYE=NP2ND
        XMIN = RHOSRF(IXXI)
        XMAX = RHOSRF(IXXE)
        YMIN = PSURF(IYYI)
        YMAX = PSURF(IYYE)
      ELSEIF (LEVGEO.EQ.2.AND.NLPOL) THEN
C  SUFFICIENT TO SEARCH ON OUTERMOST RADIAL SURFACE (BECAUSE: CONVEX)
        IXXI=MAX(1,IXXI)
        IXXE=MIN(IXXE,NR1ST)
        IF (IXXE.LE.0) IXXE=NR1ST
        IYYI=MAX(1,IYYI)
        IYYE=MIN(IYYE,NP2ND)
        IF (IYYE.LE.0) IYYE=NP2ND
        DO 5 IP = IYYI,IYYE
          XMIN = MIN(XMIN,XPOL(IXXE,IP))
          XMAX = MAX(XMAX,XPOL(IXXE,IP))
          YMIN = MIN(YMIN,YPOL(IXXE,IP))
          YMAX = MAX(YMAX,YPOL(IXXE,IP))
    5   CONTINUE
      ELSEIF (LEVGEO.EQ.3.AND.NLPOL) THEN
C  SEARCH ON WHOLE MESH
        IXXI=MAX(1,IXXI)
        IXXE=MIN(IXXE,NR1ST)
        IF (IXXE.LE.0) IXXE=NR1ST
        IYYI=MAX(1,IYYI)
        IYYE=MIN(IYYE,NP2ND)
        IF (IYYE.LE.0) IYYE=NP2ND
        DO 1 IR=IXXI,IXXE
          NP1=MAX(IYYI,NPOINT(1,1))
          NP2=MIN(IYYE,NPOINT(2,NPPLG))
          XMIN=MIN(XMIN,XPOL(IR,NP1),XPOL(IR,NP2))
          YMIN=MIN(YMIN,YPOL(IR,NP1),YPOL(IR,NP2))
          XMAX=MAX(XMAX,XPOL(IR,NP1),XPOL(IR,NP2))
          YMAX=MAX(YMAX,YPOL(IR,NP1),YPOL(IR,NP2))
    1   CONTINUE
C
        DO 4 IPART=1,NPPLG
          DO 2 IP=NPOINT(1,IPART),NPOINT(2,IPART)
            IF (IP.LT.IYYI.OR.IP.GT.IYYE) GOTO 2
            XMIN=MIN(XMIN,XPOL(IXXI,IP))
            YMIN=MIN(YMIN,YPOL(IXXI,IP))
            XMAX=MAX(XMAX,XPOL(IXXI,IP))
            YMAX=MAX(YMAX,YPOL(IXXI,IP))
    2     CONTINUE
          DO 3 IP=NPOINT(1,IPART),NPOINT(2,IPART)
            IF (IP.LT.IYYI.OR.IP.GT.IYYE) GOTO 3
            XMIN=MIN(XMIN,XPOL(IXXE,IP))
            YMIN=MIN(YMIN,YPOL(IXXE,IP))
            XMAX=MAX(XMAX,XPOL(IXXE,IP))
            YMAX=MAX(YMAX,YPOL(IXXE,IP))
    3     CONTINUE
    4   CONTINUE
      ELSEIF (LEVGEO.EQ.4) THEN
C  SEARCH ON WHOLE MESH
        DO I=1,NRKNOT
          XMIN=MIN(XMIN,XTRIAN(I))
          YMIN=MIN(YMIN,YTRIAN(I))
          XMAX=MAX(XMAX,XTRIAN(I))
          YMAX=MAX(YMAX,YTRIAN(I))
        ENDDO
      ENDIF
C
C
      CM=20.
      DX=FCABS1(IBLD) * (XMAX-XMIN)
      DY=FCABS2(IBLD) * (YMAX-YMIN)
      FAK=CM/MAX(DX,DY)
C
C  PLOT FRAME
C
      CALL EIRENE_PLNXTB (1,'VECLNE.F')

cdr  use FZJ proprietary GR plot software
      CALL GRSCLC (10.,4.,REAL(10.+DX*FAK,SP),
     .             REAL(4.+DY*FAK,SP))
      CALL GRSCLV (REAL(XMIN,SP),REAL(YMIN,SP),
     .             REAL(XMAX,SP),REAL(YMAX,SP))
      CALL GRAXS (7,'X=3,Y=3',6,'R (CM)',6,'Z (CM)')
C
C  SCALE FACTORS: USER COORDINATES TO CM:
C  X-DIRECTION:
      SCLFCX=((10.+DX*FAK)-10.)/(XMAX-XMIN)
C  Y-DIRECTION:
      SCLFCY=((4.+DY*FAK)-4.)/(YMAX-YMIN)
C
C  PLOT BOUNDARY OF MESH
C
      IF ((LEVGEO.EQ.1.OR.LEVGEO.EQ.2).AND.NLTOR) THEN
        CALL GRJMP(REAL(XMIN,SP),REAL(YMIN,SP))
        CALL GRDRW(REAL(XMIN,SP),REAL(YMAX,SP))
        CALL GRDRW(REAL(XMAX,SP),REAL(YMAX,SP))
        CALL GRDRW(REAL(XMAX,SP),REAL(YMIN,SP))
        CALL GRDRW(REAL(XMIN,SP),REAL(YMIN,SP))
      ELSEIF (LEVGEO.EQ.1.AND.NLPOL) THEN
        CALL GRJMP(REAL(XMIN,SP),REAL(YMIN,SP))
        CALL GRDRW(REAL(XMIN,SP),REAL(YMAX,SP))
        CALL GRDRW(REAL(XMAX,SP),REAL(YMAX,SP))
        CALL GRDRW(REAL(XMAX,SP),REAL(YMIN,SP))
        CALL GRDRW(REAL(XMIN,SP),REAL(YMIN,SP))
      ELSEIF (LEVGEO.EQ.2.AND.NLPOL) THEN
        DO 7 IR=IXXI,IXXE,IXXE-IXXI
          CALL GRJMP(REAL(XPOL(IR,IYYI),SP),
     .               REAL(YPOL(IR,IYYI),SP))
          DO IP = IYYI+1,IYYE
            CALL GRDRW(REAL(XPOL(IR,IP),SP),
     .                 REAL(YPOL(IR,IP),SP))
          END DO
    7   CONTINUE
      ELSEIF (LEVGEO.EQ.3.AND.NLPOL) THEN
        NP1=MAX(IYYI,NPOINT(1,1))
        NP2=MIN(IYYE,NPOINT(2,NPPLG))
        CALL GRJMP(REAL(XPOL(IXXI,NP1),SP),
     .             REAL(YPOL(IXXI,NP1),SP))
          DO 10 IR=IXXI+1,IXXE
            CALL GRDRW (REAL(XPOL(IR,NP1),SP),
     .                  REAL(YPOL(IR,NP1),SP))
   10     CONTINUE
C
        CALL GRJMP (REAL(XPOL(IXXI,NP2),SP),
     .              REAL(YPOL(IXXI,NP2),SP))
        DO 11 IR=IXXI+1,IXXE
          CALL GRDRW (REAL(XPOL(IR,NP2),SP),
     .                REAL(YPOL(IR,NP2),SP))
   11   CONTINUE
        DO 15 I=1,NPPLG
          NP1=MAX(IYYI,NPOINT(1,I))
          NP2=MIN(IYYE,NPOINT(2,I))
          CALL GRJMP (REAL(XPOL(IXXI,NP1),SP),
     .                REAL(YPOL(IXXI,NP1),SP))
          DO 12 IP=NP1,NP2
            CALL GRDRW (REAL(XPOL(IXXI,IP),SP),
     .                  REAL(YPOL(IXXI,IP),SP))
   12     CONTINUE
          CALL GRJMP (REAL(XPOL(IXXE,NPOINT(1,I)),SP),
     .                REAL(YPOL(IXXE,NPOINT(1,I)),SP))
          DO 13 IP=NP1,NP2
            CALL GRDRW (REAL(XPOL(IXXE,IP),SP),
     .                  REAL(YPOL(IXXE,IP),SP))
   13     CONTINUE
   15   CONTINUE
      ELSEIF (LEVGEO.EQ.4) THEN
        DO ITR=1,NTRII
c   .true.or. .... added, to plot all triangles, not only those with no neighbor
          IF (.true..or.NCHBAR(1,ITR) .EQ. 0) THEN
            CALL GRJMP(REAL(XTRIAN(NECKE(1,ITR)),SP),
     .                 REAL(YTRIAN(NECKE(1,ITR)),SP))
            CALL GRDRW(REAL(XTRIAN(NECKE(2,ITR)),SP),
     .                 REAL(YTRIAN(NECKE(2,ITR)),SP))
          ENDIF
          IF (.true..or.NCHBAR(2,ITR) .EQ. 0) THEN
            CALL GRJMP(REAL(XTRIAN(NECKE(3,ITR)),SP),
     .                 REAL(YTRIAN(NECKE(3,ITR)),SP))
            CALL GRDRW(REAL(XTRIAN(NECKE(2,ITR)),SP),
     .                 REAL(YTRIAN(NECKE(2,ITR)),SP))
          ENDIF
          IF (.true..or.NCHBAR(3,ITR) .EQ. 0) THEN
            CALL GRJMP(REAL(XTRIAN(NECKE(1,ITR)),SP),
     .                 REAL(YTRIAN(NECKE(1,ITR)),SP))
            CALL GRDRW(REAL(XTRIAN(NECKE(3,ITR)),SP),
     .                 REAL(YTRIAN(NECKE(3,ITR)),SP))
          ENDIF
        ENDDO
      ENDIF
C
C PLOT 2D VECTOR FIELD
C
      IF (((LEVGEO.EQ.1.OR.LEVGEO.EQ.2).AND.NLTOR).OR.
     .     (LEVGEO.EQ.1.AND.NLPOL)) THEN
C
C 1. SEARCH FOR MINIMA AND MAXIMA
          VMIN=1.D60
          VMAX=-1.D60
          DO 900 IR=IXXI,IXXE-1
           DO IP=IYYI,IYYE-1
            IRAD = IR + (IP-1)*NR1ST
            VABS = SQRT(AORIG(IRAD)**2 + BORIG(IRAD)**2)
            VMIN = MIN(VMIN,VABS)
            VMAX = MAX(VMAX,VABS)
           END DO
  900     CONTINUE
C 2. SCALING
          DO 1100 IR=IXXI,IXXE-1
           DO IP=IXXI,IYYE-1
            IRAD = IR + (IP-1)*NR1ST
            VX = (AORIG(IRAD) / VMAX)
            VY = (BORIG(IRAD) / VMAX)
C 3. PLOT VECTOR
            CALL GRNWPN(3)
            XM=XCOM(IRAD)-VX/2
            YM=YCOM(IRAD)-VY/2
            plfl=sqrt(vx*vx+vy*vy)/5.*sclfcx
            brfl=sqrt(vx*vx+vy*vy)/7.5*sclfcx
            call grarrw(REAL(xm,SP),REAL(ym,SP),
     .                  REAL(xm+vx,SP),REAL(ym+vy,SP),
     .                  REAL(PLFL,SP),REAL(BRFL,SP),1)
           END DO
 1100     CONTINUE
C
      ELSEIF ((LEVGEO.EQ.2.AND.NLPOL).OR.
     .         LEVGEO.EQ.3) THEN
C
C 1. SEARCH FOR MINIMA AND MAXIMA
          VMIN= 1.D60
          VMAX=-1.D60
          DIAMETER=1.D60
          DO 901 IR=IXXI,IXXE-1
           DO IP=IXXI,IYYE-1
            IRAD = IR + (IP-1)*NR1ST
            D=SQRT((XPOL(IR,IP)-XPOL(IR+1,IP+1))**2+
     .             (YPOL(IR,IP)-YPOL(IR+1,IP+1))**2)
            DIAMETER=MIN(DIAMETER,D)
            VABS = SQRT(AORIG(IRAD)**2 + BORIG(IRAD)**2)
            VMIN = MIN(VMIN,VABS)
            VMAX = MAX(VMAX,VABS)
C 2. SCALING
            VX =  AORIG(IRAD) / VABS * D * 0.5
            VY =  BORIG(IRAD) / VABS * D * 0.5
C 3. PLOT VECTOR
            CALL GRNWPN(3)
            XM=XCOM(IRAD)-VX/2
            YM=YCOM(IRAD)-VY/2
            plfl=sqrt(vx*vx+vy*vy)/5.*sclfcx
            brfl=sqrt(vx*vx+vy*vy)/7.5*sclfcx
            call grarrw(REAL(xm,SP),REAL(ym,SP),
     .                  REAL(xm+vx,SP),REAL(ym+vy,SP),
     .                  REAL(PLFL,SP),REAL(BRFL,SP),1)
           END DO
  901     CONTINUE
C
      ELSEIF (LEVGEO.EQ.4) THEN
C
C 1. SEARCH FOR MINIMA AND MAXIMA
          VMIN=1.D60
          VMAX=-1.D60
          DO 1500 IR=1,NTRII
            VABS = SQRT(AORIG(IR)**2 + BORIG(IR)**2)
            VMIN = MIN(VMIN,VABS)
            VMAX = MAX(VMAX,VABS)
 1500     CONTINUE
C 2. SCALING
          DO 1700 IR=1,NTRII
            VX = 10 * (AORIG(IR) / VMAX)
            VY = 10 * (BORIG(IR) / VMAX)
C 3. PLOT VECTOR
            CALL GRNWPN(3)
            XM=XCOM(IR)-VX/2
            YM=YCOM(IR)-VY/2
            plfl=sqrt(vx*vx+vy*vy)/5.*sclfcx
            brfl=sqrt(vx*vx+vy*vy)/7.5*sclfcx
            call grarrw(REAL(xm,SP),REAL(ym,SP),
     .                  REAL(xm+vx,SP),REAL(ym+vy,SP),
     .                  REAL(PLFL,SP),REAL(BRFL,SP),1)
 1700     CONTINUE
      ENDIF
C
 2000 CONTINUE
C
C     WRITE TEXT, MAXIMUM AND MINIMUM VALUE ONTO THE PLOT
C
      TEXT1='2D VECTOR FIELD'
      CALL GRNWPN (1)
      CALL GRSCLC (0.,0.,39.,28.)
      CALL GRSCLV (0.,0.,39.,28.)
      YH=27.5
      CALL GRTXT (1.,REAL(YH,SP),72,RUNID)
      YH=26.75
      CALL GRTXT (1.,REAL(YH,SP),72,HEAD)
      YH=26.00
      CALL GRTXT (1.,REAL(YH,SP),72,TXHEAD)
      YH=25.25
      CALL GRTXT (1.,REAL(YH,SP),10,'TALLY :  ')
      CALL GRTXTC (72,TEXT1)
      CALL GRTXT (1.,REAL(YH-0.5,SP),10,'SPECIES :')
      CALL GRTXTC (24,TEXT2)
      CALL GRTXT (1.,REAL(YH-1.,SP),10,'UNITS :   ')
      CALL GRTXTC (24,TEXT3)
      CALL GRTXT (1.,REAL(YH-2.,SP),10,'MAX. VALUE')
      WRITE (CH,'(1P,E10.3)') VMAX
      CALL GRTXT (1.,REAL(YH-2.5,SP),10,CH)
      CALL GRTXT (1.,REAL(YH-3.,SP),10,'MIN. VALUE')
      WRITE (CH,'(1P,E10.3)') VMIN
      CALL GRTXT (1.,REAL(YH-3.5,SP),10,CH)
C
      RETURN
      END
