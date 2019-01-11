Cdr Jan 2017:  start to synchronize with inttal. goal: should
C              become identical, if weighting function "VOL == 1"
C   some unused variables (VR, VRTP,....) intentionally retained here,
C   to facilitate the future combination of intvol and inttal into a single
C   routine.
C
C*DK INTVOL
      SUBROUTINE EIRENE_INTVOL (A,J,M,N,YINT,NX,NY,NZ,NB)
C
C   SIMILAR TO INTTAL
C   SUM UP TALLY A(J,K), K=1,N, RESULT IS YINT
C   USE 1. AS WEIGHTING. PUT TOTALS RATHER THAN MEANS ON AVERAGE-TALLY
C   LOCATIONS

C   J FIXED, (SPECIES INDEX)
C   M: LEADING DIMENSION OF A IN CALLING PROGRAM
C   A IS A TOTAL TALLY (EG. CELL VOLUME)
C
C  IN CASE NX*NY*NZ*NB NOT ZERO, A STANDARD MESH HAS BEEN DEFINED.
C                                THEN THE VARIOUS PROJECTIONS
C                                ONTO THE 1ST,2ND AND 3RD MESH ARE
C                                COMPUTED
C                                IX=1,NX-1, IY=1,NY, IZ=1,NZ-1
C
      USE EIRMOD_PRECISION
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: J, M, N, NX, NY, NZ, NB
      REAL(DP), INTENT(INOUT) :: A(M,*)
C     REAL(DP), INTENT(IN) :: VOL(*)
      REAL(DP), INTENT(OUT) :: YINT
      REAL(DP) :: VR, YR, YPR, VRTP, VRT, YRT, VTP, YTP, VT, YT, YP,
     .            VPR, VP, YRTP
      INTEGER :: IX, IY, IZ, IB, NIR, NIRTP, NIRT, NITP, NIT, NIPR, NIP,
     .           N1DEL, N2DEL, KX, KY, KZ, K, NXM, NYM, NZM, IADD,
     .           KB, NS
C
      N1DEL=0
      IF (NY.GT.1.OR.NZ.GT.1) N1DEL=NX
      N2DEL=0
      IF (NZ.GT.1) N2DEL=NY
      NXM=MAX(1,NX-1)
      NYM=MAX(1,NY-1)
      NZM=MAX(1,NZ-1)
      NS=NX*NY*NZ*NB
C
C  LOOP OVER ALL CELLS
      YINT=0.
      DO 5 KB=1,NB
        IADD=(KB-1)*NX*NY*NZ
        DO KX=1,NXM
        DO KY=1,NYM
        DO KZ=1,NZM
          K=KX+((KY-1)+(KZ-1)*N2DEL)*N1DEL+IADD
          YINT=YINT+A(J,K)
        END DO
        END DO
        END DO
    5 CONTINUE
      DO 6 K=NS+1,N
        YINT=YINT+A(J,K)
    6 CONTINUE
C
C  LOOP OVER STANDARD MESH BLOCKS
      IF (NS.EQ.0) RETURN
C
      DO 1000 IB=1,NB
      IADD=(IB-1)*NX*NY*NZ
C
C  INTEGRATE OVER RADIAL COORDINATE: YR
C  INTEGRATE OVER RADIAL AND TOROIDAL COORDINATE: YRT
C  INTEGRATE OVER ALL THREE COORDINATES: YRTP
      YRTP=0.
      VRTP=1.D-60
      NIRTP=NX+((NY-1)+(NZ-1)*N2DEL)*N1DEL+IADD
      DO 100 IY=1,NYM
        YRT=0.
        VRT=1.D-60
        NIRT=NX+((IY-1)+(NZ-1)*N2DEL)*N1DEL+IADD
        DO 101 IZ=1,NZM
          YR=0.
          VR=1.D-60
          NIR=NX+((IY-1)+(IZ-1)*N2DEL)*N1DEL+IADD
          DO 102 IX=1,NXM
            K=IX+((IY-1)+(IZ-1)*N2DEL)*N1DEL+IADD
            YR=YR+A(J,K)
            YRT=YRT+A(J,K)
            YRTP=YRTP+A(J,K)
  102     CONTINUE
          A(J,NIR)=YR
  101   CONTINUE
        A(J,NIRT)=YRT
  100 CONTINUE
      A(J,NIRTP)=YRTP
C  INTEGRATE OVER POLOIDAL COORDINATE: YP
C  INTEGRATE OVER POLOIDAL AND RADIAL COORDINATE: YPR
      DO 200 IZ=1,NZM
        YPR=0.
        VPR=1.D-60
        NIPR=NX+((NY-1)+(IZ-1)*N2DEL)*N1DEL+IADD
        DO 201 IX=1,NXM
          YP=0.
          VP=1.D-60
          NIP=IX+((NY-1)+(IZ-1)*N2DEL)*N1DEL+IADD
          DO 202 IY=1,NYM
            K=IX+((IY-1)+(IZ-1)*N2DEL)*N1DEL+IADD
            YP=YP+A(J,K)
            YPR=YPR+A(J,K)
  202     CONTINUE
          A(J,NIP)=YP
  201   CONTINUE
        A(J,NIPR)=YPR
  200 CONTINUE
C  INTEGRATE OVER TOROIDAL COORDINATE: YT
C  INTEGRATE OVER TOROIDAL AND POLOIDAL COORDINATE: YTP
      DO 300 IX=1,NXM
        YTP=0.
        VTP=1.D-60
        NITP=IX+((NY-1)+(NZ-1)*N2DEL)*N1DEL+IADD
        DO 301 IY=1,NYM
          YT=0.
          VT=1.D-60
          NIT=IX+((IY-1)+(NZ-1)*N2DEL)*N1DEL+IADD
          DO 302 IZ=1,NZM
            K=IX+((IY-1)+(IZ-1)*N2DEL)*N1DEL+IADD
            YT=YT+A(J,K)
            YTP=YTP+A(J,K)
  302     CONTINUE
          A(J,NIT)=YT
  301   CONTINUE
        A(J,NITP)=YTP
  300 CONTINUE
 1000 CONTINUE
      RETURN
      END
