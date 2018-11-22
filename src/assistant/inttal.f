Cdr Jan 2017:  start to synchronize with intvol. goal: should
C              become identical, if weighting function "VOL == 1"
C
C
C
C*DK INTTAL
      SUBROUTINE EIRENE_INTTAL (A,VOL,J,M,N,YINT,NX,NY,NZ,NB)
C
C  SIMILAR TO INTVOL
C   INTEGRATE TALLY A(J,K), K=1,N, RESULT IS YINT
C   USE VOL(K) AS WEIGHTING
C
C   J FIXED, (SPECIES INDEX)
C   M: LEADING DIMENSION OF A IN CALLING PROGRAM
C   A IS A VOLUME-AVERAGED TALLY
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
      REAL(DP), INTENT(IN) :: VOL(*)
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
          YINT=YINT+VOL(K)*A(J,K)
        END DO
        END DO
        END DO
    5 CONTINUE
      DO 6 K=NS+1,N
        YINT=YINT+VOL(K)*A(J,K)
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
            YR=YR+VOL(K)*A(J,K)
            YRT=YRT+VOL(K)*A(J,K)
            YRTP=YRTP+VOL(K)*A(J,K)
            VR=VR+VOL(K)
            VRT=VRT+VOL(K)
            VRTP=VRTP+VOL(K)
  102     CONTINUE
          A(J,NIR)=YR/VR
  101   CONTINUE
        A(J,NIRT)=YRT/VRT
  100 CONTINUE
      A(J,NIRTP)=YRTP/VRTP
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
            YP=YP+VOL(K)*A(J,K)
            YPR=YPR+VOL(K)*A(J,K)
            VP=VP+VOL(K)
            VPR=VPR+VOL(K)
  202     CONTINUE
          A(J,NIP)=YP/VP
  201   CONTINUE
        A(J,NIPR)=YPR/VPR
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
            YT=YT+VOL(K)*A(J,K)
            YTP=YTP+VOL(K)*A(J,K)
            VT=VT+VOL(K)
            VTP=VTP+VOL(K)
  302     CONTINUE
          A(J,NIT)=YT/VT
  301   CONTINUE
        A(J,NITP)=YTP/VTP
  300 CONTINUE
 1000 CONTINUE
      RETURN
      END
