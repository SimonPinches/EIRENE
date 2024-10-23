C
C
C*DK SLOPE
      FUNCTION EIRENE_SLOPE(KN,JCH,NN,EN,BU,N,M)
cdr  simple linear regression, least square estimators

C
C     CALCULATE SLOPE OF STRAIGHT LINE
C
      USE EIRMOD_PRECISION
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: KN, JCH, NN, N, M
      REAL(DP), INTENT(IN) :: EN(M), BU(N,M)
      REAL(DP) :: ZSUM1, ZSUM2, ZSUM3, ZSUM4, ZX, ZY, FKN, EIRENE_SLOPE
      INTEGER :: J
C
      ZSUM1=0.0
      ZSUM2=0.0
      ZSUM3=0.0
      ZSUM4=0.0
C
      FKN=KN
      DO 101 J=1,KN
        ZX=EN(NN+J-1)
        ZY=BU(JCH,NN+J-1)
C
        ZSUM1=ZSUM1+ZX*ZY    !     sum (x_j*y_j)
        ZSUM2=ZSUM2+ZX       !KN * X_mean
        ZSUM3=ZSUM3+ZY       !KN * Y_mean
        ZSUM4=ZSUM4+ZX*ZX    !     sum (X_j**2)
  101 CONTINUE
C
      EIRENE_SLOPE=(FKN*ZSUM1-ZSUM2*ZSUM3)/(FKN*ZSUM4-ZSUM2*ZSUM2)
C
      RETURN
      END
