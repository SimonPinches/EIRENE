C
C
C*DK SLOPE
      FUNCTION EIRENE_SLOPE(KN,JCH,NN,EN,BU,N,M,BU0,RQ)
cdr  simple linear regression, least square estimators
cdr Nov. 21: added: y axis intersection BU0
cdr          Coefficient of determination RQ
cdr          Routine moved from diagno to mathematics

C
C     CALCULATE SLOPE OF STRAIGHT LINE, fit to channel JCH: BU(JCH,J) vs. EN(J)

cdr EN     array of independent variables vector, length M, (=regressor = predictor variables)
cdr BU     array of up to N dependent variable vectors, each of length M
cdr M      dimension of independent variable array  EN (EN(M)) as in calling program
cdr N      first dimension ("channel") of array BU  (BU(N,M)) as in calling program

cdr NN     no. of first data point used. Fit range: NN to NN+KN-1
cdr KN     data points used from NN to NN+KN-1, NN+KN-1 must be less or equal M
cdr JCH    channel number in BU field. Must be less or equal N

C
cdr   USE EIRMOD_PRECISION
      IMPLICIT NONE
      INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(13)

      INTEGER, INTENT(IN) :: KN, JCH, NN, N, M
      REAL(DP), INTENT(IN) :: EN(M), BU(N,M)
      REAL(DP) :: ZSUM1, ZSUM2, ZSUM3, ZSUM4, ZSUM5, ZX, ZY, FKN,
     .            XYNOM, DENOMX, DENOMY, SLOPE
      REAL(DP), INTENT(OUT) :: BU0, RQ
      REAL(DP) :: EIRENE_SLOPE
      INTEGER :: J
C
      ZSUM1=0.0
      ZSUM2=0.0
      ZSUM3=0.0
      ZSUM4=0.0
      ZSUM5=0.0
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
        ZSUM5=ZSUM5+ZY*ZY    !     sum (Y_j**2)
  101 CONTINUE
C
cdr  simple linear least square estimate
      XYNOM= FKN*ZSUM1-ZSUM2*ZSUM3
      DENOMX=FKN*ZSUM4-ZSUM2*ZSUM2+1.0d-30
      DENOMY=FKN*ZSUM5-ZSUM3*ZSUM3+1.0d-30

      SLOPE=XYNOM/DENOMX
cdr  y axis intercept least square estimate
      BU0         =(ZSUM3-SLOPE * ZSUM2)/FKN
cdr RQ: coefficient of determination (between 0 and 1)
      RQ=XYNOM*XYNOM/(DENOMX*DENOMY)

      EIRENE_SLOPE=SLOPE
C
      RETURN
      END FUNCTION EIRENE_SLOPE
