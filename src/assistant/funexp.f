CDR,  Aug. 3rd 2017.
C
C*DK FUNEXP
      FUNCTION EIRENE_FUNEXP(X,EX)
C
C  STANDARD MONTE CARLO "SPANIER-ESTIMATOR" FUNCTION

C  FUNEXP = (1-EXP(-X))/X  =  (EXP(-X)-1)/(-X) = (EXP(Y)-1)/Y
C  WITH X >=0.,  OR, EQUVALIENTLY:  -X = Y <=0.

c  return s=funexp, as well as intermediate result ex=exp(-x)
c                   or find ex by separate evaluation in intrinsic exp -function ?
c  note: in a Taylor expansion the division 1/x cancels.
c  note: ffexp= exp(x)-1  may already be an existing intrinsic fortran function.
C
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: X
      REAL(DP), INTENT(OUT) :: EX
      REAL(DP) :: EIRENE_FUNEXP, S

      IF (X.LE.1.D-10) THEN
        EX=1.
        S=1.
      ELSEIF (X.GT.1.D2) THEN
        EX=0.D0
        S=1/X
      ELSE
        EX=EXP(-X)
        S=(1.-EX)/X
      ENDIF
c

      EIRENE_FUNEXP=S

      RETURN
      END
