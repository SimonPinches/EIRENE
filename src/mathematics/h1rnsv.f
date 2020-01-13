cdr  Fetch (and return) a seed VEC(100) from current RANMAR status,
cdr  to continue the random number sequence later at this point.
cdr  VEC(100) is obtained from Common RASET1
cdr
cdr  Corresponds to RANGET (fetch a seed of random generator, for current status).
cdr  However, for the RANMAR generator a single integer seed IJKL 
cdr  is not sufficient for all possible states of the generator. 
cdr  The full vector VEC is needed instead.
*
      SUBROUTINE H1RNSV(VEC)
*
      USE EIRMOD_PRECISION
C     USE EIRMOD_COMPRT, ONLY: IUNOUT

      IMPLICIT NONE
      REAL(DP), INTENT(OUT) :: VEC(100)
      INTEGER :: IC
      REAL(DP) :: U, C, CD, CM
      INTEGER :: I, J
      COMMON /RASET1/ U(97),C,CD,CM,I,J
*
      DO 10 IC = 1, 97
        VEC(IC) = U(IC)
   10 CONTINUE
      VEC(98) = C
      VEC(99) = REAL(I)
      VEC(100)= REAL(J)
c     WRITE(IUNOUT,*)
c    .           ' H1RNSV: H1RN (RANMAR) STATUS SAVED',
c    >           ' ON SEED ARRAY VEC(100)'

      RETURN
      END
