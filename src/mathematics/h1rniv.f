*
*
      SUBROUTINE H1RNIV(VEC)
cdr  input:  vector (seed array) VEC(100), which fully describes a status of the
cdr          RANMAR generator.
cdr  output: Common RASET1, RASET2, such that random generator is re-initialized
cdr          according to status vector VEC(100)
*
      USE EIRMOD_PRECISION
C     USE EIRMOD_COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: VEC(100)
      INTEGER :: IC
      CHARACTER*16    FLAG
      REAL(DP) :: U, C, CD, CM
      INTEGER :: IP, JP
      COMMON /RASET1/ U(97),C,CD,CM,IP,JP
      COMMON /RASET2/ FLAG
*
      DO 400 IC = 1, 97
        U(IC) = VEC(IC)
  400 CONTINUE
      C  = VEC(98)
      CD =  7654321./16777216.
      CM = 16777213./16777216.
      IP = NINT(VEC(99))
      JP = NINT(VEC(100))
*
      FLAG = 'H1RN INITIALISED'
c     WRITE(IUNOUT,*)
c    .           ' H1RNIV: H1RN (RANMAR) INITIALISED/RESTARTED WITH',
c    >           ' SEED ARRAY VEC(100)'
*
      RETURN
      END
