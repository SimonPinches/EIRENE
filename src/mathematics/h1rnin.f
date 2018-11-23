cdr  Initializes random number generator H1RN  ( = RANMAR, F. James, see below)
*
cdr   SUBROUTINE H1RNIN(IJ,KL)
cdr   IMPLICIT NONE
cdr   INTEGER, INTENT(IN,OUT) :: IJ, KL
cdr Version in eirene until April 2017: take 2 input seeds,
cdr and produce 4 smaller input seeds from them.

cdr April 2017: back to original suggestions by F. James
cdr Take only one seed IJKL, then first make 2 smaller (IJ, KL),
cdr then 4 yet smaller  (I,J,K,L) from it.
cdr IJKL must be not larger then 900.000.000, each seed then produces an
cdr      independent random sequence of length
cdr
      SUBROUTINE H1RNIN(IJKL)
cdr F. James, "A review of pseudorandom number generators",
cdr           Comp. Physics Communications 60 (1990) 329 - 344
cdr Subroutine RMARIN, p340.
cdr legal seed range 0<= IJKL<= 900.000.000 already enforced in calling
cdr routine
*
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IJKL
      INTEGER :: L, II, J, K, JJ, M, I, IJ, KL
      REAL(DP) :: S, T
cdr  next 5 lines: status of generator, initialized with IJKL
      CHARACTER*16    FLAG
      REAL(DP) :: U, C, CD, CM
      INTEGER :: IP, JP
      COMMON /RASET1/ U(97),C,CD,CM,IP,JP
      COMMON /RASET2/ FLAG
*
cdr   IJ = IABS(IJ)
cdr   KL = IABS(KL)

cdr   IJ = MOD(IJ,31329)   I do not know where this comes from.
cdr                        change back to original generator from
cdr                        F. James, CPC 60 (1990) 329 - 344, page 340

      IJ = IJKL/30082
      KL = IJKL-30082*IJ    ! = MOD(IJKL,30082)
c
      I  = MOD(IJ/177, 177) + 2
      J  = MOD(IJ, 177)     + 2
      K  = MOD(KL/169, 178) + 1
      L  = MOD(KL, 169)
*
cdr  now we have the 4 small seeds.
cdr  next: initialize RANMAR, data are transfered via common RASET1

      DO 300 II= 1, 97
         S= 0.
         T= 0.5
         DO 250 JJ= 1,24
            M = MOD(MOD(I*J,179)*K, 179)
            I = J
            J = K
            K = M
            L = MOD(53*L+1, 169)
            IF ( MOD(L*M,64) .GE. 32) S = S + T
            T = 0.5*T
  250    CONTINUE
         U(II) = S
  300 CONTINUE

      C  =   362436./16777216.
      CD =  7654321./16777216.
      CM = 16777213./16777216.
      IP = 97
      JP = 33
*
      FLAG = 'H1RN INITIALISED'

*
      RETURN
      END
