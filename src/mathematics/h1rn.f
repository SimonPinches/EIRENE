C  random number generator RANMAR, F. James, CPC 60, (1990), 329-344
C  period length: 2**144
*
      FUNCTION H1RN()
*
*#**********************************************************************
*# RANDOM NUMBER GENERATOR AS ADVOCATED BY F. JAMES FROM PROPOSAL OF   *
*# MARSAGLIA AND ZAMAN FSU-SCRI-87-50 AND MODIFIED BY F. JAMES 1988 TO *
*# PRODUCE VECTOR OF NUMBERS.                                          *
*# ENTRIES ARE:                                                        *
*#     FUNCTION    H1RN()          SINGLE RANDOM NUMBER                *
*#     SUBROUTINE  H1RNV(VEC,LEN)  VECTOR OF RANDOM NUMBERS            *
*# cdr SUBROUTINE  H1RNIN(IJ,KL)   INITIALISE WITH 2 SEEDS             *
*#     SUBROUTINE  H1RNIN(IJKL)    INITIALISE WITH 1 SEED              *
*#     SUBROUTINE  H1RNIV(VEC)     INITIALISE/RESTART WITH SEED ARRAY  *
*#     SUBROUTINE  H1RNSV(VEC)     SAVE SEED ARRAY VEC(100)            *
*#                                                                     *
*# NOTE: -H1RNIN OR H1RNIV MUST BE CALLED BEFORE GENERATING ANY        *
*#        RANDOM NUMBER(S).                                            *
*#       -H1RNSV SAVES SEED ARRAY INTO VEC(100) ONLY. THE USER HAS TO  *
*#        OUTPUT IT.                                                   *
*#                                                                     *
*# CHANGED BY: G. GRINDHAMMER AT: 90/03/14                             *
*# REASON : ?                                                           *
*# CHANGED BACK TO ORIGINAL BY D REITER AT 17/04/05                    *
*#**********************************************************************
*
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      REAL(DP) :: H1RN
cdr   INTEGER :: ISEED1, ISEED2
      INTEGER :: ISEED
      CHARACTER(16) :: CHECK

cdr  next 5 lines: status of generator, already initialized with previous call 
cdr  either to h1rnin(ijkl) or to H1RNIV(VEC), with VEC(100)
      CHARACTER(16) :: FLAG
      REAL(DP) :: U, C, CD, CM
      INTEGER :: I, J
      COMMON /RASET1/ U(97),C,CD,CM,I,J
      COMMON /RASET2/ FLAG
*
      LOGICAL, SAVE :: FIRST=.TRUE.

      DATA CHECK /'H1RN INITIALISED'/
*
      IF (FIRST) THEN
         IF (FLAG .NE. CHECK) THEN

cdr  apparently H1RN is called without initialization.
cdr  Use default Marsaglia-Zaman seeds:
cdr         WRITE(IUNOUT,*) ' H1RN (RANMAR): INITIALIZED WITH DEFAULT SEED'
cdr  changed back to single default seed, also used in seed driver routine RANSET
cdr         ISEED1      = 1802  ! = ij
cdr         ISEED2      = 9373  ! = kl
cdr         CALL H1RNIN(ISEED1,ISEED2)

cdr this single default seed produces the 4 Marsaglia-Zaman seeds.
c   ijkl= ij*30082+kl
cdr Loc.cit. F.James, CPC (1990), p340
            ISEED = 54217137  ! = ijkl
            CALL H1RNIN(ISEED)
         ENDIF
         FIRST = .FALSE.
      ENDIF
*
  100 CONTINUE
      H1RN = U(I)-U(J)
      IF(H1RN .LT. 0.) H1RN = H1RN + 1.
      U(I) = H1RN
      I = I - 1
      IF( I .EQ. 0) I=97
      J = J - 1
      IF( J .EQ. 0) J=97
      C = C - CD
      IF( C .LT. 0) C = C + CM
      H1RN = H1RN-C
      IF(H1RN .LE. 0.) H1RN = H1RN + 1.
      IF(H1RN .GE. 1.) GOTO 100
      RETURN
      END
