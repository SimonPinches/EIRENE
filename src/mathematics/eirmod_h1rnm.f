      MODULE EIRMOD_H1RNM
cym 04/2020
cym groups h1rn.f,h1rnin.f, h1rniv.f, hlrnsv.f, hlrnv.f
cym common /raset1/ /raset2/ removed
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      PRIVATE

      PUBLIC :: H1RN,H1RNIN,H1RNIV,H1RNSV,H1RNV
cym care needs to taken that some variables chaneg names in the commons
cmathematics/h1rn.f:40:      COMMON /RASET1/ U(97),C,CD,CM,I,J
cmathematics/h1rn.f:41:!$omp threadprivate(/raset1/)            
cmathematics/h1rnin.f:31:      COMMON /RASET1/ U(97),C,CD,CM,IP,JP
cmathematics/h1rnin.f:32:!$omp threadprivate(/raset1/)     
cmathematics/h1rnin.f:52:cdr  next: initialize RANMAR, data are transfered via common RASET1
cmathematics/h1rniv.f:15:      COMMON /RASET1/ U(97),C,CD,CM,IP,JP
cmathematics/h1rniv.f:16:!$omp threadprivate(/raset1/)      
cmathematics/h1rnsv.f:13:      COMMON /RASET1/ U(97),C,CD,CM,I,J
cmathematics/h1rnsv.f:14:!$omp threadprivate(/raset1/)      
cmathematics/h1rnv.f:17:      COMMON /RASET1/ U(97),C,CD,CM,I,J
cmathematics/h1rnv.f:18:!$omp threadprivate(/raset1/)

cmathematics/h1rn.f:42:      COMMON /RASET2/ FLAG
cmathematics/h1rn.f:43:!$omp threadprivate(/raset2/)      
cmathematics/h1rnin.f:33:      COMMON /RASET2/ FLAG
cmathematics/h1rnin.f:34:!$omp threadprivate(/raset2/)      
cmathematics/h1rniv.f:17:      COMMON /RASET2/ FLAG
cmathematics/h1rniv.f:18:!$omp threadprivate(/raset2/)
cmathematics/h1rnv.f:19:      COMMON /RASET2/ FLAG
cmathematics/h1rnv.f:20:!$omp threadprivate(/raset2/)


cym - former RASET1 common
      real(dp) :: U(97)
      real(dp) :: C,CD,CM
cym these where IP,JP in some places and I,J elsewhere 
      integer :: IM,JM

cym - former RASET2 common

      CHARACTER*16 :: FLAG

!$omp threadprivate(U,C,CD,CM,IM,JM,FLAG)

      LOGICAL, SAVE :: FIRST1=.TRUE.

      CONTAINS


C  random number generator RANMAR, F. James, CPC 60, (1990), 329-344
C  period length: 2**144
*
      FUNCTION H1RN(DUMMY)
*
*#**********************************************************************
*# RANDOM NUMBER GENERATOR AS ADVOCATED BY F. JAMES FROM PROPOSAL OF   *
*# MARSAGLIA AND ZAMAN FSU-SCRI-87-50 AND MODIFIED BY F. JAMES 1988 TO *
*# PRODUCE VECTOR OF NUMBERS.                                          *
*# ENTRIES ARE:                                                        *
*#     FUNCTION    H1RN(DUMMY)     SINGLE RANDOM NUMBER                *
*#     SUBROUTINE  H1RNV(VEC,LEN)  VECTOR OF RANDOM NUMBERS            *
*# cdr SUBROUTINE  H1RNIN(IJ,KL)   INITIALISE WITH SEEDS               *
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
      REAL(DP), INTENT(IN) :: DUMMY
      REAL(DP) :: H1RN
cdr   INTEGER :: ISEED1, ISEED2
      INTEGER :: ISEED
      CHARACTER(16) :: CHECK

cym cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cdr  next 5 lines: status of generator, already initialized with previous call
cdr  either to h1rnin(ijkl) or to H1RNIV(VEC), with VEC(100)
cym      CHARACTER(16) :: FLAG
cym      REAL(DP) :: U, C, CD, CM
cym      INTEGER :: I, J
cym      COMMON /RASET1/ U(97),C,CD,CM,I,J
cym!$omp threadprivate(/raset1/)            
cym      COMMON /RASET2/ FLAG
cym!$omp threadprivate(/raset2/)
cym cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
*
cym  moved to module declaration part- check - renamed FIRST1
cym      LOGICAL, SAVE :: FIRST1=.TRUE.

      DATA CHECK /'H1RN INITIALISED'/
*
      IF (FIRST1) THEN
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
         FIRST1 = .FALSE.
      ENDIF
*
  100 CONTINUE
cym cccccccccccccccccc
cym I,J->IJ,IM
cym cccccccccccccccccc

      H1RN = U(IM)-U(JM)
      IF(H1RN .LT. 0.) H1RN = H1RN + 1.
      U(IM) = H1RN
      IM = IM - 1
      IF( IM .EQ. 0) IM=97
      JM = JM - 1
      IF( JM .EQ. 0) JM=97
      C = C - CD
      IF( C .LT. 0) C = C + CM
      H1RN = H1RN-C
      IF(H1RN .LE. 0.) H1RN = H1RN + 1.
      IF(H1RN .GE. 1.) GOTO 100
      RETURN
      END FUNCTION H1RN

cdr  Initializes random number generator H1RN  ( = RANMAR, F. James, see below)
c
cdr  old code (before April 2017):
cdr   SUBROUTINE H1RNIN(IJ,KL)
cdr   IMPLICIT NONE
cdr   INTEGER, INTENT(IN,OUT) :: IJ, KL
cdr Version in eirene until April 2017: take 2 input seeds,
cdr and produce 4 smaller input seeds from them.
cdr NOTE: The seed variables can have values between: 0 <= IJ <= 31328
cdr                                                   0 <= KL <= 30081

cdr April 2017: back to original suggestions by F. James
cdr Take only one seed IJKL, then first make 2 smaller (IJ, KL) legal seeds,
cdr then 4 yet smaller  (I,J,K,L) from them.
cdr IJKL must be not larger then 900.000.000

cdr
c Default Marsaglia-Zaman-seeds:
c Use IJ = 1802 & KL = 9373 or, equivalently,
c IJKL=54217137, to test the random number generator. The
c subroutine RANMAR should be used to generate 20000 random numbers.
c Then display the next six random numbers generated multiplied by
c 4096*4096
c If the random number generator is working properly, the random numbers
c should be:
c           6533892.0  14220222.0   7275067.0
c           6172232.0   8354498.0  10633180.0
cdr


      SUBROUTINE H1RNIN(IJKL)
cdr F. James, "A review of pseudorandom number generators",
cdr           Comp. Physics Communications 60 (1990) 329 - 344
cdr Subroutine RMARIN, p340.
cdr Legal seed IJKL range 0<= IJKL<= 900.000.000 must be already enforced
cdr in calling routine. (Strictly: IJKL <= 942.408.896, see below.)
cdr Each seed then produces an
cdr independent (non-overlapping) random sequence of
cdr average length about 10**30.
c
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IJKL
      INTEGER :: L, II, J, K, JJ, M, I, IJ, KL
      REAL(DP) :: S, T
cdr  next 5 lines: status of generator, initialized with IJKL

cym ccccccccccccccccccccccccccccccccccccccccccccccc
cym      CHARACTER*16    FLAG
cym      REAL(DP) :: U, C, CD, CM
cym      INTEGER :: IP, JP
cym      COMMON /RASET1/ U(97),C,CD,CM,IP,JP
cym !$omp threadprivate(/raset1/)     
cym      COMMON /RASET2/ FLAG
cym !$omp threadprivate(/raset2/)
cym cccccccccccccccccccccccccccccccccccccccccccccccc      
*
cdr  old code
cdr  enforce legal values IJ, KL:
cdr   IJ = IABS(IJ)
cdr   KL = IABS(KL)
cdr   IJ = MOD(IJ,31329)
cdr   KL = MOD(KL,30082)

cdr                        change back to original generator from
cdr                        F. James, CPC 60 (1990) 329 - 344, page 340
cdr  the max legal 4 digit seed is:      IJKL_max=(30082*31329)-1 = 942......,
cdr  but to avoid round off errors take: IJKL_max= 30082*31328    = 942.408.896
      IJ = IJKL/30082
      KL = IJKL-30082*IJ    ! = MOD(IJKL,30082)
c
      I  = MOD(IJ/177, 177) + 2
      J  = MOD(IJ, 177)     + 2
      K  = MOD(KL/169, 178) + 1
      L  = MOD(KL, 169)
*
cdr  now we have the 4 small seeds.
cdr  next: initialize RANMAR, data are transferred via common RASET1

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
cym IP->IM
cym   IP = 97
      IM = 97
cym   JP = 33
      JM = 33
c
      FLAG = 'H1RN INITIALISED'

c
      RETURN
      END SUBROUTINE H1RNIN

c
      SUBROUTINE H1RNV(RVEC,LEN)
cdr Return a vector RVEC, of length LEN, of random numbers
cdr from H1RN (=RANMAR) generator, in a single call.
c
      USE EIRMOD_PRECISION
      USE EIRMOD_COMPRT, ONLY: IUNOUT
   
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LEN
      REAL(DP), INTENT(OUT) :: RVEC(LEN)
      REAL(DP) :: UNI
      INTEGER :: IVEC
  
      CHARACTER(16) :: CHECK

cym ccccccccccccccccccccccccccccccccccccccc
cym      CHARACTER(16) :: FLAG,CHECK
cym      REAL(DP) :: U, C, CD, CM
cym      INTEGER :: I, J
cym      COMMON /RASET1/ U(97),C,CD,CM,I,J
cym !$omp threadprivate(/raset1/)
cym      COMMON /RASET2/ FLAG
cym !$omp threadprivate(/raset2/)
cym ccccccccccccccccccccccccccccccccccccccc
*
      LOGICAL FIRST
      DATA FIRST /.TRUE./
      DATA CHECK /'H1RN INITIALISED'/
*
      IF (FIRST) THEN
         IF (FLAG .NE. CHECK) THEN
            WRITE(IUNOUT,*)
     .                 ' H1RNV (RANMAR): CALL H1RNIN OR H1RNIV BEFORE',
     >                 ' CALLING H1RN.'
            CALL EIRENE_EXIT_OWN(1)
         ELSE
            FIRST = .FALSE.
         ENDIF
      ENDIF
*
      DO 200 IVEC = 1,LEN
  190    CONTINUE

cym cccccccccccccccccc
cym I,J -> IM,JM
cym cccccccccccccccccc

         UNI = U(IM)-U(JM)
         IF(UNI .LT. 0.) UNI = UNI + 1.
         U(IM) = UNI
         IM = IM - 1
         IF( IM .EQ. 0) IM=97
         JM = JM - 1
         IF( JM .EQ. 0) JM=97
         C = C - CD
         IF( C .LT. 0) C = C + CM
         UNI = UNI-C
         IF(UNI .LE. 0.) UNI = UNI + 1.
         IF(UNI .GE. 1.) GOTO 190
         RVEC(IVEC) = UNI
  200 CONTINUE
      RETURN
      END SUBROUTINE H1RNV
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

cym ccccccccccccccccccccccccccccccccccccccc      
cym      CHARACTER*16    FLAG
cym      REAL(DP) :: U, C, CD, CM
cym      INTEGER :: IP, JP
cym      COMMON /RASET1/ U(97),C,CD,CM,IP,JP
cym      !$omp threadprivate(/raset1/)      
cym      COMMON /RASET2/ FLAG
cym      !$omp threadprivate(/raset2/)
cym cccccccccccccccccccccccccccccccccccccccc
*
      DO 400 IC = 1, 97
        U(IC) = VEC(IC)
  400 CONTINUE
      C  = VEC(98)
      CD =  7654321./16777216.
      CM = 16777213./16777216.

cym IP,JP -> IM,JM

      IM = NINT(VEC(99))
      JM = NINT(VEC(100))
*
      FLAG = 'H1RN INITIALISED'
c     WRITE(IUNOUT,*)
c    .           ' H1RNIV: H1RN (RANMAR) INITIALISED/RESTARTED WITH',
c    >           ' SEED ARRAY VEC(100)'
*
      RETURN
      END SUBROUTINE H1RNIV


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
      
cym cccccccccccccccccccccccccccccccccccccc
cym      REAL(DP) :: U, C, CD, CM
cym      INTEGER :: I, J
cym      COMMON /RASET1/ U(97),C,CD,CM,I,J
cym      !$omp threadprivate(/raset1/)
cym cccccccccccccccccccccccccccccccccccccc
      
*
      DO 10 IC = 1, 97
        VEC(IC) = U(IC)
   10 CONTINUE
      VEC(98) = C

cym I,J -> IM,JM
      VEC(99) = REAL(IM)
      VEC(100)= REAL(JM)
      RETURN
      END SUBROUTINE H1RNSV


      END MODULE EIRMOD_H1RNM
