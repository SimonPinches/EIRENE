      Subroutine RANMAR_TEST(IADD, IERR, RAN)

cdr  Nov. 2019: testing the RANMAR (=H1RN) random number generator.
c
c    Call ranmar_test at any time and from any thread, to verify
c    A)  the proper functioning of initialization, grep-status,
c    B)  the (Marsaglia) reference random sequence.
c           6533892.0  14220222.0   7275067.0
c           6172232.0   8354498.0  10633180.0
cdr
c     input:  IADD  (>= 0) advance RANMAR by IADD steps (random numbers)
c                          IADD=0: the random number sequence is not altered.
c     output: RAN   latest random number generated in this call
c                   (in case IADD=0: RAN=run_e, next random number in sequence)
c             IERR  error code: 7 digits:
c                   if any of the 6 lowest digits =1: Marsaglia sequence failed
c                   if digit 7 =1:  initialisation and/or grep-status failed.

      USE EIRMOD_PRECISION
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_H1RNM
      IMPLICIT NONE


cdr  next 5 lines: status of generator, initialized with previous call to H1RNIN
      CHARACTER(16) :: FLAG
      REAL(DP) :: U, C, CD, CM
      INTEGER :: I, J
      COMMON /RASET1/ U(97),C,CD,CM,I,J
      COMMON /RASET2/ FLAG

      INTEGER, INTENT(IN) :: IADD
      INTEGER, INTENT(OUT) :: IERR
      REAL(DP), INTENT(OUT) :: RAN
      REAL(DP) :: VEC(100), RAN_I, RAN_E , RMARSAG(6), RVEC(20000)
      real(dp) :: dum  !pb
      INTEGER :: IJKL, IRAN
!pb      real(DP), external :: H1RN

      ierr=0
      rmarsag(1)= 6533892.0
      rmarsag(2)=14220222.0
      rmarsag(3)= 7275067.0
      rmarsag(4)= 6172232.0
      rmarsag(5)= 8354498.0
      rmarsag(6)=10633180.0

c  save the current status on VEC(100)
      CALL H1RNSV(VEC)

c  fetch new random number and save it for testing further below.
!pb   RAN_I= H1RN()
      RAN_I= H1RN(dum)

c  Set the original Marsaglia-seed (converted to single integer seed):
c  Old coding: James 2-seed option:
cdr         IJ      = 1802
cdr         KL      = 9373
cdr         CALL H1RNIN(IJ,KL)

c  This single seed should produce identically the same as IJ, KL from above.
      IJKL = 54217137
      call h1rnin(ijkl)



cdr  Skip 20000 random numbers
c     do IRAN=1,20000
c       ran=H1RN()
c     enddo
c  this single call should do the same as the above loop.
      call H1RNV(RVEC,20000)

cdr  next: the next 6 random number should produce the
cdr        Marsaglia sequence exactly.
      Do IRAN=1,6
!pb     ran=H1RN()
        ran=H1RN(dum)
        ran=ran*4096*4096
        if (RAN.EQ.RMARSAG(IRAN)) cycle
        write (IUNOUT,*) 'TEST RANMAR FAILED', IRAN, RAN, RMARSAG(IRAN)
        ierr=ierr+10**(IRAN-1)
      enddo

cdr back to original status, VEC(100) --> Common
      call H1RNIV(VEC)
cdr  check, if same status as before indeed ?
!pb   ran_e=H1RN()
      ran_e=H1RN(dum)
cdr check: ran_i=ran_e ??
      if (ran_i.ne.ran_e) then
        write (IUNOUT,*) 'INIT RANMAR FAILED ', RAN_I, RAN_E
        ierr=ierr+10**6
      endif

      if (IADD.eq.0) then
c  continue previous random number sequence without any change of sequence
        ran=ran_e
        call H1RNIV(VEC)
      elseif (IADD.eq.1) then
c  random sequence is one step further, due to one call to H1RN
        ran=ran_e
      elseif (IADD.gt.1) then
cdr one step is already done, now the further IADD-1 steps
        do IRAN=1,IADD-1
!pb       RAN=H1RN()
          RAN=H1RN(dum)
        enddo
      endif

      return
      end subroutine ranmar_test
