!pb  22.11.06:   flag for shift of first parameter to rate_coeff introduced
c                rather than shifting pls directly here.
cdr  jan 2014:   comments.


      FUNCTION EIRENE_FTABEI1 (IREI,K)
C  this is the "on the fly" storage saving version to eliminate
C  pre-computed array TABEI1(irei,k) from this run
c
c  evaluate electron impact rate (1/s),  include density factor
c  for EI process no. IREI,
c         in cell no. K

c  call to fct. RATE_COEFF,

c  hard-wired: Cut-off (density) parameter for fits: 1e8 cm**-3.
c              Minimum density is set to 1e8 cm**-3
c              for evaluation of rate coefficient.
c              But then the 'true' density is used to return a 'rate'
c  May 18: sync with xstei.f, Te cut-off earlier than TVAC.
c              Cut-off (Te) parameter for fits: 0.1 eV.
c              Minimum temperature is set to 0.1 eV
c              for evaluation of rate coefficient.


      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMXS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IREI, K
      REAL(DP) :: EIRENE_FTABEI1, DEIMIN, PLS, TBEI,
     .            EIRENE_RATE_COEFF
      REAL(DP) :: TEE
      INTEGER :: KK

      TBEI=0.D0
      KK = NREAEI(IREI)   !  KK=-11:NREAC in REACDAT, IFTFLG 
!                            KK=  1:NREAC else


      DEIMIN=LOG(1.D8)
      PLS=MAX(DEIMIN,DEINL(K))

c   density parameter rescaling: now done in rate_coeff(....,1)
cdr careful, hidden link: this is only for double polynomial fit (and: MODC=3)

cdr  safety cut-off at Te= 0.1 eV. (note: TVAC=0.02)
      TEE = max(-2.3_dp,TEINL(K))

      TBEI = EIRENE_RATE_COEFF(KK,K,TEE,PLS,.TRUE.,1)*
     .       FACREI(IREI,1)
      IF (IFTFLG(KK,2) < 100) TBEI=TBEI*DEIN(K)

      EIRENE_FTABEI1 = TBEI

      RETURN
      END FUNCTION EIRENE_FTABEI1
