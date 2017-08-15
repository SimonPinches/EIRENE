!pb  22.11.06:   flag for shift of first parameter to rate_coeff introduced
c                rather than shifting pls directly here.
cdr  jan 2014:   comments.
cdr: to be done: remove erate from here (needed only for H-colrad option, move to better place)
 
 
      FUNCTION EIRENE_FTABEI1 (IREI,K)
C  this is the "on the fly" storage saving version to eliminate
C  pre-computed array TABEI1(irei,k) from this run  
c
c  evaluate electron impact rate (1/s),  include density factor  
c  for ei process no. IREI,
c         in cell no. K

c  call to fct. RATE_COEFF, 

c  hard wired: cut off (density) parameter for fits: 1e8. Minimum density is set to 1e8 cm**-3
c              for evaluation of rate coefficient. But then the 'true' density is used to return a 'rate'

 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMXS
 
      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: IREI, K
      REAL(DP) :: EIRENE_FTABEI1, DEIMIN, PLS, TBEI,
     .            EIRENE_RATE_COEFF
C      REAL(DP) :: DSUB
      INTEGER :: KK
 
      TBEI=0.D0
      KK = NREAEI(IREI)
 

      DEIMIN=LOG(1.D8)
      PLS=MAX(DEIMIN,DEINL(K))

c   density parameter rescaling: now done in rate_coeff(....,1,..)
c                                only for double polynomial fit
!pb   DSUB=LOG(1.D8)
!pb   PLS=MAX(DEIMIN,DEINL(K))-DSUB

 
      TBEI = EIRENE_RATE_COEFF(KK,K,TEINL(K),PLS,.TRUE.,1)*
     .       FACREI(IREI,1)
      IF (IFTFLG(KK,2) < 100) TBEI=TBEI*DEIN(K)
 
      EIRENE_FTABEI1 = TBEI
 
      RETURN
      END
 
