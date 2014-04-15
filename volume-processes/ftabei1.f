!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced
c              rather than shifting pls directly here.
 
 
      FUNCTION EIRENE_FTABEI1 (IREI,K)
c  evaluate electron impact rate coefficient (or rate), 
c  for ei process no. IREI,
c         in cell no. K

c  hard wired: cut off (density) parameter for fits: 1e8

 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMXS
 
      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: IREI, K
      REAL(DP) :: TBEIC(9), EIRENE_FTABEI1, DEIMIN, DSUB, PLS, TBEI,
     .            EIRENE_RATE_COEFF
      REAL(DP) :: ERATE
      INTEGER :: J, I, II, KK
 
      TBEI=0.D0
      KK = NREAEI(IREI)
 

      DEIMIN=LOG(1.D8)
      PLS=MAX(DEIMIN,DEINL(K))

c   density parameter rescaling: now done in rate_coeff(....,1,..)
c                                only for double polynomial fit
!pb   DSUB=LOG(1.D8)
!pb   PLS=MAX(DEIMIN,DEINL(K))-DSUB

 
      TBEI = EIRENE_RATE_COEFF(KK,TEINL(K),PLS,.TRUE.,1,ERATE)*
     .       FACREI(IREI,1)
      IF (IFTFLG(KK,2) < 100) TBEI=TBEI*DEIN(K)
 
      EIRENE_FTABEI1 = TBEI
 
      RETURN
      END
 
