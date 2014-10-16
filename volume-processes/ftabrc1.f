!pb  22.11.06: flag (T/F) for shift of first parameter to rate_coeff introduced.
c              this transformation of parameter PLS is now done in rate-coeff.f
 
 
      FUNCTION EIRENE_FTABRC1 (IRRC,K)
c  evaluate volume-recombination rate (1/s), also spontaneous volumetric transition rate (1/s)  
c  include density factor, if rate_coeff is in cm*3/s  (controlled by fitting flag iftflg) 
c  for rc process no. IRRC,
c         in cell no. K


c  hard wired: cut off (density) parameter for H.4 fits: 1e8
c  hard wired: density parameter in fit reduced by DSUB=1e8, done in rate_coeff.f
 
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CCONA
      USE EIRMOD_COMXS
 
      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: IRRC, K
      REAL(DP) :: DEIMIN, DSUB, EIRENE_FTABRC1, ZX, TBRC, PLS,
     .            EIRENE_RATE_COEFF, ERATE
      INTEGER :: II, KK
 
      TBRC=0.D0
      KK = NREARC(IRRC)
 
      IF (KK == 0) THEN
        ZX=EIONH/MAX(1.E-5_DP,TEIN(K))
        TBRC=1.27E-13*ZX**1.5/(ZX+0.59)*DEIN(K)
 
      ELSE
 

        DEIMIN=LOG(1.D8)
        PLS=MAX(DEIMIN,DEINL(K))

!pb     PL=MAX(DEIMIN,DEINL(K))-DSUB
!pb     DSUB=LOG(1.D8)
 
        TBRC = EIRENE_RATE_COEFF(KK,TEINL(K),PLS,.TRUE.,1,ERATE)*
     .         FACRRC(IRRC,1)
        IF (IFTFLG(KK,2) < 100) TBRC=TBRC*DEIN(K)
 
      ENDIF
 
      EIRENE_FTABRC1 = TBRC
 
      RETURN
      END
 
