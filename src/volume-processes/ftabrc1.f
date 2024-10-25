!pb  22.11.06: flag (T/F) for shift of first parameter to rate_coeff introduced.
c              this transformation of parameter PLS is now done in rate_coeff.f

cdr check dsub, cut-off etc....


      FUNCTION EIRENE_FTABRC1 (IRRC,K)
c  evaluate volume-recombination rate (1/s), also spontaneous volumetric transition rate (1/s)
c  for rc process no. IRRC,  (=== process KK=NREARC(IRRC))
c         in cell no. K
c  Include density factor, if rate_coeff is in cm*3/s  (controlled by fitting flag iftflg)


c  hard-wired: cut-off (density) parameter for H.4 fits: 1e8
c  hard-wired: density parameter in fit scaled by DSUB=1e-8, done in rate_coeff.f


      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CCONA
      USE EIRMOD_COMXS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IRRC, K
      REAL(DP) :: DEIMIN, ZX, TBRC, PLS,
     .            EIRENE_FTABRC1,
     .            EIRENE_RATE_COEFF
C      REAL(DP) :: DSUB
      INTEGER :: KK
      EXTERNAL :: EIRENE_EXIT_OWN, EIRENE_RATE_COEFF

      TBRC=0.D0
      KK = NREARC(IRRC)

c  default radiative rate coefficient, see xstrc.f
      IF (KK == -1) THEN   ! E + H+ --> H + hv
        ZX=EIONH/MAX(1.E-5_DP,TEIN(K))
        TBRC=1.27E-13*ZX**1.5/(ZX+0.59)*DEIN(K)

      ELSEIF (KK == -2) THEN  ! E + He+ --> He + hv
         ZX=EIONHE/MAX(1.E-5_DP,TEIN(K))
C  rate = [rate coeff <sig v>] times [electr. density], 1/s per ion
c    1.96e-14*sqrt(EionHe/Ry) = 3.5487E-14
         TBRC=3.5487E-14*ZX**1.5/(ZX+0.35)*DEIN(K)

      ELSEIF (KK .LE. 0 ) THEN
        WRITE (IUNOUT,*) 'INVALID KK IN FTABRC1, KK= ',KK
        CALL EIRENE_EXIT_OWN(1)

c  kk >  0
      ELSE

cdr  automatic cut-off at density 1E8: collapse fit to corona value.
        DEIMIN=LOG(1.D8)
        PLS=MAX(DEIMIN,DEINL(K))

cdr  DSUB=1e8, recaling of density parameter, done in rate coeff.
cdr  only valid for AMJUEL fits.

        TBRC = EIRENE_RATE_COEFF(KK,K,TEINL(K),PLS,.TRUE.,1)*
     .         FACRRC(IRRC,1)
        IF (IFTFLG(KK,2) < 100) TBRC=TBRC*DEIN(K)

      ENDIF

      EIRENE_FTABRC1 = TBRC

      RETURN
      END FUNCTION EIRENE_FTABRC1
