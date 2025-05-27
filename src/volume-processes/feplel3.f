!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced


      FUNCTION EIRENE_FEPLEL3 (IREL,K)

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CCONA
      USE EIRMOD_COMPRT
      USE EIRMOD_COMXS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IREL, K
cdr  tbd: change name of variable EPEL. That would be the logical name for
cdr       tally EPEL: volume sampled electron (EL) energy loss (E) from field (P),
cdr       e.g. in recombination (RC) reactions. So far that EPEL tally is missing.
      REAL(DP) :: PLS, ADD, EPEL, EIRENE_FEPLEL3,
     .            EIRENE_RATE_COEFF
      INTEGER :: KK, IPLSTI
      EXTERNAL :: EIRENE_RATE_COEFF

      EIRENE_FEPLEL3=0.D0
      KK=NELREL(IREL)
      IPLSTI=MPLSTI(IPLS)
      IF (KK < 0) THEN
        SELECT CASE (KK)
        CASE (-1)
C  DEFAULT EL MODEL
C  OUT
        CASE (-2)
C  MEAN ENERGY FROM DRIFTING MONOENERGETIC
            EIRENE_FEPLEL3=EPLEL3(IREL,1,1)
            IF (LEDRIFT) EIRENE_FEPLEL3=EIRENE_FEPLEL3+EDRIFT(IPLS,K)
        CASE (-3)
C  MEAN ENERGY FROM DRIFTING MAXWELLIAN
            EIRENE_FEPLEL3=1.5*TIIN(IPLSTI,K)
            IF (LEDRIFT) EIRENE_FEPLEL3=EIRENE_FEPLEL3+EDRIFT(IPLS,K)
        END SELECT
      ELSE
C  MEAN ENERGY FROM SINGLE PARAMETER FIT KK
        PLS=TIINL(IPLSTI,K)+ADDEL(IREL,IPLS)
        EPEL = EIRENE_RATE_COEFF(KK,K,PLS,0._DP,.FALSE.,0)
        ADD=EPLEL3(IREL,1,1)
        EIRENE_FEPLEL3=EPEL*DIIN(IPLS,K)*ADD
      END IF

      RETURN
      END
