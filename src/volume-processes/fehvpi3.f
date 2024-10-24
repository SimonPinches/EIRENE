

      FUNCTION EIRENE_FEHVPI3 (IRPI,K)
C  this is the "on the fly" storage saving version to eliminate
C  pre-computed array EHVPI3(irpi,k,..) from with run

cdr  find heavy secondary particle energy for PI process no. IRPI,  energy in eV
c    locally in cell K, for process kk= nhvrpi(irpi)
c    sum over all heavy secondaries.
c    Distribution to individual heavy secondary type and species is done later,
c    e.g. in velopi for sampling, and in update, collide,.. for scoring


      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CCONA
      USE EIRMOD_COMXS
      USE EIRMOD_COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IRPI, K
      REAL(DP) :: EIRENE_FEHVPI3, EHVPI, EIRENE_FTABPI3,
     .            EIRENE_ENERGY_RATE_COEFF
      INTEGER :: KK
      EXTERNAL :: EIRENE_EXIT_OWN, EIRENE_FTABPI3,
     .            EIRENE_ENERGY_RATE_COEFF

      EIRENE_FEHVPI3=0.D0
      KK=NHVRPI(IRPI)

      IF (KK == 0) THEN
! TAKE A CONSTANT KER "EHEAVY" FROM INPUT FILE
         EIRENE_FEHVPI3=EHVPI3(IRPI,1,1)

      ELSE IF (KK > 0) THEN
        EHVPI = EIRENE_ENERGY_RATE_COEFF(KK,K,TEINL(K),0._DP,.FALSE.,0)
        EHVPI=EXP(MAX(-100._DP,EHVPI+FACRPI(IRPI,2)))
        EIRENE_FEHVPI3=EHVPI*DEIN(K)/(EIRENE_FTABPI3(IRPI,K)+EPS60)

      ELSE
        GOTO 999
      END IF

      RETURN

  999 CONTINUE
      WRITE (IUNOUT,*) 'FEHVPI3: INVALID PARAMETER NHVRPI '
      WRITE (IUNOUT,*) 'IRPI, NHVRPI ',IRPI,NHVRPI
      CALL EIRENE_EXIT_OWN(1)
      END
