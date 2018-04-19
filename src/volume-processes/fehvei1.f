!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced
cdr  27.09.15:  renamed from fehvds1 to fehvei1, 
cdr  KK<0 cases still to be syncronized with ftabei1, feelei1, etc. 
cdr  nelrei, nreahv, 
cdr  01.01.16:  documented
cdr             unified numbering of reactions: kk= nreaei(irei), rather than nreahv(irei) (not ready)
 
      FUNCTION EIRENE_FEHVEI1 (IREI,K)
C  this is the "on the fly" storage saving version to eliminate
C  pre-computed array EHVEI1(irei,k) from with run

cdr  find heavy secondary particle energy for EI process no. IREI,  energy in eV
c    locally in cell K, for process kk= nreahv(irei)
c    sum over all heavy secondaries.
c    distribution to individual heavy secondary type and species is done later,
c    e.g. in veloei for sampling, and in update, collide for scoring 
c
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CCONA
      USE EIRMOD_COMXS
 
      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: IREI, K
      REAL(DP) :: EIRENE_FEHVEI1, EHVEI, EIRENE_FTABEI1,
     .            EIRENE_RATE_COEFF, DE_10
      INTEGER :: KK
 
      EIRENE_FEHVEI1=0.D0
      KK=NREAHV(IREI)
      IF (KK < 0) THEN
        SELECT CASE (KK)
c       CASE (-4)   ! DEFAULT PROCESS KK=-4:  H+ E --> H+ +E, no net energy transfer to H+
c          EIRENE_FEHVEI1 =0.0   
        CASE (-1)
            EIRENE_FEHVEI1=EHVEI1(IREI,1)
        CASE (-2)  ! DEFAULT PROCESS KK=-5:  H2+E --> H+H +E,  
            EIRENE_FEHVEI1=6.  ! DEFAULT PROCESS KK=-5:  H2+E --> H+H +E,    
        CASE (-3)
            EIRENE_FEHVEI1=10.0
        CASE (-4)
            EIRENE_FEHVEI1=8.6
        CASE (-5)
            EIRENE_FEHVEI1=0.5 
        CASE (-6)  ! DEFAULT PROCESS KK=-10: H2+ E --> H + H, DISS. RECOMBINATION 
C  FOR THE FACTOR -0.896... SEE: EIRENE MANUAL, INPUT BLOCK 4, EXAMPLES
            DE_10=8.964355004318D-01
            EIRENE_FEHVEI1=DE_10*TEIN(K)
c       CASE (-11)   ! DEFAULT PROCESS KK=-11:  He+ E --> He+ +E, no net energy transfer to He+!
c           EIRENE_FEHVEI1 =0.0 
        END SELECT

c  non default models, data from external databases
      ELSE IF (KK > 0) THEN
        EHVEI = EIRENE_RATE_COEFF(KK,K,TEINL(K),0._DP,.FALSE.,0)
        EHVEI=EXP(MAX(-100._DP,EHVEI+FACREI(IREI,2)))
        EIRENE_FEHVEI1=EHVEI*DEIN(K)/(EIRENE_FTABEI1(IREI,K)+EPS60)
      END IF
 
      RETURN
      END
