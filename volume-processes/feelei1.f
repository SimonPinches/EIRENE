!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced
!pb  30.11.06: DELPOT introduced

cdr  eelds1 defined twice in case of default models, here and in xsectm, xsecta, xsecti, xsecpt
cdr  done: eelds1 set in xsect... routines.
 
      FUNCTION EIRENE_FEELEI1 (IREI,K)
cdr  find electron energy loss for EI process no. IREI,
c    locally in cell K  
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CCONA
      USE EIRMOD_COMXS
 
      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: IREI, K
      REAL(DP) :: ELEIC(9), EIRENE_FEELEI1, PLS, DEIMIN, EE,
     .            EIRENE_FTABEI1, 
     .            ELEI, EIRENE_ENERGY_RATE_COEFF, DELE
      INTEGER :: J, I, KK, II
 
      EIRENE_FEELEI1=0.D0
      KK=NELREI(IREI)


      IF (KK < 0) THEN
c   electron energy losses per collision from the 10 default EI processes
        SELECT CASE (KK)
        CASE (-1)
            EIRENE_FEELEI1=-EIONHE
        CASE (-2)
            EIRENE_FEELEI1=EELDS1(IREI,1)
        CASE (-3)
            EIRENE_FEELEI1=-1.5*TEIN(K)
        CASE (-4)
            EIRENE_FEELEI1=-EIONH
        CASE (-5)
            EIRENE_FEELEI1=-10.5
        CASE (-6)
            EIRENE_FEELEI1=-25.0
        CASE (-7)
            EIRENE_FEELEI1=EELDS1(IREI,1)
        CASE (-8)
            EIRENE_FEELEI1=-10.5
        CASE (-9)
            EIRENE_FEELEI1=-15.5
        CASE (-10)
C  FOR THE FACTOR -0.88 SEE: EIRENE MANUAL, INPUT BLOCK 4, EXAMPLES
            EIRENE_FEELEI1=-0.88*TEIN(K)
        END SELECT

c  non default models, data from external databases
      ELSE IF (KK > 0) THEN
        IF (JELREI(IREI) == 1) THEN  !  Te dependence
          ELEI = EIRENE_ENERGY_RATE_COEFF(KK,TEINL(K),0._DP,.TRUE.,0)
          EIRENE_FEELEI1=-ELEI*DEIN(K)*FACREI(IREI,1)/
     .                   (EIRENE_FTABEI1(IREI,K)+EPS60)
        ELSE    !  Te, ne dependence. missing still:  EB,Te dependence
          DEIMIN=LOG(1.D8)
          PLS=MAX(DEIMIN,DEINL(K))

          ELEI = EIRENE_ENERGY_RATE_COEFF(KK,TEINL(K),PLS,.FALSE.,1)
          EE=MAX(-100._DP,ELEI+FACREI(IREI,2))
          EIRENE_FEELEI1=-EXP(EE)*DEIN(K)/(EIRENE_FTABEI1(IREI,K)+EPS60)
        END IF
        IF (DELPOT(KK).NE.0.D0) THEN
          DELE=DELPOT(KK)
          EIRENE_FEELEI1=EIRENE_FEELEI1+DELE
        END IF
      END IF
 
      RETURN
      END
