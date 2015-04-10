c  0707: new, for PI processes, copied and adapted from feelei1.f
 
      FUNCTION EIRENE_FEELPI1 (IRPI,K)

cdr  find electron energy loss for PI process no. IRPI,
c    locally in cell K  
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CCONA
      USE EIRMOD_COMXS
      USE EIRMOD_COMPRT, ONLY: IUNOUT
 
      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: IRPI, K
      REAL(DP) :: ELPIC(9), EIRENE_FEELPI1, PLS, DEIMIN, EE,
     .            EIRENE_FTABPI3,
     .            ELPI, EIRENE_ENERGY_RATE_COEFF, DELE
      INTEGER :: J, I, KK, II
 
      EIRENE_FEELPI1=0.D0
      KK=NELRPI(IRPI)


      IF (KK < 0) THEN
c   electron energy losses per collision from the default PI processes

        IF (KK == -1) THEN
          EIRENE_FEELPI1 = EPLPI3(IRPI,1,1)
        ELSE
          WRITE (IUNOUT,*) 'ERROR IN FEELPI1: KK IS LT 0. '
          WRITE (IUNOUT,*)
     .       'BUT THERE SHOULD BE NO DEFAULT PI PROCESSES '
          CALL EIRENE_EXIT_OWN(1)
        END IF

c  non default models, data from external databases
      ELSE IF (KK > 0) THEN
        IF (JELRPI(IRPI) == 1) THEN
          ELPI = EIRENE_ENERGY_RATE_COEFF(KK,TEINL(K),0._DP,.TRUE.,0)
          EIRENE_FEELPI1=-ELPI*DEIN(K)*FACRPI(IRPI,1)/
     .                   (EIRENE_FTABPI3(IRPI,K)+EPS60)
        ELSEIF(JELRPI(IRPI) == 9) THEN   !  Te, ne dependence. 
          DEIMIN=LOG(1.D8)
          PLS=MAX(DEIMIN,DEINL(K))
          ELPI = EIRENE_ENERGY_RATE_COEFF(KK,TEINL(K),PLS,.FALSE.,1)
          EE=MAX(-100._DP,ELPI+FACRPI(IRPI,2)+DEINL(K))
          EIRENE_FEELPI1=-EXP(EE)/(EIRENE_FTABPI3(IRPI,K)+EPS60)
        ELSE
CDR: missing still:  EB,Te dependence
          WRITE (IUNOUT,* ) 'ERROR IN FEELPI1, INVALID JELRPI '
          CALL EIRENE_EXIT_OWN(1)
        END IF
        IF (DELPOT(KK).NE.0.D0) THEN
          DELE=DELPOT(KK)
          EIRENE_FEELPI1=EIRENE_FEELPI1+DELE
        END IF
      END IF
 
      RETURN
      END
