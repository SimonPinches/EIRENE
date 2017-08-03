c  0707: new, for PI processes, copied and adapted from feelei1.f






CDR TO BE DONE:  when kk >0  then on the fly evaluation of rate coeff. is
cdr              repeated here. This should be avoided, by returning the energy weighted rate,
cdr              rather than the mean electron energy itself.




      FUNCTION EIRENE_FEELPI1 (IRPI,K)
C  this is the "on the fly", storage saving, version to eliminate
C  pre-computed array EELPI1(irpi,k) from this run

cdr  find electron energy loss for PI process no. IRPI,  energy in eV
c    locally in cell K, for process kk= nelrpi(irpi) 
c
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CCONA
      USE EIRMOD_COMXS
      USE EIRMOD_COMPRT, ONLY: IUNOUT
 
      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: IRPI, K
      REAL(DP) :: EIRENE_FEELPI1, PLS, DEIMIN, EE,
     .            EIRENE_FTABPI3,
     .            ELPI, EIRENE_ENERGY_RATE_COEFF, DELE
      INTEGER :: KK
 
      EIRENE_FEELPI1=0.D0
      KK=NELRPI(IRPI)


      IF (KK < 0) THEN
c   electron energy losses per collision from the default PI processes

        WRITE (IUNOUT,*) 'ERROR IN FEELPI1: KK IS LT 0. '
        WRITE (IUNOUT,*)
     .       'BUT THERE SHOULD BE NO DEFAULT PI PROCESSES '
        CALL EIRENE_EXIT_OWN(1)

c  non default models, data from external databases
      ELSE IF (KK > 0) THEN
        IF (JELRPI(IRPI) == 1) THEN  !  Te dependence
          ELPI = EIRENE_ENERGY_RATE_COEFF(KK,K,TEINL(K),0._DP,.TRUE.,0)
          EIRENE_FEELPI1=-ELPI*DEIN(K)*FACRPI(IRPI,1)/
     .                   (EIRENE_FTABPI3(IRPI,K)+EPS60)
        ELSEIF(JELRPI(IRPI) == 9) THEN   !  Te, ne dependence. 
          DEIMIN=LOG(1.D8)
          PLS=MAX(DEIMIN,DEINL(K))
          ELPI = EIRENE_ENERGY_RATE_COEFF(KK,K,TEINL(K),PLS,.FALSE.,1)
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
