c  0717: new, for PI processes.
c        currently unfinished, and unused. 

CDR TO BE DONE:  when kk >0, 
cdr              then the "on the fly evaluation" of rate coeff. 
cdr              EIRENE_FTABPI3     
cdr              is repeated here. This should be avoided,
cdr              by returning the energy-weighted rate,
cdr              rather than the mean electron energy itself,
cdr              because in calling program ftabpi3 is likely to
cdr              be available already.




      FUNCTION EIRENE_FEELPI3 (IRPI,K)
C  this is the "on the fly", storage saving, version to eliminate
C  pre-computed array EELPI3(irpi,k,1,...,9) from this run

cdr  find electron energy loss for PI process no. IRPI, energy in eV
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
      REAL(DP) :: EIRENE_FEELPI3, PLS, DEIMIN, EE,
     .            EIRENE_FTABPI3,
     .            ELPI, EIRENE_ENERGY_RATE_COEFF, DELE
      INTEGER :: KK

      EIRENE_FEELPI3=0.D0
      KK=NELRPI(IRPI)


      IF (KK < 0) THEN
c   electron energy losses per collision from the default PI processes

        WRITE (IUNOUT,*) 'ERROR IN FEELPI3: KK IS LT 0. '
        WRITE (IUNOUT,*)
     .       'BUT THERE SHOULD BE NO DEFAULT PI PROCESSES '
        CALL EIRENE_EXIT_OWN(1)

c  non-default models, data from external databases
      ELSE IF (KK > 0) THEN
        IF (JELRPI(IRPI) == 1) THEN  !  Te dependence ? Should be: Ti
          ELPI = EIRENE_ENERGY_RATE_COEFF(KK,K,TEINL(K),0._DP,.TRUE.,0)
          EIRENE_FEELPI3=-ELPI*DEIN(K)*FACRPI(IRPI,1)/
     .                   (EIRENE_FTABPI3(IRPI,K)+EPS60)
        ELSEIF(JELRPI(IRPI) == 9) THEN   !  Te, ne dependence. ?
                                         !  Should be Ti,ni.
          DEIMIN=LOG(1.D8)
          PLS=MAX(DEIMIN,DEINL(K))
          ELPI = EIRENE_ENERGY_RATE_COEFF(KK,K,TEINL(K),PLS,.FALSE.,1)
          EE=MAX(-100._DP,ELPI+FACRPI(IRPI,2)+DEINL(K))
          EIRENE_FEELPI3=-EXP(EE)/(EIRENE_FTABPI3(IRPI,K)+EPS60)
        ELSE
CDR: missing still: EB,Ti dependence
          WRITE (IUNOUT,* ) 'ERROR IN FEELPI3, INVALID JELRPI '
          CALL EIRENE_EXIT_OWN(1)
        END IF
        IF (DELPOT(KK).NE.0.D0) THEN
          DELE=DELPOT(KK)
          EIRENE_FEELPI3=EIRENE_FEELPI3+DELE
        END IF
      END IF

      RETURN
      END
