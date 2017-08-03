!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced
!pb  30.11.06: DELPOT introduced
cdr  21.09.15: default process rate coeff KK=-1 (ei on He) now moved to KK=-11, to avoid conflict
cdr            with default cross section KK=-1 (cx on H)
cdr            default process kk=-10 (diss rec of H2+) slightly changed,
cdr            to enable external database model which is truely identical to default model

CDR TO BE DONE:  when kk >0  then on the fly evaluation of rate coeff. is
cdr              repeated here. This should be avoided, by returning the energy weighted rate,
cdr              rather than the mean electron energy itself.

cdr  ARRAY eelei1 defined twice in case of default models, here and in xsectm, xsecta, xsecti, xsecpt
cdr  done: eelei1 set in xsect... routines.
!pb  APR   16: eelds -> eelei 

      FUNCTION EIRENE_FEELEI1 (IREI,K)
C  this is the "on the fly" storage saving version to eliminate
C  pre-computed array EELEI1(irei,k) from with run

cdr  find electron energy loss for EI process no. IREI,  energy in eV
c    locally in cell K, for process kk= nelrei(irei) 
c
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CCONA
      USE EIRMOD_COMXS
      USE EIRMOD_COMPRT, ONLY: IUNOUT
 
      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: IREI, K
      REAL(DP) :: EIRENE_FEELEI1, PLS, DEIMIN, EE,
     .            EIRENE_FTABEI1, 
     .            ELEI, EIRENE_ENERGY_RATE_COEFF, DELE, DE_10
      INTEGER :: KK
 
      EIRENE_FEELEI1=0.D0
      KK=NELREI(IREI)


      IF (KK < 0) THEN
c   electron energy losses per collision from the default EI processes
        SELECT CASE (KK)
        CASE (-1)
c           EIRENE_FEELEI1=-EIONHE         ! DEFAULT PROCESS KK=-1: NOT IN USE  
        CASE (-2)
c           EIRENE_FEELEI1=EELEI1(IREI,1)  ! DEFAULT PROCESS KK=-2: NOT IN USE
        CASE (-3)
c           EIRENE_FEELEI1=-1.5*TEIN(K)    ! DEFAULT PROCESS KK=-3: NOT IN USE
        CASE (-4)
            EIRENE_FEELEI1=-EIONH   !  DEFAULT PROCESS KK=-4  H+E --> H+ + 2E
        CASE (-5)    
            EIRENE_FEELEI1=-10.5  ! DEFAULT PROCESS KK=-5:  H2+E --> H+H +E,  
        CASE (-6)
            EIRENE_FEELEI1=-25.0  ! DEFAULT PROCESS KK=-6:  H2+E --> H + H+  +2E 
        CASE (-7)
            EIRENE_FEELEI1=EELEI1(IREI,1) ! DEFAULT PROCESS KK=-7: H2+E --> H2+  +2E
        CASE (-8)
            EIRENE_FEELEI1=-10.5
        CASE (-9)
            EIRENE_FEELEI1=-15.5
        CASE (-10)  ! DEFAULT PROCESS KK=-10: H2+ E --> H + H, DISS. RECOMBINATION
C  FOR THE FACTOR -0.896... SEE: EIRENE MANUAL, INPUT BLOCK 4, EXAMPLES
            DE_10=8.964355004318D-01
            EIRENE_FEELEI1=-DE_10*TEIN(K)
        CASE (-11)
            EIRENE_FEELEI1=-EIONHE   !  FORMERLY DEFAULT PROCESS KK=-1  HE+E --> HE+ +2E
        END SELECT

c  non default models, data from external databases
      ELSE IF (KK > 0) THEN
        IF (JELREI(IREI) == 1) THEN  !  Te dependence
          ELEI = EIRENE_ENERGY_RATE_COEFF(KK,K,TEINL(K),0._DP,.TRUE.,0)
          EIRENE_FEELEI1=-ELEI*DEIN(K)*FACREI(IREI,1)/
     .                   (EIRENE_FTABEI1(IREI,K)+EPS60)
        ELSEIF(JELREI(IREI) == 9) THEN    !  Te, ne dependence. 
          DEIMIN=LOG(1.D8)
          PLS=MAX(DEIMIN,DEINL(K))
          ELEI = EIRENE_ENERGY_RATE_COEFF(KK,K,TEINL(K),PLS,.FALSE.,1)
          EE=MAX(-100._DP,ELEI+FACREI(IREI,2))
          EIRENE_FEELEI1=-EXP(EE)*DEIN(K)/(EIRENE_FTABEI1(IREI,K)+EPS60)
        ELSE
CDR: missing still:  EB,Te dependence
          WRITE (IUNOUT,* ) 'ERROR IN FEELEI1, INVALID JELREI '
          CALL EIRENE_EXIT_OWN(1)
        END IF
        IF (DELPOT(KK).NE.0.D0) THEN
          DELE=DELPOT(KK)
          EIRENE_FEELEI1=EIRENE_FEELEI1+DELE
        END IF
      END IF
 
      RETURN
      END
