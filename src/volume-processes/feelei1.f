!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced
!pb  30.11.06: DELPOT introduced
cdr  21.09.15: default process rate coeff KK=-1 (ei on He) now moved to KK=-11, to avoid conflict
cdr            with default cross-section KK=-1 (cx on H)
cdr            default process kk=-10 (diss rec of H2+) slightly changed,
cdr            to enable external database model which is truly identical to default model

CDR TO BE DONE: when kk >0 then on the fly evaluation of rate coeff. is
cdr             repeated here. This should be avoided, by returning the energy-weighted rate,
cdr             rather than the mean electron energy itself.

cdr  ARRAY eelei1 defined twice in case of default models, here and in xsectm, xsecta, xsecti, xsecpt
cdr  done: eelei1 set in xsect... routines.
!pb  APR   16: eelds -> eelei
cdr  sept. 18: KK=0 options (e.g. constant electron energy cost) added

      FUNCTION EIRENE_FEELEI1 (IREI,K)
C  this is the "on the fly", storage saving, version to eliminate
C  pre-computed array EELEI1(irei,k) from this run

cdr  find electron energy loss for EI process no. IREI, energy in eV
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
C  IDENTIFY NUMBER OF PROCESS.
C  CURRENTLY KK=-11 -- KK=-4 EIRENE MINIMAL PROCESSES (FLAG JELREI IS UNUSED).
C            KK> 0  KREAD: COLLISION PROCESSES STORED ON REACDAT FROM EXTERNAL DATABASES
c                          USE JELREI FLAG FOR ELECTRON KINETIC ENERGY LOSS
C            KK= 0  SET DEFAULT ELECTRON KINETIC ENERGY LOSS FROM JELREI FLAG

      KK=NELREI(IREI)

      IF (KK < 0) THEN
c   electron energy losses per collision from the default EI processes -4 ....-11
        SELECT CASE (KK)
        CASE (-1)
          GOTO 999  !  DEFAULT PROCESS KK=-1:
                    !  NOT IN USE FOR EI PROCESSES
        CASE (-2)
          GOTO 999  !  DEFAULT PROCESS KK=-2:
                    !  NOT IN USE FOR EI PROCESSES
        CASE (-3)
          GOTO 999  !  DEFAULT PROCESS KK=-3:
                    !  NOT IN USE FOR EI PROCESSES
        CASE (-4)
            EIRENE_FEELEI1=-EIONH  ! DEFAULT PROCESS KK=-4:
                                   ! H + E --> H+ + 2E
        CASE (-5)
            EIRENE_FEELEI1=-10.5   ! DEFAULT PROCESS KK=-5:
                                   ! H2 + E --> H + H + E
        CASE (-6)
cdr changed Sept. 22:  25.0 --> 28.1, see xsectm default model kk=-6
            EIRENE_FEELEI1=-28.1   ! DEFAULT PROCESS KK=-6:
                                   ! H2 + E --> H + H+ + 2E
        CASE (-7)
            EIRENE_FEELEI1=EELEI1(IREI,1) ! DEFAULT PROCESS KK=-7:
                                   ! H2 + E --> H2+ + 2E
        CASE (-8)
            EIRENE_FEELEI1=-10.5   ! DEFAULT PROCESS KK=-8:
                                   ! H2+ + E --> H + H+ + E (DE)
        CASE (-9)
            EIRENE_FEELEI1=-15.5   ! DEFAULT PROCESS KK=-9:
                                   ! H2+ + E --> H+ + H+ + 2E (DI)
        CASE (-10)  ! DEFAULT PROCESS KK=-10:
                    ! H2+ E --> H + H, DISS. RECOMBINATION (DR)
C  FOR THE FACTOR -0.896... SEE: EIRENE MANUAL, INPUT BLOCK 4, EXAMPLES
            DE_10=8.964355004318D-01
            EIRENE_FEELEI1=-DE_10*TEIN(K)
        CASE (-11)
            EIRENE_FEELEI1=-EIONHE ! DEFAULT PROCESS KK=-11:
                                   ! HE + E --> HE+ + 2E
        END SELECT

c  non-default models, data from external databases, KK=KREAD for el. energy-weighted rates
      ELSE IF (KK > 0) THEN
        IF (JELREI(IREI) == 1) THEN  !  Te dependence only
          ELEI = EIRENE_ENERGY_RATE_COEFF(KK,K,TEINL(K),0._DP,.TRUE.,0)
          EIRENE_FEELEI1=-ELEI*DEIN(K)*FACREI(IREI,1)/
     .                   (EIRENE_FTABEI1(IREI,K)+EPS60)
        ELSEIF(JELREI(IREI) == 9) THEN    !  Te, ne dependence.
          DEIMIN=LOG(1.D8)
          PLS=MAX(DEIMIN,DEINL(K))
          ELEI = EIRENE_ENERGY_RATE_COEFF(KK,K,TEINL(K),PLS,.FALSE.,1)
          EE=MAX(-100._DP,ELEI+FACREI(IREI,2)+DEINL(K))
          EIRENE_FEELEI1=-EXP(EE)/(EIRENE_FTABEI1(IREI,K)+EPS60)
        ELSE
CDR: missing still: EB,Te dependence
          WRITE (IUNOUT,* ) 'ERROR IN FEELEI1, INVALID JELREI '
          CALL EIRENE_EXIT_OWN(1)
        END IF

        IF (DELPOT(KK).NE.0.D0) THEN
          DELE=DELPOT(KK)
          EIRENE_FEELEI1=EIRENE_FEELEI1+DELE
        END IF

      ELSE IF (KK == 0) THEN
CDR  NO EXTERNAL DATASET FOR ELECTRON COOLING RATE FOR PROCESS IREI
CDR  SET SIMPLE DEFAULTS
        IF (JELREI(IREI) == -1) THEN  !  constant electron energy loss:
          EE=EELEI1(IREI,1)
          EIRENE_FEELEI1=EE
        ELSEIF (JELREI(IREI) == -2) THEN  !  1.5 * Te
          EE=-1.5*TEIN(K)
          EIRENE_FEELEI1=EE
        ELSE
          GOTO 999
        ENDIF
      END IF

      RETURN

  999 CONTINUE
      WRITE (IUNOUT,*) 'FEELEI1: INVALID PARAMETER NELREI '
      WRITE (IUNOUT,*) 'IREI, NELREI ',IREI,NELREI
      CALL EIRENE_EXIT_OWN(1)
      END FUNCTION EIRENE_FEELEI1
