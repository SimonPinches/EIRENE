C 0406: default resonant CX for He in He+/He++ plasma added:
C       Janev (HYDHEL), 1987, reactions 5.3.1 and 6.3.1
C
C 0710: provide value of cross-section for reaction K
C       K=0 means no cross-section available for this reaction
c       k=-1,-2,-3: default (hard-wired) CX cross-sections
c 0315: increase kk>=-10 to kk>=-11 for He EI process, to fully reserve k=-1
c       for H+p CX as default process
C Nov.19: Arrhenius factor (not needed here,
cdr       but for coding consistency)
cdr Jan 20: Added proper cut-off for iftflg=3 at threshold XI.

      FUNCTION EIRENE_CROSS(AL,K,IR,FACT,TEXT)
C
C  CROSS-SECTION
C    AL=LN(ELAB), ELAB IN (EV)
C    RETURN CROSS-SECTION IN CM**2
C
C  K>0 : DATA FROM ARRAY REACDAT, I.E. FROM EXTERNAL DATABASE
C
C  K<0 : DEFAULT MODEL DEFINED IN SETUP_MINIMAL_REACTIONS, BUT NOW ALSO ON REACDAT
C
C  K=-1: H + H+ --> H+ + H   CX CROSS-SECTION, JANEV, 3.1.8
C        LINEAR EXTRAPOLATION ON LOG-LOG SCALE AT LOW ENERGY END FOR LN(SIGMA)
C        IDENTICAL TO hydhel.tex, H.1, 3.1.8
C
C  K=-2: He + He+ --> He+ + He  CX CROSS-SECTION, JANEV, 5.3.1
C        LINEAR EXTRAPOLATION AT LOW ENERGY END FOR LN(SIGMA)
C        IDENTICAL TO hydhel.tex, H.1, 5.3.1
C
C  K=-3: He + He++ --> He++ + He  CX CROSS-SECTION, JANEV, 6.3.1
C        LINEAR EXTRAPOLATION AT LOW ENERGY END FOR LN(SIGMA)
C        IDENTICAL TO hydhel.tex, H.1, 6.3.1

C  FOR -11 <= K <= -4  VARIOUS DEFAULT EI PROCESSES FOR H, H2, H2+, HE
C                      CURRENTLY ONLY RATE COEFF, SO NO DATA HERE.
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMXS
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CTRCEI, ONLY: TRCAMD

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: AL, FACT
      INTEGER, INTENT(IN) :: K, IR
      CHARACTER(LEN=*), INTENT(IN) :: TEXT
      REAL(DP) :: B(8), FP(6)
      REAL(DP) :: EIRENE_CROSS,
     .            ALMIN, ALMAX, COUMIN, COUMAX, EMIN, EMAX,
     .            RES, EIRENE_EXTRAP, E, XI,
     .            EIRENE_SNGL_POLY, EARRH0
      INTEGER :: I
      LOGICAL :: LEXP
      type(poly_data), pointer :: rpp
      type(fit_forms), pointer :: rpc
      EXTERNAL :: EIRENE_EXIT_OWN, EIRENE_EXTRAP, EIRENE_SNGL_POLY
C
      IF ((K >= -11) .AND. (K <= NREAC)) THEN

        IF (K == 0) THEN

          EIRENE_CROSS = 0._DP
          WRITE (IUNOUT,*) 'ERROR IN CROSS: K=0'
          WRITE (iunout,*) 'CALLED FROM ',TEXT
          WRITE (iunout,*) 'REACTION NO. ',IR
          WRITE (IUNOUT,*) 'NO CROSS-SECTION DATA AVAILABLE FOR ',
     .                     'THIS REACTION'
          CALL EIRENE_EXIT_OWN(1)

CDR:  SAME EVALUATION OF CROSS-SECTION FORMULAS FOR DEFAULT
CDR:  (K<0) AND NON-DEFAULT (K>0) CROSS-SECTIONS

        ELSE IF (IFTFLG(K,1) == 0) THEN

C  EVALUATE POLYNOMIAL FIT NO. K, IFTFLG=0
C  FILL CROSS-SECTION DATA, SINGLE PARAMETER POLYNOMIAL IN AL=LN(E)
          RPC => REACDAT(K)%CRS
          RPP => RPC%POLY
          FP(1:3) = RPC%FP1L
          FP(4:6) = RPC%FP1R
CDR no Arrhenius factor in case of cross-sections:
          EARRH0=0.0
          LEXP=.TRUE.
cdr  Notation: rpp%dblpol is a single parameter polynomial?
          RES = EIRENE_SNGL_POLY(RPP%DBLPOL,AL,
     .                            RPC%RC1MIN,RPC%RC1MAX,FP,
     .                            RPC%JFEX1MN,RPC%JFEX1MX,
     .                            EARRH0,TRCAMD,LEXP)

          EIRENE_CROSS = RES*FACT

        ELSE IF (IFTFLG(K,1) == 3) THEN
cdr Near threshold and high energy Born-Bethe asymptotically correct cross-section fit

cdr  XI: threshold. ie. E>XI necessarily for this fit.
          E = EXP(AL)
          XI= REACDAT(K)%CRS%POLY%DBLPOL(1,1)

cdr  careful: E and XI must relate to same mass.
          if (E .le. XI) then
            eirene_cross=0._dp
            return
          endif

C  default extrapolation ifexmn=-1 not yet available

          ALMIN=REACDAT(K)%CRS%RC1MIN
          ALMAX=REACDAT(K)%CRS%RC1MAX

C  ELAB BELOW MINIMUM ENERGY FOR FIT:
          IF (AL.LT.ALMIN) THEN
C  USE ASYMPTOTIC EXPRESSION NO. IFEXMN(K)
            EMIN = EXP(ALMIN)
            B(1:8) = REACDAT(K)%CRS%POLY%DBLPOL(2:9,1)
cdr EMIN .gt. XI is already fulfilled here.
            COUMIN = B(1)*LOG(EMIN/XI)
            DO I=1,7
              COUMIN = COUMIN + B(I+1)*(1.D0-XI/EMIN)**I
            END DO
            COUMIN = COUMIN * 1.D-13/(XI*EMIN)

            FP(1:3) = REACDAT(K)%CRS%FP1L
            RES=EIRENE_EXTRAP(AL,ALMIN,COUMIN,
     .                   REACDAT(K)%CRS%JFEX1MN,
     .                   FP(1),FP(2),FP(3))

            EIRENE_CROSS = RES*FACT

C  ELAB ABOVE MAXIMUM ENERGY FOR FIT:
          ELSEIF (AL.GT.ALMAX) THEN
C  USE ASYMPTOTIC EXPRESSION NO. IFEXMX(K,1)
            EMAX = EXP(ALMAX)
            B(1:8) = REACDAT(K)%CRS%POLY%DBLPOL(2:9,1)
            COUMAX = B(1)*LOG(EMAX/XI)
            DO I=1,7
              COUMAX = COUMAX + B(I+1)*(1.D0-XI/EMAX)**I
            END DO
            COUMAX = COUMAX * 1.D-13/(XI*EMAX)

            FP(1:3) = REACDAT(K)%CRS%FP1R
            RES=EIRENE_EXTRAP(AL,ALMAX,COUMAX,
     .                   REACDAT(K)%CRS%JFEX1MX,
     .                   FP(1),FP(2),FP(3))

            EIRENE_CROSS = RES*FACT

C  ELAB IS WITHIN VALID RANGE OF FIT EXPRESSION
C  EVALUATE FIT EXPRESSION IFTFLG=3:
          ELSE
            B(1:8) = REACDAT(K)%CRS%POLY%DBLPOL(2:9,1)
cdr E .gt. XI is already fulfilled here.
            RES = B(1)*LOG(E/XI)
            DO I=1,7
              RES = RES + B(I+1)*(1.D0-XI/E)**I
            END DO
            RES = RES * 1.D-13/(XI*E)

            EIRENE_CROSS = RES*FACT

          ENDIF

        ELSE
          WRITE (iunout,*) ' WRONG FITTING FLAG IN CROSS'
          WRITE (iunout,*) ' K = ',K,' IFTFLG = ',IFTFLG(K,1)
          WRITE (iunout,*) 'REACTION NO. ',IR
          CALL EIRENE_EXIT_OWN(1)
        END IF

      ELSE  ! K out of range?
        WRITE (iunout,*) 'ERROR IN CROSS: K= ',K
        WRITE (iunout,*) 'CALLED FROM ',TEXT
        WRITE (iunout,*) 'REACTION NO. ',IR
        CALL EIRENE_EXIT_OWN(1)
      ENDIF

      RETURN
      END FUNCTION EIRENE_CROSS
