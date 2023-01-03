!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced
cdr  27.09.15:  renamed from fehvds1 to fehvei1,
cdr  KK<0 cases still to be synchronized with ftabei1, feelei1, etc.
cdr  nelrei, nhvrei,
cdr  01.01.16:  documented
cdr             more unified numbering of reactions: kk= nreaei(irei),
cdr             for default reactions, kk<0.

      FUNCTION EIRENE_FEHVEI1 (IREI,K)
C  this is the "on the fly" storage saving version to eliminate
C  pre-computed array EHVEI1(irei,k) from with run

cdr  find heavy secondary particle energy (=KER) for EI process no. IREI,  energy in eV
c    locally in cell K, for process kk= nhvrei(irei)
c    sum over all heavy secondaries.
c    Distribution to individual heavy secondary type and species is done later,
c    e.g. in veloei for sampling, and in UPDATE, COLLIDE,.. for scoring
c
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CCONA
      USE EIRMOD_COMXS
      USE EIRMOD_COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IREI, K
      REAL(DP) :: EIRENE_FEHVEI1, EHVEI, EIRENE_FTABEI1,
     .            EIRENE_ENERGY_RATE_COEFF, DE_10
      INTEGER :: KK

      EIRENE_FEHVEI1=0.D0

C  IDENTIFY NUMBER OF PROCESS.
C  CURRENTLY KK=-11 -- KK=-4 EIRENE DEFAULT PROCESSES
C            KK> 0  KREAD: COLLISION PROCESSES STORED ON REACDAT FROM EXTERNAL DATABASES
      KK=NHVREI(IREI)


      IF (KK < 0) THEN
cdr in this case (default reactions) we have: nhvrei = nreaei
        SELECT CASE (KK)
        CASE (-4)   ! DEFAULT PROCESS KK=-4: H + E --> H+ + 2E, no net energy transfer to H+
           EIRENE_FEHVEI1 =0.0
cdr  now the 6 default reactions for H2, H2+
        CASE (-5)  ! DEFAULT PROCESS KK=-5:  H2+E --> H+H +E,
            EIRENE_FEHVEI1=6.
        CASE (-6)  ! DEFAULT PROCESS KK=-6: H2 + E --> H + p + 2E,
            EIRENE_FEHVEI1=10.0
        CASE (-7)  ! DEFAULT PROCESS KK=-7: H2 + E --> H2+ + 2E,
            EIRENE_FEHVEI1=0.0
        CASE (-8)   ! DEFAULT PROCESS KK=-8:  H2+ + E --> H+ + H +E,
            EIRENE_FEHVEI1=8.6
        CASE (-9)   ! DEFAULT PROCESS KK=-9:  H2+ + E --> H+ + H+ +E,
            EIRENE_FEHVEI1=0.8
        CASE (-10)   ! DEFAULT PROCESS KK=-10: H2+ + E --> H + H, DISS. RECOMBINATION
C  FOR THE FACTOR -0.896... SEE: EIRENE MANUAL, INPUT BLOCK 4, EXAMPLES
C  FOR THIS PROCESS: RADIATION RATE = -POTENTIAL ENERGY RATE: DP= 15.6-4.52=11.08
cdr  this DP is a bit too small, because H2+ is also vibr. excited.
cdr  additional electron cost for H2+(0) --> H2+(v) may also be radiated.
            DE_10=8.964355004318D-01
            EIRENE_FEHVEI1=DE_10*TEIN(K)
        CASE (-11)   ! DEFAULT PROCESS KK=-11:  He+ E --> He+ +E, no net energy transfer to He+!
            EIRENE_FEHVEI1 =0.0
        CASE DEFAULT
            GOTO 999
        END SELECT

c  non-default models, data from external databases
      ELSE IF (KK > 0) THEN
cdr flag jhvrei is still missing. how to we know whether to use density parameter pls ?
cdr this code corresponds to jhvrei=1
        EHVEI = EIRENE_ENERGY_RATE_COEFF(KK,K,TEINL(K),0._DP,.FALSE.,0)
        EHVEI=EXP(MAX(-100._DP,EHVEI+FACREI(IREI,2)))
        EIRENE_FEHVEI1=EHVEI*DEIN(K)/(EIRENE_FTABEI1(IREI,K)+EPS60)
c  else: jhvrei=9:   density-dependent KER, to be written

      ELSE IF (KK == 0) THEN  ! flag jhvrei is still missing.
         EIRENE_FEHVEI1 = EHVEI1(IREI,1) ! currently: only constant KER option

      END IF

      RETURN

  999 CONTINUE
      WRITE (IUNOUT,*) 'FEHVEI1: INVALID PARAMETER NHVREI '
      WRITE (IUNOUT,*) 'IREI, NHVREI ',IREI,NHVREI
      CALL EIRENE_EXIT_OWN(1)
      END  FUNCTION EIRENE_FEHVEI1
