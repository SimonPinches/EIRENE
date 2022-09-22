      module eirmod_sheath
cym 04/2020 turned into module beacuse of save 
      use eirmod_precision
      implicit none
      private

      public :: eirene_sheath

cym make sure this works fine
      integer,save :: icount=0

!$OMP THREADPRIVATE(icount)
C
      contains
C
      FUNCTION EIRENE_SHEATH(TE,DPP,VP,ZNZP,GAMMA,CUR,NP,MS)

C  JAN.93: NEW VERSION, NOW CORRECTLY ACCOUNTING
C          FOR MULTISPECIES PLASMA BACKGROUND
C  SHEATH POTENTIAL -E*DPHI (EV) FOR AN NP COMPONENT
C  PLASMA AS FUNCTION OF:
C  TE ELECTRON TEMPERATURE (EV)
C  DPP(J),J=1,NP BULK ION DENSITIES (#/CM**3)
C  VP (J),J=1,NP BULK ION DRIFT VELOCITIES (CM/SEC)
C                COMPONENT PARALLEL TO MAGNETIC FIELD B
C                (ASSUME: MOBILITY OF ELECTRONS MUCH LARGER
C                THAN OF IONS, AND ELECTRONS OBEY A
C                MAXWELL-BOLTZMANN DISTRIBUTION)
C  ZNZP(J),J+1,NP AVERAGE CHARGE OF BULK IONS (cnh 02.11.2019, integer NZP -> dble ZNZP)
C  GAMMA         SECONDARY ELECTRON EMISSION COEF.
C                SHEATH VOLTAGE
C                "SHEATH" DECREASES WITH INCREASING GAMMA
C  CUR           NET ELECTR. CURRENT DENSITY OF BULK PLASMA TO TARGET
C                CUR=(J-ION) - (J-ELECTR.) (AMP/CM**2)
C                "SHEATH" DECREASES FOR NEGATIVE CURRENTS, AND
C                         INCREASES FOR POSITIVE CURRENTS
C  MS            SURFACE NUMBER (ONLY FOR PRINTOUT,
C                                IN CASE INVALID SHEATH PARAMETERS AT THIS SURFACE)
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CCONA
      USE EIRMOD_COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NP, MS
      REAL(DP), INTENT(IN) :: DPP(NP),VP(NP), TE, GAMMA, CUR, ZNZP(NP)
      REAL(DP) :: SUM, DE, CE, EIRENE_SHEATH
      INTEGER :: J, MSS

cccccccccccccccccccccccccccccccccccc
cym      INTEGER :: J, MSS, ICOUNT
cym      SAVE ICOUNT
cym      DATA ICOUNT/0/
cccccccccccccccccccccccccccccccccccc
C
      EIRENE_SHEATH=0.
      DE=0.
      DO 10 J=1,NP
        DE=DE+DPP(J)*ZNZP(J)
   10 CONTINUE
      DE=DE+EPS60
C  UNITS OF SUM: VELOCITY (CM/SEC)
      SUM=0.
      DO 100 J=1,NP
        SUM=SUM+ZNZP(J)*DPP(J)/DE*VP(J)
  100 CONTINUE

      CE=CVEL2A*SQRT(TE/PMASSE)
      SUM=1./CE*SQRT(PI2A)/(1.-GAMMA)*(SUM-CUR/ELCHA/DE)
      IF (SUM.GT.0.D0) THEN
        EIRENE_SHEATH=-TE*LOG(SUM)
      ELSEIF (ICOUNT.LE.10 .AND. MS.NE.0) THEN
        MSS=MS
        IF (MSS.GT.NLIM) MSS=-(MSS-NLIM)
        WRITE (iunout,*) 'WARNING FROM FCT. SHEATH: INVALID ARGUMENTS '

        WRITE (iunout,*) 'SHEATH RETURNED FOR SURFACE ', MSS,': 2.8*TE'
        WRITE (iunout,*) 'NUMBER OF ION SPECIES FLOWS INTO SHEATH: ',NP
        WRITE (iunout,*) 'NUMBER,        CHARGE,'//
     &              '     ION DENS.,      ION VEL. '
        DO J=1,NP
          WRITE (iunout,60) J,ZNZP(J),DPP(J),VP(J)
        ENDDO
        CALL EIRENE_MASR1('I-CURR  ',SUM)
        ICOUNT=ICOUNT+1
        EIRENE_SHEATH = 2.8*TE
      ELSE
        EIRENE_SHEATH = 2.8*TE
      ENDIF

   60 FORMAT (1X,I6,3X,3(1PE12.4,3X))
      RETURN
      END FUNCTION EIRENE_SHEATH

      end module eirmod_sheath
