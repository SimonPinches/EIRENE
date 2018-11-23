C
C   USER-SUPPLIED SUBROUTINES
C

C
      SUBROUTINE EIRENE_PROUSR (PRO,INDX,P0,P1,P2,P3,P4,P5,PROVAC,N)

c  user defined background (plasma) profiles.
c  called from routine plasma.f, indpro=5 option
c  input:  indx  : indicates which input tally (e.g. Te, ni, vx, etc...)
c          provac: eirene default vacuum value (e.g. dvac, tvac, etc...)
c          N     : N=NSURF, standard grid without additional cell region
C
C  EXAMPLE:
C     P0 : CENTRAL VALUE
C     P1 : STARTING RADIUS FOR POLYNOMIAL
C     P2 : SWITCH FOR PHASE 1: OH-PHASE
C                           2: NI-PHASE
C     P3 : FACTOR FOR TI: TI=K*TE
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: P0, P1, P2, P3, P4, P5, PROVAC
      REAL(DP), INTENT(OUT) :: PRO(*)
      INTEGER, INTENT(IN) :: INDX, N

      PRO(1:N)=0.0_DP

      RETURN
      END
