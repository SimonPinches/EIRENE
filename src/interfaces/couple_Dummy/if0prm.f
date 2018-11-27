cdr called from find_param.f in initialization phase,
cdr when eirene is in "coupled mode": i.e. IF(NMODE.NE.0)
C
cdr Read block 14 from interfacing routines (not from eirene_input.f)
c   this version: couple_dummy, i.e. only dummy interfacing routines.
c
c   and set storage for allocatable arrays:
c   NSTEP :
c   NPTRGT:
c   NAIN  :
c   NAOT  :
c   NCPV  :
c   NKNOT :
C   NTRII :
C   NCPVI :  no. of special couple tallies

c  also set: NDX,NDY,NFL, NDXP, NDYP

      SUBROUTINE EIRENE_IF0PRM(IUNIN,IUNOUT)

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IUNIN,IUNOUT

      NCPV=0
      NAIN=0
      NPTRGT=1

      END
