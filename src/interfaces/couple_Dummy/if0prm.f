cdr called from find_param.f in initialization phase. Read block 14
c   and set storage for allocatable arrays:
c   NPTRGT:
c   NAIN  :
c   NCPV  :
c   NKNOT :
C   NTRII :
C   NCPVI :  no. of special couple tallies    

      SUBROUTINE EIRENE_IF0PRM(IUNIN,IUNOUT)
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IUNIN,IUNOUT
 
      NCPV=0
      NAIN=0
      NPTRGT=1      
 
      END
