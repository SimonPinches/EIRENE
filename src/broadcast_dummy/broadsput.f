cdr March 2016:  dimensioning of sputer parameters, second index: N2 --> 0:N2
cdr             N2 is the index for the target (wall) material
cdr             i2=0 is storage for those sputer parameters 
cdr             which are missing in the data tables and which are 
cdr             evaluated "on the fly" instead (subr. sputer.f)  
 
      SUBROUTINE EIRENE_BROADSPUT(ES,M2M1,ETF,ETH,Q,N1,N2)
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      IMPLICIT NONE
      INTEGER, INTENT(IN):: N1, N2
      REAL(DP), INTENT(IN OUT) :: ES(N1), M2M1(N1,0:N2), ETF(N1,0:N2),
     .                            ETH(N1,0:N2), Q(N1,0:N2)
      RETURN
      END
