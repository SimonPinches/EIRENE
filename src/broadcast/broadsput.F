cdr  March 2016 parameters for sputter formula, second index: N2  --> 0:N2.
cdr             N2 is the index for the target (wall) material
cdr             i2=0 is storage for those sputter parameters
cdr             which are missing in the data tables and which are
cdr             evaluated "on the fly" instead (subr. sputer.f)

      SUBROUTINE EIRENE_BROADSPUT(ES,M2M1,ETF,ETH,Q,N1,N2)
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_MPI
      IMPLICIT NONE

      INTEGER, INTENT(IN):: N1, N2
      REAL(DP), INTENT(IN OUT) :: ES(N1), M2M1(N1,0:N2), ETF(N1,0:N2),
     .                            ETH(N1,0:N2), Q(N1,0:N2)
      INTEGER :: IER

      CALL EIRENE_CHECK_EXIT
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

      CALL MPI_BCAST (ES,N1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (M2M1,N1*(N2+1),MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ETF,N1*(N2+1),MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ETH,N1*(N2+1),MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (Q,N1*(N2+1),MPI_REAL8,0,MPI_COMM_WORLD,ier)

      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

      RETURN
      END SUBROUTINE EIRENE_BROADSPUT
