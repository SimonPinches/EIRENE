      MODULE EIRMOD_SECOND_OWN
cdr  used for internal run time monitoring

      USE EIRMOD_PRECISION
#ifdef USE_OPENMP 
      use omp_lib
#endif
#ifdef USE_MPI      
      use mpi
#endif
      IMPLICIT NONE
      PRIVATE

      PUBLIC :: EIRENE_SECOND_OWN, EIRENE_RESET_SECOND

      real(sp),save :: start=0.0

      CONTAINS
c
      FUNCTION EIRENE_SECOND_OWN()
cdr returns wall clock time (in seconds) since the last call to function eirene_reset_second.
      implicit none

      real(dp) :: EIRENE_second_own
      real(dp) :: time
#ifdef USE_OPENMP 
      time=omp_get_wtime()
#elif USE_MPI
      time=mpi_wtime()
#else
      call cpu_time(time)      
#endif
      
      EIRENE_second_own=time-start
      return
      END FUNCTION EIRENE_SECOND_OWN
C
      FUNCTION EIRENE_RESET_SECOND()
cdr returns wall clock time (in seconds)
      implicit none
      real(dp) :: EIRENE_reset_second
#ifdef USE_OPENMP 
      start=omp_get_wtime()
#elif USE_MPI
      start=mpi_wtime()
#else
      call cpu_time(start)
#endif
      EIRENE_reset_second=start
      return
      END FUNCTION EIRENE_RESET_SECOND

      END MODULE EIRMOD_SECOND_OWN
