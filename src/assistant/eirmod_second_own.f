      MODULE EIRMOD_SECOND_OWN
cdr  used for internal run time monitoring

      USE EIRMOD_PRECISION
cym
#ifdef USE_OPENMP 
      use omp_lib
#endif
#ifdef USE_MPI      
      use mpi
#endif
cym
      IMPLICIT NONE
      PRIVATE

      PUBLIC :: EIRENE_SECOND_OWN, EIRENE_RESET_SECOND

      real(sp),save :: start=0.0

      CONTAINS
c
      FUNCTION EIRENE_SECOND_OWN()
      implicit none

      real(dp) :: EIRENE_second_own
cym      real(sp) :: time
cym   will need a dummy version for this
cym      call cpu_time(time)
      real(dp) :: time
#ifdef USE_OPENMP 
      time=omp_get_wtime()
#elif USE_MPI
      time=mpi_wtime()
#endif
      
      EIRENE_second_own=time-start
      END
C
      FUNCTION EIRENE_RESET_SECOND()
      implicit none
      real(dp) :: EIRENE_reset_second
cym   will need a dummy version for this
cym      call cpu_time(start)
#ifdef USE_OPENMP 
      start=omp_get_wtime()
#elif USE_MPI
      start=mpi_wtime()
#endif
      EIRENE_reset_second=start
      return
      END

      END MODULE EIRMOD_SECOND_OWN
