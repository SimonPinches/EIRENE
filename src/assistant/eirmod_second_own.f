      MODULE EIRMOD_SECOND_OWN
cdr  used for internal run time monitoring

      USE EIRMOD_PRECISION
cym
      use omp_lib
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
      
      time=omp_get_wtime()
      EIRENE_second_own=time-start
      END
C
      FUNCTION EIRENE_RESET_SECOND()
      implicit none
      real(dp) :: EIRENE_reset_second
cym   will need a dummy version for this
cym      call cpu_time(start)
      start=omp_get_wtime()
      EIRENE_reset_second=start
      return
      END

      END MODULE EIRMOD_SECOND_OWN
