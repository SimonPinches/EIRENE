
c
      FUNCTION EIRENE_SECOND_OWN()
      USE EIRMOD_PRECISION
      implicit none

      real(dp) :: EIRENE_reset_second, EIRENE_second_own
      real(sp) :: start,time
      save start
      data start /0.0/

      call cpu_time(time)
      EIRENE_second_own=time-start
      RETURN
C
      ENTRY EIRENE_RESET_SECOND

      call cpu_time(start)
      EIRENE_reset_second=start
      return
      END
