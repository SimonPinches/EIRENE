
      SUBROUTINE MPI_send (buffer,cnt,datatype,dest,tag,comm,ier)
      IMPLICIT NONE
      integer, intent(out) :: ier
      real*8, pointer, INTENT(INOUT) :: buffer
      integer, intent(in) :: cnt,datatype,dest,tag,comm
      IER=0
      RETURN
      END
