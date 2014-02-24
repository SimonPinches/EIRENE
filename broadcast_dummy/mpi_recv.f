 
      SUBROUTINE MPI_recv (buffer,cnt,datatype,source,tag,comm,st,ier)
      IMPLICIT NONE
      integer, intent(out) :: ier,st
      real*8, pointer, INTENT(INOUT) :: buffer
      integer, intent(in) :: cnt,datatype,source,tag,comm
      IER=0
      RETURN
      END
