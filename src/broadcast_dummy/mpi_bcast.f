
      SUBROUTINE MPI_BCAST (buffer,cnt,datatype,root,comm,ier)
      IMPLICIT NONE
      integer, intent(out) :: ier
      real*8, pointer, INTENT(INOUT) :: buffer
      integer, intent(in) :: cnt,datatype,root,comm
      IER=0
      RETURN
      END
