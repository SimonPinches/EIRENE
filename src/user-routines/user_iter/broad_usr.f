

      SUBROUTINE EIRENE_BROADCAST_USR(ME)

      USE EIRMOD_PARMMOD
      USE EIRMOD_CTRIG
      use eirmod_extrab25
!pbtest      USE EIRMOD_CPES
!pb     >    , ONLY : input_distribution_strategy
      use eirmod_mpi

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ME
      integer :: ier

      if (me > 0) then
         if (.not.allocated(plnxtri)) then
           allocate(plnxtri(ntrii))
           allocate(plnytri(ntrii))
           allocate(pplnxtri(ntrii))
           allocate(pplnytri(ntrii))
         end if
      end if

      CALL MPI_BCAST (plnxtri,ntrii,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (plnytri,ntrii,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (pplnxtri,ntrii,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (pplnytri,ntrii,MPI_REAL8,0,MPI_COMM_WORLD,ier)

!pbtest      CALL MPI_BCAST (input_distribution_strategy,1,MPI_INTEGER,0,
!pb     .                MPI_COMM_WORLD,ier)
c
csw 07apr2010
      call eirene_broadcast_infcop
csw

      CALL EIRENE_CHECK_EXIT
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

      RETURN
      END
