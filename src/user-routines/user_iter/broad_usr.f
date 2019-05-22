

      SUBROUTINE EIRENE_BROAD_USR

      USE EIRMOD_PARMMOD
      USE EIRMOD_CTRIG
      USE EIRMOD_CPES
      use eirmod_extrab25
      use eirmod_mpi

      IMPLICIT NONE

      integer :: ier

      if (my_pe > 0) then
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
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

      RETURN
      END
