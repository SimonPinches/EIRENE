      MODULE EIRMOD_EIRBRA

C  NEUTRAL SOURCE TERMS: SNI,SMO,SEE,SEI (EIRENE ---> BRAAMS)

      USE EIRMOD_PRECISION
      USE EIRMOD_MPI
      IMPLICIT NONE

      PRIVATE

      PUBLIC :: EIRENE_ALLOC_EIRBRA, EIRENE_DEALLOC_EIRBRA,
     P          EIRENE_INIT_EIRBRA
csw mpi
      public :: eirene_broadcast_eirbra
csw

      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R SNI(:,:,:,:), SMO(:,:,:,:),
     R SEE(:,:,:),   SEI(:,:,:)

      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R VOLSUMN(:), VOLSUMM(:), VOLSUMEI(:), VOLSUMEE(:)

      INTEGER, SAVE :: NDXD, NDYD, NFLD, NSTRAD

C AK
      REAL(DP),PUBLIC,ALLOCATABLE,SAVE ::
     .                srcstrn(:)
!pb  .               ,flxspci(:,:)
!pb  .               ,srccrfc(:,:)
c*** srcstrn: intensities of different neutral sources (fluxt)
c*** flxspci: fluxes of plasma ions to the recycling surfaces
c*** srccrfc: source correction factors (from infcop)

C AK END

      CONTAINS


      SUBROUTINE EIRENE_ALLOC_EIRBRA(NDXP, NDYP, NFL, NSTRA)

      USE EIRMOD_PARMMOD, ONLY: IUNMEM

      INTEGER, INTENT(IN) :: NDXP, NDYP, NFL, NSTRA

      IF (ALLOCATED(SNI)) RETURN

      NDXD = NDXP+1
      NDYD = NDYP+1
      NFLD = NFL
!pb  obviously SOLPS always assumes there is a time stratum
!pb  therefore add additional space in arrays
csw   NSTRAD = NSTRA
      NSTRAD = NSTRA + 1

      ALLOCATE (SNI(0:NDXD,0:NDYD,NFLD,NSTRAD))
      ALLOCATE (SMO(0:NDXD,0:NDYD,NFLD,NSTRAD))
      ALLOCATE (SEE(0:NDXD,0:NDYD,NSTRAD))
      ALLOCATE (SEI(0:NDXD,0:NDYD,NSTRAD))

      ALLOCATE (VOLSUMN(NSTRAD))
      ALLOCATE (VOLSUMM(NSTRAD))
      ALLOCATE (VOLSUMEI(NSTRAD))
      ALLOCATE (VOLSUMEE(NSTRAD))

      ALLOCATE (SRCSTRN(NSTRAD))

      WRITE (IUNMEM,'(A,T25,I15)')
     .             ' EIRBRA ',((NDXD+1)*(NDYD+1)*NSTRAD*2*NFLD+
     .                        5*NSTRAD)*8

      CALL EIRENE_INIT_EIRBRA

      RETURN
      END SUBROUTINE EIRENE_ALLOC_EIRBRA


      SUBROUTINE EIRENE_DEALLOC_EIRBRA

      IF (.NOT.ALLOCATED(SNI)) RETURN

      DEALLOCATE (SNI)
      DEALLOCATE (SMO)
      DEALLOCATE (SEE)
      DEALLOCATE (SEI)

      DEALLOCATE (VOLSUMN)
      DEALLOCATE (VOLSUMM)
      DEALLOCATE (VOLSUMEI)
      DEALLOCATE (VOLSUMEE)

      RETURN
      END SUBROUTINE EIRENE_DEALLOC_EIRBRA


      SUBROUTINE EIRENE_INIT_EIRBRA

      SNI = 0._DP
      SMO = 0._DP
      SEE = 0._DP
      SEI = 0._DP

      VOLSUMN  = 0._DP
      VOLSUMM  = 0._DP
      VOLSUMEI = 0._DP
      VOLSUMEE = 0._DP

      RETURN
      END SUBROUTINE EIRENE_INIT_EIRBRA


csw mpi
      subroutine eirene_broadcast_eirbra

      integer :: inum, ierr

      inum = (ndxd+1)*(ndyd+1)*nfld*nstrad
      call mpi_bcast (sni,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
      call mpi_bcast (smo,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)

      inum = (ndxd+1)*(ndyd+1)*nstrad
      call mpi_bcast (see,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
      call mpi_bcast (sei,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)

      inum = nstrad
      call mpi_bcast (volsumn,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
      call mpi_bcast (volsumm,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
      call mpi_bcast (volsumee,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
      call mpi_bcast (volsumei,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
      end subroutine eirene_broadcast_eirbra
csw
      END MODULE EIRMOD_EIRBRA
