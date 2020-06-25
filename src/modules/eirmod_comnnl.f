!pb  31.10.06:  definition of census arrays RPART, RPARTC, IPART, IPARTC changed
!               RPART (NPRNL,NPARTT) --> RPART (NPARTT,NPRNL)
!               IPART (MPRNL,NPARTT) --> IPART (MPARTT,NPRNL)
!               RPARTC(NPRNL,NPARTT) --> RPARTC(NPARTT,NPRNL)
!               IPARTC(MPRNL,NPARTT) --> IPARTC(MPARTT,NPRNL)
cdr:  rpart, ipart:  "true census" arrays, real and integer.
cdr:  rpartc,ipartc: copies of census arrays, needed for sampling (bootstrapping)
cdr                  from old census, while already filling the new census
cdr   rpartw:  cumulated weight from census, set at the end of timestep (TMSTEP.f)
cdr            from the weights stored on rpart, for bootstrapping in next timestep
cdr            (in subr. LOCATE.f)


      MODULE EIRMOD_COMNNL

      USE EIRMOD_PRECISION, ONLY: DP
      USE EIRMOD_PARMMOD, ONLY: IUNMEM, MPARTT, NPARTT, NPRNL, NSTRA


      IMPLICIT NONE

      PRIVATE

      PUBLIC :: EIRENE_ALLOC_COMNNL, EIRENE_DEALLOC_COMNNL,
     P          EIRENE_INIT_COMNNL, EIRENE_BROADCAST_COMNNL

      REAL(DP), PUBLIC, SAVE ::
     R  DTIMV,  DTIMVI, DTIMVN, TIME0

      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R  RPART(:,:), RPARTC(:,:), RPARTW(:)

      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I  IPART(:,:), IPARTC(:,:),
     I  NPRNLS(:)

      INTEGER, PUBLIC, SAVE ::
     I  NPRNLI, IPRNLI, IPRNLS, IPRNL,
     I  NPTST,  NTMSTP, ITMSTP


      CONTAINS


      SUBROUTINE EIRENE_ALLOC_COMNNL

      IF (ALLOCATED(RPART)) RETURN

      ALLOCATE (RPART(NPARTT,NPRNL))
      ALLOCATE (RPARTC(NPARTT,NPRNL))
      ALLOCATE (RPARTW(0:NPRNL))

      ALLOCATE (IPART(MPARTT,NPRNL))
      ALLOCATE (IPARTC(MPARTT,NPRNL))
      ALLOCATE (NPRNLS(NSTRA))

      WRITE (IUNMEM,'(A,T25,I15)')
     .      ' COMNNL ',(2*NPRNL*NPARTT+NPRNL+1)*8 +
     .                 (2*NPRNL*MPARTT+NSTRA)*4

      CALL EIRENE_INIT_COMNNL

      RETURN
      END SUBROUTINE EIRENE_ALLOC_COMNNL


      SUBROUTINE EIRENE_DEALLOC_COMNNL

      IF (.NOT.ALLOCATED(RPART)) RETURN

      DEALLOCATE (RPART)
      DEALLOCATE (RPARTC)
      DEALLOCATE (RPARTW)

      DEALLOCATE (IPART)
      DEALLOCATE (IPARTC)
      DEALLOCATE (NPRNLS)

      RETURN
      END SUBROUTINE EIRENE_DEALLOC_COMNNL


      SUBROUTINE EIRENE_INIT_COMNNL

      RPART  = 0._DP
      RPARTC = 0._DP
      RPARTW = 0._DP

      IPART  = 0
      IPARTC = 0
      NPRNLS = 0

      NPRNLI = 0
      IPRNLI = 0
      IPRNLS = 0
      IPRNL  = 0
      NPTST  = 0
      NTMSTP = 0
      ITMSTP = 0

      RETURN
      END SUBROUTINE EIRENE_INIT_COMNNL


C> \brief Broadcast quantities of module EIRMOD_COMNNL
C>
C> All quantities of the module EIRMOD_COMNNL are broadcasted here.
C> The census arrays are only broadcasted in the parallelisation mode
C> with proportional allocation. Especially, in the embarrassingly
C> parallel mode, each process should keep its own census. (Attention:
C> write out of the census of each process for the restart of a run is
C> not (yet) implemented and may be part of the plasma code interface.)
      SUBROUTINE EIRENE_BROADCAST_COMNNL(ME)

      USE EIRMOD_PARMMOD, ONLY: MPARTT, NPARTT, NPRNL, NSTRA
      USE EIRMOD_COMUSR, ONLY: NPRLL
      USE EIRMOD_MPI

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ME
      INTEGER :: IER

      IF (ME /= 0) CALL EIRENE_ALLOC_COMNNL
      
      CALL MPI_BCAST (DTIMV,1,MPI_REAL8,0,MPI_COMM_WORLD,IER)
      CALL MPI_BCAST (DTIMVI,1,MPI_REAL8,0,MPI_COMM_WORLD,IER)
      CALL MPI_BCAST (DTIMVN,1,MPI_REAL8,0,MPI_COMM_WORLD,IER)
      CALL MPI_BCAST (TIME0,1,MPI_REAL8,0,MPI_COMM_WORLD,IER)
      CALL MPI_BCAST (NPRNLI,1,MPI_INTEGER,0,MPI_COMM_WORLD,IER)
      CALL MPI_BCAST (IPRNLI,1,MPI_INTEGER,0,MPI_COMM_WORLD,IER)
      CALL MPI_BCAST (IPRNLS,1,MPI_INTEGER,0,MPI_COMM_WORLD,IER)
      CALL MPI_BCAST (IPRNL ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IER)
      CALL MPI_BCAST (NPTST ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IER)
      CALL MPI_BCAST (NTMSTP,1,MPI_INTEGER,0,MPI_COMM_WORLD,IER)
      CALL MPI_BCAST (ITMSTP,1,MPI_INTEGER,0,MPI_COMM_WORLD,IER)
      
      IF ( NPRLL == 1 ) THEN
        CALL MPI_BCAST (RPART,NPRNL*NPARTT,MPI_REAL8,
     >                  0,MPI_COMM_WORLD,IER)
        CALL MPI_BCAST (RPARTC,NPRNL*NPARTT,MPI_REAL8,
     >                  0,MPI_COMM_WORLD,IER)
        CALL MPI_BCAST (RPARTW,NPRNL+1,MPI_REAL8,0,MPI_COMM_WORLD,IER)

        CALL MPI_BCAST (IPART,NPRNL*MPARTT,MPI_INTEGER,
     >                  0,MPI_COMM_WORLD,IER)
        CALL MPI_BCAST (IPARTC,NPRNL*MPARTT,MPI_INTEGER,
     >                  0,MPI_COMM_WORLD,IER)
        CALL MPI_BCAST (NPRNLS,NSTRA,MPI_INTEGER,0,MPI_COMM_WORLD,IER)
      END IF

      RETURN

      END SUBROUTINE EIRENE_BROADCAST_COMNNL

      END MODULE EIRMOD_COMNNL
