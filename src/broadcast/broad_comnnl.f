C> \brief Broadcast quantities of module EIRMOD_COMNNL
C>
C> All quantities of the module EIRMOD_COMNNL are broadcasted here.
C> The census arrays are only broadcasted in the parallelisation mode
C> with proportional allocation. Especially, in the embarrassingly
C> parallel mode, each process should keep its own census. (Attention:
C> write out of the census of each process for the restart of a run is
C> not (yet) implemented and may be part of the plasma code interface.)
      SUBROUTINE EIRENE_BROAD_COMNNL

      USE EIRMOD_COMNNL
      USE EIRMOD_COMUSR, ONLY: NPRLL
      USE EIRMOD_MPI
      USE EIRMOD_PARMMOD, ONLY: MPARTT, NPARTT, NPRNL, NSTRA

      IMPLICIT NONE

      INTEGER :: IER

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

      END SUBROUTINE EIRENE_BROAD_COMNNL

