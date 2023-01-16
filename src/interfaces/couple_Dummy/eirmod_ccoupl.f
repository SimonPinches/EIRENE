cdr Aug 19:  added forgotten broadcast of iccpl2(:)
cdr          Perhaps obsolete code, but this way it works
cdr          also with MPI. Is NFLA really still needed in case of
cdr          dummy-interfacing routines? If not, code here
cdr          can be strongly simplified

      MODULE EIRMOD_CCOUPL

      USE EIRMOD_PARMMOD
      USE EIRMOD_MPI

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: EIRENE_ALLOC_CCOUPL, EIRENE_DEALLOC_CCOUPL,
     P          EIRENE_INIT_CCOUPL, EIRENE_BROADCAST_CCOUPL

      INTEGER, PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     I         ICCPL2(:)

      INTEGER, PUBLIC, POINTER, SAVE ::
     I NFLA

      INTEGER, PUBLIC, SAVE ::
     I MCOUPL2

      CONTAINS


      SUBROUTINE EIRENE_ALLOC_CCOUPL (ICAL)
cdr Here only deal with ICCPL2.
cdr NFLA is iccpl2(3) in corresponding non-default modules 

      INTEGER, INTENT(IN) :: ICAL

      IF (ICAL == 1) THEN

        IF (ALLOCATED(ICCPL2)) RETURN

        MCOUPL2 = 3
        ALLOCATE (ICCPL2(MCOUPL2))

        WRITE (IUNMEM,'(A,T25,I15)')
     .        ' CCOUPL ',MCOUPL2*4

        NFLA       => ICCPL2(3)


      ELSE IF (ICAL == 2) THEN

      END IF

      CALL EIRENE_INIT_CCOUPL (ICAL)

      RETURN
      END SUBROUTINE EIRENE_ALLOC_CCOUPL


      SUBROUTINE EIRENE_DEALLOC_CCOUPL

      IF (.NOT.ALLOCATED(ICCPL2)) RETURN

      DEALLOCATE (ICCPL2)

      RETURN
      END SUBROUTINE EIRENE_DEALLOC_CCOUPL


      SUBROUTINE EIRENE_INIT_CCOUPL(ICAL)

      INTEGER, INTENT(IN) :: ICAL

      IF (ICAL == 1) THEN

        ICCPL2 = 0

      ELSE IF (ICAL == 2) THEN

      END IF

      RETURN
      END SUBROUTINE EIRENE_INIT_CCOUPL


      SUBROUTINE EIRENE_BROADCAST_CCOUPL(ME)
      USE EIRMOD_COMUSR, ONLY : NMODE
      INTEGER, INTENT(IN) :: ME
      INTEGER :: IER

      CALL MPI_BCAST (NMODE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

      IF (NMODE == 0) RETURN 

      IF (ME /= 0) THEN
        CALL EIRENE_ALLOC_CCOUPL(1)
        CALL EIRENE_ALLOC_CCOUPL(2)
      END IF

      IF (ALLOCATED(ICCPL2)) THEN
        CALL MPI_BCAST (ICCPL2,MCOUPL2,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      END IF

      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

      END SUBROUTINE EIRENE_BROADCAST_CCOUPL 


      END MODULE EIRMOD_CCOUPL
