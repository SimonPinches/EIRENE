!pb  26.10.06: close file after read or write
!pb  31.10.06: definition of census arrays RPART, RPARTC, IPART, IPARTC changed
cdr:  2015
!               RPART (NPARTT,NPRNL) (now) <-- RPART (NPRNL,NPARTT) (formerly)
!               RPARTC(NPARTT,NPRNL) (now) <-- RPARTC(NPRNL,NPARTT) (formerly)
!               IPART (MPARTT,NPRNL) (now) <-- IPART (NPRNL,MPARTT) (formerly)
!               IPARTC(MPARTT,NPRNL) (now) <-- IPARTC(NPRNL,MPARTT) (formerly)
C
C
      SUBROUTINE EIRENE_WRSNAP( ISTR )
C
C  SAVE SNAPSHOT POPULATION RPARTC,IPARTC, AT END OF TIMESTEP, ON FORT.15
C
      USE EIRMOD_PARMMOD, ONLY: IFOFF, MPARTT, NPARTT
      USE EIRMOD_CTRCEI, ONLY: TRCFLE
      USE EIRMOD_COMNNL, ONLY: DTIMV, IPARTC, IPRNL, RPARTC, RPARTW
      USE EIRMOD_COMSOU, ONLY: FLUX
      USE EIRMOD_COMPRT, ONLY: IUNOUT


      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ISTR
      INTEGER :: I, J
C
      OPEN (UNIT=15+ifoff,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 15+ifoff
C
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 15: IPRNL,FLUX,DTIMV'
      WRITE (15+ifoff) IPRNL,FLUX(ISTR),DTIMV
      WRITE (15+ifoff) ((RPARTC(J,I),J=1,NPARTT),I=1,IPRNL)
      WRITE (15+ifoff)  (RPARTW(  I)            ,I=0,IPRNL)
      WRITE (15+ifoff) ((IPARTC(J,I),J=1,MPARTT),I=1,IPRNL)
      CLOSE (UNIT=15+ifoff)
C
      END SUBROUTINE EIRENE_WRSNAP
C
      SUBROUTINE EIRENE_RSNAP( ISTR )
      USE EIRMOD_PARMMOD, ONLY: IFOFF, MPARTT, NPARTT, NPRNL
      USE EIRMOD_CTRCEI, ONLY: TRCFLE
      USE EIRMOD_COMNNL, ONLY: DTIMV, IPARTC, IPRNL, RPARTC, RPARTW
      USE EIRMOD_COMSOU, ONLY: FLUX
      USE EIRMOD_COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ISTR
      INTEGER :: I, J
C
      OPEN (UNIT=15+ifoff,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 15+ifoff
      READ (15+ifoff) IPRNL,FLUX(ISTR),DTIMV

      IF (TRCFLE) WRITE (iunout,*) 'READ 15: IPRNL,FLUX,DTIMV'

cdr  tbd:
cdr  fort.15 (census) was written in a previous run.
cdr  In the present run allocation of storage for census arrays rpartc,ipartc,rpartw
cdr  is determined by input: --> nprnl
cdr  make sure that iprnl in previous run was not larger than in present run.
      IF (NPRNL.LT.IPRNL) THEN
        WRITE (IUNOUT,*)
     .     ' ERROR WHEN READING CENSUS ARRAY FOR T-DEP MODE'
        WRITE (IUNOUT,*)
     .     ' OLD CENSUS FILE CANNOT BE READ, BECAUSE NPRNL TOO SMALL'
        CALL EIRENE_MASJ2(' NPRNL, IPRNL=  ',NPRNL,IPRNL)
        CALL EIRENE_EXIT_OWN(1)
      ENDIF

      READ (15+ifoff) ((RPARTC(J,I),J=1,NPARTT),I=1,IPRNL)
      READ (15+ifoff)  (RPARTW(  I)            ,I=0,IPRNL)
      READ (15+ifoff) ((IPARTC(J,I),J=1,MPARTT),I=1,IPRNL)
      CLOSE (UNIT=15+ifoff)
C
      RETURN
      END
