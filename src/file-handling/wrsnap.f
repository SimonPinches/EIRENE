!pb  26.10.06: close file after read or write

cdr  2020
cdr: FLUX(istr), for istr=NSTRAI (scaling for stratum ISTR) is also on fort.13.
cdr: This may be a hidden link. Decouple flux(nstrai) written on fort.13 from
cdr  that written on fort.15: FLXCEN, to reduce code complexity.

C
      SUBROUTINE EIRENE_WRSNAP( ISTR )
C
C  SAVE SNAPSHOT PARAMETERS IPRNL,FLXCEN,DTIMV, AT END OF TIMESTEP, ON FORT.15
C  SAVE SNAPSHOT POPULATION RPARTC,IPARTC, AT END OF TIMESTEP, ON FORT.15
C
      USE EIRMOD_PARMMOD, ONLY: IFOFF, MPARTT, NPARTT
      USE EIRMOD_CTRCEI, ONLY: TRCFLE
      USE EIRMOD_COMNNL, ONLY: FLXCEN, DTIMV, IPRNL,
     .                         IPARTC, RPARTC, RPARTW
      USE EIRMOD_COMSOU, ONLY: FLUX
      USE EIRMOD_COMPRT, ONLY: IUNOUT


      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ISTR
      INTEGER :: I, J
c  (atomic) census flux
      FLXCEN=FLUX(ISTR)
C
      OPEN (UNIT=15+ifoff,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 15+ifoff
C
      IF (TRCFLE) WRITE (iunout,*) 'WRITE 15: IPRNL,FLXCEN,DTIMV'
      WRITE (15+ifoff) IPRNL,FLXCEN,DTIMV
      WRITE (15+ifoff) ((RPARTC(J,I),J=1,NPARTT),I=1,IPRNL)
      WRITE (15+ifoff)  (RPARTW(  I)            ,I=0,IPRNL)
      WRITE (15+ifoff) ((IPARTC(J,I),J=1,MPARTT),I=1,IPRNL)
      CLOSE (UNIT=15+ifoff)
C
      RETURN
      END SUBROUTINE EIRENE_WRSNAP
C
      SUBROUTINE EIRENE_RSNAP( ISTR )
C
C  READ SNAPSHOT PARAMETERS IPRNL,FLXCEN,DTIMV, FROM FORT.15
C  READ SNAPSHOT POPULATION RPARTC,IPARTC, FROM FORT.15
C  If fort.15 is absent, define an empty census (zero flux and
C  zero scores).
C
      USE EIRMOD_PRECISION, ONLY: DP
      USE EIRMOD_PARMMOD, ONLY: IFOFF, MPARTT, NPARTT, NPRNL
      USE EIRMOD_CTRCEI, ONLY: TRCFLE
      USE EIRMOD_COMNNL, ONLY: FLXCEN, DTIMV, IPRNL, DTIMVN,
     ,                         IPARTC, IPRNL, RPARTC, RPARTW
      USE EIRMOD_COMSOU, ONLY: FLUX
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CINIT, ONLY: FORT

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ISTR
      INTEGER :: I, J
C
      OPEN (UNIT=15+ifoff,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 15+ifoff
      IF (TRCFLE) WRITE (iunout,*) 'READ 15: IPRNL,FLXCEN,DTIMV'
      READ (15+ifoff,END=915) IPRNL,FLXCEN,DTIMV

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

      FLUX(ISTR)=FLXCEN
C
      RETURN
C
CXPB  SAFETY ADDED FOR CASE WHEN FORT.15 IS ABSENT
  915 CONTINUE
      CLOSE (UNIT=15+ifoff)
      IPRNL=0
      FLXCEN=0.0_DP
      DTIMV=DTIMVN
      RPARTW(0)=0.0_DP

      FLUX(ISTR)=FLXCEN
      WRITE (IUNOUT,*) 'FILE '//FORT//'15 WAS FOUND EMPTY'
      WRITE (IUNOUT,*) 'CREATING AN EMPTY TIME CENSUS ARRAY'
      CALL EIRENE_LEER(1)
      RETURN
      END SUBROUTINE EIRENE_RSNAP
