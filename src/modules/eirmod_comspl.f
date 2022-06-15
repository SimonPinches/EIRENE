!pb  31.10.06:  definition of RSPLST, ISPLST changed
!               RSPLST(15,NPARTC) --> RSPLST(NPARTC,15)
!               ISPLST(15,MPARTC) --> ISPLST(MPARTC,15)
!    20.06.07:  MAXLEVEL=15 defined, used in dimensioning of RSPLST and ISPLST
cdr  oct 19:
cdr             MAXLEVEL is the "depth" of splitting cascades.
cdr             Formerly: MAXLEVEL=15
cdr             now (2013) hard-coded: MAXLEVEL=300, why?
cdr             Is this intended indeed?
cpb             MAXLEVEL was increased for NLCASCAD option. If all secondary
cpb             particles of a reaction are to be traced the storage for
cpb             storing not yet followed secondaries needs to be somewhat larger.

      MODULE EIRMOD_COMSPL

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD


      IMPLICIT NONE

      PRIVATE

      PUBLIC :: EIRENE_ALLOC_COMSPL, EIRENE_DEALLOC_COMSPL,
     P          EIRENE_INIT_COMSPL, EIRENE_BROADCAST_COMSPL

      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE :: RCMSPL(:)

      REAL(DP), PUBLIC, POINTER, SAVE ::
     R WMINV,    WMINS,     WMINC,     WMINL,     SPLPAR,
     R RNUMB(:), PRMSPL(:)

      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R RSPLST(:,:)

      INTEGER, PUBLIC, PARAMETER :: MAXLEVEL=300

      INTEGER, PUBLIC, TARGET, ALLOCATABLE, SAVE :: ICMSPL(:)

      INTEGER, PUBLIC, POINTER, SAVE ::
     I MAXLEV, NLEVEL,
     I MAXRAD, MAXPOL, MAXTOR, MAXADD,
     I NODES(:), NSSPL(:)

      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I ISPLST(:,:)

      LOGICAL, PUBLIC, TARGET, ALLOCATABLE, SAVE :: LCMSPL(:)

      LOGICAL, PUBLIC, POINTER, SAVE ::
     L NLSPLT(:),  ! indicate non dev. surfacs as "splitting-rr" surfaces

     L NLPRCA(:), NLPRCM(:), NLPRCI(:), NLPRCPH(:)
cdr  NLPRCS should also become POINTER, belongs to NLPRCA; ..., cond exp. est.
      LOGICAL, PUBLIC, ALLOCATABLE, SAVE ::
     L NLPRCS(:)  ! indicate additional surfaces as attractors for cond. exp. est.

cym PRIVATE-> PUBLIV      
      INTEGER, PUBLIC, SAVE ::
     I NCMSPL, MCMSPL, KCMSPL

!$OMP  THREADPRIVATE(WMINV,WMINS,WMINC,WMINL,SPLPAR,RNUMB,PRMSPL,
!$OMP& RSPLST,MAXLEV,NLEVEL,MAXRAD,MAXPOL,MAXTOR,MAXADD,NODES,NSSPL,
!$OMP& ISPLST,
cym again not private 18/04 
cym nlsplt,nlprca,nlprcm,nlprci,nlprcph,nlprcs,
cym - not private, initialized in this module
cym!$omp& ncmspl,mcmspl,kcmspl,
cym again not private 18/04
cym!$omp& lcmspl,
!$OMP& RCMSPL,ICMSPL)

cdr  end threadprivate here

      CONTAINS


      SUBROUTINE EIRENE_ALLOC_COMSPL

      IF (ALLOCATED(RCMSPL)) RETURN

      NCMSPL=2*(N1ST+N2ND+N3RD+NLIM)+5
      MCMSPL=N1ST+N2ND+N3RD+NLIM+MAXLEVEL+6

      KCMSPL=N1ST+N2ND+N3RD+NLIM+   ! splitting surfaces
     .       NATM+NMOL+NION+NPHOT   ! cond exp. est., by species
cdr  .      +NLIMPS                 ! for NLPRCS, tbd.

      ALLOCATE (RCMSPL(NCMSPL))
      ALLOCATE (ICMSPL(MCMSPL))
      ALLOCATE (LCMSPL(KCMSPL))

      ALLOCATE (RSPLST(NPARTC,MAXLEVEL))
      ALLOCATE (ISPLST(MPARTC,MAXLEVEL))


      WRITE (IUNMEM,'(A,T25,I15)')
     .       ' COMSPL ',(NCMSPL+MAXLEVEL*NPARTC)*8 +
     .                  (MCMSPL+MAXLEVEL*MPARTC)*4 + (KCMSPL+NLIMPS+1)*4

      WMINV  => RCMSPL(1)
      WMINS  => RCMSPL(2)
      WMINC  => RCMSPL(3)
      WMINL  => RCMSPL(4)
      SPLPAR => RCMSPL(5)
      RNUMB  => RCMSPL(6:5+N1ST+N2ND+N3RD+NLIM)
      PRMSPL => RCMSPL(6+N1ST+N2ND+N3RD+NLIM : NCMSPL)

      MAXLEV => ICMSPL(1)
      NLEVEL => ICMSPL(2)
      MAXRAD => ICMSPL(3)
      MAXPOL => ICMSPL(4)
      MAXTOR => ICMSPL(5)
      MAXADD => ICMSPL(6)

cdr formerly: maxlevel=15 was hard-coded, now: maxlevel=300 ?
      NODES  => ICMSPL(7:6+MAXLEVEL)
      NSSPL  => ICMSPL(7+MAXLEVEL:MCMSPL)

cdr
c  conditional expectation estimator, if trajectory "sees" additional surface ILIM
      ALLOCATE (NLPRCS(0:NLIMPS))
cdr conditional expectation estimator for A, M, I, PH
      NLPRCA => LCMSPL(1:NATM)
      NLPRCM => LCMSPL(1+NATM:NATM+NMOL)
      NLPRCI => LCMSPL(1+NATM+NMOL:NATM+NMOL+NION)
      NLPRCPH=> LCMSPL(1+NATM+NMOL+NION:NATM+NMOL+NION+NPHOT)
cdr  NLSPLT(ISURF): surface isurf is a "splitting-rr" surface
      NLSPLT => LCMSPL(1+NATM+NMOL+NION+NPHOT:KCMSPL)

      CALL EIRENE_INIT_COMSPL

      RETURN
      END SUBROUTINE EIRENE_ALLOC_COMSPL


      SUBROUTINE EIRENE_DEALLOC_COMSPL

      IF (.NOT.ALLOCATED(RCMSPL)) RETURN

      DEALLOCATE (RCMSPL)
      DEALLOCATE (ICMSPL)
      DEALLOCATE (LCMSPL)

      DEALLOCATE (RSPLST)
      DEALLOCATE (ISPLST)
      DEALLOCATE (NLPRCS)

      RETURN
      END SUBROUTINE EIRENE_DEALLOC_COMSPL


      SUBROUTINE EIRENE_INIT_COMSPL

      RCMSPL = 0._DP
      ICMSPL = 0
      LCMSPL = .FALSE.

      RSPLST = 0._DP
      ISPLST = 0
      NLPRCS = .FALSE.

      RETURN
      END SUBROUTINE EIRENE_INIT_COMSPL


      SUBROUTINE EIRENE_BROADCAST_COMSPL(ME)
      USE EIRMOD_MPI
      INTEGER, INTENT(IN) :: ME
      INTEGER :: IER

      IF (ME /= 0) CALL EIRENE_ALLOC_COMSPL

      CALL MPI_BCAST (RCMSPL,NCMSPL,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (RSPLST,MAXLEVEL*NPARTC,MPI_REAL8,0,
     .                MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ICMSPL,MCMSPL,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ISPLST,MAXLEVEL*MPARTC,MPI_INTEGER,0,
     .                MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (LCMSPL,KCMSPL,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NLPRCS,NLIMPS+1,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)

      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

      END SUBROUTINE EIRENE_BROADCAST_COMSPL

      
      END MODULE EIRMOD_COMSPL
