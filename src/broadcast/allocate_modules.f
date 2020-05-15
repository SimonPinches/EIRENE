c     memory allocation for additional processors
       
!pb  15.05.20  allocation is now done in the broadcast routines
!pb            directly in the modules
!pb            allocations triggered here belong to modules with no
!pb            broadcasting routine     
      
!pb  03.06.09  NCHORI --> NCHOR in if-condition as NCHORI is not yet broadcasted

      SUBROUTINE EIRENE_ALLOCATE_MODULES
      USE EIRMOD_PARMMOD
      USE EIRMOD_CESTIM
      USE EIRMOD_COMUSR
      USE EIRMOD_CADGEO
      USE EIRMOD_CAI
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CINIT
      USE EIRMOD_COMSIG
      USE EIRMOD_CPOLYG
      USE EIRMOD_CGRID
      USE EIRMOD_CZT1
      USE EIRMOD_CTRCEI
      USE EIRMOD_CGEOM
      USE EIRMOD_CSDVI
      
      USE EIRMOD_CTETRA
      USE EIRMOD_COMPRT
      USE EIRMOD_CPES
      USE EIRMOD_COMNNL
      USE EIRMOD_COMSOU
      USE EIRMOD_CSTEP
      USE EIRMOD_COMSPL
      USE EIRMOD_CTEXT
      USE EIRMOD_CLGIN
      USE EIRMOD_COUTAU
      USE EIRMOD_COMXS
      USE EIRMOD_CTRIG
      USE EIRMOD_CLAST
      USE EIRMOD_CPLOT
      USE EIRMOD_CREF
      USE EIRMOD_CSPEI
      USE EIRMOD_CSPEZ
      USE EIRMOD_CUPD
      USE EIRMOD_CFPLK
      IMPLICIT NONE
      INTEGER JTRJ

      CALL EIRENE_ALLOC_COMPRT(NPRS)
      CALL EIRENE_ALLOC_CPES
      CALL EIRENE_ALLOC_CLAST
      CALL EIRENE_ALLOC_CPLOT
      CALL EIRENE_ALLOC_CSPEI
      CALL EIRENE_ALLOC_CSPEZ
      CALL EIRENE_ALLOC_CUPD
      CALL EIRENE_ALLOC_CFPLK

!  allocate and initialize storage for trajectories
      IF (.NOT.ALLOCATED(TRAJ)) THEN
        ALLOCATE (TRAJ(NCHOR+NTRJ))

        DO JTRJ = 1, NCHOR+NTRJ
          ALLOCATE(TRAJ(JTRJ)%TRJ)
          TRAJ(JTRJ)%TRJ%NCOU_CELL = 0
          NULLIFY(TRAJ(JTRJ)%TRJ%CELLS)
        END DO
      END IF

      CALL EIRENE_ALLOC_COUPLE

      RETURN
      END SUBROUTINE EIRENE_ALLOCATE_MODULES
