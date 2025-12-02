cdr  Dec. 15: added allocatable public array fnuiar(npls): species-dependent
cdr                                 collision frequency, for FP collisions

      MODULE EIRMOD_CFPLK
C   parameters for fokker-planck collision operator,
C   local, along the trajectory
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: EIRENE_ALLOC_CFPLK, EIRENE_DEALLOC_CFPLK

cdr  begin threadprivate here
      REAL(DP), PUBLIC, SAVE ::
     R E0PAR, VELPAR, VELPER, VLXPAR, VLYPAR, VLZPAR, SIGPAR,
     R BVEC(3),  !dr local B field vector, cartesian
cdr B field, full transf. tensor (dyadic)
cdr cart. B field components wrt. par,(0) perp,(1) and dia,(2) components, resp.
     R BX(0:2), BY(0:2), BZ(0:2),
     R TAUE
      LOGICAL, public, save :: LCART,COLFLAG

      REAL(DP), PUBLIC, ALLOCATABLE, SAVE :: FNUIAR(:)

#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP  THREADPRIVATE(E0PAR,VELPAR,VELPER,VLXPAR,VLYPAR,VLZPAR,SIGPAR,
!$OMP& BVEC,BX,BY,BZ,TAUE,LCART,COLFLAG,FNUIAR)
#endif

cdr  end threadprivate here

      CONTAINS


      SUBROUTINE EIRENE_ALLOC_CFPLK

      IF (ALLOCATED(FNUIAR)) RETURN

      ALLOCATE (FNUIAR(NPLS))

      CALL EIRENE_INIT_CFPLK

      RETURN
      END SUBROUTINE EIRENE_ALLOC_CFPLK


      SUBROUTINE EIRENE_DEALLOC_CFPLK

      IF (.NOT. ALLOCATED(FNUIAR)) RETURN

      DEALLOCATE(FNUIAR)

      RETURN
      END SUBROUTINE EIRENE_DEALLOC_CFPLK


      SUBROUTINE EIRENE_INIT_CFPLK

      FNUIAR = 0._DP

      RETURN
      END SUBROUTINE EIRENE_INIT_CFPLK


      END MODULE EIRMOD_CFPLK
