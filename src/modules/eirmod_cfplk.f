cdr  Dec. 15: added allocatable public array fnuiar(npls): species dependent 
cdr                                 collision frequency, for FP collisions

      MODULE EIRMOD_CFPLK
C   parameters for fokker planck collision operator
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
 
      IMPLICIT NONE
 
      PRIVATE
 
      PUBLIC :: EIRENE_ALLOC_CFPLK, EIRENE_DEALLOC_CFPLK
 
      REAL(DP), PUBLIC, SAVE ::
     R E0PAR, VELPAR, VELPER, VLXPAR, VLYPAR, VLZPAR, SIGPAR,
     R BVEC(3), BBX, BBY, BBZ,
     R TAUE
      LOGICAL, public, save :: LCART,COLFLAG

      REAL(DP), PUBLIC, ALLOCATABLE, SAVE :: FNUIAR(:)
 
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
