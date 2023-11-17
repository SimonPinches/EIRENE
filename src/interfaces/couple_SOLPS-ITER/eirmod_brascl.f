cdr Nov. 17:  prepare re-vive of short cycle. Comments.
cdr           compared to early short cycle version (1992):
cdr          --implicit correction for selected strata: NLSRON(ISTR)
cdr          In particular this allows to exclude vol.rec strata from implicit treatment
cdr             (old version: all strata or no stratum at all)
cdr          --Documentation of new version started
c
cdr          --seiod, seinw separated from seioda, seinwa. Have different meaning
cdr            (ion energy density in plasma flow, per species) !
cdr             --> missing rescaling of EI rates from atoms is now possible.
cdr            seioda(..npls) --> seioda(..natm), and seiod(..npls): additional new array
cdr          --Made alloc, dealloc, init: more symmetric code between od and nw arrays
cdr Dec. 18:  more consistent notation:
cdr           init_brascl2                   --> init_brascl
cdr           separated from dealloc_brascl: --> dealloc_rate_array
cdr           init_brascl1                   --> init_rate_array

      MODULE EIRMOD_BRASCL

C  RESCALING ARRAYS FOR NEUTRAL SOURCE TERMS IN SHORT LOOP ITERATION

c  allocate storage for short loop (semi-implicit) source term corrections
c  S..OD....: storage for source term from previous iteration.
c  Storage only for those strata which are in "short cycle loop"
c
c  Let:
c  Source terms S_i for plasma fluid equations, in iteration (TIME STEP) no. I
c
c  If the source term S_i from a stratum ISTR is derived from the previous one,
c  rather than re-calculated with new MC trajectories, then
c  (up to global rescaling factors):

c  S_i are written in a form: S_i = Sc_(i-1) + Sv_i(pl_i)
c    Sc: stays constant (frozen) between iterations i-1 and i.
c    Sv: re-evaluated from new plasma data at iteration i.
c
c...........................................................................
C   ALLOC_BRASCL,           CALLS: INIT_BRASCL            (FOR ..nw.. RATES)
C   ALLOC_RATE_ARRAY(ISTR), CALLS: INIT_RATE_ARRAY(ISTR)  (FOR ..od.. RATES)
C   DEALLOC_BRASCL                                        (FOR ..nw.. RATES)
C   DEALLOC_RATE_ARRAY (FOR ALL ISTR)                     (FOR ..od.. RATES)

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: EIRENE_ALLOC_BRASCL,     EIRENE_DEALLOC_BRASCL,
     P          EIRENE_ALLOC_RATE_ARRAY, EIRENE_DEALLOC_RATE_ARRAY,
     P          EIRENE_INIT_BRASCL,      EIRENE_INIT_RATE_ARRAY,
     P          RATE_STORE, RATES_ARRAY

      TYPE RATE_STORE
        REAL(DP), POINTER ::
     R      SEIOD(:,:),  ! needed only once, not per stratum!
     R      SPLODA(:,:,:), SEIODA(:,:), SEEODA(:,:), SMOODA(:,:),
     R      SPLODI(:,:,:), SEIODI(:,:), SEEODI(:,:), SMOODI(:,:),
     R      SPLODM(:,:,:), SEIODM(:,:), SEEODM(:,:), SMOODM(:,:)
      END TYPE RATE_STORE

      TYPE RATES_ARRAY
        TYPE(RATE_STORE), POINTER :: RTA
      END TYPE RATES_ARRAY

      TYPE(RATES_ARRAY), PUBLIC, ALLOCATABLE, SAVE :: RTS(:)
      INTEGER, PUBLIC, ALLOCATABLE, SAVE :: ITS(:), ITS_COUNT(:)


      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R SEINW(:,:),  ! needed only once, not per stratum!
     R SPLNWA(:,:,:), SEINWA(:,:), SEENWA(:,:), SMONWA(:,:),
     R SPLNWI(:,:,:), SEINWI(:,:), SEENWI(:,:), SMONWI(:,:),
     R SPLNWM(:,:,:), SEINWM(:,:), SEENWM(:,:), SMONWM(:,:)


      CONTAINS


      SUBROUTINE EIRENE_ALLOC_BRASCL

      INTEGER :: ISTR

      IF (ALLOCATED(SPLNWA)) RETURN

      ALLOCATE (RTS(NSTRA))
      ALLOCATE (ITS(NSTRA))
      ALLOCATE (ITS_COUNT(NSTRA))

C  ENERGY DENSITY IN IPLS PLASMA FLOW
      ALLOCATE (SEINW(NRAD,NPLS))
C  SOURCE TERM CORRECTIONS, CURRENT (NEW) RATES: PART,EE,EI,MOM
      ALLOCATE (SPLNWA(NRAD,NATM,NPLS))
      ALLOCATE (SEINWA(NRAD,NATM))
      ALLOCATE (SEENWA(NRAD,NATM))
      ALLOCATE (SMONWA(NRAD,NATM))
      ALLOCATE (SPLNWI(NRAD,NION,NPLS))
      ALLOCATE (SEINWI(NRAD,NION))
      ALLOCATE (SEENWI(NRAD,NION))
      ALLOCATE (SMONWI(NRAD,NION))
      ALLOCATE (SPLNWM(NRAD,NMOL,NPLS))
      ALLOCATE (SEINWM(NRAD,NMOL))
      ALLOCATE (SEENWM(NRAD,NMOL))
      ALLOCATE (SMONWM(NRAD,NMOL))

      WRITE (IUNMEM,'(A,T25,I15)')
     .      ' BRASCL ',NRAD*(NATM*(NPLS+3)+NION*(NPLS+3)+
     .                       NMOL*(NPLS+3)+NPLS)*8

      DO ISTR=1, NSTRA
        NULLIFY(RTS(ISTR)%RTA)
      END DO

      ITS = 1
      ITS_COUNT = 0

      CALL EIRENE_INIT_BRASCL

      RETURN
      END SUBROUTINE EIRENE_ALLOC_BRASCL


      SUBROUTINE EIRENE_ALLOC_RATE_ARRAY(ISTR)
C  allocate storage for semi-implicit correction, for stratum ISTR
C  This is done for all strata, for which semi-implicit corrections (short-cycling)
C  is to be done.

      INTEGER, INTENT(IN) :: ISTR

      IF (ASSOCIATED(RTS(ISTR)%RTA)) RETURN

      ALLOCATE (RTS(ISTR)%RTA)

C  ENERGY DENSITY IN IPLS PLASMA FLOW
      ALLOCATE (RTS(ISTR)%RTA%SEIOD(NRAD,NPLS))

C  SOURCE TERM CORRECTIONS, PREVIOUS (OLD) RATES: PART,EE,EI,MOM
      ALLOCATE (RTS(ISTR)%RTA%SPLODA(NRAD,NATM,NPLS))
      ALLOCATE (RTS(ISTR)%RTA%SEIODA(NRAD,NATM))
      ALLOCATE (RTS(ISTR)%RTA%SEEODA(NRAD,NATM))
      ALLOCATE (RTS(ISTR)%RTA%SMOODA(NRAD,NATM))
C
      ALLOCATE (RTS(ISTR)%RTA%SPLODI(NRAD,NION,NPLS))
      ALLOCATE (RTS(ISTR)%RTA%SEIODI(NRAD,NION))
      ALLOCATE (RTS(ISTR)%RTA%SEEODI(NRAD,NION))
      ALLOCATE (RTS(ISTR)%RTA%SMOODI(NRAD,NION))
C
      ALLOCATE (RTS(ISTR)%RTA%SPLODM(NRAD,NMOL,NPLS))
      ALLOCATE (RTS(ISTR)%RTA%SEIODM(NRAD,NMOL))
      ALLOCATE (RTS(ISTR)%RTA%SEEODM(NRAD,NMOL))
      ALLOCATE (RTS(ISTR)%RTA%SMOODM(NRAD,NMOL))

      WRITE (IUNMEM,'(A,T25,I15)')
     .      ' BRASCL ',NRAD*(NATM*(NPLS+3)+NION*(NPLS+3)+
     .                       NMOL*(NPLS+3)+NPLS)*8

      CALL EIRENE_INIT_RATE_ARRAY(ISTR)

      RETURN
      END SUBROUTINE EIRENE_ALLOC_RATE_ARRAY


      SUBROUTINE EIRENE_DEALLOC_RATE_ARRAY
cdr  dec.20
cdr  dealloc_rate_array must come before dealloc_brascl,
cdr  because the latter also deallocates RTS....
      INTEGER :: ISTR

      DO ISTR = 1, NSTRA
        IF (.NOT.ASSOCIATED(RTS(ISTR)%RTA)) CYCLE
        DEALLOCATE (RTS(ISTR)%RTA%SEIOD)
C
        DEALLOCATE (RTS(ISTR)%RTA%SPLODA)
        DEALLOCATE (RTS(ISTR)%RTA%SEIODA)
        DEALLOCATE (RTS(ISTR)%RTA%SEEODA)
        DEALLOCATE (RTS(ISTR)%RTA%SMOODA)

        DEALLOCATE (RTS(ISTR)%RTA%SPLODI)
        DEALLOCATE (RTS(ISTR)%RTA%SEIODI)
        DEALLOCATE (RTS(ISTR)%RTA%SEEODI)
        DEALLOCATE (RTS(ISTR)%RTA%SMOODI)

        DEALLOCATE (RTS(ISTR)%RTA%SPLODM)
        DEALLOCATE (RTS(ISTR)%RTA%SEIODM)
        DEALLOCATE (RTS(ISTR)%RTA%SEEODM)
        DEALLOCATE (RTS(ISTR)%RTA%SMOODM)
      END DO

      RETURN
      END SUBROUTINE EIRENE_DEALLOC_RATE_ARRAY

      SUBROUTINE EIRENE_DEALLOC_BRASCL


      IF (.NOT.ALLOCATED(SPLNWA)) RETURN

      DEALLOCATE (RTS)
      DEALLOCATE (ITS)
      DEALLOCATE (ITS_COUNT)

      DEALLOCATE (SEINW)
C
      DEALLOCATE (SPLNWA)
      DEALLOCATE (SEINWA)
      DEALLOCATE (SEENWA)
      DEALLOCATE (SMONWA)
      DEALLOCATE (SPLNWI)
      DEALLOCATE (SEINWI)
      DEALLOCATE (SEENWI)
      DEALLOCATE (SMONWI)
      DEALLOCATE (SPLNWM)
      DEALLOCATE (SEINWM)
      DEALLOCATE (SEENWM)
      DEALLOCATE (SMONWM)

      RETURN
      END SUBROUTINE EIRENE_DEALLOC_BRASCL


      SUBROUTINE EIRENE_INIT_RATE_ARRAY(ISTR)

      INTEGER, INTENT(IN) :: ISTR

      RTS(ISTR)%RTA%SEIOD = 0.D0
C
      RTS(ISTR)%RTA%SPLODA = 0.D0
      RTS(ISTR)%RTA%SEIODA = 0.D0
      RTS(ISTR)%RTA%SEEODA = 0.D0
      RTS(ISTR)%RTA%SMOODA = 0.D0
      RTS(ISTR)%RTA%SPLODI = 0.D0
      RTS(ISTR)%RTA%SEIODI = 0.D0
      RTS(ISTR)%RTA%SEEODI = 0.D0
      RTS(ISTR)%RTA%SMOODI = 0.D0
      RTS(ISTR)%RTA%SPLODM = 0.D0
      RTS(ISTR)%RTA%SEIODM = 0.D0
      RTS(ISTR)%RTA%SEEODM = 0.D0
      RTS(ISTR)%RTA%SMOODM = 0.D0

      RETURN
      END SUBROUTINE EIRENE_INIT_RATE_ARRAY


      SUBROUTINE EIRENE_INIT_BRASCL

      SEINW = 0.D0
C
      SPLNWA = 0.D0
      SEINWA = 0.D0
      SEENWA = 0.D0
      SMONWA = 0.D0
      SPLNWI = 0.D0
      SEINWI = 0.D0
      SEENWI = 0.D0
      SMONWI = 0.D0
      SPLNWM = 0.D0
      SEINWM = 0.D0
      SEENWM = 0.D0
      SMONWM = 0.D0

      RETURN
      END SUBROUTINE EIRENE_INIT_BRASCL

      END MODULE EIRMOD_BRASCL
