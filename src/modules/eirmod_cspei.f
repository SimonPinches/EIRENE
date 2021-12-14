      MODULE EIRMOD_CSPEI
cdr Handle storage for:
c   cspei: Standard deviations, co.-variances, etc.
c   plasma_bckgrnd: For interfacing with external plasma code
c                   also used for internal non-lin. iterations within eirene
c
c   Size of plasma_bckgrnd:  currently NIINTF (formerly: NIDC)
c   to be checked: is this consistent with usage of plasma_bckgrnd?
c   also in modbgk? 

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: EIRENE_ALLOC_CSPEI, EIRENE_DEALLOC_CSPEI,
     P          EIRENE_INIT_CSPEI,
     P          EIRENE_ALLOC_BCKGRND, EIRENE_DEALLOC_BCKGRND,
     P          EIRENE_INIT_BCKGRND

      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     R  SMESTV(:,:), SMESTS(:,:)

      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R  STV(:,:), STVW(:,:), STVC(:,:,:),
     R  STVS(:),  STVWS(:),  STVCS(:,:)

      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R  EE(:,:),    FF(:,:),
     R  EES(:),     FFS(:),
     R  SDVIA(:,:), SDVIAW(:,:), SDVIAC(:,:,:)

      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     R  PLASMA_BCKGRND(:,:)

      REAL(DP), PUBLIC, POINTER, SAVE ::
     .  TEINTF(:),   TIINTF(:,:), DIINTF(:,:),
     .  VXINTF(:,:), VYINTF(:,:), VZINTF(:,:),
     .  BXINTF(:),   BYINTF(:),   BZINTF(:),   BFINTF(:),
     .  VLINTF(:),   ADINTF(:,:)

      INTEGER, PUBLIC, SAVE :: IESTR

      INTEGER, PUBLIC, SAVE ::
     I  NIDC, NIDV, NIDS


      CONTAINS


      SUBROUTINE EIRENE_ALLOC_CSPEI
cdr
c  called from main routine eirene.f, allocates storage for
c  storage arrays which are needed for "sum over strata"
c  smestv, smests, smestl
c  and for intermediate storage arrays for variance per history evaluation
c  stv,sdvia,....
cdr
      INTEGER, PARAMETER :: IL = SELECTED_INT_KIND(15)
      INTEGER(IL) :: MEM

      IF (ALLOCATED(SMESTV)) RETURN  ! allocated smestv is used as indicator for: 'all fields are allocated'

      IF (NSMSTRA > 0) THEN
        NIDV=NVOLTL
        NIDS=NSRFTL
      ELSE
        NIDV=1
        NIDS=1
      END IF
CDR  storage for for sum over strata, surface and volume tallies.
      ALLOCATE (SMESTV(NIDV,NRTAL))
      ALLOCATE (SMESTS(NIDS,NLMPGS))
CDR   same for spectra, but:
cdr   ALLOCATE (SMESTL(NADSPC))   ! TO BE MOVED HERE FROM INPUT.F,  NOT POSSIBLE BECAUSE DIFFERENT DATA TYPE FOR SMESTL


cdr  these next arrays are intermediate storage array to perform variance per history calculations.
      ALLOCATE (STV(NSD,NRTAL))
      ALLOCATE (STVW(NSDW,NLIMPS))
      ALLOCATE (STVC(0:2,NCV,NRTAL))
      ALLOCATE (STVS(NSD))
      ALLOCATE (STVWS(NSDW))
      ALLOCATE (STVCS(0:2,NCV))
      ALLOCATE (EE(NSD,NRTAL))
      ALLOCATE (FF(NSDW,NLIMPS))
      ALLOCATE (EES(NSD))
      ALLOCATE (FFS(NSDW))
      ALLOCATE (SDVIA(NSD,NRTAL))
      ALLOCATE (SDVIAW(NSDW,NLIMPS))
      ALLOCATE (SDVIAC(2,NCV,NRTAL))

CDR  same again: intermediate storage for spectra, data type prevents this from having it here?
CDR BEGIN:
CDR TO BE MOVED HERE FROM INPUT.F, NOT EASILY POSSIBLE,
CDR BECAUSE, SADLY, DIFFERENT DATA TYPE FOR SSPEC
c     ALLOCATE(SSPEC)
c     ALLOCATE(SSPEC%SPC(0:NSPS+1))
c  standard deviation of spectra tallies, sum over strata intermediate storage
!     IF (NSIGI_SPC > 0) THEN
c       ALLOCATE(SSPEC%SDV(0:NSPS+1))
c       ALLOCATE(SSPEC%SGM(0:NSPS+1))
c       ALLOCATE(SSPEC%STV(0:NSPS+1))
c       ALLOCATE(SSPEC%GG(0:NSPS+1))
!     END IF
c     SSPEC = ESPEC   !  ????
c     SMESTL(J)%PSPC => SSPEC   !  ???? not sure here !!!!
CDR END

C  TOTAL ALLOCATED STORAGE IN THIS ROUTINE

      MEM = (NIDV*NRTAL+(3*NSD+5*NCV)*NRTAL +
     .                  NIDS*NLMPGS + 3*NSDW*NLIMPS + 2*NSD +
     .                  2*NSDW + 3*NCV)*8

      WRITE (IUNMEM,'(A,T25,I15)')
     .       ' CSPEI ', MEM
      
      CALL EIRENE_INIT_CSPEI

      RETURN
      END SUBROUTINE EIRENE_ALLOC_CSPEI


      SUBROUTINE EIRENE_DEALLOC_CSPEI

      IF (.NOT.ALLOCATED(SMESTV)) RETURN

      DEALLOCATE (SMESTV)
      DEALLOCATE (SMESTS)

      DEALLOCATE (STV)
      DEALLOCATE (STVW)
      DEALLOCATE (STVC)
      DEALLOCATE (STVS)
      DEALLOCATE (STVWS)
      DEALLOCATE (STVCS)
      DEALLOCATE (EE)
      DEALLOCATE (FF)
      DEALLOCATE (EES)
      DEALLOCATE (FFS)
      DEALLOCATE (SDVIA)
      DEALLOCATE (SDVIAW)
      DEALLOCATE (SDVIAC)

      RETURN
      END SUBROUTINE EIRENE_DEALLOC_CSPEI


      SUBROUTINE EIRENE_INIT_CSPEI

      SMESTV = 0._DP
      SMESTS = 0._DP

      STV    = 0._DP
      STVW   = 0._DP
      STVC   = 0._DP
      STVS   = 0._DP
      STVWS  = 0._DP
      STVCS  = 0._DP
      EE     = 0._DP
      FF     = 0._DP
      EES    = 0._DP
      FFS    = 0._DP
      SDVIA  = 0._DP
      SDVIAW = 0._DP
      SDVIAC = 0._DP

      RETURN
      END SUBROUTINE EIRENE_INIT_CSPEI


      SUBROUTINE EIRENE_ALLOC_BCKGRND
cdr  Sept 19:  NIDC renamed to NIINTF
cdr  some time in past: typo in NIDC fixed: NIDC=...+3*NPLSV
cdr                            rather than: NIDC=...+4*NPLSV
cdr  Jan 20: further comments.

cdr  Allocate storage for handling background medium,
cdr  originally only for transfer of plasma data from external code
cdr  and PROFR  profile options INDPRO=6 or INDPRO=7


!pb   NIDC=1*NPLS+NAIN+6+NPLSTI+4*NPLSV
      NIDC=1*NPLS+NAIN+6+NPLSTI+3*NPLSV

      IF (.NOT.ALLOCATED(PLASMA_BCKGRND)) THEN

        ALLOCATE(PLASMA_BCKGRND(NIDC,NRAD))

        TEINTF => PLASMA_BCKGRND(1+0+0*NPLS             ,   :)
        TIINTF => PLASMA_BCKGRND(1+1+0*NPLS        :
     .                           1+0+0*NPLS+NPLSTI,         :)
        DIINTF => PLASMA_BCKGRND(1+1+0*NPLS+NPLSTI :
     .                           1+0+1*NPLS+NPLSTI,         :)
        VXINTF => PLASMA_BCKGRND(1+1+1*NPLS+NPLSTI+0*NPLSV :
     .                           1+0+1*NPLS+NPLSTI+1*NPLSV, :)
        VYINTF => PLASMA_BCKGRND(1+1+1*NPLS+NPLSTI+1*NPLSV :
     .                           1+0+1*NPLS+NPLSTI+2*NPLSV, :)
        VZINTF => PLASMA_BCKGRND(1+1+1*NPLS+NPLSTI+2*NPLSV :
     .                           1+0+1*NPLS+NPLSTI+3*NPLSV, :)
        BXINTF => PLASMA_BCKGRND(1+1+1*NPLS+NPLSTI+3*NPLSV, :)
        BYINTF => PLASMA_BCKGRND(1+2+1*NPLS+NPLSTI+3*NPLSV, :)
        BZINTF => PLASMA_BCKGRND(1+3+1*NPLS+NPLSTI+3*NPLSV, :)
        BFINTF => PLASMA_BCKGRND(1+4+1*NPLS+NPLSTI+3*NPLSV, :)
        VLINTF => PLASMA_BCKGRND(1+5+1*NPLS+NPLSTI+3*NPLSV, :)
!pb     ADINTF => PLASMA_BCKGRND(1+6+1*NPLS+NPLSTI+3*NPLSV :
!pb  .                             6+1*NPLS+NPLSTI+4*NPLSV+NAIN, :)
        ADINTF => PLASMA_BCKGRND(1+6+1*NPLS+NPLSTI+3*NPLSV :
     .                             6+1*NPLS+NPLSTI+3*NPLSV+NAIN, :)

        CALL EIRENE_INIT_BCKGRND

      END IF

      RETURN
      END SUBROUTINE EIRENE_ALLOC_BCKGRND


      SUBROUTINE EIRENE_DEALLOC_BCKGRND

      IF (ALLOCATED(PLASMA_BCKGRND)) DEALLOCATE(PLASMA_BCKGRND)

      NULLIFY(TEINTF)
      NULLIFY(TIINTF)
      NULLIFY(DIINTF)
      NULLIFY(VXINTF)
      NULLIFY(VYINTF)
      NULLIFY(VZINTF)
      NULLIFY(BXINTF)
      NULLIFY(BYINTF)
      NULLIFY(BZINTF)
      NULLIFY(BFINTF)
      NULLIFY(VLINTF)
      NULLIFY(ADINTF)

      RETURN
      END SUBROUTINE EIRENE_DEALLOC_BCKGRND


      SUBROUTINE EIRENE_INIT_BCKGRND

      PLASMA_BCKGRND = 0._DP

      RETURN
      END SUBROUTINE EIRENE_INIT_BCKGRND

      END MODULE EIRMOD_CSPEI
