      MODULE EIRMOD_CSPEI
cdr Handle storage for:
c   CSPEI: Standard deviations, covariances, etc.
cdr and, quite unrelated, for:
c   PLASMA_BCKGRND: For interfacing with external plasma code
c                   also used for internal non-lin. iterations within eirene
c
c   Size of plasma_bckgrnd: currently NIINTF (formerly: NIDC)
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
     I  NIINTF, NIDV, NIDS


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

! allocated smestv is used as indicator for: 'all fields are allocated'
      IF (ALLOCATED(SMESTV)) RETURN

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
CDR  same for spectra, but:
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

cdr  called from:
c    input.f  (if any indpro(1:12) =6 or =7). Called TWICE
c              once after reading block 13 and before reading block 14,
c              once again after statement 4000.
c    eirsrt.f (coupling to various B2 code variants)
c    rplam.f:
c       rplam_shrt.f
c       rplam_long.f (if IFLG=10)
c    modbgk.f (transfer of modified virtual background for next iteration)
cdr

c      provide storage for
c      transfer of background (plasma) tallies
c      from an external code into eirene background tallies
c
cdr  In iterative mode, e.g. BGK iterations,
c      ALLOC_BCKGRND may also be called.
c      Then: risk of a hidden link? And unnecessary stuff to be dealt
cdr    within MODBGK (such as B field ?)
cdr    And other stuff from MODBGK may be missing:  tabel3, lgvac,....
c
c    currently: no E field information ?

cdr  size of interfacing plasma data storage
      NIINTF=6+1*NPLS+NPLSTI+3*NPLSV+NAIN

      IF (.NOT.ALLOCATED(PLASMA_BCKGRND)) THEN

        ALLOCATE(PLASMA_BCKGRND(NIINTF,NRAD))

cdr initial : final storage for ADINTF tally corrected
c  1
        TEINTF => PLASMA_BCKGRND(1+0+0*NPLS             ,   :)
c  2 -- 1+nplsti
        TIINTF => PLASMA_BCKGRND(1+1+0*NPLS        :
     .                           1+0+0*NPLS+NPLSTI,         :)
c  2+nplsti -- 1+nplsti+npls
        DIINTF => PLASMA_BCKGRND(1+1+0*NPLS+NPLSTI :
     .                           1+0+1*NPLS+NPLSTI,         :)
c  2+nplsti+npls -- 1+nplsti+npls+nplsv
        VXINTF => PLASMA_BCKGRND(1+1+1*NPLS+NPLSTI+0*NPLSV :
     .                           1+0+1*NPLS+NPLSTI+1*NPLSV, :)
        VYINTF => PLASMA_BCKGRND(1+1+1*NPLS+NPLSTI+1*NPLSV :
     .                           1+0+1*NPLS+NPLSTI+2*NPLSV, :)
        VZINTF => PLASMA_BCKGRND(1+1+1*NPLS+NPLSTI+2*NPLSV :
     .                           1+0+1*NPLS+NPLSTI+3*NPLSV, :)
c
        BXINTF => PLASMA_BCKGRND(1+1+1*NPLS+NPLSTI+3*NPLSV, :)
        BYINTF => PLASMA_BCKGRND(1+2+1*NPLS+NPLSTI+3*NPLSV, :)
        BZINTF => PLASMA_BCKGRND(1+3+1*NPLS+NPLSTI+3*NPLSV, :)
        BFINTF => PLASMA_BCKGRND(1+4+1*NPLS+NPLSTI+3*NPLSV, :)
        VLINTF => PLASMA_BCKGRND(1+5+1*NPLS+NPLSTI+3*NPLSV, :)
c  7+npls+nplsti+3nplsv -- 6+npls+nplsti+3nplsv+nain
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

      IF (ALLOCATED(PLASMA_BCKGRND)) PLASMA_BCKGRND = 0._DP

      RETURN
      END SUBROUTINE EIRENE_INIT_BCKGRND

      END MODULE EIRMOD_CSPEI
