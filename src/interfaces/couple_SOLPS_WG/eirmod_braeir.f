      MODULE EIRMOD_BRAEIR

C  PLASMA DATA: NI,TE,TI,VV,UU,PR,UP,RR,FNIX,FNIY.. (BRAAMS ---> EIRENE)

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      IMPLICIT NONE

      PRIVATE

      PUBLIC :: EIRENE_ALLOC_BRAEIR, EIRENE_DEALLOC_BRAEIR,
     P          EIRENE_INIT_BRAEIR, EIRENE_ALLOC_TARGET_DATA,
     P          targetdata

      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     &  dnib(:,:), teb(:), tib(:), tnb(:), pob(:), prb(:), neb(:),
     &  vvb(:,:), uub(:,:), wwb(:,:), upb(:,:),
     &  rrb(:), volb(:), bfeldb(:),
     &  bbxb(:), bbyb(:), bbzb(:),
     &  exb(:), eyb(:), ezb(:),
     &  fnib(:,:), feib(:), feeb(:),
     &  vxb(:,:), vyb(:,:), vzb(:,:),
     &  uudiab(:,:), vvdiab(:,:),
     &  delta_sheathb(:)

      real(DP), public, allocatable, save ::
     &  amb(:), znb(:), zab(:),
     &  zib(:,:)                    !nh spatially dependent charge

      REAL(DP), ALLOCATABLE, SAVE, PUBLIC ::
     &  PUX(:),  PUY(:),  PVX(:),  PVY(:)

      REAL(DP), PUBLIC, SAVE :: L_MACRO_B25 !nh macroscopic length scale

      type :: targetdata
        integer :: ntrgdat
        integer, allocatable :: faces(:)
        real(DP), allocatable :: fcori(:), flength(:)
        real(DP), allocatable :: 
     &    tet(:), tit(:), tnt(:), fe(:), fsh(:),
     &    flux(:,:), dnit(:,:), vxt(:,:), vyt(:,:), vzt(:,:),
     &    fi(:,:), fel(:,:), vpart(:,:), mach(:,:), usr(:,:),
     &    zit(:,:)  !nh charge
      end type targetdata

      type(targetdata), public, allocatable, save :: trgt(:)

      CONTAINS

      SUBROUTINE EIRENE_ALLOC_BRAEIR(nCv, nFc, NFLD)

      INTEGER, INTENT(IN) :: nCv, nFc, NFLD

      IF (ALLOCATED(DNIB)) RETURN

      NFL = NFLD
 
      ALLOCATE (NEB(NCV))
      ALLOCATE (TEB(NCV))
      ALLOCATE (TIB(NCV))
      ALLOCATE (TNB(NCV))
      ALLOCATE (POB(NCV))
      ALLOCATE (PRB(NCV))
      ALLOCATE (RRB(NCV))
      ALLOCATE (VOLB(NCV))
      ALLOCATE (BFELDB(NCV))
      ALLOCATE (BBXB(NCV))
      ALLOCATE (BBYB(NCV))
      ALLOCATE (BBZB(NCV))
      ALLOCATE (EXB(NCV))
      ALLOCATE (EYB(NCV))
      ALLOCATE (EZB(NCV))

      ALLOCATE (DNIB(NCV,NFL))
      ALLOCATE (VVB(NCV,NFL))
      ALLOCATE (UUB(NCV,NFL))
      ALLOCATE (WWB(NCV,NFL))
      ALLOCATE (UPB(NCV,NFL))
      ALLOCATE (UUDIAB(NCV,NFL))
      ALLOCATE (VVDIAB(NCV,NFL))

      ALLOCATE (FEIB(NFC))
      ALLOCATE (FEEB(NFC))
      ALLOCATE (FNIB(NFC,NFL))
      ALLOCATE (VXB(NFC,NFL))
      ALLOCATE (VYB(NFC,NFL))
      ALLOCATE (VZB(NFC,NFL))
      ALLOCATE (DELTA_SHEATHB(NFC))

      allocate(amb(nfl))
      allocate(znb(nfl))
      allocate(zab(nfl))
      allocate(zib(nCv,nfl))
 
      ALLOCATE (PUX(NCV))
      ALLOCATE (PUY(NCV))
      ALLOCATE (PVX(NCV))
      ALLOCATE (PVY(NCV))
 
      WRITE (IUNMEM,'(A,T25,I15)')
     .      ' BRAEIR ',NCV*(7*NFL+11)*8 + NFC*(4*NFL+2)*8 +
     .                 3*NFL*8

      CALL EIRENE_INIT_BRAEIR

      RETURN
      END SUBROUTINE EIRENE_ALLOC_BRAEIR

      subroutine eirene_alloc_target_data(ist,ndat,nfl)
      implicit none
      integer, intent(in) :: ist, ndat, nfl
      integer, save :: memr=0, memi = 0

      if (ist <= ubound(trgt,1)) then
        trgt(ist)%ntrgdat = ndat
        if (.not.allocated(trgt(ist)%faces)) then
          allocate (trgt(ist)%faces(ndat))
          allocate (trgt(ist)%fcori(ndat))
          allocate (trgt(ist)%flength(ndat))
          allocate (trgt(ist)%tet(ndat))
          allocate (trgt(ist)%tit(ndat))
          allocate (trgt(ist)%tnt(ndat))
          allocate (trgt(ist)%fe(ndat))
          allocate (trgt(ist)%fsh(ndat))
          allocate (trgt(ist)%flux(ndat,nfl))
          allocate (trgt(ist)%dnit(ndat,nfl))
          allocate (trgt(ist)%vxt(ndat,nfl))
          allocate (trgt(ist)%vyt(ndat,nfl))
          allocate (trgt(ist)%vzt(ndat,nfl))
          allocate (trgt(ist)%fi(ndat,nfl))
          allocate (trgt(ist)%fel(ndat,nfl))
          allocate (trgt(ist)%vpart(ndat,nfl))
          allocate (trgt(ist)%mach(ndat,nfl))
          allocate (trgt(ist)%usr(ndat,nfl))
          allocate (trgt(ist)%zit(ndat,nfl))
          memr = memr + ndat*(4+11*nfl)
          memi = memi + ndat + 1
        end if
        
        trgt(ist)%faces = 0
        trgt(ist)%fcori = 0.0_DP
        trgt(ist)%flength = 0.0_DP
        trgt(ist)%tet = 0.0_DP
        trgt(ist)%tit = 0.0_DP
        trgt(ist)%tnt = 0.0_DP
        trgt(ist)%fe = 0.0_DP
        trgt(ist)%fsh = 0.0_DP
        trgt(ist)%flux = 0.0_DP
        trgt(ist)%dnit = 0.0_DP
        trgt(ist)%vxt = 0.0_DP
        trgt(ist)%vyt = 0.0_DP
        trgt(ist)%vzt = 0.0_DP
        trgt(ist)%fi = 0.0_DP
        trgt(ist)%fel = 0.0_DP
        trgt(ist)%vpart = 0.0_DP
        trgt(ist)%mach = 0.0_DP
        trgt(ist)%usr = 0.0_DP
        trgt(ist)%zit = 0.0_DP
        
      end if

      if (ist == ubound(trgt,1))
     .  WRITE (IUNMEM,'(A,T25,I15)') 
     .  ' BRAEIR TARGET DATA ',MEMR*8 + MEMI*4

      return
      end subroutine eirene_alloc_target_data


      SUBROUTINE EIRENE_DEALLOC_BRAEIR
      integer :: i

      IF (.NOT.ALLOCATED(DNIB)) RETURN

      DEALLOCATE (NEB)
      DEALLOCATE (TEB)
      DEALLOCATE (TIB)
      DEALLOCATE (TNB)
      DEALLOCATE (POB)
      DEALLOCATE (PRB)
      DEALLOCATE (RRB)
      DEALLOCATE (VOLB)
      DEALLOCATE (BFELDB)
      DEALLOCATE (BBXB)
      DEALLOCATE (BBYB)
      DEALLOCATE (BBZB)
      DEALLOCATE (EXB)
      DEALLOCATE (EYB)
      DEALLOCATE (EZB)

      DEALLOCATE (DNIB)
      DEALLOCATE (VVB)
      DEALLOCATE (UUB)
      DEALLOCATE (WWB)
      DEALLOCATE (UPB)
      DEALLOCATE (UUDIAB)
      DEALLOCATE (VVDIAB)
      DEALLOCATE (DELTA_SHEATHB)

      DEALLOCATE (FEIB)
      DEALLOCATE (FEEB)
      DEALLOCATE (FNIB)
      DEALLOCATE (VXB)
      DEALLOCATE (VYB)
      DEALLOCATE (VZB)

      deallocate (amb)
      deallocate (znb)
      deallocate (zab)
      deallocate (zib)

      DEALLOCATE (PUX)
      DEALLOCATE (PUY)
      DEALLOCATE (PVX)
      DEALLOCATE (PVY)

      IF (ALLOCATED(TRGT)) THEN
        DO I=LBOUND(TRGT,1),UBOUND(TRGT,1)
          IF (ALLOCATED(TRGT(I)%FACES)) THEN
            DEALLOCATE (TRGT(I)%FACES)
            DEALLOCATE (TRGT(I)%FCORI)
            DEALLOCATE (TRGT(I)%FLENGTH)
            DEALLOCATE (TRGT(I)%TET)
            DEALLOCATE (TRGT(I)%TIT)
            DEALLOCATE (TRGT(I)%TNT)
            DEALLOCATE (TRGT(I)%FE)
            DEALLOCATE (TRGT(I)%FSH)
            DEALLOCATE (TRGT(I)%FLUX)
            DEALLOCATE (TRGT(I)%DNIT)
            DEALLOCATE (TRGT(I)%VXT)
            DEALLOCATE (TRGT(I)%VYT)
            DEALLOCATE (TRGT(I)%VZT)
            DEALLOCATE (TRGT(I)%FI)
            DEALLOCATE (TRGT(I)%FEL)
            DEALLOCATE (TRGT(I)%VPART)
            DEALLOCATE (TRGT(I)%MACH)
            DEALLOCATE (TRGT(I)%USR)
            DEALLOCATE (TRGT(I)%ZIT)
          END IF
        ENDDO
        DEALLOCATE (TRGT)
      END IF

      RETURN
      END SUBROUTINE EIRENE_DEALLOC_BRAEIR


      SUBROUTINE EIRENE_INIT_BRAEIR

      DNIB    = 0.D0
      NEB     = 0.D0
      TEB     = 0.D0
      TIB     = 0.D0
      TNB     = 0.D0
      POB     = 0.D0
      VVB     = 0.D0
      UUB     = 0.D0
      WWB     = 0.D0
      PRB     = 0.D0
      UPB     = 0.D0
      RRB     = 0.D0
      VOLB    = 0.D0
      BFELDB  = 0.D0
      BBXB    = 0.D0
      BBYB    = 0.D0
      BBZB    = 0.D0
      EXB     = 0.D0
      EYB     = 0.D0
      EZB     = 0.D0
      FNIB    = 0.D0
      FEIB    = 0.D0
      FEEB    = 0.D0

      VXB  = 0.D0
      VYB  = 0.D0
      VZB  = 0.D0
      UUDIAB = 0.D0
      VVDIAB = 0.D0
      DELTA_SHEATHB = 0.D0

      amb=0.0_DP
      znb=0.0_DP
      zab=0.0_DP
      zib=0.0_DP

      PUX  = 0.D0
      PUY  = 0.D0
      PVX  = 0.D0
      PVY  = 0.D0

      RETURN
      END SUBROUTINE EIRENE_INIT_BRAEIR

      END MODULE EIRMOD_BRAEIR
