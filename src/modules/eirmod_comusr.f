cdr  may 2017:  preparing for storage reduction by elimination of unnecessary input tallies:
cdr             commenting,
cdr             lusr, musr, nusr, nplpr1, nplpr2, nsfprm made local,
cdr             rather than public
cpb  Dec. 2017: remove type SPECT_ARRAY, not needed in Fortran 2003
cdr             remove redundant tally LGDFT (also from LUSR)

      MODULE EIRMOD_COMUSR

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: EIRENE_ALLOC_COMUSR, EIRENE_DEALLOC_COMUSR,
     P          EIRENE_INIT_COMUSR, EIRENE_ALLOC_CORNERS,
     P          EIRENE_COMUSR_REINIT

      INTEGER, SAVE :: IFIRST=0
      INTEGER, SAVE ::
     P NPLPR1, NSFPRM, NPLPR2  ! internal, not public. former storage tests in setprm are abandoned
      INTEGER, PUBLIC, SAVE ::
     P NPLPRM  ! nplprm, is also used in setprm, for a storage test.
c
      INTEGER, SAVE ::
     P NUSR,   MUSR,   LUSR             ! also only local in this module, apparently

      REAL(DP), ALLOCATABLE, PUBLIC, SAVE ::
C  NPLPRM, REAL.
C  THE FIRST NPLPR1 DATA ARE PRIMARY INPUT PROFILES, SET IN SUBROUTINE PLASMA
     R        TEIN(:),   TIIN(:,:),   DEIN(:),   DIIN(:,:),
     R        VXIN(:,:), VYIN(:,:),   VZIN(:,:),
     R        BXIN(:),   BYIN(:),     BZIN(:),   BFIN(:),
     R        ADIN(:,:), VOL(:),      WGHT(:,:),
     R        EXIN(:),   EYIN(:),     EZIN(:),   EFIN(:),
     R        POT(:),
C  NSFPRM
     R        FLXOUT(:), SAREA(:),
C  NPLPR2, REAL.
C  THIS SECOND SET OF DATA ARE DERIVED INPUT PROFILES, SET IN SUBROUTINE PLASMA_DERIV
C  (STRICTLY ALSO DEIN (ELECTRON DENSITY) FROM THE NPLPR1 BLOCK ABOVE
C   IS SUCH A DERIVED QUANTITY)
     R        TEINL(:),  TIINL(:,:),  DEINL(:),  DIINL(:,:),
     R        BVIN(:,:), PARMOM(:,:), EDRIFT(:,:),
     R        BXPERP(:), BYPERP(:),
C
     R        DIOD(:),   DATD(:),     DMLD(:),   DPLD(:),    DPHD(:),
     R        DION(:),   DATM(:),     DMOL(:),   DPLS(:),    DPHOT(:)

!  DECLARATION AS TARGET ARRAYS FOR POINTERS USED BY UNIFIED SUBROUTINES
      REAL(DP), TARGET, ALLOCATABLE, PUBLIC, SAVE ::
     R        RMASSI(:), RMASSA(:),   RMASSM(:), RMASSPH(:), RMASSP(:)

!     POINTER FOR UNIFIED SUBROUTINES
      REAL(DP), POINTER, PUBLIC, SAVE :: RMASSX

C     PLASMA PROFILES ON CELL VERTICES
      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     R        CORNER_PROFILES(:,:)

c  storage for setting tallies at cell vertices, rather than cell centres,
c  for interpolations
      REAL(DP), POINTER, PUBLIC, SAVE ::
     .        TEINCORNER(:),   TIINCORNER(:,:), DEINCORNER(:),
     .        DIINCORNER(:,:),
     .        VXINCORNER(:,:), VYINCORNER(:,:), VZINCORNER(:,:),
     .        BXINCORNER(:),   BYINCORNER(:),   BZINCORNER(:),
     .        BFINCORNER(:),   BVINCORNER(:,:),
     .        EXCORNER(:),     EYCORNER(:),     EZCORNER(:),
     .        EFCORNER(:),     POTCORNER(:)

      REAL(DP), PUBLIC, SAVE :: TVAC, DVAC, VVAC, ALLOC

      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R TEDTEDX(:), TEDTEDY(:), TEDTEDZ(:)

      CHARACTER(8), ALLOCATABLE, PUBLIC, SAVE :: TEXTS(:)

C  MUSR, INTEGER
      INTEGER, PUBLIC, SAVE ::
     I         NSPH  , NPHOTI, NPHOTIM, NPHOTI_IN,
     I         NSPA  , NATMI,  NATMIM,  NATMI_IN,
     I         NSPAM , NMOLI,  NMOLIM,  NMOLI_IN,
     I         NSPAMI, NIONI,  NIONIM,  NIONI_IN,
     I         NSPTOT, NPLSI,  NPLSIM,  NPLSI_IN,
     I         NSNVI,  NCPVI,  NADVI,   NBGVI,
     I         NALVI,  NCLVI,  NADSI,   NALSI, NAINI, NBITS
      INTEGER, ALLOCATABLE, PUBLIC, SAVE ::
     I         NMASSA(:), NCHARA(:), NFOLA(:),  NGENA(:),
     I         NMASSM(:), NCHARM(:), NFOLM(:),  NGENM(:),
     I         NMASSI(:), NCHARI(:), NCHRGI(:), NFOLI(:), NGENI(:),
     I         NMASSP(:), NCHARP(:), NCHRGP(:),
     I         NFOLPH(:), NGENPH(:),
     I         NPRT(:),   ISPEZ(:,:,:,:,:,:),   ISPEZI(:,:),
     I         MPLSTI(:), MPLSV(:)
      INTEGER, ALLOCATABLE, PUBLIC, SAVE ::
     I         ISPZ_BACK(:,:)

C  LUSR, LOGICAL
      LOGICAL, ALLOCATABLE, PUBLIC, SAVE ::
     L         LGVAC(:,:)
      LOGICAL, PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     L         LSMOPRO(:)
      LOGICAL, PUBLIC, POINTER, SAVE ::
     L         LTESMO, LTISMO, LDESMO, LDISMO,    ! ldesmo: unused?
     L         LVSMO,  LBSMO,  LESMO,  LPOTSMO

C FROM HERE ON: NO EQUIVALENCE
      INTEGER, ALLOCATABLE, PUBLIC, SAVE ::
     I         IADVE(:),  IADVS(:), IADVT(:),  IADRC(:),
     I         ICLVE(:),  ICLVS(:), ICLVT(:),  ICLRC(:),
     I         ISNVE(:),  ISNVS(:), ISNVT(:),  ISNRC(:),
     I         ICPVE(:),  ICPVS(:), ICPVT(:),  ICPRC(:),
     I         IBGVE(:),  IBGVS(:), IBGVT(:),  IBGRC(:),
     I         IADSE(:),  IADSS(:), IADST(:),  IADSC(:),
     I         NFRSTP(:), NADDP(:), NSPAN(:),  NSPEN(:),
     I         NSPANW(:), NSPENW(:)

      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I         NAINS(:), NAINT(:)

      INTEGER, PUBLIC, SAVE ::
     I         NPRLL, NMODE,  NTCPU,
     I         NFILE, NFILEN, NFILEM, NFILEL, NFILEK, NFILEJ,
     I         NITER, IITER,  NTIME,  ITIMV

!      TYPE(SPECT_ARRAY), PUBLIC, ALLOCATABLE, SAVE :: BACK_SPEC(:)
      TYPE(EIRENE_SPECTRUM), PUBLIC, ALLOCATABLE, SAVE :: BACK_SPEC(:)
      LOGICAL, PUBLIC, ALLOCATABLE, SAVE :: LSPCCLL(:)


      CONTAINS

      SUBROUTINE EIRENE_ALLOC_COMUSR (ICAL)

      INTEGER, INTENT(IN) :: ICAL

      IF (.NOT.ALLOCATED(LSMOPRO)) ALLOCATE (LSMOPRO(12))  ! used only lsmopro(1:7) ?

      IF (ICAL == 1) THEN

        IF (ALLOCATED(TEIN)) RETURN
c
        NPLPR1=(12+1*NPLS+NPLSTI+3*NPLSV)*NRAD  ! background data, set in plasma.f, 17 arrays
        NPLPRM=NPLPR1+(NAIN+NSPZMC)*NRAD        !  adin, wght,...??? adin is allocated in call with ICAL == 2
cdr BVIN: add nplsv to nplpr2 and remove npls from nplprm. tbd:  check correct dimension of bvin !
        NPLPR2=(4+3*NPLS+NPLSTI+NPLSV)*NRAD+    ! additional background data, set in plasma_deriv.f, currently 9 arrays
     .          3*(NATM+NMOL+NION+NPLS)+4+NSPZ+2*NPHOT ! species (test particle and background) related data
        NUSR=NPLPRM+NPLPR2

        MUSR=4*NATM+4*NMOL+5*NION+3*NPLS+30+NSPZ+
     .       6*(1+NPHOTP)*(1+NATMP)*(1+NMOLP)*(1+NIONP)*(1+NPLSP)+NSPZ*6
     .       +2*NPLS+NSPZ*NPLS

        LUSR=NRAD*(NPLS+2)+NRAD   ! lgvac + lspccll

C NPLPR1 + ... = NPLPRM
        ALLOCATE (TEIN(NRAD))
        ALLOCATE (TIIN(NPLSTI,NRAD))
        ALLOCATE (DEIN(NRAD))   !  ital=-3,  derived tally
        ALLOCATE (DIIN(NPLS,NRAD))
        ALLOCATE (VXIN(NPLSV,NRAD))
        ALLOCATE (VYIN(NPLSV,NRAD))
        ALLOCATE (VZIN(NPLSV,NRAD))
        ALLOCATE (BXIN(NRAD))
        ALLOCATE (BYIN(NRAD))
        ALLOCATE (BZIN(NRAD))
        ALLOCATE (BFIN(NRAD))
cdr     ALLOCATE (ADIN(NAIN,NRAD))    !  ital=-12. Not yet. done later below, ical == 2 option
        ALLOCATE (VOL(NRAD))   !  ital=-14
        ALLOCATE (WGHT(NSPZMC,NRAD))  ! check size of  nspzmc.  this "weight window" array is unused so far.
        ALLOCATE (EXIN(NRAD))
        ALLOCATE (EYIN(NRAD))
        ALLOCATE (EZIN(NRAD))
        ALLOCATE (EFIN(NRAD))
        ALLOCATE (POT(NRAD))      !  ital=-22
c NPLPR2
        ALLOCATE (TEINL(NRAD))
        ALLOCATE (TIINL(NPLSTI,NRAD))
        ALLOCATE (BVIN(NPLSV,NRAD))   ! ital=nn
        ALLOCATE (PARMOM(NPLS,NRAD))  ! ital=nn
        ALLOCATE (BXPERP(NRAD))       ! ital=-16
        ALLOCATE (BYPERP(NRAD))       ! ital=-17
        ALLOCATE (EDRIFT(NPLS,NRAD))  ! ital=-13
        ALLOCATE (DEINL(NRAD))
        ALLOCATE (DIINL(NPLS,NRAD))

        ALLOCATE (RMASSA(MAX(1,NATM)))
        ALLOCATE (RMASSM(MAX(1,NMOL)))
        ALLOCATE (RMASSI(MAX(1,NION)))
        ALLOCATE (RMASSPH(MAX(1,NPHOT)))
        ALLOCATE (RMASSP(MAX(1,NPLS)))

        ALLOCATE (DIOD(MAX(1,NION)))
        ALLOCATE (DATD(MAX(1,NATM)))
        ALLOCATE (DMLD(MAX(1,NMOL)))
        ALLOCATE (DPLD(MAX(1,NPLS)))
        ALLOCATE (DPHD(MAX(1,NPHOT)))
        ALLOCATE (DION(MAX(1,NION)))
        ALLOCATE (DATM(MAX(1,NATM)))
        ALLOCATE (DMOL(MAX(1,NMOL)))
        ALLOCATE (DPLS(MAX(1,NPLS)))
        ALLOCATE (DPHOT(MAX(1,NPHOT)))

c  3 nrtal tallies ?  only for thermal force ??  size of nrtal ??
        ALLOCATE (TEDTEDX(NRTAL))  ! ital=nn
        ALLOCATE (TEDTEDY(NRTAL))  ! ital=nn
        ALLOCATE (TEDTEDZ(NRTAL))  ! ital=nn

        ALLOCATE (TEXTS(NSPZ))
c  integer species and background tally data

        ALLOCATE (NMASSA(MAX(1,NATM)))
        ALLOCATE (NCHARA(MAX(1,NATM)))
        ALLOCATE (NFOLA(MAX(1,NATM)))
        ALLOCATE (NGENA(MAX(1,NATM)))
        ALLOCATE (NMASSM(MAX(1,NMOL)))
        ALLOCATE (NCHARM(MAX(1,NMOL)))
        ALLOCATE (NFOLM(MAX(1,NMOL)))
        ALLOCATE (NGENM(MAX(1,NMOL)))
        ALLOCATE (NMASSP(MAX(1,NPLS)))
        ALLOCATE (NCHARP(MAX(1,NPLS)))
        ALLOCATE (NCHRGP(MAX(1,NPLS)))
        ALLOCATE (NMASSI(MAX(1,NION)))
        ALLOCATE (NCHARI(MAX(1,NION)))
        ALLOCATE (NCHRGI(MAX(1,NION)))
        ALLOCATE (NFOLI(MAX(1,NION)))
        ALLOCATE (NGENI(MAX(1,NION)))
        ALLOCATE (NFOLPH(MAX(1,NPHOT)))
        ALLOCATE (NGENPH(MAX(1,NPHOT)))
        ALLOCATE (NPRT(MAX(1,NSPZ)))
        ALLOCATE (ISPEZ(-1:4,0:NPHOTP,0:NATMP,0:NMOLP,0:NIONP,0:NPLSP))
        ALLOCATE (ISPEZI(NSPZ,-1:4))

        ALLOCATE (MPLSTI(MAX(1,NPLS)))
        ALLOCATE (MPLSV(MAX(1,NPLS)))
        ALLOCATE (ISPZ_BACK(NSPZ,NPLS))
        ALLOCATE (IADVE(MAX(1,NADV)))
        ALLOCATE (IADVS(MAX(1,NADV)))
        ALLOCATE (IADVT(MAX(1,NADV)))
        ALLOCATE (IADRC(MAX(1,NADV)))
        ALLOCATE (ICLVE(MAX(1,NCLV)))
        ALLOCATE (ICLVS(MAX(1,NCLV)))
        ALLOCATE (ICLVT(MAX(1,NCLV)))
        ALLOCATE (ICLRC(MAX(1,NCLV)))
        ALLOCATE (ISNVE(MAX(1,NSNV)))
        ALLOCATE (ISNVS(MAX(1,NSNV)))
        ALLOCATE (ISNVT(MAX(1,NSNV)))
        ALLOCATE (ISNRC(MAX(1,NSNV)))
        ALLOCATE (IADSE(MAX(1,NADS)))
        ALLOCATE (IADSS(MAX(1,NADS)))
        ALLOCATE (IADST(MAX(1,NADS)))
        ALLOCATE (IADSC(MAX(1,NADS)))
        ALLOCATE (NFRSTP(NTALI))
        ALLOCATE (NADDP(NTALI))
        ALLOCATE (NSPAN(NTALV))
        ALLOCATE (NSPEN(NTALV))
        ALLOCATE (NSPANW(NTALS))
        ALLOCATE (NSPENW(NTALS))
c  logicals
        ALLOCATE (LGVAC(NRAD,0:NPLS+1))
        ALLOCATE (LSPCCLL(NRAD))

        WRITE (IUNMEM,'(A,T25,I15)')
     .        ' COMUSR(1) ',NUSR*8 + MUSR*4 + 
     .                     (LUSR+12)*4 + !  lgvac+lspccll+lsmopro
     .                      3*NRTAL*8

      ELSE IF (ICAL == 2) THEN
c  NAIN: first dimension of adin is now fixed.  correct nplprm with nain*nrad
        IF (ALLOCATED(ADIN)) RETURN

        NPLPR1=(12+1*NPLS+NPLSTI+3*NPLSV)*NRAD
        NPLPRM=NPLPR1+(NAIN+NSPZMC)*NRAD
        ALLOCATE (ADIN(NAIN,NRAD))

        ALLOCATE (NAINS(NAIN))
        ALLOCATE (NAINT(NAIN))

c  NCPV, NBGV are now set
        ALLOCATE (ICPVE(NCPV))
        ALLOCATE (ICPVS(NCPV))
        ALLOCATE (ICPVT(NCPV))
        ALLOCATE (ICPRC(NCPV))
        ALLOCATE (IBGVE(NBGV))
        ALLOCATE (IBGVS(NBGV))
        ALLOCATE (IBGVT(NBGV))
        ALLOCATE (IBGRC(NBGV))

      ELSE IF (ICAL == 3) THEN

        IF (ALLOCATED(FLXOUT)) RETURN

        NSFPRM=2*NLMPGS
        ALLOCATE (FLXOUT(NLMPGS))
        ALLOCATE (SAREA(NLMPGS))

        WRITE (IUNMEM,'(A,T25,I15)')
     .         ' COMUSR(3) ',NSFPRM*8

      END IF

      CALL EIRENE_INIT_COMUSR(ICAL)

      RETURN
      END SUBROUTINE EIRENE_ALLOC_COMUSR



      SUBROUTINE EIRENE_ALLOC_CORNERS(IUNOUT)

      INTEGER, INTENT(IN) :: IUNOUT
      INTEGER :: N1DIM(12)
      INTEGER :: NTOT, ICO


      IF (ALLOCATED(CORNER_PROFILES)) RETURN

      N1DIM = 0
      NTOT = 0
cdr  are there any FEM interpolated background tallies in this run?

      IF (LTESMO) NTOT = NTOT + 1
      IF (LTISMO) NTOT = NTOT + NPLSTI
cdr   IF (LDESMO) ....?  unused
      IF (LDISMO) NTOT = NTOT + NPLS + 1
      IF (LVSMO)  NTOT = NTOT + 4*NPLSV
      IF (LBSMO)  NTOT = NTOT + 4
      IF (LESMO)  NTOT = NTOT + 4
      IF (LPOTSMO)  NTOT = NTOT + 1

      IF (NTOT > 0) THEN
cdr  allocate storage for background tallies on cell vertices
cdr  ncorner is set in GRID.f (levgeo=4,5) or in SNEIGH.f (levgeo=1,2,3)
        ALLOCATE (CORNER_PROFILES(NCORNER,NTOT))
      ELSE
        ALLOCATE (CORNER_PROFILES(1,1))
      END IF

      ICO = 0

      IF (LTESMO) THEN
        TEINCORNER => CORNER_PROFILES(:,ICO+1)
        ICO = ICO + 1
      ELSE
        NULLIFY(TEINCORNER)
      END IF

      IF (LTISMO) THEN
        TIINCORNER => CORNER_PROFILES(:,ICO+1 : ICO+NPLSTI )
        ICO = ICO + NPLSTI
      ELSE
        NULLIFY(TIINCORNER)
      END IF

      IF (LDISMO) THEN
        DEINCORNER => CORNER_PROFILES(:,ICO+1)
        DIINCORNER => CORNER_PROFILES(:,ICO+2 : ICO+NPLS+1 )
        ICO = ICO + NPLS + 1
      ELSE
        NULLIFY(DEINCORNER)
        NULLIFY(DIINCORNER)
      END IF

      IF (LVSMO) THEN
        VXINCORNER => CORNER_PROFILES(:,ICO+1 : ICO+NPLSV )
        VYINCORNER => CORNER_PROFILES(:,ICO+1+1*NPLSV : ICO+2*NPLSV )
        VZINCORNER => CORNER_PROFILES(:,ICO+1+2*NPLSV : ICO+3*NPLSV )
        BVINCORNER => CORNER_PROFILES(:,ICO+1+3*NPLSV : ICO+4*NPLSV )
        ICO = ICO + 4*NPLSV
      ELSE
        NULLIFY(VXINCORNER)
        NULLIFY(VYINCORNER)
        NULLIFY(VZINCORNER)
        NULLIFY(BVINCORNER)
      END IF

      IF (LBSMO) THEN
        BXINCORNER => CORNER_PROFILES(:,ICO+1)
        BYINCORNER => CORNER_PROFILES(:,ICO+2)
        BZINCORNER => CORNER_PROFILES(:,ICO+3)
        BFINCORNER => CORNER_PROFILES(:,ICO+4)
        ICO = ICO + 4
      ELSE
        NULLIFY(BXINCORNER)
        NULLIFY(BYINCORNER)
        NULLIFY(BZINCORNER)
        NULLIFY(BFINCORNER)
      END IF

      IF (LESMO) THEN
        EXCORNER => CORNER_PROFILES(:,ICO+1)
        EYCORNER => CORNER_PROFILES(:,ICO+2)
        EZCORNER => CORNER_PROFILES(:,ICO+3)
        EFCORNER => CORNER_PROFILES(:,ICO+4)
        ICO = ICO + 4
      ELSE
        NULLIFY(EXCORNER)
        NULLIFY(EYCORNER)
        NULLIFY(EZCORNER)
        NULLIFY(EFCORNER)
      END IF

      IF (LPOTSMO) THEN
        POTCORNER => CORNER_PROFILES(:,ICO+1)
        ICO = ICO + 1
      ELSE
        NULLIFY(POTCORNER)
      END IF

      IF (ICO /= NTOT) THEN
        WRITE (IUNOUT,*) ' ERROR IN EIRENE_ALLOC_CORNERS '
        WRITE (IUNOUT,*) ' NTOT = ',NTOT,' /= ICO = ',ICO
        CALL EIRENE_EXIT_OWN(1)
      END IF

      CORNER_PROFILES = 0._DP

      END SUBROUTINE EIRENE_ALLOC_CORNERS



      SUBROUTINE EIRENE_DEALLOC_COMUSR
C
      IF (.NOT.ALLOCATED(TEIN)) RETURN

      DEALLOCATE (TEIN)
      DEALLOCATE (TIIN)
      DEALLOCATE (DEIN)
      DEALLOCATE (DIIN)
      DEALLOCATE (VXIN)
      DEALLOCATE (VYIN)
      DEALLOCATE (VZIN)
      DEALLOCATE (BXIN)
      DEALLOCATE (BYIN)
      DEALLOCATE (BZIN)
      DEALLOCATE (BFIN)
      DEALLOCATE (ADIN)
      DEALLOCATE (VOL)
      DEALLOCATE (WGHT)
      DEALLOCATE (BXPERP)
      DEALLOCATE (BYPERP)
      DEALLOCATE (EXIN)
      DEALLOCATE (EYIN)
      DEALLOCATE (EZIN)
      DEALLOCATE (EFIN)
      DEALLOCATE (POT)
c
      DEALLOCATE (TEINL)
      DEALLOCATE (TIINL)
      DEALLOCATE (BVIN)
      DEALLOCATE (PARMOM)
      DEALLOCATE (EDRIFT)
      DEALLOCATE (DEINL)
      DEALLOCATE (DIINL)

      DEALLOCATE (FLXOUT)
      DEALLOCATE (SAREA)


      DEALLOCATE (RMASSA)
      DEALLOCATE (RMASSM)
      DEALLOCATE (RMASSI)
      DEALLOCATE (RMASSPH)
      DEALLOCATE (RMASSP)

      DEALLOCATE (DIOD)
      DEALLOCATE (DATD)
      DEALLOCATE (DMLD)
      DEALLOCATE (DPLD)
      DEALLOCATE (DPHD)
      DEALLOCATE (DION)
      DEALLOCATE (DATM)
      DEALLOCATE (DMOL)
      DEALLOCATE (DPLS)
      DEALLOCATE (DPHOT)

      DEALLOCATE (TEDTEDX)
      DEALLOCATE (TEDTEDY)
      DEALLOCATE (TEDTEDZ)

      DEALLOCATE (TEXTS)
      DEALLOCATE (NMASSA)
      DEALLOCATE (NCHARA)
      DEALLOCATE (NFOLA)
      DEALLOCATE (NGENA)
      DEALLOCATE (NMASSM)
      DEALLOCATE (NCHARM)
      DEALLOCATE (NFOLM)
      DEALLOCATE (NGENM)
      DEALLOCATE (NMASSP)
      DEALLOCATE (NCHARP)
      DEALLOCATE (NCHRGP)
      DEALLOCATE (NMASSI)
      DEALLOCATE (NCHARI)
      DEALLOCATE (NCHRGI)
      DEALLOCATE (NFOLI)
      DEALLOCATE (NGENI)
      DEALLOCATE (NFOLPH)
      DEALLOCATE (NGENPH)
      DEALLOCATE (NPRT)
      DEALLOCATE (ISPEZ)
      DEALLOCATE (ISPEZI)
      DEALLOCATE (MPLSTI)
      DEALLOCATE (MPLSV)
      DEALLOCATE (ISPZ_BACK)
      DEALLOCATE (IADVE)
      DEALLOCATE (IADVS)
      DEALLOCATE (IADVT)
      DEALLOCATE (IADRC)
      DEALLOCATE (ICLVE)
      DEALLOCATE (ICLVS)
      DEALLOCATE (ICLVT)
      DEALLOCATE (ICLRC)
      DEALLOCATE (ISNVE)
      DEALLOCATE (ISNVS)
      DEALLOCATE (ISNVT)
      DEALLOCATE (ISNRC)
      DEALLOCATE (ICPVE)
      DEALLOCATE (ICPVS)
      DEALLOCATE (ICPVT)
      DEALLOCATE (ICPRC)
      DEALLOCATE (IBGVE)
      DEALLOCATE (IBGVS)
      DEALLOCATE (IBGVT)
      DEALLOCATE (IBGRC)
      DEALLOCATE (IADSE)
      DEALLOCATE (IADSS)
      DEALLOCATE (IADST)
      DEALLOCATE (IADSC)
      DEALLOCATE (NFRSTP)
      DEALLOCATE (NADDP)
      DEALLOCATE (NSPAN)
      DEALLOCATE (NSPEN)
      DEALLOCATE (NSPANW)
      DEALLOCATE (NSPENW)
      DEALLOCATE (NAINS)
      DEALLOCATE (NAINT)
      DEALLOCATE (LGVAC)
      DEALLOCATE (LSPCCLL)
      DEALLOCATE (LSMOPRO)

!pb      IF (NBACK_SPEC > 0) DEALLOCATE (BACK_SPEC)
      IF (ALLOCATED(BACK_SPEC)) DEALLOCATE (BACK_SPEC)

      IF (ALLOCATED(CORNER_PROFILES)) DEALLOCATE (CORNER_PROFILES)

      RETURN
      END SUBROUTINE EIRENE_DEALLOC_COMUSR


      SUBROUTINE EIRENE_INIT_COMUSR(ICAL)

      INTEGER, INTENT(IN) :: ICAL

      IF (IFIRST == 0) THEN
        LSMOPRO = .FALSE.

        LTESMO   => LSMOPRO(1)
        LTISMO   => LSMOPRO(2)
        LDISMO   => LSMOPRO(3)
        LVSMO    => LSMOPRO(4)
        LBSMO    => LSMOPRO(5)
        LESMO    => LSMOPRO(6)
        LPOTSMO  => LSMOPRO(7)

        IFIRST = 1
      ENDIF

      IF (ICAL == 1) THEN

        TEIN   = 0._DP
        TIIN   = 0._DP
        DEIN   = 0._DP
        DIIN   = 0._DP
        VXIN   = 0._DP
        VYIN   = 0._DP
        VZIN   = 0._DP
        BXIN   = 0._DP
        BYIN   = 0._DP
        BZIN   = 0._DP
        BFIN   = 0._DP
        VOL    = 0._DP       ! ital=-14
        WGHT   = 1._DP
        BXPERP = 0._DP
        BYPERP = 0._DP
        EXIN   = 0._DP
        EYIN   = 0._DP
        EZIN   = 0._DP
        EFIN   = 0._DP
        POT    = 0._DP
        TEINL  = 0._DP
        TIINL  = 0._DP
        BVIN   = 0._DP
        PARMOM = 0._DP
        EDRIFT = 0._DP
        DEINL  = 0._DP
        DIINL  = 0._DP

        RMASSA = 0._DP
        RMASSM = 0._DP
        RMASSI = 0._DP
        RMASSPH = 0._DP
        RMASSP = 0._DP

        DIOD   = 0._DP
        DATD   = 0._DP
        DMLD   = 0._DP
        DPLD   = 0._DP
        DPHD   = 0._DP
        DION   = 0._DP
        DATM   = 0._DP
        DMOL   = 0._DP
        DPLS   = 0._DP
        DPHOT  = 0._DP

        TEDTEDX = 0._DP
        TEDTEDY = 0._DP
        TEDTEDZ = 0._DP

        TEXTS  = ' '
        NMASSA = 0
        NCHARA = 0
        NFOLA  = 0
        NGENA  = 0
        NMASSM = 0
        NCHARM = 0
        NFOLM  = 0
        NGENM  = 0
        NMASSP = 0
        NCHARP = 0
        NCHRGP = 0
        NMASSI = 0
        NCHARI = 0
        NCHRGI = 0
        NFOLI  = 0
        NGENI  = 0
        NFOLPH = 0
        NGENPH = 0
        NPRT   = 0
        ISPEZ  = 0
        ISPEZI = 0
        MPLSTI = 0
        MPLSV  = 0
        ISPZ_BACK = 0
        IADVE  = 0
        IADVS  = 0
        IADVT  = 0
        IADRC  = 0
        ICLVE  = 0
        ICLVS  = 0
        ICLVT  = 0
        ICLRC  = 0
        ISNVE  = 0
        ISNVS  = 0
        ISNVT  = 0
        ISNRC  = 0
        IADSE  = 0
        IADSS  = 0
        IADST  = 0
        IADSC  = 0
        NFRSTP = 0
        NADDP  = 0
        NSPAN  = 0
        NSPEN  = 0
        NSPANW = 0
        NSPENW = 0
        LGVAC  = .FALSE.
        LSPCCLL = .FALSE.

      ELSE IF (ICAL == 2) THEN
c  at this call: first dimension of adin is known, as well as size of cop and bgk tallies
        ADIN   = 0._DP      ! ital=-12

        NAINS = 0
        NAINT = 0

        ICPVE  = 0
        ICPVS  = 0
        ICPVT  = 0
        ICPRC  = 0
        IBGVE  = 0
        IBGVS  = 0
        IBGVT  = 0
        IBGRC  = 0

      ELSE IF (ICAL == 3) THEN

        FLXOUT = 0._DP
        SAREA  = 666._DP

      END IF

      RETURN
      END SUBROUTINE EIRENE_INIT_COMUSR

      SUBROUTINE EIRENE_COMUSR_REINIT
      IMPLICIT NONE
      IFIRST = 0
      RETURN

      END SUBROUTINE EIRENE_COMUSR_REINIT

      END MODULE EIRMOD_COMUSR
