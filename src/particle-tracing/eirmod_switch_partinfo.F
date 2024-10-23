cdr feb 22: some tallies missing: pxph, exph, gen. lim. tallies, ....
cdr         Now added.
cdr 2022: reworking some of this coding, to enable condensation of
cdr       colatm,colmol,colion into the single routine collide.
cdr       (E.g. for more transparent (unified) generation limit, cascading
cdr       and species scaling implementations.
cdr       Also: accommodate photon tracing tallies now.

cdr        PXX, PXX2 pointers had excess species indices. Now removed.
cdr        LGX.. etc... pointers had excess species index. Now removed
cdr mar 22: done.
cdr
cdr         3D time_array split into 4 separate 2D time_ar, one for each type
cdr         for compiler complains from MASYR1 printout routine.

      module eirmod_switch_partinfo
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT
      USE EIRMOD_CESTIM
      USE EIRMOD_COMXS
      USE EIRMOD_CZT1
      USE EIRMOD_CSPEZ
      USE EIRMOD_CTRCEI
      USE EIRMOD_COMSOU
      USE EIRMOD_CLOGAU
      USE EIRMOD_SECOND_OWN, ONLY: eirene_second_own

      IMPLICIT NONE
cdr for trchktim option: cpu time consumption split by species
      real(dp), allocatable, save :: time_ar0(:,:),
     .                               time_ar1(:,:),
     .                               time_ar2(:,:),
     .                               time_ar3(:,:)
      integer, save :: istra_old=-1,
     .                 ityp_old=-1,
     .                 iphot_old=-1,
     .                 iatm_old=-1,
     .                 imol_old=-1,
     .                 iion_old=-1,
     .                 ipls_old=-1

!$omp  threadprivate(istra_old,ityp_old,iphot_old,iatm_old,imol_old,
!$omp&   iion_old,ipls_old)

      PRIVATE

      PUBLIC :: eirene_switch_partinfo, eirene_output_partinfo,
     .          eirene_reinit_partinfo

      CONTAINS

cdr  Jan 18: bypass this actions for photons (ityp=0). Code not ready for photon transport.
cpb:  added: cpu time statistics by particle type, species and stratum: time_array
cdr: Apr.18: testing, cleaning of time_array options (minor bug fix)
cdr          further photonic arrays added (targets,pointer)
      subroutine eirene_switch_partinfo
c  added oct. 2017:
c  this routine sets the various pointers for tallies,
c  for a unified treatment of scoring volume tallies (update),
c                                     mfp evaluation (fpath)
c  and a for scoring a number of surface tallies (locate, escape, ...)

c  It is called from: LOCATE (for primary source particles),
c          and after: COLLIDE and REFLEC (for new post-collision species)
c
c  Input:  istra, ityp, iphot, iatm, imol, iion, ipls
c  Output: ixspz, nmetoff, logphot, logatm, logmol, logion

      implicit none

      real(dp) :: tim_spent
      real(dp), save :: tim_start=0._dp, tim_end=0._dp

      if ((ityp_old == ityp) .and. (iatm_old == iatm) .and.
     .    (imol_old == imol) .and. (iion_old == iion) .and.
     .    (iphot_old == iphot) .and. (ipls_old == ipls).and.
     .    (istra_old == istra)) return

cdr  this next part: split CPU time consumption
cdr  by test particle species and type.

      if (trchktim) then
        tim_end = EIRENE_SECOND_OWN()
        tim_spent = tim_end - tim_start
        tim_start = tim_end

        select case(ityp_old)
        case(0)
!  photons
          time_ar0(iphot_old,istra_old) =
     .    time_ar0(iphot_old,istra_old) + tim_spent
        case(1)
!  atoms
          time_ar1(iatm_old,istra_old) =
     .    time_ar1(iatm_old,istra_old) + tim_spent
        case(2)
!  molecules
          time_ar2(imol_old,istra_old) =
     .    time_ar2(imol_old,istra_old) + tim_spent
        case(3)
!  test ions
          time_ar3(iion_old,istra_old) =
     .    time_ar3(iion_old,istra_old) + tim_spent
        case(4)
!  bulk ions: NOT IN USE
!         time_ar4(ipls_old,istra_old) =
!     .   time_ar4(ipls_old,istra_old) + tim_spent
        case default
cdr  initialization:
cdr  first call is with ityp_old=-1
          if (.not.allocated(time_ar0)) then
            allocate (time_ar0(0:nphot,0:nstra))
            allocate (time_ar1(0:natm,0:nstra))
            allocate (time_ar2(0:nmol,0:nstra))
            allocate (time_ar3(0:nion,0:nstra))
            time_ar0 = 0._dp
            time_ar1 = 0._dp
            time_ar2 = 0._dp
            time_ar3 = 0._dp
          end if
        end select
      end if

      if (.not.allocated(lgxcx)) then
        ALLOCATE (LGXCX(0:NRCX,0:1))
        ALLOCATE (LGXEI(0:NREI))
        ALLOCATE (LGXEL(0:NREL,0:1))
        ALLOCATE (LGXPI(0:NRPI,0:1))
      endif

C  save stratum, old type, species
      istra_old= istra
      ityp_old = ityp

      iphot_old= iphot
      iatm_old = iatm
      imol_old = imol
      iion_old = iion
      ipls_old = ipls

      NULLIFY (PDENX)
      NULLIFY (EDENX)

      NULLIFY (PXEL)
      NULLIFY (PXAT)
      NULLIFY (PXML)
      NULLIFY (PXIO)
      NULLIFY (PXPHT)
      NULLIFY (PXPL)

      NULLIFY (EXEL)
      NULLIFY (EXAT)
      NULLIFY (EXML)
      NULLIFY (EXIO)
      NULLIFY (EXPHT)
      NULLIFY (EXPL)

      NULLIFY (VXDENX)
      NULLIFY (VYDENX)
      NULLIFY (VZDENX)

      NULLIFY (PGENX)
      NULLIFY (EGENX)
      NULLIFY (VGENX)

      NULLIFY (MXPL)
      NULLIFY (RXEL)

      NULLIFY (PXX)
      NULLIFY (EXX)

      NDXX = 0

      select case(ityp)

      case(0)
!  photons

       LPDENX  => LPDENPH
       LEDENX  => LEDENPH

       LPXEL   => LPPHEL
       LPXAT   => LPPHAT
       LPXML   => LPPHML
       LPXIO   => LPPHIO
       LPXPHT  => LPPHPHT  ! same as LPXX
       LPXPL   => LPPHPL

       LEXEL   => LEPHEL
       LEXAT   => LEPHAT
       LEXML   => LEPHML
       LEXIO   => LEPHIO
       LEXPHT  => LEPHPHT  ! same as LEXX
       LEXPL   => LEPHPL

       LVXDENX => LVXDENPH
       LVYDENX => LVYDENPH
       LVZDENX => LVZDENPH

       LPGENX  => LPGENPH
       LEGENX  => LEGENPH
       LVGENX  => LVGENPH

       LMXPL   => LMPHPL
       LRXEL   => LRPHEL

       LPXX    => LPPHPHT
       LEXX    => LEPHPHT
       LSCX    => NLSPCSCL_PHOT

       IF (LPDENX)  PDENX  => PDENPH(IPHOT,:)
       IF (LEDENX)  EDENX  => EDENPH(IPHOT,:)

       IF (LPXEL)   PXEL   => PPHEL(:)
       IF (LSCX) THEN
         IF (LPXAT) THEN
           NTS_PXATA = NTS_PI+1
           NTS_PXATE = NTS_APH
           PXAT   => PPHAT(1:NATM*NPHOTP,:)
         END IF
         IF (LPXML) THEN
           NTS_PXMLA = NTS_APH+1
           NTS_PXMLE = NTS_MPH
           PXML   => PPHML(1:NMOL*NPHOTP,:)
         END IF
         IF (LPXIO) THEN
           NTS_PXIOA = NTS_MPH+1
           NTS_PXIOE = NTS_IPH
           PXIO   => PPHIO(1:NION*NPHOTP,:)
         END IF
         IF (LPXPHT) THEN
           NTS_PXPHA = NTS_IPH+1
           NTS_PXPHE = NTS_PHPH
           PXPHT  => PPHPHT(1:NPHOT*NPHOTP,:)
         END IF
         IF (LPXPL) THEN
           NTS_PXPLA = NTS_PHPH+1
           NTS_PXPLE = NTS_PPH
           IF (LPXPL) PXPL   => PPHPL(1:NPLS*NPHOTP,:)
         END IF
       ELSE
         IF (LPXAT)  PXAT   => PPHAT(1:NATMI,:)
         IF (LPXML)  PXML   => PPHML(1:NMOLI,:)
         IF (LPXIO)  PXIO   => PPHIO(1:NIONI,:)
         IF (LPXPHT) PXPHT  => PPHPHT(1:NPHOTI,:)
         IF (LPXPL)  PXPL   => PPHPL(1:NPLSI,:)
       END IF

       IF (LEXEL)   EXEL   => EPHEL(:)
       IF (LEXAT)   EXAT   => EPHAT(:)
       IF (LEXML)   EXML   => EPHML(:)
       IF (LEXIO)   EXIO   => EPHIO(:)
       IF (LEXPHT)  EXPHT  => EPHPHT(:)     ! for species I=IPHOT:
                                            ! same as EXX
       IF (LEXPL)   EXPL   => EPHPL(1:NPLSI,:)

       IF (LVXDENX) VXDENX => VXDENPH(IPHOT,:)
       IF (LVYDENX) VYDENX => VYDENPH(IPHOT,:)
       IF (LVZDENX) VZDENX => VZDENPH(IPHOT,:)

       IF (LPGENX) PGENX   => PGENPH(IPHOT,:)
       IF (LEGENX) EGENX   => EGENPH(IPHOT,:)
       IF (LVGENX) VGENX   => VGENPH(IPHOT,:)

       IF (LMXPL)   MXPL   => MPHPL(1:NPLSI,:)
       IF (LRXEL) RXEL => RPHEL(1:NPHOTI,:)

       IF (LSCX) THEN
         NDXX = NPHOT
         NDXXA = NTS_IPH+1
         NDXXE = NTS_PHPH
         IF (LPXX)  PXX    => PPHPHT(1:NPHOT*NPHOT,:)
       ELSE
         IF (LPXX)  PXX    => PPHPHT(1:NPHOTI,:)
       END IF
       IF (LEXX)    EXX    => EPHPHT(:)

       LEX     => LEPH

!       LGXCX => LGPHCX
!       LGXEI => LGPHEI
!       LGXEL => LGPHEL
!       LGXPI => LGPHPI
cdr now: allocatable, rather than pointer
!pb that means: here we do a copy
       LGXCX(0:NRCX,0:1) = LGPHCX(IPHOT,0:NRCX,0:1)
       LGXEI(0:NREI)     = LGPHEI(IPHOT,0:NREI)
       LGXEL(0:NREL,0:1) = LGPHEL(IPHOT,0:NREL,0:1)
       LGXPI(0:NRPI,0:1) = LGPHPI(IPHOT,0:NRPI,0:1)

       NXEII => NPHEII(IPHOT)
       NXCXI => NPHCXI(IPHOT)
       NXELI => NPHELI(IPHOT)
       NXPII => NPHPII(IPHOT)

       NXEIIM => NPHEIIM(IPHOT)
       NXCXIM => NPHCXIM(IPHOT)
       NXELIM => NPHELIM(IPHOT)
       NXPIIM => NPHPIIM(IPHOT)

       NPBGKX => NPBGKPH(IPHOT)
       NGENX  => NGENPH(IPHOT)

       IXSPZ = IPHOT
       NMETOFF = 0

       RMASSX => RMASSPH(IPHOT)
       CNDYNX => CNDYNPH(IPHOT)
       CVRSSX => CVRSSPH(IPHOT)

       LOGPHOT(IPHOT,ISTRA)=.TRUE.
       LOGXSPZ => LOGPHOT(IPHOT,ISTRA)

      case(1)
!  atoms

       LPDENX  => LPDENA
       LEDENX  => LEDENA

       LPXEL   => LPAEL
       LPXAT   => LPAAT  ! same as LPXX
       LPXML   => LPAML
       LPXIO   => LPAIO
       LPXPHT  => LPAPHT
       LPXPL   => LPAPL

       LEXEL   => LEAEL
       LEXAT   => LEAAT  ! same as LEXX
       LEXML   => LEAML
       LEXIO   => LEAIO
       LEXPHT  => LEAPHT
       LEXPL   => LEAPL

       LVXDENX => LVXDENA
       LVYDENX => LVYDENA
       LVZDENX => LVZDENA

       LPGENX  => LPGENA
       LEGENX  => LEGENA
       LVGENX  => LVGENA

       LMXPL   => LMAPL

       LRXEL   => LRAEL
       LPXX    => LPAAT
       LEXX    => LEAAT
       LSCX    => NLSPCSCL_ATM

       IF (LPDENX)  PDENX  => PDENA(IATM,:)
       IF (LEDENX)  EDENX  => EDENA(IATM,:)

       IF (LPXEL)   PXEL   => PAEL(:)
       IF (LSCX) THEN
         IF (LPXAT) THEN
           NTS_PXATA = NSPZTOTS+1
           NTS_PXATE = NTS_AA
           PXAT   => PAAT(1:NATM*NATMP,:)
         END IF
         IF (LPXML) THEN
           NTS_PXMLA = NTS_AA+1
           NTS_PXMLE = NTS_MA
           PXML   => PAML(1:NMOL*NATMP,:)
         END IF
         IF (LPXIO) THEN
           NTS_PXIOA = NTS_MA+1
           NTS_PXIOE = NTS_IA
           PXIO   => PAIO(1:NION*NATMP,:)
         END IF
         IF (LPXPHT) THEN
           NTS_PXPHA = NTS_IA+1
           NTS_PXPHE = NTS_PHA
           PXPHT  => PAPHT(1:NPHOT*NATMP,:)
         END IF
         IF (LPXPL) THEN
           NTS_PXPLA = NTS_PHA+1
           NTS_PXPLE = NTS_PA
           PXPL   => PAPL(1:NPLS*NATMP,:)
         END IF
       ELSE
         IF (LPXAT)  PXAT   => PAAT(1:NATMI,:)  ! for species i=iatm:
                                                ! same as PXX
         IF (LPXML)  PXML   => PAML(1:NMOLI,:)
         IF (LPXIO)  PXIO   => PAIO(1:NIONI,:)
         IF (LPXPHT) PXPHT  => PAPHT(1:NPHOTI,:)
         IF (LPXPL)  PXPL   => PAPL(1:NPLSI,:)
       END IF

       IF (LEXEL)   EXEL   => EAEL(:)
       IF (LEXAT)   EXAT   => EAAT(:)   ! for species i=iatm:
                                        ! same as EXX
       IF (LEXML)   EXML   => EAML(:)
       IF (LEXIO)   EXIO   => EAIO(:)
       IF (LEXPHT)  EXPHT  => EAPHT(:)
       IF (LEXPL)   EXPL   => EAPL(1:NPLSI,:)

       IF (LVXDENX) VXDENX => VXDENA(IATM,:)
       IF (LVYDENX) VYDENX => VYDENA(IATM,:)
       IF (LVZDENX) VZDENX => VZDENA(IATM,:)

       IF (LRXEL)   RXEL   => RAEL(1:NATMI,:)

       IF (LPGENX)  PGENX  => PGENA(IATM,:)
       IF (LEGENX)  EGENX  => EGENA(IATM,:)
       IF (LVGENX)  VGENX  => VGENA(IATM,:)

       IF (LMXPL)   MXPL   => MAPL(1:NPLSI,:)

       IF (LSCX) THEN
         NDXX = NATM
         NDXXA = NSPZTOTS+1
         NDXXE = NTS_AA
         IF (LPXX)  PXX    => PAAT(1:NATM*NATMP,:)
       ELSE
         IF (LPXX)  PXX    => PAAT(1:NATMI,:)
       END IF
       IF (LEXX)    EXX    => EAAT(:)

       LEX     => LEA

!      LGXCX => LGACX
!      LGXEI => LGAEI
!      LGXEL => LGAEL
!      LGXPI => LGAPI
cdr now: allocatable, rather than pointer
!pb that means: here we do a copy
       LGXCX(0:NRCX,0:1) = LGACX(IATM,0:NRCX,0:1)
       LGXEI(0:NREI)     = LGAEI(IATM,0:NREI)
       LGXEL(0:NREL,0:1) = LGAEL(IATM,0:NREL,0:1)
       LGXPI(0:NRPI,0:1) = LGAPI(IATM,0:NRPI,0:1)

       NXEII => NAEII(IATM)
       NXCXI => NACXI(IATM)
       NXELI => NAELI(IATM)
       NXPII => NAPII(IATM)

       NXEIIM => NAEIIM(IATM)
       NXCXIM => NACXIM(IATM)
       NXELIM => NAELIM(IATM)
       NXPIIM => NAPIIM(IATM)

       NPBGKX => NPBGKA(IATM)
       NGENX  => NGENA(IATM)

       IXSPZ = IATM
       NMETOFF = NSPH

       RMASSX => RMASSA(IATM)
       CNDYNX => CNDYNA(IATM)
       CVRSSX => CVRSSA(IATM)

       LOGATM(IATM,ISTRA)=.TRUE.
       LOGXSPZ => LOGATM(IATM,ISTRA)

      case(2)
!  molecules

       LPDENX  => LPDENM
       LEDENX  => LEDENM

       LPXEL   => LPMEL
       LPXAT   => LPMAT
       LPXML   => LPMML  ! same as LPXX
       LPXIO   => LPMIO
       LPXPHT  => LPMPHT
       LPXPL   => LPMPL

       LEXEL   => LEMEL
       LEXAT   => LEMAT
       LEXML   => LEMML  ! same as LEXX
       LEXIO   => LEMIO
       LEXPHT  => LEMPHT
       LEXPL   => LEMPL

       LVXDENX => LVXDENM
       LVYDENX => LVYDENM
       LVZDENX => LVZDENM

       LPGENX  => LPGENM
       LEGENX  => LEGENM
       LVGENX  => LVGENM

       LMXPL   => LMMPL

       LRXEL   => LRMEL

       LPXX    => LPMML
       LEXX    => LEMML
       LSCX    => NLSPCSCL_MOL

       IF (LPDENX)  PDENX  => PDENM(IMOL,:)
       IF (LEDENX)  EDENX  => EDENM(IMOL,:)

       IF (LPXEL)   PXEL   => PMEL(:)
       IF (LSCX) THEN
         IF (LPXAT) THEN
           NTS_PXATA = NTS_PA+1
           NTS_PXATE = NTS_AM
           PXAT   => PMAT(1:NATM*NMOLP,:)
         END IF
         IF (LPXML) THEN
           NTS_PXMLA = NTS_AM+1
           NTS_PXMLE = NTS_MM
           PXML   => PMML(1:NMOL*NMOLP,:)  ! same as PXX
         END IF
         IF (LPXIO) THEN
           NTS_PXIOA = NTS_MM+1
           NTS_PXIOE = NTS_IM
           PXIO   => PMIO(1:NION*NMOLP,:)
         END IF
         IF (LPXPHT) THEN
           NTS_PXPHA = NTS_IM+1
           NTS_PXPHE = NTS_PHM
           PXPHT  => PMPHT(1:NPHOT*NMOLP,:)
         END IF
         IF (LPXPL) THEN
           NTS_PXPLA = NTS_PHM+1
           NTS_PXPLE = NTS_PM
           PXPL   => PMPL(1:NPLS*NMOLP,:)
         END IF
       ELSE
         IF (LPXAT)  PXAT   => PMAT(1:NATMI,:)
         IF (LPXML)  PXML   => PMML(1:NMOLI,:)  ! for species i=imol:
                                                ! same as PXX
         IF (LPXIO)  PXIO   => PMIO(1:NIONI,:)
         IF (LPXPHT) PXPHT  => PMPHT(1:NPHOTI,:)
         IF (LPXPL)  PXPL   => PMPL(1:NPLSI,:)
       END IF

       IF (LEXEL)   EXEL   => EMEL(:)
       IF (LEXAT)   EXAT   => EMAT(:)
       IF (LEXML)   EXML   => EMML(:)   ! for species i=imol:
                                        ! same as EXX
       IF (LEXIO)   EXIO   => EMIO(:)
       IF (LEXPHT)  EXPHT  => EMPHT(:)
       IF (LEXPL)   EXPL   => EMPL(1:NPLSI,:)

       IF (LVXDENX) VXDENX => VXDENM(IMOL,:)
       IF (LVYDENX) VYDENX => VYDENM(IMOL,:)
       IF (LVZDENX) VZDENX => VZDENM(IMOL,:)

       IF (LPGENX)  PGENX  => PGENM(IMOL,:)
       IF (LEGENX)  EGENX  => EGENM(IMOL,:)
       IF (LVGENX)  VGENX  => VGENM(IMOL,:)

       IF (LMXPL)   MXPL   => MMPL(1:NPLSI,:)

       IF (LRXEL)   RXEL   => RMEL(1:NMOLI,:)

       IF (LSCX) THEN
         NDXX = NMOL
         NDXXA = NTS_AM+1
         NDXXE = NTS_MM
         IF (LPXX)  PXX    => PMML(1:NMOL*NMOLP,:)
       ELSE
         IF (LPXX)  PXX    => PMML(1:NMOLI,:)
       END IF
       IF (LEXX)    EXX    => EMML(:)

       LEX     => LEM

!      LGXCX => LGMCX
!      LGXEI => LGMEI
!      LGXEL => LGMEL
!      LGXPI => LGMPI
cdr now: allocatable, rather than pointer
!pb that means: here we do a copy
       LGXCX(0:NRCX,0:1) = LGMCX(IMOL,0:NRCX,0:1)
       LGXEI(0:NREI)     = LGMEI(IMOL,0:NREI)
       LGXEL(0:NREL,0:1) = LGMEL(IMOL,0:NREL,0:1)
       LGXPI(0:NRPI,0:1) = LGMPI(IMOL,0:NRPI,0:1)

       NXEII => NMEII(IMOL)
       NXCXI => NMCXI(IMOL)
       NXELI => NMELI(IMOL)
       NXPII => NMPII(IMOL)

       NXEIIM => NMEIIM(IMOL)
       NXCXIM => NMCXIM(IMOL)
       NXELIM => NMELIM(IMOL)
       NXPIIM => NMPIIM(IMOL)

       NPBGKX => NPBGKM(IMOL)
       NGENX  => NGENM(IMOL)

       IXSPZ = IMOL
       NMETOFF = NSPA

       RMASSX => RMASSM(IMOL)
       CNDYNX => CNDYNM(IMOL)
       CVRSSX => CVRSSM(IMOL)

       LOGMOL(IMOL,ISTRA)=.TRUE.
       LOGXSPZ => LOGMOL(IMOL,ISTRA)

      case(3)
!  test ions

       LPDENX  => LPDENI
       LEDENX  => LEDENI

       LPXEL   => LPIEL
       LPXAT   => LPIAT
       LPXML   => LPIML
       LPXIO   => LPIIO  ! same as LPXX
       LPXPHT  => LPIPHT
       LPXPL   => LPIPL

       LEXEL   => LEIEL
       LEXAT   => LEIAT
       LEXML   => LEIML
       LEXIO   => LEIIO  ! same as LEXX
       LEXPHT  => LEIPHT
       LEXPL   => LEIPL

       LVXDENX => LVXDENI
       LVYDENX => LVYDENI
       LVZDENX => LVZDENI

       LRXEL   => LRIEL

       LPGENX  => LPGENI
       LEGENX  => LEGENI
       LVGENX  => LVGENI

       LMXPL   => LMIPL

       LPXX    => LPIIO
       LEXX    => LEIIO
       LSCX    => NLSPCSCL_ION

       IF (LPDENX)  PDENX  => PDENI(IION,:)
       IF (LEDENX)  EDENX  => EDENI(IION,:)

       IF (LPXEL)   PXEL   => PIEL(:)
       IF (LSCX) THEN
         IF (LPXAT) THEN
           NTS_PXATA = NTS_PM+1
           NTS_PXATE = NTS_AI
           PXAT   => PIAT(1:NATM*NIONP,:)
         END IF
         IF (LPXML) THEN
           NTS_PXMLA = NTS_AI+1
           NTS_PXMLE = NTS_MI
           PXML   => PIML(1:NMOL*NIONP,:)
         END IF
         IF (LPXIO) THEN
           NTS_PXIOA = NTS_MI+1
           NTS_PXIOE = NTS_II
           PXIO   => PIIO(1:NION*NIONP,:)
         END IF
         IF (LPXPHT) THEN
           NTS_PXPHA = NTS_II+1
           NTS_PXPHE = NTS_PHI
           PXPHT  => PIPHT(1:NPHOT*NIONP,:)
         END IF
         IF (LPXPL) THEN
           NTS_PXPLA = NTS_PHI+1
           NTS_PXPLE = NTS_PI
           PXPL   => PIPL(1:NPLS*NIONP,:)
         END IF
       ELSE
         IF (LPXAT)  PXAT   => PIAT(1:NATMI,:)
         IF (LPXML)  PXML   => PIML(1:NMOLI,:)
         IF (LPXIO)  PXIO   => PIIO(1:NIONI,:)  ! for species i=iion:
                                                ! same as PXX
         IF (LPXPHT) PXPHT  => PIPHT(1:NPHOTI,:)
         IF (LPXPL)  PXPL   => PIPL(1:NPLSI,:)
       END IF

       IF (LEXEL)   EXEL   => EIEL(:)
       IF (LEXAT)   EXAT   => EIAT(:)
       IF (LEXML)   EXML   => EIML(:)
       IF (LEXIO)   EXIO   => EIIO(:)   ! for species i=iion:
                                        ! same as EXX
       IF (LEXPHT)  EXPHT  => EIPHT(:)
       IF (LEXPL)   EXPL   => EIPL(1:NPLSI,:)

       IF (LVXDENX) VXDENX => VXDENI(IION,:)
       IF (LVYDENX) VYDENX => VYDENI(IION,:)
       IF (LVZDENX) VZDENX => VZDENI(IION,:)

       IF (LRXEL)   RXEL   => RIEL(1:NIONI,:)

       IF (LPGENX)  PGENX  => PGENI(IION,:)
       IF (LEGENX)  EGENX  => EGENI(IION,:)
       IF (LVGENX)  VGENX  => VGENI(IION,:)

       IF (LMXPL)   MXPL   => MIPL(1:NPLSI,:)

       IF (LSCX) THEN
         NDXX = NION
         NDXXA = NTS_MI+1
         NDXXE = NTS_II
         IF (LPXX)  PXX    => PIIO(1:NION*NIONP,:)
       ELSE
         IF (LPXX)  PXX    => PIIO(1:NIONI,:)
       END IF
       IF (LEXX)    EXX    => EIIO(:)

       LEX     => LEIO

!      LGXCX => LGICX
!      LGXEI => LGIEI
!      LGXEL => LGIEL
!      LGXPI => LGIPI
cdr now: allocatable, rather than pointer
!pb that means: here we do a copy
       LGXCX(0:NRCX,0:1) = LGICX(IION,0:NRCX,0:1)
       LGXEI(0:NREI)     = LGIEI(IION,0:NREI)
       LGXEL(0:NREL,0:1) = LGIEL(IION,0:NREL,0:1)
       LGXPI(0:NRPI,0:1) = LGIPI(IION,0:NRPI,0:1)

       NXEII => NIEII(IION)
       NXCXI => NICXI(IION)
       NXELI => NIELI(IION)
       NXPII => NIPII(IION)

       NXEIIM => NIEIIM(IION)
       NXCXIM => NICXIM(IION)
       NXELIM => NIELIM(IION)
       NXPIIM => NIPIIM(IION)

       NPBGKX => NPBGKI(IION)
       NGENX  => NGENI(IION)

       IXSPZ = IION
       NMETOFF = NSPAM

       RMASSX => RMASSI(IION)
       CNDYNX => CNDYNI(IION)
       CVRSSX => CVRSSI(IION)

       LOGION(IION,ISTRA)=.TRUE.
       LOGXSPZ => LOGION(IION,ISTRA)

      case(4)
!  bulk ions: nothing to be done

      case default
         write (iunout,*) ' WRONG TYPE IN SWITCH_PARTINFO '
         write (iunout,*) ' ITYP = ',ITYP
      end select

      return

      end subroutine eirene_switch_partinfo

      subroutine eirene_output_partinfo
      implicit none
      integer :: istr, it, is

      call eirene_leer(2)

      call eirene_headng ('STATISTICS OF CPU TIME SPENT'//
     .                    ' IN FOLLOWING TRAJECTORIES ',57)

! sum over species

      do istr = 1, nstrai
        time_ar0(0,istr) = sum(time_ar0(1:,istr))
        time_ar1(0,istr) = sum(time_ar1(1:,istr))
        time_ar2(0,istr) = sum(time_ar2(1:,istr))
        time_ar3(0,istr) = sum(time_ar3(1:,istr))
      end do

! sum over strata

      do is = 0, ubound(time_ar0,1)
         time_ar0(is,0) = sum(time_ar0(is,1:))
        end do
      do is = 0, ubound(time_ar1,1)
         time_ar1(is,0) = sum(time_ar1(is,1:))
      end do
      do is = 0, ubound(time_ar2,1)
         time_ar2(is,0) = sum(time_ar2(is,1:))
      end do
      do is = 0, ubound(time_ar3,1)
         time_ar3(is,0) = sum(time_ar3(is,1:))
      end do

      do istr = 0, nstrai

        call eirene_leer(1)
        if (istr == 0) then
          write (iunout,'(1x,A,I6)') 'SUM OVER STRATA '
          write (iunout,'(1x,A,I6)') '=============== '
        else
          write (iunout,'(1x,A,I6)') 'STRATUM ',istr
          write (iunout,'(1x,A,I6)') '============== '
        end if

        if (any(time_ar0(:,istr) > 0._dp)) then
          call eirene_leer(1)
          write (iunout,*) 'TIME (SEC) SPENT IN FOLLOWING '
          CALL EIRENE_MASYR1 ('PHOTONS = ',time_ar0(0:nphot,:),
     .                LOGPHOT,ISTR,0,NPHOT,NPHOT,0,NSTRA,TEXTS(1))
          CALL EIRENE_MASAGE
     .      ('SUM OVER SPECIES                               ')
          CALL EIRENE_MASR1 ('TOTAL=  ',time_ar0(0,ISTR))
        end if

        if (any(time_ar1(:,istr) > 0._dp)) then
          call eirene_leer(1)
          write (iunout,*) 'TIME (SEC) SPENT IN FOLLOWING '
          CALL EIRENE_MASYR1 ('ATOMS =   ',time_ar1(0:natm,:),
     .                LOGATM,ISTR,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
         CALL EIRENE_MASAGE
     .      ('SUM OVER SPECIES                               ')
          CALL EIRENE_MASR1 ('TOTAL=  ',time_ar1(0,ISTR))
        end if

        if (any(time_ar2(:,istr) > 0._dp)) then
          call eirene_leer(1)
          write (iunout,*) 'TIME (SEC) SPENT IN FOLLOWING '
          CALL EIRENE_MASYR1 ('MOLECULES=',time_ar2(0:nmol,:),
     .                LOGMOL,ISTR,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
          CALL EIRENE_MASAGE
     .      ('SUM OVER SPECIES                               ')
          CALL EIRENE_MASR1 ('TOTAL=  ',time_ar2(0,ISTR))
        end if

        if (any(time_ar3(:,istr) > 0._dp)) then
          call eirene_leer(1)
          write (iunout,*) 'TIME (SEC) SPENT IN FOLLOWING '
          CALL EIRENE_MASYR1 ('TEST IONS=',time_ar3(0:nion,:),
     .                LOGION,ISTR,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
          CALL EIRENE_MASAGE
     .      ('SUM OVER SPECIES                               ')
          CALL EIRENE_MASR1 ('TOTAL=  ',time_ar3(0,ISTR))
        end if

!        if (any(time_ar4(:,istr) > 0._dp)) then
!          write (iunout,*) ' TIME (SEC) SPENT IN FOLLOWING '
!          CALL EIRENE_MASYR1 ('BULK IONS=',time_ar4(1:npls,:),
!     .                LOGPLS,ISTR,1,NPLS,NPLS,1,NSTRA,TEXTS(NSPAMI+1))
!          CALL EIRENE_MASAGE
!     .      ('SUM OVER SPECIES                               ')
!          CALL EIRENE_MASR1 ('TOTAL=  ',SUM(time_ar4(1:npls,ISTR)))
!        end if

      end do  !istra
      call eirene_leer(1)
      return
      end subroutine eirene_output_partinfo

      subroutine eirene_reinit_partinfo
      implicit none

      if (allocated(time_ar0)) then
        deallocate(time_ar0)
        deallocate(time_ar1)
        deallocate(time_ar2)
        deallocate(time_ar3)
      endif
      istra_old=-1
      ityp_old=-1
      iphot_old=-1
      iatm_old=-1
      imol_old=-1
      iion_old=-1
      ipls_old=-1

      return
      end subroutine eirene_reinit_partinfo

      end module eirmod_switch_partinfo
