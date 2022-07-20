

      module eirmod_extrab25

      use eirmod_precision
      use eirmod_parmmod
      use eirmod_COMXS
      use eirmod_COMUSR
      use eirmod_CTEXT
      use eirmod_COMSOU
      use eirmod_CGRID
      use eirmod_CESTIM
      use eirmod_CTRCEI
      use eirmod_CCONA
      use eirmod_COUTAU
      use eirmod_CPOLYG
      use eirmod_CSDVI
      use eirmod_CGEOM
      use eirmod_CSPEI
      use eirmod_CADGEO
      use eirmod_CLGIN
      use eirmod_COMPRT
      use eirmod_CLOGAU
      use eirmod_COMSIG
      use eirmod_CUPD
      use eirmod_CZT1
      use eirmod_CCOUPL
      use eirmod_comspl
      use eirmod_ctrig

      implicit none
      private

      public :: eirene_extrab25_cleanup,  eirene_extrab25_wneutrals
      public :: eirene_extrab25_wneuinit, eirene_extrab25_wneufill
      public :: eirene_extrab25_wneusave, eirene_extrab25_wneuclean
      public :: eirene_extrab25_alloc_mods
      public :: eirene_extrab25_iniusr_init
      public :: eirene_extrab25_emissivity
      public :: b2_cell

      ! eirdiag.h/eirdiag.f
      !c*** Volume data:
      !c***    srcml   :   power loss due to molecules, including
      !c***    edissml :   power loss due to molecule dissociation
      !c***    eneutrad:   power radiated due to neutral atoms
      !c***    emolrad :   power radiated due to molecules
      !c***    eionrad :   power radiated due to molecular ions
      !c*** Surface data:
      !c***    wldnek  :   heat transferred with neutrals
      !c***    wldnep  :   potential energy released by neutrals
      !c***    wldna   :   flux of atoms impinging onto the surface
      !c***    ewlda   :   their average energy
      !c***    wldnm   :   flux of molecules impinging onto the surface
      !c***    ewldm   :   their average energy
      !c***    wldra   :   flux of reflected atoms
      !c***    wldrm   :   flux of reflected molecules
      !c***    wldpp   :   flux of plasma ions impinging onto the surface
      !c***    wldpa   :   flux of resulting atoms
      !c***    wldpm   :   flux of resulting molecules
      !c***    wldpeb  :   power carried away by these atoms and molecules
      !c***    wldspt  :   flux of sputtered wall material
      !c***    isrftype:   surface type (iliin in Eirene)
      !c***    wlarea  :   areas of the surface segments from Eirene
      !c***    wlabsrp :   absorption at the surfaces (1-recyct from Eirene)
      !c***    wlpump  :   pumped flux at the surfaces
      real(DP), save, allocatable, dimension(:,:,:,:), public ::
     .  dab2,dmb2,dib2,tab2,tmb2,tib2,rfluxa,rfluxm,refluxa,refluxm,
     .  pfluxa,pfluxm,pefluxa,pefluxm,emiss,emissmol,srcml,edissml
      real(DP), save, allocatable, dimension(:,:),public ::
     .  wldnek,wldnep
      real(DP), save, allocatable, dimension (:,:,:), public ::
     .  wldna,ewlda,wldnm,ewldm,wldra,wldrm,wldpp,wldpa,wldpm
      real(DP), save, allocatable, dimension (:,:),public ::
     .  wldpeb,wldspt
      real(DP), save, allocatable, dimension (:,:,:,:), public ::
     .  eneutrad, emolrad, eionrad

      integer, save, public :: nnlimi,nnstsi,nnatmi,nnmoli,nnioni
      integer, save, public :: nnplsi,nns
      integer, save, allocatable, public :: isrftype(:)
      logical, save, public :: lhalpha=.false.,lvib=.false.

      ! wneutral globals
!c      logical, save :: hlp_pr
      integer, save :: ia0,ia1,ia2,ia3,ifirst_wneutral=0
      real(DP), save :: hlp_cnv

      ! B2.5 neutrals parameters modifications
      integer, save, public :: bn_spcsrf
      integer, save, public, allocatable :: bl_spcsrf(:), bi_spcsrf(:)
      integer, save, public, allocatable :: bj_spcsrf(:), bsps_sgrp(:)
      real(DP), save, public, allocatable :: bsps_absr(:), bsps_trno(:)
      real(DP), save, public, allocatable :: bsps_mtri(:), bsps_tmpr(:)
      real(DP), save, public, allocatable :: bsps_trni(:), bsps_spph(:)
      real(DP), save, public, allocatable :: bsps_spch(:)
      character*8, save, public, allocatable :: bsps_mtrl(:), bsps_id(:)

      ! diag2 globals and parameters
      integer, parameter :: mgwtiesx=12 ! max. number of wall segments tied to a grid edge segment

      ! remaining bits and pieces from braeir common
      real(DP), save, public :: chemical_sputter_yield, fchar_chemical
      integer, save, public :: igass_chemical, itsput_chemical,
     .                          issput_chemical

      !pb volumes of B2.5 cells
      real(DP), save, public, allocatable :: volcel(:,:)

      ! normals of B2.5 cell edges per triangle
      real(DP), save, public, allocatable :: plnxtri(:), plnytri(:)
      real(DP), save, public, allocatable :: pplnxtri(:), pplnytri(:)

      !flux_save
      real(DP), save, public, allocatable :: flux_save(:)

!pb 27012016
! flag indicating if subroutine iniusr is called from B2.5
      integer, public, save :: ini_iniusr=0

      contains

      subroutine eirene_extrab25_alloc_mods(nnx,nny)
      use eirmod_parmmod
      use eirmod_eirbra
      use eirmod_braeir
      implicit none
      integer, intent(in) :: nnx,nny
      call eirene_find_param
      call eirene_set_parmmod(1)
      call eirene_alloc_braeir(nnx,nny,nfl)
      call eirene_alloc_eirbra(nnx,nny,nfl,nstra)
      allocate (flux_save(nstra))
      flux_save=0.d0
      end subroutine eirene_extrab25_alloc_mods

      subroutine eirene_extrab25_wneutrals
      ! old version : 27.07.2000 23:07
      ! new version : 25.01.2011 (s.wiesen@fz-juelich.de) - f90 module
      !c     this subroutine produces an output file ft44 for plotting in B2
      !c     with the neutral densities, temperatures, and fluxes
      !c     (summed up for all strata)
      implicit none

      real(DP) :: dummy(0:ndxp,0:ndyp)
      !c*** label for fort.44 file
      integer, parameter  :: jvft44=20000727
      !c*** and dissociation energy of the hydrogen molecule
      real(DP), parameter  :: diss_pot_H2=4.48

      integer :: ix,iy,ir,iistra,icell,ierror,i,j,in,jatm,jmol,jion
      integer :: l,k,nred
      real(DP) :: de,te,hlp
      real(DP) :: dej,tei,tef,def
      real(DP) :: powalf,datm3,dpls3,dmol3,dion3,dnml3,sigadd
      real(DP) :: da,dpp,dm,di,dn,ratio2,ratio7
      integer :: istra_in,istra_save
      real(DP) :: vl
      real(DP) :: value
      character*36 hlp_frm
      external eirene_indmpi,eirene_neutr

!     !c======================================================================
      !c---------------------------------------------------------------------<
      entry eirene_extrab25_wneuinit
      !c      print *,'%%% wneuinit'
      !c--------------------------------------------------------------------->

      IF (IFIRST_wneutral.NE.0) return
      ifirst_wneutral=1

      allocate(dab2(0:ndxp,0:ndyp,natm,1))
      allocate(dmb2(0:ndxp,0:ndyp,nmol,1))
      allocate(dib2(0:ndxp,0:ndyp,nion,1))
      allocate(tab2(0:ndxp,0:ndyp,natm,1))
      allocate(tmb2(0:ndxp,0:ndyp,nmol,1))
      allocate(tib2(0:ndxp,0:ndyp,nion,1))
      allocate(rfluxa(0:ndxp,0:ndyp,natm,1))
      allocate(rfluxm(0:ndxp,0:ndyp,nmol,1))
      allocate(refluxa(0:ndxp,0:ndyp,natm,1))
      allocate(refluxm(0:ndxp,0:ndyp,nmol,1))
      allocate(pfluxa(0:ndxp,0:ndyp,natm,1))
      allocate(pfluxm(0:ndxp,0:ndyp,nmol,1))
      allocate(pefluxa(0:ndxp,0:ndyp,natm,1))
      allocate(pefluxm(0:ndxp,0:ndyp,nmol,1))
      allocate(emiss(0:ndxp,0:ndyp,1,1))
      allocate(emissmol(0:ndxp,0:ndyp,1,1))
      allocate(srcml(0:ndxp,0:ndyp,nmol,1))
      allocate(edissml(0:ndxp,0:ndyp,nmol,0:nstra+1))
      allocate(wldnek(nlmpgs,0:nstra+1))
      allocate(wldnep(nlmpgs,0:nstra+1))
      allocate(wldna(nlmpgs,natm,0:nstra+1))
      allocate(ewlda(nlmpgs,natm,0:nstra+1))
      allocate(wldnm(nlmpgs,nmol,0:nstra+1))
      allocate(ewldm(nlmpgs,nmol,0:nstra+1))
      allocate(wldra(nlmpgs,natm,0:nstra+1))
      allocate(wldrm(nlmpgs,nmol,0:nstra+1))
      allocate(wldpp(nlmpgs,npls,0:nstra+1))
      allocate(wldpa(nlmpgs,natm,0:nstra+1))
      allocate(wldpm(nlmpgs,nmol,0:nstra+1))
      allocate(wldpeb(nlmpgs,0:nstra+1))
      allocate(wldspt(nlmpgs,0:nstra+1))
      allocate(eneutrad(0:ndxp,0:ndyp,natm,0:nstra+1))
      allocate(emolrad(0:ndxp,0:ndyp,nmol,0:nstra+1))
      allocate(eionrad(0:ndxp,0:ndyp,nion,0:nstra+1))
      allocate(isrftype(nlmpgs))

      allocate(volcel(0:ndxp,0:ndyp))

      !tamas zero init
      wldnep = 0
      wldna  = 0
      ewlda  = 0
      wldnm  = 0
      ewldm  = 0
      wldra  = 0
      wldrm  = 0
      wldpp  = 0
      wldpa  = 0
      wldpm  = 0
!pb   wlarea = 0
!pb   wldspta = 0
!pb   wldsptm = 0

      call eirene_extraB25_wneuclean

      hlp_cnv=1./elcha
      ia0=0
      ia1=ia0+natmi+nmoli
      ia2=ia1+natmi+nmoli
      ia3=ia2+natmi+nmoli

      !c======================================================================
      !c*** fill arrays
      !c
      !c---------------------------------------------------------------------<
      return
      !c
      !c======================================================================
      entry eirene_extraB25_wneufill(istra_in)
      !c
      !c write (iunout,*),'%%% wneufill: istra,istra_in = ',istra,istra_in
      istra_save=istra
      istra=istra_in
      !c--------------------------------------------------------------------->

      !csw
      !csw 21feb2012 corrected radiation from neutrals (atoms only), taken from SOLPS4.3 (V.Kotov)
      !csw 04mar2013 shifted from wneusave to here (wneufill)
      !csw
      value=0.0
      eneutrad(:,:,:,istra) = 0.0_dp
      emolrad(:,:,:,istra) = 0.0_dp
      eionrad(:,:,:,istra) = 0.0_dp
      do icell=1,ntrii
        ix=ixtri(icell)
        iy=iytri(icell)
        if(b2_cell(ix,iy)) then
          do jatm=1,natmi
            if (lrael) eneutrad(ix,iy,jatm,istra)=
     .        eneutrad(ix,iy,jatm,istra)+rael(jatm,icell)*vol(icell)
          enddo
          do jmol=1,nmoli
            if (lrmel) emolrad(ix,iy,jmol,istra)=
     .        emolrad(ix,iy,jmol,istra)+rmel(jmol,icell)*vol(icell)
          enddo
          do jion=1,nioni
            if (lriel) eionrad(ix,iy,jion,istra)=
     .        eionrad(ix,iy,jion,istra)+riel(jion,icell)*vol(icell)
          enddo
        endif
      enddo

      volcel = 0.d0
      do in=1,ntrii
        ix=ixtri(in)
        iy=iytri(in)
        if(b2_cell(ix,iy)) volcel(ix,iy) = volcel(ix,iy) + vol(in)
      end do

      !c
      !C map 1d EIRENE neutral densities and temperatures on 2d arrays
      !c and change the units to SI for plotting in B2
      !c
      !write(iunout,*) 'istra ',istra
      !write(iunout,*) 'natmi, nmoli, nioni ',natmi,nmoli,nioni
      !crfs     IF (WTOTP(0,ISTRA).EQ.0.) GOTO 60
      do jatm=1,natmi
        do in=1,ntrii
          ix=ixtri(in)
          iy=iytri(in)
          if(b2_cell(ix,iy)) then
            vl = vol(in)/volcel(ix,iy)
            if (lpdena) dab2(ix,iy,jatm,1)=dab2(ix,iy,jatm,1)+
     .                          pdena(jatm,in)*1e6*vl
            if (laddv) then
             if(ia0+jatm.le.nadv)
     .         rfluxa(ix,iy,jatm,1)=rfluxa(ix,iy,jatm,1)+
     .                                   addv(ia0+jatm,in)*1.0e4
             if(ia2+jatm.le.nadv)
     .         pfluxa(ix,iy,jatm,1)=pfluxa(ix,iy,jatm,1)+
     .                                   addv(ia2+jatm,in)*1.0e4
             if(ia1+jatm.le.nadv)
     .         refluxa(ix,iy,jatm,1)=refluxa(ix,iy,jatm,1)+
     .                             addv(ia1+jatm,in)*1.0e4*elcha
             if(ia3+jatm.le.nadv)
     .         pefluxa(ix,iy,jatm,1)=pefluxa(ix,iy,jatm,1)+
     .                             addv(ia3+jatm,in)*1.0e4*elcha
             end if
            if (ledena) tab2(ix,iy,jatm,1)=tab2(ix,iy,jatm,1)+
     .                             edena(jatm,in)*vl
          endif
        end do
      end do
      do jmol=1,nmoli
        edissml(:,:,jmol,istra) = 0.d0
        do in=1,ntrii
          ix=ixtri(in)
          iy=iytri(in)
          if(b2_cell(ix,iy)) then
            vl = vol(in)/volcel(ix,iy)
            if (lpdenm) dmb2(ix,iy,jmol,1)=dmb2(ix,iy,jmol,1)+
     .                          pdenm(jmol,in)*1e6*vl
            if(ia0+natmi+jmol.le.nadv)
     .         rfluxm(ix,iy,jmol,1)=rfluxm(ix,iy,jmol,1)+
     .                                 addv(ia0+natmi+jmol,in)*1.0e4
            if(ia2+natmi+jmol.le.nadv)
     .         pfluxm(ix,iy,jmol,1)=pfluxm(ix,iy,jmol,1)+
     .                             addv(ia2+natmi+jmol,in)*1.0e4

            if(ia1+natmi+jmol.le.nadv .and. laddv)
     .         refluxm(ix,iy,jmol,1)=refluxm(ix,iy,jmol,1)+
     .                       addv(ia1+natmi+jmol,in)*1.0e4*elcha
            if(ia3+natmi+jmol.le.nadv .and. laddv)
     .         pefluxm(ix,iy,jmol,1)=pefluxm(ix,iy,jmol,1)+
     .                       addv(ia3+natmi+jmol,in)*1.0e4*elcha
            if (ledenm) tmb2(ix,iy,jmol,1)=tmb2(ix,iy,jmol,1)+
     .                      edenm(jmol,in)*vl
            if (lpmml) srcml(ix,iy,jmol,1)=srcml(ix,iy,jmol,1)+
     .                                     pmml(jmol,in)*vl
            !c*** Potential energy source related to hydrogen molecules
            !c*** for determining the radiation (W per cell)
            if (ncharm(jmol).eq.2 .and. lpmml)
     .         edissml(ix,iy,jmol,istra)=edissml(ix,iy,jmol,istra)+
     .            pmml(jmol,in)*vl*diss_pot_H2
          endif
        end do
      end do
      do jion=1,nioni
        do in=1,ntrii
          ix=ixtri(in)
          iy=iytri(in)
          if(b2_cell(ix,iy)) then
            vl = vol(in)/volcel(ix,iy)
            if (lpdeni) dib2(ix,iy,jion,1)=dib2(ix,iy,jion,1)+
     .                          pdeni(jion,in)*1e6*vl
            if (ledeni) tib2(ix,iy,jion,1)=tib2(ix,iy,jion,1)+
     .                          edeni(jion,in)*vl
          endif
        end do
      end do

      !c
      !c*** Rescale the surface data from A to 1/sec and average the energy
      !c
      do i=1,nlimps
        wldnek(i,istra)=0.
        wldnep(i,istra)=0.
        wldpeb(i,istra)=0.
        !c*** hlp accumulates the power taken away with reemitted particles
        hlp=0.
        do j=1,natmi
          if (leotat)  wldnek(i,istra)=wldnek(i,istra)+eotat(j,i)
          if (lerfpat) wldpeb(i,istra)=wldpeb(i,istra)+erfpat(j,i)
          if (lerfaat) hlp=hlp+erfaat(j,i)
          if (lerfmat) hlp=hlp+erfmat(j,i)
          if (lerfiat) hlp=hlp+erfiat(j,i)
          if (lpotat.and.leotat) then
            if(potat(j,i).gt.0.) then
              ewlda(i,j,istra)=eotat(j,i)/potat(j,i)
            else
              ewlda(i,j,istra)=0.
            end if
          else
            ewlda(i,j,istra)=0.
          end if
          if (lpotat) wldna(i,j,istra)=hlp_cnv*potat(j,i)
          wldra(i,j,istra)=0.
          if (lprfaat) wldra(i,j,istra)=wldra(i,j,istra)+prfaat(j,i)
          if (lprfmat) wldra(i,j,istra)=wldra(i,j,istra)+prfmat(j,i)
          if (lprfiat) wldra(i,j,istra)=wldra(i,j,istra)+prfiat(j,i)
          wldra(i,j,istra)=hlp_cnv*wldra(i,j,istra)
          if (lprfpat) wldpa(i,j,istra)=hlp_cnv*prfpat(j,i)
        end do
        do j=1,nmoli
          if (leotml)  wldnek(i,istra)=wldnek(i,istra)+eotml(j,i)
          if (lerfpml) wldpeb(i,istra)=wldpeb(i,istra)+erfpml(j,i)
          if (lerfaml) hlp=hlp+erfaml(j,i)
          if (lerfmml) hlp=hlp+erfmml(j,i)
          if (lerfiml) hlp=hlp+erfiml(j,i)
          if (lpotml.and.leotml) then
            if(potml(j,i).gt.0.) then
              ewldm(i,j,istra)=eotml(j,i)/potml(j,i)
            else
              ewldm(i,j,istra)=0.
            end if
          else
            ewldm(i,j,istra)=0.
          end if
          if (lpotml) wldnm(i,j,istra)=hlp_cnv*potml(j,i)
          wldrm(i,j,istra)=0.
          if (lprfaml) wldrm(i,j,istra)=wldrm(i,j,istra)+prfaml(j,i)
          if (lprfmml) wldrm(i,j,istra)=wldrm(i,j,istra)+prfmml(j,i)
          if (lprfiml) wldrm(i,j,istra)=wldrm(i,j,istra)+prfiml(j,i)
          wldrm(i,j,istra)=hlp_cnv*wldrm(i,j,istra)
          if (lprfpml) wldpm(i,j,istra)=hlp_cnv*prfpml(j,i)
        end do
        do j=1,nioni
          if (leotio)  wldnek(i,istra)=wldnek(i,istra)+eotio(j,i)
          if (lerfpio) wldpeb(i,istra)=wldpeb(i,istra)+erfpio(j,i)
          if (lerfaio) hlp=hlp+erfaio(j,i)
          if (lerfmio) hlp=hlp+erfmio(j,i)
          if (lerfiio) hlp=hlp+erfiio(j,i)
        end do
        do j=1,nfla
          if (lpotpl) wldpp(i,j,istra)=hlp_cnv*potpl(j,i)
        end do
        !c*** recombination energy of hydrogen molecules
        do j=1,nmoli
          if(ncharm(j).eq.2 .and. lprfaml)
     .      wldnep(i,istra)=wldnep(i,istra)+ diss_pot_H2*prfaml(j,i)
        end do
        !c*** subtract the outcoming power
        wldnek(i,istra)=wldnek(i,istra)-hlp

        !c*** calculate the flux of sputtered particles

        !CSW CHECK 26jan2011 NOT READY, sptwll missing!
        !wldspt(i,istra)=hlp_cnv*sptwll(i)
        wldspt(i,istra)=0.
      end do
      !c

      !c---------------------------------------------------------------------<
      istra=istra_save
      return
      !c
      entry eirene_extraB25_wneusave
      !c
      !c      write (iunout,*) '%%% wneusave: istra = ',istra
      !c*** Calculate the totals (stratum 0)
      !c
      do i=1,nlimps
        wldnek(i,0)=0.0_DP
        wldnep(i,0)=0.0_DP
        wldpeb(i,0)=0.0_DP
        wldspt(i,0)=0.0_DP
        do j=1,natmi
          ewlda(i,j,0)=0.0_DP
          wldna(i,j,0)=0.0_DP
          wldra(i,j,0)=0.0_DP
          wldpa(i,j,0)=0.0_DP
        end do
        do j=1,nmoli
          ewldm(i,j,0)=0.0_DP
          wldnm(i,j,0)=0.0_DP
          wldrm(i,j,0)=0.0_DP
          wldpm(i,j,0)=0.0_DP
        end do
        do j=1,nfla
          wldpp(i,j,0)=0.0_DP
        end do
        !c
        do k=1,nstrai
          wldnek(i,0)=wldnek(i,0)+wldnek(i,k)
          wldnep(i,0)=wldnep(i,0)+wldnep(i,k)
          wldpeb(i,0)=wldpeb(i,0)+wldpeb(i,k)
          wldspt(i,0)=wldspt(i,0)+wldspt(i,k)
          do j=1,natmi
            wldna(i,j,0)=wldna(i,j,0)+wldna(i,j,k)
            wldra(i,j,0)=wldra(i,j,0)+wldra(i,j,k)
            wldpa(i,j,0)=wldpa(i,j,0)+wldpa(i,j,k)
            ewlda(i,j,0)=ewlda(i,j,0)+ewlda(i,j,k)*wldna(i,j,k)
          end do
          do j=1,nmoli
            wldnm(i,j,0)=wldnm(i,j,0)+wldnm(i,j,k)
            wldrm(i,j,0)=wldrm(i,j,0)+wldrm(i,j,k)
            wldpm(i,j,0)=wldpm(i,j,0)+wldpm(i,j,k)
            ewldm(i,j,0)=ewldm(i,j,0)+ewldm(i,j,k)*wldnm(i,j,k)
          end do
          do j=1,nfla
            wldpp(i,j,0)=wldpp(i,j,0)+wldpp(i,j,k)
          end do
        end do
        do j=1,natmi
          if(wldna(i,j,0).gt.0.0_DP) then
            ewlda(i,j,0)=ewlda(i,j,0)/wldna(i,j,0)
          else
            ewlda(i,j,0)=0.0_DP
          end if
        end do
        do j=1,nmoli
          if(wldnm(i,j,0).gt.0.0_DP) then
            ewldm(i,j,0)=ewldm(i,j,0)/wldnm(i,j,0)
          else
            ewldm(i,j,0)=0.0_DP
          end if
        end do
      end do
!pb
      eneutrad(:,:,:,0) = 0.0_DP
      emolrad(:,:,:,0) = 0.0_DP
      eionrad(:,:,:,0) = 0.0_DP
      edissml(:,:,:,0) = 0.0_DP
      do ix = 1, ndxa
        do iy = 1, ndya
          do jatm=1,natmi
            eneutrad(ix,iy,jatm,0) = sum(eneutrad(ix,iy,jatm,1:nstrai))
          end do
          do jmol=1,nmoli
            edissml(ix,iy,jmol,0) = sum(edissml(ix,iy,jmol,1:nstrai))
            emolrad(ix,iy,jmol,0) = sum(emolrad(ix,iy,jmol,1:nstrai))
          end do
          do jion=1,nioni
            eionrad(ix,iy,jion,0) = sum(eionrad(ix,iy,jion,1:nstrai))
          end do
        end do
      end do
      !c--------------------------------------------------------------------->


      !c*** Surface type and properties
      do i=1,nlimps
        isrftype(i)=iliin(i)
      end do

      !c
      !c*** Calculate the temperatures
      !c
      do ix=1,ndxa
        do iy=1,ndya
          do jatm=1,natmi
            if (dab2(ix,iy,jatm,1).gt.0.) then
              tab2(ix,iy,jatm,1)=tab2(ix,iy,jatm,1)/dab2(ix,iy,jatm,1)*
     .                                                  elcha*2./3.*1.e6
            else
              tab2(ix,iy,jatm,1)=1.e-6*elcha
            end if
          end do
          do jmol=1,nmoli
            if (dmb2(ix,iy,jmol,1).gt.0.) then
              tmb2(ix,iy,jmol,1)=tmb2(ix,iy,jmol,1)/dmb2(ix,iy,jmol,1)*
     .                                                  elcha*2./3.*1.e6
            else
              tmb2(ix,iy,jmol,1)=1.e-6*elcha
            end if
          end do
          do jion=1,nioni
            if (dib2(ix,iy,jion,1).gt.0.) then
              tib2(ix,iy,jion,1)=tib2(ix,iy,jion,1)/dib2(ix,iy,jion,1)*
     .                                                  elcha*2./3.*1.e6
            else
              tib2(ix,iy,jion,1)=1.e-6*elcha
            end if
          end do
        end do
      end do

      !c*** print some neutral fluxes across the "non-default" surfaces

      write(hlp_frm,'(a,i3,a)')
     . '(/16x,2(2x,a8,1x),',3*(natmi+nmoli)+nfla,'(3x,a6,i2))'
      write(iunout,hlp_frm) '  area','  power',
     .             ('atflx',i,i=1,natmi),('mlflx',i,i=1,nmoli),
     .             ('atflxr',i,i=1,natmi),('mlflxr',i,i=1,nmoli),
     .             ('atflxp',i,i=1,natmi),('mlflxp',i,i=1,nmoli),
     .             ('plsflx',i,i=1,nfla)
      do i=1,nstsi
        j=nlim+i
        do k=1,nfla
          wldpp(j,k,0)=0.
          do l=1,nstrai
            wldpp(j,k,0)=wldpp(j,k,0)+wldpp(j,k,l)
          end do
        end do
        do k=1,natmi
          wldpa(j,k,0)=0.
          do l=1,nstrai
            wldpa(j,k,0)=wldpa(j,k,0)+wldpa(j,k,l)
          end do
        end do
        do k=1,nmoli
          wldpm(j,k,0)=0.
          do l=1,nstrai
            wldpm(j,k,0)=wldpm(j,k,0)+wldpm(j,k,l)
          end do
        end do

        write(hlp_frm,'(a,i3,a)')
     .   '(a,i4,1p,',2+3*(natmi+nmoli)+nfla,'e11.3)'
        if (nmoli.gt.0) then
          if(wldnek(j,0).ne.0. .or.
     .       wldna(j,1,0).ne.0 .or. wldnm(j,1,0).ne.0 .or.
     .       wldra(j,1,0).ne.0 .or. wldrm(j,1,0).ne.0 .or.
     .       wldpa(j,1,0).ne.0 .or. wldpm(j,1,0).ne.0 .or.
     .       wldpp(j,1,0).ne.0)
     .         write(iunout,hlp_frm)
     .           'non-def-surf ',i, 1.e-4*sarea(j), 1.e-6*wldnek(j,0),
     .               (wldna(j,k,0),k=1,natmi),(wldnm(j,k,0),k=1,nmoli),
     .               (wldra(j,k,0),k=1,natmi),(wldrm(j,k,0),k=1,nmoli),
     .               (wldpa(j,k,0),k=1,natmi),(wldpm(j,k,0),k=1,nmoli),
     .               (wldpp(j,k,0),k=1,nfla)
        else
          if(wldnek(j,0).ne.0. .or.
     .       wldna(j,1,0).ne.0 .or. wldra(j,1,0).ne.0 .or.
     .       wldpa(j,1,0).ne.0 .or. wldpp(j,1,0).ne.0)
     .         write(iunout,hlp_frm)
     .           'non-def-surf ',i, 1.e-4*sarea(j), 1.e-6*wldnek(j,0),
     .               (wldna(j,k,0),k=1,natmi),(wldra(j,k,0),k=1,natmi),
     .               (wldpa(j,k,0),k=1,natmi),(wldpp(j,k,0),k=1,nfla)
        end if
      end do
      call eirene_leer(1)
!cc%%%
!c      write (iunout,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!c      write (iunout,'(/6x,20(a8,2x))')
!c     ,         'wldnek','wldnep','wldna','ewlda','wldnm','ewldm',
!c     ,         'wldna He','ewlda He','wldna Ne','ewlda Ne','wldra H',
!c     ,         'wldra He','wldra Ne','wldrm','prfaat'
!c      do i=1,nlimi
!c      write (iunout,'(1p,i6,20e10.2)') i,wldnek(i),wldnep(i),
!c     ,       wldna(i,1),ewlda(i,1),wldnm(i,1),ewldm(i,1),
!c     ,       wldna(i,2),ewlda(i,2),wldna(i,3),ewlda(i,3),
!c     ,       wldra(i,1),wldra(i,2),wldra(i,3),wldrm(i,1),
!c     ,       prfaat(1,i)
!c      end do
!c      write (iunout,'(/6x,20(a8,2x))')
!c     ,         'wldnek','wldnep','wldna','ewlda','wldnm','ewldm',
!c     ,         'wldna He','ewlda He','wldna Ne','ewlda Ne','wldra H',
!c     ,         'wldra He','wldra Ne','wldrm','prfaat'
!c      do i=nlim+1,nlim+nstsi
!c      write (iunout,'(1p,i6,20e10.2)') i-nlim,wldnek(i),wldnep(i),
!c     ,       wldna(i,1),ewlda(i,1),wldnm(i,1),ewldm(i,1),
!c     ,       wldna(i,2),ewlda(i,2),wldna(i,3),ewlda(i,3),
!c     ,       wldra(i,1),wldra(i,2),wldra(i,3),wldrm(i,1),
!c     ,       prfaat(1,i)
!c      end do
!c      write (iunout,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!cc%%%
      write(iunout,*) 'ncutl,ncutb ',ncutl,ncutb
      write(iunout,'(1x,a,7i6)') 'ndx,ndy,natm,ndxa,ndya,nfla,n1st',
     .                            ndx,ndy,natm,ndxa,ndya,nfla,n1st
      !c
      !c*** backmapping of 2d arrays for b2
      !c
      if (ncutl.ne.ncutb) then
        call eirene_indmpi(dab2,dummy,ndx,ndy,natm,ndxa,ndya,natmi,
     .                                    ncutb,ncutl,npoint,npplg,1,1)
        call eirene_indmpi(rfluxa,dummy,ndx,ndy,natm,ndxa,ndya,natmi,
     .                                    ncutb,ncutl,npoint,npplg,1,1)
        call eirene_indmpi(refluxa,dummy,ndx,ndy,natm,ndxa,ndya,natmi,
     .                                    ncutb,ncutl,npoint,npplg,1,1)
        call eirene_indmpi(pfluxa,dummy,ndx,ndy,natm,ndxa,ndya,natmi,
     .                                    ncutb,ncutl,npoint,npplg,1,1)
        call eirene_indmpi(pefluxa,dummy,ndx,ndy,natm,ndxa,ndya,natmi,
     .                                    ncutb,ncutl,npoint,npplg,1,1)
        call eirene_indmpi(tab2,dummy,ndx,ndy,natm,ndxa,ndya,natmi,
     .                                    ncutb,ncutl,npoint,npplg,1,1)
        call eirene_indmpi(dmb2,dummy,ndx,ndy,nmol,ndxa,ndya,nmoli,
     .                                    ncutb,ncutl,npoint,npplg,1,1)
        call eirene_indmpi(rfluxm,dummy,ndx,ndy,nmol,ndxa,ndya,nmoli,
     .                                    ncutb,ncutl,npoint,npplg,1,1)
        call eirene_indmpi(refluxm,dummy,ndx,ndy,nmol,ndxa,ndya,nmoli,
     .                                    ncutb,ncutl,npoint,npplg,1,1)
        call eirene_indmpi(pfluxm,dummy,ndx,ndy,nmol,ndxa,ndya,nmoli,
     .                                    ncutb,ncutl,npoint,npplg,1,1)
        call eirene_indmpi(pefluxm,dummy,ndx,ndy,nmol,ndxa,ndya,nmoli,
     .                                    ncutb,ncutl,npoint,npplg,1,1)
        call eirene_indmpi(tmb2,dummy,ndx,ndy,nmol,ndxa,ndya,nmoli,
     .                                    ncutb,ncutl,npoint,npplg,1,1)
        call eirene_indmpi(dib2,dummy,ndx,ndy,nion,ndxa,ndya,nioni,
     .                                    ncutb,ncutl,npoint,npplg,1,1)
        call eirene_indmpi(tib2,dummy,ndx,ndy,nion,ndxa,ndya,nioni,
     .                                    ncutb,ncutl,npoint,npplg,1,1)
        call eirene_indmpi(emiss,dummy,ndx,ndy,1,ndxa,ndya,1,
     .                                    ncutb,ncutl,npoint,npplg,1,1)
        call eirene_indmpi(emissmol,dummy,ndx,ndy,1,ndxa,ndya,1,
     .                                    ncutb,ncutl,npoint,npplg,1,1)
        call eirene_indmpi(srcml,dummy,ndx,ndy,nmol,ndxa,ndya,nmoli,
     .                                    ncutb,ncutl,npoint,npplg,1,1)
!csw 04mar2013
        do iistra = 0, nstrai
        call eirene_indmpi(eneutrad,
     .                     dummy,ndx,ndy,natm,ndxa,ndya,natmi,
     .                     ncutb,ncutl,npoint,npplg,nstra+1,iistra+1)
        call eirene_indmpi(emolrad,
     .                     dummy,ndx,ndy,nmol,ndxa,ndya,nmoli,
     .                     ncutb,ncutl,npoint,npplg,nstra+1,iistra+1)
        call eirene_indmpi(eionrad,
     .                     dummy,ndx,ndy,nion,ndxa,ndya,nioni,
     .                     ncutb,ncutl,npoint,npplg,nstra+1,iistra+1)
        call eirene_indmpi(edissml,dummy,ndx,ndy,nmol,ndxa,ndya,nmoli,
     .                     ncutb,ncutl,npoint,npplg,nstra+1,iistra+1)
        end do
      end if
      !c
      !c*** writing on file ft44
      !c
      NRED=(NPPLG-1)*(NCUTL-NCUTB)
      write (iunout,*) 'nred ',nred
      OPEN (UNIT=44,ACCESS='SEQUENTIAL',FORM='FORMATTED') ! added 19980603 dpc
      rewind (44)
      WRITE(44,'(i4,2x,i4,2x,i8)') ndxa-nred,ndya,jvft44
      write(44,'(i4,2x,i4,2x,i4)') natmi,nmoli,nioni
      !cank 960623
      nnatmi=natmi
      nnmoli=nmoli
      nnplsi=nplsi
      nns=nnplsi
      !cank
      do jatm=1,natmi
        write(44,*) texts(jatm+nsph)
      end do
      do jmol=1,nmoli
        write(44,*) texts(jmol+nspa)
      end do
      do jion=1,nioni
        write(44,*) texts(jion+nspam)
      end do
      call eirene_neutr(44,ndxa-nred,ndya,natmi,dab2,ndx,ndy,natm,1,1)
      call eirene_neutr(44,ndxa-nred,ndya,natmi,tab2,ndx,ndy,natm,1,1)
      if (nmoli.gt.0) then
        call eirene_neutr(44,ndxa-nred,ndya,nmoli,dmb2,ndx,ndy,nmol,1,1)
        call eirene_neutr(44,ndxa-nred,ndya,nmoli,tmb2,ndx,ndy,nmol,1,1)
      end if
      if (nioni.gt.0) then
        call eirene_neutr(44,ndxa-nred,ndya,nioni,dib2,ndx,ndy,nion,1,1)
        call eirene_neutr(44,ndxa-nred,ndya,nioni,tib2,ndx,ndy,nion,1,1)
      end if
      call eirene_neutr(44,ndxa-nred,ndya,natmi,rfluxa,ndx,ndy,natm,1,1)
      if (nmoli.gt.0) then
        call eirene_neutr(44,ndxa-nred,ndya,nmoli,
     .                       rfluxm,ndx,ndy,nmol,1,1)
      end if
      call eirene_neutr(44,ndxa-nred,ndya,natmi,pfluxa,ndx,ndy,natm,1,1)
      if (nmoli.gt.0) then
        call eirene_neutr(44,ndxa-nred,ndya,nmoli,
     .                       pfluxm,ndx,ndy,nmol,1,1)
      end if
      call eirene_neutr(44,ndxa-nred,ndya,natmi,
     .                     refluxa,ndx,ndy,natm,1,1)
      if (nmoli.gt.0) then
        call eirene_neutr(44,ndxa-nred,ndya,nmoli,
     .                       refluxm,ndx,ndy,nmol,1,1)
      end if
      call eirene_neutr(44,ndxa-nred,ndya,natmi,
     .                     pefluxa,ndx,ndy,natm,1,1)
      if (nmoli.gt.0) then
        call eirene_neutr(44,ndxa-nred,ndya,nmoli,
     .                       pefluxm,ndx,ndy,nmol,1,1)
      end if  
      call eirene_neutr(44,ndxa-nred,ndya,1,emiss,ndx,ndy,1,1,1)
      call eirene_neutr(44,ndxa-nred,ndya,1,emissmol,ndx,ndy,1,1,1)
      !cank 960511
      !c*** save the molecule-related sources...
      if (nmoli.gt.0) then
        call eirene_neutr(44,ndxa-nred,ndya,nmoli,
     .                       srcml,ndx,ndy,nmol,1,1)
        call eirene_neutr(44,ndxa-nred,ndya,nmoli,
     .                       edissml,ndx,ndy,nmol,1,1)
      end if  
      !c*** and the data on wall loading
      !cank 960513
      nnlimi=nlimi
      nnstsi=nstsi
      write(44,'(3i6)') nlimi, nstsi, nstrai
      call neutrs(44,wldnek,1)
      call neutrs(44,wldnep,1)
      call neutrs(44,wldna,natmi)
      call neutrs(44,ewlda,natmi)
      call neutrs(44,wldnm,nmoli)
      call neutrs(44,ewldm,nmoli)
      !c*** write down the wall geometry (only valid for polygons!)
      write (44,'(8f10.4)') (0.01*p1(1,ix), 0.01*p1(2,ix),
     .               0.01*p2(1,ix), 0.01*p2(2,ix), ix=1,nlimi)
      !c*** from 960623 on:
      call neutrs(44,wldra,natmi)
      call neutrs(44,wldrm,nmoli)
      !c*** from 960727 on:
      if(nstrai.gt.1) then
        do iistra=1,nstrai
          call neutrs(44,wldnek(1,iistra),1)
          call neutrs(44,wldnep(1,iistra),1)
          call neutrs(44,wldna(1,1,iistra),natmi)
          call neutrs(44,ewlda(1,1,iistra),natmi)
          if (nmoli.gt.0) call neutrs(44,wldnm(1,1,iistra),nmoli)
          if (nmoli.gt.0) call neutrs(44,ewldm(1,1,iistra),nmoli)
          call neutrs(44,wldra(1,1,iistra),natmi)
          if (nmoli.gt.0) call neutrs(44,wldrm(1,1,iistra),nmoli)
        end do
      end if
      !c*** from 961228 on:
      call neutrs(44,wldpp,nfla)
      call neutrs(44,wldpa,natmi)
      call neutrs(44,wldpm,nmoli)
      call neutrs(44,wldpeb,1)
      call neutrs(44,wldspt,1)
      if(nstrai.gt.1) then
        do iistra=1,nstrai
          call neutrs(44,wldpp(1,1,iistra),nfla)
          call neutrs(44,wldpa(1,1,iistra),natmi)
          if (nmoli.gt.0) call neutrs(44,wldpm(1,1,iistra),nmoli)
          call neutrs(44,wldpeb(1,iistra),1)
          call neutrs(44,wldspt(1,iistra),1)
        end do
      end if
      !c*** from 20000727 on:
      write(44,'(18i4)') (isrftype(i),i=1,nnlimi)
      write(44,'(18i4)') (isrftype(nlim+i),i=1,nnstsi)
      !cank
      rewind (44)
!cc<<<
!c      write (iunout,'()')
!c      write (iunout,*) '%%% Eirene data in Eirene %%%',nlmpgs
!c      write (iunout,*) 'nnatmi,nnmoli,nnstsi,nnlimi = ',
!c     ,        nnatmi,nnmoli,nnstsi,nnlimi
!c      do k=0,nstrai !{
!c       write (iunout,'()')
!c       if (nmoli.gt.0) then
!c         write (iunout,*) 'wldna, wldnm, wldra, wldrm:',k
!c         write (iunout,'(1x,20(6x,a2,i3.2))')
!c     ,         ('na',i,i=1,natmi),('nm',i,i=1,nmoli),
!c     ,         ('ra',i,i=1,natmi),('rm',i,i=1,nmoli)
!c       else
!c         write (iunout,*) 'wldna, wldra:',k
!c         write (iunout,'(1x,20(6x,a2,i3.2))')
!c     ,         ('na',i,i=1,natmi),('ra',i,i=1,natmi)
!c       end if
!c       do j=1,nlim+nstsi !{
!c         if(j.le.nlimi .or. j.gt.nlim) then !{
!c           hlp_pr=.false.
!c           do i=1,natmi !{
!c             hlp_pr= hlp_pr .or. wldna(j,i,k).ne.0.
!c     .                                         .or. wldra(j,i,k).ne.0.
!c           end do !}
!c           do i=1,nmoli !{
!c             hlp_pr= hlp_pr .or. wldnm(j,i,k).ne.0.
!c     .                                         .or. wldrm(j,i,k).ne.0.
!c           end do !}
!c           if(hlp_pr) then !{
!c             write (iunout,'(1p,i4,20e11.4)') j,
!c     ,              (wldna(j,i,k),i=1,natmi),(wldnm(j,i,k),i=1,nmoli),
!c     ,              (wldra(j,i,k),i=1,natmi),(wldrm(j,i,k),i=1,nmoli)
!c           end if !}
!c         end if !}
!c       end do !}
!c
!c       write (iunout,'()')
!c       write (iunout,*) 'wldpp:',k
!c       do j=nlim+1,nlim+nstsi !{
!c         write (iunout,'(1p,i4,20e9.2)') j,(wldpp(j,i,k),i=1,nfla)
!c       end do !}
!c       write (iunout,'()')
!c       if (nmoli.gt.0) then
!c         write (iunout,*) 'wldpeb, wldpa, wldpm:',k
!c       else
!c         write (iunout,*) 'wldpeb, wldpa:',k
!c       end if
!c       do j=nlim+1,nlim+nstsi !{
!c         write (iunout,'(1p,i4,20e9.2)') j,wldpeb(j,k),
!c     ,               (wldpa(j,i,k),i=1,natmi),(wldpm(j,i,k),i=1,nmoli)
!c       end do !}
!c      end do !}
!cc>>>
      return
      !c======================================================================

      contains
        subroutine neutrs(kard,dummy,ldmf)
        use eirmod_parmmod
        use eirmod_cadgeo
        use eirmod_clgin
        implicit none
        integer :: kard,ldmf,is,iif
        real(DP) :: dummy(nlimps,*)
        do iif=1,ldmf
          if(nlimi.gt.0) write(kard,'(5(e16.8))')
     .                   (dummy(is,iif),is=1,nlimi)
          if(nstsi.gt.0) write(kard,'(5(e16.8))')
     .                   (dummy(is+nlim,iif),is=1,nstsi)
        enddo
        end subroutine neutrs
      end subroutine eirene_extrab25_wneutrals

      !
      ! DETERMINING IF CELL (IX,IY) BELONGS TO B2 GRID
      !
      FUNCTION B2_CELL(IX,IY)

      IMPLICIT NONE

      LOGICAL :: B2_CELL
      INTEGER,INTENT(IN) :: IX,IY

      IF(IX.GE.0.AND.IX.LE.NDXP.AND.IY.GE.0.AND.IY.LE.NDYP) THEN
        B2_CELL=.TRUE.
      ELSE
        B2_CELL=.FALSE.
      END IF

      END FUNCTION B2_CELL

      subroutine eirene_extrab25_wneuclean
      implicit none
      if(allocated(dab2)) then
        dab2=0._DP
        dmb2=0._DP
        dib2=0._DP
        tab2=0._DP
        tmb2=0._DP
        tib2=0._DP
        rfluxa=0._DP
        rfluxm=0._DP
        refluxa=0._DP
        refluxm=0._DP
        pfluxa=0._DP
        pfluxm=0._DP
        pefluxa=0._DP
        pefluxm=0._DP
        emiss=0._DP
        emissmol=0._DP
        srcml=0._DP
        edissml=0._DP

        wldnek=0._DP
        wldnep=0._DP
        wldna=0._DP
        ewlda=0._DP
        wldnm=0._DP
        ewldm=0._DP
        wldra=0._DP
        wldrm=0._DP
        wldpp=0._DP
        wldpa=0._DP
        wldpm=0._DP
        wldpeb=0._DP
        wldspt=0._DP

        eneutrad=0._DP
        emolrad=0._DP
        eionrad=0._DP
      endif
      return
      end subroutine eirene_extrab25_wneuclean


      subroutine eirene_extrab25_iniusr_init(n_spcsrf,l_spcsrf,
     .          i_spcsrf,
     .          j_spcsrf,sps_sgrp,sps_absr,sps_trno,sps_trni,
     .          sps_mtri,sps_tmpr,sps_spph,sps_spch,
     .          sps_mtrl,sps_id)
      implicit none
      integer, intent(in) :: n_spcsrf,l_spcsrf(*),
     .      i_spcsrf(*), j_spcsrf(*), sps_sgrp(*)
      real(dp), intent(in) :: sps_absr(*), sps_trno(*), sps_trni(*),
     .      sps_mtri(*), sps_tmpr(*), sps_spph(*), sps_spch(*)
      character*8, intent(in) :: sps_mtrl(*), sps_id(*)

      bn_spcsrf=n_spcsrf
      allocate(bl_spcsrf(nlimps))
      allocate(bi_spcsrf(n_spcsrf))
      allocate(bj_spcsrf(n_spcsrf))
      allocate(bsps_sgrp(n_spcsrf))
      allocate(bsps_absr(n_spcsrf))
      allocate(bsps_trno(n_spcsrf))
      allocate(bsps_trni(n_spcsrf))
      allocate(bsps_mtri(n_spcsrf))
      allocate(bsps_tmpr(n_spcsrf))
      allocate(bsps_spph(n_spcsrf))
      allocate(bsps_spch(n_spcsrf))
      allocate(bsps_mtrl(n_spcsrf))
      allocate(bsps_id(n_spcsrf))

      bl_spcsrf(1:nlimps) = l_spcsrf(1:nlimps)
      bi_spcsrf(1:n_spcsrf) = i_spcsrf(1:n_spcsrf)
      bj_spcsrf(1:n_spcsrf) = j_spcsrf(1:n_spcsrf)
      bsps_sgrp(1:n_spcsrf) = sps_sgrp(1:n_spcsrf)
      bsps_absr(1:n_spcsrf) = sps_absr(1:n_spcsrf)
      bsps_trno(1:n_spcsrf) = sps_trno(1:n_spcsrf)
      bsps_trni(1:n_spcsrf) = sps_trni(1:n_spcsrf)
      bsps_mtri(1:n_spcsrf) = sps_mtri(1:n_spcsrf)
      bsps_tmpr(1:n_spcsrf) = sps_tmpr(1:n_spcsrf)
      bsps_spph(1:n_spcsrf) = sps_spph(1:n_spcsrf)
      bsps_spch(1:n_spcsrf) = sps_spch(1:n_spcsrf)
      bsps_mtrl(1:n_spcsrf) = sps_mtrl(1:n_spcsrf)
      bsps_id(1:n_spcsrf) = sps_id(1:n_spcsrf)

!pb 27012016
! flag indicating if subroutine iniusr is called from B2.5
      ini_iniusr = 1
      end subroutine eirene_extrab25_iniusr_init

      subroutine eirene_extrab25_cleanup
      implicit none

      if(allocated(bl_spcsrf)) then
        deallocate(bl_spcsrf)
        deallocate(bi_spcsrf)
        deallocate(bj_spcsrf)
        deallocate(bsps_sgrp)
        deallocate(bsps_absr)
        deallocate(bsps_trno)
        deallocate(bsps_trni)
        deallocate(bsps_mtri)
        deallocate(bsps_tmpr)
        deallocate(bsps_spph)
        deallocate(bsps_spch)
        deallocate(bsps_mtrl)
        deallocate(bsps_id)
      endif

      if(allocated(dab2)) then
        deallocate(dab2)
        deallocate(dmb2)
        deallocate(dib2)
        deallocate(tab2)
        deallocate(tmb2)
        deallocate(tib2)
        deallocate(rfluxa)
        deallocate(rfluxm)
        deallocate(refluxa)
        deallocate(refluxm)
        deallocate(pfluxa)
        deallocate(pfluxm)
        deallocate(pefluxa)
        deallocate(pefluxm)
        deallocate(emiss)
        deallocate(emissmol)
        deallocate(srcml)
        deallocate(edissml)
        deallocate(wldnek)
        deallocate(wldnep)
        deallocate(wldna)
        deallocate(ewlda)
        deallocate(wldnm)
        deallocate(ewldm)
        deallocate(wldra)
        deallocate(wldrm)
        deallocate(wldpp)
        deallocate(wldpa)
        deallocate(wldpm)
        deallocate(wldpeb)
        deallocate(wldspt)
        deallocate(isrftype)
        deallocate(eneutrad)
        deallocate(emolrad)
        deallocate(eionrad)
      endif
      ifirst_wneutral=0

      if(allocated(flux_save)) then
        deallocate(flux_save)
      endif

      if(allocated(plnxtri)) then
        deallocate(plnxtri, plnytri, pplnxtri, pplnytri)
      endif

      end subroutine eirene_extrab25_cleanup


      subroutine eirene_extrab25_emissivity
      implicit none
      integer :: istr, i, j, iadv, icell, ncelc, ix, iy

      istr = 0

      IF (IESTR.EQ.ISTR) THEN
C  NOTHING TO BE DONE
      ELSEIF (NFILEN.EQ.1.OR.NFILEN.EQ.2) THEN
        IESTR=ISTR
        CALL EIRENE_RSTRT(ISTR,NSTRAI,NESTM1,NESTM2,NADSPC,
     .             ESTIMV,ESTIMS,ESTIML,
     .             NSDVI1,SDVI1,NSDVI2,SDVI2,
     .             NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .             NSIGI_SPC,TRCFLE)
      ELSEIF ((NFILEN.EQ.6.OR.NFILEN.EQ.7).AND.ISTR.EQ.0) THEN
        IESTR=ISTR
        CALL EIRENE_RSTRT(ISTR,NSTRAI,NESTM1,NESTM2,NADSPC,
     .             ESTIMV,ESTIMS,ESTIML,
     .             NSDVI1,SDVI1,NSDVI2,SDVI2,
     .             NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .             NSIGI_SPC,TRCFLE)
      ELSE
        WRITE (IUNOUT,*) 'ERROR IN EXTRAB25_EMISSIVITY: ' //
     .                   'DATA FOR STRATUM ISTRA= ', ISTR
        WRITE (IUNOUT,*) 
     .    'ARE NOT AVAILABLE. EXTRAB25_EMISSIVITY ABANDONED'
        RETURN
      ENDIF

      emiss(:,:,1,1) = 0._dp
      emissmol(:,:,1,1) = 0._dp
      do i = 1, num_lines
        write (iunout,*) 'in EXTRAB25_EMISSIVITY, line no. = ',i
        if (mod_addv == 0) then
c storage saving mode: 
c ADDV is overwritten when a new line comes, within a run.
c Thus re-calculate the new emissivity profile on ADDV now
          call eirene_emissivity(istr, i, i, 0)
c       else
c Sufficiently large storage on ADDV additional tally array,
c for all lines and components. No need to reset ADDV tallies.
        end if

cdr run over components
        do j = 1, emis_lines(i)%num_compo
cdr  iadv: tally number on ADDV
          iadv = emis_lines(i)%compo(j)%iadv
          write (iunout,*) 'IADV, MODADDV = ',iadv, mod_addv

!pb  This is dangerous! It is implicitely assumed that the default 
!pb  emissivity model is used.
!pb  In case of a user specific model specified in the Eirene input 
!pb  this might produce rubbish.
          if (i.eq.1.and.istr.eq.0) then ! Ba-alpha emissivity for fort.44
            do icell = 1, nsbox
              ncelc=ncltal(icell)
              if (ncelc.gt.ntrii.or.ncelc.eq.0) cycle
              ix=ixtri(ncelc)
              iy=iytri(ncelc)
              if(b2_cell(ix,iy)) then
                if (j.ge.1 .and. j.le.2) then ! Atomic components
                  emiss(ix,iy,1,1)=emiss(ix,iy,1,1)+
     .                             addv(iadv,ncelc)*1.0d6
                else if (j.ge.3 .and. j.le.6) then ! Molecular components
                  emissmol(ix,iy,1,1)=emissmol(ix,iy,1,1)+
     ,                             addv(iadv,ncelc)*1.0d6
                end if
              end if
            end do
          end if
        end do
      end do

      end subroutine eirene_extrab25_emissivity

      end module eirmod_extrab25
