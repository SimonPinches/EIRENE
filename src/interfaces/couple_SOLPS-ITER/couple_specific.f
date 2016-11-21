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

      public :: eirene_extrab25_cleanup,eirene_extrab25_wneutrals
      public :: eirene_extrab25_wneuinit, eirene_extrab25_wneufill
      public :: eirene_extrab25_wneusave, eirene_extraB25_wneuclean
      public :: eirene_extrab25_alloc_mods
      public :: eirene_extrab25_iniusr_init

      ! eirdiag.h/eirdiag.f
      !c*** Volume data:
      !c***    srcml   :   power loss due to molecules, including
      !c***    edissml :   power loss due to molecule dissociation
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
      real*8, save, allocatable, dimension(:,:,:,:),public :: 
     .  dab2,dmb2,dib2,tab2,tmb2,tib2,rfluxa,rfluxm,refluxa,refluxm,
     .  pfluxa,pfluxm,pefluxa,pefluxm,emiss,emissmol,srcml,edissml
      real*8, save, allocatable, dimension(:,:),public :: 
     .  wldnek,wldnep
      real*8, save, allocatable, dimension (:,:,:), public :: 
     .  wldna,ewlda,wldnm,ewldm,wldra,wldrm,wldpp,wldpa,wldpm
      real*8, save, allocatable, dimension (:,:),public :: 
     .  wldpeb,wldspt
      real*8, save, allocatable, dimension (:,:,:,:), public :: 
     .  eneutrad
      integer, save, public :: nnlimi,nnstsi,nnatmi,nnmoli,nnioni
      integer, save, public :: nnplsi,nns, nnstrai
      integer, save, allocatable, public :: isrftype(:)
      logical, save, public :: lhalpha=.false.,lvib=.false.

      ! wneutral globals
      real*8, save :: DA31(0:8,0:8)
      real*8, save :: DP31(0:8,0:8)
      real*8, save :: DM31(0:8,0:8)
      real*8, save :: DI31(0:8,0:8)
      real*8, save :: DN31(0:8,0:8)
      real*8, save :: RHMH2(0:8),RH2PH2(0:8,0:8)
      CHARACTER, save :: FILNAM*8,H123*4,REAC*9,CRC*3
      logical, save :: hlp_pr
      integer, save :: ia1,ia2,ia3,iindex,ifirst_wneutral=0
      real*8, save :: hlp_cnv

      ! b2.5 neutrals parameters modicifcations
      integer, save, public :: bn_spcsrf
      integer, save, public, allocatable :: bl_spcsrf(:),bi_spcsrf(:)
      integer, save, public, allocatable :: bj_spcsrf(:),bsps_sgrp(:)
      real*8, save, public, allocatable :: bsps_absr(:),bsps_trno(:)
      real*8, save, public, allocatable :: bsps_mtri(:),bsps_tmpr(:)
      real*8, save, public, allocatable :: bsps_trni(:), bsps_spph(:)
      real*8, save, public, allocatable :: bsps_spch(:)
      character*8, save, public, allocatable :: bsps_mtrl(:), bsps_id(:)

      ! diag2 globals and parameters
      integer, parameter :: mgwtiesx=12 ! max. number of wall segments tied to a grid edge segment

      ! remaining bits and pieces from braeir common
      real*8, save, public :: chemical_sputter_yield,fchar_chemical
      integer, save, public :: igass_chemical,itsput_chemical,
     .                          issput_chemical

      !pb volumes of b2.5 cells
      real*8, save, public, allocatable :: volcel(:,:)

      ! normals of B2.5 cell edges per triangle
      real*8, save, public, allocatable :: plnxtri(:), plnytri(:)
      real*8, save, public, allocatable :: pplnxtri(:), pplnytri(:)

      !flux_save
      real*8, save, public, allocatable :: flux_save(:)

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
      call eirene_alloc_braeir(nnx,nny,nfl,ifoff)
      call eirene_alloc_eirbra(nnx,nny,nfl,nstra,ifoff)
      allocate (flux_save(nstra))
      flux_save=0.d0
      end subroutine

      subroutine eirene_extrab25_wneutrals
      ! old version : 27.07.2000 23:07
      ! new version : 25.01.2011 (s.wiesen@fz-juelich.de) - f90 module
      !c     this subroutine produces an output file ft44 for plotting in B2
      !c     with the neutral densities, temperatures, and fluxes
      !c     (summed up for all strata)
      implicit none
      real*8 :: dummy(0:ndxp,0:ndyp)
      !c*** label for fort.44 file
      integer, parameter  :: jvft44=20000727
      !c*** and dissociation energy of the hydrogen molecule
      real*8, parameter  :: diss_pot_H2=4.48
      !c*** radiative transition prob. level 3-->2 (1/sec) for H-alpha calc.
      real*8, parameter :: fac32=4.410e7
      !C*** ionization potentials
      integer,parameter :: npot=20
      real(dp),save :: pot_data(npot),pot
      data pot_data /13.598, !H
     .               24.587, !He 
     .                5.392, !Li
     .                9.322, !Be
     .                8.298, !B
     .               11.260, !C
     .               14.534, !N
     .               13.618, !O
     .               17.422, !F
     .               21.564, !Ne
     .                5.139, !Na 
     .                7.646, !Mg 
     .                5.986, !Al
     .                8.151, !Si
     .               10.486, !P                                         
     .               10.360, !S                                
     .               12.967, !Cl                        
     .               15.759, !Ar                
     .                4.341, !K        
     .                6.113/ !Ca


      integer :: ix,iy,ir,ierror,i,j,in,iistra
      integer :: l,k,nred
      real*8 :: de,te,hlp,sigadd1,sigadd2,sigadd3,sigadd4,sigadd5
      real*8 :: dej,tei,tef,def,powalf1,powalf2,powalf3,powalf4,powalf5
      real*8 :: powalf,datm3,dpls3,dmol3,dion3,dnml3,sigadd
      real*8 :: da,dp,dm,di,dn,ratio2,ratio7
      integer :: istra_in,istra_save
      real*8 :: rcmin,rcmax,fp(6),vl
      integer :: jfexmn,jfexmx
      real*8 :: value
      external eirene_indmpi,eirene_neutr,eirene_slreac

!     !c======================================================================
      !c---------------------------------------------------------------------<
      entry eirene_extraB25_wneuinit
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
!pb      allocate(eneutrad(0:ndxp,0:ndyp,natm,1))
      allocate(eneutrad(0:ndxp,0:ndyp,natm,0:nstra+1))
      allocate(isrftype(nlmpgs))

      allocate(volcel(0:ndxp,0:ndyp))

      call eirene_extraB25_wneuclean

      hlp_cnv=1./elcha
      ia1=natmi+nmoli
      ia2=2*ia1
      ia3=3*ia1

      !c
      !c*** Initialise the data for H-alpha radiation
      !c
      if(lhalpha) then
          write(*,*) 'Using new SIGHA (941017)'
          IERROR=0
          rcmin=-huge(1.d0)
          rcmax= huge(1.d0)
          jfexmn=0
          jfexmx=0
          fp=0.d0
          !C
          !C  READ REDUCED POPULATION COEFFICIENT FOR HYDR. ATOMS FROM FILE AMJUEL
          !C  AND PUT THEM FROM CREAC(..,..,IR) ONTO DA,DP,DM,DI, AND DN ARRAY
          !C
          IR=NREACI
          IF (IR+7.GT.NREAC) then
            WRITE (6,*) 'FROM SUBROUTINE HALFA: '
            CALL EIRENE_MASPRM('NREAC',5,NREAC,'IR',2,IR,IERROR)
            CALL EIRENE_EXIT_OWN(1)
          end if
          !C
          FILNAM='AMJUEL  '
          H123='H.12'
          CRC='OT '
          !C
          !C  H(n=3)/H(n=1)
          REAC='2.1.5a   '
          IR=IR+1
          CALL eirene_slreac(IR,FILNAM,H123,REAC,CRC,
     .                    rcmin,rcmax,fp,jfexmn,jfexmx,'  ',0)
          do J=1,9
            do I=1,9
              !DA31(J-1,I-1)=CREAC(J,I,NREACI+1)
              DA31(J-1,I-1)=REACDAT(IR)%OTH%POLY%DBLPOL(J,I)
            end do
          end do
          !C  H(n=3)/H+
          REAC='2.1.8a   '
          IR=IR+1
          CALL eirene_slreac(IR,FILNAM,H123,REAC,CRC,
     .                    rcmin,rcmax,fp,jfexmn,jfexmx,'  ',0)
          do J=1,9
            do I=1,9
              !DP31(J-1,I-1)=CREAC(J,I,NREACI+1)
              DP31(J-1,I-1)=REACDAT(IR)%OTH%POLY%DBLPOL(J,I)
            end do
          end do
          !C  H(n=3)/H2(g)
          REAC='2.2.5a   '
          IR=IR+1
          CALL eirene_slreac(IR,FILNAM,H123,REAC,CRC,
     .                    rcmin,rcmax,fp,jfexmn,jfexmx,'  ',0)
          do J=1,9
            do I=1,9
              !DM31(J-1,I-1)=CREAC(J,I,NREACI+1)
              DM31(J-1,I-1)=REACDAT(IR)%OTH%POLY%DBLPOL(J,I)
            end do
          end do
          !C  H(n=3)/H2+(g)
          REAC='2.2.14a  '
          IR=IR+1
          CALL eirene_slreac(IR,FILNAM,H123,REAC,CRC,
     .                        rcmin,rcmax,fp,jfexmn,jfexmx,'  ',0)
          do J=1,9
            do I=1,9
              !DI31(J-1,I-1)=CREAC(J,I,NREACI+1)
              DI31(J-1,I-1)=REACDAT(IR)%OTH%POLY%DBLPOL(J,I)
            end do
          end do
          !C  H(n=3)/H-
          REAC='7.2a     '
          IR=IR+1
          CALL eirene_slreac(IR,FILNAM,H123,REAC,CRC,
     .                        rcmin,rcmax,fp,jfexmn,jfexmx,'  ',0)
          do J=1,9
            do I=1,9
              !DN31(J-1,I-1)=CREAC(J,I,NREACI+1)
              DN31(J-1,I-1)=REACDAT(IR)%OTH%POLY%DBLPOL(J,I)
            end do
          end do
          !C
          !C  NOW READ RATIO OF DENSITIES:
          !C
          !C  FIRST: H-/H2
          FILNAM='AMJUEL  '
          H123='H.11'
          !csw 28jan2011 changed from 7.0 to 7.0b, CHECK
          REAC='7.0b    '
          CRC='OT '
          IR=IR+1
          CALL eirene_slreac(IR,FILNAM,H123,REAC,CRC,
     .                        rcmin,rcmax,fp,jfexmn,jfexmx,'  ',0)
          do I=1,9
            !RHMH2(I-1)=CREAC(I,1,NREACI+1)
            RHMH2(I-1)=REACDAT(IR)%OTH%POLY%DBLPOL(I,1)
          end do

          !C  NEXT : H2+/H2
          FILNAM='AMJUEL  '
          H123='H.12'
          REAC='2.0c    '
          CRC='OT '
          IR=IR+1
          !C  2.0C INCLUDES ION CONVERION (CX) ON H2(V)
          !C  OLD VERSION (WITHOUT THIS CX) SHOULD BE RECOVERED BY
          !C  READING 2.0B INSTEAD, AND OMITTING THE H- CHANNEL 5.
          CALL eirene_slreac(IR,FILNAM,H123,REAC,CRC,
     .                        rcmin,rcmax,fp,jfexmn,jfexmx,'  ',0)
           do I=1,9
            do J=1,9
              !RH2PH2(I-1,J-1)=CREAC(I,J,NREACI+1)
              RH2PH2(I-1,J-1)=REACDAT(IR)%OTH%POLY%DBLPOL(I,J)
            end do
          end do
          !C
          write(*,*) 'NREAC,NREACI,IR     ',NREAC,NREACI,IR
          !write(*,*) 'NRCX,IRCX        ',NRCX,IRCX
          !write(*,*) 'NREL,IREL        ',NREL,IREL
          !write(*,*) 'NRII,IRII        ',NRII,IRII
          !write(*,*) 'NELI,NAELI       ',NELI,naeli
          !write(*,*) 'NREI,NMEII,NIEII ',NREI,nmeii,nieii
          !write(*,*) 'NREC,NIRCI,NPRCI ',NREC,nirci,nprci
      end if
      iindex=0
      !c======================================================================
      !c*** fill arrays
      !c
      write(*,*) 'NSTRAI,NESTIM,NSDVI',nstrai,nestim,nsdvi
      nnstrai=nstrai
      !c

      !c---------------------------------------------------------------------<
      return
      !c
      !c======================================================================
      entry eirene_extraB25_wneufill(istra_in)
      !c
      !c      print *,'%%% wneufill: istra,istra_in = ',istra,istra_in
      istra_save=istra
      istra=istra_in
      !c--------------------------------------------------------------------->

      !csw
      !csw 21feb2012 corrected radiation from neutrals (atoms only), taken from SOLPS4.3 (V.Kotov)
      !csw 04mar2013 shifted from wneusave to here (wneufill)
      !csw
      value=0.0
      eneutrad(:,:,1,istra) = 0.d0
      do ncell=1,ntrii
        ix=ixtri(ncell)
        iy=iytri(ncell)
        if(ix.gt.0) then
           eneutrad(ix,iy,1,istra)=eneutrad(ix,iy,1,istra) 
     .                      +eael(ncell)*vol(ncell)
          value=value+eael(ncell)*vol(ncell)
          do iatm=1,natmi
            if(nchara(iatm).le.npot) then
              pot=pot_data(nchara(iatm))
            else
              pot=0.
            endif
            eneutrad(ix,iy,1,istra)=eneutrad(ix,iy,1,istra) 
     .                        -paat(iatm,ncell)*pot*vol(ncell)
            value=value-paat(iatm,ncell)*pot*vol(ncell)
          enddo
        endif
      enddo
      !open(555,file='eneutrad.dat',form='formatted')
      !do ix=1,76
      !  do iy=1,28
      !    write(555,'(i6,i6,1x,e13.6)') ix,iy,eneutrad(ix,iy,1,1)
      !  enddo
      !enddo
      !close(555)
cdr    write(*,'(a,i6,1x,e13.6)') 'DBG: ISTRA, ENEUTRAD',istra,value
      !csw

      !c
      !C map 1d-EIRENE neutral densities and temperatures on 2d-arrays
      !c and change the units to SI for plotting in B2
      !c
cdr   write(6,*) 'istra ',istra
cdr   write(6,*) 'natmi, nmoli, nioni ',natmi,nmoli,nioni
      !crfs     IF (WTOTP(0,ISTRA).EQ.0.) GOTO 60
      do iatm=1,natmi
        volcel = 0.d0
        do in=1,ntrii
          ix=ixtri(in)
          iy=iytri(in)
          if(ix.gt.0) then
            volcel(ix,iy) = volcel(ix,iy) + vol(in)
          end if
        end do
        do in=1,ntrii
          ix=ixtri(in)
          iy=iytri(in)
          if(ix.gt.0) then
            vl = vol(in)/volcel(ix,iy)
            dab2(ix,iy,iatm,1)=dab2(ix,iy,iatm,1)+ 
     .                          pdena(iatm,in)*1e6*vl
            if(iindex+iatm.le.nadv) 
     .         rfluxa(ix,iy,iatm,1)=rfluxa(ix,iy,iatm,1)+ 
     .                                   addv(iindex+iatm,in)*1.0e4
            if(iindex+ia2+iatm.le.nadv) 
     .         pfluxa(ix,iy,iatm,1)=pfluxa(ix,iy,iatm,1)+ 
     .                                   addv(iindex+ia2+iatm,in)*1.0e4
            if(iindex+ia1+iatm.le.nadv) 
     .         refluxa(ix,iy,iatm,1)=refluxa(ix,iy,iatm,1)+ 
     .                             addv(iindex+ia1+iatm,in)*1.0e4*elcha
            if(iindex+ia3+iatm.le.nadv) 
     .         pefluxa(ix,iy,iatm,1)=pefluxa(ix,iy,iatm,1)+ 
     .                             addv(iindex+ia3+iatm,in)*1.0e4*elcha
            tab2(ix,iy,iatm,1)=tab2(ix,iy,iatm,1)+edena(iatm,in)*vl
          endif
        end do
      end do
      do imol=1,nmoli
        edissml(:,:,imol,istra) = 0.d0
        do in=1,ntrii
          ix=ixtri(in)
          iy=iytri(in)
          if(ix.gt.0) then
            vl = vol(in)/volcel(ix,iy)
            dmb2(ix,iy,imol,1)=dmb2(ix,iy,imol,1)+ 
     .                          pdenm(imol,in)*1e6*vl
            if(iindex+natmi+imol.le.nadv) 
     .         rfluxm(ix,iy,imol,1)=rfluxm(ix,iy,imol,1)+ 
     .                                 addv(iindex+natmi+imol,in)*1.0e4
            if(iindex+ia2+natmi+imol.le.nadv) 
     .         pfluxm(ix,iy,imol,1)=pfluxm(ix,iy,imol,1)+ 
     .                             addv(iindex+ia2+natmi+imol,in)*1.0e4
            if(iindex+ia1+natmi+imol.le.nadv) 
     .         refluxm(ix,iy,imol,1)=refluxm(ix,iy,imol,1)+ 
     .                       addv(iindex+ia1+natmi+imol,in)*1.0e4*elcha
            if(iindex+ia3+natmi+imol.le.nadv) 
     .         pefluxm(ix,iy,imol,1)=pefluxm(ix,iy,imol,1)+ 
     .                       addv(iindex+ia3+natmi+imol,in)*1.0e4*elcha
            tmb2(ix,iy,imol,1)=tmb2(ix,iy,imol,1)+edenm(imol,in)*vl
            srcml(ix,iy,imol,1)=srcml(ix,iy,imol,1)+ 
     .                                           pmml(imol,in)*vol(in)
            !c*** Potential energy source related to hydrogen molecules
            !c*** for determining the radiation (W per cell)
            if (ncharm(imol).eq.2) 
     .         edissml(ix,iy,imol,istra)=edissml(ix,iy,imol,istra)+ 
     .                               pmml(imol,in)*vol(in)*diss_pot_H2
          endif
        end do
      end do
      do iion=1,nioni
        do in=1,ntrii
          ix=ixtri(in)
          iy=iytri(in)
          if(ix.gt.0) then
            vl = vol(in)/volcel(ix,iy)
            dib2(ix,iy,iion,1)=dib2(ix,iy,iion,1)+ 
     .                          pdeni(iion,in)*1e6*vl
            tib2(ix,iy,iion,1)=tib2(ix,iy,iion,1)+edeni(iion,in)*vl
          endif
        end do
      end do
      !c
      !c*** Re-scale the surface data from A to 1/sec and average the energy
      !c
      do i=1,nlimps
        wldnek(i,istra)=0.
        wldnep(i,istra)=0.
        wldpeb(i,istra)=0.
        !c*** hlp accumulates the power taken away with re-emitted particles
        hlp=0.
        do j=1,natmi
          wldnek(i,istra)=wldnek(i,istra)+eotat(j,i)
          wldpeb(i,istra)=wldpeb(i,istra)+erfpat(j,i)
          hlp=hlp+erfaat(j,i)+erfmat(j,i)+erfiat(j,i)
          if(potat(j,i).gt.0.) then
            ewlda(i,j,istra)=eotat(j,i)/potat(j,i)
          else
            ewlda(i,j,istra)=0.
          end if
          wldna(i,j,istra)=hlp_cnv*potat(j,i)
          wldra(i,j,istra)=hlp_cnv*(prfaat(j,i)+prfmat(j,i)+ 
     .                                                    prfiat(j,i))
          wldpa(i,j,istra)=hlp_cnv*prfpat(j,i)
        end do
        do j=1,nmoli
          wldnek(i,istra)=wldnek(i,istra)+eotml(j,i)
          wldpeb(i,istra)=wldpeb(i,istra)+erfpml(j,i)
          hlp=hlp+erfaml(j,i)+erfmml(j,i)+erfiml(j,i)
          if(potml(j,i).gt.0.) then
            ewldm(i,j,istra)=eotml(j,i)/potml(j,i)
          else
            ewldm(i,j,istra)=0.
          end if
          wldnm(i,j,istra)=hlp_cnv*potml(j,i)
          wldrm(i,j,istra)=hlp_cnv*(prfaml(j,i)+prfmml(j,i)+ 
     .                                                    prfiml(j,i))
          wldpm(i,j,istra)=hlp_cnv*prfpml(j,i)
        end do
        do j=1,nioni
          wldnek(i,istra)=wldnek(i,istra)+eotio(j,i)
          wldpeb(i,istra)=wldpeb(i,istra)+erfpio(j,i)
          hlp=hlp+erfaio(j,i)+erfmio(j,i)+erfiio(j,i)
        end do
        do j=1,nplsi
          wldpp(i,j,istra)=hlp_cnv*potpl(j,i)
        end do
        !c*** recombination energy of hydrogen molecules
        do j=1,nmoli
          if(ncharm(j).eq.2) wldnep(i,istra)=wldnep(i,istra)+ 
     .                                         diss_pot_H2*prfaml(j,i)
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
      !c      print *,'%%% wneusave: istra = ',istra
      !c*** Calculate the totals (stratum 0)
      !c
      do i=1,nlimps
        wldnek(i,0)=0.
        wldnep(i,0)=0.
        wldpeb(i,0)=0.
        wldspt(i,0)=0.
        do j=1,natmi
          ewlda(i,j,0)=0.
          wldna(i,j,0)=0.
          wldra(i,j,0)=0.
          wldpa(i,j,0)=0.
        end do
        do j=1,nmoli
          ewldm(i,j,0)=0.
          wldnm(i,j,0)=0.
          wldrm(i,j,0)=0.
          wldpm(i,j,0)=0.
        end do
        do j=1,nplsi
          wldpp(i,j,0)=0.
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
          do j=1,nplsi
            wldpp(i,j,0)=wldpp(i,j,0)+wldpp(i,j,k)
          end do
        end do
        do j=1,natmi
          if(wldna(i,j,0).gt.0.) then
            ewlda(i,j,0)=ewlda(i,j,0)/wldna(i,j,0)
          else
            ewlda(i,j,0)=0.
          end if
        end do
        do j=1,nmoli
          if(wldnm(i,j,0).gt.0.) then
            ewldm(i,j,0)=ewldm(i,j,0)/wldnm(i,j,0)
          else
            ewldm(i,j,0)=0.
          end if
        end do
      end do
!pb
      eneutrad(:,:,:,0) = 0.
      edissml(:,:,:,0) = 0.
      do ix = 1, ndxa
        do iy = 1, ndya
          eneutrad(ix,iy,1,0) = sum(eneutrad(ix,iy,1,1:nstrai))
          do imol=1,nmoli
            edissml(ix,iy,imol,0) = sum(edissml(ix,iy,imol,1:nstrai)) 
          end do
        end do
      end do 
      !c--------------------------------------------------------------------->


      !c*** Surface type
      do i=1,nlimps
        isrftype(i)=iliin(i)
      end do

      !c
      !c*** Calculate the temperatures
      !c
      do ix=1,ndxa
        do iy=1,ndya
          do iatm=1,natmi
            if (dab2(ix,iy,iatm,1).gt.0.) then
              tab2(ix,iy,iatm,1)=tab2(ix,iy,iatm,1)/dab2(ix,iy,iatm,1)* 
     .                                                  elcha*2./3.*1.e6
            else
              tab2(ix,iy,iatm,1)=1.e-6*elcha
            end if
          end do
          do imol=1,nmoli
            if (dmb2(ix,iy,imol,1).gt.0.) then
              tmb2(ix,iy,imol,1)=tmb2(ix,iy,imol,1)/dmb2(ix,iy,imol,1)* 
     .                                                  elcha*2./3.*1.e6
            else
              tmb2(ix,iy,imol,1)=1.e-6*elcha
            end if
          end do
          do iion=1,nioni
            if (dib2(ix,iy,iion,1).gt.0.) then
              tib2(ix,iy,iion,1)=tib2(ix,iy,iion,1)/dib2(ix,iy,iion,1)* 
     .                                                  elcha*2./3.*1.e6
            else
              tib2(ix,iy,iion,1)=1.e-6*elcha
            end if
          end do
        end do
      end do
      !c
      !c*** Calculate H-alpha emissivity
      !c
      !C SEPT. 96: REVISE  CH. NO 4 (COUPLING TO H2+): INCLUDE ION CONVERSION ON
      !C SEPT. 96: INCLUDE CH. NO 5 (COUPLING TO H-)
      !
      !CSW 28jan2011: we are not using the official EIRENE routine ba_halpha CHECK !
      !
      if(lhalpha) then
        POWALF=0.
        POWALF1=0.
        POWALF2=0.
        POWALF3=0.
        POWALF4=0.
        POWALF5=0.
        do ncell=1,ntrii
          ix=ixtri(ncell)
          iy=iytri(ncell)
          if(ix.gt.0) then
            !C
            !C  LOCAL PLASMA DATA
            !C
            TE=TEIN(NCELL)
            DE=DEIN(NCELL)
            SIGADD1=0.
            SIGADD2=0.
            SIGADD3=0.
            SIGADD4=0.
            SIGADD5=0.
            IF (.not. LGVAC(NCELL,0)) then
              DEF=LOG(DE*1.D-8)
              TEF=LOG(TE)
              DATM3=0.
              DPLS3=0.
              DMOL3=0.
              DION3=0.
              DNML3=0.
              do J=0,8
                DEJ=DEF**J
                do I=0,8
                  TEI=TEF**I
                  DATM3=DATM3+DA31(I,J)*TEI*DEJ
                  DPLS3=DPLS3+DP31(I,J)*TEI*DEJ
                  DMOL3=DMOL3+DM31(I,J)*TEI*DEJ
                  DION3=DION3+DI31(I,J)*TEI*DEJ
                  DNML3=DNML3+DN31(I,J)*TEI*DEJ
                end do
              end do
              if(datm3.lt.500) then
                DATM3=EXP(DATM3)
              else
                write(*,*) '[DPC] Problem in wneusave: ln(datm3) = ', 
     .            datm3, ' --- exponential will overflow'
                write(*,*) '[DPC] TE, DE = ', TE, DE
                datm3=1d30
              endif
              if(dpls3.lt.500) then
                DPLS3=EXP(DPLS3)
              else
                write(*,*) '[DPC] Problem in wneusave: ln(dpls3) = ', 
     .            dpls3, ' --- exponential will overflow'
                write(*,*) '[DPC] TE, DE = ', TE, DE
                dpls3=1d30
              endif
              if(dmol3.lt.500) then
                DMOL3=EXP(DMOL3)
              else
                write(*,*) '[DPC] Problem in wneusave: ln(dmol3) = ', 
     .            dmol3, ' --- exponential will overflow'
                write(*,*) '[DPC] TE, DE = ', TE, DE
                dmol3=1d30
              endif
              if(dion3.lt.500) then
                DION3=EXP(DION3)
              else
                write(*,*) '[DPC] Problem in wneusave: ln(dion3) = ', 
     .            dion3, ' --- exponential will overflow'
                write(*,*) '[DPC] TE, DE = ', TE, DE
                dion3=1d30
              endif
              if(dnml3.lt.500) then
                DNML3=EXP(DNML3)
              else
                write(*,*) '[DPC] Problem in wneusave: ln(dnml3) = ', 
     .            dnml3, ' --- exponential will overflow'
                write(*,*) '[DPC] TE, DE = ', TE, DE
                dnml3=1d30
              endif

              !C  RATIO OF DENSITIES: H- TO H2, COLL. EQUIL. IN VIBRATION
              !C  (ONLY TE-DEPENDENT)

              RATIO7=0
              do I=0,8
                TEI=TEF**I
                RATIO7=RATIO7+RHMH2(I)*TEI
              end do
              RATIO7=EXP(RATIO7)

              !C  RATIO OF DENSITIES: H2+ TO H2, INCL. ION CONVERSION

              RATIO2=0
              do J=0,8
                DEJ=DEF**J
                do I=0,8
                  TEI=TEF**I
                  RATIO2=RATIO2+RH2PH2(I,J)*TEI*DEJ
                end do
              end do
              if(RATIO2.lt.500) then
                RATIO2=EXP(RATIO2)
              else
                write(*,*) '[DPC] Problem in wneusave: ln(ratio2) = ', 
     .            RATIO2,' --- exponential will overflow'
                RATIO2=1d30
              endif

              !C
              !C  CHANNEL 1
              !C  H ALPHA SOURCE RATE:  PHOTONS/SEC/M**3
              !c  linear in atomic density (ionisation)
              !C
              do  IATM=1,NATMI
                !C  HYDROGENIC SPECIES?
                IF (NCHARA(IATM).eq.1) then
                  !c DA=DATM3*PDENA(IATM,NCELL)
                  da=datm3*dab2(ix,iy,iatm,1)
                  SIGADD1=SIGADD1+DA*FAC32
                end if
              end do
              !C
              !C  CHANNEL 2
              !C  H ALPHA SOURCE RATE:  PHOTONS/SEC/M**3
              !c  linear in ion density (recombination)
              !C
              do  IPLS=1,NPLSI
                IF (NCHARP(IPLS).eq.1) then
                  DP=DPLS3*DIIN(IPLS,NCELL)
                  SIGADD2=SIGADD2+DP*FAC32
                end if
              end do
              sigadd2=1.e6*sigadd2
              !C
              !C  CHANNEL 3
              !C  H ALPHA SOURCE RATE:  PHOTONS/SEC/M**3
              !c  linear in molecular density (dissociation of H2)
              !C
              do  IMOL=1,NMOLI
                IF (NCHARM(IMOL).eq.2) then
                  !c DM=DMOL3*PDENM(IMOL,NCELL)
                  dm=dmol3*dmb2(ix,iy,imol,1)
                  SIGADD3=SIGADD3+DM*FAC32
                end if
              end do
              !C
              !C  CHANNEL 4
              !C  H ALPHA SOURCE RATE:  PHOTONS/SEC/M**3
              !C  LINEAR IN PDENI: (DISSOCIATION OF H2+)
              !C
              !C  REVISED: USE (PDENM * DENSITY RATIO H2+/H2) NOW, INSTEAD OF PDENI
              !C
              !C      DO 215 IION=1,NIONI
              !C        IF (NCHARI(IION).NE.2) GOTO 215
              !C        DI=DION3*PDENI(IION,NCELL)
              !C        SIGADD4=SIGADD4+DI*FAC32
              !C215    CONTINUE
              do IMOL=1,NMOLI
                IF (NCHARM(IMOL).eq.2) then
                  !c DI=DION3*PDENM(IMOL,NCELL)*RATIO2
                  di=dion3*dmb2(ix,iy,imol,1)*ratio2
                  SIGADD4=SIGADD4+DI*FAC32
                end if
              end do
              !C
              !C  CHANNEL 5
              !C  H ALPHA SOURCE RATE:  PHOTONS/SEC/M**3
              !C  LINEAR IN H- DENSITY (CHARGE-EXCHANGE RECOMBINATION)
              !C
              do IMOL=1,NMOLI
                IF (NCHARM(IMOL).eq.2) then
                  !c DN=DNML3*PDENM(IMOL,NCELL)*RATIO7
                  dn=dnml3*dmb2(ix,iy,imol,1)*ratio7
                  SIGADD5=SIGADD5+DN*FAC32
                end if
              end do
            end if
            !C
            EMISS(ix,iy,1,1)=EMISS(ix,iy,1,1)+(SIGADD1+SIGADD2)
            if(lvib) then
              EMISSMOL(ix,iy,1,1)=EMISSMOL(ix,iy,1,1)+ 
     .                                         (SIGADD3+SIGADD4+SIGADD5)
            else
              EMISSMOL(ix,iy,1,1)=EMISSMOL(ix,iy,1,1)+SIGADD3
            end if
            SIGADD=SIGADD1+SIGADD2+SIGADD3+SIGADD4+SIGADD5
            powalf=powalf+sigadd*3.028e-25*vol(ncell)
            powalf1=powalf1+sigadd1*3.028e-25*vol(ncell)
            powalf2=powalf2+sigadd2*3.028e-25*vol(ncell)
            powalf3=powalf3+sigadd3*3.028e-25*vol(ncell)
            powalf4=powalf4+sigadd4*3.028e-25*vol(ncell)
            powalf5=powalf5+sigadd5*3.028e-25*vol(ncell)
          endif
        end do
        WRITE (6,*) ' RADIATED POWER BY HALPHA:',POWALF
        WRITE (6,*) ' COUPL. TO GROUNDSTATE   :',POWALF1
        WRITE (6,*) ' COUPLING TO CONTINUUM   :',POWALF2
        WRITE (6,*) ' COUPLING TO MOLECULES   :',POWALF3
        WRITE (6,*) ' COUPLING TO MOL.IONS    :',POWALF4
        WRITE (6,*) ' COUPLING TO NEG.IONS    :',POWALF5
      end if

      !c*** print some neutral fluxes across the "non-default" surfaces

      write(6,'(/17x,2(2x,a8),100(3x,a6,i2))') '  area','  power', 
     .             ('atflx',i,i=1,natmi),('mlflx',i,i=1,nmoli),
     .             ('atflxr',i,i=1,natmi),('mlflxr',i,i=1,nmoli),
     .             ('atflxp',i,i=1,natmi),('mlflxp',i,i=1,nmoli),
     .             ('plsflx',i,i=1,nplsi)
      do i=1,nstsi
        j=nlim+i
        do k=1,nplsi
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

        if(wldnek(j,0).ne.0. .or. 
     .      wldna(j,1,0).ne.0 .or. wldnm(j,1,0).ne.0 .or.
     .      wldra(j,1,0).ne.0 .or. wldrm(j,1,0).ne.0 .or.
     .      wldpa(j,1,0).ne.0 .or. wldpm(j,1,0).ne.0 .or.
     .      wldpp(j,1,0).ne.0) 
     .         write(6,'(a,i4,1p,100e11.3)')
     .           'non-def-surf ',i, 1.e-4*sarea(j), 1.e-6*wldnek(j,0),
     .               (wldna(j,k,0),k=1,natmi),(wldnm(j,k,0),k=1,nmoli),
     .               (wldra(j,k,0),k=1,natmi),(wldrm(j,k,0),k=1,nmoli),
     .               (wldpa(j,k,0),k=1,natmi),(wldpm(j,k,0),k=1,nmoli),
     .               (wldpp(j,k,0),k=1,nplsi)
      end do
!cc%%%
!c      print *,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!c      print '(/6x,20(a8,2x))',
!c     ,         'wldnek','wldnep','wldna','ewlda','wldnm','ewldm',
!c     ,         'wldna He','ewlda He','wldna Ne','ewlda Ne','wldra H',
!c     ,         'wldra He','wldra Ne','wldrm','prfaat'
!c      do i=1,nlimi
!c      print '(1p,i6,20e10.2)',i,wldnek(i),wldnep(i),
!c     ,       wldna(i,1),ewlda(i,1),wldnm(i,1),ewldm(i,1),
!c     ,       wldna(i,2),ewlda(i,2),wldna(i,3),ewlda(i,3),
!c     ,       wldra(i,1),wldra(i,2),wldra(i,3),wldrm(i,1),
!c     ,       prfaat(1,i)
!c      end do
!c      print '(/6x,20(a8,2x))',
!c     ,         'wldnek','wldnep','wldna','ewlda','wldnm','ewldm',
!c     ,         'wldna He','ewlda He','wldna Ne','ewlda Ne','wldra H',
!c     ,         'wldra He','wldra Ne','wldrm','prfaat'
!c      do i=nlim+1,nlim+nstsi
!c      print '(1p,i6,20e10.2)',i-nlim,wldnek(i),wldnep(i),
!c     ,       wldna(i,1),ewlda(i,1),wldnm(i,1),ewldm(i,1),
!c     ,       wldna(i,2),ewlda(i,2),wldna(i,3),ewlda(i,3),
!c     ,       wldra(i,1),wldra(i,2),wldra(i,3),wldrm(i,1),
!c     ,       prfaat(1,i)
!c      end do
!c      print *,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!cc%%%
      write(6,*) 'ncutl,ncutb ',ncutl,ncutb
      write(6,*) 'ndx,ndy,natm,ndxa,ndya,nfla,n1st'
      write(6,*) ndx,ndy,natm,ndxa,ndya,nfla,n1st
      !c
      !c*** backmapping of 2d-arrays for b2
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
        do istra = 0, nstrai
        call eirene_indmpi(eneutrad,
     .                     dummy,ndx,ndy,natm,ndxa,ndya,natmi,
     .                     ncutb,ncutl,npoint,npplg,nstra+1,istra+1)
        call eirene_indmpi(edissml,dummy,ndx,ndy,nmol,ndxa,ndya,nmoli,
     .                     ncutb,ncutl,npoint,npplg,nstra+1,istra+1)
        end do
      end if
      !c
      !c*** writing on file ft44
      !c
      NRED=(NPPLG-1)*(NCUTL-NCUTB)
      write (6,*) 'nred ',nred
      OPEN (UNIT=44,ACCESS='SEQUENTIAL',FORM='FORMATTED')    ! added 19980603 dpc
      rewind (44)
      WRITE(44,'(i4,2x,i4,2x,i8)') ndxa-nred,ndya,jvft44
      write(44,'(i4,2x,i4,2x,i4)') natmi,nmoli,nioni
      !cank 960623
      nnatmi=natmi
      nnmoli=nmoli
      nnplsi=nplsi
      nns=nnplsi
      !cank
      do iatm=1,natmi
        write(44,*) texts(iatm+nsph)
      end do
      do imol=1,nmoli
        write(44,*) texts(imol+nspa)
      end do
      do iion=1,nioni
        write(44,*) texts(iion+nspam)
      end do
      call eirene_neutr(44,ndxa-nred,ndya,natmi,dab2,ndx,ndy,natm,1,1)
      call eirene_neutr(44,ndxa-nred,ndya,natmi,tab2,ndx,ndy,natm,1,1)
      call eirene_neutr(44,ndxa-nred,ndya,nmoli,dmb2,ndx,ndy,nmol,1,1)
      call eirene_neutr(44,ndxa-nred,ndya,nmoli,tmb2,ndx,ndy,nmol,1,1)
      call eirene_neutr(44,ndxa-nred,ndya,nioni,dib2,ndx,ndy,nion,1,1)
      call eirene_neutr(44,ndxa-nred,ndya,nioni,tib2,ndx,ndy,nion,1,1)
      call eirene_neutr(44,ndxa-nred,ndya,natmi,rfluxa,ndx,ndy,natm,1,1)
      call eirene_neutr(44,ndxa-nred,ndya,nmoli,rfluxm,ndx,ndy,nmol,1,1)
      call eirene_neutr(44,ndxa-nred,ndya,natmi,pfluxa,ndx,ndy,natm,1,1)
      call eirene_neutr(44,ndxa-nred,ndya,nmoli,pfluxm,ndx,ndy,nmol,1,1)
      call eirene_neutr(44,ndxa-nred,ndya,natmi,
     .                     refluxa,ndx,ndy,natm,1,1)
      call eirene_neutr(44,ndxa-nred,ndya,nmoli,
     .                     refluxm,ndx,ndy,nmol,1,1)
      call eirene_neutr(44,ndxa-nred,ndya,natmi,
     .                     pefluxa,ndx,ndy,natm,1,1)
      call eirene_neutr(44,ndxa-nred,ndya,nmoli,
     .                     pefluxm,ndx,ndy,nmol,1,1)
      call eirene_neutr(44,ndxa-nred,ndya,1,emiss,ndx,ndy,1,1,1)
      call eirene_neutr(44,ndxa-nred,ndya,1,emissmol,ndx,ndy,1,1,1)
      !cank 960511
      !c*** save the molecule-related sources...
      call eirene_neutr(44,ndxa-nred,ndya,nmoli,srcml,ndx,ndy,nmol,1,1)
      call eirene_neutr(44,ndxa-nred,ndya,nmoli,
     .                     edissml,ndx,ndy,nmol,1,1)
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
        do istra=1,nstrai
          call neutrs(44,wldnek(1,istra),1)
          call neutrs(44,wldnep(1,istra),1)
          call neutrs(44,wldna(1,1,istra),natmi)
          call neutrs(44,ewlda(1,1,istra),natmi)
          call neutrs(44,wldnm(1,1,istra),nmoli)
          call neutrs(44,ewldm(1,1,istra),nmoli)
          call neutrs(44,wldra(1,1,istra),natmi)
          call neutrs(44,wldrm(1,1,istra),nmoli)
        end do
      end if
      !c*** from 961228 on:
      call neutrs(44,wldpp,nplsi)
      call neutrs(44,wldpa,natmi)
      call neutrs(44,wldpm,nmoli)
      call neutrs(44,wldpeb,1)
      call neutrs(44,wldspt,1)
      if(nstrai.gt.1) then
        do istra=1,nstrai
          call neutrs(44,wldpp(1,1,istra),nplsi)
          call neutrs(44,wldpa(1,1,istra),natmi)
          call neutrs(44,wldpm(1,1,istra),nmoli)
          call neutrs(44,wldpeb(1,istra),1)
          call neutrs(44,wldspt(1,istra),1)
        end do
      end if
      !c*** from 20000727 on:
      write(44,'(18i4)') (isrftype(i),i=1,nnlimi)
      write(44,'(18i4)') (isrftype(nlim+i),i=1,nnstsi)
      !cank
      rewind (44)
!cc<<<
!c      print *
!c      print *,'%%% Eirene data in Eirene %%%',nlmpgs
!c      print *,'nnatmi,nnmoli,nnstsi,nnlimi = ',
!c     ,        nnatmi,nnmoli,nnstsi,nnlimi
!c      do k=0,nstrai !{
!c       print *
!c       print *,'wldna, wldnm, wldra, wldrm:',k
!c       print '(1x,20(6x,a2,i3.2))',
!c     ,        ('na',i,i=1,natmi),('nm',i,i=1,nmoli),
!c     ,        ('ra',i,i=1,natmi),('rm',i,i=1,nmoli)
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
!c             print '(1p,i4,20e11.4)',j,
!c     ,              (wldna(j,i,k),i=1,natmi),(wldnm(j,i,k),i=1,nmoli),
!c     ,              (wldra(j,i,k),i=1,natmi),(wldrm(j,i,k),i=1,nmoli)
!c           end if !}
!c         end if !}
!c       end do !}
!c
!c       print *,'wldpp:',k
!c       do j=nlim+1,nlim+nstsi !{
!c         print '(1p,i4,20e9.2)',j,(wldpp(j,i,k),i=1,nplsi)
!c       end do !}
!c       print *
!c       print *,'wldpeb, wldpa, wldpm:',k
!c       do j=nlim+1,nlim+nstsi !{
!c         print '(1p,i4,20e9.2)',j,wldpeb(j,k),
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
        real*8 :: dummy(nlimps,*)
        do iif=1,ldmf
          if(nlimi.gt.0) write(kard,'(5(e16.8))') 
     .                   (dummy(is,iif),is=1,nlimi)
          if(nstsi.gt.0) write(kard,'(5(e16.8))') 
     .                   (dummy(is+nlim,iif),is=1,nstsi)
        enddo
        end subroutine
      end subroutine

      subroutine eirene_extrab25_wneuclean
      implicit none
      if(allocated(dab2)) then
        dab2=0
        dmb2=0
        dib2=0
        tab2=0
        tmb2=0
        tib2=0
        rfluxa=0
        rfluxm=0
        refluxa=0
        refluxm=0
        pfluxa=0
        pfluxm=0
        pefluxa=0
        pefluxm=0
        emiss=0
        emissmol=0
        srcml=0
        edissml=0
        wldnek=0
        wldnep=0
        wldna=0
        ewlda=0
        wldnm=0
        ewldm=0
        wldra=0
        wldrm=0
        wldpp=0
        wldpa=0
        wldpm=0
        wldpeb=0
        wldspt=0
        eneutrad=0
      endif
      return
      end subroutine


      subroutine eirene_extrab25_iniusr_init(n_spcsrf,l_spcsrf,
     .          i_spcsrf,
     .          j_spcsrf,sps_sgrp,sps_absr,sps_trno,sps_trni,
     .          sps_mtri,sps_tmpr,sps_spph,sps_spch,
     .          sps_mtrl,sps_id)
      implicit none
      integer, intent(in) :: n_spcsrf,l_spcsrf(:),i_spcsrf(:),
     .      j_spcsrf(:),sps_sgrp(:)
      real*8, intent(in) :: sps_absr(:),sps_trno(:),sps_trni(:),
     .      sps_mtri(:), sps_tmpr(:), sps_spph(:), sps_spch(:)
      character*8, intent(in) :: sps_mtrl(:),sps_id(:)

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
      end subroutine

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
      endif
      ifirst_wneutral=0

      if(allocated(flux_save)) then
        deallocate(flux_save)
      endif

      if(allocated(plnxtri)) then
        deallocate(plnxtri, plnytri, pplnxtri, pplnytri)
      endif
      end subroutine

      end module


C
*//GEOMD//
C=======================================================================
C          S U B R O U T I N E   G E O M D
C=======================================================================
      SUBROUTINE EIRENE_GEOMD(NDXA,NDYA,NPLP,NR1ST,
     .                 PUX,PUY,PVX,PVY,itype)
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CCONA
      USE EIRMOD_CGEOM
      IMPLICIT NONE
C
      REAL(DP), INTENT(OUT) :: PUX(*),PUY(*),PVX(*),PVY(*)
      INTEGER, INTENT(INOUT) :: NDXA,NDYA,NPLP,NR1ST
      INTEGER, INTENT(IN) :: ITYPE

      IF (ITYPE == 1) THEN
        CALL EIRENE_GEOMD_SONNET(NDXA,NDYA,NPLP,NR1ST,
     .                 PUX,PUY,PVX,PVY)
C
      ELSEIF (ITYPE == 2) THEN
        CALL EIRENE_GEOMD_CARRE(NDXA,NDYA,NPLP,NR1ST,
     .                 PUX,PUY,PVX,PVY)

      ELSE
        CALL EIRENE_GEOMD_LINDA(NDXA,NDYA,NPLP,NR1ST,
     .                 PUX,PUY,PVX,PVY)
      END IF

      RETURN
*//END GEOMD//
      END

C
*//GEOMD//
C=======================================================================
C          S U B R O U T I N E   G E O M D _ C A R R E
C=======================================================================
      SUBROUTINE EIRENE_GEOMD_CARRE(NDXA,NDYA,NPLP,NR1ST,
     .                 PUX,PUY,PVX,PVY)
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CCONA
      USE EIRMOD_CGEOM
      IMPLICIT NONE
C
      REAL(DP), INTENT(OUT) :: PUX(*),PUY(*),PVX(*),PVY(*)
      INTEGER, INTENT(INOUT) :: NDXA,NDYA,NPLP,NR1ST

      character(110) :: zeile
      REAL(DP) :: br(0:ndxp,0:ndyp,4),bz(0:ndxp,0:ndyp,4)
C
C  GEOMETRY DATA: CELL VERTICES (LINDA ---> EIRENE)
      REAL(DP) ::
     R  X1(NDX),Y1(NDX),X2(NDX),Y2(NDX),X3(NDX),Y3(NDX),
     R  X4(NDX),Y4(NDX)
      REAL(DP) :: DX, DY
      INTEGER :: I0, I0E, IX, IY, I1, I2, I3, I4, IPART, I, J, NCUT
      CHARACTER(80) :: LINE
      LOGICAL :: LRDCUT

      ndxa=0
      ndya=0
      lrdcut = .false.

!pb      DO I = 1, 4
      DO 
        read (30,'(A80)') line
        if (.not.lrdcut) then
          call eirene_uppercase (line)
          i0 = index(line,'NCUT')
          if (i0 > 0) then
            i0 = index(line,'=') + 1
            i1 = i0 + verify(line(i0:),' ')-1
            i2 = i1 + scan(line(i1+1:),' ')-1
            read (line(i1:i2),*) ncut
            read (30,'(A80)') line
            i1 = index(line,'=')
            read (line(i1+1:),*) (npoint(2,j),j=1,ncut)
            do j=1,ncut
              npoint(2,j) = npoint(2,j) + j
              npoint(1,j+1) = npoint(2,j) + 1
            end do
            npoint(1,1) = 1
            ipart = ncut+1
            lrdcut = .true.
          end if 
        end if
        if (index(line,'=======') /= 0) exit
      END DO

1     continue
      read (30,'(a110)',end=99) zeile
      i0=index(zeile,'(')
      i0e=index(zeile,')')
      read (zeile(i0+1:i0e-1),*) ix,iy
      ndxa=max(ndxa,ix)
      ndya=max(ndya,iy)
      i1=index(zeile,': (')
      i2=index(zeile(i1+3:),')')+i1+2
      read (zeile(i1+3:i2-1),*) br(ix,iy,4),bz(ix,iy,4)
      i3=index(zeile(i2+1:),'(')+i2
      i4=i3+index(zeile(i3+1:),')')
      read (zeile(i3+1:i4-1),*) br(ix,iy,3),bz(ix,iy,3)

      read (30,'(a110)') zeile

      read (30,'(a110)') zeile
      i1=index(zeile,'(')
      i2=index(zeile,')')
      read (zeile(i1+1:i2-1),*) br(ix,iy,1),bz(ix,iy,1)
      i3=i2+index(zeile(i2+1:),'(')
      i4=i2+index(zeile(i2+1:),')')
      read (zeile(i3+1:i4-1),*) br(ix,iy,2),bz(ix,iy,2)

      read (30,*)
      goto 1


99    continue
      ndxa=ndxa-1
      ndya=ndya-1
C
!pb      DO 1015 IY=1,NDYA
!pb        DO 1014 IX=1,NDXA
!pb          X1(IX)=br(ix,iy,1)
!pb          Y1(IX)=bz(ix,iy,1)
!pb          X2(IX)=br(ix,iy,2)
!pb          Y2(IX)=bz(ix,iy,2)
!pb          X3(IX)=br(ix,iy,4)
!pb          Y3(IX)=bz(ix,iy,4)
!pb          X4(IX)=br(ix,iy,3)
!pb          Y4(IX)=bz(ix,iy,3)
!pb1014    CONTINUE
!pb        CALL MSHPROJ (X1,Y1,X2,Y2,X3,Y3,X4,Y4,PUX,PUY,PVX,PVY,NDXA,
!pb     .                NR1ST,IY)
!pb1015  CONTINUE
C
C SEARCH FOR THE CUTS
C
      IF (.NOT.LRDCUT) THEN
        IPART=1
        NPOINT(1,IPART)=1
        IY=1
        DO IX=1,NDXA
          DX=BR(IX+1,IY,1)-BR(IX,IY,2)
          DY=BZ(IX+1,IY,1)-BZ(IX,IY,2)
          IF (DX*DX+DY*DY.GT.EPS10) THEN
C CUT GEFUNDEN
            NPOINT(2,IPART)=IX+1+IPART-1
            IPART=IPART+1
            NPOINT(1,IPART)=IX+1+IPART-1
          ENDIF
        ENDDO
      END IF
      NPOINT(2,IPART)=NDXA+IPART
C
      NPLP=IPART
C
      DO IY=1,NDYA
        DO IPART=1,NPLP
          DO IX=NPOINT(1,IPART),NPOINT(2,IPART)-1
            XPOL(IY,IX)=BR(IX-(IPART-1),IY,1)
            YPOL(IY,IX)=BZ(IX-(IPART-1),IY,1)
          ENDDO
          XPOL(IY,NPOINT(2,IPART))=BR(NPOINT(2,IPART)-IPART,IY,2)
          YPOL(IY,NPOINT(2,IPART))=BZ(NPOINT(2,IPART)-IPART,IY,2)
        ENDDO
      ENDDO
C INTRODUCE OUTERMOST RADIAL POLYGON
      DO IPART=1,NPLP
        DO IX=NPOINT(1,IPART),NPOINT(2,IPART)-1
          XPOL(NDYA+1,IX)=BR(IX-(IPART-1),NDYA,4)
          YPOL(NDYA+1,IX)=BZ(IX-(IPART-1),NDYA,4)
        ENDDO
        XPOL(NDYA+1,NPOINT(2,IPART))=BR(NPOINT(2,IPART)-IPART,NDYA,3)
        YPOL(NDYA+1,NPOINT(2,IPART))=BZ(NPOINT(2,IPART)-IPART,NDYA,3)
      ENDDO
C
      DO J=1,NDYA+1
        DO I=1,NPOINT(2,NPLP)
          XPOL(J,I)=XPOL(J,I)*100.
          YPOL(J,I)=YPOL(J,I)*100.
        END DO
      END DO
C
      ndxa=npoint(2,nplp)-1

      DO 1015 IY=1,NDYA
        DO 1014 IX=1,NDXA
          X1(IX)=XPOL(IY,IX)
          Y1(IX)=YPOL(IY,IX)
          X2(IX)=XPOL(IY,IX+1)
          Y2(IX)=YPOL(IY,IX+1)
          X3(IX)=XPOL(IY+1,IX)
          Y3(IX)=YPOL(IY+1,IX)
          X4(IX)=XPOL(IY+1,IX+1)
          Y4(IX)=YPOL(IY+1,IX+1)
1014    CONTINUE
C
        CALL EIRENE_MSHPROJ (X1,Y1,X2,Y2,X3,Y3,X4,Y4,PUX,PUY,PVX,PVY,
     .                NDXA,NR1ST,IY)
1015  CONTINUE
c
C     do j=1,ndya+1
C       write (iunout,*)
C       write (iunout,*) 'in geomd polygon ',j
C       write (iunout,'(1p,6e12.4)') 
C    .        (xpol(j,i),ypol(j,i),i=1,npoint(2,nplp))
C     enddo
C
      RETURN
*//END GEOMD_CARRE//
      END

*//GEOMD//
C=======================================================================
C          S U B R O U T I N E   G E O M D _ L I N D A
C=======================================================================
      SUBROUTINE EIRENE_GEOMD_LINDA(NDXA,NDYA,NPLP,NR1ST,
     .                        PUX,PUY,PVX,PVY)
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CCONA
      USE EIRMOD_CPOLYG
      USE EIRMOD_CGEOM
      USE EIRMOD_COMPRT,ONLY:IUNIN,IUNOUT !VK
      IMPLICIT NONE
C
      INTEGER, INTENT(INOUT) :: NDXA, NDYA, NPLP
      INTEGER, INTENT(INOUT) :: NR1ST
      REAL(DP), INTENT(OUT) :: PUX(*),PUY(*),PVX(*),PVY(*)
C
C  GEOMETRY DATA: CELL VERTICES (LINDA ---> EIRENE)
      REAL(DP) :: 
     R  X1(NDX),Y1(NDX),X2(NDX),Y2(NDX),X3(NDX),Y3(NDX),
     R  X4(NDX),Y4(NDX)
C
      CHARACTER(80) :: LINE,LINE2
C   DIMENSIONIERUNG FUER GITTER
      INTEGER :: DIMXH,DIMYH,NNCUT,NNISO,
     1 NXCUT1(10),NXCUT2(10),NYCUT1(10),NYCUT2(10),
     2 NXISO1(10),NXISO2(10),NYISO1(10),NYISO2(10)
      INTEGER :: IX, IY, I, J, NP, NWISO

      REAL(DP) :: DUMMI(3)
      REAL(DP) :: MERK(NDY)
C  ACTUAL MESH USED IN THIS RUN
C
C      EINLESEROUTINE ANGEPASST AUF BRAAMS-OUTPUT
C   GEAENDERTE DIMENSIONIERUNG BZW. CUT-POSITION
C        MUSS PER HAND ANGEPASST WERDEN:
C        PARAMETER DIMXH,DIMYH                        RFS 14.5.1991
       NXCUT1=0
       NXCUT2=0
       NYCUT1=0
       NYCUT2=0
       NXISO1=0
       NXISO2=0
       NYISO1=0
       NYISO2=0
       nniso = -1 !pb
       OPEN (UNIT=30,ACCESS='SEQUENTIAL',FORM='FORMATTED',ERR=100) !VK
      REWIND 30
3366  FORMAT(/)
      read(30,*)
      do
        read (30,'(A80)') LINE
        i = verify(line,' ')
        if (i /= 0) then
           backspace 30
           exit
        end if
      end do

      read(30,*,ERR=100,END=100) dimxh,dimyh,nncut
      if(nncut.gt.10) stop 'Increase array sizes for cut'
      read(30,*,ERR=100,END=100) 
     r     (nxcut1(i),nxcut2(i),nycut1(i),nycut2(i),i=1,nncut)
      if (nncut.gt.2) then
         read(30,*,ERR=100,END=100) nniso
         if(nniso.gt.10) stop 'Increase array sizes for insulating cut'
         read(30,*,ERR=100,END=100) 
     r            (nxiso1(i),nxiso2(i),nyiso1(i),nyiso2(i),i=1,nniso)
      ELSE
       NNISO=0 !VK
      endif
      read(30,*,ERR=100,END=100)
C    ANZAHL DER TEILSTUECKE PRO POLYGON
      NPLP = MAX((NNCUT/2)*3,1)
C    COMPUTE WIDTH OF INSULATING CUT FOR DOUBLE NULL
      NWISO=NXISO2(1)-NXISO1(1)     
C     PRINT MESSAGE AND CHECK
      WRITE(IUNOUT,*) "GEOMD: NNCUT, NNISO, NPLP, NWISO ", 
     w                        NNCUT, NNISO, NPLP, NWISO 
      IF(NNCUT.NE.0.AND.NNCUT.NE.2.AND.NNCUT.NE.4) THEN
        WRITE(IUNOUT,*) "WARNING FROM GEOMD: UNKNOWN TOPOLOGY"
        WRITE(IUNOUT,*) " NNCUT ",NNCUT
      END IF      
      IF(NNCUT.EQ.4.AND.NNISO.NE.1) THEN
        WRITE(IUNOUT,*) "WARNING FROM GEOMD: UNKNOWN TOPOLOGY"
        WRITE(IUNOUT,*) " NNCUT, NNISO ",NNCUT,NNISO
      END IF      
      IF(NNCUT.EQ.2) THEN
       IF(NXCUT1(1).NE.NXCUT2(2)-1.OR.
     .    NXCUT1(2).NE.NXCUT2(1)-1) THEN
        WRITE(IUNOUT,*) "WARNING FROM GEOMD: CUTS DO NOT MATCH"
        WRITE(IUNOUT,*) " NXCUT1(1), NXCUT2(2) ",NXCUT1(1),NXCUT2(2)
        WRITE(IUNOUT,*) " NXCUT1(2), NXCUT2(2) ",NXCUT1(1),NXCUT2(2)
       END IF
      END IF
      IF(NNCUT.EQ.4) THEN
       IF(NXCUT1(1).NE.NXCUT2(NNCUT)-1.OR.
     .    NXCUT1(2).NE.NXCUT2(3)-1.OR.
     .    NXCUT1(3).NE.NXCUT2(2)-1.OR.
     .    NXCUT1(NNCUT).NE.NXCUT2(1)-1) THEN
        WRITE(IUNOUT,*) "WARNING FROM GEOMD: CUTS DO NOT MATCH"
        WRITE(IUNOUT,*) " NXCUT1(1), NXCUT2(4) ",NXCUT1(1),NXCUT2(4)
        WRITE(IUNOUT,*) " NXCUT1(2), NXCUT2(3) ",NXCUT1(2),NXCUT2(3)
        WRITE(IUNOUT,*) " NXCUT1(3), NXCUT2(2) ",NXCUT1(3),NXCUT2(2)
        WRITE(IUNOUT,*) " NXCUT1(4), NXCUT2(1) ",NXCUT1(4),NXCUT2(1)
       END IF
      END IF

C    READING OF POLYGON DATA
      DO 10 IX = 1, DIMXH
       IF (IX.LE.nxcut1(1)-1) THEN
        DO 12 IY = 1, DIMYH
         IF (IY.LE.dimyh-1) THEN
          READ (30,*,ERR=100,END=100) DUMMI(1),
     .                  DUMMI(2),DUMMI(3),XPOL(IY,IX)
          READ (30,*,ERR=100,END=100) DUMMI(1),
     .                  DUMMI(2),DUMMI(3),YPOL(IY,IX)
         ENDIF
         IF (IY.EQ.dimyh) THEN
          READ (30,*,ERR=100,END=100) DUMMI(1),
     .                   DUMMI(2),XPOL(dimyh+1,IX),XPOL(IY,IX)
          READ (30,*,ERR=100,END=100) DUMMI(1),
     .                   DUMMI(2),YPOL(dimyh+1,IX),YPOL(IY,IX)
         ENDIF
12      CONTINUE
       ENDIF
       IF (IX.EQ.nxcut1(1)) THEN
        DO 14 IY = 1, DIMYH
         IF (IY.LE.dimyh-1) THEN
         READ (30,*,ERR=100,END=100) XPOL(IY,nxcut2(2)),DUMMI(1),
     .                  DUMMI(2),XPOL(IY,IX)
         READ (30,*,ERR=100,END=100) YPOL(IY,nxcut2(2)),DUMMI(1),
     .                  DUMMI(2),YPOL(IY,IX)
         ENDIF
         IF (IY.EQ.dimyh) THEN
          READ (30,*,ERR=100,END=100) 
     .                  XPOL(IY,nxcut2(2)),XPOL(dimyh+1,nxcut2(2)),
     .                                XPOL(dimyh+1,IX),XPOL(IY,IX)
          READ (30,*,ERR=100,END=100) 
     .                  YPOL(IY,nxcut2(2)),YPOL(dimyh+1,nxcut2(2)),
     .                                YPOL(dimyh+1,IX),YPOL(IY,IX)
         ENDIF
14      CONTINUE
       ENDIF
       IF ((IX.GE.nxcut2(2)).AND.(IX.LE.nxcut1(2)-1)) THEN
        DO 16 IY = 1, DIMYH
         IF (IY.LE.dimyh-1) THEN
          READ (30,*,ERR=100,END=100) DUMMI(1),
     .                   DUMMI(2),DUMMI(3),XPOL(IY,IX+1)
          READ (30,*,ERR=100,END=100) DUMMI(1),
     .                   DUMMI(2),DUMMI(3),YPOL(IY,IX+1)
         ENDIF
         IF (IY.EQ.dimyh) THEN
          READ (30,*,ERR=100,END=100) DUMMI(1),DUMMI(2),
     .                   XPOL(dimyh+1,IX+1),XPOL(IY,IX+1)
          READ (30,*,ERR=100,END=100) DUMMI(1),DUMMI(2),
     .                   YPOL(dimyh+1,IX+1),YPOL(IY,IX+1)
         ENDIF
16      CONTINUE
       ENDIF
       IF (IX.EQ.nxcut1(2)) THEN
        DO 18 IY = 1, DIMYH
         IF (IY.LE.dimyh-1) THEN
         READ (30,*,ERR=100,END=100) XPOL(IY,nxcut2(1)+1),DUMMI(1),
     .                  DUMMI(2),XPOL(IY,IX+1)
         READ (30,*,ERR=100,END=100) YPOL(IY,nxcut2(1)+1),DUMMI(1),
     .                  DUMMI(2),YPOL(IY,IX+1)
         ENDIF
         IF (IY.EQ.dimyh) THEN
          READ (30,*,ERR=100,END=100) 
     .                  XPOL(IY,nxcut2(1)+1),XPOL(dimyh+1,nxcut2(1)+1),
     .                                XPOL(dimyh+1,IX+1),XPOL(IY,IX+1)
          READ (30,*,ERR=100,END=100) 
     .                  YPOL(IY,nxcut2(1)+1),YPOL(dimyh+1,nxcut2(1)+1),
     .                                YPOL(dimyh+1,IX+1),YPOL(IY,IX+1)
         ENDIF
18      CONTINUE
       ENDIF
       IF ((IX.GE.nxcut2(1)).AND.(IX.LE.dimxh-1)) THEN
        DO 22 IY = 1, DIMYH
         IF (IY.LE.dimyh-1) THEN
          READ (30,*,ERR=100,END=100) DUMMI(1),DUMMI(2),DUMMI(3),
     .                   XPOL(IY,IX+2)
          READ (30,*,ERR=100,END=100) DUMMI(1),DUMMI(2),DUMMI(3),
     .                   YPOL(IY,IX+2)
         ENDIF
         IF (IY.EQ.dimyh) THEN
          READ (30,*,ERR=100,END=100) DUMMI(1),DUMMI(2),
     .                   XPOL(dimyh+1,IX+2),XPOL(IY,IX+2)
          READ (30,*,ERR=100,END=100) DUMMI(1),DUMMI(2),
     .                   YPOL(dimyh+1,IX+2),YPOL(IY,IX+2)
         ENDIF
22      CONTINUE
       ENDIF
       IF (IX.EQ.dimxh) THEN
        DO 24 IY = 1, DIMYH
         IF (IY.LE.dimyh-1) THEN
         READ (30,*,ERR=100,END=100) XPOL(IY,dimxh+3),DUMMI(1),DUMMI(2),
     .                  XPOL(IY,IX+2)
         READ (30,*,ERR=100,END=100) YPOL(IY,dimxh+3),DUMMI(1),DUMMI(2),
     .                  YPOL(IY,IX+2)
         ENDIF
         IF (IY.EQ.dimyh) THEN
          READ (30,*,ERR=100,END=100) 
     .                  XPOL(IY,dimxh+3),XPOL(dimyh+1,dimxh+3),
     .                                XPOL(dimyh+1,IX+2),XPOL(IY,IX+2)
          READ (30,*,ERR=100,END=100) 
     .                  YPOL(IY,dimxh+3),YPOL(dimyh+1,dimxh+3),
     .                                YPOL(dimyh+1,IX+2),YPOL(IY,IX+2)
         ENDIF
24      CONTINUE
       ENDIF
10    CONTINUE

3333  FORMAT(4E15.7)

C   ANFANGSPUNKT DES ERSTEN TEILSTUECKS DES I-TEN POLYGONS
      NPOINT(1,1)=1
C   ENDPUNKT DES ERSTEN TEILSTUECKS DES I-TEN POLYGONS
      NPOINT(2,1)=nxcut1(1)+1
      IF (NNCUT.EQ.0) NPOINT(2,1)=dimxh+1
C   ANFANGSPUNKT DES ZWEITEN TEILSTUECKS DES I-TEN POLYGONS
      NPOINT(1,2)=nxcut2(nncut)+1
C   ENDPUNKT DES ZWEITEN TEILSTUECKS DES I-TEN POLYGONS
      NPOINT(2,2)=nxcut2(nncut-1)+1
C   ANFANGSPUNKT DES DRITTEN TEILSTUECKS DES I-TEN POLYGONS
      NPOINT(1,3)=nxcut2(nncut-1)+2
C   ENDPUNKT DES DRITTEN TEILSTUECKS DES I-TEN POLYGONS
      IF (NNCUT.EQ.2) NPOINT(2,3)=dimxh+3
      IF (NNCUT.EQ.4) THEN
       NPOINT(2,3)=NXISO1(1)+3
C   ANFANGSPUNKT DES VIERTEN TEILSTUECKS DES I-TEN POLYGONS
       NPOINT(1,4)=nxiso2(1)+4-NWISO
C   ENDPUNKT DES VIERTEN TEILSTUECKS DES I-TEN POLYGONS
       NPOINT(2,4)=nxcut1(3)+4-NWISO
C   ANFANGSPUNKT DES FUNFTEN TEILSTUECKS DES I-TEN POLYGONS
       NPOINT(1,5)=nxcut2(2)+4-NWISO
C   ENDPUNKT DES FUNFTEN TEILSTUECKS DES I-TEN POLYGONS
       NPOINT(2,5)=nxcut1(4)+5-NWISO
C   ANFANGSPUNKT DES SECHSTEN TEILSTUECKS DES I-TEN POLYGONS
       NPOINT(1,6)=nxcut2(1)+5-NWISO
C   ENDPUNKT DES SECHSTEN TEILSTUECKS DES I-TEN POLYGONS
       NPOINT(2,6)=dimxh+6-NWISO
      END IF
C
      DO 1015 IY=1,NDYA
        DO 1014 IX=1,NDXA
          X1(IX)=XPOL(IY,IX)
          Y1(IX)=YPOL(IY,IX)
          X2(IX)=XPOL(IY,IX+1)
          Y2(IX)=YPOL(IY,IX+1)
          X3(IX)=XPOL(IY+1,IX)
          Y3(IX)=YPOL(IY+1,IX)
          X4(IX)=XPOL(IY+1,IX+1)
          Y4(IX)=YPOL(IY+1,IX+1)
1014    CONTINUE
C
        CALL EIRENE_MSHPROJ (X1,Y1,X2,Y2,X3,Y3,X4,Y4,
     .                       PUX,PUY,PVX,PVY,NDXA,
     .                       NR1ST,IY)
1015  CONTINUE
C
      NP=NPOINT(2,NPLP)
      DO 1020 J=1,NDYA+1
        DO 1020 I=1,NP
          XPOL(J,I)=XPOL(J,I)*100.
          YPOL(J,I)=YPOL(J,I)*100.
          IF (ABS(XPOL(J,I)).LT.5.D-5) XPOL(J,I)=0.
          IF (ABS(YPOL(J,I)).LT.5.D-5) YPOL(J,I)=0.
1020  CONTINUE
      RETURN

 100  WRITE(IUNOUT,*) "COULD NOT OPEN FORT.30. ",
     w                "SKIP READING THE B2 GEOMETRY" !VK

*//END GEOMD_LINDA//
      END


C
*//GEOMD//
C=======================================================================
C          S U B R O U T I N E   G E O M D _ S O N N E T
C=======================================================================
      SUBROUTINE EIRENE_GEOMD_SONNET(NDXA,NDYA,NPLP,NR1ST,
     .                 PUX,PUY,PVX,PVY)
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CCONA
      USE EIRMOD_CGEOM
      IMPLICIT NONE
C
      REAL(DP), INTENT(OUT) :: PUX(*),PUY(*),PVX(*),PVY(*)
      INTEGER, INTENT(INOUT) :: NDXA,NDYA,NPLP,NR1ST

      character(200) :: zeile
      REAL(DP) :: br(0:ndxp,0:ndyp,4),bz(0:ndxp,0:ndyp,4)
C
C  GEOMETRY DATA: CELL VERTICES (LINDA ---> EIRENE)
      REAL(DP) ::
     R  X1(NDX),Y1(NDX),X2(NDX),Y2(NDX),X3(NDX),Y3(NDX),
     R  X4(NDX),Y4(NDX)
      REAL(DP) :: DX, DY
      INTEGER :: I0, I0E, IX, IY, I1, I2, I3, I4, IPART, I, J, NCUT
      CHARACTER(200) :: LINE
      LOGICAL :: LRDCUT

      ndxa=0
      ndya=0
      lrdcut = .false.

!pb      DO I = 1, 4
      DO 
        read (30,'(A200)') line
        if (.not.lrdcut) then
          call eirene_uppercase (line)
          i0 = index(line,' CUT')
          if (i0 > 0) then
            i0 = i0 + 4
            i1 = i0 + verify(line(i0:),' ')-1
            i2 = i1 + scan(line(i1+1:),' ')-1
            read (line(i1:i2),*) ncut
            do j=1,ncut
              i1 = i2 + verify(line(i2+1:),' ')-1
              i2 = i1 + scan(line(i1+1:),' ')-1
              read (line(i1:i2),*) npoint(2,j)
              npoint(2,j) = npoint(2,j) + j-1
              npoint(1,j+1) = npoint(2,j) + 1
            end do
            npoint(1,1) = 1
            ipart = ncut+1
            lrdcut = .true.
          end if 
        end if
        if (index(line,'=======') /= 0) exit
      END DO

1     continue
      read (30,'(a200)',end=99) zeile
      i0=index(zeile,'(')
      i0e=index(zeile,')')
      read (zeile(i0+1:i0e-1),*) ix,iy
      ndxa=max(ndxa,ix)
      ndya=max(ndya,iy)
      i1=index(zeile,': (')
      i2=index(zeile(i1+3:),')')+i1+2
      read (zeile(i1+3:i2-1),*) br(ix,iy,4),bz(ix,iy,4)
      i3=index(zeile(i2+1:),'(')+i2
      i4=i3+index(zeile(i3+1:),')')
      read (zeile(i3+1:i4-1),*) br(ix,iy,3),bz(ix,iy,3)

      read (30,'(a200)') zeile

      read (30,'(a200)') zeile
      i1=index(zeile,'(')
      i2=index(zeile,')')
      read (zeile(i1+1:i2-1),*) br(ix,iy,1),bz(ix,iy,1)
      i3=i2+index(zeile(i2+1:),'(')
      i4=i2+index(zeile(i2+1:),')')
      read (zeile(i3+1:i4-1),*) br(ix,iy,2),bz(ix,iy,2)

      read (30,*)
      goto 1


99    continue
      ndxa=ndxa-1
      ndya=ndya-1
C
!pb      DO 1015 IY=1,NDYA
!pb        DO 1014 IX=1,NDXA
!pb          X1(IX)=br(ix,iy,1)
!pb          Y1(IX)=bz(ix,iy,1)
!pb          X2(IX)=br(ix,iy,2)
!pb          Y2(IX)=bz(ix,iy,2)
!pb          X3(IX)=br(ix,iy,4)
!pb          Y3(IX)=bz(ix,iy,4)
!pb          X4(IX)=br(ix,iy,3)
!pb          Y4(IX)=bz(ix,iy,3)
!pb1014    CONTINUE
!pb        CALL MSHPROJ (X1,Y1,X2,Y2,X3,Y3,X4,Y4,PUX,PUY,PVX,PVY,NDXA,
!pb     .                NR1ST,IY)
!pb1015  CONTINUE
C
C SEARCH FOR THE CUTS
C
      IF (.NOT.LRDCUT) THEN
        IPART=1
        NPOINT(1,IPART)=1
        IY=1
        DO IX=1,NDXA
          DX=BR(IX+1,IY,1)-BR(IX,IY,2)
          DY=BZ(IX+1,IY,1)-BZ(IX,IY,2)
          IF (DX*DX+DY*DY.GT.EPS10) THEN
C CUT GEFUNDEN
            NPOINT(2,IPART)=IX+1+IPART-1
            IPART=IPART+1
            NPOINT(1,IPART)=IX+1+IPART-1
          ENDIF
        ENDDO
      END IF
      NPOINT(2,IPART)=NDXA+IPART
C
      NPLP=IPART
C
      DO IY=1,NDYA
        DO IPART=1,NPLP
          DO IX=NPOINT(1,IPART),NPOINT(2,IPART)-1
            XPOL(IY,IX)=BR(IX-(IPART-1),IY,1)
            YPOL(IY,IX)=BZ(IX-(IPART-1),IY,1)
          ENDDO
          XPOL(IY,NPOINT(2,IPART))=BR(NPOINT(2,IPART)-IPART,IY,2)
          YPOL(IY,NPOINT(2,IPART))=BZ(NPOINT(2,IPART)-IPART,IY,2)
        ENDDO
      ENDDO
C INTRODUCE OUTERMOST RADIAL POLYGON
      DO IPART=1,NPLP
        DO IX=NPOINT(1,IPART),NPOINT(2,IPART)-1
          XPOL(NDYA+1,IX)=BR(IX-(IPART-1),NDYA,4)
          YPOL(NDYA+1,IX)=BZ(IX-(IPART-1),NDYA,4)
        ENDDO
        XPOL(NDYA+1,NPOINT(2,IPART))=BR(NPOINT(2,IPART)-IPART,NDYA,3)
        YPOL(NDYA+1,NPOINT(2,IPART))=BZ(NPOINT(2,IPART)-IPART,NDYA,3)
      ENDDO
C
      DO J=1,NDYA+1
        DO I=1,NPOINT(2,NPLP)
          XPOL(J,I)=XPOL(J,I)*100.
          YPOL(J,I)=YPOL(J,I)*100.
        END DO
      END DO
C
      ndxa=npoint(2,nplp)-1

      DO 1015 IY=1,NDYA
        DO 1014 IX=1,NDXA
          X1(IX)=XPOL(IY,IX)
          Y1(IX)=YPOL(IY,IX)
          X2(IX)=XPOL(IY,IX+1)
          Y2(IX)=YPOL(IY,IX+1)
          X3(IX)=XPOL(IY+1,IX)
          Y3(IX)=YPOL(IY+1,IX)
          X4(IX)=XPOL(IY+1,IX+1)
          Y4(IX)=YPOL(IY+1,IX+1)
1014    CONTINUE
C
        CALL EIRENE_MSHPROJ (X1,Y1,X2,Y2,X3,Y3,X4,Y4,PUX,PUY,PVX,PVY,
     .                NDXA,NR1ST,IY)
1015  CONTINUE
c
C     do j=1,ndya+1
C       write (iunout,*)
C       write (iunout,*) 'in geomd polygon ',j
C       write (iunout,'(1p,6e12.4)') 
C    .        (xpol(j,i),ypol(j,i),i=1,npoint(2,nplp))
C     enddo
C
      RETURN
*//END GEOMD_SONNET//
      END


C
C
      SUBROUTINE EIRENE_MSHPROJ(X1,Y1,X2,Y2,X3,Y3,X4,Y4,PUX,PUY,PVX,PVY,
     .                   NDXA,NR1ST,IY)

      USE EIRMOD_PRECISION
      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: X1(*), Y1(*), X2(*), Y2(*),
     .                      X3(*), Y3(*), X4(*), Y4(*)
      REAL(DP), INTENT(OUT) :: PUX(*), PUY(*), PVX(*), PVY(*)
      INTEGER, INTENT(IN) :: NDXA, NR1ST, IY
      REAL(DP) :: D12, D34, D13, D24, EPS60, PUPV, PVPV, DVX, DVY,
     .          DUX, DUY
      INTEGER :: IX, IN

      EPS60 = 1.E-60_DP
C
C
      DO 1 IX=1,NDXA
C
C  CALCULATE THE NORM OF THE VECTORS (POINT2-POINT1),....
C
        D12 = SQRT((X2(IX)-X1(IX))*(X2(IX)-X1(IX))+(Y2(IX)-Y1(IX))*
     .        (Y2(IX)-Y1(IX)))+EPS60
        D34 = SQRT((X4(IX)-X3(IX))*(X4(IX)-X3(IX))+(Y4(IX)-Y3(IX))*
     .        (Y4(IX)-Y3(IX)))+EPS60
        D13 = SQRT((X3(IX)-X1(IX))*(X3(IX)-X1(IX))+(Y3(IX)-Y1(IX))*
     .        (Y3(IX)-Y1(IX)))+EPS60
        D24 = SQRT((X4(IX)-X2(IX))*(X4(IX)-X2(IX))+(Y4(IX)-Y2(IX))*
     .        (Y4(IX)-Y2(IX)))+EPS60
C
C  CALCULATE THE BISSECTING VECTORS, BUT NOT NORMALISED YET
C
        DUX = (X2(IX)-X1(IX))/D12 + (X4(IX)-X3(IX))/D34
        DUY = (Y2(IX)-Y1(IX))/D12 + (Y4(IX)-Y3(IX))/D34
        DVX = (X3(IX)-X1(IX))/D13 + (X4(IX)-X2(IX))/D24
        DVY = (Y3(IX)-Y1(IX))/D13 + (Y4(IX)-Y2(IX))/D24
C
C  CALCULATE THE COMPONENTS OF THE TWO UNIT VECTOR (= PROJECTION RATE)
C
        IN=IY+(IX-1)*NR1ST
        PUX(IN) = DUX/(SQRT(DUX*DUX+DUY*DUY)+EPS60)
        PUY(IN) = DUY/(SQRT(DUX*DUX+DUY*DUY)+EPS60)
        PVX(IN) = DVX/(SQRT(DVX*DVX+DVY*DVY)+EPS60)
        PVY(IN) = DVY/(SQRT(DVX*DVX+DVY*DVY)+EPS60)
C
C  ORTHOGONORMALIZE, CONSERVE ORIENTATION (E.SCHMIDT)
C
        PUPV=PUX(IN)*PVX(IN)+PUY(IN)*PVY(IN)
        PVX(IN)=PVX(IN)-PUPV*PUX(IN)
        PVY(IN)=PVY(IN)-PUPV*PUY(IN)
        PVPV=SQRT(PVX(IN)*PVX(IN)+PVY(IN)*PVY(IN))+EPS60
        PVX(IN)=PVX(IN)/PVPV
        PVY(IN)=PVY(IN)/PVPV
C
1     CONTINUE
      RETURN
      END


C
C
      SUBROUTINE EIRENE_INDMAP(FIELD,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .                  NCUTB,NCUTL,NPOINT,NPPLG)
C
C     INDEX MAPPING FOR BRAAMS DATA FIELDS. DATA IN DUMMY ZONES
C     (CUTS OR BOUNDARY ZONES) MAY BE NEEDED AND THUS ARE KEPT
C     AND DUBLICATED IN CASE NCUTL GT NCUTB
C
C     NCUTB= NUMBER OF CELLS IN IX DIRECTION PER CUT IN BRAAMS
C     NCUTL= NUMBER OF CELLS IN IX DIRECTION PER CUT IN LINDA (AND
C            THUS ALSO IN EIRENE) GEOMETRY

      USE EIRMOD_PRECISION
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NPOINT(2,*)
      INTEGER, INTENT(IN) :: NDX, NDY, NFL, NDXA, NDYA, NFLA, NCUTB,
     .                       NCUTL, NPPLG
      REAL(DP), INTENT(INOUT) :: FIELD(0:NDX+1,0:NDY+1,NFL),
     .                         DUMMY(0:NDX+1,0:NDY+1)
      INTEGER :: IX, IPART, IY, IF, IENDD, INB, IINID, IINIV, IENDV
C
C  LOOP FOR THE SPECIES
C
      DO 500 IF=1,NFLA
C
C  INITIALIZE DUMMY
C
        DO 10 IY=0,NDY+1
          DO 10 IX=0,NDX+1
10          DUMMY(IX,IY)=FIELD(IX,IY,IF)
C
C
C      NDX DIRECTION: IX=0: NOT MODIFIED
C                     IX=I(CUT): USE CUT VALUE
C                     IX=I(LAST X ZONE): MOVE TO NDXA+1
C
C  NPOINT(1,1)=1
C  NPOINT(2,NPPLG)=NDXA+1
C
        IF (NCUTB.LT.0) GOTO 990
        DO 211 IPART = 1,NPPLG
C  "VALID REGION"
          IINIV= NPOINT(1,IPART)
          IENDV= NPOINT(2,IPART)-1
C  "CUT REGION" AND LAST X ZONE IX = NDXA+1
          IF (IPART.LT.NPPLG) THEN
            IINID= NPOINT(2,IPART)
            IENDD= NPOINT(1,IPART+1)-1
            IF (IENDD-IINID+1.NE.NCUTL) GOTO 991
          ELSE
            IINID= NDXA+1
            IENDD= NDXA+1
          ENDIF
          DO 212 IY=0,NDYA+1
            DO 213 IX = IINIV,IENDV
              INB=IX-(IPART-1)*(NCUTL-NCUTB)
              DUMMY(IX,IY)=FIELD(INB,IY,IF)
213         CONTINUE
            DUMMY(IINID,IY) = FIELD(INB+1,IY,IF)
            IF (IENDD.NE.IINID) DUMMY(IENDD,IY) = FIELD(INB+NCUTB,IY,IF)
212       CONTINUE
211     CONTINUE
        DO 220 IY=0,NDYA+1
          DO 220 IX=0,NDXA+1
            FIELD(IX,IY,IF)=DUMMY(IX,IY)
220     CONTINUE
C
500   CONTINUE
      RETURN
C
990   CONTINUE
      WRITE (iunout,*) 'ERROR IN SUBR. INDMAP: THIS SUBR. IS VALID ONLY'
      WRITE (iunout,*) 'NCUTB>=0 BUT NCUTB = ',NCUTB
      CALL EIRENE_EXIT_OWN(1)
991   WRITE (iunout,*) 
     .  'ERROR IN SUBR. INDMAP: INCONSISTENCY IN NUMBER OF '
      WRITE (iunout,*) 'ZONES PER CUT FROM LINDA GEOMETRY DETECTED.  '
      WRITE (iunout,*) 'NCUTL = ',NCUTL, ' IENDD-IINID+1 = ',
     .                  IENDD-IINID+1
      CALL EIRENE_EXIT_OWN(1)
      END


C
C
      SUBROUTINE EIRENE_INDMPI(FIELD,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .                  NCUTB,NCUTL,NPOINT,NPPLG,NSTR,ISTR)
C
C     INDEX MAPPING: INVERS TO SUBR. INDMAP
C
      USE EIRMOD_PRECISION
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NPOINT(2,*)
      INTEGER, INTENT(IN) :: NDX, NDY, NFL, NDXA, NDYA, NFLA, NCUTB,
     .                       NCUTL, NPPLG, NSTR, ISTR
      REAL(DP), INTENT(INOUT) :: FIELD(0:NDX+1,0:NDY+1,NFL,NSTR),
     .                         DUMMY(0:NDX+1,0:NDY+1)
      INTEGER :: IX, IY, IF, IENDD, IPART, INB, IINID, IINIV, IENDV
C
C  LOOP OVER THE SPECIES
C
      DO 500 IF=1,NFLA
C
C  INITIALIZE DUMMY
C
        DO 10 IY=0,NDY+1
          DO 10 IX=0,NDX+1
10          DUMMY(IX,IY)=0.
C
C
C      NDX DIRECTION
C
C  NPOINT(1,1)=1
C  NPOINT(2,NPPLG)=NDXA+1
C
        IF (NCUTB.LT.0) GOTO 990
        DO 211 IPART = 1,NPPLG
C  "VALID REGION"
          IINIV= NPOINT(1,IPART)
          IENDV= NPOINT(2,IPART)-1
C  "CUT REGION" AND LAST X ZONE IX = NDXA+1
          IF (IPART.LT.NPPLG) THEN
            IINID= NPOINT(2,IPART)
            IENDD= NPOINT(1,IPART+1)-1
            IF (IENDD-IINID+1.NE.NCUTL) GOTO 991
          ELSE
            IINID= NDXA+1
            IENDD= NDXA+1
          ENDIF
          DO 212 IY=0,NDYA+1
            DO 213 IX = IINIV,IENDV
              INB=IX-(IPART-1)*(NCUTL-NCUTB)
              DUMMY(INB,IY)=FIELD(IX,IY,IF,ISTR)
213         CONTINUE
            DUMMY(INB+1,IY)=FIELD(IINID,IY,IF,ISTR)
            IF (IENDD.NE.IINID)
     .          DUMMY(INB+NCUTB,IY)=FIELD(IENDD,IY,IF,ISTR)
212       CONTINUE
211     CONTINUE
        DO 220 IY=0,NDYA+1
          DO 220 IX=0,NDXA+1
            FIELD(IX,IY,IF,ISTR)=DUMMY(IX,IY)
220     CONTINUE
C
500   CONTINUE
      RETURN
C
990   CONTINUE
      WRITE (iunout,*) 'ERROR IN SUBR. INDMPI: THIS SUBR. IS VALID ONLY'
      WRITE (iunout,*) 'NCUTB>=0 BUT NCUTB = ',NCUTB
      CALL EIRENE_EXIT_OWN(1)
991   WRITE (iunout,*) 
     .  'ERROR IN SUBR. INDMPI: INCONSISTENCY IN NUMBER OF'
      WRITE (iunout,*) 'ZONES PER CUT FROM LINDA GEOMETRY DETECTED. '
      WRITE (iunout,*) 'NCUTL = ',NCUTL, ' IENDD-IINID+1 = ',
     .                  IENDD-IINID+1
      CALL EIRENE_EXIT_OWN(1)
      END


C
*//PLASM//
C=======================================================================
C          S U B R O U T I N E   P L A S M
C=======================================================================
      SUBROUTINE EIRENE_PLASM(KARD,NDIMX,NDIMY,NDIMF,N,M,NF,DUMMY)

      USE EIRMOD_PRECISION
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: KARD, NDIMX, NDIMY, NDIMF, N, M, NF
      REAL(DP), INTENT(INOUT) :: DUMMY(0:N+1,0:M+1,NF)
      INTEGER :: ND1, LIM, IF, III, IX, IY, i1, i2, i3
      character(50) :: form
      character(200) :: zeile

      form = repeat(' ',50)
      read (kard,'(a200)',END=500) zeile
      i1 = index(zeile,'.')
      i2 = scan(zeile,'E,e')
      i3 = index(zeile(i2+1:),' ')
      write (form,'(A4,i0,a1,i0,a2)') '(5(E',i2+i3-1,'.',i2-i1-1,'))'
      backspace kard

      ND1 = NDIMX + 2
      LIM = (ND1/5)*5 - 4
      DUMMY(0:N+1,0:M+1,NF)=0._DP
      DO    110  IF = 1,NDIMF
      DO    110  IY = 0,NDIMY+1
      DO    100  IX = 1,LIM,5
100     READ(KARD,FORM,END=500) (DUMMY(-1+IX-1+III,IY,IF),III = 1,5)
        IF( (LIM+4).EQ.ND1 )     GOTO 110
       READ(KARD,FORM,END=500) (DUMMY(-1+IX,IY,IF),IX = LIM+5,ND1)
110   CONTINUE
500   RETURN
*//END PLASM//
      END


C
C
*//NEUTR//
C=======================================================================
C          S U B R O U T I N E   N E U T R
C=======================================================================
      SUBROUTINE EIRENE_NEUTR(KARD,NDIMX,NDIMY,NDIMF,DUMMY,LDMX,LDMY,
     .                        LDMF,LDNS,IS)

      USE EIRMOD_PRECISION
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: KARD, NDIMX, NDIMY, NDIMF, LDMX, LDMY,
     .                       LDMF, LDNS, IS
      REAL(DP), INTENT(IN) :: DUMMY(0:LDMX+1,0:LDMY+1,LDMF,LDNS)
      INTEGER :: ND1, LIM, IX, IY, III, IF
C
      ND1 = NDIMX
      LIM = (ND1/5)*5 - 4
      DO  500  IF = 1,NDIMF
        DO  110  IY = 1,NDIMY
          DO  100  IX = 1,LIM,5
  100     WRITE(KARD,910) (DUMMY(IX-1+III,IY,IF,IS),III = 1,5)
          IF( (LIM+4).EQ.ND1 )   GOTO 110
          WRITE(KARD,910) (DUMMY(IX,IY,IF,IS),IX = LIM+5,ND1)
  110   CONTINUE
  500 CONTINUE
      RETURN
csw 21jul2011 changed format, exponents with 3 digits needed
csw  910 FORMAT(5(E16.8))
  910 FORMAT(5(ES16.7E3))
*//END NEUTR//
      END




      SUBROUTINE EIRENE_SAVE_TALLIES (ISTRAI)
C
C  SAVE EIRENE TALLIES, SCALE PER UNIT FLUX (AMP), ON COMMON BRASCL
C  WTOTP IS NEGATIVE IN EIRENE (SINK FOR IONS)
C  ALL STRATA WHICH ARE NOT SPECIFIED BY INPUT BLOCK 14 (FROM
C  PLASMA CODE DATA) ARE NOT RESCALED HERE
C
csw 30nov2011 added checks LPAPL,...
csw

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_BRASPOI
      USE EIRMOD_CCOUPL
      USE EIRMOD_COUTAU
      USE EIRMOD_COMUSR
      USE EIRMOD_CGRID
      USE EIRMOD_CESTIM

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ISTRAI
      REAL(DP) :: FLXI
      INTEGER :: IN, IATM, IMOL, IPLS, IION, ICPV

      TYPE(CELLSIM), POINTER :: CPSIM
      TYPE(CELLMUL), POINTER :: CPMUL

      IF (ISTRAI.LE.NTARGI.AND.WTOTP(0,ISTRAI).NE.0.) THEN
         FLXI=-1._DP/WTOTP(0,ISTRAI)
      ELSEIF (ISTRAI.LE.NTARGI.AND.WTOTP(0,ISTRAI).EQ.0.) THEN
         RETURN
      ELSEIF (ISTRAI.GT.NTARGI) THEN
         FLXI=1._DP
      ENDIF

      CALL EIRENE_FREE_SIMARR(ISTRAI)
      CALL EIRENE_FREE_MULARR(ISTRAI)

      DO IPLS=1,NPLSI
cdr  save volumetric sources for plasma species ipls: particle, momentum, ion energy
        DO IN=1,NSBOX_TAL
          IF (LPAPL) THEN
            IF (PAPL(IPLS,IN) .NE. 0.D0) THEN
!pb            ALLOCATE(CPMUL)
              CPMUL => EIRENE_NEW_MULARR()
              CPMUL%IART = IPLS
              CPMUL%ICM = IN
              CPMUL%VALUEM = PAPL(IPLS,IN)*FLXI
              CPMUL%NXTMUL => PAPLS(ISTRAI)%PMUL
              PAPLS(ISTRAI)%PMUL => CPMUL
            ENDIF
          ENDIF
          IF (LPMPL) THEN 
            IF (PMPL(IPLS,IN) .NE. 0.D0) THEN
!PB            ALLOCATE(CPMUL)
              CPMUL => EIRENE_NEW_MULARR()
              CPMUL%IART = IPLS
              CPMUL%ICM = IN
              CPMUL%VALUEM = PMPL(IPLS,IN)*FLXI
              CPMUL%NXTMUL => PMPLS(ISTRAI)%PMUL
              PMPLS(ISTRAI)%PMUL => CPMUL
            ENDIF
          ENDIF

          IF (LPIPL) THEN
            IF (PIPL(IPLS,IN) .NE. 0.D0) THEN
!PB            ALLOCATE(CPMUL)
              CPMUL => EIRENE_NEW_MULARR()
              CPMUL%IART = IPLS
              CPMUL%ICM = IN
              CPMUL%VALUEM = PIPL(IPLS,IN)*FLXI
              CPMUL%NXTMUL => PIPLS(ISTRAI)%PMUL
              PIPLS(ISTRAI)%PMUL => CPMUL
            ENDIF
          ENDIF

          IF (LEAPL) THEN 
          IF (EAPL(IPLS,IN) .NE. 0.D0) THEN
!PB          ALLOCATE(CPMUL)
            CPMUL => EIRENE_NEW_MULARR()
            CPMUL%IART = IPLS
            CPMUL%ICM = IN
            CPMUL%VALUEM = EAPL(IPLS,IN)*FLXI
            CPMUL%NXTMUL => EAPLS(ISTRAI)%PMUL
            EAPLS(ISTRAI)%PMUL => CPMUL
          ENDIF
          ENDIF
          IF (LEMPL) THEN 
          IF (EMPL(IPLS,IN) .NE. 0.D0) THEN
!PB          ALLOCATE(CPMUL)
            CPMUL => EIRENE_NEW_MULARR()
            CPMUL%IART = IPLS
            CPMUL%ICM = IN
            CPMUL%VALUEM = EMPL(IPLS,IN)*FLXI
            CPMUL%NXTMUL => EMPLS(ISTRAI)%PMUL
            EMPLS(ISTRAI)%PMUL => CPMUL
          ENDIF
          ENDIF
          IF (LEIPL) THEN 
          IF (EIPL(IPLS,IN) .NE. 0.D0) THEN
!PB          ALLOCATE(CPMUL)
            CPMUL => EIRENE_NEW_MULARR()
            CPMUL%IART = IPLS
            CPMUL%ICM = IN
            CPMUL%VALUEM = EIPL(IPLS,IN)*FLXI
            CPMUL%NXTMUL => EIPLS(ISTRAI)%PMUL
            EIPLS(ISTRAI)%PMUL => CPMUL
          ENDIF
          ENDIF

          IF(LMAPL) THEN
            IF (MAPL(IPLS,IN) .NE. 0.D0) THEN
!pb            ALLOCATE(CPMUL)
              CPMUL => EIRENE_NEW_MULARR()
              CPMUL%IART = IPLS
              CPMUL%ICM = IN
              CPMUL%VALUEM = MAPL(IPLS,IN)*FLXI
              CPMUL%NXTMUL => MAPLS(ISTRAI)%PMUL
              MAPLS(ISTRAI)%PMUL => CPMUL
            ENDIF
          ENDIF
          IF (LMMPL) THEN
            IF (MMPL(IPLS,IN) .NE. 0.D0) THEN
!PB            ALLOCATE(CPMUL)
              CPMUL => EIRENE_NEW_MULARR()
              CPMUL%IART = IPLS
              CPMUL%ICM = IN
              CPMUL%VALUEM = MMPL(IPLS,IN)*FLXI
              CPMUL%NXTMUL => MMPLS(ISTRAI)%PMUL
              MMPLS(ISTRAI)%PMUL => CPMUL
            ENDIF
          ENDIF
          IF (LMIPL) THEN
            IF (MIPL(IPLS,IN) .NE. 0.D0) THEN
!PB            ALLOCATE(CPMUL)
              CPMUL => EIRENE_NEW_MULARR()
              CPMUL%IART = IPLS
              CPMUL%ICM = IN
              CPMUL%VALUEM = MIPL(IPLS,IN)*FLXI
              CPMUL%NXTMUL => MIPLS(ISTRAI)%PMUL
              MIPLS(ISTRAI)%PMUL => CPMUL
            ENDIF
          ENDIF
          IF (LMPHPL) THEN
            IF (MPHPL(IPLS,IN) .NE. 0.D0) THEN
!PB            ALLOCATE(CPMUL)
              CPMUL => EIRENE_NEW_MULARR()
              CPMUL%IART = IPLS
              CPMUL%ICM = IN
              CPMUL%VALUEM = MPHPL(IPLS,IN)*FLXI
              CPMUL%NXTMUL => MPHPLS(ISTRAI)%PMUL
              MPHPLS(ISTRAI)%PMUL => CPMUL
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      DO IN=1,NSBOX_TAL
    	IF (LEAEL) THEN
          IF (EAEL(IN) .NE. 0.D0) THEN
!PB          ALLOCATE(CPSIM)
            CPSIM => EIRENE_NEW_SIMARR()
            CPSIM%ICS = IN
            CPSIM%VALUES = EAEL(IN)*FLXI
            CPSIM%NXTSIM => EAELS(ISTRAI)%PSIM
            EAELS(ISTRAI)%PSIM => CPSIM
          ENDIF
        ENDIF
        IF (LEMEL) THEN
          IF (EMEL(IN) .NE. 0.D0) THEN
!PB          ALLOCATE(CPSIM)
            CPSIM => EIRENE_NEW_SIMARR()
            CPSIM%ICS = IN
            CPSIM%VALUES = EMEL(IN)*FLXI
            CPSIM%NXTSIM => EMELS(ISTRAI)%PSIM
            EMELS(ISTRAI)%PSIM => CPSIM
          ENDIF
        ENDIF
        IF (LEIEL) THEN 
          IF (EIEL(IN) .NE. 0.D0) THEN
!PB           ALLOCATE(CPSIM)
            CPSIM => EIRENE_NEW_SIMARR()
            CPSIM%ICS = IN
            CPSIM%VALUES = EIEL(IN)*FLXI
            CPSIM%NXTSIM => EIELS(ISTRAI)%PSIM
            EIELS(ISTRAI)%PSIM => CPSIM
          ENDIF
        ENDIF

      ENDDO

      DO IATM=1,NATMI
        DO IN=1,NSBOX_TAL
          IF(LPDENA) THEN
            IF (PDENA(IATM,IN) .NE. 0.D0) THEN
!PB            ALLOCATE(CPMUL)
              CPMUL => EIRENE_NEW_MULARR()
              CPMUL%IART = IATM
              CPMUL%ICM = IN
              CPMUL%VALUEM = PDENA(IATM,IN)*FLXI
              CPMUL%NXTMUL => PDENAS(ISTRAI)%PMUL
              PDENAS(ISTRAI)%PMUL => CPMUL
            ENDIF
          ENDIF
          IF (LEDENA) THEN
            IF (EDENA(IATM,IN) .NE. 0.D0) THEN
!PB            ALLOCATE(CPMUL)
              CPMUL => EIRENE_NEW_MULARR()
              CPMUL%IART = IATM
              CPMUL%ICM = IN
              CPMUL%VALUEM = EDENA(IATM,IN)*FLXI
              CPMUL%NXTMUL => EDENAS(ISTRAI)%PMUL
              EDENAS(ISTRAI)%PMUL => CPMUL
            ENDIF
	      ENDIF
        ENDDO
      ENDDO

      DO IMOL=1,NMOLI
        DO IN=1,NSBOX_TAL
	      IF (LPDENM) THEN
            IF (PDENM(IMOL,IN) .NE. 0.D0) THEN
!PB            ALLOCATE(CPMUL)
              CPMUL => EIRENE_NEW_MULARR()
              CPMUL%IART = IMOL
              CPMUL%ICM = IN
              CPMUL%VALUEM = PDENM(IMOL,IN)*FLXI
              CPMUL%NXTMUL => PDENMS(ISTRAI)%PMUL
              PDENMS(ISTRAI)%PMUL => CPMUL
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      DO IION=1,NIONI
        DO IN=1,NSBOX_TAL
	      IF (LPDENI) THEN
            IF (PDENI(IION,IN) .NE. 0.D0) THEN
!PB            ALLOCATE(CPMUL)
              CPMUL => EIRENE_NEW_MULARR()
              CPMUL%IART = IION
              CPMUL%ICM = IN
              CPMUL%VALUEM = PDENI(IION,IN)*FLXI
              CPMUL%NXTMUL => PDENIS(ISTRAI)%PMUL
              PDENIS(ISTRAI)%PMUL => CPMUL
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      DO ICPV=1,NCPVI
        DO IN=1,NSBOX_TAL
	      IF (LCOPV) THEN
            IF (COPV(ICPV,IN) .NE. 0.D0) THEN
!PB            ALLOCATE(CPMUL)
              CPMUL => EIRENE_NEW_MULARR()
              CPMUL%IART = ICPV
              CPMUL%ICM = IN
              CPMUL%VALUEM = COPV(ICPV,IN)*FLXI
              CPMUL%NXTMUL => COPVS(ISTRAI)%PMUL
              COPVS(ISTRAI)%PMUL => CPMUL
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END

C
C
C
C

