

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

      ! eirdiag.h/eirdiag.f
      !c*** Volume data:
      !c***    srcml   :   power loss due to molecules, including
      !c***    edissml :   power loss due to molecule dissociation
      !c***    eneutrad:   power radiated due to neutrals (only atoms at the moment)
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
     .  eneutrad

      integer, save, public :: nnlimi,nnstsi,nnatmi,nnmoli,nnioni
      integer, save, public :: nnplsi,nns
      integer, save, allocatable, public :: isrftype(:)
      logical, save, public :: lhalpha=.false.,lvib=.false.

      ! wneutral globals
      real(DP), save :: DA31(0:8,0:8)
      real(DP), save :: DP31(0:8,0:8)
      real(DP), save :: DM31(0:8,0:8)
      real(DP), save :: DI31(0:8,0:8)
      real(DP), save :: DN31(0:8,0:8)
      real(DP), save :: RHMH2(0:8),RH2PH2(0:8,0:8)
      CHARACTER, save :: FILNAM*8,H123*4,REAC*9,CRC*3
      logical, save :: hlp_pr
      integer, save :: ia1,ia2,ia3,iindex,ifirst_wneutral=0
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
     
      INTERFACE
        SUBROUTINE EIRENE_SLREAC (IR,FILNAM,H123,REAC,CRC,
     .             RC1MIN, RC1MAX, FP1, JFEX1MN, JFEX1MX,
     .             RC2MIN, RC2MAX, FP2, JFEX2MN, JFEX2MX,
     .             ELNAME, IZ1, IROW_ESC, ICOL_ESC, POP_ESC)
        USE EIRMOD_PRECISION
        INTEGER,      INTENT(IN) :: IR, IZ1
        INTEGER,      INTENT(IN), OPTIONAL :: IROW_ESC, ICOL_ESC
        REAL(DP),     INTENT(IN), OPTIONAL :: POP_ESC       
        CHARACTER(8), INTENT(IN) :: FILNAM
        CHARACTER(4), INTENT(IN) :: H123
        CHARACTER(LEN=*), INTENT(IN) :: REAC, ELNAME
        CHARACTER(3), INTENT(IN) :: CRC
        INTEGER,  INTENT(IN OUT) :: JFEX1MN, JFEX1MX,JFEX2MN, JFEX2MX
        REAL(DP), INTENT(IN OUT) :: RC1MIN, RC1MAX, FP1(6),
     .                              RC2MIN, RC2MAX, FP2(6)
        END SUBROUTINE EIRENE_SLREAC
      END INTERFACE

      real(DP) :: dummy(0:ndxp,0:ndyp)
      !c*** label for fort.44 file
      integer, parameter  :: jvft44=20000727
      !c*** and dissociation energy of the hydrogen molecule
      real(DP), parameter  :: diss_pot_H2=4.48
      !c*** radiative transition prob. level 3-->2 (1/sec) for H-alpha calc.
      real(DP), parameter :: fac32=4.410e7
      !C*** ionization potentials
      integer,parameter :: npot=20
      real(DP),save :: pot_data(npot),pot
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


      integer :: ix,iy,ir,ierror,i,j,in
      integer :: l,k,nred
      real(DP) :: de,te,hlp,sigadd1,sigadd2,sigadd3,sigadd4,sigadd5
      real(DP) :: dej,tei,tef,def,
     .            powalf1,powalf2,powalf3,powalf4,powalf5
      real(DP) :: powalf,datm3,dpls3,dmol3,dion3,dnml3,sigadd
      real(DP) :: da,dpp,dm,di,dn,ratio2,ratio7
      integer :: istra_in,istra_save
      real(DP) :: rc1min,rc1max,fp1(6),rc2min,rc2max,fp2(6),vl
      integer :: jfex1mn,jfex1mx,jfex2mn,jfex2mx
      real(DP) :: value
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
          write(iunout,*) 'Using new SIGHA (941017)'
          IERROR=0
          rc1min=-huge(1.d0)
          rc1max= huge(1.d0)
          jfex1mn=0
          jfex1mx=0
          fp1=0.d0
          rc2min=-huge(1.d0)
          rc2max= huge(1.d0)
          jfex2mn=0
          jfex2mx=0
          fp2=0.d0
          !C
          !C  READ REDUCED POPULATION COEFFICIENT FOR HYDR. ATOMS FROM FILE AMJUEL
          !C  AND PUT THEM FROM CREAC(..,..,IR) ONTO DA,DPP,DM,DI, AND DN ARRAY
          !C
          IR=NREACI
          IF (IR+7.GT.NREAC) then
            WRITE (IUNOUT,*) 'FROM SUBROUTINE HALFA: '
            CALL EIRENE_MASPRM('NREAC',5,NREAC,'IR',2,IR+7,IERROR)
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
     .                    rc1min,rc1max,fp1,jfex1mn,jfex1mx,
     .                    rc2min,rc2max,fp2,jfex2mn,jfex2mx,'  ',0)
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
     .                    rc1min,rc1max,fp1,jfex1mn,jfex1mx,
     .                    rc2min,rc2max,fp2,jfex2mn,jfex2mx,'  ',0)
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
     .                    rc1min,rc1max,fp1,jfex1mn,jfex1mx,
     .                    rc2min,rc2max,fp2,jfex2mn,jfex2mx,'  ',0)
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
     .                    rc1min,rc1max,fp1,jfex1mn,jfex1mx,
     .                    rc2min,rc2max,fp2,jfex2mn,jfex2mx,'  ',0)
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
     .                    rc1min,rc1max,fp1,jfex1mn,jfex1mx,
     .                    rc2min,rc2max,fp2,jfex2mn,jfex2mx,'  ',0)
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
     .                    rc1min,rc1max,fp1,jfex1mn,jfex1mx,
     .                    rc2min,rc2max,fp2,jfex2mn,jfex2mx,'  ',0)
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
          !C  2.0C INCLUDES ION CONVERSION (CX) ON H2(V)
          !C  OLD VERSION (WITHOUT THIS CX) SHOULD BE RECOVERED BY
          !C  READING 2.0B INSTEAD, AND OMITTING THE H- CHANNEL 5.
          CALL eirene_slreac(IR,FILNAM,H123,REAC,CRC,
     .                    rc1min,rc1max,fp1,jfex1mn,jfex1mx,
     .                    rc2min,rc2max,fp2,jfex2mn,jfex2mx,'  ',0)
           do I=1,9
            do J=1,9
              !RH2PH2(I-1,J-1)=CREAC(I,J,NREACI+1)
              RH2PH2(I-1,J-1)=REACDAT(IR)%OTH%POLY%DBLPOL(I,J)
            end do
          end do
          !C
          write(iunout,*) 'NREAC,NREACI,IR     ',NREAC,NREACI,IR
          !write(iunout,*) 'NRCX,IRCX        ',NRCX,IRCX
          !write(iunout,*) 'NREL,IREL        ',NREL,IREL
          !write(iunout,*) 'NRII,IRII        ',NRII,IRII
          !write(iunout,*) 'NELI,NAELI       ',NELI,naeli
          !write(iunout,*) 'NREI,NMEII,NIEII ',NREI,nmeii,nieii
          !write(iunout,*) 'NREC,NIRCI,NPRCI ',NREC,nirci,nprci
      end if
      iindex=0
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
      eneutrad(:,:,1,istra) = 0.d0
      do ncell=1,ntrii
        ix=ixtri(ncell)
        iy=iytri(ncell)
        if(ix.gt.0) then
          if (leael) eneutrad(ix,iy,1,istra)=eneutrad(ix,iy,1,istra) 
     .                      +eael(ncell)*vol(ncell)
          if (leael) value=value+eael(ncell)*vol(ncell)
          do iatm=1,natmi
            if(nchara(iatm).le.npot) then
              pot=pot_data(nchara(iatm))
            else
              pot=0.0_dp
            endif
            if (lpaat) eneutrad(ix,iy,1,istra)=eneutrad(ix,iy,1,istra) 
     .                        -paat(iatm,ncell)*pot*vol(ncell)
            if (lpaat) value=value-paat(iatm,ncell)*pot*vol(ncell)
          enddo
        endif
      enddo
      !open(555,file='eneutrad.dat',form='formatted')
      !do ix=1,76
      !  do iy=1,28
      !    write(555,'(i6,i6,1x,1p,e13.6)') ix,iy,eneutrad(ix,iy,1,istra)
      !  enddo
      !enddo
      !close(555)
      !write(iunout,'(a,i6,1x,1p,e13.6)') 'DBG: ISTRA, ENEUTRAD',istra,value
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
            if (lpdena) dab2(ix,iy,iatm,1)=dab2(ix,iy,iatm,1)+ 
     .                          pdena(iatm,in)*1e6*vl
            if (laddv) then
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
             end if
            if (ledena) tab2(ix,iy,iatm,1)=tab2(ix,iy,iatm,1)+
     .                             edena(iatm,in)*vl
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

            if(iindex+ia1+natmi+imol.le.nadv .and. laddv) 
     .         refluxm(ix,iy,imol,1)=refluxm(ix,iy,imol,1)+ 
     .                       addv(iindex+ia1+natmi+imol,in)*1.0e4*elcha
            if(iindex+ia3+natmi+imol.le.nadv .and. laddv) 
     .         pefluxm(ix,iy,imol,1)=pefluxm(ix,iy,imol,1)+ 
     .                       addv(iindex+ia3+natmi+imol,in)*1.0e4*elcha
            if (ledenm) tmb2(ix,iy,imol,1)=tmb2(ix,iy,imol,1)+
     .                      edenm(imol,in)*vl
            if (lpmml) srcml(ix,iy,imol,1)=srcml(ix,iy,imol,1)+ 
     .                                     pmml(imol,in)*vl
            !c*** Potential energy source related to hydrogen molecules
            !c*** for determining the radiation (W per cell)
            if (ncharm(imol).eq.2 .and. lpmml) 
     .         edissml(ix,iy,imol,istra)=edissml(ix,iy,imol,istra)+ 
     .            pmml(imol,in)*vl*diss_pot_H2
          endif
        end do
      end do
      do iion=1,nioni
        do in=1,ntrii
          ix=ixtri(in)
          iy=iytri(in)
          if(ix.gt.0) then
            vl = vol(in)/volcel(ix,iy)
            if (lpdeni) dib2(ix,iy,iion,1)=dib2(ix,iy,iion,1)+ 
     .                          pdeni(iion,in)*1e6*vl
            if (ledeni) tib2(ix,iy,iion,1)=tib2(ix,iy,iion,1)+
     .                          edeni(iion,in)*vl
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
        !c*** hlp accumulates the power taken away with re-emitted particles
        hlp=0.
        do j=1,natmi
          if (leotat)  wldnek(i,istra)=wldnek(i,istra)+eotat(j,i)
          if (lerfpat) wldpeb(i,istra)=wldpeb(i,istra)+erfpat(j,i)
          hlp=hlp+erfaat(j,i)+erfmat(j,i)+erfiat(j,i)
          if(potat(j,i).gt.0.) then
            ewlda(i,j,istra)=eotat(j,i)/potat(j,i)
          else
            ewlda(i,j,istra)=0.
          end if
          if (lpotat) wldna(i,j,istra)=hlp_cnv*potat(j,i)
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
          if (leotio)  wldnek(i,istra)=wldnek(i,istra)+eotio(j,i)
          if (lerfpio) wldpeb(i,istra)=wldpeb(i,istra)+erfpio(j,i)
          hlp=hlp+erfaio(j,i)+erfmio(j,i)+erfiio(j,i)
        end do
        do j=1,nplsi
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
        do j=1,nplsi
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
          do j=1,nplsi
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
      edissml(:,:,:,0) = 0.0_DP
      do ix = 1, ndxa
        do iy = 1, ndya
          eneutrad(ix,iy,1,0) = sum(eneutrad(ix,iy,1,1:nstrai))
          do imol=1,nmoli
            edissml(ix,iy,imol,0) = sum(edissml(ix,iy,imol,1:nstrai)) 
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
              if(datm3.lt.500.) then
                DATM3=EXP(DATM3)
              else
                write(iunout,*) 
     .            '[DPC] Problem in wneusave: ln(datm3) = ', 
     .             datm3, ' --- exponential will overflow'
                write(iunout,*) '[DPC] TE, DE = ', TE, DE
                datm3=1.0d30
              endif
              if(dpls3.lt.500.) then
                DPLS3=EXP(DPLS3)
              else
                write(iunout,*) 
     .           '[DPC] Problem in wneusave: ln(dpls3) = ', 
     .            dpls3, ' --- exponential will overflow'
                write(iunout,*) '[DPC] TE, DE = ', TE, DE
                dpls3=1.0d30
              endif
              if(dmol3.lt.500.) then
                DMOL3=EXP(DMOL3)
              else
                write(iunout,*) 
     .           '[DPC] Problem in wneusave: ln(dmol3) = ', 
     .            dmol3, ' --- exponential will overflow'
                write(iunout,*) '[DPC] TE, DE = ', TE, DE
                dmol3=1.0d30
              endif
              if(dion3.lt.500.) then
                DION3=EXP(DION3)
              else
                write(iunout,*) 
     .           '[DPC] Problem in wneusave: ln(dion3) = ', 
     .            dion3, ' --- exponential will overflow'
                write(iunout,*) '[DPC] TE, DE = ', TE, DE
                dion3=1.0d30
              endif
              if(dnml3.lt.500.) then
                DNML3=EXP(DNML3)
              else
                write(iunout,*) 
     .           '[DPC] Problem in wneusave: ln(dnml3) = ', 
     .            dnml3, ' --- exponential will overflow'
                write(iunout,*) '[DPC] TE, DE = ', TE, DE
                dnml3=1.0d30
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

              RATIO2=0.
              do J=0,8
                DEJ=DEF**J
                do I=0,8
                  TEI=TEF**I
                  RATIO2=RATIO2+RH2PH2(I,J)*TEI*DEJ
                end do
              end do
              if(RATIO2.lt.500.) then
                RATIO2=EXP(RATIO2)
              else
                write(iunout,*) 
     .           '[DPC] Problem in wneusave: ln(ratio2) = ', 
     .            RATIO2,' --- exponential will overflow'
                RATIO2=1.0d30
              endif

              !C
              !C  CHANNEL 1
              !C  H ALPHA SOURCE RATE:  PHOTONS/SEC/M**3
              !c  linear in atomic density (ionisation)
              !C
              do  IATM=1,NATMI
                !C  HYDROGENIC SPECIES?
                IF (NCHARA(IATM).eq.1) then
                  !c IF (LPDENA) DA=DATM3*PDENA(IATM,NCELL)
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
                  DPP=DPLS3*DIIN(IPLS,NCELL)
                  SIGADD2=SIGADD2+DPP*FAC32
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
                  !c IF (LPDENM) DM=DMOL3*PDENM(IMOL,NCELL)
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
                  !c IF (LPDENM) DI=DION3*PDENM(IMOL,NCELL)*RATIO2
                  di=dion3*dmb2(ix,iy,imol,1)*ratio2
                  SIGADD4=SIGADD4+DI*FAC32
                end if
              end do
              !C
              !C  CHANNEL 5
              !C  H ALPHA SOURCE RATE:  PHOTONS/SEC/M**3
              !C  LINEAR IN H- DENSITY (CHARGE EXCHANGE RECOMBINATION)
              !C
              do IMOL=1,NMOLI
                IF (NCHARM(IMOL).eq.2) then
                  !c IF (LPDENM) DN=DNML3*PDENM(IMOL,NCELL)*RATIO7
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
        WRITE (IUNOUT,*) ' RADIATED POWER BY HALPHA:',POWALF
        WRITE (IUNOUT,*) ' COUPL. TO GROUNDSTATE   :',POWALF1
        WRITE (IUNOUT,*) ' COUPLING TO CONTINUUM   :',POWALF2
        WRITE (IUNOUT,*) ' COUPLING TO MOLECULES   :',POWALF3
        WRITE (IUNOUT,*) ' COUPLING TO MOL.IONS    :',POWALF4
        WRITE (IUNOUT,*) ' COUPLING TO NEG.IONS    :',POWALF5
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
!c      write (iunout,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!c      write (iunout,'(/6x,20(a8,2x))')
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
!c      write (iunout,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!cc%%%
      write(iunout,*) 'ncutl,ncutb ',ncutl,ncutb
      write(iunout,'(a,7i6)') 'ndx,ndy,natm,ndxa,ndya,nfla,n1st',
     .                         ndx,ndy,natm,ndxa,ndya,nfla,n1st
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
!c       write (iunout,'()')
!c       write (iunout,*) 'wldpp:',k
!c       do j=nlim+1,nlim+nstsi !{
!c         write (iunout,'(1p,i4,20e9.2)') j,(wldpp(j,i,k),i=1,nplsi)
!c       end do !}
!c       write (iunout,'()')
!c       write (iunout,*) 'wldpeb, wldpa, wldpm:',k
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
        real(DP) :: dummy(nlimps,*)
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
      endif
      return
      end subroutine


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
