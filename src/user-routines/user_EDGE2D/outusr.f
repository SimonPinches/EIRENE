C+---------------------------------------------------------------+
C| Purpose:                                                      |
C| --------                                                      |
C| Write the eirene.tranfer file to give the EIRENE results back |
C| to EDGE2D.                                                    |
C| Write also eirene.chemFluxDep file, to store the neutral      |
C| particle fluxes on to the walls for the next EIRENE call      |
C+---------------------------------------------------------------+
C| Modifications:                                                |
C| --------------                                                |
C| 16/07/2010   D.Harting    Added writing of neutral fluxes to  |
C|                           file eirene.chemFluxDep. Added also |
C|                           two variables to eirene_user        |
C|                           namelist.                           |
C| 23/11/2010   D.Harting    If EDGE2D is used with density      |
C|                           control by puff+recycling, the      |
C|                           number of puffeing surfaces and thus|
C|                           the number of additional surfaces   |
C|                           (NLIM) may vary. Before, the neutral|
C|                           flux file eirene.chemFluxDep from   |
C|                           previous run was checked to have the|
C|                           same number of add. surfaces as the |
C|                           actual run. This forced a stopping  |
C|                           of the code. Now the actual triangle|
C|                           number and its side is checked, to  |
C|                           asure that the right neutral flux   |
C|                           is used.                            |
C| 24/11/2010   D.Harting    Added for backward compatibility    |
C|                           a version number to the neutral flux|
C|                           file eirene.chemFluxDep. If the     |
C|                           Version number in the file and in   |
C|                           code are not matching, the file is  |
C|                           not read and zero neutral flux is   |
C|                           assumed in the actual run. At the   |
C|                           end of the eirene run, a new neutral|
C|                           flux file with the current version  |
C|                           number is generated.                |
C+---------------------------------------------------------------+
      subroutine EIRENE_outusr
      use eirmod_precision
      use eirmod_parmmod
      use eirmod_comusr
      use eirmod_ccona
      use eirmod_cestim
      use eirmod_ctrcei
      use eirmod_comxs
      use eirmod_cgeom
      use eirmod_cinit
      use eirmod_ctrig
      use eirmod_csdvi
      use eirmod_clogau
      use eirmod_comsou
      use eirmod_comprt
      use eirmod_czt1
      use eirmod_coutau
      implicit none
      integer :: fp,mtri,ir,is

      real(dp) :: pa,pi,pm,pph,vvol,val,sumvol,sumval,vdenpara,bb
      integer :: j,i
      real(dp) :: c1,c2,c3,c4

csw ratecoeff.dat
      integer :: iaei, iacx, imei, imcx, iiei, iicx, iirc
      integer :: irei, ircx, irrc
      real(dp) :: de,di
      integer :: idsc,ipl,kk,nrc,ireac,k,np,nr,msg
csw 25oct07
      real(dp) :: x1,x2,y1,y2,ar,xc
      real(dp),allocatable :: sumpotpl(:)
csw
      logical :: lcxsigma

      integer, external :: EIRENE_idez

      logical, save :: ldebug

      integer, save :: eirene_nbirth,eirene_njetto
      character(len=256), save :: eirene_fbirth,eirene_ftransfer,
     &     eirene_fstoreneutflux
      real(dp) :: eirene_phi_offsets(9)
      integer :: eirene_wallFluxModel ! calculation of wall fluxes for chemical sputtering
c                = 0: no wall fluxes are used (old edge2d model)
c                = 1: only ion fluxes are used
c                = 2: ion fluxes and neutral fluxes from last eirene iteration are used
c                = 3: ion and neutral fluxes are used, and EIRENE is iterated to give
c                     converged neutral fluxes.
      logical :: eirene_use_elstepdat_bug
      logical  :: lfound
      real(dp) :: neutralFluxFileVersion
c     replicate old sputtered flux arrays sptpl, sptat, sptml, sptio, sptpht
c     for the moment these are filled with values from sptpltot, sptatot,sptmtot,sptitot,sptphtot
c     in the future it is better to pass particle resolved sputtered fluxes
c     from sptXY X=PH,I,A,M,P Y=PHT,IO,AT,ML,PL
      real(dp),dimension(npls,nlimps)  :: sptpl
      real(dp),dimension(natm,nlimps)  :: sptat
      real(dp),dimension(nmol,nlimps)  :: sptml
      real(dp),dimension(nion,nlimps)  :: sptio
      real(dp),dimension(nphot,nlimps) :: sptpht
      namelist /eirene_user/eirene_nbirth,eirene_njetto,
     .                      eirene_fbirth,eirene_ftransfer,
     .                      eirene_phi_offsets,
     .                      eirene_fstoreneutflux,
     .                      eirene_wallFluxModel,
     .                      eirene_use_elstepdat_bug

      NeutralFluxFileVersion = 1.0
      ldebug=.false.
csw
      eirene_ftransfer = 'eirene.transfer'
      eirene_njetto=0
cdmh
      eirene_fstoreneutflux = 'eirene.chemFluxDep'
      eirene_wallFluxModel = 1
cdmh
      open(unit=9998,file='eirene_user.namelist')
      read(9998,eirene_user)
      close(9998)
csw

C     fill replicated sputtered flux arrays
      sptpl(:,:) = 0.d0
      sptpl(1,1:nlimps) = sptpltot(1:nlimps)
      sptat(:,:) = 0.d0
      sptat(1,1:nlimps) = sptatot(1:nlimps)
      sptml(:,:) = 0.d0
      sptml(1,1:nlimps) = sptmtot(1:nlimps)
      sptio(:,:) = 0.d0
      sptio(1,1:nlimps) = sptitot(1:nlimps)
      sptpht(:,:) = 0.d0
      sptpht(1,1:nlimps) = sptphtot(1:nlimps)

csw 25oct07
      allocate(sumpotpl(npls))
csw
      do k=1,npls
         sumpotpl(k) = 0.
         do np=1,3
            do nr=1,ntrii
                if(  inmti(np,nr) == 2+1 .or.
     .               inmti(np,nr) == 3+1) then
                   msg = nlim+nsts+inspat(np,nr)

                   x1 = xtrian(necke(np,nr))
                   y1 = ytrian(necke(np,nr))
                   if(np == 3) then
                      x2 = xtrian(necke(1,nr))
                      y2 = ytrian(necke(1,nr))
                   else
                      x2 = xtrian(necke(np+1,nr))
                      y2 = ytrian(necke(np+1,nr))
                   endif

                   ar = sqrt( (x1-x2)**2 + (y1-y2)**2)
                   xc = (x1+x2)/2.d0
c     potpl in 1/s/cm
                   sumpotpl(k) = sumpotpl(k) +
     .                  estims(naddw(25)+k,msg)
     .                  /1.6022e-19/2.d0/pia/xc
                endif
            enddo
         enddo
      enddo
csw


      fp = 4999
      open(unit=fp,file=trim(eirene_ftransfer),access='sequential',
     .     status='replace')

      write(fp,'(a)') '* ntrii,nrad  :'
      write(fp,'(3(1x,i6))') ntrii,nrad

      write(fp,'(a)') '* nstordr:'
      write(fp,'(1x,i6)') nstordr

      write(fp,'(a)') '* natm,nmol,nion,nphot  :'
      write(fp,'(4(1x,i6))') natm,nmol,nion,nphot

      write(fp,'(a)') '* npls:'
      write(fp,'(3(1x,i6))') npls

      write(fp,'(a)') '* nlimps,nlim,nsts:'
      write(fp,'(3(1x,i6))') nlimps,nlim,nsts

c---------------------------------------
      write(fp,'(a,i6)') '* BULK SPECIES NPLS = ',npls
      do ipls=1,npls
         ityp=4
         ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
         write(fp,'(a,i6)') ' VOL.AV. IPLS = ',ipls
         do ir=1,ntrii
            c1 = mapl(ipls,ir)
            c2 = mmpl(ipls,ir)
            c3 = mipl(ipls,ir)
c            c4 = mphpl(ipls,ir)
            c4 = 0.
            write(fp,'(i6,30(1x,e14.6))') ir,
     .           papl(ipls,ir),
     .           pmpl(ipls,ir),
     .           pipl(ipls,ir),
c     .           pphpl(ipls,ir),
     .           0.,
c pspl:
     .           papl(ipls,ir)+pmpl(ipls,ir)+pipl(ipls,ir)
c     .                        +pphpl(ipls,ir),
     .                        +0.,
c mapl:
     .           c1,
c mmpl:
     .           c2,
c mipl:
     .           c3,
c mphpl:
     .           c4,
c mspl:
     .           c1+c2+c3+c4,

c not used anymore... (changed background profiles in extra file?)
     .           diin(ipls,ir),
     .           vxin(ipls,ir),
     .           vyin(ipls,ir),
     .           vzin(ipls,ir),
     .           bvin(ipls,ir),
     .           tiin(ipls,ir),
     .           edrift(ipls,ir),

c     .           eapl(ipls,ir),empl(ipls,ir),eipl(ipls,ir),ephpl(ir),
     .           eapl(ipls,ir),empl(ipls,ir),eipl(ipls,ir),0.,
c espl:
c     .           eapl(ipls,ir)+empl(ipls,ir)+eipl(ipls,ir)+ephpl(ir),
     .           eapl(ipls,ir)+empl(ipls,ir)+eipl(ipls,ir)+0.,

c     .           eael(ir),emel(ir),eiel(ir),ephel(ir),
     .           eael(ir),emel(ir),eiel(ir),0.,
c esel:
c     .           eael(ir)+emel(ir)+eiel(ir)+ephel(ir)
     .           eael(ir)+emel(ir)+eiel(ir)+0.
         enddo
         write(fp,'(a,i6)') ' SRF.AV. IPLS = ',ipls
         do is=1,nlimps
            write(fp,'(i6,20(1x,e14.6))') is,
     .           potpl(ipls,is),
     .           eotpl(ipls,is),
     .           sptpl(ipls,is),
     .           spump(ispz,is)
         enddo
         write(fp,'(20(1x,e14.6))') sumpotpl(ipls),srec(ipls,0)
      enddo

c---------------------------------------
      write(fp,'(a,i6)') '* ATOMIC SPECIES, NATM = ',natm
      do iatm=1,natm
         ityp=1
         ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
         write(fp,'(a,i6)') ' VOL.AV. IATM = ',iatm
         do ir=1,ntrii
            bb = sqrt( bxin(ir)**2 + byin(ir)**2 + bzin(ir)**2 )
            if( bb > 0.) then
               vdenpara = vxdena(iatm,ir)*bxin(ir) + vydena(iatm,ir)
     .              *byin(ir) + vzdena(iatm,ir)*bzin(ir)
               vdenpara = vdenpara / bb
            else
               vdenpara=0.
            endif
            write(fp,'(i6,20(1x,e14.6))') ir,
     .           pdena(iatm,ir),
     .           vxdena(iatm,ir),
     .           vydena(iatm,ir),
     .           vzdena(iatm,ir),
     .           vdenpara,
     .           edena(iatm,ir),
     .           edena(iatm,ir)/(pdena(iatm,ir)+1.d-60),
     .           0.,
     .           sigma(1,ir),
     .           0.,
     .           0.
         enddo
         write(fp,'(a,i6)') ' SRF.AV. IATM = ',iatm
         do is=1,nlimps
            write(fp,'(i6,20(1x,e14.6))') is,
     .           potat(iatm,is),
     .           prfaat(iatm,is),
     .           prfmat(iatm,is),
     .           prfiat(iatm,is),
c     .           prfphat(iatm,is),
     .           0.,
     .           prfpat(iatm,is),
     .           eotat(iatm,is),
     .           erfaat(iatm,is),
     .           erfmat(iatm,is),
     .           erfiat(iatm,is),
c     .           erfphat(iatm,is),
     .           0.,
     .           erfpat(iatm,is),
     .           sptat(iatm,is),
     .           spump(ispz,is)
         enddo
      enddo

c---------------------------------------
      write(fp,'(a,i6)') '* MOLECULAR SPECIES, NMOL = ',nmol
      do imol=1,nmol
         ityp=2
         ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
         write(fp,'(a,i6)') ' VOL.AV. IMOL = ',imol
         do ir=1,ntrii
            bb = sqrt( bxin(ir)**2 + byin(ir)**2 + bzin(ir)**2 )
            if( bb > 0.) then
               vdenpara = vxdenm(imol,ir)*bxin(ir) + vydenm(imol,ir)
     .              *byin(ir) + vzdenm(imol,ir)*bzin(ir)
               vdenpara = vdenpara / bb
            else
               vdenpara=0.
            endif
            write(fp,'(i6,20(1x,e14.6))') ir,
     .           pdenm(imol,ir),
     .           vxdenm(imol,ir),
     .           vydenm(imol,ir),
     .           vzdenm(imol,ir),
     .           vdenpara,
     .           edenm(imol,ir),
     .           edenm(imol,ir)/(pdenm(imol,ir)+1.d-60),
     .           0.,
     .           0.,
     .           0.
         enddo
         write(fp,'(a,i6)') ' SRF.AV. IMOL = ',imol
         do is=1,nlimps
            write(fp,'(i6,20(1x,e14.6))') is,
     .           potml(imol,is),
     .           prfaml(imol,is),
     .           prfmml(imol,is),
     .           prfiml(imol,is),
c     .           prfphml(imol,is),
     .           0.,
     .           prfpml(imol,is),
     .           eotml(imol,is),
     .           erfaml(imol,is),
     .           erfmml(imol,is),
     .           erfiml(imol,is),
c     .           erfphml(imol,is),
     .           0.,
     .           erfpml(imol,is),
     .           sptml(imol,is),
     .           spump(ispz,is)
         enddo
      enddo

c---------------------------------------
      write(fp,'(a,i6)') '* IONIC SPECIES, NION = ',nion
      do iion=1,nion
         ityp=3
         ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
         write(fp,'(a,i6)') ' VOL.AV. IION = ',iion
         do ir=1,ntrii
            bb = sqrt( bxin(ir)**2 + byin(ir)**2 + bzin(ir)**2 )
            if( bb > 0.) then
               vdenpara = vxdeni(iion,ir)*bxin(ir) + vydeni(iion,ir)
     .              *byin(ir) + vzdeni(iion,ir)*bzin(ir)
               vdenpara = vdenpara / bb
            else
               vdenpara=0.
            endif
            write(fp,'(i6,20(1x,e14.6))') ir,
     .           pdeni(iion,ir),
     .           vxdeni(iion,ir),
     .           vydeni(iion,ir),
     .           vzdeni(iion,ir),
     .           vdenpara,
     .           edeni(iion,ir),
     .           edeni(iion,ir)/(pdeni(iion,ir)+1.d-60),
     .           0.,
     .           0.,
     .           0.
         enddo
         write(fp,'(a,i6)') ' SRF.AV. IION = ',iion
         do is=1,nlimps
            write(fp,'(i6,20(1x,e14.6))') is,
     .           potio(iion,is),
     .           prfaio(iion,is),
     .           prfmio(iion,is),
     .           prfiio(iion,is),
c     .           prfphio(iion,is),
     .           0.,
     .           prfpio(iion,is),
     .           eotio(iion,is),
     .           erfaio(iion,is),
     .           erfmio(iion,is),
     .           erfiio(iion,is),
c     .           erfphio(iion,is),
     .           0.,
     .           erfpio(iion,is),
     .           sptio(iion,is),
     .           spump(ispz,is)
         enddo
      enddo

c---------------------------------------
      write(fp,'(a,i6)') '* PHOTONIC SPECIES, NPHOT = ',nphot
      do iphot=1,nphot
         ityp=0
         ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
         write(fp,'(a,i6)') ' VOL.AV. IPHOT = ',iphot
         do ir=1,ntrii
            bb = sqrt( bxin(ir)**2 + byin(ir)**2 + bzin(ir)**2 )
            if( bb > 0.) then
               vdenpara = vxdenph(iphot,ir)*bxin(ir) + vydenph(iphot,ir)
     .              *byin(ir) + vzdenph(iphot,ir)*bzin(ir)
               vdenpara = vdenpara / bb
            else
               vdenpara=0.
            endif
            write(fp,'(i6,20(1x,e14.6))') ir,
     .           pdenph(iphot,ir),
     .           vxdenph(iphot,ir),
     .           vydenph(iphot,ir),
     .           vzdenph(iphot,ir),
     .           vdenpara,
     .           edenph(iphot,ir),
     .           edenph(iphot,ir)/(pdenph(iphot,ir)+1.d-60),
     .           0.,
     .           0.,
     .           0.
         enddo
         write(fp,'(a,i6)') ' SRF.AV. IPHOT = ',iphot
         do is=1,nlimps
            write(fp,'(i6,20(1x,e14.6))') is,
     .           potpht(iphot,is),
     .           prfapht(iphot,is),
     .           prfmpht(iphot,is),
     .           prfipht(iphot,is),
     .           prfphpht(iphot,is),
     .           prfppht(iphot,is),
     .           eotpht(iphot,is),
     .           erfapht(iphot,is),
     .           erfmpht(iphot,is),
     .           erfipht(iphot,is),
     .           erfphpht(iphot,is),
     .           erfppht(iphot,is),
     .           sptpht(iphot,is),
     .           spump(ispz,is)
         enddo
      enddo

c---------------------------------------
      write(fp,'(a,i6)') '* MISC DATA'
      do ir=1,ntrii
         write(fp,'(i6,20(1x,e14.6))') ir,
     .        dble(ncltal(ir)),0.,0.,
     .        vol(ir),voltal(ir),
     .        dein(ir),tein(ir),bxin(ir),byin(ir),bzin(ir),bfin(ir)
      enddo

csw 20dec07-----------------------------
      write(fp,'(a)') '* STRATUM DATA'
      write(fp,'(i6)') nstra
      do istra=1,nstra
        do iatm=1,natm
          write(fp,'(2i6,20(1x,e14.6))') istra,iatm,
     .              wtota(iatm,istra),0.
        enddo
        do imol=1,nmol
          write(fp,'(2i6,20(1x,e14.6))') istra,imol,
     .              wtotm(imol,istra),0.
        enddo
        do iion=1,nion
          write(fp,'(2i6,20(1x,e14.6))') istra,iion,
     .              wtoti(iion,istra),0.
        enddo
        do iphot=1,nphot
          write(fp,'(2i6,20(1x,e14.6))') istra,iphot,
     .              wtotph(iphot,istra),0.
        enddo
        do ipls=1,npls
          write(fp,'(2i6,20(1x,e14.6))') istra,ipls,
     .              wtotp(ipls,istra),0.
        enddo
      enddo

csw 24jan08
c---------------------------------------
      if(eirene_njetto .gt. 0) then
        write(fp,'(a)') '* ADDITIONAL TALLY DATA FOR JETTO NEUTRALS'
        write(fp,'(6i6)') nadv
        do ir=1,ntrii
          write(fp,'(i6,20(1x,e14.6))') ir,( addv(i,ir), i=1,nadv ),
     .                                     ( ppat(i,ir), i=1,natm )
        enddo
      endif

      if(eirene_nbirth .gt. 0) then
        write(fp,'(a)') '* ADDITIONAL TALLY DATA FOR JETTO NBI'
        write(fp,'(6i6)') nstra
        do is=1,nstra
          write(fp,'(i6,20(1x,e14.6))') is, etota(is)
        enddo
      endif

      close(fp)
csw 25oct07
      deallocate(sumpotpl)
csw

cdmh 15jun10
c----------------------------------------
c     store neutral particle fluxes [A] on wall
      fp = 4999
      open(unit=fp,file=trim(eirene_fstoreneutflux),access='sequential',
     .     status='replace')
      write(fp,'(a28,f14.6)') "* Neutral flux file version:",
     &     NeutralFluxFileVersion
      write(fp,'(a,a)') '*  NLIM,   NSTS,  NGITT, NGSTAL,',
     &     '   NATM, NLMPGS,  NTRII'
      write(fp,'(7i8)') NLIM, NSTS, NGITT, NGSTAL, NATM, NLMPGS, NTRII
      do iatm=1,NATM
         write(fp,'(a,i0)') '* neutral fluxes from atom species ',iatm
         write(fp,'(a)') "*  idx,  neutral_flux, ITRIA, ISIDE, ISURF"
         do is=1,NLMPGS
c           find corresponding triangle
            lfound = .false.
            nr = 0
            np = 0
            do i=1,ntrii
               do j=1,3
                  if ((INSPAT(j,i).eq. is -(NLIM+NSTS))
     &                 .and.(INSPAT(j,i).ne.0) )then
                     if (lfound) then
                        write(iunout,*)"* EIRENE_OUTUSR:"
                        write(iunout,*)"* Edge twice found"
                        call EIRENE_exit_own(1)
                     endif
                     lfound=.true.
                     nr = i
                     np = j
                  endif
               enddo
            enddo
            if ((nr.le.0).or.(np.le.0)) then
               write(fp,'(i6,1x,e14.6,1x,i6,1x,i6,1x,i6)')
     &              is,POTAT(iatm,is),nr,np,0
            else
               write(fp,'(i6,1x,e14.6,1x,i6,1x,i6,1x,i6)')
     &              is,POTAT(iatm,is),nr,np,INMTI(np,nr)
            endif
         enddo                  !is
      enddo                     !iatm
      close(fp)

      end
