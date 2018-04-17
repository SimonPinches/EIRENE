cdr  Jan 18:  bypass this actions for photons (ityp=0). Code not ready for photon transport.
cpb:  added: cpu time statistics by particle type, species and stratum: time_array

      subroutine eirene_switch_partinfo
c  added oct. 2017:
c  this routine sets the various pointers for tallies,
c  for a unified treatment of scoring volume tallies (update), 
c                                     mfp evaluation (fpath)
c  and a for scoring a number of surface tallies (locate, escape, ...) 

c  It is called from: LOCATE (for primary source particles),
c          and after: COLLIDE and REFLEC (for new post collision species)
c
c  Input:  istra, ityp, iphot, iatm, imol, iion, ipls
c  Output: ixspz,nmetoff,logphot,logatm,logmol,logion

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
       
      implicit none

      integer, save :: istra_old=-1,
     .                 ityp_old=-1,
     .                 iphot_old=-1, 
     .                 iatm_old=-1, 
     .                 imol_old=-1, 
     .                 iion_old=-1,
     .                 ipls_old=-1
      real(dp), allocatable, save :: time_array(:,:,:)
      real(dp) :: tim_spent
      real(dp), save :: tim_start=0._dp, tim_end=0._dp 
      real(dp) :: eirene_second_own
      integer :: istr, it, is
      
      if ((ityp_old == ityp) .and. (iatm_old == iatm) .and.
     .    (imol_old == imol) .and. (iion_old == iion) .and.
     .    (iphot_old == iphot) .and. (ipls_old == ipls).and.
     .    (istra_old == istra)) return
      
      if (trchktim) then
        tim_end = EIRENE_SECOND_OWN()
        tim_spent = tim_end - tim_start
        tim_start = tim_end

        select case(ityp_old)
        case(0)
!  photons
          time_array(ityp_old,iphot_old,istra_old) = 
     .    time_array(ityp_old,iphot_old,istra_old) + tim_spent
        case(1)
!  atoms
          time_array(ityp_old,iatm_old,istra_old) = 
     .    time_array(ityp_old,iatm_old,istra_old) + tim_spent
        case(2)
!  molecules
          time_array(ityp_old,imol_old,istra_old) = 
     .    time_array(ityp_old,imol_old,istra_old) + tim_spent
        case(3)
!  test ions
          time_array(ityp_old,iion_old,istra_old) = 
     .    time_array(ityp_old,iion_old,istra_old) + tim_spent
        case(4)
!  bulk ions: NOT IN USE
!         time_array(ityp_old,ipls_old,istra_old) = 
!     .   time_array(ityp_old,ipls_old,istra_old) + tim_spent
        case default
          if (.not.allocated(time_array)) then
            allocate (time_array(0:3,
     .                     0:max(nphot,natm,nmol,nion),0:nstra))
            time_array = 0._dp
          end if
        end select
      end if

C  save stratum, old type, species
      istra_old= istra
      ityp_old = ityp

      iphot_old= iphot
      iatm_old = iatm
      imol_old = imol
      iion_old = iion
      ipls_old = ipls

      select case(ityp)

      case(0)
!  photons

       LPDENX  => LPDENPH 
       LEDENX  => LEDENPH 
       LPXEL   => LPPHEL  
       LPXAT   => LPPHAT 
       LPXML   => LPPHML 
       LPXIO   => LPPHIO
       LPXPL   => LPPHPL
       LEXEL   => LEPHEL   
       LEXAT   => LEPHAT   
       LEXML   => LEPHML
       LEXIO   => LEPHIO
       LEXPL   => LEPHPL
       LVXDENX => LVXDENPH
       LVYDENX => LVYDENPH
       LVZDENX => LVZDENPH
       LMXPL   => LMPHPL
       LPXX    => LPPHPHT 
       LEXX    => LEPHPHT 

       IF (LPDENX)  PDENX  => PDENPH(IPHOT,:) 
       IF (LEDENX)  EDENX  => EDENPH(IPHOT,:) 
       IF (LPXEL)   PXEL   => PPHEL(:)   
       IF (LPXAT)   PXAT   => PPHAT(1:NATMI,:)  
       IF (LPXML)   PXML   => PPHML(1:NMOLI,:)  
       IF (LPXIO)   PXIO   => PPHIO(1:NIONI,:) 
       IF (LPXPL)   PXPL   => PPHPL(1:NPLSI,:)
       IF (LEXEL)   EXEL   => EPHEL(:)   
       IF (LEXAT)   EXAT   => EPHAT(:)   
       IF (LEXML)   EXML   => EPHML(:)
       IF (LEXIO)   EXIO   => EPHIO(:)
       IF (LEXPL)   EXPL   => EPHPL(1:NPLSI,:)
       IF (LVXDENX) VXDENX => VXDENPH(IPHOT,:)
       IF (LVYDENX) VYDENX => VYDENPH(IPHOT,:)
       IF (LVZDENX) VZDENX => VZDENPH(IPHOT,:)
       IF (LMXPL)   MXPL   => MPHPL(1:NPLSI,:)
       IF (LPXX)    PXX    => PPHPHT(1:NPHOTI,:)
       IF (LEXX)    EXX    => EPHPHT(:)

       LEX     => LEPH

       LGXCX => LGACX
       LGXEI => LGAEI
       LGXEL => LGAEL
       LGXPI => LGAPI

       NXEII => NAEII(IATM)
       NXCXI => NACXI(IATM)
       NXELI => NAELI(IATM)
       NXPII => NAPII(IATM)

       NPBGKX => NPBGKA(IATM)

       IXSPZ = IPHOT
       NMETOFF = 0

       RMASSX => RMASSA(IATM)
       CNDYNX => CNDYNA(IATM)
 
       LOGPHOT(IPHOT,ISTRA)=.TRUE.
       LOGXSPZ => LOGPHOT(IPHOT,ISTRA)

      case(1)
!  atoms
       PDENX  => PDENA(IATM,:) 
       EDENX  => EDENA(IATM,:) 
       PXEL   => PAEL(:)   
       PXAT   => PAAT(1:NATMI,:)  
       PXML   => PAML(1:NMOLI,:)  
       PXIO   => PAIO(1:NIONI,:) 
       PXPL   => PAPL(1:NPLSI,:)
       EXEL   => EAEL(:)   
       EXAT   => EAAT(:)   
       EXML   => EAML(:)
       EXIO   => EAIO(:)
       EXPL   => EAPL(1:NPLSI,:)
       VXDENX => VXDENA(IATM,:)
       VYDENX => VYDENA(IATM,:)
       VZDENX => VZDENA(IATM,:)
       MXPL   => MAPL(1:NPLSI,:)
       PXX    => PAAT(1:NATMI,:)
       EXX    => EAAT(:)

       LPDENX  => LPDENA 
       LEDENX  => LEDENA 
       LPXEL   => LPAEL  
       LPXAT   => LPAAT 
       LPXML   => LPAML 
       LPXIO   => LPAIO
       LPXPL   => LPAPL
       LEXEL   => LEAEL   
       LEXAT   => LEAAT   
       LEXML   => LEAML
       LEXIO   => LEAIO
       LEXPL   => LEAPL
       LVXDENX => LVXDENA
       LVYDENX => LVYDENA
       LVZDENX => LVZDENA
       LMXPL   => LMAPL
       LPXX    => LPAAT 
       LEXX    => LEAAT 

       LEX     => LEA

       LGXCX => LGACX
       LGXEI => LGAEI
       LGXEL => LGAEL
       LGXPI => LGAPI

       NXEII => NAEII(IATM)
       NXCXI => NACXI(IATM)
       NXELI => NAELI(IATM)
       NXPII => NAPII(IATM)

       NPBGKX => NPBGKA(IATM)

       IXSPZ = IATM
       NMETOFF = NSPH

       RMASSX => RMASSA(IATM)
       CNDYNX => CNDYNA(IATM)

       LOGATM(IATM,ISTRA)=.TRUE.
       LOGXSPZ => LOGATM(IATM,ISTRA)

      case(2)
!  molecules

       PDENX  => PDENM(IMOL,:) 
       EDENX  => EDENM(IMOL,:) 
       PXEL   => PMEL(:)   
       PXAT   => PMAT(1:NATMI,:)  
       PXML   => PMML(1:NMOLI,:)  
       PXIO   => PMIO(1:NIONI,:) 
       PXPL   => PMPL(1:NPLSI,:)
       EXEL   => EMEL(:)   
       EXAT   => EMAT(:)   
       EXML   => EMML(:)
       EXIO   => EMIO(:)
       EXPL   => EMPL(1:NPLSI,:)
       VXDENX => VXDENM(IMOL,:)
       VYDENX => VYDENM(IMOL,:)
       VZDENX => VZDENM(IMOL,:)
       MXPL   => MMPL(1:NPLSI,:)
       PXX    => PMML(1:NMOLI,:)
       EXX    => EMML(:)

       LPDENX  => LPDENM 
       LEDENX  => LEDENM 
       LPXEL   => LPMEL  
       LPXAT   => LPMAT 
       LPXML   => LPMML 
       LPXIO   => LPMIO
       LPXPL   => LPMPL
       LEXEL   => LEMEL   
       LEXAT   => LEMAT   
       LEXML   => LEMML
       LEXIO   => LEMIO
       LEXPL   => LEMPL
       LVXDENX => LVXDENM
       LVYDENX => LVYDENM
       LVZDENX => LVZDENM
       LMXPL   => LMMPL
       LPXX    => LPMML 
       LEXX    => LEMML 

       LEX     => LEM

       LGXCX => LGMCX
       LGXEI => LGMEI
       LGXEL => LGMEL
       LGXPI => LGMPI

       NXEII => NMEII(IMOL)
       NXCXI => NMCXI(IMOL)
       NXELI => NMELI(IMOL)
       NXPII => NMPII(IMOL)

       NPBGKX => NPBGKM(IMOL)

       IXSPZ = IMOL
       NMETOFF = NSPA

       RMASSX => RMASSM(IMOL)
       CNDYNX => CNDYNM(IMOL)

       LOGMOL(IMOL,ISTRA)=.TRUE.
       LOGXSPZ => LOGMOL(IMOL,ISTRA)

      case(3)
!  test ions

       PDENX  => PDENI(IION,:) 
       EDENX  => EDENI(IION,:) 
       PXEL   => PIEL(:)   
       PXAT   => PIAT(1:NATMI,:)  
       PXML   => PIML(1:NMOLI,:)  
       PXIO   => PIIO(1:NIONI,:) 
       PXPL   => PIPL(1:NPLSI,:)
       EXEL   => EIEL(:)   
       EXAT   => EIAT(:)   
       EXML   => EIML(:)
       EXIO   => EIIO(:)
       EXPL   => EIPL(1:NPLSI,:)
       VXDENX => VXDENI(IION,:)
       VYDENX => VYDENI(IION,:)
       VZDENX => VZDENI(IION,:)
       MXPL   => MIPL(1:NPLSI,:)
       PXX    => PIIO(1:NIONI,:)
       EXX    => EIIO(:)

       LPDENX  => LPDENI 
       LEDENX  => LEDENI 
       LPXEL   => LPIEL  
       LPXAT   => LPIAT 
       LPXML   => LPIML 
       LPXIO   => LPIIO
       LPXPL   => LPIPL
       LEXEL   => LEIEL   
       LEXAT   => LEIAT   
       LEXML   => LEIML
       LEXIO   => LEIIO
       LEXPL   => LEIPL
       LVXDENX => LVXDENI
       LVYDENX => LVYDENI
       LVZDENX => LVZDENI
       LMXPL   => LMIPL
       LPXX    => LPIIO 
       LEXX    => LEIIO 

       LEX     => LEIO

       LGXCX => LGICX
       LGXEI => LGIEI
       LGXEL => LGIEL
       LGXPI => LGIPI

       NXEII => NIEII(IION)
       NXCXI => NICXI(IION)
       NXELI => NIELI(IION)
       NXPII => NIPII(IION)

       NPBGKX => NPBGKI(IION)

       IXSPZ = IION
       NMETOFF = NSPAM

       RMASSX => RMASSI(IION)
       CNDYNX => CNDYNI(IION)

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

