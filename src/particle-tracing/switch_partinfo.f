cdr  Jan 18:  bypass this actions for photons (ityp=0). Code not ready for photon transport.

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
       
      implicit none

      integer, save :: istra_old=-1,
     .                 ityp_old=-1,
     .                 iphot_old=-1, 
     .                 iatm_old=-1, 
     .                 imol_old=-1, 
     .                 iion_old=-1,
     .                 ipls_old=-1
      
      if ((ityp_old == ityp) .and. (iatm_old == iatm) .and.
     .    (imol_old == imol) .and. (iion_old == iion) .and.
     .    (iphot_old == iphot) .and. (ipls_old == ipls).and.
     .    (istra_old == istra)) return

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
!  photons: 
cdr: not ready.
cdr  currently: by-pass this code-section for photons (ityp=0),
cdr  as long as update, collide, fpath for photons are still kept as separate routines.

       return   ! for the time being....


       PDENX  => PDENPH(IPHOT,:) 
       EDENX  => EDENPH(IPHOT,:) 
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

