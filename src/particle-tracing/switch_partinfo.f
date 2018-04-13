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
!  photons: not ready
          time_array(ityp_old,iphot_old,istra_old) = 
     .      time_array(ityp_old,iphot_old,istra_old) + tim_spent
        case(1)
!  atoms
          time_array(ityp_old,iatm_old,istra_old) = 
     .      time_array(ityp_old,iatm_old,istra_old) + tim_spent
        case(2)
!  molecules
          time_array(ityp_old,imol_old,istra_old) = 
     .      time_array(ityp_old,imol_old,istra_old) + tim_spent
        case(3)
!  test ions
          time_array(ityp_old,iion_old,istra_old) = 
     .      time_array(ityp_old,iion_old,istra_old) + tim_spent
        case(4)
!  bulk ions
!          time_array(ityp_old,ipls_old,istra_old) = 
!     .      time_array(ityp_old,ipls_old,istra_old) + tim_spent
        case default
          if (.not.allocated(time_array)) then
            allocate (time_array(0:3,
     .                     0:max(nphot,natm,nmol,nion),0:nstra))
            time_array = 0._dp
          end if
        end select
      end if

C  save stratum, old typ, species
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
      NULLIFY (PXPL)
      NULLIFY (EXEL)
      NULLIFY (EXAT)
      NULLIFY (EXML)
      NULLIFY (EXIO)
      NULLIFY (EXPL)
      NULLIFY (VXDENX)
      NULLIFY (VYDENX)
      NULLIFY (VZDENX)
      NULLIFY (MXPL)
      NULLIFY (PXX)
      NULLIFY (EXX)
      
      select case(ityp)

      case(0)
!  photons: not ready

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

       IF (LPDENX)  PDENX  => PDENA(IATM,:) 
       IF (LEDENX)  EDENX  => EDENA(IATM,:) 
       IF (LPXEL)   PXEL   => PAEL(:)   
       IF (LPXAT)   PXAT   => PAAT(1:NATMI,:)  
       IF (LPXML)   PXML   => PAML(1:NMOLI,:)  
       IF (LPXIO)   PXIO   => PAIO(1:NIONI,:) 
       IF (LPXPL)   PXPL   => PAPL(1:NPLSI,:)
       IF (LEXEL)   EXEL   => EAEL(:)   
       IF (LEXAT)   EXAT   => EAAT(:)   
       IF (LEXML)   EXML   => EAML(:)
       IF (LEXIO)   EXIO   => EAIO(:)
       IF (LEXPL)   EXPL   => EAPL(1:NPLSI,:)
       IF (LVXDENX) VXDENX => VXDENA(IATM,:)
       IF (LVYDENX) VYDENX => VYDENA(IATM,:)
       IF (LVZDENX) VZDENX => VZDENA(IATM,:)
       IF (LMXPL)   MXPL   => MAPL(1:NPLSI,:)
       IF (LPXX)    PXX    => PAAT(1:NATMI,:)
       IF (LEXX)    EXX    => EAAT(:)

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

       IF (LPDENX)  PDENX  => PDENM(IMOL,:) 
       IF (LEDENX)  EDENX  => EDENM(IMOL,:) 
       IF (LPXEL)   PXEL   => PMEL(:)   
       IF (LPXAT)   PXAT   => PMAT(1:NATMI,:)  
       IF (LPXML)   PXML   => PMML(1:NMOLI,:)  
       IF (LPXIO)   PXIO   => PMIO(1:NIONI,:) 
       IF (LPXPL)   PXPL   => PMPL(1:NPLSI,:)
       IF (LEXEL)   EXEL   => EMEL(:)   
       IF (LEXAT)   EXAT   => EMAT(:)   
       IF (LEXML)   EXML   => EMML(:)
       IF (LEXIO)   EXIO   => EMIO(:)
       IF (LEXPL)   EXPL   => EMPL(1:NPLSI,:)
       IF (LVXDENX) VXDENX => VXDENM(IMOL,:)
       IF (LVYDENX) VYDENX => VYDENM(IMOL,:)
       IF (LVZDENX) VZDENX => VZDENM(IMOL,:)
       IF (LMXPL)   MXPL   => MMPL(1:NPLSI,:)
       IF (LPXX)    PXX    => PMML(1:NMOLI,:)
       IF (LEXX)    EXX    => EMML(:)

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

       IF (LPDENX)  PDENX  => PDENI(IION,:) 
       IF (LEDENX)  EDENX  => EDENI(IION,:) 
       IF (LPXEL)   PXEL   => PIEL(:)   
       IF (LPXAT)   PXAT   => PIAT(1:NATMI,:)  
       IF (LPXML)   PXML   => PIML(1:NMOLI,:)  
       IF (LPXIO)   PXIO   => PIIO(1:NIONI,:) 
       IF (LPXPL)   PXPL   => PIPL(1:NPLSI,:)
       IF (LEXEL)   EXEL   => EIEL(:)   
       IF (LEXAT)   EXAT   => EIAT(:)   
       IF (LEXML)   EXML   => EIML(:)
       IF (LEXIO)   EXIO   => EIIO(:)
       IF (LEXPL)   EXPL   => EIPL(1:NPLSI,:)
       IF (LVXDENX) VXDENX => VXDENI(IION,:)
       IF (LVYDENX) VYDENX => VYDENI(IION,:)
       IF (LVZDENX) VZDENX => VZDENI(IION,:)
       IF (LMXPL)   MXPL   => MIPL(1:NPLSI,:)
       IF (LPXX)    PXX    => PIIO(1:NIONI,:)
       IF (LEXX)    EXX    => EIIO(:)

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

      entry eirene_output_partinfo

      call eirene_leer(2)

      call eirene_headng ('STATISTICS OVER CPU TIME SPENT'//
     .                    ' IN FOLLOWING TRAJECTORIES ',57)

! sum over species
!      time_array(:,0,1:nstra) = sum(time_array,2)
      do it =0, 3
        do is = 1, nstrai
          time_array(it,0,is) = sum(time_array(it,1:,is))
        end do
      end do
! sum over strata
!      time_array(:,:,0) = sum(time_array,3)
      do it = 0, 3
        do is = 0, ubound(time_array,2)
           time_array(it,is,0) = sum(time_array(it,is,1:))
        end do
      end do

      write (0,*) ' vor do'
      do istr = 0, nstrai

        call eirene_leer(1)
        if (istr == 0) then
          write (iunout,'(1x,A,I6)') 'SUM OVER STRATA '
          write (iunout,'(1x,A,I6)') '=============== '
        else
          write (iunout,'(1x,A,I6)') 'STRATUM ',istr
          write (iunout,'(1x,A,I6)') '============== '
        end if

        if (any(time_array(0,:,:) > 0._dp)) then
          call eirene_leer(1)
          write (iunout,*) 'TIME (SEC) SPENT IN FOLLOWING '
          CALL EIRENE_MASYR1 ('PHOTONS = ',time_array(0,0:nphot,:),
     .                LOGPHOT,ISTR,0,NPHOT,0,NSTRA,TEXTS(1))
          CALL EIRENE_MASAGE
     .      ('SUM OVER SPECIES                               ')
          CALL EIRENE_MASR1 ('TOTAL=  ',SUM(time_array(0,1:nphot,ISTR)))
        end if

        if (any(time_array(1,:,:) > 0._dp)) then
          call eirene_leer(1)
          write (iunout,*) 'TIME (SEC) SPENT IN FOLLOWING '
          CALL EIRENE_MASYR1 ('ATOMS =   ',time_array(1,0:natm,:),
     .                LOGATM,ISTR,0,NATM,0,NSTRA,TEXTS(NSPH+1))
         CALL EIRENE_MASAGE
     .      ('SUM OVER SPECIES                               ')
          CALL EIRENE_MASR1 ('TOTAL=  ',SUM(time_array(1,1:natm,ISTR)))
        end if

        if (any(time_array(2,:,:) > 0._dp)) then
          call eirene_leer(1)
          write (iunout,*) 'TIME (SEC) SPENT IN FOLLOWING '
          CALL EIRENE_MASYR1 ('MOLECULES=',time_array(2,0:nmol,:),
     .                LOGMOL,ISTR,0,NMOL,0,NSTRA,TEXTS(NSPA+1))
          CALL EIRENE_MASAGE
     .      ('SUM OVER SPECIES                               ')
          CALL EIRENE_MASR1 ('TOTAL=  ',SUM(time_array(2,1:nmol,ISTR)))
        end if

        if (any(time_array(3,:,:) > 0._dp)) then
          call eirene_leer(1)
          write (iunout,*) 'TIME (SEC) SPENT IN FOLLOWING '
          CALL EIRENE_MASYR1 ('TEST IONS=',time_array(3,0:nion,:),
     .                LOGION,ISTR,0,NION,0,NSTRA,TEXTS(NSPAM+1))
          CALL EIRENE_MASAGE
     .      ('SUM OVER SPECIES                               ')
          CALL EIRENE_MASR1 ('TOTAL=  ',SUM(time_array(3,1:nion,ISTR)))
        end if

!        if (any(time_array(4,:,:) > 0._dp)) then
!          write (iunout,*) ' TIME (SEC) SPENT IN FOLLOWING '
!          CALL EIRENE_MASYR1 ('BULK IONS=',time_array(4,1:npls,:),
!     .                LOGPLS,ISTR,1,NPLS,1,NSTRA,TEXTS(NSPAMI+1))
!          CALL EIRENE_MASAGE
!     .      ('SUM OVER SPECIES                               ')
!          CALL EIRENE_MASR1 ('TOTAL=  ',SUM(time_array(4,1:npls,ISTR)))
!        end if

      end do
      return

      entry eirene_reinit_partinfo

      if (allocated(time_array)) deallocate(time_array)
      istra_old=-1
      ityp_old=-1
      iphot_old=-1 
      iatm_old=-1 
      imol_old=-1 
      iion_old=-1
      ipls_old=-1

      return      
      end subroutine eirene_switch_partinfo

