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
      USE EIRMOD_CTRCEI
      USE EIRMOD_COMSOU
      USE EIRMOD_CPES
      USE EIRMOD_MPI
       
      implicit none

      integer, save :: istra_old=-1,
     .                 ityp_old=-1,
     .                 iphot_old=-1, 
     .                 iatm_old=-1, 
     .                 imol_old=-1, 
     .                 iion_old=-1,
     .                 ipls_old=-1
      real(dp), allocatable, save :: time_array(:,:,:), helpt(:,:,:)
      real(dp) :: tim_spent
      real(dp), save :: tim_start=0._dp, tim_end=0._dp 
      real(dp) :: eirene_second_own
      integer :: istr, nn, ier, itp, isp
      
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

      entry eirene_output_partinfo

      if (nprs > 0) then
        nn = max(nphot,natm,nmol,nion)
        allocate (helpt(0:3,0:nn,0:nstra))
        helpt = 0._dp
        call mpi_reduce(time_array,helpt,4*(nn+1)*(nstra+1),
     .       mpi_double_precision,mpi_sum,0,mpi_comm_world,ier)
        time_array = helpt
        deallocate (helpt)
        if (my_pe /= 0) return
      end if

      call eirene_leer(2)

      call eirene_headng ('STATISTICS OVER CPU TIME SPENT'//
     .                    ' IN FOLLOWING TRAJECTORIES ',57)

! sum over species
      do istr = 1, nstrai
        time_array(0,0,istr) = sum(time_array(0,1:,istr))
        time_array(1,0,istr) = sum(time_array(1,1:,istr))
        time_array(2,0,istr) = sum(time_array(2,1:,istr))
        time_array(3,0,istr) = sum(time_array(3,1:,istr))
      end do
! sum over strata
      do itp = 0,3
        do isp = 0,nn
          time_array(itp,isp,0) = sum(time_array(itp,isp,1:nstrai))
        end do
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

        if (any(time_array(0,:,istr) > 0._dp)) then
          call eirene_leer(1)
          write (iunout,*) 'TIME (SEC) SPENT IN FOLLOWING '
          CALL EIRENE_MASYR1 ('PHOTONS = ',time_array(0,0:nphot,:),
     .                LOGPHOT,ISTR,0,NPHOT,0,NSTRA,TEXTS(1))
          CALL EIRENE_MASAGE
     .      ('SUM OVER SPECIES                               ')
          CALL EIRENE_MASR1 ('TOTAL=  ',SUM(time_array(0,1:nphot,ISTR)))
        end if

        if (any(time_array(1,:,istr) > 0._dp)) then
          call eirene_leer(1)
          write (iunout,*) 'TIME (SEC) SPENT IN FOLLOWING '
          CALL EIRENE_MASYR1 ('ATOMS =   ',time_array(1,0:natm,:),
     .                LOGATM,ISTR,0,NATM,0,NSTRA,TEXTS(NSPH+1))
         CALL EIRENE_MASAGE
     .      ('SUM OVER SPECIES                               ')
          CALL EIRENE_MASR1 ('TOTAL=  ',SUM(time_array(1,1:natm,ISTR)))
        end if

        if (any(time_array(2,:,istr) > 0._dp)) then
          call eirene_leer(1)
          write (iunout,*) 'TIME (SEC) SPENT IN FOLLOWING '
          CALL EIRENE_MASYR1 ('MOLECULES=',time_array(2,0:nmol,:),
     .                LOGMOL,ISTR,0,NMOL,0,NSTRA,TEXTS(NSPA+1))
          CALL EIRENE_MASAGE
     .      ('SUM OVER SPECIES                               ')
          CALL EIRENE_MASR1 ('TOTAL=  ',SUM(time_array(2,1:nmol,ISTR)))
        end if

        if (any(time_array(3,:,istr) > 0._dp)) then
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

