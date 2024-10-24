!pb  This subroutine has been moved from the SOLPS-ITER interface into the main EIRENE code
!pb  as EION is needed for calculation of radiation tallies

      !c
      !c*** Obtain ionization potentials for consistency with B2.5
      !c
      subroutine eirene_init_eion
      use eirmod_precision
      use eirmod_parmmod
      use eirmod_comusr
      use eirmod_ccona
      use eirmod_mpi
!pb
      use eirmod_comprt, only : IUNOUT
      use eirmod_cinit , only : MASTER_PATH
!pb
      implicit none
      integer iss, jatm
      character*256 filename
      logical found
      logical, save :: eion_set

      !C*** ionization potentials
      integer, save :: npot=20
      real(dp), allocatable, save :: pot_data(:)

      data eion_set /.false./

      if (eion_set) return

      filename='ionization_potentials'
      inquire(file=filename,exist=found)
      if(.not.found) then
        filename='../'//trim(filename)
        inquire(file=filename,exist=found)
      endif
      if (found) then
        open(UNIT=99,FILE=trim(filename))
      else
        inquire (FILE=trim(MASTER_PATH)//
     .   '/data.local/ionization_potentials',exist=found)
        if (found) then
          open(UNIT=99,FILE=trim(MASTER_PATH)//
     .     '/data.local/ionization_potentials')
        else
          inquire (FILE=trim(MASTER_PATH)//
     .     '/modules/B2.5/Database/ionization_potentials',exist=found)
          if (found) open(UNIT=99,FILE=trim(MASTER_PATH)//
     .     '/modules/B2.5/Database/ionization_potentials')
        endif
      endif

      if (found) then

        do jatm=1,natmi
          rewind(99)
          do iss=1,nchara(jatm)-1
            read(99,*)
          enddo
          read(99,*) Eion(jatm)
          if(nchara(jatm).eq.1 .and. nmassa(jatm).ne.1)
     .     Eion(jatm)=Eion(jatm)*
     .      (nmassa(jatm)*(pmassa+pmasse))/(pmasse+nmassa(jatm)*pmassa)
        enddo
        close(99)

      else

        write (iunout,*) 'in INIT_EION, no ionization file '
        allocate(pot_data(npot))
        pot_data = (/13.598_DP,  ! H
     .               24.587_DP,  ! He
     .                5.392_DP,  ! Li
     .                9.322_DP,  ! Be
     .                8.298_DP,  ! B
     .               11.260_DP,  ! C
     .               14.534_DP,  ! N
     .               13.618_DP,  ! O
     .               17.422_DP,  ! F
     .               21.564_DP,  ! Ne
     .                5.139_DP,  ! Na
     .                7.646_DP,  ! Mg
     .                5.986_DP,  ! Al
     .                8.151_DP,  ! Si
     .               10.486_DP,  ! P
     .               10.360_DP,  ! S
     .               12.967_DP,  ! Cl
     .               15.759_DP,  ! Ar
     .                4.341_DP,  ! K
     .                6.113_DP/) ! Ca

        write (iunout,'(a6,a12)') '  JATM','        Eion'
        do jatm = 1, natm
          if(nchara(jatm).eq.1) then
             Eion(jatm)=EionH
             if (nmassa(jatm).ne.1) Eion(jatm)=Eion(jatm)*
     .         (nmassa(jatm)*(pmassa+pmasse))/
     .         (pmasse+nmassa(jatm)*pmassa)
          elseif (nchara(jatm).eq.2) then
            Eion(jatm)=EionHe
          elseif (nchara(jatm).le.npot) then
            Eion(jatm)=pot_data(nchara(jatm))
          else
            Eion(jatm)=0.0_dp
          endif
          write (iunout,'(i6,es12.4)') jatm, eion(jatm)
        enddo

        deallocate (pot_data)
      endif

      eion_set=.true.
      return

      end subroutine eirene_init_eion
