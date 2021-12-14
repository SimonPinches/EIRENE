c  14.5.06:  bug fix: 1 line added: if nchtal.ne.1 and. nchtal.ne.3:  cycle
C  oct.14.  variance tallies corrected
cpb  Dec. 2017: remove type SPECT_ARRAY, not needed in Fortran 2003
cdr  May 18: Probably here we use the data structure traj(i..) for storing the
cdr  line of sight, e.g for scoring spectra along lines of sight?
cdr  This same data structure is (or was) probably also used
cdr  for an unfinished correlated sampling option.
cdr  In either case it may not be complete any more.

cdr: tbd: Try to document status and purpose.  started....

cdr  called from: ... if ...


      subroutine EIRENE_setup_chord_spectra

      use EIRMOD_precision
      use EIRMOD_parmmod
      use EIRMOD_cestim
      use EIRMOD_comsig
      use EIRMOD_comprt
      use EIRMOD_cupd
      use EIRMOD_ccona
      use EIRMOD_LININT, ONLY: EIRENE_LININT

      implicit none

      real(dp) :: c1(3), c2(3), PSIG(0:NSPZ+10)
      real(dp) :: ze, timax
      integer :: ichori, ifirst, ichrd, ipvot, nbc2, nac2,
     .           ispc, ntot_cell, ntotsp, iprtyp
      type(eirene_spectrum), allocatable :: svestiml(:), svsmestl(:)
      TYPE(EIRENE_SPECTRUM), POINTER :: ESPEC, SSPEC
      TYPE(CELL_INFO), POINTER :: FIRST, CUR

!  FIND CELLS INTERSECTED BY CHORDS

      ze = 1._dp

      ifirst = -1
      ntot_cell = 0

      do ichori = 1, nchori

        IF (.NOT.NLSTCHR(ICHORI)) CYCLE
        IF ((NCHTAL(ICHORI) /= 1) .AND. (NCHTAL(ICHORI) /= 3) .AND.
     .      (NCHTAL(ICHORI) /= 4) ) CYCLE

        IPVOT=IPIVOT(ICHORI)
        C1(1)=XPIVOT(ICHORI)
        C1(2)=YPIVOT(ICHORI)
        C1(3)=ZPIVOT(ICHORI)
C
        ICHRD=ICHORD(ICHORI)
        C2(1)=XCHORD(ICHORI)
        C2(2)=YCHORD(ICHORI)
        C2(3)=ZCHORD(ICHORI)
C
        NBC2=NSPBLC(ICHORI)
        NAC2=NSPADD(ICHORI)

        ALLOCATE(TRAJ(ICHORI)%TRJ)
        TRAJ(ICHORI)%TRJ%P1 = C1
        TRAJ(ICHORI)%TRJ%P2 = C2
        TRAJ(ICHORI)%TRJ%NCOU_CELL = 0
        NULLIFY(TRAJ(ICHORI)%TRJ%CELLS)

        CALL EIRENE_LININT
     .  (IFIRST,ICHORI,C1,C2,ICHRD,IPVOT,NBC2,NAC2,ZE,
     .               PSIG,TIMAX,1,1,1,IABS(NCHENI))

cdr probably: this call to linint provides ncou_cell,
cdr           the total number of cells visited by chord no. ICHORI
        ntot_cell = ntot_cell + traj(ichori)%trj%ncou_cell

      end do

      IF (NTOT_CELL == 0) RETURN

!dr  there are 'NTOT_CELL' directional CELL-BASED SPECTRA TO BE ADDED TO SPECTRUM TALLIES

!  SAVE SPECTRA SPECIFIED VIA INPUT

      IF (NADSPC > 0) THEN
C  SAVE ESTIML, SMESTL,...
        IF (ALLOCATED(ESTIML)) THEN
          ALLOCATE(SVESTIML(NADSPC))
          DO ISPC = 1, NADSPC

            SVESTIML(ISPC) = ESTIML(ISPC)

            DEALLOCATE (ESTIML(ISPC)%SPC)
            IF (ASSOCIATED(ESTIML(ISPC)%SDV)) THEN
              DEALLOCATE (ESTIML(ISPC)%SDV)
              DEALLOCATE (ESTIML(ISPC)%SGM)
              DEALLOCATE (ESTIML(ISPC)%STV)
              DEALLOCATE (ESTIML(ISPC)%GG)
            END IF
          END DO
          DEALLOCATE(ESTIML)
        END IF

        IF (ALLOCATED(SMESTL)) THEN
          ALLOCATE(SVSMESTL(NADSPC))
          DO ISPC = 1, NADSPC

            SVSMESTL(ISPC) = SMESTL(ISPC)

            DEALLOCATE (SMESTL(ISPC)%SPC)
            IF (ASSOCIATED(SMESTL(ISPC)%SDV)) THEN
              DEALLOCATE (SMESTL(ISPC)%SDV)
              DEALLOCATE (SMESTL(ISPC)%SGM)
              DEALLOCATE (SMESTL(ISPC)%STV)
              DEALLOCATE (SMESTL(ISPC)%GG)
            END IF
          END DO
          DEALLOCATE(SMESTL)
        END IF

      END IF

cdr  set up (allocate) additional tally arrays for directional cell-based spectra
cdr  SMESTL


      NTOTSP = NADSPC + NTOT_CELL

      ALLOCATE(ESTIML(NTOTSP))
      DO ISPC = 1, NADSPC
        ESTIML(ISPC) = SVESTIML(ISPC)
      END DO

      IF (ALLOCATED(SVSMESTL).or.NSMSTRA.GT.0) THEN
        ALLOCATE(SMESTL(NTOTSP))
        DO ISPC = 1, NADSPC
          SMESTL(ISPC) = SVSMESTL(ISPC)
        END DO
      END IF

!  add spectra for cells along chords

      ispc = nadspc
      do ichori = 1,nchori

        IF ((NCHTAL(ICHORI) /= 1) .AND. (NCHTAL(ICHORI) /= 3) .AND.
     .      (NCHTAL(ICHORI) /= 4) ) CYCLE

cdr spectrally resolved lines of sight tallies for options 1,3 and 4 ??


         if (.not.associated(traj(ichori)%trj%cells)) cycle
         first => traj(ichori)%trj%cells
         cur => first

         if (nchtal(ichori) == 1) then
           iprtyp = 1
         else if (nchtal(ichori) == 3) then
           iprtyp = 0
         else if (nchtal(ichori) == 4) then
           iprtyp = 1
         else
           write (iunout,*) ' wrong chord type for spectrum '
           write (iunout,*) ' no spectrum set up for cell',cur%no_cell
           cycle
         end if

!  loop over all cells along trajectory
         do

           ispc = ispc + 1
           espec => estiml(ispc)

           espec%isrfcll = 2 ! SURFACE OR CELL-BASED OR DIRECTIONAL CELL-BASED
           espec%ispcsrf = cur%no_cell
           espec%iprtyp = iprtyp   !TYP
           espec%iprsp = nspspz(ichori)  !SPECIES
           espec%ispctyp = 1  ! TYPE OF SPECTRUM AMP/EV, OR WATT/EV,  ETC...
           espec%nspc = abs(ncheni)
           espec%imetsp = 0
           espec%idirec = 1
           if (ncheni > 0) then
             espec%spcmin = emin1(ichori)
             espec%spcmax = emax1(ichori)
             espec%log = .false.
           else
             espec%spcmin = log10(emin1(ichori))
             espec%spcmax = log10(emax1(ichori))
             espec%log = .true.
           end if
           espec%esp_00 = 0._dp
           espec%spc_xplt = 0._dp
           espec%spc_yplt = 0._dp
           espec%spc_same = 0._dp
           espec%spcvx = traj(ichori)%trj%vx
           espec%spcvy = traj(ichori)%trj%vy
           espec%spcvz = traj(ichori)%trj%vz
           espec%esp_min =1.e30_dp
           espec%esp_max = -1.e30_dp
           espec%spcdel=(espec%spcmax-espec%spcmin)/real(espec%nspc,dp)
           espec%spcdeli = 1._dp / (espec%spcdel+eps60)

           allocate(espec%spc(0:espec%nspc+1))

           allocate(espec%sdv(0:espec%nspc+1))
           allocate(espec%sgm(0:espec%nspc+1))
           allocate(espec%stv(0:espec%nspc+1))
           allocate(espec%gg(0:espec%nspc+1))

           espec%spc(0:espec%nspc+1) = 0


           if (allocated(smestl)) then
C SUM OVER STRATA SPECTRA TALLIES

             sspec => smestl(ispc)

             allocate(sspec%spc(0:espec%nspc+1))

             allocate(sspec%sdv(0:espec%nspc+1))
             allocate(sspec%sgm(0:espec%nspc+1))
             allocate(sspec%stv(0:espec%nspc+1))
             allocate(sspec%gg(0:espec%nspc+1))

             smestl(ispc) = estiml(ispc)
           end if

           cur%no_spect = ispc
           cur => cur%nextc
!pb associated with two arguments tests if both arguments point to the same target
           if (associated(cur,first)) exit
         end do
      end do

      NADSPC = ISPC

      end subroutine EIRENE_setup_chord_spectra
