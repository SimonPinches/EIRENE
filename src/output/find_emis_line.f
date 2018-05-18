      subroutine eirene_find_emis_line (ist, ichori, ener, lno)

      use eirmod_precision
      use eirmod_parmmod
      use eirmod_comsig
      use eirmod_ccona
      use eirmod_comsou
      use eirmod_comusr
      use eirmod_cgeom
      use eirmod_cgrid
      use eirmod_cspei
      use eirmod_ctrcei
      USE EIRMOD_CESTIM
      USE EIRMOD_CSDVI
      USE EIRMOD_CSDVI_BGK
      USE EIRMOD_CSDVI_COP
      USE EIRMOD_COMPRT

      implicit none
      
      integer, intent(in) :: ist, ichori
      real(dp), intent(in) :: ener
      integer, intent(out) :: lno
      integer :: i
      character(len=:), allocatable :: ctest1, ctest2
      logical :: found

      lno = 0
      found = .false.

      if (len_trim(ch_line_name(ichori)) > 0) then
! find corresponding line from line names
        ctest1 = adjustl(trim(ch_line_name(ichori)))

        do i = 1, num_lines
           ctest2 = adjustl(trim(emis_lines(i)%line_name))
           if (ctest1 == ctest2) then
             found = .true.
             lno = i
             exit
           end if
        end do
      end if

      if (.not.found) then
! check energies
        do i = 1, num_lines
          if (abs((ener-emis_lines(i)%energy)/emis_lines(i)%energy) <= 
     .        eps5) then
            found = .true.
            lno = i
            exit
          end if
        end do
      end if

      if (.not.found) then
        write (iunout,*) ' NO MATCHING EMISSION LINE FOUND FOR CHORD ',
     .                   ICHORD
        return
      end if

      IF (IESTR.EQ.IST) THEN
C  NOTHING TO BE DONE
      ELSEIF (NFILEN.EQ.1.OR.NFILEN.EQ.2) THEN
        IESTR=IST
        CALL EIRENE_RSTRT(IST,NSTRAI,NESTM1,NESTM2,NADSPC,
     .             ESTIMV,ESTIMS,ESTIML,
     .             NSDVI1,SDVI1,NSDVI2,SDVI2,
     .             NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .             NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .             NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .             NSIGI_SPC,TRCFLE)
      ELSEIF ((NFILEN.EQ.6.OR.NFILEN.EQ.7).AND.IST.EQ.0) THEN
        IESTR=IST
        CALL EIRENE_RSTRT(IST,NSTRAI,NESTM1,NESTM2,NADSPC,
     .             ESTIMV,ESTIMS,ESTIML,
     .             NSDVI1,SDVI1,NSDVI2,SDVI2,
     .             NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .             NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .             NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .             NSIGI_SPC,TRCFLE)
      ELSE
        WRITE (IUNOUT,*) 'ERROR IN EMIS_PROFILES: ' // 
     .                   'DATA FOR STRATUM ISTRA= ', IST
        WRITE (IUNOUT,*) 'ARE NOT AVAILABLE. EMIS_PROFILES ABANDONNED'
        RETURN
      ENDIF
C      
      if (mod_addv == 0) then
! ADDV is always overwritten thus recalculate the emissivity profile
         call eirene_emissivity(ist, lno, lno)
      end if

      return
      end subroutine eirene_find_emis_line
