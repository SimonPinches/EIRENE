      subroutine eirene_find_emis_line (istr, ichori, ener, lno)

cdr documentation ? comments ?

cdr may 18:  try to identify the line LNO,
cdr          als specified by input flags 
cdr                                       ICHORI  (line of sight number)
cdr                                       ENER    (flag for selecting a particular line)

cdr          this is done my trying to find a match of 'ch_line_name(ichori)'
cdr          read from block 12 for chord ICHORI
cdr          with 'emis_lines(i)%line_name'

cdr          Then fill the appropriate additional tallies ADDV
cdr          with the needed volumetric line emissivities, for stratum ISTR,
cdr          by calling  EIRENE_EMISSIVITY(...)

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
      
      integer, intent(in) :: istr, ichori
      real(dp), intent(in) :: ener
      integer, intent(out) :: lno
      real(dp) :: ener_il
      integer :: iline
      character(len=:), allocatable :: ctest1, ctest2
      logical :: found

      lno = 0
      found = .false.

      if (len_trim(ch_line_name(ichori)) > 0) then
! find corresponding line from line names
        ctest1 = adjustl(trim(ch_line_name(ichori)))

        do iline = 1, num_lines
           ctest2 = adjustl(trim(emis_lines(iline)%line_name))
           if (ctest1 == ctest2) then
             found = .true.
             lno = iline
             exit
           end if
        end do
      end if

      if (.not.found) then
! check energies, for backward compatibility with old input block 12.
        do iline = 1, num_lines
          ener_il = emis_lines(iline)%energy
          if (abs((ener-ener_il)/ener_il) <= eps5) then
            found = .true.
            lno = iline
            exit
          end if
        end do
      end if

      if (.not.found) then
        write (iunout,*) ' NO MATCHING EMISSION LINE FOUND FOR CHORD ',
     .                   ICHORD
        return
      end if

      IF (IESTR.EQ.ISTR) THEN
C  NOTHING TO BE DONE
      ELSEIF (NFILEN.EQ.1.OR.NFILEN.EQ.2) THEN
        IESTR=ISTR
        CALL EIRENE_RSTRT(ISTR,NSTRAI,NESTM1,NESTM2,NADSPC,
     .             ESTIMV,ESTIMS,ESTIML,
     .             NSDVI1,SDVI1,NSDVI2,SDVI2,
     .             NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .             NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .             NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .             NSIGI_SPC,TRCFLE)
      ELSEIF ((NFILEN.EQ.6.OR.NFILEN.EQ.7).AND.ISTR.EQ.0) THEN
        IESTR=ISTR
        CALL EIRENE_RSTRT(ISTR,NSTRAI,NESTM1,NESTM2,NADSPC,
     .             ESTIMV,ESTIMS,ESTIML,
     .             NSDVI1,SDVI1,NSDVI2,SDVI2,
     .             NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .             NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .             NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .             NSIGI_SPC,TRCFLE)
      ELSE
        WRITE (IUNOUT,*) 'ERROR IN FIND_EMIS_LINE: ' // 
     .                   'DATA FOR STRATUM ISTRA= ', ISTR
        WRITE (IUNOUT,*) 'ARE NOT AVAILABLE. FIND_EMIS_LINE ABANDONNED'
        RETURN
      ENDIF
C      
      if (mod_addv == 0) then
c ADDV is overwritten when a new line comes, within a run.
c Thus recalculate the new emissivity profile on ADDV
         call eirene_emissivity(istr, lno, lno)
c     else
c Sufficiently large storage on ADDV additional tally array, 
c for all lines and components. No need to reset ADDV tallies.    
      end if

      return
      end subroutine eirene_find_emis_line
