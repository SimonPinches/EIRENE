      subroutine eirene_find_emis_line (istr, ichori, ener_ch, lno)

cdr documentation ? comments ?  ongoing....
cdr
c  input:   istr  :  stratum number
c           ichori:  chord number of line of sight
c           ener_ch  :  energy parameter for identifying a particular 
C                                  emission profile by "trans_en"
c           ch_line_name(ichori)): optional: identify a particular
c                                  emission profile by "name of line" 
c
c  output:  lno   :  "line number", i.e. the volumetric emission profile.

cdr 
cdr currently called from SIGLINE in diagno block. 
cdr    
cdr SIGLINE is called only for chords ICHORI, for which NCHTAL(ichori)=2.
cdr Calls are whenever either a change in stratum number ISTR, 
cdr                    or                 transition energy ENER,
cdr                    or                 a new internal iteration 
cdr  (time stepping, nonlinear BGK iterations)

cdr In storage saving mode (mod_addv = 0) from here we call EIRENE_EMISSIVITY.F
cdr   for the present stratum and present line, 
cdr   to fill the proper tallies ADDV with emissivities.
cdr    (similar to former Ba_alpha.f,...etc.)
cdr In full storage mode (mod_addv > 0) the relevant ADDV tallies have already
cdr   been set in calls to EIRENE_EMISSIVITY from the stratum-loop in MCARLO.f
cdr 
cdr May 18:  try to identify the line LNO,
cdr          or (old options) as specified by input flags
cdr                                       ICHORI  (line of sight number)
cdr                                       ENER    (flag for selecting a particular line)

cdr          this is done by trying to find a match of 'ch_line_name(ichori)'
cdr          read from block 12 for chord ICHORI
cdr          with 'emis_lines(i)%line_name'
cdr          If this is not successful, then we try to use the ENER identifier,
cdr          which may have been read in input for emission profiles
cdr          or set from old default hydrogenic line models (setup_default_emissivity).
cdr          The latter is done whenever there are chords and no emissivities
cdr          have been read from external files.

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
      USE EIRMOD_COMPRT

      implicit none

      integer, intent(in) :: istr, ichori
      real(dp), intent(in) :: ener_ch
      integer, intent(out) :: lno
      real(dp) :: ener_il
      integer :: iline
      character(len=80) :: ctest1, ctest2
      logical :: found

      lno = 0
      found = .false.

      if (ichori.gt. 0) then
cdr  for ichori <= 0: bypass ichori,
cdr  to set volumetric emission profiles on ADDV array, also without any chords.

      if (len_trim(ch_line_name(ichori)) > 0) then
! find corresponding line from line names
        ctest1 = ch_line_name(ichori)

        do iline = 1, num_lines
           ctest2 = emis_lines(iline)%line_name
           if (ctest1 == ctest2) then
             found = .true.
             lno = iline
             exit
           end if
        end do
CDR
        IF (FOUND) THEN
          WRITE (IUNOUT,*) 'EMISSION LINE identified by NAME,',
     .                     ' iline=',LNO
          WRITE (IUNOUT,*) 'name: ',trim(CTEST2)   ! =CTEST1
        ELSE
          WRITE (IUNOUT,*) 'NO EMISSION LINE identified by NAME'
          WRITE (IUNOUT,*) 'name: ',trim(CTEST1)
        ENDIF
      else
        WRITE (IUNOUT,*) 'NO EMISSION LINE identified by NAME'
        WRITE (IUNOUT,*) 'name: ', 'not specified'

      end if  ! ch_line_name(ichori) defined

      endif   !  ichori.gt.0

      if (.not.found) then
! check transition energies, for backward compatibility with old input block 12.
        do iline = 1, num_lines
          ener_il = emis_lines(iline)%trans_en
          if (abs((ener_ch-ener_il)/(ener_il+eps10)) <= eps5) then
            found = .true.
            lno = iline
            ctest1 = adjustl(trim(emis_lines(iline)%line_name))
            exit
          end if
        end do
CDR
        IF  (FOUND) THEN
          WRITE (IUNOUT,*) 'EMISSION LINE identified by ENERGY,',
     .                     ' iline=',LNO
          CALL EIRENE_MASR1('ENERGY: ',ener_il)
          WRITE (IUNOUT,*) 'name: ',trim(ctest1)   ! =cname_ch
        ELSE
          WRITE (IUNOUT,*) 'NO EMISSION LINE identified by ENERGY'
          CALL EIRENE_MASR1('ENERGY: ',ENER_CH)
        endif

      end if

      if (.not.found) then
        write (iunout,*) ' NO MATCHING EMISSION LINE FOUND FOR CHORD ',
     .                     ICHORI
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
     .             NSIGI_SPC,TRCFLE)
      ELSEIF ((NFILEN.EQ.6.OR.NFILEN.EQ.7).AND.ISTR.EQ.0) THEN
        IESTR=ISTR
        CALL EIRENE_RSTRT(ISTR,NSTRAI,NESTM1,NESTM2,NADSPC,
     .             ESTIMV,ESTIMS,ESTIML,
     .             NSDVI1,SDVI1,NSDVI2,SDVI2,
     .             NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .             NSIGI_SPC,TRCFLE)
      ELSE
        WRITE (IUNOUT,*) 'ERROR IN FIND_EMIS_LINE: ' //
     .                   'DATA FOR STRATUM ISTRA= ', ISTR
        WRITE (IUNOUT,*) 'ARE NOT AVAILABLE. FIND_EMIS_LINE ABANDONED'
        RETURN
      ENDIF
C
      if (mod_addv == 0) then
c storage saving mode: 
c ADDV is overwritten when a new line comes, within a run.
c Thus re-calculate the new emissivity profile on ADDV now
         call eirene_emissivity(istr, lno, lno, 1)
c     else
c Sufficiently large storage on ADDV additional tally array,
c for all lines and components. No need to reset ADDV tallies.
      end if

      return
      end subroutine eirene_find_emis_line
