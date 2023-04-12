      subroutine post_couple(nx_b2,ny_b2,ns_b2,nstrat_b2)
      use eirmod_PARMMOD
      use eirmod_CGRID
      use eirmod_CPOLYG
      use eirmod_COMUSR
      use eirmod_COMSOU
      use eirmod_CCOUPL
      use eirmod_COMNNL
      use eirmod_COMPRT
      use eirmod_wneutrals
      implicit none
      integer nx_b2,ny_b2,ns_b2,nstrat_b2
      logical error
      error=.false.
      if (levgeo.le.3) then
        if(nx_b2.ne.np2nd-npplg) then
          write(iunout,*) 'NX(B2) <> NX(EIRENE) ',nx_b2,np2nd-npplg
          error=.true.
        endif
        if(ny_b2.ne.nr1st-1) then
          write(iunout,*) 'NY(B2) <> NY(EIRENE) ',ny_b2,nr1st-1
          error=.true.
        endif
      elseif (levgeo.eq.4) then
        if(nx_b2.ne.np2tal_save-npplg) then
          write(iunout,*) 'NX(B2) <> NX(EIRENE) ',nx_b2,np2tal_save-npplg
          error=.true.
        endif
        if(ny_b2.ne.nr1tal_save-1) then
          write(iunout,*) 'NY(B2) <> NY(EIRENE) ',ny_b2,nr1tal_save-1
          error=.true.
        endif
      else
        write(iunout,*) 'UNEXPECTED VALUE OF LEVGEO: ', levgeo
      endif
      if(ns_b2.ne.natma+nflb) then
        write(iunout,*) 'NS(B2) <> NS(EIRENE) ',ns_b2,natma+nflb
        error=.true.
      endif
      if(nstrat_b2.ne.nstrai) then
        write(iunout,*) 'NSTRA(B2) <> NSTRA(EIRENE) ',nstrat_b2,nstrai
        error=.true.
      endif
      if(error) then
        write(iunout,*) 'ERROR IN POST-COUPLING PARAMETERS'
        write(iunout,*) 'POST-COUPLING TEST FAILED'
        write(iunout,*) 'THE ACTUALLY USED ARRAY SIZES DO NOT MATCH'
        write(iunout,*) 'CHECK INPUT FILES'
!PB     call xerrab ('Error in post-coupling parameters')
        write(iunout,*) 'Error in post-coupling parameters'
        call eirene_exit_own(1)
      else
        write(iunout,*) 'POST-COUPLING TEST PASSED'
      endif
      return
      end

!!!Local Variables:
!!! mode: f90
!!! End:
