      subroutine post_couple(ncv_b2,nfc_b2,ns_b2,nstrai_b2,nplsh)
      use eirmod_PARMMOD
      use eirmod_CGRID
      use eirmod_CPOLYG
      use eirmod_COMUSR
      use eirmod_COMSOU
      use eirmod_CCOUPL
      use eirmod_COMNNL
      use eirmod_COMPRT
      implicit none
      integer ncv_b2,nfc_b2,ns_b2,nstrai_b2,nplsh
      logical error
      external eirene_exit_own

      error=.false.
      if (levgeo.eq.4) then
        if(ncv_b2.ne.ncvp) then
          write(iunout,*) 'NCV(B2) <> NCV(EIRENE) ',ncv_b2,ncvp
          error=.true.
        endif
        if(nfc_b2.ne.nfcp) then
          write(iunout,*) 'NFC(B2) <> NFC(EIRENE) ',nfc_b2,nfcp
          error=.true.
        endif
      else
        write(iunout,*) 'UNEXPECTED VALUE OF LEVGEO: ', levgeo
      endif
      if(ns_b2.ne.natma+nflb-nplsh) then
        write(iunout,*) 'NS(B2) <> NS(EIRENE) '
        write(iunout,*) 'NS(B2)', ns_b2
        write(iunout,*) 'NS(EIRENE) (standard run)', natma+nflb
        write(iunout,*) 'NS(EIRENE) (spatially hybrid)', natma+nflb-nplsh
        error=.true.
      endif
      if(ntime.gt.0) then
        if(nstrai_b2+1.ne.nstrai) then
          write(iunout,*) 'NSTRA(B2) <> NSTRA(EIRENE) ',nstrai_b2+1,nstrai
          write(iunout,*) 'INCLUDING CENSUS STRATUM'
          error=.true.
        endif
      else
        if(nstrai_b2.ne.nstrai) then
          write(iunout,*) 'NSTRA(B2) <> NSTRA(EIRENE) ',nstrai_b2,nstrai
          write(iunout,*) 'NO CENSUS STRATUM DETECTED'
          error=.true.
        endif
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
      end subroutine post_couple

!!!Local Variables:
!!! mode: f90
!!! End:
