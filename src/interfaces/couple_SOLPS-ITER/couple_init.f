cpg called from eirene.f in initialization phase,
cpg when eirene is in "coupled mode": i.e. IF(NMODE.NE.0)
C

      subroutine eirene_couple_init_output

       use EIRMOD_CPES, only: MY_PE, NPRS
       use EIRMOD_COMPRT, only: IUNOUT
       use EIRMOD_PARMMOD, only: LOUTAPP
      
       implicit none
       

      character(20) :: outname
      character(6) :: outpos
      
#ifndef B25_EIRENE
      if (NPRS > 1) IUNOUT = 7
#endif  
         if (my_pe.ne.0) then
     
            OUTNAME='output.'
            write (OUTNAME(8:),'(I4.4)') MY_PE
#ifndef NAGFOR
            if ( LOUTAPP ) THEN
              OUTPOS='APPEND'
            else
              OUTPOS='ASIS'
            end if
            
            open (unit=IUNOUT, file=OUTNAME, access='SEQUENTIAL',
     .        form='FORMATTED', position=OUTPOS)
#else
            open (unit=IUNOUT, file=OUTNAME, access='SEQUENTIAL',
     .        FORM='FORMATTED')
#endif
cpg #ifdef WINDOWS
cpg DIR$attributes c, alias: 'ioflush_' :: ioflush
cpg #endif
            call ioflush_usr
cpg end
   
#ifdef B25_EIRENE 
          end if
#endif      
      end if
       return   
     
      end subroutine eirene_couple_init_output
      
      
      subroutine eirene_couple_alloc
        use eirmod_parmmod, only : nmol,nion,npls
        use eirmod_extrab25
      
!pb        CALL EIRENE_EXTRAB25_INIT_EION
        call eirene_extrab25_alloc_mods
        call eirene_extrab25_eirpbls_init(nmol,nion,npls)
      
      end subroutine eirene_couple_alloc
      
      subroutine eirene_couple_post_input
        use eirmod_precision
        use eirmod_ccona, only : ELCHA
        use eirmod_cpes, only : my_pe
        use eirmod_comsou, only : FLUX,nlvol,nlcns,NSTRAI
        use eirmod_comprt, only : ISTRA
        use eirmod_ccoupl, only : NTARGI
        use eirmod_extrab25, only : FLUX_SAVE
        
        integer :: ISTRAI
csw 24oct2011      
        if(my_pe==0) then
          if (nmode /= 0) then
            DO ISTRAI=NTARGI+1,NSTRAI
              ISTRA = ISTRAI
              if(.not.nlvol(istra) .and. .not. nlcns(istra)) then
                if (allocated(flux_save)) then !wd avoid issue with non-allocated flux_save in standalone runs
                  IF (FLUX_SAVE(ISTRA).NE.0._DP ) THEN
                    FLUX(ISTRA)=FLUX_SAVE(ISTRA)*ELCHA
                  else
                    FLUX(ISTRA)=1._DP
                  end if
                ELSE
                  FLUX(ISTRA)=1._DP
                ENDIF
              endif
            ENDDO
          endif
        endif

        return
      end subroutine eirene_couple_post_input
