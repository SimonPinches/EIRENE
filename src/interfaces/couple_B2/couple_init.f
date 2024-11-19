cpg called from eirene.f in initialization phase,
cpg when eirene is in "coupled mode": i.e. IF(NMODE.NE.0)
C

      subroutine eirene_couple_init_output

      use EIRMOD_CPES, only: MY_PE, NPRS
      use EIRMOD_COMPRT, only: IUNOUT
      use EIRMOD_PARMMOD, only: LOUTAPP
      use EIRMOD_OPENFILE, only: EIRENE_OPENFILE

      implicit none

      character(20) :: outname
      character(6) :: outpos

      if (NPRS > 1) IUNOUT = 7

      if (my_pe.ne.0) then

        OUTNAME='output.'
        WRITE (OUTNAME(8:),'(I4.4)') MY_PE
        IF ( LOUTAPP ) THEN
          OUTPOS='APPEND'
        ELSE
          OUTPOS='ASIS'
        END IF

        CALL EIRENE_OPENFILE (IUNIT=IUNOUT, FILE=OUTNAME,
     .        ACCESS='SEQUENTIAL', FORM='FORMATTED', POSITION=OUTPOS)
cpg #ifdef WINDOWS
cpg DIR$attributes c, alias: 'ioflush_' :: ioflush
cpg #endif
        call ioflush_usr
cpg end

      end if
      return

      end subroutine eirene_couple_init_output


      subroutine eirene_couple_alloc
      return
      end subroutine eirene_couple_alloc

      subroutine eirene_couple_post_input
      return
      end subroutine eirene_couple_post_input
