!     Function to retrieve the $SOLPSTOP environment
!     variable, or the global default

      function get_solpstop() result(solpstop)
      use eirmod_comprt
      implicit none

!     Code taken from B2.5 open_file.F

      ! internal
      character*256 solpstop
#ifndef NO_GETENV
      integer lenval, ierror
#ifndef USE_PXFGETENV
      intrinsic get_environment_variable
#endif
#endif
      logical first
      save first
      data first /.true./


!      if(first) then
        solpstop=' '
#ifndef NO_GETENV
#ifdef USE_PXFGETENV
        CALL PXFGETENV ('SOLPSTOP', 0, solpstop, lenval, ierror)
#else
        call get_environment_variable('SOLPSTOP',
     .   status=ierror, length=lenval)
        if (ierror.eq.0) call get_environment_variable('SOLPSTOP',
     .   value=solpstop)
#endif
#endif
        if(solpstop.eq.' ') then
          if (first) 
     >     write(IUNOUT,'(A,A)') 'get_solpstop: no SOLPSTOP set'
        else
          if (first) write(IUNOUT,'(A,A)')
     1     'get_solpstop: SOLPSTOP = ',trim(solpstop)
        endif
        first=.false.
!      endif

      return
      end function get_solpstop
