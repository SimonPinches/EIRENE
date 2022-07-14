!     Function to retrieve the $SOLPSTOP environment
!     variable, or the global default

      function get_solpstop() result(solpstop)
      use eirmod_comprt
      implicit none

!     Code taken from B2.5 open_file.F

      ! internal
      character*256 solpstop
#ifdef USE_PXFGETENV
      integer lenval, ierror
#else
#ifdef NAGFOR
      integer lenval, ierror
#endif
#endif
      logical first
      save first
      data first /.true./


!      if(first) then
#ifdef NO_GETENV
        solpstop=' '
#else
#ifdef NAGFOR
        call get_environment_variable('SOLPSTOP', 
     .   status=ierror, length=lenval)
        if (ierror.eq.0) then
          call get_environment_variable('SOLPSTOP',value=solpstop)
        elseif (ierror.eq.1) then
          solpstop=' '
        endif
#else
#ifdef USE_PXFGETENV
        CALL PXFGETENV ('SOLPSTOP', 0, solpstop, lenval, ierror)
#else
        call getenv ('SOLPSTOP', solpstop)
#endif
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


      end 
