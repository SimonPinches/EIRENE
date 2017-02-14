      SUBROUTINE EIR_SHOW_GIT_INFO(LOUT)
      implicit none
      integer, intent(in) :: LOUT

C................. Write out GIT specific configuration ................
C

      WRITE (LOUT,*) ""
      WRITE (LOUT,*) "EIRENE REPOSITORY STATUS DURING COMPILATION :-"
      WRITE (LOUT,*) "----------------------------------------------"

      WRITE(LOUT,"(a,a)") " Current GIT repository:  ",
     &__GITREPOSITORY
      WRITE(LOUT,"(a,a)") " Current GIT release tag: ",
     &__GITRELEASETAG
      WRITE(LOUT,"(a,a)") " Current GIT branch:      ",
     &__GITBRANCHNAME
      WRITE(LOUT,"(a,a)") " Last commit SHA1-key:    ",
     &__GITSHAKEY
      WRITE(LOUT,"(a,a)") " Hostname at compilation: ",
     &__HOSTNAME
      if (__GITSTATUSSTRING.eq."modyf") then
         WRITE (LOUT,*) ""
         WRITE (LOUT,*) "********************************************"
         WRITE (LOUT,*) "***  WARNING:                            ***"
         WRITE (LOUT,*) "***  There were uncommitted changes in   ***"
         WRITE (LOUT,*) "***  the repository during compilation.  ***"
         WRITE (LOUT,*) "***  This executable may not be          ***"
         WRITE (LOUT,*) "***  reproducible by it's SHA1-key!      ***"
         WRITE (LOUT,*) "********************************************"
      endif

      END SUBROUTINE EIR_SHOW_GIT_INFO
C
      SUBROUTINE EIR_SHOW_GIT_INFO_SHORT(LOUT)
C***********************************************************************
C   AUTHOR  :  Derek Harting (d.harting@fz-juelich.de)
C
C   PURPOSE :  Print out short format GIT repository status like repository
C              name, last commit SHA1-key and if the repsository was up to date.
C              The information is passed during compilation time via
C              the preprocessor in the __GIT variables.
C
C   INPUT   :   LOUT    - Unit number for output
C
C***********************************************************************
      implicit none
      INTEGER, intent(in) :: LOUT
C
      INTEGER :: LENSTR
      CHARACTER :: CSTR*70

C
C................. Write out GIT specific configuration ................
C
C
      WRITE(LOUT,"(a,a)") " EIRENE GIT repository : ",
     &  __GITREPOSITORY
      CSTR=__GITSHAKEY
      if (__GITSTATUSSTRING.eq."modyf") then
        CSTR=CSTR(1:LENSTR(CSTR))//" ( + uncommitted changes !! )"
      ENDIF
      WRITE(LOUT,"(a,a)") " EIRENE SHA1-key       : ",CSTR

      END  SUBROUTINE EIR_SHOW_GIT_INFO_SHORT
