c
c
      function ranset_eirene(ise)
      USE EIRMOD_PRECISION
      USE EIRMOD_CLOGAU
      implicit none
      integer, intent(in) :: ise
      integer :: iseed1, iseed2
 
      integer :: iseed
      common /cmem/ iseed
      real(dp) :: ranset_eirene, ranset_eirene_reinit
      integer, save :: ifirst=0
 
      IF (NLOLDRAN) THEN
         iseed1=ise
         iseed2=2000000-ise
CPB   iseed2=iseed+1
         call h1rnin(iseed1,iseed2)
c     ranset=ranfinit(iseed)      !VK for the new random number gen.
         ranset_eirene=0
      ELSE
         if (ise <= 0) then
            if (ifirst == 0) then
               iseed = 9876543
!     else  seed schon gesetzt, uebernimmt dies
            end if
         else
            iseed = abs(ise)
         end if
 
         ifirst = 1
         ranset_eirene=0
      END IF
      return
 
C     The following ENTRY is for reinitialization of EIRENE
 
      ENTRY ranset_eirene_reinit()
      ifirst = 0
      ranset_eirene_reinit = 0.D0
      return
      end
