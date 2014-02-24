c
c
      function ranget_eirene(ise)
      USE EIRMOD_PRECISION
      USE EIRMOD_CLOGAU
      implicit none
      integer, intent(inout) :: ise
 
      integer :: iseed
      common /cmem/ iseed
      real(dp) :: ran, ranf_eirene
      integer :: ranget_eirene, ifirst, large
      data ifirst/0/
      save large

      IF (NLOLDRAN) THEN
         if (ifirst.eq.0) then
            large=HUGE(large)   !vk 1->large
            ifirst=1
         endif
         ran=ranf_eirene()
         ise=ran*large
         ranget_eirene=ise
      ELSE
         ranget_eirene=iseed
         ise = iseed
      END IF
      return
      end
