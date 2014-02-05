      subroutine eirene_bfield (ind, x, y, z, bx, by, bz, bf)

      use eirmod_precision
      use eirmod_parmmod
      use eirmod_comusr
      use eirmod_cinit

      implicit none

      integer, intent(in) :: ind
      real(dp), intent(in) :: x, y, z
      real(dp), intent(out) :: bx, by, bz, bf
      real(dp) :: eirene_femint, bni
      integer :: eirene_idez

      IF (INDPRO(5) == 8) THEN

         CALL EIRENE_VECUSR (1,BX,BY,BZ,1)
         bf = 1.

      ELSE IF (LBSMO) THEN

         bx = eirene_Femint(bxincorner, ind, x, y, z, .false.)
         by = eirene_Femint(byincorner, ind, x, y, z, .true.)
         bz = eirene_Femint(bzincorner, ind, x, y, z, .true.)
         bf = eirene_Femint(bfincorner, ind, x, y, z, .true.)
         bni = 1._dp/sqrt(bx*bx + by*by + bz*bz)
         bx = bx * bni
         by = by * bni
         bz = bz * bni

      ELSE

         BX=BXIN(IND)
         BY=BYIN(IND)
         BZ=BZIN(IND)
         BF=BFIN(IND)

      END IF
      
      return
      end subroutine eirene_bfield
