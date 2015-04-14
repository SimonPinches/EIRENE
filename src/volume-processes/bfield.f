cdr  sept 2014:  comments added
c    provide carthesian local magnetic field unit vector bx,by,bz,
c    as well as b-field strength bf (Tesla)
c    at point x,y,z, in cell icell

c  current options:
c  default    :  use input background tallies. B-field is constant per cell
c  indpro(5)=8:  user provided B-field
c  LBSMO (?)  :  apparently: only in case of levgeo=4, interpolation in triangles

      subroutine eirene_bfield (icell, x, y, z, bx, by, bz, bf)

      use eirmod_precision
      use eirmod_parmmod
      use eirmod_comusr
      use eirmod_cinit

      implicit none

      integer, intent(in) :: icell
      real(dp), intent(in) :: x, y, z
      real(dp), intent(out) :: bx, by, bz, bf
      real(dp) :: eirene_femint, bni
      integer :: eirene_idez

      IF (INDPRO(5) == 8) THEN

         CALL EIRENE_VECUSR (1,BX,BY,BZ,1)
         bf = 1.

      ELSE IF (LBSMO) THEN

c  if not levgeo.eq.4  error exit

         bx = eirene_Femint(bxincorner, icell, x, y, z, .false.)
         by = eirene_Femint(byincorner, icell, x, y, z, .true.)
         bz = eirene_Femint(bzincorner, icell, x, y, z, .true.)
         bf = eirene_Femint(bfincorner, icell, x, y, z, .true.)
         bni = 1._dp/sqrt(bx*bx + by*by + bz*bz)
         bx = bx * bni
         by = by * bni
         bz = bz * bni

      ELSE

         BX=BXIN(ICELL)
         BY=BYIN(ICELL)
         BZ=BZIN(ICELL)
         BF=BFIN(ICELL)

      END IF
      
      return
      end subroutine eirene_bfield
