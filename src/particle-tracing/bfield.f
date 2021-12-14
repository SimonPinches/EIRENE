cdr  jan. 2020:  internal consistency enforced for vector components,
cdr              and their smoothing/interpolation.
cdr              New: LBIN  flag:  all B field tallies exist.
cdr              NEW  LBSMO flag:  all B field tallies smoothed (interpolated)
cdr  aug. 2015:  logical flag L added. position x,y,z known (L=true)
cdr                           or else: use COM of cell ICELL
cdr  sept 2014:  comments added
c    provide cartesian local magnetic field unit vector bx,by,bz,
c    as well as B field strength bf (Tesla)
c    at point x,y,z, in cell icell

c  current options:
c  default    :  use input background tallies. B field is constant per cell
c  indpro(5)=8:  user-provided B field
c  LBSMO (?)  :  apparently: only in case of levgeo=4,5 interpolation in triangles, tetrahedra
cpb              switch used for interpolation of magnetic field input tally
cdr              from cell vertices to a local x,y,z point inside a cell.
cdr              Not fully available for all levgeo=1,2,3 optins. Check FEMINT.f

      subroutine eirene_bfield (icell, x, y, z, bx, by, bz, bf,l)

      use eirmod_precision, only: dp
      use eirmod_comusr, only: BXIN, BYIN, BZIN, BFIN, BXINCORNER, 
     >                         BYINCORNER, BZINCORNER, BFINCORNER, 
     >                         LBSMO, LBXIN, LBYIN, LBZIN, LBFIN
      use eirmod_cinit, only: INDPRO

      implicit none

      integer, intent(in) :: icell
      logical             :: l
      real(dp), intent(in) :: x, y, z
      real(dp), intent(out) :: bx, by, bz, bf
      real(dp) :: eirene_femint, bni

      IF (INDPRO(5) == 8) THEN
cdr  user defined B field. Units of Bx, By, Bz?
cdr  BF=1. ?  BF should be in Tesla.
cdr  L=true : spatial coordinates x,y,z are known here, 
cdr           i.e. call vecusr with L=true
cdr  or else:  x,y,z are unknown here.
cdr           Then VECUSR returns B field at COM (center of mass)
         CALL EIRENE_VECUSR (1,ICELL,X,Y,Z,BX,BY,BZ,1,L)
         bf = 1.

      ELSE IF (LBSMO) THEN
cdr    use logical input flag L:  return either value at COM or local value at x,y,z

         bx = eirene_Femint(bxincorner, icell, x, y, z, .false.)
         by = eirene_Femint(byincorner, icell, x, y, z, .true.)
         bz = eirene_Femint(bzincorner, icell, x, y, z, .true.)
         bf = eirene_Femint(bfincorner, icell, x, y, z, .true.)
         bni = 1._dp/sqrt(bx*bx + by*by + bz*bz)
         bx = bx * bni
         by = by * bni
         bz = bz * bni

      ELSE
cdr if no BFIELD input tallies: 
cdr use default B field: 1 [T] in z-direction
         BX = 0._DP
         BY = 0._DP
         BZ = 1._DP
         BF = 1._DP
         IF (LBXIN) BX=BXIN(ICELL)
         IF (LBYIN) BY=BYIN(ICELL)
         IF (LBZIN) BZ=BZIN(ICELL)
         IF (LBFIN) BF=BFIN(ICELL)

      END IF

      RETURN
      END SUBROUTINE EIRENE_BFIELD
