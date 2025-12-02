cdr  jan. 2020:  internal consistency enforced for vector components,
cdr              and their smoothing/interpolation.
cdr              New: LBIN  flag:  all B field tallies exist.
cdr              NEW  LBSMO flag:  all B field tallies smoothed (interpolated)
cdr  aug. 2015:  logical flag L added. position x,y,z known (L=true)
cdr                           or else: use COM of cell ICELL
cdr              This flag is only used in case indpro(5)=8, i.e.
cdr              if B field is provided by user routine VECUSR
cdr  sept 2014:  comments added
c    provide cartesian local magnetic field unit vector bx,by,bz,
c    as well as B field strength bf (Tesla)
c    at point x,y,z, in cell icell

cdr Jan 22: ...and also the perp and diamagn components:
cdr         bx,by,bz (=parallel B)       --> bx(0),by(0),bz(0)
cdr         perp.    (=grad Psi)         --> bx(1),by(1),bz(1)
cdr         diamag.  (=B cross grad Psi) --> bx(2),by(2),bz(2)

c  current options:
c  default    :  use input background tallies. B field is constant per cell
c  indpro(5)=8:  user-provided B field
c  LBSMO (?)  :  apparently: only in case of levgeo=4,5 interpolation in triangles, tetrahedra
cpb              switch used for interpolation of magnetic field input tally
cdr              from cell vertices to a local x,y,z point inside a cell.
cdr              Not fully available for all levgeo=1,2,3 optins. Check FEMINT.f

      subroutine eirene_bfield (icell, x, y, z, bx, by, bz, bf, l)

      use eirmod_precision, only: dp
      use eirmod_comusr, only: BXIN, BYIN, BZIN, BFIN, BXINCORNER,
     >                         BYINCORNER, BZINCORNER, BFINCORNER,
     >                         bxperp, byperp,
     >                         LBSMO, LBIN, lbxperp, lbyperp
      use eirmod_cinit, only: INDPRO

      implicit none

      integer, intent(in) :: icell
      logical             :: l, lsame
      real(dp), intent(in) :: x, y, z
      real(dp), intent(out) :: bx(0:2), by(0:2), bz(0:2), bf
      real(dp) :: eirene_femint, bni
      external :: eirene_vecusr, eirene_femint

cdr  default "parallel" BX(0),BY(0),BZ(0) and BF are set below.
      BX(1:2) = 0.0
      BY(1:2) = 0.0
      BZ(1:2) = 0.0

      IF (INDPRO(5) == 8) THEN
cdr  user defined B field. Units of Bx, By, Bz?
cdr  BF=1. ?  BF should be in Tesla.
cdr  L=true : spatial coordinates x,y,z are known here,
cdr           i.e. call vecusr with L=true
cdr  L=false: position x,y,z is unknown here.
cdr           Then VECUSR returns B field at COM (center of mass) in cell ICELL
        CALL EIRENE_VECUSR (1,ICELL,X,Y,Z,BX(0),BY(0),BZ(0),1,L)

cdr  unfinished coding here? BF in Tesla?
        bf = 1.

      ELSE IF (LBSMO) THEN
cdr
        lsame=.false.
        bx(0) = eirene_Femint(bxincorner, icell, x, y, z, lsame)
        lsame=.true.   ! next calls to FEMINT at same position x,y,z
        by(0) = eirene_Femint(byincorner, icell, x, y, z, lsame)
        bz(0) = eirene_Femint(bzincorner, icell, x, y, z, lsame)
cdr  re-normalize unit vector, after interpolation.
        bni = 1._dp/sqrt(bx(0)*bx(0) + by(0)*by(0) + bz(0)*bz(0))
        bx(0) = bx(0) * bni
        by(0) = by(0) * bni
        bz(0) = bz(0) * bni
        bf = eirene_Femint(bfincorner, icell, x, y, z, lsame)

      ELSEIF (LBIN) THEN
cdr cartesian unit vector
        BX(0)=BXIN(ICELL)
        BY(0)=BYIN(ICELL)
        BZ(0)=BZIN(ICELL)
cdr  perpendicular (grad Psi)
cdr  so far: for a 2d (tor. symmetric) B field only
        if (lbxperp) bx(1)=bxperp(icell)
        if (lbyperp) by(1)=byperp(icell)
                     bz(1)=0.
cdr  diamagnetic (bi-normal)
        bx(2) = by(0)*bz(1) - bz(0)*by(1)
        by(2) =-bz(0)*bx(1) + bx(0)*bz(1)
        bz(2) = bx(0)*by(1) - by(0)*bx(1)

cdr modulus of B field, T
        BF=BFIN(ICELL)

      ELSE
cdr Default: B field input tallies BXIN,BYIN,BZIN,BFIN
cdr IF no BFIELD input tallies: use default B field: 1 [T] in z-direction
cdr
        BX(0) = 0._DP
        BY(0) = 0._DP
        BZ(0) = 1._DP
cdr  modulus of B field, T
        BF = 1._DP

      END IF

      RETURN
      END SUBROUTINE EIRENE_BFIELD
