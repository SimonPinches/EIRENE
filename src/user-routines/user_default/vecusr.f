
      SUBROUTINE EIRENE_VECUSR (I,VEC_X,VEC_Y,VEC_Z,IPLS)

c  This call is at position x0,y0,z0, in cell NCELL   (module: comprt.f).
c  I=1:  return local B-field vector as vec_x,vec_y,vec_z
c  I=2:  return local plasma drift velocity vector of eirene background species "ipls" 
c               as vec_x,vec_y,vec_z
c  note: ipls is the true background species, 
c        not to be confused with the mapped species index mplsv(ipls)
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMPRT, only: NCELL, X0,Y0,Z0 
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I, IPLS
      REAL(DP), INTENT(OUT) :: VEC_X,VEC_Y,VEC_Z
c
      VEC_X=0.0
      VEC_Y=0.0
      VEC_Z=1.0
      RETURN
      END
