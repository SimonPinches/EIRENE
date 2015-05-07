
      SUBROUTINE EIRENE_VECUSR (I,VEC_X,VEC_Y,VEC_Z,IPLS)
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: I, IPLS
      REAL(DP), INTENT(OUT) :: VEC_X,VEC_Y,VEC_Z
c  I=1:  return local B-field vector as vec_x,vec_y,vec_z
c  I=2:  return local plasma drift velocity vector of eirene background species "ipls" 
c               as vec_x,vec_y,vec_z
c  note: ipls is the true background species, 
c        not to be confused with the mapped species index mplsv(ipls)
      RETURN
      END 
 

