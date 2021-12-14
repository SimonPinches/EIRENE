      SUBROUTINE EIRENE_SAMVOL_USR(n,x,y,z)
cdr  added Nov. 2018:
CDR  for user geometry option LEVGEO= 10, 
cdr  provide a particle birth point x,y,z in cell N,
cdr  to be used in samvol.f (volume sources).
cdr  corresponding call addded by HF, to samvol, in case levgeo=10
cdr  
      USE EIRMOD_PRECISION

      IMPLICIT NONE
    
      INTEGER, INTENT(IN) :: N
      REAL(DP), INTENT(OUT) :: X,Y,Z
      x=0._DP
      y=0._DP
      z=0._DP
      RETURN
      END SUBROUTINE EIRENE_SAMVOL_USR
