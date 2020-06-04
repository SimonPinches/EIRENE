CDK USER
C
C   USER SUPPLIED SUBROUTINES
C
C
      SUBROUTINE EIRENE_PROUSR (PRO,INDX,P0,P1,P2,P3,P4,P5,PROVAC,N)

c  user defined backgound (plasma) profiles.
c  called from routine plasma.f, indpro=5 option
c  input:  indx  : indicates which input tally (e.g. Te, ni, vx, etc...)
c          provac: eirene default vacuum value (e.g. dvac, tvac, etc...)
c          N     : N=NSURF, standard grid without additional cell region
C 
C  EXAMPLE:
c  mache exp. dichteprofil in y richtung, und dichte in x richtung konstant.
c     tbd. : genauer:  indx anschliesen...

C     P0 : CENTRAL VALUE  (at psurf(1))
C     P1 : outer   value  (at psurf(np2nd)) , profile: exponential decay
c    indx not used, da eh nur ein aufruf: indpro=5 genau nur fuer dichteprofil
c                   im anwendungsfall hier.

C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: P0, P1, P2, P3, P4, P5, PROVAC
      REAL(DP), INTENT(OUT) :: PRO(*)
      INTEGER, INTENT(IN) :: INDX, N
      integer  :: i,j
      real(dp) :: E, densj
      

      PRO(1:N)=0.0

c   find decay length E in y direction
      E=-(psurf(np2nd)-psurf(1))/(log(p1)-log(p0))
           
c set density constant in x direction and exp profile in y direction (2nd grid)
      do i=1,nr1st-1
        do j=1,np2nd-1
          densj= p0*exp(-(phzone(j)-psurf(1))/e)
          pro(i+nr1st*(j-1))= densj
        enddo
      enddo      
 
      RETURN
      END
