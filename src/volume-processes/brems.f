      FUNCTION EIRENE_BREMS(Te,ne,zi)
     .         result(res)

cdr   bremsstrahlung, for ions of charge Zi,
cdr   in a bath of electrons at Te, ne.

c  input:
c  Te, electron temperature, eV
c  Ti, ion temperature, eV  (not used)
c  Zi, charge number of ions.
c  ne, electron density,  #/cm^3

      USE EIRMOD_PRECISION
      implicit none

      REAL(DP), INTENT(IN) :: TE,NE,ZI
      REAL(DP) :: EIRENE_NGFFMH_B, brems, RES


c  bremsstrahlung in W, per ion, based on free-free Gaunt factors
c  formula from ADAS, see function ngffmh_b
      if (zi.ne.0._DP) then
        BREMS = 1.54E-32_DP * TE**0.5 * ZI**2 *
     .         eirene_ngffmh_B(ZI**2 * 13.6_DP/TE) *ne
      else
        BREMS = 0._DP
      endif

      RES=BREMS

      return
      end function eirene_brems


c  free-free Gaunt factor routine, obtained from Martin O'Mullane in 2007
c  slightly adapted to use eirene precision convention (eirmod_precision)

       FUNCTION EIRENE_NGFFMH_B(GAM2)
       USE EIRMOD_PRECISION
       IMPLICIT NONE
C-----------------------------------------------------------------------
C
C  ************** FORTRAN77 SUBROUTINE: NGFFMH  ************************
C
C  VERSION:  1.0
C
C  PURPOSE:
C
C  EVALUATES ELECTRON TEMPERATURE- AND FREQUENCY-AVERAGED HYDROGENIC
C  FREE-FREE GAUNT FACTOR.
C  OBTAINED FROM INTERPOLATION OF KARZAS & LATTER (1959) FIG.6
C  FOR -3<LOG10(Z0*Z0*IH/KTE)<1. OUTSIDE THIS RANGE A VERY APPROXIMATE
C  EXTRAPOLATION IS PERFORMED WITH GFFMH=1 IN THE INFINITE LIMITS.
C
C  INPUT:
C       GAM2=Z0*Z0*IH/KTE
C  OUTPUT:
C       NGFFMH=MAXWELL- AND FREQUENCY-AVERAGED FREE-FREE GAUNT FACTOR.
C  ********* H.P.SUMMERS, JET         12 JAN 1987    ******************
C-----------------------------------------------------------------------
C
C UNIX-IDL CONVERSION:
C
C VERSION: 1.1                          DATE: 22-08-96
C MODIFIED: WILLIAM OSBORN
C               - FIRST CONVERTED. NO CHANGES.
C
C VERSION: 1.2                          DATE: 13-09-99
C MODIFIED: Martin O'Mullane
C               - Define function name as real*8.
C
C-----------------------------------------------------------------------
       REAL(DP), INTENT(IN) :: GAM2
       REAL(DP) :: EIRENE_NGFFMH_B, GAM2L
       INTEGER :: K

!PB    DIMENSION GAM2LA(17),GA(17)
       REAL(DP) :: GAM2LA(17),GA(17)
       DATA GAM2LA/-3.0D0,-2.75D0,-2.50D0,-2.25D0,-2.00D0,-1.75D0,
     &-1.50D0,-1.25D0,-1.00D0,-0.75D0,-0.50D0,-0.25D0,0.00D0,0.25D0,
     &0.50D0,0.75D0,1.00D0/
       DATA GA/1.139D0,1.151D0,1.167D0,1.189D0,1.215D0,1.248D0,1.283D0,
     &1.326D0,1.370D0,1.411D0,1.431D0,1.436D0,1.433D0,1.415D0,1.379D0,
     &1.338D0,1.296D0/

       GAM2L=DLOG10(GAM2)
       IF(GAM2L.LE.-3.0D0)GO TO 30
       IF(GAM2L.GE.1.0D0)GO TO 40
C  INTERPOLATION REGION, LOCATE GAM2L IN GAM2LA ARRAY
       K=0
   20  K=K+1
       IF(GAM2L.GT.GAM2LA(K)) GO TO 20
       K=K-1
       IF (K.EQ.16) K=15
C  QUADRATIC INTERPOLATION
       EIRENE_NGFFMH_B=(GAM2L-GAM2LA(K+1))*(GAM2L-GAM2LA(K+2))/
     &((GAM2LA(K)-GAM2LA(K+1))*(GAM2LA(K)-GAM2LA(K+2)))*GA(K)+
     &(GAM2L-GAM2LA(K))*(GAM2L-GAM2LA(K+2))/
     &((GAM2LA(K+1)-GAM2LA(K))*(GAM2LA(K+1)-GAM2LA(K+2)))*GA(K+1)+
     &(GAM2L-GAM2LA(K))*(GAM2L-GAM2LA(K+1))/
     &((GAM2LA(K+2)-GAM2LA(K))*(GAM2LA(K+2)-GAM2LA(K+1)))*GA(K+2)
       RETURN
C  EXTRAPOLATION FOR LOW GAM2
   30  EIRENE_NGFFMH_B=1.0D0-0.417D0/GAM2L
       RETURN
   40  EIRENE_NGFFMH_B=1.0D0+0.296D0/GAM2L
       RETURN
       END FUNCTION EIRENE_NGFFMH_B
