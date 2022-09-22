cdr  Nov. 2019: add functionality for Arrhenius factors EXP(-DE/T)
cdr             excluded from rest of fit.
cdr             Should replace the need for
cdr             low temp asymptotics in H.2, H.5, and H.8 polynomial data.
cdr  Aug. 2016:  generalized (ifexmx<0 enabled), two new parameters in list for fct. extrap

cdr  This function evaluates the standard single parameter 8th-order polynomial
cdr  fits for cross-section and rate coefficients, used in the
cdr  HYDHEL  (Janev, Langer et al, Springer, 1987)
cdr  METHANE (Ehrhardt, Langer et al, PPPL report)
cxb  AMMONX  (Mougenot, Touchard et al., LSPM-CNRS, 2017)
cdr  databases. See references in online manual.
cdr  The same fit format is also used most of the time in the eirene home-made
cdr  databases AMJUEL, H2VIBR.

      function EIRENE_sngl_poly (cf, al, rcmin, rcmax, fpp,
     .                                   ifexmn, ifexmx, 
     .                                   earrh0, trc, lexp)
     .                   result(cou)
c  input:
c  cf      : fit coefficients for fit POLY=f(parm=)=sum_1^9 (cf(i) log(parm)^(i-1))
c  al      : argument of fit, al=log(parm)
c  rcmin   : left boundary of valid range of PARM
c  rcmax   : right boundary of valid range of PARM
c  fpp(1:6): parameters for extrapolation from valid range
c            1:3  left (low PARM values) extrapolation.
c            4:6  right (high PARM values) extrapolation
c  ifexmn: flag for choice of left (low end) extrapolation expression
c  ifexmx: flag for choice of right (high end) extrapolation expression
c  earrh0: Arrhenius factor exp(-earrh0/T) separated from fit
cdr           ifex=0:  constant extrapolation
cdr           ifex<0:  determine extrapolation parameters here
cdr                    (linear extrapolation),
cdr                    and call extrap.f with model IFEX=3
cdr           ifex>0:  evaluate fit at corresponding boundary,
cdr                    and call extrap.f with model IFEX
c  trc     : flag for diagnostic print output
c  lexp    : return exp(POLY), i.e. the fit polynomial POLY is the
c                                   Logarithm of the requested result
c            a cutoff at COU == -100 is introduced to avoid floating
c            point underflow exceptions 
c output:
c  lexp=false:  return POLY
c  lexp=true:   return exp(POLY)
      
      use EIRMOD_precision
      USE EIRMOD_COMPRT, ONLY: IUNOUT

      implicit none

      real(dp), intent(in) :: cf(9), fpp(6)
      real(dp), intent(in) :: al, rcmin, rcmax, earrh0
      integer, intent(in) :: ifexmn, ifexmx
      logical, intent(in) :: trc, lexp
      real(dp) :: p1, ep1, cou, fp(6), s01, s02, ds12, expo1, expo2,
     .            ccxm1, ccxm2, almin, almax, coumin, coumax,
     .            EIRENE_extrap
      integer :: ii, if8, ifex

      p1=al        ! fit parameterfor log-log fit: log(T), log(E),...
c     ep1=exp(p1)  ! physical parameter T, E

      if (p1 < rcmin) then

C  PARM BELOW MINIMUM PARAMETER FOR POLYNOMIAL FIT:


C  USE ASYMPTOTIC EXPRESSION NO. IFEXMN
        IF (IFEXMN.LT.0) THEN
C  DETERMINE EXTRAPOLATION COEFFICIENTS FOR LINEAR EXTRAP. OF LOG(FIT) IN LN(parm)
          S01=RCMIN
          S02=LOG(1.25_DP)+RCMIN    ! use parm=exp(rcmin) and 1.25*parm for extrapolation
          DS12=S02-S01
          EXPO1=CF(9)
          EXPO2=CF(9)
          DO 1 II=1,8
            IF8=9-II
            EXPO1=EXPO1*S01+CF(IF8)  !  evaluate fit at left boundary -->EXPO1=log(fit)
            EXPO2=EXPO2*S02+CF(IF8)  !  evaluate fit at PARM=(1.25*parleft), parleft= left boundary
    1     CONTINUE
          CCXM1=EXPO1
          CCXM2=EXPO2
          FP(1)=CCXM1+(CCXM2-CCXM1)/DS12*(-S01)
          FP(2)=      (CCXM2-CCXM1)/DS12
          FP(3)=0.D0
C  Linear extrapolation on log-log scale
          IFEX=3
          ALMIN =RCMIN
          COUMIN=EXP(EXPO1)

        ELSEIF (IFEXMN.GT.0) THEN
C  USE ASYMPTOTIC EXPRESSION NO. IFEXMN
C  IFEXMN IS .GT. 0, use pre-programmed extrapolation scheme no. ifexmn

          coumin = cf(9)
          do ii = 8, 1, -1
            coumin = coumin * RCMIN + cf(ii)
          end do
          ifex = ifexmn
C  determine parameter and fit value at left boundary. May be needed by fct. extrap
          ALMIN= RCMIN
          COUMIN=EXP(COUMIN)
          FP(1:3) = FPP(1:3)

        ELSE

C  AL IS OUT OF RANGE, BUT NO EXTRAPOLATION SCHEME SPECIFIED (IFEXMN=0)
C  Continue with a constant: poly evalutated at RCMIN
          P1=RCMIN
          if (trc) write (iunout,*) 'unclear extrapolation in sngl_poly'
          GOTO 100

        ENDIF

        COU=EIRENE_EXTRAP(P1,ALMIN,COUMIN,IFEX,FP(1),FP(2),FP(3))
        cou = log(cou)
        goto 1000

      elseif (p1 > rcmax) then

C  PARM IS ABOVE MAXIMUM VALID PARAMETER FOR FIT:

C  USE ASYMPTOTIC EXPRESSION NO. IFEXMX
        IF (IFEXMX.LT.0) THEN
C  DETERMINE EXTRAPOLATION COEFFICIENTS FOR LINEAR EXTRAP. OF LOG(FIT) IN LN(parm)
          S01=RCMAX
          S02=LOG(0.75_DP)+RCMAX   ! use parm=exp(rcmin) and 0.75*parm for extrapolation
          DS12=S02-S01
          EXPO1=CF(9)
          EXPO2=CF(9)
          DO 2 II=1,8
            IF8=9-II
            EXPO1=EXPO1*S01+CF(IF8)  !  evaluate fit at right boundary -->EXPO1=log(fit)
            EXPO2=EXPO2*S02+CF(IF8)  !  evaluate fit at PARM=(0.75*parright), parright= right boundary
    2     CONTINUE
          CCXM1=EXPO1
          CCXM2=EXPO2
          FP(4)=CCXM1+(CCXM2-CCXM1)/DS12*(-S01)
          FP(5)=      (CCXM2-CCXM1)/DS12
          FP(6)=0.D0
C  Linear extrapolation on log-log scale
          IFEX=3
          ALMAX =RCMAX
          COUMAX=EXP(EXPO1)

        ELSEIF (IFEXMX.GT.0) THEN
C  USE ASYMPTOTIC EXPRESSION NO. IFEXMN
C  IFEXMX IS .GT. 0, use pre-programmed extrapolation scheme no. ifexmx

          COUMAX = cf(9)
          do ii = 8, 1, -1
            coumax = coumax * RCMAX + cf(ii)
          end do
          ifex = ifexmx
C  determine parameter and fit value at right boundary. may be needed by fct. extrap
          ALMAX= RCMAX
          COUMAX=EXP(COUMAX)
          FP(4:6) = FPP(4:6)
          
        ELSE

C  AL IS OUT OF RANGE, BUT NO EXTRAPOLATION SCHEME SPECIFIED (IFEXMX=0)
C  Continue with a constant: poly evalutated at RCMAX
          P1=RCMAX
          if (trc) write (iunout,*) 'unclear extrapolation in sngl_poly'
          GOTO 100

        ENDIF

        COU=EIRENE_EXTRAP(P1,ALMAX,COUMAX,IFEX,FP(4),FP(5),FP(6))
        cou = log(cou)
        goto 1000

      ENDIF

C  PARAMETER "P1=AL" IS WITHIN VALID RANGE OF FIT:

  100 cou = cf(9)

      do ii = 8, 1, -1
        cou = cou * p1 + cf(ii)
      end do

c  Arrhenius factor exp(-earrh0/T), here: add log thereof to the fit POLY.
c  ie. add -earrh0/parm
      if (earrh0.gt.0.0) then
        ep1=exp(p1)
        cou=cou-earrh0/ep1
      endif

 1000 continue
cdr  return cou=POLY, or cou=exp(POLY):
      if (lexp) cou = exp(max(-100._dp, cou))
      
      return
      end function EIRENE_sngl_poly
