cdr  Aug. 2016:  generalized (ifexmx<0 enabled), two new parameters in list for fct. extrap

cdr  this function evaluates the standard single parameter 8th order polynomial
cdr  fits for cross section and rate coefficients, used in the 
cdr  HYDHEL  (Janev, Langer et al, Springer, 1987)
cdr  METHANE (Ehrhardt, Langer et al, PPPL report) 
cdr  databases. See references in online manual.
cdr  the same fit format is also used most of the time in the eirene-home
cdr  databases amjuel, h2vibr,

      function EIRENE_sngl_poly (cf, al, rcmin, rcmax, fpp, 
     .                                   ifexmn, ifexmx)
     .                   result(cou)
c  input:
c  cf    : fit coefficients for fit f(parm=)=sum_1^9 (cf(i) log(parm)^(i-1)) 
c  al    : argument of fit, log(parm)
c  rcmin : left boundary of valid range
c  rcmax : right boundary of valid range
c  fpp   : parameters for extrapolation from valid range
c  ifexmn: flag for choice of left (low end) extrapolation expression
c  ifexmx: flag for choice of right (high end) extrapolation expression 
 
      use EIRMOD_precision
 
      implicit none
 
      real(dp), intent(in) :: cf(9), fpp(6)
      real(dp), intent(in) :: al, rcmin, rcmax
      integer, intent(in) :: ifexmn, ifexmx
      real(dp) :: cou, fp(6), s01, s02, ds12, expo1, expo2, ccxm1,
     .            ccxm2, almin,almax,coumin,coumax, 
     .            EIRENE_extrap
      integer :: ii, if8, ifex
 
      if (al < rcmin) then
 
C  PARM BELOW MINIMUM PARAMETER FOR POLYNOM FIT:
 
        FP = FPP

C  USE ASYMPTOTIC EXPRESSION NO. IFEXMN
        IF (IFEXMN.LT.0) THEN
C  DETERMINE EXTRAPOLATION COEFFICIENTS FOR LINEAR EXTRAP. OF LOG(FIT) IN LN(SIGMA)
          S01=RCMIN
          S02=LOG(1.25_DP)+RCMIN    ! use parm=exp(rcmin) and 1.25*parm for extrapolation
          DS12=S02-S01
          EXPO1=CF(9)
          EXPO2=CF(9)
          DO 1 II=1,8
            IF8=9-II
            EXPO1=EXPO1*S01+CF(IF8)  !  evaluate fit at left boundary -->EXPO1=log(fit)
            EXPO2=EXPO2*S02+CF(IF8)  !  evaluate fit at PARM=(1.25*parleft), parleft= left boundary
 1        CONTINUE
          CCXM1=EXPO1
          CCXM2=EXPO2
          FP(1)=CCXM1+(CCXM2-CCXM1)/DS12*(-S01)
          FP(2)=      (CCXM2-CCXM1)/DS12
          FP(3)=0.D0
C
          IFEX=5
          ALMIN =RCMIN
          COUMIN=EXP(EXPO1)

        ELSEIF  (IFEXMN.GT.0) THEN
C  USE ASYMPTOTIC EXPRESSION NO. IFEXMN  
C  IFEXMN IS .GT. 0,  use preprogrammed extrapolation scheme no. ifexmn
          
          coumin = cf(9)
          do ii = 8, 1, -1
            coumin = coumin * RCMIN + cf(ii)
          end do
          ifex = ifexmn
C  determine parameter and fit value at left boundary. may be needed by fct. extrap
          ALMIN= RCMIN
          COUMIN=EXP(COUMIN)

        ELSE

C  AL IS OUT OF RANGE, BUT NO EXTRAPOLATION SCHEME SPECIFIED
C  WHAT TO WE DO NOW ???
          GOTO 100

        ENDIF
 
        COU=EIRENE_EXTRAP(AL,ALMIN,COUMIN,IFEX,FP(1),FP(2),FP(3))
        cou = log(cou)
        return
 
      elseif (al > rcmax) then
 
C  PARM IS ABOVE MAXIMUM VALID PARAMETER FOR FIT:
 
        FP = FPP
C  USE ASYMPTOTIC EXPRESSION NO. IFEXMX
        IF (IFEXMX.LT.0) THEN
C  DETERMINE EXTRAPOLATION COEFFICIENTS FOR LINEAR EXTRAP. OF LOG(FIT) IN LN(SIGMA)
          S01=RCMAX
          S02=LOG(0.75_DP)+RCMAX   ! use parm=exp(rcmin) and 0.75*parm for extrapolation 
          DS12=S02-S01
          EXPO1=CF(9)
          EXPO2=CF(9)
          DO 2 II=1,8
            IF8=9-II
            EXPO1=EXPO1*S01+CF(IF8)  !  evaluate fit at right boundary -->EXPO1=log(fit)
            EXPO2=EXPO2*S02+CF(IF8)  !  evaluate fit at PARM=(0.75*parleft), parleft= left boundary
2         CONTINUE
          CCXM1=EXPO1
          CCXM2=EXPO2
          FP(4)=CCXM1+(CCXM2-CCXM1)/DS12*(-S01)
          FP(5)=      (CCXM2-CCXM1)/DS12
          FP(6)=0.D0
C
          IFEX=5
          ALMAX =RCMAX
          COUMAX=EXP(EXPO1)

        ELSEIF (IFEXMX.GT.0) THEN 
C  USE ASYMPTOTIC EXPRESSION NO. IFEXMN  
C  IFEXMX IS .GT. 0,  use preprogrammed extrapolation scheme no. ifexmx   
          
          COUMAX = cf(9)
          do ii = 8, 1, -1
            coumax = coumax * RCMAX + cf(ii)
          end do
          ifex = ifexmx
C  determine parameter and fit value at right boundary. may be needed by fct. extrap
          ALMAX= RCMAX
          COUMAX=EXP(COUMAX)

        ELSE

C  AL IS OUT OF RANGE, BUT NO EXTRAPOLATION SCHEME SPECIFIED
C  WHAT TO WE DO NOW ???
          GOTO 100

        ENDIF

        COU=EIRENE_EXTRAP(AL,ALMAX,COUMAX,IFEX,FP(4),FP(5),FP(6))
        cou = log(cou)
        return
 
      ENDIF
 
C  PARAMETER "AL" IS WITHIN VALID RANGE OF FIT:
 
100   cou = cf(9)
 
      do ii = 8, 1, -1
        cou = cou * al + cf(ii)
      end do
 

      return
      end function EIRENE_sngl_poly
