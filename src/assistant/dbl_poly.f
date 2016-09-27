cdr: aug 16:  comments, bug fixes re.   ji,je=1,1 option
cdr           finally: remove ji:je parameters, set back to defaults: 1:9 
cdr: sept.16: prepare extrapolation options for double poly. fits.
cdr:          to be written: 2nd parameter out of range   (2 options)
cdr:          to be written: both parameters out of range (4 options)

c  called from: rate_coeff.f
c               energy_rate_coeff.f
c               sgnal.f

      subroutine EIRENE_dbl_poly (cf, al1, al2, cou, dum, 
     .                      rc1min, rc1max, fp1, ifex1mn, ifex1mx,
     .                      rc2min, rc2max, fp2, ifex2mn, ifex2mx)

c  used for evaluating double polynomial fit. 
c  return        COU=fit2(AL1,AL2)
C  INTERNAL INTERPRETATION: AL1=LOG(PARM1), 
C                            AL2=LOG(PARM2), WITH PHYSICAL PARAMETERS PARM1, PARM2
c  and return also:
c  return reduced fit coefficients dum(1,9) for AL2 dependency 
c  at fixed first parameter AL1. Then COU=fit1(AL2)=sum_1^9 [dum(j) log(parm2)^(j-1)]
 
c  input:
c  cf    : fit coefficients for 2D polyn. fit. We require data: cf(1:9,1:9) 
c          f(parm1,parm2)=sum_1^9 [sum_1^9 (cf(i,j) log(parm1)^(i-1)] log(parm2)^(j-1)
c  al1   : 1st argument of fit, log(parm1)  e.g. Te, Ti, ...
c  al2   : 2nd argument of fit, log(parm2)  e.g. ne, E0, ...

c  rc1min : left boundary of valid range for first parameter al1
c  rc1max : right boundary of valid range for first parameter al2
c  fpp    : parameters for extrapolation from valid range
c  ifex1mn: flag for choice of left (low end) extrapolation expression
c  ifex1mx: flag for choice of right (high end) extrapolation expression 
c  output:
c  cou   : log of value of fit(al1,al2)
c  dum   : dum(1,9), summation over 1st parameter is done, 
c          dum provides 1D fit coefficients for 2nd parameter dependency (e.g. ne, E0)
c          for reduced 1D fit evaluation "on the fly", at fixed parameter AL1.     
   
      use EIRMOD_precision
 
      implicit none
 
      real(dp), intent(in) :: cf(9,9), fp1(6), fp2(6)
      real(dp), intent(in) :: al1, al2, rc1min, rc1max, rc2min, rc2max
      real(dp), intent(out) :: cou, dum(9)
      integer, intent(in) :: ifex1mn, ifex1mx, ifex2mn, ifex2mx
      real(dp) :: s01, s02, ds12, expo1, expo2, ccxm1, ccxm2,
     .            fpar1, fpar2, fpar3,
     .            al1min,al1max,cou1min,cou1max, 
     .            eirene_extrap
      integer :: kk, jj, j, i, ii, ifex
 
      dum = 0._dp
 
      if (al1 < rc1min) then
        IF (IFEX1MN.LT.0) THEN
C  DETERMINE EXTRAPOLATION COEFFICIENTS FOR AL1 - LINEAR EXTRAP. IN LN(FIT) AT FIXED AL2
C  EVALUATED AT AL2, WHICH MUST BE INSIDE VALID RANGE
          S01=RC1MIN
          S02=LOG(1.25_DP)+RC1MIN
          DS12=S02-S01
          EXPO1=0.
          EXPO2=0.
c  evaluate double parameter fit at s01,al2,  and at s02,al2
          DO 1 J=1,9
            JJ=J-1
            DO 1 I=1,9
              II=I-1
              EXPO1=EXPO1+S01**II*AL2**JJ*CF(I,J)
              EXPO2=EXPO2+S02**II*AL2**JJ*CF(I,J)
 1        CONTINUE
          CCXM1=EXPO1
          CCXM2=EXPO2
          FPAR1=CCXM1+(CCXM2-CCXM1)/DS12*(-S01)
          FPAR2=      (CCXM2-CCXM1)/DS12
          FPAR3=0.D0
C
          IFEX=5
          AL1MIN=RC1MIN
          COU1MIN=EXP(EXPO1)
          
        ELSEIF (IFEX1MN.GT.0) THEN
C  IFEX1MN IS GT 0, USE ONE OF THE PREPROGRAMMED EXTRAPOLATION SCHEMES 
          IFEX=IFEX1MN

        ELSE
C  FIT OUT OF VALID RANGE, AL1 < RC1MIN, BUT NO EXTRAPOLATION SCHEME AVAILABLE
          GOTO 100
        ENDIF

        COU=EIRENE_EXTRAP(AL1,AL1MIN,COU1MIN,IFEX,FPAR1,FPAR2,FPAR3)
        cou=log(cou)
 
      elseif (al1 > rc1max) then
        IF (IFEX1MX.LT.0) THEN     
C  DETERMINE EXTRAPOLATION COEFFICIENTS FOR AL1 - LINEAR EXTRAP. IN LN(FIT) AT FIXED AL2
C  EVALUATED AT AL2, WHICH MUST BE INSIDE VALID RANGE
          S01=RC1MAX
          S02=log(0.75_dp)+RC1MAX
          DS12=S02-S01
          EXPO1=0.
          EXPO2=0.
c  evaluate double parameter fit at s01,al2,  and at s02,al2
          DO J=1,9
            JJ=J-1
            DO I=1,9
              II=I-1
              EXPO1=EXPO1+S01**II*AL2**JJ*CF(I,J)
              EXPO2=EXPO2+S02**II*AL2**JJ*CF(I,J)
            enddo
          enddo
          CCXM1=EXPO1
          CCXM2=EXPO2
          FPAR1=CCXM1+(CCXM2-CCXM1)/DS12*(-S01)
          FPAR2=      (CCXM2-CCXM1)/DS12
          FPAR3=0.D0
C
          IFEX=5
          AL1MAX=RC1MAX
          COU1MAX=EXP(EXPO1)

        ELSEIF  (IFEX1MN.GT.0) THEN 
C  IFEX1MX IS GT 0, USE ONE OF THE PREPROGRAMMED EXTRAPOLATION SCHEMES
          IFEX=IFEX1MX

        ELSE
C  FIT OUT OF VALID RANGE, AL1 > RC1MAX, BUT NO EXTRAPOLATION SCHEME AVAILABLE
          GOTO 100
        ENDIF

        COU=EIRENE_EXTRAP(AL1,AL1MAX,COU1MAX,IFEX,FPAR1,FPAR2,FPAR3)
        cou=log(cou)

      elseif (al2 < rc2min) then

cdr   to be written
        goto 100

      elseif (al2 > rc2max) then

cdr   to be written
        goto 100        

      else

cdr  both parameters al1 and al2 are within valid range of fit
cdr  evaluate double polynomial fit.  
cdr  result = cou
cdr  the reduced single parameter fit coefficients at fixed AL1 are returned on dum(1:9)
        
        goto 100

      endif


100   CONTINUE
 
      do jj = 9, 1, -1   
        dum(jj) = cf(9,jj)
        do kk = 8, 1, -1
          dum(jj) = dum(jj) * al1 + cf(kk,jj)
        end do
      end do
 
      cou = dum(9)          
      do jj = 8, 1, -1      
        cou = cou * al2 + dum(jj)
      end do
 
   
 
      return
      end subroutine EIRENE_dbl_poly
