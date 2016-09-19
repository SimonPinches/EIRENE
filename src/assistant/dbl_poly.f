cdr: aug 16:  comments, bug fixes re.   ji,je=1,1 option
cdr           finally: remove ji:je parameters, set back to defaults: 1:9 
cdr: sept.16: prepare extrapolation options for double poly. fits.
c  called from: rate_coeff.f
c               energy_rate_coeff.f
c               sgnal.f

      subroutine EIRENE_dbl_poly (cf, al1, al2, cou, dum, 
     .                      rcmin, rcmax, fp, ifexmn, ifexmx )

c  used for evaluating double polynomial fit. 
c  return        COU=fit2(AL1,AL2)
C  INTERNAL INTERPRETATION: AL1=LOG(PARM1), 
C                            AL2=LOG(PARM2), WITH PHYSICAL PARAMETERS PARM1, PARM2
c  and return also:
c  return reduced fit coefficients dum(1,9) for AL2 dependency 
c  at fixed first parameter AL1. Then COU=fit1(AL2)=sum_1^9 [dum(j) log(parm2)^(j-1)] 
 
c  input:
c  cf    : fit coefficients for 2D polyn. fit. We require data: cf(1:9,1:9) 
c          f(parm=)=sum_1^9 [sum_1^9 (cf(i,j) log(parm1)^(i-1)] log(parm2)^(j-1)
c  al1   : 1st argument of fit, log(parm1)  e.g. Te, Ti, ...
c  al2   : 2nd argument of fit, log(parm2)  e.g. ne, E0, ...

c  rcmin : left boundary of valid range for first parameter al1
c  rcmax : right boundary of valid range for first parameter al2
c  fpp   : parameters for extrapolation from valid range
c  ifexmn: flag for choice of left (low end) extrapolation expression
c  ifexmx: flag for choice of right (high end) extrapolation expression 
c  output:
c  cou   : log of value of fit(al1,al2)
c  dum   : dum(1,9) ,summation over 1st parameter is done, 
c          dum provides 1D fit coefficients for 2nd parameter dependency (e.g. ne, E0)
c          for reduced 1D fit evaluation "on the fly", at fixed parameter al1.     
 

      use EIRMOD_precision
 
      implicit none
 
      real(dp), intent(in) :: cf(9,9), fp(6)
      real(dp), intent(in) :: al1, al2, rcmin, rcmax
      real(dp), intent(out) :: cou, dum(9)
      integer, intent(in) :: ifexmn, ifexmx
      real(dp) :: s01, s02, ds12, expo1, expo2, ccxm1, ccxm2,
     .            fpar1, fpar2, fpar3,
     .            almin,almax,coumin,coumax, 
     .            eirene_extrap
      integer :: kk, jj, j, i, ii, ifex
 
      dum = 0._dp
 
      if ((ifexmn .ne. 0) .and. (al1 < rcmin)) then
        IF (IFEXMN.LT.0) THEN
C  DETERMINE EXTRAPOLATION COEFFICIENTS FOR AL1 - LINEAR EXTRAP. IN LN(FIT) AT FIXED AL2
C  EVALUATED AT AL2, WHICH MUST BE INSIDE VALID RANGE
          S01=RCMIN
          S02=LOG(1.25_DP)+RCMIN
          DS12=S02-S01
          EXPO1=0.
          EXPO2=0.
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
          ALMIN=RCMIN
          COUMIN=EXP(EXPO1)
          
        ELSE  !  IFEXMN IS GT 0
          IFEX=IFEXMN
        ENDIF

        COU=EIRENE_EXTRAP(AL1,ALMIN,COUMIN,IFEX,FPAR1,FPAR2,FPAR3)
        cou=log(cou)
 
      elseif ((ifexmx .ne. 0) .and. (al1 > rcmax)) then
C  DETERMINE EXTRAPOLATION COEFFICIENTS FOR AL1 - LINEAR EXTRAP. IN LN(FIT) AT FIXED AL2
C  EVALUATED AT AL2, WHICH MUST BE INSIDE VALID RANGE
        IF (IFEXMX.LT.0) THEN     
          S01=RCMAX
          S02=log(0.75_dp)+RCMAX
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
          ALMAX=RCMAX
          COUMAX=EXP(EXPO1)

        ELSE  !  IFEXMX IS GT 0
          IFEX=IFEXMX
        ENDIF

        COU=EIRENE_EXTRAP(AL1,ALMAX,COUMAX,IFEX,FPAR1,FPAR2,FPAR3)
        cou=log(cou)

      else

cdr  both parameters al1 and al2 are within valid range of fit
cdr  evaluate double polynomial fit.  
cdr  result = cou
cdr  the reduced single parameter fit coefficients at fixed AL1 are returned on dum(1:9)
 
        do jj = 9, 1, -1   ! 
          dum(jj) = cf(9,jj)
          do kk = 8, 1, -1
            dum(jj) = dum(jj) * al1 + cf(kk,jj)
          end do
        end do
 
        cou = dum(9)          
        do jj = 8, 1, -1      
          cou = cou * al2 + dum(jj)
        end do
 
      end if
 
      return
      end subroutine EIRENE_dbl_poly
