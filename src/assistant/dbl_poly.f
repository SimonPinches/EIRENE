cdr: aug 16:  comments, bug fixes re.   ji,je=1,1 option
cdr           finally: remove ji:je parameters, set back to defaults: 1:9 
c  called from: rate_coeff.f
c               energy_rate_coeff.f
c               prep_rtcs.f
c               sgnal.f

      subroutine EIRENE_dbl_poly (cf, al, pl, cou, dum, 
     .                      rcmin, rcmax, fp, ifexmn, ifexmx )

c  used for evaluating double polynomial fit. 
c  return        COU=fit2(AL,PL)
c  and also:
c  return reduced fit coefficients dum(1,9) for PL dependency 
c  at fixed first parameter AL. Then COU=fit1(PL)=sum_1^9 [dum(j) log(parm2)^(j-1)] 
 
c  input:
c  cf    : fit coefficients for 2D polyn. fit. We require data: cf(1:9,1:9) 
c          f(parm=)=sum_1^9 [sum_1^9 (cf(i,j) log(parm1)^(i-1)] log(parm2)^(j-1)
c  al    : 1st argument of fit, log(parm1)  e.g. Te, Ti, ...
c  pl    : 2nd argument of fit, log(parm2)  e.g. ne, E0, ...

c  rcmin : left boundary of valid range for first parameter al
c  rcmax : right boundary of valid range for first parameter al
c  fpp   : parameters for extrapolation from valid range
c  ifexmn: flag for choice of left (low end) extrapolation expression
c  ifexmx: flag for choice of right (high end) extrapolation expression 
c  output:
c  cou   : log of value of fit(al,pl)
c  dum   : dum(1,9) ,summation over 1nd parameter is done, 
c          dum provides 1D fit coefficients for 2nd parameter dependency (e.g. E0)
c          for reduced 1D fit evaluation "on the fly", at fixed parameter al.     
 

      use EIRMOD_precision
 
      implicit none
 
      real(dp), intent(in) :: cf(9,9), fp(6)
      real(dp), intent(in) :: al, pl, rcmin, rcmax
      real(dp), intent(out) :: cou, dum(9)
      integer, intent(in) :: ifexmn, ifexmx
      real(dp) :: s01, s02, ds12, expo1, expo2, ccxm1, ccxm2,
     .            fpar1, fpar2, fpar3,
     .            almin,almax,coumin,coumax, 
     .            eirene_extrap
      integer :: kk, jj, j, i, ii, ifex
 
      dum = 0._dp
 
      if ((ifexmn .ne. 0) .and. (al < rcmin)) then
        IF (IFEXMN.LT.0) THEN
C  DETERMINE EXTRAPOLATION COEFFICIENTS FOR AL - LINEAR EXTRAP. IN LN(FIT)
          S01=RCMIN
          S02=LOG(1.25_DP)+RCMIN
          DS12=S02-S01
          EXPO1=0.
          EXPO2=0.
          DO 1 J=1,9
            JJ=J-1
            DO 1 I=1,9
              II=I-1
              EXPO1=EXPO1+S01**II*PL**JJ*CF(I,J)
              EXPO2=EXPO2+S02**II*PL**JJ*CF(I,J)
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

        COU=EIRENE_EXTRAP(AL,ALMIN,COUMIN,IFEX,FPAR1,FPAR2,FPAR3)
        cou=log(cou)
 
      elseif ((ifexmx .ne. 0) .and. (al > rcmax)) then
C  DETERMINE EXTRAPOLATION COEFFICIENTS FOR AL - LINEAR EXTRAP. IN LN(FIT)
C  EVALUATED AT PL, WHICH MUST BE INSIDE VALID RANGE
        IF (IFEXMX.LT.0) THEN     
          S01=RCMAX
          S02=log(0.75_dp)+RCMAX
          DS12=S02-S01
          EXPO1=0.
          EXPO2=0.
c  evaluate double parameter fit at s01,pl,  and at s02,pl
          DO J=1,9
            JJ=J-1
            DO I=1,9
              II=I-1
              EXPO1=EXPO1+S01**II*PL**JJ*CF(I,J)
              EXPO2=EXPO2+S02**II*PL**JJ*CF(I,J)
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

        COU=EIRENE_EXTRAP(AL,ALMAX,COUMAX,IFEX,FPAR1,FPAR2,FPAR3)
        cou=log(cou)

      else

cdr  both parameters al and pl are within valid range of fit
cdr  evaluate double polynomial fit.  result = cou
cdr  the reduced single parameter fit coefficients is on dum(1:9)
 
        do jj = 9, 1, -1   ! 
          dum(jj) = cf(9,jj)
          do kk = 8, 1, -1
            dum(jj) = dum(jj) * al + cf(kk,jj)
          end do
        end do
 
        cou = dum(9)          
        do jj = 8, 1, -1      
          cou = cou * pl + dum(jj)
        end do
 
      end if
 
      return
      end subroutine EIRENE_dbl_poly
