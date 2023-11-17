cdr: aug 16:  comments, bug fixes re. ji,je=1,1 option
cdr           finally: remove ji:je parameters, set back to defaults: 1:9
cdr: sept.16: prepare extrapolation options for double poly. fits.
cdr:          to be written: 2nd parameter out of range   (2 options)
cdr:          to be written: both parameters out of range (4 options)
cdr: may 17 : speedup possible, if only dum(..) but not cou is needed. Tbd.
cdr: july 17: trc: print warning in case of extrapolation
cdr: aprl 18: extrapolation wrt. 2nd parameter added.
cdr           ifex=0: constant extrapolation
cdr           ifex<0: find extrapolation parameters here, and call extrap.f
cdr           ifex>0: find parameters boundary, and call extrap.f
cdr  aug. 20: code safeties from ITER branch: activate transfer of fp1, fp2
cdr           to routine extrap.f
cdr           Still missing: what if both al1 AND al2 are out of range?
cdr           Still missing: Arrhenius factors (if specified)

c  called from: rate_coeff.f
c               energy_rate_coeff.f
c               other_rate_coeff.f

      subroutine EIRENE_dbl_poly (cf, al1, al2, cou, dum,
     .                      rc1min, rc1max, fp1, ifex1mn, ifex1mx,
     .                      rc2min, rc2max, fp2, ifex2mn, ifex2mx,
     .                      trc)

c  used for evaluating double polynomial fit.
c  return        COU=fit2(AL1,AL2)
C  INTERNAL INTERPRETATION: AL1=LOG(PARM1),
C                           AL2=LOG(PARM2), WITH PHYSICAL PARAMETERS PARM1, PARM2
c  and return also:
c  reduced fit coefficients dum(1,9) for AL2 dependency at fixed first parameter AL1.
c  Then COU=fit1(AL2)=sum_1^9 [dum(j) log(parm2)^(j-1)]

c  input:
c  nd1, nd2: orders of double fit polynomial: nd1-1,nd2-1 (here still: nd1=nd2=9, merging tbd.)
c  cf    : fit coefficients for 2D polyn. fit. We require data: cf(1:9,1:9)
c          f(parm1,parm2)=sum_1^9 [sum_1^9 (cf(i,j) log(parm1)^(i-1)] log(parm2)^(j-1)
c  al1   : 1st argument of fit, log(parm1)  e.g. Te, Ti, ...
c  al2   : 2nd argument of fit, log(parm2)  e.g. ne, E0, ...

c  rc1min : left boundary of valid range for first parameter al1
c  rc1max : right boundary of valid range for first parameter al1
c  fp1    : parameters for extrapolation from valid range
c  ifex1mn: flag for choice of left (low end) extrapolation expression
c  ifex1mx: flag for choice of right (high end) extrapolation expression

c  rc2min : bottom boundary of valid range for second parameter al2
c  rc2max : top    boundary of valid range for second parameter al2
c  fp2    : parameters for extrapolation from valid range
c  ifex2mn: flag for choice of bottom (low end) extrapolation expression
c  ifex2mx: flag for choice of top (high end) extrapolation expression

c  trc    : print warnings in case of extrapolation beyond specified range

c  output:
c  cou   : log of value of fit(al1,al2)
c  dum   : dum(1,9), summation over 1st parameter is done,
c          dum provides 1D fit coefficients for 2nd parameter dependency (e.g. ne, E0)
c          for reduced 1D fit evaluation "on the fly", at fixed parameter AL1.
c          Also: dum(1) is the fit value at al2=0, e.g. at density=1e8 for amjuel.

      use EIRMOD_precision
      USE EIRMOD_COMPRT, ONLY: IUNOUT

      implicit none

      real(dp), intent(in) :: cf(9,9), fp1(6), fp2(6)
      real(dp), intent(in) :: al1, al2, rc1min, rc1max, rc2min, rc2max
      real(dp), intent(out) :: cou, dum(9)
      integer, intent(in) :: ifex1mn, ifex1mx, ifex2mn, ifex2mx

      real(dp) :: p1, p2, s01, s02, ds12, expo1, expo2, ccxm1, ccxm2,
     .            fpar1, fpar2, fpar3,
     .            al1min,al1max,cou1min,cou1max,
     .            al2min,al2max,cou2min,cou2max,
     .            eirene_extrap
      integer :: kk, jj, j, i, ii, ifex
      logical :: trc

      dum = 0._dp
      p1=al1
      p2=al2

      if (al1 < rc1min) then

        IF (IFEX1MN.EQ.0) THEN
C  FIT OUT OF VALID RANGE, AL1 < RC1MIN
C  DEFAULT: TAKE THE FIT AT AL1=RC1MIN
c  same as IFEXMN1=4 option further below
          p1=rc1min
          if (trc) write (iunout,*) 'extrap option 0 dbl_pol',al1,rc1min
          GOTO 100

        ELSEIF (IFEX1MN.LT.0) THEN
C  DETERMINE EXTRAPOLATION COEFFICIENTS FOR AL1 - LINEAR EXTRAP. IN LN(FIT) AT FIXED AL2
C  EVALUATED AT AL2, WHICH MUST BE INSIDE VALID RANGE
          S01=RC1MIN
          S02=LOG(1.25_DP)+RC1MIN
          DS12=S02-S01
          EXPO1=0.
          EXPO2=0.
c  evaluate double parameter fit at s01,al2, and at s02,al2
          DO J=1,9
            JJ=J-1
            DO I=1,9
              II=I-1
              EXPO1=EXPO1+S01**II*AL2**JJ*CF(I,J)
              EXPO2=EXPO2+S02**II*AL2**JJ*CF(I,J)
            ENDDO
          ENDDO
          CCXM1=EXPO1
          CCXM2=EXPO2
          FPAR1=CCXM1+(CCXM2-CCXM1)/DS12*(-S01)
          FPAR2=      (CCXM2-CCXM1)/DS12
          FPAR3=0.D0
C
          IFEX=3
          AL1MIN=RC1MIN
          COU1MIN=EXP(EXPO1)

        ELSEIF (IFEX1MN.GT.0) THEN
C  IFEX1MN IS GT 0, USE ONE OF THE PRE-PROGRAMMED EXTRAPOLATION SCHEMES
C  PROVIDE FIT EXPRESSION AT BOUNDARY RC1MIN
C  EVALUATED AT SECOND PARAMETER AL2, WHICH MUST BE INSIDE ITS VALID RANGE
          S01=RC1MIN
          EXPO1=0.
c  evaluate double parameter fit at s01,al2
          DO J=1,9
            JJ=J-1
            DO I=1,9
              II=I-1
              EXPO1=EXPO1+S01**II*AL2**JJ*CF(I,J)
            enddo
          enddo
          IFEX=IFEX1MN
          AL1MIN=RC1MIN
          COU1MIN=EXP(EXPO1)
          FPAR1=FP1(1)
          FPAR2=FP1(2)
          FPAR3=FP1(3)
        ENDIF

        COU=EIRENE_EXTRAP(AL1,AL1MIN,COU1MIN,IFEX,FPAR1,FPAR2,FPAR3)
        cou=log(cou)
C       DUM(1:9)= ???
        return

      elseif (al1 > rc1max) then
       IF (IFEX1MX.EQ.0) THEN
C  FIT OUT OF VALID RANGE, AL1 > RC1MAX
C  DEFAULT: TAKE THE FIT AT AL1=RC1MAX
c  same as IFEXMX1=4 option further below
          p1=rc1max
          if (trc) write (iunout,*) 'extrap option 0 dbl_pol',al1,rc1max
          GOTO 100

        ELSEIF (IFEX1MX.LT.0) THEN
C  DETERMINE EXTRAPOLATION COEFFICIENTS FOR AL1 - LINEAR EXTRAP. IN LN(FIT) AT FIXED AL2
C  EVALUATED AT AL2, WHICH MUST BE INSIDE VALID RANGE
          S01=RC1MAX
          S02=log(0.75_dp)+RC1MAX
          DS12=S02-S01
          EXPO1=0.
          EXPO2=0.
c  evaluate double parameter fit at s01,al2, and at s02,al2
cdr  ne or E dependence:
          DO J=1,9
            JJ=J-1
cdr  T-dependence:
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
          IFEX=3
          AL1MAX=RC1MAX
          COU1MAX=EXP(EXPO1)

        ELSEIF (IFEX1MX.GT.0) THEN
C  IFEX1MX IS GT 0, USE ONE OF THE PREPROGRAMMED EXTRAPOLATION SCHEMES
C  PROVIDE FIT EXPRESSION AT BOUNDARY RC1MAX
C  EVALUATED AT SECOND PARAMETER AL2, WHICH MUST BE INSIDE ITS VALID RANGE
          S01=RC1MAX
          EXPO1=0.
c  evaluate double parameter fit at s01,al2
          DO J=1,9
            JJ=J-1
            DO I=1,9
              II=I-1
              EXPO1=EXPO1+S01**II*AL2**JJ*CF(I,J)
             enddo
          enddo
          IFEX=IFEX1MX
          AL1MAX=RC1MAX        !  HERE: AL1MAX=AL1MAX(AL2)
          COU1MAX=EXP(EXPO1)   !  HERE: COU1MAX=COUMAX(AL2)
          FPAR1=FP1(4)
          FPAR2=FP1(5)
          FPAR3=FP1(6)
        ENDIF

        COU=EIRENE_EXTRAP(AL1,AL1MAX,COU1MAX,IFEX,FPAR1,FPAR2,FPAR3)
        cou=log(cou)
C       DUM(1:9)= ???
        return

      elseif (al2 < rc2min) then
        IF (IFEX2MN.EQ.0) THEN
C  FIT OUT OF VALID RANGE, AL1 < RC1MIN
C  DEFAULT: TAKE THE FIT AT AL1=RC1MIN
c  same as IFEXMN2=4 option further below
          p2=rc2min
          if (trc) write (iunout,*) 'extrap option 0 dbl_pol',al2,rc2min
          GOTO 100
        ELSEIF (IFEX2MN.LT.0) THEN
C  DETERMINE EXTRAPOLATION COEFFICIENTS FOR AL2 - LINEAR EXTRAP. IN LN(FIT) AT FIXED AL1
C  EVALUATED AT AL1, WHICH MUST BE INSIDE VALID RANGE
          S01=RC2MIN
          S02=LOG(1.25_DP)+RC2MIN
          DS12=S02-S01
          EXPO1=0.
          EXPO2=0.
c  evaluate double parameter fit at s01,al1, and at s02,al1
          DO J=1,9
            JJ=J-1
            DO I=1,9
              II=I-1
              EXPO1=EXPO1+AL1**II*S01**JJ*CF(I,J)
              EXPO2=EXPO2+AL1**II*S02**JJ*CF(I,J)
            ENDDO
          ENDDO
          CCXM1=EXPO1
          CCXM2=EXPO2
          FPAR1=CCXM1+(CCXM2-CCXM1)/DS12*(-S01)
          FPAR2=      (CCXM2-CCXM1)/DS12
          FPAR3=0.D0
C
          IFEX=3
          AL2MIN=RC2MIN
          COU2MIN=EXP(EXPO1)

        ELSEIF (IFEX2MN.GT.0) THEN
C  IFEX2MN IS GT 0, USE ONE OF THE PREPROGRAMMED EXTRAPOLATION SCHEMES
C  PROVIDE FIT EXPRESSION AT BOUNDARY RC2MIN
C  EVALUATED AT SECOND PARAMETER AL1, WHICH MUST BE INSIDE ITS VALID RANGE
          S01=RC2MIN
          EXPO1=0.
c  evaluate double parameter fit at s01,al1
          DO J=1,9
            JJ=J-1
            DO I=1,9
              II=I-1
              EXPO1=EXPO1+AL1**II*S01**JJ*CF(I,J)
            enddo
          enddo
          IFEX=IFEX2MN
          AL2MIN=RC2MIN       !  HERE: AL2MIN=AL2MIN(AL1)
          COU2MIN=EXP(EXPO1)  !  HERE: COU2MIN=COU2MIN(AL1)
          FPAR1=FP2(1)
          FPAR2=FP2(2)
          FPAR3=FP2(3)
        ENDIF


        COU=EIRENE_EXTRAP(AL2,AL2MIN,COU2MIN,IFEX,FPAR1,FPAR2,FPAR3)
        cou=log(cou)
C       DUM(1:9)= ???
        return

      elseif (al2 > rc2max) then
       IF (IFEX2MX.EQ.0) THEN
C  FIT OUT OF VALID RANGE, AL2 > RC2MAX
C  DEFAULT: TAKE THE FIT AT AL2=RC2MAX
c  same as IFEXMX2=4 option further below
          p2=rc2max
          if (trc) write (iunout,*) 'extrap option 0 dbl_pol',al2,rc2max
          GOTO 100
        ELSEIF (IFEX2MX.LT.0) THEN
C  DETERMINE EXTRAPOLATION COEFFICIENTS FOR AL2 - LINEAR EXTRAP. IN LN(FIT) AT FIXED AL1
C  EVALUATED AT AL1, WHICH MUST BE INSIDE VALID RANGE
          S01=RC2MAX
          S02=log(0.75_dp)+RC2MAX
          DS12=S02-S01
          EXPO1=0.
          EXPO2=0.
c  evaluate double parameter fit at s01,al1,  and at s02,al1
          DO J=1,9
            JJ=J-1
            DO I=1,9
              II=I-1
              EXPO1=EXPO1+AL1**II*S01**JJ*CF(I,J)
              EXPO2=EXPO2+AL1**II*S02**JJ*CF(I,J)
            enddo
          enddo
          CCXM1=EXPO1
          CCXM2=EXPO2
          FPAR1=CCXM1+(CCXM2-CCXM1)/DS12*(-S01)
          FPAR2=      (CCXM2-CCXM1)/DS12
          FPAR3=0.D0
C
          IFEX=3
          AL2MAX=RC2MAX
          COU2MAX=EXP(EXPO1)

        ELSEIF (IFEX2MX.GT.0) THEN
C  IFEX2MX IS GT 0, USE ONE OF THE PREPROGRAMMED EXTRAPOLATION SCHEMES
C  PROVIDE FIT EXPRESSION AT BOUNDARY RC2MAX
C  EVALUATED AT FIRST PARAMETER AL1, WHICH MUST BE INSIDE ITS VALID RANGE
          S01=RC2MAX
          EXPO1=0.
c  evaluate double parameter fit at s01,al1
          DO J=1,9
            JJ=J-1
            DO I=1,9
              II=I-1
              EXPO1=EXPO1+AL1**II*S01**JJ*CF(I,J)
            enddo
          enddo
          IFEX=IFEX2MX
          AL2MAX=RC2MAX        !  HERE: AL2MAX=AL1MAX(AL1)
          COU2MAX=EXP(EXPO1)   !  HERE: COU2MAX=COUMAX(AL1)
          FPAR1=FP2(4)
          FPAR2=FP2(5)
          FPAR3=FP2(6)
        ENDIF

        COU=EIRENE_EXTRAP(AL2,AL2MAX,COU2MAX,IFEX,FPAR1,FPAR2,FPAR3)
        cou=log(cou)
C       DUM(1:9)= ???
        return

      else

cdr  both parameters al1 and al2 are within valid range of fit
cdr  evaluate double polynomial fit.
cdr  result = cou
cdr  the reduced single parameter fit coefficients at fixed AL1 are returned on dum(1:9)

        goto 100

      endif


  100 CONTINUE

cdr  split the fit evaluation into two steps.
cdr  first:  p1=al1 dependence, to provide collapsed fit dum(..)
cdr  second: p2=al2 dependence using dum(..)

cdr  each of the 9 coefficients dum(jj) for the al2 dependence
cdr  is obtained by summing over all 9 terms
      do jj = 9, 1, -1
        dum(jj) = cf(9,jj)
        do kk = 8, 1, -1
          dum(jj) = dum(jj) * p1 + cf(kk,jj)
        end do
      end do

cdr  this second evaluation may not be needed, if only a collapsed fit is wanted,
cdr  or if al2=0.0 (as in H.4, H.10, AMJUEL fits at 1e8 density), for automatic corona limit.

cdr  H.4, H.10, H.12 fits from AMJUEL: corona at p2 <= log(ne/10**8) = rc2min = 0.0
cdr  if p2.le.0.0, just return cou=dum(1) = fit2(AL1,AL2)

      cou = dum(9)
      do jj = 8, 1, -1
        cou = cou * p2 + dum(jj)
      end do



      return
      end subroutine EIRENE_dbl_poly
