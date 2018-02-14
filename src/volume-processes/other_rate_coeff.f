!pb  01.06.2017  copied from rate_coeff.f

      function EIRENE_other_rate_coeff(ir, ic, p1, p2, lexp, iprshft)
     .                     result (rate)

!  evaluate rate coefficient for other reactions 
!  and return this as "rate"

!  currently 5 different options controlled by 'reacdat(ir)%rtc%ifit'
!  Only ifit=2 and ifit=3 tested so far. Caution!
!  ifit=1:   single polynom fit, use P1, (e.g. HYDHEL, H.2)
!  ifit=2:   double polynom fit, use P1, P2, (e.g. HYDHEL, H.3, AMJUEL, H.4,...)
!  ifit=3:   interpolation in 2-parameter table (e.g. ADAS)
!  ifit=4:   interpolation in single parameter table (e.g. open ADAS, HYDKIN,....)
!  ifit=5:   use internal eirene collision radiative code. To be generalized
!            (currently here also energy rates, erate  for this particular option. 
!            More logical if the latter are moved
!            to routine "eirene_energy-rate-coeff"

!   input:
!   ir:        reaction number, as stored in eirene arrays.
!   p1:        first parameter (usually:  log_e temperature,...)
!   p2:        second parameter  (if any, e.g.  log_e (density),...,log_e(test particle energy),...) 
!   lexp:      return rate=rate coefficient in cm**3/sec
!   not lexp:  return rate=log_e(rate coefficient) with rate-coefficient in cm**3/sec
!   iprshft:   >0: carry out shift in parameter p2, 
!              currently hard wired: 1e-8. 
!              Currently : only for ifit=2, polynomial fits vs. ne, T, ne in units 1e8 *cm**-3.
!              emissivity.f relies on the current use of iprshft in the tested cases!

! to be done:  lexp option for ifit=4, ifit=5 not written.
!              remove erate in case of ifit=5 and generalize to more cr models.
!              iprshft option: currently hard wired only for ifit=2 and shift = 1e-8
!              what happens if later call with other shift ?  coding to be reconsidered !
!              remove ifirst and ifsub conditions and set the data once, and save. 
 
      use EIRMOD_precision
      use EIRMOD_parmmod
      use EIRMOD_comxs
      use EIRMOD_comprt, only: iunout
 
      implicit none
 
      integer, intent(in) :: ir, iprshft, ic
      real(dp), intent(in) :: p1, p2
      logical, intent(in) :: lexp
      real(dp) :: rate, EIRENE_sngl_poly, dum(9), rc1min, rc1max,
     .            fp1(6), fp2(6), pp1, pp2, rc2min, rc2max,
     .            rrc2min, rrc2max, SCR
      real(dp), save :: xlog10e =  4.34294482d-01,      !1./ln(10) = log10(e)
     .                  xln10   =  2.30258509299_dp,    !ln(10) 
     .                  dsub    = 18.420680744_dp       !ln(1e8)

      integer :: jfex1mn, jfex1mx,jfex2mn, jfex2mx
      integer :: ip1, ip2, iflavor, ivar           
 
      interface
        function EIRENE_intp_adas (ad,p1,p2,ip1,ip2) result(res)
          use EIRMOD_precision
          use EIRMOD_comxs, only: adas_data
          type(adas_data), pointer :: ad
          real(dp), intent(in) :: p1, p2
          integer, intent(out) :: ip1,ip2
          real(dp) :: res
        end function EIRENE_intp_adas
 
        function EIRENE_intp_table (tb,p1,ip1) result(res)
          use EIRMOD_precision
          use EIRMOD_comxs, only: hydkin_data
          type(hydkin_data), pointer :: tb
          real(dp), intent(in) :: p1
          integer, intent(out) :: ip1
          real(dp) :: res
        end function EIRENE_intp_table
      end interface
 
 
      if (.not.reacdat(ir)%loth) then
        write (iunout,*) ' no data for other reaction available',
     .                   ' for reaction ',ir
        call EIRENE_exit_own(1)
      end if
 
      rate = 0._dp

c.............................................................

 
      if (mod(iftflg(ir,2),100) == 10) then

!  SET A CONSTANT RATE 
        rate = reacdat(ir)%oth%poly%dblpol(1,1)

cdr   lexp missing

c.............................................................
 
      elseif (reacdat(ir)%oth%ifit == 1) then

!  SINGLE POLYNOMIAL FIT VS. P1 (TEMPERATURE), FOR LN OF ENERGY WEIGHTED RATE
 
c  extrapolation data:  for 1d polynomial fits 
        rc1min  = reacdat(ir)%oth%rc1min
        rc1max  = reacdat(ir)%oth%rc1max
        fp1(1:3)= reacdat(ir)%oth%fp1l
        fp1(4:6)= reacdat(ir)%oth%fp1r
        jfex1mn = reacdat(ir)%oth%jfex1mn
        jfex1mx = reacdat(ir)%oth%jfex1mx
 
        rate = eirene_sngl_poly(reacdat(ir)%oth%poly%dblpol(1:9,1),
     .                   p1, rc1min, rc1max, fp1, jfex1mn, jfex1mx)

C       if (.not. lexp)  rate=rate
        if (lexp)        rate = exp(max(-100._dp,rate))

c..............................................................

 
      else if (reacdat(ir)%oth%ifit == 2) then

!  DOUBLE POLYNOMIAL FIT VS. P1 (TEMPERATURE) AND P2,  FOR LN OF RATE 

c  extrapolation data:  for 2d polynomial fits 
        rc1min  = reacdat(ir)%oth%rc1min
        rc1max  = reacdat(ir)%oth%rc1max
        rc2min  = reacdat(ir)%oth%rc2min
        rc2max  = reacdat(ir)%oth%rc2max
        fp1(1:3)= reacdat(ir)%oth%fp1l
        fp1(4:6)= reacdat(ir)%oth%fp1r
        fp2(1:3)= reacdat(ir)%oth%fp2b
        fp2(4:6)= reacdat(ir)%oth%fp2t
        jfex1mn = reacdat(ir)%oth%jfex1mn
        jfex1mx = reacdat(ir)%oth%jfex1mx
        jfex2mn = reacdat(ir)%oth%jfex2mn
        jfex2mx = reacdat(ir)%oth%jfex2mx
 

c  rescale parameter p2  (currently only by 1e-8 for density):  pp2 
        pp2 = p2
        if (iprshft > 0) then
          pp2 = pp2 - dsub
          rrc2min=rc2min - dsub
          rrc2max=rc2max - dsub
        endif
cdr     write (6,*) 'particle rate '
 
        call EIRENE_dbl_poly
     .       (reacdat(ir)%oth%poly%dblpol,p1,pp2,rate,dum,
     .        rc1min, rc1max, fp1, jfex1mn, jfex1mx,
     .        rc2min, rc2max, fp2, jfex2mn, jfex2mx)

C       if (.not. lexp)  rate=rate
        if (lexp)        rate = exp(max(-100._dp,rate))

c.............................................................. 

      else if (reacdat(ir)%oth%ifit == 3) then

! 2D TABULAR INPUT,  FOR LOG10 OF RATE,  cm^3/s
! E.G.: ADAS adf11 ACD and SCD FILES
cdr  extrapolation data: for 2d tabulated data, option not ready
cdr  to be added here

!  currently hard wired:  input parameters pp1, pp2 and table coefficients are log10

c  convert parameters p1 and p2 from ln to log10:  pp1,pp2 
        pp1 = xlog10e*p1
        pp2 = xlog10e*p2
C  assume here: tabulated data are log10  (to be generalized)
        rate = eirene_intp_adas(reacdat(ir)%oth%adas,pp1,pp2,ip1,ip2)
 
        if (lexp) then
          rate=10._dp**rate
        else
          rate = xln10*rate     !    convert from log10(rate) to ln(rate)

        end if

cdr this unit conversion must be wrong in case lexp !!


c..............................................................
  
      else if (reacdat(ir)%oth%ifit == 4) then
 
! SINGLE PARAMETER TABLE  (E.G. HYDKIN)
cdr  extrapolation data: for 1d tabulated data:  option not ready (only CxHy data ?) 
cdr  to be added here
! currently hard wired:  input parameters q1 and table coefficients are neither ln nor log10
 
        pp1 = exp(p1)
C  assume here: tabulated data are neither ln nor log10  (to be generalized)
        rate = eirene_intp_table(reacdat(ir)%oth%hyd,pp1,ip1)

!  lexp option not connected here !

c..............................................................
 
      else if (reacdat(ir)%oth%ifit == 5) then

! INTERNAL COLLISION RADIATIVE CODE
 
c  convert parameters p1, p2 to exp(p1), exp(p2):  PP1,PP2
        PP1 = EXP(P1)
        PP2 = EXP(P2)
        iflavor = reacdat(ir)%oth%crm%iflav
        ivar = reacdat(ir)%oth%crm%ivarst

        CALL EIRENE_COLRAD(IR, IC, IFLAVOR, IVAR, PP1, PP2, SCR)

!  lexp option was not connected here, but used in xstei.f ! corrected, Oct. 28th 2015

        rate=scr 
        if (.not.lexp) rate = log(scr)
         
      end if
 
      return
 
      end function EIRENE_other_rate_coeff
