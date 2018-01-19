!pb  01.06.2017  copied from rate_coeff.f

c  to be done: h_colrad called twice per cell ??
c              re-use erate from previous call to rate_coeff

cdr  19.02.14: COMMENTS
cdr sept.15:  lexp not fully written, in case of adas 2d tables
cdr           ifit=4 option was missing (1D tables). added, but not checked.





      function EIRENE_other_rate_coeff (ir, p1, p2, lexp, iprshft)
     .                     result (orate)

!  evaluate other atomic data:  mostly: population coefficients, density ratios
!                 e.g. from AMJUEL H.11, H.12 sections
!  and return this as "orate"

!  currently 5 different options controlled by 'reacdat(ir)%rtc%ifit'
!  ifit=1:   single polynom fit, use P1, (e.g. AMJUEL, H.11)
!  ifit=2:   double polynom fit, use P1, P2, (e.g. AMJUEL, H.12,...)
!  ifit=3:   interpolation in 2-parameter table (e.g. ADAS)
!  ifit=4:   interpolation in single parameter table (e.g. open ADAS, HYDKIN,....)
!  ifit=5:   use internal eirene collision radiative code. To be generalized
!            (currently here also other rates, orate  for this particular option.
!            More logical if the latter are moved
!            to routine "eirene_energy-rate-coeff"

!   input:
!   ir:        reaction number, as stored in eirene arrays.
!   p1:        first parameter (usually:  log_e temperature,...)
!   p2:        second parameter  (if any, e.g.  log_e (density),...,log_e(test particle energy),...)
!   lexp:      return orate=rate coefficient in ... units
!   not lexp:  return orate=log_e(rate coefficient) with rate-coefficient in ...units
!   iprshft:   >0: carry out shift in parameter p2 for fit expression evaluation,
!              currently hard wired: 1e-8.
!             (currently : only for ifit=2, polynomial fits vs. ne, T, ne in units 1e8 *cm**-3)

! to be done:  lexp option for ifit=4, ifit=5 not written.
!              remove orate in case of ifit=5 and generalize to more cr models.
!              iprshft option: currently hard wired only for ifit=2 and shift = 1e-8
!              what happens if later call with other shift ?  coding to be reconsidered !

!              remove ifirst and ifsub conditions and set the data once, and save.  DONE (Nov. 15)

      use EIRMOD_precision
      use EIRMOD_parmmod
      use EIRMOD_comxs
      use EIRMOD_ccona
      use EIRMOD_ctrcei, only: trcamd
      use EIRMOD_comprt, only: iunout

      implicit none

      integer, intent(in) :: ir, iprshft
      real(dp), intent(in) :: p1, p2
      logical, intent(in) :: lexp

      real(dp) :: orate, EIRENE_sngl_poly, dum(9),
     .            pp1, rc1min,  rc1max,fp1(6),
     .            pp2, rc2min,  rc2max,fp2(6),
     .                 rrc2min, rrc2max,
     .            ALPCR, SCR, SCR_EXT, E_ALPCR, E_SCR, E_SCR_EXT,
     .            E_ALPCR_T, E_SCR_T, E_SCR_EXT_T
      real(dp), save :: xlog10e =  4.34294482d-01,      !1./ln(10) = log10(e)
     .                  xln10   =  2.30258509299_dp,    !ln(10)
     .                  dsub    = 18.420680744_dp       !ln(1e8)

      real(dp), allocatable, save :: pop0(:), pop1(:), pop2(:), q_ext(:)
      integer :: jfex1mn, jfex1mx,jfex2mn, jfex2mx
      integer :: ic,ip1,ip2

      interface
        function EIRENE_intp_tab2d (ad,p1,p2,ip1,ip2) result(res)
          use EIRMOD_precision
          use EIRMOD_comxs, only: adas_data
          type(adas_data), pointer :: ad
          real(dp), intent(in) :: p1, p2
          integer, intent(out) :: ip1,ip2
          real(dp) :: res
        end function EIRENE_intp_tab2d

        function EIRENE_intp_tab1d (tb,p1,ip1) result(res)
          use EIRMOD_precision
          use EIRMOD_comxs, only: hydkin_data
          type(hydkin_data), pointer :: tb
          real(dp), intent(in) :: p1
          integer, intent(out) :: ip1
          real(dp) :: res
        end function EIRENE_intp_tab1d
      end interface


      if (.not.reacdat(ir)%loth) then
        write (iunout,*) ' no data for other reaction available',
     .                   ' for reaction ',ir
        call EIRENE_exit_own(1)
      end if

      orate = 0._dp

c.............................................................


      if (mod(iftflg(ir,2),100) == 10) then

!  SET A CONSTANT RATE
        orate = reacdat(ir)%oth%poly%dblpol(1,1)

cdr   missing: iftflg < 100:  multiply density,  else: not
cdr   lexp missing

c.............................................................

      elseif (reacdat(ir)%oth%ifit == 1) then

!  SINGLE POLYNOMIAL FIT VS. P1 =LN(TEMPERATURE), FOR LN(OTHER RATE)

c  extrapolation data:  for 1d polynomial fits
        rc1min  = reacdat(ir)%oth%rc1min
        rc1max  = reacdat(ir)%oth%rc1max
        fp1(1:3)= reacdat(ir)%oth%fp1l
        fp1(4:6)= reacdat(ir)%oth%fp1r
        jfex1mn = reacdat(ir)%oth%jfex1mn
        jfex1mx = reacdat(ir)%oth%jfex1mx

        orate = eirene_sngl_poly(reacdat(ir)%oth%poly%dblpol(1:9,1),
     .                   p1,rc1min,rc1max,fp1,jfex1mn,jfex1mx,
     .                   trcamd)

C       if (.not. lexp)  orate=orate
        if (lexp)        orate = exp(max(-100._dp,orate))

c..............................................................


      else if (reacdat(ir)%oth%ifit == 2) then

!  DOUBLE POLYNOMIAL FIT VS. P1 =LN(TEMPERATURE) AND P2,  FOR LN OF RATE

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
cdr     write (6,*) 'other rate '

        call EIRENE_dbl_poly
     .       (reacdat(ir)%oth%poly%dblpol,p1,pp2,orate,dum,
     .        rc1min,  rc1max,  fp1, jfex1mn, jfex1mx,
     .        rrc2min, rrc2max, fp2, jfex2mn, jfex2mx,
     .        trcamd)

C       if (.not. lexp)  orate=orate
        if (lexp)        orate = exp(max(-100._dp,orate))

c..............................................................

      else if (reacdat(ir)%oth%ifit == 3) then

! 2D TABULAR INPUT,  FOR LOG10 OF RATE,  cm^3/s
! E.G.: ADAS FILES
cdr  extrapolation data: for 2d tabulated data, option not ready
cdr  to be added here

!  currently hard wired:  input parameters pp1, pp2 and table coefficients are log10

c  convert parameters p1 and p2 from ln to log10:  pp1,pp2
        pp1 = xlog10e*p1
        pp2 = xlog10e*p2
C  assume here: tabulated data are log10  (to be generalized)
        orate = eirene_intp_tab2d(reacdat(ir)%oth%adas,pp1,pp2,ip1,ip2)

        if (lexp) then
          orate=10._dp**orate
        else
          orate = xln10*orate     !    convert from log10(orate) to ln(orate)

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
        orate = eirene_intp_tab1d(reacdat(ir)%oth%hyd,pp1,ip1)

!  lexp option not connected here !

c..............................................................

      else if (reacdat(ir)%oth%ifit == 5) then

! INTERNAL COLLISION RADIATIVE CODE

        if (.not.allocated(pop0)) then
          allocate(pop0(40))
          allocate(pop1(40))
          allocate(pop2(40))

          allocate(q_ext(40))    !   e.g. photo excitation rate for H*(n)
        end if

        Q_EXT = 0._DP

c  convert parameters p1, p2 to exp(p1), exp(p2):  PP1,PP2
        PP1 = EXP(P1)
        PP2 = EXP(P2)
        CALL EIRENE_H_COLRAD(IC, PP1, PP2 ,Q_EXT,POP0,POP1,POP2,
     .                ALPCR,    SCR,    SCR_EXT,
     .                E_ALPCR,  E_SCR,  E_SCR_EXT,
     .                E_ALPCR_T,E_SCR_T,E_SCR_EXT_T)

!  lexp option was not connected here, but used in xstei.f ! corrected, Oct. 28th 2015

cdr     IF (.NOT.LEXP) orate = log(-o_scr)
cdr     IF (LEXP)      orate = -O_SCR
        orate=1.0

      end if


      return

      end function EIRENE_other_rate_coeff
