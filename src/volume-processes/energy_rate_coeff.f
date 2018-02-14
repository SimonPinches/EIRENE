!pb  22.11.06: flag iprshft for shift of second parameter to rate_coeff introduced
!pb  24.11.06: get extrapolation parameters for polynomial fit only
!pb  30.11.06: divide energy weighted rate coefficient by ELCHA to get correct units

c  to be done: h_colrad called twice per cell ??
c              re-use erate from previous call to rate_coeff

cdr  19.02.14: COMMENTS
cdr sept.15:  lexp not fully written, in case of adas 2d tables
cdr           also unit conversion incorrect in that case. --> tbd
cdr           ifit=4 option was missing (1D tables). added, but not checked.

cdr  16.11.15: bug fix: error in arguments in call to H_colrad
cdr            added: ifit=4 option (1d table interpolation)
cdr            additional parameters in calls to 1d and 2d table interpolation
cdr  26.11.15: additional parameter IC in call to H_colrad,
cdr            for later use to identify "visited cells"
cdr  sept. 16: started to add extrapolation options. not ready....

      function EIRENE_energy_rate_coeff (ir, p1, p2, lexp, iprshft)
     .                            result (erate)

!  evaluate energy weighted rate coefficient, eV/s per incident particle,
!  and return this as "erate"

!  currently 5 different options controlled by 'reacdat(ir)%rtcew%ifit'
!  ifit=1:   single polynom fit, use P1, (e.g. HYDHEL, AMJUEL, H.8)
!  ifit=2:   double polynom fit, use P1, P2, (e.g. HYDHEL, H.9, AMJUEL, H.10,...)
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
!   lexp:      return erate=energy weighted rate coefficient in eV*cm**3/sec
!   not lexp:  return erate=log_e(erate coefficient) with rate-coefficient in cm**3/sec
!   iprshft:   >0: carry out shift in parameter p2 for fit expression evaluation,
!                  currently hard wired: 1e-8. 
!                 (currently : only for ifit=2, polynomial fits vs. ne, T, ne in units 1e8 *cm**-3)

! to be done:  lexp option for ifit=4, ifit=5 not written.
!              erate in case of ifit=5 hard wired to E_scr. 
!              What happens in case of recombination ?
!              and generalize to more cr models.
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

      real(dp) :: erate, EIRENE_sngl_poly, dum(9),
     .            pp1, rc1min,  rc1max, fp1(6),
     .            pp2, rc2min,  rc2max, fp2(6),
     .                 rrc2min, rrc2max,
     .            ALPCR, SCR, SCR_EXT, E_ALPCR, E_SCR, E_SCR_EXT,
     .            E_ALPCR_T, E_SCR_T, E_SCR_EXT_T
      real(dp), save :: xlog10e =  4.34294482d-01,      !1./ln(10) = log10(e)
     .                  xln10   =  2.30258509299_dp,    !ln(10) 
     .                  dsub    = 18.420680744_dp,      !ln(1e8), hard wired. But should come from database
     .                  xlnelch =-43.2777390821         !ln(elcha) 
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
 
      if (.not.reacdat(ir)%lrtcew) then
        write (iunout,*) ' no data for energy weighted rate',
     .                   ' coefficient available for reaction ',ir
        call EIRENE_exit_own(1)
      end if
 
      erate = 0._dp

c.............................................................

 
      if (mod(iftflg(ir,4),100) == 10) then

!  SET A CONSTANT RATE 
        erate = reacdat(ir)%rtcew%poly%dblpol(1,1)

cdr missing: iftflg < 100:  multiply density,  else: not
cdr   lexp missing

c.............................................................

      elseif (reacdat(ir)%rtcew%ifit == 1) then

!  SINGLE POLYNOMIAL FIT VS. P1 =LN(TEMPERATURE), FOR LN(ENERGY WEIGHTED RATE)
 
c  extrapolation data: for 1d polynomial fits 
        rc1min  = reacdat(ir)%rtcew%rc1min
        rc1max  = reacdat(ir)%rtcew%rc1max
        fp1(1:3)= reacdat(ir)%rtcew%fp1l
        fp1(4:6)= reacdat(ir)%rtcew%fp1r
        jfex1mn = reacdat(ir)%rtcew%jfex1mn
        jfex1mx = reacdat(ir)%rtcew%jfex1mx
 
        erate = eirene_sngl_poly(reacdat(ir)%rtcew%poly%dblpol(1:9,1),
     .                           p1,rc1min,rc1max,fp1,jfex1mn,jfex1mx,
     .                           trcamd)

C       if (.not. lexp)  erate=erate
        if (lexp)        erate = exp(max(-100._dp,erate))

c..............................................................

 
      else if (reacdat(ir)%rtcew%ifit == 2) then

!  DOUBLE POLYNOMIAL FIT VS. P1 =LN(TEMPERATURE) AND P2,  FOR LN(ENERGY WEIGHTED RATE)

c  extrapolation data:  for 2d polynomial fits 
        rc1min  = reacdat(ir)%rtcew%rc1min
        rc1max  = reacdat(ir)%rtcew%rc1max
        rc2min  = reacdat(ir)%rtcew%rc2min
        rc2max  = reacdat(ir)%rtcew%rc2max
        fp1(1:3)= reacdat(ir)%rtcew%fp1l
        fp1(4:6)= reacdat(ir)%rtcew%fp1r
        fp2(1:3)= reacdat(ir)%rtcew%fp2b
        fp2(4:6)= reacdat(ir)%rtcew%fp2t
        jfex1mn = reacdat(ir)%rtcew%jfex1mn
        jfex1mx = reacdat(ir)%rtcew%jfex1mx
        jfex2mn = reacdat(ir)%rtcew%jfex2mn
        jfex2mx = reacdat(ir)%rtcew%jfex2mx


c  rescale parameter p2  (currently only by 1e-8 for density):  pp2 
        pp2 = p2
        if (iprshft > 0) then
          pp2 = pp2 - dsub
          rrc2min=rc2min - dsub
          rrc2max=rc2max - dsub
        endif
cdr     write (6,*) 'energy rate '

        call EIRENE_dbl_poly
     .       (reacdat(ir)%rtcew%poly%dblpol,p1,pp2,erate,dum,
     .        rc1min,  rc1max,  fp1, jfex1mn, jfex1mx,
     .        rrc2min, rrc2max, fp2, jfex2mn, jfex2mx,
     .        trcamd)

C       if (.not. lexp)  erate=erate
        if (lexp)        erate = exp(max(-100._dp,erate))

c..............................................................

      else if (reacdat(ir)%rtcew%ifit == 3) then

! 2D TABULAR INPUT,  FOR LOG10 OF ENERGY WEIGHTED RATE,  joule*cm^3/s
! E.G.: ADAS adf11 PLT and PRB FILES
cdr  extrapolation data: for 2d tabulated data, option not ready
cdr  to be added here

!  currently hard wired:  input parameters pp1, pp2 and table coefficients are log10

c  convert parameters p1 and p2 from ln to log10:  pp1,pp2 
        pp1 = xlog10e*p1
        pp2 = xlog10e*p2
C  assume here: tabulated data are log10  (to be generalized)
        erate=eirene_intp_tab2d(reacdat(ir)%rtcew%adas,pp1,pp2,ip1,ip2)
 
        if (lexp) then
          erate=10._dp**erate
        else
          erate = xln10*erate     !    convert from log10(erate) to ln(erate)
          erate = erate - xlnelch !    convert ln(erate) from joule cm^3/s to eV cm^3/s
        end if
cdr this unit conversion must be wrong in case lexp !!
     

c..............................................................
  
      else if (reacdat(ir)%rtcew%ifit == 4) then
 
! SINGLE PARAMETER TABLE  (E.G. HYDKIN)
cdr  extrapolation data: for 1d tabulated data:  option not ready (only CxHy data ?) 
cdr  to be added here

! currently hard wired:  input parameters q1 and table coefficients are neither ln nor log10
 
        pp1 = exp(p1)
C  assume here: tabulated data are neither ln nor log10  (to be generalized)
        erate = eirene_intp_tab1d(reacdat(ir)%rtcew%hyd,pp1,ip1)

!  lexp option not connected here !

c..............................................................

      else if (reacdat(ir)%rtcew%ifit == 5) then
 
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
        CALL EIRENE_H_COLRAD(IC,PP1, PP2 ,Q_EXT,POP0,POP1,POP2,
     .                ALPCR,    SCR,    SCR_EXT,
     .                E_ALPCR,  E_SCR,  E_SCR_EXT,
     .                E_ALPCR_T,E_SCR_T,E_SCR_EXT_T)

!  lexp option was not connected here, but used in xstei.f ! corrected, Oct. 28th 2015
 
        IF (.NOT.LEXP) erate = log(-e_scr)
        IF (LEXP)      erate = -E_SCR
 
      end if
 
 
      return
 
      end function EIRENE_energy_rate_coeff
 
