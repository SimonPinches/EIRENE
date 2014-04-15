!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced
!pb  24.11.06: get extrapolation parameters for polynomial fit only
!pb  07.12.06: double declaration of dsub removed
!dr  19.02.14: COMMENTS
 
      function EIRENE_rate_coeff (ir, p1, p2, lexp, iprshft, erate)
     .                     result (rate)
!  evaluate reaction rate coefficient (cm^3/s), and return this as "rate"

!  currently 5 different options controlled by :reacdat(ir)%rtc%ifit
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
!   p1:        first parameter (usually:  log temperature,...)
!   p2:        second parameter  (if any, e.g.  log (density),...,log(test particle energy),...) 
!   lexp:      return rate=rate coefficient in cm**3/sec
!   not lexp:  return rate=ln(rate coefficient) with rate-coefficient in cm**3/sec
!   iprshft:   >0: carry out shift in parameter p2, currently hard wired: 1e-8. (currently : only for ifit=2)

! to be done:  lexp option for ifit=4, ifit=5 not written.
!              remove erate in case of ifit=5 and generalize to more cr models.
!              iprshft option: currently hard wired only for ifit=2 and shift = 1e-8
 
      use EIRMOD_precision
      use EIRMOD_parmmod
      use EIRMOD_comxs
      use EIRMOD_comprt, only: iunout
 
      implicit none
 
      integer, intent(in) :: ir, iprshft
      real(dp), intent(in) :: p1, p2
      logical, intent(in) :: lexp
      real(dp), intent(out) :: erate
      real(dp) :: rate, eirene_sngl_poly, dum(9), rcmin, rcmax,
     .            fp(6), q1, q2,
     .            ALPCR, SCR, SCRRAD, E_ALPCR, E_SCR, E_SCRRAD,
     .            E_ALPCR_T, E_SCR_T, E_SCRRAD_T
      real(dp), save :: xlog10e, xln10, dsub
      real(dp), allocatable, save :: pop0(:), pop1(:), pop2(:), qcol2(:)
      integer :: jfexmn, jfexmx
      integer, save :: ifirst=0, ifsub=0
 
      interface
        function EIRENE_intp_adas (ad,p1,p2) result(res)
          use EIRMOD_precision
          use EIRMOD_comxs, only: adas_data
          type(adas_data), pointer :: ad
          real(dp), intent(in) :: p1, p2
          real(dp) :: res
        end function EIRENE_intp_adas
 
        function EIRENE_intp_table (tb,p1,p2) result(res)
          use EIRMOD_precision
          use EIRMOD_comxs, only: hydkin_data
          type(hydkin_data), pointer :: tb
          real(dp), intent(in) :: p1, p2
          real(dp) :: res
        end function EIRENE_intp_table
      end interface
 
 
      if (.not.reacdat(ir)%lrtc) then
        write (iunout,*) ' no data for rate coefficient available',
     .                    ' for reaction ',ir
        call EIRENE_exit_own(1)
      end if
 
      rate = 0._dp
      erate = 0._dp

c  extrapolation data: currently only for polynomial fits 
      if ((reacdat(ir)%rtc%ifit == 1) .or.
     .    (reacdat(ir)%rtc%ifit == 2)) then
        rcmin  = reacdat(ir)%rtc%poly%rcmn
        rcmax  = reacdat(ir)%rtc%poly%rcmx
        fp     = reacdat(ir)%rtc%poly%fparm
        jfexmn = reacdat(ir)%rtc%poly%ifexmn
        jfexmx = reacdat(ir)%rtc%poly%ifexmx
      end if
 
      if (mod(iftflg(ir,2),100) == 10) then

!  SET A CONSTANT RATE 
        rate = reacdat(ir)%rtc%poly%dblpol(1,1)
 
      elseif (reacdat(ir)%rtc%ifit == 1) then

!  SINGLE POLYNOMIAL FIT 
        rate = eirene_sngl_poly(reacdat(ir)%rtc%poly%dblpol(1:9,1),p1,
     .                   rcmin, rcmax, fp, jfexmn, jfexmx)
        if (lexp) rate = exp(max(-100._dp,rate))
 
      else if (reacdat(ir)%rtc%ifit == 2) then

! DOUPLE POLYNOMIAL FIT
        if (ifsub == 0) then
          ifsub = 1
          dsub = log(1.e8_dp)
        end if
c  rescale parameter p2  (currently only by 1e-8):  q2 
        q2 = p2
        if (iprshft > 0) q2 = q2 - dsub
 
        call EIRENE_dbl_poly
     .       (reacdat(ir)%rtc%poly%dblpol,p1,q2,rate,dum,1,9,
     .                rcmin, rcmax, fp, jfexmn, jfexmx)
        if (lexp) rate = exp(max(-100._dp,rate))
 
      else if (reacdat(ir)%rtc%ifit == 3) then
 
! DOUBLE PARAMETER TABLE  (E.G. ADAS)

!  currently hard wired:  input parameters q1, q2 and table coefficients are log10
        if (ifirst == 0) then
          ifirst = 1
          xln10 = log(10._dp)
          xlog10e = 1._dp/xln10
        end if
c  convert parameters p1 and p2 from ln to log10:  q1,q2 
        q1 = xlog10e*p1
        q2 = xlog10e*p2
C  assume here: tabulated data are log10  (to be generalized)
        rate = eirene_intp_adas(reacdat(ir)%rtc%adas,q1,q2)
 
        if (lexp) then
          rate=10._dp**rate
        else
          rate = xln10*rate
        end if
 
      else if (reacdat(ir)%rtc%ifit == 4) then
 
! SINGLE PARAMETER TABLE  (E.G. HYDKIN)
!  currently hard wired:  input parameters q1 and table coefficients are neither ln nor log10
 
        q1 = exp(p1)
C  assume here: tabulated data are not ln nor log10  (to be generalized)
        rate = eirene_intp_table(reacdat(ir)%rtc%hyd,q1,p2)

!  lexp option not connected here !
 
      else if (reacdat(ir)%rtc%ifit == 5) then
 
! H-colrad   RATE AND ENERGY LOSS RATE IN ONE SIGNLE STEP
 
        if (.not.allocated(pop0)) then
          allocate(pop0(40))
          allocate(pop1(40))
          allocate(pop2(40))
          allocate(qcol2(40))
        end if
 
        QCOL2 = 0._DP
c  convert parameters p1, p2 to exp(p1), exp(p2):  q1,q2
        Q1 = EXP(P1)
        Q2 = EXP(P2)
        CALL EIRENE_H_COLRAD(Q1, Q2 ,QCOL2,POP0,POP1,POP2,
     .                ALPCR,    SCR,    SCRRAD,
     .                E_ALPCR,  E_SCR,  E_SCRRAD,
     .                E_ALPCR_T,E_SCR_T,E_SCRRAD_T)

!  lexp option not connected here !
 
        rate = log(scr)
        erate = log(-e_scr)
 
      end if
 
      return
 
      end function EIRENE_rate_coeff
 
