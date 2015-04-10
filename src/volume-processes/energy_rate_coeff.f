!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced
!pb  24.11.06: get extrapolation parameters for polynomial fit only
!pb  30.11.06: divide energy rate coefficient by ELCHA to get correct units

c  to be done: h_colrad called twice per cell ??
c              re-use erate from previous call to rate_coeff
 
      function EIRENE_energy_rate_coeff (ir, p1, p2, lexp, iprshft)
     .                           result (rate)

cdr  return an energy weighted rate, eV/s per incident particle  

      use EIRMOD_precision
      use EIRMOD_parmmod
      use EIRMOD_comxs
      use EIRMOD_ccona
      use EIRMOD_comprt, only: iunout
 
      implicit none
 
      integer, intent(in) :: ir, iprshft
      real(dp), intent(in) :: p1, p2
      logical, intent(in) :: lexp
      real(dp) :: rate, EIRENE_sngl_poly, dum(9), rcmin, rcmax, fp(6), 
     .            q1, q2,
     .            ALPCR, SCR, SCRRAD, E_ALPCR, E_SCR, E_SCRRAD,
     .            E_ALPCR_T, E_SCR_T, E_SCRRAD_T
      real(dp), save :: xlog10e =  4.34294482d-01,       !1./ln(10) = log10(e) 
     .                  xln10   =  2.30258509299_dp,     !ln(10)  
     .                  dsub    = 18.420680744_dp,       !ln(1e8) 
     .                  xlnelch =-43.2777390821          !ln(elcha) 
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
      end interface
 
      if (.not.reacdat(ir)%lrtcew) then
        write (iunout,*) ' no data for energy weighted rate',
     .                   ' coefficient available for reaction ',ir
        call EIRENE_exit_own(1)
      end if
 
      rate = 0._dp
 
      if ((reacdat(ir)%rtcew%ifit == 1) .or.
     .    (reacdat(ir)%rtcew%ifit == 2)) then
        rcmin  = reacdat(ir)%rtcew%poly%rcmn
        rcmax  = reacdat(ir)%rtcew%poly%rcmx
        fp     = reacdat(ir)%rtcew%poly%fparm
        jfexmn = reacdat(ir)%rtcew%poly%ifexmn
        jfexmx = reacdat(ir)%rtcew%poly%ifexmx
      end if
 
      if (mod(iftflg(ir,2),100) == 10) then
 
        rate = reacdat(ir)%rtcew%poly%dblpol(1,1)
 
      elseif (reacdat(ir)%rtcew%ifit == 1) then
 
        rate = EIRENE_sngl_poly(reacdat(ir)%rtcew%poly%dblpol(1:9,1),p1,
     .                   rcmin, rcmax, fp, jfexmn, jfexmx)
        if (lexp) rate = exp(max(-100._dp,rate))
 
      else if (reacdat(ir)%rtcew%ifit == 2) then
 
 
        q2 = p2
        if (iprshft > 0) q2 = q2 - dsub
 
        call EIRENE_dbl_poly
     .                (reacdat(ir)%rtcew%poly%dblpol,p1,q2,rate,dum,1,9,
     .                rcmin, rcmax, fp, jfexmn, jfexmx)
        if (lexp) rate = exp(max(-100._dp,rate))
 
      else if (reacdat(ir)%rtcew%ifit == 3) then
 
! ADAS
 
        q1 = xlog10e*p1
        q2 = xlog10e*p2
        rate = EIRENE_intp_adas(reacdat(ir)%rtcew%adas,q1,q2)
 
        if (lexp) then
          rate=10._dp**rate
        else
          rate = xln10*rate
        end if
        rate = rate - xlnelch
 
      else if (reacdat(ir)%rtcew%ifit == 5) then
 
! H-colrad
 
        if (.not.allocated(pop0)) then
          allocate(pop0(40))
          allocate(pop1(40))
          allocate(pop2(40))
          allocate(qcol2(40))
        end if
 
        QCOL2 = 0._DP
        q1 = exp(p1)
        q2 = exp(p2)
        CALL EIRENE_H_COLRAD(P1, P2 ,Q2,POP0,POP1,POP2,
     .                ALPCR,    SCR,    SCRRAD,
     .                E_ALPCR,  E_SCR,  E_SCRRAD,
     .                E_ALPCR_T,E_SCR_T,E_SCRRAD_T)
 
        rate = log(e_scr)
 
      end if
 
 
      return
 
      end function EIRENE_energy_rate_coeff
 
