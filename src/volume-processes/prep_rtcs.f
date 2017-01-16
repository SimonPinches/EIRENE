cdr: Aug. 2016  started commenting,documenting
cdr: Sept.  16  this routine is not used any longer. calls to dbl_poly now are
cdr             direct from xstcx, xstel and xstpi without this intermediate step


      subroutine EIRENE_prep_rtcs (ir, iflg, al, dum)
c  called from: xstcx,xstel,xstpi
c  prepare rate coefficients, originally given with double polynomial fit
c  (two independent parameters p1, p2)
c  1) find pointer to coefficients "rp" for reaction ir, from data structure reacdat(ir)
c  2) then evaluate this fit (call to dbl.poly.f) with log(p1)=al, log(p2)=0.0
c     and return the reduced (1D) fit coefficients DUM(1:9) for the 2nd parameter
c     dependency, evaluated at the fixed first parameter p1 

c  input:
c  ir   :  process number on data structure reacdat
c  iflg :  =3:  reaction rate coefficient                   (H.3, H.4)
c  iflg :  =4:  momentum weighted reaction rate coefficient (H.6, H.7)
c  iflg :  =5:  energy weighted reaction rate coefficient   (H.9, H.10)
c  iflg :  =6:  ???  not in use ???
c  al   :  log of 1'st parameter of double. polyn. fit, e.g. Te, Ti

c  output:
c  cf   :  fit parameters for p2 dependency, at fixed p1: al=log(p1)
c          i.e. for the fit=sum_1^9 dum(i) log(p2)^(i-1)
c          

 
      use EIRMOD_precision
      use EIRMOD_parmmod
      use EIRMOD_comxs
      use EIRMOD_comprt, only: iunout
 
      implicit none
 
      integer, intent(in) :: ir, iflg
      real(dp), intent(in) :: al
      real(dp), intent(out) :: dum(9)
      real(dp) :: cou
      real(dp) :: fp1(6),fp2(6)
      type(poly_data), pointer :: rp
      type(fit_forms), pointer :: rt
 
      select case (iflg)
      case (3)
        if (.not.reacdat(ir)%lrtc) then
          WRITE (IUNOUT,*) ' NO DATA AVAILABLE FOR RATE COEFFICIENT',ir
          CALL EIRENE_EXIT_OWN(1)
        END IF
        rt => reacdat(ir)%rtc
        rp => reacdat(ir)%rtc%poly
 
      case (4)
        if (.not.reacdat(ir)%lrtcmw) then
          WRITE (IUNOUT,*) ' NO DATA AVAILABLE FOR',
     .                     ' MOMENTUM-WEIGHTED RATE COEFFICIENT',ir
          CALL EIRENE_EXIT_OWN(1)
        END IF
        rt => reacdat(ir)%rtcmw
        rp => reacdat(ir)%rtcmw%poly
 
      case (5)
        if (.not.reacdat(ir)%lrtcew) then
          WRITE (IUNOUT,*) ' NO DATA AVAILABLE FOR',
     .                     ' ENERGY-WEIGHTED RATE COEFFICIENT',ir
          CALL EIRENE_EXIT_OWN(1)
        END IF
        rt => reacdat(ir)%rtcew
        rp => reacdat(ir)%rtcew%poly
 
      case (6)
        if (.not.reacdat(ir)%loth) then
          WRITE (IUNOUT,*) ' NO DATA AVAILABLE FOR',
     .                     ' OTHER REACTION',ir
          CALL EIRENE_EXIT_OWN(1)
        END IF
        rt => reacdat(ir)%oth
        rp => reacdat(ir)%oth%poly
 
      case default
        write (iunout,*)
     .  ' call EIRENE_to prep_rtcs with wrong flag, iflg = ',
     .                     iflg
        write (iunout,*) ' 3 <= iflg <= 6 assumed '
        call EIRENE_exit_own(1)
      end select
 
      fp1(1:3) = rt%fp1l
      fp1(4:6) = rt%fp1r
      fp2(1:3) = rt%fp2b
      fp2(4:6) = rt%fp2t
      call EIRENE_dbl_poly (rp%dblpol,al,0._dp,cou,dum,
     .     rt%rc1min, rt%rc1max, fp1, rt%jfex1mn, rt%jfex1mx,
     .     rt%rc2min, rt%rc2max, fp2, rt%jfex2mn, rt%jfex2mx)

      return
      end subroutine EIRENE_prep_rtcs
