module mod_nimbus_xs
  implicit none

  public :: nimxs_cxsig,nimxs_hmsig
  private

  !COMMON/DQHCOM/IDQH,KES,AR,TS,C1,C2,C3,C4
  real*8 :: xtfct_c1,xtfct_c2,xtfct_c3,xtfct_c4,xtfct_ts,xtfct_ar
  integer :: xtfct_idqh

  logical,save, public :: lflag_icx=.false.


  !C============================== DATA FOR MOLECULAR PROCESSES =========
  !C
  !C  .................................................................
  !C  .  MAXWELLIAN RATE COEFFICIENTS (E.M.JONES - CLM-R 175 (1977))  .
  !C  .................................................................
  !C
  !C     1)FAST DISSOCIATION:  H2-->H + H (TWO ATOMS AT 3EV)
  !C
  real*8, parameter,dimension(9) :: RH2E2H = &
      (/ -26.65989, &
      9.312341, &
      -5.315732, &
      2.551588, &
      -9.873073E-1, &
      2.490946E-1, &
      -3.733209E-2, &
      3.016028E-3, &
      -1.015095E-4/)
  !C
  !C     2)IONIZATION: H2-->2*H+ FOLLOWED BY -->(H+) + H (ONE ATOM AT 3EV)
  !C
  real*8, parameter, dimension(9) :: rh2= &
       (/ -34.23880, &
       16.14829, &
       -7.646064, &
       2.372854, &
       -0.4906446, &
       6.507808E-2, &
       -5.235834E-3, &
       2.302402E-4, &
       -4.222022E-6/)
  !C
  !C     3)DISSOCIATIVE IONIZATION: H2-->(H+) + H (ONE ATOM AT 3EV)
  !C
  real*8, parameter, dimension(9) :: rh2eh = &
       (/ -115.6823, &
       141.0502, &
       -95.33784, &
       37.05283, &
       -8.853047, &
       1.320344, &
       -1.198239E-1, &
       6.053489E-3, &
       -1.305276E-4 /)
  !C
  !C     4)SLOW DISSOCIATION: H2-->H + H  (TWO ATOMS AT 0.3EV)
  !C                                 (JANUARY 1988)
  real*8, parameter, dimension(7) :: rh2ehs = &
       (/ -34.906693, &
       14.250277, &
       -5.706695, &
       1.314436, &
       -1.776308E-1, &
       1.290029E-2, &
       -3.870771E-4/)
  !C
  !C  ...................................................................
  !C  .  MOLECULAR RATE COEFFICIENTS (JANEV - 'ELEMENTARY PROCESSES IN  .
  !C  .                                        HYDROGEN-HELIUM PLASMA', .
  !C  .                                        SPRINGER (1987))         .
  !C  .................................................................
  !C
  !C     1) DISSOCIATION:  H2 + E ---> H + H + E
  !C
  !C        REACTION 2.2.5 (P.44 WITH CORRECTIONS ACCORDING TO D.REITER)
  !C
  real*8, parameter, dimension(9) :: r225 = &
       (/-2.7872175E+01, &
       1.0522527E+01, &
       -4.9732123E+00, &
       1.4511982E+00, &
       -3.0627906E-01, &
       4.4333795E-02, &
       -4.0963442E-03, &
       2.1596703E-04, &
       -4.9285453E-06/)
  !C
  !C     2) IONIZATION: H2 ---> (H2+) + E + E
  !C                              |
  !C                              ---> H + (H+)
  !C
  !C        REACTION 2.2.9 (P.52)
  !C
  real*8, parameter, dimension(9) :: r229 = &
       (/-3.568640293666E+01, &
       1.733468989961E+01, &
       -7.767469363538E+00, &
       2.211579405415E+00, &
       -4.169840174384E-01, &
       5.088289820867E-02, &
       -3.832737518325E-03, &
       1.612863120371E-04, &
       -2.893391904431E-06/)
  !C
  !C     3) DISSOCIATIVE IONIZATION: H2 + E ---> H + (H+) + E + E
  !C
  !C        REACTION 2.2.10 (P.54)
  !C
  real*8, parameter, dimension(9) :: r2210 = &
      (/-3.834597006782E+01, &
      1.426322356722E+01, &
      -5.826468569506E+00, &
      1.727940947913E+00, &
      -3.598120866343E-01, &
      4.822199350494E-02, &
      -3.909402993006E-03, &
      1.738776657690E-04, &
      -3.252844486351E-06/)

  !C                    MICROSCOPIC CROSS SECTION OF REACTION
  !C                    (H+)+H-->H+(H+) IS GIVEN
  real*8, parameter, dimension(9) :: r318 = &
       (/ -3.274123792568E+01, &
       -8.916456579806E-02, &
       -3.016990732025E-02,&
       9.205482406462E-03,&
       2.400266568315E-03,&
       -1.927122311323E-03,&
       3.654750340106E-04,&
       -2.788866460622E-05,&
       7.422296363524E-07/)


contains

  subroutine nimxs_hmsig(iz,a,te,de,svmdf,svmds,svmi,sigmdf,sigmds,sigmi,isig,ixtype)
    implicit none
    real*8, intent(out) :: svmdf,svmi,svmds,sigmdf,sigmds,sigmi
    real*8, intent(in) :: a,te,de
    integer, intent(in) :: iz,isig,ixtype
    real*8 :: rc1,rc2,rc3,rc4,tran

    if(iz.ne.1) then
       write(*,'(a,i3)') ' error: nimxs_hmsig: non-hydrogenic species requested, called for z=',iz
       stop
    endif

    svmdf = 0.0
    svmi = 0.0
    svmds = 0.0
    sigmdf =0.0
    sigmds =0.0
    sigmi =0.0

    if(te.le.0.0 .or. de.le.0.0) return
    !c
    !c              rate coefficients  rc(cm**3/sec) v*sigma(v)
    !c                     rc1: dissociation
    !c                     rc2: ionization
    !c                     rc3: dissociative ionization
    !c                     rc4: slow dissociation
    !c
    if( ixtype.eq.1 ) then
       !c         janev
       call pnfit(9,r225(1) ,te,rc1)
       call pnfit(9,r229(1) ,te,rc2)
       call pnfit(9,r2210(1),te,rc3)
       rc4=0.0
    else
       !c         jones
       call pnfit(9,rh2e2h(1),te,rc1)
       call pnfit(9,rh2(1)   ,te,rc2)
       call pnfit(9,rh2eh(1) ,te,rc3)
       call pnfit(7,rh2ehs(1),te,rc4)
       !c           correct dissociation
       !c           --------------------
       !c           coefficient 0.3 was introduced on 4 february 1987
       !c           upon request of m.harrison, according to a suggestion of
       !c           janev and langer to heifetz (see heifetz letter 14/10/86
       !c           to de matteis about intor benchmark).
       !c           reference:usa contr. to the 14-th workshop meeting intor
       !c           phase ii dec.1986 brussels/vienna eurfubru/xii-52/86/edv 22
       !c           pag.154(heifetz)
       rc1=0.3*rc1
       rc4=0.3*rc4
    end if

    svmdf=rc1
    svmi= rc2+rc3
    svmds=rc4

    if(isig.eq.0) return
    !c                     from erg to ev
    tran=sqrt(9.57967e+11/a)
    !c        9.57967e+11 = 2.0*1.6022e-12/(2.0*a*1.6725e-24)
    !c                     pure dissociation (fast, 3 ev)
    sigmdf=rc1/tran
    !c                     pure dissociation (slow, 3 ev)
    sigmds=rc4/tran
    !c                     half dissociation and half ionization
    sigmi=(rc2+rc3)/tran
    !c              (to compute x-sections for h2-like molecules,
    !c               divide sigmdf,sigmds and sigmi by sqrt(e),
    !c               where e is the energy in ev of the molecule
    !c               composed of two atoms of mass a each)
    return
  end subroutine nimxs_hmsig


  subroutine nimxs_cxsig(iz,a,e,wx,wy,wz,ti,te,fmach,atarg,cdx,cdy,cdz,ixtype,sigcx,svcx)
    implicit none
    integer, intent(in) :: iz,ixtype
    real*8, intent(in) :: a,e,wx,wy,wz,ti,te,fmach,atarg,cdx,cdy,cdz
    real*8, intent(out) :: sigcx,svcx

    real*8 :: s,sx,sy,sz,r,rv2,es,ens,ests,rrate,er
    integer :: ip

    !C I   IZ       = ATOMIC NUMBER
    !C I   A        = MASS OF NEUTRAL (AMU)
    !C I   E        = NEUTRAL ENERGY (EV)
    !C I   WX/Y/Z   = NEUTRAL DIRECTION COSINES
    !C I   TI       = ION TEMPERATURE (EV)
    !C I   TE       = ELECTRON TEMPERATURE
    !C I   FMACH    = MACH NUMBER OF PLASMA FLOW
    !C I   ATARG    = ION MASS (AMU)
    !C I   CDX/Y/Z  = FLOW DIRECTION COSINES
    !C I   IXTYPE   = C.X. CROSS SECTION MODEL
    !C O   SIGCX    = C.X. CROSS SECTION (EFFECTIVE)
    !C O   SVCX     = C.X. REACTION RATE
    !C
    !C     RETURNS MICROSCOPIC CHARGE EXCHANGE X-SECTION  BETWEEN H-ISOTOPES
    !C     AVERAGED OVER ALL TARGET VELOCITIES (MAXWELL + FLOW) FOR NEUTRALS
    !C     IN A PLASMA REGION CONTAINING H+,D+,DT+,T+.
    !C
    IF(IZ.NE.1) then
       WRITE(*,'(a,i3)') ' ***** ERROR: NIMBXS NON-HYDROGENIC SPECIES REQUESTED, CALLED FOR Z=',iz
       stop
    endif



    !C                    VELOCITY OF THE NEUTRAL
    S=SQRT(E/A)*1.3841E+06
    SX=S*WX
    SY=S*WY
    SZ=S*WZ


    !C                       MEAN RELATIVE VELOCITY BETWEEN FLOW AND NEUTRAL
    IF(FMACH.LE.0.0) THEN
       RV2=S**2
    ELSE
       R=SQRT( (TE+TI)/ATARG )*1.0E+06 * FMACH
       !C                      0.978729E+06
       RV2=(R*CDX-SX)**2+(R*CDY-SY)**2+(R*CDZ-SZ)**2
    ENDIF
    !C                          SCALED SHIFTED NEUTRAL ENERGY (EV)
    ES=0.52197E-12*RV2
    !C                          SHIFTED NEUTRAL ENERGY
    ENS=ES*A
    !C                          SCALED ION TEMPERATURE
    xtfct_TS=TI/ATARG
    !C
    ESTS=ES/xtfct_TS

    ip=0
    if(es >= 40000.0) then
       if(ests < 100.0) ip=3
    elseif(es >=5000.) then
       if(ests < 26.)   ip=2
    else
       if(ests < 5.)    ip=1
    endif

    if(ip > 0) then
       rrate = nimxs_cxsig_integrate(atarg,a,ens,ti,ip)
       SIGCX=RRATE/S
    else
       !C                    MEAN RELATIVE NEUTRAL/ION ENERGY
       ER= 1.5*A*xtfct_TS+ENS
       CALL HXHP(ER/A,SIGCX)
       !C                    CONSERVE THE REACTION RATE
       SIGCX=SIGCX*SQRT(ER/E)
    endif

    !csw
    !sigcx=1.e-30
    !csw

    SVCX=SIGCX * S
    RETURN

  end subroutine nimxs_cxsig

  real*8 function nimxs_cxsig_integrate(atarg,a,ens,ti,ip) result(rrate)
    implicit none
    integer, intent(in) :: ip
    real*8, intent(in) :: atarg,a,ens,ti
    real*8 :: rm,c
!    real*8, external:: xtfct

    C=SQRT(ATARG/A)
    xtfct_C1=SQRT(ENS/TI)
    xtfct_C2=C*xtfct_C1
    xtfct_C3=0.0
    IF(xtfct_C2**2.LT.174.0) xtfct_C3=EXP(-xtfct_C2**2)
    xtfct_C4=2.0*xtfct_C2
    xtfct_IDQH=1

    select case(ip)
    case(1)
       CALL DQH04P(XTFCT,RM)
    case(2)
       CALL DQH12P(XTFCT,RM)
    case(3)
       CALL DQH32P(XTFCT,RM)
    end select
    RRATE=0.780939E+6*RM*TI/ATARG*SQRT(A/ENS)
  end function nimxs_cxsig_integrate

  SUBROUTINE HXHP(ER,SIGMA)
    implicit none
    real*8, intent(in) :: er
    real*8, intent(out) :: sigma
    !C
    !C     RETURNS MICROSCOPIC CROSS SECTION SIGMA FOR CHARGE-EXCHANGE
    !C     REACTION BETWEEN HYDROGEN ATOMS/IONS
    !C
    !C     ER=RELATIVE ENERGY
    !C
    !C     ICXTYP=1:JANEV - 'ELEMENTARY PROCESSES IN HYDROGEN-HELIUM PLASMA',
    !C                       SPRINGER (1987))
    real*8 :: escale
    ESCALE=ER
    IF(ESCALE.LT.0.1) ESCALE=0.1
    IF( ESCALE.LT.2.E6 )THEN
       CALL PNFIT(9,R318(1) ,ESCALE,SIGMA)
    ELSE
       SIGMA=0.0
    ENDIF
    RETURN
  END SUBROUTINE HXHP


  SUBROUTINE PNFIT(NJ,A,T,SVJ)
    implicit none
    !C
    !C     POLYNOMIAL FORM USED FOR MOLECULAR REACTIONS WITH ELECTRONS
    !C     AND JANEV C.X. X.S.
    !C
    !C     SVJ=RATE COEFFICIENT (CM**3/SEC)
    !C     T=TEMPERATURE (EV)
    !C     ALOG(SVJ)=A(1)+A(2)*ALOG(T)+...+A(NJ)*ALOG(T)**(NJ-1)
    !C              =A(1)+X*(A(2)+X*(A(3)+X*(...+X*(A(N-1)+A(NJ)*X)...)))
    !C
    integer, intent(in) :: nj
    real*8, intent(in) , dimension(nj) :: a
    real*8, intent(in) :: t
    real*8, intent(out) :: svj
    real*8 :: x,tt,sv
    integer :: i
    TT=T
    X=DLOG(TT)
    SV=A(NJ)
    DO I=NJ-1,1,-1
       SV=SV*X+A(I)
    enddo
    SVJ=EXP(SV)
    RETURN
  END SUBROUTINE PNFIT


  real*8 FUNCTION XTFCT(X) result(res)
    implicit none
    real*8, intent(in) :: x
    real*8 :: D,D1,D2,F,er,sigma

    ER=xtfct_TS*X*X
    select case(xtfct_idqh)
    case(1)
       CALL HXHP(ER,SIGMA)
    case(2)
       !CALL HXHEP(ER,SIGMA)
    case(3)
       !CALL XIPP(ER,SIGMA)
    case(4)
       !CALL ELXS(KES,ER*AR,SIGMA)
    end select


    if(xtfct_c3 /= 0.) then
       D=xtfct_C4*X
       IF(D.GT.174.0) then
          D1=xtfct_C2*(2.0D0*X-xtfct_C2)
          D2=-xtfct_C2*(2.0D0*X+xtfct_C2)
          F=DEXP(D1)-DEXP(D2)
       else
          D=DEXP(D)
          F=(D-1.0D0/D)*xtfct_C3
       endif
    else
       D1=xtfct_C2*(2.0D0*X-xtfct_C2)
       D2=-xtfct_C2*(2.0D0*X+xtfct_C2)
       F=DEXP(D1)-DEXP(D2)
    endif

    res=X*X*F*SIGMA
    RETURN
  END FUNCTION XTFCT



  SUBROUTINE DQH04P(FCT,Y)
    implicit none
    REAL*8 :: X,Y
    real*8, external :: FCT
    X=.29306374202572440D1
    Y=.19960407221136762D-3*FCT(X)
    X=.19816567566958429D1
    Y=Y+.17077983007413475D-1*FCT(X)
    X=.11571937124467802D1
    Y=Y+.20780232581489188D0*FCT(X)
    X=.38118699020732212D0
    Y=Y+.66114701255824129D0*FCT(X)
    RETURN
  END SUBROUTINE DQH04P

  SUBROUTINE DQH12P(FCT,Y)
    implicit none
    REAL*8 :: X,Y
    real*8, external :: FCT
    X=.60159255614257397D1
    Y=.16643684964891089D-15*FCT(X)
    X=.52593829276680444D1
    Y=Y+.65846202430781701D-12*FCT(X)
    X=.46256627564237873D1
    Y=Y+.30462542699875639D-9*FCT(X)
    X=.40536644024481495D1
    Y=Y+.40189711749414297D-7*FCT(X)
    X=.35200068130345247D1
    Y=Y+.21582457049023336D-5*FCT(X)
    X=.30125461375655648D1
    Y=Y+.56886916364043798D-4*FCT(X)
    X=.25238810170114270D1
    Y=Y+.8236924826884175D-3*FCT(X)
    X=.20490035736616989D1
    Y=Y+.70483558100726710D-2*FCT(X)
    X=.15842500109616941D1
    Y=Y+.37445470503230746D-1*FCT(X)
    X=.11267608176112451D1
    Y=Y+.12773962178455916D0*FCT(X)
    X=.67417110703721224D0
    Y=Y+.28617953534644302D0*FCT(X)
    X=.22441454747251559D0
    Y=Y+.42693116386869925D0*FCT(X)
    RETURN
  END SUBROUTINE DQH12P

  SUBROUTINE DQH32P(FCT,Y)
    implicit none
    REAL*8 ::  X,Y
    real*8, external :: FCT
    X=.10526123167960546D2
    Y=.55357065358569428D-48*FCT(X)
    X=.9895287586829539D1
    Y=Y+.16797479901081592D-42*FCT(X)
    X=.9373159549646721D1
    Y=Y+.34211380112557405D-38*FCT(X)
    X=.8907249099964770D1
    Y=Y+.15573906246297638D-34*FCT(X)
    X=.8477529083379863D1
    Y=Y+.25496608991129993D-31*FCT(X)
    X=.8073687285010225D1
    Y=Y+.19291035954649669D-28*FCT(X)
    X=.7689540164040497D1
    Y=Y+.7861797788925910D-26*FCT(X)
    X=.7321013032780949D1
    Y=Y+.19117068833006428D-23*FCT(X)
    X=.69652411205511075D1
    Y=Y+.29828627842798512D-21*FCT(X)
    X=.66201122626360274D1
    Y=Y+.31522545665037814D-19*FCT(X)
    X=.62840112287748282D1
    Y=Y+.23518847106758191D-17*FCT(X)
    X=.59556663267994860D1
    Y=Y+.12800933913224380D-15*FCT(X)
    X=.56340521643499721D1
    Y=Y+.52186237265908475D-14*FCT(X)
    X=.53183252246332709D1
    Y=Y+.16283407307097204D-12*FCT(X)
    X=.50077796021987682D1
    Y=Y+.39591777669477239D-11*FCT(X)
    X=.47018156474074998D1
    Y=Y+.7615217250145451D-10*FCT(X)
    X=.43999171682281376D1
    Y=Y+.11736167423215493D-8*FCT(X)
    X=.41016344745666567D1
    Y=Y+.14651253164761094D-7*FCT(X)
    X=.38065715139453605D1
    Y=Y+.14955329367272471D-6*FCT(X)
    X=.35143759357409062D1
    Y=Y+.12583402510311846D-5*FCT(X)
    X=.32247312919920357D1
    Y=Y+.8788499230850359D-5*FCT(X)
    X=.29373508230046218D1
    Y=Y+.51259291357862747D-4*FCT(X)
    X=.26519724354306350D1
    Y=Y+.25098369851306249D-3*FCT(X)
    X=.23683545886324014D1
    Y=Y+.10363290995075777D-2*FCT(X)
    X=.20862728798817620D1
    Y=Y+.36225869785344588D-2*FCT(X)
    X=.18055171714655449D1
    Y=Y+.10756040509879137D-1*FCT(X)
    X=.15258891402098637D1
    Y=Y+.27203128953688918D-1*FCT(X)
    X=.12472001569431179D1
    Y=Y+.58739981964099435D-1*FCT(X)
    X=.9692694230711780D0
    Y=Y+.10849834930618684D0*FCT(X)
    X=.69192230581004458D0
    Y=Y+.17168584234908370D0*FCT(X)
    X=.41498882412107868D0
    Y=Y+.23299478606267805D0*FCT(X)
    X=.13830224498700972D0
    Y=Y+.27137742494130398D0*FCT(X)
    RETURN
  END SUBROUTINE DQH32P

end module mod_nimbus_xs


