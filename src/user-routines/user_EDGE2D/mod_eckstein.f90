module mod_eckstein
  implicit none

  public :: rf0eck, rf1eck, eck_didimo, eck_recycle, eck_isect2d

  private
  real*8, save :: pia,pia2
  integer, save ::  rand=1,mulcor=13163521
  real*8, save :: table(24,11)
  integer, save :: itargtab(300)
  integer, save :: iptab(4,30), istab(30), ittab(30), iatab(30)
  real*8, save :: recmat(30,30)
  character(len=8) :: elemstr(30)
  real*8, save :: eatmd
  integer, save :: ntest,iscred

  logical, parameter :: ldebug=.false.
  integer, parameter :: ifeck = 598

  ! molecule emission?
  logical, save :: lmolec=.true.
  integer, save :: iaemis = 0

contains
  subroutine rf0eck
    implicit none
    integer :: i
    logical :: lex

    ! constants
    pia = 4.d0*atan(1.d0)
    pia2 = 2.d0*pia

    ! define table / iproj/itarg
    !C     (1)  SYMBOL,
    !C     (2)  CHARGE,
    !C     (3)  A.M.U.,
    !C     (4)  ATOMIC NUMBER Z,
    !C     (5)  THE ROW INDEX OF ITS ASSUMED FUNDAMENTAL ISOTOPE,
    !C     (6)  NEUTRAL SYMBOL,
    !C     (7)  MOLECULE DISSOCIATION TEMPERATURE,
    !C     (8)  MOLECULE SYMBOL,
    !C     (9)  AVERAGE IONIZATION CHARGE,
    !C     (10) DENSITY OF STRUCTURAL MATERIALS(GR/CM3).
    !C     (11) GAS(0) OR METAL(1)? (METALS STICK IN THE WALL)
    table=0.
    itargtab=0
    !C                              E-
    table(1,2) = -1.
    table(1,3) = 5.447e-4
    table(1,4) = 0.
    table(1,5) = 1.
    elemstr(1) = 'E-'
    !C                              H+
    table(2,2) = 1.
    table(2,3) = 1.
    table(2,4) = 1.
    table(2,5) = 2.
    table(2,7) = 3.0
    elemstr(2) = 'H(H+)'
    !C                              D+
    table(3,2) = 1.
    table(3,3) = 2.
    table(3,4) = 1.
    table(3,5) = 2.
    table(3,7) = 3.0
    elemstr(3) = 'D(D+)'
    !C                              DT+
    table(4,2) = 1.
    table(4,3) = 2.5
    table(4,4) = 1.
    table(4,5) = 2.
    table(4,7) = 3.0
    elemstr(4) = 'DT(DT+)'
    !C                              T+
    table(5,2) = 1.
    table(5,3) = 3.
    table(5,4) = 1.
    table(5,5) = 2.
    table(5,7) = 3.0
    elemstr(5) = 'T(T+)'
    !C                              HE+
    table(6,2) = 1.
    table(6,3) = 4.
    table(6,4) = 2.
    table(6,5) = 6.
    elemstr(6) = 'He(He+)'
    !C                              HE++
    table(7,2) = 2.
    table(7,3) = 4.
    table(7,4) = 2.
    table(7,5) = 6.
    elemstr(7) = 'He/(He++)'
    !C                              FE
    table(8,3) = 56.
    itargtab(56) = 8
    table(8,4) = 26.
    table(8,5) = 8.
    table(8,9) = 5.
    table(8,10)= 7.86
    table(8,11)= 1.
    elemstr(8) = 'Fe'
    !C                              CU
    table(9,3) = 63.5
    itargtab(63) = 9
    table(9,4) = 29.
    table(9,5) = 9.
    table(9,9) = 0.
    table(9,10)= 8.96
    table(9,11)= 1.
    elemstr(9) = 'Cu'
    !C                              C
    table(10,3) = 12.
    itargtab(12) = 10
    table(10,4) = 6.
    table(10,5) = 10.
    table(10,9) = 5.
    table(10,10)= 2.26
    table(10,11)= 1.
    elemstr(10) = 'C'
    !C                              MO
    table(11,3) = 96.
    itargtab(96) = 11
    table(11,4) = 42.
    table(11,5) = 11.
    table(11,9) = 0.
    table(11,10)= 10.2
    table(11,11)= 1.
    elemstr(11) = 'Mo'
    !C                              NI
    table(12,3) = 59.
    itargtab(59) = 12
    table(12,4) = 28.
    table(12,5) = 12.
    table(12,9) = 0.
    table(12,10)= 8.9
    table(12,11)= 1.
    elemstr(12) = 'Ni'
    !C                              W
    table(13,3) = 184.
    itargtab(184) = 13
    table(13,4) = 74.
    table(13,5) = 13.
    table(13,9) = 4.
    table(13,10)= 19.3
    table(13,11)= 1.
    elemstr(13) = 'W'
    !C                              AL
    table(14,3) = 27.
    itargtab(27) = 14
    table(14,4) = 13.
    table(14,5) = 14.
    table(14,9) = 0.
    table(14,10)= 2.7
    table(14,11)= 1.
    elemstr(14) = 'Al'
    !C                              AU
    table(15,3) = 197.
    itargtab(197) = 15
    table(15,4) = 79.
    table(15,5) = 15.
    table(15,9) = 0.
    table(15,10)= 19.3
    table(15,11)= 1.
    elemstr(15) = 'Au'
    !C                              BE
    table(16,3) = 9.
    itargtab(9) = 16
    table(16,4) = 4.
    table(16,5) = 16.
    table(16,9) = 0.
    table(16,10)= 1.85
    table(16,11)= 1.
    elemstr(16) = 'Be'
    !C                              SI
    table(17,3) = 28.
    itargtab(28) = 17
    table(17,4) = 14.
    table(17,5) = 17.
    table(17,9) = 0.
    table(17,10)= 2.33
    table(17,11)= 1.
    elemstr(17) = 'Si'
    !C                              TA
    table(18,3) = 181.
    itargtab(181) = 18
    table(18,4) = 73.
    table(18,5) = 18.
    table(18,9) = 0.
    table(18,10)= 16.6
    table(18,11)= 1.
    elemstr(18) = 'Ta'
    !C                              TI
    table(19,3) = 48.
    itargtab(48) = 19
    table(19,4) = 22.
    table(19,5) = 19.
    table(19,9) = 0.
    table(19,10)= 4.51
    table(19,11)= 1.
    elemstr(19) = 'Ta'
    !C                              V
    table(20,3) = 51.
    itargtab(51) = 20
    table(20,4) = 23.
    table(20,5) = 20.
    table(20,9) = 0.
    table(20,10)= 6.1
    table(20,11)= 1.
    elemstr(20) = 'V'
    !C                              ZR
    table(21,3) = 91.
    itargtab(91) = 21
    table(21,4) = 40.
    table(21,5) = 21.
    table(21,9) = 0.
    table(21,10)= 6.49
    table(21,11)= 1.
    elemstr(21) = 'Zr'
    !C                              NE
    table(22,3) = 20.2
    itargtab(20) = 22
    table(22,4) = 10.
    table(22,5) = 22.
    table(22,9) = 0.
    table(22,10)= 0.
    table(22,11)= 0.
    elemstr(22) = 'Ne'
    !C                              HE IMPURITY
    table(23,3) = 4.
    table(23,4) = 2.
    table(23,5) = 6.
    elemstr(23) = 'HeImp'
    !C                              N
    table(24,3) = 14.
    table(24,4) = 7.
    table(24,5) = 24.
    elemstr(24) = 'N'


    ! init random generator
    rand=666

    ! energy of particles re-emitted as atoms
    ! <0. : use dissociation energy from table(:,7)
    eatmd = 0.025d0
    eatmd = 0.03d0

    ! ntest test-particles
    ntest=3

    !
    ! particle tables EIRENE --> NIMBUS
    iptab=0

    ! D
    iptab(1,1) = 3
    ! D2
    iptab(2,1) = 3
    ! D+
    iptab(4,1) = 3

    ! C
    iptab(1,2) = 10
    iptab(4,2) = 10
    iptab(4,3) = 10
    iptab(4,4) = 10
    iptab(4,5) = 10
    iptab(4,6) = 10
    iptab(4,7) = 10

    !
    ! particle tables NIMBUS(ncom) --> EIRENE
    istab=0
    ittab=0
    iatab=0

    ! D
    istab(1)   = 1
    ittab(1)   = 1
    iatab(3)   = istab(1)

    ! D2
    istab(2)   = 1
    ittab(2)   = 2

    ! carbon
    istab(3)   = 2
    ittab(3)   = 1
    iatab(10)  = istab(3)


    open(unit=4999, file='iaemis.tag',form='formatted')
    read(4999,'(i5)') iaemis
    close(4999)
    if(iaemis .eq. 1) lmolec=.false.

    if(ldebug) then
       write(ifeck,*) 'RF0ECK: LMOLEC=',lmolec
    endif

    ! define default recmat for all surfaces
    recmat = 0.
    if(lmolec) then
       recmat(2,3) = 1.
    else
       recmat(1,3) = 1.
    endif
    ! carbon impurity, sticky or not?
    recmat(3,10) = 1.


  end subroutine rf0eck

  subroutine rf1eck(xmw,xcw,xmp,xcp,igasf,igast,zcos,zsin,expi,rprob,e0term, &
       velx,vely,velz, e0, ityp,isp, weight, crtx,crty,crtz)
    implicit none
    real*8, intent(in) :: xmw,xcw,xmp,xcp,zcos,zsin,expi,rprob,e0term,crtx,crty,crtz
    integer, intent(in) :: igasf,igast
    real*8, intent(inout) :: e0,velx,vely,velz
    integer, intent(inout) :: ityp,isp
    real*8, intent(inout) :: weight

    real*8 :: a1,a2,z1,z2,polar,cospol,den,epskev,aloge,e,coskt
    real*8 :: sinp,cosp,sint,cost,sinpn,cospn,sintn,costn,rn,re,ewscal
    real*8 :: c,totp,t,vn,cx,cy,cz,cc
    integer :: ia,iproj,itarg,ip,i,in,is,ity


    !e0=.1

    a1 = xmp
    z1 = xcp
    a2 = xmw
    z2 = xcw
    e  = e0

    if(ityp == 2) then
       e =e *0.5
       a1=a1*0.5
    endif

    ia = int(a1)
    select case(ia)
       case(1)
          ! H/H+
          iproj = 2
       case(2)
          ! D/D+
          iproj = 3
       case(5)
          ! DT/DT+
          iproj = 4
          a1=a1/2.
       case(3)
          ! T/T+
          iproj = 5
       case(4)
          ! He+/He
          ! check charge! --> iproj 6 or 7 ?
          stop
       case(12)
          ! carbon
          iproj = 10
       case default
          write(*,*) 'mod_eckstein, a1 not recognised...'
          stop
    end select

    ia = int(a2)
    itarg = itargtab(ia)

    if(ldebug) then
       write(ifeck,'(a)') '------------------------------------'
       write(ifeck,'(4a)') 'RF1ECK: ',elemstr(iproj),' on ',elemstr(itarg)
       write(ifeck,'(a)') 'before:'
       write(ifeck,'(a,2(1x,i5))')  'ityp,isp = ',ityp,isp
       write(ifeck,'(a,4(1x,e14.7))') 'e0,e:  ',e0,e
       write(ifeck,'(a,4(1x,e14.7))') 'vel: ',velx,vely,velz,velx**2+vely**2+velz**2
    endif

    ! new direction
    cospn = -crtx
    sinpn = -crty
    !costn = -crtz
    !sintn = sqrt(1.0-costn**2)
    costn = 0.
    sintn = 1.

    ! old direction
    cost = velz
    sint = sqrt(1.0-cost**2)
    if(abs(sint) > 0.) then
       cosp = velx/sint
       sinp = vely/sint
    else
       cosp=0.
       sinp=0.
    endif

    cospol = -(( cosp*cospn + sinp*sinpn)* sint*sintn+ cost*costn)
    if(ldebug) then
       write(ifeck,*) ' COSPOL=',cospol
    endif

    !if(cospol < -.999) cospol = -.999
    !if(cospol > +.999) cospol = +.999

    den = sqrt(z1**0.6666+z2**0.6666)*z1*z2*(a1+a2)
    epskev = e*(0.0325*a2/den)
    aloge = log10(epskev)
    polar=abs(asin(sqrt(1.0-cospol**2)))

    ! re-deposition? sticking? for test-particles
    !iscred=0 ! NO RE-DEP YET
    !if(table(iproj,11) /= 0. .and. ityp /= 4) then
    !   if(iscred == 0) then
    !      weight=0.
    !      return
    !   endif
    !endif

    if(iproj <= 7) then
       !
       ! H/D/T/DT/He projectile
       !
       if(itarg == 8 .or. itarg == 12) then
          !
          ! Fe/Ni wall
          !
          if (e > 1.5 .or. iproj > 5) then
             ! rn = ((1.0+3.2116*epskev**0.34334)**1.5+(1.3288*epskev**1.5)**1.5)**(-0.6666667)
             ! following formula equivalent to previous up to epskev=4 (He or heavier nuclei,
             ! any energy):
             rn = 0.1926-0.2037*aloge
          else
             ! low energy scaling rules (Eckstein 1984)
             select case(iproj)
                case(2)
                   ewscal = e
                case(3)
                   ewscal = e-0.10
                case(4)
                   ewscal = e-0.15
                case(5)
                   ewscal = e-0.20
             end select
             if(ewscal < 0.) ewscal = 0.001
             rn = 1.125*ewscal-0.225*ewscal**2.78-0.2
          endif

          if(e>1.5) then
             ! re = ((1.0+7.1172*epskev**0.3525)**1.5+(5.2757*epskev**1.5)**1.5)**(-0.666667)
             ! the following formular is equivalent to the previous up to epskev=4
             ! (nuclei heavier than He, any energy)
             re = 0.10-0.1667*aloge
          elseif(iproj <=5) then
             ! low energy, H/D/DT/T
             re = 0.6*rn
          else
             ! low energy, He
             re = 0.6
          endif


       elseif(itarg == 10) then
          !
          ! C wall
          !
          if(e>3.0) then
             ! high energy, Eckstein-Verbeek
             if(e <= 1000.0) then
                select case(iproj)
                case(1,2)
                   rn=e**(-0.3)
                   re=0.78*e**(-0.417)
                case(3)
                   rn=0.758*e**(-0.3)
                   re=0.522*e**(-0.417)
                case(4)
                   rn=0.578*e**(-0.3)
                   re=0.392*e**(-0.417)
                case(5,6,7)
                   rn=0.459*e**(-0.3)
                   re=0.313*e**(-0.417)
                end select
             else
                select case(iproj)
                case(1,2)
                   rn=575.4*e**(-1.22)
                   re=1178.0*e**(-1.476)
                case(3)
                   rn=411.0*e**(-1.22)
                   re=743.0*e**(-1.476)
                case(4)
                   rn=320.0*e**(-1.22)
                   re=589.0*e**(-1.476)
                case(5,6,7)
                   rn=274.0*e**(-1.22)
                   re=473.0*e**(-1.476)
                end select
             endif
          else
             ! low energy scaling rules for H (Eckstein 1984)
             ! D/T missing
             rn = 0.348*e-0.00146*e**4.777-0.156
             if(rn < 0.) then
                rn=0.001
             endif
             if(rn > 1.) rn=1.0
             re = 0.3*rn
          endif

       elseif(itarg == 13) then
          !
          ! W wall
          !
          if(e > 1.5) then
             ! high energy
             rn = 0.1885-0.2265*aloge
             if(epskev > 0.1) then
                re = 0.07-0.18*aloge
             else
                re = -0.25*aloge
             endif
          else
             ! low energy scaling rules for H (Eckstein 1984)
             ! D/T missing
             rn = 0.85-0.332*(1.5-e)**3
             re = 0.6*rn
          endif

       elseif(itarg == 13) then
          !
          ! Au wall
          !
          if(e > 0.0) then
             ! high energy
             rn = 0.306-0.1318*aloge
             re = 0.061-0.138*aloge
          else
             rn=1.0
             re=1.0
          endif

       else
          !
          ! all other walls
          !
          rn = 0.1885-0.2265*aloge
          if(epskev > 0.1) then
             re = 0.07-0.18*aloge
          else
             re = -0.25*aloge
          endif

          ! special check for H on Cu
          ! Oen/Robinson, Nucl.Instr and Method 132 pp647 (1976)
          if(iproj ==2 .and. itarg== 9) then
             if(e <= 2000.) then
                rn=rn*(0.597*polar**2-0.007*polar+1.0)
                re=re*(0.405*polar**2+0.446*polar+1.0)
             else
                rn=rn*(1.807*polar**2-0.250*polar+1.0)
                re=re*(5.320*polar**2-1.531*polar+1.0)
             endif
          endif

       endif


    else
       !
       ! all other projectiles
       !
       if(itarg == 8 .or. itarg == 12) then
          ! Fe/Ni wall
          ! rn = ((1.0+3.2116*epskev**0.34334)**1.5+(1.3288*epskev**1.5)**1.5)**(-0.6666667)
          ! following formula equivalent to previous up to epskev=4 (He or heavier nuclei,
          ! any energy):
          rn = 0.1926-0.2037*aloge

          ! re = ((1.0+7.1172*epskev**0.3525)**1.5+(5.2757*epskev**1.5)**1.5)**(-0.666667)
          ! the following formular is equivalent to the previous up to epskev=4
          ! (nuclei heavier than He, any energy)
          re = 0.10-0.1667*aloge
       else
          ! all other walls
          rn = 0.1885-0.2265*aloge

          if(epskev > 0.1) then
             re = 0.07-0.18*aloge
          else
             re = -0.25*aloge
          endif
       endif

    endif

    if(rn < 0.) rn=0.001
    if(rn > 1.) rn=1.0
    if(re < 0.) re=0.001
    if(re > 1.) re=1.0

    ! backscattering/reflection or re-emission ?
    c = myranf()

    if(ityp == 2) then
       ! molecules?, then force re-emission (RECYCF == 0)
       c=1.1
    endif

    !c=0.

    if(c > rn ) then
       !
       ! re-emission
       ! particle type may change
       !
       if(ldebug) WRITE(ifeck,*) 'RE-EMISSION'

       ! sample new particle type from recmat
       ip = iptab(ityp,isp)
       ! cumulative
       totp = 0.
       do i=1,ntest
          totp=totp+recmat(i,ip)
       enddo

       !
       if(totp <= 0.) then
          weight = 0.
          return
       endif

       ! sample
       t=0.
       c = myranf()
       do i=1,ntest
          in = i
          t=t+recmat(i,ip)/totp
          if(c <= t) exit
       enddo
       is  = istab(in)
       ity = ittab(in)

       !if(ity == 2 .and. ityp/=2) then
       !   weight=weight*0.5
       !endif

       ! new energy
       select case(ity)
       case(1)
          !atoms
          if(iproj <= 7) then
             ! H/D/DT/He projectiles on standard wall
             if(eatmd < 0.) then
                ! dissoc. energy
                e = table(iproj,7)
             else
                ! prescribed energy
                e = eatmd
             endif
             e=e
          else
             ! all other projectiles
             e=0.
          endif

          if(e <= 0.) then
             ! if not specified, new energy=wall-temperature
             e=e0term
          endif
       case(2)
          !molecules, thermal energy according to wall temp
          !e=1.5*e0term
          ! here, factor 1.5 already included in inputfile
          e=e0term
       case default
          write(*,*) 'mod_eckstein: error, ityp=',ityp,' not available in re-emission'
          stop
       end select

    else
       !
       ! backscattering/reflection
       !
       if(ldebug) WRITE(ifeck,*) 'BACKSCATTERING'

       ! new energy
       e=e*re/rn

       ! same particle type
       ity = ityp
       is  = isp

       ! hack, force atoms (no molecules) when plasma-ion hits surface
       if(ityp == 4) then
          ip  = iptab(ityp,isp)
          ity = 1
          is  = iatab(ip)
       endif

       if(ity == 2) e=e*2.0

    endif

    ! sample new direction, cosine law for polar axis
    call coslaw(coskt,2)
    cosp=cospn
    sinp=sinpn
    cost=costn !0
    sint=sintn !1
    ! turn polar axis & sample azimuthal angle
    call eck_didimo(coskt,cosp,sinp,cost,sint)

    ! cartesian
    velx = cosp*sint
    vely = sinp*sint
    velz = cost
    !velz=0.

    vn = sqrt(velx**2+vely**2+velz**2)
    velx=velx/vn
    vely=vely/vn
    velz=velz/vn

    e0   = e
    ityp = ity
    isp  = is

    if(ldebug) then
       write(ifeck,'(a)') 'after:'
       write(ifeck,'(a,2(1x,i5))')  'ityp,isp = ',ityp,isp
       write(ifeck,'(a,4(1x,e14.7))') 'e0:  ',e0
       write(ifeck,'(a,4(1x,e14.7))') 'vel: ',velx,vely,velz,velx**2+vely**2+velz**2
    endif

  end subroutine rf1eck


  subroutine coslaw(coskt,mod)
    !c     source cosine with respect to polar axis between 0 and 90 deg.
    !c     mod=1  uniform distribution
    !c     mod=2  cosine distribution ( g(coskt)=2*coskt , 0<=coskt<=1 )
    implicit none
    integer, intent(in) :: mod
    real*8, intent(out) :: coskt
    real*8 :: c
    c=myranf()
    select case(mod)
    case(1)
       coskt=c
    case(2)
       coskt=sqrt(c)
    end select
  end subroutine coslaw

  subroutine eck_didimo(coskt, cosp,sinp,cost,sint)
    !c                given a polar axis through cosp,sinp,cost
    !c                and the cosine coskt of a cone around this axis,
    !c                s.r.didimo returns a direction laying on the cone,
    !c                with azimuthal angle randomly uniform (old polar axis
    !c                direction is lost, new direction is given through
    !c                cosp,sinp,cost)
    !c                (didu,didv,didw)=(cosp*sint,sinp*sint,cost)=dir.cosines
    implicit none
    real*8, intent(in) :: coskt
    real*8, intent(inout) :: cosp,sinp,cost,sint
    real*8 :: c,sinkt,sinkp,coskp,didu1,didv1,didw1
    real*8 :: didbc,didbcw,didbd,didu,didv,temp

    sinkt=sqrt(1.0d0-coskt**2)
    c=myranf()
    coskp=cos(pia2*c)

    sinkp=sqrt(1.0d0-coskp**2)
    c=myranf()
    if(c-0.5d0 < 0.d0) then
       sinkp=-sinkp
    endif

    if(sint-.1e-4  <= 0.) then
       didu1=sinkt*coskp
       didv1=sinkt*sinkp
       didw1=coskt*cost
    else
       didbc=sinkt*coskp
       didbcw=didbc*cost
       didbd=sinkt*sinkp
       didu=sint*cosp
       didv=sint*sinp
       didu1=(didbcw*didu-didbd*didv)/sint+coskt*didu
       didv1=(didbcw*didv+didbd*didu)/sint+coskt*didv
       didw1=-didbc*sint+coskt*cost
    endif

    cost=didw1


    !if(abs(cost) >= 0.999d0) then
    !   cost=sign(0.999d0,cost)
    !elseif(abs(cost) <= 0.001) then
    !   cost=sign(0.001d0,cost)
    !endif

    sint=sqrt(1.0d0-cost**2)
    !tant=sint/cost
    if(abs(sint) > 0.) then
       cosp=didu1/sint
       sinp=didv1/sint

       temp=sqrt(sinp**2+cosp**2)
       sinp=sinp/temp
       cosp=cosp/temp
    else
       cosp=0.
       sinp=0.
    endif



    !if(abs(cosp) >= 0.999d0) then
    !   cosp=sign(0.999d0,cosp)
    !elseif(abs(cosp) <= 0.001d0) then
    !   cosp=sign(0.001d0,cosp)
    !endif
    !sinp=sign(sqrt(1.0d0-cosp**2),sinp)

    return
  end subroutine eck_didimo


  subroutine eck_recycle(eshet,tiwd,tewd,nmass,ee)
    implicit none
    real*8, intent(in) :: eshet,tiwd,tewd
    integer, intent(in) :: nmass
    real*8, intent(out) :: ee
    real*8 :: e,et,avsou2
    real*8, parameter :: coefmx=1.66666d0

    call sqxexp(e)
    e=e*coefmx
    e=e*tiwd

    et=e

    avsou2 = (1.d0*tewd+tiwd)/2.d0
    e=et+.5*nmass*avsou2

    e=e+eshet


    ee=e
  end subroutine eck_recycle

  subroutine eck_isect2d(x0,y0,vx,vy,x1,y1,x2,y2, flag)
    implicit none
    real*8, intent(in) :: vx,vy,x1,y1,x2,y2
    real*8, intent(inout) :: x0,y0
    logical, intent(out) :: flag
    real*8 :: rx,ry,ax,ay,v,t

    flag=.false.
    rx = x2-x1
    ry = y2-y1
    ax = x1 - x0
    ay = y1 - y0
    v  = (ax*vy - ay*vx) / (ry*vx-rx*vy+1.d-60)
    if(v >= 0. .and. v <= 1.) then
       if(abs(vx) > abs(vy)) then
          t = (ax+v*rx)/vx
       else
          t = (ay+v*ry)/vy
       endif

       if(t > 0.) then
          x0 = x0+t*vx
          y0 = y0+t*vy
          flag=.true.
       endif
    endif
    return
  end subroutine eck_isect2d

  SUBROUTINE  SQXEXP(xX)
    !C----- SCEGLIE Xx DA SQRT(Xx)*EXP(-Xx)    (N(E) MAXWELLIANA)
    !C----- METODO MISTO INTERVALLI EQUIPROBABILI - REIEZIONE
    !C     SI TRASCURA LA CODA X>10 (ERRORE RELATIVO 0.0003)
    !C
    real*8, intent(out) :: xx
    real*4, parameter :: fntab=128.0
    real*4, dimension(129), parameter :: tab = &
         (/ 0.0       ,4.8526E-02,7.7934E-02,1.0315E-01,1.2609E-01,1.4756E-01,&
         1.6797E-01,1.8758E-01,2.0658E-01,2.2507E-01,2.4315E-01,2.6090E-01,&
         2.7835E-01,2.9558E-01,3.1259E-01,3.2943E-01,3.4613E-01,3.6271E-01,&
         3.7918E-01,3.9556E-01,4.1188E-01,4.2813E-01,4.4435E-01,4.6053E-01,&
         4.7668E-01,4.9283E-01,5.0897E-01,5.2512E-01,5.4127E-01,5.5745E-01,&
         5.7365E-01,5.8990E-01,6.0617E-01,6.2250E-01,6.3888E-01,6.5531E-01,&
         6.7182E-01,6.8838E-01,7.0503E-01,7.2176E-01,7.3857E-01,7.5547E-01,&
         7.7248E-01,7.8958E-01,8.0680E-01,8.2412E-01,8.4157E-01,8.5914E-01,&
         8.7684E-01,8.9467E-01,9.1265E-01,9.3078E-01,9.4904E-01,9.6748E-01,&
         9.8609E-01,1.0049E+00,1.0238E+00,1.0430E+00,1.0623E+00,1.0818E+00,&
         1.1016E+00,1.1215E+00,1.1417E+00,1.1621E+00,1.1827E+00,1.2036E+00,&
         1.2248E+00,1.2462E+00,1.2679E+00,1.2899E+00,1.3122E+00,1.3348E+00,&
         1.3578E+00,1.3810E+00,1.4046E+00,1.4286E+00,1.4530E+00,1.4777E+00,&
         1.5029E+00,1.5285E+00,1.5546E+00,1.5810E+00,1.6081E+00,1.6356E+00,&
         1.6636E+00,1.6923E+00,1.7214E+00,1.7513E+00,1.7817E+00,1.8129E+00,&
         1.8447E+00,1.8774E+00,1.9108E+00,1.9451E+00,1.9803E+00,2.0164E+00,&
         2.0535E+00,2.0917E+00,2.1310E+00,2.1716E+00,2.2136E+00,2.2568E+00,&
         2.3016E+00,2.3480E+00,2.3962E+00,2.4463E+00,2.4984E+00,2.5529E+00,&
         2.6098E+00,2.6695E+00,2.7322E+00,2.7983E+00,2.8682E+00,2.9424E+00,&
         3.0215E+00,3.1062E+00,3.1973E+00,3.2961E+00,3.4038E+00,3.5225E+00,&
         3.6546E+00,3.8037E+00,3.9754E+00,4.1769E+00,4.4228E+00,4.7374E+00,&
         5.1755E+00,5.9152E+00,1.0000E+01 /)


    real*4, dimension(128), parameter :: alt = &
         (/ 2.0985E-01,2.5824E-01,2.8969E-01,3.1303E-01,3.3144E-01,3.4647E-01,&
         3.5903E-01,3.6968E-01,3.7880E-01,3.8667E-01,3.9349E-01,3.9940E-01,&
         4.0455E-01,4.0901E-01,4.1287E-01,4.1620E-01,4.1904E-01,4.2145E-01,&
         4.2346E-01,4.2512E-01,4.2644E-01,4.2745E-01,4.2818E-01,4.2864E-01,&
         4.2886E-01,4.2888E-01,4.2885E-01,4.2862E-01,4.2819E-01,4.2757E-01,&
         4.2677E-01,4.2579E-01,4.2466E-01,4.2337E-01,4.2194E-01,4.2036E-01,&
         4.1866E-01,4.1683E-01,4.1487E-01,4.1280E-01,4.1062E-01,4.0833E-01,&
         4.0594E-01,4.0345E-01,4.0086E-01,3.9818E-01,3.9542E-01,3.9257E-01,&
         3.8963E-01,3.8662E-01,3.8352E-01,3.8036E-01,3.7712E-01,3.7381E-01,&
         3.7043E-01,3.6698E-01,3.6347E-01,3.5990E-01,3.5627E-01,3.5257E-01,&
         3.4882E-01,3.4501E-01,3.4115E-01,3.3723E-01,3.3326E-01,3.2924E-01,&
         3.2517E-01,3.2105E-01,3.1688E-01,3.1266E-01,3.0840E-01,3.0409E-01,&
         2.9974E-01,2.9535E-01,2.9091E-01,2.8643E-01,2.8191E-01,2.7735E-01,&
         2.7275E-01,2.6811E-01,2.6343E-01,2.5872E-01,2.5397E-01,2.4918E-01,&
         2.4435E-01,2.3949E-01,2.3461E-01,2.2967E-01,2.2471E-01,2.1972E-01,&
         2.1469E-01,2.0962E-01,2.0453E-01,1.9940E-01,1.9424E-01,1.8906E-01,&
         1.8383E-01,1.7858E-01,1.7330E-01,1.6798E-01,1.6263E-01,1.5726E-01,&
         1.5186E-01,1.4643E-01,1.4096E-01,1.3547E-01,1.2995E-01,1.2440E-01,&
         1.1882E-01,1.1321E-01,1.0756E-01,1.0190E-01,9.6196E-02,9.0464E-02,&
         8.4704E-02,7.8904E-02,7.3083E-02,6.7221E-02,6.1341E-02,5.5414E-02,&
         4.9461E-02,4.3470E-02,3.7429E-02,3.1364E-02,2.5237E-02,1.9070E-02,&
         1.2862E-02,6.5623E-03/)

    integer, parameter :: modo = 2

    real*4 :: c,fint,x,h,xmin,delta,c1,f,y
    integer :: k

    C=sngl(myRANF())
    select case(modo)
    case(1)
       FINT=C*FNTAB
       K=FINT
       FINT=FINT-FLOAT(K)
       K=K+1
       X=TAB(K)+(TAB(K+1)-TAB(K))*FINT
    case(2)
       !C LINEAR INTERPOLATION
       !C INTERPOLATION AND REJECTION
       K=C*FNTAB
       K=K+1
       H=ALT(K)
       XMIN=TAB(K)
       DELTA=TAB(K+1)-XMIN
       do
          C=sngl(myRANF())
          C1=sngl(myRANF())
          X=XMIN+C*DELTA
          Y=C1*H
          F=SQRT(X)*EXP(-X)
          if( (y-f) <= 0.) exit
       enddo
    end select
    xx = dble(x)
  END SUBROUTINE SQXEXP

  real*8 function myranf() result(res)
    implicit none
    real*4 :: rand
    real*8, external :: EIR_ranf_eirene
    do
       !call myrandom( rand )
       !res = dble(rand)
       res= EIR_ranf_eirene()
       if(res> 0.d0 .and. res < 1.d0) return
    enddo
    return
  end function myranf

  subroutine myrandom(fl)
    implicit none
    real*4, intent(out) :: fl

    integer, parameter :: maxint = 2147483647
    integer, parameter :: modulo =  268435456
    integer, parameter :: moltip =   41475557
    integer, parameter :: period =   67108864

    rand = rand*moltip
    rand = ibclr(rand,31)
    rand = mod(rand,modulo)
    fl=rand
    fl=fl*3.72528e-09
    return
  end subroutine myrandom


end module mod_eckstein
