cdr: nov 2013: comments added
cdr  oct. 14 : formerly: routine df_xyz
cdr: nov 2015: further comments
cdr  oct.  18: warning: dfdz=0 in 2D (x,y) cases is only correct if
cdr            z coordinate is straght and orthogonal to X,Y.
cdr            This is not the case for NLTRA and NLTRT options.

      subroutine  eirene_fem_differentiate (fecken, icell, x, y, z,
     .                            dfdx, dfdy, dfdz)
c
c  AFEM:  course "Advanced Finate Element Methods", 
c         Department of Aerospace Enginerring Sciences, 
c         University of Colorado at Boulder
c         https://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/
c  IFEM:  course "Introduction to Finite Element Methods"
c         Department of Aerospace Enginerring Sciences, 
c         University of Colorado at Boulder
c         https://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/Home.html
c

c  return partial derivatives of function fecken, at internal point x,y,z, 
c         which is known to be located in grid cell icell

c  input:  fecken: values of function f on cell vertices
c          function fecken must be defined already, e.g. from an earlier call to 'cell-to-corner.f'
c
c  for speed-up, and overhead reduction:
c  fill array 'visited' to indicate, which cells have been visitied in earlier calls
c  currently: array 'visited' is only set for levgeo=4, 
C  to be done for levgeo=5
C  TO BE DONE: deallocate 'visited(icell)' at the end of a run.


c  not done here: check if x,y,z really inside cell icell ?
c  to be done:    range test for local coordinates r,s,t,u

      use eirmod_precision
      use eirmod_parmmod
      use eirmod_ctrig
      use eirmod_ctetra
      use eirmod_cgrid
      use eirmod_cgeom
      use eirmod_clogau
      use eirmod_ccona
      USE EIRMOD_COMPRT, ONLY: IUNOUT

      implicit none

      real(dp), intent(in) :: fecken(:)
      real(dp), intent(in) :: x, y, z
      real(dp), intent(out) :: dfdx, dfdy, dfdz
      integer, intent(in) :: icell

      real(dp) :: x1, x2, x3, x4, y1, y2, y3, y4, f1, f2, f3, f4, 
     .            z1, z2, z3, z4, det, eirene_deter4x4, 
     .            deti, r, s, t, u, dxdr, dxds, dydr, dyds,  
     .            drdx, dsdx, drdy, dsdy
      real(dp) :: a(4,4), ad(4,4), am1(4,4), 
     .            dndr(4), dnds(4), j(2,2), jm1(2,2)
      integer :: ir, ip, it, ia, ib

      logical, allocatable, save :: visited(:)
      real(dp), allocatable, save :: x32(:), x13(:), x21(:), 
     .                               y23(:), y31(:), y12(:), twoai(:)
      
cdr   integer, save :: icount=0  !  for test-output only
      
      real(dp) :: dummy

c  2d grid, quadrangels, x-y plane. z: ignorable.

c  levgeo=1,  2d regular orthogonal carthesian grid, x-y plane. z: ignorable.
c  levgeo=2:  quadrangels obtained from second order alg. surfaces,
c  levgeo=3:  general polygon grid, but convex cell.


      if (((levgeo == 1) .and. nlrad .and. nlpol) .or.
     .    ((levgeo == 2) .and. nlpol).or. 
     .     (levgeo == 3)) then

cdr in 2D this next statement is only correct for scalar tallies
cdr       or, for vector components, only if NLTRZ (ignorable z coordinate) is straight, cartesian
        dfdz = 0._dp

c  find 2D grid indices ir, ip, from 1d cell number icell

        call  eirene_ncelln(icell,ir,ip,it,ia,ib,nr1st,np2nd,nt3rd,
     .              nbmlt,nlrad,nlpol,nltor)

        if (levgeo == 1 .or. (levgeo == 2 .and. nlcrc)) then
cdr  levgeo=2 and nlcrc:  metric coefficients missing here ?
          x1=rsurf(ir)
          x2=rsurf(ir+1)
          x3=rsurf(ir+1)
          x4=rsurf(ir)
          y1=psurf(ip)
          y2=psurf(ip)
          y3=psurf(ip+1)
          y4=psurf(ip+1)
        else if ((levgeo == 2) .or. (levgeo == 3)) then
          x1=xpol(ir+1,ip)
          x2=xpol(ir,ip)
          x3=xpol(ir,ip+1)
          x4=xpol(ir+1,ip+1)
          y1=ypol(ir+1,ip)
          y2=ypol(ir,ip)
          y3=ypol(ir,ip+1)
          y4=ypol(ir+1,ip+1)
        end if

        f1=fecken(INDPOINT(IR+1,IP))
        f2=fecken(INDPOINT(IR,IP))
        f3=fecken(INDPOINT(IR,IP+1))
        f4=fecken(INDPOINT(IR+1,IP+1))

c  find local coordinates r,s,t,u of point x,y,z within cell icell

        call eirene_xyz_to_rst(icell, x1, y1, 0._dp, x2, y2, 0._dp,
     .                  x3, y3, 0._dp, x4, y4, 0._dp,
     .                  x, y, z, r, s, t, u)

c: to be done: what happens if, by error, x,y,z, are not in cell icell?
c: to be done: check for valid range of r,s,t,u ?

! partial derivatives of shape functions
        dndr(1) = -0.25_dp * (1._dp - s)
        dndr(2) =  0.25_dp * (1._dp - s)
        dndr(3) =  0.25_dp * (1._dp + s)
        dndr(4) = -0.25_dp * (1._dp + s)

        dnds(1) = -0.25_dp * (1._dp - r)
        dnds(2) = -0.25_dp * (1._dp + r)
        dnds(3) =  0.25_dp * (1._dp + r)
        dnds(4) =  0.25_dp * (1._dp - r)

        dxdr = x1*dndr(1) + x2*dndr(2) + x3*dndr(3) + x4*dndr(4)
        dxds = x1*dnds(1) + x2*dnds(2) + x3*dnds(3) + x4*dnds(4)
        dydr = y1*dndr(1) + y2*dndr(2) + y3*dndr(3) + y4*dndr(4)
        dyds = y1*dnds(1) + y2*dnds(2) + y3*dnds(3) + y4*dnds(4)

        j(1,1) = dxdr
        j(1,2) = dydr
        j(2,1) = dxds
        j(2,2) = dyds

        jm1(1,1) = j(2,2)
        jm1(1,2) = -j(1,2)
        jm1(2,1) = -j(2,1)
        jm1(2,2) = j(1,1)


        dummy = (j(1,1)*j(2,2) - j(1,2)*j(2,1))

        if ( abs(dummy) > eps30 ) then
           jm1 = jm1 / dummy
	else
	   jm1 = 0._DP
	endif

        drdx = jm1(1,1)
        dsdx = jm1(1,2)
        drdy = jm1(2,1)
        dsdy = jm1(2,2)

c  return partial derivaties wrt. cartesian coordinates, at point x,y,z

        dfdx =   f1 * (dndr(1)*drdx + dnds(1)*dsdx)
     .         + f2 * (dndr(2)*drdx + dnds(2)*dsdx)
     .         + f3 * (dndr(3)*drdx + dnds(3)*dsdx)
     .         + f4 * (dndr(4)*drdx + dnds(4)*dsdx)
        


        dfdy =   f1 * (dndr(1)*drdy + dnds(1)*dsdy)
     .         + f2 * (dndr(2)*drdy + dnds(2)*dsdy)
     .         + f3 * (dndr(3)*drdy + dnds(3)*dsdy)
     .         + f4 * (dndr(4)*drdy + dnds(4)*dsdy)

c  triangles, 2D grid
      elseif (levgeo == 4) then

c  precompute some parameters in cell icell, unless done so on earlier call
        if (.not.allocated(visited)) then
          allocate (x32(0:nrad))
          allocate (x13(0:nrad))
          allocate (x21(0:nrad))
          allocate (y23(0:nrad))
          allocate (y31(0:nrad))
          allocate (y12(0:nrad))

          allocate (twoai(0:nrad))
          allocate (visited(0:nrad))
          visited = .false.
        end if

        if (.not.visited(icell)) then
c  pre-compute some variables per cell (which are independent of point x,y,z)
          x1=xtrian(necke(1,icell)) 
          x2=xtrian(necke(2,icell)) 
          x3=xtrian(necke(3,icell)) 
          y1=ytrian(necke(1,icell)) 
          y2=ytrian(necke(2,icell)) 
          y3=ytrian(necke(3,icell))
          x32(icell)=x3-x2
          x13(icell)=x1-x3
          x21(icell)=x2-x1
          y23(icell)=y2-y3
          y31(icell)=y3-y1
          y12(icell)=y1-y2
          twoai(icell)=1._dp / 
     .         (x1*y23(icell) + x2*y31(icell) + x3*y12(icell))
        end if

        f1=fecken(necke(1,icell)) 
        f2=fecken(necke(2,icell)) 
        f3=fecken(necke(3,icell))
        dfdx = twoai(icell) *
     .         (f1*y23(icell) + f2*y31(icell) + f3*y12(icell))
        dfdy = twoai(icell) *
     .         (f1*x32(icell) + f2*x13(icell) + f3*x21(icell))

cdr in 2D this next statement is only correct for scalar tallies
cdr       or, for vector components, only if NLTRZ (ignorable z coordinate) is straight, cartesian
        dfdz = 0._dp

cdr     icount=icount+1
cdr     if (icount <= 1000) then
cdr        write (56,*) ' icount = ',icount, ' icell = ',icell
cdr        write (56,*) ' x, y, z ', x, y, z
cdr        write (56,*) ' dfdx,dfdy,dfdz ', dfdx, dfdy, dfdz 
cdr        write (56,*)
cdr     end if

        visited(icell) = .true.
        visited(0) = .false.    ! reset cell 0 for cell outside mesh

c  tetrahedra, 3d grid.

      else if (levgeo == 5) then

c  setting of array 'visited(icell)':  to be done

        x1=xtetra(nteck(1,icell))
        x2=xtetra(nteck(2,icell))
        x3=xtetra(nteck(3,icell))
        x4=xtetra(nteck(4,icell))
        y1=ytetra(nteck(1,icell))
        y2=ytetra(nteck(2,icell))
        y3=ytetra(nteck(3,icell))
        y4=ytetra(nteck(4,icell))
        z1=ztetra(nteck(1,icell))
        z2=ztetra(nteck(2,icell))
        z3=ztetra(nteck(3,icell))
        z4=ztetra(nteck(4,icell))
        f1=fecken(nteck(1,icell))
        f2=fecken(nteck(2,icell))
        f3=fecken(nteck(3,icell))
        f4=fecken(nteck(4,icell))

        a(1,1:4) = (/ 1._dp, 1._dp, 1._dp, 1._dp /)
        a(2,1:4) = (/ x1, x2, x3, x4 /)
        a(3,1:4) = (/ y1, y2, y3, y4 /)
        a(4,1:4) = (/ z1, z2, z3, z4 /)
        det = eirene_deter4x4(a)
        deti = 1._dp / det

        ad(1,1) = x2*(y3*z4-y4*z3)+x3*(y4*z2-y2*z4)+x4*(y2*z3-y3*z2)
        ad(1,2) = x1*(y4*z3-y3*z4)+x3*(y1*z4-y4*z1)+x4*(y3*z1-y1*z3)
        ad(1,3) = x1*(y2*z4-y4*z2)+x2*(y4*z1-y1*z4)+x4*(y1*z2-y2*z1)
        ad(1,4) = x1*(y3*z2-y2*z3)+x2*(y1*z3-y3*z1)+x3*(y2*z1-y1*z2)

        ad(2,1) = y2*(z4-z3) + y3*(z2-z4) + y4*(z3-z2)
        ad(2,2) = y1*(z3-z4) + y3*(z4-z1) + y4*(z1-z3)
        ad(2,3) = y1*(z4-z2) + y2*(z1-z4) + y4*(z2-z1)
        ad(2,4) = y1*(z2-z3) + y2*(z3-z1) + y4*(z1-z2)

        ad(3,1) = x2*(z3-z4) + x3*(z4-z2) + x4*(z2-z3)
        ad(3,2) = x1*(z4-z3) + x3*(z1-z4) + x4*(z3-z1)
        ad(3,3) = x1*(z2-z4) + x2*(z4-z1) + x4*(z1-z2)
        ad(3,4) = x1*(z3-z2) + x2*(z1-z3) + x3*(z2-z1)

        ad(4,1) = x2*(y4-y3) + x3*(y2-y4) + x4*(y3-y2)
        ad(4,2) = x1*(y3-y4) + x3*(y4-y1) + x4*(y1-y3)
        ad(4,3) = x1*(y4-y2) + x2*(y1-y4) + x4*(y2-y1)
        ad(4,4) = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)

        am1 = transpose(ad) * deti

        dfdx = f1*am1(1,2) + f2*am1(2,2) +
     .         f3*am1(3,2) + f4*am1(4,2)
        dfdy = f1*am1(1,3) + f2*am1(2,3) +
     .         f3*am1(3,3) + f4*am1(4,3)
        dfdz = f1*am1(1,4) + f2*am1(2,4) +
     .         f3*am1(3,4) + f4*am1(4,4)

      else
        write (iunout,*) ' levgeo = ',levgeo,' to be written in',
     .              ' subroutine derivative: fem_differentiate'
!pb         write (iunout,*) ' calculation abandonned '
!pb         call eirene_exit_own(1)
        dfdx = 0._dp
        dfdy = 0._dp
        dfdz = 0._dp
      end if

      return
      end subroutine eirene_fem_differentiate

      
