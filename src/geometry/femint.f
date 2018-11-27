cdr order of points 3 and 4 in quadrangles (xpol,ypol) got changed in 2015. IFEM documents
cdr See IFEM documents
cdr june 17:  comments
cdr nov. 17: comments.  Some unfinished options re levgeo= 1 and levgeo = 2



c

c                     fem_interpolate  (this routine)
c  related routines:  fem_differentiate
c                     fem_local-coord
c                     fem_cell-corner
c
c  AFEM:  course "Advanced Finite Element Methods",
c         Department of Aerospace Engineering Sciences,
c         University of Colorado at Boulder
c         https://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/
c  IFEM:  course "Introduction to Finite Element Methods"
c         Department of Aerospace Engineering Sciences,
c         University of Colorado at Boulder
c         https://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/Home.html
c

      function eirene_femint (fecken, icell, x, y, z, lsame)
     .         result(res)

c  former function femint.f :  (-->  fem_interpolate.f)
c  interpolate a given function fecken, defined on cell vertices of grid cell no. icell,
c  using certain FEM-shape functions.

c  input:
c  lsame:   call with same coordinates x,y,z as in previous call, just another function 'fecken'
c           if lsame    : local coordinates r,s,t are taken from previous call
c           if not lsame: local coordinates are calculated here (call fem_local-coord)
c  output:
c  res : interpolated function fecken, evaluated at x,y,z

c  programmed for levgeo=4,5,
c  as well as partially for 2D (x,y) grids in case of levgeo=1,2,3


      use eirmod_precision
      use eirmod_parmmod
      use eirmod_clogau
      use eirmod_cgrid
      use eirmod_cgeom
      use eirmod_cpolyg
      use eirmod_ctrig
      use eirmod_ctetra
      USE EIRMOD_COMPRT, ONLY: IUNOUT


      implicit none

      real(dp), intent(in) :: fecken(*)
      real(dp), intent(in) :: x, y, z
      integer, intent(in) :: icell
      logical, intent(in) :: lsame

      real(dp) :: x1, x2, x3, x4, y1, y2, y3, y4, f1, f2, f3, f4,
     .            z1, z2, z3, z4, res
      real(dp), save :: r, s, t, u
      integer, save :: ir, ip, it, ia, ib

      res = 0._dp

c  2D computational grid, quadrangles
      if (((levgeo == 1) .and. nlrad .and. nlpol) .or.
     .    ((levgeo == 2 .and. .not.nlcrc) .and. nlpol) .or.
     .    (levgeo == 3)) then

        if (.not.lsame) then

c  new point x,y in icell. Find local coord. r,s

          call eirene_ncelln(icell,ir,ip,it,ia,ib,
     .                       nr1st,np2nd,nt3rd,nbmlt,
     .                       nlrad,nlpol,nltor)

          if (levgeo == 1) then
cdr   not ready for 2D x-z grids, nor for 3d x-y-z- grids
            x1=rsurf(ir)
            x2=rsurf(ir+1)
            x3=rsurf(ir+1)
            x4=rsurf(ir)
            y1=psurf(ip)
            y2=psurf(ip)
            y3=psurf(ip+1)
            y4=psurf(ip+1)

          else if ((levgeo == 2) .or. (levgeo == 3)) then
cdr  not ready: in case levgeo=2 and nlcrc:  xpol and ypol are not set.
cdr
c  warning, possibly a hidden link:
c  here we use xpol,ypol in case of levgeo=2.  These arrays may not be defined at all ??

            x1=xpol(ir+1,ip)
            x2=xpol(ir,ip)
            x3=xpol(ir+1,ip+1)
            x4=xpol(ir,ip+1)
            y1=ypol(ir+1,ip)
            y2=ypol(ir,ip)
            y3=ypol(ir+1,ip+1)
            y4=ypol(ir,ip+1)
          end if

          call eirene_xyz_to_rst(icell, x1, y1, 0._dp, x2, y2, 0._dp,
     .                           x3, y3, 0._dp, x4, y4, 0._dp,
     .                           x, y, z, r, s, t, u)
c       elseif (lsame)
c  same cell no. icell, x,y,z same as in previous call.
c  ir,ip, and local coordinates r,s are already set in previous call

        end if

        f1=fecken(INDPOINT(IR+1,IP))
        f2=fecken(INDPOINT(IR,IP))
        f3=fecken(INDPOINT(IR+1,IP+1))
        f4=fecken(INDPOINT(IR,IP+1))


        res = f1 * 0.25_dp * (1._dp - r) * (1._dp - s)
     .      + f2 * 0.25_dp * (1._dp + r) * (1._dp - s)
     .      + f3 * 0.25_dp * (1._dp + r) * (1._dp + s)
     .      + f4 * 0.25_dp * (1._dp - r) * (1._dp + s)

c...................................................................

c  2D computational grid, triangles
      else if (levgeo == 4) then

        if (.not.lsame) then
c  new point x,y in icell. Find local coord. r,s,t

          x1=xtrian(necke(1,icell))
          x2=xtrian(necke(2,icell))
          x3=xtrian(necke(3,icell))
          y1=ytrian(necke(1,icell))
          y2=ytrian(necke(2,icell))
          y3=ytrian(necke(3,icell))

          call eirene_xyz_to_rst(icell, x1, y1, 0._dp, x2, y2, 0._dp,
     .                           x3, y3, 0._dp, x4, y4, 0._dp,
     .                           x, y, z, r, s, t, u)
c       elseif (lsame)
c  same cell no icell, x,y as in previous call.
c  Local coordinates r,s,t  are already set in previous call

        end if

        f1=fecken(necke(1,icell))
        f2=fecken(necke(2,icell))
        f3=fecken(necke(3,icell))

        res = f1*r + f2*s + f3*t
c.........................................................................
c  3D computational grid, tetrahedra
      else if (levgeo == 5) then

        if (.not.lsame) then
c  new point x,y,z in icell. Find local coord. r,s,t,u

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

          call eirene_xyz_to_rst(icell, x1, y1, z1, x2, y2, z2,
     .                           x3, y3, z3, x4, y4, z4,
     .                           x, y, z, r, s, t, u)
c       elseif (lsame)
c  same cell no. icell, x,y,z, as in previous call.
c  Local coordinates r,s,t,u  are already set in previous call
        end if

        f1=fecken(nteck(1,icell))
        f2=fecken(nteck(2,icell))
        f3=fecken(nteck(3,icell))
        f4=fecken(nteck(4,icell))

        res = f1*r + f2*s + f3*t + f4*u

      else

        write (iunout,*) 'error in femint, un-written levgeo option'
        call eirene_exit_own(1)

      end if

      return
      end function eirene_femint
