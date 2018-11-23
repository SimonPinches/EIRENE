!pb  070109  ff(2,2) introduced to avoid array temporaries

cdr  oct.15  comments. Parameters ip1,ip2 added, to indicate extrapolation
cdr          tbd: generalize to general 2d (2 indep., parameters) tabular data
cdr               not just adas data.

      function EIRENE_intp_tab2d(ad,p1,p2,ip1,ip2) result(res)
c  bilinear interpolation in 2-parameter table ad
c  linear extrapolation at all 4 parameter boundaries

c input:
c  table: ad(..,..)
c  parameters:  p1, p2

c output:
c  res,  interpolated value at (p1,p2)
c  ip1 = 0:  parameter p1 is in valid range of table
c  ip1 =-1:  parameter p1 is below left boundary of tabulated range
c  ip1 =+1:  parameter p1 is above right boundary of tabulated range

c  ip2    :  same as ip1, for parameter p2.


      use EIRMOD_precision
      use EIRMOD_comxs, only: adas_data

      implicit none

      type(adas_data), pointer :: ad
      real(dp), intent(in) :: p1, p2
      integer, intent(out) :: ip1,ip2
      real(dp) :: res, rx, ry, ff(2,2)
      integer :: ide, ite

c  first parameter, find cell ite, ite+1
      if (p1 <= ad%temp(1)) then
        ite = 1
        ip1 = -1
      else if (p1 >= ad%temp(ad%ntemp)) then
        ite = ad%ntemp-1
        ip1 = 1
      else
        call EIRENE_binsearch_2 (ad%temp, ad%ntemp, p1, ite)
        ip1= 0
      end if
c  second parameter, find cell ide, ide+1
      if (p2 <= ad%dens(1)) then
        ide = 1
        ip2 = -1
      else if (p2 >= ad%dens(ad%ndens)) then
        ide = ad%ndens-1
        ip2 = 1
      else
        call EIRENE_binsearch_2 (ad%dens, ad%ndens, p2, ide)
        ip2 = 0
      end if

c  both cell indices ite and ide are found.
c  the parameter (p1,p2) is either in cell [ite, ite+1, ide, ide+1],
c  or outside the parameter range of this table.


      rx = (ad%temp(ite+1) - p1) * ad%dte(ite)
      ry = (ad%dens(ide+1) - p2) * ad%dde(ide)
c  fill ff: intermediate 2 by 2 matrix of tabulated values at cell vertices
      ff = ad%tab2d(ite:ite+1,ide:ide+1)

      call EIRENE_bilinear_int (ff, rx, ry, res)

      return
      end function EIRENE_intp_tab2d
