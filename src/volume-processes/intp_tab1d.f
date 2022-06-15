      function EIRENE_intp_tab1d(tb,p1,ip1) result(res)

cdr  (linear) interpolate in 1d table tb at argument p1.
cdr  binary search in table. linear extrapolation outside range of table

cdr  Nov. 15: argument p2 removed, added: ip1:  indicator for extrapolation or interpolation
cdr           ip1: as in intp_tab2d (2d tables)
cdr      tbd:  generalize to general 1d tables, not just hydkin_data
cdr  Aug. 19: started. Old data structure hydkin rebuild, for 1D tables

      use EIRMOD_precision
      use EIRMOD_comxs, only: tab1d_data

      implicit none

      type(tab1d_data), pointer :: tb
      real(dp), intent(in) :: p1
      integer, intent(out) :: ip1
      real(dp) :: res
      integer :: ite

      if (p1 <= tb%temps(1)) then
        ite = 1
        ip1 = -1
      else if (p1 >= tb%temps(tb%ntemps)) then
        ite = tb%ntemps-1
        ip1 = 1
      else
        call EIRENE_binsearch_2 (tb%temps, tb%ntemps, p1, ite)
        ip1 = 0
      end if

      res = tb%tab1d(ite) + (p1-tb%temps(ite))*tb%diffquot(ite)

      return
      end function EIRENE_intp_tab1d
