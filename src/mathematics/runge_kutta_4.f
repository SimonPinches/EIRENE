!***************************************************************************
!*        SOLVING DIFFERENTIAL EQUATIONS WITH 1 VARIABLE OF ORDER 1        *
!*                       of type dy/dx = f(x,y)                            *
!*        BY RUNGE-KUTTA METHOD OF ORDER 4                                 *
!* ----------------------------------------------------------------------- *
!*  INPUTS:                                                                *
!*    fp        Given equation to integrate (see test program)             *
!*              must be declared EXTERNAL in calling program               *
!*    xi, xf    Begin, end values of variable x                            *
!*    yi        Begin Value of y at x=xi                                   *
!*    m         Number of points to calculate                              *
!*    fi        finesse (number of intermediate points)                    *
!*                                                                         *
!*  OUTPUTS:                                                               *
!*    xb      : real vector storing m nodes for function y                 *
!*    t       : real vector storing m results for function y               *
!***************************************************************************
      Subroutine Runge_Kutta_4(fp,xi,xf,yi,m,fi,xb,t)
      USE EIRMOD_PRECISION
      implicit none
      integer fi,i,j,m,ni
      real(DP), intent(out) :: xb(m+1),t(m+1)
      real(DP) :: xi,xf,yi,a,b,c,d,h,x,y,fp

      if (fi < 1) return

      h  = (xf - xi) / fi / (m-1)
      y  = yi
      xb(1) = xi
      t(1) = yi

      do i=1, m
        ni = (i - 1) * fi - 1
        do j=1, fi
          x = xi + h * (ni + j)
          a = h * fp(x,y)
          b = h * fp(x+h/2.d0,y+a/2.d0)
          c = h * fp(x+h/2.d0,y+b/2.d0)
          x = x + h
          d = h * fp(x,y+c)
          y = y + (a + b + b + c + c + d) / 6.d0
        end do
        xb(i+1) = x
        t(i+1) = y
      end do

      return
      end
