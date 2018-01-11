!********************************************************************
!*   Solve Y' = F(X,Y) with initial conditions using the Adams-     *
!*   Moulton Prediction-Correction Method                           *
!* ---------------------------------------------------------------- *
!*  INPUTS:                                                                *
!*    s1000     Given equation to integrate (see test program)             *
!*              must be declared EXTERNAL in calling program               *
!*    xi, xf    Begin, end values of variable x                            *
!*    yi        Begin Value of y at x=xi                                   *
!*    m         Number of points to calculate                              *
!*                                                                         *
!*  OUTPUTS:                                                               *
!*    xb      : real vector storing m nodes for function y                 *
!*    yb      : real vector storing m results for function y               *
!***************************************************************************

      subroutine adam_moulton (fp, xi, xf, yi, ec, m, xb, yb)

      real*8, intent(in) :: xi, xf, yi, ec
      integer, intent(in) :: m
      real*8, intent(out) :: xb(m+1), yb(m+1)
      real*8 :: fp

      real*8 X(0:3), Y(0:3)
      real*8 H,C1,C2,C3,C4,XX,YY,YC,YP
      integer K, L, ms
      
      H = (xf - xi) / real(m-1,kind(1.d0))  !integration step
      X(0) = xi
      Y(0) = yi               !Initial conditions

      ms = 1
      xb(ms) = xi
      yb(ms) = yi


!Start with Runge-Kutta
      DO K = 0, 1
        XX = X(K) 
        YY = Y(K) 
        C1 = FP(XX,YY)

        XX = X(K) + H / 2.d0
        YY = Y(K) + H / 2.d0 * C1 
        C2 = FP(XX,YY)

        YY = Y(K) + H / 2.d0 * C2
        C3 = FP(XX,YY)

        X(K + 1) = X(K) + H
        XX = X(K + 1)
        YY = Y(K) + H * C3 
        C4 = FP(XX,YY)

        Y(K + 1) = Y(K) + H * (C1 + 2 * C2 + 2 * C3 + C4) / 6.d0
        XX = X(K + 1) 
        
        ms = ms + 1
        xb(ms) = x(k+1)
        yb(ms) = y(k+1)

      END DO

! explicit precdictor step Adams-Bashford
 100  K = 2
      XX = X(K) 
      YY = Y(K) 
      C1 = FP(XX,YY)

      XX = X(K - 1)
      YY = Y(K - 1)
      C2 = FP(XX,YY)

      XX = X(K - 2) 
      YY = Y(K - 2) 
      C3 = FP(XX,YY)

      X(K + 1) = X(K) + H
      YP = Y(K) + H / 12.d0 * (23.d0 * C1 - 16.d0 * C2 + 5.d0 * C3)
      
! implicit corrector step Adams-Moulton
      L = 0
 200  XX = X(K + 1) 
      YY = YP 
      C1 = FP(XX,YY)

      XX = X(K) 
      YY = Y(K) 
      C2 = FP(XX,YY)

      XX = X(K - 1)
      YY = Y(K - 1) 
      C3 = FP(XX,YY)

      YC = Y(K) + H / 12.d0 * (5.d0 * C1 + 8.d0 * C2 - C3)
      
      IF (DABS(YP - YC) > EC) THEN
         YP = YC 
         L = L + 1 
         GOTO 200
      END IF
      
      Y(K + 1) = YC 

      ms = ms + 1
      xb(ms) = x(k+1)
      yb(ms) = y(k+1)
     
      DO K = 0, 2
        X(K) = X(K + 1) 
        Y(K) = Y(K + 1)
      END DO

      IF (X(2) < xf)  GOTO 100
           
      END
