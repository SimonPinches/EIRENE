chk  March 24: Evaluate the nth Legendre polynomial at x
              
      SUBROUTINE EIRENE_LEGENDRE(n,x,l)

      USE EIRMOD_PRECISION
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      REAL(DP), INTENT(IN) :: x
      REAL(DP), INTENT(OUT) :: l
      REAL(dp) :: l0, l1
      INTEGER :: i 
      SELECT CASE(n)
chk Degrees 10 and lower pre-computed for speed
        case(0)
          l=1
        case(1)
          l=x
        case(2)
          l=1.5*x*x - 0.5
        case(3)
          l=2.5*x**3 - 1.5*x
        case(4)
          l=4.375*x**4 - 3.75*x*x + 0.375
        case(5)
          l=7.875*x**5 - 8.75*x**3 + 1.875*x
        case(6)
          l=(231*x**6 - 315*x**4 + 105*x**2 - 5)/16.0
        case(7)
          l=(429*x**7 - 693*x**5 + 315*x**3 - 35*x)/16.0
        case(8)
          l=(6435*x**8-12012*x**6+6930*x**4-1260*x**2+35)/128.0
        case(9)
          l=(12155*x**9-25740*x**7+18018*x**5-4620*x**3+315*x)/128.0
        case(10)
          l=(46189*x**10-109395*x**8+90090*x**6-30030*x**4+3465*x**2-63)
     .       / 256.0
        case default
chk Compute higher degree polynomials with a loop
          l0  = 1
          l1  = x
          DO i = 2, n
            l = (2 * i-1) * x * l1 / i - (i-1) * l0 / i
            l0 = l1
            l1 = l
          END DO
        END SELECT
        RETURN
      END

