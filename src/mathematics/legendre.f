chk  March 24: Evaluate the nth Legendre polynomial at x

      SUBROUTINE EIRENE_LEGENDRE(n,x,l)

      USE EIRMOD_PRECISION
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      REAL(DP), INTENT(IN) :: x
      REAL(DP), INTENT(OUT) :: l
      INTEGER :: i
      REAL(DP) :: l0, l1

      SELECT CASE(n)
chk Degrees 10 and lower pre-computed for speed
      case(0)
        l=1.0_DP
      case(1)
        l=x
      case(2)
        l=1.5_DP*x*x - 0.5_DP
      case(3)
        l=2.5_DP*x**3 - 1.5_DP*x
      case(4)
        l=4.375_DP*x**4 - 3.75_DP*x*x + 0.375_DP
      case(5)
        l=7.875_DP*x**5 - 8.75_DP*x**3 + 1.875_DP*x
      case(6)
        l=(231*x**6 - 315*x**4 + 105*x**2 - 5.0_DP)/16.0_DP
      case(7)
        l=(429*x**7 - 693*x**5 + 315*x**3 - 35*x)/16.0_DP
      case(8)
        l=(6435*x**8-12012*x**6+6930*x**4-1260*x**2+35.0_DP)/128.0_DP
      case(9)
        l=(12155*x**9-25740*x**7+18018*x**5-4620*x**3+315*x)/128.0_DP
      case(10)
        l=(46189*x**10-109395*x**8+90090*x**6-30030*x**4+
     .      3465*x**2-63.0_DP)/ 256.0_DP
      case default
chk Compute higher degree polynomials with a loop
        l0  = 1.0_DP
        l1  = x
        DO i = 2, n
          l = (2 * i-1) * x * l1 / real(i,DP)
     .              - (i-1) * l0 / real(i,DP)
          l0 = l1
          l1 = l
        END DO
      END SELECT
      RETURN
      END SUBROUTINE EIRENE_LEGENDRE
