      function EIRENE_deter4x4 (a)
cdr  Determinant of (4,4) matrix A, Laplace expansion from first row A(1,1:4)
cdr  (better: find a row or column with as many as possible zeros and expand from there)
      use EIRMOD_precision
      implicit none

      real(dp), intent(in) :: a(4,4)
      real(dp) :: EIRENE_deter4x4, d11, d12, d13, d14, EIRENE_deter3x3
      external :: EIRENE_deter3x3

      d11 = EIRENE_deter3x3(a(2,2),a(3,2),a(4,2),
     .            a(2,3),a(3,3),a(4,3),
     .            a(2,4),a(3,4),a(4,4))

      d12 = EIRENE_deter3x3(a(2,1),a(3,1),a(4,1),
     .            a(2,3),a(3,3),a(4,3),
     .            a(2,4),a(3,4),a(4,4))

      d13 = EIRENE_deter3x3(a(2,1),a(3,1),a(4,1),
     .            a(2,2),a(3,2),a(4,2),
     .            a(2,4),a(3,4),a(4,4))

      d14 = EIRENE_deter3x3(a(2,1),a(3,1),a(4,1),
     .            a(2,2),a(3,2),a(4,2),
     .            a(2,3),a(3,3),a(4,3))

      EIRENE_deter4x4 = a(1,1)*d11 - a(1,2)*d12 +
     .                  a(1,3)*d13 - a(1,4)*d14

      return
      end function EIRENE_deter4x4
