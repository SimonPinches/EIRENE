C
C
C*DK DETER
      FUNCTION EIRENE_DETER3x3(A11,A21,A31,A12,A22,A32,A13,A23,A33)
C
C  DETERMINANT OF (3x3) MATRIX A, Laplace expansion, first row a11,a12,a13.
cdr   better find row or column with zeros first, and expand from there...
C
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: A11, A21, A31, A12, A22, A32, A13, A23,
     .                        A33
      REAL(DP) :: S, EIRENE_DETER3x3

      S=A11*(A22*A33-A23*A32)
      S=S-A12*(A21*A33-A23*A31)
      S=S+A13*(A21*A32-A22*A31)
      EIRENE_DETER3x3=S
      RETURN
      END FUNCTION EIRENE_DETER3x3
