cdr  june 16: comments added
cdr  NS=N1=1 always, except: call from STEP.f
cdr  to be done: binary search
C
C*DK LEARCA
      FUNCTION EIRENE_LEARCA (X,R,N1,N,NS,TEXT)
C
C   THIS FUNCTION COMPUTES THE INDEX OF THE SMALLER MESHPOINT OF THE
C   INTERVAL CONTAINING THE POINT X IN THE 2ND COORDINATE OF A 2D MESH R(NS,I),I=1,N
C   N1 IS THE LEADING DIMENSION OF THE FIELD R, AS SPECIFIED IN THE
C   CALLING PROGRAM AND NS IS A FIXED INDEX, TO REDUCE A 2D FIELD R(J,I)
C   TO A 1D STRUCTURE R(J=NS,I)

C   return learca=1 or learca=N, if x out of bounds at left or right end, respectively
C
      USE EIRMOD_PRECISION
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
 
      CHARACTER(*), INTENT(IN) :: TEXT
      INTEGER, INTENT(IN) :: N1, N, NS
      REAL(DP), INTENT(IN) :: R(N1,*)
      REAL(DP), INTENT(IN) :: X
      INTEGER :: NNN, I, J, EIRENE_LEARCA
 
      NNN=1
      IF (X.LT.R(NS,1)-1.D-12) GOTO 20

cdr  this loop should be replaced with a binary search
13    DO 10 J=2,N
        I=J
        IF (X-R(NS,J).LE.0.0) GOTO 15
10    CONTINUE


      NNN=N
      IF (X.GT.R(NS,N)+1.D-12) GOTO 20


15    EIRENE_LEARCA=I-1
      RETURN
C
20    WRITE (iunout,*) 'X OUT OF RANGE IN LEARCA'
      WRITE (iunout,*)  X,NNN,R(NS,NNN)
      WRITE (iunout,*) 'LEARCA= ',NNN,' RETURNED TO SUBR. ',TEXT
      EIRENE_LEARCA=NNN
C
      RETURN
      END
