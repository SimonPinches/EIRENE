cdr  june 16: comments added
cdr  same as learca2, but here: NS=N1=1 always, except: learca2 call from STEP.f
cdr  to be done: binary search
C
C*DK LEARCA1
      FUNCTION EIRENE_LEARCA1 (X,R,N,TEXT)
C
C   THIS FUNCTION COMPUTES THE INDEX OF THE SMALLER MESHPOINT OF THE
C   INTERVAL CONTAINING THE POINT X IN THE 1D MESH R(I),I=1,N
C

C   return learca1=1 or learca1=N, if x out of bounds at left or right end, respectively
C
      USE EIRMOD_PRECISION
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: TEXT
      INTEGER, INTENT(IN) :: N
      REAL(DP), INTENT(IN) :: R(*)
      REAL(DP), INTENT(IN) :: X
      INTEGER :: NNN, I, J, EIRENE_LEARCA1

      NNN=1
      IF (X.LT.R(1)-1.D-12) GOTO 20

cdr  this loop should be replaced with a binary search
   13 DO 10 J=2,N
        I=J
        IF (X-R(J).LE.0.0) GOTO 15
   10 CONTINUE


      NNN=N
      IF (X.GT.R(N)+1.D-12) GOTO 20


   15 EIRENE_LEARCA1=I-1
      RETURN
C
   20 WRITE (iunout,*) 'X OUT OF RANGE IN LEARCA1'
      WRITE (iunout,*)  X,NNN,R(NNN)
      WRITE (iunout,*) 'LEARCA1= ',NNN,' RETURNED TO SUBR. ',TEXT
      EIRENE_LEARCA1=NNN
C
      RETURN
      END
