C
cdr  aug. 20: code safeties from ITER branch
C
      SUBROUTINE EIRENE_MASIR2 (A,NLFLD,NFIRST0,NFIRST1,N0,N1,M)
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: A
      INTEGER, INTENT(IN) :: NFIRST0, NFIRST1, N0, N1, M
      INTEGER, INTENT(IN) :: NLFLD(NFIRST0:NFIRST1,*)
      INTEGER :: I, J, K, JA, IZ, JE, NCH
C
      NCH=110
      IZ=M/NCH
      IF (MOD(M,NCH).NE.0) IZ=IZ+1
C
      WRITE (iunout,'(///1X,A)') A
      DO 1 I=1,IZ
         JA=(I-1)*NCH+1
         JE=MIN(I*NCH,M)
         WRITE (iunout,'(/7X,A1,12(1X,10I1))') 'I',(MOD(J,10),J=JA,JE)
         WRITE (iunout,'(1X,140A1)') ('-',J=JA,JE+18)
         DO 2 K=N0,N1
            WRITE (iunout,'(1X,I5,1X,A1,12(1X,10I1))') K,'I',
     .            (NLFLD(K,J),J=JA,JE)
    2    CONTINUE
    1 CONTINUE
      RETURN
      END SUBROUTINE EIRENE_MASIR2
