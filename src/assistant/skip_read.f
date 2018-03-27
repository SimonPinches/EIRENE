      subroutine eirene_skip_read_comment(IREAD,IUNIN,ZEILE)

      IMPLICIT NONE

c  skip optional comment lines, starting with *, 
c  after an input block starts *** .....
c
c  input:  IREAD  =0: next input card is to be read from IUNIN
c          IREAD  =1: next input card is already given on ZEILE
c  output: IREAD  =1: always
c          ZEILE    : contains content of next input card 
c                     which is not a comment (not starting with "*").
c
c
      INTEGER :: IREAD,IUNIN
      CHARACTER(72) :: ZEILE, ZEILE_IN, ZEILE_OUT

331   CONTINUE
      IF (IREAD.EQ.0) THEN 
        READ (IUNIN,'(A72)') ZEILE_OUT
      ELSE
        ZEILE_IN=ZEILE
        READ (ZEILE_IN,'(A72)') ZEILE_OUT
      ENDIF
      IREAD=1
      IF (ZEILE_OUT(1:1).EQ.'*') THEN
        IREAD=0
C       WRITE (iunout,......)
        GOTO 331
      ENDIF
     
      ZEILE=ZEILE_OUT

c  now on exit:  iread=1:
c  next read must come from ZEILE:  READ(ZEILE, FORMAT) ....
c 
       return
       end
