      SUBROUTINE EIRENE_PLTTLY
     .  (X,Y,VBAR,YMN,YMX,IR1,IR2,IRS,NKURV,TXTTAL,
     .                   TXTSPC,TXTUNT,TXTRUN,TXHEAD,
     .                   LBAR,XMI,XMA,YMNLG,YMXLG,LPLOT,LHIST,IERR,
     .                   N1BAR,N1DIM,L_SAME)

      USE EIRMOD_PRECISION

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: N1BAR, N1DIM
      REAL(DP), INTENT(IN) :: X(*), VBAR(N1BAR,*),
     .                      YMN(*), YMX(*), YMNLG(*), YMXLG(*)
      REAL(DP), INTENT(IN) :: XMI, XMA
      REAL(DP), INTENT(INOUT) :: Y(N1DIM,*)
      INTEGER, INTENT(IN) :: IR1(*), IR2(*), IRS(*), NKURV
      INTEGER, INTENT(OUT) :: IERR
      LOGICAL, INTENT(IN) :: LBAR(*), LPLOT(*), L_SAME
      LOGICAL, INTENT(IN) :: LHIST
      CHARACTER(LEN=*), INTENT(IN) :: TXTTAL(*),TXTSPC(*),TXTUNT(*),
     .                                TXTRUN, TXHEAD

      IERR = 0

      RETURN
      END
