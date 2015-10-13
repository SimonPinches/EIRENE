C
C
      SUBROUTINE EIRENE_OUTUSR
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD

      USE EIRMOD_CPLOT
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_COUTAU
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      USE EIRMOD_CSPEZ
      IMPLICIT NONE
      real(dp) :: dummy(nrtal)
      integer :: iadv
      integer :: i,ios,iun
      character*200 :: filename

1000  format(t5, a)
2000  format(20(1pe14.5E2))

      iun = 55

!...  write bfield out for plotting
      filename='bfield_for_plot'
      open(iun, file=filename, status='REPLACE', iostat=ios)
      IF (ios /= 0) THEN
         write(*,1000) 'Could NOT write ', filename, ' !'
         STOP
      END IF
      DO i=1,NRAD
         write(iun,2000) XCOM(i), YCOM(i), 0., BXIN(i), BYIN(i),
     >                   BZIN(i), BFIN(i), TEIN(i), TIIN(1,i),
     >                   DEIN(i), DIIN(1,i)
      END DO
      close(iun)

      RETURN
      END
