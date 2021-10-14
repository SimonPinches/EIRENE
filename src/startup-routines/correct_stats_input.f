
      SUBROUTINE EIRENE_CORRECT_STATS_INPUT
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CSDVI
      USE EIRMOD_CESTIM
      USE EIRMOD_COMPRT, ONLY : IUNOUT
      IMPLICIT NONE
      INTEGER :: J

      IF (NSIGVI > 0) THEN
        J=1
        DO WHILE (J <= NSIGVI)
          IF (LMISTALV(IIH(J))) THEN
            WRITE (iunout,*)
     .        ' NO STATISTICS IS DONE FOR VOLUME TALLY NO. ',IIH(J)
            WRITE (iunout,*) ' BECAUSE TALLY HAS BEEN SWITCHED OFF'
            IGH(J)=IGH(NSIGVI)
            IIH(J)=IIH(NSIGVI)
            NSIGVI=NSIGVI-1
          ELSE
            J=J+1
          END IF
        END DO
      END IF

      IF (NSIGSI > 0) THEN
        J=1
        DO WHILE (J <= NSIGSI)
          IF (LMISTALS(IIHW(J))) THEN
            WRITE (iunout,*)
     .        ' NO STATISTICS IS DONE FOR SURFACE TALLY NO. ',IIHW(J)
            WRITE (iunout,*) ' BECAUSE TALLY HAS BEEN SWITCHED OFF'
            IGHW(J)=IGHW(NSIGSI)
            IIHW(J)=IIHW(NSIGSI)
            NSIGSI=NSIGSI-1
          ELSE
            J=J+1
          END IF
        END DO
      END IF

      IF (NSIGCI > 0) THEN
        J=1
        DO WHILE (J <= NSIGCI)
          IF (LMISTALV(IIHC(1,J)) .OR. LMISTALV(IIHC(2,J))) THEN
            WRITE (iunout,*)
     .        ' NO CORRELATION COEFFICIENT IS CALCULATED',
     .        ' BETWEEN TALLIES ',IIHC(1,J),' AND ',IIHC(2,J)
            WRITE (iunout,*) ' BECAUSE TALLIES HAVE BEEN SWITCHED OFF'
            IGHC(1,J)=IGHC(1,NSIGCI)
            IIHC(1,J)=IIHC(1,NSIGCI)
            IGHC(2,J)=IGHC(2,NSIGCI)
            IIHC(2,J)=IIHC(2,NSIGCI)
            NSIGCI=NSIGCI-1
          ELSE
            J=J+1
          END IF
        END DO
      END IF
C
      NSIGI=NSIGVI+NSIGSI+NSIGCI
      
      RETURN
      END SUBROUTINE EIRENE_CORRECT_STATS_INPUT
      
