C
*//GEOMD//
C=======================================================================
C          S U B R O U T I N E   G E O M D
C=======================================================================
      SUBROUTINE EIRENE_GEOMD(NDXA,NDYA,NPLP,NR1ST,
     .                 PUX,PUY,PVX,PVY,itype)
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CCONA
      USE EIRMOD_CGEOM
      IMPLICIT NONE
C
      REAL(DP), INTENT(OUT) :: PUX(*),PUY(*),PVX(*),PVY(*)
      INTEGER, INTENT(INOUT) :: NDXA,NDYA,NPLP,NR1ST
      INTEGER, INTENT(IN) :: ITYPE

      IF (ITYPE == 1) THEN
        CALL EIRENE_GEOMD_SONNET(NDXA,NDYA,NPLP,NR1ST,
     .                 PUX,PUY,PVX,PVY)
C
      ELSEIF (ITYPE == 2) THEN
        CALL EIRENE_GEOMD_CARRE(NDXA,NDYA,NPLP,NR1ST,
     .                 PUX,PUY,PVX,PVY)

      ELSE
        CALL EIRENE_GEOMD_LINDA(NDXA,NDYA,NPLP,NR1ST,
     .                 PUX,PUY,PVX,PVY)
      END IF

      RETURN
*//END GEOMD//
      END
