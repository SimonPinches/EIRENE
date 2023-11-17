
      SUBROUTINE EIRENE_SETUP_INMP
C
C  SET NON-DEFAULT STANDARD SURFACE IDENTIFIERS INMP...
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CLGIN
      USE EIRMOD_CGRID
      USE EIRMOD_CTRIG

      IMPLICIT NONE

      INTEGER :: ISTS, IR, J, K, JP, I, KT

      INMP1I=0
      INMP2I=0
      INMP3I=0

      DO 2019 ISTS=1,NSTSI

C  RADIAL SURFACE
        DO 2014 IR=1,NR1ST
          IF (IR.EQ.INUMP(ISTS,1)) THEN
            INMP1I(IR,0,0)=ISTS
            DO J=IRPTA(ISTS,2),IRPTE(ISTS,2)-1
              INMP1I(IR,J,0)=ISTS
              DO K=IRPTA(ISTS,3),IRPTE(ISTS,3)-1
                INMP1I(IR,0,K)=ISTS
                INMP1I(IR,J,K)=ISTS
              END DO
            END DO
          ENDIF
 2014   CONTINUE

C  POLOIDAL SURFACE
        DO 2016 JP=1,NP2ND
          IF (JP.EQ.INUMP(ISTS,2)) THEN
            INMP2I(0,JP,0)=ISTS
            DO I=IRPTA(ISTS,1),IRPTE(ISTS,1)-1
              INMP2I(I,JP,0)=ISTS
              DO K=IRPTA(ISTS,3),IRPTE(ISTS,3)-1
                INMP2I(0,JP,K)=ISTS
                INMP2I(I,JP,K)=ISTS
              END DO
            END DO
          ENDIF
 2016   CONTINUE

C  TOROIDAL SURFACE
        DO 2018 KT=1,NT3RD
          IF (KT.EQ.INUMP(ISTS,3)) THEN
            INMP3I(0,0,KT)=ISTS
            DO I=IRPTA(ISTS,1),IRPTE(ISTS,1)-1
              INMP3I(I,0,KT)=ISTS
              DO J=IRPTA(ISTS,2),IRPTE(ISTS,2)-1
                INMP3I(0,J,KT)=ISTS
                INMP3I(I,J,KT)=ISTS
              END DO
            END DO
            IF (LEVGEO == 4) THEN
              DO I=IRPTA(ISTS,1),IRPTE(ISTS,1)-1
                INMTI3(I,KT)=ISTS
              END DO
            END IF
          ENDIF
 2018   CONTINUE
C
 2019 CONTINUE

      RETURN
      END SUBROUTINE EIRENE_SETUP_INMP
