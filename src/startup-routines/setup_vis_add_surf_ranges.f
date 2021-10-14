
      SUBROUTINE EIRENE_SETUP_VIS_ADD_SURF_RANGES
      
C  SET 'VISIBLE ADDITIONAL SURFACES' RANGES nlimii(j),nlimie(j), for each grid cell j
C  FROM INFORMATION ON IGJUM3
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CLGIN
      USE EIRMOD_CGRID
      USE EIRMOD_CADGEO
      USE EIRMOD_COMUSR, ONLY : NBITS
      
      IMPLICIT NONE
      INTEGER :: J, I, IIN, IEN , NO, IGO, IBEND, IB, ILA, NSOPT,
     .           EIRENE_ILLZ
      LOGICAL :: LHELP(NLIMPS)
      
C  DEFAULT
      NLIMII=1
      NLIMIE=NLIMI
C
      NSOPT=MIN(NOPTIM,NSBOX)
      IF (NLIMPB >= NLIMPS) THEN

        DO J=1,NSOPT
            DO 8005 I=1,NLIMI
              LHELP(I) = IGJUM3(J,I)==0
 8005       CONTINUE
            IIN=EIRENE_ILLZ(NLIMI,LHELP,1)+1
            IEN=NLIMI-EIRENE_ILLZ(NLIMI,LHELP,-1)
            NLIMII(J)=IIN
            NLIMIE(J)=IEN
          ENDDO   ! NSOPT LOOP, GRID CELLS


        ELSEIF (NLIMPB < NLIMPS) THEN
C  now try the same thing but with bit arithmetic, in case of storage
c  saving mode (igjum3 array stored in single bit integer format)

          DO J=1,NSOPT
            IIN = 1
            IEN = NLIMI
! NO='1111....111'B ALL BITS SET TO 1
            NO=NOT(0)
! NUMBER OF INTEGERS USED TO STORE SURFACE INFORMATION
            IGO=NLIMI/NBITS
            IF (MOD(NLIMI,NBITS) > 0) IGO = IGO + 1
! CHECK FOR FIRST ACTIVE SURFACE
            DO I=1,IGO
              IF (IAND(IGJUM3(J,I),NO) /= NO) THEN
                IBEND = NBITS-1
                IF (I == IGO) IBEND = NLIMI-(I-1)*NBITS - 1
                DO IB=0,IBEND
                  IF (.NOT.BTEST(IGJUM3(J,I),IB)) THEN
                    ILA=IB
                    EXIT
                  END IF
                END DO
                IIN = (I-1)*NBITS+ILA+1
                IF (IB <= IBEND) EXIT
              END IF
            END DO
! CHECK FOR LAST ACTIVE SURFACE
            DO I=IGO,1,-1
              IF (IAND(IGJUM3(J,I),NO) /= NO) THEN
                IBEND = NBITS-1
                IF (I == IGO) IBEND = NLIMI-(I-1)*NBITS - 1
                DO IB=IBEND,0,-1
                  IF (.NOT.BTEST(IGJUM3(J,I),IB)) THEN
                    ILA=IB
                    EXIT
                  END IF
                END DO
                IEN = (I-1)*NBITS+ILA+1
                IF (IB >= 0) EXIT
              END IF
            END DO  !
            NLIMII(J)=IIN
            NLIMIE(J)=IEN
          ENDDO   ! NSOPT LOOP, GRID CELLS

        ENDIF

      RETURN
      END SUBROUTINE EIRENE_SETUP_VIS_ADD_SURF_RANGES
