cdr  12.5.2015:  move general interface driver routine "EIRSRT" up, own routine.
cdr:             check: is eirsrt universal, then: move even further up to "main routines".
cdr  09.02.2016:  done ! synchronization of eirsrt.f started, but not completed fully
c  jan 2017: synchronisation with corresponding version in couple_SOLPS-ITER,
c            re. reading polygon data in geomd_linda from fort.30
c            added: species index in eapl,empl,eipl tallies
C
C  ASSISTANT ROUTINES, SPECIFIC TO A PARTICULAR EDGE CODE INTERFACE
C  DATA STRUCTURES (grid, plasma data, etc...)
c
cdr:  GEOMD : DRIVER FOR DIFFERENT VERSIONS OF GEOMD_..ROUTINES, diff. geometry file formats
C     GEOMD_CARRE
C     GEOMD_LINDA
C     GEOMD_SONNET
c
cdr   MESHPROJ
cdr   INDMAP
cdr   INDMPI
cdr   PLASM
cdr   NEUTR
cdr   SAVE_TALLIES
C
C
C
C
      SUBROUTINE EIRENE_MSHPROJ(X1,Y1,X2,Y2,X3,Y3,X4,Y4,PUX,PUY,PVX,PVY,
     .                   NDXA,NR1ST,IY)

      USE EIRMOD_PRECISION
      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: X1(*), Y1(*), X2(*), Y2(*),
     .                      X3(*), Y3(*), X4(*), Y4(*)
      REAL(DP), INTENT(OUT) :: PUX(*), PUY(*), PVX(*), PVY(*)
      INTEGER, INTENT(IN) :: NDXA, NR1ST, IY
      REAL(DP) :: D12, D34, D13, D24, EPS60, PUPV, PVPV, DVX, DVY,
     .          DUX, DUY
      INTEGER :: IX, IN

      EPS60 = 1.E-60_DP
C
C
      DO 1 IX=1,NDXA
C
C  CALCULATE THE NORM OF THE VECTORS (POINT2-POINT1),....
C
        D12 = SQRT((X2(IX)-X1(IX))*(X2(IX)-X1(IX))+(Y2(IX)-Y1(IX))*
     .        (Y2(IX)-Y1(IX)))+EPS60
        D34 = SQRT((X4(IX)-X3(IX))*(X4(IX)-X3(IX))+(Y4(IX)-Y3(IX))*
     .        (Y4(IX)-Y3(IX)))+EPS60
        D13 = SQRT((X3(IX)-X1(IX))*(X3(IX)-X1(IX))+(Y3(IX)-Y1(IX))*
     .        (Y3(IX)-Y1(IX)))+EPS60
        D24 = SQRT((X4(IX)-X2(IX))*(X4(IX)-X2(IX))+(Y4(IX)-Y2(IX))*
     .        (Y4(IX)-Y2(IX)))+EPS60
C
C  CALCULATE THE BISSECTING VECTORS, BUT NOT NORMALISED YET
C
        DUX = (X2(IX)-X1(IX))/D12 + (X4(IX)-X3(IX))/D34
        DUY = (Y2(IX)-Y1(IX))/D12 + (Y4(IX)-Y3(IX))/D34
        DVX = (X3(IX)-X1(IX))/D13 + (X4(IX)-X2(IX))/D24
        DVY = (Y3(IX)-Y1(IX))/D13 + (Y4(IX)-Y2(IX))/D24
C
C  CALCULATE THE COMPONENTS OF THE TWO UNIT VECTOR (= PROJECTION RATE)
C
        IN=IY+(IX-1)*NR1ST
        PUX(IN) = DUX/(SQRT(DUX*DUX+DUY*DUY)+EPS60)
        PUY(IN) = DUY/(SQRT(DUX*DUX+DUY*DUY)+EPS60)
        PVX(IN) = DVX/(SQRT(DVX*DVX+DVY*DVY)+EPS60)
        PVY(IN) = DVY/(SQRT(DVX*DVX+DVY*DVY)+EPS60)
C
C  ORTHOGONORMALIZE, CONSERVE ORIENTATION (E.SCHMIDT)
C
        PUPV=PUX(IN)*PVX(IN)+PUY(IN)*PVY(IN)
        PVX(IN)=PVX(IN)-PUPV*PUX(IN)
        PVY(IN)=PVY(IN)-PUPV*PUY(IN)
        PVPV=SQRT(PVX(IN)*PVX(IN)+PVY(IN)*PVY(IN))+EPS60
        PVX(IN)=PVX(IN)/PVPV
        PVY(IN)=PVY(IN)/PVPV
C
    1 CONTINUE
      RETURN
      END


C
C
      SUBROUTINE EIRENE_INDMPI(FIELD,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .                  NCUTB,NCUTL,NPOINT,NPPLG,NSTR,ISTR)
C
C     INDEX MAPPING: INVERSE TO SUBR. INDMAP
C
      USE EIRMOD_PRECISION
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NPOINT(2,*)
      INTEGER, INTENT(IN) :: NDX, NDY, NFL, NDXA, NDYA, NFLA, NCUTB,
     .                       NCUTL, NPPLG, NSTR, ISTR
      REAL(DP), INTENT(INOUT) :: FIELD(0:NDX+1,0:NDY+1,NFL,NSTR),
     .                         DUMMY(0:NDX+1,0:NDY+1)

      RETURN
C
      END
