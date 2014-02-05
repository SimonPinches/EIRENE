*//GEOMD//
C=======================================================================
C          S U B R O U T I N E   G E O M D _ L I N D A
C=======================================================================
      SUBROUTINE EIRENE_GEOMD_LINDA(NDXA,NDYA,NPLP,NR1ST,
     .                        PUX,PUY,PVX,PVY)
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CCONA
      USE EIRMOD_CPOLYG
      USE EIRMOD_CGEOM
      USE EIRMOD_COMPRT,ONLY:IUNIN,IUNOUT !VK
      IMPLICIT NONE
C
      INTEGER, INTENT(INOUT) :: NDXA, NDYA, NPLP
      INTEGER, INTENT(INOUT) :: NR1ST
      REAL(DP), INTENT(OUT) :: PUX(*),PUY(*),PVX(*),PVY(*)
C
C  GEOMETRY DATA: CELL VERTICES (LINDA ---> EIRENE)
      REAL(DP) :: 
     R  X1(NDX),Y1(NDX),X2(NDX),Y2(NDX),X3(NDX),Y3(NDX),
     R  X4(NDX),Y4(NDX)
C
      CHARACTER(80) :: LINE,LINE2
C   DIMENSIONIERUNG FUER GITTER
      INTEGER :: DIMXH,DIMYH,NNCUT,NNISO,
     1 NXCUT1(10),NXCUT2(10),NYCUT1(10),NYCUT2(10),
     2 NXISO1(10),NXISO2(10),NYISO1(10),NYISO2(10)
      INTEGER :: IX, IY, I, J, NP, NWISO

      REAL(DP) :: DUMMI(3)
      REAL(DP) :: MERK(NDY)
C  ACTUAL MESH USED IN THIS RUN
C
C      EINLESEROUTINE ANGEPASST AUF BRAAMS-OUTPUT
C   GEAENDERTE DIMENSIONIERUNG BZW. CUT-POSITION
C        MUSS PER HAND ANGEPASST WERDEN:
C        PARAMETER DIMXH,DIMYH                        RFS 14.5.1991
       NXCUT1=0
       NXCUT2=0
       NYCUT1=0
       NYCUT2=0
       NXISO1=0
       NXISO2=0
       NYISO1=0
       NYISO2=0
       nniso = -1 !pb
       OPEN (UNIT=30,ACCESS='SEQUENTIAL',FORM='FORMATTED',ERR=100) !VK
      REWIND 30
3366  FORMAT(/)
      read(30,*)
      do
        read (30,'(A80)') LINE
        i = verify(line,' ')
        if (i /= 0) then
           backspace 30
           exit
        end if
      end do

      read(30,*,ERR=100,END=100) dimxh,dimyh,nncut
      if(nncut.gt.10) stop 'Increase array sizes for cut'
      read(30,*,ERR=100,END=100) 
     r     (nxcut1(i),nxcut2(i),nycut1(i),nycut2(i),i=1,nncut)
      if (nncut.gt.2) then
         read(30,*,ERR=100,END=100) nniso
         if(nniso.gt.10) stop 'Increase array sizes for insulating cut'
         read(30,*,ERR=100,END=100) 
     r            (nxiso1(i),nxiso2(i),nyiso1(i),nyiso2(i),i=1,nniso)
      ELSE
       NNISO=0 !VK
      endif
      read(30,*,ERR=100,END=100)
C    ANZAHL DER TEILSTUECKE PRO POLYGON
      NPLP = MAX((NNCUT/2)*3,1)
C    COMPUTE WIDTH OF INSULATING CUT FOR DOUBLE NULL
      NWISO=NXISO2(1)-NXISO1(1)     
C     PRINT MESSAGE AND CHECK
      WRITE(IUNOUT,*) "GEOMD: NNCUT, NNISO, NPLP, NWISO ", 
     w                        NNCUT, NNISO, NPLP, NWISO 
      IF(NNCUT.NE.0.AND.NNCUT.NE.2.AND.NNCUT.NE.4) THEN
        WRITE(IUNOUT,*) "WARNING FROM GEOMD: UNKNOWN TOPOLOGY"
        WRITE(IUNOUT,*) " NNCUT ",NNCUT
      END IF      
      IF(NNCUT.EQ.4.AND.NNISO.NE.1) THEN
        WRITE(IUNOUT,*) "WARNING FROM GEOMD: UNKNOWN TOPOLOGY"
        WRITE(IUNOUT,*) " NNCUT, NNISO ",NNCUT,NNISO
      END IF      
      IF(NNCUT.EQ.2) THEN
       IF(NXCUT1(1).NE.NXCUT2(2)-1.OR.
     .    NXCUT1(2).NE.NXCUT2(1)-1) THEN
        WRITE(IUNOUT,*) "WARNING FROM GEOMD: CUTS DO NOT MATCH"
        WRITE(IUNOUT,*) " NXCUT1(1), NXCUT2(2) ",NXCUT1(1),NXCUT2(2)
        WRITE(IUNOUT,*) " NXCUT1(2), NXCUT2(2) ",NXCUT1(1),NXCUT2(2)
       END IF
      END IF
      IF(NNCUT.EQ.4) THEN
       IF(NXCUT1(1).NE.NXCUT2(NNCUT)-1.OR.
     .    NXCUT1(2).NE.NXCUT2(3)-1.OR.
     .    NXCUT1(3).NE.NXCUT2(2)-1.OR.
     .    NXCUT1(NNCUT).NE.NXCUT2(1)-1) THEN
        WRITE(IUNOUT,*) "WARNING FROM GEOMD: CUTS DO NOT MATCH"
        WRITE(IUNOUT,*) " NXCUT1(1), NXCUT2(4) ",NXCUT1(1),NXCUT2(4)
        WRITE(IUNOUT,*) " NXCUT1(2), NXCUT2(3) ",NXCUT1(2),NXCUT2(3)
        WRITE(IUNOUT,*) " NXCUT1(3), NXCUT2(2) ",NXCUT1(3),NXCUT2(2)
        WRITE(IUNOUT,*) " NXCUT1(4), NXCUT2(1) ",NXCUT1(4),NXCUT2(1)
       END IF
      END IF

C    READING OF POLYGON DATA
      DO 10 IX = 1, DIMXH
       IF (IX.LE.nxcut1(1)-1) THEN
        DO 12 IY = 1, DIMYH
         IF (IY.LE.dimyh-1) THEN
          READ (30,3333) DUMMI(1),
     .                  DUMMI(2),DUMMI(3),XPOL(IY,IX)
          READ (30,3333) DUMMI(1),
     .                  DUMMI(2),DUMMI(3),YPOL(IY,IX)
         ENDIF
         IF (IY.EQ.dimyh) THEN
          READ (30,3333) DUMMI(1),
     .                   DUMMI(2),XPOL(dimyh+1,IX),XPOL(IY,IX)
          READ (30,3333) DUMMI(1),
     .                   DUMMI(2),YPOL(dimyh+1,IX),YPOL(IY,IX)
         ENDIF
12      CONTINUE
       ENDIF
       IF (IX.EQ.nxcut1(1)) THEN
        DO 14 IY = 1, DIMYH
         IF (IY.LE.dimyh-1) THEN
         READ (30,3333) XPOL(IY,nxcut2(2)),DUMMI(1),
     .                  DUMMI(2),XPOL(IY,IX)
         READ (30,3333) YPOL(IY,nxcut2(2)),DUMMI(1),
     .                  DUMMI(2),YPOL(IY,IX)
         ENDIF
         IF (IY.EQ.dimyh) THEN
          READ (30,3333) XPOL(IY,nxcut2(2)),XPOL(dimyh+1,nxcut2(2)),
     .                                XPOL(dimyh+1,IX),XPOL(IY,IX)
          READ (30,3333) YPOL(IY,nxcut2(2)),YPOL(dimyh+1,nxcut2(2)),
     .                                YPOL(dimyh+1,IX),YPOL(IY,IX)
         ENDIF
14      CONTINUE
       ENDIF
       IF ((IX.GE.nxcut2(2)).AND.(IX.LE.nxcut1(2)-1)) THEN
        DO 16 IY = 1, DIMYH
         IF (IY.LE.dimyh-1) THEN
          READ (30,3333) DUMMI(1),
     .                   DUMMI(2),DUMMI(3),XPOL(IY,IX+1)
          READ (30,3333) DUMMI(1),
     .                   DUMMI(2),DUMMI(3),YPOL(IY,IX+1)
         ENDIF
         IF (IY.EQ.dimyh) THEN
          READ (30,3333) DUMMI(1),DUMMI(2),
     .                   XPOL(dimyh+1,IX+1),XPOL(IY,IX+1)
          READ (30,3333) DUMMI(1),DUMMI(2),
     .                   YPOL(dimyh+1,IX+1),YPOL(IY,IX+1)
         ENDIF
16      CONTINUE
       ENDIF
       IF (IX.EQ.nxcut1(2)) THEN
        DO 18 IY = 1, DIMYH
         IF (IY.LE.dimyh-1) THEN
         READ (30,3333) XPOL(IY,nxcut2(1)+1),DUMMI(1),
     .                  DUMMI(2),XPOL(IY,IX+1)
         READ (30,3333) YPOL(IY,nxcut2(1)+1),DUMMI(1),
     .                  DUMMI(2),YPOL(IY,IX+1)
         ENDIF
         IF (IY.EQ.dimyh) THEN
          READ (30,3333) XPOL(IY,nxcut2(1)+1),XPOL(dimyh+1,nxcut2(1)+1),
     .                                XPOL(dimyh+1,IX+1),XPOL(IY,IX+1)
          READ (30,3333) YPOL(IY,nxcut2(1)+1),YPOL(dimyh+1,nxcut2(1)+1),
     .                                YPOL(dimyh+1,IX+1),YPOL(IY,IX+1)
         ENDIF
18      CONTINUE
       ENDIF
       IF ((IX.GE.nxcut2(1)).AND.(IX.LE.dimxh-1)) THEN
        DO 22 IY = 1, DIMYH
         IF (IY.LE.dimyh-1) THEN
          READ (30,3333) DUMMI(1),DUMMI(2),DUMMI(3),
     .                   XPOL(IY,IX+2)
          READ (30,3333) DUMMI(1),DUMMI(2),DUMMI(3),
     .                   YPOL(IY,IX+2)
         ENDIF
         IF (IY.EQ.dimyh) THEN
          READ (30,3333) DUMMI(1),DUMMI(2),
     .                   XPOL(dimyh+1,IX+2),XPOL(IY,IX+2)
          READ (30,3333) DUMMI(1),DUMMI(2),
     .                   YPOL(dimyh+1,IX+2),YPOL(IY,IX+2)
         ENDIF
22      CONTINUE
       ENDIF
       IF (IX.EQ.dimxh) THEN
        DO 24 IY = 1, DIMYH
         IF (IY.LE.dimyh-1) THEN
         READ (30,3333) XPOL(IY,dimxh+3),DUMMI(1),DUMMI(2),
     .                  XPOL(IY,IX+2)
         READ (30,3333) YPOL(IY,dimxh+3),DUMMI(1),DUMMI(2),
     .                  YPOL(IY,IX+2)
         ENDIF
         IF (IY.EQ.dimyh) THEN
          READ (30,3333) XPOL(IY,dimxh+3),XPOL(dimyh+1,dimxh+3),
     .                                XPOL(dimyh+1,IX+2),XPOL(IY,IX+2)
          READ (30,3333) YPOL(IY,dimxh+3),YPOL(dimyh+1,dimxh+3),
     .                                YPOL(dimyh+1,IX+2),YPOL(IY,IX+2)
         ENDIF
24      CONTINUE
       ENDIF
10    CONTINUE

3333  FORMAT(4E15.7)

C   ANFANGSPUNKT DES ERSTEN TEILSTUECKS DES I-TEN POLYGONS
      NPOINT(1,1)=1
C   ENDPUNKT DES ERSTEN TEILSTUECKS DES I-TEN POLYGONS
      NPOINT(2,1)=nxcut1(1)+1
      IF (NNCUT.EQ.0) NPOINT(2,1)=dimxh+1
C   ANFANGSPUNKT DES ZWEITEN TEILSTUECKS DES I-TEN POLYGONS
      NPOINT(1,2)=nxcut2(nncut)+1
C   ENDPUNKT DES ZWEITEN TEILSTUECKS DES I-TEN POLYGONS
      NPOINT(2,2)=nxcut2(nncut-1)+1
C   ANFANGSPUNKT DES DRITTEN TEILSTUECKS DES I-TEN POLYGONS
      NPOINT(1,3)=nxcut2(nncut-1)+2
C   ENDPUNKT DES DRITTEN TEILSTUECKS DES I-TEN POLYGONS
      IF (NNCUT.EQ.2) NPOINT(2,3)=dimxh+3
      IF (NNCUT.EQ.4) THEN
       NPOINT(2,3)=NXISO1(1)+3
C   ANFANGSPUNKT DES VIERTEN TEILSTUECKS DES I-TEN POLYGONS
       NPOINT(1,4)=nxiso2(1)+4-NWISO
C   ENDPUNKT DES VIERTEN TEILSTUECKS DES I-TEN POLYGONS
       NPOINT(2,4)=nxcut1(3)+4-NWISO
C   ANFANGSPUNKT DES FUNFTEN TEILSTUECKS DES I-TEN POLYGONS
       NPOINT(1,5)=nxcut2(2)+4-NWISO
C   ENDPUNKT DES FUNFTEN TEILSTUECKS DES I-TEN POLYGONS
       NPOINT(2,5)=nxcut1(4)+5-NWISO
C   ANFANGSPUNKT DES SECHSTEN TEILSTUECKS DES I-TEN POLYGONS
       NPOINT(1,6)=nxcut2(1)+5-NWISO
C   ENDPUNKT DES SECHSTEN TEILSTUECKS DES I-TEN POLYGONS
       NPOINT(2,6)=dimxh+6-NWISO
      END IF
C
      DO 1015 IY=1,NDYA
        DO 1014 IX=1,NDXA
          X1(IX)=XPOL(IY,IX)
          Y1(IX)=YPOL(IY,IX)
          X2(IX)=XPOL(IY,IX+1)
          Y2(IX)=YPOL(IY,IX+1)
          X3(IX)=XPOL(IY+1,IX)
          Y3(IX)=YPOL(IY+1,IX)
          X4(IX)=XPOL(IY+1,IX+1)
          Y4(IX)=YPOL(IY+1,IX+1)
1014    CONTINUE
C
        CALL EIRENE_MSHPROJ (X1,Y1,X2,Y2,X3,Y3,X4,Y4,
     .                       PUX,PUY,PVX,PVY,NDXA,
     .                       NR1ST,IY)
1015  CONTINUE
C
      NP=NPOINT(2,NPLP)
      DO 1020 J=1,NDYA+1
        DO 1020 I=1,NP
          XPOL(J,I)=XPOL(J,I)*100.
          YPOL(J,I)=YPOL(J,I)*100.
          IF (ABS(XPOL(J,I)).LT.5.D-5) XPOL(J,I)=0.
          IF (ABS(YPOL(J,I)).LT.5.D-5) YPOL(J,I)=0.
1020  CONTINUE
      RETURN

 100  WRITE(IUNOUT,*) "COULD NOT OPEN FORT.30. ",
     w                "SKIP READING THE B2 GEOMETRY" !VK

*//END GEOMD_LINDA//
      END
