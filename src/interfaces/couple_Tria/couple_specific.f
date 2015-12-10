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

C
*//GEOMD//
C=======================================================================
C          S U B R O U T I N E   G E O M D _ C A R R E
C=======================================================================
      SUBROUTINE EIRENE_GEOMD_CARRE(NDXA,NDYA,NPLP,NR1ST,
     .                 PUX,PUY,PVX,PVY)
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CCONA
      USE EIRMOD_CGEOM
      IMPLICIT NONE
C
      REAL(DP), INTENT(OUT) :: PUX(*),PUY(*),PVX(*),PVY(*)
      INTEGER, INTENT(INOUT) :: NDXA,NDYA,NPLP,NR1ST

      character(110) :: zeile
      REAL(DP) :: br(0:ndxp,0:ndyp,4),bz(0:ndxp,0:ndyp,4)
C
C  GEOMETRY DATA: CELL VERTICES (LINDA ---> EIRENE)
      REAL(DP) ::
     R  X1(NDX),Y1(NDX),X2(NDX),Y2(NDX),X3(NDX),Y3(NDX),
     R  X4(NDX),Y4(NDX)
      REAL(DP) :: DX, DY
      INTEGER :: I0, I0E, IX, IY, I1, I2, I3, I4, IPART, I, J, NCUT
      CHARACTER(80) :: LINE
      LOGICAL :: LRDCUT

      ndxa=0
      ndya=0
      lrdcut = .false.

!pb      DO I = 1, 4
      DO 
        read (30,'(A80)') line
        if (.not.lrdcut) then
          call eirene_uppercase (line)
          i0 = index(line,'NCUT')
          if (i0 > 0) then
            i0 = index(line,'=') + 1
            i1 = i0 + verify(line(i0:),' ')-1
            i2 = i1 + scan(line(i1+1:),' ')-1
            read (line(i1:i2),*) ncut
            read (30,'(A80)') line
            i1 = index(line,'=')
            read (line(i1+1:),*) (npoint(2,j),j=1,ncut)
            do j=1,ncut
              npoint(2,j) = npoint(2,j) + j
              npoint(1,j+1) = npoint(2,j) + 1
            end do
            npoint(1,1) = 1
            ipart = ncut+1
            lrdcut = .true.
          end if 
        end if
        if (index(line,'=======') /= 0) exit
      END DO

1     continue
      read (30,'(a110)',end=99) zeile
      i0=index(zeile,'(')
      i0e=index(zeile,')')
      read (zeile(i0+1:i0e-1),*) ix,iy
      ndxa=max(ndxa,ix)
      ndya=max(ndya,iy)
      i1=index(zeile,': (')
      i2=index(zeile(i1+3:),')')+i1+2
      read (zeile(i1+3:i2-1),*) br(ix,iy,4),bz(ix,iy,4)
      i3=index(zeile(i2+1:),'(')+i2
      i4=i3+index(zeile(i3+1:),')')
      read (zeile(i3+1:i4-1),*) br(ix,iy,3),bz(ix,iy,3)

      read (30,'(a110)') zeile

      read (30,'(a110)') zeile
      i1=index(zeile,'(')
      i2=index(zeile,')')
      read (zeile(i1+1:i2-1),*) br(ix,iy,1),bz(ix,iy,1)
      i3=i2+index(zeile(i2+1:),'(')
      i4=i2+index(zeile(i2+1:),')')
      read (zeile(i3+1:i4-1),*) br(ix,iy,2),bz(ix,iy,2)

      read (30,*)
      goto 1


99    continue
      ndxa=ndxa-1
      ndya=ndya-1
C
!pb      DO 1015 IY=1,NDYA
!pb        DO 1014 IX=1,NDXA
!pb          X1(IX)=br(ix,iy,1)
!pb          Y1(IX)=bz(ix,iy,1)
!pb          X2(IX)=br(ix,iy,2)
!pb          Y2(IX)=bz(ix,iy,2)
!pb          X3(IX)=br(ix,iy,4)
!pb          Y3(IX)=bz(ix,iy,4)
!pb          X4(IX)=br(ix,iy,3)
!pb          Y4(IX)=bz(ix,iy,3)
!pb1014    CONTINUE
!pb        CALL MSHPROJ (X1,Y1,X2,Y2,X3,Y3,X4,Y4,PUX,PUY,PVX,PVY,NDXA,
!pb     .                NR1ST,IY)
!pb1015  CONTINUE
C
C SEARCH FOR THE CUTS
C
      IF (.NOT.LRDCUT) THEN
        IPART=1
        NPOINT(1,IPART)=1
        IY=1
        DO IX=1,NDXA
          DX=BR(IX+1,IY,1)-BR(IX,IY,2)
          DY=BZ(IX+1,IY,1)-BZ(IX,IY,2)
          IF (DX*DX+DY*DY.GT.EPS10) THEN
C CUT GEFUNDEN
            NPOINT(2,IPART)=IX+1+IPART-1
            IPART=IPART+1
            NPOINT(1,IPART)=IX+1+IPART-1
          ENDIF
        ENDDO
      END IF
      NPOINT(2,IPART)=NDXA+IPART
C
      NPLP=IPART
C
      DO IY=1,NDYA
        DO IPART=1,NPLP
          DO IX=NPOINT(1,IPART),NPOINT(2,IPART)-1
            XPOL(IY,IX)=BR(IX-(IPART-1),IY,1)
            YPOL(IY,IX)=BZ(IX-(IPART-1),IY,1)
          ENDDO
          XPOL(IY,NPOINT(2,IPART))=BR(NPOINT(2,IPART)-IPART,IY,2)
          YPOL(IY,NPOINT(2,IPART))=BZ(NPOINT(2,IPART)-IPART,IY,2)
        ENDDO
      ENDDO
C INTRODUCE OUTERMOST RADIAL POLYGON
      DO IPART=1,NPLP
        DO IX=NPOINT(1,IPART),NPOINT(2,IPART)-1
          XPOL(NDYA+1,IX)=BR(IX-(IPART-1),NDYA,4)
          YPOL(NDYA+1,IX)=BZ(IX-(IPART-1),NDYA,4)
        ENDDO
        XPOL(NDYA+1,NPOINT(2,IPART))=BR(NPOINT(2,IPART)-IPART,NDYA,3)
        YPOL(NDYA+1,NPOINT(2,IPART))=BZ(NPOINT(2,IPART)-IPART,NDYA,3)
      ENDDO
C
      DO J=1,NDYA+1
        DO I=1,NPOINT(2,NPLP)
          XPOL(J,I)=XPOL(J,I)*100.
          YPOL(J,I)=YPOL(J,I)*100.
        END DO
      END DO
C
      ndxa=npoint(2,nplp)-1

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
        CALL EIRENE_MSHPROJ (X1,Y1,X2,Y2,X3,Y3,X4,Y4,PUX,PUY,PVX,PVY,
     .                NDXA,NR1ST,IY)
1015  CONTINUE
c
C     do j=1,ndya+1
C       write (iunout,*)
C       write (iunout,*) 'in geomd polygon ',j
C       write (iunout,'(1p,6e12.4)') 
C    .        (xpol(j,i),ypol(j,i),i=1,npoint(2,nplp))
C     enddo
C
      RETURN
*//END GEOMD_CARRE//
      END

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
          READ (30,*) DUMMI(1),
     .                  DUMMI(2),DUMMI(3),XPOL(IY,IX)
          READ (30,*) DUMMI(1),
     .                  DUMMI(2),DUMMI(3),YPOL(IY,IX)
         ENDIF
         IF (IY.EQ.dimyh) THEN
          READ (30,*) DUMMI(1),
     .                   DUMMI(2),XPOL(dimyh+1,IX),XPOL(IY,IX)
          READ (30,*) DUMMI(1),
     .                   DUMMI(2),YPOL(dimyh+1,IX),YPOL(IY,IX)
         ENDIF
12      CONTINUE
       ENDIF
       IF (IX.EQ.nxcut1(1)) THEN
        DO 14 IY = 1, DIMYH
         IF (IY.LE.dimyh-1) THEN
         READ (30,*) XPOL(IY,nxcut2(2)),DUMMI(1),
     .                  DUMMI(2),XPOL(IY,IX)
         READ (30,*) YPOL(IY,nxcut2(2)),DUMMI(1),
     .                  DUMMI(2),YPOL(IY,IX)
         ENDIF
         IF (IY.EQ.dimyh) THEN
          READ (30,*) XPOL(IY,nxcut2(2)),XPOL(dimyh+1,nxcut2(2)),
     .                                XPOL(dimyh+1,IX),XPOL(IY,IX)
          READ (30,*) YPOL(IY,nxcut2(2)),YPOL(dimyh+1,nxcut2(2)),
     .                                YPOL(dimyh+1,IX),YPOL(IY,IX)
         ENDIF
14      CONTINUE
       ENDIF
       IF ((IX.GE.nxcut2(2)).AND.(IX.LE.nxcut1(2)-1)) THEN
        DO 16 IY = 1, DIMYH
         IF (IY.LE.dimyh-1) THEN
          READ (30,*) DUMMI(1),
     .                   DUMMI(2),DUMMI(3),XPOL(IY,IX+1)
          READ (30,*) DUMMI(1),
     .                   DUMMI(2),DUMMI(3),YPOL(IY,IX+1)
         ENDIF
         IF (IY.EQ.dimyh) THEN
          READ (30,*) DUMMI(1),DUMMI(2),
     .                   XPOL(dimyh+1,IX+1),XPOL(IY,IX+1)
          READ (30,*) DUMMI(1),DUMMI(2),
     .                   YPOL(dimyh+1,IX+1),YPOL(IY,IX+1)
         ENDIF
16      CONTINUE
       ENDIF
       IF (IX.EQ.nxcut1(2)) THEN
        DO 18 IY = 1, DIMYH
         IF (IY.LE.dimyh-1) THEN
         READ (30,*) XPOL(IY,nxcut2(1)+1),DUMMI(1),
     .                  DUMMI(2),XPOL(IY,IX+1)
         READ (30,*) YPOL(IY,nxcut2(1)+1),DUMMI(1),
     .                  DUMMI(2),YPOL(IY,IX+1)
         ENDIF
         IF (IY.EQ.dimyh) THEN
          READ (30,*) XPOL(IY,nxcut2(1)+1),XPOL(dimyh+1,nxcut2(1)+1),
     .                                XPOL(dimyh+1,IX+1),XPOL(IY,IX+1)
          READ (30,*) YPOL(IY,nxcut2(1)+1),YPOL(dimyh+1,nxcut2(1)+1),
     .                                YPOL(dimyh+1,IX+1),YPOL(IY,IX+1)
         ENDIF
18      CONTINUE
       ENDIF
       IF ((IX.GE.nxcut2(1)).AND.(IX.LE.dimxh-1)) THEN
        DO 22 IY = 1, DIMYH
         IF (IY.LE.dimyh-1) THEN
          READ (30,*) DUMMI(1),DUMMI(2),DUMMI(3),
     .                   XPOL(IY,IX+2)
          READ (30,*) DUMMI(1),DUMMI(2),DUMMI(3),
     .                   YPOL(IY,IX+2)
         ENDIF
         IF (IY.EQ.dimyh) THEN
          READ (30,*) DUMMI(1),DUMMI(2),
     .                   XPOL(dimyh+1,IX+2),XPOL(IY,IX+2)
          READ (30,*) DUMMI(1),DUMMI(2),
     .                   YPOL(dimyh+1,IX+2),YPOL(IY,IX+2)
         ENDIF
22      CONTINUE
       ENDIF
       IF (IX.EQ.dimxh) THEN
        DO 24 IY = 1, DIMYH
         IF (IY.LE.dimyh-1) THEN
         READ (30,*) XPOL(IY,dimxh+3),DUMMI(1),DUMMI(2),
     .                  XPOL(IY,IX+2)
         READ (30,*) YPOL(IY,dimxh+3),DUMMI(1),DUMMI(2),
     .                  YPOL(IY,IX+2)
         ENDIF
         IF (IY.EQ.dimyh) THEN
          READ (30,*) XPOL(IY,dimxh+3),XPOL(dimyh+1,dimxh+3),
     .                                XPOL(dimyh+1,IX+2),XPOL(IY,IX+2)
          READ (30,*) YPOL(IY,dimxh+3),YPOL(dimyh+1,dimxh+3),
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


C
*//GEOMD//
C=======================================================================
C          S U B R O U T I N E   G E O M D _ S O N N E T
C=======================================================================
      SUBROUTINE EIRENE_GEOMD_SONNET(NDXA,NDYA,NPLP,NR1ST,
     .                 PUX,PUY,PVX,PVY)
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CCONA
      USE EIRMOD_CGEOM
      IMPLICIT NONE
C
      REAL(DP), INTENT(OUT) :: PUX(*),PUY(*),PVX(*),PVY(*)
      INTEGER, INTENT(INOUT) :: NDXA,NDYA,NPLP,NR1ST

      character(200) :: zeile
      REAL(DP) :: br(0:ndxp,0:ndyp,4),bz(0:ndxp,0:ndyp,4)
C
C  GEOMETRY DATA: CELL VERTICES (LINDA ---> EIRENE)
      REAL(DP) ::
     R  X1(NDX),Y1(NDX),X2(NDX),Y2(NDX),X3(NDX),Y3(NDX),
     R  X4(NDX),Y4(NDX)
      REAL(DP) :: DX, DY
      INTEGER :: I0, I0E, IX, IY, I1, I2, I3, I4, IPART, I, J, NCUT
      CHARACTER(200) :: LINE
      LOGICAL :: LRDCUT

      ndxa=0
      ndya=0
      lrdcut = .false.

!pb      DO I = 1, 4
      DO 
        read (30,'(A200)') line
        if (.not.lrdcut) then
          call eirene_uppercase (line)
          i0 = index(line,' CUT')
          if (i0 > 0) then
            i0 = i0 + 4
            i1 = i0 + verify(line(i0:),' ')-1
            i2 = i1 + scan(line(i1+1:),' ')-1
            read (line(i1:i2),*) ncut
            do j=1,ncut
              i1 = i2 + verify(line(i2+1:),' ')-1
              i2 = i1 + scan(line(i1+1:),' ')-1
              read (line(i1:i2),*) npoint(2,j)
              npoint(2,j) = npoint(2,j) + j-1
              npoint(1,j+1) = npoint(2,j) + 1
            end do
            npoint(1,1) = 1
            ipart = ncut+1
            lrdcut = .true.
          end if 
        end if
        if (index(line,'=======') /= 0) exit
      END DO

1     continue
      read (30,'(a200)',end=99) zeile
      i0=index(zeile,'(')
      i0e=index(zeile,')')
      read (zeile(i0+1:i0e-1),*) ix,iy
      ndxa=max(ndxa,ix)
      ndya=max(ndya,iy)
      i1=index(zeile,': (')
      i2=index(zeile(i1+3:),')')+i1+2
      read (zeile(i1+3:i2-1),*) br(ix,iy,4),bz(ix,iy,4)
      i3=index(zeile(i2+1:),'(')+i2
      i4=i3+index(zeile(i3+1:),')')
      read (zeile(i3+1:i4-1),*) br(ix,iy,3),bz(ix,iy,3)

      read (30,'(a200)') zeile

      read (30,'(a200)') zeile
      i1=index(zeile,'(')
      i2=index(zeile,')')
      read (zeile(i1+1:i2-1),*) br(ix,iy,1),bz(ix,iy,1)
      i3=i2+index(zeile(i2+1:),'(')
      i4=i2+index(zeile(i2+1:),')')
      read (zeile(i3+1:i4-1),*) br(ix,iy,2),bz(ix,iy,2)

      read (30,*)
      goto 1


99    continue
      ndxa=ndxa-1
      ndya=ndya-1
C
!pb      DO 1015 IY=1,NDYA
!pb        DO 1014 IX=1,NDXA
!pb          X1(IX)=br(ix,iy,1)
!pb          Y1(IX)=bz(ix,iy,1)
!pb          X2(IX)=br(ix,iy,2)
!pb          Y2(IX)=bz(ix,iy,2)
!pb          X3(IX)=br(ix,iy,4)
!pb          Y3(IX)=bz(ix,iy,4)
!pb          X4(IX)=br(ix,iy,3)
!pb          Y4(IX)=bz(ix,iy,3)
!pb1014    CONTINUE
!pb        CALL MSHPROJ (X1,Y1,X2,Y2,X3,Y3,X4,Y4,PUX,PUY,PVX,PVY,NDXA,
!pb     .                NR1ST,IY)
!pb1015  CONTINUE
C
C SEARCH FOR THE CUTS
C
      IF (.NOT.LRDCUT) THEN
        IPART=1
        NPOINT(1,IPART)=1
        IY=1
        DO IX=1,NDXA
          DX=BR(IX+1,IY,1)-BR(IX,IY,2)
          DY=BZ(IX+1,IY,1)-BZ(IX,IY,2)
          IF (DX*DX+DY*DY.GT.EPS10) THEN
C CUT GEFUNDEN
            NPOINT(2,IPART)=IX+1+IPART-1
            IPART=IPART+1
            NPOINT(1,IPART)=IX+1+IPART-1
          ENDIF
        ENDDO
      END IF
      NPOINT(2,IPART)=NDXA+IPART
C
      NPLP=IPART
C
      DO IY=1,NDYA
        DO IPART=1,NPLP
          DO IX=NPOINT(1,IPART),NPOINT(2,IPART)-1
            XPOL(IY,IX)=BR(IX-(IPART-1),IY,1)
            YPOL(IY,IX)=BZ(IX-(IPART-1),IY,1)
          ENDDO
          XPOL(IY,NPOINT(2,IPART))=BR(NPOINT(2,IPART)-IPART,IY,2)
          YPOL(IY,NPOINT(2,IPART))=BZ(NPOINT(2,IPART)-IPART,IY,2)
        ENDDO
      ENDDO
C INTRODUCE OUTERMOST RADIAL POLYGON
      DO IPART=1,NPLP
        DO IX=NPOINT(1,IPART),NPOINT(2,IPART)-1
          XPOL(NDYA+1,IX)=BR(IX-(IPART-1),NDYA,4)
          YPOL(NDYA+1,IX)=BZ(IX-(IPART-1),NDYA,4)
        ENDDO
        XPOL(NDYA+1,NPOINT(2,IPART))=BR(NPOINT(2,IPART)-IPART,NDYA,3)
        YPOL(NDYA+1,NPOINT(2,IPART))=BZ(NPOINT(2,IPART)-IPART,NDYA,3)
      ENDDO
C
      DO J=1,NDYA+1
        DO I=1,NPOINT(2,NPLP)
          XPOL(J,I)=XPOL(J,I)*100.
          YPOL(J,I)=YPOL(J,I)*100.
        END DO
      END DO
C
      ndxa=npoint(2,nplp)-1

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
        CALL EIRENE_MSHPROJ (X1,Y1,X2,Y2,X3,Y3,X4,Y4,PUX,PUY,PVX,PVY,
     .                NDXA,NR1ST,IY)
1015  CONTINUE
c
C     do j=1,ndya+1
C       write (iunout,*)
C       write (iunout,*) 'in geomd polygon ',j
C       write (iunout,'(1p,6e12.4)') 
C    .        (xpol(j,i),ypol(j,i),i=1,npoint(2,nplp))
C     enddo
C
      RETURN
*//END GEOMD_SONNET//
      END


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
1     CONTINUE
      RETURN
      END


C
      SUBROUTINE EIRENE_INDMAP(FIELD,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .                  NCUTB,NCUTL,NPOINT,NPPLG)
C
C     INDEX MAPPING FOR BRAAMS DATA FIELDS. DATA IN DUMMY ZONES
C     (CUTS OR BOUNDARY ZONES) MAY BE NEEDED AND THUS ARE KEPT
C     AND DUBLICATED IN CASE NCUTL GT NCUTB
C
C     NCUTB= NUMBER OF CELLS IN IX DIRECTION PER CUT IN BRAAMS
C     NCUTL= NUMBER OF CELLS IN IX DIRECTION PER CUT IN LINDA (AND
C            THUS ALSO IN EIRENE) GEOMETRY

      USE EIRMOD_PRECISION
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NPOINT(2,*)
      INTEGER, INTENT(IN) :: NDX, NDY, NFL, NDXA, NDYA, NFLA, NCUTB,
     .                       NCUTL, NPPLG
      REAL(DP), INTENT(INOUT) :: FIELD(0:NDX+1,0:NDY+1,NFL),
     .                         DUMMY(0:NDX+1,0:NDY+1)
      INTEGER :: IX, IPART, IY, IF, IENDD, INB, IINID, IINIV, IENDV
C
C  LOOP FOR THE SPECIES
C
      DO 500 IF=1,NFLA
C
C  INITIALIZE DUMMY
C
        DO 10 IY=0,NDY+1
          DO 10 IX=0,NDX+1
10          DUMMY(IX,IY)=FIELD(IX,IY,IF)
C
C
C      NDX DIRECTION: IX=0: NOT MODIFIED
C                     IX=I(CUT): USE CUT VALUE
C                     IX=I(LAST X ZONE): MOVE TO NDXA+1
C
C  NPOINT(1,1)=1
C  NPOINT(2,NPPLG)=NDXA+1
C
        IF (NCUTB.LT.0) GOTO 990
        DO 211 IPART = 1,NPPLG
C  "VALID REGION"
          IINIV= NPOINT(1,IPART)
          IENDV= NPOINT(2,IPART)-1
C  "CUT REGION" AND LAST X ZONE IX = NDXA+1
          IF (IPART.LT.NPPLG) THEN
            IINID= NPOINT(2,IPART)
            IENDD= NPOINT(1,IPART+1)-1
            IF (IENDD-IINID+1.NE.NCUTL) GOTO 991
          ELSE
            IINID= NDXA+1
            IENDD= NDXA+1
          ENDIF
          DO 212 IY=0,NDYA+1
            DO 213 IX = IINIV,IENDV
              INB=IX-(IPART-1)*(NCUTL-NCUTB)
              DUMMY(IX,IY)=FIELD(INB,IY,IF)
213         CONTINUE
            DUMMY(IINID,IY) = FIELD(INB+1,IY,IF)
            IF (IENDD.NE.IINID) DUMMY(IENDD,IY) = FIELD(INB+NCUTB,IY,IF)
212       CONTINUE
211     CONTINUE
        DO 220 IY=0,NDYA+1
          DO 220 IX=0,NDXA+1
            FIELD(IX,IY,IF)=DUMMY(IX,IY)
220     CONTINUE
C
500   CONTINUE
      RETURN
C
990   CONTINUE
      WRITE (6,*) 'ERROR IN SUBR. INDMAP: THIS SUBR. IS VALID ONLY'
      WRITE (6,*) 'NCUTB>=0 BUT NCUTB = ',NCUTB
      CALL EIRENE_EXIT_OWN(1)
991   WRITE (6,*) 'ERROR IN SUBR. INDMAP: INCONSISTENCY IN NUMBER OF '
      WRITE (6,*) 'ZONES PER CUT FROM LINDA GEOMETRY DETECTED.  '
      WRITE (6,*) 'NCUTL = ',NCUTL, ' IENDD-IINID+1 = ',IENDD-IINID+1
      CALL EIRENE_EXIT_OWN(1)
      END


C
C
      SUBROUTINE EIRENE_INDMPI(FIELD,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .                  NCUTB,NCUTL,NPOINT,NPPLG,NSTR,ISTR)
C
C     INDEX MAPPING: INVERS TO SUBR. INDMAP
C
      USE EIRMOD_PRECISION
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NPOINT(2,*)
      INTEGER, INTENT(IN) :: NDX, NDY, NFL, NDXA, NDYA, NFLA, NCUTB,
     .                       NCUTL, NPPLG, NSTR, ISTR
      REAL(DP), INTENT(INOUT) :: FIELD(0:NDX+1,0:NDY+1,NFL,NSTR),
     .                         DUMMY(0:NDX+1,0:NDY+1)
      INTEGER :: IX, IY, IF, IENDD, IPART, INB, IINID, IINIV, IENDV
C
C  LOOP OVER THE SPECIES
C
      DO 500 IF=1,NFLA
C
C  INITIALIZE DUMMY
C
        DO 10 IY=0,NDY+1
          DO 10 IX=0,NDX+1
10          DUMMY(IX,IY)=0.
C
C
C      NDX DIRECTION
C
C  NPOINT(1,1)=1
C  NPOINT(2,NPPLG)=NDXA+1
C
        IF (NCUTB.LT.0) GOTO 990
        DO 211 IPART = 1,NPPLG
C  "VALID REGION"
          IINIV= NPOINT(1,IPART)
          IENDV= NPOINT(2,IPART)-1
C  "CUT REGION" AND LAST X ZONE IX = NDXA+1
          IF (IPART.LT.NPPLG) THEN
            IINID= NPOINT(2,IPART)
            IENDD= NPOINT(1,IPART+1)-1
            IF (IENDD-IINID+1.NE.NCUTL) GOTO 991
          ELSE
            IINID= NDXA+1
            IENDD= NDXA+1
          ENDIF
          DO 212 IY=0,NDYA+1
            DO 213 IX = IINIV,IENDV
              INB=IX-(IPART-1)*(NCUTL-NCUTB)
              DUMMY(INB,IY)=FIELD(IX,IY,IF,ISTR)
213         CONTINUE
            DUMMY(INB+1,IY)=FIELD(IINID,IY,IF,ISTR)
            IF (IENDD.NE.IINID)
     .          DUMMY(INB+NCUTB,IY)=FIELD(IENDD,IY,IF,ISTR)
212       CONTINUE
211     CONTINUE
        DO 220 IY=0,NDYA+1
          DO 220 IX=0,NDXA+1
            FIELD(IX,IY,IF,ISTR)=DUMMY(IX,IY)
220     CONTINUE
C
500   CONTINUE
      RETURN
C
990   CONTINUE
      WRITE (6,*) 'ERROR IN SUBR. INDMPI: THIS SUBR. IS VALID ONLY'
      WRITE (6,*) 'NCUTB>=0 BUT NCUTB = ',NCUTB
      CALL EIRENE_EXIT_OWN(1)
991   WRITE (6,*) 'ERROR IN SUBR. INDMPI: INCONSISTENCY IN NUMBER OF'
      WRITE (6,*) 'ZONES PER CUT FROM LINDA GEOMETRY DETECTED. '
      WRITE (6,*) 'NCUTL = ',NCUTL, ' IENDD-IINID+1 = ',IENDD-IINID+1
      CALL EIRENE_EXIT_OWN(1)
      END


C
*//PLASM//
C=======================================================================
C          S U B R O U T I N E   P L A S M
C=======================================================================
      SUBROUTINE EIRENE_PLASM(KARD,NDIMX,NDIMY,NDIMF,N,M,NF,DUMMY)

      USE EIRMOD_PRECISION
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: KARD, NDIMX, NDIMY, NDIMF, N, M, NF
      REAL(DP), INTENT(INOUT) :: DUMMY(0:N+1,0:M+1,NF)
      INTEGER :: ND1, LIM, IF, III, IX, IY

      ND1 = NDIMX + 2
      LIM = (ND1/5)*5 - 4
      DO    110  IF = 1,NDIMF
      DO    110  IY = 0,NDIMY+1
      DO    100  IX = 1,LIM,5
100     READ(KARD,910,END=500) (DUMMY(-1+IX-1+III,IY,IF),III = 1,5)
        IF( (LIM+4).EQ.ND1 )     GOTO 110
        READ(KARD,910,END=500) (DUMMY(-1+IX,IY,IF),IX = LIM+5,ND1)
110   CONTINUE
500   RETURN
910   FORMAT(5(E16.8))
!910   FORMAT(5(E23.16))
*//END PLASM//
      END


C
C
*//NEUTR//
C=======================================================================
C          S U B R O U T I N E   N E U T R
C=======================================================================
      SUBROUTINE EIRENE_NEUTR(KARD,NDIMX,NDIMY,NDIMF,DUMMY,LDMX,LDMY,
     .                        LDMF,LDNS,IS)

      USE EIRMOD_PRECISION
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: KARD, NDIMX, NDIMY, NDIMF, LDMX, LDMY,
     .                       LDMF, LDNS, IS
      REAL(DP), INTENT(IN) :: DUMMY(0:LDMX+1,0:LDMY+1,LDMF,LDNS)
      INTEGER :: ND1, LIM, IX, IY, III, IF
C
      ND1 = NDIMX
      LIM = (ND1/5)*5 - 4
      DO  500  IF = 1,NDIMF
        DO  110  IY = 1,NDIMY
          DO  100  IX = 1,LIM,5
  100     WRITE(KARD,910) (DUMMY(IX-1+III,IY,IF,IS),III = 1,5)
          IF( (LIM+4).EQ.ND1 )   GOTO 110
          WRITE(KARD,910) (DUMMY(IX,IY,IF,IS),IX = LIM+5,ND1)
  110   CONTINUE
  500 CONTINUE
      RETURN
  910 FORMAT(5(E16.8))
*//END NEUTR//
      END


!pb  121206  check if tallies are available before using them

      SUBROUTINE EIRENE_SAVE_TALLIES (ISTRAI)
C
C  SAVE EIRENE TALLIES, SCALE PER UNIT FLUX (AMP), ON COMMON BRASCL
C  WTOTP IS NEGATIVE IN EIRENE (SINK FOR IONS)
C  ALL STRATA WHICH ARE NOT SPECIFIED BY INPUT BLOCK 14 (FROM
C  PLASMA CODE DATA) ARE NOT RESCALED HERE
C

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_BRASPOI
      USE EIRMOD_CCOUPL
      USE EIRMOD_COUTAU
      USE EIRMOD_COMUSR
      USE EIRMOD_CGRID
      USE EIRMOD_CESTIM

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ISTRAI
      REAL(DP) :: FLXI
      INTEGER :: IN, IATM, IMOL, IPLS, IION, ICPV

      TYPE(CELLSIM), POINTER :: CPSIM
      TYPE(CELLMUL), POINTER :: CPMUL

      IF (ISTRAI.LE.NTARGI.AND.WTOTP(0,ISTRAI).NE.0.) THEN
         FLXI=-1._DP/WTOTP(0,ISTRAI)
      ELSEIF (ISTRAI.LE.NTARGI.AND.WTOTP(0,ISTRAI).EQ.0.) THEN
         RETURN
      ELSEIF (ISTRAI.GT.NTARGI) THEN
         FLXI=1._DP
      ENDIF

      CALL EIRENE_FREE_SIMARR(ISTRAI)
      CALL EIRENE_FREE_MULARR(ISTRAI)

      DO IPLS=1,NPLSI
        DO IN=1,NSBOX_TAL
          IF (LPAPL) THEN 
	  IF (PAPL(IPLS,IN) .NE. 0.D0) THEN
!pb            ALLOCATE(CPMUL)
            CPMUL => EIRENE_NEW_MULARR()
            CPMUL%IART = IPLS
            CPMUL%ICM = IN
            CPMUL%VALUEM = PAPL(IPLS,IN)*FLXI
            CPMUL%NXTMUL => PAPLS(ISTRAI)%PMUL
            PAPLS(ISTRAI)%PMUL => CPMUL
          ENDIF
          ENDIF
          IF (LPMPL) THEN 
          IF (PMPL(IPLS,IN) .NE. 0.D0) THEN
!PB            ALLOCATE(CPMUL)
            CPMUL => EIRENE_NEW_MULARR()
            CPMUL%IART = IPLS
            CPMUL%ICM = IN
            CPMUL%VALUEM = PMPL(IPLS,IN)*FLXI
            CPMUL%NXTMUL => PMPLS(ISTRAI)%PMUL
            PMPLS(ISTRAI)%PMUL => CPMUL
          ENDIF
          ENDIF
          IF (LPIPL) THEN 
          IF (PIPL(IPLS,IN) .NE. 0.D0) THEN
!PB            ALLOCATE(CPMUL)
            CPMUL => EIRENE_NEW_MULARR()
            CPMUL%IART = IPLS
            CPMUL%ICM = IN
            CPMUL%VALUEM = PIPL(IPLS,IN)*FLXI
            CPMUL%NXTMUL => PIPLS(ISTRAI)%PMUL
            PIPLS(ISTRAI)%PMUL => CPMUL
          ENDIF
          ENDIF
          IF (LMAPL) THEN 
          IF (MAPL(IPLS,IN) .NE. 0.D0) THEN
!pb            ALLOCATE(CPMUL)
            CPMUL => EIRENE_NEW_MULARR()
            CPMUL%IART = IPLS
            CPMUL%ICM = IN
            CPMUL%VALUEM = MAPL(IPLS,IN)*FLXI
            CPMUL%NXTMUL => MAPLS(ISTRAI)%PMUL
            MAPLS(ISTRAI)%PMUL => CPMUL
          ENDIF
          ENDIF
          IF (LMMPL) THEN 
          IF (MMPL(IPLS,IN) .NE. 0.D0) THEN
!PB            ALLOCATE(CPMUL)
            CPMUL => EIRENE_NEW_MULARR()
            CPMUL%IART = IPLS
            CPMUL%ICM = IN
            CPMUL%VALUEM = MMPL(IPLS,IN)*FLXI
            CPMUL%NXTMUL => MMPLS(ISTRAI)%PMUL
            MMPLS(ISTRAI)%PMUL => CPMUL
          ENDIF
          ENDIF
          IF (LMIPL) THEN 
          IF (MIPL(IPLS,IN) .NE. 0.D0) THEN
!PB            ALLOCATE(CPMUL)
            CPMUL => EIRENE_NEW_MULARR()
            CPMUL%IART = IPLS
            CPMUL%ICM = IN
            CPMUL%VALUEM = MIPL(IPLS,IN)*FLXI
            CPMUL%NXTMUL => MIPLS(ISTRAI)%PMUL
            MIPLS(ISTRAI)%PMUL => CPMUL
          ENDIF
          ENDIF
          IF (LMPHPL) THEN 
          IF (MPHPL(IPLS,IN) .NE. 0.D0) THEN
!PB            ALLOCATE(CPMUL)
            CPMUL => EIRENE_NEW_MULARR()
            CPMUL%IART = IPLS
            CPMUL%ICM = IN
            CPMUL%VALUEM = MPHPL(IPLS,IN)*FLXI
            CPMUL%NXTMUL => MPHPLS(ISTRAI)%PMUL
            MPHPLS(ISTRAI)%PMUL => CPMUL
          ENDIF
          ENDIF
        ENDDO
      ENDDO

      DO IN=1,NSBOX_TAL
	IF (LEAEL) THEN
        IF (EAEL(IN) .NE. 0.D0) THEN
!PB          ALLOCATE(CPSIM)
          CPSIM => EIRENE_NEW_SIMARR()
          CPSIM%ICS = IN
          CPSIM%VALUES = EAEL(IN)*FLXI
          CPSIM%NXTSIM => EAELS(ISTRAI)%PSIM
          EAELS(ISTRAI)%PSIM => CPSIM
        ENDIF
        ENDIF
        IF (LEMEL) THEN 
        IF (EMEL(IN) .NE. 0.D0) THEN
!PB          ALLOCATE(CPSIM)
          CPSIM => EIRENE_NEW_SIMARR()
          CPSIM%ICS = IN
          CPSIM%VALUES = EMEL(IN)*FLXI
          CPSIM%NXTSIM => EMELS(ISTRAI)%PSIM
          EMELS(ISTRAI)%PSIM => CPSIM
        ENDIF
        ENDIF
        IF (LEIEL) THEN 
        IF (EIEL(IN) .NE. 0.D0) THEN
!PB          ALLOCATE(CPSIM)
          CPSIM => EIRENE_NEW_SIMARR()
          CPSIM%ICS = IN
          CPSIM%VALUES = EIEL(IN)*FLXI
          CPSIM%NXTSIM => EIELS(ISTRAI)%PSIM
          EIELS(ISTRAI)%PSIM => CPSIM
        ENDIF
        ENDIF
        IF (LEAPL) THEN 
        IF (EAPL(IN) .NE. 0.D0) THEN
!PB          ALLOCATE(CPSIM)
          CPSIM => EIRENE_NEW_SIMARR()
          CPSIM%ICS = IN
          CPSIM%VALUES = EAPL(IN)*FLXI
          CPSIM%NXTSIM => EAPLS(ISTRAI)%PSIM
          EAPLS(ISTRAI)%PSIM => CPSIM
        ENDIF
        ENDIF
        IF (LEMPL) THEN 
        IF (EMPL(IN) .NE. 0.D0) THEN
!PB          ALLOCATE(CPSIM)
          CPSIM => EIRENE_NEW_SIMARR()
          CPSIM%ICS = IN
          CPSIM%VALUES = EMPL(IN)*FLXI
          CPSIM%NXTSIM => EMPLS(ISTRAI)%PSIM
          EMPLS(ISTRAI)%PSIM => CPSIM
        ENDIF
        ENDIF
        IF (LEIPL) THEN 
        IF (EIPL(IN) .NE. 0.D0) THEN
!PB          ALLOCATE(CPSIM)
          CPSIM => EIRENE_NEW_SIMARR()
          CPSIM%ICS = IN
          CPSIM%VALUES = EIPL(IN)*FLXI
          CPSIM%NXTSIM => EIPLS(ISTRAI)%PSIM
          EIPLS(ISTRAI)%PSIM => CPSIM
        ENDIF
	ENDIF
      ENDDO

      DO IATM=1,NATMI
        DO IN=1,NSBOX_TAL
	  IF (LPDENA) THEN
          IF (PDENA(IATM,IN) .NE. 0.D0) THEN
!PB            ALLOCATE(CPMUL)
            CPMUL => EIRENE_NEW_MULARR()
            CPMUL%IART = IATM
            CPMUL%ICM = IN
            CPMUL%VALUEM = PDENA(IATM,IN)*FLXI
            CPMUL%NXTMUL => PDENAS(ISTRAI)%PMUL
            PDENAS(ISTRAI)%PMUL => CPMUL
          ENDIF
          ENDIF
          IF (LEDENA) THEN 
          IF (EDENA(IATM,IN) .NE. 0.D0) THEN
!PB            ALLOCATE(CPMUL)
            CPMUL => EIRENE_NEW_MULARR()
            CPMUL%IART = IATM
            CPMUL%ICM = IN
            CPMUL%VALUEM = EDENA(IATM,IN)*FLXI
            CPMUL%NXTMUL => EDENAS(ISTRAI)%PMUL
            EDENAS(ISTRAI)%PMUL => CPMUL
          ENDIF
	  ENDIF
        ENDDO
      ENDDO

      DO IMOL=1,NMOLI
        DO IN=1,NSBOX_TAL
	  IF (LPDENM) THEN
          IF (PDENM(IMOL,IN) .NE. 0.D0) THEN
!PB            ALLOCATE(CPMUL)
            CPMUL => EIRENE_NEW_MULARR()
            CPMUL%IART = IMOL
            CPMUL%ICM = IN
            CPMUL%VALUEM = PDENM(IMOL,IN)*FLXI
            CPMUL%NXTMUL => PDENMS(ISTRAI)%PMUL
            PDENMS(ISTRAI)%PMUL => CPMUL
          ENDIF
          ENDIF
        ENDDO
      ENDDO

      DO IION=1,NIONI
        DO IN=1,NSBOX_TAL
	  IF (LPDENI) THEN
          IF (PDENI(IION,IN) .NE. 0.D0) THEN
!PB            ALLOCATE(CPMUL)
            CPMUL => EIRENE_NEW_MULARR()
            CPMUL%IART = IION
            CPMUL%ICM = IN
            CPMUL%VALUEM = PDENI(IION,IN)*FLXI
            CPMUL%NXTMUL => PDENIS(ISTRAI)%PMUL
            PDENIS(ISTRAI)%PMUL => CPMUL
          ENDIF
          ENDIF
        ENDDO
      ENDDO

      DO ICPV=1,NCPVI
        DO IN=1,NSBOX_TAL
	  IF (LCOPV) THEN
          IF (COPV(ICPV,IN) .NE. 0.D0) THEN
!PB            ALLOCATE(CPMUL)
            CPMUL => EIRENE_NEW_MULARR()
            CPMUL%IART = ICPV
            CPMUL%ICM = IN
            CPMUL%VALUEM = COPV(ICPV,IN)*FLXI
            CPMUL%NXTMUL => COPVS(ISTRAI)%PMUL
            COPVS(ISTRAI)%PMUL => CPMUL
          ENDIF
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END


C
C
C
C
      SUBROUTINE EIRENE_EIRSRT(LSTOP,LTIME,DELTAT,FLUXES,
     .                  B2BRM,B2RD,B2Q,B2VP)

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_BRASPOI
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CINIT
      USE EIRMOD_CPOLYG
      USE EIRMOD_CGRID
      USE EIRMOD_CSPEZ
      USE EIRMOD_CZT1
      USE EIRMOD_CTRCEI
      USE EIRMOD_CCOUPL
      USE EIRMOD_CGEOM
      USE EIRMOD_CSDVI
      USE EIRMOD_CSDVI_BGK
      USE EIRMOD_CSDVI_COP
      USE EIRMOD_COMPRT
      USE EIRMOD_COMNNL
      USE EIRMOD_COMSOU
      USE EIRMOD_COUTAU
      USE EIRMOD_COMXS
      USE EIRMOD_CSPEI
      USE EIRMOD_BRASCL
      USE EIRMOD_CTRIG

      IMPLICIT NONE
C
      REAL(DP), INTENT(IN) :: FLUXES(*)
      REAL(DP), INTENT(IN) :: DELTAT, B2BRM, B2RD, B2Q, B2VP
      LOGICAL, INTENT(IN) :: LSTOP, LTIME

      REAL(DP), ALLOCATABLE, SAVE :: FLUXS(:)
      REAL(DP) :: EIRENE_FTABEI1, EIRENE_FEELEI1, FLXI, ESIG, 
     .            EIRENE_RESET_SECOND, DUMMY,
     .          EIRENE_SECOND_OWN, DTIMVO
      INTEGER :: IN, IAEI, IRDS, IIDS, ICPV, IMDS, IFIRST, K, JC, NDXY,
     .           J, IRC, NREC10, NREC11, ITNR, IPLSTI, IST_RATE, IST,
     .           IFRSTR, ISTH, ISTNEW, ISTIN
      REAL(DP), ALLOCATABLE :: OUTAU(:)
      INTEGER, ALLOCATABLE :: IHELP(:)
      LOGICAL :: LSTP, LLST, LPLASM, NLSRON_SAVE(NSTRA)
C
      TYPE(CELLSIM), POINTER :: CPSIM
      TYPE(CELLMUL), POINTER :: CPMUL
      TYPE(RATE_STORE), POINTER :: RTIS
C
C
      SAVE
      DATA IFIRST/0/
C
      IF (LTIME) THEN
C
        B2BREM=B2BRM
        B2RAD=B2RD
        B2QIE=B2Q
        B2VDP=B2VP
        DUMMY=EIRENE_RESET_SECOND()
        IF(IFIRST.EQ.0) THEN
C
          CALL GRSTRT(35,8)
C
C  READ FORMATTED INPUT FILE IUNIN
C  AND RUN EIRENE FOR ONE TIME-CYCLE: ITIMV=1
C  WITH OR WITHOUT INITIAL DISTRIBUTION ON FILE FT15 (NFILE-J FLAG)
C  AS FINAL STRATUM
C  EXPECT PLASMA DATA ON FORT.31 (NLPLAS=.FALSE.)
C
          CALL EIRENE_EIRENE(DELTAT,.FALSE.,.FALSE.,1,.TRUE.)
C
C  EIRENE RUN DONE. CENSUS ARRAY WRITTEN
C  NOW ITIMV=ITIMV+1, NLPLAS=.TRUE.
C
          IF (.NOT.NLPLAS) THEN
            WRITE (6,*) 'INCONSISTENT COUPLING '
            WRITE (6,*) 'LTIME=TRUE, BUT NTIME = ', NTIME
            CALL EIRENE_EXIT_OWN(1)
          ENDIF
          IF (.NOT.ALLOCATED(FLUXS)) ALLOCATE (FLUXS(NSTRA))
          DO 3 ISTRA=1,NSTRAI
            FLUXS(ISTRA)=FLUX(ISTRA)
3         CONTINUE
          IFIRST=1
        ELSE
C
C  NOW: NLPLAS=.TRUE., I.E., PLASMA DATA EXPECTED ON BRAEIR
C  NOW: ITIMV=ITIMV+1
C  BUT: COMMON BRAEIR REDONE IN EXTERNAL CODE.
C  REACTIVATE INDEX MAPPING, EVEN WITHOUT READING INPUT BLOCK 14 AGAIN
          NCUTB_SAVE=NCUTB
C
          DTIMVO=DTIMV
          DTIMVN=DELTAT
C
C-----------------------------------------------------------------------
C
C  STRATA 1 TO NTARGI ARE SCALED IN PLASMA CODE  (RECYCLING STRATA)
C
C     RETURN TO PLASMA CODE THE PROFILES PER UNIT SOURCE STRENGTH
C     IE. THE PROFILES ARE SCALED BY 1./FLUX(ISTRA) BEFORE RETURN
C
C  STRATA NTARGI+1 TO NSTRAI-1  ARE SCALED BY EIRENE
C
C     (EG. GAS PUFF, VOLUME RECOMBINATION, ETC.)
C     THEY MAY BE RESCALED BY PLASMA CODE FACTORS: FLUXES(ISTRA)
C     RETURN TO PLASMA CODE THE PROFILES SCALED WITH
C     SOURCE STRENGTH: FLUX(ISTRA) (AMP)
C
C  STRATUM NSTRAI IS RESCALED WITH RATIO OF OLD TO NEW TIMESTEP
C
C     RETURN TO PLASMA CODE THE PROFILES WITH FLUX(ISTRA) (AMP)
C
          DO ISTRA=NTARGI+1,NSTRAI-1
            IF (FLUXES(ISTRA).NE.0.) THEN
              FLUX(ISTRA)=FLUXS(ISTRA)*FLUXES(ISTRA)*ELCHA
            ELSE
              FLUX(ISTRA)=FLUXS(ISTRA)
            ENDIF
          ENDDO
C
          IF (DTIMVN.NE.DTIMVO) THEN
            FLUX(NSTRAI)=FLUX(NSTRAI)*DTIMVO/DTIMVN
C
            WRITE (6,*) 'FLUX IS RESCALED BY DTIMV_OLD/DTIMV_NEW '
            CALL EIRENE_MASR1('FLUX    ',FLUX(NSTRAI))
            CALL EIRENE_LEER(1)
          ENDIF
C
C-----------------------------------------------------------------------
C
          DTIMV=DTIMVN
C
C  RUN EIRENE ON TIMESTEP DTIMV
C  THEN CALL INTERFACING ROUTINE AT ENTRY IF3COP (FROM EIRENE MAIN)
C
          IITER=1
          IPRNLI=0
          NLSRON=.TRUE.
          CALL EIRENE_EIRENE_COUPLE (LSTOP,1,.TRUE.)
          IF (LSTOP) THEN
            CALL GREND
          ENDIF
        ENDIF
        CALL EIRENE_LEER(2)
        WRITE(*,*) 'EIRENE USED ',EIRENE_SECOND_OWN(),' CPU SECONDS'
        CALL EIRENE_LEER(2)
C
        RETURN
C
      ELSEIF (.NOT.LTIME) THEN

!swpb for multiprocessor calculation
        DUMMY=EIRENE_RESET_SECOND()
C
        IF (IFIRST.GE.1) GOTO 10000

        CALL GRSTRT(35,8)
C
C  READ FORMATTED INPUT FILE IUNIN
C  AND RUN EIRENE FOR ONE TIME-CYCLE: ITIMV=1
C  WITH OR WITHOUT INITIAL DISTRIBUTION ON FILE FT15 (NFILE-J FLAG)
C  AS FINAL STRATUM
C  EXPECT PLASMA DATA ON FORT.31 (NLPLAS=.FALSE.)
C
        B2BREM=B2BRM
        B2RAD=B2RD
        B2QIE=B2Q
        B2VDP=B2VP
        LPLASM=.FALSE.
        LLST=LSTOP
        ITNR=1

 10     CONTINUE

        CALL EIRENE_EIRENE(DELTAT,LPLASM,LLST,ITNR,.TRUE.)
C
C  IN THIS CALL TO EIRENE ALREADY IF3COP IS CALLED FOR EACH STRATUM
C  THOSE WITH NLSRON = TRUE  HAVE BEEN RECOMPUTED BY EIRENE
C  THOSE WITH NLSRON = FALSE HAVE BEEN SHORT-CYCLED
C  AT IFIRST   =0: ALL NLSRON=TRUE
C  AT IFIRST.GE.1: FIRST A SHORT CYCLE TEST IS DONE, AND NLSRON IS FOUND
C
        IF (.NOT.LLST) THEN

        IF (IFIRST.GE.1) NLSRON = NLSRON_SAVE

csw 12apr2011       NDXY=(NDXA-1)*NR1ST+NDYA
        NDXY=NTRII
C
        CALL EIRENE_ALLOC_BRASCL

! find new calculated stratum with smallest number
        DO IST = 1, NSTRAI
          IF (NLSRON(IST)) THEN
            IFRSTR = IST
            EXIT
          END IF
        END DO

! determine index of rate storage which has been used in the last iteration 
        IST_RATE = ITS(IFRSTR)

! reduce counters of rate storages for all new calculated strata
        DO IST = 1, NSTRAI
          IF (NLSRON(IST)) THEN
            ISTIN = ITS(IST)
            ITS_COUNT(ISTIN) = ITS_COUNT(ISTIN) - 1
          END IF
        END DO

! check how often storage IST_RATE is still used
        ISTH = 0
        IF (IST_RATE > 0) ISTH= ITS_COUNT(IST_RATE)
        
        IF (ISTH < 1) THEN
! rate storage can be used again
        ELSE
! rate storage still in use, look for an empty slot
          ISTNEW = MINLOC(ITS_COUNT,DIM=1)
          IF (ITS_COUNT(ISTNEW) > 0) THEN
            WRITE (IUNOUT,*) ' PROBLEM IN EIRSRT '
            WRITE (IUNOUT,*) ' ITS_COUNT > 0 '
            WRITE (IUNOUT,*) ' ITS_COUNT ',ITS_COUNT
            CALL EIRENE_EXIT_OWN(1)
          END IF
          IST_RATE = ISTNEW
        END IF
        
        WHERE (NLSRON) 
          ITS = IST_RATE
        END WHERE

        ITS_COUNT(IST_RATE) = COUNT(NLSRON)
C
        CALL EIRENE_ALLOC_RATE_ARRAY(IST_RATE)
        CALL EIRENE_INIT_BRASCL1(IST_RATE)

        RTIS => RTS(IST_RATE)%RTA
C
C  INITIAL: ATOMS, EI-PROCESSES
C
        DO 21 IATM=1,NATMI
        DO 21 IPLS=1,NPLSI
        DO 21 IAEI=1,NAEII(IATM)
          IRDS=LGAEI(IATM,IAEI)
          IF (PPLDS(IRDS,IPLS).EQ.0.) GOTO 21
          DO 22 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              RTIS%SPLODA(IN,IATM,IPLS)=RTIS%SPLODA(IN,IATM,IPLS)+
     .                        TABDS1(IRDS,IN)*PPLDS(IRDS,IPLS)
            ELSE
              RTIS%SPLODA(IN,IATM,IPLS)=RTIS%SPLODA(IN,IATM,IPLS)+
     .                        EIRENE_FTABEI1(IRDS,IN)*PPLDS(IRDS,IPLS)
            END IF
22        CONTINUE
21      CONTINUE
        DO 23 IPLS=1,NPLSI
          IPLSTI = MPLSTI(IPLS)
          DO 24 IN=1,NDXY
            RTIS%SEIODA(IN,IPLS)=DIIN(IPLS,IN)*
     .                      (1.5*TIIN(IPLSTI,IN)+EDRIFT(IPLS,IN))
24        CONTINUE
23      CONTINUE
C
        DO 25 IATM=1,NATMI
        DO 25 IAEI=1,NAEII(IATM)
          IRDS=LGAEI(IATM,IAEI)
          DO 25 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              RTIS%SEEODA(IN,IATM)=RTIS%SEEODA(IN,IATM)+EELDS1(IRDS,IN)*
     .                                        TABDS1(IRDS,IN)
            ELSE
              RTIS%SEEODA(IN,IATM)=RTIS%SEEODA(IN,IATM)+
     .                                        EIRENE_FEELEI1(IRDS,IN)*
     .                                        EIRENE_FTABEI1(IRDS,IN)
            END IF
25      CONTINUE
C
C  INITIAL: TEST IONS, EI-PROCESSES
C
        DO 26 IION=1,NIONI
        DO 26 IIDS=1,NIDSI(IION)
          IRDS=LGIEI(IION,IIDS)
          DO 26 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              RTIS%SEEODI(IN,IION)=RTIS%SEEODI(IN,IION)+EELDS1(IRDS,IN)*
     .                                        TABDS1(IRDS,IN)
            ELSE
              RTIS%SEEODI(IN,IION)=RTIS%SEEODI(IN,IION)+
     .                                        EIRENE_FEELEI1(IRDS,IN)*
     .                                        EIRENE_FTABEI1(IRDS,IN)
            END IF
26      CONTINUE
C
        DO 27 IION=1,NIONI
        DO 27 IPLS=1,NPLSI
        DO 27 IIDS=1,NIDSI(IION)
          IRDS=LGIEI(IION,IIDS)
          IF (PPLDS(IRDS,IPLS).EQ.0.) GOTO 27
          DO 28 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              RTIS%SPLODI(IN,IION,IPLS)=RTIS%SPLODI(IN,IION,IPLS)+
     .                             TABDS1(IRDS,IN)*PPLDS(IRDS,IPLS)
            ELSE
              RTIS%SPLODI(IN,IION,IPLS)=RTIS%SPLODI(IN,IION,IPLS)+
     .                          EIRENE_FTABEI1(IRDS,IN)*PPLDS(IRDS,IPLS)
            ENDIF
28        CONTINUE
27      CONTINUE
C
        DO 29 IION=1,NIONI
        DO 29 IPLS=1,NPLSI
        DO 29 IIDS=1,NIDSI(IION)
          IRDS=LGIEI(IION,IIDS)
          ESIG=EPLDS(IRDS,2)
          DO 30 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              RTIS%SEIODI(IN,IION)=RTIS%SEIODI(IN,IION)+
     .                        TABDS1(IRDS,IN)*ESIG
            ELSE
              RTIS%SEIODI(IN,IION)=RTIS%SEIODI(IN,IION)+
     .                        EIRENE_FTABEI1(IRDS,IN)*ESIG
            END IF
30        CONTINUE
29      CONTINUE
C
C
C  INITIAL: MOLECULES, EI-PROCESSES
C
        DO 35 IMOL=1,NMOLI
        DO 35 IMDS=1,NMDSI(IMOL)
          IRDS=LGMEI(IMOL,IMDS)
          DO 35 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              RTIS%SEEODM(IN,IMOL)=RTIS%SEEODM(IN,IMOL)+EELDS1(IRDS,IN)*
     .                                        TABDS1(IRDS,IN)
            ELSE
              RTIS%SEEODM(IN,IMOL)=RTIS%SEEODM(IN,IMOL)+
     .                                        EIRENE_FEELEI1(IRDS,IN)*
     .                                        EIRENE_FTABEI1(IRDS,IN)
            ENDIF
35      CONTINUE
C
        DO 47 IMOL=1,NMOLI
        DO 47 IPLS=1,NPLSI
        DO 47 IMDS=1,NMDSI(IMOL)
          IRDS=LGMEI(IMOL,IMDS)
          IF (PPLDS(IRDS,IPLS).EQ.0.) GOTO 47
          DO 48 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              RTIS%SPLODM(IN,IMOL,IPLS)=RTIS%SPLODM(IN,IMOL,IPLS)+
     .                             TABDS1(IRDS,IN)*PPLDS(IRDS,IPLS)
            ELSE
              RTIS%SPLODM(IN,IMOL,IPLS)=RTIS%SPLODM(IN,IMOL,IPLS)+
     .                  EIRENE_FTABEI1(IRDS,IN)*PPLDS(IRDS,IPLS)
            END IF
48        CONTINUE
47      CONTINUE
C
        DO 49 IMOL=1,NMOLI
        DO 49 IPLS=1,NPLSI
        DO 49 IMDS=1,NMDSI(IMOL)
          IRDS=LGMEI(IMOL,IMDS)
          ESIG=EPLDS(IRDS,2)
          DO 50 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              RTIS%SEIODM(IN,IMOL)=RTIS%SEIODM(IN,IMOL)+
     .                        TABDS1(IRDS,IN)*ESIG
            ELSE
              RTIS%SEIODM(IN,IMOL)=RTIS%SEIODM(IN,IMOL)+
     .                        EIRENE_FTABEI1(IRDS,IN)*ESIG
            END IF
50        CONTINUE
49      CONTINUE

        END IF
C
C
        B2BREM=B2BRM
        B2RAD=B2RD
        B2QIE=B2Q
        B2VDP=B2VP
C
        IFIRST=IFIRST+1

        IF (LSTOP) THEN
          CALL EIRENE_DEALLOC_COMUSR
          CALL EIRENE_DEALLOC_CESTIM
          CALL EIRENE_DEALLOC_BRASCL
          CALL EIRENE_DEALLOC_BRASPOI

          CALL GREND
        END IF

        RETURN
C
C  NOT THE FIRST CALL IN THIS CYCLE: CHECK: SHORT LOOP CORRECTION
C                                           OR FULL EIRENE, FOR EACH
C                                           STRATUM INDIVIDUALLY
10000   CONTINUE
C
        LSTP = LSTOP
        NCUTB_SAVE=NCUTB

        CALL EIRENE_ALLOC_BCKGRND

        CALL EIRENE_INTER1
C
        CALL EIRENE_PLASMA
C
        CALL EIRENE_PLASMA_DERIV(0)
C
        CALL EIRENE_SETAMD(2)
C
C  IN PLASMA_DERIV THE BACKGROUND PLASMA STATE HAS BEEN 
C  WRITTEN TO FORT.13
C  NFILEL HAS BEEN CHANGED TO NFILEL = 3 OR 9
C  ==> PLASMA AND REACTION DATA ARE READ IN SUBR. INPUT
C  NOW SAVE REACTION DATA AS WELL IN ORDER TO HAVE A 
C  CONSISTENT PLASMA STATE ON FORT.13
C  
      IF ((NFILEL >=1) .AND. (NFILEL <=5)) THEN
         NFILEL=3
         CALL EIRENE_WRPLAM(TRCFLE,0)
      ELSE IF (NFILEL > 5) THEN
         NFILEL=9
         CALL EIRENE_WRPLAM_XDR(TRCFLE,0)
      END IF

C
        CALL EIRENE_ALLOC_BRASCL
        CALL EIRENE_INIT_BRASCL2
C
C  NEW: ATOMS, EI PROCESSES
C
        DO 101 IATM=1,NATMI
        DO 101 IPLS=1,NPLSI
          DO 102 IAEI=1,NAEII(IATM)
            IRDS=LGAEI(IATM,IAEI)
            IF (PPLDS(IRDS,IPLS).EQ.0.) GOTO 101
            DO 102 IN=1,NDXY
              IF (NSTORDR >= NRAD) THEN
                SPLNWA(IN,IATM,IPLS)=SPLNWA(IN,IATM,IPLS)+
     .                          TABDS1(IRDS,IN)*PPLDS(IRDS,IPLS)
              ELSE
                SPLNWA(IN,IATM,IPLS)=SPLNWA(IN,IATM,IPLS)+
     .                          EIRENE_FTABEI1(IRDS,IN)*PPLDS(IRDS,IPLS)
              END IF
102       CONTINUE
101     CONTINUE
C
        DO 103 IPLS=1,NPLSI
          IPLSTI = MPLSTI(IPLS)
          DO 104 IN=1,NDXY
            SEINWA(IN,IPLS)=DIIN(IPLS,IN)*
     .                      (1.5*TIIN(IPLSTI,IN)+EDRIFT(IPLS,IN))
104       CONTINUE
103     CONTINUE
C
        DO 105 IATM=1,NATMI
          DO 105 IAEI=1,NAEII(IATM)
            IRDS=LGAEI(IATM,IAEI)
            DO 105 IN=1,NDXY
              IF (NSTORDR >= NRAD) THEN
                SEENWA(IN,IATM)=SEENWA(IN,IATM)+EELDS1(IRDS,IN)*
     .                                          TABDS1(IRDS,IN)
              ELSE
                SEENWA(IN,IATM)=SEENWA(IN,IATM)+EIRENE_FEELEI1(IRDS,IN)*
     .                                          EIRENE_FTABEI1(IRDS,IN)
              END IF
105     CONTINUE
C
C  NEW: TEST IONS, EI PROCESSES
C
        DO 106 IION=1,NIONI
          DO 106 IIDS=1,NIDSI(IION)
            IRDS=LGIEI(IION,IIDS)
            DO 106 IN=1,NDXY
              IF (NSTORDR >= NRAD) THEN
                SEENWI(IN,IION)=SEENWI(IN,IION)+EELDS1(IRDS,IN)*
     .                                          TABDS1(IRDS,IN)
              ELSE
                SEENWI(IN,IION)=SEENWI(IN,IION)+EIRENE_FEELEI1(IRDS,IN)*
     .                                          EIRENE_FTABEI1(IRDS,IN)
              END IF
106     CONTINUE
C
        DO 107 IION=1,NIONI
        DO 107 IPLS=1,NPLSI
        DO 107 IIDS=1,NIDSI(IION)
          IRDS=LGIEI(IION,IIDS)
          IF (PPLDS(IRDS,IPLS).EQ.0.) GOTO 107
          DO 108 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SPLNWI(IN,IION,IPLS)=SPLNWI(IN,IION,IPLS)+
     .                             TABDS1(IRDS,IN)*PPLDS(IRDS,IPLS)
            ELSE
              SPLNWI(IN,IION,IPLS)=SPLNWI(IN,IION,IPLS)+
     .                  EIRENE_FTABEI1(IRDS,IN)*PPLDS(IRDS,IPLS)
            END IF
108       CONTINUE
107     CONTINUE
C
        DO 109 IION=1,NIONI
        DO 109 IPLS=1,NPLSI
        DO 109 IIDS=1,NIDSI(IION)
          IRDS=LGIEI(IION,IIDS)
          ESIG=EPLDS(IRDS,2)
          DO 110 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SEINWI(IN,IION)=SEINWI(IN,IION)+TABDS1(IRDS,IN)*ESIG
            ELSE
              SEINWI(IN,IION)=SEINWI(IN,IION)+
     .                        EIRENE_FTABEI1(IRDS,IN)*ESIG
            END IF
110       CONTINUE
109     CONTINUE
C
C  NEW: MOLECULES, EI PROCESSES
C
        DO 115 IMOL=1,NMOLI
        DO 115 IMDS=1,NMDSI(IMOL)
          IRDS=LGMEI(IMOL,IMDS)
          DO 116 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SEENWM(IN,IMOL)=SEENWM(IN,IMOL)+EELDS1(IRDS,IN)*
     .                                        TABDS1(IRDS,IN)
            ELSE
              SEENWM(IN,IMOL)=SEENWM(IN,IMOL)+EIRENE_FEELEI1(IRDS,IN)*
     .                                        EIRENE_FTABEI1(IRDS,IN)
            END IF
116       CONTINUE
115     CONTINUE
C
        DO 117 IMOL=1,NMOLI
        DO 117 IPLS=1,NPLSI
        DO 117 IMDS=1,NMDSI(IMOL)
          IRDS=LGMEI(IMOL,IMDS)
          IF (PPLDS(IRDS,IPLS).EQ.0.) GOTO 117
          DO 118 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SPLNWM(IN,IMOL,IPLS)=SPLNWM(IN,IMOL,IPLS)+
     .                             TABDS1(IRDS,IN)*PPLDS(IRDS,IPLS)
            ELSE
              SPLNWM(IN,IMOL,IPLS)=SPLNWM(IN,IMOL,IPLS)+
     .                  EIRENE_FTABEI1(IRDS,IN)*PPLDS(IRDS,IPLS)
            END IF
118       CONTINUE
117     CONTINUE
C
        DO 119 IMOL=1,NMOLI
        DO 119 IPLS=1,NPLSI
        DO 119 IMDS=1,NMDSI(IMOL)
          IRDS=LGMEI(IMOL,IMDS)
          ESIG=EPLDS(IRDS,2)
          DO 120 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SEINWM(IN,IMOL)=SEINWM(IN,IMOL)+TABDS1(IRDS,IN)*ESIG
            ELSE
              SEINWM(IN,IMOL)=SEINWM(IN,IMOL)+
     .                        EIRENE_FTABEI1(IRDS,IN)*ESIG
            END IF
120       CONTINUE
119     CONTINUE
C
        B2BREM=B2BRM
        B2RAD=B2RD
        B2QIE=B2Q
        B2VDP=B2VP
        CALL EIRENE_INTER3(LSTP,IFIRST,1,NSTRAI,0)

        NLSRON_SAVE = NLSRON
        
        IF (ANY(NLSRON(1:NSTRAI))) THEN
!pb           IFIRST=0
           LPLASM=.TRUE.
           LSTP=LSTOP
           ITNR=ITNR+1
           GOTO 10
        END IF
C
        IFIRST=IFIRST+1

        IF (LSTOP) THEN
          CALL EIRENE_DEALLOC_COMUSR
          CALL EIRENE_DEALLOC_CESTIM
          CALL EIRENE_DEALLOC_BRASCL
          CALL EIRENE_DEALLOC_BRASPOI

          CALL GREND
        END IF

        RETURN
C
      ENDIF

      END

