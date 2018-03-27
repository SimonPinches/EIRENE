cdr  12.5.2015:  move general interface driver routine "EIRSRT" up, own routine.
cdr:             check: is eirsrt universal, then: move even further up to "main routines".
cdr  09.02.2016:  done ! syncronization of eirsrt.f started, but not completed fully
c  jan 2017: syncronisation with corresponding version in couple_SOLPS-ITER,
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
      USE EIRMOD_COMPRT,ONLY:IUNOUT
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
      CHARACTER(80) :: LINE
C   DIMENSIONIERUNG FUER GITTER
      INTEGER :: DIMXH,DIMYH,NNCUT,NNISO,
     1 NXCUT1(10),NXCUT2(10),NYCUT1(10),NYCUT2(10),
     2 NXISO1(10),NXISO2(10),NYISO1(10),NYISO2(10)
      INTEGER :: IX, IY, I, J, NP, NWISO

      REAL(DP) :: DUMMI(3)

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
          READ (30,*,ERR=100,END=100) DUMMI(1),
     .                  DUMMI(2),DUMMI(3),XPOL(IY,IX)
          READ (30,*,ERR=100,END=100) DUMMI(1),
     .                  DUMMI(2),DUMMI(3),YPOL(IY,IX)
         ENDIF
         IF (IY.EQ.dimyh) THEN
          READ (30,*,ERR=100,END=100) DUMMI(1),
     .                   DUMMI(2),XPOL(dimyh+1,IX),XPOL(IY,IX)
          READ (30,*,ERR=100,END=100) DUMMI(1),
     .                   DUMMI(2),YPOL(dimyh+1,IX),YPOL(IY,IX)
         ENDIF
12      CONTINUE
       ENDIF
       IF (IX.EQ.nxcut1(1)) THEN
        DO 14 IY = 1, DIMYH
         IF (IY.LE.dimyh-1) THEN
         READ (30,*,ERR=100,END=100) XPOL(IY,nxcut2(2)),DUMMI(1),
     .                  DUMMI(2),XPOL(IY,IX)
         READ (30,*,ERR=100,END=100) YPOL(IY,nxcut2(2)),DUMMI(1),
     .                  DUMMI(2),YPOL(IY,IX)
         ENDIF
         IF (IY.EQ.dimyh) THEN
          READ (30,*,ERR=100,END=100) 
     .                  XPOL(IY,nxcut2(2)),XPOL(dimyh+1,nxcut2(2)),
     .                                XPOL(dimyh+1,IX),XPOL(IY,IX)
          READ (30,*,ERR=100,END=100) 
     .                  YPOL(IY,nxcut2(2)),YPOL(dimyh+1,nxcut2(2)),
     .                                YPOL(dimyh+1,IX),YPOL(IY,IX)
         ENDIF
14      CONTINUE
       ENDIF
       IF ((IX.GE.nxcut2(2)).AND.(IX.LE.nxcut1(2)-1)) THEN
        DO 16 IY = 1, DIMYH
         IF (IY.LE.dimyh-1) THEN
          READ (30,*,ERR=100,END=100) DUMMI(1),
     .                   DUMMI(2),DUMMI(3),XPOL(IY,IX+1)
          READ (30,*,ERR=100,END=100) DUMMI(1),
     .                   DUMMI(2),DUMMI(3),YPOL(IY,IX+1)
         ENDIF
         IF (IY.EQ.dimyh) THEN
          READ (30,*,ERR=100,END=100) DUMMI(1),DUMMI(2),
     .                   XPOL(dimyh+1,IX+1),XPOL(IY,IX+1)
          READ (30,*,ERR=100,END=100) DUMMI(1),DUMMI(2),
     .                   YPOL(dimyh+1,IX+1),YPOL(IY,IX+1)
         ENDIF
16      CONTINUE
       ENDIF
       IF (IX.EQ.nxcut1(2)) THEN
        DO 18 IY = 1, DIMYH
         IF (IY.LE.dimyh-1) THEN
         READ (30,*,ERR=100,END=100) XPOL(IY,nxcut2(1)+1),DUMMI(1),
     .                  DUMMI(2),XPOL(IY,IX+1)
         READ (30,*,ERR=100,END=100) YPOL(IY,nxcut2(1)+1),DUMMI(1),
     .                  DUMMI(2),YPOL(IY,IX+1)
         ENDIF
         IF (IY.EQ.dimyh) THEN
          READ (30,*,ERR=100,END=100) 
     .                  XPOL(IY,nxcut2(1)+1),XPOL(dimyh+1,nxcut2(1)+1),
     .                                XPOL(dimyh+1,IX+1),XPOL(IY,IX+1)
          READ (30,*,ERR=100,END=100) 
     .                  YPOL(IY,nxcut2(1)+1),YPOL(dimyh+1,nxcut2(1)+1),
     .                                YPOL(dimyh+1,IX+1),YPOL(IY,IX+1)
         ENDIF
18      CONTINUE
       ENDIF
       IF ((IX.GE.nxcut2(1)).AND.(IX.LE.dimxh-1)) THEN
        DO 22 IY = 1, DIMYH
         IF (IY.LE.dimyh-1) THEN
          READ (30,*,ERR=100,END=100) DUMMI(1),DUMMI(2),DUMMI(3),
     .                   XPOL(IY,IX+2)
          READ (30,*,ERR=100,END=100) DUMMI(1),DUMMI(2),DUMMI(3),
     .                   YPOL(IY,IX+2)
         ENDIF
         IF (IY.EQ.dimyh) THEN
          READ (30,*,ERR=100,END=100) DUMMI(1),DUMMI(2),
     .                   XPOL(dimyh+1,IX+2),XPOL(IY,IX+2)
          READ (30,*,ERR=100,END=100) DUMMI(1),DUMMI(2),
     .                   YPOL(dimyh+1,IX+2),YPOL(IY,IX+2)
         ENDIF
22      CONTINUE
       ENDIF
       IF (IX.EQ.dimxh) THEN
        DO 24 IY = 1, DIMYH
         IF (IY.LE.dimyh-1) THEN
         READ (30,*,ERR=100,END=100) XPOL(IY,dimxh+3),DUMMI(1),DUMMI(2),
     .                  XPOL(IY,IX+2)
         READ (30,*,ERR=100,END=100) YPOL(IY,dimxh+3),DUMMI(1),DUMMI(2),
     .                  YPOL(IY,IX+2)
         ENDIF
         IF (IY.EQ.dimyh) THEN
          READ (30,*,ERR=100,END=100) 
     .                  XPOL(IY,dimxh+3),XPOL(dimyh+1,dimxh+3),
     .                                XPOL(dimyh+1,IX+2),XPOL(IY,IX+2)
          READ (30,*,ERR=100,END=100) 
     .                  YPOL(IY,dimxh+3),YPOL(dimyh+1,dimxh+3),
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
      USE EIRMOD_COMPRT, ONLY: IUNOUT
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
C  INITIALISE DUMMY
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
      WRITE (iunout,*) 'ERROR IN SUBR. INDMAP: THIS SUBR. IS VALID ONLY'
      WRITE (iunout,*) 'NCUTB>=0 BUT NCUTB = ',NCUTB
      CALL EIRENE_EXIT_OWN(1)
991   WRITE (iunout,*) 
     .  'ERROR IN SUBR. INDMAP: INCONSISTENCY IN NUMBER OF '
      WRITE (iunout,*) 'ZONES PER CUT FROM LINDA GEOMETRY DETECTED.  '
      WRITE (iunout,*) 'NCUTL = ',NCUTL, ' IENDD-IINID+1 = ',
     .                  IENDD-IINID+1
      CALL EIRENE_EXIT_OWN(1)
      END


C
C
      SUBROUTINE EIRENE_INDMPI(FIELD,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .                  NCUTB,NCUTL,NPOINT,NPPLG,NSTR,ISTR)
C
C     INDEX MAPPING: INVERSE TO SUBR. INDMAP
C
      USE EIRMOD_PRECISION
      USE EIRMOD_COMPRT, ONLY: IUNOUT
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
C  INITIALISE DUMMY
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
      WRITE (iunout,*) 'ERROR IN SUBR. INDMPI: THIS SUBR. IS VALID ONLY'
      WRITE (iunout,*) 'NCUTB>=0 BUT NCUTB = ',NCUTB
      CALL EIRENE_EXIT_OWN(1)
991   WRITE (iunout,*) 
     .  'ERROR IN SUBR. INDMPI: INCONSISTENCY IN NUMBER OF'
      WRITE (iunout,*) 'ZONES PER CUT FROM LINDA GEOMETRY DETECTED. '
      WRITE (iunout,*) 'NCUTL = ',NCUTL, ' IENDD-IINID+1 = ',
     .                  IENDD-IINID+1
      CALL EIRENE_EXIT_OWN(1)
      END


C
*//PLASM//
C=======================================================================
C          S U B R O U T I N E   P L A S M
C=======================================================================
!pb  dec. 2015 automatic detection of data format
      SUBROUTINE EIRENE_PLASM(KARD,NDIMX,NDIMY,NDIMF,N,M,NF,DUMMY)

      USE EIRMOD_PRECISION
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: KARD, NDIMX, NDIMY, NDIMF, N, M, NF
      REAL(DP), INTENT(INOUT) :: DUMMY(0:N+1,0:M+1,NF)
      INTEGER :: ND1, LIM, IF, III, IX, IY, i1, i2, i3
      character(50) :: form
      character(200) :: zeile

cdr  construct the proper format for reading from file FORT(KARD)  
cdr  character string FORM replaces old card: 910  FORMAT(5(E16.8))
c
      form = repeat(' ',50)
      read (kard,'(a200)',END=500) zeile
      i1 = index(zeile,'.')
      i2 = scan(zeile,'E,e')
      i3 = index(zeile(i2+1:),' ')
      write (form,'(A4,i0,a1,i0,a2)') '(5(E',i2+i3-1,'.',i2-i1-1,'))'
      backspace kard
c     write (6,*) 'plasm: detected format ', form

      ND1 = NDIMX + 2
      LIM = (ND1/5)*5 - 4
      DUMMY(0:N+1,0:M+1,NF)=0._DP
      DO    110  IF = 1,NDIMF
      DO    110  IY = 0,NDIMY+1
      DO    100  IX = 1,LIM,5
100     READ(KARD,FORM,END=500) (DUMMY(-1+IX-1+III,IY,IF),III = 1,5)
        IF( (LIM+4).EQ.ND1 )     GOTO 110
        READ(KARD,FORM,END=500) (DUMMY(-1+IX,IY,IF),IX = LIM+5,ND1)
110   CONTINUE
500   RETURN
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
  910 FORMAT(5(ES16.7E3))
*//END NEUTR//
      END




      SUBROUTINE EIRENE_SAVE_TALLIES (ISTRAI)
C
C  SAVE EIRENE TALLIES, SCALE PER UNIT FLUX (AMP), ON COMMON BRASCL
C  WTOTP IS NEGATIVE IN EIRENE (SINK FOR IONS).
C
C  STRATA WHICH ARE SPECIFIED BY INPUT BLOCK 14 
C     (SOURCES DEFINED FROM PLASMA CODE DATA DIRECTLY) 
C     ARE RESCALED HERE TO UNIT SOURCE STRENGTH (FLXI)

c  added in Nov. 15: ipls resolved ion energy sources eapl,empl,eipl
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
cdr  save volumetric sources for plasma species ipls: particle, momentum, ion energy
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
!PB           ALLOCATE(CPMUL)
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
!PB           ALLOCATE(CPMUL)
              CPMUL => EIRENE_NEW_MULARR()
              CPMUL%IART = IPLS
              CPMUL%ICM = IN
              CPMUL%VALUEM = PIPL(IPLS,IN)*FLXI
              CPMUL%NXTMUL => PIPLS(ISTRAI)%PMUL
              PIPLS(ISTRAI)%PMUL => CPMUL
            ENDIF
          ENDIF

          IF (LEAPL) THEN 
          IF (EAPL(IPLS,IN) .NE. 0.D0) THEN
!PB         ALLOCATE(CPMUL)
            CPMUL => EIRENE_NEW_MULARR()
            CPMUL%IART = IPLS
            CPMUL%ICM = IN
            CPMUL%VALUEM = EAPL(IPLS,IN)*FLXI
            CPMUL%NXTMUL => EAPLS(ISTRAI)%PMUL
            EAPLS(ISTRAI)%PMUL => CPMUL
          ENDIF
          ENDIF
          IF (LEMPL) THEN 
          IF (EMPL(IPLS,IN) .NE. 0.D0) THEN
!PB         ALLOCATE(CPMUL)
            CPMUL => EIRENE_NEW_MULARR()
            CPMUL%IART = IPLS
            CPMUL%ICM = IN
            CPMUL%VALUEM = EMPL(IPLS,IN)*FLXI
            CPMUL%NXTMUL => EMPLS(ISTRAI)%PMUL
            EMPLS(ISTRAI)%PMUL => CPMUL
          ENDIF
          ENDIF

          IF (LEIPL) THEN 
          IF (EIPL(IPLS,IN) .NE. 0.D0) THEN
!PB         ALLOCATE(CPMUL)
            CPMUL => EIRENE_NEW_MULARR()
            CPMUL%IART = IPLS
            CPMUL%ICM = IN
            CPMUL%VALUEM = EIPL(IPLS,IN)*FLXI
            CPMUL%NXTMUL => EIPLS(ISTRAI)%PMUL
            EIPLS(ISTRAI)%PMUL => CPMUL
          ENDIF
          ENDIF

          IF (LMAPL) THEN
            IF (MAPL(IPLS,IN) .NE. 0.D0) THEN
!pb           ALLOCATE(CPMUL)
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
!PB           ALLOCATE(CPMUL)
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
!PB           ALLOCATE(CPMUL)
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
!PB           ALLOCATE(CPMUL)
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
!PB           ALLOCATE(CPSIM)
            CPSIM => EIRENE_NEW_SIMARR()
            CPSIM%ICS = IN
            CPSIM%VALUES = EIEL(IN)*FLXI
            CPSIM%NXTSIM => EIELS(ISTRAI)%PSIM
            EIELS(ISTRAI)%PSIM => CPSIM
          ENDIF
        ENDIF

      ENDDO

      DO IATM=1,NATMI
        DO IN=1,NSBOX_TAL
          IF (LPDENA) THEN
            IF (PDENA(IATM,IN) .NE. 0.D0) THEN
!PB           ALLOCATE(CPMUL)
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

