cdr  12.5.2015:  move general interface driver routine "EIRSRT" up, own routine.
cdr:             check: is eirsrt universal, then: move even further up to "main routines".
cdr  09.02.2016:  done ! synchronization of eirsrt.f started, but not completed fully
c  jan 2017: synchronisation with corresponding version in other couple_...
c            re. reading polygon data in geomd_linda from fort.30
c            added: species index in eapl,empl,eipl tallies
cdr july 21: added: additional test re. 2nd dimension (np2nd), from forks/iter...
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
      SUBROUTINE EIRENE_GEOMD(NDXA,NDYA,NPLP,NR1ST,NP2ND,
     .                 PUX,PUY,PVX,PVY,itype)
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CCONA
      USE EIRMOD_CGEOM
      IMPLICIT NONE
C
      REAL(DP), INTENT(OUT) :: PUX(*),PUY(*),PVX(*),PVY(*)
      INTEGER, INTENT(INOUT) :: NDXA,NDYA,NPLP,NR1ST,NP2ND
      INTEGER, INTENT(IN) :: ITYPE

      TARGINDEX = 3
      IF (ITYPE == 1) THEN
        CALL EIRENE_GEOMD_SONNET(NDXA,NDYA,NPLP,NR1ST,
     .                 PUX,PUY,PVX,PVY)
C
      ELSEIF (ITYPE == 2) THEN
        CALL EIRENE_GEOMD_CARRE(NDXA,NDYA,NPLP,NR1ST,
     .                 PUX,PUY,PVX,PVY)

      ELSE
        CALL EIRENE_GEOMD_LINDA(NDXA,NDYA,NPLP,NR1ST,NP2ND,
     .                 PUX,PUY,PVX,PVY)
      END IF

      RETURN
*//END GEOMD//
      END SUBROUTINE EIRENE_GEOMD

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

    1 continue
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


   99 continue
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
!pb 1014   CONTINUE
!pb        CALL MSHPROJ (X1,Y1,X2,Y2,X3,Y3,X4,Y4,PUX,PUY,PVX,PVY,NDXA,
!pb     .                NR1ST,IY)
!pb 1015 CONTINUE
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
 1014   CONTINUE
C
        CALL EIRENE_MSHPROJ (X1,Y1,X2,Y2,X3,Y3,X4,Y4,PUX,PUY,PVX,PVY,
     .                NDXA,NR1ST,IY)
 1015 CONTINUE
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
      END SUBROUTINE EIRENE_GEOMD_CARRE

*//GEOMD//
C=======================================================================
C          S U B R O U T I N E   G E O M D _ L I N D A
C=======================================================================
      SUBROUTINE EIRENE_GEOMD_LINDA(NDXA,NDYA,NPLP,NR1ST,NP2ND,
     .                        PUX,PUY,PVX,PVY)
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CCONA
      USE EIRMOD_CPOLYG
      USE EIRMOD_CGEOM
      USE EIRMOD_CINIT
      USE EIRMOD_COMPRT,ONLY:IUNOUT
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: NDXA, NDYA, NR1ST, NP2ND
      INTEGER, INTENT(OUT) :: NPLP
      REAL(DP), INTENT(OUT) :: PUX(*),PUY(*),PVX(*),PVY(*)
C
C  GEOMETRY DATA: CELL VERTICES (LINDA ---> EIRENE)
      REAL(DP) ::
     R  X1(NDX),Y1(NDX),X2(NDX),Y2(NDX),X3(NDX),Y3(NDX),
     R  X4(NDX),Y4(NDX)
C
      CHARACTER(80) :: LINE
C   DIMENSIONIERUNG FUER GITTER
      INTEGER :: DIMXH,DIMYH,NNCUT,NNISO
      INTEGER, ALLOCATABLE ::
     1 NXCUT1(:),NXCUT2(:),NYCUT1(:),NYCUT2(:),
     2 NXISO1(:),NXISO2(:),NYISO1(:),NYISO2(:)
      INTEGER :: IX, IY, I, J, NP, NWISO

      character*54 format_string(2)
      character*7 filename
      integer, save :: new_format, exp_location

      data new_format/1/
      data format_string/
     . '(T2,1E15.7,T17,1E15.7,T32,1E15.7,T47,1E15.7)',
     . '(T2,1E16.8,T18,1E16.8,T34,1E16.8,T50,1E16.8)' /

      REAL(DP) :: DUMMI(4)

C  ACTUAL MESH USED IN THIS RUN
C
C      EINLESEROUTINE ANGEPASST AUF BRAAMS-OUTPUT
C   GEAENDERTE DIMENSIONIERUNG BZW. CUT-POSITION
C        MUSS PER HAND ANGEPASST WERDEN:
C        PARAMETER DIMXH,DIMYH                        RFS 14.5.1991
      filename = fort_lc//'30'
      OPEN (UNIT=30,ACCESS='SEQUENTIAL',FORM='FORMATTED',ERR=100) !VK
      REWIND 30
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
      allocate(nxcut1(0:nncut),nxcut2(0:nncut),
     .         nycut1(0:nncut),nycut2(0:nncut))
      read(30,*,ERR=100,END=100)
     r     (nxcut1(i),nxcut2(i),nycut1(i),nycut2(i),i=1,nncut)
      if (nncut.gt.2) then
         read(30,*,ERR=100,END=100) nniso
         allocate(nxiso1(0:nniso),nxiso2(0:nniso),
     .            nyiso1(0:nniso),nyiso2(0:nniso))
         read(30,*,ERR=100,END=100)
     r            (nxiso1(i),nxiso2(i),nyiso1(i),nyiso2(i),i=1,nniso)
      ELSE
       NNISO=0 !VK
      endif
      read(30,*,ERR=100,END=100)
C    ANZAHL DER TEILSTUECKE PRO POLYGON
      NPLP = MAX((NNCUT/2)*3,1)
C    COMPUTE WIDTH OF INSULATING CUT FOR DOUBLE NULL
      IF (NNISO.GE.1) THEN
        NWISO=NXISO2(1)-NXISO1(1)
      ELSE
        NWISO=0
      ENDIF
C     PRINT MESSAGE AND CHECK
      CALL EIRENE_LEER (1)
      CALL EIRENE_MASJ4 ('GEOMD: NNCUT, NNISO, NPLP, NWISO',
     w                           NNCUT, NNISO, NPLP, NWISO)
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
      read (30,'(A80)') line
      exp_location = 81
      if (index(line,'E').ne.0)
     . exp_location = min(exp_location,index(line,'E'))
      if (index(line,'e').ne.0)
     . exp_location = min(exp_location,index(line,'e'))
      if (index(line,'d').ne.0)
     . exp_location = min(exp_location,index(line,'d'))
      if (index(line,'D').ne.0)
     . exp_location = min(exp_location,index(line,'D'))
      if (exp_location.eq.12) then
        new_format = 1
      else if (exp_location.eq.13) then
        new_format = 2
      else
        call eirene_masage
     .   ('Unrecognized format in '//trim(filename)//' file')
        call eirene_exit_own(1)
      endif
      backspace(30)
      write(iunout,'(a,a)') 'Detected '//trim(filename)//' is using ',
     .  format_string(new_format)

      DO 10 IX = 1, DIMXH
        IF (NNCUT.GE.1) THEN
          IF (IX.LE.nxcut1(1)-1) THEN
            DO 12 IY = 1, DIMYH
              IF (IY.LE.dimyh-1) THEN
                READ (30,format_string(new_format)) DUMMI(1),
     .                  DUMMI(2),DUMMI(3),XPOL(IY,IX)
                READ (30,format_string(new_format)) DUMMI(1),
     .                  DUMMI(2),DUMMI(3),YPOL(IY,IX)
              ENDIF
              IF (IY.EQ.dimyh) THEN
                READ (30,format_string(new_format)) DUMMI(1),
     .                   DUMMI(2),XPOL(dimyh+1,IX),XPOL(IY,IX)
                READ (30,format_string(new_format)) DUMMI(1),
     .                   DUMMI(2),YPOL(dimyh+1,IX),YPOL(IY,IX)
              ENDIF
   12       CONTINUE
          ENDIF
          IF (IX.EQ.nxcut1(1)) THEN
            DO 14 IY = 1, DIMYH
              IF (IY.LE.dimyh-1) THEN
                READ (30,format_string(new_format))
     .                         XPOL(IY,nxcut2(NNCUT)),DUMMI(1),
     .                         DUMMI(2),XPOL(IY,IX)
                READ (30,format_string(new_format))
     .                         YPOL(IY,nxcut2(NNCUT)),DUMMI(1),
     .                         DUMMI(2),YPOL(IY,IX)
              ENDIF
              IF (IY.EQ.dimyh) THEN
                READ (30,format_string(new_format))
     .                         XPOL(IY,nxcut2(NNCUT)),
     .                         XPOL(dimyh+1,nxcut2(NNCUT)),
     .                         XPOL(dimyh+1,IX),XPOL(IY,IX)
                READ (30,format_string(new_format))
     .                         YPOL(IY,nxcut2(NNCUT)),
     .                         YPOL(dimyh+1,nxcut2(NNCUT)),
     .                         YPOL(dimyh+1,IX),YPOL(IY,IX)
              ENDIF
 14         CONTINUE
          ENDIF
        ENDIF
C    READING OF SECOND POLOIDAL PATCH
        IF (NNCUT.GE.2) THEN
C    Snowflake Case
         IF (nxcut2(1).LE.nxcut1(2)) THEN
          IF ((IX.GE.nxcut2(NNCUT)).AND.(IX.LE.nxcut1(NNCUT)-1)) THEN
            DO 19 IY = 1, DIMYH
              IF (IY.LE.dimyh-1) THEN
                READ (30,format_string(new_format)) DUMMI(1),
     .                       DUMMI(2),DUMMI(3),XPOL(IY,IX+1)
                READ (30,format_string(new_format)) DUMMI(1),
     .                       DUMMI(2),DUMMI(3),YPOL(IY,IX+1)
              ENDIF
              IF (IY.EQ.dimyh) THEN
                READ (30,format_string(new_format))
     .                       DUMMI(1),DUMMI(2),
     .                       XPOL(dimyh+1,IX+1),XPOL(IY,IX+1)
                READ (30,format_string(new_format))
     .                       DUMMI(1),DUMMI(2),
     .                       YPOL(dimyh+1,IX+1),YPOL(IY,IX+1)
              ENDIF
   19       CONTINUE
          ENDIF
          IF (IX.EQ.nxcut1(NNCUT)) THEN
            DO 20 IY = 1, DIMYH
              IF (IY.LE.dimyh-1) THEN
                READ (30,format_string(new_format))
     .                       XPOL(IY,IX+2),DUMMI(1),
     .                       DUMMI(2),XPOL(IY,IX+1)
                READ (30,format_string(new_format))
     .                       YPOL(IY,IX+2),DUMMI(1),
     .                       DUMMI(2),YPOL(IY,IX+1)
              ENDIF
              IF (IY.EQ.dimyh) THEN
                READ (30,format_string(new_format)) XPOL(IY,IX+2),
     .                       XPOL(dimyh+1,IX+2),
     .                       XPOL(dimyh+1,IX+1),XPOL(IY,IX+1)
                READ (30,format_string(new_format)) YPOL(IY,IX+2),
     .                       YPOL(dimyh+1,IX+2),
     .                       YPOL(dimyh+1,IX+1),YPOL(IY,IX+1)
              ENDIF
   20       CONTINUE
          ENDIF
C   DDN/SN Case
         ELSE
          IF ((IX.GE.nxcut2(NNCUT)).AND.(IX.LE.nxcut1(2)-1)) THEN
            DO 190 IY = 1, DIMYH
              IF (IY.LE.dimyh-1) THEN
                READ (30,format_string(new_format)) DUMMI(1),
     .                         DUMMI(2),DUMMI(3),XPOL(IY,IX+1)
                READ (30,format_string(new_format)) DUMMI(1),
     .                         DUMMI(2),DUMMI(3),YPOL(IY,IX+1)
              ENDIF
              IF (IY.EQ.dimyh) THEN
                READ (30,format_string(new_format))
     .                         DUMMI(1),DUMMI(2),
     .                         XPOL(dimyh+1,IX+1),XPOL(IY,IX+1)
                READ (30,format_string(new_format))
     .                         DUMMI(1),DUMMI(2),
     .                         YPOL(dimyh+1,IX+1),YPOL(IY,IX+1)
              ENDIF
 190        CONTINUE
          ENDIF
          IF (IX.EQ.nxcut1(2)) THEN
            DO 200 IY = 1, DIMYH
              IF (IY.LE.dimyh-1) THEN
                READ (30,format_string(new_format))
     .                         XPOL(IY,IX+2),DUMMI(1),
     .                         DUMMI(2),XPOL(IY,IX+1)
                READ (30,format_string(new_format))
     .                         YPOL(IY,IX+2),DUMMI(1),
     .                         DUMMI(2),YPOL(IY,IX+1)
              ENDIF
              IF (IY.EQ.dimyh) THEN
                READ (30,format_string(new_format)) XPOL(IY,IX+2),
     .                         XPOL(dimyh+1,IX+2),
     .                         XPOL(dimyh+1,IX+1),XPOL(IY,IX+1)
                READ (30,format_string(new_format)) YPOL(IY,IX+2),
     .                         YPOL(dimyh+1,IX+2),
     .                         YPOL(dimyh+1,IX+1),YPOL(IY,IX+1)
              ENDIF
 200        CONTINUE
          ENDIF
         ENDIF
        ENDIF
C    READING OF THIRD AND FOURTH POLOIDAL PATCHES
        IF (NNCUT.GE.4 .AND. NNISO.GE.1) THEN
C    Snowflake Case
         IF (nxcut2(1).LE.nxcut1(2)) THEN
          IF ((IX.GE.nxcut1(4)+1).AND.(IX.LE.nxcut1(2)-1)) THEN
            DO 26 IY = 1, DIMYH
              IF (IY.LE.dimyh-1) THEN
                READ (30,format_string(new_format)) DUMMI(1),
     .                         DUMMI(2),DUMMI(3),XPOL(IY,IX+2)
                READ (30,format_string(new_format)) DUMMI(1),
     .                         DUMMI(2),DUMMI(3),YPOL(IY,IX+2)
              ENDIF
              IF (IY.EQ.dimyh) THEN
                READ (30,format_string(new_format)) DUMMI(1),DUMMI(2),
     .                         XPOL(dimyh+1,IX+2),XPOL(IY,IX+2)
                READ (30,format_string(new_format)) DUMMI(1),DUMMI(2),
     .                         YPOL(dimyh+1,IX+2),YPOL(IY,IX+2)
              ENDIF
   26       CONTINUE
          ENDIF
          IF (IX.EQ.nxcut1(2)) THEN
            DO 36 IY = 1, DIMYH
              IF (IY.LE.dimyh-1) THEN
                READ (30,format_string(new_format))
     .                         XPOL(IY,IX+3),DUMMI(1),
     .                         DUMMI(2),XPOL(IY,IX+2)
                READ (30,format_string(new_format))
     .                         YPOL(IY,IX+3),DUMMI(1),
     .                         DUMMI(2),YPOL(IY,IX+2)
              ENDIF
              IF (IY.EQ.dimyh) THEN
                READ (30,format_string(new_format))
     .                         XPOL(IY,IX+3),XPOL(IY+1,IX+3),
     .                         XPOL(IY+1,IX+2),XPOL(IY,IX+2)
                READ (30,format_string(new_format))
     .                         YPOL(IY,IX+3),YPOL(IY+1,IX+3),
     .                         YPOL(IY+1,IX+2),YPOL(IY,IX+2)
              ENDIF
   36       CONTINUE
          ENDIF
          IF ((IX.GE.nxcut1(2)+1).AND.(IX.LE.nxiso1(1)-1)) THEN
            DO 38 IY = 1, DIMYH
              IF (IY.LE.dimyh-1) THEN
                READ (30,format_string(new_format)) DUMMI(1),
     .                         DUMMI(2),DUMMI(3),XPOL(IY,IX+3-NWISO)
                READ (30,format_string(new_format)) DUMMI(1),
     .                         DUMMI(2),DUMMI(3),YPOL(IY,IX+3-NWISO)
              ENDIF
              IF (IY.EQ.dimyh) THEN
                READ (30,format_string(new_format)) DUMMI(1),DUMMI(2),
     .                         XPOL(dimyh+1,IX+3-NWISO),
     .                         XPOL(IY,IX+3-NWISO)
                READ (30,format_string(new_format)) DUMMI(1),DUMMI(2),
     .                         YPOL(dimyh+1,IX+3-NWISO),
     .                         YPOL(IY,IX+3-NWISO)
              ENDIF
   38       CONTINUE
          ENDIF
          IF (IX.EQ.nxiso1(1)) THEN
            DO 28 IY = 1, DIMYH
              IF (IY.LE.dimyh-1) THEN
                READ (30,format_string(new_format))
     .                         XPOL(IY,IX+4-NWISO),DUMMI(1),
     .                         DUMMI(2),XPOL(IY,IX+3-NWISO)
                READ (30,format_string(new_format))
     .                         YPOL(IY,IX+4-NWISO),DUMMI(1),
     .                         DUMMI(2),YPOL(IY,IX+3-NWISO)
              ENDIF
              IF (IY.EQ.dimyh) THEN
                READ (30,format_string(new_format))
     .                         XPOL(IY,IX+4-NWISO),
     .                         XPOL(dimyh+1,IX+4-NWISO),
     .                         XPOL(dimyh+1,IX+3-NWISO),
     .                         XPOL(IY,IX+3-NWISO)
                READ (30,format_string(new_format))
     .                         YPOL(IY,IX+4-NWISO),
     .                         YPOL(dimyh+1,IX+4-NWISO),
     .                         YPOL(dimyh+1,IX+3-NWISO),
     .                         YPOL(IY,IX+3-NWISO)
              ENDIF
   28       CONTINUE
          ENDIF
          IF ((IX.GE.nxiso2(1)+1).AND.(IX.LE.nxcut1(3)-1)) THEN
            DO 32 IY = 1, DIMYH
              IF (IY.LE.dimyh-1) THEN
                READ (30,format_string(new_format)) DUMMI(1),DUMMI(2),
     .                         DUMMI(3),
     .                         XPOL(IY,IX+4-NWISO)
                READ (30,format_string(new_format)) DUMMI(1),DUMMI(2),
     .                         DUMMI(3),
     .                         YPOL(IY,IX+4-NWISO)
              ENDIF
              IF (IY.EQ.dimyh) THEN
                READ (30,format_string(new_format)) DUMMI(1),DUMMI(2),
     .                         XPOL(dimyh+1,IX+4-NWISO),
     .                         XPOL(IY,IX+4-NWISO)
                READ (30,format_string(new_format)) DUMMI(1),DUMMI(2),
     .                         YPOL(dimyh+1,IX+4-NWISO),
     .                         YPOL(IY,IX+4-NWISO)
              ENDIF
   32       CONTINUE
          ENDIF
          IF (IX.EQ.nxcut1(3)) THEN
            DO 34 IY = 1, DIMYH
              IF (IY.LE.dimyh-1) THEN
                READ (30,format_string(new_format))
     .                         XPOL(IY,IX+5-NWISO),DUMMI(1),
     .                         DUMMI(2),XPOL(IY,IX+4-NWISO)
                READ (30,format_string(new_format))
     .                         YPOL(IY,IX+5-NWISO),DUMMI(1),
     .                         DUMMI(2),YPOL(IY,IX+4-NWISO)
              ENDIF
              IF (IY.EQ.dimyh) THEN
                READ (30,format_string(new_format))
     .                         XPOL(IY,IX+5-NWISO),
     .                         XPOL(dimyh+1,IX+5-NWISO),
     .                         XPOL(dimyh+1,IX+4-NWISO),
     .                         XPOL(IY,IX+4-NWISO)
                READ (30,format_string(new_format))
     .                         YPOL(IY,IX+5-NWISO),
     .                         YPOL(dimyh+1,IX+5-NWISO),
     .                         YPOL(dimyh+1,IX+4-NWISO),
     .                         YPOL(IY,IX+4-NWISO)
              ENDIF
   34       CONTINUE
          ENDIF
          ELSE
C    DDN case
           IF ((IX.GE.nxcut2(3)).AND.(IX.LE.nxiso1(1)-1)) THEN
             DO 260 IY = 1, DIMYH
               IF (IY.LE.dimyh-1) THEN
                 READ (30,format_string(new_format)) DUMMI(1),
     .                          DUMMI(2),DUMMI(3),XPOL(IY,IX+2)
                 READ (30,format_string(new_format)) DUMMI(1),
     .                          DUMMI(2),DUMMI(3),YPOL(IY,IX+2)
               ENDIF
               IF (IY.EQ.dimyh) THEN
                 READ (30,format_string(new_format)) DUMMI(1),DUMMI(2),
     .                          XPOL(dimyh+1,IX+2),XPOL(IY,IX+2)
                 READ (30,format_string(new_format)) DUMMI(1),DUMMI(2),
     .                          YPOL(dimyh+1,IX+2),YPOL(IY,IX+2)
               ENDIF
  260        CONTINUE
           ENDIF
           IF (IX.EQ.nxiso1(1)) THEN
             DO 360 IY = 1, DIMYH
               IF (IY.LE.dimyh-1) THEN
                 READ (30,format_string(new_format))
     .                          XPOL(IY,IX+3),DUMMI(1),
     .                          DUMMI(2),XPOL(IY,IX+2)
                 READ (30,format_string(new_format))
     .                          YPOL(IY,IX+3),DUMMI(1),
     .                          DUMMI(2),YPOL(IY,IX+2)
               ENDIF
               IF (IY.EQ.dimyh) THEN
                 READ (30,format_string(new_format))
     .                          XPOL(IY,IX+3),XPOL(IY+1,IX+3),
     .                          XPOL(IY+1,IX+2),XPOL(IY,IX+2)
                 READ (30,format_string(new_format))
     .                          YPOL(IY,IX+3),YPOL(IY+1,IX+3),
     .                          YPOL(IY+1,IX+2),YPOL(IY,IX+2)
               ENDIF
  360        CONTINUE
           ENDIF
           IF ((IX.GE.nxiso2(1)+1).AND.(IX.LE.nxcut1(3)-1)) THEN
             DO 380 IY = 1, DIMYH
               IF (IY.LE.dimyh-1) THEN
                 READ (30,format_string(new_format)) DUMMI(1),
     .                          DUMMI(2),DUMMI(3),XPOL(IY,IX+3-NWISO)
                 READ (30,format_string(new_format)) DUMMI(1),
     .                          DUMMI(2),DUMMI(3),YPOL(IY,IX+3-NWISO)
               ENDIF
               IF (IY.EQ.dimyh) THEN
                 READ (30,format_string(new_format)) DUMMI(1),DUMMI(2),
     .                          XPOL(dimyh+1,IX+3-NWISO),
     .                          XPOL(IY,IX+3-NWISO)
                 READ (30,format_string(new_format)) DUMMI(1),DUMMI(2),
     .                          YPOL(dimyh+1,IX+3-NWISO),
     .                          YPOL(IY,IX+3-NWISO)
               ENDIF
  380        CONTINUE
           ENDIF
           IF (IX.EQ.nxcut1(3)) THEN
             DO 280 IY = 1, DIMYH
               IF (IY.LE.dimyh-1) THEN
                 READ (30,format_string(new_format))
     .                          XPOL(IY,IX+4-NWISO),DUMMI(1),
     .                          DUMMI(2),XPOL(IY,IX+3-NWISO)
                 READ (30,format_string(new_format))
     .                          YPOL(IY,IX+4-NWISO),DUMMI(1),
     .                          DUMMI(2),YPOL(IY,IX+3-NWISO)
               ENDIF
               IF (IY.EQ.dimyh) THEN
                 READ (30,format_string(new_format))
     .                          XPOL(IY,IX+4-NWISO),
     .                          XPOL(dimyh+1,IX+4-NWISO),
     .                          XPOL(dimyh+1,IX+3-NWISO),
     .                          XPOL(IY,IX+3-NWISO)
                 READ (30,format_string(new_format))
     .                          YPOL(IY,IX+4-NWISO),
     .                          YPOL(dimyh+1,IX+4-NWISO),
     .                          YPOL(dimyh+1,IX+3-NWISO),
     .                          YPOL(IY,IX+3-NWISO)
               ENDIF
  280        CONTINUE
           ENDIF
           IF ((IX.GE.nxcut2(2)).AND.(IX.LE.nxcut1(4)-1)) THEN
             DO 320 IY = 1, DIMYH
               IF (IY.LE.dimyh-1) THEN
                 READ (30,format_string(new_format)) DUMMI(1),DUMMI(2),
     .                          DUMMI(3),
     .                          XPOL(IY,IX+4-NWISO)
                 READ (30,format_string(new_format)) DUMMI(1),DUMMI(2),
     .                          DUMMI(3),
     .                          YPOL(IY,IX+4-NWISO)
               ENDIF
               IF (IY.EQ.dimyh) THEN
                 READ (30,format_string(new_format)) DUMMI(1),DUMMI(2),
     .                          XPOL(dimyh+1,IX+4-NWISO),
     .                          XPOL(IY,IX+4-NWISO)
                 READ (30,format_string(new_format)) DUMMI(1),DUMMI(2),
     .                          YPOL(dimyh+1,IX+4-NWISO),
     .                          YPOL(IY,IX+4-NWISO)
               ENDIF
  320        CONTINUE
           ENDIF
           IF (IX.EQ.nxcut1(4)) THEN
             DO 340 IY = 1, DIMYH
               IF (IY.LE.dimyh-1) THEN
                 READ (30,format_string(new_format))
     .                          XPOL(IY,IX+5-NWISO),DUMMI(1),
     .                          DUMMI(2),XPOL(IY,IX+4-NWISO)
                 READ (30,format_string(new_format))
     .                          YPOL(IY,IX+5-NWISO),DUMMI(1),
     .                          DUMMI(2),YPOL(IY,IX+4-NWISO)
               ENDIF
               IF (IY.EQ.dimyh) THEN
                 READ (30,format_string(new_format))
     .                          XPOL(IY,IX+5-NWISO),
     .                          XPOL(dimyh+1,IX+5-NWISO),
     .                          XPOL(dimyh+1,IX+4-NWISO),
     .                          XPOL(IY,IX+4-NWISO)
                 READ (30,format_string(new_format))
     .                          YPOL(IY,IX+5-NWISO),
     .                          YPOL(dimyh+1,IX+5-NWISO),
     .                          YPOL(dimyh+1,IX+4-NWISO),
     .                          YPOL(IY,IX+4-NWISO)
               ENDIF
  340        CONTINUE
           ENDIF
          ENDIF
        ENDIF
        IF (NNCUT.GE.1) THEN
C    snowflake/XPT case
         IF (nxcut2(1).LE.nxcut1(2)) THEN
          IF ((IX.GE.nxcut1(3)+1).AND.(IX.LE.dimxh-1)) THEN
            DO 22 IY = 1, DIMYH
              IF (IY.LE.dimyh-1) THEN
                READ (30,format_string(new_format)) DUMMI(1),DUMMI(2),
     .                         DUMMI(3),XPOL(IY,IX+NPLP-1-NWISO)
                READ (30,format_string(new_format)) DUMMI(1),DUMMI(2),
     .                         DUMMI(3),YPOL(IY,IX+NPLP-1-NWISO)
              ENDIF
              IF (IY.EQ.dimyh) THEN
                READ (30,format_string(new_format)) DUMMI(1),DUMMI(2),
     .                         XPOL(dimyh+1,IX+NPLP-1-NWISO),
     .                         XPOL(IY,IX+NPLP-1-NWISO)
                READ (30,format_string(new_format)) DUMMI(1),DUMMI(2),
     .                         YPOL(dimyh+1,IX+NPLP-1-NWISO),
     .                         YPOL(IY,IX+NPLP-1-NWISO)
              ENDIF
 22         CONTINUE
          ENDIF
         ELSE
C   DDN/SN cases
          IF ((IX.GE.nxcut2(1)).AND.(IX.LE.dimxh-1)) THEN
            DO 220 IY = 1, DIMYH
              IF (IY.LE.dimyh-1) THEN
                READ (30,format_string(new_format)) DUMMI(1),DUMMI(2),
     .                        DUMMI(3),XPOL(IY,IX+NPLP-1-NWISO)
                READ (30,format_string(new_format)) DUMMI(1),DUMMI(2),
     .                        DUMMI(3),YPOL(IY,IX+NPLP-1-NWISO)
              ENDIF
              IF (IY.EQ.dimyh) THEN
                READ (30,format_string(new_format)) DUMMI(1),DUMMI(2),
     .                         XPOL(dimyh+1,IX+NPLP-1-NWISO),
     .                         XPOL(IY,IX+NPLP-1-NWISO)
                READ (30,format_string(new_format)) DUMMI(1),DUMMI(2),
     .                         YPOL(dimyh+1,IX+NPLP-1-NWISO),
     .                         YPOL(IY,IX+NPLP-1-NWISO)
              ENDIF
  220       CONTINUE
          ENDIF
         ENDIF
        ELSEIF (NNCUT.EQ.0) THEN
          IF (IX.LE.dimxh-1) THEN
            DO IY = 1, DIMYH
              IF (IY.LE.dimyh-1) THEN
                READ (30,format_string(new_format)) DUMMI(1),DUMMI(2),
     .                         DUMMI(3),XPOL(IY,IX+NPLP-1-NWISO)
                READ (30,format_string(new_format)) DUMMI(1),DUMMI(2),
     .                         DUMMI(3),YPOL(IY,IX+NPLP-1-NWISO)
              ENDIF
              IF (IY.EQ.dimyh) THEN
                READ (30,format_string(new_format)) DUMMI(1),DUMMI(2),
     .                         XPOL(dimyh+1,IX+NPLP-1-NWISO),
     .                         XPOL(IY,IX+NPLP-1-NWISO)
                READ (30,format_string(new_format)) DUMMI(1),DUMMI(2),
     .                         YPOL(dimyh+1,IX+NPLP-1-NWISO),
     .                         YPOL(IY,IX+NPLP-1-NWISO)
              ENDIF
            ENDDO
          ENDIF
        ENDIF
        IF (IX.EQ.dimxh) THEN
          DO 24 IY = 1, DIMYH
            IF (IY.LE.dimyh-1) THEN
              READ (30,format_string(new_format))
     .                         XPOL(IY,dimxh+NPLP-NWISO),DUMMI(1),
     .                         DUMMI(2),XPOL(IY,IX+NPLP-1-NWISO)
              READ (30,format_string(new_format))
     .                         YPOL(IY,dimxh+NPLP-NWISO),DUMMI(1),
     .                         DUMMI(2),YPOL(IY,IX+NPLP-1-NWISO)
            ENDIF
            IF (IY.EQ.dimyh) THEN
              READ (30,format_string(new_format))
     .                         XPOL(IY,dimxh+NPLP-NWISO),
     .                         XPOL(dimyh+1,dimxh+NPLP-NWISO),
     .                         XPOL(dimyh+1,IX+NPLP-1-NWISO),
     .                         XPOL(IY,IX+NPLP-1-NWISO)
              READ (30,format_string(new_format))
     .                         YPOL(IY,dimxh+NPLP-NWISO),
     .                         YPOL(dimyh+1,dimxh+NPLP-NWISO),
     .                         YPOL(dimyh+1,IX+NPLP-1-NWISO),
     .                         YPOL(IY,IX+NPLP-1-NWISO)
            ENDIF
   24     CONTINUE
        ENDIF
   10 CONTINUE

C   ANFANGSPUNKT DES ERSTEN TEILSTUECKS DES I-TEN POLYGONS
C   START OF FIRST POLOIDAL SECTION OF THE GRID
      NPOINT(1,1)=1
C   ENDPUNKT DES ERSTEN TEILSTUECKS DES I-TEN POLYGONS
C   Fill array of cuts in the case of a Low field Snowflake/XPT
      IF (nxcut2(1).LE.nxcut1(2)) THEN
       IF (NNCUT.LE.1) THEN
         NPOINT(2,1)=dimxh+1
       ELSE
        NPOINT(2,1)=nxcut1(1)+1
        NPOINT(1,2)=nxcut1(1)+2
        NPOINT(2,2)=nxcut1(4)+2
        NPOINT(1,3)=nxcut1(4)+3
        IF (NNCUT.EQ.2) NPOINT(2,3)=dimxh+3
        IF (NNCUT.EQ.4) THEN
          NPOINT(2,3)=nxcut1(2)+3
          NPOINT(1,4)=nxcut1(2)+4
          NPOINT(2,4)=nxiso1(1)+4
          NPOINT(1,5)=nxiso2(1)+5-NWISO
          NPOINT(2,5)=nxcut1(3)+5-NWISO
          NPOINT(1,6)=nxcut1(3)+6-NWISO
          NPOINT(2,6)=dimxh+6-NWISO
        ENDIF
       ENDIF
       TARGINDEX = 4
C   Fill array of cuts in the case of a ddn
      ELSE
       IF (NNCUT.LE.1) THEN
C   ENDPUNKT DES ERSTEN TEILSTUECKS DES I-TEN POLYGONS
C   END OF FIRST POLOIDAL SECTION OF THE GRID
        NPOINT(2,1)=dimxh+1
       ELSE
        NPOINT(2,1)=nxcut1(1)+1
C   ANFANGSPUNKT DES ZWEITEN TEILSTUECKS DES I-TEN POLYGONS
C   START OF SECOND POLOIDAL SECTION OF THE GRID
        NPOINT(1,2)=nxcut2(nncut)+1
C   ENDPUNKT DES ZWEITEN TEILSTUECKS DES I-TEN POLYGONS
C   END OF SECOND POLOIDAL SECTION OF THE GRID
        NPOINT(2,2)=nxcut2(nncut-1)+1
C   ANFANGSPUNKT DES DRITTEN TEILSTUECKS DES I-TEN POLYGONS
C   START OF THIRD POLOIDAL SECTION OF THE GRID
        NPOINT(1,3)=nxcut2(nncut-1)+2
C   ENDPUNKT DES DRITTEN TEILSTUECKS DES I-TEN POLYGONS
        IF (NNCUT.EQ.2) NPOINT(2,3)=dimxh+3
        IF (NNCUT.EQ.4) THEN
C   ENDPUNKT DES DRITTEN TEILSTUECKS DES I-TEN POLYGONS
C   END OF THIRD POLOIDAL SECTION OF THE GRID
          NPOINT(2,3)=NXISO1(1)+3
C   ANFANGSPUNKT DES VIERTEN TEILSTUECKS DES I-TEN POLYGONS
C   START OF FOURTH POLOIDAL SECTION OF THE GRID
          NPOINT(1,4)=nxiso2(1)+4-NWISO
C   ENDPUNKT DES VIERTEN TEILSTUECKS DES I-TEN POLYGONS
C   END OF FOURTH POLOIDAL SECTION OF THE GRID
          NPOINT(2,4)=nxcut1(3)+4-NWISO
C   ANFANGSPUNKT DES FUNFTEN TEILSTUECKS DES I-TEN POLYGONS
C   START OF FIFTH POLOIDAL SECTION OF THE GRID
          NPOINT(1,5)=nxcut2(2)+4-NWISO
C   ENDPUNKT DES FUNFTEN TEILSTUECKS DES I-TEN POLYGONS
C   END OF FIFTH POLOIDAL SECTION OF THE GRID
          NPOINT(2,5)=nxcut1(4)+5-NWISO
C   ANFANGSPUNKT DES SECHSTEN TEILSTUECKS DES I-TEN POLYGONS
C   START OF SIXTH POLOIDAL SECTION OF THE GRID
          NPOINT(1,6)=nxcut2(1)+5-NWISO
C   ENDPUNKT DES SECHSTEN TEILSTUECKS DES I-TEN POLYGONS
C   END OF SIXTH POLOIDAL SECTION OF THE GRID
          NPOINT(2,6)=dimxh+6-NWISO
        ENDIF
      ENDIF
       TARGINDEX = 3
      ENDIF
C
C   CHECK THAT ARRAY SIZES DEDUCED FROM fort.30 MATCH THOSE FROM INPUT FILE
C
      IF (NPPLG.NE.NPLP) THEN
        WRITE(IUNOUT,*)
     .   'MISMATCH BETWEEN INPUT FILE AND '//trim(filename)
        WRITE(IUNOUT,*) 'FOR THE NUMBER OF SEGMENTS : '
        WRITE(IUNOUT,*) 'INPUT    : NPPLG = ', NPPLG
        WRITE(IUNOUT,*) filename//'  : NNCUT = ', NNCUT
        WRITE(IUNOUT,*) 'WE MUST HAVE NPPLG = MAX((NNCUT/2)*3,1)'
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
      IF (NR1ST.NE.dimyh+1) THEN
        WRITE(IUNOUT,*)
     .   'MISMATCH BETWEEN INPUT FILE AND '//trim(filename)
        WRITE(IUNOUT,*) 'FOR THE NUMBER OF RADIAL CELLS : '
        WRITE(IUNOUT,*) 'INPUT    : NR1ST = ', NR1ST
        WRITE(IUNOUT,*) filename//'  : NY    = ', dimyh
        WRITE(IUNOUT,*) 'WE MUST HAVE NR1ST = NY+1'
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
      IF (NP2ND.NE.NPOINT(2,NPLP)) THEN
        WRITE(IUNOUT,*)
     .   'MISMATCH BETWEEN INPUT FILE AND '//trim(filename)
        WRITE(IUNOUT,*) 'FOR THE NUMBER OF POLOIDAL CELLS : '
        WRITE(IUNOUT,*) 'INPUT    : NP2ND = ', NP2ND
        WRITE(IUNOUT,*) filename//'  : NX    = ', dimxh
        WRITE(IUNOUT,*) filename//'  : NNCUT = ', NNCUT
        WRITE(IUNOUT,*) filename//'  : NNISO = ', NNISO
        WRITE(IUNOUT,*) filename//'  : NWISO = ', NWISO
        WRITE(IUNOUT,*) 
     .   'WE MUST HAVE NP2ND = NX+MAX(3*(NNCUT/2),1)+NNISO*NWISO'
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
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
 1014   CONTINUE
C
        CALL EIRENE_MSHPROJ (X1,Y1,X2,Y2,X3,Y3,X4,Y4,
     .                       PUX,PUY,PVX,PVY,NDXA,
     .                       NR1ST,IY)
 1015 CONTINUE
C
      NP=NPOINT(2,NPLP)
      DO 1020 J=1,NDYA+1
        DO I=1,NP
          XPOL(J,I)=XPOL(J,I)*100.
          YPOL(J,I)=YPOL(J,I)*100.
cdr nov 2020: in some versions the next two lines are found
cdr           I do not understand, it seems to be in absolute
cdr           rather than relative units. So I remove these two lines.
cdr       IF (ABS(XPOL(J,I)).LT.5.D-5) XPOL(J,I)=0.
cdr       IF (ABS(YPOL(J,I)).LT.5.D-5) YPOL(J,I)=0.
        END DO
 1020 CONTINUE
      RETURN

  100 WRITE(IUNOUT,*) "COULD NOT OPEN "//trim(filename)//". ",
     w                "SKIP READING THE B2 GEOMETRY" !VK

*//END GEOMD_LINDA//
      END SUBROUTINE EIRENE_GEOMD_LINDA


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
C  GEOMETRY DATA: CELL VERTICES (SONNET ---> EIRENE)
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

    1 continue
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

   99 continue
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
 1014   CONTINUE
C
        CALL EIRENE_MSHPROJ (X1,Y1,X2,Y2,X3,Y3,X4,Y4,PUX,PUY,PVX,PVY,
     .                NDXA,NR1ST,IY)
 1015 CONTINUE
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
      END SUBROUTINE EIRENE_GEOMD_SONNET


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
      END SUBROUTINE EIRENE_MSHPROJ


C
C
      SUBROUTINE EIRENE_INDMAP(FIELD,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .                  NCUTB,NCUTL,NPOINT,NPPLG)
C
C     INDEX MAPPING FOR BRAAMS DATA FIELDS. DATA IN DUMMY ZONES
C     (CUTS OR BOUNDARY ZONES) MAY BE NEEDED AND THUS ARE KEPT
C     AND DUPLICATED IN CASE NCUTL GT NCUTB
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
      INTEGER :: IX, IPART, IY, IFL, IENDD, INB, IINID, IINIV, IENDV
C
C  LOOP FOR THE SPECIES
C
      DO 500 IFL=1,NFLA
C
C  INITIALIZE DUMMY
C
        DO IY=0,NDY+1
          DO IX=0,NDX+1
            DUMMY(IX,IY)=FIELD(IX,IY,IFL)
          END DO
        END DO
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
              DUMMY(IX,IY)=FIELD(INB,IY,IFL)
  213       CONTINUE
            IF (MOD(IPART,3).EQ.0 .AND. IPART.NE.NPPLG) THEN
              DUMMY(IINID,IY) = FIELD(INB  ,IY,IFL)
            ELSE
              DUMMY(IINID,IY) = FIELD(INB+1,IY,IFL)
            ENDIF
            IF (IENDD.NE.IINID)
     >       DUMMY(IENDD,IY) = FIELD(INB+NCUTB,IY,IFL)
  212     CONTINUE
  211   CONTINUE
        DO IY=0,NDYA+1
          DO IX=0,NDXA+1
            FIELD(IX,IY,IFL)=DUMMY(IX,IY)
          END DO
        END DO
C
  500 CONTINUE
      RETURN
C
  990 CONTINUE
      WRITE (iunout,*) 'ERROR IN SUBR. INDMAP: THIS SUBR. IS VALID ONLY'
      WRITE (iunout,*) 'NCUTB>=0 BUT NCUTB = ',NCUTB
      CALL EIRENE_EXIT_OWN(1)
  991 WRITE (iunout,*)
     .  'ERROR IN SUBR. INDMAP: INCONSISTENCY IN NUMBER OF '
      WRITE (iunout,*) 'ZONES PER CUT FROM LINDA GEOMETRY DETECTED.  '
      WRITE (iunout,*) 'NCUTL = ',NCUTL, ' IENDD-IINID+1 = ',
     .                  IENDD-IINID+1
      CALL EIRENE_EXIT_OWN(1)
      END SUBROUTINE EIRENE_INDMAP


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
      INTEGER :: IX, IY, IFL, IENDD, IPART, INB, IINID, IINIV, IENDV
C
C  LOOP OVER THE SPECIES
C
      DO 500 IFL=1,NFLA
C
C  INITIALIZE DUMMY
C
        DO IY=0,NDY+1
          DO IX=0,NDX+1
            DUMMY(IX,IY)=0.
          END DO
        END DO
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
              DUMMY(INB,IY)=FIELD(IX,IY,IFL,ISTR)
  213       CONTINUE
            IF (MOD(IPART,3).EQ.0 .AND. IPART.NE.NPPLG) THEN
              DUMMY(INB  ,IY)=FIELD(IINID,IY,IFL,ISTR)
            ELSE
              DUMMY(INB+1,IY)=FIELD(IINID,IY,IFL,ISTR)
            ENDIF
            IF (IENDD.NE.IINID)
     .          DUMMY(INB+NCUTB,IY)=FIELD(IENDD,IY,IFL,ISTR)
  212     CONTINUE
  211   CONTINUE
        DO IY=0,NDYA+1
          DO IX=0,NDXA+1
            FIELD(IX,IY,IFL,ISTR)=DUMMY(IX,IY)
          END DO
        END DO
C
  500 CONTINUE
      RETURN
C
  990 CONTINUE
      WRITE (iunout,*) 'ERROR IN SUBR. INDMPI: THIS SUBR. IS VALID ONLY'
      WRITE (iunout,*) 'NCUTB>=0 BUT NCUTB = ',NCUTB
      CALL EIRENE_EXIT_OWN(1)
  991 WRITE (iunout,*)
     .  'ERROR IN SUBR. INDMPI: INCONSISTENCY IN NUMBER OF'
      WRITE (iunout,*) 'ZONES PER CUT FROM LINDA GEOMETRY DETECTED. '
      WRITE (iunout,*) 'NCUTL = ',NCUTL, ' IENDD-IINID+1 = ',
     .                  IENDD-IINID+1
      CALL EIRENE_EXIT_OWN(1)
      END SUBROUTINE EIRENE_INDMPI


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
      INTEGER :: ND1, LIM, IF, I1, I2, I3, III, IX, IY
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
!      write (6,*) 'plasm: detected format ', form

      ND1 = NDIMX + 2
      LIM = (ND1/5)*5 - 4
      DUMMY(0:N+1,0:M+1,NF)=0._DP
      DO IF = 1,NDIMF
        DO IY = 0,NDIMY+1
          DO IX = 1,LIM,5
           READ(KARD,FORM,END=500) (DUMMY(-1+IX-1+III,IY,IF),III = 1,5)
          END DO
          IF( (LIM+4).EQ.ND1 )     GOTO 110
          READ(KARD,FORM,END=500) (DUMMY(-1+IX,IY,IF),IX = LIM+5,ND1)
  110   CONTINUE
       END DO
      END DO

      RETURN
 500  backspace kard
      RETURN
*//END PLASM//
      END SUBROUTINE EIRENE_PLASM


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
      INTEGER :: ND1, LIM, IX, IY, III, IFL
C
      ND1 = NDIMX
      LIM = (ND1/5)*5 - 4
      DO  500  IFL = 1,NDIMF
        DO  110  IY = 1,NDIMY
          DO  100  IX = 1,LIM,5
            WRITE(KARD,910) (DUMMY(IX-1+III,IY,IFL,IS),III = 1,5)
  100     CONTINUE
          IF( (LIM+4).EQ.ND1 )   GOTO 110
          WRITE(KARD,910) (DUMMY(IX,IY,IFL,IS),IX = LIM+5,ND1)
  110   CONTINUE
  500 CONTINUE
      RETURN
  910 FORMAT(5(ES16.7E3))
*//END NEUTR//
      END SUBROUTINE EIRENE_NEUTR




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
      USE EIRMOD_CLOGAU

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
              CPMUL => EIRENE_NEW_MULARR()
              CPMUL%IART = IPLS
              CPMUL%ICM = IN
              CPMUL%VALUEM = PAPL(IPLS,IN)*FLXI
              IF (NLSPCSCL_ATM) THEN
                PAPL2(1:NPLS,0:NATM) => PAPL(:,IN)
                CPMUL%VALUAM(1:NATMI) = PAPL2(IPLS,1:NATMI)*FLXI
              END IF
              CPMUL%NXTMUL => PAPLS(ISTRAI)%PMUL
              PAPLS(ISTRAI)%PMUL => CPMUL
            ENDIF
          ENDIF
          IF (LPMPL) THEN
            IF (PMPL(IPLS,IN) .NE. 0.D0) THEN
              CPMUL => EIRENE_NEW_MULARR()
              CPMUL%IART = IPLS
              CPMUL%ICM = IN
              CPMUL%VALUEM = PMPL(IPLS,IN)*FLXI
              IF (NLSPCSCL_MOL) THEN
                PMPL2(1:NPLS,0:NMOL) => PMPL(:,IN)
                CPMUL%VALUMM(1:NMOL) = PMPL2(IPLS,1:NMOL)*FLXI
              END IF
              CPMUL%NXTMUL => PMPLS(ISTRAI)%PMUL
              PMPLS(ISTRAI)%PMUL => CPMUL
            ENDIF
          ENDIF
          IF (LPIPL) THEN
            IF (PIPL(IPLS,IN) .NE. 0.D0) THEN
              CPMUL => EIRENE_NEW_MULARR()
              CPMUL%IART = IPLS
              CPMUL%ICM = IN
              CPMUL%VALUEM = PIPL(IPLS,IN)*FLXI
              IF (NLSPCSCL_ION) THEN
                PIPL2(1:NPLS,0:NION) => PIPL(:,IN)
                CPMUL%VALUIM(1:NION) = PIPL2(IPLS,1:NION)*FLXI
              END IF
              CPMUL%NXTMUL => PIPLS(ISTRAI)%PMUL
              PIPLS(ISTRAI)%PMUL => CPMUL
            ENDIF
          ENDIF

          IF (LEAPL) THEN
            IF (EAPL(IPLS,IN) .NE. 0.D0) THEN
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
            CPSIM => EIRENE_NEW_SIMARR()
            CPSIM%ICS = IN
            CPSIM%VALUES = EAEL(IN)*FLXI
            CPSIM%NXTSIM => EAELS(ISTRAI)%PSIM
            EAELS(ISTRAI)%PSIM => CPSIM
          ENDIF
        ENDIF
        IF (LEMEL) THEN
          IF (EMEL(IN) .NE. 0.D0) THEN
            CPSIM => EIRENE_NEW_SIMARR()
            CPSIM%ICS = IN
            CPSIM%VALUES = EMEL(IN)*FLXI
            CPSIM%NXTSIM => EMELS(ISTRAI)%PSIM
            EMELS(ISTRAI)%PSIM => CPSIM
          ENDIF
        ENDIF
        IF (LEIEL) THEN
          IF (EIEL(IN) .NE. 0.D0) THEN
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
      END SUBROUTINE EIRENE_SAVE_TALLIES

C
C
C
C
