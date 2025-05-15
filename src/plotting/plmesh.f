c  nov 03:  use relative distances to find neighbor segment,
c           otherwise sometimes problems with non-closing polygons encountered.
cdr june 17:  separate WRMESH (WRITING) and PLMESH (PLOTTING).

      SUBROUTINE EIRENE_PLMESH
c  create closed polygonal contours, from the eirene standard and
c  additional surfaces
c  use ILPLG(isurf) flag, from input blocks 3A LEVGEO=3 OR LEVGEO=4,
C                         or certain additional surfaces, input block 3B,
c                         0<RLB<2.
c  This set of closed contours, together with their orientation, can be used
c  in 2D grid generators to produce multiply connected triangular grids.
c  The orientation indicates whether the inner or outer part of a closed contour
c  is a valid computational volume for triangular grid generation (not needed
c  for plotting)


c  EIRENE_WRMESH: Write closed contours onto output stream 78+ifoff.
c  EIRENE_PLMESH: plots these contours, using GR plot software.


      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CADGEO
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CPLOT
      USE EIRMOD_CPOLYG
      USE EIRMOD_CGEOM
      USE EIRMOD_CLGIN
      USE EIRMOD_CGRPTL
      USE EIRMOD_CTRIG
      USE EIRMOD_CGRID
      USE EIRMOD_CTRCEI
      IMPLICIT NONE

      REAL(DP), ALLOCATABLE :: partcont(:,:,:)
      INTEGER, SAVE :: MAXPOIN=2000
      REAL(DP) :: XPE, YPE, HELP
      REAL(DP) :: DISTQI, DISTQJ1, DISTQJ2
      INTEGER  :: ICONT, IPOIN, I, J,
     .            IN, IS, IS1, ITRI, INBS, INBT, IFOUND, NCONT
      REAL(SP) :: xmin,xmax,ymin,ymax,deltax,deltay,delta,xcm,ycm
      REAL(SP) :: XP,YP
      LOGICAL  :: LCLOSED, LFOUND
      LOGICAL, ALLOCATABLE :: FOUND(:,:)
      EXTERNAL :: EIRENE_LEER
      EXTERNAL :: GRDRW, GRJMP, GRNXTF, GRNWPN, GRSCLC, GRSCLV

C INITIALISIERUNG DER PLOTDATEN
      xmin = REAL(CH2X0-CH2MX,SP)
      ymin = REAL(CH2Y0-CH2MY,SP)
      xmax = REAL(CH2X0+CH2MX,SP)
      ymax = REAL(CH2Y0+CH2MY,SP)
      deltax = abs(xmax-xmin)
      deltay = abs(ymax-ymin)
      delta = max(deltax,deltay)
      xcm = 24._SP * deltax/delta
      ycm = 24._SP * deltay/delta

C ANZAHL DER KONTOUREN BESTIMMEN
C ILPLG WIRD IM INPUT BLOCK 3 EINGELESEN
      CALL EIRENE_LEER(2)
      WRITE (iunout,*) 'SUBROUTINE PLMESH CALLED'
      CALL EIRENE_LEER(1)

      ALLOCATE (partcont(maxpoin,2,2))

      NCONT = 0
      DO I=1,NLIMI
        NCONT = MAX(NCONT,ABS(ILPLG(I)))
      ENDDO
      DO I=NLIM+1,NLIM+NSTSI
        NCONT = MAX(NCONT,ABS(ILPLG(I)))
      ENDDO

      if (ncont == 0) return

      call grnxtf
      call grsclc(3._SP,3._SP,3._SP+xcm,3._SP+ycm)
      call grsclv(xmin,ymin,xmax,ymax)

      DO ICONT = 1,NCONT
        IPOIN = 0
C AKTUELLE KONTOUR BESTIMMEN, STUECKE MIT ILPLG=ICONT GEHOEREN ZUR
C AKTUELLEN KONTOUR, ANFANGS UND ENDPUNKT DIESES STUECKES WERDEN AUF
C PARTCONT GESPEICHERT

        select case (LEVGEO)
        case (3)
c  ADDITIONAL SURFACES
          DO I=1,NLIMI
            IF (ABS(ILPLG(I)) .EQ. ICONT) THEN
C 0 < RLB(I) < 2
C 2-PUNKT OPTION WIRD IM TIMEA0 AUF RLB=1 ZURUECKGEFUEHRT
              IF ((RLB(I) .GT. 0._DP) .AND. (RLB(I) .LT. 2._DP) .AND.
     >            (P3(1,I) .EQ. 1.D55 .OR. P3(2,I) .EQ. 1.D55
     >            .OR. P3(3,I) .EQ. 1.D55)) THEN
                IPOIN = IPOIN + 1
                IF (IPOIN.GT.MAXPOIN) THEN
                  WRITE(IUNOUT,*)
     .             'INSUFFICIENT NUMBER OF POINTS FOR CONTOUR ',
     .              ICONT
                  WRITE(IUNOUT,*) 'INCREASE VALUE OF MAXPOIN'
                  CALL EIRENE_EXTEND_ARRAY
                  WRITE(IUNOUT,*) 'MAXPOIN SET TO = ', MAXPOIN
                ENDIF
                IF (A3LM(I) .EQ. 0._DP) THEN
C               X,Y-KOORDINATEN
                  PARTCONT(IPOIN,1,1) = P1(1,I)
                  PARTCONT(IPOIN,1,2) = P1(2,I)
                  PARTCONT(IPOIN,2,1) = P2(1,I)
                  PARTCONT(IPOIN,2,2) = P2(2,I)
                ELSEIF (A2LM(I) .EQ. 0._DP) THEN
C               X,Z-KOORDINATEN
                  PARTCONT(IPOIN,1,1) = P1(1,I)
                  PARTCONT(IPOIN,1,2) = P1(3,I)
                  PARTCONT(IPOIN,2,1) = P2(1,I)
                  PARTCONT(IPOIN,2,2) = P2(3,I)
                ELSEIF (A1LM(I) .EQ. 0._DP) THEN
C               Y,Z-KOORDINATEN
                  PARTCONT(IPOIN,1,1) = P1(2,I)
                  PARTCONT(IPOIN,1,2) = P1(3,I)
                  PARTCONT(IPOIN,2,1) = P2(2,I)
                  PARTCONT(IPOIN,2,2) = P2(3,I)
                ENDIF
              ELSE
C  ERROR
                WRITE(iunout,'(a,f11.4,2i4)')
     >           'FALSCHE ANGABE FUER RLB, RLB = ',RLB(I),ILPLG(I),I
              ENDIF
            ENDIF
          ENDDO

          DO I=1,NSTSI
            IF (ABS(ILPLG(NLIM+I)) .EQ. ICONT) THEN
              IF (INUMP(I,2) .NE. 0) THEN
C  POLOIDAL SURFACES
                DO J=IRPTA(I,1),IRPTE(I,1)-1
                  IF ((XPOL(J,INUMP(I,2)).NE.XPOL(J+1,INUMP(I,2))) .OR.
     >                (YPOL(J,INUMP(I,2)).NE.YPOL(J+1,INUMP(I,2)))) THEN
                    IPOIN = IPOIN + 1
                    IF (IPOIN.GT.MAXPOIN) THEN
                      WRITE(IUNOUT,*)
     .                 'INSUFFICIENT NUMBER OF POINTS FOR CONTOUR ',
     .                  ICONT
                      WRITE(IUNOUT,*) 'INCREASE VALUE OF MAXPOIN'
                      CALL EIRENE_EXTEND_ARRAY
                      WRITE(IUNOUT,*) 'MAXPOIN SET TO = ', MAXPOIN
                    ENDIF
                    PARTCONT(IPOIN,1,1) = XPOL(J,INUMP(I,2))
                    PARTCONT(IPOIN,1,2) = YPOL(J,INUMP(I,2))
                    PARTCONT(IPOIN,2,1) = XPOL(J+1,INUMP(I,2))
                    PARTCONT(IPOIN,2,2) = YPOL(J+1,INUMP(I,2))
                  ENDIF
                ENDDO
              ELSEIF (INUMP(I,1) .NE. 0) THEN
C  RADIAL SURFACES
                DO J=IRPTA(I,2),IRPTE(I,2)-1
                  IF ((XPOL(INUMP(I,1),J).NE.XPOL(INUMP(I,1),J+1)) .OR.
     >                (YPOL(INUMP(I,1),J).NE.YPOL(INUMP(I,1),J+1))) THEN
                    IPOIN = IPOIN + 1
                    IF (IPOIN.GT.MAXPOIN) THEN
                      WRITE(IUNOUT,*)
     .                 'INSUFFICIENT NUMBER OF POINTS FOR CONTOUR ',
     .                  ICONT
                      WRITE(IUNOUT,*) 'INCREASE VALUE OF MAXPOIN'
                      CALL EIRENE_EXTEND_ARRAY
                      WRITE(IUNOUT,*) 'MAXPOIN SET TO = ', MAXPOIN
                    ENDIF
                    PARTCONT(IPOIN,1,1) = XPOL(INUMP(I,1),J)
                    PARTCONT(IPOIN,1,2) = YPOL(INUMP(I,1),J)
                    PARTCONT(IPOIN,2,1) = XPOL(INUMP(I,1),J+1)
                    PARTCONT(IPOIN,2,2) = YPOL(INUMP(I,1),J+1)
                  ENDIF
                ENDDO
              ELSE
C  ERROR
                WRITE(iunout,*) 'CASE NOT FORESEEN: INUMP: ',
     >                       (INUMP(I,J),J=1,3)
              ENDIF
            ENDIF
          ENDDO

        case (4)
C  TRIANGLE SIDES
          ALLOCATE (FOUND(1:3,1:NTRII))
          FOUND = .FALSE.
          DO ITRI = 1, NTRII
            DO IS = 1, 3
              IN=INMTI(IS,ITRI)
              LFOUND=FOUND(IS,ITRI)
              IF (IN /= 0 .AND. .NOT.LFOUND) THEN
                IF (ABS(ILPLG(IN)) == ICONT) THEN
                  IS1 = IS+1
                  IF (IS1 > 3) IS1=1
                  IPOIN = IPOIN + 1
                  IF (IPOIN.GT.MAXPOIN) THEN
                    WRITE(IUNOUT,*)
     .               'INSUFFICIENT NUMBER OF POINTS FOR CONTOUR ',
     .                ICONT
                    WRITE(IUNOUT,*) 'INCREASE VALUE OF MAXPOIN'
                    CALL EIRENE_EXTEND_ARRAY
                    WRITE(IUNOUT,*) 'MAXPOIN SET TO = ', MAXPOIN
                  ENDIF
                  PARTCONT(IPOIN,1,1) = XTRIAN(NECKE(IS,ITRI))
                  PARTCONT(IPOIN,1,2) = YTRIAN(NECKE(IS,ITRI))
                  PARTCONT(IPOIN,2,1) = XTRIAN(NECKE(IS1,ITRI))
                  PARTCONT(IPOIN,2,2) = YTRIAN(NECKE(IS1,ITRI))
                  FOUND(IS,ITRI)=.TRUE.
                  INBT=NCHBAR(IS,ITRI)
cwdk Make sure the corresponding side of the neighboring triangle is not found again
                  IF (INBT.GT.0) THEN
                    INBS=NSEITE(IS,ITRI)
                    FOUND(INBS,INBT)=.TRUE.
                  ENDIF
                END IF
              END IF
            END DO
          END DO
          DEALLOCATE (FOUND)
        end select

        IF (IPOIN.LE.0) THEN
          WRITE(iunout,*) 'CONTOUR ',ICONT,' NOT FOUND'
          GOTO 1000
        ENDIF

C STUECKE DER AKTUELLEN KONTOUR WERDEN SORTIERT
        DO I=1,IPOIN-1
          XPE = PARTCONT(I,2,1)
          YPE = PARTCONT(I,2,2)
          DISTQI=(PARTCONT(I,2,1)-PARTCONT(I,1,1))**2+
     .           (PARTCONT(I,2,2)-PARTCONT(I,1,2))**2
          IFOUND=0
          DO J=I+1,IPOIN
            IF (IFOUND.EQ.1) CYCLE
            DISTQJ1=(XPE-PARTCONT(J,1,1))**2+
     .              (YPE-PARTCONT(J,1,2))**2
            DISTQJ2=(XPE-PARTCONT(J,2,1))**2+
     .              (YPE-PARTCONT(J,2,2))**2
            IF (DISTQJ1/DISTQI.LE.1.D-10) THEN
              IFOUND=1
              HELP = PARTCONT(I+1,1,1)
              PARTCONT(I+1,1,1) = PARTCONT(J,1,1)
              PARTCONT(J,1,1) = HELP
              HELP = PARTCONT(I+1,1,2)
              PARTCONT(I+1,1,2) = PARTCONT(J,1,2)
              PARTCONT(J,1,2) = HELP

              HELP = PARTCONT(I+1,2,1)
              PARTCONT(I+1,2,1) = PARTCONT(J,2,1)
              PARTCONT(J,2,1) = HELP
              HELP = PARTCONT(I+1,2,2)
              PARTCONT(I+1,2,2) = PARTCONT(J,2,2)
              PARTCONT(J,2,2) = HELP

            ELSEIF (DISTQJ2/DISTQI.LE.1.D-10) THEN

              IFOUND=1
              HELP = PARTCONT(J,1,1)
              PARTCONT(J,1,1) = PARTCONT(J,2,1)
              PARTCONT(J,2,1) = HELP
              HELP = PARTCONT(J,1,2)
              PARTCONT(J,1,2) = PARTCONT(J,2,2)
              PARTCONT(J,2,2) = HELP

              HELP = PARTCONT(I+1,1,1)
              PARTCONT(I+1,1,1) = PARTCONT(J,1,1)
              PARTCONT(J,1,1) = HELP
              HELP = PARTCONT(I+1,1,2)
              PARTCONT(I+1,1,2) = PARTCONT(J,1,2)
              PARTCONT(J,1,2) = HELP

              HELP = PARTCONT(I+1,2,1)
              PARTCONT(I+1,2,1) = PARTCONT(J,2,1)
              PARTCONT(J,2,1) = HELP
              HELP = PARTCONT(I+1,2,2)
              PARTCONT(I+1,2,2) = PARTCONT(J,2,2)
              PARTCONT(J,2,2) = HELP

            ENDIF
          ENDDO

        ENDDO

        IF ((PARTCONT(1,1,1) .NE. PARTCONT(IPOIN,2,1)) .OR.
     >      (PARTCONT(1,1,2) .NE. PARTCONT(IPOIN,2,2))) THEN
          WRITE(iunout,*) 'CONTOUR ',ICONT,' IS NOT CLOSED'
          LCLOSED = .FALSE.
        ELSE
          WRITE(iunout,*) 'CLOSED CONTOUR ',ICONT
          LCLOSED = .TRUE.
        ENDIF


c  PLOT CONTOUR ICONT

        call grnwpn(icont)
c   first point on contour
        XP = REAL(PARTCONT(1,1,1),SP)
        YP = REAL(PARTCONT(1,1,2),SP)
        call grjmp(XP,YP)
        DO I=2,IPOIN
          XP = REAL(PARTCONT(I,1,1),SP)
          YP = REAL(PARTCONT(I,1,2),SP)
          call grdrw(XP,YP)
        ENDDO
c  last point on contour
        IF (LCLOSED) THEN
          XP = REAL(PARTCONT(1,1,1),SP)
          YP = REAL(PARTCONT(1,1,2),SP)
          call grdrw(XP,YP)
        ELSE
          XP = REAL(PARTCONT(IPOIN,2,1),SP)
          YP = REAL(PARTCONT(IPOIN,2,2),SP)
          call grdrw(XP,YP)
        END IF

 1000   CONTINUE
      ENDDO    ! END OF DO ICONT.... LOOP

c  re-initialize gr plot software for next picture
      call eirene_leer(1)
      call grnwpn(1)
      call grnxtf

      DEALLOCATE (partcont)

      RETURN

      CONTAINS

      SUBROUTINE EIRENE_EXTEND_ARRAY

      IMPLICIT NONE
      REAL(DP), ALLOCATABLE :: pc(:,:,:)
      INTEGER :: NEWPOIN

      ALLOCATE(PC(MAXPOIN,2,2))
      PC = PARTCONT

      DEALLOCATE (PARTCONT)

      NEWPOIN = MAXPOIN + 2000
      ALLOCATE (partcont(newpoin,2,2))

      partcont(1:maxpoin,:,:) = pc(1:maxpoin,:,:)

      maxpoin = newpoin

      deallocate(pc)

      return
      END SUBROUTINE EIRENE_EXTEND_ARRAY

      END SUBROUTINE EIRENE_PLMESH
