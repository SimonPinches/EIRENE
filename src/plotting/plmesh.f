c  nov 03:  use relative distances to find neighbor segment,
c           otherwise sometimes problems with non-closing polygons encountered.
cdr june 17:  separate WRMESH (WRITING) and PLMESH (PLOTTING).

      SUBROUTINE EIRENE_PLMESH
c  create closed polygonal contours, from the eirene standard and additional surfaces
c  use ILPLG(isurf) flag, from input blocks 3A LEVGEO=3 OR LEVGEO=4,
C                         or certain additional surfaces, input block 3B,
c                         0<RLB<2.
c  This set of closed contours, together with their orientation, can be used
c  in 2D grid generators to produce multiply connected triangular grids.
c  The orientation indicates whether the inner or outer part of a closed contour
c  is a valid computational volume for triangular grid generation (not needed for plotting)


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


      INTEGER, PARAMETER :: MAXPOIN=2000
      REAL(DP) :: partcont(maxpoin,2,2), maxlen
      REAL(DP) :: XPE, YPE, HELP, XT, YT, PHI1, X1, X2, Y1, Y2, PHI2
      REAL(DP) :: DISTQI, DISTQJ1, DISTQJ2, YMN
      INTEGER  :: ICONT, IPOIN, IWST, IWEN, IWL, IWP,
     .            IWAN, IMN, I, NCONT, J, IUHR, ISTORE, IP, IH, IFOUND,
     .            ICO, IPO, IN, IS, IS1, ITRI, INBT, INBS
      INTEGER  :: IDIAG(MAXPOIN),irip(maxpoin,2)
      REAL(SP) :: xmin,xmax,ymin,ymax,deltax,deltay,delta,xcm,ycm
      REAL(SP) :: XP,YP
      LOGICAL  :: LCLOSED, LFOUND
      LOGICAL, ALLOCATABLE :: FOUND(:,:)

C INITIALISIERUNG DER PLOTDATEN
      xmin = CH2X0-CH2MX
      ymin = CH2Y0-CH2MY
      xmax = CH2X0+CH2MX
      ymax = CH2Y0+CH2MY
      deltax = abs(xmax-xmin)
      deltay = abs(ymax-ymin)
      delta = max(deltax,deltay)
      xcm = 24. * deltax/delta
      ycm = 24. * deltay/delta

C ANZAHL DER KONTOUREN BESTIMMEN
C ILPLG WIRD IM INPUT BLOCK 3 EINGELESEN
      CALL EIRENE_LEER(2)
      WRITE (iunout,*) 'SUBROUTINE PLMESH CALLED'
      CALL EIRENE_LEER(1)

      NCONT = 0
      DO I=1,NLIMI
        NCONT = MAX(NCONT,ABS(ILPLG(I)))
      ENDDO
      DO I=NLIM+1,NLIM+NSTSI
        NCONT = MAX(NCONT,ABS(ILPLG(I)))
      ENDDO

      if (ncont == 0) return
      
      IF (.NOT.ALLOCATED(NCONPOINT)) THEN
        ALLOCATE (NCONPOINT(NCONT))
        ALLOCATE (XCONTOUR(MAXPOIN,NCONT))
        ALLOCATE (YCONTOUR(MAXPOIN,NCONT))
      ENDIF
      NCONPOINT = 0
      ICO = 0

      call grnxtf
      call grsclc(3.,3.,3.+real(xcm,kind(1.e0)),3.+real(ycm,kind(1.e0)))
      call grsclv(real(xmin,kind(1.e0)),real(ymin,kind(1.e0)),
     .            real(xmax,kind(1.e0)),real(ymax,kind(1.e0)))

      DO ICONT = 1,NCONT
        IPOIN = 0
        MAXLEN = 0.
        irip=0
C AKTUELLE KONTOUR BESTIMMEN, STUECKE MIT ILPLG=ICONT GEHOEREN ZUR
C AKTUELLEN CONTOUR, ANFANGS UND ENDPUNKT DIESES STUECKES WERDEN AUF
C PARTCONT GESPEICHERT

c  ADDITIONAL SURFACES
        DO I=1,NLIMI
          IF (ABS(ILPLG(I)) .EQ. ICONT) THEN
            IUHR=ILPLG(I)
C 0 < RLB(I) < 2
C 2-PUNKT OPTION WIRD IM TIMEA0 AUF RLB=1 ZURUECKGEFUEHRT
            IF ((RLB(I) .GT. 0.) .AND. (RLB(I) .LT. 2.) .AND.
     >          (P3(1,I) .EQ. 1.D55 .OR. P3(2,I) .EQ. 1.D55
     >          .OR. P3(3,I) .EQ. 1.D55)) THEN
              IPOIN = IPOIN + 1
              IF (IPOIN.GT.MAXPOIN) THEN
                WRITE(IUNOUT,*)
     .           'INSUFFICIENT NUMBER OF POINTS FOR CONTOUR ',
     .            ICONT
                WRITE(IUNOUT,*)
     .           'INCREASE VALUE OF MAXPOIN IN wrmesh.F'
                WRITE(IUNOUT,*)
     .           'CURRENTLY MAXPOIN = ', MAXPOIN
                CALL EIRENE_EXIT_OWN(1)
              ENDIF
              IF (A3LM(I) .EQ. 0.) THEN
C               X,Y-KOORDINATEN
                PARTCONT(IPOIN,1,1) = P1(1,I)
                PARTCONT(IPOIN,1,2) = P1(2,I)
                PARTCONT(IPOIN,2,1) = P2(1,I)
                PARTCONT(IPOIN,2,2) = P2(2,I)
                idiag(ipoin)=i
              ELSEIF (A2LM(I) .EQ. 0.) THEN
C               X,Z-KOORDINATEN
                PARTCONT(IPOIN,1,1) = P1(1,I)
                PARTCONT(IPOIN,1,2) = P1(3,I)
                PARTCONT(IPOIN,2,1) = P2(1,I)
                PARTCONT(IPOIN,2,2) = P2(3,I)
                idiag(ipoin)=i
              ELSEIF (A1LM(I) .EQ. 0.) THEN
C               Y,Z-KOORDINATEN
                PARTCONT(IPOIN,1,1) = P1(2,I)
                PARTCONT(IPOIN,1,2) = P1(3,I)
                PARTCONT(IPOIN,2,1) = P2(2,I)
                PARTCONT(IPOIN,2,2) = P2(3,I)
                idiag(ipoin)=i
              ENDIF
              maxlen = maxlen +
     >               sqrt((partcont(ipoin,1,1)-partcont(ipoin,2,1))**2
     >                   +(partcont(ipoin,1,2)-partcont(ipoin,2,2))**2)
            ELSE
C  ERROR
              WRITE(iunout,'(a,f11.4,2i4)')
     >         'FALSCHE ANGABE FUER RLB, RLB = ',RLB(I),ILPLG(I),I
            ENDIF
          ENDIF
        ENDDO

        select case (LEVGEO)
        case (3)
        DO I=1,NSTSI
          IF (ABS(ILPLG(NLIM+I)) .EQ. ICONT) THEN
            IUHR=ILPLG(NLIM+I)
            IF (INUMP(I,2) .NE. 0) THEN
C  POLOIDAL SURFACES
              DO J=IRPTA(I,1),IRPTE(I,1)-1
                IF ((XPOL(J,INUMP(I,2)) .NE. XPOL(J+1,INUMP(I,2))) .OR.
     >              (YPOL(J,INUMP(I,2)) .NE. YPOL(J+1,INUMP(I,2)))) THEN
                  IPOIN = IPOIN + 1
                  IF (IPOIN.GT.MAXPOIN) THEN
                    WRITE(IUNOUT,*)
     .               'INSUFFICIENT NUMBER OF POINTS FOR CONTOUR ',
     .                ICONT
                    WRITE(IUNOUT,*)
     .               'INCREASE VALUE OF MAXPOIN IN wrmesh.F'
                    WRITE(IUNOUT,*)
     .               'CURRENTLY MAXPOIN = ', MAXPOIN
                    CALL EIRENE_EXIT_OWN(1)
                  ENDIF
                  PARTCONT(IPOIN,1,1) = XPOL(J,INUMP(I,2))
                  PARTCONT(IPOIN,1,2) = YPOL(J,INUMP(I,2))
                  PARTCONT(IPOIN,2,1) = XPOL(J+1,INUMP(I,2))
                  PARTCONT(IPOIN,2,2) = YPOL(J+1,INUMP(I,2))
                idiag(ipoin)=-i
                irip(ipoin,1)=j
                irip(ipoin,2)=INUMP(I,2)
              maxlen = maxlen +
     >               sqrt((partcont(ipoin,1,1)-partcont(ipoin,2,1))**2
     >                   +(partcont(ipoin,1,2)-partcont(ipoin,2,2))**2)
                ENDIF
              ENDDO
            ELSEIF (INUMP(I,1) .NE. 0) THEN
C  RADIAL SURFACES
              DO J=IRPTA(I,2),IRPTE(I,2)-1
                IF ((XPOL(INUMP(I,1),J) .NE. XPOL(INUMP(I,1),J+1)) .OR.
     >              (YPOL(INUMP(I,1),J) .NE. YPOL(INUMP(I,1),J+1))) THEN
                  IPOIN = IPOIN + 1
                  IF (IPOIN.GT.MAXPOIN) THEN
                    WRITE(IUNOUT,*)
     .               'INSUFFICIENT NUMBER OF POINTS FOR CONTOUR ',
     .                ICONT
                    WRITE(IUNOUT,*)
     .               'INCREASE VALUE OF MAXPOIN IN wrmesh.F'
                    WRITE(IUNOUT,*)
     .               'CURRENTLY MAXPOIN = ', MAXPOIN
                    CALL EIRENE_EXIT_OWN(1)
                  ENDIF
                  PARTCONT(IPOIN,1,1) = XPOL(INUMP(I,1),J)
                  PARTCONT(IPOIN,1,2) = YPOL(INUMP(I,1),J)
                  PARTCONT(IPOIN,2,1) = XPOL(INUMP(I,1),J+1)
                  PARTCONT(IPOIN,2,2) = YPOL(INUMP(I,1),J+1)
                  idiag(ipoin)=-i
                  irip(ipoin,1)=INUMP(I,1)
                  irip(ipoin,2)=j
                  maxlen = maxlen +
     >               sqrt((partcont(ipoin,1,1)-partcont(ipoin,2,1))**2
     >                   +(partcont(ipoin,1,2)-partcont(ipoin,2,2))**2)
                ENDIF
              ENDDO
            ELSE
C  ERROR
              WRITE(iunout,*) 'CASE NOT FORESEEN: INUMP: ',
     >                     (INUMP(I,J),J=1,3)
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
                  IUHR=ILPLG(IN)
                  IS1 = IS+1
                  IF (IS1 > 3) IS1=1
                  IPOIN = IPOIN + 1
                  IF (IPOIN.GT.MAXPOIN) THEN
                    WRITE(IUNOUT,*)
     .               'INSUFFICIENT NUMBER OF POINTS FOR CONTOUR ',
     .                ICONT
                    WRITE(IUNOUT,*)
     .               'INCREASE VALUE OF MAXPOIN IN wrmesh.F'
                    WRITE(IUNOUT,*)
     .               'CURRENTLY MAXPOIN = ', MAXPOIN
                    CALL EIRENE_EXIT_OWN(1)
                  ENDIF
                  PARTCONT(IPOIN,1,1) = XTRIAN(NECKE(IS,ITRI))
                  PARTCONT(IPOIN,1,2) = YTRIAN(NECKE(IS,ITRI))
                  PARTCONT(IPOIN,2,1) = XTRIAN(NECKE(IS1,ITRI))
                  PARTCONT(IPOIN,2,2) = YTRIAN(NECKE(IS1,ITRI))
                  idiag(ipoin)=IN
                  irip(ipoin,1)=itri
                  irip(ipoin,2)=is
                  maxlen = maxlen +
     >               sqrt((partcont(ipoin,1,1)-partcont(ipoin,2,1))**2
     >                   +(partcont(ipoin,1,2)-partcont(ipoin,2,2))**2)
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

              ih=idiag(i+1)
              idiag(i+1)=idiag(j)
              idiag(j)=ih

              ih=irip(i+1,1)
              irip(i+1,1)=irip(j,1)
              irip(j,1)=ih
              ih=irip(i+1,2)
              irip(i+1,2)=irip(j,2)
              irip(j,2)=ih
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

              ih=idiag(i+1)
              idiag(i+1)=idiag(j)
              idiag(j)=ih

              ih=irip(i+1,1)
              irip(i+1,1)=irip(j,1)
              irip(j,1)=ih
              ih=irip(i+1,2)
              irip(i+1,2)=irip(j,2)
              irip(j,2)=ih
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
        XP = PARTCONT(1,1,1)
        YP = PARTCONT(1,1,2)
        call grjmp(REAL(XP,KIND(1.E0)),REAL(YP,KIND(1.E0)))
        DO I=2,IPOIN
          XP = PARTCONT(I,1,1)
          YP = PARTCONT(I,1,2)
          call grdrw(REAL(XP,KIND(1.E0)),REAL(YP,KIND(1.E0)))
        ENDDO
c  last point on contour
        IF (LCLOSED) THEN
          XP = PARTCONT(1,1,1)
          YP = PARTCONT(1,1,2)
          call grDRW(REAL(XP,KIND(1.E0)),REAL(YP,KIND(1.E0)))
        ELSE
          XP = PARTCONT(IPOIN,2,1)
          YP = PARTCONT(IPOIN,2,2)
          call grDRW(REAL(XP,KIND(1.E0)),REAL(YP,KIND(1.E0)))
        END IF

1000    CONTINUE
      ENDDO    ! END OF DO ICONT.... LOOP


c  re-initialize gr plot software for next picture
      call grnwpn(1)
      call grnxtf

cdr
      if (allocated(nconpoint)) then
        DEALLOCATE (NCONPOINT)
        DEALLOCATE (XCONTOUR)
        DEALLOCATE (YCONTOUR)
      endif

      return
      END
