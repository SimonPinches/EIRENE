C
c sept 2014 : some indexing unified, comments added.
c             ngitt (storage for spatial resolution on single surface) now evaluated in calling program input.f
c
!pb  3.12.06: allow INDGRD /= 6 for NLTET option
!pb  3.12.06: specify NGITT in case of NLTET
!pb           initialize XDIFF=0
!pb 22.03.07: LEVGEO=6 --> LEVGEO=10
!pb 28.06.10: consistency check for LEVGEO=3 or LEVGEO=4 introduced
!             stop run if cell side is transparent but no neighbor cell
!             is defined
!pb   ??      use nrplg rather than np2nd for 1D radial polygons,
!             to allow 1D levgeo=3 runs
!dr 17.01.14  test printout removed, some comments added
cdr 19.07.16  arguments corrected in call to grid_1 for levgeo=1,
cdr           to allow also there an outer vacuum (void) zone.
cdr           tested, ok.
cdr july 17:  initialize NCORNER also in case LEVGEO=10.
cdr           For several other LEVGEO options NCORNER is set here,
cdr           for some other is routine sneigh.f (geometry module)

      SUBROUTINE EIRENE_GRID (IND)

C  SET STANDARD GRIDS AND RELATED DATA
C    INPUT:   IND
C    OUTPUT:  IN MODULES
C    IND=1:  1ST GRID, X OR RADIAL COORDINATE, AS WELL AS TRIANGULAR (LEVGEO=4) AND TETRAHEDRON (LEVGEO=5) GRIDS.
C    IND=2:  2ND GRID, Y OR POLOIDAL COORDINATE
C    IND=3:  3RD GRID, Z OR TOROIDAL COORDINATE

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CADGEO
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CINIT
      USE EIRMOD_CUPD
      USE EIRMOD_CPOLYG
      USE EIRMOD_CGRID
      USE EIRMOD_CTRCEI
      USE EIRMOD_CGEOM
      USE EIRMOD_CTETRA
      USE EIRMOD_CLGIN
      USE EIRMOD_CTRIG

      IMPLICIT NONE

      REAL(DP) :: PC1(3), EDGELEN(6)
      REAL(DP) :: ELPARM, X1, X2, SY, Y1, Y2, SX, AELL, X3, Y3, X4, Y4,
     .          VPXX, RN, RRN, FN, VPYY, PLABS, XNORM, VPX, VPY, QUOTI,
     .          GESFL, FRING, CONST, RRR, FL, FR, RL, RR, RRL, XD,
     .          PLEN, XDIFF, RORIG, XS3, PLABS2, PLABS3, XD1, YD,
     .          XS, PLABS1, YD1, XS2, XD3, YD3, XS1, XD2, YD2, R, PIN,
     .          POUT, EX1, XX1, XX2, YY1, YY2, DSD, COM, S, SQ,
     .          DP1, DS1, DS2
      REAL(DP), EXTERNAL :: EIRENE_ARTRI3
      INTEGER :: ITSIDE(3,4)
      INTEGER :: IP, IRP, IPP, IT, KDN, KUP, NCELL, IR, I, K, IUP, IDN,
     .           J, IND, NLOCAL, ND, IC3, IC4, ITET, IC1, IC2, IC, NLJ,
     .           IFLAG, ISTS, IECKE2, MSURFG, IS, IT1, NCELL1, NSRFTR,
     .           NT, IS1, IM, IMP, JLIM
      LOGICAL :: LERROR
      LOGICAL, ALLOCATABLE :: VISITED(:,:)
!pb
      TYPE(TRI_ELEM), POINTER :: CUR

      DATA ITSIDE /1,2,3,
     .             1,4,2,
     .             2,4,3,
     .             3,4,1/
C STATEMENT FUNCTION FOR GRID PARAMETERS FOR LEVGEO=2 OPTION
      ELPARM(R,PIN,POUT,EX1)=(PIN-POUT)*(1.-R**EX1)**1.+POUT
C
      SELECT CASE (IND)
C
C   RADIAL GRID
C
      CASE (1) ! IND
C
      IF (NR1ST.LT.2) RETURN
C
      SELECT CASE (LEVGEO)
      CASE (1) ! LEVGEO
C
C  GRID DATA GENERATION FOR LEVGEO.EQ.1
C
        IF (INDGRD(IND).LE.4) THEN
          IF (RRA.GT.RAA) THEN
            NLOCAL=NR1STM
          ELSE
            NLOCAL=NR1ST
          ENDIF
C** USE ONE OF THE EIRENE DEFAULT GRID OPTIONS
          CALL EIRENE_GRID_1(RSURF,NR1ST,NRSEP,NLOCAL,RIA,RGA,RAA,RRA,1)
        ELSEIF (INDGRD(IND).EQ.5) THEN
C** TAKE RADIAL GRID DATA FROM USER-SUPPLIED SUBROUTINE
C         CALL PROUSR (RSURF,2+4*NPLS+3,0._DP,0._DP,0._DP,0._DP,
C    .                 0._DP,0._DP,0._DP,NR1ST)
        ELSEIF (INDGRD(IND).EQ.6) THEN
C** TAKE RADIAL GRID DATA FROM INTERFACE
C         CALL PROFR (RSURF,6+5*NPLS+NAIN,1,1,NR1ST)
        ENDIF
C  SET DERIVED GRID DATA FOR LEVGEO = 1 OPTION
C
C  NOTHING TO BE DONE HERE
C
        IF (TRCGRD) THEN
          CALL EIRENE_LEER(1)
          WRITE (iunout,*) 'GRIDPOINTS IN X DIRECTION '
          CALL EIRENE_LEER(1)
          CALL EIRENE_MASRR1('  N, RSURF ',RSURF,NR1ST,3)
          CALL EIRENE_LEER(2)
        ENDIF
C
      CASE (2) ! LEVGEO
C
C  GRID DATA GENERATION FOR LEVGEO.EQ.2
C
        IF (INDGRD(IND).LE.4) THEN
C** USE ONE OF THE EIRENE DEFAULT GRID OPTIONS
          IF (RRA.GT.RAA) THEN
            NLOCAL=NR1STM
            RSURF(NR1ST)=RRA
            ELL(NR1ST)=ELLCH
            EP1(NR1ST)=EP1CH
            TRI(NR1ST)=TRICH
          ELSE
            NLOCAL=NR1ST
          ENDIF
C
          IF (INDGRD(IND).EQ.1) THEN
C** EQUIDISTANT IN RADIUS
            ND=NRSEP
            DO 105 J=1,ND
              RSURF(J)=RIA+DBLE(J-1)/DBLE(ND-1)*(RGA-RIA)
  105       CONTINUE
            DO 106 J=ND+1,NLOCAL
              RSURF(J)=RGA+DBLE(J-ND)/DBLE(NLOCAL-ND)*(RAA-RGA)
  106       CONTINUE
          ELSEIF (INDGRD(IND).EQ.2) THEN
C** RADIAL GRID WITH CONSTANT AREA
            IF (NLCRC) THEN
              GESFL=(RAA*RAA-RIA*RIA)*PIA
              FRING=GESFL/(NLOCAL-1)
              RSURF(1)=RIA
              RSURF(NLOCAL)=RAA
              DO 103 J=2,NLOCAL-1
                RSURF(J)=RIA+SQRT((J-1)*FRING*PIAI)
  103         CONTINUE
            ELSEIF (NLELL) THEN
              GESFL=(RAA*RAA*ELLOT-RIA*RIA*ELLIN)*PIA
              FRING=GESFL/(NLOCAL-1)
              RSURF(1)=RIA
              RSURF(NLOCAL)=RAA
              DO 104 J=2,NLOCAL-1
C  SOLVE RSURF**2*ELL(RSURF)-(J-1)*FRING/PIA=0., RSURF(J-1)<RSURF<RAA
                CONST=(J-1)*FRING*PIAI
                RL=RSURF(J-1)+EPS30
                RR=RAA
  108           RRL=(RL-RSURF(1))/(RSURF(NLOCAL)-RSURF(1))
  109           RRR=(RR-RSURF(1))/(RSURF(NLOCAL)-RSURF(1))
                FL=RL**2*ELPARM(RRL,ELLIN,ELLOT,EXELL)-CONST
                FR=RR**2*ELPARM(RRR,ELLIN,ELLOT,EXELL)-CONST
                QUOTI=(RR-RL)/(FR-FL)
                RN=-FL*QUOTI+RL
                RRN=(RN-RSURF(1))/(RSURF(NLOCAL)-RSURF(1))
                FN=RN**2*ELPARM(RRN,ELLIN,ELLOT,EXELL)-CONST
                IF (ABS(FN/CONST).LT.EPS6) THEN
                  RSURF(J)=RN
                  GOTO 104
                ELSEIF (FN.LT.0.D0) THEN
                  RL=RN
                  GOTO 108
                ELSEIF (FN.GT.0.D0) THEN
                  RR=RN
                  GOTO 109
                ENDIF
  104         CONTINUE
            ELSEIF (NLTRI) THEN
              GESFL=((RAA*RAA-2.*TRIOT*TRIOT)*ELLOT-
     .               (RIA*RIA-2.*TRIIN*TRIIN)*ELLIN)*PIA
              FRING=GESFL/(NLOCAL-1)
              RSURF(1)=RIA
              RSURF(NLOCAL)=RAA
              WRITE (iunout,*) 'NLTRI NOT READY IN GRID'
              CALL EIRENE_EXIT_OWN(1)
            ENDIF
          ELSE
            WRITE (iunout,*) 'INVALID OPTION ENCOUNTERED IN SUBR. GRID'
            WRITE (iunout,*) 'INDGRD(1),LEVGEO ',INDGRD(IND),LEVGEO
            WRITE (iunout,*) 'EXIT CALLED FROM SUBR. GRID'
            CALL EIRENE_EXIT_OWN(1)
          ENDIF
        ELSEIF (INDGRD(IND).EQ.5) THEN
C** TAKE RADIAL GRID DATA FROM USER-SUPPLIED SUBROUTINE
          CALL EIRENE_PROUSR (RSURF,2+4*NPLS+3,0._DP,0._DP,0._DP,0._DP,
     .                 0._DP,0._DP,0._DP,NR1ST)
        ELSEIF (INDGRD(IND).EQ.6) THEN
C** TAKE RADIAL GRID DATA FROM INTERFACE
C         CALL PROFR (RSURF,2+4*NPLS+3,1,1,NR1ST)
        ENDIF
C
C  SET DATA FOR GRID OF ELLIPSES
C
        IF (INDGRD(IND).LE.2) THEN
          IF (NLCRC) THEN
            DO 111 J=1,NR1ST
              EP1(J)=0.
              ELL(J)=1.
              TRI(J)=0.
  111       CONTINUE
          ELSEIF (NLELL) THEN
            EP1(1)=EP1IN
            ELL(1)=ELLIN
            DO 112 J=2,NLOCAL
              RR=(RSURF(J)-RSURF(1))/(RSURF(NLOCAL)-RSURF(1))
              EP1(J)=ELPARM(RR,EP1IN,EP1OT,EXEP1)
              ELL(J)=ELPARM(RR,ELLIN,ELLOT,EXELL)
  112       CONTINUE
            DO 113 J=1,NR1ST
              TRI(J)=0.
  113       CONTINUE
          ELSEIF (NLTRI) THEN
            EP1(1)=EP1IN
            ELL(1)=ELLIN
            TRI(1)=TRIIN
            IF (ABS(TRI(1)/(RSURF(1)+EPS30)).GT.EPS30) THEN
              WRITE (iunout,*) 'FROM SUBR. GRID: '
              WRITE (iunout,*) 'ERROR IN TRIANGULARITY PARAMETERS '
              CALL EIRENE_EXIT_OWN(1)
            ENDIF
            DO 114 J=2,NLOCAL
              RR=(RSURF(J)-RSURF(1))/(RSURF(NLOCAL)-RSURF(1))
              EP1(J)=ELPARM(RR,EP1IN,EP1OT,EXEP1)
              ELL(J)=ELPARM(RR,ELLIN,ELLOT,EXELL)
              TRI(J)=ELPARM(RR,TRIIN,TRIOT,EXTRI)
  114       CONTINUE
          ENDIF
        ELSEIF (INDGRD(IND).EQ.5) THEN
C** TAKE ELLIP. GRID DATA FROM USER-SUPPLIED SUBROUTINE
C         CALL PROUSR (EP1,2+4*NPLS+?,0._DP,0._DP,0._DP,0._DP,0._DP,0._DP,0._DP,
C         CALL PROUSR (ELL,2+4*NPLS+?,0._DP,0._DP,0._DP,0._DP,0._DP,0._DP,0._DP,
C         CALL PROUSR (TRI,2+4*NPLS+?,0._DP,0._DP,0._DP,0._DP,0._DP,0._DP,0._DP,
        ELSEIF (INDGRD(IND).EQ.6) THEN
C** TAKE ELLIP. GRID DATA FROM INTERFACE
C         CALL PROFR (EP1,2+4*NPLS+?,1,1,NR1ST)
C         CALL PROFR (ELL,2+4*NPLS+?,1,1,NR1ST)
C         CALL PROFR (TRI,2+4*NPLS+?,1,1,NR1ST)
        ENDIF
C
C  SET DERIVED GRID DATA FOR LEVGEO = 2 OPTION
C  (SAME FOR ALL INDGRD OPTIONS)
C
        DO 115 J=1,NR1ST
          RQ(J)=RSURF(J)*RSURF(J)
          ELLQ(J)=ELL(J)*ELL(J)
  115   CONTINUE
C
        IF (TRCGRD) THEN
          CALL EIRENE_MASRR4('  N, RSURF,EP1,ELL,TRI',
     .                      RSURF,EP1,ELL,TRI,NR1ST)
          CALL EIRENE_LEER(2)
        ENDIF
C
      case (3)
C
C  GRID DATA GENERATION FOR LEVGEO.EQ.3
C
        IF (INDGRD(IND).LE.4) THEN
C  ALL POLYGON DATA HAVE BEEN READ FROM INPUT FILE, NOTHING ELSE
C  TO BE DONE HERE
        ELSEIF (INDGRD(IND).EQ.5) THEN
C*** POLYGON DATA NOT YET AVAILABLE FROM PROUSR (INDGRD.EQ.5 OPTION)
        ELSEIF (INDGRD(IND).EQ.6) THEN
C*** GEOMETRICAL DATA ARE SET IN IF0COP NOTHING TO BE DONE HERE
        ENDIF
C
        IF (PLREFL.GT.0._DP) THEN
C  SET ONE ADDITIONAL OUTERMOST POLYGON, DISTANCE PLREFL (CM) FROM THE
C  OUTERMOST POLYGON SPECIFIED SO FAR
          I=NR1STM
          DO 133 J=1,NPPLG
            DO 132 K=NPOINT(1,J),NPOINT(2,J)
              VPXX=XPOL(I,K)-XPOL(I-1,K)
              VPYY=YPOL(I,K)-YPOL(I-1,K)
              XNORM=SQRT(VPXX**2+VPYY**2)
              VPX=VPXX/XNORM*PLREFL
              VPY=VPYY/XNORM*PLREFL
              XPOL(I+1,K)=XPOL(I,K)+VPX
              YPOL(I+1,K)=YPOL(I,K)+VPY
  132       CONTINUE
  133     CONTINUE
        ENDIF
C
        IF (XPCOR.NE.0._DP) THEN
C  SHIFT WHOLE POLYGON MESH IN X DIRECTION BY XPCOR (CM)
C
          DO 135 I=1,NR1ST
            DO J=1,NPPLG
              DO K=NPOINT(1,J),NPOINT(2,J)
                XPOL(I,K)=XPOL(I,K)+XPCOR
              END DO
            END DO
  135     CONTINUE
        ENDIF
C
        IF (YPCOR.NE.0._DP) THEN
C  SHIFT WHOLE POLYGON MESH IN Y DIRECTION BY YPCOR (CM)
C
          DO 136 I=1,NR1ST
            DO J=1,NPPLG
              DO K=NPOINT(1,J),NPOINT(2,J)
                YPOL(I,K)=YPOL(I,K)+YPCOR
              END DO
            END DO
  136     CONTINUE
        ENDIF
C
C  SET DERIVED GRID DATA FOR LEVGEO = 3 OPTION
C  (SAME FOR ALL INDGRD OPTIONS)
C
        DO 140 I=1,NR1ST
          DO J=1,NPPLG
            DO K=NPOINT(1,J),NPOINT(2,J)-1
              VPLX(I,K)=XPOL(I,K+1)-XPOL(I,K)
              VPLY(I,K)=YPOL(I,K+1)-YPOL(I,K)
            END DO
          END DO
  140   CONTINUE
        DO 141 I=1,NR1ST
          DO J=1,NPPLG
            IF (J.EQ.1) THEN
              DO 142 K=1,NPOINT(1,1)
                BGL(I,K)=0.
  142         CONTINUE
            ELSE
              DO 143 K=NPOINT(2,J-1),NPOINT(1,J)
                BGL(I,K)=BGL(I,NPOINT(2,J-1))
  143         CONTINUE
            ENDIF
            DO K=NPOINT(1,J)+1,NPOINT(2,J)
              BGL(I,K)=BGL(I,K-1)+SQRT(VPLX(I,K-1)**2+VPLY(I,K-1)**2)
            END DO
          END DO
  141   CONTINUE
C
C   CALCULATE THE OUTER NORMALS OF POLYGONS
C
        DO 144 I=1,NR1ST
          DO J=1,NPPLG
            DO K=NPOINT(1,J),NPOINT(2,J)-1
              PLABS=SQRT(VPLX(I,K)**2+VPLY(I,K)**2)
              PLNX(I,K)=VPLY(I,K)/(PLABS+EPS60)
              PLNY(I,K)=-VPLX(I,K)/(PLABS+EPS60)
            END DO
          END DO
  144   CONTINUE
C
        DO 147 I=1,NR1ST
          IUP=I+1
          IDN=I
          IF (IUP.GT.NR1ST) THEN
            IUP=I
            IDN=I-1
          ENDIF
          DO J=1,NPPLG
           DO K=NPOINT(1,J),NPOINT(2,J)-1
  146       XD=XPOL(IUP,K+1)-XPOL(IDN,K+1)
            YD=YPOL(IUP,K+1)-YPOL(IDN,K+1)
            IF (XD*XD+YD*YD.LT.EPS30) THEN
              IF (IUP.LT.NR1ST) THEN
                IUP=IUP+1
              ELSE
                IDN=IDN-1
              ENDIF
              GOTO 146
            ENDIF
            XS=SIGN(1._DP,XD*PLNX(I,K)+YD*PLNY(I,K))
            PLNX(I,K)=PLNX(I,K)*XS
            PLNY(I,K)=PLNY(I,K)*XS
           END DO
          END DO
  147   CONTINUE
C
        IF (TRCGRD) THEN
          WRITE (iunout,*) ' NO. OF VALID PARTS = ',NPPLG
          DO 155 J=1,NR1ST
            WRITE (iunout,*) ' POLYGON NO. J = ',J
            DO 156 K=1,NPPLG
              WRITE (iunout,*) 'IA = ',NPOINT(1,K),' IE = ',NPOINT(2,K)
              WRITE (iunout,'(/1X,1P,6E12.4)') (XPOL(J,I),YPOL(J,I),
     .                                   I=NPOINT(1,K),NPOINT(2,K))
  156       CONTINUE
  155     CONTINUE
          CALL EIRENE_LEER(2)
          WRITE (iunout,*)
     .      'ARCLENGTH BGL(I,K) OF RADIAL SURFACES AT Z=0.'
          DO 153 I=1,NR1ST
            WRITE (iunout,*) 'I = ',I
            WRITE (iunout,'(/1X,1P,6E12.4)') (BGL(I,K),K=1,NRPLG)
            CALL EIRENE_LEER(1)
  153     CONTINUE
        ENDIF
C
        CALL EIRENE_SNEIGH
C
      CASE (4) ! LEVGEO
C
C  GRID DATA GENERATION FOR LEVGEO.EQ.4
C
C
C  SET DERIVED GRID DATA FOR LEVGEO = 4 OPTION
C  (SAME FOR ALL INDGRD OPTIONS)
C

cdr  ncorner: number of cell vertices in case levgeo=4, for FEM interpolation
        NCORNER = NKNOT

!pb initialize list of triangles per gridpoint
        ALLOCATE (COORTRI(NKNOT))
        DO I=1,NKNOT
          NULLIFY(COORTRI(I)%PTRI)
        END DO

        DO 165 I=1,NTRII
          VTRIX(1,I)=XTRIAN(NECKE(2,I))-XTRIAN(NECKE(1,I))
          VTRIY(1,I)=YTRIAN(NECKE(2,I))-YTRIAN(NECKE(1,I))
          VTRIX(2,I)=XTRIAN(NECKE(3,I))-XTRIAN(NECKE(2,I))
          VTRIY(2,I)=YTRIAN(NECKE(3,I))-YTRIAN(NECKE(2,I))
          VTRIX(3,I)=XTRIAN(NECKE(1,I))-XTRIAN(NECKE(3,I))
          VTRIY(3,I)=YTRIAN(NECKE(1,I))-YTRIAN(NECKE(3,I))

          ALLOCATE (CUR)
          CUR%NOTRI = I
          CUR%NEXT_TRI => COORTRI(NECKE(1,I))%PTRI
          COORTRI(NECKE(1,I))%PTRI => CUR

          ALLOCATE (CUR)
          CUR%NOTRI = I
          CUR%NEXT_TRI => COORTRI(NECKE(2,I))%PTRI
          COORTRI(NECKE(2,I))%PTRI => CUR

          ALLOCATE (CUR)
          CUR%NOTRI = I
          CUR%NEXT_TRI => COORTRI(NECKE(3,I))%PTRI
          COORTRI(NECKE(3,I))%PTRI => CUR
  165   CONTINUE
C
C
C   CALCULATE THE OUTER NORMALS OF TRIANGLES
C
        DO 161 I=1,NTRII
          PLABS1=SQRT(VTRIX(1,I)**2+VTRIY(1,I)**2)
          PLABS2=SQRT(VTRIX(2,I)**2+VTRIY(2,I)**2)
          PLABS3=SQRT(VTRIX(3,I)**2+VTRIY(3,I)**2)
          PTRIX(1,I)=VTRIY(1,I)/(PLABS1+EPS60)
          PTRIX(2,I)=VTRIY(2,I)/(PLABS2+EPS60)
          PTRIX(3,I)=VTRIY(3,I)/(PLABS3+EPS60)
          PTRIY(1,I)=-VTRIX(1,I)/(PLABS1+EPS60)
          PTRIY(2,I)=-VTRIX(2,I)/(PLABS2+EPS60)
          PTRIY(3,I)=-VTRIX(3,I)/(PLABS3+EPS60)
  161   CONTINUE
C
C PTRIX/PTRIY POINT OUT OF TRIANGE FOR A MATHEMATICAL POSITVE
C             ORIENTATION OF TRIANGLE (1-2-3-1: COUNTER-CLOCKWISE,
C                                               AS IT MUST BE)
C PTRIX/PTRIY POINT INTO TRIANGE FOR A MATHEMATICAL NEGATIVE
C             ORIENTATION OF TRIANGLE (1-2-3-1: CLOCKWISE,
C                                               ERROR MESSAGE FROM LEARC1)
C
        DO 162 I=1,NTRII
          XD1=VTRIX(1,I)+PTRIX(1,I)
          YD1=VTRIY(1,I)+PTRIY(1,I)
          XS1=SIGN(1._DP,XD1*VTRIX(1,I)+YD1*VTRIY(1,I))
          PTRIX(1,I)=PTRIX(1,I)*XS1
          PTRIY(1,I)=PTRIY(1,I)*XS1
          XD2=VTRIX(2,I)+PTRIX(2,I)
          YD2=VTRIY(2,I)+PTRIY(2,I)
          XS2=SIGN(1._DP,XD2*VTRIX(2,I)+YD2*VTRIY(2,I))
          PTRIX(2,I)=PTRIX(2,I)*XS2
          PTRIY(2,I)=PTRIY(2,I)*XS2
          XD3=VTRIX(3,I)+PTRIX(3,I)
          YD3=VTRIY(3,I)+PTRIY(3,I)
          XS3=SIGN(1._DP,XD3*VTRIX(3,I)+YD3*VTRIY(3,I))
          PTRIX(3,I)=PTRIX(3,I)*XS3
          PTRIY(3,I)=PTRIY(3,I)*XS3
  162   CONTINUE
C
C  BUILD LIST OF TRIANGLE SIDES COMPRISING A NON-DEFAULT STANDARD SURFACE
C                               OR AN ADDITIONAL SURFACE
C
        IF (MAXVAL(INMTI(1:3,1:NTRII)) > NSTSI+NLIM) THEN
          WRITE (iunout,*)
     .      ' WRONG INDEX OF REFLECTION MODEL SPECIFIED '
          WRITE (iunout,*) ' CHECK DEFINITION OF TRIANGLES ',
     .                     ' AND THEIR MODEL-FLAGS FOR SIDES '
          CALL EIRENE_EXIT_OWN(1)
        END IF

C  INDEX J:            EIRENE SURFACE  (NON-DEFAULT, OR ADDITIONAL)
C  INDEX IT, OR ITRI:  TRIANGLE
C  INDEX IS:           TRIANGLE SIDE

        DO J = 1, NLIMPS
          NSRFTR = COUNT(INMTI(1:3,1:NTRII) .EQ. J)
          IF (NSRFTR > 0) THEN
            ALLOCATE (SURF_TRIAN(J)%ITRIAS(NSRFTR))
            ALLOCATE (SURF_TRIAN(J)%ITRISI(NSRFTR))
            ALLOCATE (SURF_TRIAN(J)%BGLT(NSRFTR+1))
            SURF_TRIAN(J)%BGLT(1) = 0._DP
          END IF
        END DO

        LERROR = .FALSE.
        DO IT=1,NTRII
          DO IS = 1, 3
            J = INMTI(IS,IT)
            IF ( J .NE. 0) THEN
              IF ((NCHBAR(IS,IT) == 0) .AND. (ILIIN(J) <= 0)) THEN
                IF ((ILIIN(J) < 0) .AND. (ILSWCH(J) < 1000)) THEN
c  transparent surfaces, that switch into additional cells, are legal.
c  all other transparent cell faces must either have a neighbor, or a surface boundary condition.
                WRITE (iunout,*) 'SIDE',IS,' OF TRIANGLE ',IT,
     .              ' IS TRANSPARENT BUT HAS NO NEIGHBOR '
                LERROR = .TRUE.
                END IF
              END IF
              SURF_TRIAN(J)%NUMTR = SURF_TRIAN(J)%NUMTR + 1
              SURF_TRIAN(J)%ITRIAS(SURF_TRIAN(J)%NUMTR) = IT
              SURF_TRIAN(J)%ITRISI(SURF_TRIAN(J)%NUMTR) = IS
              SURF_TRIAN(J)%BGLT(SURF_TRIAN(J)%NUMTR+1) =
     .          SURF_TRIAN(J)%BGLT(SURF_TRIAN(J)%NUMTR) +
     .          SQRT(VTRIX(IS,IT)**2+VTRIY(IS,IT)**2)
            END IF
          END DO
        END DO

        IF (LERROR) THEN
          WRITE (iunout,*) ' CALCULATION STOPPED DUE TO ERRORS',
     .                     ' LISTED ABOVE '
          CALL EIRENE_EXIT_OWN(1)
        END IF

! ADDITIONAL SURFACES
        DO J=1, NLIMI
          IF (ILPLG(J) /= 0) THEN
            NT = SURF_TRIAN(J)%NUMTR
!pb31072012            IF (NT == 1) CYCLE
            IF (NT <= 1) CYCLE
            DP1 = HUGE(1._DP)
! FIND TRIANGLE STARTING AT FIRST POINT OF SURFACE J
            DO I = 1, NT
              IT = SURF_TRIAN(J)%ITRIAS(I)
              IS = SURF_TRIAN(J)%ITRISI(I)
              IS1 = IS + 1
              IF (IS1 > 3) IS1 = 1
              DS1 = SQRT((P1(1,J)-XTRIAN(NECKE(IS,IT)))**2 +
     .                   (P1(2,J)-YTRIAN(NECKE(IS,IT)))**2)
              DS2 = SQRT((P1(1,J)-XTRIAN(NECKE(IS1,IT)))**2 +
     .                   (P1(2,J)-YTRIAN(NECKE(IS1,IT)))**2)
              IF (DS1 < DP1) THEN
                DP1 = DS1
                IM = I
                IMP = NECKE(IS1,IT)
              END IF
              IF (DS2 < DP1) THEN
                DP1 = DS2
                IM = I
                IMP = NECKE(IS,IT)
              END IF
            END DO
            IF ((DP1/SURF_TRIAN(J)%BGLT(NT+1) > 1.D-2)
     .          .OR. (IM*IMP == 0)) THEN
               WRITE (IUNOUT,*) ' PROBLEM FINDING STARTING POINT',
     .                          ' FOR SORTING OF TRIANGLES '
               WRITE (IUNOUT,*)
     .            ' SORTING ABANDONED FOR SURFACE NUMBER ',I
               CYCLE
            END IF

! SET FIRST TRIANGLE ALONG SURFACE
            IF (IM /= 1) THEN
              IT = SURF_TRIAN(J)%ITRIAS(1)
              IS = SURF_TRIAN(J)%ITRISI(1)
              SURF_TRIAN(J)%ITRIAS(1) = SURF_TRIAN(J)%ITRIAS(IM)
              SURF_TRIAN(J)%ITRISI(1) = SURF_TRIAN(J)%ITRISI(IM)
              SURF_TRIAN(J)%ITRIAS(IM) = IT
              SURF_TRIAN(J)%ITRISI(IM) = IS
              IT = SURF_TRIAN(J)%ITRIAS(1)
              IS = SURF_TRIAN(J)%ITRISI(1)
              SURF_TRIAN(J)%BGLT(2) =
     .             SURF_TRIAN(J)%BGLT(1) +
     .             SQRT(VTRIX(IS,IT)**2+VTRIY(IS,IT)**2)
            END IF

! SORT THE TRIANGLES
            DO I = 2, NT-1
              DO K = I, NT
                IT = SURF_TRIAN(J)%ITRIAS(K)
                IS = SURF_TRIAN(J)%ITRISI(K)
                IS1 = IS + 1
                IF (IS1 > 3) IS1 = 1
                IF (IMP == NECKE(IS,IT)) THEN
                  IMP = NECKE(IS1,IT)
                  EXIT
                END IF
                IF (IMP == NECKE(IS1,IT)) THEN
                  IMP = NECKE(IS,IT)
                  EXIT
                END IF
              END DO
              IF (K > NT) THEN
                WRITE (IUNOUT,*) ' NO MATCHING TRIANGLE FOUND '
              ELSE
                IF (I /= K) THEN
                   IT = SURF_TRIAN(J)%ITRIAS(I)
                   IS = SURF_TRIAN(J)%ITRISI(I)
                   SURF_TRIAN(J)%ITRIAS(I) = SURF_TRIAN(J)%ITRIAS(K)
                   SURF_TRIAN(J)%ITRISI(I) = SURF_TRIAN(J)%ITRISI(K)
                   SURF_TRIAN(J)%ITRIAS(K) = IT
                   SURF_TRIAN(J)%ITRISI(K) = IS
                END IF
                IT = SURF_TRIAN(J)%ITRIAS(I)
                IS = SURF_TRIAN(J)%ITRISI(I)
                SURF_TRIAN(J)%BGLT(I+1) =
     .               SURF_TRIAN(J)%BGLT(I) +
     .               SQRT(VTRIX(IS,IT)**2+VTRIY(IS,IT)**2)
              END IF
            END DO
          END IF
        END DO
C
        IF (TRCGRD) THEN
          WRITE (iunout,*) ' NUMBER OF TRIANGLES = ',NTRII
          WRITE (iunout,*) ' I,(XTRIAN(J),YTRIAN(J),J=1,3) '
          CALL EIRENE_LEER(1)
          DO 163 I=1,NTRII
            WRITE (iunout,'(1X,I4,1X,1P,6E12.4)')
     .                               I,(XTRIAN(NECKE(J,I)),
     .                                  YTRIAN(NECKE(J,I)),J=1,3)
  163     CONTINUE
          CALL EIRENE_LEER(2)
          WRITE (iunout,*) ' NGITT SET TO ',NGITT

          DO J=1, NLIMPS
            IF (J.LE.NLIM) THEN
              JLIM = J
            ELSE
              JLIM = -(J-NLIM)
            ENDIF
            WRITE (IUNOUT,*)
            WRITE (IUNOUT,*) ' SURFACE NO. ',JLIM
            WRITE (IUNOUT,'(5A6,A12)') 'J','ITRI','ISIDE',
     .                                 'IP1','IP2','BLGT'
            DO I=1, SURF_TRIAN(J)%NUMTR
              IT = SURF_TRIAN(J)%ITRIAS(I)
              IS = SURF_TRIAN(J)%ITRISI(I)
              IS1 = IS + 1
              IF (IS1 > 3) IS1 = 1
              WRITE (IUNOUT,'(5I6,ES12.4)') I, IT, IS, NECKE(IS,IT),
     .                                         NECKE(IS1,IT),
     .                                         SURF_TRIAN(J)%BGLT(I+1)
            END DO
          END DO
        ENDIF
C
C
      CASE (5) ! LEVGEO
C
C  GRID DATA GENERATION FOR LEVGEO.EQ.5
C
C  GRID DATA FOR TETRAHEDRA ARE SET IN COUPLING ROUTINE
C  NOTHING TO BE DONE HERE
C
C  SET DERIVED GRID DATA FOR LEVGEO = 5 OPTION
C  (SAME FOR ALL INDGRD OPTIONS)
C
cdr  ncorner: number of cell vertices in case levgeo=5, for FEM interpolation
        NCORNER = NCOORD

        DO ITET=1,NTET
          IC1 = NTECK(1,ITET)
          IC2 = NTECK(2,ITET)
          IC3 = NTECK(3,ITET)
          IC4 = NTECK(4,ITET)
          RINCRC(1:4,ITET) = 1._DP
C  CALCULATE DIRECTIONS OF EDGES
C  EDGE  1-2
          VTETX(1,ITET) = XTETRA(IC2) - XTETRA(IC1)
          VTETY(1,ITET) = YTETRA(IC2) - YTETRA(IC1)
          VTETZ(1,ITET) = ZTETRA(IC2) - ZTETRA(IC1)
C  EDGE  2-3
          VTETX(2,ITET) = XTETRA(IC3) - XTETRA(IC2)
          VTETY(2,ITET) = YTETRA(IC3) - YTETRA(IC2)
          VTETZ(2,ITET) = ZTETRA(IC3) - ZTETRA(IC2)
C  EDGE  3-1
          VTETX(3,ITET) = XTETRA(IC1) - XTETRA(IC3)
          VTETY(3,ITET) = YTETRA(IC1) - YTETRA(IC3)
          VTETZ(3,ITET) = ZTETRA(IC1) - ZTETRA(IC3)
C  EDGE  1-4
          VTETX(4,ITET) = XTETRA(IC4) - XTETRA(IC1)
          VTETY(4,ITET) = YTETRA(IC4) - YTETRA(IC1)
          VTETZ(4,ITET) = ZTETRA(IC4) - ZTETRA(IC1)
C  EDGE  2-4
          VTETX(5,ITET) = XTETRA(IC4) - XTETRA(IC2)
          VTETY(5,ITET) = YTETRA(IC4) - YTETRA(IC2)
          VTETZ(5,ITET) = ZTETRA(IC4) - ZTETRA(IC2)
C  EDGE  3-4
          VTETX(6,ITET) = XTETRA(IC4) - XTETRA(IC3)
          VTETY(6,ITET) = YTETRA(IC4) - YTETRA(IC3)
          VTETZ(6,ITET) = ZTETRA(IC4) - ZTETRA(IC3)

          EDGELEN(1:6) = SQRT(VTETX(1:6,ITET)**2 +
     .                        VTETY(1:6,ITET)**2 +
     .                        VTETZ(1:6,ITET)**2)
C  CALCULATE THE OUTER NORMALS OF TETRAHEDRA
C  SIDE 1-2-3
          PTETX(1,ITET) = VTETY(3,ITET)*VTETZ(1,ITET) -
     .                    VTETZ(3,ITET)*VTETY(1,ITET)
          PTETY(1,ITET) = VTETZ(3,ITET)*VTETX(1,ITET) -
     .                    VTETX(3,ITET)*VTETZ(1,ITET)
          PTETZ(1,ITET) = VTETX(3,ITET)*VTETY(1,ITET) -
     .                    VTETY(3,ITET)*VTETX(1,ITET)
          S = 0.5_DP *( EDGELEN(1) + EDGELEN(2) + EDGELEN(3) )
          SQ = SQRT( (S-EDGELEN(1)) *
     .               (S-EDGELEN(2)) * (S-EDGELEN(3)) / S)
          IF (SQ > EPS30) RINCRC(1,ITET) = 1._DP / SQ
C  SIDE 1-2-4
          PTETX(2,ITET) = VTETY(4,ITET)*VTETZ(1,ITET) -
     .                    VTETZ(4,ITET)*VTETY(1,ITET)
          PTETY(2,ITET) = VTETZ(4,ITET)*VTETX(1,ITET) -
     .                    VTETX(4,ITET)*VTETZ(1,ITET)
          PTETZ(2,ITET) = VTETX(4,ITET)*VTETY(1,ITET) -
     .                    VTETY(4,ITET)*VTETX(1,ITET)
          S =  0.5_DP *( EDGELEN(1) + EDGELEN(5) + EDGELEN(4) )
          SQ = SQRT( (S-EDGELEN(1)) *
     .               (S-EDGELEN(5)) * (S-EDGELEN(4)) / S)
          IF (SQ > EPS30) RINCRC(2,ITET) = 1._DP / SQ
C  SIDE 2-3-4
          PTETX(3,ITET) = VTETY(5,ITET)*VTETZ(2,ITET) -
     .                    VTETZ(5,ITET)*VTETY(2,ITET)
          PTETY(3,ITET) = VTETZ(5,ITET)*VTETX(2,ITET) -
     .                    VTETX(5,ITET)*VTETZ(2,ITET)
          PTETZ(3,ITET) = VTETX(5,ITET)*VTETY(2,ITET) -
     .                    VTETY(5,ITET)*VTETX(2,ITET)
          S =  0.5_DP *( EDGELEN(2) + EDGELEN(6) + EDGELEN(5) )
          SQ = SQRT( (S-EDGELEN(2)) *
     .               (S-EDGELEN(6)) * (S-EDGELEN(5)) / S)
          IF (SQ > EPS30) RINCRC(3,ITET) = 1._DP / SQ
C  SIDE 3-1-4
          PTETX(4,ITET) = VTETY(4,ITET)*VTETZ(3,ITET) -
     .                    VTETZ(4,ITET)*VTETY(3,ITET)
          PTETY(4,ITET) = VTETZ(4,ITET)*VTETX(3,ITET) -
     .                    VTETX(4,ITET)*VTETZ(3,ITET)
          PTETZ(4,ITET) = VTETX(4,ITET)*VTETY(3,ITET) -
     .                    VTETY(4,ITET)*VTETX(3,ITET)
          S =  0.5_DP *( EDGELEN(3) + EDGELEN(4) + EDGELEN(6) )
          SQ = SQRT( (S-EDGELEN(3)) *
     .               (S-EDGELEN(4)) * (S-EDGELEN(6)) / S)
          IF (SQ > EPS30) RINCRC(4,ITET) = 1._DP / SQ

          DO J=1,4
            PLEN=SQRT(PTETX(J,ITET)**2+PTETY(J,ITET)**2+
     .                PTETZ(J,ITET)**2)+EPS60
            PTETX(J,ITET)=PTETX(J,ITET)/PLEN
            PTETY(J,ITET)=PTETY(J,ITET)/PLEN
            PTETZ(J,ITET)=PTETZ(J,ITET)/PLEN
          END DO
        END DO


        CALL EIRENE_SUCHE_NACHBARN

        IC=0
        NTET_COLLAPS=0
        LERROR = .FALSE.
        DO ITET=1,NTET
          DO IS=1,4
            IF ((NTBAR(IS,ITET) == 0) .AND. (INMTIT(IS,ITET) == 0)) THEN
              IC=IC+1
              WRITE (iunout,*) ' TETRAHEDRON WITH NO NEIGHBORS AND NO ',
     .                    'REFLECTION MODEL FOUND'
              WRITE (iunout,*) ' ITET = ',ITET,' ISIDE = ',IS
              write (iunout,*) nteck(itside(1,is),itet),
     .                    nteck(itside(2,is),itet),
     .                    nteck(itside(3,is),itet)
            END IF
            IF (NTBAR(IS,ITET) < 0) THEN
              IF (SUM(NTBAR(1:4,ITET)) > -4) THEN
                WRITE (iunout,*)
     .            ' TETRAHEDRON WITH NEIGHBOR -1 DETECTED'
                WRITE (iunout,*) ' ITET = ',ITET,' ISIDE = ',IS
                WRITE (iunout,*) ' NTBAR(ITET) = ',NTBAR(1:4,ITET)
                IC=IC+1
              ELSE
                NTET_COLLAPS=NTET_COLLAPS+1
!pb             WRITE(iunout,*) ' COLLAPSED TETRAHEDRON ITET = ',ITET
                EXIT
              END IF
            END IF
            J = INMTIT(IS,ITET)
            IF ((NTBAR(IS,ITET) == 0) .AND. (ILIIN(J) <= 0)) THEN
              IF ((ILIIN(J) < 0) .AND. (ILSWCH(J) < 1000)) THEN
c  transparent surfaces, that switch into additional cells, are legal.
c  all other transparent cell faces must either have a neighbor, or a surface boundary condition.
                WRITE (iunout,*) 'SIDE',IS,' OF TETRAHEDRON ',ITET,
     .              ' IS TRANSPARENT BUT HAS NO NEIGHBOR'
                LERROR = .TRUE.
              END IF
            END IF
          END DO
        END DO

        IF ((IC > 0) .OR. LERROR) THEN
          WRITE (iunout,*) ' CALCULATION STOPPED DUE TO ERRORS',
     .                     ' LISTED ABOVE '
          CALL EIRENE_EXIT_OWN(1)
        END IF
C
        IF (TRCGRD) THEN
          WRITE (iunout,*) ' NUMBER OF COORDINATES = ',NCOOR
          WRITE (iunout,*) ' I,(XETRA(J),YTETRA(J),ZTETRA(J),J=1,3) '
          DO I=1,NCOOR
            WRITE (iunout,'(1X,I4,1X,6ES12.4)')
     .             I,XTETRA(I),YTETRA(I),ZTETRA(I)
          END DO
          CALL EIRENE_LEER(2)

          WRITE (iunout,*) ' NUMBER OF TETRAHEDRA = ',NTET
          DO ITET=1,NTET
            WRITE (iunout,*)
            WRITE (iunout,*) ' TETRAHEDRON ',ITET
            DO J=1,4
              IC=NTECK(J,ITET)
              WRITE (iunout,'(1X,I6,3ES12.4)')
     .               IC,XTETRA(IC),YTETRA(IC),ZTETRA(IC)
            END DO
            DO J=1,4
              WRITE (iunout,'(A,I3,A,3I6,4X,A,I6,A,I6)')
     .              ' SIDE ',J,': ',
     .                NTECK(ITSIDE(1,J),ITET),
     .                NTECK(ITSIDE(2,J),ITET),
     .                NTECK(ITSIDE(3,J),ITET),
     .              ' NEIGHBOR ',NTBAR(J,ITET),
     .              ' SIDE ',NTSEITE(J,ITET)
            END DO
            WRITE (iunout,*) 'OUTER NORMALS '
            DO J=1,4
              WRITE (iunout,'(7X,3ES12.4)')
     .              PTETX(J,ITET),PTETY(J,ITET),PTETZ(J,ITET)
            END DO
          END DO

C  CHECK OUTER NORMALS
          DO ITET = 1,NTET
            DO J=1,4
             IF (NTBAR(J,ITET) /= 0) THEN
               PC1(1:3) = (/
     .         ABS(PTETX(J,ITET)+PTETX(NTSEITE(J,ITET),NTBAR(J,ITET))),
     .         ABS(PTETY(J,ITET)+PTETY(NTSEITE(J,ITET),NTBAR(J,ITET))),
     .         ABS(PTETZ(J,ITET)+PTETZ(NTSEITE(J,ITET),NTBAR(J,ITET)))/)
               IF (ANY(PC1 > 1.E-6)) THEN
                WRITE(iunout,*) ' PROBLEM WITH TETRAHEDRON ',ITET,
     .                          ' SIDE ',J
                WRITE(iunout,*) ' NORMALS DO NOT MATCH '
                WRITE(iunout,*) PTETX(J,ITET),PTETY(J,ITET),
     .                          PTETZ(J,ITET)
                WRITE(iunout,*) PTETX(NTSEITE(J,ITET),NTBAR(J,ITET)),
     .                     PTETY(NTSEITE(J,ITET),NTBAR(J,ITET)),
     .                     PTETZ(NTSEITE(J,ITET),NTBAR(J,ITET))
               END IF
             END IF
            END DO
          END DO

        ENDIF
C
      CASE (10) ! LEVGEO
C
C  GENERAL GEOMETRY OPTION: NOTHING TO DONE HERE
        NCORNER=0
C
      END SELECT
C
C  SET GEOMETRICAL CONSTANTS FOR IGNORABLE Y OR POLOIDAL COORDINATE
C  THESE MAY BE REVISED IF A 2ND (Y- OR POL.) GRID IS DEFINED BELOW
C
      IF (LEVGEO.EQ.1) THEN
        YDF=YAA-YIA
        PSURF(1)=YIA
      ELSEIF (LEVGEO.EQ.2.AND..NOT.NLPOL) THEN
        YDF=1.
        PSURF(1)=YIA
      ELSEIF (LEVGEO.EQ.2.AND.NLPOL) THEN
        YDF=(YAA-YIA)*DEGRAD
        PSURF(1)=YIA
      ELSEIF (LEVGEO.EQ.3) THEN
        YDF=1.
        PSURF(1)=YIA
      ELSEIF (LEVGEO.EQ.4) THEN
        YDF=1.
        PSURF(1)=YIA
      ELSEIF (LEVGEO.EQ.5) THEN
        YDF=1.
      ELSEIF (LEVGEO.EQ.10) THEN
        YDF=1.
C
C  GENERAL GEOMETRY OPTION: NOTHING TO BE DONE HERE
C
      ENDIF
      IF (YDF.LE.0._DP) GOTO 992
C
      IF (TRCGRD.AND..NOT.NLPOL) THEN
        CALL EIRENE_LEER(2)
        WRITE (iunout,*) 'CONSTANTS FOR POLOIDAL OR Y DIRECTION'
        CALL EIRENE_MASR3('YDF,YIA,YAA=            ',YDF,YIA,YAA)
        CALL EIRENE_LEER(1)
      ENDIF
C
C  SET GEOMETRICAL CONSTANTS FOR IGNORABLE Z COORDINATE
C  THESE MAY BE REVISED IF A 3RD (Z- OR TOR.) GRID IS DEFINED BELOW
C
C  A) IN TOROIDAL APPROXIMATION:
C
      IF (NLTRA) THEN

        IF (NTTRA.LE.3.OR.ROA.LT.0._DP) GOTO 991
        select case (LEVGEO)
        case (1,3:5)
          XDIFF=0.
        case (2)
          XDIFF=EP1OT
C
C  GENERAL GEOMETRY OPTION: NOTHING TO BE DONE HERE
C
        end select
C
C  TOROIDAL ANGLE, INPUT IS IN DEGREES, CONVERT TO RADIANS
        ZDF=(ZAA-ZIA)*DEGRAD
C
C  ALPHA: HALF OF THE ANGLE INCREMENT IN EQUIDISTANT TOROIDAL ANGLE GRID
        ALPHA=0.5*(ZDF/DBLE(NTTRAM))
        TANAL=TAN(ALPHA)
        SINAL=SIN(2.*ALPHA)
        COSAL=COS(2.*ALPHA)
C
        DPHI=1./(2.*ALPHA)
C
C  ROA IS THE LARGE RADIUS OF THE TORUS
C  (RMTOR,0,0) IS THE ORIGIN OF LOCAL COORDINATE SYSTEM IN
C              EACH TOROIDAL CELL
C  RMTOR SUCH THAT VOLUME OF TORUS = VOLUME OF THE NTTRAM SEGMENTS
C  AT PRESENT: FULLFILLED AT SURFACE DEFINED BY (RAA,EP1OT,ELLOT)
C              OR AT A POLYGON WITH XDIFF=0. (IF THERE IS ONE)
        RMTOR=(ROA+XDIFF)*ALPHA/TANAL-XDIFF
        ZHALF=ALPHA
        ZFULL=ZHALF*2.
C  SET ZSURF EVEN IF NLTOR=FALSE, FOR 3D GEOMETRY PLOTS
        DO 170 J=1,NTTRA
          ZSURF(J)=ZIA*DEGRAD+(J-1)/DPHI
  170   CONTINUE
        DO 172 J=1,NTTRAM
          ZZONE(J)=0.5*(ZSURF(J)+ZSURF(J+1))
  172   CONTINUE
        RORIG=RMTOR
C
C  B) IN CYLIND. APPROXIMATION:
C     ROA AND RMTOR ARE IRRELEVANT IN THIS CASE, AND ARE NOT DEFINED
C
      ELSEIF (NLTRZ) THEN

        ZDF=ZAA-ZIA
        IF (ZDF.LE.0._DP) GOTO 991
        ZSURF(1)=ZIA
        ZZONE(1)=(ZAA+ZIA)*0.5
        DPHI=1./ZDF
C
        RORIG=0.
C
C  C) IN TORUS COORDINATES
C
      ELSEIF (NLTRT) THEN
        ZDF=(ZAA-ZIA)*DEGRAD
        IF (ZDF.LE.0._DP) GOTO 991
        ZSURF(1)=ZIA
        ZZONE(1)=(ZAA+ZIA)*0.5
        DPHI=1./ZDF
C
        RORIG=0.
C
      ELSE
        WRITE (iunout,*) ' ERROR IN INPUT DATA! '
        WRITE (iunout,*) ' NLTRA OR NLTRZ OR NLTRT MUST BE .TRUE. '
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
C
      IF (TRCGRD) THEN
        CALL EIRENE_LEER(2)
        IF (.NOT.NLTOR) THEN
          WRITE (iunout,*) 'CONSTANTS FOR TOROIDAL OR Z DIRECTION'
          CALL EIRENE_MASR3('ZDF,ZIA,ZAA=            ',ZDF,ZIA,ZAA)
        ENDIF
        IF (NLTRA) THEN
          CALL EIRENE_MASR2('ROA,RMTOR=      ',ROA,RMTOR)
          IF (.NOT.NLTOR) THEN
            CALL EIRENE_MASRR1 (' N,  ZSURF ',ZSURF,NTTRA,3)
            CALL EIRENE_MASRR1 (' N,  ZZONE ',ZZONE,NTTRAM,3)
          ENDIF
          CALL EIRENE_LEER(2)
        ENDIF
        CALL EIRENE_LEER(1)
      ENDIF
C
C  SET SURFACE AREA OF NON-DEFAULT STANDARD SURFACES
C
      SELECT CASE (LEVGEO)
      CASE (1) ! LEVGEO
        DO 180 ISTS=1,NSTSI
          IF (INUMP(ISTS,1).NE.0) THEN
            IR=INUMP(ISTS,1)
            NLJ=NLIM+ISTS
            IF (NLTRZ) THEN
              SAREA(NLJ)=YDF*ZDF
            ELSEIF (NLTRA) THEN
              SAREA(NLJ)=YDF*ZDF*(RSURF(IR)+RMTOR)
            ELSEIF (NLTRT) THEN
            ENDIF
          ENDIF
  180   CONTINUE

      CASE (4) ! LEVGEO
        DO ISTS=1,NSTSI
          NLJ=NLIM+ISTS
          SAREA(NLJ)=0.
        ENDDO
!        DO I=1,NTRII
!          DO J=1,3
!            IECKE2 = J+1
!            IF (IECKE2 .EQ. 4) IECKE2 = 1
!            ISTS=ABS(INMTI(J,I))
!            IF (ISTS .GT. NLIM) THEN
C  SEITE J VON DREIECK I GEHOERT ZUM RAND ISTS
!              XX1=XTRIAN(NECKE(J,I))
!              YY1=YTRIAN(NECKE(J,I))
!              XX2=XTRIAN(NECKE(IECKE2,I))
!              YY2=YTRIAN(NECKE(IECKE2,I))
!              DSD=((XX1-XX2)**2+(YY1-YY2)**2)**0.5
!              IF (NLTRA) THEN
!                XX1=XX1+RMTOR
!                XX2=XX2+RMTOR
!                COM=0.5*(XX1+XX2)
!                DSD=DSD*COM*2.*PIA
!              ELSE
!                DSD=DSD*ZDF
!              ENDIF
!              IF (NLMPGS.GT.NLIMPS) THEN
!                MSURFG=NLIM+NSTS+INSPAT(J,I)
!                SAREA(MSURFG)=DSD
!              END IF
!              SAREA(ISTS)=SAREA(ISTS)+DSD
!            ENDIF
!          ENDDO
!        ENDDO

        ALLOCATE (VISITED(0:3,0:NTRI))
        VISITED = .FALSE.

        DO ISTS = 1, NLIMPS

cdr       write (iunout,*) ' area for surface ',ists
          DO J= 1, SURF_TRIAN(ISTS)%NUMTR
            I = SURF_TRIAN(ISTS)%ITRIAS(J)
            IS = SURF_TRIAN(ISTS)%ITRISI(J)
            IF (VISITED(IS,I)) CYCLE
            IF (NCHBAR(IS,I) > 0) THEN
              IF (VISITED(NSEITE(IS,I),NCHBAR(IS,I))) CYCLE
            END IF
            IECKE2 = IS+1
            IF (IECKE2 .EQ. 4) IECKE2 = 1
            XX1=XTRIAN(NECKE(IS,I))
            YY1=YTRIAN(NECKE(IS,I))
            XX2=XTRIAN(NECKE(IECKE2,I))
            YY2=YTRIAN(NECKE(IECKE2,I))
            DSD=((XX1-XX2)**2+(YY1-YY2)**2)**0.5
            VISITED(IS,I) = .TRUE.
            IF (NCHBAR(IS,I) > 0) THEN
              VISITED(NSEITE(IS,I),NCHBAR(IS,I)) = .TRUE.
            END IF
            IF (NLTRA) THEN
              XX1=XX1+RMTOR
              XX2=XX2+RMTOR
              COM=0.5*(XX1+XX2)
cdr           write (iunout,'(3i6,2es12.4)')
cdr  .            i, ixtri(i), iytri(i), dsd, com

              DSD=DSD*COM*TANAL/ALPHA*PI2A
            ELSE
              DSD=DSD*ZDF
            ENDIF
            IF (NLMPGS.GT.NLIMPS) THEN
              MSURFG=NLIM+NSTS+INSPAT(IS,I)
              SAREA(MSURFG)=DSD
              IF (NCHBAR(IS,I) > 0) THEN
                MSURFG=NLIM+NSTS+INSPAT(NSEITE(IS,I),NCHBAR(IS,I))
                SAREA(MSURFG)=DSD
              END IF
            END IF
            SAREA(ISTS)=SAREA(ISTS)+DSD
          END DO
cdr       write (iunout,*) ' sarea ', ists, SAREA(ISTS)
        END DO

        DEALLOCATE (VISITED)

      case (5)
        DO ISTS=1,NSTSI
          NLJ=NLIM+ISTS
          SAREA(NLJ)=0.
        ENDDO
        DO ITET=1,NTET
          DO J=1,4
            ISTS=ABS(INMTIT(J,ITET))
            IF (ISTS .GT. 0) THEN
              IC1=NTECK(ITSIDE(1,J),ITET)
              IC2=NTECK(ITSIDE(2,J),ITET)
              IC3=NTECK(ITSIDE(3,J),ITET)
!pb              NLJ=NLIM+ISTS   ! INMTIT CHANGED IN INFCOP
              SAREA(ISTS)=SAREA(ISTS)+
     .            EIRENE_ARTRI3(XTETRA(IC1),YTETRA(IC1),ZTETRA(IC1),
     .                          XTETRA(IC2),YTETRA(IC2),ZTETRA(IC2),
     .                          XTETRA(IC3),YTETRA(IC3),ZTETRA(IC3))
            END IF
          END DO
        END DO
C     ELSEIF (LEVGEO.EQ....) THEN
      end select

!  SET NSTGRD FOR AVERAGING CELLS

      IR = NR1ST
      DO IT = 1, NT3RD
        DO IP = 1, NP2ND
          NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
          NSTGRD(NCELL) = 3
        END DO
      END DO

C
C   POLOIDAL OR Y-GRID
C
      CASE (2) ! IND
C
C  IF NLSYMP, Y-GRID MUST BE SYMMETRIC: PSURF(I)=PSURF(NP2ND-I+1)
C
      select case (LEVGEO)
      case (1) !  CARTESIAN, Y-DIRECTION
C   Y-GRID
        IF (INDGRD(IND).LE.4) THEN
          CALL EIRENE_GRID_1(PSURF,NP2ND,NPSEP,NPPLA,YIA,YGA,YAA,YYA,2)
C       ELSEIF (INDGRD(IND).EQ.5) THEN
C** TAKE Y GRID DATA FROM USER-SUPPLIED SUBROUTINE
C TO BE WRITTEN
C         CALL PROUSR (PSURF,2+4*NPLS+3,0._DP,0._DP,0._DP,0._DP,0._DP,0._DP,0._D
C       ELSEIF (INDGRD(IND).EQ.6) THEN
C** TAKE Y GRID DATA FROM INTERFACE
C TO BE WRITTEN
C         CALL PROFR (PSURF,2+4*NPLS+5,1,1,NP2ND)
        ENDIF
C
        DO 210 J=1,NP2NDM
          PHZONE(J)=(PSURF(J)+PSURF(J+1))/2.
  210   CONTINUE
C
        IF (TRCGRD) THEN
          CALL EIRENE_LEER(1)
          WRITE (iunout,*) 'GRIDPOINTS IN Y DIRECTION '
          CALL EIRENE_LEER(1)
          CALL EIRENE_MASRR1('  N, PSURF ',PSURF,NP2ND,3)
          CALL EIRENE_LEER(2)
        ENDIF
C
        CALL EIRENE_SNEIGH
C
      case (2) ! POLAR ANGLE GRID, "THETA"-GRID, IN RADIANS
C
        IF (INDGRD(IND).EQ.1) THEN
          ND=NPSEP
          DO 211 J=1,ND
            PSURF(J)=(YIA+DBLE((J-1))/DBLE(ND-1)*(YGA-YIA))*DEGRAD
  211     CONTINUE
          DO 212 J=ND+1,NP2ND
            PSURF(J)=(YGA+DBLE(J-ND)/DBLE(NP2ND-ND)*(YAA-YGA))*DEGRAD
  212     CONTINUE
        ELSEIF (INDGRD(IND).EQ.2) THEN
          DO 213 J=1,NP2ND
            PSURF(J)=(YIA+(J-1)/DBLE(NP2NDM)*(YAA-YIA))*DEGRAD
  213     CONTINUE
        ELSE
          WRITE (iunout,*) 'INVALID OPTION ENCOUNTERED IN SUBR. GRID'
          WRITE (iunout,*) 'INDGRD(2),LEVGEO ',INDGRD(IND),LEVGEO
          WRITE (iunout,*) 'EXIT CALLED FROM SUBR. GRID'
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
        DO 215 J=1,NP2ND
          COSPH(J)=COS(PSURF(J))
          SINPH(J)=SIN(PSURF(J))
  215   CONTINUE
C
        IF (TRCGRD) THEN
          CALL EIRENE_LEER(1)
          WRITE (iunout,*) 'GRIDPOINTS IN POLOIDAL DIRECTION '
          CALL EIRENE_LEER(1)
          CALL EIRENE_MASRR1('  N, PSURF ',PSURF,NP2ND,3)
          CALL EIRENE_LEER(2)
        ENDIF
C
C
        NPPLG=1
        NRPLG=NP2ND
        NPOINT(1,1)=1
        NPOINT(2,1)=NP2ND
        IFLAG=2
        DO 1240 IR=1,NR1STM
          IRP=IR+1
          DO 1250 IP=1,NRPLG-1
            IPP=IP+1
            CALL EIRENE_ARELLP(EP1(IRP),EP1(IR),ELL(IRP),ELL(IR),
     .                  TRI(IRP),TRI(IR),
     .                  RSURF(IRP),RSURF(IR),PSURF(IPP),PSURF(IP),IFLAG,
     .                  AELL,SX,SY,X1,Y1,X2,Y2,X3,Y3,X4,Y4)
C
            XPOL(IRP,IPP)=X1
            YPOL(IRP,IPP)=Y1
            XPOL(IR,IPP)=X2
            YPOL(IR,IPP)=Y2
            XPOL(IRP,IP)=X3
            YPOL(IRP,IP)=Y3
            XPOL(IR,IP)=X4
            YPOL(IR,IP)=Y4
 1250     CONTINUE
 1240   CONTINUE
C
        CALL EIRENE_SNEIGH

      end select
C
      select case (LEVGEO)
      case (2:3)
C
        IF (TRCGRD) THEN
          DO 219 I=1,NRPLG
            WRITE (iunout,*) ' PERP. POLYGON NO. I = ',I
            WRITE (iunout,*) ' JA = ',1,' JE = ',NR1ST
            WRITE (iunout,'(/1X,1P,6E12.4)') (XPOL(K,I),YPOL(K,I),
     .             K=1,NR1ST)
  219     CONTINUE
        ENDIF
C
        DO 220 K=1,NRPLG
          DO I=1,NR1STM
            VVTX(I,K)=XPOL(I+1,K)-XPOL(I,K)
            VVTY(I,K)=YPOL(I+1,K)-YPOL(I,K)
          END DO
  220   CONTINUE
C
        DO 221 K=1,NRPLG
          BGLP(1,K)=0.
          DO 222 I=1,NR1STM
            BGLP(I+1,K)=BGLP(I,K)+SQRT(VVTX(I,K)**2+VVTY(I,K)**2)
  222     CONTINUE
  221   CONTINUE
C
        IF (TRCGRD) THEN
          CALL EIRENE_LEER(2)
          WRITE (iunout,*)
     .      'ARCLENGTH BGLP(I,K) OF POLOIDAL SURFACES AT Z=0.'
          DO 223 K=1,NRPLG
            WRITE (iunout,*) 'K = ',K
            WRITE (iunout,'(/1X,1P,6E12.4)') (BGLP(I,K),I=1,NR1ST)
            CALL EIRENE_LEER(1)
  223     CONTINUE
        ENDIF
C
C   CALCULATE THE OUTER NORMALS OF POLYGONS
C
        DO 224 K=1,NRPLG
          DO 225 I=1,NR1STM
            IF (ABS(VVTY(I,K)).LT.EPS12) THEN
              PPLNX(I,K)=0.
              PPLNY(I,K)=1.
            ELSE
              PPLNX(I,K)=1.
              PPLNY(I,K)=-VVTX(I,K)/VVTY(I,K)
            ENDIF
  225     CONTINUE
          DO I=1,NR1STM
            PLABS=SQRT(PPLNX(I,K)**2+PPLNY(I,K)**2)
            PPLNX(I,K)=PPLNX(I,K)/PLABS
            PPLNY(I,K)=PPLNY(I,K)/PLABS
          END DO
  224   CONTINUE
C
        DO 227 I=1,NR1STM
         DO J=1,NPPLG
          DO K=NPOINT(1,J),NPOINT(2,J)
            KUP=K+1
            KDN=K
            IF (KUP.GT.NPOINT(2,J)) THEN
              KUP=K
              KDN=K-1
            ENDIF
  226       XD=XPOL(I+1,KUP)-XPOL(I+1,KDN)
            YD=YPOL(I+1,KUP)-YPOL(I+1,KDN)
            IF (XD*XD+YD*YD.LT.EPS30) THEN
              IF (KUP.LT.NPOINT(2,J)) THEN
                KUP=KUP+1
              ELSE
                KDN=KDN-1
              ENDIF
              GOTO 226
            ENDIF
            XS=SIGN(1._DP,XD*PPLNX(I,K)+YD*PPLNY(I,K))
            PPLNX(I,K)=PPLNX(I,K)*XS
            PPLNY(I,K)=PPLNY(I,K)*XS
          END DO
         END DO
  227   CONTINUE
C
C  IDENTIFY DEAD CELLS IN GRID CUTS
C
        IT=1
        DO 229 IR=1,NR1ST-1
         DO J=1,NPPLG-1
          DO IP=NPOINT(2,J),NPOINT(1,J+1)-1
            NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
            NSTGRD(NCELL)=2
          END DO
         END DO
  229   CONTINUE
C
      case (4:)
C
        WRITE (iunout,*) 'ERROR EXIT FROM GRID. NLPOL ',LEVGEO
      end select
C
C  1ST AND 2ND GRID DEFINED
C  SET SURFACE AREA OF NON-DEFAULT STANDARD SURFACES
C
C  RADIAL (1ST GRID) SURFACES. OVERWRITE EARLIER VALUES FROM
C  CALL EIRENE_GRID(1)
C
      select case (LEVGEO)
      case (1)
        DO 280 ISTS=1,NSTSI
          IF (INUMP(ISTS,1).NE.0) THEN
            IR=INUMP(ISTS,1)
            NLJ=NLIM+ISTS
            IF (NLTRZ) THEN
              SAREA(NLJ)=PSURF(IRPTE(ISTS,2))-PSURF(IRPTA(ISTS,2))
              SAREA(NLJ)=SAREA(NLJ)*ZDF
            ELSEIF (NLTRA) THEN
              YDF=PSURF(IRPTE(ISTS,2))-PSURF(IRPTA(ISTS,2))
              SAREA(NLJ)=YDF*ZDF*(RSURF(IR)+RMTOR)
            ELSEIF (NLTRT) THEN
            ENDIF
          ENDIF
  280   CONTINUE
      case (2:3)
        DO 290 ISTS=1,NSTSI
          IF (INUMP(ISTS,1).NE.0) THEN
            IR=INUMP(ISTS,1)
            NLJ=NLIM+ISTS
cdr         write (iunout,*) ' area for surface ',nlj
            IF (NLTRZ) THEN
              SAREA(NLJ)=BGL(IR,IRPTE(ISTS,2))-BGL(IR,IRPTA(ISTS,2))
              SAREA(NLJ)=SAREA(NLJ)*ZDF
            ELSEIF (NLTRA) THEN
              SAREA(NLJ)=0.
              DO 291 IP=IRPTA(ISTS,2),IRPTE(ISTS,2)-1
                XS=((XPOL(IR,IP+1)+XPOL(IR,IP))*0.5)+RMTOR
                SAREA(NLJ)=SAREA(NLJ)+(BGL(IR,IP+1)-BGL(IR,IP))*XS
cdr             write (iunout,'(2i6,2es12.4)')
cdr  .             ir, ip, BGL(IR,IP+1)-BGL(IR,IP), xs
  291         CONTINUE
              SAREA(NLJ)=SAREA(NLJ)*TANAL/ALPHA*PI2A
cdr           write (iunout,*) 'sarea ', nlj, SAREA(NLJ)
            ENDIF
          ENDIF
  290   CONTINUE
      case default
C TO BE WRITTEN
      end select
C
C  POLOIDAL (2ND GRID) SURFACES.
C
      select case (LEVGEO)
      case (1)
        DO 285 ISTS=1,NSTSI
          IF (INUMP(ISTS,2).NE.0) THEN
            IP=INUMP(ISTS,2)
            NLJ=NLIM+ISTS
            IF (NLTRZ) THEN
              SAREA(NLJ)=RSURF(IRPTE(ISTS,1))-RSURF(IRPTA(ISTS,1))
              SAREA(NLJ)=SAREA(NLJ)*ZDF
            ELSEIF (NLTRA) THEN
              SAREA(NLJ)=0.
              DO 286 IR=IRPTA(ISTS,1),IRPTE(ISTS,1)-1
                XS=(RSURF(IR+1)+RSURF(IR))*0.5+RMTOR
                SAREA(NLJ)=SAREA(NLJ)+(RSURF(IR+1)-RSURF(IR))*XS
  286         CONTINUE
              SAREA(NLJ)=SAREA(NLJ)*TANAL/ALPHA*PI2A
            ELSEIF (NLTRT) THEN
              SAREA(NLJ)=0.
              DO 287 IR=IRPTA(ISTS,1),IRPTE(ISTS,1)-1
                XS=(RSURF(IR+1)+RSURF(IR))*0.5*2.*PIA
                SAREA(NLJ)=SAREA(NLJ)+(RSURF(IR+1)-RSURF(IR))*XS
  287         CONTINUE
            ENDIF
          ENDIF
  285   CONTINUE
      case (2:3)
        DO 295 ISTS=1,NSTSI
          IF (INUMP(ISTS,2).NE.0) THEN
            IP=INUMP(ISTS,2)
            NLJ=NLIM+ISTS
            IF (NLTRZ) THEN
              SAREA(NLJ)=BGLP(IRPTE(ISTS,1),IP)-BGLP(IRPTA(ISTS,1),IP)
              SAREA(NLJ)=SAREA(NLJ)*ZDF
            ELSEIF (NLTRA) THEN
              SAREA(NLJ)=0.
              DO 296 IR=IRPTA(ISTS,1),IRPTE(ISTS,1)-1
                XS=(XPOL(IR+1,IP)+XPOL(IR,IP))*0.5+RMTOR
                SAREA(NLJ)=SAREA(NLJ)+(BGLP(IR+1,IP)-BGLP(IR,IP))*XS
  296         CONTINUE
              SAREA(NLJ)=SAREA(NLJ)*TANAL/ALPHA*PI2A
            ELSEIF (NLTRT) THEN
              SAREA(NLJ)=0.
              DO 297 IR=IRPTA(ISTS,1),IRPTE(ISTS,1)-1
                XS=(XPOL(IR+1,IP)+XPOL(IR,IP))*0.5*2.*PIA
                SAREA(NLJ)=SAREA(NLJ)+(BGLP(IR+1,IP)-BGLP(IR,IP))*XS
  297         CONTINUE
            ENDIF
          ENDIF
  295   CONTINUE
      case default
C TO BE WRITTEN
      end select

!  SET NSTGRD FOR AVERAGING CELLS

      IP = NP2ND
      DO IT = 1, NT3RD
        DO IR = 1, NR1ST
          NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
          NSTGRD(NCELL) = 3
        END DO
      END DO

C
C   TOROIDAL OR Z-GRID
C
      CASE (3) ! IND
C
C  IF NLSYMT, Z-GRID MUST BE SYMMETRIC
C
      IF (NLTRZ) THEN
C   Z-GRID, CARTESIAN
        IF (INDGRD(IND).LE.4) THEN
          CALL EIRENE_GRID_1(ZSURF,NT3RD,NTSEP,NTTRA,ZIA,ZGA,ZAA,ZZA,3)
C       ELSEIF (INDGRD(IND).EQ.5) THEN
C** TAKE Z GRID DATA FROM USER-SUPPLIED SUBROUTINE
C TO BE WRITTEN
C         CALL PROUSR (ZSURF,2+4*NPLS+3,0._DP,0._DP,0._DP,0._DP,0._DP,0._DP,0._D
C       ELSEIF (INDGRD(IND).EQ.6) THEN
C** TAKE Z GRID DATA FROM INTERFACE
C TO BE WRITTEN
C         CALL PROFR (ZSURF,2+4*NPLS+5,1,1,NT3RD)
        ENDIF
C
        DO 310 J=1,NT3RDM
          ZZONE(J)=(ZSURF(J)+ZSURF(J+1))/2.
  310   CONTINUE
C
C     ELSEIF (NLTRA) THEN
C   GRID FOR TOROIDAL APPROXIMATION OF CYLINDER: ALREADY DONE IN CALL EIRENE_GRID(1)
C   THERE: ZSURF, ZZONE HAVE BEEN DEFINED IN RADIANS
C
      ENDIF
C
      IF (TRCGRD) THEN
        CALL EIRENE_LEER(1)
        IF (NLTRZ)
     .   WRITE (iunout,*) 'GRIDPOINTS IN Z DIRECTION'
        IF (NLTRA)
     .   WRITE (iunout,*) 'GRIDPOINTS IN TOROIDAL DIRECTION, IN RADIANS'
        CALL EIRENE_LEER(1)
        CALL EIRENE_MASRR1 (' N,  ZSURF ',ZSURF,NT3RD,3)
        CALL EIRENE_LEER(2)
      ENDIF


!  SET NSTGRD FOR CELLS IT=NT3RD CONTAINING 2-DIMENSIONAL AVERAGES

      IT = NT3RD
      DO IR = 1, NR1ST
        DO IP = 1, NP2ND
          NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
          NSTGRD(NCELL) = 3
        END DO
      END DO

!  COPY SWITCHING OFF OF DEAD CELLS FOR TOROIDAL CELL 1 TO ALL OTHERS
!  i.e.:
!  SET NSTGRD FOR CELLS IT CONTAINING 0D AND 1D AVERAGES OF 2D PROJECTION

      IT1=1
      DO IR=1,NR1ST
        DO IP=1,NP2ND
          NCELL1=IR+((IP-1)+(IT1-1)*NP2T3)*NR1P2
          DO IT=2,NT3RD-1
            NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
            NSTGRD(NCELL)=NSTGRD(NCELL1)
          END DO
        END DO
      END DO
C
      END SELECT ! IND
      RETURN
C
  991 CONTINUE
      WRITE (iunout,*) 'GRID DATA INCONSISTENCY: 3RD GRID.  ZAA > ZIA ?'
      WRITE (iunout,*) 'ZIA,ZAA,NTTRA,ROA= ',ZIA,ZAA,NTTRA,ROA
      CALL EIRENE_EXIT_OWN(1)
  992 CONTINUE
      WRITE (iunout,*) 'GRID DATA INCONSISTENCY: 2ND GRID.  YAA > YIA ?'
      WRITE (iunout,*) 'YIA,YAA = ',YIA,YAA
      CALL EIRENE_EXIT_OWN(1)
      END
