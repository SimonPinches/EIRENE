      MODULE EIRMOD_TIMEA
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_CADGEO
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CGRID
      USE EIRMOD_CLGIN
      USE EIRMOD_CTRIG
      USE EIRMOD_COMSPL
      USE EIRMOD_CPLOT
      USE EIRMOD_COMPRT, ONLY: IUNOUT, IVTKOUT
      USE EIRMOD_CTRCEI, only: trcoct
      USE EIRMOD_CPES

c     in case: NLOCTREE...
c     using special octree stuff, to speed up additional surface handling
      USE EIRMOD_OCTREE
      IMPLICIT NONE
      PRIVATE

      PUBLIC :: EIRENE_TIMEA0, EIRENE_TIMEA3_VTKOUT,
     .          EIRENE_TIMEA1, EIRENE_TIMEA0_OC,
     .          EIRENE_DEALLOC_TIMEA

c     the octree
      TYPE(octree), POINTER, SAVE :: TREE
c     active surface ids of surfs not in octree saved in lookup table
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: NOTOCSURFS
c     how many surfaces are not in the octree (NSURFNOT)
      INTEGER, SAVE :: NSURFNOT = 0

      CONTAINS

Cdr:  at some point (after 2010) this routine was extended by an "octree option" to
c     to speed up calculation for intersections of rays with add. surfaces
c     implemented as triangles in 3d space (rlb >=3 options>).

cdr Jan 2017:  the octree option NLOCTREE is part of a proprietary development branch
c              at FZJ, and should not be used by a 3rd party. It is likely not finished.
cdr see: TODO comment below. Further: is the sequence of additional surfaces changed
c        by doing the RLB< 3 surfaces first and the RLB >=3 surfaces next?
c        This may invalidate some of the surface options which depend on the
c        order of additional surfaces, see manual, input block 3b
cdr  Feb 2017:  some documentation, and clearer separation of NLOCTREE options
cdr             from the default (.NOT.NLOCTREE) procedure.  More to do.....
C
      SUBROUTINE EIRENE_TIMEA0
C
C   1 ST INTERSECTION OF THE RAY X+T*VX,Y+T*VY,Z+T*VZ WITH ONE OF THE NLIM
C   ADDITIONAL SURFACES, DEFINED BY 1.ST OR 2.ND ORDER EQUATIONS IN INPUT BLOCK 3B.
C   IT IS ALSO CHECKED, WHETHER THIS INTERSECTION TAKES PLACE INSIDE THE
C   SPECIFIED BOUNDARIES OF THOSE SURFACES
C
C  ENTRY TIMEA0: INITIALIZATION, PRE-COMPUTING SOME EXPRESSIONS, CONSISTENCY CHECKS, ETC..
C                IN CASE: NLOCTREE: SET UP OCTREE, FOR SPEED UP OF GEOMETRICAL CALCULATIONS
C  ENTRY TIMEA1: FIND CLOSEST LEGAL INTERSECTION POINT. MORE COMMENTS: SEE AT ENTRY TIMEA1.
C
C

      IMPLICIT NONE

C
      REAL(DP) :: XB(3), XC(3), XD(3), XE(3), XF(3), XG(3)
      REAL(DP) :: DZ3, XS, DZ, TEST4, TEST5, DY, DL, DX1, DX,
     .          B, XK5, V445, XN5, V335, XK4,
     .          V334, V3V45, XN4, AT, XNORM, V224, CG1, CG2,
     .          CG3, AR1, C, D, E, F, G, XK3, HELP, V2V34, XN3,
     .          AR2, V113, V1V23, V223, AR3, DY2
C      REAL(DP) :: XLS1, XLS2, XLS3, XMS1, XMS2, XMS3

      INTEGER :: IPLGN, ITRII,
     .           ISTS, I1000, J, I,
     .           ISPZ, IAB, JUM
      INTEGER, EXTERNAL :: EIRENE_IDEZ

      SAVE
C
C
C  INITIALIZATION: SET SOME SWITCHES, PRECOMPUTE ALGEBRAIC EXPRESSIONS
C                  TO SPEED UP LATER COMPUTATION OF INTERSECTION POINTS
C                  IN MONTE CARLO LOOP
C
      IF (NLIMI.LT.1) RETURN
C
C
      DO J=1,NLIMI
        IF (IGJUM0(J).NE.0) THEN
          IF (NLIMPB >= NLIMPS) THEN
            DO I=0,NLIMI
              IGJUM1(I,J)=1
            END DO
          ELSE
            DO I=0,NLIMI
              CALL EIRENE_BITSET (IGJUM1,0,NLIMPS,I,J,1,NBITS)
            END DO
          END IF
        END IF
C
        ISWICH(1,J)=EIRENE_IDEZ(ILSWCH(J),1,6)
        IF (ISWICH(1,J).EQ.1) ISWICH(1,J)=-1
        IF (ISWICH(1,J).EQ.2) ISWICH(1,J)=1
        ISWICH(2,J)=EIRENE_IDEZ(ILSWCH(J),2,6)
        IF (ISWICH(2,J).EQ.1) ISWICH(2,J)=-1
        IF (ISWICH(2,J).EQ.2) ISWICH(2,J)=1
        ISWICH(3,J)=EIRENE_IDEZ(ILSWCH(J),3,6)
        IF (ISWICH(3,J).EQ.1) ISWICH(3,J)=-1
        IF (ISWICH(3,J).EQ.2) ISWICH(3,J)=1
        ISWICH(4,J)=EIRENE_IDEZ(ILSWCH(J),4,6)
        IF (ISWICH(4,J).EQ.1) ISWICH(4,J)=-1
        IF (ISWICH(4,J).EQ.2) ISWICH(4,J)=1
        ISWICH(5,J)=EIRENE_IDEZ(ILSWCH(J),5,6)
        IF (ISWICH(5,J).EQ.1) ISWICH(5,J)=-1
        IF (ISWICH(5,J).EQ.2) ISWICH(5,J)=1
        ISWICH(6,J)=EIRENE_IDEZ(ILSWCH(J),6,6)
        IF (ISWICH(6,J).EQ.1) ISWICH(6,J)=-1
        IF (ISWICH(6,J).EQ.2) ISWICH(6,J)=1
C
        IF (ISWICH(4,J).NE.0.OR.ISWICH(5,J).NE.0.OR.ISWICH(6,J).NE.0)
     .  THEN
          ILBLCK(J)=EIRENE_IDEZ(ILCELL(J),4,4)
          I1000=1000*ILBLCK(J)
          ILACLL(J)=ILCELL(J)-I1000
        ENDIF
      END DO

C     iterate over all additional surfaces
      DO 97 J=1,NLIMI
        IF (IGJUM0(J).NE.0) THEN
          IF (LEVGEO.EQ.4) THEN
            DO ITRII=1,NTRII
              DO IPLGN=1,3
                ISTS=ABS(INMTI(IPLGN,ITRII))
                IF (J.EQ.ISTS) GOTO 85
              ENDDO
            ENDDO
          ENDIF
          GOTO 97
        ENDIF

   85   IF (RLB(J).LT.2.0) THEN
C
C   SURFACE COEFFICIENTS ARE INPUT
C
C  INDICATE INDEPENDENCE OF SURFACE EQUATION FROM X, Y, Z, RESP.
          IF (A1LM(J).EQ.0..AND.A4LM(J).EQ.0..AND.
     .        A7LM(J).EQ.0..AND.A8LM(J).EQ.0.)
     .    P3(1,J)=1.D55
          IF (A2LM(J).EQ.0..AND.A5LM(J).EQ.0..AND.
     .        A7LM(J).EQ.0..AND.A9LM(J).EQ.0.)
     .    P3(2,J)=1.D55
          IF (A3LM(J).EQ.0..AND.A6LM(J).EQ.0..AND.
     .        A8LM(J).EQ.0..AND.A9LM(J).EQ.0.)
     .    P3(3,J)=1.D55
          GOTO 90
        ELSEIF (RLB(J).GE.2) THEN
C
C  PLANE SURFACE, VERTICES ARE INPUT, SET SURFACE COEFFICIENTS:
C
          A4LM(J)=0.
          A5LM(J)=0.
          A6LM(J)=0.
          A7LM(J)=0.
          A8LM(J)=0.
          A9LM(J)=0.
        ENDIF
C
        IF (RLB(J).GE.2.6.AND.RLB(J).LT.2.9) THEN
C  2 POINTS, PLANE SURFACE WITH IGNORABLE X COORDINATE
          IF (RLB(J).LT.2.75) RLB(J)=1.
          IF (RLB(J).GE.2.75) RLB(J)=1.5
          XLIMS1(1,J)=P1(1,J)
          XLIMS2(1,J)=P2(1,J)
          DZ3=(P1(3,J)-P2(3,J))
          DY2=(P1(2,J)-P2(2,J))
!pb  avoid a0=a1=a2=a3=0
          IF (ABS(DZ3)+ABS(DY2) .LT. EPS12) GOTO 98
          A0LM(J)=DZ3*P1(2,J)-DY2*P1(3,J)
          A1LM(J)=0.
          A2LM(J)=-DZ3
          A3LM(J)=DY2
          YLIMS1(1,J)=MIN(P1(2,J),P2(2,J))
          YLIMS2(1,J)=MAX(P1(2,J),P2(2,J))
          ZLIMS1(1,J)=MIN(P1(3,J),P2(3,J))
          ZLIMS2(1,J)=MAX(P1(3,J),P2(3,J))
!
!          write (iunout,'(a,6es12.4)') 'xylims1,2 ',
!     .      xlims1(1,j),ylims1(1,j),zlims1(1,j),
!     .      xlims2(1,j),ylims2(1,j),zlims2(1,j)
!
          IF (P1(2,J).EQ.P2(2,J)) THEN
            YLIMS1(1,J)=YLIMS1(1,J)-0.1
            YLIMS2(1,J)=YLIMS2(1,J)+0.1
          ENDIF
          IF (P1(3,J).EQ.P2(3,J)) THEN
            ZLIMS1(1,J)=ZLIMS1(1,J)-0.1
            ZLIMS2(1,J)=ZLIMS2(1,J)+0.1
          ENDIF
C  INDICATE 2-POINT OPTION (BECAUSE RLB IS OVERWRITTEN)
C  AND X-INDEPENDENCE OF SURFACE EQUATION
          P3(1,J)=1.D55
          DL=SQRT(DY2**2+DZ3**2)
          DX=XLIMS2(1,J)-XLIMS1(1,J)
          IF (DX.GT.1.D20) DX=XDF
          SAREA(J)=DL*DX
          GOTO 90
        ELSEIF (RLB(J).GE.2.3.AND.RLB(J).LT.2.6) THEN
C  2 POINTS, PLANE SURFACE WITH IGNORABLE Y COORDINATE
          IF (RLB(J).LT.2.45) RLB(J)=1.
          IF (RLB(J).GE.2.45) RLB(J)=1.5
          YLIMS1(1,J)=P1(2,J)
          YLIMS2(1,J)=P2(2,J)
          DZ3=(P1(3,J)-P2(3,J))
          DX1=(P1(1,J)-P2(1,J))
!pb  avoid a0=a1=a2=a3=0
          IF (ABS(DZ3)+ABS(DX1) .LT. EPS12) GOTO 98
          A0LM(J)=DZ3*P1(1,J)-DX1*P1(3,J)
          A1LM(J)=-DZ3
          A2LM(J)=0.
          A3LM(J)=DX1
          XLIMS1(1,J)=MIN(P1(1,J),P2(1,J))
          XLIMS2(1,J)=MAX(P1(1,J),P2(1,J))
          ZLIMS1(1,J)=MIN(P1(3,J),P2(3,J))
          ZLIMS2(1,J)=MAX(P1(3,J),P2(3,J))
!
!          write (iunout,'(a,6es12.4)') 'xylims1,2 ',
!     .      xlims1(1,j),ylims1(1,j),zlims1(1,j),
!     .      xlims2(1,j),ylims2(1,j),zlims2(1,j)
!
          IF (P1(1,J).EQ.P2(1,J)) THEN
            XLIMS1(1,J)=XLIMS1(1,J)-0.1
            XLIMS2(1,J)=XLIMS2(1,J)+0.1
          ENDIF
          IF (P1(3,J).EQ.P2(3,J)) THEN
            ZLIMS1(1,J)=ZLIMS1(1,J)-0.1
            ZLIMS2(1,J)=ZLIMS2(1,J)+0.1
          ENDIF
C  INDICATE 2-POINT OPTION (BECAUSE RLB IS OVERWRITTEN)
C  AND Y-INDEPENDENCE OF SURFACE EQUATION
          P3(2,J)=1.D55
          DL=SQRT(DX1**2+DZ3**2)
          DY=YLIMS2(1,J)-YLIMS1(1,J)
          IF (DY.GT.1.E20) DY=YDF
          SAREA(J)=DL*DY
          GOTO 90
        ELSEIF (RLB(J).GE.2.0.AND.RLB(J).LT.2.3) THEN
C  2 POINTS, PLANE SURFACE WITH IGNORABLE Z COORDINATE
          IF (RLB(J).LT.2.15) RLB(J)=1.
          IF (RLB(J).GE.2.15) RLB(J)=1.5
          ZLIMS1(1,J)=P1(3,J)
          ZLIMS2(1,J)=P2(3,J)
          DY2=(P1(2,J)-P2(2,J))
          DX1=(P1(1,J)-P2(1,J))
!pb  avoid a0=a1=a2=a3=0
          IF (ABS(DY2)+ABS(DX1) .LT. EPS12) GOTO 98
          A0LM(J)=DY2*P1(1,J)-DX1*P1(2,J)
          A1LM(J)=-DY2
          A2LM(J)=DX1
          A3LM(J)=0.
          XLIMS1(1,J)=MIN(P1(1,J),P2(1,J))
          XLIMS2(1,J)=MAX(P1(1,J),P2(1,J))
          YLIMS1(1,J)=MIN(P1(2,J),P2(2,J))
          YLIMS2(1,J)=MAX(P1(2,J),P2(2,J))
!
!          write (iunout,'(a,6es12.4)') 'xylims1,2 ',
!     .      xlims1(1,j),ylims1(1,j),zlims1(1,j),
!     .      xlims2(1,j),ylims2(1,j),zlims2(1,j)
!
          IF (P1(1,J).EQ.P2(1,J)) THEN
            XLIMS1(1,J)=XLIMS1(1,J)-0.1
            XLIMS2(1,J)=XLIMS2(1,J)+0.1
          ENDIF
          IF (P1(2,J).EQ.P2(2,J)) THEN
            YLIMS1(1,J)=YLIMS1(1,J)-0.1
            YLIMS2(1,J)=YLIMS2(1,J)+0.1
          ENDIF

C  INDICATE 2-POINT OPTION (BECAUSE RLB IS OVERWRITTEN)
C  AND Z-INDEPENDENCE OF SURFACE EQUATION
          P3(3,J)=1.D55
          DL=SQRT(DX1**2+DY2**2)
          IF (NLTRZ) THEN
            DZ=ZLIMS2(1,J)-ZLIMS1(1,J)
            IF (DZ.GT.1.D20) DZ=ZDF
          ELSEIF (NLTRA) THEN
            XS=(P1(1,J)+P2(1,J))*0.5+RMTOR
            IF (ILTOR(J).GT.0) THEN
              DZ=ZLIMS2(1,J)-ZLIMS1(1,J)
C             IF (DZ.GT.1.D20) ????
            ELSEIF (ILTOR(J).EQ.0) THEN
              DZ=ZLIMS2(1,J)-ZLIMS1(1,J)
              IF (DZ.GT.1.D20) DZ=XS*TANAL/ALPHA*PI2A
            ENDIF
          ENDIF
          SAREA(J)=DL*DZ
          GOTO 90
        ELSEIF (RLB(J).GE.3.AND.RLB(J).LT.4.) THEN
C  1 TRIANGLE
          P4(1,J)=P3(1,J)
          P4(2,J)=P3(2,J)
          P4(3,J)=P3(3,J)
          P5(1,J)=P4(1,J)
          P5(2,J)=P4(2,J)
          P5(3,J)=P4(3,J)
        ELSEIF (RLB(J).GE.4..AND.RLB(J).LT.5) THEN
C  1 QUADRANGLE = 2 TRIANGLES
          P5(1,J)=P4(1,J)
          P5(2,J)=P4(2,J)
          P5(3,J)=P4(3,J)
        ENDIF
C
        TEST4=0.
        TEST5=0.
C
C  SET PLANE SURFACE FROM COORDINATES OF THE FIRST 3 VERTICES
C  RLB.GE.3.0 AT THIS POINT
C
c     vektoren bilden (alle die in einem fuenfeck vorkommen koennen)
c     -> bei dreiecken/vierecken sind halt welche gleich...
        DO I=1,3
          XB(I)=P1(I,J)-P3(I,J)
          XC(I)=P2(I,J)-P3(I,J)
          XD(I)=P2(I,J)-P4(I,J)
          XE(I)=P3(I,J)-P4(I,J)
          XF(I)=P3(I,J)-P5(I,J)
          XG(I)=P4(I,J)-P5(I,J)
        END DO
C     coordinaten darstellung der ebene
        A1LM(J)=XB(2)*XC(3)-XB(3)*XC(2)
        A2LM(J)=XB(3)*XC(1)-XB(1)*XC(3)
        A3LM(J)=XB(1)*XC(2)-XB(2)*XC(1)
        A0LM(J)=-(A1LM(J)*P1(1,J)+A2LM(J)*P1(2,J)+A3LM(J)*P1(3,J))
C
C   SURFACE AREA: SAREA (flaeche berechnen)
        B=SQRT(XB(1)**2+XB(2)**2+XB(3)**2+EPS60)
        C=SQRT(XC(1)**2+XC(2)**2+XC(3)**2+EPS60)
        D=SQRT(XD(1)**2+XD(2)**2+XD(3)**2+EPS60)
        E=SQRT(XE(1)**2+XE(2)**2+XE(3)**2+EPS60)
        F=SQRT(XF(1)**2+XF(2)**2+XF(3)**2+EPS60)
        G=SQRT(XG(1)**2+XG(2)**2+XG(3)**2+EPS60)
        CG1=(XB(1)*XC(1)+XB(2)*XC(2)+XB(3)*XC(3))/B/C
        CG2=(XD(1)*XE(1)+XD(2)*XE(2)+XD(3)*XE(3))/D/E
        CG3=(XF(1)*XG(1)+XF(2)*XG(2)+XF(3)*XG(3))/F/G
        AR1=0.5*B*C*SQRT(1.-CG1*CG1+EPS60)
        AR2=0.5*D*E*SQRT(1.-CG2*CG2+EPS60)
        AR3=0.5*F*G*SQRT(1.-CG3*CG3+EPS60)
        SAREA(J)=AR1+AR2+AR3
C
C   TEST, IF 4TH  AND 5TH POINT ON SURFACE
        TEST4=A0LM(J)+A1LM(J)*P4(1,J)+A2LM(J)*P4(2,J)+A3LM(J)*P4(3,J)
        TEST5=A0LM(J)+A1LM(J)*P5(1,J)+A2LM(J)*P5(2,J)+A3LM(J)*P5(3,J)
C
C  NEW COORDINATE SYSTEM IN P3,P1,P2 ; P4,P2,P3 ; P5,P3,P4
C    XS = P3 + XLS1*(P1-P3) + XMS1*(P2-P3)
C    XS = P4 + XLS2*(P2-P4) + XMS2*(P3-P4)
C    XS = P5 + XLS3*(P3-P5) + XMS3*(P4-P5)
C  PREPARE ARRAYS FOR COMPUTATION OF XLS1,...XMS3
        V1V23=XB(1)*XC(1)+XB(2)*XC(2)+XB(3)*XC(3)
        V113=XB(1)*XB(1)+XB(2)*XB(2)+XB(3)*XB(3)
        V223=XC(1)*XC(1)+XC(2)*XC(2)+XC(3)*XC(3)

C
        IF (ABS(V113).LE.EPS12) GOTO 98
        HELP=V1V23*V1V23/V113-V223
        IF (ABS(HELP).LE.EPS12) GOTO 98
        XK3=V1V23/V113
        XN3=1./HELP
        PS13(1,J)=(XB(1)*XK3-XC(1))*XN3
        PS13(2,J)=(XB(2)*XK3-XC(2))*XN3
        PS13(3,J)=(XB(3)*XK3-XC(3))*XN3
        P1A(J)=-P3(1,J)*PS13(1,J)-P3(2,J)*PS13(2,J)-P3(3,J)*PS13(3,J)
C
        IF (ABS(V223).LE.EPS12) GOTO 98
        HELP=V1V23*V1V23/V223-V113
        IF (ABS(HELP).LE.EPS12) GOTO 98
        XK3=V1V23/V223
        XN3=1./HELP
        PS23(1,J)=(XC(1)*XK3-XB(1))*XN3
        PS23(2,J)=(XC(2)*XK3-XB(2))*XN3
        PS23(3,J)=(XC(3)*XK3-XB(3))*XN3
        P2A(J)=-P3(1,J)*PS23(1,J)-P3(2,J)*PS23(2,J)-P3(3,J)*PS23(3,J)
C
        IF (RLB(J).GE.4) THEN
          V2V34=XD(1)*XE(1)+XD(2)*XE(2)+XD(3)*XE(3)
          V224=XD(1)*XD(1)+XD(2)*XD(2)+XD(3)*XD(3)
          V334=XE(1)*XE(1)+XE(2)*XE(2)+XE(3)*XE(3)
C
          IF (ABS(V224).LE.EPS12) GOTO 98
          HELP=V2V34*V2V34/V224-V334
          IF (ABS(HELP).LE.EPS12) GOTO 98
          XK4=V2V34/V224
          XN4=1./HELP
          PS24(1,J)=(XD(1)*XK4-XE(1))*XN4
          PS24(2,J)=(XD(2)*XK4-XE(2))*XN4
          PS24(3,J)=(XD(3)*XK4-XE(3))*XN4
          P1B(J)=-P4(1,J)*PS24(1,J)-P4(2,J)*PS24(2,J)-P4(3,J)*PS24(3,J)
C
          IF (ABS(V334).LE.EPS12) GOTO 98
          HELP=V2V34*V2V34/V334-V224
          IF (ABS(HELP).LE.EPS12) GOTO 98
          XK4=V2V34/V334
          XN4=1./HELP
          PS34(1,J)=(XE(1)*XK4-XD(1))*XN4
          PS34(2,J)=(XE(2)*XK4-XD(2))*XN4
          PS34(3,J)=(XE(3)*XK4-XD(3))*XN4
          P2B(J)=-P4(1,J)*PS34(1,J)-P4(2,J)*PS34(2,J)-P4(3,J)*PS34(3,J)
C
          IF (RLB(J).GE.5) THEN
            V3V45=XF(1)*XG(1)+XF(2)*XG(2)+XF(3)*XG(3)
            V335=XF(1)*XF(1)+XF(2)*XF(2)+XF(3)*XF(3)
            V445=XG(1)*XG(1)+XG(2)*XG(2)+XG(3)*XG(3)
C
            IF (ABS(V335).LE.EPS12) GOTO 98
            HELP=V3V45*V3V45/V335-V445
            IF (ABS(HELP).LE.EPS12) GOTO 98
            XK5=V3V45/V335
            XN5=1./HELP
            PS35(1,J)=(XF(1)*XK5-XG(1))*XN5
            PS35(2,J)=(XF(2)*XK5-XG(2))*XN5
            PS35(3,J)=(XF(3)*XK5-XG(3))*XN5
            P1C(J)=-P5(1,J)*PS35(1,J)-P5(2,J)*PS35(2,J)-
     -              P5(3,J)*PS35(3,J)
C
            IF (ABS(V445).LE.EPS12) GOTO 98
            HELP=V3V45*V3V45/V445-V335
            IF (ABS(HELP).LE.EPS12) GOTO 98
            XK5=V3V45/V445
            XN5=1./HELP
            PS45(1,J)=(XG(1)*XK5-XF(1))*XN5
            PS45(2,J)=(XG(2)*XK5-XF(2))*XN5
            PS45(3,J)=(XG(3)*XK5-XF(3))*XN5
            P2C(J)=-P5(1,J)*PS45(1,J)-P5(2,J)*PS45(2,J)-
     -              P5(3,J)*PS45(3,J)
          ENDIF
        ENDIF
C
        IF (ABS(TEST4).GT.EPS10) THEN
          WRITE (iunout,*) 'WARNING FROM TIMEA0, TEST FOR 4.TH POINT:'
          WRITE (iunout,*) 'SF. NUMBER, TEST= ',J,TEST4
        ELSEIF (ABS(TEST5).GT.EPS10) THEN
          WRITE (iunout,*) 'WARNING FROM TIMEA0, TEST FOR 5.TH POINT:'
          WRITE (iunout,*) 'SF. NUMBER, TEST= ',J,TEST5
        ENDIF
C
C  SET SOME MORE ASSISTENT SURFACE DATA AND CHECK CONSISTENCY
C  ST. NO. 90 --- 99
C
   90   CONTINUE
C
        RLBNOT(J)=RLB(J).EQ.1.5.OR.RLB(J).EQ.2.5.OR.RLB(J).EQ.3.5.OR.
     .            RLB(J).EQ.4.5.OR.RLB(J).EQ.5.5
C
        IF (ILSWCH(J).NE.0.AND.ILIIN(J).EQ.0) THEN
          WRITE (iunout,*) 'EXIT FROM TIMEA0: SURFACE NO J IS OPERATING'
          WRITE (iunout,*) 'A SWITCH BUT IS NOT SEEN BY HISTORY'
          WRITE (iunout,*) 'J= ',J
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
C
        IF (ILSWCH(J).NE.0.AND.ILIIN(J).GT.0) THEN
          DO ISPZ=1,NSPTOT
            IF (TRANSP(ISPZ,1,J).NE.0.D0.OR.TRANSP(ISPZ,2,J).NE.0.D0)
     .      GOTO 95
          ENDDO
          GOTO 96
   95     WRITE (iunout,*) 'EXIT FROM TIMEA0: SURFACE NO J IS OPERATING'
          WRITE (iunout,*) 'A SWITCH BUT IS SOMETIMES TRANSPARENT AND '
          WRITE (iunout,*)
     .      'SOMETIMES REFLECTING (SEMI-TRANSPARENCY OPTION)'
          WRITE (iunout,*)
     .      'POSSIBLE FIXES: MANUAL, CHAPTER 2, SECTION 6 '
          WRITE (iunout,*) 'J= ',J
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
C
   96   GOTO 94
C
   98   WRITE (iunout,*)
     .    'WARNING FROM TIMEA0, SURFACE NO J IS TURNED OFF'
        WRITE (iunout,*) 'BECAUSE ILL-DEFINED '
        WRITE (iunout,*) 'J= ',J
        IGJUM0(J)=1
        IF (NLIMPB >= NLIMPS) THEN
          DO IAB=0,NLIMI
            IGJUM1(IAB,J)=1
          END DO
        ELSE
          DO IAB=0,NLIMI
            CALL EIRENE_BITSET (IGJUM1,0,NLIMPS,IAB,J,1,NBITS)
          END DO
        END IF
        GOTO 97
C
   94   CONTINUE
C
        JUMLIM(J)=0
        ALM(J)=2.*A4LM(J)
        BLM(J)=2.*A5LM(J)
        CLM(J)=2.*A6LM(J)
c       if first order surface:
        IF (A4LM(J).EQ.0..AND.A5LM(J).EQ.0..AND.A6LM(J).EQ.0..AND.
     .      A7LM(J).EQ.0..AND.A8LM(J).EQ.0..AND.A9LM(J).EQ.0.) THEN
          IF (NLIMPB >= NLIMPS) THEN
            IGJUM1(J,J)=1
          ELSE
            CALL EIRENE_BITSET (IGJUM1,0,NLIMPS,J,J,1,NBITS)
          END IF
          AT=MAX(ABS(A1LM(J)),ABS(A2LM(J)),ABS(A3LM(J)))
          IF (ABS(A1LM(J)).EQ.AT) JUMLIM(J)=1
          IF (ABS(A2LM(J)).EQ.AT) JUMLIM(J)=2
          IF (ABS(A3LM(J)).EQ.AT) JUMLIM(J)=3

c  create HNF coefficients (normalized vector and distance)
          XNORM=SQRT(A1LM(J)*A1LM(J)+A2LM(J)*A2LM(J)+A3LM(J)*A3LM(J))
          A0LM(J)=A0LM(J)/XNORM
          A1LM(J)=A1LM(J)/XNORM
          A2LM(J)=A2LM(J)/XNORM
          A3LM(J)=A3LM(J)/XNORM

          JUM=JUMLIM(J)
          SELECT CASE (JUM)
          CASE (1)
            ALM(J)=-A0LM(J)/A1LM(J)
            BLM(J)=-A2LM(J)/A1LM(J)
            CLM(J)=-A3LM(J)/A1LM(J)
          CASE (2)
            ALM(J)=-A0LM(J)/A2LM(J)
            BLM(J)=-A1LM(J)/A2LM(J)
            CLM(J)=-A3LM(J)/A2LM(J)
          CASE (3)
            ALM(J)=-A0LM(J)/A3LM(J)
            BLM(J)=-A1LM(J)/A3LM(J)
            CLM(J)=-A2LM(J)/A3LM(J)
          END SELECT
        ENDIF
   97 CONTINUE

c.....................................................

!pb 04122019
!pb   clean up in case NLIMI changed from a prior iteration      
      call eirene_dealloc_timea
      IF (NLOCTREE) then
c  calling an internal subroutine for preparing octree procedure:
c  set flag NOTOCSURFS of those additional surfaces, which are not
c  handled by octree optimization.
c
        CALL EIRENE_TIMEA0_BUILDOC()

        CALL EIRENE_LEER(2)

      ELSE  !   NOT NLOCTREE
c  old default: ALL surfaces are outside the octree procedure,
c               i.e. OCTREE optimization in in-active
        if (.not.allocated(NOTOCSURFS)) ALLOCATE(NOTOCSURFS(NLIMI))
        DO J=1,NLIMI
          NOTOCSURFS(J) = J
        END DO
        NSURFNOT = NLIMI
      ENDIF


      END SUBROUTINE EIRENE_TIMEA0

c  write the octree to vtk tool output
      SUBROUTINE EIRENE_TIMEA3_VTKOUT
      IMPLICIT NONE
        if (associated(tree) .and. NLOCTREE) then
c  call octree function to plot the tree recursive
          call OCTREE_PrintVTK(tree, IVTKOUT)

          open(200,file='graphviz.out')
          call OCTREE_PrintGraphviz(tree, 200)
          close(200)
        end if
      END SUBROUTINE EIRENE_TIMEA3_VTKOUT

C
      SUBROUTINE EIRENE_TIMEA1
     . (MSURF,NCELL,NLI,NLE,NTCELL,IPERID,XX,YY,ZZ,TMT,
     .  VXX,VYY,VZZ,VV,MASURF,XR,YR,ZR,SG,TL,NLTRC,LCNDEXP)
      IMPLICIT NONE

      REAL(DP) :: TMIN, TL
      REAL(DP) :: XX, YY, ZZ, TMT, VXX, VYY, VZZ, VV, X_S, Y_S, Z_S,
     .          SG_S, XR, YR, ZR, SG
     .          
      INTEGER MSURF, NNTCL, NTCELL, IPERID, MASURF, NCELL, 
     .        NLI, NLE, MASURF_S
      LOGICAL, INTENT(OUT) :: LCNDEXP
      LOGICAL :: NLTRC, LCNDEXP_S
c     some vars for checking on intersection/where we are
      TYPE(ocnode), POINTER :: block
      REAL(DP), DIMENSION(3) :: ip, start, direction
      REAL(DP) :: runlength, norm
      LOGICAL :: status


      SAVE 

      INTERFACE
        SUBROUTINE EIRENE_TIMEA_CheckInter(
     .             MSURF,NCELL,NLI,NLE,NNTCL,XX,YY,ZZ,TMT,
     .             VXX,VYY,VZZ,VV,
     .             MASURF,XR,YR,ZR,SG,TL,NLTRC,LCNDEXP,
     .             ACTSURFS, NACTSURFS)
          USE EIRMOD_PRECISION
          INTEGER, INTENT(in) :: MSURF,NCELL,NLI,NLE,NNTCL,NACTSURFS
          INTEGER, DIMENSION(:), INTENT(in) :: ACTSURFS
          INTEGER, INTENT(out) :: MASURF
          LOGICAL, INTENT(out) :: LCNDEXP
          LOGICAL, INTENT(in) :: NLTRC
          REAL(DP), INTENT(in) :: XX, YY, ZZ, VXX, VYY, VZZ, VV, TMT
          REAL(DP), INTENT(out) :: XR, YR, ZR, SG, TL
        END SUBROUTINE
      END INTERFACE

!pb 08.01.14
      LCNDEXP = .FALSE.

C  PARTICLE ON STANDARD SURFACE?
      IF (MSURF.GT.NLIM) MSURF=0
C  FIND LOCAL COORDINATE SYSTEM IN CASE OF NLTRA
      NNTCL=1
      IF (NLTRA.AND.NLTOR) THEN
        NNTCL=NTCELL
      ELSEIF (NLTRA.AND..NOT.NLTOR) THEN
        NNTCL=IPERID
      ENDIF


#ifdef TRACE
      IF (NLTRC) THEN
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) 'TIMEA ,X,Y,Z,T ',XX,YY,ZZ,TMT
        WRITE (iunout,*) '      VX,VY,VZ,V ',VXX,VYY,VZZ,VV
        IF (NLTRA) WRITE (iunout,*) 'MSURF,NNTCL ',MSURF,NNTCL
        IF (.NOT.NLTRA) WRITE (iunout,*) 'MSURF ',MSURF
      ENDIF
#endif

!trc  if(msurf .eq. 0 .and. pladd) then
!trc    WRITE(trcnum,*) 'RESTART', XX, YY, ZZ, VXX, VYY, VZZ
!trc  end if

      TMIN=1.D30
      TL=1.D30
      MASURF=0

      if (NLOCTREE .and. trcoct) then
        WRITE(iunout,*)
        WRITE(iunout,*) "PROCESSING ADDITIONAL SURFACES WITH RLB < 3"
      endif

c   call timea2 for all those SURFACES which ARE no plane triangle
c   -> we will get a time minimum out of this, if any of these
c   surfaces are hit.
      CALL EIRENE_TIMEA_CheckInter(MSURF,NCELL,NLI,NLE,NNTCL,
     .                             XX,YY,ZZ,TMT,VXX,VYY,VZZ,VV,
     .                             MASURF_S,X_S,Y_S,Z_S,SG_S,TMIN,
     .                             NLTRC,LCNDEXP_S,
     .                             NOTOCSURFS, NSURFNOT)

c   if we actually found a nearest valid intersection on second order surfs,
c   etc, save these values for later comparison with octree surfs values
      if (masurf_s .gt. 0) then
        if (NLOCTREE .and. trcoct) then
          WRITE(iunout,*) "found intersection with surface", masurf_s
          WRITE(iunout,*) "-> continuing with this candidate in octree"
        end if
        XR = X_S
        YR = Y_S
        ZR = Z_S
        TL = TMIN
        SG = SG_S
        MASURF = MASURF_S
        LCNDEXP = LCNDEXP_S
      end if

c  if we do not have a tree, skip this part, all done, return to calling program.
      if(.not. associated(tree) .or. .not. NLOCTREE) then
        return
      end if

      if (trcoct) then
        WRITE(iunout,*)
        WRITE(iunout,*) "PROCESSING ADDITIONAL SURFACES WITH RLB >= 3"
      end if

c  set working coords to initial coords (xx,yy,zz)
      start = (/XX, YY, ZZ/)
      direction = (/VXX, VYY, VZZ/)
c  give status an initial "true" - if we do not hit the block
c  (and if outside), it will be set to false
      status = .TRUE.

c  we want to find out which tri/quad/quintangles are possible candidates
c  for intersection, by analysing the octree.

c  first: are we within octree or outside of it?
c     -> transform point to octree space
      if (.not. OCTREE_CheckVolume(start, tree%root, .true.)) then
c  second: if we are outside, check if we could ever hit our octree-block
        call OCTREE_CheckBlock(start, direction,
     .                         tree%root, status, ip, runlength)
        if (status) then
c  TODO: make a check if we hit the space within proper time... -> PETRA?
c       -> for now presume that this is the case
c  third: get a new starting point by moving the start to the
c         intersection point with the block (then we are inside)
          if (trcoct) then
            WRITE(iunout,*) "moving ray into octree space first:"
            WRITE(iunout,*) "start: ", start
            WRITE(iunout,*) "new start: ", ip
          end if
          start = ip
        end if
      end if

c     if we are inside the block (or if we moved there as we hit the
c     convex hull while following the path...), continue with octree
c     processing. if we are not within (even with the check if we intersect),
c     just continue with the other add. surfaces not in our octree...
      if(status) then
        if(trcoct) then
          WRITE(iunout,*) "starting trace @",start, "in direction",
     .                    (/VXX, VYY, VZZ/)
        end if
c       calc the norm of the direction vector, we need this for traversal
        norm = sqrt(dot_product(direction,direction))
      end if
c     step through the octree until we have found an intersection
c     or we leave the octree space
      do while(status)
c       now find out where the hell we are in the octree space...
c       -> get the pointer to our leaf-block containing the IP
        block => OCTREE_GetLeafchild(start, tree)
        if (trcoct) then
          WRITE(iunout,*) "searching in block:",block%number,
     .                    " on layer", block%layer, " testing",
     .                    block%nsurfaces, " surfaces"
        end if

c       if we have no surfaces to check in this block, continue to next
        if(block%nsurfaces .gt. 0) then
c         now get all surfaces in this block and search with these
c         in CheckInter
          CALL EIRENE_TIMEA_CheckInter(MSURF,NCELL,NLI,NLE,NNTCL,
     .                             XX,YY,ZZ,TMT,VXX,VYY,VZZ,VV,
     .                             MASURF_S,X_S,Y_S,Z_S,SG_S,TMIN,
     .                             NLTRC,LCNDEXP_S,
     .                             block%surfaces, block%nsurfaces)

c         if we found a valid intersection and the time is less than TL
c         save these values for later comparison with octree surfs values
          if (masurf_s .gt. 0 .and. tmin .lt. TL) then
c            WRITE(iunout,*) "found intersection with surface", masurf_s
            XR = X_S
            YR = Y_S
            ZR = Z_S
            IP = (/XR, YR, ZR/)
            TL = TMIN
            SG = SG_S
            MASURF = MASURF_S
            LCNDEXP = LCNDEXP_S
c           if we got a intersection WITHIN our octet, we can stop here
c           (there is no one with a smaller time reachable) else continue...
            if(OCTREE_CheckVolume(ip, block, .true.)) exit
          end if
        end if

c       if we did not end the search before, we need to traverse
c       through this block and start over in the neighbor.
c       -> if we get out of octree space, start%p will be =-1
        start = OCTREE_Traverse(tree, block, start,
     .                          direction, norm)
c       check if we are inside of the octree space anymore...
        status = OCTREE_CheckVolume(start, tree%root, .true.)
        if(.not.status.and.trcoct) WRITE(iunout,*)'left octree space...'
      end do

!trc      if(MASURF .gt. 0 .and. pladd) then
!trc        WRITE(TRCNUM,*) MASURF, XX, YY, ZZ, XR, YR, ZR, VXX, VYY, VZZ
!trc      end if

      END SUBROUTINE EIRENE_TIMEA1

c     subroutine to call the octree building internal subroutine
      SUBROUTINE EIRENE_TIMEA0_OC
      IMPLICIT NONE
        CALL EIRENE_TIMEA0_BUILDOC()
      RETURN
      END SUBROUTINE EIRENE_TIMEA0_OC

c................................................................

c     defining the internal subroutine for the octree build

      SUBROUTINE EIRENE_TIMEA0_BUILDOC()
      IMPLICIT NONE
      INTEGER J, TRCNUM

      SAVE

      INTERFACE
        FUNCTION EIRENE_TIMEA_BuildOctree() RESULT(tree)
          USE EIRMOD_OCTREE
          TYPE(octree), POINTER :: tree
        END FUNCTION
      END INTERFACE

cdr  this routine is only called if NLOCTREE
c
c    It returns NOTOCSURFS(1:NSURFNOT), a list of additional surfaces not
c       covered by the octree-procedure
c
c    unit number for the tracing debug ouput
      if (pladd) then
        TRCNUM = 27
!trc    OPEN(unit=TRCNUM,file='timea-trace.out')
      end if

C   as the code above is pretty much with gotos, simply put
C   the octree stuff in here...
c
c   build octree
      tree => EIRENE_TIMEA_BuildOctree()

c   count the surfaces that will not be in octree to get the
c   right dimension for the array NOTOCSURFS
c
      DO J=1,NLIMI
c  only add surfaces that are valid and only surfaces
c  that are no tri/quad/quintangles
        IF(RLB(J) .lt. 3 .and. IGJUM0(J) .eq. 0) THEN
           NSURFNOT = NSURFNOT +1
        END IF
      END DO

c   allocate enough space for lookup table which we will use for
c   surfaces that we do not have in the octree (like 2nd order surfaces, etc.)
      if (.not.allocated(NOTOCSURFS)) ALLOCATE(NOTOCSURFS(NSURFNOT))
cdr
      NSURFNOT = 0
c   add the second order surfaces, etc into this
      DO J=1,NLIMI
c  only add surfaces that are valid and only surfaces
c  that are no tri/quad/quintangles
        IF(RLB(J) .lt. 3 .and. IGJUM0(J) .eq. 0) THEN
           NSURFNOT = NSURFNOT +1
           NOTOCSURFS(NSURFNOT) = J
        END IF
      END DO


      END SUBROUTINE EIRENE_TIMEA0_BUILDOC

      
      SUBROUTINE EIRENE_DEALLOC_TIMEA

      if (ALLOCATED(NOTOCSURFS)) DEALLOCATE (NOTOCSURFS)
      if (ASSOCIATED(tree)) call OCTREE_DeleteTree(tree)
      
      END SUBROUTINE EIRENE_DEALLOC_TIMEA
      
      END MODULE EIRMOD_TIMEA
