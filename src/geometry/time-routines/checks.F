cdr  Jan 2020: Below here is also hidden the former timea1.f, which is
cdr                 always called, not only in case of octree.f
cdr            Adding and correcting comments.
cdr   tbd:
cdr            Identify: ACTSURFS, NACTSURFS, when no octree is present?
cdr                      and how they may be related to
cdr                      lgjum1,2,3(..).
cdr

cdr  aug. 20: code cleanup from ITER branch
c*********************************************
c     This file contains routines for checking on intersections
c     between blocks and planes (triangles, quadrangles, ...)
c     It is only used by octree.f.   cdr:  No: the timea_CheckInter(...) below
c                                    cdr   is always used !
c     -------------------------------
c     Author: Oliver Schmidt <o.schmidt@fz-juelich.de>
c     Revision: 0.1
c     Date: 25.11.2011
c     -------------------------------
c*********************************************

c     --- EIRENE_TIMEA_FirstCheck ---
c     first check: is at least one point of this surface in the block?
c     sid    : index number of the surface to test against block
c     child  : pointer on the child node = block to test against
c     RETURNS: logical result (true/false) if inside or not
      FUNCTION EIRENE_TIMEA_FirstCheck(SID, CHILD) RESULT(STATUS)
        USE EIRMOD_PRECISION
        USE EIRMOD_CADGEO
        USE EIRMOD_OCTREE
        USE EIRMOD_COMPRT
        IMPLICIT NONE

        TYPE(ocNode), POINTER, INTENT(IN) :: CHILD
        INTEGER, INTENT(IN) :: SID
        LOGICAL :: P1IN, P2IN, P3IN, P4IN, P5IN, STATUS

c       test if the points are within the block by checking
c       if their components in xyz are between the edges
        P1IN = OCTREE_CheckVolume(P1(:,SID), CHILD, .false.)
        P2IN = OCTREE_CheckVolume(P2(:,SID), CHILD, .false.)
        P3IN = OCTREE_CheckVolume(P3(:,SID), CHILD, .false.)
        P4IN = OCTREE_CheckVolume(P4(:,SID), CHILD, .false.)
        P5IN = OCTREE_CheckVolume(P5(:,SID), CHILD, .false.)

c       only one has to be inside, so we use or here
        STATUS = P1IN .or. P2IN .or. P3IN .or. P4IN .or. P5IN
      END FUNCTION EIRENE_TIMEA_FirstCheck

c     --- EIRENE_TIMEA_ThirdCheck ---
c     third check: make intersection tests with block and surface
c     sid    : index number of the surface to test against block
c     child  : pointer on the child node = block to test against
c     RETURNS: logical result (true/false) if surface hits block or not
      FUNCTION EIRENE_TIMEA_ThirdCheck(SID, CHILD) RESULT(STATUS)
        USE EIRMOD_PRECISION
        USE EIRMOD_COMPRT
        USE EIRMOD_CADGEO
        USE EIRMOD_OCTREE
        USE EIRMOD_CTRCEI
        IMPLICIT NONE

        TYPE(ocNode), POINTER, INTENT(IN) :: CHILD
        INTEGER, INTENT(IN) :: SID
        INTEGER :: I, T, IND
        REAL(DP), DIMENSION(3) :: A, B, C, O, d, RES, ip, h
        LOGICAL :: STATUS
        REAL(DP) :: S

c       lookup tables for the block-edge to line conversion in the
c       checks performed in TIMEA
c       (i,1) = x, (i,2) = y, (i,3) = z, only 4 needed (reuse with function)
        INTEGER, DIMENSION(3,4), SAVE :: POINTS = reshape(
     .                                      (/1, 1, 1,
     .                                        2, 2, 1,
     .                                        2, 1, 2,
     .                                        1, 2, 2/),
     .                                      (/3, 4/))
c       saving the diagonal entries for the 4 points.
c       every point has 3 vectors for the edge direction.
c       to generate these we use the matrix
c       x 0 0
c       0 x 0
c       0 0 x
c       where we have the x as the number below, returning for every block
        INTEGER, DIMENSION(12), SAVE :: NORMVEC = (/1,  1,  1,
     .                                             -1, -1,  1,
     .                                             -1,  1, -1,
     .                                              1, -1, -1/)

c       make sure we are looking at a tria/quadr/quintangle
        IF (RLB(SID).ge.3) THEN
c         for every triangle we may get (triangle = 1, quadrangle = 2, quintangle = 3)
c         make a test run! (wow... that is much work to do...)
          DO T=1,nint(RLB(SID)-2.)
c           get us these points of the triangle we look at
            SELECT CASE(T)
              CASE(1)
                A=P1(:,SID)
                B=P2(:,SID)
                C=P3(:,SID)
              CASE(2)
                A=P2(:,SID)
                B=P3(:,SID)
                C=P4(:,SID)
              CASE(3)
                A=P3(:,SID)
                B=P4(:,SID)
                C=P5(:,SID)
            END SELECT
c           a) check if any of the 12 edges hits the surface (check plane
c              and if inside the triangles, if so, intersection point in block?)
c           iterate over the 12 edges, using I=0..11 to be able to use modulo
            DO I=0,11
c             define the aufpunkt of the line (this cannot be written
c             in compact form as we need the entries of POINTS and take
c             a few things out of matrix B, this cannot be done in range
c             mode)
              if (i < 3) then
                ind = 1
              else if (i < 6) then
                ind = 2
              else if (i < 9) then
                ind = 3
              else
                ind = 4
              end if
              O(1) = child%B(1,POINTS(1,ind))
              O(2) = child%B(2,POINTS(2,ind))
              O(3) = child%B(3,POINTS(3,ind))

c             define the direction vector
              d = (/0, 0, 0/)
c             replace the correct element with the corresponding entry
c             to get the normal vector of the edge we check
              d(mod(I,3)+1) = NORMVEC(I+1)

c             solve the LES via cramer
              RES = OCTREE_Cramer(A, B, C, O, d)
c             check if there is a intersection...
c             -> check first via barycentric coordinates algorithm (see also cramer.f)
c             -> the RES(1) is always 0 if cramer unsolvable or
c                if beta<=0 / gamma<=0 / t<=0 / beta+gamma>=1
              IF (RES(1) .gt. 0) THEN
c               we are here, so now test if the intersection point
c               we just found is inside the block (right now it is only
c               inside the triangle!)
c               -> calculate the intersection point by inserting t(RES(3))
c                  in the line equation
                ip = O+RES(3)*d

c               if the intersection point is within the block volume,
c               return true
                STATUS = OCTREE_CheckVolume(IP, CHILD, .TRUE.)
                IF (STATUS) RETURN
              END IF
            END DO

c           b) check if the edges of the triangle go through the block
c           -> check them as lines in R3, because they cannot stop
c              in the middle of the block (P1,P2,P3,... is not in the block!)
c           -> this check is the same as the check for a particle trace,
c              use this as simplification of code here
c           -> ONLY use this if abs(S) <= 1 ==> only use the edge one time,
c               do not make it longer as it actually is!
c           line from A to B:
            H = B-A
            CALL OCTREE_CheckBlock(A, H, CHILD, STATUS, IP, S)
            IF (STATUS .and. abs(S) .le. 1) THEN
              RETURN
            END IF
c           line from B to C
            H = C-B
            CALL OCTREE_CheckBlock(B, H, CHILD, STATUS, IP, S)
            IF (STATUS .and. abs(S) .le. 1) THEN
              RETURN
            END IF
c           line from C to A
            H = A-C
            CALL OCTREE_CheckBlock(C, H, CHILD, STATUS, IP, S)
            IF (STATUS .and. abs(S) .le. 1) THEN
              RETURN
            END IF
          END DO
        END IF
c       we found nothing, so we have no intersection -> return false
        STATUS = .FALSE.
      END FUNCTION EIRENE_TIMEA_ThirdCheck

cdr  Here starts the original timea1.f: check additional surfaces.
cdr  It is now called: EIRENE_TIMEA_CheckInter
cdr  For comparison with older versions: then compare with eirene_timea1.
cdr  this part is always used, not only in case nloctree

c     --- CHECK FOR INTERSECTIONS ON SURFACES (FORMER EIRENE_TIMEA1) ---
c
C         THE RAY STARTS AT "XX,YY,ZZ", IN DIRECTION "VXX,VYY,VZZ"
C                 STARTING POINT IS IN CELL "NCELL", AND TOROIDAL PERIODICITY
C                 SEGMENT "IPERID, NTCELL".
C                 SEARCH FOR LEGAL INTERSECTIONS WITH ADDITIONAL SURFACES
C                 MASURF, IN THE RANGE "NLI LE. MASURF LE. NLE". (DO LOOP: 100)
C                 IN DO 100 LOOP: TMIN IS THE TIME TO INTERSECTION ALREADY FOUND
C                                 FOR SURFACES ALREADY CHECKED.
C                                 LATER SURFACES ARE ONLY CHECKED IF
C                                 SMALLER TIMES TL ARE POSSIBLE. OTHERWISE
C                                 GOTO 100, ALREADY WITHOUT FULL EVALUATION.
C                                 HENCE: TRY TO FIND OUT FOR EACH SURFACE
C                                        IF THERE IS A CANDIDATE INTERSECTION
C                                        WITH TL LT TMIN. OTHERWISE: NEXT SURFACE.
!    NEW OPTION SINCE 2005: IF NLPRCS(J) FOR SURFACE "J": CARRY OUT FULL CHECK,
!                           EVEN IF TL > TMIN, TO IDENTIFY POSSIBLE LEGAL INTERSECTION.
!                           THIS INTERSECTION IS NOT NECESSARILY RETURNED AS VALID
!                           INTERSECTION WITH SMALLEST DISTANCE,
!                           BUT THE FLAG LCNDEXP=TRUE, INDICATING
!                           A POSSIBLE VALID INTERSECTION WITH AT LEAST ONE OF THE SURFACES J
!                           FOR WHICH NLPRCS(J)=.TRUE.
!                           THIS OPTION IS USED FOR CONDITIONAL EXPECTATION ESTIMATORS.
C                 RETURN THE SURFACE NUMBER "MASURF" WITH CLOSEST LEGAL INTERSECTION
C                 RETURN MASURF=0 IF NO LEGAL INTERSECTION WAS FOUND
C                 IF MASURF.GT.0: ALSO RETURN POINT OF INTERSECTION XR,YR,ZR,
C                 THE ORIENTATION "SG" OF RAY RELATIVE TO SURFACE NORMAL
C                 AND THE DISTANCE "TL" ALONG RAY UNTIL INTERSECTION POINT.

      SUBROUTINE EIRENE_TIMEA_CheckInter(
     .             MSURF,NCELL,NLI,NLE,NNTCL,XX,YY,ZZ,TMT,
     .             VXX,VYY,VZZ,VV,
     .             MASURF,XR,YR,ZR,SG,TL,NLTRC,LCNDEXP,
     .             ACTSURFS, NACTSURFS)
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
#ifdef TRACE
        USE EIRMOD_COMPRT, ONLY: IUNOUT
#endif
        IMPLICIT NONE

        INTERFACE
         LOGICAL FUNCTION EIRENE_BITGET(IBIT,N1LOW,N1UP,IROW,JCOL,NBITS)
            INTEGER, INTENT(IN) :: N1LOW, N1UP, IROW, JCOL, NBITS
            INTEGER, INTENT(IN) :: IBIT(N1LOW:N1UP,*)
          END FUNCTION

          SUBROUTINE EIRENE_FZRTOR(X,Z,NOLD,XNEW,PH,NNEW,LTEST,NTEST)
            USE EIRMOD_PRECISION
            REAL(DP), INTENT(IN) :: X, Z
            REAL(DP), INTENT(OUT) :: XNEW, PH
            INTEGER, INTENT(IN) :: NOLD, NTEST
            INTEGER, INTENT(OUT) :: NNEW
            LOGICAL, INTENT(IN) :: LTEST
          END SUBROUTINE
        END INTERFACE

        INTEGER, INTENT(in) :: MSURF, NCELL, NLI, NLE, NNTCL, NACTSURFS
        INTEGER, DIMENSION(:), INTENT(in) :: ACTSURFS
        INTEGER, INTENT(out) :: MASURF
        INTEGER :: NLLLI, NNTCLS, NN, I, J, K, ICOUNT, NNR,
     .             NTNEW

        LOGICAL, INTENT(out) :: LCNDEXP
        LOGICAL, INTENT(in) :: NLTRC
        LOGICAL :: LGJ, LTSTCXP
c        LOGICAL :: LMTSRF(NLIMPS)

        REAL(DP), INTENT(in) :: XX, YY, ZZ, VXX, VYY, VZZ, VV, TMT
        REAL(DP), INTENT(out) :: XR, YR, ZR, SG, TL
        REAL(DP) :: XS, YS, ZS, VXS, VYS, VZS, VVS, TS,
     .              X,  Y,  Z,  VX,  VY,  VZ,  V, T, TADD, TMIN,
     .              VSAVE, VXJ, VZJ, WR, XN, YN, ZN,
     .              XLS1, XLS2, XLS3, XMS1, XMS2, XMS3,
     .              A1, A2, A3, F, G, TMI, TMA, TMX, ROT, TST, TUP,
     .              XXR, PPP

c      LMTSRF=.FALSE.

C  TENTATIVELY ASSUME: NO LEGAL INTERSECTION WITH ANY OF THE
C                      SURFACES J FOR WHICH NLPRCS(J)=TRUE
      LCNDEXP=.FALSE.

C  SAVE INITIAL COORDINATES
      XS=XX
      YS=YY
      ZS=ZZ
      TS=TMT
      VXS=VXX
      VYS=VYY
      VZS=VZZ
      VVS=VV
      NNTCLS=NNTCL
C  WORKING COORDINATES
      X=XX
      Y=YY
      Z=ZZ
      T=TMT
      VX=VXX
      VY=VYY
      VZ=VZZ
      V=VV
      NN=NNTCL

      NLLLI=MSURF
!dr PARTICLE IS ON SURFACE THAT HAS BEEN SWITCHED OFF
      IF(IGJUM0(MSURF).NE.0) NLLLI=0     !VK!!!!!!!!!
      TADD=0.
C
 1000 TMIN=1.D30
      TL=1.D30
      MASURF=0
C
C  LOOP OVER SURFACE NUMBER, DO 100
C  --> loop only over the surfaces we have in our lookup table,
C      respect offset, if given
      DO 100 K=1,NACTSURFS

c       get J for the former processing out of the running index K
c       using an offset (we can jump in the array, if driver timea1 wants this)
        J = ACTSURFS(K)

C  FIRST ELIMINATE ALL SURFACES WHICH ARE KNOWN A PRIORI TO BE
C        NOT POSSIBLE CANDIATES  (IGJUM FLAGS)
C
        IF (IGJUM0(J).NE.0) GOTO 100
        IF (NLIMPB >= NLIMPS) THEN
          IF (IGJUM1(NLLLI,J) .NE. 0) GOTO 100
        ELSE
          IF (EIRENE_BITGET(IGJUM1,0,NLIMPS,NLLLI,J,NBITS)) GOTO 100
        END IF
        IF (NCELL.LE.NOPTIM) THEN
          IF (NLIMPB >= NLIMPS) THEN
            IF (IGJUM3(NCELL,J).NE.0) GOTO 100
          ELSE
            IF (EIRENE_BITGET(IGJUM3,0,NOPTIM,NCELL,J,NBITS)) GOTO 100
          END IF
        ENDIF
C
C  FOR NLTRA OPTION ONLY:
C  X,Z,VX AND VZ ARE GIVEN IN TOROIDAL CELL NN,
C  TRANSFORM COORDINATES FOR THIS TRACK FROM LOCAL SYSTEM NN
C  TO THE LOCAL SYSTEM ILTOR(J), IN WHICH SURFACE J IS GIVEN
C  IF (ILTOR(J).LE.0) THIS SURFACE HAS TOROIDAL SYMMETRY
C
        IF (NLTRA) THEN
          IF (ILTOR(J).GT.0.AND.ILTOR(J).NE.NN) THEN
            CALL EIRENE_FZRTOR (X,Z,NN,XXR,PPP,NTNEW,.FALSE.,0)
            CALL EIRENE_FZRTRI (X,Z,ILTOR(J),XXR,PPP,NTNEW)
            ROT=2.*(NN-ILTOR(J))*ALPHA
            VSAVE=VX
            VX=COS(ROT)*VSAVE-SIN(ROT)*VZ
            VZ=SIN(ROT)*VSAVE+COS(ROT)*VZ
            NN=ILTOR(J)
C#ifdef TRACE
C           WRITE (iunout,*) 'J,ILTOR(J),X,Z,VX,VZ ',
C    .                        J,ILTOR(J),X,Z,VX,VZ
C#endif
C  X,Z,VX AND VZ ARE NOW GIVEN IN CELL NN=ILTOR(J). SO ARE THE COEFFICIENTS
C  OF SURFACE NO. J. FIND INTERSECTION IN THIS LOCAL SYSTEM
          ELSEIF (ILTOR(J).EQ.0) THEN
C  TOROIDALLY SYMMETRIC SURFACE, SURFACE COEFFICIENTS ARE THE SAME
C  IN EACH TOROIDAL CELL, THUS ESPECIALLY IN CELL NN
            NN=NNTCLS
            X=XS
            Z=ZS
            VX=VXS
            VZ=VZS
C#ifdef TRACE
C           WRITE (iunout,*) 'J,ILTOR(J),X,Z,VX,VZ ',
C    .                        J,ILTOR(J),X,Z,VX,VZ
C#endif
          ENDIF
        ENDIF
C
C  FIND INTERSECTION TIME TMX WITH BOUNDARY NO. J
C
C  TENTATIVELY ASSUME: NO LEGAL INTERSECTION WITH SURFACE J WITH DISTANCE
C                      LARGER THAN TMIN
        LTSTCXP=.FALSE.
        A1=0.

C  FIRST ORDER SURFACE ?
        GOTO (60,63,66),JUMLIM(J)

C  NO. THIS IS A SECOND ORDER SURFACE

C  A1*TMX*TMX+A2*TMX+A3=0
        A1=(A4LM(J)*VX+A7LM(J)*VY+A8LM(J)*VZ)*VX+
     .     (A5LM(J)*VY+A9LM(J)*VZ)*VY+A6LM(J)*VZ*VZ
        A2=(A1LM(J)+ALM(J)*X)*VX+(A2LM(J)+BLM(J)*Y)*VY+
     .     (A3LM(J)+CLM(J)*Z)*VZ+
     .      A7LM(J)*(VX*Y+VY*X)+A8LM(J)*(VX*Z+VZ*X)+A9LM(J)*(VY*Z+VZ*Y)
C  A3 =  0. ?
        IF (NLIMPB >= NLIMPS) THEN
          IF (IGJUM2(NLLLI,J).NE.0) GOTO 40
        ELSE
          IF (EIRENE_BITGET(IGJUM2,0,NLIMPS,NLLLI,J,NBITS)) GOTO 40
        END IF
C  NO
        A3=A0LM(J)+(A1LM(J)+A4LM(J)*X+A7LM(J)*Y+A8LM(J)*Z)*X
     .            +(A2LM(J)+A5LM(J)*Y+A9LM(J)*Z)*Y
     .            +(A3LM(J)+A6LM(J)*Z)*Z
        IF (A1.EQ.0.) GOTO 50
        F=-A2/(A1+A1)
        G=F*F-A3/A1
        IF (G.LT.0.) GOTO 100
        G=SQRT(G)
        TMA=F+G
        TMI=F-G
#ifdef TRACE
        IF (NLTRC) WRITE (iunout,*) 'TIMEA, J,TMI,TMA ',J,TMI,TMA
#endif
        IF (TMA.LE.EPS12) GOTO 100
        IF (TMI.GT.TMIN) THEN
          IF (.NOT.NLPRCS(J)) THEN
!  NO SWITCHING SURFACE FOR CONDITIONAL EXPECTATION ESTIMATOR
            GOTO 100
          ELSE
!  SWITCHING SURFACE FOR CONDITIONAL EXPECTATION ESTIMATOR
!  SURFACE J IS A CANDIDATE, BUT NOT ONE WITH THE SMALLEST DISTANCE
            LTSTCXP=.TRUE.
          END IF
        END IF

!  AT THIS POINT: EITHER TMI.LE.TMIN, OR LTSTCXP=TRUE

        IF (RLB(J).LT.0.) GOTO 21
        IF (RLB(J).GT.0.) GOTO 31
C
C  RLB(J) .EQ. 0.  NO FURTHER CONSTRAINS FOR THIS SURFACE. STATEMENT 11---20
C
        IF (TMI.LE.EPS12) GOTO 12
        IF (LTSTCXP) THEN
          LCNDEXP =.TRUE.
          GOTO 100
        END IF
        WR=A2+A1*(TMI+TMI)
        XN=X+TMI*VX
        YN=Y+TMI*VY
        ZN=Z+TMI*VZ
        TMX=TMI
        GOTO 70
C
   12   IF (TMA.GT.TMIN) THEN
          IF (LTSTCXP) LCNDEXP =.TRUE.
          GOTO 100
        END IF
        WR=A2+A1*(TMA+TMA)
   16   XN=X+TMA*VX
        YN=Y+TMA*VY
        ZN=Z+TMA*VZ
   18   TMX=TMA
        GOTO 70
C
C  CHECK BOUNDARY INEQUALITIES OF SURFACE
C  LGJ=.TRUE.: INTERSECTION POINT IS STILL VALID
C  LGJ=.FALSE.: INTERSECTION POINT IS OUTSIDE THE SPECIFIED AREA
C
C  RLB(J) .LT. 0.  STATEMENT 21---30
C
   21   IF (TMI.LE.EPS12) GOTO 22
        WR=A2+A1*(TMI+TMI)
        XN=X+TMI*VX
        YN=Y+TMI*VY
        ZN=Z+TMI*VZ
        TUP=TMI
        ICOUNT=1
        GOTO 28
C
   22   IF (TMA.GT.TMIN) THEN
          IF (.NOT.NLPRCS(J)) THEN
            GOTO 100
          ELSE
            LTSTCXP=.TRUE.
          END IF
        END IF
        WR=A2+A1*(TMA+TMA)
   26   XN=X+TMA*VX
        YN=Y+TMA*VY
        ZN=Z+TMA*VZ
        TUP=TMA
        ICOUNT=2
C
   28   LGJ=.TRUE.
        IF (ILIN(J).GT.0) THEN
          I=0
   27     I=1+I
          TST=ALIMS(I,J)+XLIMS(I,J)*XN+YLIMS(I,J)*YN+ZLIMS(I,J)*ZN
          LGJ=TST.LE.0.
          IF (LGJ.AND.I.LT.ILIN(J)) GOTO 27
        ENDIF
        IF (LGJ.AND.ISCN(J).GT.0) THEN
          I=0
   29     I=1+I
          TST=ALIMS0(I,J)+
     .        XN*(XLIMS1(I,J)+XN*XLIMS2(I,J)+YN*XLIMS3(I,J))+
     .        YN*(YLIMS1(I,J)+YN*YLIMS2(I,J)+ZN*ZLIMS3(I,J))+
     .        ZN*(ZLIMS1(I,J)+ZN*ZLIMS2(I,J)+XN*YLIMS3(I,J))
          LGJ=TST.LE.0.
          IF (LGJ.AND.I.LT.ISCN(J)) GOTO 29
        ENDIF
        IF (NLPRCS(J).AND.LGJ) LCNDEXP=.TRUE.
        IF (LGJ) THEN
          IF (LTSTCXP) GOTO 100
          TMX=TUP
          GOTO 70
        ELSEIF (ICOUNT.EQ.1) THEN
          GOTO 22
        ENDIF
        GOTO 100
C
C   RLB(J) .GT. 0.  STATEMENT 31---40
C
   31   IF (TMI.LE.EPS12) GOTO 32
        WR=A2+A1*(TMI+TMI)
        XN=X+TMI*VX
        YN=Y+TMI*VY
        ZN=Z+TMI*VZ
        LGJ=XN.LE.XLIMS2(1,J).AND.XN.GE.XLIMS1(1,J).AND.
     .      YN.LE.YLIMS2(1,J).AND.YN.GE.YLIMS1(1,J).AND.
     .      ZN.LE.ZLIMS2(1,J).AND.ZN.GE.ZLIMS1(1,J)
        IF (RLBNOT(J)) LGJ=.NOT.LGJ
        IF (NLPRCS(J).AND.LGJ) LCNDEXP=.TRUE.

C  IF DISTANCE TOO LARGE: LOSE INTEREST IN THIS SURFACE
        IF (LTSTCXP) GOTO 100

        TMX=TMI
        IF (LGJ) GOTO 70
C
   32   IF (TMA.GT.TMIN) THEN
          IF (.NOT.NLPRCS(J)) THEN
            GOTO 100
          ELSE
            LTSTCXP=.TRUE.
          END IF
        END IF
        WR=A2+A1*(TMA+TMA)
   36   XN=X+TMA*VX
        YN=Y+TMA*VY
        ZN=Z+TMA*VZ
   38   CONTINUE
        IF (RLB(J).LT.2.) THEN
          LGJ=XN.LE.XLIMS2(1,J).AND.XN.GE.XLIMS1(1,J).AND.
     .        YN.LE.YLIMS2(1,J).AND.YN.GE.YLIMS1(1,J).AND.
     .        ZN.LE.ZLIMS2(1,J).AND.ZN.GE.ZLIMS1(1,J)
        ELSE
          XMS1=XN*PS13(1,J)+YN*PS13(2,J)+ZN*PS13(3,J)+P1A(J)
          XLS1=XN*PS23(1,J)+YN*PS23(2,J)+ZN*PS23(3,J)+P2A(J)
          LGJ=XMS1.GE.0..AND.XLS1.GE.0..AND.XMS1+XLS1.LE.1.
          IF (RLB(J).GE.4.AND..NOT.LGJ) THEN
            XMS2=XN*PS24(1,J)+YN*PS24(2,J)+ZN*PS24(3,J)+P1B(J)
            XLS2=XN*PS34(1,J)+YN*PS34(2,J)+ZN*PS34(3,J)+P2B(J)
            LGJ=XMS2.GE.0..AND.XLS2.GE.0..AND.XMS2+XLS2.LE.1.
            IF (RLB(J).GE.5.AND..NOT.LGJ) THEN
              XMS3=XN*PS35(1,J)+YN*PS35(2,J)+ZN*PS35(3,J)+P1C(J)
              XLS3=XN*PS45(1,J)+YN*PS45(2,J)+ZN*PS45(3,J)+P2C(J)
              LGJ=XMS3.GE.0..AND.XLS3.GE.0..AND.XMS3+XLS3.LE.1.
            ENDIF
          ENDIF
        ENDIF
        IF (RLBNOT(J)) LGJ=.NOT.LGJ
        IF (NLPRCS(J).AND.LGJ) LCNDEXP=.TRUE.

C  IF DISTANCE TOO LARGE: LOSE INTEREST IN THIS SURFACE
        IF (LTSTCXP) GOTO 100

        TMX=TMA
        IF (LGJ) GOTO 70
        GOTO 100
C
C   A1*TMX+A2=0
C
   40   TMA=-A2/A1
#ifdef TRACE
        IF (NLTRC)
     .    WRITE (iunout,*) 'TIMEA AT 40, J,TMA,A1,A2 ',J,TMA,A1,A2
#endif
        IF (TMA.LE.EPS12) GOTO 100
        IF (TMA.GT.TMIN) THEN
          IF (.NOT.NLPRCS(J)) THEN
            GOTO 100
          ELSE
            LTSTCXP=.TRUE.
          END IF
        END IF
        WR=-A2
        IF (RLB(J).LT.0.) THEN
          GOTO 26
        ELSE IF (RLB(J).GT.0.) THEN
          GOTO 36
        ELSE
          GOTO 16
        END IF
C
C   A2*TMX+A3=0
C
   50   IF (A2.EQ.0.) GOTO 100
        TMA=-A3/A2
#ifdef TRACE
        IF (NLTRC) WRITE (iunout,*) 'TIMEA AT 50, J,TMA,A2,A3 ',
     .                                            J,TMA,A2,A3
#endif
        IF (TMA.LE.EPS12) GOTO 100
        IF (TMA.GT.TMIN) THEN
          IF (.NOT.NLPRCS(J)) THEN
            GOTO 100
          ELSE
            LTSTCXP=.TRUE.
          END IF
        END IF
        WR=A2
        IF (RLB(J).LT.0.) THEN
          GOTO 26
        ELSE IF (RLB(J).GT.0.) THEN
          GOTO 36
        ELSE
          GOTO 16
        END IF
C
C  A1LM(J).NE.0
C
   60   CONTINUE
        A2=A1LM(J)*VX+A2LM(J)*VY+A3LM(J)*VZ
        A3=A0LM(J)+A1LM(J)*X+A2LM(J)*Y+A3LM(J)*Z
        IF (A2.EQ.0.) GOTO 100
        TMA=-A3/A2
        IF (TMA.LE.EPS12) GOTO 100
        IF (TMA.GT.TMIN) THEN
          IF (.NOT.NLPRCS(J)) THEN
            GOTO 100
          ELSE
            LTSTCXP=.TRUE.
          END IF
        END IF
        WR=A2
        YN=Y+TMA*VY
        ZN=Z+TMA*VZ
        XN=ALM(J)+BLM(J)*YN+CLM(J)*ZN
        TUP=TMA
        ICOUNT=2
        IF (RLB(J).LT.0.) THEN
          GOTO 28
        ELSE IF (RLB(J).GT.0.) THEN
          GOTO 38
        ELSE
          GOTO 18
        END IF
C
C  A2LM(J).NE.0
C
   63   CONTINUE
        A2=A1LM(J)*VX+A2LM(J)*VY+A3LM(J)*VZ
        A3=A0LM(J)+A1LM(J)*X+A2LM(J)*Y+A3LM(J)*Z
        IF (A2.EQ.0.) GOTO 100
        TMA=-A3/A2
        IF (TMA.LE.EPS12) GOTO 100
        IF (TMA.GT.TMIN) THEN
          IF (.NOT.NLPRCS(J)) THEN
            GOTO 100
          ELSE
            LTSTCXP=.TRUE.
          END IF
        END IF
        WR=A2
        XN=X+TMA*VX
        ZN=Z+TMA*VZ
        YN=ALM(J)+BLM(J)*XN+CLM(J)*ZN
        TUP=TMA
        ICOUNT=2
        IF (RLB(J).LT.0.) THEN
          GOTO 28
        ELSE IF (RLB(J).GT.0.) THEN
          GOTO 38
        ELSE
          GOTO 18
        END IF
C
C  A3LM(J).NE.0
C
   66   CONTINUE
        A2=A1LM(J)*VX+A2LM(J)*VY+A3LM(J)*VZ
        A3=A0LM(J)+A1LM(J)*X+A2LM(J)*Y+A3LM(J)*Z
        IF (A2.EQ.0.) GOTO 100
        TMA=-A3/A2
        IF (TMA.LE.EPS12) GOTO 100
        IF (TMA.GT.TMIN) THEN
          IF (.NOT.NLPRCS(J)) THEN
            GOTO 100
          ELSE
            LTSTCXP=.TRUE.
          END IF
        END IF
        WR=A2
        XN=X+TMA*VX
        YN=Y+TMA*VY
        ZN=ALM(J)+BLM(J)*XN+CLM(J)*YN
        TUP=TMA
        ICOUNT=2
        IF (RLB(J).LT.0.) THEN
          GOTO 28
        ELSE IF (RLB(J).GT.0.) THEN
          GOTO 38
        ELSE
          GOTO 18
        END IF
C
C
C   DATA RETURNED TO CALLING PROGRAM
C   TENTATIVELY FOR SURFACE NO. J
C
   70   TL=TMX+TADD
        TMIN=TMX
        XR=XN
        YR=YN
        ZR=ZN
        NNR=NN
        VXJ=VX
        VZJ=VZ
        MASURF=J
        SG=SIGN(1._DP,WR)
#ifdef TRACE
        IF (NLTRC) THEN
          WRITE (iunout,*) 'TIMEA, TL,XR,YR,ZR,NNR,MASURF,SG '
          WRITE (iunout,*)         TL,XR,YR,ZR,NNR,MASURF,SG
        ENDIF
#endif
C
C  LOOP OVER SURFACE-INDEX FINISHED
C
  100 CONTINUE
C
C **********************************************************************
C
C  IF NO INTERSECTION FOUND, RETURN
      IF (MASURF.EQ.0) RETURN
C  INTERSECTION AT SURFACE NO. MASURF
C  IF NOT TRANSPARENT, RETURN
      IF (ILIIN(MASURF).GT.0) RETURN
C  IF TRANSPARENT BUT WRONG SIDE, RETURN
      IF (ILSIDE(MASURF)*SG.LT.0) RETURN
C
#ifdef TRACE
      IF (NLTRC) WRITE (iunout,*) 'NOT RETURNED FROM TIMEA, OTHER LOOP '
#endif
C  ILIIN=0, CONTINUE WITH ANOTHER LOOP IN SUBR. TIMEA
C  E.G.: THIS SURFACE MASURF IS A HOLE IN ANOTHER SURFACE.
C  SET STARTING POINT OF RAY TO THIS INTERSECTION AND REPEAT
C  SEARCH FROM THIS NEW POINT.
C
      IF (ILIIN(MASURF).EQ.0) THEN
        X=XR
        XS=X
        Y=YR
        YS=Y
        Z=ZR
        ZS=Z
C
        TADD=TL
        NLLLI=MASURF
        NN=NNR
        NNTCLS=NN
        VX=VXJ
        VXS=VXJ
        VZ=VZJ
        VZS=VZJ
        GOTO 1000
      ELSEIF (ILIIN(MASURF).LT.0) THEN
C  TRANSPARENT, BUT SWITCH AND/OR SURFACE TALLIES
        RETURN
      ENDIF

      RETURN
      END SUBROUTINE EIRENE_TIMEA_CheckInter
