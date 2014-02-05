c*********************************************
c     OCTREE USAGE FOR TIMEA-SPEEDUP
c     This file contains all necessary routines to use the octree
c     data structure defined in EIRENE_OCTREE modul.
c     It is only used by timea.f.
c     -------------------------------
c     Author: Oliver Schmidt <o.schmidt@fz-juelich.de>
c     Revision: 0.1
c     Date: 25.11.2011
c     -------------------------------
c*********************************************

c     --- EIRENE_TIMEA_BuildOctree ---
c     This function builds a new octree with the additional surfaces.
c     It creates first a convex hull around all of them; right now
c     only triangles, quadrangles or five-edge-objects are used.
c     The tree is built recursivly after this.
c     RETURNS: pointer to the freshly built tree
      FUNCTION EIRENE_TIMEA_BuildOctree() RESULT(TREE)
c       this imports all points (p1-p5), rlb, igjum0 stuff, which
c       has been worked on before this function is called for the first time
        USE EIRMOD_PRECISION
        USE EIRMOD_CADGEO
        USE EIRMOD_CLGIN
        USE EIRMOD_OCTREE
        USE EIRMOD_COMPRT
        USE EIRMOD_CTRCEI
        IMPLICIT NONE
c       define an explicit interface to the block building function
        INTERFACE
          RECURSIVE SUBROUTINE EIRENE_TIMEA_BuildBlocks(TREE, PARENT)
            USE EIRMOD_OCTREE
            TYPE(ocNode), POINTER, INTENT(IN) :: PARENT
            TYPE(ocTree), POINTER, INTENT(IN) :: TREE
          END SUBROUTINE EIRENE_TIMEA_BuildBlocks
        END INTERFACE

        TYPE(octree), POINTER :: TREE
        TYPE(ocNode), POINTER :: ROOT
        REAL(DP), DIMENSION(2) :: XLIM, YLIM, ZLIM
        REAL(DP) :: EIRENE_SECOND_OWN, start
        INTEGER :: J, LAYERS, NSURFACES
        
        if (trcoc) then
          WRITE (iunout,*)
          WRITE (iunout,*) 'STARTING ON BUILDING OCTREE...'
          start = EIRENE_SECOND_OWN()
        end if

c       init the maxima and minima with a very low value for max and very
c       high value for min, the first real point will change this to
c       actual values...
        XLIM(1)=1.D55
        XLIM(2)=-1.D55
        YLIM(1)=1.D55
        YLIM(2)=-1.D55
        ZLIM(1)=1.D55
        ZLIM(2)=-1.D55

        nsurfaces = 0
C       iterate over all add. surfaces
        DO J=1,NLIMI
C         we now get the min/max from our coords for building a konvex hull
C         over our add. surfaces, ignore higher order surfaces, ignore two point lines
c         !we also ignore all surfaces which have been chosen to be ignored!
          if (RLB(J) .GE. 3 .AND. IGJUM0(J).EQ.0) then
            XLIM(1)=MIN(P1(1,J),P2(1,J),P3(1,J),P4(1,J),P5(1,J),XLIM(1))
            XLIM(2)=MAX(P1(1,J),P2(1,J),P3(1,J),P4(1,J),P5(1,J),XLIM(2))
            YLIM(1)=MIN(P1(2,J),P2(2,J),P3(2,J),P4(2,J),P5(2,J),YLIM(1))
            YLIM(2)=MAX(P1(2,J),P2(2,J),P3(2,J),P4(2,J),P5(2,J),YLIM(2))
            ZLIM(1)=MIN(P1(3,J),P2(3,J),P3(3,J),P4(3,J),P5(3,J),ZLIM(1))
            ZLIM(2)=MAX(P1(3,J),P2(3,J),P3(3,J),P4(3,J),P5(3,J),ZLIM(2))

c           find out how many surfaces we will get as a maximum
            nsurfaces = nsurfaces + 1
          end if
        END DO
        
c       if we do not have any surfaces we could add, do not build a tree
        if(nsurfaces == 0) then
          if(trcoc) WRITE(iunout,*)'NO USABLE SURFACES FOUND, SKIPPING'
          tree => null()
          return
        end if

c       calc how many layers we will get by the amount of surfaces
c       we have: int(log10(# surface)) + 1 (add 0.0 to have a real value for
c       log10()-function)
        LAYERS = floor(log10(nsurfaces+0.0))+1
c        LAYERS=32
        
c       for only 10 surfaces we do not split the root...
c       -> this shall be the same criteria we have for deciding of
c          splitting up child nodes
        IF(nsurfaces .le. 10) LAYERS=1

        if (trcoc) then
          WRITE (iunout,*) 'CONVEX HULL: X,   Y,   Z'
          WRITE (iunout,*) 'MIN:', xlim(1), ylim(1), zlim(1)
          WRITE (iunout,*) 'MAX:', xlim(2), ylim(2), zlim(2)
          WRITE (iunout,*)
          WRITE (iunout,*) 'TREE WILL HAVE', layers, 'LAYER(S)'
          WRITE (iunout,*) 'SURFACES:', nsurfaces, 'OUT OF', nlimi,
     .                    'WILL BE USED'
          WRITE (iunout,*)
        end if

c       now build a new octree with a new root which is konvex hull
c       with the xyzlims
        TREE => OCTREE_NewTree(XLIM, YLIM, ZLIM, LAYERS, nsurfaces)
        ROOT => tree%root

c       associate all surfaces with the root element
        DO J=1,NLIMI
C         ignore higher order surfaces, ignore two point lines
c         !we also ignore all surfaces which have been chosen to be ignored!
          if (RLB(J) .GE. 3 .AND. IGJUM0(J).EQ.0) then
            CALL OCTREE_AddSurface(J, ROOT)
          end if
        END DO

c       now we build our real tree by adding children to the root node
c       and associating them with the correct surfaces
c       -> only build children if we want more than 1 layer (=> NLIMI >= 10)
c          and if there actually are surfaces matching our criterias,
c          which actually have been associated with the root node
        IF(root%nsurfaces .gt. 0 .and. layers .gt. 1)THEN
          CALL EIRENE_TIMEA_BuildBlocks(tree, root)
        END IF
        
        if (trcoc) then
          WRITE(iunout,*)
          WRITE(iunout,*) 'OCTREE BUILD PHASE FINISHED'
          WRITE(iunout,*)'CPU TIME USED IN BUILD:',
     .                   EIRENE_SECOND_OWN()-start
          WRITE(iunout,*)
        end if
      END FUNCTION EIRENE_TIMEA_BuildOctree

c     --- EIRENE_TIMEA_BuildBlocks ---
c     Via this recursive function we are able to build all these nasty
c     little blocks with their children...
c     RETURNS: number of surfaces that have been assigned to a node
      RECURSIVE SUBROUTINE EIRENE_TIMEA_BuildBlocks(TREE, PARENT)
c       this imports all points (p1-p5), rlb, igjum0 stuff, which
c       has been worked on before this function is called for the first time
        USE EIRMOD_PRECISION
        USE EIRMOD_COMPRT
        USE EIRMOD_CADGEO
        USE EIRMOD_CLGIN
        USE EIRMOD_CCONA
        USE EIRMOD_OCTREE
        USE EIRMOD_CTRCEI
        IMPLICIT NONE

c       definie explicit interfaces to the checking routines
        INTERFACE
          LOGICAL FUNCTION EIRENE_TIMEA_FirstCheck(SID,CHILD)
            USE EIRMOD_OCTREE
            TYPE(ocNode), POINTER, INTENT(IN) :: CHILD
            INTEGER, INTENT(IN) :: SID
          END FUNCTION
          LOGICAL FUNCTION EIRENE_TIMEA_ThirdCheck(SID,CHILD)
            USE EIRMOD_OCTREE
            TYPE(ocNode), POINTER, INTENT(IN) :: CHILD
            INTEGER, INTENT(IN) :: SID
          END FUNCTION
        END INTERFACE

        TYPE(ocNode), POINTER, INTENT(IN) :: PARENT
        TYPE(ocTree), POINTER, INTENT(IN) :: TREE
        TYPE(ocNode), POINTER :: CHILD
        INTEGER :: I, J, SID
        REAL(DP) :: DISTANCE

c       give birth to childs of the parental node
        CALL OCTREE_CreateChildren(TREE, PARENT)

c       now iterate over all children and give them the love... eh... the
c       surfaces they need
        DO I=1,8
c         which child do we operate on?
          CHILD => parent%children(i)%node

c         now iterate over all surface of the parent and check if the
c         child has this surface...
          DO J=1,parent%nsurfaces
c           fetch real surface number
            SID = parent%surfaces(J)

c           only for nice 3,4 or 5 edge objects
            IF (RLB(SID) .ge. 3) THEN
c             first check: is at least one point of this surface
c                          in the block?
              IF (EIRENE_TIMEA_FirstCheck(SID, CHILD)) THEN
                if(trcoc) WRITE(iunout,*) 'FIRST CHECK: ASSOC. SURFACE',
     .                             SID, 'WITH NODE', child%number,
     .                             'ON LAYER', child%layer
                  CALL OCTREE_AddSurface(SID, CHILD)
                  CYCLE
              END IF

c             second check: is the surface at a distance from center less
c                           then the radius of the convex ball?
c                           (=distance point <-> plane)
c             -> A0LM, A1LM, A2LM, A3LM are the coord.version of the plane
c             -> x3,y3,z3 is the center point of the block
c             first calc the numerator, then the denominator
              distance=abs(A0LM(SID) + A1LM(SID)*child%b(1,3) +
     .                    A2LM(SID)*child%b(2,3)+A3LM(SID)*child%b(3,3))
              distance = distance /
     .                   sqrt(A1LM(SID)**2+A2LM(SID)**2+A3LM(SID)**2)
c             is this distance smaller than the radius + very tiny epsilon?
              IF (distance .le. child%radius+EPS10) THEN
c               third check: now we have to check on real intersections.
c                           if we intersect, add this surface to the block
                IF (EIRENE_TIMEA_ThirdCheck(SID, CHILD)) THEN
                  if(trcoc)WRITE(iunout,*)'THIRD CHECK: ASSOC. SURFACE',
     .                             SID, 'WITH NODE', child%number,
     .                             'ON LAYER', child%layer
                  CALL OCTREE_AddSurface(SID, CHILD)
                END IF
              END IF
            END IF
          END DO

c         if we have not reached our lowest level (=0), do another recursive step
c         and let there be children, don't do so if we have less then 4 surfaces added
c         to this block!
          IF(child%layer .gt. 0 .and. child%nsurfaces .gt. 10)THEN
            CALL EIRENE_TIMEA_BuildBlocks(TREE, child)
          END IF
        END DO

      END SUBROUTINE EIRENE_TIMEA_BuildBlocks


