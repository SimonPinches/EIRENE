!****************************************
!* Description:
!*   Implementation of an OCTREE
!*   in Fortran
!* Author:
!*   Oliver Schmidt <o.schmidt@fz-juelich.de>
!*   IEK-4 FZJ
!*   WS 2011
!****************************************
 
      MODULE EIRMOD_OCTREE
 
      USE EIRMOD_PRECISION 
!real precision parameter
      USE EIRMOD_CCONA  
!constants like EPS12, etc
      USE EIRMOD_COMPRT, only: iunout !common printing (unit nums)
      USE EIRMOD_CTRCEI, only: trcoct !tracing switches
 
      IMPLICIT NONE
 
      PRIVATE
 
      PUBLIC :: OCTREE_NewTree,
     .          OCTREE_CreateChildren, OCTREE_AddSurface,
     .          OCTREE_GetLeafchild, OCTREE_Traverse,
     .          OCTREE_CheckBlock, OCTREE_CheckVolume,
     .          OCTREE_Cramer, OCTREE_PrintVTK, OCTREE_PrintGraphviz

      INTEGER, PARAMETER :: less = -1,
     .                      equal = 0,
     .                      more = 1

c     lookup table for the correct combination of coords in
c     OCTREE_CreateChildren
      INTEGER, DIMENSION(8,2) :: XTAB = reshape(
     .                                    (/1, 3, 1, 3, 1, 3, 1, 3,
     .                                      3, 2, 3, 2, 3, 2, 3, 2/),
     .                                    (/8, 2/))
      INTEGER, DIMENSION(8,2) :: YTAB = reshape(
     .                                    (/1, 1, 3, 3, 1, 1, 3, 3,
     .                                      3, 3, 2, 2, 3, 3, 2, 2/),
     .                                    (/8, 2/))
      INTEGER, DIMENSION(8,2) :: ZTAB = reshape(
     .                                    (/1, 1, 1, 1, 3, 3, 3, 3,
     .                                      3, 3, 3, 3, 2, 2, 2, 2/),
     .                                    (/8, 2/))

c     --- DERIVED POINTER TYPE FOR ARRAY ---
      TYPE, PUBLIC :: pOcNode
        TYPE(ocNode), POINTER :: node
      END TYPE pOcNode

c     --- OCTREE TYPE, CONTAINING THE ROOT-NODE ---
      TYPE, PUBLIC :: octree
        TYPE(ocNode), POINTER :: root
c       how many layers will this tree have?
        INTEGER :: layers, maxvalue
        REAL(DP), DIMENSION(3) :: length
c       what is the shortest edge length of the tree, already divided
c       by the factor below!
        REAL(DP) :: shortest
c       by what factor do we divide the length to add to the "ray"?
        REAL(DP) :: factor
c       bounds are the B1 and B2 point of the convex hull with a
c       small amount more to have a slightly larger volume. we 
c       divide by this length, so we don't have problems with floating
c       point arithmetics anymore
        REAL(DP), DIMENSION(3,2) :: bounds
      END TYPE octree

c     --- OCTREE NODE TYPE ---
      TYPE, PUBLIC :: ocNode
c       node id in the directions, coded as integers, using the bits
        INTEGER, DIMENSION(3) :: number
c       remember on which layer in the octree we are, default to layer 1
        INTEGER :: layer = 1
c       as we take cubes as our octree-nodes, define two points
c       => (:,1) left front corner
c          (:,2) right back corner
c          (:,3) center point of the cube
        REAL(DP), DIMENSION(3,3) :: B
c       this is the radius of the convex ball of this block
        REAL(DP) :: RADIUS

c       pointer to parent node
        TYPE(ocNode), POINTER :: parent => NULL()
c       field of 8 child nodes
        TYPE(pOcNode), DIMENSION(8), ALLOCATABLE :: children(:)

c       associated objects with this node
c       -> save the index-number of the add. surface here
        INTEGER, DIMENSION(:), ALLOCATABLE :: surfaces
c       save how many surface we actually have (surfaces array is bigger ;) )
        INTEGER :: nsurfaces = 0
      END TYPE

      CONTAINS
 

      FUNCTION OCTREE_NewTree (X, Y, Z, LAYERS, MAXNSURF) result(tree)

c     --- BUILD A NEW TREE ---
c     x,y,z   : 2d-arrays with the lower left and upper right points
c              of the convex hull of all surfaces we want in our tree
c     layers  : number of layers the tree will have. 0 = lowest level (leafs)
c     maxnsurf : maximum number of surfaces we will have in our tree
c     RETURNS: pointer to the new tree

        IMPLICIT NONE
        TYPE(octree), POINTER :: tree
        TYPE(ocNode), POINTER :: preroot
        REAL(DP), DIMENSION(2), INTENT(IN) :: x, y, z
        REAL(DP), DIMENSION(2) :: xx, yy, zz
        INTEGER, INTENT(IN) :: MAXNSURF, LAYERS

        INTEGER, DIMENSION(3) :: rootnumber = (/0, 0, 0/)

        if(trcoct) WRITE (iunout,*) 'ALLOCATING NEW TREE OBJECT'

c       strech the konvex hull a bit, so we have a closed intervall
c       at the right ends of all directions, as we only check with .lt.
c       for these ends... (otherwise a point exactly on the edge or vertex
c       on the right ends would not be in the block)
c       9.2.12: -> also add a little 
        XX=X
        YY=Y
        ZZ=Z

c       allocate the tree object
        ALLOCATE (tree)
        tree%layers = layers
c       we want #LAYERS bits for the location codes, so we have to
c       start at 2^(layers-1) as 2^0 is the last bit in the integer!
        tree%maxvalue = 2**(layers-1)
        tree%bounds = reshape((/X(1)-EPS12, Y(1)-EPS12, Z(1)-EPS12,
     .                          X(2)+EPS12, Y(2)+EPS12, Z(2)+EPS12/),
     .                        (/3, 2/))
        tree%length = tree%bounds(:,2)-tree%bounds(:,1)
        tree%shortest = 1.D30
        tree%factor = 1.D3

c       build a pre-root node-object so we can use this stuff
c       in the node building routine for the real root
        ALLOCATE(preroot)
c       root node will have layer LAYERS-1, all leafs have level 0
        preroot%layer=LAYERS
        preroot%nsurfaces=MAXNSURF

c       build a new root
        tree%root => OCTREE_NewNode(tree, XX,YY,ZZ, rootnumber, preroot)
c       delete the preroot parent, we don't need this anymore
        tree%root%parent => NULL()
        DEALLOCATE(preroot)
      END FUNCTION OCTREE_NewTree


      FUNCTION OCTREE_NewNode(TREE, X, Y, Z, NUMBER,PARENT) RESULT(NODE)

c     --- CREATE A NEW NODE (AS LEAF) ---
c     x,y,z    : 2d-arrays with the lower left and upper right points
c     number   : the index numbers in the 3 directions (bits are used)
c     parent   : pointer to the parent node
c     RETURNS: pointer to the new node

        IMPLICIT NONE
        REAL(DP), DIMENSION(2), INTENT(IN) :: X, Y, Z
        INTEGER, DIMENSION(3), INTENT(IN) :: NUMBER
        TYPE(ocNode), POINTER, INTENT(IN) :: PARENT
        TYPE(ocTree), POINTER, INTENT(IN) :: TREE

        TYPE(ocNode), POINTER :: NODE
c       the x, y, z coordinates for building the new node
c       -> these coords save the lower and upper edge...
        REAL(DP), DIMENSION(3,3) :: B
c       shortest edge length, divided by factor in tree
        REAL(DP) :: short

c       create the node object
        ALLOCATE (node)

c        NULLIFY(node%children)
        node%parent => PARENT
        node%number = NUMBER
        if (associated(parent)) then
          node%layer = parent%layer-1
c         allocate enough space for our surfaces (max will be
c         number of parents surfaces, so use this)
          ALLOCATE(node%surfaces(parent%nsurfaces))
        else
          node%layer = 1
        end if
        node%nsurfaces = 0

c       calc the center point of our new block
c       ->x3, y3, z3 are in the middle between x1..x2, ...
c       ->as we may assume x2 > x1, y2 > y1, z2 > z1 (max/min in
c         convex hull, this will do...
        B(1,1:2) = x
        B(1,3) = x(1) + abs(x(2)-x(1))*0.5
        B(2,1:2) = y
        B(2,3) = y(1) + abs(y(2)-y(1))*0.5
        B(3,1:2) = z
        B(3,3) = z(1) + abs(z(2)-z(1))*0.5

        if (trcoct) then
          WRITE (iunout,*) 'ALLOCATING NEW NODE OBJECT', NUMBER,
     .                     'ON LAYER', parent%layer-1
          WRITE (iunout,*) 'B1: ', B(:,1)
          WRITE (iunout,*) 'B2: ', B(:,2)
        end if

c       put the coords into the node object
        node%b = B

c       calc the radius of the convex ball of this block (equals the
c       half length of the norm of the vector between the edges)
        node%radius = 0.5 * sqrt((x(2)-x(1))**2 +
     .                           (y(2)-y(1))**2 +
     .                           (z(2)-z(1))**2)
     
c       get the length of the nodes edges, compare with shortest of tree
c       and replace if shorter
        short = minval(B(:,2)-B(:,1))/tree%factor
        if (short .lt. tree%shortest) then
          tree%shortest = short
        end if
      END FUNCTION OCTREE_NewNode


      SUBROUTINE OCTREE_CreateChildren(TREE, PARENT)

c     --- CREATE CHILDREN OF A NODE ---
c     parent : pointer to the parent node (whose children we create ;) )

        IMPLICIT NONE
        TYPE(ocNode), POINTER, INTENT(IN) :: PARENT
        TYPE(ocTree), POINTER, INTENT(IN) :: TREE
        TYPE(ocNode), POINTER :: TMPNODE
        INTEGER, DIMENSION(3) :: TMPNUMBER
        INTEGER :: I, LAYER
c       the x, y, z coordinates for building the new nodes
c       -> these coords save the lower and upper edge...
        REAL(DP), DIMENSION(3,3) :: B
        REAL(DP), DIMENSION(2) :: C1, C2, C3
        
        if (trcoct) then
          WRITE (iunout,*)'WE CREATE CHILDREN FOR BLOCK', 
     .                     parent%number,'ON LAYER', parent%layer
        end if
        
c       new layer!
        layer = parent%layer-1
        
c       fill our coordinates arrays - we already calculated
c       x3, y3, z3 as we created the parent node -> center point!
        B = parent%B

c       allocate the space for the new pointers
        ALLOCATE(parent%children(8))

c       we always create all 8 children at once (anything else makes
c       no sense)
        DO I=1,8
          tmpnumber = parent%number
c         calculate the numbers by setting the appropriate bits
c         -> the schema for setting these follows the coord combi,
c         so reuse this tables! => 1 = 0, 3 = 1
          if (xtab(I,1) .eq. 3) then
            tmpnumber(1) = ibset(tmpnumber(1), layer)
          end if
          if (ytab(I,1) .eq. 3) then
            tmpnumber(2) = ibset(tmpnumber(2), layer)
          end if
          if (ztab(I,1) .eq. 3) then
            tmpnumber(3) = ibset(tmpnumber(3), layer)
          end if
c         create the node
!          tmpnode=>OCTREE_NewNode(tree,
!     .                           (/B(1,xtab(I,1)), B(1,xtab(I,2))/),
!     .                           (/B(2,ytab(I,1)), B(2,ytab(I,2))/),
!     .                           (/B(3,ztab(I,1)), B(3,ztab(I,2))/),
!     .                           tmpnumber, parent)
          c1 = (/B(1,xtab(I,1)), B(1,xtab(I,2))/)
          c2 = (/B(2,ytab(I,1)), B(2,ytab(I,2))/)
          c3 = (/B(3,ztab(I,1)), B(3,ztab(I,2))/)
          tmpnode=>OCTREE_NewNode(tree, c1, c2, c3,
     .                           tmpnumber, parent)

c         put the child into the parents basket
          parent%children(I)%node => tmpnode
        END DO
      END SUBROUTINE OCTREE_CreateChildren

c     --- ASSOCIATE SURFACE TO A NODE ---
c     num  : index number of the surface that is added
c     node : pointer to the node to which we add the surface
      SUBROUTINE OCTREE_AddSurface(NUM, NODE)
        IMPLICIT NONE
        TYPE(ocNode), POINTER, INTENT(IN) :: NODE
        INTEGER, INTENT(IN) :: NUM

        node%nsurfaces = node%nsurfaces + 1
        node%surfaces(node%nsurfaces) = num
      END SUBROUTINE OCTREE_AddSurface


      SUBROUTINE OCTREE_CheckBlock(O, d, BLOCK, STATUS, IP, S)

c     --- CHECK BLOCK INTERSECTION WITH LINE ---
c     This subroutine checks if a line, given by the vector O and the direction
c     in vector d is going through the block given by the pointer BLOCK
c     O      : aufpunkt of the line
c     d      : direction vector of the line
c     BLOCK  : pointer on the block we want to check
c     RETURNS: status=>true/false if the line intersects with the block
c              ip    =>array with the intersection point
c              s     =>the parameter s of the line (s*direction vector)

        IMPLICIT NONE
        TYPE(ocNode), POINTER, INTENT(IN) :: BLOCK
        REAL(DP), DIMENSION(3), INTENT(IN) :: O, d
        REAL(DP), DIMENSION(3), INTENT(OUT) :: IP
        REAL(DP), INTENT(OUT) :: S
        LOGICAL, INTENT(OUT) :: STATUS

        INTEGER :: I, REST
        REAL(DP) :: DOTP, DISTANCE

c       the normal vectors of the 6 planes that form the sides of the block
c       all oriented TO THE INNER!
c       x=1 -> left, -1 -> right
c       y=1 -> front, -1 -> back
c       z=1 -> bottom, -1 -> top
        REAL(DP),DIMENSION(3,6),SAVE ::NORMAL=reshape((/1,0,0,-1, 0, 0,
     .                                                  0,1,0, 0,-1, 0,
     .                                                  0,0,1, 0, 0,-1/)
     .                                                 ,(/3, 6/))

c       first check which sides of the block face the line.
c       -> check what dot product of all normal vectors with d are > 0
c       => only inspect these layers!
        DO I=1,6
c         dot product of side norm vec and d
          DOTP = dot_product(d,NORMAL(:,I))
          IF (DOTP > EPS12) THEN
c           now check on the intersection point. we dont need cramer here,
c           as we have the hesse normal form of the plane:
c           s = -(plane distance + <normal vec , aufpunkt>)/<normal vec, direction>
c           ==> we need to calc the distance of the plane from (0,0,0) first:
c               do this by getting the dot prod. from the norm vec & B1 or B2
c               -> all planes with even I have B2 (x2,y2,z2), uneven I has B1 (x1,y1,z1)
c               -> this simple calc is possible as the planes are parallel to the axes!
c           for even numbered planes this is zero, else 1
            rest = mod(I,2)
            distance = dot_product(NORMAL(:,I),block%B(:,2-rest))
            s = (distance-dot_product(NORMAL(:,I),O))/DOTP
c           put s in our line equation, get our intersection point
            IP = O+s*d
c           is our intersection within our block?
            STATUS = OCTREE_CheckVolume(IP, BLOCK, .TRUE.)
            IF (STATUS) RETURN
          END IF
        END DO
        STATUS = .FALSE.
      END SUBROUTINE OCTREE_CheckBlock



      FUNCTION OCTREE_CheckVolume(point, block, incl) RESULT(yesno)

c     --- CHECK POINT IN BLOCK-VOLUME ---
c     with this we check if a given point is in the volume of the block
c     point   : the point we want to check
c     block   : pointer to the block which volume we check on
c     RETURNS : true or false, depending on in or out of block

        IMPLICIT NONE
        TYPE(ocNode), POINTER, INTENT(IN) :: block
        REAL(DP), DIMENSION(3), INTENT(IN) :: point
        LOGICAL, INTENT(IN) :: incl
        LOGICAL :: yesno

        yesno = .FALSE.
c       we define always an closed intervall at B1 and an open intervall
c       at B2 (the exact point is normally a B1 in another block)
        yesno = ALL(point .ge. (block%b(:,1)-EPS12)) .and.
     .          (ALL(point .lt. (block%b(:,2)+EPS12)) .or.
     .          (incl .and. ALL(point .le. (block%b(:,2)+EPS12))))
      END FUNCTION OCTREE_CheckVolume

c     --- CHECK POINT IN BLOCK-VOLUME INCLUDING B2 COORDS---
c     with this we check if a given point is in the volume of the block
c     point   : the point we want to check
c     block   : pointer to the block which volume we check on
c     RETURNS : true or false, depending on in or out of block
c      FUNCTION OCTREE_CheckVolumeIncl(point, block) RESULT(yesno)
c        IMPLICIT NONE
c        TYPE(ocNode), POINTER, INTENT(IN) :: block
c        REAL(DP), DIMENSION(3), INTENT(IN) :: point
c        LOGICAL :: yesno

c        yesno = .FALSE.
c       we define always an closed intervall at B1 and an open intervall
c       at B2 (the exact point is normally a B1 in another block)
c        IF(ALL(point.ge.block%b(:,1)) .and.
c     .     ALL(point.le.block%b(:,2))) yesno = .TRUE.
c      END FUNCTION OCTREE_CheckVolumeIncl

c     --- GET A LEAF-CHILD FROM THE OCTREE BY A LOCATION CODE ---
c     we get a pointer to the block, where the given point (given as
c     location code in the tree) lies.
c     root    : pointer to the pseudo (or real) root where we start to search
c     location: the location code we are locking for
c     RETURNS : a pointer to the block, in which the location code lies
      FUNCTION OCTREE_GetLeafchild(point, tree) RESULT(child)
        TYPE(ocTree), POINTER, INTENT(IN) :: tree
        REAL(DP), DIMENSION(3), INTENT(IN) :: point
        TYPE(ocnode), POINTER :: child
        INTEGER, DIMENSION(3) :: location
        INTEGER :: nextlevel, branchbit, indx
        
c       convert the given point to a location code we can search for
        location = OCTREE_PointToLocation(tree, point)
c       if any of the location code is < 0 we are NOT in the octree!
        if(ANY(location .lt. 0)) then
          nullify(child)
          return
        end if
        
c       start search at the given root node
        child => tree%root
c       do not stop until we found a leaf
        do while(allocated(child%children))
          nextlevel = child%layer-1
          branchBit = ISHFT(1,nextlevel)
          indx = 4 * ISHFT(IAND(location(3),branchbit),-nextlevel) +
     .           2 * ISHFT(IAND(location(2),branchbit),-nextlevel) +
     .           ISHFT(IAND(location(1),branchbit),-nextlevel) + 1
          child => child%children(indx)%node
        end do
      END FUNCTION OCTREE_GetLeafchild

c     --- TRAVERSE THE TREE ---
c     this function returns a point at which we land, if we "fly"
c     through the given block. this point lies in the correct neighbor 
c     (if there is any).
c     tree    : the actual tree in which we traverse
c     block   : the block through which we "fly"
c     point   : where do we take off?
c     dir     : in which direction?
c     norm    : and how long is that flying vector?
c     RETURNS : the point where the traversel lead us (is E R^3)
      FUNCTION OCTREE_Traverse(tree, block,point,dir,norm) RESULT(travp)
        TYPE(ocTree), POINTER, INTENT(IN) :: tree
        TYPE(ocNode), POINTER, INTENT(IN) :: block
        REAL(DP), DIMENSION(3), INTENT(IN) :: point, dir
        REAL(DP), INTENT(IN) :: norm
        REAL(DP), DIMENSION(3) :: travp
        INTEGER :: I, REST
        REAL(DP) :: DOTP, DISTANCE, S, S_MIN

c       the normal vectors of the 6 planes that form the sides of the block
c       all oriented TO THE OUTER!
c       x=-1 -> left, 1 -> right
c       y=-1 -> front, 1 -> back
c       z=-1 -> bottom, 1 -> top
        REAL(DP),DIMENSION(3,6),SAVE ::NORMAL=reshape((/-1,0,0, 1,0,0,
     .                                                  0,-1,0, 0,1,0,
     .                                                  0,0,-1, 0,0,1/)
     .                                                 ,(/3, 6/))
     
c       we want to traverse our "ray" through the block we are in,
c       so we first get the correct intersection point of ray and block
c       walls...

c       first check which sides of the block face the line.
c       -> check what dot product of all normal vectors with dir are > 0
c       => only inspect these layers!
        S_MIN = 1.D30
        DO I=1,6
c         dot product of side norm vec and direction dir
          DOTP = dot_product(dir,NORMAL(:,I))
          IF (DOTP > EPS12) THEN
c           now check on the intersection point. we dont need cramer here,
c           as we have the hesse normal form of the plane:
c           s = -(plane distance + <normal vec , aufpunkt>)/<normal vec, direction>
c           ==> we need to calc the distance of the plane from (0,0,0) first:
c               do this by getting the dot prod. from the norm vec & B1 or B2
c               -> all planes with even I have B2 (x2,y2,z2), uneven I has B1 (x1,y1,z1)
c               -> this simple calc is possible as the planes are parallel to the axes!
c           for even numbered planes this is zero, else 1
            rest = mod(I,2)
            distance = dot_product(NORMAL(:,I),block%B(:,2-rest))
            s = (distance-dot_product(NORMAL(:,I),point))/DOTP
c           if this runlength s is shorter than a former, replace it
c           -> S is always positive because we do not check on surfaces in
c           our back!
            if(s .lt. s_min ) s_min = s
          END IF
        END DO
        
c       get the resulting intersection point out of the search...
        travp = point + s_min * dir
        
c       now as we have the point where we wanted to be, add an extra
c       to be sure we land in the neighbor block
        travp = travp + tree%shortest/norm*dir
      END FUNCTION OCTREE_Traverse


c     --- 3D POINT TO LOCATION CODE ---
c     this converts a point from R^3 to the octree space
c     location code
c     tree   : pointer to the tree which space we use
c     point  : vector with the r^3 coords
c     RETURNS: the location code...
      FUNCTION OCTREE_PointToLocation(tree, point) RESULT(convl)
        TYPE(octree), POINTER, INTENT(IN) :: tree
        REAL(DP), DIMENSION(3), INTENT(IN) :: point
        REAL(DP), DIMENSION(3) :: convp
        INTEGER, DIMENSION(3) :: convl
c       convert the point out of the normal space into the octree space
c       [0,1]x[0,1]x[0,1]
        convp = (point-tree%bounds(:,1))/tree%length
c       if the point is not in the space, there are values <0 or >1
c       so return a -1,-1,-1 for the user
        IF(ANY(convp .gt. 1 .or. convp .lt. 0))
     .    convp = (/-1.0_DP, -1.0_DP, -1.0_DP/)
c       if this point is in the octree space, convert it to a location code
        convl = int(convp * tree%maxvalue)
      END FUNCTION OCTREE_PointToLocation

c     --- SOLVE BARYCENTRIC COORDINATE LES WITH CRAMER ---
c     This function calculates the linear equation system
c     (A-B A-C d)*(beta,gamma,t)^T=(A-O)
c     which is the calculation of the intersection point and
c     if the point is inside of the triangle defined by the
c     three points A,B,C.
c     The intersection point is between the plane (defined by triangle)
c     and the line g: O+t*d
c     A, B, C, 0, d: double array with xyz containing the points/vectors
c                    mentioned above
c     RETURNS: res = (0,0,0) if no solution or if beta|gamma|t <= 0, else
c                    the solution of LES
      FUNCTION OCTREE_Cramer(A, B, C, O, d) RESULT(RES)
        REAL(DP), DIMENSION(3), INTENT(IN) :: A,B,C,O,d
        REAL(DP), DIMENSION(3) :: RES, BB
        REAL(DP), DIMENSION(3,3) :: AA
        REAL(DP) :: detA, detA1, detA2, detA3

        AA(:,1) = A-B
        AA(:,2) = A-C
        AA(:,3) = d
        BB = A-O

        detA  = AA(1,1)*AA(2,2)*AA(3,3)+AA(1,2)*AA(2,3)*AA(3,1)+AA(1,3)*
     .          AA(2,1)*AA(3,2)-AA(1,3)*AA(2,2)*AA(3,1)-AA(1,2)*AA(2,1)*
     .          AA(3,3)-AA(1,1)*AA(2,3)*AA(3,2)
c       system not solvable, return 000
        IF (abs(detA).lt.EPS12) THEN
          RES = (/0, 0, 0/)
          RETURN
        END IF

c       calc x1 (=beta)
        detA1 = BB(1)*AA(2,2)*AA(3,3)+AA(1,2)*AA(2,3)*BB(3)+AA(1,3)*
     .          BB(2)*AA(3,2)-AA(1,3)*AA(2,2)*BB(3)-AA(1,2)*BB(2)*
     .          AA(3,3)-BB(1)*AA(2,3)*AA(3,2)
        RES(1) = detA1/detA
c       if the parameter is not >0, we can stop here
        IF (RES(1).le.0) THEN
          RES = (/0, 0, 0/)
          RETURN
        END IF

c       calc x2 (=gamma)
        detA2 = AA(1,1)*BB(2)*AA(3,3)+BB(1)*AA(2,3)*AA(3,1)+AA(1,3)*
     .          AA(2,1)*BB(3)-AA(1,3)*BB(2)*AA(3,1)-BB(1)*AA(2,1)*
     .          AA(3,3)-AA(1,1)*AA(2,3)*BB(3)
        RES(2) = detA2/detA
c       if the parameter is not >0 or if beta + gamma >=1, we can stop here
        IF (RES(2).le.0 .or. (RES(1)+RES(2)) .ge. 1 ) THEN
          RES = (/0, 0, 0/)
          RETURN
        END IF

c       calc x3 (=t =line-equation parameter)
        detA3 = AA(1,1)*AA(2,2)*BB(3)+AA(1,2)*BB(2)*AA(3,1)+BB(1)*
     .          AA(2,1)*AA(3,2)-BB(1)*AA(2,2)*AA(3,1)-AA(1,2)*AA(2,1)*
     .          BB(3)-AA(1,1)*BB(2)*AA(3,2)
        RES(3) = detA3/detA
c       if the parameter is not >0, we can stop here
        IF (RES(3).le.0) RES = (/0, 0, 0/)
      END FUNCTION OCTREE_Cramer
      
c     --- PRINTING OCTREE IN XML FOR VTK TOOL ---
c     this function prints the octree as cubes for the vtk tool EGview3D
c     tree : pointer to the tree that we want to be printed
c     unum : unit number of the ALREADY OPEN file we write to
      SUBROUTINE OCTREE_PrintVTK(tree, unum)
        TYPE(ocTree), POINTER, INTENT(IN) :: tree
        INTEGER, INTENT(IN) :: unum
        
        write(unum,*) '<cubes>'
c       recursive tour through the octree, printing leaf childs only
        call PrintCube(tree%root, unum)
        write(unum,*) '</cubes>'
        
      END SUBROUTINE OCTREE_PrintVTK
      
      RECURSIVE SUBROUTINE PrintCube(block, unum)
        TYPE(ocNode), POINTER, INTENT(IN) :: block
        INTEGER, INTENT(IN) :: unum
        INTEGER :: i
        
c       if this is a leaf, print it
        if(.not. allocated(block%children)) then
          write(unum,*) '<cube>'
          write(unum,*) '<point>',block%B(:,1),'</point>'
          write(unum,*) '<point>',block%B(:,2),'</point>'
          write(unum,*) '</cube>'
c       else go one level deeper
        else
          do i=1,8
            call PrintCube(block%children(I)%node, unum)
          end do
        end if
      END SUBROUTINE PrintCube
      
      SUBROUTINE OCTREE_PrintGraphviz(tree, unum)
        TYPE(ocTree), POINTER, INTENT(IN) :: tree
        INTEGER, INTENT(IN) :: unum
        
        write(unum,*) 'digraph G {'
c       recursive tour through the octree, printing leaf childs only
        call PrintDot(tree%root, unum)
        write(unum,*) '}'
      END SUBROUTINE OCTREE_PrintGraphviz
      
      RECURSIVE SUBROUTINE PrintDot(block, unum)
        TYPE(ocNode), POINTER, INTENT(IN) :: block
        INTEGER, INTENT(IN) :: unum
        INTEGER :: i
        CHARACTER(40) :: me, kidname
        TYPE(ocNode), POINTER :: child
        
        WRITE(me,"(4(i0))") block%layer,block%number
        
c       first print ourself
        write(unum,"('  ',a,a,a,':',i0,a)")
     .        trim(adjustl(me)),' [label="',trim(adjustl(me)),
     .        block%nsurfaces,'"];'

c       now lets print our children and call with them on this again
c       to print their labels
        if(allocated(block%children)) then
          do i=1,8
            child => block%children(I)%node
            WRITE(kidname,"(4(i0))") child%layer,child%number
            
            write(unum,*)'  ',trim(adjustl(me)),
     .                   ' -> ',trim(adjustl(kidname)),';'
            call PrintDot(child, unum)
          end do
        end if
      END SUBROUTINE PrintDot
      

      END MODULE EIRMOD_OCTREE
 
 
