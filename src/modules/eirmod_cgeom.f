      MODULE EIRMOD_CGEOM
cdr  precomputed grid quantities: cell volumes, center of mass,
cdr  typical cell diameter

cdr  cleanup needed
cdr  comments re: cell_elem, cell_list needed

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: EIRENE_ALLOC_CGEOM, EIRENE_DEALLOC_CGEOM,
     P          EIRENE_INIT_CGEOM, EIRENE_BROADCAST_CGEOM,
     P          CELL_ELEM, CELL_LIST

      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     R        RCGM1(:), RCGM2(:,:),   AREAG(:)

      REAL(DP), PUBLIC, POINTER, SAVE ::
     R VOLADD(:), VOLCOR(:), VOLG(:), VOLTAL(:), VOLTOT,
     R AREA(:),   CELDIA(:), XCOM(:), YCOM(:),
     R XPOINT(:), YPOINT(:)

cdr  only for polygons ?
      REAL(DP), PUBLIC, POINTER, SAVE ::
     R XPOL(:,:), YPOL(:,:)

      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I NPOINT(:,:),   NSTGRD(:),
     I NGHPLS(:,:,:),  ! neighboring polygon surface
     I NGHPOL(:,:,:),  ! neighboring polygon cell
     I NCLTAL(:),      ! index mapping: fine - coarse grid
     I INDPOINT(:,:), NOPNT(:)

      INTEGER, PUBLIC, SAVE :: NCGM1, NCGM2, NNODES

cdr  damaged cell ?
      LOGICAL, PUBLIC, ALLOCATABLE, SAVE ::
     L LDAMCEL(:)

cdr  comments ??
      TYPE :: CELL_ELEM
        INTEGER :: NOCELL
        TYPE(CELL_ELEM), POINTER :: NEXT_CELL
      END TYPE CELL_ELEM

      TYPE :: CELL_LIST
        TYPE(CELL_ELEM), POINTER :: PCELL
      END TYPE CELL_LIST

      TYPE(CELL_LIST), ALLOCATABLE, SAVE, PUBLIC :: COORCELL(:)


      CONTAINS


      SUBROUTINE EIRENE_ALLOC_CGEOM(ICAL)
c  allocate storage for geometrical arrays.
c  parameter NRAD: size of grid. NRADS: size of grid after possible elimination
c                                       of storage, in module parmmod.f (FLAG: NGEOM_USR)

      INTEGER, INTENT(IN) :: ICAL

      IF (ICAL == 1) THEN
         IF (ALLOCATED(RCGM1)) RETURN

         NCGM1 = NADD+NBMAX+7*NRAD+NRTAL+1
         NCGM2 = 2*N1STS*N2NDPLGS
c  real arrays
         ALLOCATE (RCGM1(NCGM1))
         ALLOCATE (RCGM2(N1STS,2*N2NDPLGS))
c  integer arrays
         ALLOCATE (NPOINT(2,NPPART))
         ALLOCATE (NSTGRD(NRAD))
         ALLOCATE (NGHPLS(4,N1STS,N2NDPLGS))
         ALLOCATE (NGHPOL(4,N1STS,N2NDPLGS))
         ALLOCATE (NCLTAL(NRAD))
         ALLOCATE (INDPOINT(N1STS,N2NDPLGS))
         ALLOCATE (NOPNT(NRAD))

         ALLOCATE (COORCELL(NRAD))

         ALLOCATE (LDAMCEL(NRAD))

         WRITE (IUNMEM,'(A,T25,I15)')
     .        ' CGEOM(1) ',(NCGM1+2*N1STS*N2NDPLGS)*8 +
     .        (2*NPPART+2*NRAD+8*N1STS*N2NDPLGS)*4 + nrad*4

         VOLADD => RCGM1(1 : NADD)
         VOLTAL => RCGM1(1+NADD : NADD+NRTAL)
         VOLCOR => RCGM1(1+NADD+NRTAL : NADD+NRTAL+NBMAX)
         VOLG   => RCGM1(1+NADD+NRTAL+NBMAX : NADD+NRTAL+NBMAX+NRAD)
         AREA   => RCGM1(1+NADD+NRTAL+NBMAX+NRAD :
     .                     NADD+NRTAL+NBMAX+2*NRAD)
         CELDIA => RCGM1(1+NADD+NRTAL+NBMAX+2*NRAD :
     .                     NADD+NRTAL+NBMAX+3*NRAD)
         XCOM   => RCGM1(1+NADD+NRTAL+NBMAX+3*NRAD :
     .                     NADD+NRTAL+NBMAX+4*NRAD)
         YCOM   => RCGM1(1+NADD+NRTAL+NBMAX+4*NRAD :
     .                     NADD+NRTAL+NBMAX+5*NRAD)
         XPOINT => RCGM1(1+NADD+NRTAL+NBMAX+5*NRAD :
     .                     NADD+NRTAL+NBMAX+6*NRAD)
         YPOINT => RCGM1(1+NADD+NRTAL+NBMAX+6*NRAD :
     .                     NADD+NRTAL+NBMAX+7*NRAD)
         VOLTOT => RCGM1(1+NADD+NRTAL+NBMAX+7*NRAD)

         XPOL => RCGM2(:,1:N2NDPLGS)
         YPOL => RCGM2(:,1+N2NDPLGS:2*N2NDPLGS)

      ELSE IF (ICAL == 2) THEN

         IF (ALLOCATED(AREAG)) RETURN
         ALLOCATE (AREAG(NLMPGS))
         WRITE (IUNMEM,'(A,T25,I15)') ' CGEOM(2) ',NLMPGS*8

      END IF


      CALL EIRENE_INIT_CGEOM (ICAL)

      RETURN
      END SUBROUTINE EIRENE_ALLOC_CGEOM


      SUBROUTINE EIRENE_DEALLOC_CGEOM

      IF (.NOT.ALLOCATED(RCGM1)) RETURN

      DEALLOCATE (RCGM1)
      DEALLOCATE (RCGM2)

      DEALLOCATE (NPOINT)
      DEALLOCATE (NSTGRD)
      DEALLOCATE (NGHPLS)
      DEALLOCATE (NGHPOL)
      DEALLOCATE (NCLTAL)
      DEALLOCATE (INDPOINT)
      DEALLOCATE (NOPNT)

      DEALLOCATE (COORCELL)

      DEALLOCATE (LDAMCEL)

      DEALLOCATE (AREAG)

      RETURN
      END SUBROUTINE EIRENE_DEALLOC_CGEOM


      SUBROUTINE EIRENE_INIT_CGEOM(ICAL)

      INTEGER, INTENT(IN) :: ICAL
      INTEGER :: I

      IF (ICAL == 1) THEN

         RCGM1    = 0._DP
         RCGM2    = 0._DP

         NPOINT   = 0
         NSTGRD   = 0
         NGHPLS   = 0
         NGHPOL   = 0
         NCLTAL   = 0
         INDPOINT = 0
         NOPNT = 0

         DO I=1,NRAD
            NULLIFY (COORCELL(I)%PCELL)
         END DO

         LDAMCEL = .FALSE.

      ELSE IF (ICAL == 2) THEN

         AREAG = 0._DP

      END IF

      RETURN
      END SUBROUTINE EIRENE_INIT_CGEOM

      
      SUBROUTINE EIRENE_BROADCAST_CGEOM(ME)
      USE EIRMOD_MPI
      INTEGER, INTENT(IN) :: ME
      INTEGER :: IER

      IF (ME /= 0) THEN
        CALL EIRENE_ALLOC_CGEOM(1)
        CALL EIRENE_ALLOC_CGEOM(2)
      END IF

      CALL MPI_BCAST (RCGM1,NCGM1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (RCGM2,NCGM2,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (AREAG,NLMPGS,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NPOINT,2*NPPART,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NSTGRD,NRAD,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NGHPLS,4*N1STS*N2NDPLGS,
     .                MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NGHPOL,4*N1STS*N2NDPLGS,
     .                MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NCLTAL,NRAD,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NNODES,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (LDAMCEL,NRAD,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
      
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

      END SUBROUTINE EIRENE_BROADCAST_CGEOM
      

      END MODULE EIRMOD_CGEOM
