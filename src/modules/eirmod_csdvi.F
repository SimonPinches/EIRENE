cdr july 17:  nsnv already contained in nspztot, due to previous fix.
cdr nov. 17:  nspztotw introduced, in analogy with nspztot

      MODULE EIRMOD_CSDVI

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: EIRENE_ALLOC_CSDVI, EIRENE_DEALLOC_CSDVI,
     P          EIRENE_INIT_CSDVI, EIRENE_BROADCAST_CSDVI

      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     R        SDVI1(:,:), SDVI2(:,:)

      REAL(DP), PUBLIC, POINTER, SAVE ::
     R SIGMA(:,:),    SGMS(:),
     R SIGMAW(:,:),   SGMWS(:)

      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R SIGMAC(:,:,:), SGMCS(:,:)

      INTEGER, PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     I         ISDVI(:)

      INTEGER, PUBLIC, POINTER, SAVE ::
     I IIH(:),    IGH(:),
     I IIHW(:),   IGHW(:),
     I ICOV(:),
     I NSIGI,     NSIGVI,   NSIGSI, NSIGCI,
     I NSIGI_SPC
C SPEED UP OF SUBROUTINE STATIS

cdr This list stores information along history.
cdr It is refreshed in MCARLO.f after each completed history.
cdr Currently stored per history:
c    cell visited:    ICELL=ICLMT(IC),   IC=1,..,NCLMTS
c    cell score:      SCELL=SCLMT(IS,IC),IC=1,..,NCLMTS, IS=1,.. tbd.
c    surface visited: ISURF=IWLMT(ICW), ICW=1,..,NWLMTS

c    species scoring in volume, somewhere:  ISPZ=LMETSP (...)
c    species scoring at surface, somewhere: ISPZ=LMETSPW(...)

      INTEGER, PUBLIC, POINTER, SAVE ::
     I IMETCL(:), ICLMT(:), NCLMT, NCLMTS,
     I IMETWL(:), IWLMT(:), NWLMT, NWLMTS

      LOGICAL, PUBLIC, ALLOCATABLE, TARGET, SAVE ::
     L LMETSP(:), LMETSPW(:)

      LOGICAL, PUBLIC, POINTER, SAVE ::
     l LMETSP2(:,:), LMETSPW2(:,:)

#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP  THREADPRIVATE(ISDVI,IMETCL,ICLMT,NCLMT,NCLMTS,IMETWL,IWLMT,
!$OMP& NWLMT,NWLMTS)
#endif

      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I IIHC(:,:), IGHC(:,:)

      INTEGER, PUBLIC, SAVE ::
     I NSDVI1, NSDVI2, NSDVC1, NSDVC2, NSDVI, MSDVI


      CONTAINS


      SUBROUTINE EIRENE_ALLOC_CSDVI (ICAL)

      INTEGER, INTENT(IN) :: ICAL

      IF (ICAL == 1) THEN

        IF (ALLOCATED(SDVI1)) RETURN

        NSDVI1 = NSD*(NRTAL+1)
        NSDVI2 = NSDW*(NLIMPS+1)
        NSDVC1 = 3*NCV*NRTAL
        NSDVC2 = 3*NCV
        NSDVI  = NSDVI1+NSDVI2+NSDVC1+NSDVC2
        MSDVI  = NSD*2+NSDW*2+NCV+7+2*NRTAL+2+2*NLIMPS+2

        ALLOCATE (SDVI1(NSD,NRTAL+1))
        ALLOCATE (SDVI2(NSDW,NLIMPS+1))
        ALLOCATE (SIGMAC(0:2,NCV,NRTAL))
        ALLOCATE (SGMCS(0:2,NCV))

        ALLOCATE (ISDVI(MSDVI))
        ALLOCATE (IIHC(2,NCV))
        ALLOCATE (IGHC(2,NCV))

        WRITE (IUNMEM,'(A,T25,I15)')
     .        ' CSDVI ',(NSD*(NRTAL+1) + NSDW*(NLIMPS+1) +
     .                   (NRTAL+1)*3*NCV)*8 +
     .                  (MSDVI+4*NCV)*4

        SIGMA  => SDVI1(:,1:NRTAL)
        SGMS   => SDVI1(:,NRTAL+1)

        SIGMAW => SDVI2(:,1:NLIMPS)
        SGMWS  => SDVI2(:,NLIMPS+1)

        NSIGI     => ISDVI(1)
        NSIGVI    => ISDVI(2)
        NSIGSI    => ISDVI(3)
        NSIGCI    => ISDVI(4)
c   out           => ISDVI(5)
c   out           => ISDVI(6)
        NSIGI_SPC => ISDVI(7)
        NCLMT     => ISDVI(8)
        NCLMTS    => ISDVI(9)
        NWLMT     => ISDVI(10)
        NWLMTS    => ISDVI(11)
        IIH       => ISDVI(12                  : 11+  NSD)
        IGH       => ISDVI(12+  NSD            : 11+2*NSD)
        IIHW      => ISDVI(12+2*NSD            : 11+2*NSD+  NSDW)
        IGHW      => ISDVI(12+2*NSD+  NSDW     : 11+2*NSD+2*NSDW)
        ICOV      => ISDVI(12+2*NSD+2*NSDW     : 11+2*NSD+2*NSDW+NCV)
        IMETCL    => ISDVI(12+2*NSD+2*NSDW+NCV :
     .                     11+2*NSD+2*NSDW+NCV+NRTAL)
        ICLMT     => ISDVI(12+2*NSD+2*NSDW+NCV+NRTAL :
     .                     11+2*NSD+2*NSDW+NCV+2*NRTAL)
        IMETWL    => ISDVI(12+2*NSD+2*NSDW+NCV+2*NRTAL :
     .                     11+2*NSD+2*NSDW+NCV+2*NRTAL+NLIMPS)
        IWLMT    => ISDVI(12+2*NSD+2*NSDW+NCV+2*NRTAL+NLIMPS : MSDVI)

      ELSE IF (ICAL == 2) THEN

        IF (ALLOCATED(LMETSP)) RETURN

c  arrays updated on the fly, to speed up statistical variance calculations in statis.f
c  lmetsp(i1) : logical, indicator, whether for any particular trajectory a
c               volume-averaged tally has been scored with species (first) index "i1"
c  lmetspw(i1): ditto, for surface-averaged tallies
        ALLOCATE (LMETSP(NSPZTOT))
        ALLOCATE (LMETSPW(NSPZTOTW))

      END IF

      CALL EIRENE_INIT_CSDVI (ICAL)

      RETURN
      END SUBROUTINE EIRENE_ALLOC_CSDVI


      SUBROUTINE EIRENE_DEALLOC_CSDVI

      IF (.NOT.ALLOCATED(SDVI1)) RETURN

      DEALLOCATE (SDVI1)
      DEALLOCATE (SDVI2)
      DEALLOCATE (SIGMAC)
      DEALLOCATE (SGMCS)

      DEALLOCATE (ISDVI)
      DEALLOCATE (IIHC)
      DEALLOCATE (IGHC)

      DEALLOCATE (LMETSP)
      DEALLOCATE (LMETSPW)

      RETURN
      END SUBROUTINE EIRENE_DEALLOC_CSDVI


      SUBROUTINE EIRENE_INIT_CSDVI (ICAL)

      INTEGER, INTENT(IN) :: ICAL

      IF (ICAL == 1) THEN

        SDVI1  = 0._DP
        SDVI2  = 0._DP
        SIGMAC = 0._DP
        SGMCS  = 0._DP

        ISDVI  = 0
        IIHC   = 0
        IGHC   = 0

      ELSE IF (ICAL == 2) THEN

        LMETSP = .FALSE.
        LMETSPW = .FALSE.

      END IF

      RETURN
      END SUBROUTINE EIRENE_INIT_CSDVI


      SUBROUTINE EIRENE_BROADCAST_CSDVI(ME)
      USE EIRMOD_MPI
      INTEGER, INTENT(IN) :: ME
      INTEGER :: IER
#if ( defined(USE_MPI) && !defined(GFORTRAN) && !defined(CRAY) )
      EXTERNAL :: MPI_BCAST
#endif

      IF (ME /= 0) THEN
        CALL EIRENE_ALLOC_CSDVI(1)
        CALL EIRENE_ALLOC_CSDVI(2)
      END IF

      CALL MPI_BCAST (ISDVI,MSDVI,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IIHC,2*NCV,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IGHC,2*NCV,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

      RETURN
      END SUBROUTINE EIRENE_BROADCAST_CSDVI

      END MODULE EIRMOD_CSDVI
