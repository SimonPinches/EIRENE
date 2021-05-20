!pb  30.10.06:  XNUE removed
cdr  sept. 2015: npartt=11, rather than 12 (xgener not stored on census)
cdr  april 2017: some cleanup carried over from solps_iter branch
cdr  nov.17    : dead flag: nlstor, (and call store...) now removed
cdr  oct.19    :  started to remove unused variables, ..._mean, stemis, etc.

c.........................................................................
c
c  COMPRT contains particle coordinates along track.
c  The full information for restart after splitting is contained in the
c  npartc (real) and mpartc (integer) variables.
c
c  A reduced set for restart from a census array (initial condition in time)
c  is contained in the smaller
c  NPARTT (real) and MPARTT (integer) variables.

c  mpartc, npartc and mpartt, npartt are set in eirmod_parmmod


      MODULE EIRMOD_COMPRT

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: EIRENE_ALLOC_COMPRT, EIRENE_DEALLOC_COMPRT,
     P          EIRENE_INIT_COMPRT

cdr..............................................
cdr   event-type is currently unused (still set in locate). Formerly: characterize
cdr   previous event, e.g. for modified scoring.
cdr   For example: removal of emission and
cdr   absorption within same cell from scoring.

      PUBLIC :: EVENT_TYPE

      TYPE :: EVENT_TYPE
        INTEGER :: NCELL, ITYP, ISPEZ, IFLAG
        REAL(DP) :: E0, WEIGHT
      END TYPE EVENT_TYPE

      TYPE(EVENT_TYPE), PUBLIC, SAVE :: LAST_EVENT
cdr   later use: last_event%e0, last_event%weight ...
cdr..............................................

      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE :: RPST(:)

      REAL(DP), PUBLIC, POINTER, SAVE :: RPSTT(:)

C NPARTT PARTICLE COORDINATES FOR CENSUS ARRAY
C NPARTC PARTICLE COORDINATES, REAL, (E.G.: SPLITTING)
C NPARTT AND NPARTC ARE SET IN EIRMOD_PARMMOD, CURRENTLY:
C NPARTT=11
C NPARTC=12
      REAL(DP), PUBLIC, POINTER, SAVE ::
     R X0,     Y0,     Z0,
     R VEL,    VELX,   VELY,   VELZ,
     R E0,     WEIGHT, TIME,   PHI,  ! UP TO HERE: STORE ON CENSUS
     R XGENER  ! UP TO HERE: STORE FULL PARTICLE INFORMATION

C  SOME FURTHER REAL VARIABLES USED ALONG PARTICLE TRAJECTORY

      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R TIMINT(:), TIMPOL(:,:)

      REAL(DP), PUBLIC, SAVE ::
     R TL,     TT,     TS,     TF,     ZT,         ZDT1,
     R CRTX,   CRTY,   CRTZ,   SCOS,   SCOS_SAVE,  WGHTSP, WGHTSC,
     R CRTXG,  CRTYG,  CRTZG

c............................................................
cdr the following variables are unused. Left over from former, mostly
cdr unfinished options. To be removed?

      REAL(DP), PUBLIC, SAVE :: DE0_RAYL, DE0_RAYR
      REAL(DP), PUBLIC, ALLOCATABLE, SAVE :: E0_RAY(:)

c.............................................................
      INTEGER, PUBLIC, TARGET, ALLOCATABLE, SAVE :: IPSTD(:)

      INTEGER, PUBLIC, POINTER, SAVE :: IPST(:), IPSTT(:)

C MPARTT PARTICLE COORDINATES, REDUCED SET FOR CENSUS ARRAY
C MPARTC PARTICLE COORDINATES, FULL SET, INTEGER, (E.G.: SPLITTING)
C MPARTT AND MPARTC ARE SET IN EIRMOD_PARMMOD, CURRENTLY:
C MPARTT=10
C MPARTC=14

      INTEGER, PUBLIC, POINTER, SAVE ::
     I NPANU,
     I IPOLG,  IPERID, NCELL,
     I ITIME,  IFPATH, IUPDTE,
     I ISTRA,
     I ISPZ,  ! UP TO HERE: STORE ON CENSUS
     I MRSURF, MPSURF, MTSURF, MASURF, MSURF,  ! UP TO HERE: STORE FULL PARTICLE INFORMATION
     I MSURFG
!$OMP THREADPRIVATE(ISTRA)
C  SOME FURTHER INTEGER VARIABLES USED ALONG PARTICLE TRAJECTORY

      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I NTIM(:), IIMPOL(:,:), IIMINT(:)

      INTEGER, PUBLIC, SAVE ::
     I NRCELL, NPCELL, NTCELL, NACELL, NBLOCK, NBLCKA, NSTCLL,
     I IC_NEUT, IC_ION,
     I ITYP,   IATM,   IMOL,   IION,   IPLS,   IPHOT,
     I ICOL,   IPOLGN, NINCX,  NINCY,  NINCZ,  NINCA,  NJUMP,
     I NIMINT, ITRJ


!     VARIABLES FOR UNIFIED SUBROUTINES
      INTEGER, PUBLIC, SAVE ::
     I IXSPZ, NMETOFF

      LOGICAL, PUBLIC, SAVE ::
     L LGPART, LGLAST, LGTIME,
     L NLSRFX, NLSRFY, NLSRFZ, NLSRFA,
     L NLTRC,  NLTRJ

!$OMP  THREADPRIVATE(RPSTT,X0,Y0,Z0,VEL,VELX,VELY,VELZ,E0,WEIGHT,
!$OMP& IPSTD,RPST,
!$OMP& TIME,PHI,XGENER,TIMINT,TIMPOL,TL,TT,TS,TF,ZT,ZDT1,CRTX,CRTY,
!$OMP& CRTZ,SCOS,SCOS_SAVE,WGHTSP,WGHTSC,CRTXG,CRTYG,CRTZG,
!$OMP& DE0_RAYL,
!$OMP& DE0_RAYR,E0_RAY,IPST,IPSTT,NPANU,IPOLG,IPERID,NCELL,
!$OMP& ITIME,IFPATH,IUPDTE,ISPZ,MRSURF,MPSURF,MTSURF,MASURF,MSURF,
!$OMP& MSURFG,NTIM,IIMPOL,IIMINT,NRCELL,NPCELL,NTCELL,NACELL,NBLOCK,
!$OMP& NBLCKA,NSTCLL,IC_NEUT,IC_ION,ITYP,IATM,IMOL,IION,IPLS,IPHOT,ICOL,
!$OMP& IPOLGN,NINCX,NINCY,NINCZ,NINCA,NJUMP,NIMINT,ITRJ,IXSPZ,NMETOFF,
!$OMP& LGPART,LGLAST,LGTIME,NLSRFX,NLSRFY,NLSRFZ,NLSRFA,NLTRC,NLTRJ)


    

c  unrelated to particle trajectories:  IO streams
      INTEGER, PUBLIC, SAVE ::
     I IUNIN,  IUNOUT, IVTKOUT
cym
!!!$OMP THREADPRIVATE(IUNOUT)

      DATA IUNIN / 1 /  ! must be known already during compile time.
c                       ! better: move iunin, iunout, etc.. to parmmod ??




      CONTAINS


      SUBROUTINE EIRENE_ALLOC_COMPRT(NPRS)

      INTEGER, INTENT(IN) :: NPRS

      IF (ALLOCATED(RPST)) RETURN

      ALLOCATE (RPST(NPARTC))
      ALLOCATE (IPSTD(MPARTC+1))

      ALLOCATE (TIMINT(NRADS))
      ALLOCATE (TIMPOL(N1STS,N2NDPLGS))
      ALLOCATE (NTIM(NRADS))
      ALLOCATE (IIMPOL(N1STS,N2NDPLGS))
      ALLOCATE (IIMINT(NRADS))

      WRITE (IUNMEM,'(A,T25,I15)')
     .      ' COMPRT ',(NPARTC+NRADS+N1STS*N2NDPLGS)*8 +
     .                 (MPARTC+1+NRADS+N1STS*N2NDPLGS+NRADS)*4

      RPSTT => RPST      !  full (1: npartc) particle information, real

      X0     => RPST( 1)
      Y0     => RPST( 2)
      Z0     => RPST( 3)
      VEL    => RPST( 4)
      VELX   => RPST( 5)
      VELY   => RPST( 6)
      VELZ   => RPST( 7)
      E0     => RPST( 8)
      WEIGHT => RPST( 9)
      TIME   => RPST(10)
      PHI    => RPST(11)
c  up to here: for census, npartt
      XGENER => RPST(12)
c  up to here: for splitting, npartc

      IPST  => IPSTD(2:MPARTC+1)  !  full (2:mpartc+1) particle information, integer

      IPSTT => IPSTD(1:MPARTT)    !  reduced (1:mpartt), for census

      NPANU  => IPSTD( 1)
      IPOLG  => IPSTD( 2)
      IPERID => IPSTD( 3)
      NCELL  => IPSTD( 4)
      ITIME  => IPSTD( 5)
      IFPATH => IPSTD( 6)
      IUPDTE => IPSTD( 7)
      ISTRA  => IPSTD( 8)
      ISPZ   => IPSTD( 9)
      MRSURF => IPSTD(10)
c  up to here: for census, mpartt
      MPSURF => IPSTD(11)
      MTSURF => IPSTD(12)
      MASURF => IPSTD(13)
      MSURF  => IPSTD(14)
c  up to here: for splitting, mpartc
      MSURFG => IPSTD(15)

      CALL EIRENE_INIT_COMPRT (NPRS)

      RETURN
      END SUBROUTINE EIRENE_ALLOC_COMPRT


      SUBROUTINE EIRENE_DEALLOC_COMPRT

      IF (.NOT.ALLOCATED(RPST)) RETURN

      DEALLOCATE (RPST)
      DEALLOCATE (IPSTD)

      DEALLOCATE (TIMINT)
      DEALLOCATE (TIMPOL)
      DEALLOCATE (NTIM)
      DEALLOCATE (IIMPOL)
      DEALLOCATE (IIMINT)

      RETURN
      END SUBROUTINE EIRENE_DEALLOC_COMPRT


      SUBROUTINE EIRENE_INIT_COMPRT (NPRS)
cdr:  called from ??
cdr   purpose ??
C NPRS: NUMBER OF COMPUTE THREADS IN MPI PARALLEL MODE
C     needed for IUNOUT

      INTEGER, INTENT(IN) :: NPRS

      RPST   = 0._DP
      IPSTD  = 0

      TIMINT = 0._DP
      TIMPOL = 0._DP
      NTIM   = 0
      IIMPOL = 0
      IIMINT = 0

      TL     = 0._DP
      TT     = 0._DP
      TS     = 0._DP
      TF     = 0._DP
      ZT     = 0._DP
      ZDT1   = 0._DP
      CRTX   = 0._DP
      CRTY   = 0._DP
      CRTZ   = 0._DP
      SCOS   = 0._DP
      SCOS_SAVE  = 0._DP
      WGHTSP = 0._DP
      WGHTSC = 0._DP
      CRTXG  = 0._DP
      CRTYG  = 0._DP
      CRTZG  = 0._DP

      NRCELL = 0
      NPCELL = 0
      NTCELL = 0
      NACELL = 0
      NBLOCK = 0
      NBLCKA = 0
      NSTCLL = 0
      IC_NEUT= 0
      IC_ION = 0
      ITYP   = 0
      IATM   = 0
      IMOL   = 0
      IION   = 0
      IPLS   = 0
      IPHOT  = 0
      ICOL   = 0
      IPOLGN = 0
      NINCX  = 0
      NINCY  = 0
      NINCZ  = 0
      NINCA  = 0
      NJUMP  = 0

      NIMINT = 0
      ITRJ   = 0

      LGPART = .FALSE.
      LGLAST = .FALSE.
      LGTIME = .FALSE.
      NLSRFX = .FALSE.
      NLSRFY = .FALSE.
      NLSRFZ = .FALSE.
      NLSRFA = .FALSE.
      NLTRC  = .FALSE.
      NLTRJ  = .FALSE.

      LAST_EVENT%IFLAG  = 0
      LAST_EVENT%NCELL  = 0
      LAST_EVENT%ITYP   = 0
      LAST_EVENT%ISPEZ  = 0
      LAST_EVENT%E0     = 0._DP
      LAST_EVENT%WEIGHT = 0._DP

      DE0_RAYL = 0._DP
      DE0_RAYR = 0._DP

c  io files
cdr same code as in subr. EIRENE
      IUNOUT = 6 + IFOFF
      IF (NPRS > 1) IUNOUT = 7 + IFOFF

      IVTKOUT= 28


      RETURN
      END SUBROUTINE EIRENE_INIT_COMPRT

      END MODULE EIRMOD_COMPRT
