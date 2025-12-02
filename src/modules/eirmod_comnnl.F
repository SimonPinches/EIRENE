!pb  31.10.06:  definition of census arrays RPART, RPARTC, IPART, IPARTC changed
!               RPART (NPRNL,NPARTT) --> RPART (NPARTT,NPRNL)
!               IPART (MPRNL,NPARTT) --> IPART (MPARTT,NPRNL)
!               RPARTC(NPRNL,NPARTT) --> RPARTC(NPARTT,NPRNL)
!               IPARTC(MPRNL,NPARTT) --> IPARTC(MPARTT,NPRNL)
cdr:  rpart, ipart: "true census" arrays, real and integer.
cdr:  rpartc,ipartc: copies of census arrays, needed for sampling (bootstrapping)
cdr                  from old census, while already filling the new census
cdr   rpartw: cumulated weight from census, set at the end of timestep (TMSTEP.f)
cdr           from the weights stored on rpart, for bootstrapping in next timestep
cdr           (in subr. LOCATE.f)
cdr Jan 2020:  fix nonlinear iteration problems.
cdr            In case of short version of FORT.13:
cdr              Other data than just the target PLASMA_BCKGRND input tallies
cdr              for virtual species IPLS: NPLS_FIX+1:NPLSI are needed:
cdr              Add: reaction rates, tracklength estimator rates, etc.
cdr              derived tallies: lgvac(ipls,.): ipls=npls_fix+1,npls
cdr              and maybe more.
cdr              Same for "semi-bgk" photon transport iterations.
cdr              But we can remove the B field, E field, etc...
cdr            In case of long version of FORT.13:
cdr              we could run cross-collisions with "wrong" DI, TI,
cdr              but correct LGVAC, TBEL3,...
cdr              or the other way round, but not both.


      MODULE EIRMOD_COMNNL

cdr  Some variables must become threadprivate: e.g. itmstp,
cdr  and the count no. of scores on census?

      USE EIRMOD_PRECISION, ONLY: DP
      USE EIRMOD_PARMMOD, ONLY: IUNMEM, MPARTT, NPARTT, NPRNL, NSTRA,
     .                          NREI, NRCX, NREL, NRPI, NREC, NPLS, NRAD


      IMPLICIT NONE

      PRIVATE

      PUBLIC :: EIRENE_ALLOC_COMNNL, EIRENE_DEALLOC_COMNNL,
     P          EIRENE_INIT_COMNNL, EIRENE_BROADCAST_COMNNL

cdr  data for iterative mode, e.g. bgk, photon transport, ...
      INTEGER, PUBLIC, SAVE ::
     I  NFLA_VIRT, NREA_VIRT, NREA_MAX
      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I  NFLA_IPLS(:), NFLA_ISWR(:), NFLA_IR(:)

      LOGICAL, PUBLIC, ALLOCATABLE, SAVE ::
     L  LG_STORE(:,:)
      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R  TAB_STORE(:,:,:), E_STORE(:,:,:)

cdr  data for time dep. mode, e.g. census flux, time stepping, census arrays
      REAL(DP), PUBLIC, SAVE ::
     R  FLXCEN, DTIMV,  DTIMVI, DTIMVN, TIME0
cdr  census arrays, containing the full state vectors of all scores at census (=end of timestep)
cdr  real(dp) coordinates
      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R  RPART(:,:), RPARTC(:,:),
cdr  derived: cumulative distribution of indices for re-sampling (bootstrapping) from census
     R  RPARTW(:)

      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
cdr  census arrays, containing the full state vectors of all scores at census (=end of timestep)
cdr  integer coordinates
     I  IPART(:,:), IPARTC(:,:),
     I  NPRNLS(:)

      INTEGER, PUBLIC, SAVE ::
     I  NPRNLI, IPRNLI, IPRNLS, IPRNL,
     I  NPTST,  NTMSTP, ITMSTP

      INTEGER, PUBLIC, SAVE ::
     I  NPRNLI_IN, NINITL_READ, NPRMUL

      CONTAINS


      SUBROUTINE EIRENE_ALLOC_COMNNL(IPRT)
      INTEGER, INTENT(IN) :: IPRT
cdr Jan 2020: added: part 1, for bgk-type and photon iteration data structures
cdr Part 1: nonlinear bgk (and photon) iterations
cdr Part 2: time-dep mode: census arrays

cdr 1.) deal with nonlinear (BGK) iteration data structures
cdr     (also: nonlinear coupling to radiation field)
cpb can only be called after call to input
      IF (IPRT == 1) THEN
         IF (.NOT. ALLOCATED(NFLA_ISWR)) THEN
            ALLOCATE (NFLA_IPLS(NPLS))
            NREA_MAX=NREI+NRCX+NREL+NRPI+NREC
            ALLOCATE (NFLA_ISWR(NREA_MAX))
            ALLOCATE (NFLA_IR(NREA_MAX))

            WRITE (IUNMEM,'(A,T25,I15)')
     .           ' COMNNL(1) ',(NPLS+2*NREA_MAX)*4

            CALL EIRENE_INIT_COMNNL(1)
         ENDIF
      END IF

cdr 2.) deal with census array (t-dep mode) data structures
      IF (IPRT == 2) THEN
        IF (.NOT.ALLOCATED(RPART)) THEN

          ALLOCATE (RPART(NPARTT,NPRNL))
          ALLOCATE (RPARTC(NPARTT,NPRNL))
          ALLOCATE (RPARTW(0:NPRNL))

          ALLOCATE (IPART(MPARTT,NPRNL))
          ALLOCATE (IPARTC(MPARTT,NPRNL))
          ALLOCATE (NPRNLS(NSTRA))

          WRITE (IUNMEM,'(A,T25,I15)')
     .      ' COMNNL(2) ',(2*NPRNL*NPARTT+NPRNL+1)*8 +
     .                    (2*NPRNL*MPARTT+NSTRA)*4

          CALL EIRENE_INIT_COMNNL(2)
         ENDIF
      END IF

      RETURN
      END SUBROUTINE EIRENE_ALLOC_COMNNL


      SUBROUTINE EIRENE_DEALLOC_COMNNL
cdr  part 1: nonlinear iterations, bgk, photons,...

      IF (ALLOCATED(NFLA_ISWR)) THEN
        DEALLOCATE (NFLA_IPLS)
        DEALLOCATE (NFLA_ISWR)
        DEALLOCATE (NFLA_IR)
      ENDIF

      IF (ALLOCATED(LG_STORE)) DEALLOCATE(LG_STORE)
      IF (ALLOCATED(TAB_STORE)) DEALLOCATE(TAB_STORE)
      IF (ALLOCATED(E_STORE)) DEALLOCATE(E_STORE)

cdr  part 2: time dep mode, census:
      IF (ALLOCATED(RPART)) THEN
        DEALLOCATE (RPART)
        DEALLOCATE (RPARTC)
        DEALLOCATE (RPARTW)

        DEALLOCATE (IPART)
        DEALLOCATE (IPARTC)
        DEALLOCATE (NPRNLS)
      ENDIF

      RETURN
      END SUBROUTINE EIRENE_DEALLOC_COMNNL


      SUBROUTINE EIRENE_INIT_COMNNL(iprt)
      INTEGER, INTENT(IN) :: IPRT

      if (iprt.eq.1) then
cdr  Part 1: nonlinear BGK (and photon) iterations
cdr           Modified (iterated) model parameters, e.g. reaction rates, etc.
        NFLA_VIRT = 0

        NFLA_IPLS = 0
        NFLA_ISWR = 0
        NFLA_IR = 0

      elseif (iprt.eq.2) then
cdr  Part 2: time dep mode, census arrays, etc.
        RPART  = 0._DP
        RPARTC = 0._DP
        RPARTW = 0._DP
        FLXCEN = 0._DP

        IPART  = 0
        IPARTC = 0

        NPRNLI = 0
        NPRNLS = 0
        IPRNLI = 0
        IPRNLS = 0
        IPRNL  = 0
        NPTST  = 0
        NTMSTP = 0
        ITMSTP = 0

      endif

      RETURN
      END SUBROUTINE EIRENE_INIT_COMNNL


C> \brief Broadcast quantities of module EIRMOD_COMNNL
C>
C> All quantities of the module EIRMOD_COMNNL are broadcasted here.
C> The census arrays are only broadcasted in the parallelisation mode
C> with proportional allocation. Especially, in the embarrassingly
C> parallel mode, each process should keep its own census. (Attention:
C> write out of the census of each process for the restart of a run is
C> not (yet) implemented and may be part of the plasma code interface.)
      SUBROUTINE EIRENE_BROADCAST_COMNNL(ME)

      USE EIRMOD_PARMMOD, ONLY: MPARTT, NPARTT, NPRNL, NSTRA
      USE EIRMOD_MPI

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ME
      INTEGER :: IER
#if ( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) && !defined(CRAY) )
      EXTERNAL :: MPI_BCAST
#endif

      IF (ME /= 0) THEN
        CALL EIRENE_ALLOC_COMNNL(1)
        CALL EIRENE_ALLOC_COMNNL(2)
      END IF

      CALL MPI_BCAST (NFLA_VIRT,1,MPI_INTEGER,0,MPI_COMM_WORLD,IER)
      CALL MPI_BCAST (NREA_VIRT,1,MPI_INTEGER,0,MPI_COMM_WORLD,IER)
      CALL MPI_BCAST (NREA_MAX,1,MPI_INTEGER,0,MPI_COMM_WORLD,IER)
      CALL MPI_BCAST (NFLA_IPLS,NPLS,MPI_INTEGER,0,MPI_COMM_WORLD,IER)
      CALL MPI_BCAST (NFLA_ISWR,NREA_MAX,MPI_INTEGER,
     >                                           0,MPI_COMM_WORLD,IER)
      CALL MPI_BCAST (NFLA_IR,NREA_MAX,MPI_INTEGER,0,MPI_COMM_WORLD,IER)

      IF (ME /= 0) THEN
        IF (NREA_VIRT.GT.0 .AND. .NOT.ALLOCATED(TAB_STORE) ) THEN
          ALLOCATE (TAB_STORE(NREA_VIRT,NRAD,9))
          ALLOCATE (E_STORE(NREA_VIRT,NRAD,9))
          TAB_STORE = 0._DP
          E_STORE = 0._DP
        END IF
        IF (NFLA_VIRT.GT.0 .AND. .NOT.ALLOCATED(LG_STORE) ) THEN
          ALLOCATE (LG_STORE(NRAD,NFLA_VIRT))
          LG_STORE = .FALSE.
        END IF
      END IF

      IF (NFLA_VIRT.GT.0) THEN
        CALL MPI_BCAST (LG_STORE,NRAD*NFLA_VIRT,MPI_LOGICAL,
     >                                           0,MPI_COMM_WORLD,IER)
      END IF
      IF (NREA_VIRT.GT.0) THEN
        CALL MPI_BCAST (E_STORE,NREA_VIRT*NRAD*9,MPI_REAL8,
     >                                           0,MPI_COMM_WORLD,IER)
        CALL MPI_BCAST (TAB_STORE,NREA_VIRT*NRAD*9,MPI_REAL8,
     >                                           0,MPI_COMM_WORLD,IER)
      END IF

      CALL MPI_BCAST (FLXCEN,1,MPI_REAL8,0,MPI_COMM_WORLD,IER)
      CALL MPI_BCAST (DTIMV,1,MPI_REAL8,0,MPI_COMM_WORLD,IER)
      CALL MPI_BCAST (DTIMVI,1,MPI_REAL8,0,MPI_COMM_WORLD,IER)
      CALL MPI_BCAST (DTIMVN,1,MPI_REAL8,0,MPI_COMM_WORLD,IER)
      CALL MPI_BCAST (TIME0,1,MPI_REAL8,0,MPI_COMM_WORLD,IER)
      CALL MPI_BCAST (NPRNLI,1,MPI_INTEGER,0,MPI_COMM_WORLD,IER)
      CALL MPI_BCAST (IPRNLI,1,MPI_INTEGER,0,MPI_COMM_WORLD,IER)
      CALL MPI_BCAST (IPRNLS,1,MPI_INTEGER,0,MPI_COMM_WORLD,IER)
      CALL MPI_BCAST (IPRNL ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IER)
      CALL MPI_BCAST (NPTST ,1,MPI_INTEGER,0,MPI_COMM_WORLD,IER)
      CALL MPI_BCAST (NTMSTP,1,MPI_INTEGER,0,MPI_COMM_WORLD,IER)
      CALL MPI_BCAST (ITMSTP,1,MPI_INTEGER,0,MPI_COMM_WORLD,IER)

      IF ( NPRNL > 0 ) THEN
        CALL MPI_BCAST (RPART,NPRNL*NPARTT,MPI_REAL8,
     >                  0,MPI_COMM_WORLD,IER)
        CALL MPI_BCAST (RPARTC,NPRNL*NPARTT,MPI_REAL8,
     >                  0,MPI_COMM_WORLD,IER)
        CALL MPI_BCAST (RPARTW,NPRNL+1,MPI_REAL8,0,MPI_COMM_WORLD,IER)

        CALL MPI_BCAST (IPART,NPRNL*MPARTT,MPI_INTEGER,
     >                  0,MPI_COMM_WORLD,IER)
        CALL MPI_BCAST (IPARTC,NPRNL*MPARTT,MPI_INTEGER,
     >                  0,MPI_COMM_WORLD,IER)
        CALL MPI_BCAST (NPRNLS,NSTRA,MPI_INTEGER,0,MPI_COMM_WORLD,IER)
      END IF

      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
      RETURN

      END SUBROUTINE EIRENE_BROADCAST_COMNNL

      END MODULE EIRMOD_COMNNL
