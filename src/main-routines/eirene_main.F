Cdr  june 17: gr-cleanup: call grstrt, grend --> call eirene_plstrt, eirene_plend
cdr  comments
cd   oct 20: modified to obtain total cpu time, sum over all iterations.
C
C     EIRENE VERSION SVN ....  (Jan.2014)... MOVED TO GIT REPOSITORY
C
C
      PROGRAM EIRENE_MAIN

cdr  Main program, to run EIRENE as standalone code.
c
c     a)  initialize graphics routines
c     b)  set default global run parameters NLM, DT, NLL, ITNR, MPI_INIT
c         for standalone runs.
c     c)  call EIRENE(DT,NLM,NLL,ITNR,MPI_INIT)
c     d)  close graphics routines

c    NLM     :
c    NLL     :
c    DT      :
c    ITNR    :
c    MPI_INIT:

cdr  As an alternative, there are other entry points into EIRENE,
cdr  to run EIRENE from within other codes, e.g. in iterative mode.
cdr  These calls are then to SUBR. EIRENE (rather than: PROGR. MAIN)
cdr  or to the entry point EIRENE_COUPLE in SUBR. EIRENE.

cdr  For example:
cdr  CALL EIRENE_EIRENE and CALL_EIRENE_COUPLE
cdr  are preprogrammed in the main interfacing routine EIRSRT.f,
cdr  for some frequently used coupled applications (with B2, B2.5, etc..)
cdr  The parameters DT, NLM, NLL, ITNR, MPI_INIT... are then
cdr  set from the external code, or in EIRSRT, and are problem-specific.

      USE EIRMOD_PRECISION
      USE EIRMOD_SECOND_OWN, ONLY: EIRENE_SECOND_OWN,
     .                             EIRENE_RESET_SECOND
      USE EIRMOD_EIRENE, ONLY: EIRENE_EIRENE
      USE EIRMOD_CPES, ONLY: MY_PE
#ifdef USE_MPI
      USE MPI, ONLY: MPI_FINALIZE    ! IGNORE
      USE EIRMOD_MPI, ONLY: EIRENE_MPI_INIT
#else
      USE EIRMOD_MPI, ONLY: EIRENE_MPI_INIT, MPI_FINALIZE
#endif

      IMPLICIT NONE
      REAL(DP) :: DT
      REAL(DP) :: TIMI, TIMEND
      REAL(DP), SAVE :: TIMI0
      INTEGER :: ITNR, IER
      LOGICAL :: NLM, NLL, MPI_INIT
      EXTERNAL :: EIRENE_PLSTRT, EIRENE_PLEND

C     Initialize MPI
      CALL EIRENE_MPI_INIT(IER)
C
      TIMI=EIRENE_SECOND_OWN()
      TIMI0=TIMI
c     write (*,*) 'timi ',timi0
C
      CALL EIRENE_PLSTRT
C
C  call eirene, to carry out a "standalone eirene run"
      NLM=.FALSE.
C  time step set later internally from input file
      DT=0._DP
c  last call, deallocate arrays at the end of run
      NLL=.TRUE.
c  iteration number for iterations with external code
      ITNR=1
c  initialize MPI routines
      MPI_INIT=.FALSE.
!!! debug / commented until -cpp option added
!!!      write(6,*) "EIRENE_MAIN: starting"
!!!#ifdef USE_EXT_OPENMP
!!!      write(6,*) "Entering parallel region replicating external region"
!!!$OMP PARALLEL DEFAULT(SHARED)
!!!      write(6,*) "Entered parallel region replicating external region"
!!!#endif
      CALL EIRENE_EIRENE(DT,NLM,NLL,ITNR,MPI_INIT)
!!!#ifdef USE_EXT_OPENMP
!!!!$OMP END PARALLEL
!!!#endif
C
      CALL EIRENE_PLEND
C
      TIMEND=EIRENE_RESET_SECOND()
      IF (MY_PE == 0) WRITE (*,*) 'TOTAL CPU_TIME (SEC) OF THIS RUN: ',
     .                  TIMEND-TIMI0
C     Finalize MPI
      CALL MPI_FINALIZE(IER)
C
      STOP
      END PROGRAM EIRENE_MAIN
