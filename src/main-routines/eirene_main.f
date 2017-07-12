Cdr  june 17:  gr-cleanup: call grstrt, grend --> call eirene_plstrt, eirene_plend
C
C     EIRENE VERSION SVN ....  (Jan.2014)... MOVED TO GIT REPOSITORY
C
C
      PROGRAM EIRENE_MAIN

cdr  Main program, to run eirene as stand alone code.

cdr  As an alternative, there are other entry points,
cdr  to run eirene from within other codes, e.g. in iterative mode.
cdr  These calls to eirene_eirene (or entry point eirene_eirene_couple)
cdr  are from main driving interfacing routine EIRSRT.f 
cdr  and the parameters DT, NLM, NLL, ITNR, MPI_INIT... are
cdr  external code and problem specific.


      USE EIRMOD_PRECISION
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      REAL(DP) :: DT
      REAL(DP) :: EIRENE_SECOND_OWN, TIMI, TIMEND
      INTEGER :: ITNR
      LOGICAL :: NLM,NLL,MPI_INIT
C
      TIMI=EIRENE_SECOND_OWN()
C
      CALL EIRENE_PLSTRT
C
C  call eirene, to carry out a "stand alone eirene run"
      NLM=.FALSE.
C  time step set internally from input file
      DT=0._DP
c  last call, deallocate arrays at the end of run
      NLL=.TRUE.
c  iteration number for iterations with external code
      ITNR=1
c  initialize MPI routines
      MPI_INIT=.TRUE.

      CALL EIRENE_EIRENE(DT,NLM,NLL,ITNR,MPI_INIT)
C
      CALL EIRENE_PLEND
C
      TIMEND=EIRENE_SECOND_OWN()
      WRITE (IUNOUT,*) 'TOTAL CPU_TIME OF THIS RUN: ',TIMEND-TIMI
C
      STOP
      END
