Cdr  june 17:  gr-cleanup: call grstrt, grend --> call eirene_plstrt, eirene_plend
cdr  comments
C
C     EIRENE VERSION SVN ....  (Jan.2014)... MOVED TO GIT REPOSITORY
C
C
      PROGRAM EIRENE_MAIN

cdr  Main program, to run eirene as stand alone code.
c
c     a)  initialize graphics routines
c     b)  set default global run parameters NLM, DT, NLL, ITNR, MPI_INIT
c         for stand alone runs.
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
cdr  are preprogrammed in the main interfacing routine EIRSRT.f ,
cdr  for some frequently used coupled applications (with B2, B2.5, etc..)
cdr  The parameters DT, NLM, NLL, ITNR, MPI_INIT... are then
cdr  set from the external code, or in EIRSRT, and are problem specific.


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
