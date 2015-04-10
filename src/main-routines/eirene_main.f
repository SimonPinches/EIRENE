C
C
C     EIRENE VERSION SVN ....  (Jan.2014)
C
C
      PROGRAM EIRENE_MAIN
      USE EIRMOD_PRECISION
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      REAL(DP) :: DT
      REAL(DP) :: EIRENE_SECOND_OWN, TIMI, TIMEND
      INTEGER :: ITNR
      LOGICAL :: NLM,NLL
C
      TIMI=EIRENE_SECOND_OWN()
C
      CALL GRSTRT(35,8)
C 
C  call eirene, to carry out a "stand alone eirene run"
      NLM=.FALSE.
C  time step set internally from input file
      DT=0._DP
c  last call, deallocate arrays at the end of run
      NLL=.TRUE.
c  iteration number for iterations with external code
      ITNR=1
      CALL EIRENE_EIRENE(DT,NLM,NLL,ITNR,.TRUE.)
C
      CALL GREND
C
      TIMEND=EIRENE_SECOND_OWN()
      WRITE (IUNOUT,*) 'TOTAL CPU_TIME OF THIS RUN: ',TIMEND-TIMI
C
      STOP
      END
