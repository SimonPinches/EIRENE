cdr: This routine strictly should be a third party routine,
cdr  because allocation of cpu time to strata may be done according
cdr  different criteria. (load balancing, variance minimization via stratification....)
cdr  The present version of PEDIST is aiming at "proportional allocation",
cdr  See EIRENE manual, "stratified source sampling". 

cdr  currently it is ruled out that one processor deals
cdr  with more than one stratum, except in the serial case (only one processor)
cdr  To generalize this, some coding in MCARLO.f and perhaps elsewhere 
cdr  may need to be adjusted...

!pb  18.12.06: COMPUTATION TIME PER PROCESSOR IS SET TO THE MAXIMUM TIME
!pb            THAT IS AVAILABLE
C
      SUBROUTINE EIRENE_PEDIST (XTIM,XX1)
C  PURPOSE:
C  SET:  PROCFORSTRA(ISTRA,IPE):   IF TRUE: PROCESSOR IPE WORKS ON STRATUM ISTRA
C
C   IF THERE ARE MORE PROCESSORS THAN STRATA:
C   SUBROUTINE PEDIST CALCULATES THE ASSIGNMENT OF PROCESSORS TO
C   STRATA, ACCORDING TO CERTAIN CRITERIA.
C
C  PRESENT VERSION:
C   DISTRIBUTION OF PE'S IS DONE ACCORDING TO THE DISTRIBUTION OF
C   COMPUTATION TIME.

C   IF THERE ARE FEWER PROCESSORS THAN STRATA:
C   CASE A: ONLY ONE PROCESSOR:  ALL STRATA TO THIS SINGLE PROCESSOR
C   CASE B: SEVERAL PROCESSORS:  ASIGN A PROCESSOR TO EACH STRATUM. SOME PROCESSORS
C                                MAY RECEIVE MORE THAN ONE STRATUM.
C                                DO NOT ASSIGN SEVERAL PROCESSORS TO ONE STRATUM
cdr June 17: the last criterion may  be too restricitve 
cdr          and perhaps not be needed either
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CCONA
      USE EIRMOD_CPES
      USE EIRMOD_COMSOU
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CTRCEI, ONLY: TRCCEN
csw 18mar2013
      use EIRMOD_COUTAU
 
      IMPLICIT NONE
 
      REAL(DP), INTENT(INOUT) :: XTIM(0:NSTRA)
      REAL(DP), INTENT(IN) :: XX1
      REAL(DP) :: TIMPE(0:NSTRA), TSTRPE(NSTRA,0:NPRS-1)
      REAL(DP) :: FACP, DELT, SUMTIM, TMEAN, TPE
      INTEGER :: IPE, K, I, ISTRA, NPRS_FREE, NPRS_OPT,n
 
      PROCFORSTRA = .FALSE.

      IF (NPRS == 1) THEN

! 1 PROCESSOR: ALL STRATA ARE DONE BY PROCESSOR 0
!              XTIM REMAINS UNCHANGED
        
        PROCFORSTRA(1:NSTRAI,0) = NLSRON(1:NSTRAI)

      ELSE IF (NPRS <= COUNT(NLSRON(1:NSTRAI))) THEN

! FEWER PROCESSORS THAN STRATA
! ROUND ROBIN DISTRIBUTION OF PROCESSORS
! EACH PROCESSOR CAN CALCULATE SEVERAL STRATA
! BUT EACH STRATUM IS CALCULATED BY EXACTLY ONE PROCESSOR
! ADJUST XTIM TO OPTIMIZE USE OF AVAILABLE CPU TIME       
        TSTRPE = 0._DP
        IPE = -1
        DO ISTRA = 1, NSTRAI
          IF (NLSRON(ISTRA)) THEN
            IPE = IPE + 1
            IF (IPE >= NPRS) IPE = 0
            PROCFORSTRA(ISTRA,IPE) = .TRUE.
            TSTRPE(ISTRA,IPE) = XTIM(ISTRA)
          END IF
        END DO

        sumtim=xtim(0)
        DO IPE = 0, NPRS-1
          TPE = SUM(TSTRPE(1:NSTRAI,IPE))
          FACP = SUMTIM / TPE
          TSTRPE(1:NSTRAI,IPE) = TSTRPE(1:NSTRAI,IPE) * FACP
        END DO

        IPE = -1
        DO ISTRA = 1, NSTRAI
          IF (NLSRON(ISTRA)) THEN
            IPE = IPE + 1
            IF (IPE >= NPRS) IPE = 0
            XTIM(ISTRA) = TSTRPE(ISTRA,IPE)
          ELSE
            XTIM(ISTRA) = 0._DP
          END IF
        END DO

        xtim(0) = sum(xtim(1:nstrai))
        CALL EIRENE_MASAGE
     .    ('REDEFINED CPU TIME ASSIGNED TO STRATA (SEC) :')
        DO ISTRA=1,NSTRAI
          CALL EIRENE_MASJ1R ('STRATUM, TIME   ',ISTRA,XTIM(ISTRA))
        END DO

      ELSE

! calculate mean cpu time per stratum
        sumtim=xtim(0)
        TMEAN=SUMTIM/FLOAT(NPRS)
 
        WRITE (iunout,*) ' SUMTIM = ',SUMTIM,' MEAN TIME = ',TMEAN
 
        NPRS_OPT=0
        NPRS_FREE=NPRS


csw 18mar2013 added branch to test xmct from previous run
C XMCT not stored on fort.11 any more (better place fort.14)
C Without activating fort.11:
C 1st iteration, XMCT == 0
C 2nd iteration, XMCT value of 1st iteration 
C etc.
C When XMCT still was stored on fort.11 and read from fort.11 was 
C active the situation was as follows:
C 1st iteration, XMCT == 0, at the end of MCARLO XMCT was read from 
C   some old fort.11
C 2nd iteration, XMCT used from the "some old fort.11", definitly not 
C   the last
C   => feature broken anyway... 
        if(xmct(0) <= 0.0 ) then
          DO ISTRA=1,NSTRAI
            delt=xtim(istra)
            IF (delt/tmean.GE.1.E-5) THEN
! a stratum that has got computation time gets at least 1 processor
              NPESTR(ISTRA)=1
              NPRS_FREE=NPRS_FREE-1
            ELSE
              NPESTR(ISTRA)=0
            ENDIF
! calculate the optimal number of additional processors according to
! distribution of cpu time done in mcarlo (according to number of particles
! and source strength specified in the input)
            TIMPE(ISTRA)=MAX(delt-TMEAN,0._DP)/TMEAN
            NPRS_OPT=NPRS_OPT+int(TIMPE(ISTRA))
          ENDDO
          WRITE (iunout,*) ' ISTRA, TIMPE '
          DO ISTRA=1,NSTRAI
            WRITE (iunout,*) ISTRA,TIMPE(ISTRA)
          ENDDO
 
          WRITE (iunout,*) ' NPRS_FREE ',NPRS_FREE
        
! distribute free processors to strata by their optimal number of processors
          FACP=MIN(1.D0,REAL(NPRS_FREE,KIND(1.D0))/
     .               (REAL(NPRS_OPT,KIND(1.D0))+eps30))
          write (iunout,*) ' facp ',facp
          NPESTR(0)=NPRS
          DO ISTRA=1,NSTRAI
            NPESTR(ISTRA)=NPESTR(ISTRA)+int(TIMPE(ISTRA)*FACP)
            NPRS_FREE=NPRS_FREE-int(TIMPE(ISTRA)*FACP)
          ENDDO
          WRITE (iunout,*) ' NPESTR ',(NPESTR(ISTRA),ISTRA=1,NSTRAI)
          WRITE (iunout,*) ' NPRS_FREE ',NPRS_FREE


        else

csw attempting better work load balancing           
          npestr(0)=nprs
          tmean=xtim(0)/dble(nprs)
          do istra=1,nstrai
            timpe(istra) = max(xtim(istra)-tmean,0.d0)/tmean
            if(xtim(istra) > 0.) then
             facp=max(1.0, dble(nprs)*xmct(istra)/xmct(0))
             n=int(facp)
             npestr(istra)=n
             nprs_free=nprs_free-n
            else
             npestr(istra)=0
            endif
          enddo

          do istra=1,nstrai
            write(iunout,'(a,2i6,2(1x,e13.6))') 
     .              'XMCT ',istra,npestr(istra),xmct(istra),xmcp(istra)
          enddo
        endif
 
csw 14jul2011
        do while (nprs_free < 0) 
          WRITE (iunout,*) ' NPRS_FREE ',NPRS_FREE
          i=maxloc(npestr(1:nstrai),dim=1)
          npestr(i)=npestr(i)-1
          nprs_free=nprs_free+1
        enddo
csw
 
! if there are still free processors left, distribute them to all
! strata with more than tmean cpu time assigned to them using a
! daisy chain mechanism
        ISTRA=0
        DO WHILE (NPRS_FREE.GT.0)
          ISTRA=ISTRA+1
          IF (ISTRA.GT.NSTRAI) ISTRA=1
          IF (TIMPE(ISTRA).GT.1.E-10) THEN
            NPESTR(ISTRA)=NPESTR(ISTRA)+1
            NPRS_FREE=NPRS_FREE-1
          ENDIF
        ENDDO
        WRITE (iunout,*) ' NPESTR '
        WRITE (iunout,'(12I6)') (NPESTR(ISTRA),ISTRA=1,NSTRAI)
        WRITE (iunout,*) ' NPRS_FREE ',NPRS_FREE

csw 14jul2011
        if(sum(npestr(1:nstrai)) /= npestr(0) ) then
          write(iunout,*) 'pedist: wrong number of processors in npestr'
          call eirene_exit_own(1)
        endif
csw
 
! assign each processor the numbers ISTRA of the strata it shall work on
        IPE=0
        DO ISTRA=1,NSTRAI
          DO K=1,NPESTR(ISTRA)          
            NSTRPE(IPE)=ISTRA
            PROCFORSTRA(ISTRA,IPE) = .TRUE.
            IPE=IPE+1
          ENDDO
        ENDDO
        WRITE (iunout,*) 'pedist:  proc. IPE works on stratum ISTRA '
        WRITE (iunout,*) ' IPE, ISTRA '
        WRITE (iunout,'(12I6)') (I,NSTRPE(I),I=0,NPRS-1)
 
! for each stratum define the number of the first processor NPESTA
! NPESTA(istra) is the "Master processor" for stratum no. ISTRA.

! It does calculations for this stratum.
! This is used to determine the groups of further processors in the
! accumulation of the results for one stratum
        NPESTA(0)=0
        NPESTA(1)=0
        DO ISTRA=2,NSTRAI
          NPESTA(ISTRA)=NPESTA(ISTRA-1)+NPESTR(ISTRA-1)
        ENDDO
        WRITE (iunout,*) ' MASTER PROCESSOR FOR STRATUM '
        WRITE (iunout,*) ' ISTRA, NPESTA '
        WRITE (iunout,'(12I6)') (I,NPESTA(I),I=0,NSTRAI)
 
        XTIM(1:NSTRAI) = XX1
        CALL EIRENE_MASAGE
     .    ('REDEFINED CPU TIME ASSIGNED TO STRATA (SEC) :')
        DO ISTRA=1,NSTRAI
          CALL EIRENE_MASJ1R ('STRATUM, TIME   ',ISTRA,XTIM(ISTRA))
        END DO

      END IF  
 
      RETURN
      END
