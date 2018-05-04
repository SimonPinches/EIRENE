!pb  18.12.06: COMPUTATION TIME PER PROCESSOR IS SET TO THE MAXIMUM TIME
!pb            THAT IS AVAILABLE
C
C> \brief Allocation of MPI processes to strata
C>
C> Allocates the available MPI processes to strata.
C> Allocation of CPU time to strata may be done according different 
C> criteria (load balancing, variance minimization via stratification, 
C> ...). Two standard techniques are implemented in EIRENE and can be 
C> controlled via the block 1 input parameter NPRLL.
C> The default set-up is a simple "embarrassingly parallel" schema as 
C> typical for Monte Carlo codes. An more advanced method using a 
C> proportional allocation (NPRLL == 1) to attempting load balancing 
C> when applying stratification is also available.
C> Furthermore, a user defined set-up (subroutine EIRENE_PEDIST_USR) 
C> can be used (NPRLL == -1).
C>
C> Within this subroutine three arrays are set that define the entrie 
C> parallelisation of EIRENE.
C> - PROCFORSTRA(ISTRA,IPE):   if .TRUE.: process IPE works on stratum ISTRA
C> - NPESTR(ISTRA): number of processes calculating stratum ISTRA
C> - NPESTA(ISTRA): master process for stratum ISTRA
      SUBROUTINE EIRENE_PEDIST (XTIM,XX1)
      USE EIRMOD_PRECISION, ONLY: DP
      USE EIRMOD_COMUSR, ONLY: NPRLL
      USE EIRMOD_PARMMOD, ONLY: NSTRA

      IMPLICIT NONE

      REAL(DP), INTENT(INOUT) :: XTIM(0:NSTRA) !< time allocated for stratum
      REAL(DP), INTENT(IN) :: XX1 !< remaining CPU time

      SELECT CASE( NPRLL )
        CASE( -1 )
          CALL EIRENE_PEDIST_USR( XTIM, XX1 )
        CASE( 0 )
          CALL EIRENE_PEDIST_EMBPARALL
        CASE( 1 )
          CALL EIRENE_PEDIST_PROPALLOC( XTIM, XX1 )
      END SELECT

      CONTAINS

C> \brief "Embarrassingly parallel" schema.
C>
C> Here, all strata are calculated by all processes. XTIM remains
C> unchanged. 
      SUBROUTINE EIRENE_PEDIST_EMBPARALL
      USE EIRMOD_COMSOU, ONLY: NLSRON
      USE EIRMOD_CPES, ONLY: NPESTA, NPESTR, NPRS, PROCFORSTRA
      USE EIRMOD_PARMMOD, ONLY: NSTRA

      IMPLICIT NONE

      INTEGER :: IPE

      PROCFORSTRA = .FALSE.

      DO IPE = 0, NPRS-1
        PROCFORSTRA(1:NSTRA,IPE) = NLSRON
      END DO
      NPESTA = 0
      NPESTR = NPRS

      RETURN
      END SUBROUTINE EIRENE_PEDIST_EMBPARALL

C> \brief Proportional allocation schema.
C>
C> Here, the aim is a "proportional allocation", see EIRENE manual, 
C> "stratified source sampling". 
C>
C> As long as there are more processes than strata the distribution of 
C> processes to strata is done according to the distribution of 
C> computation time.
C>
C> If there are fewer processes than strata:
C> - Case A: only one process: all strata to this single process
C> - Case B: several processes: asign a process to each stratum. Some 
C>   processes may receive more than one stratum. Do not assign several 
C>   processes to one stratum.
      SUBROUTINE EIRENE_PEDIST_PROPALLOC( XTIM, XX1 )
      USE EIRMOD_PRECISION, ONLY: DP
      USE EIRMOD_PARMMOD, ONLY: NSTRA
      USE EIRMOD_CCONA, ONLY: EPS30
      USE EIRMOD_CPES, ONLY: NPESTA, NPESTR, NPRS, PROCFORSTRA
      USE EIRMOD_COMSOU, ONLY: NLSRON, NPTS
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_COUTAU, ONLY: XMCT, XMCP
 
      IMPLICIT NONE
 
      REAL(DP), INTENT(INOUT) :: XTIM(0:NSTRA) !< time allocated for stratum
      REAL(DP), INTENT(IN) :: XX1 !< remaining CPU time
      REAL(DP) :: TIMPE(0:NSTRA), TSTRPE(NSTRA,0:NPRS-1)
      REAL(DP) :: FACP, DELT, SUMTIM, TMEAN, TPE
      INTEGER :: IPE, K, I, ISTRA, NPRS_FREE, NPRS_OPT,n
      INTEGER, DIMENSION(1) :: NSTRPE(0:NPRS-1)
 
      PROCFORSTRA = .FALSE.

      IF (NPRS == 1) THEN

! 1 PROCESSOR: ALL STRATA ARE DONE BY PROCESSOR 0
!              XTIM REMAINS UNCHANGED
        
        PROCFORSTRA(:,0) = NLSRON
        NPESTA = 0
        NPESTR = 1

      ELSE IF (NPRS <= COUNT(NLSRON)) THEN

! FEWER PROCESSORS THAN STRATA
! ROUND ROBIN DISTRIBUTION OF PROCESSORS
! EACH PROCESSOR CAN CALCULATE SEVERAL STRATA
! BUT EACH STRATUM IS CALCULATED BY EXACTLY ONE PROCESSOR
! ADJUST XTIM TO OPTIMIZE USE OF AVAILABLE CPU TIME       
        NPESTR = 1
        TSTRPE = 0._DP
        IPE = -1
        DO ISTRA = 1, NSTRA
          IF (NLSRON(ISTRA)) THEN
            IPE = IPE + 1
            IF (IPE >= NPRS) IPE = 0
            PROCFORSTRA(ISTRA,IPE) = .TRUE.
            NPESTA(ISTRA) = IPE
            TSTRPE(ISTRA,IPE) = XTIM(ISTRA)
          END IF
        END DO

        sumtim=xtim(0)
        DO IPE = 0, NPRS-1
          TPE = SUM(TSTRPE(:,IPE))
          FACP = SUMTIM / TPE
          TSTRPE(:,IPE) = TSTRPE(:,IPE) * FACP
        END DO

        IPE = -1
        DO ISTRA = 1, NSTRA
          IF (NLSRON(ISTRA)) THEN
            IPE = IPE + 1
            IF (IPE >= NPRS) IPE = 0
            XTIM(ISTRA) = TSTRPE(ISTRA,IPE)
          ELSE
            XTIM(ISTRA) = 0._DP
          END IF
        END DO

        xtim(0) = sum(xtim(1:nstra))
        CALL EIRENE_MASAGE
     .    ('REDEFINED CPU TIME ASSIGNED TO STRATA (SEC) :')
        DO ISTRA=1,NSTRA
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
        if(xmct(0) <= 0.0 ) then
          DO ISTRA=1,NSTRA
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
          DO ISTRA=1,NSTRA
            WRITE (iunout,*) ISTRA,TIMPE(ISTRA)
          ENDDO
 
          WRITE (iunout,*) ' NPRS_FREE ',NPRS_FREE
        
! distribute free processors to strata by their optimal number of processors
          FACP=MIN(1.D0,REAL(NPRS_FREE,KIND(1.D0))/
     .               (REAL(NPRS_OPT,KIND(1.D0))+eps30))
          write (iunout,*) ' facp ',facp
          DO ISTRA=1,NSTRA
            NPESTR(ISTRA)=NPESTR(ISTRA)+int(TIMPE(ISTRA)*FACP)
            NPRS_FREE=NPRS_FREE-int(TIMPE(ISTRA)*FACP)
          ENDDO
          WRITE (iunout,*) ' NPESTR ',(NPESTR(ISTRA),ISTRA=1,NSTRA)
          WRITE (iunout,*) ' NPRS_FREE ',NPRS_FREE


        else

csw attempting better work load balancing           
          tmean=xtim(0)/dble(nprs)
          do istra=1,nstra
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

          do istra=1,nstra
            write(iunout,'(a,2i6,2(1x,e13.6))') 
     .              'XMCT ',istra,npestr(istra),xmct(istra),xmcp(istra)
          enddo
        endif
 
csw 14jul2011
        do while (nprs_free < 0) 
          WRITE (iunout,*) ' NPRS_FREE ',NPRS_FREE
          i=maxloc(npestr,dim=1)
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
          IF (ISTRA.GT.NSTRA) ISTRA=1
          IF (TIMPE(ISTRA).GT.1.E-10) THEN
            NPESTR(ISTRA)=NPESTR(ISTRA)+1
            NPRS_FREE=NPRS_FREE-1
          ENDIF
        ENDDO
        WRITE (iunout,*) ' NPESTR '
        WRITE (iunout,'(12I6)') (NPESTR(ISTRA),ISTRA=1,NSTRA)
        WRITE (iunout,*) ' NPRS_FREE ',NPRS_FREE

csw 14jul2011
        if(sum(npestr) /= nprs ) then
          write(iunout,*) 'pedist: wrong number of processors in npestr'
          call eirene_exit_own(1)
        endif
csw
 
! assign each processor the numbers ISTRA of the strata it shall work on
        IPE=0
        DO ISTRA=1,NSTRA
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
        NPESTA(1)=0
        DO ISTRA=2,NSTRA
          NPESTA(ISTRA)=NPESTA(ISTRA-1)+NPESTR(ISTRA-1)
        ENDDO
        WRITE (iunout,*) ' MASTER PROCESSOR FOR STRATUM '
        WRITE (iunout,*) ' ISTRA, NPESTA '
        WRITE (iunout,'(12I6)') (I,NPESTA(I),I=1,NSTRA)
 
        XTIM(1:NSTRA) = XX1
        CALL EIRENE_MASAGE
     .    ('REDEFINED CPU TIME ASSIGNED TO STRATA (SEC) :')
        DO ISTRA=1,NSTRA
          CALL EIRENE_MASJ1R ('STRATUM, TIME   ',ISTRA,XTIM(ISTRA))
        END DO

C Rescaling of particles per stratum, to keep total particle number 
C independent of parallelisation (strong scaling appraoch):
        NPTS = NPTS / NPESTR

      END IF  
 
      RETURN
      END SUBROUTINE EIRENE_PEDIST_PROPALLOC

      END SUBROUTINE EIRENE_PEDIST
