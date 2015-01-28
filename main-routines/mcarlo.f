c  nov.16th 2005: npts_save = npts always, not only for nlmovie option
c                 because otherwise in iterative mode a stratum cannot be
c                 re-activated, once it was de-activated in a particlar iteration.
c                 v.kotov
c  19.12.05:  bug: no printout of surface tally std. dev., for sum over strata
c             bug fix: here in macrlo.f: sigmaw = stvw and sgmws=stvws added
c             also needed for this bug fix: clear_sumostra, stat_sumostra

!PB 02.03.06: storing of trajectories
!pb 08.11.06: definition of splitting arrays changed
!             RSPLST(NLEVEL,1:NPARTT) --> RSPLST(1:NPARTT,NLEVEL)
!pb 01.12.06: open and close of fort.10 moved to WRSTRT
!pb 05.12.06: COLLECT_CENSUS introduced to allow for time dependent mode in
!             parallel calculation
!pb 15.12.06: COLSUM replaced by COLLECT_COUTAU
!pb 18.12.06: call to PEDIST only done by processor 0 to avoid trouble because
!             of inaccuracies, call to broad_pedist needed to distribute
!             informations calculated in pedist
!pb 05.02.07: copy ALGV to help array before integrating to avoid array bound violation
!pb 22.05.07: option introduced to set the random number seed after a specified
!             number of particles (used to check parallelization)
!pb 15.11.07: Meaning of NTCPU changed:
!             now NTCPU is the amount of cpu time used for particle tracing
!             times used for initialization and integration of result is not
!             taken into account
!   21.07.09: Meaning of XTIM changed: now XTIM is the time allocated for each stratum
!             no longer the end time
!dr 10.05.10: LOCAT0 might also turn off a stratum. Then: skip this is MCARLO, added after call to LOCAT0
c
      SUBROUTINE EIRENE_MCARLO
C
C  MONTE CARLO CALCULATION
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CTSURF
      USE EIRMOD_CESTIM
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CRAND
      USE EIRMOD_CGRID
      USE EIRMOD_CSPEZ
      USE EIRMOD_CZT1
      USE EIRMOD_CTRCEI
      USE EIRMOD_CGEOM
      USE EIRMOD_CSDVI
      USE EIRMOD_CSDVI_BGK
      USE EIRMOD_CSDVI_COP
      USE EIRMOD_COMPRT
      USE EIRMOD_CPES
      USE EIRMOD_COMNNL
      USE EIRMOD_COMSOU
      USE EIRMOD_COMSPL
      USE EIRMOD_COMSIG
      USE EIRMOD_CLGIN
      USE EIRMOD_COUTAU
      USE EIRMOD_CSPEI
      USE EIRMOD_CUPD
      USE EIRMOD_PHOTON

      IMPLICIT NONE

      INCLUDE 'mpif.h'
C
      CHARACTER(6) :: CIS
      CHARACTER(10) :: CDATE, CTIME

      REAL(DP), ALLOCATABLE :: OUTAU(:)
!      REAL(DP) :: DUMMY(NRTAL)
!      REAL(DP) :: ZVOLIN(NRTAL), ZVOLIW(NRTAL),
!     .          XTIM(0:NSTRA), SCLTAL(N1MX,NTALV), DXTIM(0:NSTRA)
      REAL(DP), ALLOCATABLE, SAVE :: DUMMY(:),
     .                               ZVOLIN(:),ZVOLIW(:),SCLTAL(:,:)
      REAL(DP) :: XTIM(0:NSTRA), DXTIM(0:NSTRA)
      REAL(DP) :: ST, FFF, DELT, XFL1,
     .          XPRNLS, XFACT, OVER_ACC, XPRNLI, STW, STWS,
     .          TIMI, EIRENE_SECOND_OWN, XPT, XX1, XPT1, XFL, SECND, XX,
     .          FLX, VAL, ZW, ZWW, VALUE, ZVOLWT, ZVOLNT, FSIG, ZFLUX,
     .          SECND2, OVER, SECND1, WTT, SECDEL, DUMRAN, timan, timen,
     .          tim1, tim2
      REAL(DP), EXTERNAL :: RANF_EIRENE, RANSET_EIRENE

      INTEGER :: NPTS_SAVE(NSTRA), NINITL_SAVE(NSTRA)
      INTEGER :: ITAL, ISDV, IALS, ISTRAA, ISTRAE, ICELL,
     .           IGFFT, IALV, IDV, I, K, IER, IRC, IBGV, NMX, NINIST,
     .           IPANU, ISEED, ISTR, NPTTOT, NREC11, IB, N2,
     .           IC, IR, IGFF, IADD, INDX, ICLV, IADV, ICPV, ISNV,
     .           INODES, J, ISEE, IPTSI, I1, I2, I3, IA, IT, IMCP,
     .           ISUM, NPX, IS, NEW_ITER, ISPC, IN
csw
!pb 03122013      real(dp) :: timstart,timend,timused
!pb 03122013      real(dp), external :: mpi_wtime
      real(dp) :: timused
      integer :: itimstart, itimend, itimrate
csw
      INTEGER, EXTERNAL :: RANGET_EIRENE
C
      LOGICAL :: LGSTOP, NLPOLS, NLTORS
C  OVERHEAD FOR POST PROCESING (SECONDS)
      DATA N2/2/
C
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
      TIMI=EIRENE_SECOND_OWN()
      timan=timi
      tim1 = timi
C
      IF (NFILEN.NE.0) THEN
        NREC11=NOUTAU
        OPEN (UNIT=11+ifoff,ACCESS='DIRECT',FORM='UNFORMATTED',
     .        RECL=8*NREC11)
      ENDIF

      IF (.NOT.ALLOCATED(DUMMY)) THEN
        ALLOCATE ( DUMMY(NRTAL),
     .             ZVOLIN(NRTAL),
     .             ZVOLIW(NRTAL),
     .             SCLTAL(N1MX,NTALV))
      END IF

C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
C-------------------------------------------------------------------
C
C** INITIALIZE SOME DATA AND SUBROUTINES (ONCE FOR ALL STRATA) *****
C
C  SCLTAL: FLAG FOR SCALING OF VOLUME AVERAGED TALLY
C  SCLTAL =0  1.
C         =1  ZVOLIN(ICELL)
C         =2  ZW
C         =3  ZVOLIW(ICELL)
C         =4  ZWW
C
      SCLTAL=0.D0
      SCLTAL(1,1:8)=1
      SCLTAL(1,9:56)=3
C
C  INITIALIZE RANDOM NUMBER ARRAYS
      INIV1=0
      INIV2=0
      INIV3=0
      INIV4=0
C  DETERMINE MAXIMAL INTEGER (DEPENDING ON MACHINE)
      IF (NLCRR) THEN
        INTMAX=HUGE(1)
      ENDIF
C
C  IRNDVC MUST BE EVEN AND NOT LARGER THEN 64 (COMMON CRAND)
C  IRNDVC IS THE NUMBER OF RANDOM VECTORS PRODUCED IN ONE CALL TO
C  TO RANDOM SAMPLING ROUTINES
      IF (NLCRR) THEN
        IRNDVC=2
      ELSE
        IRNDVC=64
      ENDIF
      IRNDVH=IRNDVC/2

      TIMen=EIRENE_SECOND_OWN()
cdr      write (iunout,*) 'cpu time for init of mcarlo ', timen-tim1
      tim1 = timen
C
C  INITIALIZE SUBR. STATIS
C
      CALL EIRENE_LEER(1)
      CALL EIRENE_STATS0
      TIM2=EIRENE_SECOND_OWN()
cdr      write (iunout,*) 'cpu time for stats0 ', tim2-tim1
      tim1 = tim2
      CALL EIRENE_STATS0_BGK
      TIM2=EIRENE_SECOND_OWN()
cdr      write (iunout,*) 'cpu time for stats0_bgk ', tim2-tim1
      tim1 = tim2
      CALL EIRENE_STATS0_COP
      TIM2=EIRENE_SECOND_OWN()
cdr      write (iunout,*) 'cpu time for stats0_cop ', tim2-tim1
      tim1 = tim2
      CALL EIRENE_STATS0_SPC
      TIM2=EIRENE_SECOND_OWN()
cdr      write (iunout,*) 'cpu time for stats0_spc ', tim2-tim1
      tim1 = tim2
C  INITIALIZE SUBR. REFLEC AND SPUTER
      CALL EIRENE_REFLC0
      TIM2=EIRENE_SECOND_OWN()
cdr      write (iunout,*) 'cpu time for reflec0 ', tim2-tim1
      tim1 = tim2
      IF (NPHOT > 0) THEN
        CALL EIRENE_REFLC0_PHOTON
        CALL EIRENE_LINE_CUTOFF
        TIM2=EIRENE_SECOND_OWN()
cdr        write (iunout,*) 'cpu time for reflc0_photon ', tim2-tim1
        tim1 = tim2
      END IF
      CALL EIRENE_SPUTR0
      TIM2=EIRENE_SECOND_OWN()
cdr      write (iunout,*) 'cpu time for sputr0 ', tim2-tim1
      tim1 = tim2
C  INITIALIZE SUBR. SAMVOL
      CALL EIRENE_SAMVL0
      TIM2=EIRENE_SECOND_OWN()
cdr      write (iunout,*) 'cpu time for samvl0 ', tim2-tim1
      tim1 = tim2
C  INITIALIZE SUBR. SAMSRF
      CALL EIRENE_SAMSF0
      TIM2=EIRENE_SECOND_OWN()
cdr      write (iunout,*) 'cpu time for samsf0 ', tim2-tim1
      tim1 = tim2
C
C
      IESTR=-1
      IF (NFILEN.EQ.2.OR.NFILEN.EQ.7) GOTO 2000
C
C**** CLEAR WORK AREA FOR SUM OVER STRATA ****************************
C
      CALL EIRENE_CLEAR_SUMOSTRA

!pb      LOGATM=.FALSE.
!pb      LOGION=.FALSE.
!pb      LOGMOL=.FALSE.
!pb      LOGPLS=.FALSE.
!pb      LOGPHOT=.FALSE.
C
C
C   MAXIMAL CALCULATION TIME ALLOWED FOR EACH STRATUM,
C   PROPORTIONAL TO NPTS(ISTRA), OR FLUX(ISTRA) (INPUT)
C   OR LINEAR COMBINATION THEREOF
C   THEREFORE NUMBER OF TEST PARTICLES MAY BE LESS THAN NPTS
C   BUT DO AT LEAST 2 PARTICLES, IN CASE NPTS(ISTRA).GE.2
C
!pb      CALL TRMAIN(XX,NTCPU)
      XX = NTCPU
CVKMPI      XTIM(0)=EIRENE_SECOND_OWN()
CVKMPI      SECND=XTIM(0)
      SECND=EIRENE_SECOND_OWN()

      NPTS_SAVE=NPTS
      NINITL_SAVE = NINITL

      timan=secnd
C
C  REMAINING CPU TIME, SUBSTRACT N2 SECONDS FOR PRINTOUT AND PLOTS
!pb   XX1=XX-N2

C  CHANGED:  use XX=NTCPU seconds of cpu-time for calculation of trajectories
      XX1 = XX

      XPT=0.
      XFL=0.
      DO 7 ISTRA=1,NSTRAI
        IF (NPTS(ISTRA).LE.0.AND.FLUX(ISTRA).GT.0.D0) THEN
          FLUX(ISTRA)=0.D0
          NLSRON(ISTRA) = .FALSE.
          WRITE (iunout,*) 'STRATUM ISTRA= ',ISTRA,
     .                     ' TURNED OFF, BECAUSE NPTS=0'
          CALL EIRENE_LEER(1)
        ENDIF
        IF (NPTS(ISTRA).GT.0.AND.FLUX(ISTRA).LE.0.D0) THEN
          NPTS(ISTRA)=0
          NLSRON(ISTRA) = .FALSE.
          WRITE (iunout,*) 'STRATUM ISTRA= ',ISTRA,
     .                     ' TURNED OFF, BECAUSE FLUX=0.0'
          CALL EIRENE_LEER(1)
        ENDIF
        IF (NPTS(ISTRA).LE.0.OR.FLUX(ISTRA).LE.0.D0)
     .     NLSRON(ISTRA) = .FALSE.
        XPT=XPT+FLOAT(NPTS(ISTRA))
        XFL=XFL+FLUX(ISTRA)
7     CONTINUE
      XPT1=0.
      XFL1=0.
      nsteff=0
      DO 8 ISTRA=1,NSTRAI
        if (npts(istra) .gt. 0) then
CVKMPI          XPT1=XPT1+NPTS(ISTRA)
CVKMPI          XFL1=XFL1+FLUX(ISTRA)
CVKMPI          XTIM(ISTRA)=XTIM(0)+XX1*((1.-ALLOC)*XPT1/(XPT+EPS60)+
CVKMPI     +                             (   ALLOC)*XFL1/(XFL+EPS60))
          XPT1=NPTS(ISTRA) !VKMPI
          XFL1=FLUX(ISTRA) !VKMPI
          XTIM(ISTRA)=XX1*((1.-ALLOC)*XPT1/(XPT+EPS60)+
     +                     (   ALLOC)*XFL1/(XFL+EPS60)) !VKMPI
          nsteff=nsteff+1
        else
CVKMPI          xtim(istra)=xtim(istra-1)
          xtim(istra)=0.0 !VKMPI
        end if
8     CONTINUE

C  REDISTRIBUTE XTIM IN CASE THAT SOURCES ARE SWITCHED OFF (SHORT CYCLE)
CVKMPI      DO ISTRA=1,NSTRAI
CVKMPI        DXTIM(ISTRA)=XTIM(ISTRA)-XTIM(ISTRA-1)
CVKMPI        IF (.NOT.NLSRON(ISTRA)) DXTIM(ISTRA)=0._DP
CVKMPI      END DO

CVKMPI      DO ISTRA=1,NSTRAI
CVKMPI        XTIM(ISTRA)=XTIM(ISTRA-1)+DXTIM(ISTRA)
CVKMPI      END DO

      XTIM(0)=SUM(XTIM(1:NSTRA))
C

      TIMen=EIRENE_SECOND_OWN()

      CALL EIRENE_LEER(2)
      CALL EIRENE_MASAGE
     .  ('LOOP OVER STRATA STARTS AT CPU TIME(SEC):    ')
CVKMPI       CALL EIRENE_MASR1 ('STARTTIM',XTIM(0))
      CALL EIRENE_MASR1 ('STARTTIM',SECND) !VKMPI
      CALL EIRENE_MASAGE
     .  ('CPU TIME ASSIGNED TO STRATA (SEC) :          ')
      IF (ALLOC.EQ.0.D0) THEN
        CALL EIRENE_MASAGE
     .  ('PROPORTIONAL NPTS(ISTRA)                     ')
      ELSEIF (ALLOC.EQ.1.) THEN
        CALL EIRENE_MASAGE
     .  ('PROPORTIONAL FLUX(ISTRA)                     ')
      ELSE
        CALL EIRENE_MASAGE
     .  ('WEIGHTED ALLOCATION BETWEEN NPTS AND FLUX    ')
      ENDIF
      DO 9 ISTRA=1,NSTRAI
CVKMPI        DELT=XTIM(ISTRA)-XTIM(ISTRA-1)
CVKMPI        CALL EIRENE_MASJ1R ('STRATUM, TIME   ',ISTRA,DELT)
        CALL EIRENE_MASJ1R ('STRATUM, TIME   ',ISTRA,XTIM(ISTRA)) !VKMPI
9     CONTINUE
      CALL EIRENE_LEER(2)
C
C  ASSIGN NUMBER OF PARTICLES TO BE STORED ON CENSUS, PROPORTIONAL
C  TO CPU TIME ASSIGNED TO EACH STRATUM
C
      IF (NPRNLI.GT.0) THEN
        WRITE(iunout,*)
     .    'MAXIMUM NUMBER OF PARTICLES THAT WILL BE SAVED '
        WRITE(iunout,*) 'FOR SNAPSHOT ESTIMATORS: PROPORTIONAL TO CPU-'
        WRITE(iunout,*) 'TIME ALLOCATED FOR EACH STRATUM'
        DO  ISTRA=1,NSTRAI
CVKMPI          XFACT=(XTIM(ISTRA)-XTIM(ISTRA-1))/XX1
          XFACT=XTIM(ISTRA)/XX1 !VKMPI
          XPRNLS       =NPRNLI*XFACT+0.5
          NPRNLS(ISTRA)=XPRNLS
        ENDDO
10      ISUM=SUM(NPRNLS(1:NSTRAI))
        IF (ISUM.NE.NPRNLI) THEN
C  ROUND OFF ERRORS
          WRITE (iunout,*) 'ISUM,NPRNLI ',ISUM,NPRNLI
          NMX=0
          NPX=-1
          DO ISTRA=1,NSTRAI
            IF (NPRNLS(ISTRA).GT.NPX) THEN
              NMX=ISTRA
              NPX=NPRNLS(ISTRA)
            ENDIF
          ENDDO
          IS=ISIGN(1,ISUM-NPRNLI)
          NPRNLS(NMX)=NPRNLS(NMX)-IS
          GOTO 10
        ENDIF
        DO  ISTRA=1,NSTRAI
          CALL EIRENE_MASJ2 ('STRATUM, NUMBER ',ISTRA,NPRNLS(ISTRA))
        ENDDO
      ENDIF
C
C  ASSIGN PE'S TO STRATA
C
!pb      IF ((NSTEFF > 0) .AND. (NPRS.GT.nsteff)) THEN
      if (my_pe == 0) CALL EIRENE_PEDIST(XTIM,XX1)
      if (nprs > 1) then
!pb021213        call EIRENE_broad_pedist(xtim,npts,nminpts,trcdbgmpi)
        call EIRENE_broad_pedist(xtim)
        if (.not.nlident) then
          do istra=1,nstrai
            ninitl(istra)=ninitl(istra)+my_pe*10000
          enddo
        endif
        NPRNLS=NPRNLI
      ENDIF
C
C**** INITIALIZE COMMONS COUTAU AND CSPEZ
C
csw 19mar2013 moved to here after call to pedist (xmct/xmcp)
      CALL EIRENE_INIT_COUTAU(NLSRON)
      FASCL(0)=1.
      FMSCL(0)=1.
      FISCL(0)=1.
      FPHSCL(0)=1.
C
C
C**** STRATA LOOP ****************************************************
C
      NPANU=0
      OVER_ACC=0.D0
      NEW_ITER=0
      DO 1000 ISTR=1,NSTRAI

        timan=EIRENE_second_own()

        ISTRA=ISTR
        IF (.NOT.NLSRON(ISTRA)) CYCLE
        IF (PROCFORSTRA(ISTRA,MY_PE)) THEN

C  SPECIAL TREATMENT FOR MOVIE OPTION:
          IF (NLMOVIE) THEN
!pb            ISTRA=NSTRAI-ISTR+1
!pb            IF (ISTRA.EQ.NSTRAI-1) THEN
            IF (ISTR.EQ.1) THEN
C  TOTAL NUMBER OF PARTICLES TO BE LAUNCHED FROM ALL NON-CENSUS STRATA
              NPTTOT=NPRNLI-NPANU
C  REDEFINE NPTS ACCORDING TO XTIM(ISTRA)
              CALL EIRENE_LEER(2)
              WRITE (iunout,*)
     .          'REDEFINE NPTS(ISTRA) BECAUSE OF NLMOVIE OPTION'
              ISUM=0
              DO IS=1,NSTRAI-1
CVKMPI                XFACT=(XTIM(IS)-XTIM(IS-1))/XTIM(NSTRAI-1)
!pb                XFACT=XTIM(ISTRA)/(XTIM(0)-XTIM(NSTRA))  !VKMPI
                XFACT=XTIM(IS)/XTIM(0)  !PB
                XPRNLI=NPTTOT*XFACT+0.5
                NPTS(IS)=XPRNLI
                ISUM=ISUM+NPTS(IS)
                WRITE(iunout,*) 'ISTRA, NPTS = ',IS,NPTS(IS)
              ENDDO
            ENDIF
          ENDIF

C  MOVIE OPTION (NLMOVIE):  DONE

          CALL EIRENE_LEER(2)
          IF (NPTS(ISTRA).GT.0) THEN
            WRITE (iunout,*) 'BEGIN TO WORK ON STRATUM NO. ',ISTRA
          ELSEIF (NPTS(ISTRA).LE.0) THEN
            WRITE (iunout,*) 'STRATUM NO. ',ISTRA,' ABANDONED'
          ENDIF
          CALL EIRENE_LEER(2)
          XMCP(ISTRA)=0.
csw 19feb2013
          XMCT(ISTRA)=0.
csw
!pb        if( ((nprs.le.nsteff).and.(mod(ISTRA-1,nprs).eq.my_pe)) .or.
!pb     .      ((nprs.gt.nsteff).and.(nstrpe(my_pe).eq.istra)) ) then
          IPANU=0
C
C  INITIALIZE RANDOM NUMBER GENERATOR FOR STRATUM ISTRA
        IF (NINITL(ISTRA).GT.0) THEN
          NINIST=NINITL(ISTRA)
          dumran=ranset_eirene(ninist)
          iseed=ranget_eirene(isee)
          ISEEDR=ISEED*0.3D0
          INIV1=0
          INIV2=0
          INIV3=0
          INIV4=0
        ELSEIF (NINITL(ISTRA).LT.0) THEN
          CALL DATE_AND_TIME(CDATE,CTIME)
          READ(CTIME(1:6),*) NINITL(ISTRA)
          WRITE (iunout,*) 'NINITL(ISTRA) SET TO ',NINITL(ISTRA)
          NINIST=NINITL(ISTRA)
          dumran=ranset_eirene(ninist)
          iseed=ranget_eirene(isee)
          ISEEDR=ISEED*0.3D0
          INIV1=0
          INIV2=0
          INIV3=0
          INIV4=0
C       ELSEIF (NINITL(ISTRA).EQ.0) THEN
C  DON'T INITIALIZE FOR THIS STRATUM, NOTHING TO BE DONE HERE
        ENDIF
C
        FASCL(ISTRA)=1.
        FMSCL(ISTRA)=1.
        FISCL(ISTRA)=1.
        FPHSCL(ISTRA)=1.

        LOGATM(:,ISTRA)=.FALSE.
        LOGION(:,ISTRA)=.FALSE.
        LOGMOL(:,ISTRA)=.FALSE.
        LOGPLS(:,ISTRA)=.FALSE.
        LOGPHOT(:,ISTRA)=.FALSE.

        timen=EIRENE_second_own()
C
C  CLEAR WORK AREA FOR THIS STRATUM
C
        CALL EIRENE_CLEAR_STRATUM
C
C  ENFORCE TOROIDAL OR POLOIDAL SYMMETRY FOR THIS STRATUM
C
        IF (NLAVRP(ISTRA)) THEN
          NLPOLS=NLPOL
          NLPOL=.FALSE.
        ENDIF
C
        IF (NLAVRT(ISTRA)) THEN
          NLTORS=NLTOR
          NLTOR=.FALSE.
        ENDIF
C
        IPRNLS=0
C
        IF (NPTS(ISTRA).LE.0) GOTO 1000
C
C  INITIALIZE SUBR. LOCATE
C
        CALL EIRENE_LOCAT0
        IF (NPTS(ISTRA).LE.0) GOTO 1000 ! LOCAT0 might also turn off a stratum        
C
C  LOCATE AND FOLLOW MC-PARTICLES
C
        CALL EIRENE_FTCRI(ISTRA,CIS)
        CALL EIRENE_MASBOX
     .  ('LAUNCH PARTICLES FOR STRATUM NUMBER ISTRA='//CIS)
        OVER=EIRENE_SECOND_OWN()-SECND
C  ACCUMULATED OVERHEAD BETWEEN STRATA
        OVER_ACC=OVER_ACC+OVER
CVKMPI        CALL EIRENE_MASR1 ('OVERHEAD',OVER)
CVKMPI        XTIM(ISTRA)=XTIM(ISTRA)+OVER_ACC
        WRITE (iunout,*) 'XTIM(ISTRA)= ',XTIM(ISTRA)

        TIMI=EIRENE_SECOND_OWN()            !VKMPI
        XTIM(ISTRA)=XTIM(ISTRA)+TIMI !VKMPI
C
        LGLAST=.FALSE.
        LGSTOP=.FALSE.
C
C
csw 19feb2013
!pb 03122013        timstart=mpi_wtime()
        call system_clock (itimstart, itimrate)
csw

csw 18oct2012 TEST
csw        DO 100 IPTSI=1,NPTS(ISTRA)
        DO 100 IPTSI=1,NPTS(ISTRA)/max(1,npestr(istra))
csw
C
C  RESET INDEX-ARRAY
          NCLMT = 0
          DO I=1,NCLMTS
            IN=ICLMT(I)
            IMETCL(IN) = 0
          END DO
          LMETSP=.FALSE.
          NCLMTS = 0

          DO I=1,NWLMT
            IN=IWLMT(I)
            IMETWL(IN) = 0
          END DO
          LMETSPW=.FALSE.
          NWLMT = 0
          NWLMTS = 0

          IF (NADSPC > 0) THEN
            DO ISPC=1,NADSPC
              ESTIML(ISPC)%PSPC%IMETSP = 0
            END DO
          END IF
C
          IF (LGLAST.AND.LGSTOP) THEN
            CALL EIRENE_LEER(1)
            WRITE (iunout,*)
     .        'NO FURTHER COMP.TIME AVAIL. FOR THIS STRATUM'
            WRITE (iunout,*)
     .        'M.C. HISTORIES FOLLOWED UNTIL THAT TIME FOR'
            WRITE (iunout,*) 'THIS STRATUM'
cdr            CALL EIRENE_MASJ2 ('ISTRA,IPANU=    ',ISTRA,IPANU)
            call system_clock (itimend, itimrate)
            timused=real(itimend-itimstart,DP)/REAL(itimrate,DP)
            CALL EIRENE_MASJ2R('ISTRA,IPANU,TIMUSED     ',
     .                          ISTRA,IPANU,TIMUSED)
            IF (NPRNLI.GT.0) THEN
              WRITE (iunout,*) 'M.C. HISTORIES THAT SCORED AT CENSUS'
              CALL EIRENE_MASJ1 ('IPRNLS= ',IPRNLS)
            ENDIF
            IF (TRCLST) CALL EIRENE_OUTLST
            GOTO 101
          ELSEIF (LGLAST.AND..NOT.LGSTOP) THEN
            CALL EIRENE_LEER(1)
            WRITE (iunout,*) 'CENSUS ARRAYS FILLED FOR THIS STRATUM'
            WRITE (iunout,*)
     .        'M.C. HISTORIES FOLLOWED UNTIL THAT TIME FOR'
            WRITE (iunout,*) 'THIS STRATUM'
            call system_clock (itimend, itimrate)
            timused=real(itimend-itimstart,DP)/REAL(itimrate,DP)
            CALL EIRENE_MASJ2R('ISTRA,IPANU,TIMUSED     ',
     .                          ISTRA,IPANU,TIMUSED)
cdr            CALL EIRENE_MASJ2 ('ISTRA,IPANU=    ',ISTRA,IPANU)
            WRITE (iunout,*) 'M.C. HISTORIES THAT SCORED AT CENSUS'
            CALL EIRENE_MASJ1 ('IPRNLS= ',IPRNLS)
            IF (TRCLST) CALL EIRENE_OUTLST
            GOTO 101
          ENDIF
          SECND1=EIRENE_SECOND_OWN()
          LGLAST = IPTSI.EQ.NPTS(ISTRA)
          LGLAST = LGLAST.OR.(SECND1.GT.XTIM(ISTRA).AND.
     .                        IPTSI.GE.NMINPTS(ISTRA).AND.
     .                        .NOT.NLMOVIE)
          LGSTOP = LGLAST
C  NEXT MONTE CARLO HISTORY
          IF (NLCRR) THEN
C  INITIALIZE RANDOM NUMBERS FOR EACH PARTICLE, TO GENERATE CORRELATION
C           Call RANSET_eirene(ISEED)
            dumran=ranset_eirene(iseed)
            DUMRAN=RANF_EIRENE( )
            iseed=ranget_eirene(isee)
            ISEED=INTMAX-ISEED
            INIV1=0
            INIV2=0
            INIV3=0
            INIV4=0
          ELSE IF ((MY_PE == 0) .AND. (NPTSDEL(ISTRA).GT.0)) THEN
            IF (MOD(IPTSI-1,NPTSDEL(ISTRA)) == 0) THEN
              NINIST=NINITL(ISTRA)+IPTSI/NPTSDEL(ISTRA)*10000
              dumran=ranset_eirene(ninist)
              iseed=ranget_eirene(isee)
              ISEEDR=ISEED*0.3D0
              INIV1=0
              INIV2=0
              INIV3=0
              INIV4=0
            END IF
          ENDIF
          XMCP(ISTRA)=XMCP(ISTRA)+1.
          NPANU=NPANU+1
          IPANU=IPANU+1
          ITRJ = NCHORI + MOD(IPANU,NTRJ) + 1
          NLEVEL=0
          CALL EIRENE_LOCAT1(IPANU)
C  IS BIRTH PROCESS SURVIVED?
          IF (.NOT.LGPART) GOTO 110
C
102       CONTINUE
C  FOLLOW NEUTRAL PARTICLE
          IF (ITYP.EQ.0.OR.ITYP.EQ.1.OR.ITYP.EQ.2) THEN
            CALL EIRENE_FOLNEUT
C  FOLLOW TEST ION
          ELSEIF (ITYP.EQ.3) THEN
            CALL EIRENE_FOLION
          ENDIF
C  NEXT GENERATION ?
          IF (LGPART) GOTO 102
C
110       CONTINUE

          IF (NLRAY(ISTRA)) THEN
            CALL EIRENE_CLEAR_TRAJECTORY (ITRJ)
          END IF
C  NUMBER OF REMAINING NODES AND NUMBER OF LEVELS AT NEXT NODE
          IF (NLEVEL.GT.0) THEN
104         INODES=NODES(NLEVEL)-1
            NODES(NLEVEL)=INODES
            IF(INODES.LE.0) GO TO 103
C  RESTORE VARIABLES AND START NEW TRACK
            DO 105 J=1,NPARTC
              RPST(J)=RSPLST(J,NLEVEL)
105         CONTINUE
            DO 106 J=1,MPARTC
              IPST(J)=ISPLST(J,NLEVEL)
106         CONTINUE
            ITYP=ISPEZI(ISPZ,-1)
            IPHOT=ISPEZI(ISPZ,0)
            IATM=ISPEZI(ISPZ,1)
            IMOL=ISPEZI(ISPZ,2)
            IION=ISPEZI(ISPZ,3)
            IPLS=ISPEZI(ISPZ,4)
            CALL EIRENE_NCELLN(NCELL,NRCELL,NPCELL,NTCELL,NACELL,NBLOCK,
     .                  NR1ST,NP2ND,NT3RD,NBMLT,NLRAD,NLPOL,NLTOR)
            NBLCKA=NSTRD*(NBLOCK-1)+NACELL
            NLSRFX=MRSURF.GT.0
            NLSRFY=MPSURF.GT.0
            NLSRFZ=MTSURF.GT.0
            NLSRFA=MASURF.GT.0
            IF (NLTRC) CALL EIRENE_CHCTRC(X0,Y0,Z0,0,12)
            IF (NLSTOR) CALL EIRENE_STORE(200)
            IC_NEUT=0
            IC_ION=0
            GOTO 102
C  RETURN TO PREVIOUS LEVEL
103         CONTINUE
            NLEVEL=NLEVEL-1
            IF(NLEVEL.GT.0) GOTO 104
          ENDIF
C  HISTORY HAS ENDED
C
C  IN CASE NLERG: EITHER LOGATM(1,ISTRA) OR LOGMOL(1,ISTRA)
C  ACTIVATE CORRESPONDING STANDARD DEVIATION ESTIMATOR
C
          IF (NLERG.AND.IPTSI.EQ.1) THEN
            IF (NMOLI >= 1) THEN
              IF (LOGMOL(1,ISTRA)) THEN
                IIH(1)=2
                CALL EIRENE_STATS0
              ENDIF
            ENDIF
          ENDIF

!pb 16012013
!pb  update couple tally after finishing trajektory
          call eirene_upfcop

C
C   MEAN SQUARE
          IF (NSIGI.GT.0) CALL EIRENE_STATS1
     .  (NSBOX_TAL,NR1TAL,NP2TAL,
     .                                     NT3TAL,NLIMPS,
     .                                     NLSYMP(ISTRA),NLSYMT(ISTRA))
          IF (NSIGI_BGK.GT.0) CALL EIRENE_STATS1_BGK
     .  (NSBOX_TAL,NR1TAL,NP2TAL,
     .                                     NT3TAL,NLIMPS,
     .                                     NLSYMP(ISTRA),NLSYMT(ISTRA))
          IF (NSIGI_COP.GT.0) CALL EIRENE_STATS1_COP
     .  (NSBOX_TAL,NR1TAL,NP2TAL,
     .                                     NT3TAL,NLIMPS,
     .                                     NLSYMP(ISTRA),NLSYMT(ISTRA))
          IF (NSIGI_SPC.GT.0) CALL EIRENE_STATS1_SPC
     .  (NSBOX_TAL,NR1TAL,NP2TAL,
     .                                     NT3TAL,NLIMPS,
     .                                     NLSYMP(ISTRA),NLSYMT(ISTRA))
C
          IF (TRCTIM) THEN
            SECND2=EIRENE_SECOND_OWN( )
            SECDEL=SECND2-SECND1
            CALL EIRENE_MASJ1R('PART., CPU TIME ',NPANU,SECDEL)
          ENDIF
100     CONTINUE
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) 'ALL REQUESTED TRAJECTORIES COMPLETED'
        WRITE (iunout,*) 'M.C. HISTORIES FOLLOWED UNTIL THAT TIME FOR'
        WRITE (iunout,*) 'THIS STRATUM'
CDR        CALL EIRENE_MASJ2 ('ISTRA,IPANU=    ',ISTRA,IPANU)

csw 19feb2019
!pb 0312 2013        timend=mpi_wtime()
        call system_clock (itimend, itimrate)
        timused=real(itimend-itimstart,DP)/REAL(itimrate,DP)
        CALL EIRENE_MASJ2R('ISTRA,IPANU,TIMUSED     ',
     .                      ISTRA,IPANU,TIMUSED)
CDR     write(iunout,'(a,2i8,e13.6)') 'TIMUSED: ',istra,ipanu,timused

        IF (NPRNLI.GT.0) THEN
          WRITE (iunout,*) 'M.C. HISTORIES THAT SURVIVED TO CENSUS'
          CALL EIRENE_MASJ1 ('IPRNLS= ',IPRNLS)
        ENDIF
        IF (TRCLST) CALL EIRENE_OUTLST
C       GOTO 101
101     CONTINUE
C
C

        XMCT(istra)=timused
csw
        SECND=EIRENE_SECOND_OWN()
C
C**** PARTICLE TRACING FOR THIS STRATUM FINISHED **********************
C
c
c     collect data for one stratum from all pe's performing calculations
c     for this stratum
c
!pb       if ((nprs.gt.nsteff).and.(nstrpe(my_pe).eq.istra))
!csw        if (count(procforstra(istra,0:nprs-1)) > 1)
       if (nprs.gt.nsteff)
     .  call EIRENE_calstr
C
C  UPDATE AND CHECK LOGICALS FOR TALLIES
C
      DO 120  IMOL=1,NMOLI
        LOGMOL(0,ISTRA)=LOGMOL(0,ISTRA).OR.LOGMOL(IMOL,ISTRA)
        LOGMOL(IMOL,0)=LOGMOL(IMOL,0).OR.LOGMOL(IMOL,ISTRA)
120   CONTINUE
      DO 130  IATM=1,NATMI
        LOGATM(IATM,0)=LOGATM(IATM,0).OR.LOGATM(IATM,ISTRA)
        LOGATM(0,ISTRA)=LOGATM(0,ISTRA).OR.LOGATM(IATM,ISTRA)
130   CONTINUE
      DO 133  IION=1,NIONI
        LOGION(IION,0)=LOGION(IION,0).OR.LOGION(IION,ISTRA)
        LOGION(0,ISTRA)=LOGION(0,ISTRA).OR.LOGION(IION,ISTRA)
133   CONTINUE
      DO 135  IPLS=1,NPLSI
        LOGPLS(IPLS,0)=LOGPLS(IPLS,0).OR.LOGPLS(IPLS,ISTRA)
        LOGPLS(0,ISTRA)=LOGPLS(0,ISTRA).OR.LOGPLS(IPLS,ISTRA)
135   CONTINUE
      DO IPHOT=1,NPHOTI
        LOGPHOT(IPHOT,0)=LOGPHOT(IPHOT,0).OR.LOGPHOT(IPHOT,ISTRA)
        LOGPHOT(0,ISTRA)=LOGPHOT(0,ISTRA).OR.LOGPHOT(IPHOT,ISTRA)
      END DO
      LOGMOL(0,0)=LOGMOL(0,0).OR.LOGMOL(0,ISTRA)
      LOGION(0,0)=LOGION(0,0).OR.LOGION(0,ISTRA)
      LOGATM(0,0)=LOGATM(0,0).OR.LOGATM(0,ISTRA)
      LOGPLS(0,0)=LOGPLS(0,0).OR.LOGPLS(0,ISTRA)
      LOGPHOT(0,0)=LOGPHOT(0,0).OR.LOGPHOT(0,ISTRA)
C
C  NUMBER OF LOCATED M.C. HISTORIES FOR THIS STRATUM: XMCP(ISTRA)
C
csw      if ((nsteff.ge.nprs).or.(npesta(istra).eq.my_pe)) then
      if ((nsteff.ge.nprs).or. procforstra(istra,my_pe)) then

      IF(XMCP(ISTRA).LT.1.) GOTO 1111
C
      WTT=0.
      DO IPHOT=1,NPHOTI
        WTOTPH(0,ISTRA)=WTOTPH(0,ISTRA)+WTOTPH(IPHOT,ISTRA)
        WTT=WTT+WTOTPH(IPHOT,ISTRA)*NPRT(IPHOT)
      END DO
      DO 200 IATM=1,NATMI
        WTOTA(0,ISTRA)=WTOTA(0,ISTRA)+WTOTA(IATM,ISTRA)
        WTT=WTT+WTOTA(IATM,ISTRA)*NPRT(NSPH+IATM)
200   CONTINUE
      DO 201 IMOL=1,NMOLI
        WTOTM(0,ISTRA)=WTOTM(0,ISTRA)+WTOTM(IMOL,ISTRA)
        WTT=WTT+WTOTM(IMOL,ISTRA)*NPRT(NSPA+IMOL)
 201  CONTINUE
      DO 202 IION=1,NIONI
        WTOTI(0,ISTRA)=WTOTI(0,ISTRA)+WTOTI(IION,ISTRA)
        WTT=WTT+WTOTI(IION,ISTRA)*NPRT(NSPAM+IION)
 202  CONTINUE
      WTOTE(ISTRA)=0._DP
      DO 203 IPLS=1,NPLSI
        WTOTP(0,ISTRA)=WTOTP(0,ISTRA)+WTOTP(IPLS,ISTRA)
        WTOTE(ISTRA)=WTOTE(ISTRA)+WTOTP(IPLS,ISTRA)*NCHRGP(IPLS)
        WTT=WTT-WTOTP(IPLS,ISTRA)*NPRT(NSPAMI+IPLS)
 203  CONTINUE
      CALL EIRENE_LEER(2)
      WRITE (iunout,*) 'TOTAL WEIGHT OF PRIMARY SOURCE PARTICLES '
      WRITE (iunout,*) 'BULK IONS, ATOMS, MOLECULES, TEST IONS '
      CALL EIRENE_MASR5 ('WTPLS,WTATM,WTMOL,WTION,WTPHOT          ',
     .     WTOTP(0,ISTRA),WTOTA(0,ISTRA),WTOTM(0,ISTRA),WTOTI(0,ISTRA),
     .     WTOTPH(0,ISTRA))
      WRITE (iunout,*) 'TOTAL NUMBER OF MONTE CARLO HISTORIES'
      IMCP=XMCP(ISTRA)
      CALL EIRENE_MASJ1 ('NPART   ',IMCP)
C
C
C  SET SOME SCALING CONSTANTS
C
      CALL EIRENE_SET_SCAL_CONST (ISTRA, WTT,ZWW, ZW, ZVOLNT, ZVOLWT,
     .                     ZVOLIN, ZVOLIW, SCLTAL, N1MX)
C
C   STATISTICS , IF REQUESTED
C
      IF (XMCP(ISTRA).LE.1.) GOTO 219
C
C  FACTORS FOR STANDARD DEVIATION
      ZFLUX=FLXFAC(ISTRA)*XMCP(ISTRA)
      FSIG=SQRT(XMCP(ISTRA)/(XMCP(ISTRA)-1.))
C
      IF (NSIGI.GT.0) THEN
        CALL EIRENE_STATS2(XMCP(ISTRA),FSIG,ZFLUX)
C  CONVERT TO %
        SDVI1=MAX(0._DP,SDVI1-EPS6)*100.D0
        SDVI2=MAX(0._DP,SDVI2-EPS6)*100.D0
      ENDIF
      IF (NSIGI_BGK.GT.0) THEN
        CALL EIRENE_STATS2_BGK(XMCP(ISTRA),FSIG,ZFLUX)
C  CONVERT TO %
        DO 211 IB=1,NBGVI_STAT
          SGMS_BGK(IB)=MAX(0._DP,SGMS_BGK(IB)-EPS6)*100.D0
          DO 212 J=1,NSBOX_TAL
            SIGMA_BGK(IB,J)=MAX(0._DP,SIGMA_BGK(IB,J)-EPS6)*100.D0
212       CONTINUE
211     CONTINUE
      ENDIF
      IF (NSIGI_COP.GT.0) THEN
        CALL EIRENE_STATS2_COP(XMCP(ISTRA),FSIG,ZFLUX)
C  CONVERT TO %
        DO 213 IC=1,NCPVI_STAT
          SGMS_COP(IC)=MAX(0._DP,SGMS_COP(IC)-EPS6)*100.D0
          DO 214 J=1,NSBOX_TAL
            SIGMA_COP(IC,J)=MAX(0._DP,SIGMA_COP(IC,J)-EPS6)*100.D0
214       CONTINUE
213     CONTINUE
      ENDIF
      IF (NSIGI_SPC.GT.0) THEN
        CALL EIRENE_STATS2_SPC(XMCP(ISTRA),FSIG,ZFLUX)
C  CONVERT TO %
        DO ISPC=1,NADSPC
          ESTIML(ISPC)%PSPC%SGMS=MAX(0._DP,ESTIML(ISPC)%PSPC%SGMS-EPS6)*
     .                           100.D0
          DO J=0,ESTIML(ISPC)%PSPC%NSPC+1
            ESTIML(ISPC)%PSPC%SGM(J)=MAX(0._DP,ESTIML(ISPC)%PSPC%SGM(J)-
     .                               EPS6)*100.D0
          END DO
        END DO
      ENDIF
C
219   CONTINUE

      CALL EIRENE_SCAL_VOLAV_TALLIES (ISTRA, ZWW, ZW,
     .                         ZVOLIN, ZVOLIW, SCLTAL, N1MX)
C
C   REPLACE DEFAULT TALLIES BY USER SUPPLIED
C   COLLISION ESTIMATED TALLIES
C   THIS IS DONE BEFORE VOLUME INTEGRATION THUS THERE IS THE RISK TO
C   DESTROY TERMS NEEDED FOR GLOBAL BALANCES
C
      IF (LCOLV) THEN
        DO 285 ICLV=1,NCLVI
          IS=ICLVS(ICLV)
          IT=ICLVT(ICLV)
          IF (IT.LE.0.OR.IT.GE.NTALA) GOTO 285
          IGFFT=NFSTVI(IT)
          IGFF=NFIRST(IT)
          IF (IS.LE.0.OR.IS.GT.IGFFT) GOTO 285
          IADD=NADDV(IT)
          DO 286 J=1,NSBOX_TAL
            INDX=IADD+(J-1)*IGFF+IS
            ESTIMV(IADD+IS,J)=COLV(ICLV,J)
 286      CONTINUE
 285    CONTINUE
      END IF
C
C   REPLACE DEFAULT TALLIES BY USER SUPPLIED
C   TRACKLENGTH ESTIMATED TALLIES
C   THIS IS DONE BEFORE VOLUME INTEGRATION THUS THERE IS THE RISK TO
C   DESTROY TERMS NEEDED FOR GLOBAL BALANCES
C
      IF (LADDV) THEN
        DO 290 IADV=1,NADVI
          IS=IADVS(IADV)
          IT=IADVT(IADV)
          IF (IT.LE.0.OR.IT.GE.NTALA) GOTO 290
          IGFFT=NFSTVI(IT)
          IGFF=NFIRST(IT)
          IF (IS.LE.0.OR.IS.GT.IGFFT) GOTO 290
          IADD=NADDV(IT)
          DO 295 J=1,NSBOX_TAL
            INDX=IADD+(J-1)*IGFF+IS
            ESTIMV(IADD+IS,J)=ADDV(IADV,J)
 295      CONTINUE
 290    CONTINUE
      END IF
C
C
C   INTEGRATE VOLUME AVERAGED PROFILES   450 --- 459
C
      CALL EIRENE_INTEGRATE_TALLIES (ISTRA)
C
C   SYMMETRISE VOLUME AVERAGED TALLIES?
      IF (NLSYMP(ISTRA).OR.NLSYMT(ISTRA)) THEN
!pb        CALL EIRENE_SYMET(ESTIMV,NTALV,NRTAL,NR1TAL,NP2TAL,NT3TAL,
!pb     .             NLSYMP(ISTRA),NLSYMT(ISTRA))
        CALL EIRENE_SYMET(ESTIMV,NVOLTL,NRTAL,NR1TAL,NP2TAL,NT3TAL,
     .             NLSYMP(ISTRA),NLSYMT(ISTRA))
      ENDIF
C
C  WORK WITH VOLUME AVERAGED TALLIES FOR THIS STRATUM FINISHED
C
C  SCALE SURFACE AVERAGED ESTIMATORS AND OTHER FLUXES 600 - 630
C
      CALL EIRENE_SCAL_SURF_TALLIES (ISTRA)
C
C
C   SUM OVER SURFACE INDEX
C   IN THE SURFACE AVERAGED ESTIMATORS
C
C
C  SUM OVER SPECIES INDEX FOR INTEGRATED VOLUME AVERAGED TALLIES
C                         AND INTEGRATED SURFACE AVERAGED TALLIES
C
      CALL EIRENE_SUM_AVERAGE (ISTRA)
C
      CALL EIRENE_SCALE_TALLIES (ISTRA)
C
C   ALGEBRAIC EXPRESSION IN TALLIES 801--900
C
      IF (NALVI.GT.0.OR.NALSI.GT.0) THEN
C
        CALL EIRENE_ALGTAL
C
        IF (LALGV) THEN
          DO 830 IALV=1,NALVI
            DUMMY(1:NSBOX_TAL) = ALGV(IALV,1:NSBOX_TAL)
            CALL EIRENE_INTTAL (DUMMY,VOLTAL,1,1,
     .                   NSBOX_TAL,ALGVI(IALV,ISTRA),
     .                   NR1TAL,NP2TAL,NT3TAL,NBMLT)
            ALGV(IALV,1:NSBOX_TAL) = DUMMY(1:NSBOX_TAL)
830       CONTINUE
        END IF
C
        IF (LALGS) THEN
          DO 832 IALS=1,NALSI
            ALGSI(IALS,ISTRA)=0.
            DO 831 J=1,NLIMPS
              ALGSI(IALS,ISTRA)=ALGSI(IALS,ISTRA)+ALGS(IALS,J)
831         CONTINUE
832       CONTINUE
        END IF
C
      ENDIF
C
C  SCALE STANDARD DEVIATIONS, WHICH ARE NOT GIVEN IN % REL.ERROR
C  1/XMCP IS INCLUDED IN ZVOLIN,ZW,ZWW,... FOR TALLY AVERAGING
C  THEREFORE IT MUST BE MULTIPLIED HERE BECAUSE ONLY FLUX SCALING
C
      IF (XMCP(ISTRA).LE.1.D0) GOTO 950
C
      CALL EIRENE_SCALE_DEVIATION(ISTRA, ZWW, ZW, ZVOLNT, ZVOLWT,
     .                     ZVOLIN, ZVOLIW, SCLTAL, N1MX)
C
950   CONTINUE
C
C  CALL INTERFACE TO OTHER CODES TO RETURN DATA. STRATUM ISTRA
csw 13mar2013 ONLY WHEN RUN IN NON-PARALLEL MODE 
C
      IF (NMODE.GT.0) THEN
        if(nprs == 1) then
          IESTR=ISTRA
          ISTRAA=ISTRA
          ISTRAE=ISTRA
          CALL EIRENE_IF3COP(ISTRAA,ISTRAE,NEW_ITER)
          NEW_ITER=1
        endif
      ENDIF
C
C  WRITE RESULTS FOR THIS STRATUM ON TEMP. FILE
C
cpara  hier muss fuer den fall nprs > nstrai noch was getan werden!!!
cpara  csw 08mar2013: hat sich jetzt erledigt..
      IESTR=ISTRA
      IF (NFILEN.EQ.1) THEN
csw 18jul2011
csw 08mar2013 added check nprs < nstrai
        if(nprs==1.or.(nprs > 1 .and. npesta(istra)==my_pe)
     .            .or.(nprs > 1 .and. nprs < nstrai) ) then
        CALL EIRENE_WRSTRT(ISTRA,NSTRAI,NESTM1,NESTM2,NADSPC,
     .              ESTIMV,ESTIMS,ESTIML,
     .              NSDVI1,SDVI1,NSDVI2,SDVI2,
     .              NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .              NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .              NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .              NSIGI_SPC,TRCFLE)
        endif
      ENDIF
C
C  UPDATE TALLIES FOR  "SUM OVER STRATA"
C
      IF (NSTRAI.EQ.1) GOTO 1111
C
      CALL EIRENE_SUMOSTRA (ISTRA)
C
C
      IF (NSMSTRA == 1) THEN
        do idv=1,nidv
          smestv(idv,1:nrtal) = smestv(idv,1:nrtal) +
     .                          estimv(idv,1:nrtal)
        end do
        SMESTS = SMESTS + ESTIMS
        DO ISPC=1,NADSPC
          SMESTL(ISPC)%PSPC%SPC = SMESTL(ISPC)%PSPC%SPC +
     .                            ESTIML(ISPC)%PSPC%SPC
          SMESTL(ISPC)%PSPC%SPCINT = SMESTL(ISPC)%PSPC%SPCINT +
     .                               ESTIML(ISPC)%PSPC%SPCINT
        END DO
      END IF
C
      DO 1170 ISDV=1,NSIGCI
        DO 1172 ICELL=1,NSBOX_TAL
          STVC(0,ISDV,ICELL)=STVC(0,ISDV,ICELL)+SIGMAC(0,ISDV,ICELL)
          STVC(1,ISDV,ICELL)=STVC(1,ISDV,ICELL)+SIGMAC(1,ISDV,ICELL)**2
          STVC(2,ISDV,ICELL)=STVC(2,ISDV,ICELL)+SIGMAC(2,ISDV,ICELL)**2
1172    CONTINUE
        STVCS(0,ISDV)=STVCS(0,ISDV)+SGMCS(0,ISDV)
        STVCS(1,ISDV)=STVCS(1,ISDV)+SGMCS(1,ISDV)**2
        STVCS(2,ISDV)=STVCS(2,ISDV)+SGMCS(2,ISDV)**2
1170  CONTINUE
C
1111  CONTINUE
        WRITE(iunout,*)
     .  'CPU TIME USED UNTIL END OF STRATUM ISTRA '
        WRITE(iunout,*) 'ISTRA, CPU(S) ',ISTRA,EIRENE_SECOND_OWN()
        CALL EIRENE_LEER(2)
        endif  ! nprs > nsteff ...
      end if  ! nprs < nsteff ... or  nprs > nstef ...
1000  CONTINUE
C
C*** STRATA LOOP FINISHED *******************************************
C

      NPTS=NPTS_SAVE
      NINITL = NINITL_SAVE
C
      IF (NPRS > 1) THEN
        call EIRENE_collect_coutau
        IF (NPRNLI > 0) CALL EIRENE_COLLECT_CENSUS
        CALL EIRENE_COLLECT_USRDATA
      END IF

C
C  CALL INTERFACE TO OTHER CODES TO RETURN DATA. STRATUM ISTRA
C
csw 08mar2013 shifted behind STRATA LOOP, do all strata in one go
csw 13mar2013 do it here iff in parallel mode
      IF (NMODE.GT.0) THEN
        if(nprs > 1) then
          iestr=istra
          istraa=1
          istrae=nstrai
          CALL EIRENE_IF3COP(ISTRAA,ISTRAE,NEW_ITER)
          NEW_ITER=1
        endif
      ENDIF
csw

      IF ((MY_PE .EQ. 0) .AND. (NSTRAI.EQ.1)) THEN
C
C  WRITE RESULTS FOR SUM OVER STRATA ON TEMP. FILE
C  USE THE DATA FOR STRATUM NO. 1 FOR THIS, RATHER THEN DOING
C  A USELESS SUMMATION
C
C  INDICATE: DATA FOR ISTRA=1 ARE ON CESTIM, BUT WRITE AS SUM OVER
C  STRATA
        IESTR=1
        IF (NFILEN.EQ.1.OR.NFILEN.EQ.6) THEN
          CALL EIRENE_WRSTRT(0,NSTRAI,NESTM1,NESTM2,NADSPC,
     .              ESTIMV,ESTIMS,ESTIML,
     .              NSDVI1,SDVI1,NSDVI2,SDVI2,
     .              NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .              NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .              NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .              NSIGI_SPC,TRCFLE)
        ENDIF
        GOTO 2000
      ENDIF
      IF (XMCP(0).LE.1) GOTO 2000

C SEQUENTIAL REGION

      IF(MY_PE .EQ. 0) THEN

C
C    STATISTICS, SUM OVER STRATA
C
      CALL EIRENE_STAT_SUMOSTRA
C
C  PUT SUM OVER STRATA BACK ONTO CESTIM, CSDVI
C
      IF (NSMSTRA == 1) THEN
        ESTIMV(1:NIDV,1:NRTAL) = SMESTV(1:NIDV,1:NRTAL)
        ESTIMS = SMESTS
        DO ISPC=1,NADSPC
          ESTIML(ISPC)%PSPC%SPC = SMESTL(ISPC)%PSPC%SPC
          ESTIML(ISPC)%PSPC%SPCINT = SMESTL(ISPC)%PSPC%SPCINT
          IF (NSIGI_SPC > 0) THEN
            ESTIML(ISPC)%PSPC%SGM = SMESTL(ISPC)%PSPC%SGM
            ESTIML(ISPC)%PSPC%SGMS = SMESTL(ISPC)%PSPC%STVS
          END IF
        END DO


        SIGMA  = STV
        SGMS   = STVS
        SIGMAW = STVW
        SGMWS  = STVWS

        IF (NSIGI_BGK.GT.0) THEN
          DO 1271 IB=1,NBGVI_STAT
            SGMS_BGK(IB)=STVS_BGK(IB)
            DO 1272 J=1,NSBOX_TAL
              SIGMA_BGK(IB,J)=STV_BGK(IB,J)
1272        CONTINUE
1271      CONTINUE
        ENDIF
        IF (NSIGI_COP.GT.0) THEN
          DO 1273 IC=1,NCPVI_STAT
            SGMS_COP(IC)=STVS_COP(IC)
            DO 1274 J=1,NSBOX_TAL
              SIGMA_COP(IC,J)=STV_COP(IC,J)
1274        CONTINUE
1273      CONTINUE
        ENDIF
C
C   ALGEBRAIC EXPRESSION IN TALLIES, SUM OVER STRATA  1271--1279
C
        IF (NALVI.GT.0.OR.NALSI.GT.0) THEN
C
          CALL EIRENE_ALGTAL
C
          DO 1571 IALV=1,NALVI
            DUMMY(1:NSBOX_TAL) = ALGV(IALV,1:NSBOX_TAL)
            CALL EIRENE_INTTAL (DUMMY,VOLTAL,1,1,
     .                   NSBOX_TAL,ALGVI(IALV,0),
     .                   NR1TAL,NP2TAL,NT3TAL,NBMLT)
            ALGV(IALV,1:NSBOX_TAL) = DUMMY(1:NSBOX_TAL)
1571      CONTINUE
C
          DO 1572 IALS=1,NALSI
            ALGSI(IALS,0)=0.
            DO 1573 J=1,NLIMPS
              ALGSI(IALS,0)=ALGSI(IALS,0)+ALGS(IALS,J)
1573        CONTINUE
1572      CONTINUE
C
        ENDIF
C
C  WRITE RESULTS FOR SUM OVER STRATA ON TEMP. FILE
C
        IESTR=0
        IF (NFILEN.EQ.1.OR.NFILEN.EQ.6) THEN
          CALL EIRENE_WRSTRT(0,NSTRAI,NESTM1,NESTM2,NADSPC,
     .                ESTIMV,ESTIMS,ESTIML,
     .                NSDVI1,SDVI1,NSDVI2,SDVI2,
     .                NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .                NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .                NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .                NSIGI_SPC,TRCFLE)
        ENDIF
      ENDIF
C
      ENDIF
C
2000  CONTINUE

      CALL EIRENE_BROAD_IESTR(IESTR)

      IF(MY_PE .EQ. 0) THEN
C
C  SAVE OR RESTORE SOME DATA FOR "EIRENE RECALL OPTION NFILE.NE.0"
C  FROM FILE "FT11"
C  NOTE: RECORD IRC=3 MAY BE USED IN INTERFACING ROUTINE INFCOP
C
      IF (NFILEN.EQ.1.OR.NFILEN.EQ.6) THEN
        IF (TRCFLE) WRITE (iunout,*) 'WRITE DATA FOR RECALL OPTION '
        IRC=1
        WRITE (11+ifoff,REC=IRC) LOGATM,LOGION,LOGMOL,LOGPLS,LOGPHOT
        IF (TRCFLE)   WRITE (iunout,*) 'WRITE 11  IRC= ',IRC
        IRC=2
        ALLOCATE (OUTAU(NOUTAU))
        CALL EIRENE_WRITE_COUTAU (OUTAU, IUNOUT)
        WRITE (11+ifoff,REC=IRC) OUTAU
        DEALLOCATE (OUTAU)
        IF (TRCFLE)   WRITE (iunout,*) 'WRITE 11  IRC= ',IRC
      ELSEIF (NFILEN.EQ.2.OR.NFILEN.EQ.7) THEN
        IF (TRCFLE) WRITE (iunout,*) 'READ DATA FOR RECALL OPTION'
        IRC=1
        READ (11+ifoff,REC=IRC) LOGATM,LOGION,LOGMOL,LOGPLS,LOGPHOT
        IF (TRCFLE)   WRITE (iunout,*) 'READ 11  IRC= ',IRC
        IRC=2
        ALLOCATE (OUTAU(NOUTAU))
        READ (11+ifoff,REC=IRC) OUTAU
        CALL EIRENE_READ_COUTAU (OUTAU, IUNOUT)
        DEALLOCATE (OUTAU)
        IF (TRCFLE)   WRITE (iunout,*) 'READ 11  IRC= ',IRC
      ENDIF

C END SEQUENTIAL REGION
      ENDIF

!pb 30012013
      call eirene_reset_upfcop

      CALL MPI_BARRIER (MPI_COMM_WORLD,IER)

C
      RETURN

      ENTRY MCARLO2

      IF (ALLOCATED(DUMMY)) THEN
         DEALLOCATE (DUMMY,ZVOLIN,ZVOLIW,SCLTAL)
      END IF

      RETURN
      END
