C
C
C
C
      SUBROUTINE EIRENE_EIRSRT(LSTOP_in,LTIME_in,DELTAT_in,FLUXES_in,
     .                  B2BRM,B2RD,B2Q,B2VP,STEP_CPU)

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_BRASPOI
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CINIT
      USE EIRMOD_CPOLYG
      USE EIRMOD_CGRID
      USE EIRMOD_CSPEZ
      USE EIRMOD_CZT1
      USE EIRMOD_CTRCEI
      USE EIRMOD_CCOUPL
      USE EIRMOD_CGEOM
      USE EIRMOD_CSDVI
      USE EIRMOD_CSDVI_BGK
      USE EIRMOD_CSDVI_COP
      USE EIRMOD_COMPRT
      USE EIRMOD_COMNNL
      USE EIRMOD_COMSOU
      USE EIRMOD_COUTAU
      USE EIRMOD_COMXS
      USE EIRMOD_CSPEI
      USE EIRMOD_BRASCL
      use eirmod_eirbra

      IMPLICIT NONE
c
      include 'mpif.h'
C
      REAL(DP), INTENT(IN) :: FLUXES_in(NSTRA)
      REAL(DP), INTENT(IN) :: DELTAT_in, B2BRM, B2RD, B2Q, B2VP,STEP_CPU
      LOGICAL, INTENT(IN) :: LSTOP_in, LTIME_in

      integer :: rank_mpi,ierr_mpi,size_mpi
      REAL(DP) :: FLUXES(NSTRA)
      REAL(DP) :: DELTAT
      LOGICAL :: LSTOP, LTIME
      logical :: ltrigger

      REAL(DP) :: FLUXS(NSTRA)
      REAL(DP) :: EIRENE_FTABEI1, EIRENE_FEELEI1, FLXI, ESIG, 
     .            EIRENE_RESET_SECOND, DUMMY,
     .            EIRENE_SECOND_OWN, DTIMVO
      INTEGER :: IN, IAEI, IRDS, IIDS, ICPV, IMDS, IFIRST, K, JC, NDXY,
     .           J, IRC, NREC10, NREC11, ITNR, IPLSTI, IST_RATE, IST,
     .           IFRSTR, ISTH, ISTNEW, ISTIN, ICOSTP
      REAL(DP), ALLOCATABLE :: OUTAU(:)
      INTEGER, ALLOCATABLE :: IHELP(:)
      LOGICAL :: LSTP, LLST, LPLASM, NLSRON_SAVE(NSTRA)
      LOGICAL, ALLOCATABLE, SAVE :: NLSRON_INIT(:)
C
      TYPE(CELLSIM), POINTER :: CPSIM
      TYPE(CELLMUL), POINTER :: CPMUL
      TYPE(RATE_STORE), POINTER :: RTIS
C
C
      SAVE
      DATA IFIRST/0/, ICOSTP/0/

      call mpi_comm_rank(MPI_COMM_WORLD,rank_mpi,ierr_mpi)
      call mpi_comm_size(MPI_COMM_WORLD,size_mpi,ierr_mpi)
      if(rank_mpi == 0) then
        lstop=lstop_in
        ltime=ltime_in
        deltat=deltat_in
        fluxes(1:nstra) = fluxes_in(1:nstra)
        ltrigger=.true.
        call mpi_bcast(ltrigger,1,MPI_LOGICAL,
     .                 0,MPI_COMM_WORLD,ierr_mpi)
      endif    
      call mpi_bcast(lstop,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr_mpi)
      call mpi_bcast(ltime,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr_mpi)
      call mpi_bcast(deltat,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     .               ierr_mpi)
      call mpi_bcast(fluxes,nstra,MPI_DOUBLE_PRECISION, 
     .               0,MPI_COMM_WORLD,ierr_mpi)

! only needed on pe 0
      B2BREM=B2BRM
      B2RAD=B2RD
      B2QIE=B2Q
      B2VDP=B2VP
C
      IF (LTIME) THEN
        if (size_mpi > 1) then
          stop ' MPI for LTIME=.true. not yet available '
        end if
C
!out        B2BREM=B2BRM
!out        B2RAD=B2RD
!out        B2QIE=B2Q
!out        B2VDP=B2VP
        DUMMY=EIRENE_RESET_SECOND()
        IF(IFIRST.EQ.0) THEN
C
          CALL GRSTRT(35,8)
C
C  READ FORMATTED INPUT FILE IUNIN
C  AND RUN EIRENE FOR ONE TIME-CYCLE: ITIMV=1
C  WITH OR WITHOUT INITIAL DISTRIBUTION ON FILE FT15 (NFILE-J FLAG)
C  AS FINAL STRATUM
C  EXPECT PLASMA DATA ON FORT.31 (NLPLAS=.FALSE.)
C
!pb          CALL EIRENE_EIRENE(DELTAT,.FALSE.,.FALSE.,1,.TRUE.)
csw --> MPI_INIT=.false., NLPLAS=.TRUE.
          CALL EIRENE_EIRENE(DELTAT,.TRUE.,.FALSE.,1,.FALSE.)
C
C  EIRENE RUN DONE. CENSUS ARRAY WRITTEN
C  NOW ITIMV=ITIMV+1, NLPLAS=.TRUE.
C
          IF (.NOT.NLPLAS) THEN
            WRITE (iunout,*) 'INCONSISTENT COUPLING '
            WRITE (iunout,*) 'LTIME=TRUE, BUT NTIME = ', NTIME
            CALL EIRENE_EXIT_OWN(1)
          ENDIF
          DO 3 ISTRA=1,NSTRAI
            FLUXS(ISTRA)=FLUX(ISTRA)
3         CONTINUE
          IFIRST=1
        ELSE
C
C  NOW: NLPLAS=.TRUE., I.E., PLASMA DATA EXPECTED ON BRAEIR
C  NOW: ITIMV=ITIMV+1
C  BUT: COMMON BRAEIR REDONE IN EXTERNAL CODE.
C  REACTIVATE INDEX MAPPING, EVEN WITHOUT READING INPUT BLOCK 14 AGAIN
          NCUTB_SAVE=NCUTB
C
          DTIMVO=DTIMV
          DTIMVN=DELTAT
C
C-----------------------------------------------------------------------
C
C  STRATA 1 TO NTARGI ARE SCALED IN PLASMA CODE  (RECYCLING STRATA)
C
C     RETURN TO PLASMA CODE THE PROFILES PER UNIT SOURCE STRENGTH
C     IE. THE PROFILES ARE SCALED BY 1./FLUX(ISTRA) BEFORE RETURN
C
C  STRATA NTARGI+1 TO NSTRAI-1  ARE SCALED BY EIRENE
C
C     (EG. GAS PUFF, VOLUME RECOMBINATION, ETC.)
C     THEY MAY BE RESCALED BY PLASMA CODE FACTORS: FLUXES(ISTRA)
C     RETURN TO PLASMA CODE THE PROFILES SCALED WITH
C     SOURCE STRENGTH: FLUX(ISTRA) (AMP)
C
C  STRATUM NSTRAI IS RESCALED WITH RATIO OF OLD TO NEW TIMESTEP
C
C     RETURN TO PLASMA CODE THE PROFILES WITH FLUX(ISTRA) (AMP)
C
          DO ISTRA=NTARGI+1,NSTRAI-1
            IF (FLUXES(ISTRA).NE.0.) THEN
              FLUX(ISTRA)=FLUXS(ISTRA)*FLUXES(ISTRA)*ELCHA
            ELSE
              FLUX(ISTRA)=FLUXS(ISTRA)
            ENDIF
          ENDDO
C
          IF (DTIMVN.NE.DTIMVO) THEN
            FLUX(NSTRAI)=FLUX(NSTRAI)*DTIMVO/DTIMVN
C
            WRITE (iunout,*) 'FLUX IS RESCALED BY DTIMV_OLD/DTIMV_NEW '
            CALL EIRENE_MASR1('FLUX    ',FLUX(NSTRAI))
            CALL EIRENE_LEER(1)
          ENDIF
C
C-----------------------------------------------------------------------
C
          DTIMV=DTIMVN
C
C  RUN EIRENE ON TIMESTEP DTIMV
C  THEN CALL INTERFACING ROUTINE AT ENTRY IF3COP (FROM EIRENE MAIN)
C
          IITER=1
          IPRNLI=0
          NLSRON=.TRUE.
          CALL EIRENE_EIRENE_COUPLE (LSTOP,1,.TRUE.)
          IF (LSTOP) THEN
            CALL GREND
          ENDIF
        ENDIF
        CALL EIRENE_LEER(2)
        WRITE(*,*) 'EIRENE USED ',EIRENE_SECOND_OWN(),' CPU SECONDS'
        CALL EIRENE_LEER(2)
C
        RETURN
C
      ELSEIF (.NOT.LTIME) THEN

!swpb for multiprocessor calculation
        DUMMY=EIRENE_RESET_SECOND()
C
        IF (IFIRST.GE.1) GOTO 10000

        if(rank_mpi .eq. 0) then

          CALL GRSTRT(35,8)
C
C  READ FORMATTED INPUT FILE IUNIN
C  AND RUN EIRENE FOR ONE TIME-CYCLE: ITIMV=1
C  WITH OR WITHOUT INITIAL DISTRIBUTION ON FILE FT15 (NFILE-J FLAG)
C  AS FINAL STRATUM
C  EXPECT PLASMA DATA ON FORT.31 (NLPLAS=.FALSE.)
C
! already done at top of routine
!out        B2BREM=B2BRM
!out        B2RAD=B2RD
!out        B2QIE=B2Q
!out        B2VDP=B2VP
          LPLASM=.FALSE.
          LLST=LSTOP
          ITNR=1

        end if

 10     CONTINUE

!pb        CALL EIRENE_EIRENE(DELTAT,LPLASM,LLST,ITNR,.TRUE.)
        CALL EIRENE_EIRENE(DELTAT,LPLASM,LLST,ITNR,.FALSE.)
        
        IF (.NOT. ALLOCATED(NLSRON_INIT)) THEN
          ALLOCATE (NLSRON_INIT(NSTRA))
          NLSRON_INIT = NLSRON
          IF (NFULL == 0) NFULL = HUGE(1)
          WRITE (IUNOUT,*) ' SHORT CYCLE STARTED '
          WRITE (IUNOUT,*) 
     .      ' PERFORM FULL EIRENE RUN EVERY NFULL TIMESTEPS '
          WRITE (IUNOUT,*) 'NFULL = ',NFULL
        END IF
C
C  IN THIS CALL TO EIRENE ALREADY IF3COP IS CALLED FOR EACH STRATUM
C  THOSE WITH NLSRON = TRUE  HAVE BEEN RECOMPUTED BY EIRENE
C  THOSE WITH NLSRON = FALSE HAVE BEEN SHORT-CYCLED
C  AT IFIRST   =0: ALL NLSRON=TRUE
C  AT IFIRST.GE.1: FIRST A SHORT CYCLE TEST IS DONE, AND NLSRON IS FOUND
C
        if(rank_mpi .eq. 0) then

        IF (.NOT.LLST) THEN

        IF (IFIRST.GE.1) NLSRON = NLSRON_SAVE

        NDXY=(NDXA-1)*NR1ST+NDYA

        CALL EIRENE_ALLOC_BRASCL

! find new calculated stratum with smallest number
        DO IST = 1, NSTRAI
          IF (NLSRON(IST)) THEN
            IFRSTR = IST
            EXIT
          END IF
        END DO

! determine index of rate storage which has been used in the last iteration 
        IST_RATE = ITS(IFRSTR)

! reduce counters of rate storages for all new calculated strata
        DO IST = 1, NSTRAI
          IF (NLSRON(IST)) THEN
            ISTIN = ITS(IST)
            ITS_COUNT(ISTIN) = ITS_COUNT(ISTIN) - 1
          END IF
        END DO

! check how often storage IST_RATE is still used
        ISTH = 0
        IF (IST_RATE > 0) ISTH= ITS_COUNT(IST_RATE)
        
        IF (ISTH < 1) THEN
! rate storage can be used again
        ELSE
! rate storage still in use, look for an empty slot
          ISTNEW = MINLOC(ITS_COUNT,DIM=1)
          IF (ITS_COUNT(ISTNEW) > 0) THEN
            WRITE (IUNOUT,*) ' PROBLEM IN EIRSRT '
            WRITE (IUNOUT,*) ' ITS_COUNT > 0 '
            WRITE (IUNOUT,*) ' ITS_COUNT ',ITS_COUNT
            CALL EIRENE_EXIT_OWN(1)
          END IF
          IST_RATE = ISTNEW
        END IF
        
        WHERE (NLSRON) 
          ITS = IST_RATE
        END WHERE

        ITS_COUNT(IST_RATE) = COUNT(NLSRON)
C
        CALL EIRENE_ALLOC_RATE_ARRAY(IST_RATE)
        CALL EIRENE_INIT_BRASCL1(IST_RATE)

        RTIS => RTS(IST_RATE)%RTA
C
C  INITIAL: ATOMS, EI-PROCESSES
C
        DO 21 IATM=1,NATMI
        DO 21 IPLS=1,NPLSI
        DO 21 IAEI=1,NAEII(IATM)
          IRDS=LGAEI(IATM,IAEI)
          IF (PPLDS(IRDS,IPLS).EQ.0.) GOTO 21
          DO 22 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              RTIS%SPLODA(IN,IATM,IPLS)=RTIS%SPLODA(IN,IATM,IPLS)+
     .                        TABDS1(IRDS,IN)*PPLDS(IRDS,IPLS)
            ELSE
              RTIS%SPLODA(IN,IATM,IPLS)=RTIS%SPLODA(IN,IATM,IPLS)+
     .                        EIRENE_FTABEI1(IRDS,IN)*PPLDS(IRDS,IPLS)
            END IF
22        CONTINUE
21      CONTINUE
        DO 23 IPLS=1,NPLSI
          IPLSTI= MPLSTI(IPLS)
          DO 24 IN=1,NDXY
            RTIS%SEIODA(IN,IPLS)=DIIN(IPLS,IN)*
     .                      (1.5*TIIN(IPLSTI,IN)+EDRIFT(IPLS,IN))
24        CONTINUE
23      CONTINUE
C
        DO 25 IATM=1,NATMI
        DO 25 IAEI=1,NAEII(IATM)
          IRDS=LGAEI(IATM,IAEI)
          DO 25 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              RTIS%SEEODA(IN,IATM)=RTIS%SEEODA(IN,IATM)+EELDS1(IRDS,IN)*
     .                                        TABDS1(IRDS,IN)
            ELSE
              RTIS%SEEODA(IN,IATM)=RTIS%SEEODA(IN,IATM)+
     .                                        EIRENE_FEELEI1(IRDS,IN)*
     .                                        EIRENE_FTABEI1(IRDS,IN)
            END IF
25      CONTINUE
C
C  INITIAL: TEST IONS, EI-PROCESSES
C
        DO 26 IION=1,NIONI
        DO 26 IIDS=1,NIDSI(IION)
          IRDS=LGIEI(IION,IIDS)
          DO 26 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              RTIS%SEEODI(IN,IION)=RTIS%SEEODI(IN,IION)+EELDS1(IRDS,IN)*
     .                                        TABDS1(IRDS,IN)
            ELSE
              RTIS%SEEODI(IN,IION)=RTIS%SEEODI(IN,IION)+
     .                                        EIRENE_FEELEI1(IRDS,IN)*
     .                                        EIRENE_FTABEI1(IRDS,IN)
            END IF
26      CONTINUE
C
        DO 27 IION=1,NIONI
        DO 27 IPLS=1,NPLSI
        DO 27 IIDS=1,NIDSI(IION)
          IRDS=LGIEI(IION,IIDS)
          IF (PPLDS(IRDS,IPLS).EQ.0.) GOTO 27
          DO 28 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              RTIS%SPLODI(IN,IION,IPLS)=RTIS%SPLODI(IN,IION,IPLS)+
     .                             TABDS1(IRDS,IN)*PPLDS(IRDS,IPLS)
            ELSE
              RTIS%SPLODI(IN,IION,IPLS)=RTIS%SPLODI(IN,IION,IPLS)+
     .                         EIRENE_FTABEI1(IRDS,IN)*PPLDS(IRDS,IPLS)
            ENDIF
28        CONTINUE
27      CONTINUE
C
        DO 29 IION=1,NIONI
        DO 29 IPLS=1,NPLSI
        DO 29 IIDS=1,NIDSI(IION)
          IRDS=LGIEI(IION,IIDS)
          ESIG=EPLDS(IRDS,2)
          DO 30 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              RTIS%SEIODI(IN,IION)=RTIS%SEIODI(IN,IION)+
     .                        TABDS1(IRDS,IN)*ESIG
            ELSE
              RTIS%SEIODI(IN,IION)=RTIS%SEIODI(IN,IION)+
     .                        EIRENE_FTABEI1(IRDS,IN)*ESIG
            END IF
30        CONTINUE
29      CONTINUE
C
C
C  INITIAL: MOLECULES, EI-PROCESSES
C
        DO 35 IMOL=1,NMOLI
        DO 35 IMDS=1,NMDSI(IMOL)
          IRDS=LGMEI(IMOL,IMDS)
          DO 35 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              RTIS%SEEODM(IN,IMOL)=RTIS%SEEODM(IN,IMOL)+EELDS1(IRDS,IN)*
     .                                        TABDS1(IRDS,IN)
            ELSE
              RTIS%SEEODM(IN,IMOL)=RTIS%SEEODM(IN,IMOL)+
     .                                        EIRENE_FEELEI1(IRDS,IN)*
     .                                        EIRENE_FTABEI1(IRDS,IN)
            ENDIF
35      CONTINUE
C
        DO 47 IMOL=1,NMOLI
        DO 47 IPLS=1,NPLSI
        DO 47 IMDS=1,NMDSI(IMOL)
          IRDS=LGMEI(IMOL,IMDS)
          IF (PPLDS(IRDS,IPLS).EQ.0.) GOTO 47
          DO 48 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              RTIS%SPLODM(IN,IMOL,IPLS)=RTIS%SPLODM(IN,IMOL,IPLS)+
     .                             TABDS1(IRDS,IN)*PPLDS(IRDS,IPLS)
            ELSE
              RTIS%SPLODM(IN,IMOL,IPLS)=RTIS%SPLODM(IN,IMOL,IPLS)+
     .                         EIRENE_FTABEI1(IRDS,IN)*PPLDS(IRDS,IPLS)
            END IF
48        CONTINUE
47      CONTINUE
C
        DO 49 IMOL=1,NMOLI
        DO 49 IPLS=1,NPLSI
        DO 49 IMDS=1,NMDSI(IMOL)
          IRDS=LGMEI(IMOL,IMDS)
          ESIG=EPLDS(IRDS,2)
          DO 50 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              RTIS%SEIODM(IN,IMOL)=RTIS%SEIODM(IN,IMOL)+
     .                        TABDS1(IRDS,IN)*ESIG
            ELSE
              RTIS%SEIODM(IN,IMOL)=RTIS%SEIODM(IN,IMOL)+
     .                        EIRENE_FTABEI1(IRDS,IN)*ESIG
            END IF
50        CONTINUE
49      CONTINUE

        END IF
C
C
! already done at top of routine
!out        B2BREM=B2BRM
!out        B2RAD=B2RD
!out        B2QIE=B2Q
!out        B2VDP=B2VP
C
        end if ! rank_mpi == 0

        IFIRST=IFIRST+1
        ICOSTP = ICOSTP + 1
        write (iunout,*) 'icostp ',icostp

        IF (LSTOP) THEN
          CALL EIRENE_DEALLOC_COMUSR
          CALL EIRENE_DEALLOC_CESTIM
          CALL EIRENE_DEALLOC_BRASCL
          CALL EIRENE_DEALLOC_BRASPOI

          CALL GREND
        END IF

        RETURN
C
C  NOT THE FIRST CALL IN THIS CYCLE: CHECK: SHORT LOOP CORRECTION
C                                           OR FULL EIRENE, FOR EACH
C                                           STRATUM INDIVIDUALLY
10000   CONTINUE
C
        if(rank_mpi .eq. 0) then

        LSTP = LSTOP
        NCUTB_SAVE=NCUTB

        CALL EIRENE_ALLOC_BCKGRND

        CALL EIRENE_INTER1
C
        CALL EIRENE_PLASMA
C
        CALL EIRENE_PLASMA_DERIV(0)
C
        CALL EIRENE_SETAMD(2)
C
C  IN PLASMA_DERIV THE BACKGROUND PLASMA STATE HAS BEEN 
C  WRITTEN TO FORT.13
C  NFILEL HAS BEEN CHANGED TO NFILEL = 3 OR 9
C  ==> PLASMA AND REACTION DATA ARE READ IN SUBR. INPUT
C  NOW SAVE REACTION DATA AS WELL IN ORDER TO HAVE A 
C  CONSISTENT PLASMA STATE ON FORT.13
C  
      IF ((NFILEL >=1) .AND. (NFILEL <=5)) THEN
         NFILEL=3
         CALL EIRENE_WRPLAM(TRCFLE,0)
      ELSE IF (NFILEL > 5) THEN
         NFILEL=9
         CALL EIRENE_WRPLAM_XDR(TRCFLE,0)
      END IF

C
        CALL EIRENE_ALLOC_BRASCL
        CALL EIRENE_INIT_BRASCL2
C
C  NEW: ATOMS, EI PROCESSES
C
        DO 101 IATM=1,NATMI
        DO 101 IPLS=1,NPLSI
          DO 102 IAEI=1,NAEII(IATM)
            IRDS=LGAEI(IATM,IAEI)
            IF (PPLDS(IRDS,IPLS).EQ.0.) GOTO 101
            DO 102 IN=1,NDXY
              IF (NSTORDR >= NRAD) THEN
                SPLNWA(IN,IATM,IPLS)=SPLNWA(IN,IATM,IPLS)+
     .                          TABDS1(IRDS,IN)*PPLDS(IRDS,IPLS)
              ELSE
                SPLNWA(IN,IATM,IPLS)=SPLNWA(IN,IATM,IPLS)+
     .                          EIRENE_FTABEI1(IRDS,IN)*PPLDS(IRDS,IPLS)
              END IF
102       CONTINUE
101     CONTINUE
C
        DO 103 IPLS=1,NPLSI
          IPLSTI= MPLSTI(IPLS)
          DO 104 IN=1,NDXY
            SEINWA(IN,IPLS)=DIIN(IPLS,IN)*
     .                      (1.5*TIIN(IPLSTI,IN)+EDRIFT(IPLS,IN))
104       CONTINUE
103     CONTINUE
C
        DO 105 IATM=1,NATMI
          DO 105 IAEI=1,NAEII(IATM)
            IRDS=LGAEI(IATM,IAEI)
            DO 105 IN=1,NDXY
              IF (NSTORDR >= NRAD) THEN
                SEENWA(IN,IATM)=SEENWA(IN,IATM)+EELDS1(IRDS,IN)*
     .                                          TABDS1(IRDS,IN)
              ELSE
                SEENWA(IN,IATM)=SEENWA(IN,IATM)+EIRENE_FEELEI1(IRDS,IN)*
     .                                          EIRENE_FTABEI1(IRDS,IN)
              END IF
105     CONTINUE
C
C  NEW: TEST IONS, EI PROCESSES
C
        DO 106 IION=1,NIONI
          DO 106 IIDS=1,NIDSI(IION)
            IRDS=LGIEI(IION,IIDS)
            DO 106 IN=1,NDXY
              IF (NSTORDR >= NRAD) THEN
                SEENWI(IN,IION)=SEENWI(IN,IION)+EELDS1(IRDS,IN)*
     .                                          TABDS1(IRDS,IN)
              ELSE
                SEENWI(IN,IION)=SEENWI(IN,IION)+EIRENE_FEELEI1(IRDS,IN)*
     .                                          EIRENE_FTABEI1(IRDS,IN)
              END IF
106     CONTINUE
C
        DO 107 IION=1,NIONI
        DO 107 IPLS=1,NPLSI
        DO 107 IIDS=1,NIDSI(IION)
          IRDS=LGIEI(IION,IIDS)
          IF (PPLDS(IRDS,IPLS).EQ.0.) GOTO 107
          DO 108 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SPLNWI(IN,IION,IPLS)=SPLNWI(IN,IION,IPLS)+
     .                             TABDS1(IRDS,IN)*PPLDS(IRDS,IPLS)
            ELSE
              SPLNWI(IN,IION,IPLS)=SPLNWI(IN,IION,IPLS)+
     .                         EIRENE_FTABEI1(IRDS,IN)*PPLDS(IRDS,IPLS)
            END IF
108       CONTINUE
107     CONTINUE
C
        DO 109 IION=1,NIONI
        DO 109 IPLS=1,NPLSI
        DO 109 IIDS=1,NIDSI(IION)
          IRDS=LGIEI(IION,IIDS)
          ESIG=EPLDS(IRDS,2)
          DO 110 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SEINWI(IN,IION)=SEINWI(IN,IION)+TABDS1(IRDS,IN)*ESIG
            ELSE
              SEINWI(IN,IION)=SEINWI(IN,IION)+
     .                        EIRENE_FTABEI1(IRDS,IN)*ESIG
            END IF
110       CONTINUE
109     CONTINUE
C
C  NEW: MOLECULES, EI PROCESSES
C
        DO 115 IMOL=1,NMOLI
        DO 115 IMDS=1,NMDSI(IMOL)
          IRDS=LGMEI(IMOL,IMDS)
          DO 116 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SEENWM(IN,IMOL)=SEENWM(IN,IMOL)+EELDS1(IRDS,IN)*
     .                                        TABDS1(IRDS,IN)
            ELSE
              SEENWM(IN,IMOL)=SEENWM(IN,IMOL)+EIRENE_FEELEI1(IRDS,IN)*
     .                                        EIRENE_FTABEI1(IRDS,IN)
            END IF
116       CONTINUE
115     CONTINUE
C
        DO 117 IMOL=1,NMOLI
        DO 117 IPLS=1,NPLSI
        DO 117 IMDS=1,NMDSI(IMOL)
          IRDS=LGMEI(IMOL,IMDS)
          IF (PPLDS(IRDS,IPLS).EQ.0.) GOTO 117
          DO 118 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SPLNWM(IN,IMOL,IPLS)=SPLNWM(IN,IMOL,IPLS)+
     .                             TABDS1(IRDS,IN)*PPLDS(IRDS,IPLS)
            ELSE
              SPLNWM(IN,IMOL,IPLS)=SPLNWM(IN,IMOL,IPLS)+
     .                         EIRENE_FTABEI1(IRDS,IN)*PPLDS(IRDS,IPLS)
            END IF
118       CONTINUE
117     CONTINUE
C
        DO 119 IMOL=1,NMOLI
        DO 119 IPLS=1,NPLSI
        DO 119 IMDS=1,NMDSI(IMOL)
          IRDS=LGMEI(IMOL,IMDS)
          ESIG=EPLDS(IRDS,2)
          DO 120 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SEINWM(IN,IMOL)=SEINWM(IN,IMOL)+TABDS1(IRDS,IN)*ESIG
            ELSE
              SEINWM(IN,IMOL)=SEINWM(IN,IMOL)+
     .                        EIRENE_FTABEI1(IRDS,IN)*ESIG
            END IF
120       CONTINUE
119     CONTINUE
C
! already done at top of routine
!out        B2BREM=B2BRM
!out        B2RAD=B2RD
!out        B2QIE=B2Q
!out        B2VDP=B2VP
        CALL EIRENE_INTER3(LSTP,IFIRST,1,NSTRAI,0)

        IF (ICOSTP == NFULL) THEN
          NLSRON = NLSRON_INIT
          ICOSTP = 0
          WRITE (IUNOUT,*) ' FULL CALL TO EIRENE '
        END IF

        NLSRON_SAVE = NLSRON

        endif ! rank_mpi == 0

        call mpi_bcast(nlsron,nstrai,MPI_LOGICAL,
     .                 0,MPI_COMM_WORLD,ierr_mpi)
        call eirene_broadcast_eirbra

        IF (ANY(NLSRON(1:NSTRAI))) THEN
!pb           IFIRST=0
           LPLASM=.TRUE.
           LSTP=LSTOP
           ITNR=ITNR+1
           GOTO 10
        END IF
C
        IFIRST=IFIRST+1
        ICOSTP = ICOSTP + 1
        write (iunout,*) 'icostp ',icostp

        IF (LSTOP) THEN
          CALL EIRENE_DEALLOC_COMUSR
          CALL EIRENE_DEALLOC_CESTIM
          CALL EIRENE_DEALLOC_BRASCL
          CALL EIRENE_DEALLOC_BRASPOI

          CALL GREND
        END IF

        RETURN
C
      ENDIF
      return

      entry eirene_eirsrt_broad
      call mpi_bcast(nlsron,nstrai,MPI_LOGICAL,
     .                 0,MPI_COMM_WORLD,ierr_mpi)
      call eirene_broadcast_eirbra
      return

      END
