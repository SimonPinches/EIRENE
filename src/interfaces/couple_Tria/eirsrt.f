!pb APR   16:   pplds -> pplei
!pb APR   16:   eelds -> eelei
!pb MAY   16:   tabds1 -> tabei1
cdr Nov   16    finalizing notational syncronisation (..DS.. (legacy) --> ..EI..)

C  MAIN INTERFACING ROUTINE FOR COUPLED CFD-PLASMA - EIRENE APPLICATIONS

C  This routine is called from CFD PLASMA CODE and provides the entry point into EIRENE.
C
C   SPECIAL TREATMENT OF FIRST CALL TO EIRENE IN THIS (COUPLED) RUN: 

C      CALL EIRENE(..)     (main-routines)
C
C   LATER CALLS:
C
C      CALL EIRENE_COUPLE  (entry to EIRENE  main-routines, bypassing some  initialization stuff)
c

C
      SUBROUTINE EIRENE_EIRSRT(LSTOP,LTIME,DELTAT,FLUXES,
     .                  B2BRM,B2RD,B2Q,B2VP)

C   INPUT:
C     LSTOP: 
C     LTIME: TIME DEPENDENT MODE. PREPARE TIME DEPENDENT OPTIONS,
C            AND THEN CALL EIRENE
C     DELTAT: TIME STEP  (IRRELEVANT IN CASE LTIME=.FALSE.)
C     
C   ONLY FOR EIRENE ENERGY BALANCE DIAGNOSTICS:
C     B2BRM:  TOTAL BREMSSTAHLUNG LOSS IN PREVIOUS B2 STEP
C     B2RD :  TOTAL (LINE) RADIATION LOSS IN PREVIOUS B2 STEP
C     B2Q  :  VOLUMETRIC ENERGY EXCHANGE (ELECTRONS-IONS) DUE TO COULOMB INTERACTION
C     B2VP :  VOLUMETRIC ENERGY EXCHANGE (ELECTRONS-IONS) DUE TO WORK DONE BY ELECTRIC FIELD

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
      USE EIRMOD_CTRIG

      IMPLICIT NONE
C
      REAL(DP), INTENT(IN) :: FLUXES(*)
      REAL(DP), INTENT(IN) :: DELTAT, B2BRM, B2RD, B2Q, B2VP
      LOGICAL, INTENT(IN) :: LSTOP, LTIME

      REAL(DP), ALLOCATABLE, SAVE :: FLUXS(:)
      REAL(DP) :: EIRENE_FTABEI1, EIRENE_FEELEI1, FLXI, ESIG, 
     .            EIRENE_RESET_SECOND, DUMMY,
     .          EIRENE_SECOND_OWN, DTIMVO
      INTEGER :: IN, IAEI, IMEI, IIEI, IREI, ICPV, IFIRST, K, JC, NDXY,
     .           J, IRC, NREC10, NREC11, ITNR, IPLSTI, IST_RATE, IST,
     .           IFRSTR, ISTH, ISTNEW, ISTIN
      REAL(DP), ALLOCATABLE :: OUTAU(:)
      INTEGER, ALLOCATABLE :: IHELP(:)
      LOGICAL :: LSTP, LLST, LPLASM, NLSRON_SAVE(NSTRA)
C
      TYPE(CELLSIM), POINTER :: CPSIM
      TYPE(CELLMUL), POINTER :: CPMUL
      TYPE(RATE_STORE), POINTER :: RTIS
C
C
      SAVE
      DATA IFIRST/0/
C
      IF (LTIME) THEN
C
        B2BREM=B2BRM
        B2RAD=B2RD
        B2QIE=B2Q
        B2VDP=B2VP
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
          LPLASM=.FALSE.
C
          CALL EIRENE_EIRENE(DELTAT,LPLASM,.FALSE.,1,.TRUE.)
C
C  EIRENE RUN DONE. CENSUS ARRAY WRITTEN
C  NOW ITIMV=ITIMV+1, NLPLAS=.TRUE.
C
          IF (.NOT.NLPLAS) THEN
            WRITE (iunout,*) 'INCONSISTENT COUPLING '
            WRITE (iunout,*) 'LTIME=TRUE, BUT NTIME = ', NTIME
            CALL EIRENE_EXIT_OWN(1)
          ENDIF
          IF (.NOT.ALLOCATED(FLUXS)) ALLOCATE (FLUXS(NSTRA))
          DO 3 ISTRA=1,NSTRAI
            FLUXS(ISTRA)=FLUX(ISTRA)
3         CONTINUE
          IFIRST=1


        ELSE  !(IFIRST.GE.1)
C  THIS IS NOT THE FIRST CALL TO EIRENE
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
C  RESCALE CENSUS ARRAY FLUX, DUE TO DIFFERENT TIME STEPS IN PREVIOUS AND CURRENT EIRENE STEP 
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
C  MAIN ENTRY POINT FROM B2 INTO EIRENE, IN CASE OF TIME-INDEPENDENT RUNS  
C
      ELSEIF (.NOT.LTIME) THEN

!swpb for multiprocessor calculation
        DUMMY=EIRENE_RESET_SECOND()
C
        IF (IFIRST.GE.1) GOTO 10000
C
C  FIRST CALL IN PRESENT RUN. 
C  1) INITIALIZE EIRENE
C  2) CALL EIRENE
C  3) PREPARE ARRAYS FOR SEMI-IMPLICIT "SHORT CYCLE" CORRECTION. 
C             STORE SOME A&M RATES FROM PRESENT STEP, FOR NEXT STEP

        CALL GRSTRT(35,8)
C
C  READ FORMATTED INPUT FILE IUNIN
C  AND RUN EIRENE FOR ONE TIME-CYCLE: ITIMV=1
C  WITH OR WITHOUT INITIAL DISTRIBUTION ON FILE FT15 (NFILE-J FLAG)
C  AS FINAL STRATUM
C  EXPECT PLASMA DATA ON FORT.31 (NLPLAS=.FALSE.)
C
        B2BREM=B2BRM
        B2RAD=B2RD
        B2QIE=B2Q
        B2VDP=B2VP
        LPLASM=.FALSE.
        LLST=LSTOP
        ITNR=1

 10     CONTINUE

        CALL EIRENE_EIRENE(DELTAT,LPLASM,LLST,ITNR,.TRUE.)
C
C  IN THIS CALL TO EIRENE ALREADY IF3COP IS CALLED FOR EACH STRATUM
C  THOSE WITH NLSRON = TRUE  HAVE BEEN RECOMPUTED BY EIRENE
C  THOSE WITH NLSRON = FALSE HAVE BEEN SHORT-CYCLED
C  AT IFIRST   =0: ALL NLSRON=TRUE
C  AT IFIRST.GE.1: FIRST A SHORT CYCLE TEST IS DONE, AND NLSRON IS FOUND
C
        IF (.NOT.LLST) THEN

        IF (IFIRST.GE.1) NLSRON = NLSRON_SAVE

csw 12apr2011       NDXY=(NDXA-1)*NR1ST+NDYA
        NDXY=NTRII
C
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
C.....................................................................................

cdr
cdr  now start to store rates from present cycle, for future short cycle corrections
cdr
cdr  to be done
cdr  all these "short cycle data" should only be computed if short cycle is turned on at all

C
C  CURRENT RUN: ION ENERGY DENSITY: FOR ALL IPLS, BUT TIIN(IPLS) MAY BE THE SAME FOR ALL IPLS 
C                                       
C                                       
          DO IPLS=1,NPLSI
            IPLSTI= MPLSTI(IPLS)
            DO IN=1,NDXY
              RTIS%SEIODA(IN,IPLS)=DIIN(IPLS,IN)*
     .                        (1.5*TIIN(IPLSTI,IN)+EDRIFT(IPLS,IN))
            ENDDO
          ENDDO
C
C  NEXT: ATOMS, EI RATES: SPLODA, SEEODA
C               MISSING:  SEIDOA
cdr  correct energy exchange with bulk ions: e0* eplei(irei,ipls,1)+ eheavy* eplei(irei,ipls,2)
cdr  sum over ipls:                          e0* eplei(irei,0,1)   + eheavy *eplei(irei,0,2)
cdr  e0 is taken as center of mass (COM)energy (as appropriate in ei processes, but not in pi processes)
cdr  and eheavy is the kinetic energy release (KER) in reaction irei
cdr  the present short cycle correction only accounts for the KER (=0 for atoms), not for the COM part
C
C  CURRENT RUN: PARTICLE RATE: ATOMS, EI-PROCESSES, FROM IATM TO IPLS,
C                                     SUM OVER ALL EI PROCESSES
C
        DO 21 IATM=1,NATMI
        DO 21 IPLS=1,NPLSI
        DO 21 IAEI=1,NAEII(IATM)
          IREI=LGAEI(IATM,IAEI)
          IF (PPLEI(IREI,IPLS).EQ.0.) GOTO 21
          DO 22 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              RTIS%SPLODA(IN,IATM,IPLS)=RTIS%SPLODA(IN,IATM,IPLS)+
     .                        TABEI1(IREI,IN)*PPLEI(IREI,IPLS)
            ELSE
              RTIS%SPLODA(IN,IATM,IPLS)=RTIS%SPLODA(IN,IATM,IPLS)+
     .                        EIRENE_FTABEI1(IREI,IN)*PPLEI(IREI,IPLS)
            END IF
22        CONTINUE
21      CONTINUE
C
C
C  CURRENT RUN: ELECTRON COOLING RATE: ATOMS, EI-PROCESSES, FROM IATM,
C                                      SUM OVER ALL EI PROCESSES
C
        DO 25 IATM=1,NATMI
        DO 25 IAEI=1,NAEII(IATM)
          IREI=LGAEI(IATM,IAEI)
          DO 25 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
                  RTIS%SEEODA(IN,IATM)=RTIS%SEEODA(IN,IATM)+
     .                                        EELEI1(IREI,IN)*
     .                                        TABEI1(IREI,IN)
            ELSE
              RTIS%SEEODA(IN,IATM)=RTIS%SEEODA(IN,IATM)+
     .                                        EIRENE_FEELEI1(IREI,IN)*
     .                                        EIRENE_FTABEI1(IREI,IN)
            END IF
25      CONTINUE
C
C
C  NEXT: TEST IONS, EI RATES: SPLODI, SEEODI, SEIODI

C
C  CURRENT RUN: PARTICLE RATE: TEST IONS, EI-PROCESSES, FROM IION TO IPLS,
C                                     SUM OVER ALL EI PROCESSES
C
        DO 27 IION=1,NIONI
        DO 27 IPLS=1,NPLSI
        DO 27 IIEI=1,NIEII(IION)
          IREI=LGIEI(IION,IIEI)
          IF (PPLEI(IREI,IPLS).EQ.0.) GOTO 27
          DO 28 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              RTIS%SPLODI(IN,IION,IPLS)=RTIS%SPLODI(IN,IION,IPLS)+
     .                             TABEI1(IREI,IN)*PPLEI(IREI,IPLS)
            ELSE
              RTIS%SPLODI(IN,IION,IPLS)=RTIS%SPLODI(IN,IION,IPLS)+
     .                         EIRENE_FTABEI1(IREI,IN)*PPLEI(IREI,IPLS)
            ENDIF
28        CONTINUE
27      CONTINUE
C
C  CURRENT RUN: ELECTRON COOLING RATE: TEST IONS, EI-PROCESSES, FROM IION,
C                                      SUM OVER ALL EI PROCESSES
          DO 26 IION=1,NIONI
          DO 26 IIEI=1,NIEII(IION)
            IREI=LGIEI(IION,IIEI)
            DO 26 IN=1,NDXY
              IF (NSTORDR >= NRAD) THEN
                RTIS%SEEODI(IN,IION)=RTIS%SEEODI(IN,IION)+
     .                                        EELEI1(IREI,IN)*
     .                                        TABEI1(IREI,IN)
              ELSE
                RTIS%SEEODI(IN,IION)=RTIS%SEEODI(IN,IION)+
     .                                        EIRENE_FEELEI1(IREI,IN)*
     .                                        EIRENE_FTABEI1(IREI,IN)
              END IF
26        CONTINUE
C
C  CURRENT RUN: ION ENERGY EXCHANGE RATE: TEST IONS, EI-PROCESSES, FROM IION
C                                         SUM OVER ALL EI PROCESSES
C                                         SUM OVER ALL IPLS
C
        DO 29 IION=1,NIONI
        DO 29 IPLS=1,NPLSI
        DO 29 IIEI=1,NIEII(IION)
          IREI=LGIEI(IION,IIEI)
!pb 09022016            ESIG=EPLDS(IREI,2)  this was incorrect, 
cdr                     because it was already summed over ipls
          ESIG=EPLEI(IREI,IPLS,2)  ! only KER -part is corrected in short cycle
          DO 30 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              RTIS%SEIODI(IN,IION)=RTIS%SEIODI(IN,IION)+
     .                        TABEI1(IREI,IN)*ESIG
            ELSE
              RTIS%SEIODI(IN,IION)=RTIS%SEIODI(IN,IION)+
     .                        EIRENE_FTABEI1(IREI,IN)*ESIG
            END IF
30        CONTINUE
29      CONTINUE
C
C
C  NEXT: MOLECULES, EI RATES: SPLODM, SEEODM, SEIODM
C
        DO 47 IMOL=1,NMOLI
        DO 47 IPLS=1,NPLSI
        DO 47 IMEI=1,NMEII(IMOL)
            IREI=LGMEI(IMOL,IMEI)
            IF (PPLEI(IREI,IPLS).EQ.0.) GOTO 47
          DO 48 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              RTIS%SPLODM(IN,IMOL,IPLS)=RTIS%SPLODM(IN,IMOL,IPLS)+
     .                             TABEI1(IREI,IN)*PPLEI(IREI,IPLS)
            ELSE
              RTIS%SPLODM(IN,IMOL,IPLS)=RTIS%SPLODM(IN,IMOL,IPLS)+
     .                         EIRENE_FTABEI1(IREI,IN)*PPLEI(IREI,IPLS)
            END IF
48        CONTINUE
47      CONTINUE
C
C
C  CURRENT RUN: ELECTRON COOLING RATE: MOLECULES, EI-PROCESSES, FROM IMOL,
C                                      SUM OVER ALL EI PROCESSES
C
          DO 35 IMOL=1,NMOLI
          DO 35 IMEI=1,NMEII(IMOL)
          IREI=LGMEI(IMOL,IMEI)
            DO 35 IN=1,NDXY
              IF (NSTORDR >= NRAD) THEN
                  RTIS%SEEODM(IN,IMOL)=RTIS%SEEODM(IN,IMOL)+
     .                                          EELEI1(IREI,IN)*
     .                                          TABEI1(IREI,IN)
              ELSE
                RTIS%SEEODM(IN,IMOL)=RTIS%SEEODM(IN,IMOL)+
     .                                        EIRENE_FEELEI1(IREI,IN)*
     .                                        EIRENE_FTABEI1(IREI,IN)
              ENDIF
35        CONTINUE
C
C  CURRENT RUN: ION ENERGY EXCHANGE RATE: TEST IONS, EI-PROCESSES, FROM IION
C                                         SUM OVER ALL EI PROCESSES
C                                         SUM OVER ALL IPLS
        DO 49 IMOL=1,NMOLI
        DO 49 IPLS=1,NPLSI
        DO 49 IMEI=1,NMEII(IMOL)
            IREI=LGMEI(IMOL,IMEI)
!pb 09022106         ESIG=EPLDS(IREI,2)  this was incorrect,
cdr                     because it was already summed over ipls 
            ESIG=EPLEI(IREI,IPLS,2) ! only KER -part is corrected in short cycle
          DO 50 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              RTIS%SEIODM(IN,IMOL)=RTIS%SEIODM(IN,IMOL)+
     .                        TABEI1(IREI,IN)*ESIG
            ELSE
              RTIS%SEIODM(IN,IMOL)=RTIS%SEIODM(IN,IMOL)+
     .                        EIRENE_FTABEI1(IREI,IN)*ESIG
            END IF
50        CONTINUE
49      CONTINUE

        END IF
C
C
        B2BREM=B2BRM
        B2RAD=B2RD
        B2QIE=B2Q
        B2VDP=B2VP
C
        IFIRST=IFIRST+1

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
10000   CONTINUE  ! IFIRST.GE 1

C  PREPARE ARRAY FOR SEMI-IMPLICIT "SHORT CYCLE" CORRECTIONS FOR NEW STEP.
C
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
cdr
cdr  now start to store rates FOR NEXT cycle, for short cycle corrections
cdr
cdr  to be done
cdr  all these "short cycle data" should only be computed if short cycle is turned on at all
C
C  NEW RUN: ION ENERGY DENSITY: FOR ALL IPLS, BUT TIIN(IPLS) MAY BE THE SAME FOR ALL IPLS 
C
        DO IPLS=1,NPLSI
          IPLSTI= MPLSTI(IPLS)
          DO IN=1,NDXY
            SEINWA(IN,IPLS)=DIIN(IPLS,IN)*
     .                      (1.5*TIIN(IPLSTI,IN)+EDRIFT(IPLS,IN))
          ENDDO
        ENDDO
C
C
C  NEXT: ATOMS, EI RATES: SPLNWA, SEENWA
C               MISSING:  SEINWA
cdr  correct energy exchange with bulk ions: e0* eplei(IREI,ipls,1)+ eplei(IREI,ipls,2)
cdr  sum over ipls:                          e0* eplei(IREI,0,1)   + eplei(IREI,0,2)
cdr  the present short cycle correction only accounts for the KER (=0 for atoms)
C
C  NEXT RUN: PARTICLE RATE: ATOMS, EI-PROCESSES, FROM IATM TO IPLS,
C                                     SUM OVER ALL EI PROCESSES
        DO 101 IATM=1,NATMI
        DO 101 IPLS=1,NPLSI
          DO 102 IAEI=1,NAEII(IATM)
            IREI=LGAEI(IATM,IAEI)
            IF (PPLEI(IREI,IPLS).EQ.0.) GOTO 101
            DO 102 IN=1,NDXY
              IF (NSTORDR >= NRAD) THEN
                SPLNWA(IN,IATM,IPLS)=SPLNWA(IN,IATM,IPLS)+
     .                          TABEI1(IREI,IN)*PPLEI(IREI,IPLS)
              ELSE
                SPLNWA(IN,IATM,IPLS)=SPLNWA(IN,IATM,IPLS)+
     .                          EIRENE_FTABEI1(IREI,IN)*PPLEI(IREI,IPLS)
              END IF
102       CONTINUE
101     CONTINUE

C
        DO 105 IATM=1,NATMI
          DO 105 IAEI=1,NAEII(IATM)
            IREI=LGAEI(IATM,IAEI)
            DO 105 IN=1,NDXY
              IF (NSTORDR >= NRAD) THEN
                SEENWA(IN,IATM)=SEENWA(IN,IATM)+EELEI1(IREI,IN)*
     .                                          TABEI1(IREI,IN)
              ELSE
                SEENWA(IN,IATM)=SEENWA(IN,IATM)+EIRENE_FEELEI1(IREI,IN)*
     .                                          EIRENE_FTABEI1(IREI,IN)
              END IF
105     CONTINUE

cdr  no seinwa, because only KER part is in short cycle correction for EI processes
cdr             and for atoms this is identical == 0.0
C
C  NEW: TEST IONS, EI PROCESSES
C
        DO 107 IION=1,NIONI
        DO 107 IPLS=1,NPLSI
        DO 107 IIEI=1,NIEII(IION)
          IREI=LGIEI(IION,IIEI)
          IF (PPLEI(IREI,IPLS).EQ.0.) GOTO 107
          DO 108 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SPLNWI(IN,IION,IPLS)=SPLNWI(IN,IION,IPLS)+
     .                             TABEI1(IREI,IN)*PPLEI(IREI,IPLS)
            ELSE
              SPLNWI(IN,IION,IPLS)=SPLNWI(IN,IION,IPLS)+
     .                         EIRENE_FTABEI1(IREI,IN)*PPLEI(IREI,IPLS)
            END IF
108       CONTINUE
107     CONTINUE
C
C
        DO 106 IION=1,NIONI
          DO 106 IIEI=1,NIEII(IION)
            IREI=LGIEI(IION,IIEI)
            DO 106 IN=1,NDXY
              IF (NSTORDR >= NRAD) THEN
                SEENWI(IN,IION)=SEENWI(IN,IION)+EELEI1(IREI,IN)*
     .                                          TABEI1(IREI,IN)
              ELSE
                SEENWI(IN,IION)=SEENWI(IN,IION)+EIRENE_FEELEI1(IREI,IN)*
     .                                          EIRENE_FTABEI1(IREI,IN)
              END IF
106     CONTINUE
C

        DO 109 IION=1,NIONI
        DO 109 IPLS=1,NPLSI
        DO 109 IIEI=1,NIEII(IION)
          IREI=LGIEI(IION,IIEI)
!pb 09022016          ESIG=EPLDS(IREI,2)
cdr  correct energy exchange with bulk ions: e0* eplei(IREI,ipls,1)+ eplei(IREI,ipls,2)
cdr  sum over ipls:                          e0* eplei(IREI,0,1)   + eplei	(IREI,0,2)
cdr  the present short cycle correction only accounts for the KER (=0 for atoms)
          ESIG=EPLEI(IREI,IPLS,2)
          DO 110 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SEINWI(IN,IION)=SEINWI(IN,IION)+TABEI1(IREI,IN)*ESIG
            ELSE
              SEINWI(IN,IION)=SEINWI(IN,IION)+
     .                        EIRENE_FTABEI1(IREI,IN)*ESIG
            END IF
110       CONTINUE
109     CONTINUE
C
C  NEW: MOLECULES, EI PROCESSES
C
        DO 117 IMOL=1,NMOLI
        DO 117 IPLS=1,NPLSI
        DO 117 IMEI=1,NMEII(IMOL)
          IREI=LGMEI(IMOL,IMEI)
          IF (PPLEI(IREI,IPLS).EQ.0.) GOTO 117
          DO 118 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SPLNWM(IN,IMOL,IPLS)=SPLNWM(IN,IMOL,IPLS)+
     .                             TABEI1(IREI,IN)*PPLEI(IREI,IPLS)
            ELSE
              SPLNWM(IN,IMOL,IPLS)=SPLNWM(IN,IMOL,IPLS)+
     .                         EIRENE_FTABEI1(IREI,IN)*PPLEI(IREI,IPLS)
            END IF
118       CONTINUE
117     CONTINUE
C
C
        DO 115 IMOL=1,NMOLI
        DO 115 IMEI=1,NMEII(IMOL)
          IREI=LGMEI(IMOL,IMEI)
          DO 116 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SEENWM(IN,IMOL)=SEENWM(IN,IMOL)+EELEI1(IREI,IN)*
     .                                        TABEI1(IREI,IN)
            ELSE
              SEENWM(IN,IMOL)=SEENWM(IN,IMOL)+EIRENE_FEELEI1(IREI,IN)*
     .                                        EIRENE_FTABEI1(IREI,IN)
            END IF
116       CONTINUE
115     CONTINUE

        DO 119 IMOL=1,NMOLI
        DO 119 IPLS=1,NPLSI
        DO 119 IMEI=1,NMEII(IMOL)
          IREI=LGMEI(IMOL,IMEI)
!pb 09022016         ESIG=EPLDS(IREI,2)
          ESIG=EPLEI(IREI,IPLS,2)
          DO 120 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SEINWM(IN,IMOL)=SEINWM(IN,IMOL)+TABEI1(IREI,IN)*ESIG
            ELSE
              SEINWM(IN,IMOL)=SEINWM(IN,IMOL)+
     .                        EIRENE_FTABEI1(IREI,IN)*ESIG
            END IF
120       CONTINUE
119     CONTINUE
C
        B2BREM=B2BRM
        B2RAD=B2RD
        B2QIE=B2Q
        B2VDP=B2VP
        CALL EIRENE_INTER3(LSTP,IFIRST,1,NSTRAI,0)

        NLSRON_SAVE = NLSRON

        IF (ANY(NLSRON(1:NSTRAI))) THEN
!pb           IFIRST=0
           LPLASM=.TRUE.
           LSTP=LSTOP
           ITNR=ITNR+1
           GOTO 10
        END IF
C
        IFIRST=IFIRST+1

        IF (LSTOP) THEN
          CALL EIRENE_DEALLOC_COMUSR
          CALL EIRENE_DEALLOC_CESTIM
          CALL EIRENE_DEALLOC_BRASCL
          CALL EIRENE_DEALLOC_BRASPOI

          CALL GREND
        END IF

        RETURN
C
      ENDIF  !(LTIME)

      END

