cdr  sept.2015:  NLSCL scaling option for rpartc (weights of census scores),
cdr              but not for census flux FLUX(NSTRAI) ???  to be done ??
cdr              added census fluxes resolved wrt. species and stratum
cdr              currently only for diagnostic printout, but should be used
cdr              also for stratifying re-sampling to preserve species specific
cdr              fluxes exactly
cdr  aug. 2016:  NLSCL corrections are also not on partw, i.e. not
cdr              accounted for during bootstrapping (re-sampling) from census
cdr  jul. 2020:  statement 300 continue moved up a bit.
cdr              This ensures that fort.15 (census array) is written
cdr              even in case of zero flux to census (empty census then). To
cdr              facilitate continuation in time-dep runs even
cdr              if "zeroth time step" (census initialization) was too large.
cdr  oct. 2021:  Remove target-pointer structures RPSTT, IPSTT. Only use
cdr              the really necessary components of state vectors.
cdr  jun. 2022:  include particle balance scaling factors in the
cdr              census fluxes and weights.
cdr  jan. 2023:  separate counting of Scores and Histories on census,
cdr              due to new cascading options

C
      SUBROUTINE EIRENE_MOD_TMSTEP
C
C  THIS SUBROUTINE IS CALLED AFTER EACH TIME CYCLE. IT ALLOWS TO
C  MODIFY SOME PLASMA BACKGROUND AND PRIMARY SOURCE DATA (I.E., STRATA
C  ISTRA=1,NSTRAI-1) ACCORDING TO INPUT SPECIFICATIONS IN BLOCK 13.
C
C  IT THEN DEFINES THE ADDITIONAL "CENSUS-STRATUM" ISTRA=NSTRAI, I.E., THE SOURCE DUE TO
C  THE INITIAL CONDITION AT THE BEGINNING OF THE NEXT TIMESTEP.
C  FLUX(NSTRAI)    : ATOMIC FLUX
C  RPARTW(1:IPNRL) : CUMULATIVE DISTRIBUTION FOR SAMPLING INDEX I ON CENSUS
C                    NOT NORMALIZED, AND WITHOUT NPRT FLUX FACTORS
C
      USE EIRMOD_PRECISION, ONLY: DP
      USE EIRMOD_PARMMOD, ONLY: MPARTT, NATM, NION, NMOL, NPHOT, NPARTT,
     >                          NSTRA
      USE EIRMOD_COMUSR, ONLY: ISPEZI, ITIMV, NATMI, NIONI, NMOLI,
     >                         NPHOTI, NFILEJ, NPRT, NSPA, NSPAM, NSPH,
     >                         NTIME, TEXTS
      USE EIRMOD_CCONA, ONLY: EPS10, EPS60
      USE EIRMOD_CLOGAU, ONLY: NLMOVIE, NLPLAS
      USE EIRMOD_CPLOT, ONLY: PLTSRC
      USE EIRMOD_CTRCEI, ONLY: TRCPLT, TRCGRD, TRCCEN
      USE EIRMOD_COMPRT, ONLY: IATM, IION, IMOL, IPHOT,
     >                         ITYP, IUNOUT, WEIGHT
cdr   state vector: real(dp) and integer parts.
cdr  >                         IPSTT, RPSTT   !dr  targets for full set of particle "coordinates"
cdr   Oct.21:  removed. The required state vector coordinates are now made explicit:
cdr            NPANUS, ISTRAS, ISPZS, WEIGHTS
      USE EIRMOD_COMNNL, ONLY: FLXCEN, DTIMV, TIME0,
     >                         IPRNL, IPRNLI,
     >                         NPTST, RPARTW,
     >                         IPART, RPART,   !dr  census from
                                               !    this present run.
     >                         IPARTC, RPARTC  !dr  initial condition
                                               !    next timestep
      USE EIRMOD_COMSOU, ONLY: FLUX, NLSRON, NMINPTS, NPTS, NSTRAI,
     >                         NSRFSI, SORWGT
      USE EIRMOD_COUTAU, ONLY: FASCL, FISCL, FMSCL, FLXFAC, FPHSCL, XMCP

      IMPLICIT NONE

      REAL(DP) :: SGMTOT(NSTRA),FLX(NSTRA),
     .            ADDA(0:NATM,0:NSTRA),ADDM(0:NMOL,0:NSTRA),
     .            ADDI(0:NION,0:NSTRA),ADDPH(0:NPHOT,0:NSTRA)
      REAL(DP) :: FLXQ, FCT, SGMREL, SGMTQN, ADDS, ADD, ADDP, SGMTQ1,
     .            SUMM
cdr some state vector components from census score no. I
      REAL(DP) :: WEIGHTS
      INTEGER :: NPANUS, ISTRAS, ISPZS
c
      INTEGER :: ISTRAO, NPANUO,   !dr  from first state vector,
                                   ! for initializing stand. deviations
     .           ISTRAI, JATM, JMOL, JION, JPHOT,
     .           I
      EXTERNAL :: EIRENE_LEER, EIRENE_MASBOX, EIRENE_MASJ2,
     .            EIRENE_MASR1, EIRENE_MASR2,
     .            EIRENE_TMSUSR, EIRENE_WRSNAP
C
C  STEP 1
C
C  REDUCE REDUNDANT PRINTOUT IN NEXT TIME STEP
      TRCPLT=.FALSE.
      TRCGRD=.FALSE.
      PLTSRC(NSTRAI)=.FALSE.
      DO 120 ISTRAI=1,NSTRAI-1
        PLTSRC(NSTRAI)=PLTSRC(NSTRAI).OR.PLTSRC(ISTRAI)
  120 CONTINUE
C
C  SPEED UP GEOMETRY. USE INFO FROM PREVIOUS TIME-STEP
C  tbd.
C
C  STEP 2
C
C  MODIFY BACKGROUND DATA AS COMPARED TO PREVIOUS ITERATION
C  FOR TIME DEP. MODE
C
C  STARTING TIME FOR NEXT TIMESTEP
C
      TIME0=TIME0+DTIMV
      NLSRON(NSTRAI) =.TRUE.
C
C  CALL USER-SUPPLIED ROUTINE TMSUSR
C  E.G.: FILL COMMON BRAEIR WITH NEW PLASMA, IF NMODE.NE.0
C  BE CAREFUL: NO INDEX MAPPING IS DONE, UNLESS
C  NCUTL.NE.NCUTB
      CALL EIRENE_TMSUSR(TIME0)
      NLPLAS=.TRUE.
C
C  STEP 3
C
C  SET SOURCE DUE TO INITIAL CONDITION FOR NEXT TIME CYCLE
C  THERE HAVE BEEN IPRNL SCORES ON THE CENSUS IN THIS PRESENT RUN.
C
C  SOURCE STRENGTH OF INITIAL DISTRIBUTION IN NEW TIME CYCLE

C  A: set NPTS for time stratum NSTRAI IN NEXT TIME STEP
      IPRNL=IPRNLI
      IPRNLI=0
      IF (NPTST.EQ.0) THEN
        NPTS(NSTRAI)=IPRNL
        IF (NLMOVIE) THEN
          WRITE (IUNOUT,*) 'NLMOVIE TURNED OFF, BECAUSE NPTST.EQ.0'
          NLMOVIE=.FALSE.
        ENDIF
      ELSEIF (NPTST.GT.0) THEN
        NPTS(NSTRAI)=NPTST
        IF (NLMOVIE) THEN
          WRITE (IUNOUT,*) 'NLMOVIE TURNED OFF, BECAUSE NPTST.GT.0'
          NLMOVIE=.FALSE.
        ENDIF
      ELSEIF (NPTST.LT.0.OR.NLMOVIE) THEN
C  ONE BY ONE RELAUNCH FROM OLD CENSUS
C  OLD CENSUS CONTAINS IPRNL ENTRIES.
cdr oct.21: I am not sure if that option still works properly
        NPTS(NSTRAI)=IPRNL
        NMINPTS(NSTRAI)=IPRNL  !ENFORCE: FULL RE-LOCATION
                               !OF ALL PARTICLES FROM OLD CENSUS
      ENDIF


C  B: set FLUX for time stratum NSTRAI IN NEXT TIME STEP
      FLUX(NSTRAI)=0.
      FLXCEN=0.
      RPARTW(0)=0.0
      DO 130 ISTRAI=1,NSTRAI
        SGMTOT(ISTRAI)=0.0
        FLX(ISTRAI)=0.0
  130 CONTINUE

C  PREPARE VARIANCE OF CENSUS DATA
      SGMREL=0.0
C
cdr  empty census?
      IF (IPRNL.EQ.0) GOTO 300

C  FIRST SCORE ON CENSUS
cdr what if iprnl=0? are these next arrays properly allocated?

cdr     IPART(1,I)  ! particle number of score no. I on census = NPANU
cdr     IPART(8,I)  ! stratum  number of score no. I on census = ISTRA

cdr  needed for initializing the std. deviation estimates. First score to census.
      NPANUO=IPART(1,1)
      ISTRAO=IPART(8,1)

      ADDS=0.

      ADDPH=0.
      ADDA =0.
      ADDM =0.
      ADDI =0.
C
C  SET "ATOMIC" FLUXES ONTO CENSUS ARRAY
C  IF NLSCL, APPLY PART. BALANCE CORRECTION SCALING, FOR THE TOTAL CENSUS FLUX

C  THIS IS ALREADY DONE ON "TIME SURFACE" TALLY ELSEWHERE,
C  AS IT IS ON ANY OTHER SURFACE.
C  BUT IS DONE HERE ADDITIONALLY ON THE PARTICULAR "CENSUS ARRAYS" RPARTW (RESAMPLING), RPART(9,..),
C  AND THE SPECIES TYPE-RESOLVED FLUXES ADDPH, ADDA, ADDM, ADDI

      DO 140 I=1,IPRNL

cdr  Oct. 21: remove reference to full state vector here. We only need
cdr           npanu, ispz, istra and weight

cdr     IPART(1,I)  ! particle number of score no. I on census = NPANU
cdr     IPART(8,I)  ! stratum  number of score no. I on census = ISTRA
cdr     IPART(9,I)  ! species index of score no. I on census   = ISPZ
cdr     RPART(9,I)  ! particle weight at score to census       = WEIGHT

        WEIGHTS = RPART(9,I)
        NPANUS  = IPART(1,I)
        ISPZS   = IPART(9,I)
        ISTRAS  = IPART(8,I)

        ITYP=ISPEZI(ISPZS,-1)

        IF (ITYP.EQ.0) THEN
          IPHOT=ISPEZI(ISPZS,0)
          WEIGHT=WEIGHTS*FPHSCL(IPHOT,ISTRAS)
          ADDP=WEIGHTS*FLXFAC(ISTRAS)
          ADD=ADDP*NPRT(IPHOT)
          ADDPH(IPHOT,ISTRAS)=ADDPH(IPHOT,ISTRAS)+ADDP
        ELSEIF (ITYP.EQ.1) THEN
          IATM=ISPEZI(ISPZS,1)
          WEIGHT=WEIGHTS*FASCL(IATM,ISTRAS)
          ADDP=WEIGHTS*FLXFAC(ISTRAS)
          ADD=ADDP*NPRT(NSPH+IATM)
          ADDA(IATM,ISTRAS)=ADDA(IATM,ISTRAS)+ADDP
        ELSEIF (ITYP.EQ.2) THEN
          IMOL=ISPEZI(ISPZS,2)
          WEIGHT=WEIGHTS*FMSCL(IMOL,ISTRAS)
          ADDP=WEIGHTS*FLXFAC(ISTRAS)
          ADD=ADDP*NPRT(NSPA+IMOL)
          ADDM(IMOL,ISTRAS)=ADDM(IMOL,ISTRAS)+ADDP
        ELSEIF (ITYP.EQ.3) THEN
          IION=ISPEZI(ISPZS,3)
          WEIGHT=WEIGHTS*FISCL(IION,ISTRAS)
          ADDP=WEIGHTS*FLXFAC(ISTRAS)
          ADD=ADDP*NPRT(NSPAM+IION)
          ADDI(IION,ISTRAS)=ADDI(IION,ISTRAS)+ADDP
        ENDIF

C  SET DISCRETE CUMULATIVE CENSUS FLUX DISTRIBUTION FOR RESAMPLING OF SCORE-INDEX "I"
C  DO NOT INCLUDE NPRT FACTORS, BECAUSE THIS WEIGHT-FACTOŔ WILL ALREADY BE CARRIED
C  BY RELAUNCHED PARTICLE
C  ALSO THE BALANCE CORRECTION FACTORS FATM,... ARE NOT INCLUDED HERE
cdr     RPARTW(I)=RPARTW(I-1)+ADDP  ! same as next line

        RPARTW(I)=RPARTW(I-1)+WEIGHTS*FLXFAC(ISTRAS)
C
C  ACCUMULATE CONTRIBUTION FROM TEST FLIGHT NO. NPANUO
C  FOR ESTIMATION OF STATISTICAL VARIANCE OF CENSUS FLUX
C  NPANUO MAY HAVE SCORED AT CENSUS SEVERAL TIMES
        IF (NPANUS.EQ.NPANUO) THEN
          ADDS=ADDS+ADD
        ENDIF

C  NPANU IS A NEW PARTICLE ?
        IF (NPANUS.NE.NPANUO) THEN
C  ADD PREVIOUS CENSUS SCORE NPANUO TO FLUX, SGMTOT,...
C  FLUX IS TOTAL "ATOMIC FLUX" ON CENSUS
          FLUX(NSTRAI)=FLUX(NSTRAI)+ADDS
          FLX(ISTRAO)=FLX(ISTRAO)+ADDS
          SGMTOT(ISTRAO)=SGMTOT(ISTRAO)+ADDS*ADDS
          NPANUO=NPANUS
          ISTRAO=ISTRAS
C   PUT WEIGHT OF CURRENT CENSUS SCORE NPANU ONTO ADDS
          ADDS=ADD
        ENDIF

C   WEIGHT may have been altered above. So: redefine this component of state vector.
        RPART(9,I)=WEIGHT
  140 CONTINUE   ! IPRNL

C  CONTRIBUTION FROM LAST CENSUS SCORE NO. IPRNL

      FLUX(NSTRAI)=FLUX(NSTRAI)+ADDS
      FLX(ISTRAO)=FLX(ISTRAO)+ADDS
      SGMTOT(ISTRAO)=SGMTOT(ISTRAO)+ADDS*ADDS
C
C  SUM OVER SPECIES, CENSUS FLUXES
      DO ISTRAI=1,NSTRAI
        ADDA(0,ISTRAI) = SUM(ADDA(1:NATMI,ISTRAI))
        ADDM(0,ISTRAI) = SUM(ADDM(1:NMOLI,ISTRAI))
        ADDI(0,ISTRAI) = SUM(ADDI(1:NIONI,ISTRAI))
        ADDPH(0,ISTRAI)= SUM(ADDPH(1:NPHOTI,ISTRAI))
      ENDDO
C  SUM OVER STRATA
      DO JATM=1,NATMI
        ADDA(JATM,0) =SUM(ADDA(JATM,1:NSTRAI))
      ENDDO
      DO JMOL=1,NMOLI
        ADDM(JMOL,0) =SUM(ADDM(JMOL,1:NSTRAI))
      ENDDO
      DO JION=1,NIONI
        ADDI(JION,0) =SUM(ADDI(JION,1:NSTRAI))
      ENDDO
      DO JPHOT=1,NPHOTI
        ADDPH(JPHOT,0)=SUM(ADDPH(JPHOT,1:NSTRAI))
      ENDDO
C  SUM OVER SPECIES AND STRATA
      ADDA(0,0) =SUM(ADDA(0,1:NSTRAI))
      ADDM(0,0) =SUM(ADDM(0,1:NSTRAI))
      ADDI(0,0) =SUM(ADDI(0,1:NSTRAI))
      ADDPH(0,0)=SUM(ADDPH(0,1:NSTRAI))
C
C  VARIANCE OF CENSUS FLUX: ACCOUNT FOR SOURCE STRATIFICATION
C
      SGMTQ1=0.
      SGMREL=0.
      DO 150 ISTRAI=1,NSTRAI
        IF (XMCP(ISTRAI).GT.1.) THEN
          FCT=XMCP(ISTRAI)/(XMCP(ISTRAI)-1.)
          FLXQ=FLX(ISTRAI)*FLX(ISTRAI)
          SGMTQ1=SGMTQ1+(SGMTOT(ISTRAI)-FLXQ/XMCP(ISTRAI))*FCT
        ENDIF
  150 CONTINUE
      IF (SGMTQ1.GT.0.D0) THEN
        SGMTQN=SQRT(SGMTQ1)
        SGMREL=SGMTQN/(FLUX(NSTRAI)+EPS60)
        SGMREL=MAX(0._DP,SGMREL-EPS10)*100.
      ENDIF
C
      IF (FLUX(NSTRAI).GT.0.D0) THEN
        FLXCEN=FLUX(NSTRAI)
        NSRFSI(NSTRAI)=1
        SORWGT(1,NSTRAI)=1.D0
      ENDIF

C  SET ARRAYS OF CENSUS PARTICLE COORDINATES
C  FOR RE-SAMPLING (SUBR. LOCATE) IN NEXT TIMESTEP (=initial condition).
C  WEIGHT OF RE-SAMPLED PARTICLE = 1.0, SINCE ORIGINAL CENSUS PARTICLE
C  WEIGHT IS ALREADY ACCOUNTED FOR IN RE-SAMPLING DISTRIBUTION "RPARTW" SET ABOVE:
      RPARTC(1:NPARTT,1:IPRNL)=RPART(1:NPARTT,1:IPRNL)
      IPARTC(1:MPARTT,1:IPRNL)=IPART(1:MPARTT,1:IPRNL)
C
cdr probably better not to be done here. tbd: move to locate.f, or input.f.
      ISTRAS=NSTRAI
      IPARTC(8,1:IPRNL)=ISTRAS

C
C  CALL WRSNAP TO WRITE SNAPSHOT POPULATION ON FT 15
C  TO BE USED AS INITIAL CONDITION IN NEXT RUN (TIMESTEP)
C
  300 CONTINUE
      IF ((NFILEJ.EQ.1.OR.NFILEJ.EQ.3).AND.ITIMV.GE.NTIME) THEN
cdr  Even in case iprnl=0 (no flux on census):
cdr  Open and write fort.15, at least the first line.
        CALL EIRENE_WRSNAP(NSTRAI)
        WRITE (iunout,*) 'CENSUS ARRAY, FLUX AND TOTAL TIMESTEP STORED'
      ENDIF
C
cdr   300 CONTINUE  ! moved up to ensure fort.15 is written.

      CALL EIRENE_LEER(2)
      WRITE (iunout,*) 'TIME CYCLE COMPLETED, NEXT TIME CYCLE PREPARED'
      WRITE (iunout,*) 'NEXT TIME CYCLE RUNS FROM TIM1 TO TIM2:'
      CALL EIRENE_MASR2('TIM1, TIM2      ',TIME0,TIME0+DTIMV)
      CALL EIRENE_MASJ1('IPRNL   ',IPRNL)
      WRITE (IUNOUT,*) '"ATOMIC" FLUX AT CENSUS (AMP):'
      CALL EIRENE_MASR1('FLUX    ',FLXCEN)
      CALL EIRENE_MASR1('+-%     ',SGMREL)

      CALL EIRENE_LEER(2)
C
      IF (TRCCEN) THEN
        CALL EIRENE_MASBOX
     .          ('DETAILED CENSUS FLUXES, PER SPECIES, STRATUM')

        DO 950 ISTRAI=1,NSTRAI
          CALL EIRENE_LEER(1)
          WRITE (IUNOUT,*) 'STRATUM NO. ISTRA=', ISTRAI
          CALL EIRENE_LEER(1)
          SUMM=ADDA(0,ISTRAI)+ADDM(0,ISTRAI)+ADDI(0,ISTRAI)+
     .         ADDPH(0,ISTRAI)
          IF (SUMM.EQ.0.0) THEN
            WRITE (IUNOUT,*) 'NO FLUXES FROM THIS STRATUM ON CENSUS'
            GOTO 950
          ENDIF

C   PRINT ADDA(IATM) FOR ISTRAI
          IF  (ADDA(0,ISTRAI).GT.0.0) THEN
            WRITE (IUNOUT,*) 'ATOM FLUX AT CENSUS (AMP):'
            DO JATM = 1, NATMI
              WRITE (IUNOUT,'(A10,ES12.4)')
     .               TEXTS(NSPH+JATM), ADDA(JATM,ISTRAI)
            END DO
            WRITE (IUNOUT,'(A10,ES12.4)')
     .               'TOTAL     ',ADDA(0,ISTRAI)
            CALL EIRENE_LEER(1)
          ENDIF

C   PRINT ADDM(IMOL) FOR ISTRAI
          IF  (ADDM(0,ISTRAI).GT.0.0) THEN
            WRITE (IUNOUT,*) 'MOLECULE FLUX AT CENSUS (AMP):'
            DO JMOL = 1, NMOLI
              WRITE (IUNOUT,'(A10,ES12.4)')
     .               TEXTS(NSPA+JMOL), ADDM(JMOL,ISTRAI)
            END DO
            WRITE (IUNOUT,'(A10,ES12.4)')
     .               'TOTAL     ',ADDM(0,ISTRAI)
            CALL EIRENE_LEER(1)
          ENDIF

C   PRINT ADDI(IION) FOR ISTRAI
          IF  (ADDI(0,ISTRAI).GT.0.0) THEN
            WRITE (IUNOUT,*) 'TEST ION FLUX AT CENSUS (AMP):'
            DO JION = 1, NIONI
              WRITE (IUNOUT,'(A10,ES12.4)')
     .               TEXTS(NSPAM+JION), ADDI(JION,ISTRAI)
            END DO
            WRITE (IUNOUT,'(A10,ES12.4)')
     .               'TOTAL     ',ADDI(0,ISTRAI)
            CALL EIRENE_LEER(1)
          ENDIF

C   PRINT ADDPH(IPHOT) FOR ISTRAI
          IF  (ADDPH(0,ISTRAI).GT.0.0) THEN
            WRITE (IUNOUT,*) 'PHOTON FLUX AT CENSUS (AMP):'
            DO JPHOT = 1, NPHOTI
              WRITE (IUNOUT,'(A10,ES12.4)')
     .               TEXTS(JPHOT), ADDPH(JPHOT,ISTRAI)
            END DO
            WRITE (IUNOUT,'(A10,ES12.4)')
     .             'TOTAL     ',ADDPH(0,ISTRAI)
            CALL EIRENE_LEER(1)
          ENDIF

  950   CONTINUE
        CALL EIRENE_LEER(2)

        WRITE (IUNOUT,*) 'SUM OVER STRATA'

        SUMM=ADDA(0,0)+ADDM(0,0)+ADDI(0,0)+ADDPH(0,0)
        IF (SUMM.EQ.0.0) THEN
          WRITE (IUNOUT,*) 'NO FLUXES ON CENSUS '
          GOTO 960
        ENDIF

C   PRINT ADDA(IATM) FOR ISTRAI=0
        IF  (ADDA(0,0).GT.0.0) THEN
          WRITE (IUNOUT,*) 'ATOM FLUX AT CENSUS (AMP):'
          DO JATM = 1, NATMI
            WRITE (IUNOUT,'(A10,ES12.4)')
     .             TEXTS(NSPH+JATM), ADDA(JATM,0)
          END DO
          WRITE (IUNOUT,'(A10,ES12.4)')
     .             'TOTAL     ',ADDA(0,0)
          CALL EIRENE_LEER(1)
        ENDIF

C   PRINT ADDM(IMOL) FOR ISTRAI=0
        IF  (ADDM(0,0).GT.0.0) THEN
          WRITE (IUNOUT,*) 'MOLECULE FLUX AT CENSUS (AMP):'
          DO JMOL = 1, NMOLI
            WRITE (IUNOUT,'(A10,ES12.4)')
     .             TEXTS(NSPA+JMOL), ADDM(JMOL,0)
          END DO
          WRITE (IUNOUT,'(A10,ES12.4)')
     .             'TOTAL     ',ADDM(0,0)
          CALL EIRENE_LEER(1)
        ENDIF

C   PRINT ADDI(IION) FOR ISTRAI=0
        IF  (ADDI(0,0).GT.0.0) THEN
          WRITE (IUNOUT,*) 'TEST ION FLUX AT CENSUS (AMP):'
          DO JION = 1, NIONI
            WRITE (IUNOUT,'(A10,ES12.4)')
     .             TEXTS(NSPAM+JION), ADDI(JION,0)
          END DO
          WRITE (IUNOUT,'(A10,ES12.4)')
     .             'TOTAL     ',ADDI(0,0)
          CALL EIRENE_LEER(1)
        ENDIF

C   PRINT ADDPH(IPHOT) FOR ISTRAI=0
        IF  (ADDPH(0,0).GT.0.0) THEN
          WRITE (IUNOUT,*) 'PHOTON FLUX AT CENSUS (AMP):'
          DO JPHOT = 1, NPHOTI
            WRITE (IUNOUT,'(A10,ES12.4)')
     .             TEXTS(JPHOT), ADDPH(JPHOT,0)
          END DO
          WRITE (IUNOUT,'(A10,ES12.4)')
     .             'TOTAL     ',ADDPH(0,0)
          CALL EIRENE_LEER(1)
        ENDIF

  960   CONTINUE

      ENDIF  !TRCCEN

      WRITE (IUNOUT,*) '...............................................'
C
      RETURN
      END SUBROUTINE EIRENE_MOD_TMSTEP
