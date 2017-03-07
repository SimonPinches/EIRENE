C 27.6.05 : colphot: iadd removed
C 27.6.05 : colatm : mode, il removed (lgaot(..2),lgaot(..5)
C 27.6.05 : colphot: mode, il removed (lgphot(..2),lgphot(..5)
C 15.12.05: irds --> irei, iids --> iiei
C 16.12.05: wminv re-connected to "ei"-processes. suppress
C           reactions with zero test particle secondaries
C           now connected for: colatm, colmol, colion
C           still to be done: include other processes, and colphot
C 2.2.06:  wghtO set at suppression of absorption, for collision estimators.
C 2.2.06:  REMOVED: OT PROCESSES FOR ATOMS
C          GENERATION LIMIT FOR POST COLLISION ATOMS FROM PHOTONS: REMOVED
C 10.3.06: bug fix: LGEI_RED(NREI) --> LGEI_RED(0:NREI)
C          (some compilers had been unhappy with this)
C 20.3.07: PI reactions revised

cdr oct 14.14 some hard wired additional tallies ADDV removed again
cdr oct.21.14 evaluate v-parallel of incident particle only in case of need
c             i.e.  momentum collision estimators, or generation limit
c             otherwise: avoid calls to bfield.f
c
cdr  5. 8.15: ARGUMENTS ADDED TO VECUSR
cdr 20.10.15: arguments in chctrc: type of collision process: corrected for PI and OT
cdr 24.11.15:  bug fix re coll est for pi processes, in colion: eiml --> eiio
cdr Dec.15  :  bug fix pi reaction and cascading was wrong: 
cdr            irei, rather than irpi, and p2nd 
cdr            rather than p2np, were used also for PI reactions. now corrected

cdr         :  further: collision estimators for PI processes, e§pl and e§el tallies: activated
cdr         :  see also corresponding corrections/changes in update for tracklength estimators
cdr DEC. 15 :  bulk ion energy estimators: species resolved.
cdr            not ready: esigei(4, ...), esigpi(4,...) must be species resolved.

cdr            tbd:  check setting of iestm..flags for collision estimators. 
cdr                  probably not correct (outdated).


!pb  APR  16:  ipplds -> ipplei, pplds -> pplei
!pb  APR  16:  patds -> patei
!pb  APR  16:  pmlds -> pmlei
!pb  APR  16:  piods -> pioei
!pb  MAY  16:  nrds  -> nrei
cdr  sept 16:  nmdsi -> nmeii, nidsi -> nieii


cdr Aug 16:    bug fix: IPPLEI --> IPPLPI at one instance
cdr Nov 16:
cdr analog cascading NLCASCAD: started to document, 
cdr        syncronize and re-activate option, not ready !!
c   this version: prepare cascading at collisions, 
c   e.g. for antithetic variate sampling to reduce stochastic cancellation
c   start to clean up splitting, for analogue game and for anticorrelated momentum estimators
c   started for colatm, and ei processes.
c   not sure if ispz is known, NOW
cdr tbd: 
c   cascading with EI: nlevel =nlevel+ptot-1 (because one particle continues)
c   cascading with CX: define analogue PTOT
c   cascading with PI: identical to EI ?? 

cdr Nov. 16:   cflag(7,3) --> cflag(7,mstor0) 
cdr            (was already corrected much earlier in SOLPS_4.3 by VK,
cdr             then correction somehow lost in more recent EIRENE branches)
cdr Jan. 17:    started to separate more clearly the (unfinished) NLCASCAD option from active code
C               Done for COLATM and EI processes.
C               wminv activated in colmol for ei processes (analog to colatm)




      SUBROUTINE EIRENE_COLLIDE
C
C  SAMPLE FROM COLLISION KERNEL C
C
C  INPUT:  COMPRT, COMMON BLOCK, CONTAINING ACTUAL PARTICLE PARAMETERS
C          CFLAG,  FLAG FOR POST COLLISION KINETICS
C  OUTPUT: COMPRT, MODIFIED TO POST COLLISION PARTICLE PARAMETERS
C          COLTYP, FLAG: =1 CONTINUE IN CALLING ROUTINE
C                           (FOLNEUT OR FOLION)
C                        =2 EXIT FROM CALLING ROUTINE
C                           EITHER ABSORBTION, OR
C                           TRANSITION NEUTRAL-->ION (IF CALLED
C                           BY FOLNEUT), OR
C                           TRANSITION ION-->NEUTRAL (IF CALLED
C                           BY FOLION)
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_CADGEO
      USE EIRMOD_CCONA
      USE EIRMOD_CRAND
      USE EIRMOD_CINIT
      USE EIRMOD_CZT1
      USE EIRMOD_CGEOM
      USE EIRMOD_COMPRT
      USE EIRMOD_CSDVI
      USE EIRMOD_COUTAU
      USE EIRMOD_COMXS
      USE EIRMOD_COMSPL
      USE EIRMOD_CLOGAU
      USE EIRMOD_CSPEZ
      USE EIRMOD_PHOTON
 
      IMPLICIT NONE
 
      REAL(DP), INTENT(IN) :: CFLAG(7,MSTOR0), DIST
      REAL(DP), INTENT(OUT) :: COLTYP
      REAL(DP) :: DUMT(3), DUMV(3)
      REAL(DP) :: ZEP1, SIGSUM, WGHTO, FRSTP, PTOT, E0O, VELXO,
     .          VELYO, VELZO, BX, BY, BZ, V0_PARBO, VELO, SCNDP,
     .          EDEL, VDEL, SIG, V0_PARB, FP, FLTEST, ZEP3, VELQ, VX,
     .          VY, VZ, VPLASP, RMAIO, RMMIO, RMIIO, BF, ZEP
      REAL(DP) :: SIG_ELIM, SIG_TOT_N, SIG_TOT_O, SIG_TEST
      INTEGER :: IICX, IIEI, IMEL, IOLD, NOLD, IACX, IRCX, IAEI, IREI,
     .           IBGK, IAEL, IREL, IP, IMEI, IMCX, IAPI, II, NFLAG,
     .           IATMN, IPLSN, IRPI, NCLLO, IPLSV, IMPI, IIPI, I, J, IPL
      INTEGER :: NEII_RED,LGEI_RED(0:NREI)

Cdr  additional arrays for  ANALOG CASCADE and SPLITTING AT COLLISIONS. 
Cdr (should be set in initialization phase, not here)
CDR  check: are the corresponding arrays PATEI,PMLEI, PIOEI real or integer (1/2 particle possible?)
      INTEGER, ALLOCATABLE, SAVE :: NAMIEI(:),NAMIPI(:)
 
 
csw add n 2lines
      INTEGER :: iaot,irot,kk,updf,t1,t2
      real(dp):: sump
csw external
      real(dp), external :: ranf_eirene
 
      SAVE
C
      ENTRY EIRENE_COLATM(CFLAG,COLTYP,DIST)
C
C  INCIDENT SPECIES: IOLD
      VELXO=VELX
      VELYO=VELY
      VELZO=VELZ
      VELO=VEL
      NCLLO = NCELL
      NCELL = NCLTAL(NCLLO)

C  PARALLEL MOMENTUM OF TEST PARTICLE INCIDENT TO COLLISION
      IF (LMAPL.OR.NGENA(IATM).NE.0) THEN
        CALL EIRENE_BFIELD (NCLLO, X0, Y0, Z0, BX, BY, BZ, BF,.TRUE.)
        V0_PARBO=VEL*(VELX*BX+VELY*BY+VELZ*BZ)
        V0_PARBO=V0_PARBO*AMUA*RMASSA(IATM)
      ENDIF

      E0O=E0
      WGHTO=WEIGHT
      IOLD=IATM
      NOLD=NSPH+IATM

      IF (IMETCL(NCELL) == 0) THEN
        NCLMT = NCLMT+1
        ICLMT(NCLMT) = NCELL
        IMETCL(NCELL) = NCLMT
      END IF
C
C  ABSORPTION BIASSING: CURRENTLY ONLY IMPLEMENTED FOR "EI-TYPE" (ELECTRON IMPACT) PROCESSES 

C  SUPPRESS THOSE IREI PROCESSES WITH ZERO
C                      TEST PARTICLE SECONDARIES
 
      SIG_ELIM=0.
      SIG_TOT_N=SIGTOT
      SIG_TOT_O=SIGTOT
      NEII_RED=0

      IF (WEIGHT.LT.WMINV) THEN
C  WEIGHT ALREADY TOO SMALL, NO SUPPRESION OF ABSORPTION
        NEII_RED=NAEII(IOLD)
        LGEI_RED(:)=LGAEI(IOLD,:)
      ELSE
C  TRY TO SUPPRESS ABSORPTION. IDENTIFY POSSIBLE EI PROCESSES
C                              WITH ZERO TEST PARTICLE SECONDARIES

        DO IAEI=1,NAEII(IOLD)
          IREI=LGAEI(IOLD,IAEI)
C  WHILE BEING IN THIS LOOP WEIGHT MAY BE REPEATEDLY REDUCED, FOR EARLIER (LOWER) IAEI
          IF (WEIGHT.GT.WMINV) THEN
C  REMAINING RATE AFTER POSSIBLE ELIMINATION OF IREI
C  SIG_TEST=0 WOULD VIOLATE RADON-NYKODYM CONDITION OF WEIGHTING
            SIG_TEST=SIG_TOT_N-SIGVEI(IREI)
            PTOT=P2NDS(IREI)
            IF (PTOT.EQ.0..AND.SIG_TEST.GT.0.) THEN
C  IREI IS A PURELY ABSORBING PROCESS
C  ELIMINATE THIS PROCESS IREI FROM ALL NAEII POSSIBLE EI PROCESSES
C  REDUCE WEIGHT ACCORDINGLY
              SIG_ELIM=SIG_ELIM+SIGVEI(IREI)
              SIG_TOT_N=SIG_TEST
              WEIGHT=WEIGHT*SIG_TOT_N/SIG_TOT_O
              WGHTO=WEIGHT
              SIG_TOT_O=SIG_TOT_N
              IF (IESTEI(IREI,1).NE.0) GOTO 990
              IF (IESTEI(IREI,2).NE.0) GOTO 990
              IF (IESTEI(IREI,3).NE.0) GOTO 990
            ELSE
C  NO, THIS PROCESS REMAINS ACTIVE, BECAUSE THERE ARE TEST-PARTICLE SECONDARIES
              NEII_RED=NEII_RED+1
              LGEI_RED(NEII_RED)=IREI
            ENDIF
          ELSE
C  WEIGHT TOO SMALL COMPARED TO WMINV. ANALOG GAME
            NEII_RED=NEII_RED+1
            LGEI_RED(NEII_RED)=IREI
          ENDIF
        ENDDO
      ENDIF  ! SUPPRESSION OF ABSORPTION AT EI PROCESSES: DONE.

C  WEIGHT MAY HAVE BEEN REDUCED NOW, AND ALSO THE NUMBER OF ACTIVE EI PROCESSES.
C  
C
C  FIRST DECIDE: ELECTRON IMPACT (COLLISION TYPE: EI) OR OTHER PROCESS
C
      ZEP1=SIG_ELIM+RANF_EIRENE( )*SIG_TOT_N
      SIGSUM=SIG_ELIM
C
      IF (ZEP1.LE.SIGEIT) THEN
C
C  AT THIS POINT: NEII.GE.1, FOR OTHERWISE ZEP1 COULD NOT HAVE
C                 POINTED TO EI-PROCESSES
C
C  ELECTRON IMPACT COLLISION:
C
        IF (NLTRC) CALL EIRENE_CHCTRC(X0,Y0,Z0,16,2)
        IF (NLSTOR) CALL EIRENE_STORE(2)
C
C  FIND TYP OF ELECTR. IMPACT COLLISION PROCESS: IREI
        DO 240 IAEI=1,NEII_RED-1
          IREI=LGEI_RED(IAEI)
          SIGSUM=SIGSUM+SIGVEI(IREI)
          IF (ZEP1.LE.SIGSUM) GOTO 245
240     CONTINUE
        IREI=LGEI_RED(NEII_RED)
245     CONTINUE
C
C  CALCULATE WEIGHT OF THE NEXT GENERATION PARTICLE
C  ONLY ONE ATOM, MOLECULE OR TEST-ION HISTORY WITH MODIFIED WEIGHT
C  IS FOLLOWED. 
C  PTOT IS THE (INTEGER) NUMBER OF ANALOG NEXT GENERATION TEST PARTICLES
C
        PTOT=P2NDS(IREI)
C       PTOTAL=PTOT+PPLEI(IREI,0)
C  ABSORBED WEIGHT: WEIABS
C       WEIABS=WEIGHT*PPLEI(IREI,0)
C
C  PRE-COLLISION ESTIMATOR FOR EAAT, 
C  PRE- AND POST COLLISION ESTIMATOR FOR EAPL AND EAEL
        IF (IESTEI(IREI,3).NE.0) THEN
C  score loss of incoming test particle energy
          IF (LEAAT) EAAT(NCELL)=EAAT(NCELL)-WEIGHT*E0

cdr EAPL, EAEL       :  SCORE NET CHANGES HERE.
cdr EAAT, EAML, EAIO :  SCORE EXACT GAINS LATER. 
          IF (LEAPL) THEN
            DO IP=1,IPPLEI(IREI,0)
cdr: this is incorrect. esigei must be split into ipl secondaries
cdr  it only happens to be correct if the post collision bulk species are all the same (=ipl),
cdr  because then esigei is the total for this species.
              IPL=IPPLEI(IREI,IP)
              LOGPLS(IPL,ISTRA)=.TRUE.
              EAPL(IPL,NCELL)=EAPL(IPL,NCELL)+WEIGHT*ESIGEI(IREI,4)
              LMETSP(NSPAMI+IPL)=.TRUE.
            END DO
          END IF
          IF (LEAEL) EAEL(NCELL)=EAEL(NCELL)+WEIGHT*ESIGEI(IREI,5)
        ENDIF
C
C  ABSORBTION (INTO BULK SPECIES) IS SUPPRESSED
C  STRICTLY: A WMINV CRITERION MAY BE USED HERE FOR THIS PROCESS IREI AGAIN
C            WHEN THIS PROCESS RESULTS IN BOTH: TEST AND BULK SECONDARIES
C            ABOVE: ONLY PROCESSES IREI WITH ZERO TEST SECONDARIES MIGHT HAVE
C            BEEN SUPRESSED.
        WEIGHT=WEIGHT*PTOT
C
C  ARE THERE TEST PARTICLE SECONDARIES AT ALL?
        IF (WEIGHT.LE.EPS30) THEN
C  NO !
          LGPART=.FALSE.
          ITYP=4
          COLTYP=2
          NCELL=NCLLO
          RETURN
        ENDIF

Cdr  PTOT=0,1,2,etc..., = integer number of next generation particles

c.......................................................................
        IF (NLCASCAD .AND. (NLEVEL+PTOT <= MAXLEVEL)) THEN  !  EI PROCESS CASCADING  ATM
cdr  ANALOGUE SAMPLING, I.E. SPLITTING, IN CASE OF MORE THAN ONE SECONDARY.

          IF (.NOT.ALLOCATED(NAMIEI)) THEN
            ALLOCATE(NAMIEI(NSPAMI))
          END IF
cdr  build one single distribution of secondary test particle species, all types, include photons
cdr  this should not be done here, but instead only once, in preproc. phase !!
cdr  this NAMIEI is the underlying discrete pdf, which led to the normalized cummulative p2nd(IREI) ?
          NAMIEI = 0

          NAMIEI(1:NSPH)         = 0    !  PPHEI(IREI,1:NPHOTI) IS NOT YET SET IN XSTEI.F
          NAMIEI(NSPH+1:NSPA)    = PATEI(IREI,1:NATMI)
          NAMIEI(NSPA+1:NSPAM)   = PMLEI(IREI,1:NMOLI)
          NAMIEI(NSPAM+1:NSPAMI) = PIOEI(IREI,1:NIONI)

!  RESET WEIGHT BACK TO ORIGINAL VALUE
          WEIGHT=WEIGHT / PTOT

cdr  generate secondaries, one by one, call veloei, and store them on splitting arrays

          DO I = NSPAMI, NSPH+1, -1  ! LOOP OVER ALL POTENTIAL SECONDARY SPECIES 'I'
            DO J=1, NAMIEI(I)   ! THERE ARE NAMIEI(I) COPIES OF THIS SECONDARY 'I'
C  FIND A "RANDOM NUMBER" TO ENFORCE "SAMPLING" OF THIS PARTICULAR SPECIES 'I' IN VELOEI
              ZEP = 0.5_DP * (P2ND(IREI,I-1)+P2ND(IREI,I))
              CALL EIRENE_VELOEI(NCLLO,IREI,VELXO,VELYO,VELZO,VELO,ZEP)
              ISPZ = ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
C
C  SPLITTING: EACH SECONDARY IS A NEW SPLITTING LEVEL.
C
              NLEVEL=NLEVEL+1
C  SAVE LOCATION, WEIGHT AND OTHER PARAMETERS AT CURRENT LEVEL
              RSPLST(1:NPARTC,NLEVEL)=RPST(1:NPARTC)
              ISPLST(1:MPARTC,NLEVEL)=IPST(1:MPARTC)
C  NUMBER OF NODES AT THIS LEVEL
              NODES(NLEVEL)=2  !  CDR  ONE PARTICLE SCORE IN EACH LEVEL

              IF (NLTRC) THEN 
                WRITE (IUNOUT,*) 'SPLITTING IN COLATM, EI PROCESS '
                WRITE (IUNOUT,*) 'STORE ', TEXTS(ISPZ)
              ENDIF
            END DO
          END DO  ! LOOP OVER ALL SECONDARIES DONE
C  FOR ALL SECONDARIES WE HAVE CALLED VELOEL, AND STORED POST COLLISION PARAMETERS
C  ON SPLITTING ARRAYS.

!  REMOVE LAST PARTICLE FROM STORAGE AS IT'S TRAJECTORY IS CONTINUED
          NLEVEL = NLEVEL - 1
          IF (NLTRC) WRITE(IUNOUT,*) 'REMOVE FROM STORAGE ', TEXTS(ISPZ)

CDR:   VELOEI FOR THIS CONTINUED PARTICLE HAS ALREADY BEEN CALLED

          ELSE  ! NOT ENOUGH STORAGE FOR CASCADING

          IF (NLCASCAD) THEN
            WRITE (iunout,*) 
     .        'ANALOG CALCULATION ABANDONED FOR PART. NO. ',NPANU
            WRITE (iunout,*) 'CASCADE OVERFLOW: NEVEL: ',NLEVEL
          ENDIF

CDR:  (NORMAL) NON-ANALOG GAME AT EI PROCESSES

          CALL EIRENE_VELOEI(NCLLO,IREI,VELXO,VELYO,VELZO,VELO,-1._DP)

        END IF

        XGENER=0.D0
C
C  UPDATE POST-COLLISION ESTIMATORS CONTRIBUTION TO EAAT;EAML;EAIO
C         ACCOUNT FOR POST COLLISION CONTRIBUTIONS
        IF (ITYP.EQ.1) THEN
          IF (IESTEI(IREI,3).NE.0) THEN
            IF (LEAAT) EAAT(NCELL)=EAAT(NCELL)+WEIGHT*E0
          ENDIF
          COLTYP=1
        ELSEIF (ITYP.EQ.2) THEN
          IF (IESTEI(IREI,3).NE.0) THEN
            IF (LEAML) EAML(NCELL)=EAML(NCELL)+WEIGHT*E0
          ENDIF
          COLTYP=1
        ELSEIF (ITYP.EQ.3) THEN
          IF (IESTEI(IREI,3).NE.0) THEN
            IF (LEAIO) EAIO(NCELL)=EAIO(NCELL)+WEIGHT*E0
          ENDIF
          COLTYP=2
        ENDIF
        NCELL = NCLLO
        RETURN
C
      ELSEIF (ZEP1.LE.SIGEIT+SIGCXT) THEN
C
C  CHARGE EXCHANGE:
C
        IF (NLTRC) CALL EIRENE_CHCTRC(X0,Y0,Z0,16,6)
C
C   FIND PROCESS IRCX AND SPECIES INDEX IPLS OF INCIDENT BULK ION
        SIGSUM=SIGEIT
        DO 271 IACX=1,NACXIM(IATM)
          IRCX=LGACX(IATM,IACX,0)
          IPLS=LGACX(IATM,IACX,1)
          SIGSUM=SIGSUM+SIGVCX(IRCX)
          IF (ZEP1.LT.SIGSUM) GOTO 272
271     CONTINUE
        IRCX=LGACX(IATM,NACXI(IATM),0)
        IPLS=LGACX(IATM,NACXI(IATM),1)
272     CONTINUE
C
C  ARE THERE SECONDARY TEST PARTICLES AT ALL?
        FRSTP=N1STX(IRCX,3)
        SCNDP=N2NDX(IRCX,3)
 
        IPLSV=MPLSV(IPLS)
C
C  ARE THERE SECONDARY TEST PARTICLES AT ALL?

        IF (SCNDP.LE.EPS30) THEN
C  POST COLLISION ESTIMATOR FOR PAPL,EAPL,MAPL: TO BE WRITTEN
C  E.G. FOR CX RECOMBINATION
          LGPART=.FALSE.
          IF (IESTCX(IRCX,1).NE.0) GOTO 999
          IF (IESTCX(IRCX,2).NE.0) GOTO 999
          IF (IESTCX(IRCX,3).NE.0) GOTO 999
          ITYP=4
          COLTYP=2
          NCELL = NCLLO
          RETURN
        ENDIF

        IF (NLCASCAD .AND. NLEVEL < MAXLEVEL) THEN  ! CX PROCESS CASCADING ATM
! JUST OPPOSITE TO EI CASE:
CDR IN EI CASE: LAST TEST SECONDARY WAS FOLLOWED, ALL OTHERS STORED ON SPLITTING ARRAY.
CDR IN CX CASE: OPPOSITE.   TRY TO UNIFY !!
C NORMALLY THERE IS ONLY ONE SECONDARY, AND WE FOLLOW (THIS ONLY) TEST SECONDAY
C HERE FIRST AND LAST HAVE A SPECIAL MEANING (EXCHANGE OF IDENTITY, VELOCITIES, ETC..)


C  STORE 2ND SECONDARY, SPEED OF PREVIOUS TEST PARTICLE,
c  (i.e. scattering angle = PI), energy may have changed.

          ITYP=N2NDX(IRCX,1)
          NCELL = NCLLO
          XGENER=0.D0

          IF (ITYP /= 4) THEN
            SELECT CASE (ITYP)
C
            CASE (1)
              IATM=N2NDX(IRCX,2)
              E0=CVRSSA(IATM)*VELO*VELO
 
            CASE (2)
              IMOL=N2NDX(IRCX,2)
              E0=CVRSSM(IMOL)*VELO*VELO
 
            CASE(3)
              IION=N2NDX(IRCX,2)
              E0=CVRSSI(IION)*VELO*VELO
 
            CASE DEFAULT
              WRITE (iunout,*) ' ITYP ',ITYP,' AS 2ND SECONDARY IS NOT',
     .                    ' FORESEEN IN COLLIDE '
            END SELECT

            ISPZ = ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
C
C.....................................................................
C  SPLITTING
C
            NLEVEL=NLEVEL+1
C  SAVE LOCATION, WEIGHT AND OTHER PARAMETERS AT CURRENT LEVEL
            RSPLST(1:NPARTC,NLEVEL)=RPST(1:NPARTC)
            ISPLST(1:MPARTC,NLEVEL)=IPST(1:MPARTC)
C  NUMBER OF NODES AT THIS LEVEL
            NODES(NLEVEL)=2

            IF (NLTRC) 
     .        WRITE (IUNOUT,*) 'CX CASCADING: STORE ', TEXTS(ISPZ)

          ENDIF  !  splitting done.

C  FOLLOW 1ST SECONDARY
          
          ZEP3 = 0.5*FRSTP

        ELSE ! NOT ENOUGH STORAGE FOR CASCADING

          IF (NLCASCAD) THEN
            WRITE (iunout,*) 
     .        'ANALOG CALCULATION ABANDONED FOR PART. NO. ',NPANU
            WRITE (iunout,*) 'CX CASCADE OVERFLOW: NEVEL: ',NLEVEL
          ENDIF

C  FROM HERE: OLD GAME, NO SPLITTING

C  SUPPRESSION OF ABSORBTION AT CX
C  SUPPRESSION OF CASCADING  AT CX
C  NO RANDOM DECISION BETWEEN BULK AND TEST SECONDARIES, BUT WEIGHTING
          WEIGHT=WEIGHT*SCNDP
          ZEP3=RANF_EIRENE( )*SCNDP

        END IF
C
C  NEW SPECIES TYPE, INDEX AND ENERGY


        IF (ZEP3.LE.FRSTP) THEN
C  FOLLOW FIRST SECONDARY, SPEED FROM BULK POPULATION
          ITYP=N1STX(IRCX,1)
          NFLAG=CFLAG(3,IRCX)
          CALL EIRENE_VELOCX
     .         (NCLLO,VELXO,VELYO,VELZO,VELO,IOLD,NOLD,VELQ,
     .          NFLAG,IRCX,DUMT,DUMV)
 
          SELECT CASE(ITYP)
C
          CASE(1)
C  1ST SECONDARY IS ATOM: IATM
            IATM=N1STX(IRCX,2)
            E0=CVRSSA(IATM)*VELQ
C

C  CX GENERATION LIMIT, ATOMS, IATM
            IF (NGENA(IATM).GT.0) THEN
              IF (IATM.EQ.IOLD) THEN
                XGENER=XGENER+1.D0
              ELSE
                XGENER=0.D0
              ENDIF
              IF (XGENER.GE.NGENA(IATM)) THEN
C  UPDATE GENERATION LIMIT TALLIES, THEN STOP TRAJECTORY
C  USE POST COLLISION WEIGHT, VELOCITY AND ENERGY (NOT: PRE COLLISION DATA)
C  SHOULD MAKE NO DIFFERENCE ON AVERAGE, IF GENERATION LIMIT IS VALID.
C  IF NOT, ONLY THIS FORM OF ABSORPTION ESTIMATOR GIVES CORRECT BALANCES.
c               write (iunout,*) 'particle stopped at generation limit ',
c    .                            npanu,xgener
                IF (LPGENA) PGENA(IATM,NCELL)=PGENA(IATM,NCELL)-WEIGHT
                IF (LEGENA)
     .            EGENA(IATM,NCELL)=EGENA(IATM,NCELL)-WEIGHT*E0
                IF (LVGENA) THEN
C  FIND POST COLLISION PARALLEL VELOCITY. LOCAL BFIELD: KNOWN ALREADY FORM INITIALISATION
                  V0_PARB=VEL*(VELX*BX+VELY*BY+VELZ*BZ)
                  V0_PARB=V0_PARB*AMUA*RMASSA(IATM)
                  VGENA(IATM,NCELL)=VGENA(IATM,NCELL)-WEIGHT*V0_PARB
                END IF
                LGPART=.FALSE.
                IF (LPGENA.OR.LEGENA.OR.LVGENA) LMETSP(NSPH+IATM)=.TRUE.
                IF (NLTRC) CALL EIRENE_CHCTRC(X0,Y0,Z0,16,16)
                ITYP=4
                COLTYP=2
                NCELL = NCLLO
                RETURN
              ENDIF
            ENDIF
C
C  FLUID LIMIT
            IF (NGENA(IATM).LT.0) THEN
              FP=VELO/SIGVCX(IRCX)  !mfp
              FLTEST=FP/DIST
              IF (FLTEST.LT.FDLMCX(IRCX)) THEN
C  UPDATE FLUID LIMIT TALLIES
C  USE POST COLLISION WEIGHT, VELOCITY AND ENERGY
C  SHOULD MAKE NO DIFFERENCE ON AVERAGE, IF GENERATION LIMIT IS VALID.
C  IF NOT, ONLY THIS GIVES CORRECT BALANCES.
                IF (LPGENA) PGENA(IATM,NCELL)=PGENA(IATM,NCELL)-WEIGHT
                IF (LEGENA)
     .            EGENA(IATM,NCELL)=EGENA(IATM,NCELL)-WEIGHT*E0
                IF (LVGENA) THEN
                  V0_PARB=VEL*(VELX*BX+VELY*BY+VELZ*BZ)
                  V0_PARB=V0_PARB*AMUA*RMASSA(IATM)
                  VGENA(IATM,NCELL)=VGENA(IATM,NCELL)-WEIGHT*V0_PARB
                END IF
                LGPART=.FALSE.
                IF (LPGENA.OR.LEGENA.OR.LVGENA) LMETSP(NSPH+IATM)=.TRUE.
                IF (NLTRC) CALL EIRENE_CHCTRC(X0,Y0,Z0,16,17)
                ITYP=4
                COLTYP=2
                NCELL = NCLLO
                RETURN
              ENDIF
            ENDIF
C
C  NEXT LINES: COLLISION ESTIMATOR FOR CHARGE EXCHANGE NO. IRCX
C  CONSERVE CHARGE IN EACH COLLISION, NOT ONLY ON AVERAGE
C
            IF (IESTCX(IRCX,1).NE.0) THEN
C  IATMN: ATOM SPECIES AFTER CX
              IATMN=IATM
              IF (LPAAT) THEN
                PAAT(IOLD,NCELL) =PAAT(IOLD,NCELL)-WGHTO
                PAAT(IATMN,NCELL)=PAAT(IATMN,NCELL)+WEIGHT
                LMETSP(NSPH+IOLD)=.TRUE.
                LMETSP(NSPH+IATMN)=.TRUE.
              END IF
              IF (LPAPL) THEN
                PAPL(IPLS,NCELL) =PAPL(IPLS,NCELL)-WEIGHT
                LMETSP(NSPAMI+IPLS)=.TRUE.
              END IF
              IF (LPAEL) PAEL(NCELL)      =PAEL(NCELL)-WEIGHT
              IF (N2NDX(IRCX,1).EQ.4) THEN
C  IPLSN: ION SPECIES AFTER CX
                IPLSN=N2NDX(IRCX,2)
                IF (LPAPL) THEN
                  PAPL(IPLSN,NCELL)=PAPL(IPLSN,NCELL)+WGHTO
                  LMETSP(NSPAMI+IPLSN)=.TRUE.
                END IF
                IF (LPAEL) PAEL(NCELL)      =PAEL(NCELL)+WGHTO
              ELSEIF (N2NDX(IRCX,1).NE.4) THEN
                GOTO 999
              ENDIF
            ENDIF
c  UPDATE collision estimator for CX energy exchange tallies
            IF (IESTCX(IRCX,3).NE.0) THEN
              IF (LEAAT) EAAT(NCELL)=EAAT(NCELL)-E0O*WGHTO
              IF (LEAAT) EAAT(NCELL)=EAAT(NCELL)+E0*WEIGHT
              IF (LEAPL) THEN
                EAPL(IPLS,NCELL)=EAPL(IPLS,NCELL)-E0*WEIGHT
                LMETSP(NSPAMI+IPLS)=.TRUE.
              END IF
              IF (N2NDX(IRCX,1).EQ.4) THEN
                IF (LEAPL) THEN
                  IPLSN=N2NDX(IRCX,2)
                  EAPL(IPLSN,NCELL)=EAPL(IPLSN,NCELL)+E0O*WGHTO
                  LMETSP(NSPAMI+IPLSN)=.TRUE.
                ENDIF
              ELSE
                GOTO 999
              ENDIF
            ENDIF
C  UPDATE COLLISION ESTIMATOR CONTRIBUTION TO MAPL (FORMERLY: COPV)
            IF (IESTCX(IRCX,2).NE.0) THEN
              IF (LMAPL) THEN
C  SET THE POST COLLISION NEUTRAL PARALLEL VELOCITY = OLD PRE COLLISION BULK (ION) VELOCITY
                V0_PARB=VEL*(VELX*BX+VELY*BY+VELZ*BZ)
                V0_PARB=V0_PARB*AMUA*RMASSA(IATM)
                IF (INDPRO(4) == 8) THEN
                  CALL EIRENE_VECUSR(2,NCELL,X0,Y0,Z0,VX,VY,VZ,IPLS,
     .                               .TRUE.)
                  VPLASP=VX*BX+VY*BY+VZ*BZ
                ELSE
                  VPLASP = BVIN(IPLSV,NCLLO)
                ENDIF
                SIG=SIGN(1._DP,VPLASP)
C  ASSUME: OLD (INCIDENT) ION MOMENTUM IS EQUAL TO NEW ATOM MOMENTUM
                MAPL(IPLS,NCELL)=MAPL(IPLS,NCELL)-WEIGHT*V0_PARB*SIG
                LMETSP(NSPAMI+IPLS)=.TRUE.
              END IF

              IF (N2NDX(IRCX,1).EQ.4) THEN
                IF (LMAPL) THEN
C  IPLSN ION SPECIES AFTER CX
                  IPLSN=N2NDX(IRCX,2)
C  ASSUME: NEW ION MOMENTUM IS EQUAL TO INCIDENT ATOM MOMENTUM
                  MAPL(IPLSN,NCELL)=MAPL(IPLSN,NCELL)+
     .                              WGHTO*V0_PARBO*SIG
                  LMETSP(NSPAMI+IPLSN)=.TRUE.
                END IF
              ELSEIF (N2NDX(IRCX,1).NE.4) THEN
                GOTO 999
              ENDIF
            ENDIF
            COLTYP=1
            NCELL=NCLLO
            RETURN

          CASE(2)
C  1ST SECONDARY IS MOLECULE
            IMOL=N1STX(IRCX,2)
            E0=CVRSSM(IMOL)*VELQ
            XGENER=0.D0
C
            IF (IESTCX(IRCX,1).NE.0) GOTO 999
            IF (IESTCX(IRCX,2).NE.0) GOTO 999
            IF (IESTCX(IRCX,3).NE.0) GOTO 999
            COLTYP=1
            NCELL = NCLLO
            RETURN

          CASE(3)
C  1ST SECONDARY IS TEST ION
            IION=N1STX(IRCX,2)
            E0=CVRSSI(IION)*VELQ
            XGENER=0.D0
C
            IF (IESTCX(IRCX,1).NE.0) GOTO 999
            IF (IESTCX(IRCX,2).NE.0) GOTO 999
            IF (IESTCX(IRCX,3).NE.0) GOTO 999
            COLTYP=2
            NCELL = NCLLO
            RETURN

          CASE DEFAULT
            WRITE (iunout,*) ' ITYP = ',ITYP,' AS FIRST SECONDARY IS',
     .                  ' NOT FORESEEN IN COLLIDE '
          END SELECT

        ELSE
 
C  FOLLOW 2ND SECONDARY, SPEED OF PREVIOUS TEST PARTICLE
          ITYP=N2NDX(IRCX,1)

          SELECT CASE(ITYP)
C
          CASE(1)
            IATM=N2NDX(IRCX,2)
            XGENER= 0.D0
C
            E0=CVRSSA(IATM)*VELO*VELO
            IF (IESTCX(IRCX,1).NE.0) GOTO 999
            IF (IESTCX(IRCX,2).NE.0) GOTO 999
            IF (IESTCX(IRCX,3).NE.0) GOTO 999
            COLTYP=1
            NCELL = NCLLO
            RETURN
C
          CASE(2)
            IMOL=N2NDX(IRCX,2)
            XGENER= 0.D0
C
            E0=CVRSSM(IMOL)*VELO*VELO
            IF (IESTCX(IRCX,1).NE.0) GOTO 999
            IF (IESTCX(IRCX,2).NE.0) GOTO 999
            IF (IESTCX(IRCX,3).NE.0) GOTO 999
            COLTYP=1
            NCELL = NCLLO
            RETURN
C
          CASE(3)
            IION=N2NDX(IRCX,2)
            XGENER=0.D0
C
            E0=CVRSSI(IION)*VELO*VELO
            IF (IESTCX(IRCX,1).NE.0) GOTO 999
            IF (IESTCX(IRCX,2).NE.0) GOTO 999
            IF (IESTCX(IRCX,3).NE.0) GOTO 999
            COLTYP=2
            NCELL = NCLLO
            RETURN

          CASE DEFAULT
            WRITE (iunout,*) ' ITYP = ',ITYP,' AS SECOND SECONDARY IS',
     .                  ' NOT FORESEEN IN COLLIDE '
          END SELECT

        ENDIF
C
C  ELASTIC COLLISION
C
      ELSEIF (ZEP1.LE.SIGEIT+SIGCXT+SIGELT) THEN
C
        IF (NLTRC) CALL EIRENE_CHCTRC(X0,Y0,Z0,16,5)
C   FIND SPECIES INDEX OF BULK (ION) COLLISION PARTNER
        SIGSUM=SIGEIT+SIGCXT
        DO 281 IAEL=1,NAELIM(IATM)
          IREL=LGAEL(IATM,IAEL,0)
          IPLS=LGAEL(IATM,IAEL,1)
          SIGSUM=SIGSUM+SIGVEL(IREL)
          IF (ZEP1.LT.SIGSUM) GOTO 282
281     CONTINUE
        IREL=LGAEL(IATM,NAELI(IATM),0)
        IPLS=LGAEL(IATM,NAELI(IATM),1)
282     CONTINUE
 
        IPLSV=MPLSV(IPLS)
C
C  NEW SPECIES INDEX AND ENERGY
C       WEIGHT=WEIGHT*1.
C  FOLLOW SECONDARY, NEW SPEED FROM SUBROUTINE VELOEL
C       ITYP=1
        NFLAG=CFLAG(5,IREL)
        RMAIO=RMASSA(IOLD)
        CALL EIRENE_VELOEL(NCLLO,VELXO,VELYO,VELZO,VELO,IOLD,NOLD,VELQ,
     .              NFLAG,IREL,RMAIO)
C
        IATM=IOLD
C  NOT: WEIGHT=WGHTO, BECAUSE WEIGHT MAY HAVE CHANGED DUE TO NON-ANALOGUE SAMPLING IN VELOEL
        E0=CVRSSA(IATM)*VELQ


C  DO NOT UPDATE BGK TALLIES HERE
        IBGK=NPBGKP(IPLS,1)
        IF (IBGK.NE.0) GOTO 300

C  UPDATE COLLISION ESTIMATOR CONTRIBUTION
C  ASSUME, AS BEFORE, NO CHANGE IN SPECIES/TYP
        IF (IESTEL(IREL,1).NE.0) THEN
          IF (LPAAT) THEN
            PAAT(IOLD,NCELL) =PAAT(IOLD,NCELL)-WGHTO
            PAAT(IATM,NCELL) =PAAT(IATM,NCELL)+WEIGHT
            LMETSP(NSPH+IOLD)=.TRUE.
            LMETSP(NSPH+IATM)=.TRUE.
          END IF
        ENDIF
c  UPDATE collision estimator for EL energy exchange tallies
        IF (IESTEL(IREL,3).NE.0) THEN
          EDEL=E0O*WGHTO-E0*WEIGHT
          IF (LEAAT) EAAT(NCELL)      =EAAT(NCELL)-EDEL
          IF (LEAPL) THEN
            EAPL(IPLS,NCELL) =EAPL(IPLS,NCELL)+EDEL
            LMETSP(NSPAMI+IPLS)=.TRUE.
          END IF
        ENDIF
C  UPDATE COLLISION ESTIMATOR CONTRIBUTION TO MAPL (FORMERLY: COPV)
        IF (IESTEL(IREL,2).NE.0) THEN
          IF (LMAPL) THEN
C  SET THE POST COLLISION NEUTRAL PARALLEL VELOCITY 
            V0_PARB=VEL*(VELX*BX+VELY*BY+VELZ*BZ)
            V0_PARB=V0_PARB*AMUA*RMASSA(IATM)
C
            VDEL=V0_PARBO*WGHTO-V0_PARB*WEIGHT
            IF (INDPRO(4) == 8) THEN
              CALL EIRENE_VECUSR(2,NCELL,X0,Y0,Z0,VX,VY,VZ,IPLS,
     .                           .TRUE.)
              VPLASP=VX*BX+VY*BY+VZ*BZ
            ELSE
              VPLASP=BVIN(IPLSV,NCLLO)
            ENDIF
            SIG=SIGN(1._DP,VPLASP)
            MAPL(IPLS,NCELL)=MAPL(IPLS,NCELL)+VDEL*SIG
            LMETSP(NSPAMI+IPLS)=.TRUE.
          END IF
        ENDIF
300     CONTINUE
        COLTYP=1
        NCELL = NCLLO
        RETURN
C
C  GENERAL ION IMPACT COLLISION: PI-PROCESSES. NOT READY
C
      ELSEIF (ZEP1.LE.SIGEIT+SIGCXT+SIGELT+SIGPIT) THEN
C
        IF (NLTRC) CALL EIRENE_CHCTRC(X0,Y0,Z0,16,3)
        SIGSUM=SIGEIT+SIGCXT+SIGELT
        DO 261 IAPI=1,NAPIIM(IATM)
C   FIND INDEX OF THAT ION IMPACT COLLISION
          IRPI=LGAPI(IATM,IAPI,0)
          IPLS=LGAPI(IATM,IAPI,1)
          SIGSUM=SIGSUM+SIGVPI(IRPI)
          IF (ZEP1.LT.SIGSUM) GOTO 262
261     CONTINUE
        IRPI=LGAPI(IATM,NAPII(IATM),0)
        IPLS=LGAPI(IATM,NAPII(IATM),1)
262     CONTINUE
C
C  CALCULATE WEIGHT OF THE NEXT GENERATION PARTICLE
C  ONLY ONE ATOM, MOLECULE OR TEST-ION HISTORY WITH MODIFIED WEIGHT
C  IS FOLLOWED
C
        PTOT=P2NPI(IRPI)
C       PTOTAL=PTOT+PPLPI(IRPI,0)
C  ABSORBED WEIGHT: WEIABS
C       WEIABS=WEIGHT*PPLPI(IRPI,0)
C
C  PRE- COLLISION ESTIMATOR FOR EAAT,
C  PRE- AND POST COLLISION ESTIMATOR FOR EAPL AND EAEL
        IF (IESTPI(IRPI,3).NE.0) THEN
C  score loss of incoming test particle energy
          IF (LEAAT) EAAT(NCELL)=EAAT(NCELL)-WEIGHT*E0

cdr EAPL, EAEL       :  SCORE NET CHANGES HERE.
cdr EAAT, EAML, EAIO :  SCORE EXACT GAINS LATER. 
          IF (LEAPL) THEN
            DO IP=1,IPPLPI(IRPI,0)
cdr:  this is incorrect. esigpi must be split into ipl secondaries
              IPL=IPPLPI(IRPI,IP)
              LOGPLS(IPL,ISTRA)=.TRUE.
              EAPL(IPL,NCELL)=EAPL(IPL,NCELL)+WEIGHT*ESIGPI(IRPI,4)
              LMETSP(NSPAMI+IPL)=.TRUE.
            END DO
          END IF
          IF (LEAEL) EAEL(NCELL)=EAEL(NCELL)+WEIGHT*ESIGPI(IRPI,5)
        ENDIF
C
C  ABSORBTION (INTO BULK SPECIES) IS SUPPRESSED
        WEIGHT=WEIGHT*PTOT
C
C  ARE THERE TEST PARTICLE SECONDARIES AT ALL?
        IF (WEIGHT.LE.EPS30) THEN
          LGPART=.FALSE.
          ITYP=4
          COLTYP=2
          NCELL=NCLLO
          RETURN
        ENDIF
C
        NFLAG=CFLAG(4,IRPI)
        RMAIO=RMASSA(IOLD)

Cdr  PTOT=0,1,2,etc..., = integer,  number of next generation particles

        IF (NLCASCAD .AND. (NLEVEL+PTOT <= MAXLEVEL)) THEN  ! PI PROCESS CASCADING ATM

          IF (.NOT.ALLOCATED(NAMIPI)) THEN
            ALLOCATE(NAMIPI(NSPAMI))
          END IF
          NAMIPI = 0
          NAMIPI(NSPH+1:NSPA) = PATPI(IRPI,1:NATMI)
          NAMIPI(NSPA+1:NSPAM) = PMLPI(IRPI,1:NMOLI)
          NAMIPI(NSPAM+1:NSPAMI) = PIOPI(IRPI,1:NIONI)

!  RESET WEIGHT TO ORIGINAL VALUE
          WEIGHT=WEIGHT / PTOT

          DO I = NSPAMI, NSPH+1, -1
            DO J=1, NAMIPI(I)
              ZEP = 0.5_DP * (P2NP(IRPI,I-1)+P2NP(IRPI,I))
              CALL EIRENE_VELOPI(NCLLO,VELXO,VELYO,VELZO,VELO,IOLD,
     .                           NOLD,VELQ,NFLAG,IRPI,RMAIO,ZEP)
              ISPZ = ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
C
C.....................................................................
C  SPLITTING
C
              NLEVEL=NLEVEL+1
C  SAVE LOCATION, WEIGHT AND OTHER PARAMETERS AT CURRENT LEVEL
              RSPLST(1:NPARTC,NLEVEL)=RPST(1:NPARTC)
              ISPLST(1:MPARTC,NLEVEL)=IPST(1:MPARTC)
C  NUMBER OF NODES AT THIS LEVEL
              NODES(NLEVEL)=2

              IF (NLTRC) WRITE (IUNOUT,*) 'STORE ', TEXTS(ISPZ)
              
            END DO
          END DO

!  REMOVE LAST PARTICLE FROM STORAGE AS IT'S TRAJECTORY IS CONTINUED
          NLEVEL = NLEVEL - 1
          IF (NLTRC) WRITE(IUNOUT,*) 'REMOVE FROM STORAGE ', TEXTS(ISPZ)

        ELSE  ! NOT ENOUGH STORAGE FOR CASCADING

          IF (NLCASCAD) THEN
            WRITE (iunout,*) 
     .        'ANALOG CALCULATION ABANDONED FOR PART. NO. ',NPANU
            WRITE (iunout,*) 'CASCADE OVERFLOW: NEVEL: ',NLEVEL
          ENDIF

          CALL EIRENE_VELOPI(NCLLO,VELXO,VELYO,VELZO,VELO,IOLD,
     .                       NOLD,VELQ,NFLAG,IRPI,RMAIO,-1._DP)

        END IF

        XGENER=0.D0
C
C  UPDATE COLLISION ESTIMATORS CONTRIBUTION TO EAAT;EAML;EAIO
        IF (ITYP.EQ.1) THEN
          IF (IESTPI(IRPI,3).NE.0) THEN
            IF (LEAAT) EAAT(NCELL)=EAAT(NCELL)+WEIGHT*E0
          ENDIF
          COLTYP=1
        ELSEIF (ITYP.EQ.2) THEN
          IF (IESTPI(IRPI,3).NE.0) THEN
            IF (LEAML) EAML(NCELL)=EAML(NCELL)+WEIGHT*E0
          ENDIF
          COLTYP=1
        ELSEIF (ITYP.EQ.3) THEN
          IF (IESTPI(IRPI,3).NE.0) THEN
            IF (LEAIO) EAIO(NCELL)=EAIO(NCELL)+WEIGHT*E0
          ENDIF
          COLTYP=2
        ENDIF
        NCELL = NCLLO
        RETURN
C
C
      ELSE
C
C
        WRITE (iunout,*) 'ERROR IN COLATM, UNKNOWN TYPE OF COLLISION '
        CALL EIRENE_EXIT_OWN(1)
C
C
      ENDIF
      GOTO 999
C
      ENTRY EIRENE_COLMOL(CFLAG,COLTYP,DIST)
C
C  INCIDENT SPECIES: IOLD

      VELXO=VELX
      VELYO=VELY
      VELZO=VELZ
      VELO=VEL
      NCLLO = NCELL
      NCELL = NCLTAL(NCLLO)

C  PARALLEL MOMENTUM OF TEST PARTICLE INCIDENT TO COLLISION
      IF (LMMPL.OR.NGENM(IMOL).NE.0) THEN
        CALL EIRENE_BFIELD (NCLLO, X0, Y0, Z0, BX, BY, BZ, BF,.TRUE.)
        V0_PARBO=VEL*(VELX*BX+VELY*BY+VELZ*BZ)
        V0_PARBO=V0_PARBO*AMUA*RMASSM(IMOL)
      ENDIF

      E0O=E0
      WGHTO=WEIGHT
      IOLD=IMOL
      NOLD=NSPA+IMOL

      IF (IMETCL(NCELL) == 0) THEN
        NCLMT = NCLMT+1
        ICLMT(NCLMT) = NCELL
        IMETCL(NCELL) = NCLMT
      END IF
C
C  ABSORPTION BIASSING: CURRENTLY ONLY IMPLEMENTED FOR "EI-TYPE" (ELECTRON IMPACT) PROCESSES 

C  SUPPRESS THOSE IREI PROCESSES WITH ZERO
C                      TEST PARTICLE SECONDARIES: TO BE DONE, SEE ATOM PART.
 
      SIG_ELIM=0.
      SIG_TOT_N=SIGTOT
      SIG_TOT_O=SIGTOT
      NEII_RED=0
 
      DO IMEI=1,NMEII(IOLD)
        IREI=LGMEI(IOLD,IMEI)
        IF (WEIGHT.GT.WMINV) THEN
C  REMAINING RATE AFTER POSSIBLE ELIMINATION OF IREI
C  SIG_TEST=0 WOULD VIOLATE RADON-NYKODYM CONDITION OF WEIGHTING
          SIG_TEST=SIG_TOT_N-SIGVEI(IREI)
          PTOT=P2NDS(IREI)
          IF (PTOT.EQ.0..AND.SIG_TEST.GT.0.) THEN
C  IREI IS A PURELY ABSORBING PROCESS
C  ELIMINATE PROCESS IREI FROM ALL NMEII POSSIBLE PROCESSES
C  REDUCE WEIGHT ACCORDINGLY
            SIG_ELIM=SIG_ELIM+SIGVEI(IREI)
            SIG_TOT_N=SIG_TEST
            WEIGHT=WEIGHT*SIG_TOT_N/SIG_TOT_O
            WGHTO=WEIGHT
            SIG_TOT_O=SIG_TOT_N
            IF (IESTEI(IREI,1).NE.0) GOTO 990
            IF (IESTEI(IREI,2).NE.0) GOTO 990
            IF (IESTEI(IREI,3).NE.0) GOTO 990
          ELSE
C  NO, THIS PROCESS REMAINS ACTIVE, BECAUSE THERE ARE TEST-PARTICLE SECONDARIES
            NEII_RED=NEII_RED+1
            LGEI_RED(NEII_RED)=IREI
          ENDIF
        ELSE
C  WEIGHT TOO SMALL COMPARED TO WMINV. ANALOG GAME
          NEII_RED=NEII_RED+1
          LGEI_RED(NEII_RED)=IREI
        ENDIF
      ENDDO
C
C  FIRST DECIDE: ELECTRON IMPACT (COLLISION TYPE: EI) OR OTHER PROCESS
C
      ZEP1=SIG_ELIM+RANF_EIRENE( )*SIG_TOT_N
      SIGSUM=SIG_ELIM
C
      IF (ZEP1.LE.SIGEIT) THEN
C
C  AT THIS POINT: NEII.GE.1, FOR OTHERWISE ZEP1 COULD NOT HAVE
C                 POINTED TO EI-PROCESSES
C
C  ELECTRON IMPACT COLLISION:
C
        IF (NLTRC) CALL EIRENE_CHCTRC(X0,Y0,Z0,16,2)
        IF (NLSTOR) CALL EIRENE_STORE(2)
C
C  FIND TYP OF ELECTR. IMPACT COLLISION PROCESS: IREI
        DO 340 IMEI=1,NEII_RED-1
          IREI=LGEI_RED(IMEI)
          SIGSUM=SIGSUM+SIGVEI(IREI)
          IF (ZEP1.LE.SIGSUM) GOTO 345
340     CONTINUE
        IREI=LGEI_RED(NEII_RED)
345     CONTINUE
C
C  CALCULATE WEIGHT OF THE NEXT GENERATION PARTICLE
C  ONLY ONE ATOM, MOLECULE OR TEST-ION HISTORY WITH MODIFIED WEIGHT
C  IS FOLLOWED.
C  PTOT IS THE (INTEGER) NUMBER OF ANALOG NEXT GENERATION TEST PARTICLES
C
        PTOT=P2NDS(IREI)
C       PTOTAL=PTOT+PPLEI(IREI,0)
C  ABSORBED WEIGHT: WEIABS
C       WEIABS=WEIGHT*PPLEI(IREI,0)
C
C  PRE- COLLISION ESTIMATOR FOR EMML,
C  PRE- AND POST COLLISION ESTIMATOR FOR EMPL AND EMEL
        IF (IESTEI(IREI,3).NE.0) THEN
C  score loss of incoming test particle energy
          IF (LEMML) EMML(NCELL)=EMML(NCELL)-WEIGHT*E0

cdr EMPL, EMEL       :  SCORE NET CHANGES HERE.
cdr EMAT, EMML, EMIO :  SCORE EXACT GAINS LATER. 
          IF (LEMPL) THEN
            DO IP=1,IPPLEI(IREI,0)
cdr:  this is incorrect. esigei must be split into ipl secondaries
cdr  it only happens to be correct if the post collision bulk species are all the same (=ipl),
cdr  because then esigei is the total for this species.
              IPL=IPPLEI(IREI,IP)
              LOGPLS(IPL,ISTRA)=.TRUE.
              EMPL(IPL,NCELL)=EMPL(IPL,NCELL)+WEIGHT*ESIGEI(IREI,4)
              LMETSP(NSPAMI+IPL)=.TRUE.
            END DO
         END IF
          IF (LEMEL) EMEL(NCELL)=EMEL(NCELL)+WEIGHT*ESIGEI(IREI,5)
        ENDIF
C
C  ABSORBTION (INTO BULK SPECIES) IS SUPPRESSED
        WEIGHT=WEIGHT*PTOT
C
C  ARE THERE TEST PARTICLE SECONDARIES AT ALL?
        IF (WEIGHT.LE.EPS30) THEN
C  NO !
          LGPART=.FALSE.
          ITYP=4
          COLTYP=2
          NCELL = NCLLO
          RETURN
        ENDIF

Cdr  PTOT=0,1,2,etc..., = integer number of next generation particles

c.......................................................................
        IF (NLCASCAD .AND. (NLEVEL+PTOT <= MAXLEVEL)) THEN  ! EI PROCESS CASCADING MOL
cdr  ANALOG SAMPLING, I.E. SPLITTING, IN CASE OF MORE THAN ONE SECONDARY.

          IF (.NOT.ALLOCATED(NAMIEI)) THEN
            ALLOCATE(NAMIEI(NSPAMI))
          END IF
cdr  build one single distribution of secondary test particle species, all types, include photons
cdr  this should not be done here, but instead only once, in preproc. phase !!
cdr  this NAMIEI is the underlying discrete pdf, which led to the normalized cummulative p2nd(IREI) ?
          NAMIEI = 0

          NAMIEI(1:NSPH)         = 0    !  PPHEI(IREI,1:NPHOTI) IS NOT YET SET IN XSTEI.F
          NAMIEI(NSPH+1:NSPA)    = PATEI(IREI,1:NATMI)
          NAMIEI(NSPA+1:NSPAM)   = PMLEI(IREI,1:NMOLI)
          NAMIEI(NSPAM+1:NSPAMI) = PIOEI(IREI,1:NIONI)

!  RESET WEIGHT BACK TO ORIGINAL VALUE
          WEIGHT=WEIGHT / PTOT

cdr  generate secondaries, one by one, call veloei, and store them on splitting arrays

          DO I = NSPAMI, NSPH+1, -1  ! LOOP OVER ALL POTENTIAL SECONDARY SPECIES 'I'
            DO J=1, NAMIEI(I)   ! THERE ARE NAMIEI(I) COPIES OF THIS SECONDARY 'I'
C  FIND A "RANDOM NUMBER" TO ENFORCE "SAMPLING" OF THIS PARTICULAR SPECIES 'I' IN VELOEI
              ZEP = 0.5_DP * (P2ND(IREI,I-1)+P2ND(IREI,I))
              CALL EIRENE_VELOEI(NCLLO,IREI,VELXO,VELYO,VELZO,VELO,ZEP)
              ISPZ = ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
C
C.....................................................................
C  SPLITTING: EACH SECONDARY IS A NEW SPLITTING LEVEL.
C
              NLEVEL=NLEVEL+1
C  SAVE LOCATION, WEIGHT AND OTHER PARAMETERS AT CURRENT LEVEL
              RSPLST(1:NPARTC,NLEVEL)=RPST(1:NPARTC)
              ISPLST(1:MPARTC,NLEVEL)=IPST(1:MPARTC)
C  NUMBER OF NODES AT THIS LEVEL
              NODES(NLEVEL)=2

              IF (NLTRC) WRITE (IUNOUT,*) 'STORE ', TEXTS(ISPZ)
              
            END DO
          END DO

!  REMOVE LAST PARTICLE FROM STORAGE AS IT'S TRAJECTORY IS CONTINUED
          NLEVEL = NLEVEL - 1
          IF (NLTRC) WRITE(IUNOUT,*) 'REMOVE FROM STORAGE ', TEXTS(ISPZ)

CDR:   VELOEI FOR THIS CONTINUED PARTICLE HAS ALREADY BEEN CALLED

        ELSE  ! NOT ENOUGH STORAGE FOR CASCADING

          IF (NLCASCAD) THEN
            WRITE (iunout,*) 
     .        'ANALOG CALCULATION ABANDONED FOR PART. NO. ',NPANU
            WRITE (iunout,*) 'CASCADE OVERFLOW: NEVEL: ',NLEVEL
          ENDIF

CDR:  (NORMAL) NON-ANALOG GAME AT EI PROCESSES

          CALL EIRENE_VELOEI(NCLLO,IREI,VELXO,VELYO,VELZO,VELO,-1._DP)

        END IF

        XGENER=0.D0
C
C  UPDATE POST-COLLISION ESTIMATORS CONTRIBUTION TO EMAT;EMML;EMIO
C         ACCOUNT FOR POST COLLISION CONTRIBUTIONS
        IF (ITYP.EQ.1) THEN
          IF (IESTEI(IREI,3).NE.0) THEN
            IF (LEMAT) EMAT(NCELL)=EMAT(NCELL)+WEIGHT*E0
          ENDIF
          COLTYP=1
        ELSEIF (ITYP.EQ.2) THEN
          IF (IESTEI(IREI,3).NE.0) THEN
            IF (LEMML) EMML(NCELL)=EMML(NCELL)+WEIGHT*E0
          ENDIF
          COLTYP=1
        ELSEIF (ITYP.EQ.3) THEN
          IF (IESTEI(IREI,3).NE.0) THEN
            IF (LEMIO) EMIO(NCELL)=EMIO(NCELL)+WEIGHT*E0
          ENDIF
          COLTYP=2
        ENDIF
        NCELL = NCLLO
        RETURN
C
      ELSEIF (ZEP1.LE.SIGEIT+SIGCXT) THEN
C
C  CHARGE EXCHANGE:
C
        IF (NLTRC) CALL EIRENE_CHCTRC(X0,Y0,Z0,16,6)
C
C   FIND PROCESS IRCX AND SPECIES INDEX IPLS OF INCIDENT BULK ION
        SIGSUM=SIGEIT
        DO 371 IMCX=1,NMCXIM(IMOL)
          IRCX=LGMCX(IMOL,IMCX,0)
          IPLS=LGMCX(IMOL,IMCX,1)
          SIGSUM=SIGSUM+SIGVCX(IRCX)
          IF (ZEP1.LT.SIGSUM) GOTO 372
371     CONTINUE
        IRCX=LGMCX(IMOL,NMCXI(IMOL),0)
        IPLS=LGMCX(IMOL,NMCXI(IMOL),1)
372     CONTINUE
C
C  ARE THERE SECONDARY TEST PARTICLES AT ALL?
        FRSTP=N1STX(IRCX,3)
        SCNDP=N2NDX(IRCX,3)
 
        IPLSV=MPLSV(IPLS)
C
C  ARE THERE SECONDARY TEST PARTICLES AT ALL?

        IF (SCNDP.LE.EPS30) THEN
C  POST COLLISION ESTIMATOR FOR PMPL,EMPL,MMPL: TO BE WRITTEN
C  E.G. FOR CX RECOMBINATION
          LGPART=.FALSE.
          IF (IESTCX(IRCX,1).NE.0) GOTO 999
          IF (IESTCX(IRCX,2).NE.0) GOTO 999
          IF (IESTCX(IRCX,3).NE.0) GOTO 999
          ITYP=4
          COLTYP=2
          NCELL = NCLLO
          RETURN
        ENDIF

        IF (NLCASCAD .AND. NLEVEL < MAXLEVEL) THEN  ! CX PROCESS CASCADING MOL
! JUST OPPOSITE TO EI CASE:
CDR IN EI CASE: LAST TEST SECONDARY WAS FOLLOWED, ALL OTHERS STORED ON SPLITTING ARRAY.
CDR IN CX CASE: OPPOSITE.   TRY TO UNIFY !!
C NORMALLY ONLY ONE SECONDARY, AND FOLLOW (ONLY) TEST SECONDAY
C HERE FIRST AND LAST HAVE A SPECIAL MEANING (EXCHANGE OF IDENTITY, VELOCITIES, ETC..)


C  STORE 2ND SECONDARY, SPEED OF PREVIOUS TEST PARTICLE,
c  (i.e. scattering angle = PI), energy may have changed.

          ITYP=N2NDX(IRCX,1)
          NCELL = NCLLO
          XGENER=0.D0

          IF (ITYP /= 4) THEN
            SELECT CASE (ITYP)
C
            CASE(1)
              IATM=N2NDX(IRCX,2)
              E0=CVRSSA(IATM)*VELO*VELO
 
            CASE(2)
              IMOL=N2NDX(IRCX,2)
              E0=CVRSSM(IMOL)*VELO*VELO
 
            CASE(3)
              IION=N2NDX(IRCX,2)
              E0=CVRSSI(IION)*VELO*VELO
 
            CASE DEFAULT
              WRITE (iunout,*) ' ITYP ',ITYP,' AS 2ND SECONDARY IS NOT',
     .                    ' FORESEEN IN COLLIDE '
            END SELECT

            ISPZ = ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
C
C.....................................................................
C  SPLITTING
C
            NLEVEL=NLEVEL+1
C  SAVE LOCATION, WEIGHT AND OTHER PARAMETERS AT CURRENT LEVEL
            RSPLST(1:NPARTC,NLEVEL)=RPST(1:NPARTC)
            ISPLST(1:MPARTC,NLEVEL)=IPST(1:MPARTC)
C  NUMBER OF NODES AT THIS LEVEL
            NODES(NLEVEL)=2

            IF (NLTRC) 
     .        WRITE (IUNOUT,*) 'CX CASCADING: STORE ', TEXTS(ISPZ)

          ENDIF  !  splitting done.

C  FOLLOW 1ST SECONDARY
          
          ZEP3 = 0.5*FRSTP

        ELSE  ! NOT ENOUGH STORAGE FOR CASCADING

          IF (NLCASCAD) THEN
            WRITE (iunout,*) 
     .        'ANALOG CALCULATION ABANDONED FOR PART. NO. ',NPANU
            WRITE (iunout,*) 'CX CASCADE OVERFLOW: NEVEL: ',NLEVEL
          ENDIF

          WEIGHT=WEIGHT*SCNDP
          ZEP3=RANF_EIRENE( )*SCNDP

        END IF

C
C  NEW SPECIES TYPE, INDEX AND ENERGY
C  SUPPRESSION OF ABSORBTION AT CX
C  I.E., NO RANDOM DECISION BETWEEN BULK AND TEST SECONDARIES
        IF (ZEP3.LE.FRSTP) THEN
C  FOLLOW FIRST SECONDARY, SPEED FROM BULK POPULATION
          ITYP=N1STX(IRCX,1)
          NFLAG=CFLAG(3,IRCX)
          CALL EIRENE_VELOCX
     .         (NCLLO,VELXO,VELYO,VELZO,VELO,IOLD,NOLD,VELQ,
     .          NFLAG,IRCX,DUMT,DUMV)
 
          SELECT CASE(ITYP)
C
          CASE(1)
C  1ST SECONDARY IS ATOM: IATM
            IATM=N1STX(IRCX,2)
            E0=CVRSSA(IATM)*VELQ
            XGENER=0.D0
C
C  NEXT LINES: COLLISION ESTIMATOR FOR CHARGE EXCHANGE NO. IRCX
C  CONSERVE CHARGE IN EACH COLLISION, NOT ONLY ON AVERAGE
C
            IF (IESTCX(IRCX,1).NE.0) THEN
C  IATMN: ATOM SPECIES AFTER CX
              IATMN=IATM

              IF (LPMML) THEN
                PMML(IOLD,NCELL) =PMML(IOLD,NCELL)-WGHTO
                LMETSP(NSPA+IOLD)=.TRUE.
              END IF
              IF (LPMAT) THEN
                PMAT(IATMN,NCELL)=PMAT(IATMN,NCELL)+WEIGHT
                LMETSP(NSPH+IATMN)=.TRUE.
              END IF
              IF (LPMPL) THEN
                PMPL(IPLS,NCELL) =PMPL(IPLS,NCELL)-WEIGHT
                LMETSP(NSPAMI+IPLS)=.TRUE.
              END IF
              IF (LPMEL) PMEL(NCELL)      =PMEL(NCELL)-WEIGHT
              IF (N2NDX(IRCX,1).EQ.4) THEN
C  IPLSN: ION SPECIES AFTER CX
                IPLSN=N2NDX(IRCX,2)
                IF (LPMPL) THEN
                  PMPL(IPLSN,NCELL)=PMPL(IPLSN,NCELL)+WGHTO
                  LMETSP(NSPAMI+IPLSN)=.TRUE.
                END IF
                IF (LPMEL) PMEL(NCELL)      =PMEL(NCELL)+WGHTO
              ELSEIF (N2NDX(IRCX,1).NE.4) THEN
                GOTO 999
              ENDIF
            ENDIF
c  UPDATE collision estimator for CX energy exchange tallies
            IF (IESTCX(IRCX,3).NE.0) THEN
              IF (LEMML) EMML(NCELL)=EMML(NCELL)-E0O*WGHTO
              IF (LEMAT) EMAT(NCELL)=EMAT(NCELL)+E0*WEIGHT
              IF (LEMPL) THEN
                EMPL(IPLS,NCELL)=EMPL(IPLS,NCELL)-E0*WEIGHT
                LMETSP(NSPAMI+IPLS)=.TRUE.
              END IF
              IF (N2NDX(IRCX,1).EQ.4) THEN
                IF (LEMPL) THEN
                  IPLSN=N2NDX(IRCX,2)
                  EMPL(IPLSN,NCELL)=EMPL(IPLSN,NCELL)+E0O*WGHTO
                  LMETSP(NSPAMI+IPLSN)=.TRUE.
                ENDIF
              ELSE
                GOTO 999
              ENDIF
            ENDIF
C  UPDATE COLLISION ESTIMATOR CONTRIBUTION TO MMPL (FORMERLY: COPV)
            IF (IESTCX(IRCX,2).NE.0) THEN
              IF (LMMPL) THEN
C  SET THE POST COLLISION NEUTRAL PARALLEL VELOCITY = OLD PRE COLLISION ION VELOCITY
                V0_PARB=VEL*(VELX*BX+VELY*BY+VELZ*BZ)
                V0_PARB=V0_PARB*AMUA*RMASSM(IMOL)
                IPLSV=MPLSV(IPLS)
                IF (INDPRO(4) == 8) THEN
                  CALL EIRENE_VECUSR(2,NCELL,X0,Y0,Z0,VX,VY,VZ,IPLS,
     .                               .TRUE.)
                  VPLASP=VX*BX+VY*BY+VZ*BZ
                ELSE
                  VPLASP=BVIN(IPLSV,NCLLO)
                ENDIF
                SIG=SIGN(1._DP,VPLASP)
C ASSUME: OLD (INCIDENT) ION MOMENTUM IS EQUAL TO NEW MOLECULE MOMENTUM
                MMPL(IPLS,NCELL)=MMPL(IPLS,NCELL)-WEIGHT*V0_PARB*SIG
                LMETSP(NSPAMI+IPLS)=.TRUE.
              END IF

              IF (N2NDX(IRCX,1).EQ.4) THEN
                IF (LMMPL) THEN
C  IPLSN ION SPECIES AFTER CX
                  IPLSN=N2NDX(IRCX,2)
C  ASSUME: NEW ION MOMENTUM IS EQUAL TO INCIDENT MOLECULE MOMENTUM
                  MMPL(IPLSN,NCELL)=MMPL(IPLSN,NCELL)+
     .                              WGHTO*V0_PARBO*SIG
                  LMETSP(NSPAMI+IPLSN)=.TRUE.
                END IF
              ELSEIF (N2NDX(IRCX,1).NE.4) THEN
                GOTO 999
              ENDIF
            ENDIF
            COLTYP=1
            NCELL = NCLLO
            RETURN

          CASE(2)
C  1ST SECONDARY IS MOLECULE
            IMOL=N1STX(IRCX,2)
            E0=CVRSSM(IMOL)*VELQ
C
            IF (NGENM(IMOL).GT.0) THEN
              IF (IMOL.EQ.IOLD) THEN
                XGENER=XGENER+1.D0
              ELSE
                XGENER=0.D0
              ENDIF
              IF (XGENER.GE.NGENM(IMOL)) THEN
C  UPDATE GENERATION LIMIT TALLIES
                IF (LPGENM) PGENM(IMOL,NCELL)=PGENM(IMOL,NCELL)-WEIGHT
                IF (LEGENM)
     .            EGENM(IMOL,NCELL)=EGENM(IMOL,NCELL)-WEIGHT*E0
                IF (LVGENM) THEN
                  V0_PARB=VEL*(VELX*BX+VELY*BY+VELZ*BZ)
                  V0_PARB=V0_PARB*AMUA*RMASSM(IMOL)
                  VGENM(IMOL,NCELL)=VGENM(IMOL,NCELL)-WEIGHT*V0_PARB
                END IF
                IF (LPGENM.OR.LEGENM.OR.LVGENM) LMETSP(NSPA+IMOL)=.TRUE.
                LGPART=.FALSE.
                IF (NLTRC) CALL EIRENE_CHCTRC(X0,Y0,Z0,16,16)
                ITYP=4
                COLTYP=2
                NCELL = NCLLO
                RETURN
              ENDIF
            ENDIF
C
            IF (IESTCX(IRCX,1).NE.0) GOTO 999
            IF (IESTCX(IRCX,2).NE.0) GOTO 999
            IF (IESTCX(IRCX,3).NE.0) GOTO 999
            COLTYP=1
            NCELL = NCLLO
            RETURN

          CASE(3)
C  1ST SECONDARY IS TEST ION
            IION=N1STX(IRCX,2)
            E0=CVRSSI(IION)*VELQ
            XGENER=0.D0
C
            IF (IESTCX(IRCX,1).NE.0) GOTO 999
            IF (IESTCX(IRCX,2).NE.0) GOTO 999
            IF (IESTCX(IRCX,3).NE.0) GOTO 999
            COLTYP=2
            NCELL = NCLLO
            RETURN

          CASE DEFAULT
            WRITE (iunout,*) ' ITYP = ',ITYP,' AS FIRST SECONDARY IS',
     .                  ' NOT FORESEEN IN COLLIDE '
          END SELECT

        ELSE
 
C  FOLLOW 2ND SECONDARY, SPEED OF PREVIOUS TEST PARTICLE
          ITYP=N2NDX(IRCX,1)

          SELECT CASE(ITYP)
C
          CASE(1)
            IATM=N2NDX(IRCX,2)
            XGENER=0.D0
C
            E0=CVRSSA(IATM)*VELO*VELO
            IF (IESTCX(IRCX,1).NE.0) GOTO 999
            IF (IESTCX(IRCX,2).NE.0) GOTO 999
            IF (IESTCX(IRCX,3).NE.0) GOTO 999
            COLTYP=1
            NCELL = NCLLO
            RETURN
C
          CASE(2)
            IMOL=N2NDX(IRCX,2)
            XGENER=0.D0
C
            E0=CVRSSM(IMOL)*VELO*VELO
            IF (IESTCX(IRCX,1).NE.0) GOTO 999
            IF (IESTCX(IRCX,2).NE.0) GOTO 999
            IF (IESTCX(IRCX,3).NE.0) GOTO 999
            COLTYP=1
            NCELL = NCLLO
            RETURN
C
          CASE(3)
            IION=N2NDX(IRCX,2)
            XGENER=0.D0
C
            E0=CVRSSI(IION)*VELO*VELO
            IF (IESTCX(IRCX,1).NE.0) GOTO 999
            IF (IESTCX(IRCX,2).NE.0) GOTO 999
            IF (IESTCX(IRCX,3).NE.0) GOTO 999
            COLTYP=2
            NCELL = NCLLO
            RETURN

          CASE DEFAULT
            WRITE (iunout,*) ' ITYP = ',ITYP,' AS SECOND SECONDARY IS',
     .                  ' NOT FORESEEN IN COLLIDE '
          END SELECT

        ENDIF
C
C  ELASTIC COLLISION
C
      ELSEIF (ZEP1.LE.SIGEIT+SIGCXT+SIGELT) THEN
C
        IF (NLTRC) CALL EIRENE_CHCTRC(X0,Y0,Z0,16,5)
C   FIND SPECIES INDEX OF BULK (ION) COLLISION PARTNER
        SIGSUM=SIGEIT+SIGCXT
        DO 398 IMEL=1,NMELIM(IMOL)
          IREL=LGMEL(IMOL,IMEL,0)
          IPLS=LGMEL(IMOL,IMEL,1)
          SIGSUM=SIGSUM+SIGVEL(IREL)
          IF (ZEP1.LT.SIGSUM) GOTO 399
398     CONTINUE
        IREL=LGMEL(IMOL,NMELI(IMOL),0)
        IPLS=LGMEL(IMOL,NMELI(IMOL),1)
399     CONTINUE
C
C  NEW SPECIES INDEX AND ENERGY
C       WEIGHT=WEIGHT*1.
C  FOLLOW SECONDARY, NEW SPEED FROM SUBROUTINE VELOEL
C       ITYP=2
        NFLAG=CFLAG(5,IREL)
        RMMIO=RMASSM(IOLD)
        CALL EIRENE_VELOEL(NCLLO,VELXO,VELYO,VELZO,VELO,IOLD,NOLD,VELQ,
     .              NFLAG,IREL,RMMIO)
C
        IMOL=IOLD
C  NOT: WEIGHT=WGHTO, BECAUSE WEIGHT MAY HAVE CHANGED DUE TO NON-ANALOGUE SAMPLING IN VELOEL
        E0=CVRSSM(IMOL)*VELQ
C  DO NOT UPDATE BGK TALLIES HERE
        IBGK=NPBGKP(IPLS,1)
        IF (IBGK.NE.0) GOTO 400
C  UPDATE COLLISION ESTIMATOR CONTRIBUTION
C  ASSUME, AS BEFORE, NO CHANGE IN SPECIES/TYP
        IF (IESTEL(IREL,1).NE.0) THEN
          IF (LPMML) THEN
            PMML(IOLD,NCELL) =PMML(IOLD,NCELL)-WGHTO
            PMML(IMOL,NCELL) =PMML(IMOL,NCELL)+WEIGHT
            LMETSP(NSPA+IOLD)=.TRUE.
            LMETSP(NSPA+IMOL)=.TRUE.
          END IF
        ENDIF
c  UPDATE collision estimator for EL energy exchange tallies
        IF (IESTEL(IREL,3).NE.0) THEN
          EDEL=E0O*WGHTO-E0*WEIGHT
          IF (LEMML) EMML(NCELL)      =EMML(NCELL)-EDEL
          IF (LEMPL) THEN
            EMPL(IPLS,NCELL) =EMPL(IPLS,NCELL)+EDEL
            LMETSP(NSPAMI+IPLS)=.TRUE.
          END IF
        ENDIF
C  UPDATE COLLISION ESTIMATOR CONTRIBUTION TO MMPL (FORMERLY: COPV)
        IF (IESTEL(IREL,2).NE.0) THEN
          IF (LMMPL) THEN
C  SET THE POST COLLISION NEUTRAL PARALLEL VELOCITY 
            V0_PARB=VEL*(VELX*BX+VELY*BY+VELZ*BZ)
            V0_PARB=V0_PARB*AMUA*RMASSM(IMOL)
            VDEL=V0_PARBO*WGHTO-V0_PARB*WEIGHT
            IPLSV=MPLSV(IPLS)
            IF (INDPRO(4) == 8) THEN
              CALL EIRENE_VECUSR(2,NCELL,X0,Y0,Z0,VX,VY,VZ,IPLS,
     .                           .TRUE.)
              VPLASP=VX*BX+VY*BY+VZ*BZ
            ELSE
              VPLASP=BVIN(IPLSV,NCLLO)
            ENDIF
            SIG=SIGN(1._DP,VPLASP)
            MMPL(IPLS,NCELL)=MMPL(IPLS,NCELL)+VDEL*SIG
            LMETSP(NSPAMI+IPLS)=.TRUE.
          END IF
        ENDIF
400     CONTINUE
        COLTYP=1
        NCELL = NCLLO
        RETURN
C
C  GENERAL ION IMPACT COLLISION: PI-PROCESSES. NOT READY
C
      ELSEIF (ZEP1.LE.SIGEIT+SIGCXT+SIGELT+SIGPIT) THEN
C
        IF (NLTRC) CALL EIRENE_CHCTRC(X0,Y0,Z0,16,3)
        SIGSUM=SIGEIT+SIGCXT+SIGELT
        DO 461 IMPI=1,NMPIIM(IMOL)
C   FIND INDEX OF THAT ION IMPACT COLLISION
          IRPI=LGMPI(IMOL,IMPI,0)
          IPLS=LGMPI(IMOL,IMPI,1)
          SIGSUM=SIGSUM+SIGVPI(IRPI)
          IF (ZEP1.LT.SIGSUM) GOTO 462
461     CONTINUE
        IRPI=LGMPI(IMOL,NMPII(IMOL),0)
        IPLS=LGMPI(IMOL,NMPII(IMOL),1)
462     CONTINUE
C
C  CALCULATE WEIGHT OF THE NEXT GENERATION PARTICLE
C  ONLY ONE ATOM, MOLECULE OR TEST-ION HISTORY WITH MODIFIED WEIGHT
C  IS FOLLOWED
C
        PTOT=P2NPI(IRPI)
C       PTOTAL=PTOT+PPLPI(IRPI,0)
C  ABSORBED WEIGHT: WEIABS
C       WEIABS=WEIGHT*PPLPI(IRPI,0)
C
C  PRE- COLLISION ESTIMATOR FOR EMML,
C  PRE- AND POST COLLISION ESTIMATOR FOR EMPL AND EMEL
        IF (IESTPI(IRPI,3).NE.0) THEN
C  score loss of incoming test particle energy
          IF (LEMML) EMML(NCELL)=EMML(NCELL)-WEIGHT*E0

cdr EMPL, EMEL       :  SCORE NET CHANGES HERE.
cdr EMAT, EMML, EMIO :  SCORE EXACT GAINS LATER. 
          IF (LEMPL) THEN
            DO IP=1,IPPLPI(IRPI,0)
cdr:  this is incorrect. esigpi must be split into ipl secondaries
              IPL=IPPLPI(IRPI,IP)
              LOGPLS(IPL,ISTRA)=.TRUE.
              EMPL(IPL,NCELL)=EMPL(IPL,NCELL)+WEIGHT*ESIGPI(IRPI,4)
              LMETSP(NSPAMI+IPL)=.TRUE.
            END DO
          END IF
          IF (LEMEL) EMEL(NCELL)=EMEL(NCELL)+WEIGHT*ESIGPI(IRPI,5)
        ENDIF
C
C  ABSORBTION (INTO BULK SPECIES) IS SUPPRESSED
        WEIGHT=WEIGHT*PTOT
C
C  ARE THERE TEST PARTICLE SECONDARIES AT ALL?
        IF (WEIGHT.LE.EPS30) THEN
          LGPART=.FALSE.
          ITYP=4
          COLTYP=2
          NCELL=NCLLO
          RETURN
        ENDIF
C
        NFLAG=CFLAG(4,IRPI)
        RMMIO=RMASSM(IOLD)

        IF (NLCASCAD .AND. (NLEVEL+PTOT <= MAXLEVEL)) THEN  ! PI PROCESS CASCADING MOL

          IF (.NOT.ALLOCATED(NAMIPI)) THEN
            ALLOCATE(NAMIPI(NSPAMI))
          END IF
          NAMIPI = 0
          NAMIPI(NSPH+1:NSPA) = PATPI(IRPI,1:NATMI)
          NAMIPI(NSPA+1:NSPAM) = PMLPI(IRPI,1:NMOLI)
          NAMIPI(NSPAM+1:NSPAMI) = PIOPI(IRPI,1:NIONI)

!  RESET WEIGHT TO ORIGINAL VALUE
          WEIGHT=WEIGHT / PTOT

          DO I = NSPAMI, NSPH+1, -1
            DO J=1, NAMIPI(I)
              ZEP = 0.5_DP * (P2NP(IRPI,I-1)+P2NP(IRPI,I))
              CALL EIRENE_VELOPI(NCLLO,VELXO,VELYO,VELZO,VELO,IOLD,
     .                           NOLD,VELQ,NFLAG,IRPI,RMMIO,ZEP)
              ISPZ = ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
C
C.....................................................................
C  SPLITTING
C
              NLEVEL=NLEVEL+1
C  SAVE LOCATION, WEIGHT AND OTHER PARAMETERS AT CURRENT LEVEL
              RSPLST(1:NPARTC,NLEVEL)=RPST(1:NPARTC)
              ISPLST(1:MPARTC,NLEVEL)=IPST(1:MPARTC)
C  NUMBER OF NODES AT THIS LEVEL
              NODES(NLEVEL)=2

              IF (NLTRC) WRITE (IUNOUT,*) 'STORE ', TEXTS(ISPZ)
              
            END DO
          END DO

!  REMOVE LAST PARTICLE FROM STORAGE AS IT'S TRAJECTORY IS CONTINUED
          NLEVEL = NLEVEL - 1
          IF (NLTRC) WRITE(IUNOUT,*) 'REMOVE FROM STORAGE ', TEXTS(ISPZ)

        ELSE  ! NOT ENOUGH STORAGE FOR CASCADING

          IF (NLCASCAD) THEN
            WRITE (iunout,*) 
     .        'ANALOG CALCULATION ABANDONED FOR PART. NO. ',NPANU
            WRITE (iunout,*) 'CASCADE OVERFLOW: NEVEL: ',NLEVEL
          ENDIF

          CALL EIRENE_VELOPI(NCLLO,VELXO,VELYO,VELZO,VELO,IOLD,
     .                       NOLD,VELQ,NFLAG,IRPI,RMMIO,-1._DP)

        END IF

        XGENER=0.D0
C
C  UPDATE COLLISION ESTIMATORS CONTRIBUTION TO EAAT;EAML;EAIO
        IF (ITYP.EQ.1) THEN
          IF (IESTPI(IRPI,3).NE.0) THEN
            IF (LEMAT) EMAT(NCELL)=EMAT(NCELL)+WEIGHT*E0
          ENDIF
          COLTYP=1
        ELSEIF (ITYP.EQ.2) THEN
          IF (IESTPI(IRPI,3).NE.0) THEN
            IF (LEMML) EMML(NCELL)=EMML(NCELL)+WEIGHT*E0
          ENDIF
          COLTYP=1
        ELSEIF (ITYP.EQ.3) THEN
          IF (IESTPI(IRPI,3).NE.0) THEN
            IF (LEMIO) EMIO(NCELL)=EMIO(NCELL)+WEIGHT*E0
          ENDIF
          COLTYP=2
        ENDIF
        NCELL = NCLLO
        RETURN
C
C
      ELSE
C
C
        WRITE (iunout,*) 'ERROR IN COLMOL, UNKNOWN TYPE OF COLLISION '
        CALL EIRENE_EXIT_OWN(1)
C
C
      ENDIF
      GOTO 999
C
      ENTRY EIRENE_COLION(CFLAG,COLTYP,DIST)
C
C  INCIDENT SPECIES: IOLD
      VELXO=VELX
      VELYO=VELY
      VELZO=VELZ
      VELO=VEL
      NCLLO = NCELL
      NCELL = NCLTAL(NCLLO)
  
      IF (LMIPL.OR.NGENI(IION).NE.0) THEN
        CALL EIRENE_BFIELD (NCLLO, X0, Y0, Z0, BX, BY, BZ, BF,.TRUE.)
        V0_PARBO=VEL*(VELX*BX+VELY*BY+VELZ*BZ)
        V0_PARBO=V0_PARBO*AMUA*RMASSI(IION)
      ENDIF

      E0O=E0
      WGHTO=WEIGHT
      IOLD=IION
      NOLD=NSPAM+IION

      IF (IMETCL(NCELL) == 0) THEN
        NCLMT = NCLMT+1
        ICLMT(NCLMT) = NCELL
        IMETCL(NCELL) = NCLMT
      END IF
C
C  ABSORBTION BIASSING: SUPPRESS IREI PROCESSES WITH ZERO
C                       TEST PARTICLE SECONDARIES
 
      SIG_ELIM=0.
      SIG_TOT_N=SIGTOT
      SIG_TOT_O=SIGTOT
      NEII_RED=0
 
      DO IIEI=1,NIEII(IOLD)
        IREI=LGIEI(IOLD,IIEI)
        IF (WEIGHT.GT.WMINV) THEN
C  REMAINING RATE AFTER POSSIBLE ELIMINATION OF IREI
C  SIG_TEST=0 WOULD VIOLATE RADON-NYKODYM CONDITION OF WEIGHTING
          SIG_TEST=SIG_TOT_N-SIGVEI(IREI)
          PTOT=P2NDS(IREI)
          IF (PTOT.EQ.0..AND.SIG_TEST.GT.0.) THEN
C  ELIMINATE PROCESS IREI FROM ALL POSSIBLE PROCESSES
C  REDUCE WEIGHT ACCORDINGLY
            SIG_ELIM=SIG_ELIM+SIGVEI(IREI)
            SIG_TOT_N=SIG_TEST
            WEIGHT=WEIGHT*SIG_TOT_N/SIG_TOT_O
            WGHTO=WEIGHT
            SIG_TOT_O=SIG_TOT_N
            IF (IESTEI(IREI,1).NE.0) GOTO 990
            IF (IESTEI(IREI,2).NE.0) GOTO 990
            IF (IESTEI(IREI,3).NE.0) GOTO 990
          ELSE
C  NO, THIS PROCESS REMAINS ACTIVE
            NEII_RED=NEII_RED+1
            LGEI_RED(NEII_RED)=IREI
          ENDIF
        ELSE
C  WEIGHT TOO SMALL COMPARED TO WMINV. ANALOG GAME
          NEII_RED=NEII_RED+1
          LGEI_RED(NEII_RED)=IREI
        ENDIF
      ENDDO
C
C  FIRST DECIDE: ELECTRON IMPACT OR ION IMPACT
C
      ZEP1=SIG_ELIM+RANF_EIRENE( )*SIG_TOT_N
      SIGSUM=SIG_ELIM
      
C
      IF (ZEP1.LE.SIGEIT) THEN
C
C  AT THIS POINT: NEII.GE.1, FOR OTHERWISE ZEP1 COULD NOT HAVE
C                 POINTED TO EI-PROCESSES
C
C  ELECTRON IMPACT COLLISION:
C
        IF (NLTRC) CALL EIRENE_CHCTRC(X0,Y0,Z0,16,2)
        IF (NLSTOR) CALL EIRENE_STORE(2)
C
C  FIND TYP OF ELECTR. IMPACT COLLISION PROCESS: IREI
        DO 440 IIEI=1,NEII_RED-1
          IREI=LGEI_RED(IIEI)
          SIGSUM=SIGSUM+SIGVEI(IREI)
          IF (ZEP1.LE.SIGSUM) GOTO 445
440     CONTINUE
        IREI=LGEI_RED(NEII_RED)
445     CONTINUE
C
C  CALCULATE WEIGHT OF THE NEXT GENERATION PARTICLE
C  ONLY ONE ATOM, MOLECULE OR TEST-ION HISTORY WITH MODIFIED WEIGHT
C  IS FOLLOWED
C
        PTOT=P2NDS(IREI)
C       PTOTAL=PTOT+PPLEI(IREI,0)
C  ABSORBED WEIGHT: WEIABS
C       WEIABS=WEIGHT*PPLEI(IREI,0)
C
C  COLLISION ESTIMATOR FOR EIIO, EIPL AND EIEL
        IF (IESTEI(IREI,3).NE.0) THEN
C  score loss of incoming test particle energy
          IF (LEIIO) EIIO(NCELL)=EIIO(NCELL)-WEIGHT*E0

cdr EIPL, EIEL       :  SCORE NET CHANGES HERE.
cdr EIAT, EIML, EIIO :  SCORE EXACT GAINS LATER. 
          IF (LEIPL)  THEN
            DO IP=1,IPPLEI(IREI,0)
cdr:  this is incorrect. esigei must be split into ipl secondaries
              IPL=IPPLEI(IREI,IP)
              LOGPLS(IPL,ISTRA)=.TRUE.
              EIPL(IPL,NCELL)=EIPL(IPL,NCELL)+WEIGHT*ESIGEI(IREI,4)
              LMETSP(NSPAMI+IPL)=.TRUE.
            END DO
          END IF
          IF (LEIEL) EIEL(NCELL)=EIEL(NCELL)+WEIGHT*ESIGEI(IREI,5)
        ENDIF
C
C  ABSORBTION (INTO BULK SPECIES) IS SUPPRESSED
        WEIGHT=WEIGHT*PTOT
C
C  ARE THERE TEST PARTICLE SECONDARIES AT ALL?
        IF (WEIGHT.LE.EPS30) THEN
          LGPART=.FALSE.
          ITYP=4
          COLTYP=2
          NCELL = NCLLO
          RETURN
        ENDIF
C
        IF (NLCASCAD .AND. (NLEVEL+PTOT <= MAXLEVEL)) THEN ! EI PROCESS CASCADING ION

          IF (.NOT.ALLOCATED(NAMIEI)) THEN
            ALLOCATE(NAMIEI(NSPAMI))
          END IF
          NAMIEI = 0
          NAMIEI(NSPH+1:NSPA) = PATEI(IREI,1:NATMI)
          NAMIEI(NSPA+1:NSPAM) = PMLEI(IREI,1:NMOLI)
          NAMIEI(NSPAM+1:NSPAMI) = PIOEI(IREI,1:NIONI)

!  RESET WEIGHT TO ORIGINAL VALUE
          WEIGHT=WEIGHT / PTOT

          DO I = NSPAMI, NSPH+1, -1
            DO J=1, NAMIEI(I)
              ZEP = 0.5_DP * (P2ND(IREI,I-1)+P2ND(IREI,I))
              CALL EIRENE_VELOEI(NCLLO,IREI,VELXO,VELYO,VELZO,VELO,ZEP)
              ISPZ = ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
C
C.....................................................................
C  SPLITTING
C
              NLEVEL=NLEVEL+1
C  SAVE LOCATION, WEIGHT AND OTHER PARAMETERS AT CURRENT LEVEL
              RSPLST(1:NPARTC,NLEVEL)=RPST(1:NPARTC)
              ISPLST(1:MPARTC,NLEVEL)=IPST(1:MPARTC)
C  NUMBER OF NODES AT THIS LEVEL
              NODES(NLEVEL)=2
              
              IF (NLTRC) WRITE (IUNOUT,*) ' STORE ', TEXTS(ISPZ)

            END DO
          END DO

!  REMOVE LAST PARTICLE FROM STORAGE AS IT'S TRAJECTORY IS CONTINUED
          NLEVEL = NLEVEL - 1
          IF (NLTRC) WRITE(IUNOUT,*) 'REMOVE FROM STORAGE ', TEXTS(ISPZ)

        ELSE  ! NOT ENOUGH STORAGE FOR CASCADING

          IF (NLCASCAD) THEN
            WRITE (iunout,*) 
     .        'ANALOG CALCULATION ABANDONED FOR PART. NO. ',NPANU
            WRITE (iunout,*) 'CASCADE OVERFLOW: NEVEL: ',NLEVEL
          ENDIF

          CALL EIRENE_VELOEI(NCLLO,IREI,VELXO,VELYO,VELZO,VELO,-1._DP)

        END IF
C
C  UPDATE COLLISION ESTIMATORS CONTRIBUTION TO EIAT;EIML;EIIO
        IF (ITYP.EQ.1) THEN
          IF (IESTEI(IREI,3).NE.0) THEN
            IF (LEIAT) EIAT(NCELL)=EIAT(NCELL)+WEIGHT*E0
          ENDIF
          COLTYP=2
        ELSEIF (ITYP.EQ.2) THEN
          IF (IESTEI(IREI,3).NE.0) THEN
            IF (LEIML) EIML(NCELL)=EIML(NCELL)+WEIGHT*E0
          ENDIF
          COLTYP=2
        ELSEIF (ITYP.EQ.3) THEN
          IF (IESTEI(IREI,3).NE.0) THEN
            IF (LEIIO) EIIO(NCELL)=EIIO(NCELL)+WEIGHT*E0
          ENDIF
          COLTYP=1
        ENDIF
        NCELL = NCLLO
        RETURN
C
      ELSEIF (ZEP1.LE.SIGEIT+SIGCXT) THEN
C
C  CHARGE EXCHANGE:
C
        IF (NLTRC) CALL EIRENE_CHCTRC(X0,Y0,Z0,16,6)
C
C   FIND SPECIES INDEX OF CHARGE EXCHANGING BULK ION
        SIGSUM=SIGEIT
        DO 490 IICX=1,NICXIM(IION)
          IRCX=LGICX(IION,IICX,0)
          IPLS=LGICX(IION,IICX,1)
          SIGSUM=SIGSUM+SIGVCX(IRCX)
          IF (ZEP1.LT.SIGSUM) GOTO 491
490     CONTINUE
        IRCX=LGICX(IION,NICXI(IION),0)
        IPLS=LGICX(IION,NICXI(IION),1)
491     CONTINUE
C
C  ARE THERE SECONDARY TEST PARTICLES AT ALL?
        FRSTP=N1STX(IRCX,3)
        SCNDP=N2NDX(IRCX,3)
        IF (SCNDP.LE.EPS30) THEN
          LGPART=.FALSE.
          ITYP=4
          COLTYP=2
          NCELL = NCLLO
          RETURN
        ENDIF

        IF (NLCASCAD .AND. NLEVEL < MAXLEVEL) THEN  !  CX PROCESS CASCADING ION

C  STORE 2ND SECONDARY, SPEED OF PREVIOUS TEST PARTICLE
          ITYP=N2NDX(IRCX,1)
          NCELL = NCLLO
          XGENER=0.D0

          IF (ITYP /= 4) THEN
            SELECT CASE (ITYP)
C
            CASE(1)
              IATM=N2NDX(IRCX,2)
              E0=CVRSSA(IATM)*VELO*VELO
 
            CASE(2)
              IMOL=N2NDX(IRCX,2)
              E0=CVRSSM(IMOL)*VELO*VELO
 
            CASE(3)
              IION=N2NDX(IRCX,2)
              E0=CVRSSI(IION)*VELO*VELO
 
            CASE DEFAULT
              WRITE (iunout,*) ' ITYP ',ITYP,' AS 2ND SECONDARY IS NOT',
     .                    ' FORESEEN IN COLLIDE '
            END SELECT

            ISPZ = ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
C
C.....................................................................
C  SPLITTING
C
            NLEVEL=NLEVEL+1
C  SAVE LOCATION, WEIGHT AND OTHER PARAMETERS AT CURRENT LEVEL
            RSPLST(1:NPARTC,NLEVEL)=RPST(1:NPARTC)
            ISPLST(1:MPARTC,NLEVEL)=IPST(1:MPARTC)
C  NUMBER OF NODES AT THIS LEVEL
            NODES(NLEVEL)=2

            IF (NLTRC) WRITE (IUNOUT,*) ' STORE ', TEXTS(ISPZ)

          END IF

C  FOLLOW 1ST SECONDARY
          
          ZEP3 = 0.5*FRSTP

        ELSE  ! NOT ENOUGH STORAGE FOR CASCADING

          IF (NLCASCAD) THEN
            WRITE (iunout,*) 
     .        'ANALOG CALCULATION ABANDONED FOR PART. NO. ',NPANU
            WRITE (iunout,*) 'CASCADE OVERFLOW: NEVEL: ',NLEVEL
          ENDIF

          WEIGHT=WEIGHT*SCNDP
          ZEP3=RANF_EIRENE( )*SCNDP

        END IF
C
C  NEW SPECIES TYPE, INDEX AND ENERGY
C  SUPPRESSION OF ABSORBTION AT CX
C  I.E., NO RANDOM DECISION BETWEEN BULK AND TEST SECONDARIES
        IF (ZEP3.LE.FRSTP) THEN
C  FOLLOW FIRST SECONDARY, SPEED FROM BULK POPULATION
          ITYP=N1STX(IRCX,1)
          NFLAG=CFLAG(3,IRCX)
          CALL
     .    EIRENE_VELOCX(NCLLO,VELXO,VELYO,VELZO,VELO,IOLD,NOLD,VELQ,
     .                  NFLAG,IRCX,DUMT,DUMV)
 
          SELECT CASE (ITYP)
C
          CASE(1)
C  1ST SECONDARY IS ATOM
            IATM=N1STX(IRCX,2)
            E0=CVRSSA(IATM)*VELQ
            XGENER=0.D0
C
C  NEXT LINES: COLLISION ESTIMATOR FOR CHARGE EXCHANGE NO. IRCX
C  CONSERVE CHARGE IN EACH COLLISION, NOT ONLY ON AVERAGE
            IF (IESTCX(IRCX,1).NE.0) THEN
C  IATMN ATOM SPECIES AFTER CX
              IATMN=IATM
              IF (LPIIO) THEN
                PIIO(IOLD,NCELL) =PIIO(IOLD,NCELL)-WGHTO
                LMETSP(NSPAM+IOLD)=.TRUE.
              END IF
              IF (LPIAT) THEN
                PIAT(IATMN,NCELL)=PIAT(IATMN,NCELL)+WEIGHT
                LMETSP(NSPH+IATMN)=.TRUE.
              END IF
              IF (LPIPL) THEN
                PIPL(IPLS,NCELL) =PIPL(IPLS,NCELL)-WEIGHT
                LMETSP(NSPAMI+IPLS)=.TRUE.
              END IF
              IF (LPIEL) PIEL(NCELL)      =PIEL(NCELL)-WEIGHT
              IF (N2NDX(IRCX,1).EQ.4) THEN
C  IPLSN ION SPECIES AFTER CX
                IPLSN=N2NDX(IRCX,2)
                IF (LPIPL) THEN
                  PIPL(IPLSN,NCELL)=PIPL(IPLSN,NCELL)+WGHTO
                  LMETSP(NSPAMI+IPLSN)=.TRUE.
                END IF
                IF (LPIEL) PIEL(NCELL)      =PIEL(NCELL)+WGHTO
              ELSEIF (N2NDX(IRCX,1).NE.4) THEN
                GOTO 999
              ENDIF
            ENDIF
            IF (IESTCX(IRCX,3).NE.0) THEN
              IF (LEIIO) EIIO(NCELL)=EIIO(NCELL)-E0O*WGHTO
              IF (LEIAT) EIAT(NCELL)=EIAT(NCELL)+E0*WEIGHT
              IF (LEIPL) THEN
                EIPL(IPLS,NCELL)=EIPL(IPLS,NCELL)-E0*WEIGHT
                LMETSP(NSPAMI+IPLS)=.TRUE.
              END IF
              IF (N2NDX(IRCX,1).EQ.4) THEN
                IF (LEIPL) THEN
                  IPLSN=N2NDX(IRCX,2)
                  EIPL(IPLSN,NCELL)=EIPL(IPLSN,NCELL)+E0O*WGHTO
                  LMETSP(NSPAMI+IPLSN)=.TRUE.
                END IF
              ELSE
                GOTO 999
              ENDIF
            ENDIF
C  UPDATE COLLISION ESTIMATOR CONTRIBUTION TO MIPL (COPV)
            IF (IESTCX(IRCX,2).NE.0) THEN
              IF (LMIPL) THEN
                V0_PARB=VEL*(VELX*BX+VELY*BY+VELZ*BZ)
                V0_PARB=V0_PARB*AMUA*RMASSI(IION)
                IPLSV=MPLSV(IPLS)
                IF (INDPRO(4) == 8) THEN
                  CALL EIRENE_VECUSR(2,NCELL,X0,Y0,Z0,VX,VY,VZ,IPLS,
     .                               .TRUE.)
                  VPLASP=VX*BX+VY*BY+VZ*BZ
                ELSE
                  VPLASP=BVIN(IPLSV,NCLLO)
                ENDIF
                SIG=SIGN(1._DP,VPLASP)
C ASSUME: OLD (INCIDENT) ION MOMENTUM IS EQUAL TO NEW ATOM MOMENTUM
                MIPL(IPLS,NCELL)=MIPL(IPLS,NCELL)-WEIGHT*V0_PARB*SIG
                LMETSP(NSPAMI+IPLS)=.TRUE.
              END IF
              IF (N2NDX(IRCX,1).EQ.4) THEN
                IF (LMIPL) THEN
C  IPLSN = BULK ION SPECIES AFTER CX
                  IPLSN=N2NDX(IRCX,2)
C ASSUME: NEW ION MOMENTUM IS EQUAL TO INCIDENT TEST ION MOMENTUM
                  MIPL(IPLSN,NCELL)=MIPL(IPLSN,NCELL)+
     .                              WGHTO*V0_PARBO*SIG
                  LMETSP(NSPAMI+IPLSN)=.TRUE.
                END IF
              ELSEIF (N2NDX(IRCX,1).NE.4) THEN
                GOTO 999
              ENDIF
            ENDIF
            COLTYP=2
            NCELL = NCLLO
            RETURN

          CASE(2)
C  1ST SECONDARY IS MOLECULE
            IMOL=N1STX(IRCX,2)
            E0=CVRSSM(IMOL)*VELQ
            XGENER=0.D0
C
            IF (IESTCX(IRCX,1).NE.0) GOTO 999
            IF (IESTCX(IRCX,2).NE.0) GOTO 999
            IF (IESTCX(IRCX,3).NE.0) GOTO 999
            COLTYP=2
            NCELL = NCLLO
            RETURN

          CASE(3)
C  1ST SECONDARY IS TEST ION
            IION=N1STX(IRCX,2)
            E0=CVRSSI(IION)*VELQ
            XGENER=0.D0
C
            IF (IESTCX(IRCX,1).NE.0) GOTO 999
            IF (IESTCX(IRCX,2).NE.0) GOTO 999
            IF (IESTCX(IRCX,3).NE.0) GOTO 999
            COLTYP=1
            NCELL = NCLLO
            RETURN

          CASE DEFAULT
            WRITE (iunout,*) ' ITYP ',ITYP,' AS 1ST SECONDARY IS NOT ',
     .                  ' FORESEEN IN COLLIDE '
          END SELECT

        ELSE
C  FOLLOW 2ND SECONDARY, SPEED OF PREVIOUS TEST PARTICLE
          ITYP=N2NDX(IRCX,1)

          SELECT CASE (ITYP)
C
          CASE (1)
            IATM=N2NDX(IRCX,2)
            E0=CVRSSA(IATM)*VELO*VELO
            XGENER=0.D0
C
            IF (IESTCX(IRCX,1).NE.0) GOTO 999
            IF (IESTCX(IRCX,2).NE.0) GOTO 999
            IF (IESTCX(IRCX,3).NE.0) GOTO 999
            COLTYP=2
            NCELL = NCLLO
            RETURN

          CASE (2)
            IMOL=N2NDX(IRCX,2)
            E0=CVRSSM(IMOL)*VELO*VELO
            XGENER=0.D0
C
            IF (IESTCX(IRCX,1).NE.0) GOTO 999
            IF (IESTCX(IRCX,2).NE.0) GOTO 999
            IF (IESTCX(IRCX,3).NE.0) GOTO 999
            COLTYP=2
            NCELL = NCLLO
            RETURN

          CASE (3)
            IION=N2NDX(IRCX,2)
            E0=CVRSSI(IION)*VELO*VELO
            XGENER=0.D0

            IF (IESTCX(IRCX,1).NE.0) GOTO 999
            IF (IESTCX(IRCX,2).NE.0) GOTO 999
            IF (IESTCX(IRCX,3).NE.0) GOTO 999
            COLTYP=1
            NCELL = NCLLO
            RETURN

          CASE DEFAULT
            WRITE (iunout,*) ' ITYP ',ITYP,' AS 2ND SECONDARY IS NOT ',
     .                  ' FORESEEN IN COLLIDE '
          END SELECT
        ENDIF
C
cdr:  at this place to be done: elastic collisions of test ions
cdr   in particular: fokker planck (velocity space diffusion) approximation
cdr:  currently still somewhere in folion. To be moved here, 
cdr   build on analogy with other elastic collisions

C
      ELSEIF (ZEP1.LE.SIGEIT+SIGCXT+SIGPIT) THEN
C
C  GENERAL ION IMPACT COLLISION: PI-PROCESSES
C
        IF (NLTRC) CALL EIRENE_CHCTRC(X0,Y0,Z0,16,3)
        SIGSUM=SIGEIT+SIGCXT
        DO 561 IIPI=1,NIPIIM(IION)
C   FIND INDEX OF THAT ION IMPACT COLLISION
          IRPI=LGIPI(IION,IIPI,0)
          IPLS=LGIPI(IION,IIPI,1)
          SIGSUM=SIGSUM+SIGVPI(IRPI)
          IF (ZEP1.LT.SIGSUM) GOTO 562
561     CONTINUE
        IRPI=LGIPI(IION,NIPII(IION),0)
        IPLS=LGIPI(IION,NIPII(IION),1)
562     CONTINUE
C
C  CALCULATE WEIGHT OF THE NEXT GENERATION PARTICLE
C  ONLY ONE ATOM, MOLECULE OR TEST-ION HISTORY WITH MODIFIED WEIGHT
C  IS FOLLOWED
C
        PTOT=P2NPI(IRPI)
C       PTOTAL=PTOT+PPLPI(IRPI,0)
C  ABSORBED WEIGHT: WEIABS
C       WEIABS=WEIGHT*PPLPI(IRPI,0)
C
C  PRE- COLLISION ESTIMATOR FOR EIIO,
C  PRE- AND POST COLLISION ESTIMATOR FOR EIPL AND EIEL
        IF (IESTPI(IRPI,3).NE.0) THEN
C  score loss of incoming test particle energy
          IF (LEIIO) EIIO(NCELL)=EIIO(NCELL)-WEIGHT*E0

cdr EIPL, EIEL       :  SCORE NET CHANGES HERE.
cdr EIAT, EIML, EIIO :  SCORE EXACT GAINS LATER. 
          IF (LEIPL)  THEN
            DO IP=1,IPPLPI(IRPI,0)
cdr:  this is incorrect. esigpi must be split into ipl secondaries
              IPL=IPPLPI(IRPI,IP)
              LOGPLS(IPL,ISTRA)=.TRUE.
              EIPL(IPL,NCELL)=EIPL(IPL,NCELL)+WEIGHT*ESIGPI(IRPI,4)
              LMETSP(NSPAMI+IPL)=.TRUE.
            END DO
          END IF
          IF (LEIEL) EIEL(NCELL)=EIEL(NCELL)+WEIGHT*ESIGPI(IRPI,5)
        ENDIF
C
C  ABSORBTION (INTO BULK SPECIES) IS SUPPRESSED
        WEIGHT=WEIGHT*PTOT
C
C  ARE THERE TEST PARTICLE SECONDARIES AT ALL?
        IF (WEIGHT.LE.EPS30) THEN
          LGPART=.FALSE.
          ITYP=4
          COLTYP=2
          NCELL=NCLLO
          RETURN
        ENDIF
C
        NFLAG=CFLAG(4,IRPI)
        RMIIO=RMASSI(IOLD)

        IF (NLCASCAD .AND. (NLEVEL+PTOT <= MAXLEVEL)) THEN !  PI PROCESS CASCADING ION

          IF (.NOT.ALLOCATED(NAMIPI)) THEN
            ALLOCATE(NAMIPI(NSPAMI))
          END IF
          NAMIPI = 0
          NAMIPI(NSPH+1:NSPA) = PATPI(IRPI,1:NATMI)
          NAMIPI(NSPA+1:NSPAM) = PMLPI(IRPI,1:NMOLI)
          NAMIPI(NSPAM+1:NSPAMI) = PIOPI(IRPI,1:NIONI)

!  RESET WEIGHT TO ORIGINAL VALUE
          WEIGHT=WEIGHT / PTOT

          DO I = NSPAMI, NSPH+1, -1
            DO J=1, NAMIPI(I)
              ZEP = 0.5_DP * (P2NP(IRPI,I-1)+P2NP(IRPI,I))
              CALL EIRENE_VELOPI(NCLLO,VELXO,VELYO,VELZO,VELO,IOLD,
     .                           NOLD,VELQ,NFLAG,IRPI,RMIIO,ZEP)
              ISPZ = ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
C
C.....................................................................
C  SPLITTING
C
              NLEVEL=NLEVEL+1
C  SAVE LOCATION, WEIGHT AND OTHER PARAMETERS AT CURRENT LEVEL
              RSPLST(1:NPARTC,NLEVEL)=RPST(1:NPARTC)
              ISPLST(1:MPARTC,NLEVEL)=IPST(1:MPARTC)
C  NUMBER OF NODES AT THIS LEVEL
              NODES(NLEVEL)=2

              IF (NLTRC) WRITE (IUNOUT,*) 'STORE ', TEXTS(ISPZ)
              
            END DO
          END DO

!  REMOVE LAST PARTICLE FROM STORAGE AS IT'S TRAJECTORY IS CONTINUED
          NLEVEL = NLEVEL - 1
          IF (NLTRC) WRITE(IUNOUT,*) 'REMOVE FROM STORAGE ', TEXTS(ISPZ)

        ELSE  ! NOT ENOUGH STORAGE FOR CASCADING

          IF (NLCASCAD) THEN
            WRITE (iunout,*) 
     .        'ANALOG CALCULATION ABANDONED FOR PART. NO. ',NPANU
            WRITE (iunout,*) 'CASCADE OVERFLOW: NEVEL: ',NLEVEL
          ENDIF

          CALL EIRENE_VELOPI(NCLLO,VELXO,VELYO,VELZO,VELO,IOLD,
     .                       NOLD,VELQ,NFLAG,IRPI,RMIIO,-1._DP)
        END IF

        XGENER=0.D0
C
C  UPDATE COLLISION ESTIMATORS CONTRIBUTION TO EIAT;EIML;EIIO
        IF (ITYP.EQ.1) THEN
          IF (IESTPI(IRPI,3).NE.0) THEN
            IF (LEIAT) EIAT(NCELL)=EIAT(NCELL)+WEIGHT*E0
          ENDIF
          COLTYP=1
        ELSEIF (ITYP.EQ.2) THEN
          IF (IESTPI(IRPI,3).NE.0) THEN
            IF (LEIML) EIML(NCELL)=EIML(NCELL)+WEIGHT*E0
          ENDIF
          COLTYP=1
        ELSEIF (ITYP.EQ.3) THEN
          IF (IESTPI(IRPI,3).NE.0) THEN
            IF (LEIIO) EIIO(NCELL)=EIIO(NCELL)+WEIGHT*E0
          ENDIF
          COLTYP=2
        ENDIF
        NCELL = NCLLO
        RETURN
C
C  ELASTIC COLLISION
C
      ELSE
C
C
        WRITE (iunout,*) 'ERROR IN COLION '
        CALL EIRENE_EXIT_OWN(1)
C
C
      ENDIF
      GOTO 999
C
      ENTRY EIRENE_COLPHOT(CFLAG,COLTYP,DIST)
C
C  INCIDENT SPECIES: IOLD
      VELXO=VELX
      VELYO=VELY
      VELZO=VELZ
      VELO=VEL
      NCLLO = NCELL
      NCELL = NCLTAL(NCLLO)

C  parallel momentum of photon:  not ready
C     CALL EIRENE_BFIELD (NCLLO, X0, Y0, Z0, BX, BY, BZ, BF,.TRUE.)
C     V0_PARBO=VEL*(VELX*BX+VELY*BY+VELZ*BZ)
c     V0_PARBO=V0_PARBO*AMUA*RMASSA(IATM)

      E0O=E0
      WGHTO=WEIGHT
      IOLD=IPHOT
      NOLD=0+IPHOT

      IF (IMETCL(NCELL) == 0) THEN
        NCLMT = NCLMT+1
        ICLMT(NCLMT) = NCELL
        IMETCL(NCELL) = NCLMT
      END IF
C
C  FIRST DECIDE: ELECTRON IMPACT OR ION IMPACT
C
      ZEP1=RANF_EIRENE( )*SIGTOT
      SIGSUM=0.
C
      IF (ZEP1.LE.SIGEIT) THEN
C
C  ELECTRON IMPACT COLLISION:
C
        goto 999
      ELSEIF (ZEP1.LE.SIGEIT+SIGCXT) THEN
C
C  CHARGE EXCHANGE:
C
        goto 999
      ELSEIF (ZEP1.LE.SIGEIT+SIGCXT+SIGELT) THEN
C
C  ELASTIC COLLISION
C
        goto 999
      ELSEIF (ZEP1.LE.SIGEIT+SIGCXT+SIGELT+SIGOTT) THEN
C
C  PHOTON (OT) COLLISION (analog to cx in colatm)
C
        IF (NLTRC) CALL EIRENE_CHCTRC(X0,Y0,Z0,16,4)
C
C   FIND SPECIES INDEX OF BULK COLLISION PARTNER
        SIGSUM=SIGEIT+SIGCXT+SIGELT
        DO IAOT=1,PHV_NPHOTI(IPHOT)-1
          IROT=PHV_LGPHOT(IPHOT,IAOT,0)
          IPLS=PHV_LGPHOT(IPHOT,IAOT,1)
          KK  =PHV_LGPHOT(IPHOT,IAOT,3)
          UPDF=PHV_LGPHOT(IPHOT,IAOT,4)
          SIGSUM=SIGSUM+SIGVOT(IROT)
          IF (ZEP1.LT.SIGSUM) GOTO 1272
        enddo
        IROT=PHV_LGPHOT(IPHOT,PHV_NPHOTI(IPHOT),0)
        IPLS=PHV_LGPHOT(IPHOT,PHV_NPHOTI(IPHOT),1)
        KK  =PHV_LGPHOT(IPHOT,PHV_NPHOTI(IPHOT),3)
        UPDF=PHV_LGPHOT(IPHOT,PHV_NPHOTI(IPHOT),4)
1272    CONTINUE
C
C  ARE THERE SECONDARY TEST PARTICLES AT ALL?
        FRSTP=dble(PHV_N1STOTph(iphot,IROT,3))
        SCNDP=dble(PHV_N2NDOTph(iphot,IROT,3))
        SUMP=frstp+scndp
        IF (SUMP.LE.EPS30) THEN
          LGPART=.FALSE.
          ITYP=4
          COLTYP=2
          NCELL = NCLLO
          RETURN
        ENDIF

csw check type first secondary
         if(phv_n1stotph(iphot,irot,1) == 4) then
            t1=phv_n2ndotph(iphot,irot,1)
            select case(t1)
            case(0)
csw check ipl
               if(phv_n2ndotph(iphot,irot,2) == 0) then
                  lgpart=.false.
                  ityp=4
                  coltyp=2
                  ncell=ncllo
                  return
               endif
            case(4)
               lgpart=.false.
               ityp=4
               coltyp=2
               ncell=ncllo
               return
            end select
         endif
csw check type second secondary
         if(phv_n2ndotph(iphot,irot,1) == 4) then
            t1=phv_n1stotph(iphot,irot,1)
            select case(t1)
            case(0)
csw check ipl
               if(phv_n1stotph(iphot,irot,2) == 0) then
                  lgpart=.false.
                  ityp=4
                  coltyp=2
                  ncell=ncllo
                  return
               endif
            case(4)
               lgpart=.false.
               ityp=4
               coltyp=2
               ncell=ncllo
               return
            end select
         endif


C
C  NEW SPECIES TYPE, INDEX AND ENERGY

C  I.E., NO RANDOM DECISION BETWEEN BULK AND TEST SECONDARIES
        WEIGHT=WEIGHT*SUMP
        ZEP3=RANF_EIRENE( )*SUMP
        IF (ZEP3.LE.FRSTP) THEN
C  FOLLOW FIRST SECONDARY, SPEED FROM BULK POPULATION
csw no  coll.estim.
            IF (PHV_IESTOTph(iphot,IROT,1).NE.0) goto 999
            IF (PHV_IESTOTph(iphot,IROT,2).NE.0) goto 999
            IF (PHV_IESTOTph(iphot,IROT,3).NE.0) goto 999
            ITYP=PHV_N1STOTph(iphot,IROT,1)
            write (iunout,*)
     .        'ot not ready for photons. exit from collide '
            call EIRENE_exit_own(1)
c           call PH_POST_ENERGY(ncllo,kk,mode,il,
c    .           iold,0,velxo,velyo,velzo,velo,e0o,ityp)
 
          SELECT CASE(ITYP)
C
          CASE(0)
C  1ST secondary is PHOTON, E0 set by PH_POST_ENERGY
            IPHOT=PHV_N1STOTph(iphot,IROT,2)
c  implement generation limit in the future...?
            COLTYP=1
            NCELL=NCLLO
            RETURN
C
          CASE(1)
C  1ST SECONDARY IS ATOM
            IATM=PHV_N1STOTph(iphot,IROT,2)
            E0=CVRSSA(IATM)*VEL*VEL
C
C
            COLTYP=1
            NCELL=NCLLO
            RETURN

          CASE(2)
C  1ST SECONDARY IS MOLECULE
            IMOL=PHV_N1STOTph(iphot,IROT,2)
            E0=CVRSSM(IMOL)*VEL*VEL
            XGENER=0.D0
C
            COLTYP=1
            NCELL = NCLLO
            RETURN

          CASE(3)
C  1ST SECONDARY IS TEST ION
            IION=PHV_N1STOTph(iphot,IROT,2)
            E0=CVRSSI(IION)*VEL*VEL
            XGENER=0.D0
C
            COLTYP=2
            NCELL = NCLLO
            RETURN
c
          case(4)
            lgpart=.false.
            ipls=phv_n1stotph(iphot,irot,2)
            e0=cvrssp(ipls)*vel*vel

            coltyp=2
            ncell=ncllo
            return

          CASE DEFAULT
            WRITE (iunout,*) ' ITYP = ',ITYP,' AS FIRST SECONDARY IS',
     .                  ' NOT FORESEEN IN COLLIDE '
          END SELECT

        ELSE
C  FOLLOW 2ND SECONDARY, SPEED OF PREVIOUS TEST PARTICLE
csw no coll.estim.
            IF (PHV_IESTOTph(iphot,IROT,1).NE.0) GOTO 999
            IF (PHV_IESTOTph(iphot,IROT,2).NE.0) GOTO 999
            IF (PHV_IESTOTph(iphot,IROT,3).NE.0) GOTO 999
          ITYP=PHV_N2NDOTph(iphot,IROT,1)
          vel=velo
            write (iunout,*)
     .        'ot not ready for photons. exit from collide '
            call EIRENE_exit_own(1)
c new energy?
c         call PH_POST_ENERGY(ncllo,kk,mode,il,
c    .           iold,0,velxo,velyo,velzo,velo,e0o,ityp)
 
          SELECT CASE(ITYP)
C
          CASE(0)
            IPHOT=PHV_N2NDOTph(iphot,IROT,2)
            XGENER=0.D0
csw e0 set by ph_post energy
            coltyp=1
            ncell=ncllo
            return

C
            CASE(1)
              IATM=PHV_N2NDOTph(iphot,IROT,2)
              XGENER= 0.D0
C
              E0=CVRSSA(IATM)*VEL*VEL
              COLTYP=1
              NCELL = NCLLO
              RETURN
C
            CASE(2)
              IMOL=PHV_N2NDOTph(iphot,IROT,2)
              XGENER= 0.D0
C
              E0=CVRSSM(IMOL)*VEL*VEL
              COLTYP=1
              NCELL = NCLLO
              RETURN
C
            CASE(3)
              IION=PHV_N2NDOTph(iphot,IROT,2)
              XGENER=0.D0
C
              E0=CVRSSI(IION)*VEL*VEL
              COLTYP=2
              NCELL = NCLLO
              RETURN
c
            case(4)
              ipls=phv_n2ndotph(iphot,irot,2)
              lgpart=.false.
              ityp=4
              e0=cvrssp(ipls)*vel*vel
              coltyp=2
              ncell=ncllo
              return
            CASE DEFAULT
              WRITE (iunout,*) ' ITYP = ',ITYP,
     .                    ' AS SECOND SECONDARY IS',
     .                    ' NOT FORESEEN IN COLLIDE '
            END SELECT

         ENDIF

      ELSE
C     GENERAL IMPACT COLLISION: NOT READY
      ENDIF

      GOTO 999
C
990   WRITE (iunout,*) 'ERROR IN COLLIDE '
      WRITE (iunout,*) 'IREI=  ',IREI,' IS SUPPRESSED, BUT'
      WRITE (iunout,*) 'COLLISION ESTIMATOR WAS SELECTED  '
      WRITE (iunout,*)
     .  'SET WMINV = INFTY, OR USE TRACKLENGTH ESTIM. '
      CALL EIRENE_EXIT_OWN(1)
C
999   WRITE (iunout,*) 'ERROR IN COLLIDE '
      WRITE (iunout,*) 'ITYP ',ITYP,IPHOT,IATM,IMOL,IION,IPLS
      CALL EIRENE_EXIT_OWN(1)
      END
