C 27.6.05 : colphot: iadd removed
C 27.6.05 : colatm : mode, il removed (lgaot(..2),lgaot(..5)
C 27.6.05 : colphot: mode, il removed (lgphot(..2),lgphot(..5)
C 15.12.05: irds --> irei, iids --> iiei
C 16.12.05: wminv re-connected to "ei"-processes. suppress
C           reactions with zero test particle secondaries
C           now connected for: colatm, colmol, colion
C           still to be done: include other processes, and colphot
C 2.2.06:  wghtO set at suppression of absorption, for collision estimators.
C 2.2.06:  REMOVED: PH PROCESSES FOR ATOMS
C          GENERATION LIMIT FOR POST-COLLISION ATOMS FROM PHOTONS: REMOVED
C 10.3.06: bug fix: LGEI_RED(NREI) --> LGEI_RED(0:NREI)
C          (some compilers had been unhappy with this)
C 20.3.07: PI reactions revised

cdr oct 14.14 some hard-wired additional tallies ADDV removed again
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
cdr DEC. 15 :  bulk ion energy estimators: species-resolved.
cdr            not ready: esigei(4, ...), esigpi(4,...) must be species-resolved.

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
cdr        synchronize and re-activate option, not ready !!
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
cdr May 17: some spelling error corrections in comments adopted from ITER branch
c            AE: analog, --> BE: analogue, etc..
cdr Nov. 17: remove call to subr.store  (flag NLSTOR: out)
cdr          comments for further unification of colatm,colmol,colion routines
cdr          P2NDS --> P2NEI




      SUBROUTINE EIRENE_COLION(CFLAG,COLTYP,DIST)
C
C  SAMPLE FROM COLLISION KERNEL C
C
C  INPUT:  COMPRT, COMMON BLOCK, CONTAINING ACTUAL PARTICLE PARAMETERS
C          CFLAG,  FLAG FOR POST-COLLISION KINETICS
C  OUTPUT: COMPRT, MODIFIED TO POST-COLLISION PARTICLE PARAMETERS
C          COLTYP, FLAG: =1 CONTINUE IN CALLING ROUTINE
C                           (FOLNEUT OR FOLION)
C                        =2 EXIT FROM CALLING ROUTINE
C                           EITHER ABSORPTION, OR
C                           TRANSITION NEUTRAL-->ION (IF CALLED
C                           BY FOLNEUT), OR
C                           TRANSITION ION-->NEUTRAL (IF CALLED
C                           BY FOLION)
C  LGPART: TRUE,  TRAJECTORY CONTINUES, AT LEAST FOR POST-COLL. SCORING.
C  LGPART: FALSE, TRAJECTORY STOPS, NO FURTHER SCORING
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
cdr  .         ,ss,ssr  ! for consistency test only. Now de-activated
      REAL(DP) :: SIG_ELIM, SIG_TOT_N, SIG_TOT_O, SIG_TEST
      INTEGER ::
     .           IICX, IIEI, IIPI, IIEL,
c    .           IMCX, IMEI, IMPI, IMEL,
c    .           IACX, IAEI, IAPI, IAEL, IAPH,
c    .
     .           IOLD, NOLD,
     .           IRCX, IREI, IRPI, IREL, IRPH,
     .           IBGK, IP, NFLAG,
     .           IATMN, IPLSN, NCLLO, IPLSV,  I, J, IPL
      INTEGER :: NEII_RED,LGEI_RED(0:NREI)

Cdr  additional arrays for  ANALOG CASCADE and SPLITTING AT COLLISIONS.
Cdr (should be set in initialization phase, not here)
CDR  check: are the corresponding arrays PATEI,PMLEI, PIOEI real
CDR         or integer (1/2 particle possible?)
      INTEGER, ALLOCATABLE, SAVE :: NAMIEI(:),NAMIPI(:)


csw add n 2lines
cdr   INTEGER :: kk,updf,t1
cdr   real(dp):: sump
csw external
      real(dp), external :: ranf_eirene

      SAVE

C  INCIDENT SPECIES: IOLD
      VELXO=VELX
      VELYO=VELY
      VELZO=VELZ
      VELO=VEL
      NCLLO = NCELL
      NCELL = NCLTAL(NCLLO)

C  PARALLEL MOMENTUM OF TEST PARTICLE INCIDENT TO COLLISION
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
C  ABSORPTION BIASING: CURRENTLY ONLY IMPLEMENTED FOR "EI-TYPE" (ELECTRON IMPACT) PROCESSES

C  SUPPRESS THOSE IREI PROCESSES WITH ZERO
C                      TEST PARTICLE SECONDARIES

      SIG_ELIM=0.
      SIG_TOT_N=SIGTOT
      SIG_TOT_O=SIGTOT
      NEII_RED=0

      IF (WEIGHT.LT.WMINV) THEN
C  WEIGHT ALREADY TOO SMALL, NO SUPPRESSION OF ABSORPTION
        NEII_RED=NIEII(IOLD)
        LGEI_RED(:)=LGIEI(IOLD,:)
      ELSE
C  TRY TO SUPPRESS ABSORPTION. IDENTIFY POSSIBLE EI PROCESSES
C                              WITH ZERO TEST PARTICLE SECONDARIES
cdr     ss=0.
        DO IIEI=1,NIEII(IOLD)
          IREI=LGIEI(IOLD,IIEI)
cdr       ss=ss+SIGVEI(IREI)
C  WHILE BEING IN THIS LOOP WEIGHT MAY BE REPEATEDLY REDUCED, FOR EARLIER (LOWER) IIEI
          IF (WEIGHT.GT.WMINV) THEN
C  REMAINING RATE AFTER POSSIBLE ELIMINATION OF IREI
C  SIG_TEST=0 WOULD VIOLATE RADON-NYKODYM CONDITION OF WEIGHTING
          SIG_TEST=SIG_TOT_N-SIGVEI(IREI)
            PTOT=P2NEI(IREI)
          IF (PTOT.EQ.0..AND.SIG_TEST.GT.0.) THEN
C  IREI IS A PURELY ABSORBING EI PROCESS, but other EI processes exist.
C  ELIMINATE THIS PROCESS IREI FROM ALL NIEII POSSIBLE EI PROCESSES
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
C  NO, THIS PROCESS REMAINS ACTIVE, BECAUSE THERE ARE TEST PARTICLE SECONDARIES
            NEII_RED=NEII_RED+1
            LGEI_RED(NEII_RED)=IREI
          ENDIF
        ELSE
C  WEIGHT TOO SMALL COMPARED TO WMINV. ANALOGUE GAME
          NEII_RED=NEII_RED+1
          LGEI_RED(NEII_RED)=IREI
        ENDIF
      ENDDO
cdr  test: ss=sigeit ?
cdr     ssr=abs(ss-sigeit)/(ss+eps6)
cdr     if (ssr.gt.1e-5) write (iunout,*) 'colatm', ss, sigeit,iold,
cdr  .                                              naeii(iold),
cdr  .                                      sigvei(1:naeii(iold))
      ENDIF  ! SUPPRESSION OF ABSORPTION AT EI PROCESSES: DONE.

C  WEIGHT MAY HAVE BEEN REDUCED NOW, AND ALSO THE NUMBER OF ACTIVE EI PROCESSES.
C  similar weight reduction (wminv-criterion) also to be done for PI and CX
C
C
C  FIRST DECIDE: ELECTRON IMPACT (COLLISION TYPE: EI) OR OTHER PROCESS
C
      ZEP1=SIG_ELIM+RANF_EIRENE( )*SIG_TOT_N
      SIGSUM=SIG_ELIM
C
      IF (ZEP1.LE.SIGEIT) THEN
C
C  AT THIS POINT: NEII_RED.GE.1, FOR OTHERWISE ZEP1 COULD NOT HAVE
C                 POINTED TO EI-PROCESSES
C
C  ELECTRON IMPACT COLLISION:
C
        IF (NLTRC) CALL EIRENE_CHCTRC(X0,Y0,Z0,16,2)
C
C  FIND TYPE OF ELECTR. IMPACT COLLISION PROCESS: IREI
        DO 440 IIEI=1,NEII_RED-1
          IREI=LGEI_RED(IIEI)
          SIGSUM=SIGSUM+SIGVEI(IREI)
          IF (ZEP1.LE.SIGSUM) GOTO 445
  440   CONTINUE
        IREI=LGEI_RED(NEII_RED)
  445   CONTINUE
C
C  CALCULATE WEIGHT OF THE NEXT GENERATION PARTICLE
C  ONLY ONE ATOM, MOLECULE OR TEST ION HISTORY WITH MODIFIED WEIGHT
C  IS FOLLOWED.
C  PTOT IS THE (INTEGER) NUMBER OF ANALOGUE NEXT GENERATION TEST PARTICLES
C
        PTOT=P2NEI(IREI)
C       PTOTAL=PTOT+PPLEI(IREI,0)
C  ABSORBED WEIGHT: WEIABS
C       WEIABS=WEIGHT*PPLEI(IREI,0)
C
C  PRE-COLLISION ESTIMATOR FOR EIIO,
C  PRE- AND POST-COLLISION ESTIMATOR FOR EIPL AND EIEL
        IF (IESTEI(IREI,3).NE.0) THEN
C  score loss of incoming test particle energy
          IF (LEIIO) EIIO(NCELL)=EIIO(NCELL)-WEIGHT*E0

cdr EIPL, EIEL       :  SCORE NET CHANGES HERE.
cdr EIAT, EIML, EIIO :  SCORE EXACT GAINS LATER.
          IF (LEIPL) THEN
            DO IP=1,IPPLEI(IREI,0)

cdr: this is incorrect. esigei must be split into ipl secondaries
cdr  it only happens to be correct if the post-collision bulk species are all the same (=ipl),
cdr  because then esigei is the total for this species.
              IPL=IPPLEI(IREI,IP)
              LOGPLS(IPL,ISTRA)=.TRUE.
              EIPL(IPL,NCELL)=EIPL(IPL,NCELL)+WEIGHT*ESIGEI(IREI,4)
              LMETSP(NSPAMI+IPL)=.TRUE.
            END DO
          END IF
          IF (LEIEL) EIEL(NCELL)=EIEL(NCELL)+WEIGHT*ESIGEI(IREI,5)
        ENDIF
C
C  ABSORPTION (INTO BULK SPECIES) IS SUPPRESSED
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
          NCELL = NCLLO
          RETURN
        ENDIF

Cdr  PTOT=0,1,2,etc..., = integer, number of next generation test particles

CC.......................................................................
        IF (.NOT.NLCASCAD) GOTO 451  !  EI PROCESS CASCADING ION
cdr
c    splitting of post-collision particles, i.e. create a true cascade

cdr  ANALOGUE SAMPLING, I.E. SPLITTING, IN CASE OF MORE THAN ONE SECONDARY.
        IF (NLEVEL+PTOT <= MAXLEVEL) THEN   ! there is still storage for splitting

cdr
          IF (.NOT.ALLOCATED(NAMIEI)) THEN
            ALLOCATE(NAMIEI(NSPAMI))
          END IF
cdr  build one single distribution of secondary test particle species, all types, include photons
cdr  this should not be done here, but instead only once, in preproc. phase !!
cdr  this NAMIEI is the underlying discrete pdf, which led to the normalized cumulative p2nd(IREI) ?
          NAMIEI = 0

          NAMIEI(1:NSPH)         = 0    !  PPHEI(IREI,1:NPHOTI) IS NOT YET SET IN XSTEI.F
          NAMIEI(NSPH+1:NSPA) = PATEI(IREI,1:NATMI)
          NAMIEI(NSPA+1:NSPAM) = PMLEI(IREI,1:NMOLI)
          NAMIEI(NSPAM+1:NSPAMI) = PIOEI(IREI,1:NIONI)

!  RESET WEIGHT BACK TO ORIGINAL VALUE
          WEIGHT=WEIGHT / PTOT

cdr  generate secondaries, one by one, call veloei, and store them on splitting arrays

          DO I = NSPAMI, NSPH+1, -1  ! LOOP OVER ALL POTENTIAL SECONDARY SPECIES 'I'
            DO J=1, NAMIEI(I)   ! THERE ARE NAMIEI(I) COPIES OF THIS SECONDARY 'I'
C  FIND A "RANDOM NUMBER" TO ENFORCE "SAMPLING" OF THIS PARTICULAR SPECIES 'I' IN VELOEI
cdr
cdr die drei zeilen hier vor: ggfls. sehr lange do loop, meist aber nur 1 oder hoechstens 2 treffer
cdr (1 oder 2 test folgeteilchen). grund in der der naechsten zeile soll ggfls 2 mal das gleiche
cdr teilchen durch zep ausgewaehlt werden.
cdr
cdr alternative: p2nei folgeteilchen gibt es. anstatt zep zu setzen: nur loop ueber diese, deren
cdr ispz dann fest mitgeben, und in veloel nicht auswürfeln

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
              NODES(NLEVEL)=2  !  ONE PARTICLE SCORE IN EACH LEVEL

              IF (NLTRC) THEN
                WRITE (IUNOUT,*) 'SPLITTING IN COLION, EI PROCESS '
                WRITE (IUNOUT,*) 'STORE ', TEXTS(ISPZ)
              ENDIF
            END DO
          END DO  ! LOOP OVER ALL POTENTIAL SECONDARIES DONE
C  FOR ALL SECONDARIES WE HAVE CALLED VELOEI, AND STORED POST-COLLISION PARAMETERS
C  ON SPLITTING ARRAYS.

!  REMOVE LAST PARTICLE FROM STORAGE AS ITS TRAJECTORY WILL BE CONTINUED NOW
          NLEVEL = NLEVEL - 1
          IF (NLTRC) WRITE(IUNOUT,*) 'REMOVE FROM STORAGE ', TEXTS(ISPZ)

CDR:   VELOEI FOR THIS CONTINUED PARTICLE HAS ALREADY BEEN CALLED
          GOTO 450

        ELSE  ! NOT ENOUGH STORAGE FOR CASCADING

          WRITE (iunout,*)
     .      'ANALOGUE CALCULATION ABANDONED FOR PART. NO. ',NPANU
          WRITE (iunout,*) 'CASCADE OVERFLOW: NEVEL: ',NLEVEL

          GOTO 451
        ENDIF  !  DONE WITH NLCASCAD OPTION

CC................................................................................
CDR:  (NORMAL) NON-CASCADING GAME AT EI PROCESSES

  451   CALL EIRENE_VELOEI(NCLLO,IREI,VELXO,VELYO,VELZO,VELO,-1._DP)

  450   CONTINUE
        XGENER=0.D0
C
C  UPDATE POST-COLLISION ESTIMATORS CONTRIBUTION TO EIAT;EIML;EIIO
C         ACCOUNT FOR POST-COLLISION CONTRIBUTIONS
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
C  CHARGE-EXCHANGE:
C
        IF (NLTRC) CALL EIRENE_CHCTRC(X0,Y0,Z0,16,6)
C
C   FIND PROCESS IRCX AND SPECIES INDEX IPLS OF INCIDENT BULK ION
        SIGSUM=SIGEIT
        DO 490 IICX=1,NICXIM(IION)
          IRCX=LGICX(IION,IICX,0)
          IPLS=LGICX(IION,IICX,1)
          SIGSUM=SIGSUM+SIGVCX(IRCX)
          IF (ZEP1.LT.SIGSUM) GOTO 491
  490   CONTINUE
        IRCX=LGICX(IION,NICXI(IION),0)
        IPLS=LGICX(IION,NICXI(IION),1)
  491   CONTINUE
C
C  ARE THERE SECONDARY TEST PARTICLES AT ALL?
        FRSTP=N1STX(IRCX,3)
        SCNDP=N2NDX(IRCX,3)

        IPLSV=MPLSV(IPLS)
C
C  ARE THERE SECONDARY TEST PARTICLES AT ALL?

        IF (SCNDP.LE.EPS30) THEN
C  POST-COLLISION ESTIMATOR FOR PAPL,EAPL,MAPL: TO BE WRITTEN
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

        IF (NLCASCAD .AND. NLEVEL < MAXLEVEL) THEN  ! CX PROCESS CASCADING ION
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

        ELSE ! NOT ENOUGH STORAGE FOR CASCADING

          IF (NLCASCAD) THEN
            WRITE (iunout,*)
     .        'ANALOG CALCULATION ABANDONED FOR PART. NO. ',NPANU
            WRITE (iunout,*) 'CX CASCADE OVERFLOW: NEVEL: ',NLEVEL
          ENDIF

C  FROM HERE: OLD GAME, NO SPLITTING

C  SUPPRESSION OF ABSORPTION AT CX
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

          SELECT CASE (ITYP)
C
          CASE(1)

C  1ST SECONDARY IS ATOM: IATM
            IATM=N1STX(IRCX,2)
            E0=CVRSSA(IATM)*VELQ
            XGENER=0.D0
C
C  NEXT LINES: COLLISION ESTIMATOR FOR CHARGE-EXCHANGE NO. IRCX
C  CONSERVE CHARGE IN EACH COLLISION, NOT ONLY ON AVERAGE
C
            IF (IESTCX(IRCX,1).NE.0) THEN
C  IATMN: ATOM SPECIES AFTER CX
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
C  IPLSN: ION SPECIES AFTER CX
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
c  UPDATE collision estimator for CX energy exchange tallies
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
C  UPDATE COLLISION ESTIMATOR CONTRIBUTION TO MIPL (FORMERLY: COPV)
            IF (IESTCX(IRCX,2).NE.0) THEN
              IF (LMIPL) THEN
C  SET THE POST-COLLISION TEST PARTICLE PARALLEL VELOCITY = OLD PRE-COLLISION BULK (ION) VELOCITY
                V0_PARB=VEL*(VELX*BX+VELY*BY+VELZ*BZ)
                V0_PARB=V0_PARB*AMUA*RMASSI(IION)
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
C  IPLSN: BULK ION SPECIES AFTER CX
                  IPLSN=N2NDX(IRCX,2)
C  ASSUME: NEW ION MOMENTUM IS EQUAL TO INCIDENT TEST ION MOMENTUM
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
C
            IF (NGENI(IION).GT.0) THEN
              IF (IION.EQ.IOLD) THEN
                XGENER=XGENER+1.D0
            ELSE
              XGENER=0.D0
            ENDIF
            IF (XGENER.GE.NGENI(IION)) THEN
C  UPDATE GENERATION LIMIT TALLIES
                IF (LPGENI) PGENI(IION,NCELL)=PGENI(IION,NCELL)-WEIGHT
                IF (LEGENI)
     .            EGENI(IION,NCELL)=EGENI(IION,NCELL)-WEIGHT*E0
                IF (LVGENI) THEN
                  V0_PARB=VEL*(VELX*BX+VELY*BY+VELZ*BZ)
                  V0_PARB=V0_PARB*AMUA*RMASSM(IMOL)
                  VGENI(IION,NCELL)=VGENI(IION,NCELL)-WEIGHT*V0_PARB
                END IF
                IF (LPGENI.OR.LEGENI.OR.LVGENI)
     .            LMETSP(NSPAM+IION)=.TRUE.
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

          CASE DEFAULT
            WRITE (iunout,*) ' ITYP = ',ITYP,' AS FIRST SECONDARY IS',
     .                  ' NOT FORESEEN IN COLLIDE '
          END SELECT

        ELSE

C  FOLLOW 2ND SECONDARY, SPEED OF PREVIOUS TEST PARTICLE
          ITYP=N2NDX(IRCX,1)

          SELECT CASE (ITYP)
C
          CASE (1)
            IATM=N2NDX(IRCX,2)
            XGENER=0.D0
C
            E0=CVRSSA(IATM)*VELO*VELO
            IF (IESTCX(IRCX,1).NE.0) GOTO 999
            IF (IESTCX(IRCX,2).NE.0) GOTO 999
            IF (IESTCX(IRCX,3).NE.0) GOTO 999
            COLTYP=2
            NCELL = NCLLO
            RETURN
C
          CASE (2)
            IMOL=N2NDX(IRCX,2)
            XGENER=0.D0
C
            E0=CVRSSM(IMOL)*VELO*VELO
            IF (IESTCX(IRCX,1).NE.0) GOTO 999
            IF (IESTCX(IRCX,2).NE.0) GOTO 999
            IF (IESTCX(IRCX,3).NE.0) GOTO 999
            COLTYP=2
            NCELL = NCLLO
            RETURN
C
          CASE (3)
            IION=N2NDX(IRCX,2)
            XGENER=0.D0
C
            E0=CVRSSI(IION)*VELO*VELO
            IF (IESTCX(IRCX,1).NE.0) GOTO 999
            IF (IESTCX(IRCX,2).NE.0) GOTO 999
            IF (IESTCX(IRCX,3).NE.0) GOTO 999
            COLTYP=1
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
cdr:  at this place to be done: elastic collisions of test ions
cdr   in particular: Fokker-Planck (velocity space diffusion--> TAU approximation?)
cdr:  currently still somewhere in folion. To be moved here,
cdr   build on analogy with other elastic collisions
C
        IF (NLTRC) CALL EIRENE_CHCTRC(X0,Y0,Z0,16,5)
C   FIND IREL, AND SPECIES INDEX IPLS OF BULK (ION) COLLISION PARTNER
        SIGSUM=SIGEIT+SIGCXT
        DO 281 IIEL=1,NIELIM(IION)
          IREL=LGIEL(IION,IIEL,0)
          IPLS=LGIEL(IION,IIEL,1)
          SIGSUM=SIGSUM+SIGVEL(IREL)
          IF (ZEP1.LT.SIGSUM) GOTO 282
  281   CONTINUE
        IREL=LGIEL(IION,NIELI(IION),0)
        IPLS=LGIEL(IION,NIELI(IION),1)
  282   CONTINUE

        IPLSV=MPLSV(IPLS)
C
C  NEW SPECIES INDEX AND ENERGY
C       WEIGHT=WEIGHT*1.
C  FOLLOW SECONDARY, NEW SPEED FROM SUBROUTINE VELOEL
C       ITYP=1
        NFLAG=CFLAG(5,IREL)
        RMAIO=RMASSI(IOLD)
        CALL EIRENE_VELOEL(NCLLO,VELXO,VELYO,VELZO,VELO,IOLD,NOLD,VELQ,
     .              NFLAG,IREL,RMAIO)
C
        IION=IOLD
C  NOT: WEIGHT=WGHTO, BECAUSE WEIGHT MAY HAVE CHANGED DUE TO NON-ANALOGUE SAMPLING IN VELOEL
        E0=CVRSSI(IION)*VELQ


C  DO NOT UPDATE BGK TALLIES HERE
        IBGK=NPBGKP(IPLS,1)
        IF (IBGK.NE.0) GOTO 300

C  UPDATE COLLISION ESTIMATOR CONTRIBUTION
C  ASSUME, AS BEFORE, NO CHANGE IN SPECIES/TYP
        IF (IESTEL(IREL,1).NE.0) THEN
          IF (LPIIO) THEN
            PIIO(IOLD,NCELL) =PIIO(IOLD,NCELL)-WGHTO
            PIIO(IION,NCELL) =PIIO(IION,NCELL)+WEIGHT
            LMETSP(NSPAMI+IOLD)=.TRUE.
            LMETSP(NSPAMI+IION)=.TRUE.
          END IF
        ENDIF
c  UPDATE collision estimator for EL energy exchange tallies
        IF (IESTEL(IREL,3).NE.0) THEN
          EDEL=E0O*WGHTO-E0*WEIGHT
          IF (LEIIO) EIIO(NCELL)      =EIIO(NCELL)-EDEL
          IF (LEIPL) THEN
            EIPL(IPLS,NCELL) =EIPL(IPLS,NCELL)+EDEL
            LMETSP(NSPAMI+IPLS)=.TRUE.
          END IF
        ENDIF
C  UPDATE COLLISION ESTIMATOR CONTRIBUTION TO MIPL (FORMERLY: COPV)
        IF (IESTEL(IREL,2).NE.0) THEN
          IF (LMIPL) THEN
C  SET THE POST-COLLISION TEST PARTICLE PARALLEL VELOCITY
            V0_PARB=VEL*(VELX*BX+VELY*BY+VELZ*BZ)
            V0_PARB=V0_PARB*AMUA*RMASSI(IION)
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
            MIPL(IPLS,NCELL)=MIPL(IPLS,NCELL)+VDEL*SIG
            LMETSP(NSPAMI+IPLS)=.TRUE.
          END IF
        ENDIF
  300   CONTINUE
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
        DO 561 IIPI=1,NIPIIM(IION)
C   FIND INDEX OF THAT ION IMPACT COLLISION
          IRPI=LGIPI(IION,IIPI,0)
          IPLS=LGIPI(IION,IIPI,1)
          SIGSUM=SIGSUM+SIGVPI(IRPI)
          IF (ZEP1.LT.SIGSUM) GOTO 562
  561   CONTINUE
        IRPI=LGIPI(IION,NIPII(IION),0)
        IPLS=LGIPI(IION,NIPII(IION),1)
  562   CONTINUE
C
C  CALCULATE WEIGHT OF THE NEXT GENERATION PARTICLE
C  ONLY ONE ATOM, MOLECULE OR TEST ION HISTORY WITH MODIFIED WEIGHT
C  IS FOLLOWED
C
        PTOT=P2NPI(IRPI)
C       PTOTAL=PTOT+PPLPI(IRPI,0)
C  ABSORBED WEIGHT: WEIABS
C       WEIABS=WEIGHT*PPLPI(IRPI,0)
C
C  PRE- COLLISION ESTIMATOR FOR EIIO,
C  PRE- AND POST-COLLISION ESTIMATOR FOR EIPL AND EIEL
        IF (IESTPI(IRPI,3).NE.0) THEN
C  score loss of incoming test particle energy
          IF (LEIIO) EIIO(NCELL)=EIIO(NCELL)-WEIGHT*E0

cdr EIPL, EIEL       :  SCORE NET CHANGES HERE.
cdr EIAT, EIML, EIIO :  SCORE EXACT GAINS LATER.
          IF (LEIPL) THEN
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
C  ABSORPTION (INTO BULK SPECIES) IS SUPPRESSED
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

Cdr  PTOT=0,1,2,etc..., = integer,  number of next generation particles

        IF (NLCASCAD .AND. (NLEVEL+PTOT <= MAXLEVEL)) THEN  ! PI PROCESS CASCADING ION

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

!  REMOVE LAST PARTICLE FROM STORAGE AS ITS TRAJECTORY IS CONTINUED
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
C
      ELSE
C
C
        WRITE (iunout,*) 'ERROR IN COLION UNKNOWN TYPE OF COLLISION '
        CALL EIRENE_EXIT_OWN(1)
C
C
      ENDIF

      GOTO 999
C

C
  990 WRITE (iunout,*) 'ERROR IN COLLIDE '
      WRITE (iunout,*) 'IREI=  ',IREI,' IS SUPPRESSED, BUT'
      WRITE (iunout,*) 'COLLISION ESTIMATOR WAS SELECTED  '
      WRITE (iunout,*)
     .  'SET WMINV = INFTY, OR USE TRACKLENGTH ESTIM. '
      CALL EIRENE_EXIT_OWN(1)
C
  999 WRITE (iunout,*) 'ERROR IN COLLIDE '
      WRITE (iunout,*) 'ITYP ',ITYP,IPHOT,IATM,IMOL,IION,IPLS
      CALL EIRENE_EXIT_OWN(1)
      END
