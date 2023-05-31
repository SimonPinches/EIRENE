C 16.12.05: wminv re-connected to EI processes. suppress
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
cdr 20.10.15: arguments in chctrc: type of collision process: corrected for PI and PH
cdr 24.11.15:  bug fix re coll est for PI processes, in colion: eiml --> eiio
cdr Dec.15  :  bug fix PI reaction and cascading was wrong:
cdr            irei, rather than irpi, and p2nd
cdr            rather than p2np, were used also for PI reactions. now corrected

cdr         :  further: collision estimators for PI processes, e§pl and e§el tallies: activated
cdr         :  see also corresponding corrections/changes in update for tracklength estimators
cdr DEC. 15 :  bulk ion energy estimators: species-resolved.
cdr            not ready: esigei(4, ...), esigpi(4,...) must be species-resolved.

cdr            tbd:  check setting of iestm..flags for collision estimators.
cdr                  probably not correct (outdated).

cdr Aug 16:    bug fix: IPPLEI --> IPPLPI at one instance
cdr Nov 16:
cdr analog cascading NLCASCAD: started to document,
cdr        synchronize and re-activate option, not ready !!
c   this version: prepare cascading at collisions,
c   e.g. for antithetic variate sampling to reduce stochastic cancellation
c   start to clean up splitting, for analogue game and for anti-correlated momentum estimators
c   started for colatm, and EI processes.
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
C            wminv activated in colmol for EI processes (analog to colatm)
cdr May 17: some spelling error corrections in comments adopted from ITER branch
c            AE: analog, --> BE: analogue, etc..
cdr Nov. 17: remove call to subr.store  (flag NLSTOR: out)
cdr          comments for further unification of colatm,colmol,colion routines
cdr          P2NDS --> P2NEI
c  unify: iold: nxeii, lgxei, rmassx, LEXEL, EXEL, EXPL, etc.
cdr 2021   : former routines colatm, colmol, colion unified into collide.f
cdr 2022   : cascading at collisions enabled: currently for CX 1st secondary
cdr 2023   : absorption biassing generalized, no also for CX.
cdr          tbd: for PI, EL, processes

      MODULE EIRMOD_COLLIDE


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
      USE EIRMOD_RANF, ONLY: RANF_EIRENE
      USE EIRMOD_VELOPI, ONLY: EIRENE_VELOPI
      USE EIRMOD_VELOEL, ONLY: EIRENE_VELOEL
      USE EIRMOD_VELOCX, ONLY: EIRENE_VELOCX
      USE EIRMOD_PLT2D, ONLY: EIRENE_CHCTRC

      IMPLICIT NONE

      PRIVATE
      
      PUBLIC :: EIRENE_COLLIDE

      REAL(DP) :: DUMT(3), DUMV(3)
      REAL(DP) :: SIGSUM, WGHTO, FRSTP, SCNDP, PTOT,
     .          VELXO, VELYO, VELZO, VELO, E0O,
     .          BX, BY, BZ, V0_PARBO,
     .          BXN(0:2), BYN(0:2), BZN(0:2),
     .          V0_O(0:2), M0_O(0:2), VP_O(0:2), MP_O(0:2),
     .          EDEL, VDEL, SIGNUM,
     .          V0_PARB,
     .          V0_N(0:2), M0_N(0:2), VP_N(0:2), MP_N(0:2),
     .          FP, FLTEST, ZEP3, VELQ,
     .          VX, VY, VZ, VPLASP, RMXIO, BF, ZEP
      REAL(DP) :: SIG_ELIM, SIG_TOT_N, SIG_TOT_O, SIG_TEST
      INTEGER ::
     .           IXCX, IXEI, IXPI, IXEL, !IXPH,
     .           IOLD, NOLD, INEW, NNEW, ITYPO, ITYPN, IPLSO,
     .           IRCX, IREI, IRPI, IREL, !IRPH,
     .           IBGK, IP, NFLAG,
     .           IATMN, IMOLN, IIONN, IPLSN, 
     .           NCLLO, IPLSV, IPL,
     .           IPTYPO, IPTYPN

Cdr  additional arrays for ANALOG CASCADE and SPLITTING AT COLLISIONS.
Cdr (should be set in initialization phase, not here)
CDR  check: are the corresponding arrays PATEI,PMLEI, PIOEI real
CDR         or integer (1/2 particle possible?)
      INTEGER, ALLOCATABLE :: NAMIEI(:),NAMIPI(:)
      CHARACTER(50) :: CCOLEST

      SAVE  

cym IAPH, RMMIO, RMIIO and IRPH removed during merge
!$OMP THREADPRIVATE (dumt,dumv,
!$OMP& SIGSUM, WGHTO, FRSTP, SCNDP, PTOT, VELXO,
!$OMP& VELYO, VELZO, VELO, E0O, BX, BY, BZ, V0_PARBO,
!$OMP& BXN, BYN, BZN, V0_O, M0_O,
!$OMP& EDEL, VDEL, SIGNUM, V0_PARB, V0_N, M0_N, VP_N, MP_N,
!$OMP& FP, FLTEST, ZEP3, VELQ,
!$OMP& VX, VY, VZ, VPLASP, RMXIO, BF, ZEP,
!$OMP& SIG_ELIM, SIG_TOT_N, SIG_TOT_O, SIG_TEST,
!$OMP& IXCX,IXEI,IXPI,IXEL,IOLD,NOLD,INEW,NNEW,ITYPO,ITYPN,IPLSO,
!$OMP& IRCX,IREI,IRPI,IREL,
!$OMP& IBGK, IP, NFLAG, IATMN, IMOLN, IIONN, IPLSN, NCLLO, IPLSV, IPL,
!$OMP& NAMIEI,NAMIPI,IPTYPO,IPTYPN,CCOLEST)

      contains

      SUBROUTINE EIRENE_COLLIDE(CFLAG,COLTYP,DIST,KK)
C
C  SAMPLE FROM COLLISION KERNEL C
C
C  INPUT:  COMPRT, COMMON BLOCK, CONTAINING ACTUAL PARTICLE PARAMETERS
C          CFLAG,  FLAG FOR POST-COLLISION KINETICS
C  OUTPUT: COMPRT, MODIFIED TO POST-COLLISION PARTICLE PARAMETERS
C          COLTYP, FLAG: =1 CONTINUE IN CALLING ROUTINE
C                           (EITHER FOLNEUT OR FOLION)
C                        =2 EXIT FROM CALLING ROUTINE
C                           EITHER ABSORPTION, OR
C                           TRANSITION NEUTRAL-->ION (IF CALLED
C                           BY FOLNEUT), OR
C                           TRANSITION ION-->NEUTRAL (IF CALLED
C                           BY FOLION)
C          KK    , GLOBAL REACTION NUMBER
C  LGPART: TRUE,  TRAJECTORY CONTINUES, AT LEAST FOR POST-COLL. SCORING.
C  LGPART: FALSE, TRAJECTORY STOPS, NO FURTHER SCORING
C
      REAL(DP), INTENT(IN) :: CFLAG(7,MSTOR0),DIST
      INTEGER, INTENT(OUT) :: COLTYP
      INTEGER, INTENT(OUT) :: KK
      INTEGER :: NEII_RED,LGEI_RED(0:NREI)
      REAL(DP) :: ZEP1
      INTEGER :: I,J
      REAL(DP), POINTER :: PXX2(:,:), PXPL2(:,:), PXAT2(:,:),
     .                     PXML2(:,:), PXIO2(:,:)

C  INIT
      KK = 0

C  INCIDENT SPECIES: IOLD
      VELXO=VELX
      VELYO=VELY
      VELZO=VELZ
      VELO=VEL
      NCLLO = NCELL
      NCELL = NCLTAL(NCLLO)

C  VELOCITY AND MOMENTUM OF TEST PARTICLE INCIDENT TO COLLISION,
C  RELATIVE TO B FIELD
cdr  needed only in case if(any(iestab(:,2) .ne. 0)), "ab" stands for ei,cx,pi,el

      IF (LMXPL.OR.NGENX.NE.0) THEN
        CALL EIRENE_BFIELD (NCLLO, X0, Y0, Z0, BX, BY, BZ, BF,.TRUE.)
        V0_PARBO=VEL*(VELX*BX+VELY*BY+VELZ*BZ)
        V0_PARBO=V0_PARBO*AMUA*RMASSX
      ENDIF

      E0O=E0
      WGHTO=WEIGHT

      IOLD=IXSPZ
      NOLD=NMETOFF+IXSPZ
      ITYPO=ITYP

      IF (ITYPO.EQ.3) THEN
cdr  called from folion 
        IPTYPO=1
      ELSE
cdr  called from folneut
        IPTYPO=0
      ENDIF      

      IF (IMETCL(NCELL) == 0) THEN
        NCLMT = NCLMT+1
        ICLMT(NCLMT) = NCELL
        IMETCL(NCELL) = NCLMT
      END IF
C
C  ABSORPTION BIASING: CURRENTLY ONLY IMPLEMENTED FOR "EI-TYPE" (ELECTRON IMPACT) PROCESSES

C  SUPPRESS ALL IREI PROCESSES WITH ZERO
C                    TEST PARTICLE SECONDARIES

      SIG_ELIM=0.
      SIG_TOT_N=SIGTOT
      SIG_TOT_O=SIGTOT
      NEII_RED=0

      IF (WEIGHT.LT.WMINV) THEN
C  WEIGHT ALREADY TOO SMALL, NO SUPPRESSION OF ABSORPTION
        NEII_RED=NXEII
        LGEI_RED(:)=LGXEI(:)
      ELSE
C  TRY TO SUPPRESS ABSORPTION. IDENTIFY POSSIBLE EI PROCESSES
C                              WITH ZERO TEST PARTICLE SECONDARIES
        DO IXEI=1,NXEII
          IREI=LGXEI(IXEI)
C  WHILE BEING IN THIS LOOP WEIGHT MAY BE REPEATEDLY REDUCED, FOR EARLIER (LOWER) IXEI
          IF (WEIGHT.GT.WMINV) THEN
C  REMAINING RATE AFTER POSSIBLE ELIMINATION OF IREI
C  SIG_TEST=0 WOULD VIOLATE RADON-NYKODYM CONDITION OF WEIGHTING
            SIG_TEST=SIG_TOT_N-SIGVEI(IREI)
            PTOT=P2NEI(IREI)
            IF (PTOT.EQ.0..AND.SIG_TEST.GT.0.) THEN
C  IREI IS A PURELY ABSORBING EI PROCESS, but other EI processes exist.
C  ELIMINATE THIS PROCESS IREI FROM ALL NXEII POSSIBLE EI PROCESSES
C  REDUCE WEIGHT ACCORDINGLY
              SIG_ELIM=SIG_ELIM+SIGVEI(IREI)
              SIG_TOT_N=SIG_TEST
              WEIGHT=WEIGHT*SIG_TOT_N/SIG_TOT_O
              WGHTO=WEIGHT
              SIG_TOT_O=SIG_TOT_N
              IF (IESTEI(IREI,1).NE.0) GOTO 997
              IF (IESTEI(IREI,2).NE.0) GOTO 997
              IF (IESTEI(IREI,3).NE.0) GOTO 997
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
C                 POINTED TO EI PROCESSES
C
C  ELECTRON IMPACT COLLISION:
C
        IF (NLTRC) THEN
!$OMP CRITICAL
          CALL EIRENE_CHCTRC(X0,Y0,Z0,16,2)
!$OMP END CRITICAL
        ENDIF
C
C  FIND TYPE OF ELECTR. IMPACT COLLISION PROCESS: IREI
        DO 240 IXEI=1,NEII_RED-1
          IREI=LGEI_RED(IXEI)
          SIGSUM=SIGSUM+SIGVEI(IREI)
          IF (ZEP1.LE.SIGSUM) GOTO 245
  240   CONTINUE
        IREI=LGEI_RED(NEII_RED)
  245   CONTINUE
C       GET GLOBAL REACTION NUMBER          
        KK = NREAEI(IREI)
C
C  CALCULATE WEIGHT OF THE NEXT GENERATION PARTICLE FOR PROCESS IREI
C  ONLY ONE ATOM, MOLECULE OR TEST ION HISTORY WITH MODIFIED WEIGHT
C  IS FOLLOWED.
C  PTOT IS THE (INTEGER) NUMBER OF ANALOGUE NEXT GENERATION TEST PARTICLES
C
        PTOT=P2NEI(IREI)
C       PTOTAL=PTOT+PPLEI(IREI,0)
C  ABSORBED WEIGHT: WEIABS
C       WEIABS=WEIGHT*PPLEI(IREI,0)
C
C  PRE-COLLISION ESTIMATOR FOR EXX,
C  NET PRE- AND POST-COLLISION ESTIMATOR FOR EXPL AND EXEL
        IF (IESTEI(IREI,3).NE.0) THEN
C  score loss of incoming test particle energy
          IF (LEXX) THEN
!$OMP ATOMIC
             EXX(NCELL)=EXX(NCELL)-WEIGHT*E0
          ENDIF

cdr EXPL, EXEL       :  SCORE NET CHANGES HERE.
cdr EXAT, EXML, EXIO :  SCORE EXACT POST-COLLISION GAINS LATER.
          IF (LEXPL) THEN
            DO IP=1,IPPLEI(IREI,0)
cdr: This is incorrect. esigei must be split into ipl secondaries.
cdr  It only happens to be correct if the post-collision bulk species are all the same (=ipl),
cdr  because then esigei is the total for this species.
cdr  For atomic test particles this error should not matter.
              IPL=IPPLEI(IREI,IP)
              LOGPLS(IPL,ISTRA)=.TRUE.
!$OMP ATOMIC
              EXPL(IPL,NCELL)=EXPL(IPL,NCELL)+WEIGHT*ESIGEI(IREI,4)
              LMETSP(NSPAMI+IPL)=.TRUE.
            END DO
          END IF
          IF (LEXEL) THEN
!$OMP ATOMIC
            EXEL(NCELL)=EXEL(NCELL)+WEIGHT*ESIGEI(IREI,5)
          ENDIF
        ENDIF

        IF (IESTEI(IREI,1).NE.0) THEN
          CCOLEST='PRE COL. PARTICLE RATE, EI PROCESS'
          GOTO 998
        ENDIF
        IF (IESTEI(IREI,2).NE.0) THEN
          CCOLEST='PRE COL. MOMENTUM RATE, EI PROCESS'
          GOTO 998
        ENDIF
C
C  ABSORPTION (INTO BULK SPECIES) IS SUPPRESSED
C  STRICTLY: A WMINV CRITERION MAY BE USED HERE FOR THIS PROCESS IREI AGAIN
C            WHEN THIS PROCESS RESULTS IN BOTH: TEST AND BULK SECONDARIES
C            ABOVE: ONLY PROCESSES IREI WITH ZERO TEST SECONDARIES MIGHT HAVE
C            BEEN SUPRESSED.
C
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

Cdr  PTOT=0,1,2,etc..., = integer, number of next generation test particles

CC.......................................................................
        IF (.NOT.NLCASCAD) GOTO 251  !  EI PROCESS CASCADING
cdr
c    splitting of post-collision particles, i.e. create a true cascade

cdr  ANALOGUE SAMPLING, I.E. SPLITTING, IN CASE OF MORE THAN ONE SECONDARY.
        IF (NLEVEL+PTOT <= MAXLEV) THEN   ! there is still storage for splitting

cdr
          IF (.NOT.ALLOCATED(NAMIEI)) THEN
            ALLOCATE(NAMIEI(NSPAMI))
          END IF
cdr  build one single distribution of secondary test particle species, all types, include photons
cdr  this should not be done here, but instead only once, in preproc. phase !!
cdr  this NAMIEI is the underlying discrete pdf, which led to the normalized cumulative p2nd(IREI)
cdr  NAMIEI is not normalized. The entries are the number of secondaries,
          NAMIEI = 0

          NAMIEI(1:NSPH)         = 0    !  PPHEI(IREI,1:NPHOTI) IS NOT YET SET IN XSTEI.F
          NAMIEI(NSPH+1:NSPA)    = INT(PATEI(IREI,1:NATMI))
          NAMIEI(NSPA+1:NSPAM)   = INT(PMLEI(IREI,1:NMOLI))
          NAMIEI(NSPAM+1:NSPAMI) = INT(PIOEI(IREI,1:NIONI))

!  RESET WEIGHT BACK TO ORIGINAL VALUE
          WEIGHT=WEIGHT / PTOT

cdr  generate secondaries, one by one, call veloei, and store them on splitting arrays

          DO I = NSPAMI, NSPH+1, -1  ! LOOP OVER ALL POTENTIAL SECONDARY SPECIES 'I'
            DO J=1, NAMIEI(I)   ! THERE ARE NAMIEI(I) COPIES OF THIS SECONDARY 'I'
C  FIND A "RANDOM NUMBER" TO ENFORCE "SAMPLING" OF THIS PARTICULAR SPECIES 'I' IN VELOEI
cdr
cdr WIP: unclear code here. Still not unravelled.
cdr die drei zeilen hier vor: ggfls. sehr lange do loop, meist aber nur 1 oder hoechstens 2 treffer
cdr (1 oder 2 test folgeteilchen). Grund in der naechsten zeile soll ggfls 2 mal das gleiche
cdr teilchen durch zep ausgewaehlt werden.
cdr
cdr alternative: p2nei folgeteilchen gibt es. anstatt zep zu setzen: nur loop ueber diese, deren
cdr ispz dann fest mitgeben, und in veloel nicht mehr auswürfeln

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
                WRITE (IUNOUT,*) 'SPLITTING IN COLLIDE, EI PROCESS '
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
          GOTO 250

        ELSE  ! NOT ENOUGH STORAGE FOR CASCADING

          WRITE (iunout,*)
     .      'ANALOGUE CALCULATION ABANDONED FOR PART. NO. ',NPANU
          WRITE (iunout,*) 'CASCADE OVERFLOW: NEVEL: ',NLEVEL

          GOTO 251
        ENDIF  !  DONE WITH NLCASCAD OPTION

CC................................................................................
CDR:  (NORMAL) NON-CASCADING GAME AT EI PROCESSES

  251   CALL EIRENE_VELOEI(NCLLO,IREI,VELXO,VELYO,VELZO,VELO,-1._DP)

  250   CONTINUE
        XGENER=0.D0

        ITYPN=ITYP
        IF (ITYP.EQ.3) THEN
cdr  return to folion 
          IPTYPN=1
        ELSE
cdr  return to folneut
          IPTYPN=0
        ENDIF

        if (iptypo .eq. iptypn) then
          coltyp=1
        else
          coltyp=2
        endif

C
C  UPDATE POST-COLLISION ESTIMATORS CONTRIBUTION TO EXAT;EXML;EXIO
C         ACCOUNT FOR POST-COLLISION CONTRIBUTIONS
        IF (IESTEI(IREI,3).NE.0) THEN
          IF (ITYP.EQ.1) THEN
            IF (LEXAT) THEN
!$OMP ATOMIC
              EXAT(NCELL)=EXAT(NCELL)+WEIGHT*E0
            ENDIF
          ELSEIF (ITYP.EQ.2) THEN
            IF (LEXML) THEN
!$OMP ATOMIC
              EXML(NCELL)=EXML(NCELL)+WEIGHT*E0
            ENDIF
          ELSEIF (ITYP.EQ.3) THEN
            IF (LEXIO) THEN
!$OMP ATOMIC
              EXIO(NCELL)=EXIO(NCELL)+WEIGHT*E0
            ENDIF
          ENDIF
        ENDIF

C  UPDATE POST-COLLISION ESTIMATORS CONTRIBUTION TO MXPL_VEC
C         ACCOUNT FOR POST-COLLISION CONTRIBUTIONS
        IF (IESTEI(IREI,2).NE.0) THEN
          CCOLEST='POST COL. MOMENTUM RATE, EI PROCESS'
          GOTO 998
        ENDIF

        NCELL = NCLLO
        RETURN
C
      ELSEIF (ZEP1.LE.SIGEIT+SIGCXT) THEN
C
C  CHARGE-EXCHANGE:
C
        IF (NLTRC) THEN
!$OMP CRITICAL              
          CALL EIRENE_CHCTRC(X0,Y0,Z0,16,6)
!$OMP END CRITICAL
        ENDIF
C
C   FIND CX PROCESS IRCX AND SPECIES INDEX IPLS OF INCIDENT BULK ION
        SIGSUM=SIGEIT
        DO 271 IXCX=1,NXCXIM
          IRCX=LGXCX(IXCX,0)
          IPLS=LGXCX(IXCX,1)
          SIGSUM=SIGSUM+SIGVCX(IRCX)
          IF (ZEP1.LT.SIGSUM) GOTO 272
  271   CONTINUE
        IRCX=LGXCX(NXCXI,0)
        IPLS=LGXCX(NXCXI,1)
  272   CONTINUE
C       GET GLOBAL REACTION NUMBER          
        KK = NREACX(IRCX)

        IPLSO=IPLS
        IPLSV=MPLSV(IPLS)
C
c  secondary test particles for the selected CX reaction IRCX
        FRSTP=N1STX(IRCX,3)
        SCNDP=N2NDX(IRCX,3)

C
C  ARE THERE SECONDARY TEST PARTICLES AT ALL?

        IF (SCNDP.LE.EPS30) THEN
C  COLLISION ESTIMATOR FOR PXPL,EXPL,MXPL_VEC: TO BE WRITTEN
C  E.G. FOR CX RECOMBINATION
          LGPART=.FALSE.
          IF (IESTCX(IRCX,1).NE.0) GOTO 999
          IF (IESTCX(IRCX,2).NE.0) GOTO 999
          IF (IESTCX(IRCX,3).NE.0) GOTO 999
          ITYP=4
          COLTYP=2
          NCELL=NCLLO
          RETURN
        ENDIF

        IF (NLCASCAD .AND. NLEVEL < MAXLEV) THEN  ! CX PROCESS CASCADING
! JUST OPPOSITE TO EI CASE:
CDR IN EI CASE: LAST TEST SECONDARY WAS FOLLOWED, ALL OTHERS STORED ON SPLITTING ARRAY.
CDR IN CX CASE: OPPOSITE.   TRY TO UNIFY !!
C NORMALLY THERE IS ONLY ONE SECONDARY, AND WE FOLLOW (THIS ONLY) TEST SECONDARY
C HERE FIRST AND LAST HAVE A SPECIAL MEANING (EXCHANGE OF IDENTITY, VELOCITIES, ETC..)

C  STORE 2ND SECONDARY, SPEED OF PREVIOUS TEST PARTICLE,
c  (i.e. scattering angle = PI), energy may have changed.

          ITYP=N2NDX(IRCX,1)
          NCELL = NCLLO
          XGENER=0.D0

          IF (ITYP /= 4) THEN

cdr 2nd secondary is a test particle
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

          ITYPN=ITYP
          IF (ITYPN.EQ.3) THEN
cdr  return to folion 
            IPTYPN=1
          ELSE
cdr  return to folneut
            IPTYPN=0
          ENDIF

          NFLAG=NINT(CFLAG(3,IRCX))
cdr  uses incident test particle velocity, reaction index IRCX, IPLSO
cdr  and returns a sampled field ion velocity as new velx,vely velz, and weight
          CALL EIRENE_VELOCX
     .         (NCLLO,VELXO,VELYO,VELZO,VELO,IOLD,NOLD,VELQ,
     .          NFLAG,IRCX,DUMT,DUMV)

          SELECT CASE(ITYP)
C
          CASE(1)

C  1ST SECONDARY IS ATOM: IATM, WEIGHT
            IATM=N1STX(IRCX,2)
            IATMN=IATM
            NNEW=NSPH+IATM
            E0=CVRSSA(IATM)*VELQ

C  GENERATION LIMIT, CX, AND SAME SPECIES
            IF (NGENX.GT.0) THEN
              IF (NNEW.EQ.NOLD) THEN
                XGENER=XGENER+1.D0
              ELSE
                XGENER=0.D0
              ENDIF
              IF (XGENER.GE.NGENX) THEN
                CALL EIRENE_GENLIM
                RETURN
              ENDIF
            ENDIF
C
C  FLUID LIMIT, CX, AND SAME SPECIES
            IF (NGENX.LT.0) THEN
              IF (NNEW.EQ.NOLD) THEN
                FP=VELO/SIGVCX(IRCX)  !mfp
                FLTEST=FP/DIST
                IF (FLTEST.LT.FDLMCX(IRCX)) THEN
                  CALL EIRENE_GENLIM
                  RETURN
                ENDIF
              ENDIF
            ENDIF
C
C  NEXT LINES: COLLISION ESTIMATOR FOR CHARGE-EXCHANGE NO. IRCX
C  CONSERVE CHARGE IN EACH COLLISION, NOT ONLY ON AVERAGE
C  "X TO AT AND PL"
C
            IF (IESTCX(IRCX,1).NE.0) THEN
C  pre-collision estimator, iold, nold, iplso
              IF (LPXX) THEN
!$OMP ATOMIC
                PXX(IOLD,NCELL) =PXX(IOLD,NCELL)-WGHTO
                LMETSP(NOLD)=.TRUE.
                IF (LSCX) THEN
                  PXX2(1:NDXX,0:NDXX) => PXX(:,NCELL)
!$OMP ATOMIC
                  PXX2(IOLD,IOLD)=PXX2(IOLD,IOLD)-WGHTO
                  LMETSP2(1:NDXX,0:NDXX) => LMETSP(NDXXA:NDXXE)
                  LMETSP2(IOLD,0) = .TRUE.
                  LMETSP2(IOLD,IOLD) = .TRUE.
                END IF
              END IF
              IF (LPXPL) THEN
!$OMP ATOMIC
                PXPL(IPLSO,NCELL) =PXPL(IPLSO,NCELL)-WEIGHT
                LMETSP(NSPAMI+IPLSO)=.TRUE.
                IF (LSCX) THEN
                  PXPL2(1:NPLS,0:NDXX) => PXPL(:,NCELL)
!$OMP ATOMIC
                  PXPL2(IPLSO,IOLD)=PXPL2(IPLSO,IOLD)-WEIGHT
                  LMETSP2(1:NPLS,0:NDXX) => LMETSP(NTS_PXPLA:NTS_PXPLE)
                  LMETSP2(IPLSO,0) = .TRUE.
                  LMETSP2(IPLSO,IOLD) = .TRUE.
                END IF
              END IF
              IF (LPXEL) THEN
!$OMP ATOMIC
                PXEL(NCELL)      =PXEL(NCELL)-WEIGHT
              ENDIF

C  post-collision estimator
C  IATMN: ATOM SPECIES AFTER CX
              IF (LPXAT) THEN
!$OMP ATOMIC
                PXAT(IATMN,NCELL)=PXAT(IATMN,NCELL)+WEIGHT
                LMETSP(NSPH+IATMN)=.TRUE.
                IF (LSCX) THEN
                  PXAT2(1:NATM,0:NDXX) => PXAT(:,NCELL)
!$OMP ATOMIC
                  PXAT2(IATMN,IOLD)=PXAT2(IATMN,IOLD)+WEIGHT
                  LMETSP2(1:NATM,0:NDXX) => LMETSP(NTS_PXATA:NTS_PXATE)
                  LMETSP2(IATMN,0) = .TRUE.
                  LMETSP2(IATMN,IOLD) = .TRUE.
                END IF
              END IF
              IF (N2NDX(IRCX,1).EQ.4) THEN
C  IPLSN: ION SPECIES AFTER CX
                IPLSN=N2NDX(IRCX,2)
                IF (LPXPL) THEN
!$OMP ATOMIC
                  PXPL(IPLSN,NCELL)=PXPL(IPLSN,NCELL)+WGHTO
                  LMETSP(NSPAMI+IPLSN)=.TRUE.
                  IF (LSCX) THEN
                    PXPL2(1:NPLS,0:NDXX) => PXPL(:,NCELL)
!$OMP ATOMIC
                    PXPL2(IPLSN,IOLD)=PXPL2(IPLSN,IOLD)+WGHTO
                    LMETSP2(1:NPLS,0:NDXX) => 
     .                      LMETSP(NTS_PXPLA:NTS_PXPLE)
                    LMETSP2(IPLSN,0) = .TRUE.
                    LMETSP2(IPLSN,IOLD) = .TRUE.
                  END IF
                END IF
                IF (LPXEL) THEN
!$OMP ATOMIC
                  PXEL(NCELL)      =PXEL(NCELL)+WGHTO
                ENDIF

              ELSEIF (N2NDX(IRCX,1).NE.4) THEN
cdr  col estim. particle tally, but second secondary is not a bulk ion
                GOTO 999
              ENDIF
            ENDIF

c  collision estimator for CX energy exchange tallies
            IF (IESTCX(IRCX,3).NE.0) THEN
              IF (LEXX) THEN
!$OMP ATOMIC
                EXX(NCELL)=EXX(NCELL)-E0O*WGHTO
              ENDIF
              IF (LEXAT) THEN
!$OMP ATOMIC
                EXAT(NCELL)=EXAT(NCELL)+E0*WEIGHT
              ENDIF
              IF (LEXPL) THEN
!$OMP ATOMIC
                EXPL(IPLSO,NCELL)=EXPL(IPLSO,NCELL)-E0*WEIGHT
                LMETSP(NSPAMI+IPLSO)=.TRUE.
              END IF
              IF (N2NDX(IRCX,1).EQ.4) THEN
                IF (LEXPL) THEN
                  IPLSN=N2NDX(IRCX,2)
!$OMP ATOMIC
                  EXPL(IPLSN,NCELL)=EXPL(IPLSN,NCELL)+E0O*WGHTO
                  LMETSP(NSPAMI+IPLSN)=.TRUE.
                ENDIF
              ELSE
cdr  col estim. particle tally, but second secondary is not a bulk ion
                GOTO 999
              ENDIF
            ENDIF

C  COLLISION ESTIMATOR CONTRIBUTION TO MXPL_VEC (FORMERLY: COPV)
            IF (IESTCX(IRCX,2).NE.0) THEN
              IF (LMXPL) THEN
C  SET THE POST-COLLISION TEST PARTICLE PARALLEL VELOCITY = OLD PRE-COLLISION BULK (ION) VELOCITY
                V0_PARB=VEL*(VELX*BX+VELY*BY+VELZ*BZ)
                V0_PARB=V0_PARB*AMUA*RMASSX
                IF (INDPRO(4) == 8) THEN
                  CALL EIRENE_VECUSR(2,NCELL,X0,Y0,Z0,VX,VY,VZ,IPLS,
     .                               .TRUE.)
                  VPLASP=VX*BX+VY*BY+VZ*BZ
                  SIGNUM=SIGN(1._DP,VPLASP)
                ELSE
                  SIGNUM=1._DP
                  IF (LBVIN) SIGNUM =SIGN(1._DP,BVIN(IPLSV,NCLLO))
                ENDIF

C  ASSUME: OLD (INCIDENT) FIELD PARTICE VELOCITY WAS EQUAL TO NEW TEST PARTICLE VELOCITY
C          Scattering angle = PI in COM.
C  remove its momentum from field particles
!$OMP ATOMIC
                MXPL(IPLS,NCELL)=MXPL(IPLS,NCELL)-
     .                           WEIGHT*V0_PARB
     .                           *SIGNUM
                LMETSP(NSPAMI+IPLS)=.TRUE.
              END IF

              IF (N2NDX(IRCX,1).EQ.4) THEN
                IF (LMXPL) THEN
C  IPLSN: BULK ION SPECIES AFTER CX
                  IPLSN=N2NDX(IRCX,2)
C  ASSUME: NEW FIELD PARTICLE VELOCITY IS EQUAL TO (OLD) INCIDENT TEST PARTICLE VELOCITY
C          Scattering angle = PI in COM. Exchange of identity
!$OMP ATOMIC
                  MXPL(IPLSN,NCELL)=MXPL(IPLSN,NCELL)+
     .                              WGHTO*V0_PARBO
     .                              *SIGNUM
                  LMETSP(NSPAMI+IPLSN)=.TRUE.
                END IF
              ELSEIF (N2NDX(IRCX,1).NE.4) THEN
cdr  col estim. momentum tally, but second secondary is not a bulk ion
                GOTO 999
              ENDIF
            ENDIF

            if (iptypo .eq. iptypn) then
cdr continue particle tracing in calling routine
              coltyp=1
            else
cdr also exit from calling routine
              coltyp=2
            endif

            NCELL=NCLLO
            RETURN

          CASE(2)
C  1ST SECONDARY IS MOLECULE IMOL
            IMOL=N1STX(IRCX,2)
            IMOLN=IMOL
            E0=CVRSSM(IMOL)*VELQ
            NNEW=NSPA+IMOL
C
C  CX GENERATION LIMIT, CX, AND SAME SPECIES
            IF (NGENX.GT.0) THEN
              IF (NNEW.EQ.NOLD) THEN
                XGENER=XGENER+1.D0
              ELSE
                XGENER=0.D0
              ENDIF
              IF (XGENER.GE.NGENX) THEN
                CALL EIRENE_GENLIM
                RETURN
              ENDIF
            ENDIF
C
C  FLUID LIMIT, CX, AND SAME SPECIES
            IF (NGENX.LT.0) THEN
              IF (NNEW.EQ.NOLD) THEN
                FP=VELO/SIGVCX(IRCX)  !mfp
                FLTEST=FP/DIST
                IF (FLTEST.LT.FDLMCX(IRCX)) THEN
                  CALL EIRENE_GENLIM
                  RETURN
                ENDIF
              ENDIF
            ENDIF
C
cdr  col estim. particle tally
            IF (IESTCX(IRCX,1).NE.0) GOTO 999
cdr  col estim. momentum tally
            IF (IESTCX(IRCX,2).NE.0) GOTO 999
cdr  col estim. energy tally
            IF (IESTCX(IRCX,3).NE.0) GOTO 999

            if (iptypo .eq. iptypn) then
              coltyp=1
            else
              coltyp=2
            endif

            NCELL = NCLLO
            RETURN

          CASE(3)
C  1ST SECONDARY IS TEST ION
            IION=N1STX(IRCX,2)
            IIONN=IION
            E0=CVRSSI(IION)*VELQ
            NNEW=NSPAM+IION
C
C  CX GENERATION LIMIT, CX, AND SAME SPECIES
            IF (NGENX.GT.0) THEN
              IF (NNEW.EQ.NOLD) THEN
                XGENER=XGENER+1.D0
              ELSE
                XGENER=0.D0
              ENDIF
              IF (XGENER.GE.NGENX) THEN
                CALL EIRENE_GENLIM
                RETURN
              ENDIF
            ENDIF
C
C  FLUID LIMIT, CX, AND SAME SPECIES
            IF (NGENX.LT.0) THEN
              IF (NNEW.EQ.NOLD) THEN
                FP=VELO/SIGVCX(IRCX)  !mfp
                FLTEST=FP/DIST
                IF (FLTEST.LT.FDLMCX(IRCX)) THEN
                  CALL EIRENE_GENLIM
                  RETURN
                ENDIF
              ENDIF
            ENDIF
C
cdr  col estim. particle tally
            IF (IESTCX(IRCX,1).NE.0) GOTO 999
cdr  col estim. momentum tally
            IF (IESTCX(IRCX,2).NE.0) GOTO 999
cdr  col estim. energy tally
            IF (IESTCX(IRCX,3).NE.0) GOTO 999

            if (iptypo .eq. iptypn) then
              coltyp=1
            else
              coltyp=2
            endif

            NCELL = NCLLO
            RETURN

          CASE DEFAULT
            WRITE (iunout,*) ' ITYP = ',ITYP,' AS FIRST SECONDARY IS',
     .                  ' NOT FORESEEN IN COLLIDE '
          END SELECT

        ELSE

C  FOLLOW 2ND SECONDARY, SPEED OF PREVIOUS TEST PARTICLE
          ITYP=N2NDX(IRCX,1)

          ITYPN=ITYP
          IF (ITYPN.EQ.3) THEN
cdr  return to folion 
            IPTYPN=1
          ELSE
cdr  return to folneut
            IPTYPN=0
          ENDIF

          SELECT CASE(ITYP)
C
          CASE(1)
            IATM=N2NDX(IRCX,2)
            XGENER=0.D0
C
            E0=CVRSSA(IATM)*VELO*VELO
cdr  col estim. tally, but second secondary is an atom
            IF (IESTCX(IRCX,1).NE.0) GOTO 999
            IF (IESTCX(IRCX,2).NE.0) GOTO 999
            IF (IESTCX(IRCX,3).NE.0) GOTO 999

            if (iptypo .eq. iptypn) then
              coltyp=1
            else
              coltyp=2
            endif

            NCELL = NCLLO
            RETURN
C
          CASE(2)
            IMOL=N2NDX(IRCX,2)
            XGENER=0.D0
C
            E0=CVRSSM(IMOL)*VELO*VELO
cdr  col estim. tally, but second secondary is a molecule
            IF (IESTCX(IRCX,1).NE.0) GOTO 999
            IF (IESTCX(IRCX,2).NE.0) GOTO 999
            IF (IESTCX(IRCX,3).NE.0) GOTO 999

            if (iptypo .eq. iptypn) then
              coltyp=1
            else
              coltyp=2
            endif

            NCELL = NCLLO
            RETURN
C
          CASE(3)
            IION=N2NDX(IRCX,2)
            XGENER=0.D0
C
            E0=CVRSSI(IION)*VELO*VELO
cdr  col estim. tally, but second secondary is a test ion
            IF (IESTCX(IRCX,1).NE.0) GOTO 999
            IF (IESTCX(IRCX,2).NE.0) GOTO 999
            IF (IESTCX(IRCX,3).NE.0) GOTO 999

            if (iptypo .eq. iptypn) then
              coltyp=1
            else
              coltyp=2
            endif

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
cdr:  at this place to be done: elastic collisions of test ions with field ions.
cdr   in particular: Fokker-Planck (velocity space diffusion--> TAU approximation?)
cdr:  currently still somewhere in folion. To be moved here,
cdr   build on analogy with other elastic collisions
C
        IF (NLTRC) THEN
!$OMP CRITICAL
          CALL EIRENE_CHCTRC(X0,Y0,Z0,16,5)
!$OMP END CRITICAL
        ENDIF

C   FIND IREL, AND SPECIES INDEX IPLS OF BULK (ION) COLLISION PARTNER
        SIGSUM=SIGEIT+SIGCXT
        DO 281 IXEL=1,NXELIM
          IREL=LGXEL(IXEL,0)
          IPLS=LGXEL(IXEL,1)
          SIGSUM=SIGSUM+SIGVEL(IREL)
          IF (ZEP1.LT.SIGSUM) GOTO 282
  281   CONTINUE
        IREL=LGXEL(NXELI,0)
        IPLS=LGXEL(NXELI,1)
  282   CONTINUE
C       GET GLOBAL REACTION NUMBER          
        KK = NREAEL(IREL)

        IPLSV=MPLSV(IPLS)
C
C  NEW SPECIES INDEX AND ENERGY
C       WEIGHT=WEIGHT*1.
C  FOLLOW SECONDARY, NEW SPEED FROM SUBROUTINE VELOEL
C       ITYP=2
        NFLAG=NINT(CFLAG(5,IREL))
        RMXIO=RMASSX
        CALL EIRENE_VELOEL(NCLLO,VELXO,VELYO,VELZO,VELO,IOLD,NOLD,VELQ,
     .              NFLAG,IREL,RMXIO)
C
        INEW=IOLD
C  NOTE: WEIGHT .NE. WGHTO IS POSSIBLE HERE,BECAUSE WEIGHT MAY HAVE CHANGED
C       DUE TO NON-ANALOGUE SAMPLING IN VELOEL
        E0=CVRSSX*VELQ


C  DO NOT UPDATE BGK TALLIES HERE
        IBGK=NPBGKP(IPLS,1)
        IF (IBGK.NE.0) GOTO 300

C  UPDATE COLLISION ESTIMATOR CONTRIBUTION
C  ASSUME, AS BEFORE, NO CHANGE IN SPECIES/TYPE
        IF (IESTEL(IREL,1).NE.0) THEN
          IF (LPXX) THEN
!$OMP ATOMIC
            PXX(IOLD,NCELL) =PXX(IOLD,NCELL)-WGHTO
c  IOLD=INEW for EL processes
!$OMP ATOMIC
            PXX(INEW,NCELL) =PXX(INEW,NCELL)+WEIGHT
            LMETSP(NOLD)=.TRUE.
c  NOLD=NNEW for EL processes
            LMETSP(NNEW)=.TRUE.
            IF (LSCX) THEN
              PXX2(1:NDXX,0:NDXX) => PXX(:,NCELL)
!$OMP ATOMIC
              PXX2(IOLD,IOLD)=PXX2(IOLD,IOLD)-WGHTO
!$OMP ATOMIC
              PXX2(INEW,IOLD)=PXX2(INEW,IOLD)+WEIGHT
              LMETSP2(1:NDXX,0:NDXX) => LMETSP(NDXXA:NDXXE)
              LMETSP2(IOLD,0) = .TRUE.
              LMETSP2(INEW,0) = .TRUE.
              LMETSP2(IOLD,IOLD) = .TRUE.
              LMETSP2(INEW,IOLD) = .TRUE.
            END IF
          END IF
        ENDIF
c  UPDATE NET collision estimator for EL energy exchange tallies
        IF (IESTEL(IREL,3).NE.0) THEN
          EDEL=E0O*WGHTO-E0*WEIGHT
          IF (LEXX) THEN
!$OMP ATOMIC
            EXX(NCELL)      =EXX(NCELL)-EDEL
          ENDIF
          IF (LEXPL) THEN
!$OMP ATOMIC
            EXPL(IPLS,NCELL) =EXPL(IPLS,NCELL)+EDEL
            LMETSP(NSPAMI+IPLS)=.TRUE.
          END IF
        ENDIF

C  UPDATE COLLISION ESTIMATOR CONTRIBUTION TO MXPL_VEC (FORMERLY: COPV)
        IF (IESTEL(IREL,2).NE.0) THEN
          IF (LMXPL) THEN
C  SET THE POST-COLLISION TEST PARTICLE PARALLEL VELOCITY WRT. B FIELD
            V0_PARB=VEL*(VELX*BX+VELY*BY+VELZ*BZ)
            V0_PARB=V0_PARB*AMUA*RMASSX
C
            VDEL=V0_PARBO*WGHTO-V0_PARB*WEIGHT
            IF (INDPRO(4) == 8) THEN
              CALL EIRENE_VECUSR(2,NCELL,X0,Y0,Z0,VX,VY,VZ,IPLS,
     .                           .TRUE.)
              VPLASP=VX*BX+VY*BY+VZ*BZ
              SIGNUM=SIGN(1._DP,VPLASP)
            ELSE
              SIGNUM=1._DP
              IF (LBVIN) SIGNUM=SIGN(1._DP,BVIN(IPLSV,NCLLO))
            ENDIF
!$OMP ATOMIC
            MXPL(IPLS,NCELL)=MXPL(IPLS,NCELL)+VDEL*SIGNUM
            LMETSP(NSPAMI+IPLS)=.TRUE.
          END IF
        ENDIF
  300   CONTINUE
        COLTYP=1
        NCELL = NCLLO
        RETURN
C
C  GENERAL ION IMPACT COLLISION: PI PROCESSES. NOT READY
C
      ELSEIF (ZEP1.LE.SIGEIT+SIGCXT+SIGELT+SIGPIT) THEN
C
        IF (NLTRC) THEN
!$OMP CRITICAL
          CALL EIRENE_CHCTRC(X0,Y0,Z0,16,3)
!$OMP END CRITICAL
        ENDIF

        SIGSUM=SIGEIT+SIGCXT+SIGELT
        DO 261 IXPI=1,NXPIIM
C   FIND INDEX OF THAT ION IMPACT COLLISION
          IRPI=LGXPI(IXPI,0)
          IPLS=LGXPI(IXPI,1)
          SIGSUM=SIGSUM+SIGVPI(IRPI)
          IF (ZEP1.LT.SIGSUM) GOTO 262
  261   CONTINUE
        IRPI=LGXPI(NXPII,0)
        IPLS=LGXPI(NXPII,1)
  262   CONTINUE
C       GET GLOBAL REACTION NUMBER          
        KK = NREAPI(IRPI)
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
C  PRE- COLLISION ESTIMATOR FOR EXX,
C  PRE- AND POST-COLLISION ESTIMATOR FOR EXPL AND EXEL
        IF (IESTPI(IRPI,3).NE.0) THEN
C  score loss of incoming test particle energy
          IF (LEXX) THEN
!$OMP ATOMIC
            EXX(NCELL)=EXX(NCELL)-WEIGHT*E0
          ENDIF

cdr EXPL, EXEL       :  SCORE NET CHANGES HERE.
cdr EXAT, EXML, EXIO :  SCORE EXACT GAINS LATER.
          IF (LEXPL) THEN
            DO IP=1,IPPLPI(IRPI,0)
cdr:  this is incorrect. esigpi must be split into ipl secondaries
              IPL=IPPLPI(IRPI,IP)
              LOGPLS(IPL,ISTRA)=.TRUE.
!$OMP ATOMIC
              EXPL(IPL,NCELL)=EXPL(IPL,NCELL)+WEIGHT*ESIGPI(IRPI,4)
              LMETSP(NSPAMI+IPL)=.TRUE.
            END DO
          END IF
          IF (LEXEL) THEN
!$OMP ATOMIC
            EXEL(NCELL)=EXEL(NCELL)+WEIGHT*ESIGPI(IRPI,5)
          ENDIF
        ENDIF

        IF (IESTPI(IRPI,1).NE.0) THEN
          CCOLEST='PRE COL. PARTICLE RATE, PI PROCESS'
          GOTO 998
        ENDIF
        IF (IESTPI(IRPI,2).NE.0) THEN
          CCOLEST='PRE COL. MOMENTUM RATE, PI PROCESS'
          GOTO 998
        ENDIF
C
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
        NFLAG=NINT(CFLAG(4,IRPI))
        RMXIO=RMASSX

Cdr  PTOT=0,1,2,etc..., = integer,  number of next generation particles

        IF (NLCASCAD .AND. (NLEVEL+PTOT <= MAXLEV)) THEN  ! PI PROCESS CASCADING

          IF (.NOT.ALLOCATED(NAMIPI)) THEN
            ALLOCATE(NAMIPI(NSPAMI))
          END IF
          NAMIPI = 0
          NAMIPI(NSPH+1:NSPA) = INT(PATPI(IRPI,1:NATMI))
          NAMIPI(NSPA+1:NSPAM) = INT(PMLPI(IRPI,1:NMOLI))
          NAMIPI(NSPAM+1:NSPAMI) = INT(PIOPI(IRPI,1:NIONI))

!  RESET WEIGHT TO ORIGINAL VALUE
          WEIGHT=WEIGHT / PTOT

          DO I = NSPAMI, NSPH+1, -1
            DO J=1, NAMIPI(I)
              ZEP = 0.5_DP * (P2NP(IRPI,I-1)+P2NP(IRPI,I))
              CALL EIRENE_VELOPI(NCLLO,VELXO,VELYO,VELZO,VELO,IOLD,
     .                           NOLD,VELQ,NFLAG,IRPI,RMXIO,ZEP)
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
     .                       NOLD,VELQ,NFLAG,IRPI,RMXIO,-1._DP)

        END IF

        XGENER=0.D0

        ITYPN=ITYP
        IF (ITYPN.EQ.3) THEN
cdr  return to folion 
          IPTYPN=1
        ELSE
cdr  return to folneut
          IPTYPN=0
        ENDIF

C
C  UPDATE POST-COLLISION ESTIMATORS CONTRIBUTION TO EXAT, EXML, EXIO
C  NEW TYP: ITYP
        IF (IESTPI(IRPI,3).NE.0) THEN
          IF (ITYP.EQ.1) THEN
            IF (LEXAT) THEN
!$OMP ATOMIC
              EXAT(NCELL)=EXAT(NCELL)+WEIGHT*E0
            ENDIF
          ELSEIF (ITYP.EQ.2) THEN
            IF (LEXML) THEN
!$OMP ATOMIC
              EXML(NCELL)=EXML(NCELL)+WEIGHT*E0
            ENDIF
          ELSEIF (ITYP.EQ.3) THEN
            IF (LEXIO) THEN
!$OMP ATOMIC
              EXIO(NCELL)=EXIO(NCELL)+WEIGHT*E0
            ENDIF
          ENDIF
        ENDIF

        if (iptypo .eq. iptypn) then
          coltyp=1
        else
          coltyp=2
        endif

        NCELL = NCLLO
        RETURN
C
C
      ELSE
C
C
        WRITE (iunout,*) 'ERROR IN COLLIDE, UNKNOWN TYPE OF COLLISION '
        CALL EIRENE_EXIT_OWN(1)
C
C
      ENDIF

      GOTO 999
C

C
  997 WRITE (iunout,*) 'ERROR IN COLLIDE '
      WRITE (iunout,*) 'IREI=  ',IREI,' IS SUPPRESSED, BUT'
      WRITE (iunout,*) 'COLLISION ESTIMATOR WAS SELECTED  '
      WRITE (iunout,*)
     .  'SET WMINV = INFINITY, OR USE TRACKLENGTH ESTIM. '
      CALL EIRENE_EXIT_OWN(1)

  998 WRITE (iunout,*) 'ERROR IN COLLIDE '
      WRITE (iunout,*) 'COLLISION ESTIMATOR WAS SELECTED  '
      WRITE (iunout,*) 'BUT IS NOT READY IN SUBR. COLLIDE'
      WRITE (iunout,*) 'TYPE ',trim(ccolest)
      CALL EIRENE_EXIT_OWN(1)

C
  999 WRITE (iunout,*) 'ERROR IN COLLIDE '
      WRITE (iunout,*) 'ITYP ',ITYP,IPHOT,IATM,IMOL,IION,IPLS
      CALL EIRENE_EXIT_OWN(1)

      CONTAINS

      SUBROUTINE EIRENE_GENLIM
C  UPDATE GENERATION LIMIT TALLIES, THEN STOP TRAJECTORY
C  USE POST-COLLISION WEIGHT, VELOCITY AND ENERGY (NOT: PRE-COLLISION DATA)
C  SHOULD MAKE NO DIFFERENCE ON AVERAGE, IF GENERATION LIMIT IS VALID.
C  IF NOT, ONLY THIS FORM OF ABSORPTION ESTIMATOR GIVES CORRECT BALANCES.

      IF (LPGENX) THEN
!$OMP ATOMIC
        PGENX(NCELL)=PGENX(NCELL)-WEIGHT
      ENDIF
      IF (LEGENX) THEN
!$OMP ATOMIC
        EGENX(NCELL)=EGENX(NCELL)-WEIGHT*E0
      ENDIF
      IF (LVGENX) THEN
C  FIND POST-COLLISION PARALLEL VELOCITY.
C  THE LOCAL B FIELD: KNOWN ALREADY FROM INITIALISATION
        V0_N(0)=VEL*(VELX*BXN(0)+VELY*BYN(0)+VELZ*BZN(0))
c  Mass of post-collision test particle is the same as pre-collision mass RMASSX
        M0_N(0)=V0_N(0)*AMUA*RMASSX
!$OMP ATOMIC
        VGENX(NCELL)=VGENX(NCELL)-WEIGHT*M0_N(0)
      END IF

      IF (LPGENX.OR.LEGENX.OR.LVGENX) LMETSP(NOLD)=.TRUE.

      IF (NLTRC) THEN
!$OMP CRITICAL                
        CALL EIRENE_CHCTRC(X0,Y0,Z0,16,16)
!$OMP END CRITICAL
      ENDIF
cdr  also exit the calling routine
      LGPART=.FALSE.
      ITYP=4
      COLTYP=2
      NCELL = NCLLO
cdr   write (iunout,*) 'genlim, npanu ',npanu
      RETURN
      END SUBROUTINE EIRENE_GENLIM

      END SUBROUTINE EIRENE_COLLIDE

      END MODULE EIRMOD_COLLIDE
