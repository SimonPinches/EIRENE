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
cdr May 17: some spelling error corrections in comments adopted from ITER branch
c            AE: analog, --> BE: analogue, etc..




      SUBROUTINE EIRENE_COLPHOT(CFLAG,COLTYP,DIST)
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
     .           IBGK, IAEL, IREL, IP, IMEI, IMCX, IAPI, NFLAG,
     .           IATMN, IPLSN, IRPI, NCLLO, IPLSV, IMPI, IIPI, I, J, IPL
      INTEGER :: NEII_RED,LGEI_RED(0:NREI)

Cdr  additional arrays for  ANALOG CASCADE and SPLITTING AT COLLISIONS. 
Cdr (should be set in initialization phase, not here)
CDR  check: are the corresponding arrays PATEI,PMLEI, PIOEI real or integer (1/2 particle possible?)
      INTEGER, ALLOCATABLE, SAVE :: NAMIEI(:),NAMIPI(:)
 
 
csw add n 2lines
      INTEGER :: iaot,irot,kk,updf,t1
      real(dp):: sump
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
C  CHARGE-EXCHANGE:
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
