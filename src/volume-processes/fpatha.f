c  25.11.05: option modcol(3,4...)=3 added
c            (adopted from fpatha)
c            cx rate option 4 added (adopted from fpatha)


C               added: jcou,ncou
!pb  30.08.06:  data structure for reaction data redefined
!pb  12.10.06:  modcol revised
!pb  22.11.06:  flag for shift of first parameter to rate_coeff introduced
!pb  28.11.06:  initialization of XSTOR reactivated because of trouble in
!pb             BGK iteration
!pb  22.03.07:  PI reactions revised
cdr  oct.14  :  ftabcx3 added. Full tests still to be done
cdr  oct.14  :  syncronized with fpathm, fpathi

cdr 31.10.14 :  speedup of final cut off evaluations

cdr note:       sgnl_poly evaluations are just the 8th order polynom, 
cdr             plus rcmin,rcmax consideration.
cdr             unless rcmin,rcmax are set (as it is the case currently here), 
cdr             there is no need to call  --> move to in-line 
cdr 06.08.15 :  arguments added to vecusr

cdr dec. 15:    missing: ftabel3
cdr jan. 16:    call to ftabcx3 added and tested for modcol=1 option 



!pb APR  16:    eatds -> eatei
!pb APR  16:    emlds -> emlei
!pb APR  16:    eiods -> eioei
!pb APR  16:    eelds -> eelei
!pb MAY  16:    tabds1 -> tabei1
!pb JUL  16:    ehvds1 -> ehvei1
cdr Nov. 16:    cflag(7,mstor0) rather than cflag(6,3), see comments

cdr aug. 16:    bug fix re EXPO in PI branch
cdr sept.16:    pi process: use v0/vth >> 1. to switch to beam-rate coeff
cdr             ei process: started to check for H.3, H.1 options for EI processes
cdr                         according to v0/vth >> 1. criteria

C
      FUNCTION EIRENE_FPATHA (K,CFLAG,JCOU,NCOU)
C
C   CALCULATE MEAN FREE PATH AND REACTION RATES FOR NEUTRAL
C   "BEAM ATOMS" , SPECIES IATM, OF VELOCITY VEL IN DRIFTING MAXWELLIAN PLASMA-BACKGROUND
C   IN CELL K

C
C   INPUT:
C   IATM      :  ATOM SPECIES INDEX (INPUT VIA COMMON)
C   K         :  CURRENT GRID CELL
C   JCOU, NCOU:  THERE WILL BE NCOU CALLS TO FPATH, FOR SAME TEST PARTICLE
C                COORDINATES. THIS CURRENT CALL IS CALL NO. JCOU.
 
C   OUTPUT: COMMON COMLCA
C           CFLAG: FLAG FOR SAMPLING OF POST COLLISION STATES
C           CFLAG(1,...): EI
C           CFLAG(2,...): NOT IN USE, was DS process class in very old versions
C           CFLAG(3,...): CX
C           CFLAG(4,...): PI
C           CFLAG(5,...): EL
C           CFLAG(6,...): RC
c           CFLAG(7,...): OT
C
C   FLAG FOR POST COLLISION DISTRIBUTION IN VELOCITY SPACE
C  CFLAG(...,IRCL),  IRCL: IREI,..., IRCX,IRPI,IREL,IRRC,IROT
C      =0:   VI: DELTA COLLISION IN VELOCITY SPACE (BUT DIFFERENT
C                                                   SPECIES ALLOWED)
C      =1:   VI: MONOENERGETIC AND ISOTROPIC IN FRAME MOVING WITH BULK SPECIES
C      =2:   VI: DRIFTING MAXWELLIAN
C      =3:   VI: SIGMA-V-WEIGHTED MAXWELLIAN IN FRAME MOVING WITH BULK SPECIES
C      =X    VI: DELTA COLLISION IN VELOCITY SPACE: VI=V0 (BUT DIFFERENT SPECIES ALLOWED)
C                TO BE WRITTEN
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CINIT
      USE EIRMOD_CZT1
      USE EIRMOD_COMPRT
      USE EIRMOD_COMXS
      USE EIRMOD_CESTIM , ONLY: LEA
 
      IMPLICIT NONE
 
      REAL(DP), INTENT(OUT) :: CFLAG(7,MSTOR0)
      INTEGER, INTENT(IN) :: K,JCOU,NCOU
 
      REAL(DP) :: DENIO(NPLS), ZTI(NPLS)
      REAL(DP) :: PVELQ(NPLSV)
      REAL(DP) :: TBPI3(9), TBCX3(9), TBEL3(9), FP(6)
      REAL(DP) :: EPPI3(9), EPCX3(9), EPEL3(9)
      REAL(DP) :: EIRENE_FPATHA,
     .          EIRENE_CROSS, 
     .          EIRENE_RATE_COEFF, EIRENE_SNGL_POLY,
     .          EIRENE_ENERGY_RATE_COEFF, 
     .          CEL, RMN, RLMS, ER, RMI, RMSI, CXS, VEFFQ, TBCX, VEFF,
     .          SIG,  TBEL,
     .          ELTHDUM, CTCHDUM, SIGMAX,  EHEAVY,
     .          DENEL, VX, VY, VZ, PVELQ0, ELAB,
     .          VRELQ, VREL, XC,YC,ZC,
     .          CII, ELB,TII,V0_REL,VI_TH,TEE,VE_TH,
     .          PLS, TBPI, EXPO,
cdr  functions for 'on the fly' evaluation of a&m data
     .          EIRENE_FEELEI1, EIRENE_FEELPI1,
     .          EIRENE_FEHVEI1, EIRENE_FEHVPI3,
     .          EIRENE_FEPLCX3, EIRENE_FEPLPI3, EIRENE_FEPLEL3,
     .          EIRENE_FTABCX3, EIRENE_FTABPI3, 
     .          EIRENE_FTABEI1,

     .          RCMIN, RCMAX,
     .          ERATE
      INTEGER :: IBGK, IAEL, IREL, IAEI, IREI, IAPI,
     .           IRPI, IACX, IRCX, 
     .           II, IF8, JAN, J, I1, I2, KK, IPLSTI,
     .           IPL, IAT, IPLSV, IREAC
C
C  SET DEFAULTS: NO REACTIONS
C
      XSTORV=0.D0
!pb   IF (NCOU.GT.1) THEN
        XSTOR=0.D0
!pb   ENDIF
      EIRENE_FPATHA=1.D10
      SIGMAX=0.D0
C
      IF (LGVAC(K,0)) RETURN
C
C   LOCAL PLASMA PARAMETERS
C
      DENEL=DEIN(K)
      
      DO 2 IPLS=1,NPLSI
        ZTI(IPLS)=ZT1(IPLS,K)
2       DENIO(IPLS)=DIIN(IPLS,K)
C
C  TRANSFORM TEST PARTICLE VELOCITY TO FRAME MOVING WITH BULK SPECIES IPLS
C            PVELQ(IPLS) IS SQUARED THE ATOM VELOCITY IN THESE FRAMES 
C
      PVELQ0=VEL*VEL
      DO 3 IPLS=1,NPLSV
        IF (NLDRFT) THEN
          IF (INDPRO(4) == 8) THEN
            XC=0.
            YC=0.
            ZC=0.
            CALL EIRENE_VECUSR (2,K,XC,YC,ZC,VX,VY,VZ,IPLS,.FALSE.)
          ELSE
            VX=VXIN(IPLS,K)
            VY=VYIN(IPLS,K)
            VZ=VZIN(IPLS,K)
          END IF
          PVELQ(IPLS)=(VELX*VEL-VX)**2+
     .                (VELY*VEL-VY)**2+
     .                (VELZ*VEL-VZ)**2
        ELSE
          PVELQ(IPLS)=PVELQ0
        ENDIF
3     CONTINUE
C
C
C  ELECTRON IMPACT COLLISION - RATE - COEFFICIENT
C  NO MASS SCALING NEEDED FOR BULK ELECTRONS
C
20    IF (LGAEI(IATM,0).EQ.0.OR.LGVAC(K,NPLS+1)) GOTO 30
      DO 10 IAEI=1,NAEII(IATM)
        IREI=LGAEI(IATM,IAEI)
        IF (MODCOL(1,2,IREI).EQ.1) THEN
          IF (NSTORDR >= NRAD) THEN
            SIGVEI(IREI)=TABEI1(IREI,K)
          ELSE
            SIGVEI(IREI)=EIRENE_FTABEI1(IREI,K)
          END IF

        ELSEIF (MODCOL(1,2,IREI).EQ.2) THEN 

          WRITE (IUNOUT,*) 'UNFINISHED OPTION IN FPATHA, EXIT '
          CALL EIRENE_EXIT_OWN(1)
! scale log collision energy to projectile energy for proper isotope, for rate coefficient, i.e. use neutral particle mass
C set hard wired MINIMUM PROJECTILE ENERGY: 0.1 EV
C         ELB=MAX(-2.3_DP,LOG(PVELQ0)+ ???, TO CONVERT TO LOG ENERGY)
          V0_REL=SQRT(PVELQ0)
! scale log temperature to target temperature for proper isotope, for rate coefficient, i.e. use charged particle mass
C         TEE=LOG(TEIN(K))
c thermal velocity at Te
          VE_TH=CVELAA*SQRT(TEIN(K)/PMASSE)
C rather than elb>>tii, one should compare V0_REL and the thermal velocity, then: also ok. for ei processes 
          IF (TEIN(K).LT.TVAC .OR. (V0_REL/VE_TH).GT.10.) THEN

c         IF ((ELB-TII).GT.4.6) THEN
C  HERE: T_E IS SO LOW, THAT ALL ENERGY IS IN TEST PARTICLE MOTION.
c        use cross section times v0, rather than rate coefficient
cdr         WRITE (IUNOUT,*) 'K,EI',K, V0_REL/VE_TH

          ELSEIF ((V0_REL/VE_TH).GE.0.1.AND.(V0_REL/VE_TH).LE.10.) THEN

c  USE H.3 RATE COEFFICIENT, needs tabei3, to be written..... 

          ELSE  ! NORMAL CASE FOR EI COLLISIONS: V0 << VTH
c  use original H.2 rate coefficient, for test particle at rest relative to electron speed
          ENDIF
 

        ELSE
          GOTO 990
        ENDIF
C
        IF (NSTORDR >= NRAD) THEN
          ESIGEI(IREI,5)=EELEI1(IREI,K)
          EHEAVY        =EHVEI1(IREI,K)
        ELSE
          ESIGEI(IREI,5)=EIRENE_FEELEI1(IREI,K)
          EHEAVY        =EIRENE_FEHVEI1(IREI,K)
        ENDIF
C
        ESIGEI(IREI,1)=EATEI(IREI,0,1)*E0+EATEI(IREI,0,2)*EHEAVY
        ESIGEI(IREI,2)=EMLEI(IREI,0,1)*E0+EMLEI(IREI,0,2)*EHEAVY
        ESIGEI(IREI,3)=EIOEI(IREI,0,1)*E0+EIOEI(IREI,0,2)*EHEAVY

        ESIGEI(IREI,4)=EPLEI(IREI,0,1)*E0+EPLEI(IREI,0,2)*EHEAVY
C
        SIGMAX=MAX(SIGMAX,SIGVEI(IREI))
        SIGEIT=SIGEIT+SIGVEI(IREI)
10    CONTINUE
C
C  GENERAL ION IMPACT ON ATOM IATM, BULK ION SPEZIES IPLS=1,NPLSI
C  30--->40
C
30    IF (LGAPI(IATM,0,0).EQ.0) GOTO 40
      DO 36 IAPI=1,NAPII(IATM)
        IRPI=LGAPI(IATM,IAPI,0)
        IPLS=LGAPI(IATM,IAPI,1)
        IPLSV=MPLSV(IPLS)
        IPLSTI=MPLSTI(IPLS)
        IF (LGVAC(K,IPLS)) GOTO 36
C
C  1.) RATE COEFFICIENT
C
        IF (MODCOL(4,2,IRPI).EQ.1) THEN
C  MAXWELL, AT FIXED BEAM ENERGY, MOSTLY E0=0.0
          IF (NSTORDR >= NRAD) THEN
            SIGVPI(IRPI)=TABPI3(IRPI,K,1)
          ELSE
CDR  SIGVPI(IRPI)=FTABPI3 : NOT READY
!pb         KK=NREAPI(IRPI)
!pb         PLS=TIINL(IPLSTI,K)+ADDPI(IRPI,IPLS)
!pb         TBPI = EIRENE_RATE_COEFF(KK,PLS,0._DP,.TRUE.,0,ERATE)*DIIN(IPLS,K)
!pb         SIGVPI(IRPI)=TBPI
c
            SIGVPI(IRPI)=EIRENE_FTABPI3(IRPI,K)
          END IF
        ELSEIF (MODCOL(4,2,IRPI).EQ.2) THEN

C  MODEL 2:
C  BEAM - MAXWELLIAN RATE IN PLASMA FRAME

! scale log collision energy to projectile energy for proper isotope, for rate coefficient, i.e. use neutral particle mass
C set hard wired MINIMUM PROJECTILE ENERGY: 0.1 EV
          ELB=MAX(-2.3_DP,LOG(PVELQ(IPLSV))+EEFPI(IRPI))
          V0_REL=SQRT(PVELQ(IPLSV))
! scale log temperature to target temperature for proper isotope, for rate coefficient, i.e. use charged particle mass
          TII=TIINL(IPLSTI,K)+ADDPI(IRPI,IPLS)
c thermal velocity at Ti
          VI_TH=CVELAA*SQRT(TIIN(IPLSTI,K)/RMASSP(IPLS))
C rather than elb>>tii, one should compare V0_REL and the thermal velocity, then: also ok. for ei processes 
          IF (TIIN(IPLSTI,K).LT.TVAC .OR. (V0_REL/VI_TH).GT.10.) THEN
c         IF ((ELB-TII).GT.4.6) THEN
C  HERE: T_I IS SO LOW, THAT ALL ION ENERGY IS IN DRIFT MOTION.
c           WRITE (IUNOUT,*) 'K,PI',K, exp(ELB-TII),V0_REL/VI_TH
C           HENCE: USE BEAM-BEAM RATE INSTEAD.
            VRELQ=PVELQ(IPLSV)
            VREL=SQRT(VRELQ)
! scale collision energy to proper isotope, for cross section, i.e. use charged particle mass
            ELAB=LOG(VRELQ)+DEFPI(IRPI)  
            IREAC=MODCOL(4,1,IRPI)
            CII=EIRENE_CROSS(ELAB,IREAC,IRPI,FACRPI(IRPI,1),
     .                       'FPATHA PI1')
            SIGVPI(IRPI)=CII*VREL*DENIO(IPLS)
          ELSE
            IF (NSTORDR >= NRAD) THEN
              TBPI3(1:NSTORDT) = TABPI3(IRPI,K,1:NSTORDT)
              FP = 0._DP
              RCMIN = -HUGE(1._DP)
              RCMAX = HUGE(1._DP)
              EXPO = EIRENE_SNGL_POLY(TBPI3,ELB,RCMIN,RCMAX,FP,0,0)
            ELSE
! CALCULATE RATE-COEFFICIENT "ON THE FLY"
              KK=NREAPI(IRPI)

              EXPO = EIRENE_RATE_COEFF(KK,TII,ELB,.FALSE.,0,ERATE)
     .             + DIINL(IPLS,K) + FACRPI(IRPI,2)
            ENDIF
            SIGVPI(IRPI)=EXP(EXPO)
          END IF
      
        ELSEIF (MODCOL(4,2,IRPI).EQ.3) THEN
C  BEAM - BEAM, BUT WITH EFFECTIVE INTERACTION ENERGY
          VRELQ=ZTI(IPLS)+PVELQ(IPLSV)
          VREL=SQRT(VRELQ)
          ELAB=LOG(VRELQ)+DEFPI(IRPI)
          IREAC=MODCOL(4,1,IRPI)
          CII=EIRENE_CROSS(ELAB,IREAC,IRPI,FACRPI(IRPI,1),'FPATHA II')
          SIGVPI(IRPI)=CII*VREL*DENIO(IPLS)
        ELSE
          GOTO 991
        ENDIF
C
C  2.A ELECTRON ENERGY LOSS PER COLLISION (EV)
C
        IF (NSTORDR >= NRAD) THEN
          ESIGPI(IRPI,5)=EELPI1(IRPI,K)
          EHEAVY        =EHVPI3(IRPI,K,1)
        ELSE
          ESIGEI(IREI,5)=EIRENE_FEELPI1(IRPI,K)
          EHEAVY        =EIRENE_FEHVPI3(IRPI,K)
        ENDIF
C
C  2.B SECONDARY PARTICLE ENERGY LOSS/GAIN PER COLLISION (EV)
C
        ESIGPI(IRPI,1)=EATPI(IRPI,0,1)*E0+EATPI(IRPI,0,2)*EHEAVY
        ESIGPI(IRPI,2)=EMLPI(IRPI,0,1)*E0+EMLPI(IRPI,0,2)*EHEAVY
        ESIGPI(IRPI,3)=EIOPI(IRPI,0,1)*E0+EIOPI(IRPI,0,2)*EHEAVY

        ESIGPI(IRPI,4)=EPLPI(IRPI,0,1)*E0+EPLPI(IRPI,0,2)*EHEAVY
C
        SIGMAX=MAX(SIGMAX,SIGVPI(IRPI))
        SIGPIT=SIGPIT+SIGVPI(IRPI)
C
C  2.C BULK ION ENERGY LOSS PER COLLISION (EV)
C
cdr     IF (NSTORDR >= NRAD) THEN
cdr       ESIGPI(IRPI,4)=EPLPI3(IRPI,K,1)
cdr     ELSE
cdr       ESIGPI(IRPI,4)=EIRENE_FEPLPI3(IRPI,K)
cdr     END IF

        CFLAG(4,IRPI)=1
c
36    CONTINUE
C
C  CHARGE EXCHANGE RATE COEFFICIENT OF ATOMS IATM
C  WITH BULK IONS OF SPEZIES IPLS=1,NPLSI
C  40--->50
C
40    CONTINUE
      IF (LGACX(IATM,0,0).EQ.0.OR.LGVAC(K,0)) GOTO 50
      DO 41 IACX=1,NACXI(IATM)
        IRCX=LGACX(IATM,IACX,0)
        IPLS=LGACX(IATM,IACX,1)
        IPLSTI=MPLSTI(IPLS)
        IPLSV=MPLSV(IPLS)
        IF (LGVAC(K,IPLS)) GOTO 41
C
C  1.) RATE COEFFICIENT
C
        IF (MODCOL(3,2,IRCX).EQ.1) THEN
C  MODEL 1:
C  MAXWELLIAN RATE, IGNORE NEUTRAL VELOCITY
          IF (NSTORDR >= NRAD) THEN
            SIGVCX(IRCX)=TABCX3(IRCX,K,1)
          ELSE
            SIGVCX(IRCX)=EIRENE_FTABCX3(IRCX,K)
          END IF

        ELSEIF (MODCOL(3,2,IRCX).EQ.2) THEN
C  MODEL 2:
C  BEAM - MAXWELLIAN RATE IN PLASMA FRAME
          IF (TIIN(IPLSTI,K).LT.TVAC) THEN
C  HERE: T_I IS SO LOW, THAT ALL ION ENERGY IS IN DRIFT MOTION.
C           HENCE: USE BEAM-BEAM RATE INSTEAD.
            VRELQ=PVELQ(IPLSV)
            VREL=SQRT(VRELQ)
C   PMASS FOR CROSS SECTION RELATIVE VELOCITY
            ELAB=LOG(VRELQ)+DEFCX(IRCX)
            IREAC=MODCOL(3,1,IRCX)
            CXS=EIRENE_CROSS(ELAB,IREAC,IRCX,FACRCX(IRCX,1),
     .                       'FPATHA CX1')
            SIGVCX(IRCX)=CXS*VREL*DENIO(IPLS)
          ELSE
C  MINIMUM PROJECTILE ENERGY: 0.1 EV
cdr         if (LOG(PVELQ(IPLSV))+EEFCX(IRCX).le.-2.3) then
cdr           elb=LOG(PVELQ(IPLSV))+EEFCX(IRCX)
cdr           write (6,*) 'elb in fpatha-1 ',elb, exp(elb)
cdr         endif
C   TMASS FOR RATE COEFF. BEAM VELOCITY
            ELB=MAX(-2.3_DP,LOG(PVELQ(IPLSV))+EEFCX(IRCX))
            IF (NSTORDR >= NRAD) THEN
! DOUBLE POLYNOMIAL FIT REDUCED TO SINGLE POLYNOMIAL FIT BY
! PRECALCULATING TEMPERATURE DEPENDENCIES
              TBCX3(1:NSTORDT) = TABCX3(IRCX,K,1:NSTORDT)
              FP = 0._DP
              RCMIN = -HUGE(1._DP)
              RCMAX = HUGE(1._DP)
              EXPO = EIRENE_SNGL_POLY(TBCX3,ELB,RCMIN,RCMAX,FP,0,0)
            ELSE
! CALCULATE RATE-COEFFICIENT
CDR  THIS SHOULD BE DONE IN FTABCX3.  NOT READY
              KK=NREACX(IRCX)
              TII=TIINL(IPLSTI,K)+ADDCX(IRCX,IPLS)
              EXPO = EIRENE_RATE_COEFF(KK,TII,ELB,.FALSE.,0,ERATE)
     .               + DIINL(IPLS,K) + FACRCX(IRCX,2)
            END IF
            SIGVCX(IRCX)=EXP(EXPO)
          ENDIF
        ELSEIF (MODCOL(3,2,IRCX).EQ.3) THEN
C  MODEL 3:  (ALSO:  DEFAULT CX MODEL, ONLY CROSS SECTION IS USED, NO RATE COEFFICIENTS)
C  BEAM - BEAM RATE, BUT WITH EFFECTIVE INTERACTION ENERGY
          VEFFQ=ZTI(IPLS)+PVELQ(IPLSV)
          VEFF=SQRT(VEFFQ)
          ELAB=LOG(VEFFQ)+DEFCX(IRCX)
          IREAC=MODCOL(3,1,IRCX)
          CXS=EIRENE_CROSS(ELAB,IREAC,IRCX,FACRCX(IRCX,1),'FPATHA CX2')
          SIGVCX(IRCX)=CXS*VEFF*DENIO(IPLS)
        ELSEIF (MODCOL(3,2,IRCX).EQ.4) THEN
C  MODEL 4
C  BEAM - BEAM RATE, IGNORE THERMAL ION ENERGY
          VRELQ=PVELQ(IPLSV)
          VREL=SQRT(VRELQ)
          ELAB=LOG(VRELQ)+DEFCX(IRCX)
          IREAC=MODCOL(3,1,IRCX)
          CXS=EIRENE_CROSS(ELAB,IREAC,IRCX,FACRCX(IRCX,1),'FPATHA CX3')
          SIGVCX(IRCX)=CXS*VREL*DENIO(IPLS)
        ELSE
          GOTO 992
        ENDIF
C
        SIGMAX=MAX(SIGMAX,SIGVCX(IRCX))
        SIGCXT=SIGCXT+SIGVCX(IRCX)
C
C  2.) BULK ION ENERGY LOSS RATE:
C
        IF (MODCOL(3,4,IRCX).EQ.1) THEN
C  MODEL 1:
C  MEAN ENERGY FROM DRIFTING MAXWELLIAN
C  (ONLY NEEDED FOR TRACKLENGTH ESTIMATOR)
C  ION SAMPLING FROM MAXWELLIAN
          IF (LEA.AND.(IESTCX(IRCX,3).EQ.0)) THEN  ! for tracklength estimator only
            IF (NSTORDR >= NRAD) THEN
              ESIGCX(IRCX,1)=EPLCX3(IRCX,K,1)
            ELSE
              ESIGCX(IRCX,1)=EIRENE_FEPLCX3(IRCX,K)
            END IF
          END IF  ! this was for tracklength estimator only
          CFLAG(3,IRCX)=2
        ELSEIF (MODCOL(3,4,IRCX).EQ.2) THEN
C  MODEL 2:
C  MEAN ENERGY FROM CROSS SECTION WEIGHTED DRIFTING MAXWELLIAN
C  (ONLY NEEDED FOR TRACKLENGTH ESTIMATOR)
C  ION SAMPLING FROM WEIGHTED DRIFTING MAXWELLIAN (E.G., BY REJECTION)
          IF (LEA.AND.(IESTCX(IRCX,3).EQ.0)) THEN  ! for tracklength estimator only
C  MINIMUM PROJECTILE ENERGY: 0.1 EV
cdr         if (LOG(PVELQ(IPLSV))+EEFCX(IRCX).le.-2.3) then
cdr           elb=LOG(PVELQ(IPLSV))+EEFCX(IRCX)
cdr           write (6,*) 'elb in fpatha-2 ',elb, exp(elb)
cdr         endif
            ELB=MAX(-2.3_DP,LOG(PVELQ(IPLSV))+EEFCX(IRCX))
            IF (NSTORDR >= NRAD) THEN
! DOUBLE POLYNOMIAL FIT REDUCED TO SINGLE POLYNOMIAL FIT BY
! PRECALCULATING TEMPERATURE DEPENDENCIES
              EPCX3(1:NSTORDT) = EPLCX3(IRCX,K,1:NSTORDT)
              FP = 0._DP
              RCMIN = -HUGE(1._DP)
              RCMAX = HUGE(1._DP)
              EXPO = EIRENE_SNGL_POLY(EPCX3,ELB,RCMIN,RCMAX,FP,0,0)
            ELSE
! CALCULATE ENERGY-WEIGHTED RATE-COEFFICIENT ON THE FLY
              KK=NELRCX(IRCX)
              TII=TIINL(IPLSTI,K)+ADDCX(IRCX,IPLS)
              EXPO = EIRENE_ENERGY_RATE_COEFF(KK,TII,ELB,.FALSE.,0)
     .               + DIINL(IPLS,K) + FACRCX(IRCX,2)
            END IF
            ESIGCX(IRCX,1)=EXP(EXPO)/SIGVCX(IRCX)
            ESIGCX(IRCX,1)=ESIGCX(IRCX,1)+EDRIFT(IPLS,K)
          ENDIF  ! this was for tracklength estimator only
          CFLAG(3,IRCX)=3
        ELSEIF (MODCOL(3,4,IRCX).EQ.3) THEN
C  MODEL 3:
C  MEAN ENERGY FROM DRIFTING ISOTROPIC ONE SPEED DISTRIBUTION
C  (ONLY NEEDED FOR TRACKLENGTH ESTIMATOR)
C  ION SAMPLING FROM WEIGHTED DRIFTING ISOTROPIC ONE SPEED DISTRIBUTION
          IF (LEA.AND.(IESTCX(IRCX,3).EQ.0)) THEN  ! for tracklength estimator only
            IF (NSTORDR >= NRAD) THEN
              ESIGCX(IRCX,1)=EPLCX3(IRCX,K,1)
            ELSE
              ESIGCX(IRCX,1)=EIRENE_FEPLCX3(IRCX,K)
            END IF
          ENDIF  ! this was for tracklength estimator only
          CFLAG(3,IRCX)=1
        ELSE
          GOTO 992
        ENDIF
41    CONTINUE
C
C  ELASTIC COLLISIONS OF ATOMS IATM  WITH IONS OF SPEZIES IPLS=1,NPLSI
C  50--->60
C
50    CONTINUE
      IF (LGAEL(IATM,0,0).EQ.0.OR.LGVAC(K,0)) GOTO 60
      DO 51 IAEL=1,NAELI(IATM)
        IREL=LGAEL(IATM,IAEL,0)
        IPLS=LGAEL(IATM,IAEL,1)
        IPLSTI=MPLSTI(IPLS)
        IPLSV=MPLSV(IPLS)
        IBGK=NPBGKP(IPLS,1)
        IF (LGVAC(K,IPLS)) GOTO 51
C
C  1.) RATE COEFFICIENT
C
        IF (MODCOL(5,2,IREL).EQ.1) THEN
C  MODEL 1:
C  MAXWELLIAN RATE, IGNORE ATOM VELOCITY
          IF (NSTORDR >= NRAD) THEN
            SIGVEL(IREL)=TABEL3(IREL,K,1)
          ELSE
cdr  here should be call to ftabel3,  to be done
            KK=NREAEL(IREL)
            TII=TIINL(IPLSTI,K)+ADDEL(IREL,IPLS)
            TBEL = EIRENE_RATE_COEFF(KK,TII,0._DP,.TRUE.,0,ERATE)*
     .             DIIN(IPLS,K)
            SIGVEL(IREL)=TBEL
          END IF
        ELSEIF (MODCOL(5,2,IREL).EQ.2) THEN
C  BEAM - MAXWELL
          IF (TIIN(IPLSTI,K).LT.TVAC) THEN
C  TEMPERATURE TOO LOW, USE: BEAM_ATOM - BEAM_DRIFT RATECOEFF.
            VRELQ=PVELQ(IPLSV)
            VREL=SQRT(VRELQ)
            ELAB=LOG(VRELQ)+DEFEL(IREL)
            IREAC=MODCOL(5,1,IREL)
            CEL=EIRENE_CROSS(ELAB,IREAC,IREL,FACREL(IREL,1),
     .                       'FPATHA EL1')
            SIGVEL(IREL)=CEL*VREL*DENIO(IPLS)
          ELSE
C  MINIMUM PROJECTILE ENERGY: 0.1 EV
            ELB=MAX(-2.3_DP,LOG(PVELQ(IPLSV))+EEFEL(IREL))
            IF (NSTORDR >= NRAD) THEN
! DOUBLE POLYNOMIAL FIT REDUCED TO SINGLE POLYNOMIAL FIT BY
! PRECALCULATING TEMPERATURE DEPENDENCIES
              TBEL3(1:NSTORDT) = TABEL3(IREL,K,1:NSTORDT)
              FP = 0._DP
              RCMIN = -HUGE(1._DP)
              RCMAX = HUGE(1._DP)
              EXPO = EIRENE_SNGL_POLY(TBEL3,ELB,RCMIN,RCMAX,FP,0,0)
            ELSE
! CALCULATE RATE-COEFFICIENT
              KK=NREAEL(IREL)
              TII=TIINL(IPLSTI,K)+ADDEL(IREL,IPLS)
              EXPO = EIRENE_RATE_COEFF(KK,TII,ELB,.FALSE.,0,ERATE)
     .               + DIINL(IPLS,K) + FACREL(IREL,2)
            END IF
            SIGVEL(IREL)=EXP(EXPO)
          ENDIF
        ELSEIF (MODCOL(5,2,IREL).EQ.3) THEN
C  BEAM - BEAM RATE, BUT WITH EFFECTIVE INTERACTION ENERGY
          VEFFQ=ZTI(IPLS)+PVELQ(IPLSV)
          VEFF=SQRT(VEFFQ)
          ELAB=LOG(VEFFQ)+DEFEL(IREL)
          IREAC=MODCOL(5,1,IREL)
C  FIND SIGMA FROM OAK RIDGE "ELASTIC" DATA TABLES
          IF (LHABER) THEN
            RMN=RMASSA(IATM)
            RMI=RMASSP(IPLS)
            RMSI=1./(RMN+RMI)
            RLMS=RMN*RMI*RMSI
            ER=RLMS*VEFFQ*CVELI2
cdr  flag -1.0_DP: only sigma(ER), but no scattering angle evaluated
            CALL EIRENE_SCATANG (ER,-1.0_DP,ELTHDUM,CTCHDUM,SIG)
            CEL= SIG*AU_TO_CM2
          ELSE
C  FIND SIGMA FROM AMJUEL DATA TABLES (BACHMANN ET AL.)
            CEL=EIRENE_CROSS(ELAB,IREAC,IREL,FACREL(IREL,1),
     .                       'FPATHA EL2')
          END IF
          SIGVEL(IREL)=CEL*VEFF*DENIO(IPLS)
        ELSEIF (MODCOL(5,2,IREL).EQ.4) THEN
C  MODEL 4
C  BEAM - BEAM RATE, IGNORE THERMAL ION ENERGY
          VRELQ=PVELQ(IPLSV)
          VREL=SQRT(VRELQ)
          ELAB=LOG(VRELQ)+DEFEL(IREL)
          IREAC=MODCOL(5,1,IREL)
          CEL=EIRENE_CROSS(ELAB,IREAC,IREL,FACREL(IREL,1),'FPATHA EL3')
          SIGVEL(IREL)=CEL*VREL*DENIO(IPLS)
        ELSE
          GOTO 995
        ENDIF
 
        SIGMAX=MAX(SIGMAX,SIGVEL(IREL))
        SIGELT=SIGELT+SIGVEL(IREL)
 
        IF (IBGK.NE.0) SIGBGK=SIGBGK+SIGVEL(IREL)
C
C  2.) BULK ION ENERGY LOSS RATE:
C
        IF (MODCOL(5,4,IREL).EQ.1) THEN
C  MODEL 1:
C  MEAN ENERGY FROM DRIFTING MAXWELLIAN
C  (ONLY NEEDED FOR TRACKLENGTH ESTIMATOR)
C  ION SAMPLING FROM MAXWELLIAN
          IF (NSTORDR >= NRAD) THEN
            ESIGEL(IREL,1)=EPLEL3(IREL,K,1)
          ELSE
            ESIGEL(IREL,1)=EIRENE_FEPLEL3(IREL,K)
          END IF
          CFLAG(5,IREL)=2
        ELSEIF (MODCOL(5,4,IREL).EQ.2) THEN
C  MODEL 2:
C  MEAN ENERGY FROM CROSS SECTION WEIGHTED DRIFTING MAXWELLIAN
C  (ONLY NEEDED FOR TRACKLENGTH ESTIMATOR)
C  ION SAMPLING FROM WEIGHTED DRIFTING MAXWELLIAN (E.G., BY REJECTION)
          IF (IESTEL(IREL,3).EQ.0) THEN  ! for tracklength estimator only
C  MINIMUM PROJECTILE ENERGY: 0.1 EV
            ELB=MAX(-2.3_DP,LOG(PVELQ(IPLSV))+EEFEL(IREL))
            IF (NSTORDR >= NRAD) THEN
! DOUBLE POLYNOMIAL FIT REDUCED TO SINGLE POLYNOMIAL FIT BY
! PRECALCULATING TEMPERATURE DEPENDENCIES
              EPEL3(1:NSTORDT) = EPLEL3(IREL,K,1:NSTORDT)
              FP = 0._DP
              RCMIN = -HUGE(1._DP)
              RCMAX = HUGE(1._DP)
              EXPO = EIRENE_SNGL_POLY(EPEL3,ELB,RCMIN,RCMAX,FP,0,0)
            ELSE
! CALCULATE ENERGY-WEIGHTED RATE-COEFFICIENT ON THE FLY
              KK=NELREL(IREL)
              TII=TIINL(IPLSTI,K)+ADDEL(IREL,IPLS)
              EXPO = EIRENE_ENERGY_RATE_COEFF(KK,TII,ELB,.FALSE.,0)
     .               + DIINL(IPLS,K) + FACREL(IREL,2)
            END IF
            ESIGEL(IREL,1)=EXP(EXPO)/SIGVEL(IREL)
            ESIGEL(IREL,1)=ESIGEL(IREL,1)+EDRIFT(IPLS,K)
          ENDIF  ! this was for tracklength estimator only
          CFLAG(5,IREL)=3
        ELSEIF (MODCOL(5,4,IREL).EQ.3) THEN
C  MODEL 3:
C  MEAN ENERGY FROM DRIFTING ISOTROPIC ONE SPEED DISTRIBUTION
C  (ONLY NEEDED FOR TRACKLENGTH ESTIMATOR)
C  ION SAMPLING FROM WEIGHTED DRIFTING ISOTROPIC ONE SPEED DISTRIBUTION
          IF (NSTORDR >= NRAD) THEN
            ESIGEL(IREL,1)=EPLEL3(IREL,K,1)
          ELSE
            ESIGEL(IREL,1)=EIRENE_FEPLEL3(IREL,K)
          END IF
          CFLAG(5,IREL)=1
        ELSE
          GOTO 995
        ENDIF
51    CONTINUE
C
60    CONTINUE
C
C     TOTAL
C
100   CONTINUE

C
C  CUT OF RESIDUAL RATES, WHICH SHOULD STRICTLY BE ZERO
C  TO AVOID SPURIOUS ENTRIES TO COLLISION RATE TALLIES
C  CURRENTLY: CUT OFF AT 1E-10 TIMES SIGMAX
C
      IF (SIGEIT.GT.0._DP) THEN
        DO IAEI=1,NAEII(IATM)
          IREI=LGAEI(IATM,IAEI)
          IF (SIGVEI(IREI) .LE. SIGMAX*1.D-10) THEN
            SIGEIT=SIGEIT-SIGVEI(IREI)
            SIGVEI(IREI) = 0.D0
          END IF
        END DO
      END IF
 
      IF (SIGPIT.GT.0._DP) THEN
        DO IAPI=1,NAPII(IATM)
          IRPI=LGAPI(IATM,IAPI,0)
          IF (SIGVPI(IRPI) .LE. SIGMAX*1.D-10) THEN
            SIGPIT=SIGPIT-SIGVPI(IRPI)
            SIGVPI(IRPI) = 0.D0
          END IF
        END DO
      END IF
 
      IF (SIGCXT.GT.0._DP) THEN
        DO IACX=1,NACXI(IATM)
          IRCX=LGACX(IATM,IACX,0)
          IF (SIGVCX(IRCX) .LE. SIGMAX*1.D-10) THEN
            SIGCXT=SIGCXT-SIGVCX(IRCX)
            SIGVCX(IRCX) = 0.D0
          END IF
        END DO
      END IF
C
      IF (SIGELT.GT.0._DP) THEN
        DO IAEL=1,NAELI(IATM)
          IREL=LGAEL(IATM,IAEL,0)
          IF (SIGVEL(IREL) .LE. SIGMAX*1.D-10) THEN
            SIGELT=SIGELT-SIGVEL(IREL)
            SIGVEL(IREL) = 0.D0
          END IF
        END DO
      END IF
C
      SIGTOT=SIGEIT+SIGPIT+SIGCXT+SIGELT
      IF (SIGTOT.GT.1.D-20) THEN
        EIRENE_FPATHA=VEL/SIGTOT
        ZMFPI=1./EIRENE_FPATHA
      ENDIF
C
      RETURN
990   CONTINUE
      WRITE (iunout,*) 'ERROR IN FPATHA: INCONSISTENT ELEC. IMP. DATA'
      WRITE (iunout,*) 'IATM,IREI,MODCOL(1,J=1,4,IREI) '
      WRITE (iunout,*) IATM,IREI,(MODCOL(1,J,IREI),J=1,4)
      CALL EIRENE_EXIT_OWN(1)
991   CONTINUE
      WRITE (iunout,*) 'ERROR IN FPATHA: INCONSISTENT ION IMP. DATA'
      WRITE (iunout,*) 'IATM,IRPI,MODCOL(4,J,IRPI) '
      WRITE (iunout,*) IATM,IRPI,(MODCOL(4,J,IRPI),J=1,4)
      CALL EIRENE_EXIT_OWN(1)
992   CONTINUE
      WRITE (iunout,*)
     .  'ERROR IN FPATHA: INCONSISTENT CHARGE EXCHANGE DATA'
      WRITE (iunout,*) 'IATM,IRCX,(MODCOL(3,J,IRCX),J=1,4 '
      WRITE (iunout,*) IATM,IRCX,(MODCOL(3,J,IRCX),J=1,4)
      CALL EIRENE_EXIT_OWN(1)
995   CONTINUE
      WRITE (iunout,*)
     .  'ERROR IN FPATHA: INCONSISTENT ELASTIC COLL. DATA'
      WRITE (iunout,*) 'IATM,IREL,(MODCOL(5,J,IREL),J=1,4) '
      WRITE (iunout,*) IATM,IREL,(MODCOL(5,J,IREL),J=1,4)
      CALL EIRENE_EXIT_OWN(1)
      END
