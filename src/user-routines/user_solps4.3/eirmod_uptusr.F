! 23.08.06: VPX, VPY, VRX, VRY changed to ALLOCATABLE, SAVE to speed up
!           subroutine call (save time in storage allocation)
cdr Jan 17: remove local allocatable cndyn.. arrays. These are now
cdr         set in code initialisation phase
C
      MODULE EIRMOD_UPTUSR

      PUBLIC

      CONTAINS

C
      SUBROUTINE EIRENE_UPTUSR(XSTOR2,XSTORV2,WV,IFLAG)
C
C  USER-SUPPLIED TRACKLENGTH ESTIMATOR, VOLUME-AVERAGED

C  ALSO: SUMMED OVER IRCX (ALL CX PROCESSES)

C  THIS VERSION:

CCC   1     PARTICLE CX RATE (ONLY TOTAL, ZERO PARTICLE SOURCE WITH THIS PROCESS, IN SINGLE SPECIES RUN)
CCC   2-4   ENERGY   CX RATES

C    ADDV(IATM,ICELL)         : VOLUMETRIC CX RATE (REACTIONS/S/CM**3), ATOMS

C    ADDV(  NATMI+IATM,ICELL) : VOLUMETRIC INCIDENT ION ENERGY CX RATE (EV/S/CM**3), ATOMS
C    ADDV(2*NATMI+IATM,ICELL) : VOLUMETRIC INCIDENT NEUTRAL ENERGY CX RATE (EV/S/CM**3), ATOMS
C    ADDV(3*NATMI+IATM,ICELL) : VOLUMETRIC NET ION ENERGY CX RATE (EV/S/CM**3), ATOMS
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CESTIM
      USE EIRMOD_CSDVI
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT
      USE EIRMOD_CUPD
      USE EIRMOD_COMXS
      USE EIRMOD_CSPEZ
      USE EIRMOD_CGRID
      USE EIRMOD_CLOGAU
      USE EIRMOD_COMSIG
      USE EIRMOD_CGEOM
      USE EIRMOD_CZT1

      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: XSTOR2(MSTOR1,MSTOR2,N2ND+N3RD),
     .                        XSTORV2(NSTORV,N2ND+N3RD),
     .                        WV
      INTEGER, INTENT(IN) :: IFLAG
      INTEGER :: ICOU,K,IRD,IACX,IRCX,IRDO,nti,nte,ia,
     .           IMCX,IMEL,IREL,IFIRST,I
      REAL(DP) :: DIST,WTR,WTRSIG
      REAL(DP) :: VSIG_PARB(NPLS), VAL_PARB(NPLS),
     .            VSIG_PERP(NPLS), VAL_PERP(NPLS),
     .            V0_PARB,PARMOM_0,
     .            V0_PERP,PERPMOM_0
     .           ,DD, VR, VP
      REAL(DP), ALLOCATABLE, SAVE :: VPX(:),VPY(:),VRX(:),VRY(:)
      DATA IFIRST/0/

      IF (IFIRST.EQ.0) THEN
        IFIRST=1
C
CDR
CDR  PROVIDE A RADIAL UNIT VECTOR PER CELL
CDR  VPX,VPY, NEEDED FOR PROJECTING PARTICLE VELOCITIES
CDR  SAME FOR POLOIDAL UNIT VECTOR VRX,VRY
C
        if(allocated(vpx)) deallocate(vpx,vpy,vrx,vry)
        ALLOCATE (VPX(NRAD))
        ALLOCATE (VPY(NRAD))
        ALLOCATE (VRX(NRAD))
        ALLOCATE (VRY(NRAD))
        VPX=0.
        VPY=0.
        VRX=0.
        VRY=0.
        IF (LBXIN.AND.LBYIN.AND.LBXPERP.AND.LBYPERP) THEN
!          DO I=1,ntrii
          DO I=1,nrad
cdr  probably incorrect coding in infcop. pln.tri, ppln.tri: do not use.
!           VPX(I)=PLNXTRI(i)    ! radial unit vector
!           VPY(I)=PLNYTRI(i)    ! => bxperp, byperp
!           VRX(I)=PPLNXTRI(i)   ! poloidal unit vector
!           VRY(I)=PPLNYTRI(i)   ! => BXIN, BYIN TO BE NORMALIZED
cdr tbd: check liftali: bx,by,,bxperp,byperp
            VPX(I)=BXPERP(i)     ! radial unit vector <= bxperp, byperp
            VPY(I)=BYPERP(i)     !
            DD = SQRT(BXIN(I)**2 + BYIN(I)**2)
            VRX(I)=BXIN(i)/DD    ! poloidal unit vector <= BXIN, BYIN TO BE NORMALIZED
            VRY(I)=BYIN(i)/DD    !
          END DO
        END IF
      ENDIF

C
C  ON INPUT: WV=WEIGHT/VEL


C
C  ATOMS, CX ENERGY
      IF (ITYP.NE.1) GOTO 999

C  CHECK: STORAGE FOR AT LEAST 12 ADDITIONAL TRACKLENGTH-ESTIMATED TALLIES?
!pb      IF (NADV.LT.MOD_ADDV+12*NATMI) THEN
      IF (NADV.LT.12*NATMI) THEN
        GOTO 9999
      ELSE
C  THIS ROUTINE: SCORE 12 ADDITIONAL TALLIES ADDV
c  initial species index for ADDV tally: nspan
        nti=NSPAN(ntala)
        nte=Nti+12*NATM-1
        LMETSP(nti:nte)=.TRUE.
      ENDIF
C
      DO 200 ICOU=1,NCOU
          DIST=CLPD(ICOU)
          WTR=WV*DIST
          IRDO=NRCELL+NUPC(ICOU)*NR1P2+NBLCKA
          IRD=NCLTAL(IRDO)
c  set parallel plasma flow parameters
c  assume here: bvin, parmom are set in plasma_deriv.
c               In case of other options (indpro): see update.f
cdr only signum needed: default tbd: signum=1
          IF (LBVIN) THEN
            VAL_PARB(1:NPLSI) = BVIN(MPLSV(1:NPLSI),IRDO)
          ELSE
            VAL_PARB(1:NPLSI) = 0._DP
          END IF
          IF (LPARMOM) THEN
            VSIG_PARB(1:NPLSI) = PARMOM(1:NPLSI,IRDO)
          ELSE
            VSIG_PARB(1:NPLSI) = 0._DP
          END IF
C  for the time being: no perpendicular plasma flow.
         VAL_PERP(1:NPLSI) =0.0
         VSIG_PERP(1:NPLSI)=0.0


C
cdr       IF (LGVAC(IRDO,0)) GOTO 200
C
          if (ncou.gt.1) then
            XSTOR(:,:) = XSTOR2(:,:,ICOU)
            XSTORV(:)  = XSTORV2(:,ICOU)
          endif

C
C
          IF (LGACX(IATM,0,0).EQ.0) GOTO 590
            DO 560 IACX=1,NACXI(IATM)
!pb  MOD_ADDV is no incremental value. It is a flag indicating whether all the
!pb  rates used for emissivity lines are to be stored or whether storage saving
!pb  mode ist to be used, only storing the rates for the latest used line
!pb           IA=MOD_ADDV   !  increment for addv tally 1st index: cx energy, atoms
              IA=0
              IRCX=LGACX(IATM,IACX,0)
              IPLS=LGACX(IATM,IACX,1)
              IF (LGVAC(IRDO,IPLS)) GOTO 560
c  volumetric charge exchange rate
              WTRSIG=WTR*SIGVCX(IRCX)
              ADDV(IA+IATM,IRD)=ADDV(IA+IATM,IRD)+WTRSIG
c  note: in case of collision estimator for energy tallies in ircx process:
c        esigcx is not set (== 0.) in fpatha.f
c  volumetric incident ion energy loss rate due to charge exchange with iatm.
c  Gain for neutrals
              ADDV(IA+1*NATM+IATM,IRD)=ADDV(IA+1*NATM+IATM,IRD)+
     .                            WTRSIG*ESIGCX(IRCX,1)
c  volumetric incident neutral energy loss rate due to charge exchange
c  Loss for neutrals
              ADDV(IA+2*NATM+IATM,IRD)=ADDV(IA+2*NATM+IATM,IRD)+
     .                              WTRSIG*(-E0)
c  volumetric net ion energy loss rate due to charge exchange with iatm.
c  sign: for neutrals.
              ADDV(IA+3*NATM+IATM,IRD)=ADDV(IA+3*NATM+IATM,IRD)+
     .                              WTRSIG*(ESIGCX(IRCX,1)-E0)

ccc
ccc  next: parallel momentum exchange rates due to CX, ATOMS
ccc

             IF (LBXIN.AND.LBYIN.AND.LBZIN) THEN
               V0_PARB=VEL*
     .               (VELX*BXIN(IRDO)+VELY*BYIN(IRDO)+VELZ*BZIN(IRDO))
             ELSE
               V0_PARB=0._DP
             END IF
             PARMOM_0=V0_PARB*CNDYNA(IATM)

c   next: PERP NEUTRAL MOMENTUM = neutral momentum - par-neutral momentum

             V0_PERP=VEL*
     .               (VELX*(1.0-BXIN(IRDO))+
     .                VELY*(1.0-BYIN(IRDO))+
     .                VELZ*(1.0-BZIN(IRDO)))
             PERPMOM_0=V0_PERP*CNDYNA(IATM)

C  CHARGE EXCHANGE CONTRIBUTION FROM ATOMS
C

C
C  PRESENTLY: PARALLEL COMPONENT OF VSIGCX(IRCX) IS NOT AVAILABLE
C             FROM FUNCTION FPATHA
C  this uptusr: single species, ipls=iatm=1=natm=npls
C
            WTRSIG=WTR*SIGVCX(IRCX)

c  tally addv(ia+4*natm+...): currently free.

C  PREVIOUS INCIDENT BULK ION IPLS, NOW LOST FOR BULK
            ADDV(IA+5*NATM+IPLS,IRD)=ADDV(IA+5*NATM+IPLS,IRD)-
     .                               WTRSIG*VSIG_PARB(IPLS)
c  volumetric incident neutral par. momentum rate due to charge exchange
c  loss for neutrals, GAIN FOR IPLS OR VICE VERSA, DEPENDS ON FLOW AND NEUTRAL VELOCITY DIRECTION.
            ADDV(IA+6*NATM+IPLS,IRD)=ADDV(IA+6*NATM+IPLS,IRD)+
     .                      WTRSIG*PARMOM_0*
     .                      SIGN(1._DP,VAL_PARB(IPLS))
c  volumetric net ion parallel momentum loss/gain rate due to charge exchange
            ADDV(IA+7*NATM+IPLS,IRD)=ADDV(IA+7*NATM+IPLS,IRD)+
     .                      WTRSIG*(-VSIG_PARB(IPLS)+PARMOM_0*
     .                      SIGN(1._DP,VAL_PARB(IPLS)))

  560       CONTINUE
  590     CONTINUE

C  RADIAL GESCHWINDIGKEITSKOMPONENTE (CM/SEC)
          VR=(VELX*VPX(IRDO)+VELY*VPY(IRDO))*VEL
          ADDV(IA+8*NATM+IATM,IRD)=ADDV(IA+8*NATM+IATM,IRD)+WTR*VR
          ADDV(IA+9*NATM+IATM,IRD)=ADDV(IA+9*NATM+IATM,IRD)+WTR*VR*E0
C  POLOIDALE GESCHWINDIGKEITSKOMPONENTE (CM/SEC)
          VP=(VELX*VRX(IRDO)+VELY*VRY(IRDO))*VEL
          ADDV(IA+10*NATM+IATM,IRD)=ADDV(IA+10*NATM+IATM,IRD)+WTR*VP
          ADDV(IA+11*NATM+IATM,IRD)=ADDV(IA+11*NATM+IATM,IRD)+WTR*VP*E0

  200 CONTINUE  ! NCOU LOOP
C
C
C
      RETURN

  999 CONTINUE


 1000 CONTINUE
C  MOLECULES, CX ENERGY
      IF (ITYP.NE.2) GOTO 1999

C  CHECK: STORAGE FOR AT LEAST 8 MORE ADDITIONAL TRACKLENGTH-ESTIMATED TALLIES?
!pb   IF (NADV.LT.MOD_ADDV+12*NATMI+8*NMOLI) THEN
      IF (NADV.LT.12*NATMI+8*NMOLI) THEN
        GOTO 9999
      ELSE
C  THIS ROUTINE: SCORE 4 ADDITIONAL TALLIES ADDV
c  initial species index for ADDV tally: nspan
        nti=NSPAN(ntala)+12*NATM
        nte=Nti+8*NMOL-1
        LMETSP(nti:nte)=.TRUE.
      ENDIF
C
      DO 1200 ICOU=1,NCOU
          DIST=CLPD(ICOU)
          WTR=WV*DIST
          IRDO=NRCELL+NUPC(ICOU)*NR1P2+NBLCKA
          IRD=NCLTAL(IRDO)
c  assume here: bvin, parmom are set in plasma_deriv.
c               In case of other options (indpro): see update.f
cdr only signum needed: default tbd: signum=1
          IF (LBVIN) THEN
            VAL_PARB(1:NPLSI) = BVIN(MPLSV(1:NPLSI),IRDO)
          ELSE
            VAL_PARB(1:NPLSI) = 0._DP
          END IF
          IF (LPARMOM) THEN
            VSIG_PARB(1:NPLSI) = PARMOM(1:NPLSI,IRDO)
          ELSE
            VSIG_PARB(1:NPLSI) = 0._DP
          END IF
C
          IF (LGVAC(IRDO,0)) GOTO 1200
C
          if (ncou.gt.1) then
            XSTOR(:,:) = XSTOR2(:,:,ICOU)
            XSTORV(:)  = XSTORV2(:,ICOU)
          endif

C
C
          IF (LGMCX(IMOL,0,0).EQ.0) GOTO 1590
            DO 1560 IMCX=1,NMCXI(IMOL)
!pb           IA=MOD_ADDV+12*NATM  !  increment for addv tally 1st index: cx energy, molecules
              IA = 12*NATM
              IRCX=LGMCX(IMOL,IMCX,0)
              IPLS=LGMCX(IMOL,IMCX,1)
              IF (LGVAC(IRDO,IPLS)) GOTO 1560
c  volumetric charge exchange rate
              WTRSIG=WTR*SIGVCX(IRCX)
              ADDV(IA+IMOL,IRD)=ADDV(IA+IMOL,IRD)+WTRSIG
c  note: in case of collision estimator for energy tallies in ircx process:
c        esigcx is not set (== 0.) in fpatha.f
c  volumetric incident ion energy loss rate due to charge exchange
              ADDV(IA+1*NMOL+IMOL,IRD)=ADDV(IA+1*NMOL+IMOL,IRD)+
     .                            WTRSIG*ESIGCX(IRCX,1)
c  volumetric incident neutral energy loss rate due to charge exchange
              ADDV(IA+2*NMOL+IMOL,IRD)=ADDV(IA+2*NMOL+IMOL,IRD)+
     .                              WTRSIG*(-E0)
c  volumetric net ion energy loss rate due to charge exchange
              ADDV(IA+3*NMOL+IMOL,IRD)=ADDV(IA+3*NMOL+IMOL,IRD)+
     .                              WTRSIG*(ESIGCX(IRCX,1)-E0)
ccc
ccc  next: parallel momentum exchange rates due to CX
ccc

 1560       CONTINUE
 1590     CONTINUE

C  RADIAL GESCHWINDIGKEITSKOMPONENTE (CM/SEC)
          VR=(VELX*VPX(IRDO)+VELY*VPY(IRDO))*VEL
          ADDV(IA+4*NMOL+IMOL,IRD)=ADDV(IA+4*NMOL+IMOL,IRD)+WTR*VR
          ADDV(IA+5*NMOL+IMOL,IRD)=ADDV(IA+5*NMOL+IMOL,IRD)+WTR*VR*E0
C  POLOIDALE GESCHWINDIGKEITSKOMPONENTE (CM/SEC)
          VP=(VELX*VRX(IRDO)+VELY*VRY(IRDO))*VEL
          ADDV(IA+6*NMOL+IMOL,IRD)=ADDV(IA+6*NMOL+IMOL,IRD)+WTR*VP
          ADDV(IA+7*NMOL+IMOL,IRD)=ADDV(IA+7*NMOL+IMOL,IRD)+WTR*VP*E0

 1200 CONTINUE
C
C
C
      RETURN

 1999 CONTINUE

 2000 CONTINUE
C  MOLECULES, EL ENERGY
      IF (ITYP.NE.2) GOTO 2999

C  CHECK: STORAGE FOR AT LEAST 8 MORE ADDITIONAL TRACKLENGTH-ESTIMATED TALLIES?
!pb   IF (NADV.LT.MOD_ADDV+12*NATMI+12*NMOLI) THEN
      IF (NADV.LT.12*NATMI+12*NMOLI) THEN
        GOTO 9999
      ELSE
C  THIS ROUTINE: SCORE 4 ADDITIONAL TALLIES ADDV
c  initial species index for ADDV tally: nspan
        nti=NSPAN(ntala)+12*NATM+8*NMOL
        nte=Nti+4*NMOL-1
        LMETSP(nti:nte)=.TRUE.
      ENDIF
C
      DO 2200 ICOU=1,NCOU
          DIST=CLPD(ICOU)
          WTR=WV*DIST
          IRDO=NRCELL+NUPC(ICOU)*NR1P2+NBLCKA
          IRD=NCLTAL(IRDO)
C
          IF (LGVAC(IRDO,0)) GOTO 2200
C
          if (ncou.gt.1) then
            XSTOR(:,:) = XSTOR2(:,:,ICOU)
            XSTORV(:)  = XSTORV2(:,ICOU)
          endif

C
C
          IF (LGMEL(IMOL,0,0).EQ.0) GOTO 2590
            DO 2560 IMEL=1,NMELI(IMOL)
!pb           IA=MOD_ADDV+12*NATM+8*NMOL  !  increment for addv tally 1st index: el energy, molecules
              IA=12*NATM+8*NMOL
              IREL=LGMEL(IMOL,IMEL,0)
              IPLS=LGMEL(IMOL,IMEL,1)
              IF (LGVAC(IRDO,IPLS)) GOTO 2560
c  volumetric charge exchange rate
              WTRSIG=WTR*SIGVEL(IREL)
              ADDV(IA+IMOL,IRD)=ADDV(IA+IMOL,IRD)+WTRSIG
c  note: in case of collision estimator for energy tallies in irel process:
c        esigel is not set (== 0.) in fpathm.f
c  volumetric incident ion energy loss rate due to charge exchange
              ADDV(IA+1*NMOL+IMOL,IRD)=ADDV(IA+1*NMOL+IMOL,IRD)+
     .                            WTRSIG*ESIGEL(IREL,1)
c  volumetric incident neutral energy loss rate due to charge exchange
              ADDV(IA+2*NMOL+IMOL,IRD)=ADDV(IA+2*NMOL+IMOL,IRD)+
     .                              WTRSIG*(-E0)
c  volumetric net ion energy loss rate due to charge exchange
              ADDV(IA+3*NMOL+IMOL,IRD)=ADDV(IA+3*NMOL+IMOL,IRD)+
     .                              WTRSIG*(ESIGEL(IREL,1)-E0)
ccc
ccc  next: parallel momentum exchange rates due to EL
ccc

 2560       CONTINUE
 2590     CONTINUE
 2200 CONTINUE
C
C
C
      RETURN

 2999 CONTINUE



 9999 CONTINUE

C     WRITE (IUNOUT,*) 'NOTHING DONE IN UPTUSR'
      RETURN
      END SUBROUTINE EIRENE_UPTUSR

      END MODULE EIRMOD_UPTUSR
