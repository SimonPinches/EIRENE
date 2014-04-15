! 23.08.06: VPX, VPY, VRX, VRY changed to ALLOCATABLE, SAVE to speed up
!           subroutine call (save time in storage allocation)
C
C
      SUBROUTINE EIRENE_UPTUSR(XSTOR2,XSTORV2,WV,IFLAG)
C
C  USER SUPPLIED TRACKLENGTH ESTIMATOR, VOLUME AVERAGED
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CESTIM
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT
      USE EIRMOD_CUPD
      USE EIRMOD_COMXS
      USE EIRMOD_CSPEZ
      USE EIRMOD_CGRID
      USE EIRMOD_CLOGAU
      USE EIRMOD_CCONA
      USE EIRMOD_CPOLYG
      USE EIRMOD_CZT1
      IMPLICIT NONE
      REAL(DP), INTENT(INOUT) :: XSTOR2(MSTOR1,MSTOR2,N2ND+N3RD),
     .                         XSTORV2(NSTORV,N2ND+N3RD), WV
      INTEGER, INTENT(IN) :: IFLAG
      REAL(DP), ALLOCATABLE, SAVE :: CNDYNA(:),CNDYNP(:)
CDR
      REAL(DP), ALLOCATABLE, SAVE :: VPX(:),VPY(:),VRX(:),VRY(:)
CDR
      INTEGER :: IAT, IPL, I, IR, IP, IRD
      INTEGER, SAVE :: IFIRST, IA1, IA2, IA3, NA4, INDEXM, INDEXF
      integer :: icou
      real(dp) :: wtr,vr,vp,dist
      DATA IFIRST/0/
 
      IF (IFIRST.EQ.0) THEN
        IFIRST=1
        ALLOCATE (CNDYNA(NATM))
        ALLOCATE (CNDYNP(NPLS))
        DO IAT=1,NATMI
          CNDYNA(IAT)=1.D3*AMUA*RMASSA(IAT)
        END DO
        DO IPL=1,NPLSI
          CNDYNP(IPL)=1.D3*AMUA*RMASSP(IPL)
        END DO
C
CDR
CDR  PROVIDE A RADIAL UNIT VECTOR PER CELL
CDR  VPX,VPY,  NEEDED FOR PROJECTING PARTICLE VELOCITIES
CDR  SAME FOR POLOIDAL UNIT VECTOR VRX,VRY
C
        ALLOCATE (VPX(NRAD))
        ALLOCATE (VPY(NRAD))
        ALLOCATE (VRX(NRAD))
        ALLOCATE (VRY(NRAD))
        DO I=1,NRAD
          VPX(I)=0.
          VPY(I)=0.
          VRX(I)=0.
          VRY(I)=0.
        END DO
        DO IR=1,NR1STM
          DO IP=1,NP2NDM
            IRD=IR+(IP-1)*NR1P2
            VPX(IRD)=PLNX(IR,IP)
            VPY(IRD)=PLNY(IR,IP)
            VRX(IRD)=PPLNX(IR,IP)
            VRY(IRD)=PPLNY(IR,IP)
          END DO
        END DO
        IA1=NATMI+NMOLI
        IA2=2*IA1
        IA3=3*IA1
        NA4=4*IA1
        INDEXM=NPLSI
        INDEXF=2*NPLSI
      ENDIF
csw
csw added for B2.5 coupling (ank_mods in eirene_mc.F of B2.5)
csw so called 'STANDARD' option in Vlad's SOLPS4.3 user/uptusr.f
csw 21oct2011
csw
C
C  WV=WEIGHT/VEL
C
C  ATOMS
      IF (ITYP.EQ.1) THEN
        DO 20 ICOU=1,NCOU
          DIST=CLPD(ICOU)
          WTR=WV*DIST
          IRD=NRCELL+NUPC(ICOU)*NR1ST+NBLCKA
C
          IF (LGVAC(IRD,0)) GOTO 20
C
CDR
C  RADIAL GESCHWINDIGKEITSKOMPONENTE  (CM/SEC)
          VR=(VELX*VPX(IRD)+VELY*VPY(IRD))*VEL
          if(iatm.gt.nadv) goto 20
          ADDV(IATM,IRD)=ADDV(IATM,IRD)+WTR*VR
          if(ia1+iatm.gt.nadv) goto 20
          ADDV(IA1+IATM,IRD)=ADDV(IA1+IATM,IRD)+WTR*VR*E0
C  POLOIDALE GESCHWINDIGKEITSKOMPONENTE (CM/SEC)
          VP=(VELX*VRX(IRD)+VELY*VRY(IRD))*VEL
          if(ia2+iatm.gt.nadv) goto 20
          ADDV(IA2+IATM,IRD)=ADDV(IA2+IATM,IRD)+WTR*VP
          if(ia3+iatm.gt.nadv) goto 20
          ADDV(IA3+IATM,IRD)=ADDV(IA3+IATM,IRD)+WTR*VP*E0
CDR
20      CONTINUE
C
C  MOLECULES
      ELSEIF (ITYP.EQ.2) THEN
CDR
        DO 200 ICOU=1,NCOU
          DIST=CLPD(ICOU)
          WTR=WV*DIST
          IRD=NRCELL+NUPC(ICOU)*NR1ST+NBLCKA
C
          IF (LGVAC(IRD,0)) GOTO 200
C  RADIAL GESCHWINDIGKEITSKOMPONENTE (CM/SEC)
          VR=(VELX*VPX(IRD)+VELY*VPY(IRD))*VEL
          if(natmi+imol.gt.nadv) goto 200
          ADDV(NATMI+IMOL,IRD)=ADDV(NATMI+IMOL,IRD)+WTR*VR
          if(ia1+natmi+imol.gt.nadv) goto 200
          ADDV(IA1+NATMI+IMOL,IRD)=ADDV(IA1+NATMI+IMOL,IRD)+WTR*VR*E0
C  POLOIDALE GESCHWINDIGKEITSKOMPONENTE (CM/SEC)
          VP=(VELX*VRX(IRD)+VELY*VRY(IRD))*VEL
          if(ia2+natmi+imol.gt.nadv) goto 200
          ADDV(IA2+NATMI+IMOL,IRD)=ADDV(IA2+NATMI+IMOL,IRD)+WTR*VP
          if(ia3+natmi+imol.gt.nadv) goto 200
          ADDV(IA3+NATMI+IMOL,IRD)=ADDV(IA3+NATMI+IMOL,IRD)+WTR*VP*E0
C
200     CONTINUE
CDR
C
C  TEST IONS
      ELSEIF (ITYP.EQ.3) THEN
C  TO BE WRITTEN
C
      ENDIF
      RETURN
      END
 
 
 
