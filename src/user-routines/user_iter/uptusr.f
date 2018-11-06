! 23.08.06: VPX, VPY, VRX, VRY changed to ALLOCATABLE, SAVE to speed up
!           subroutine call (save time in storage allocation)
cdr Jan 17: remove local allocatable cndyn.. arrays. These are now
cdr         set in code initialisation phase
cdr may 18: revised, particle currents, particle flux,...., comments..
cdr         I am not sure that the rad and pol normal vectors are correct.
cdr         In solps5.0 we use the underlying polygon grid. 
C
C
      SUBROUTINE EIRENE_UPTUSR(XSTOR2,XSTORV2,WV,IFLAG)
C
C  USER-SUPPLIED TRACKLENGTH ESTIMATOR, VOLUME-AVERAGED
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
      USE EIRMOD_CTRIG
      use eirmod_extrab25

      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: XSTOR2(MSTOR1,MSTOR2,N2ND+N3RD),
     .                        XSTORV2(NSTORV,N2ND+N3RD), WV
      INTEGER, INTENT(IN) :: IFLAG
CDR
      REAL(DP), ALLOCATABLE, SAVE :: VPX(:),VPY(:),VRX(:),VRY(:)
CDR
      INTEGER :: IAT, IPL, I, IR, IP, IRD
      INTEGER, SAVE :: IFIRST, IA0, IA1, IA2, IA3, IA4
      integer :: icou
      real(dp) :: wtr,vr,vp,dist
      DATA IFIRST/0/
 
      IF (IFIRST.EQ.0) THEN
        IFIRST=1
C
CDR
CDR  PROVIDE A RADIAL UNIT VECTOR PER CELL
CDR  VPX,VPY NEEDED FOR PROJECTING PARTICLE VELOCITIES
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
        DO I=1,NTRII
            VPX(I)=PLNXTRI(i)
            VPY(I)=PLNYTRI(i)
            VRX(I)=PPLNXTRI(i)
            VRY(I)=PPLNYTRI(i)
        END DO

cdr  increments for tally number iadv
        IA0=0               !  RADIAL CURRENT
        IA1=NATMI+NMOLI     !  RADIAL ENREGY FLUX
        IA2=2*IA1           !  POLOIDAL CURRENT
        IA3=3*IA1           !  POLOIDAL ENERGY FLUX
        IA4=4*IA1           !  FLUX (ANGLE-AVERAGED)
      ENDIF

csw
csw added for B2.5 coupling (ank_mods in eirene_mc.F of B2.5)
csw so-called 'STANDARD' option in SOLPS4.3 user/uptusr.f from Vlad
csw 21oct2011
csw
C
C  WV=WEIGHT/VEL
C
C  ATOMS
      IF (ITYP.EQ.1) THEN
        DO 20 ICOU=1,NCOU

          if (ncou.gt.1) then
            XSTOR(:,:) = XSTOR2(:,:,ICOU)
            XSTORV(:)  = XSTORV2(:,ICOU)
          endif

          DIST=CLPD(ICOU)
          WTR=WV*DIST
          IRD=NRCELL+NUPC(ICOU)*NR1ST+NBLCKA
C
cdr       IF (LGVAC(IRD,0)) GOTO 20   ! score neutral fluxes also in Vac. region !
C
CDR
C  particle current, radial component  (CM/SEC)
          VR=(VELX*VPX(IRD)+VELY*VPY(IRD))*VEL
          if(ia0+iatm.gt.nadv) goto 20
          ADDV(IATM,IRD)=ADDV(IATM,IRD)+WTR*VR
          if(ia1+iatm.gt.nadv) goto 20
          ADDV(IA1+IATM,IRD)=ADDV(IA1+IATM,IRD)+WTR*VR*E0
C  particle current, poloidal component (CM/SEC)
          VP=(VELX*VRX(IRD)+VELY*VRY(IRD))*VEL
          if(ia2+iatm.gt.nadv) goto 20
          ADDV(IA2+IATM,IRD)=ADDV(IA2+IATM,IRD)+WTR*VP
          if(ia3+iatm.gt.nadv) goto 20
          ADDV(IA3+IATM,IRD)=ADDV(IA3+IATM,IRD)+WTR*VP*E0

c  particle current: toroidal component (cm/sec)  
c    can be found from default tallies: vden_xzy scalar product bxin,....bzin
c    note   particle current, cartesian, vden_xyz is now a default tally.

C  particle flux, integrated over all directions
          if(ia4+iatm.gt.nadv) goto 20
          ADDV(IA4+IATM,IRD)=ADDV(IA4+IATM,IRD)+WTR*VEL
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
C  particle current, radial component  (CM/SEC)
          VR=(VELX*VPX(IRD)+VELY*VPY(IRD))*VEL
          if(ia0+natmi+imol.gt.nadv) goto 200
          ADDV(NATMI+IMOL,IRD)=ADDV(NATMI+IMOL,IRD)+WTR*VR
          if(ia1+natmi+imol.gt.nadv) goto 200
          ADDV(IA1+NATMI+IMOL,IRD)=ADDV(IA1+NATMI+IMOL,IRD)+WTR*VR*E0
C  particle current, poloidal component (CM/SEC)
          VP=(VELX*VRX(IRD)+VELY*VRY(IRD))*VEL
          if(ia2+natmi+imol.gt.nadv) goto 200
          ADDV(IA2+NATMI+IMOL,IRD)=ADDV(IA2+NATMI+IMOL,IRD)+WTR*VP
          if(ia3+natmi+imol.gt.nadv) goto 200
          ADDV(IA3+NATMI+IMOL,IRD)=ADDV(IA3+NATMI+IMOL,IRD)+WTR*VP*E0

c  particle current: toroidal component (cm/sec)  
c    can be found from default tallies: vden_xzy scalar product bxin,....bzin
c    note   particle current, cartesian, vden_xyz is now a default tally.

c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C  particle flux, integrated over all directions

          if(ia4+NATMI+IMOL.gt.nadv) goto 200
          ADDV(IA4+NATMI+IMOL,IRD)=ADDV(IA4+NATMI+IMOL,IRD)+WTR*VEL
C
200     CONTINUE
CDR
C
C  TEST IONS
      ELSEIF (ITYP.EQ.3) THEN
C  TO BE WRITTEN
C
      ENDIF
C
      RETURN
      END
