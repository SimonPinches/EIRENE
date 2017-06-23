!pb  100107: ENTRY VELOCX_REINIT added for reinitialization of EIRENE
!pb  110311: avoid relative velocity VREL=0
!pb  110311: ensure ELMIN <= ELAB <= ELMAX
!DR  250311: ensure ELMIN <= ELAB <= ELMAX disabled again: would lead
!DR          to wrong cross sections, e.g. for beam penetration
CDR  5.8.15: ARGUMENTS ADDED TO VECUSR
cdr  aug.16: some test output, re asymptotic, rejection sampling. commented out.
C
      SUBROUTINE EIRENE_VELOCX(K,VXO,VYO,VZO,VLO,IOLD,NOLD,VELQ,NFLAG,
     .                  IRCX,DUMT,DUMV)
C
C  THIS SUBROUTINE CARRIES OUT A CHARGE EXCHANGE COLLISION OF A TEST PARTICLE
C  WITH A BULK PARTICLE.
C  IT RETURNS THE POST COLLISION VELOCITY VECTOR.
C
C  NFLAG= 1:       SAMPLING FROM MONOENERGETIC DISTRIBUTION
C                  OF ION SPEED IN 3D, X,Y,Z DIRECTION
C                  (I.E., DELTA FUNCTION IN ENERGY SPACE)
C                  E=M/2 V_M^2 =3/2 KT
C                  to be generalized to E=ESIGCX(IRCX,1)
C  NFLAG= 2:       SAMPLING FROM SHIFTED MAXWELLIAN
C                  "FMAXW" AT TI AND V-DRIFT IN CELL K
C  NFLAG= 3:       SAMPLING FROM SHIFTED MAXWELLIAN + WEIGHT CORRECTION
C                  FACTOR = SIGMA*VREL*FMAXW/<SIGMA*VREL>
C                  OR ALTERNATIVELY: REJECTION
C
 
C  K   : .NE.0 :CELL INDEX FOR LOCAL BULK ION TI AND V_DRIFT
C  note: ti has already been converted into thermal velocity units: zrg(ipls,k)
 
C  K   : .EQ.0 :TX,TY,TZ,V-DRIFT_X,Y,Z ARE NOT FROM LOCAL BULK ION
C               SPECIES IPLS, BUT EXPLICITLY IN THE PARAMETERS DUMT AND DUMV.
C               RESPECTIVELY.
c  note: here dumt must also be in thermal velocity units
 
C  VXO : X COMPONENT OF SPEED UNIT VECTOR OF TEST PARTICLE BEFORE EVENT
C  VYO : Y COMPONENT OF SPEED UNIT VECTOR OF TEST PARTICLE BEFORE EVENT
C  VZO : Z COMPONENT OF SPEED UNIT VECTOR OF TEST PARTICLE BEFORE EVENT
C  VLO : VELOCITY OF TEST PARTICLE BEFORE EVENT
C  IOLD: SPECIES INDEX OF THE TEST PARTICLE BEFORE THE EVENT
C  NOLD: DITO, IN MODCOL-ARRAY
C  IPLS: SPECIES INDEX FOR THE THERMAL PLASMA ION VELOCITY
C        AND FOR THE PLASMA DRIFT VELOCITY TO BE USED AS
C        SHIFT VECTOR   (IPLS IN COMMON COMUSR)
C  IRCX: LABEL FOR CX-REACTION, E.G., FOR SIGVCX(IRCX)
C        NOT NEEDED FOR NFLAG=2, THEN SET E.G.: IRCX=1
 
C  USED E.G. FOR VOLUME RECOMBINATION SOURCE (NFLAG=2)
C  OR TO  FETCH A NEW  VELOCITY FOR A NEUTRAL ATOM "IATM",
C  A NEUTRAL MOLECULE "IMOL" OR A TEST ION "IION"
C  AFTER CX-EVENT WITH BULK ION "IPLS" IN CELL NO. K FROM A SHIFTED
C  MAXWELLIAN (NFLAG=2), WEIGHTED BY SIGMA*VREL (NFLAG=3)
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CRAND
      USE EIRMOD_CINIT
      USE EIRMOD_CZT1
      USE EIRMOD_CTRCEI
      USE EIRMOD_COMPRT
      USE EIRMOD_COMXS
      USE EIRMOD_CLAST
 
      IMPLICIT NONE
 
      REAL(DP), INTENT(IN) :: DUMT(3), DUMV(3)
      REAL(DP), INTENT(IN) :: VXO, VYO, VZO, VLO
      REAL(DP), INTENT(OUT) :: VELQ
      INTEGER, INTENT(IN) :: K, IOLD, NOLD, NFLAG, IRCX

      REAL(DP) :: VXN, VYN, VZN, VX,VY,VZ, VN, ZARGX, ZARGY, ZARGZ,
     .          VXDR, VYDR, VZDR, VRELQ, E0MAX, TIMAX, SIGS, VRELS,
     .          WRMEAN, TEST, VREL, WRAT, WO, ELAB, ELB, CXS,
     .          VR, VRQ, EIRENE_CROSS, ELMAX, ELMIN
      REAL(DP), EXTERNAL :: RANF_EIRENE
 
      INTEGER :: ICOUNT, J, JJ, IRL, IREAC
      INTEGER :: IFIRST = 0
 
      SAVE
C
      IF (IFIRST.EQ.0) THEN
        IFIRST=1
        DO IRL=1,NRCXI
          IFLRCX(IRL)=0
          NCMEAN(IRL)=0
          XCMEAN(IRL)=0.D0
        ENDDO
      ENDIF
C
      IF (IFLRCX(IRCX).EQ.0.AND.NFLAG.NE.2) THEN
        IFLRCX(IRCX)=-1
C  PREPARE REJECTION SAMPLING OF INCIDENT ION VELOCITY
C  IS CROSS SECTION AVAILABLE?
        IREAC=MODCOL(3,1,IRCX)
        IF (IREAC.EQ.0) GOTO 1
C CURRENTLY: HARD WIRED SEARCH RANGE
        elmin=log(0.1_dp)
        elmax=log(1.e4_dp)
        SGCVMX(IRCX)=-1.D60
        JJ=1
        do j=1,1000
c  elab:  here ln(E), with E from 0.1 to 1e4 eV
          elab=elmin+(j-1)/999._dp*(elmax-elmin)
          CXS=EIRENE_CROSS(ELAB,IREAC,IRCX,FACRCX(IRCX,1),'VELOCX 1')
          vrq=exp(elab-defCX(IRCX))
          vr=sqrt(vrq)
          if (cXS*vr.gt.SGCVMX(IRCX)) then
            JJ=J
            SGCVMX(IRCX)=cXS*vr
          endif
        enddo
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) 'FIRST CALL TO VELOCX FOR IRCX= ',IRCX
        WRITE (iunout,*) 'PREPARE REJECTION TECHNIQUE '
        WRITE (iunout,*) 'FIND MAX. "SGCVMX" OF SIGMA(VEL) * VEL '
        CALL EIRENE_MASJ1R('JJ, SGCVMX      ',JJ, SGCVMX(IRCX))
        IF (JJ.NE.1.AND.JJ.NE.1000) THEN
          elab=elmin+(JJ-1)/999.*(elmax-elmin)
          ELAB=EXP(ELAB)
          WRITE (iunout,*) 'TRUE MAXIMUM FOUND AT ELAB(EV) = ',ELAB
          IFLRCX(IRCX)=1
        ELSE
          WRITE (iunout,*) 'NO TRUE MAXIMUM FOUND, USE WEIGHTING '
        ENDIF
        CALL EIRENE_LEER(1)
      ENDIF
1     CONTINUE
C
C  INITIALIZE COUNTER FOR REJECTION SAMPLING OF INCIDENT BULK PARTICLE
C
      ICOUNT=1

C  NEXT: STEP 1
C
C    set parameters for random sampling in cell icell=K

C
      IF (K.GT.0) THEN  ! K is the grid cell number. Use local bulk medium parameters
c  scaled 1d temperatures, per degree of fredom
        ZARGX=ZRG(IPLS,K)
        ZARGY=ZRG(IPLS,K)
        ZARGZ=ZRG(IPLS,K)
c  drift velocity, cm/s
        IF (NLDRFT) THEN
          IF (INDPRO(4) == 8) THEN
            CALL EIRENE_VECUSR (2,K,X0,Y0,Z0,VXDR,VYDR,VZDR,IPLS,
     .                          .TRUE.)
          ELSE
            VXDR=VXIN(IPLS,K)
            VYDR=VYIN(IPLS,K)
            VZDR=VZIN(IPLS,K)
          END IF
        ELSE
          VXDR=0.D0
          VYDR=0.D0
          VZDR=0.D0
        ENDIF
      ELSE  !  K=0, USE ARGUMENTS DUMT AND DUMV AS PARAMETERS FOR DRIFTING MAXWELLIAN
        IF (NFLAG.NE.2) GOTO 999
        ZARGX=DUMT(1)
        ZARGY=DUMT(2)
        ZARGZ=DUMT(3)
        VXDR=DUMV(1)
        VYDR=DUMV(2)
        VZDR=DUMV(3)
      ENDIF
C

c   start random sampling here

123   CONTINUE
      IF (INIV2.LE.0) CALL EIRENE_FGAUSS
C
C  SAMPLE FROM 3D NORMALIZED MAXWELLIAN
      VXN=FG1(INIV2)
      VYN=FG2(INIV2)
      VZN=FG3(INIV2)
      INIV2=INIV2-1
C
      IF (NFLAG.EQ.1) THEN
C  DRIFTING, MONOENERGETIC ISOTROPIC DISTRIBUTION
C  ZT1 CORRESPONDS TO MEAN SQUARE VELOCITY AT TIIN(IPLS,K)
        VEL=SQRT(ZT1(IPLS,K))
        VN=VEL/SQRT(VXN*VXN+VYN*VYN+VZN*VZN)
        VXN=VXN*VN+VXDR
        VYN=VYN*VN+VYDR
        VZN=VZN*VN+VZDR
C  ALL OTHER CASES: MAXWELLIAN AT LOCAL TEMPERATURE AND DRIFT
      ELSE
        VXN=VXN*ZARGX+VXDR
        VYN=VYN*ZARGY+VYDR
        VZN=VZN*ZARGZ+VZDR
      ENDIF
C
C  DRIFTING MAXWELLIAN DISTRIBUTION (FOR MAXWELL-1/r^4-POTENTIAL: SIGMA*V = CONST.)
C
      IF (NFLAG.EQ.2) THEN
C

        VELQ=VXN*VXN+VYN*VYN+VZN*VZN
        VEL=SQRT(VELQ)
        VN=1./VEL
        VELX=VXN*VN
        VELY=VYN*VN
        VELZ=VZN*VN


C   NOTHING MORE TO BE DONE
C
        RETURN
C
      ELSE
C
C  SAVE  INCIDENT TEST PARTICLE VELOCITY
        VX=VXO*VLO
        VY=VYO*VLO
        VZ=VZO*VLO
C
C   ALL OTHER DISTRIBUTIONS
C
C   WEIGHT CORRECTION DUE TO ENERGY DEPENDENCE IN CROSS SECTION
C   OR: REJECTION     DUE TO ENERGY DEPENDENCE IN CROSS SECTION
C   PRESENT VERSION: REJECTION
        VRELQ=MAX((VXN-VX)**2+(VYN-VY)**2+(VZN-VZ)**2, EPS30)
        VREL=SQRT(VRELQ)
        ELAB=LOG(VRELQ)+DEFCX(IRCX)
        IREAC=MODCOL(3,1,IRCX)

cdr.........................................................  
cdr  test output only
c       elb=exp(elab)
c       if (elb.le.0.1) then
c         write (6,*) 'elb velocx ',elab,elb
c       endif
cdr.........................................................  

        CXS=EIRENE_CROSS(ELAB,IREAC,IRCX,FACRCX(IRCX,1),'VELOCX 2')
C
C       IF (NLREJC) THEN    !  REJECTION IS NOW DEFAULT OPTION
C
        IF (IFLRCX(IRCX).GT.0) THEN
          TEST=RANF_EIRENE()*SGCVMX(IRCX)
          if (test.gt.cxs*vrel) then
c  reject
            icount=icount+1
            if (icount.lt.500) goto 123  ! fetch a new bulk ion velocity
c  rejection loop failed, too many attempts.
            write (iunout,*)
     .        'icount too large ( > 500) IN VELOCX. ACCEPT SAMPLE '
cdr............................................................   
cdr  test output only
cdr         ELB=EXP(ELAB)
cdr         write (iunout,*) 'npanu, ireac, ircx, ELAB(EV),icell ',
cdr  .                        npanu, ireac, ircx, ELB,  K
cdr............................................................
          else
c  accept
            xcmean(ircx)=xcmean(ircx)+icount
            ncmean(ircx)=ncmean(ircx)+1
          endif
C       ELSEIF (NLWEIGHT) THEN
 
        ELSE
C  FOR SOME REASON SGCVMX COULD NOT BE FOUND.
C  SO USE WEIGHTING RATHER THAN REJECTION
          WEIGHT=WEIGHT*CXS*VREL*DIIN(IPLS,K)/SIGVCX(IRCX)
        ENDIF
C
C
        VELQ=VXN*VXN+VYN*VYN+VZN*VZN
        VEL=SQRT(VELQ)
        VN=1./VEL
        VELX=VXN*VN
        VELY=VYN*VN
        VELZ=VZN*VN
C

      ENDIF
C
      RETURN
C
999   CONTINUE
      WRITE (iunout,*)
     .  'PARAMETER ERROR IN SUBR. VELOCX. EXIT CALLED'
      CALL EIRENE_EXIT_OWN(1)
 
C  the following ENTRY is for reinitialization of EIRENE 
      ENTRY EIRENE_VELOCX_REINIT
      IFIRST = 0
      return
      END
