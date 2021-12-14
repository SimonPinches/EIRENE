      MODULE EIRMOD_VELOCX
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
      USE EIRMOD_RANF, ONLY: RANF_EIRENE
      IMPLICIT NONE
      PRIVATE

      PUBLIC :: EIRENE_VELOCX, EIRENE_VELOCX_REINIT

      INTEGER, SAVE :: IFIRST = 0

      REAL(DP) :: VXN, VYN, VZN, VX,VY,VZ, VN, VXI, VYI, VZI,
     .          ZARGX, ZARGY, ZARGZ,
     .          VXDR, VYDR, VZDR, VRELQ,
     .          TEST, VREL, ELAB, CXS,
     .          VR, VRQ, ELMAX, ELMIN
C      REAL(DP) :: ELB
ctk      REAL(DP), EXTERNAL :: RANF_EIRENE

      INTEGER :: ICOUNT, J, JJ, IRL, IREAC

      SAVE

!$omp threadprivate(ifirst,vxn,vyn,vzn,vx,vy,vz,vn,vxi,vyi,vzi,
!$OMP& zargx,zargy,zargz,vxdr,vydr,vzdr,vrelq,test,vrel,elab,
!$OMP& cxs,vr,vrq,elmax,elmin,icount,j,jj,irl,ireac)


      CONTAINS

!pb  100107: SUBROUTINE VELOCX_REINIT added for reinitialization of EIRENE
!pb  110311: avoid relative velocity VREL=0
!pb  110311: ensure ELMIN <= ELAB <= ELMAX
!DR  250311: ensure ELMIN <= ELAB <= ELMAX disabled again: would lead
!DR          to wrong cross-sections, e.g. for beam penetration
CDR  5.8.15: ARGUMENTS ADDED TO VECUSR
cdr  aug.16: some test output, re asymptotic, rejection sampling. Commented out.
cdr  sep.17: sync with veloel. Prepare bgk relaxation. perhaps ready: nflag=2
cdr  jan.18: comments, cleanup. Sync with veloel, velopi, for incident ion sampling
cdr          then here: only relaxation, Delta_E=0. Scattering angle= Pi in COM.
cdr          but exchange of masses also allowed (distrinct from EL processes).
C
      SUBROUTINE EIRENE_VELOCX(K,VXO,VYO,VZO,VLO,IOLD,NOLD,VELQ,NFLAG,
     .                  IRCX,DUMT,DUMV)
C
C  THIS SUBROUTINE CARRIES OUT A CHARGE EXCHANGE COLLISION OF A TEST PARTICLE
C  WITH A BULK PARTICLE.
C  IT RETURNS THE POST-COLLISION VELOCITY VECTOR.
C
C  NFLAG= 1:       SAMPLING FROM MONOENERGETIC DISTRIBUTION
C                  OF ION SPEED IN 3D, X,Y,Z DIRECTION
C                  (I.E., DELTA FUNCTION IN ENERGY SPACE)
C                  E=M/2 V_M^2 =3/2 KT, IN REST FRAME OF IPLS
C                  USE WEIGHT CORRECTION OR REJECTION
C                  to be generalized to E=ESIGCX(IRCX,1)
C  NFLAG= 2:       SAMPLING FROM SHIFTED MAXWELLIAN
C                  "FMAXW" AT TI AND V-DRIFT IN CELL K
C  NFLAG= 3:       SAMPLING FROM SHIFTED MAXWELLIAN + WEIGHT CORRECTION
C                  FACTOR = SIGMA*VREL*FMAXW/<SIGMA*VREL>
C                  OR ALTERNATIVELY: REJECTION
C
C  K   : CELL INDEX

C  K   : .NE.0 :CELL INDEX FOR LOCAL BULK ION TI AND V_DRIFT
C  note: Ti has already been converted into thermal velocity units: zrg(ipls,k) in [cm/s]

C  K   : .EQ.0 :TX,TY,TZ,V-DRIFT_X,Y,Z ARE NOT FROM LOCAL BULK ION
C               SPECIES IPLS PARAMETERS, BUT EXPLICITLY DEFINED IN THE
C               PARAMETERS DUMT AND DUMV, RESPECTIVELY.
c  note: here dumt must also be in thermal velocity units

C  VXO : X COMPONENT OF SPEED UNIT VECTOR OF TEST PARTICLE BEFORE EVENT
C  VYO : Y COMPONENT OF SPEED UNIT VECTOR OF TEST PARTICLE BEFORE EVENT
C  VZO : Z COMPONENT OF SPEED UNIT VECTOR OF TEST PARTICLE BEFORE EVENT
C  VLO : VELOCITY OF TEST PARTICLE BEFORE EVENT
C  IOLD: SPECIES INDEX OF THE TEST PARTICLE BEFORE THE EVENT
C  NOLD: DITTO, IN MODCOL ARRAY
C  IPLS: SPECIES INDEX FOR THE THERMAL PLASMA ION VELOCITY
C        AND FOR THE PLASMA DRIFT VELOCITY TO BE USED AS
C        SHIFT VECTOR (IPLS IN COMMON COMUSR)
C  IRCX: LABEL FOR CX-REACTION, E.G., FOR SIGVCX(IRCX)
C        NOT NEEDED FOR NFLAG=2, THEN SET E.G.: IRCX=1


C  OR TO  FETCH A NEW  VELOCITY FOR A NEUTRAL ATOM "IATM",
C  A NEUTRAL MOLECULE "IMOL" OR A TEST ION "IION"
C  AFTER CX-EVENT WITH BULK ION "IPLS" IN CELL NO. K
C  FROM A SHIFTED
C  MAXWELLIAN (NFLAG=2), WEIGHTED BY SIGMA*VREL (NFLAG=3)

C  ADDITIONALLY:
C  USED E.G. FOR VOLUME RECOMBINATION SOURCE (NFLAG=2)
C      
      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: DUMT(3), DUMV(3)
      REAL(DP), INTENT(IN) :: VXO, VYO, VZO, VLO
      REAL(DP), INTENT(OUT) :: VELQ
      INTEGER, INTENT(IN) :: K, IOLD, NOLD, NFLAG, IRCX
      REAL(DP) ::EIRENE_CROSS

cpg      SAVE
C
c initialize arrays for "on the fly" rejection efficiency estimates
C IFLAG=1 AND IFLAG=3 OPTIONS
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
C  IS CROSS-SECTION AVAILABLE?
        IREAC=MODCOL(3,1,IRCX)
        IF (IREAC.EQ.0) GOTO 1
C CURRENTLY: HARD-WIRED SEARCH RANGE
        elmin=log(0.1_dp)
        elmax=log(1.e4_dp)
        SGCVMX(IRCX)=-1.D60
        JJ=1
        do j=1,1000
c  elab:  here ln(E), with E from 0.1 to 1e4 eV
          elab=elmin+(j-1)/999._dp*(elmax-elmin)

c  find cross-section at ENERGY ELAB from a fit or table.
          CXS=EIRENE_CROSS(ELAB,IREAC,IRCX,FACRCX(IRCX,1),'VELOCX 1')
c
          vrq=exp(elab-defCX(IRCX))
          vr=sqrt(vrq)
          if (cxs*vr.gt.SGCVMX(IRCX)) then
            JJ=J
            SGCVMX(IRCX)=cxs*vr
          endif
        enddo

!$OMP CRITICAL
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
!$OMP END CRITICAL

      ENDIF
    1 CONTINUE

c  preparations for process IRCX done. Start sampling procedure here.
C
C  INITIALIZE COUNTER FOR REJECTION SAMPLING OF INCIDENT BULK PARTICLE
C
      ICOUNT=1
C
C  NEXT: STEP 1
C
C    set parameters for random sampling in cell icell=K
C
      IF (K.GT.0.AND.K.LE.NRAD) THEN  ! K is the grid cell number. Use local bulk medium parameters
c  scaled 1d temperatures, per degree of freedom
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
      ELSEIF (K.EQ.0) THEN  !  K=0, USE ARGUMENTS DUMT AND DUMV AS PARAMETERS FOR DRIFTING MAXWELLIAN
        IF (NFLAG.NE.2) GOTO 999
        ZARGX=DUMT(1)
        ZARGY=DUMT(2)
        ZARGZ=DUMT(3)
        VXDR=DUMV(1)
        VYDR=DUMV(2)
        VZDR=DUMV(3)
      ELSE
        GOTO 999
      ENDIF
C
C
C  SAVE VELOCITY VECTOR (CM/S) OF INCIDENT TEST PARTICLE
      VX=VXO*VLO
      VY=VYO*VLO
      VZ=VZO*VLO
C

c   start random sampling here

  123 CONTINUE
      IF (INIV2.LE.0) CALL EIRENE_FGAUSS
C
C  SAMPLE FROM 3D NORMALIZED MAXWELLIAN (m=0;s=1)
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
C  ALL OTHER CASES: MAXWELLIAN AT LOCAL TEMPERATURE TIIN AND DRIFT VDR
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
        VXI=VXN   ! INCIDENT ION VELOCITY; CM/S.
        VYI=VYN
        VZI=VZN

C   NOTHING MORE TO BE DONE

C
      ELSE  !   NFLAG.NE.2, ALL OTHER OPTIONS
C
C   ALL OTHER DISTRIBUTIONS
C
C   WEIGHT CORRECTION DUE TO ENERGY DEPENDENCE IN CROSS-SECTION
C   OR: REJECTION     DUE TO ENERGY DEPENDENCE IN CROSS-SECTION
C   PRESENT VERSION: REJECTION
        VRELQ=MAX((VXN-VX)**2+(VYN-VY)**2+(VZN-VZ)**2, EPS30)
        VREL=SQRT(VRELQ)
        ELAB=LOG(VRELQ)+DEFCX(IRCX)
        IREAC=MODCOL(3,1,IRCX)
        CXS=EIRENE_CROSS(ELAB,IREAC,IRCX,FACRCX(IRCX,1),'VELOCX 2')
C
c...........................................................
cdr  test output only
c       elb=exp(elab)
c       if (elb.le.0.1) then
c         write (iunout,*) 'elb velocx ',elab,elb
c       endif
cdr
c.....................................................................

C
C       IF (NLREJC) THEN    !  REJECTION IS NOW DEFAULT OPTION
C
        IF (IFLRCX(IRCX).GT.0) THEN
          TEST=RANF_EIRENE()*SGCVMX(IRCX)
          IF (TEST.GT.CXS*VREL) THEN
C  REJECT
            ICOUNT=ICOUNT+1
            IF (ICOUNT.LT.500) GOTO 123  ! fetch a new bulk ion velocity
c  rejection loop failed, too many attempts.
            WRITE (iunout,*)
     .        'ICOUNT TOO LARGE ( > 500) IN VELOCX. ACCEPT SAMPLE '
cdr............................................................
cdr  test output only
cdr         ELB=EXP(ELAB)
cdr         write (iunout,*) 'npanu, ireac, ircx, ELAB(EV),icell ',
cdr  .                        npanu, ireac, ircx, ELB,  K
cdr............................................................
          ELSE
C  ACCEPT
            XCMEAN(IRCX)=XCMEAN(IRCX)+ICOUNT
            NCMEAN(IRCX)=NCMEAN(IRCX)+1
          ENDIF
C       ELSEIF (NLWEIGHT) THEN

        ELSE
C  FOR SOME REASON SGCVMX COULD NOT BE FOUND, or rejection is too inefficient.
C  SO USE WEIGHTING RATHER THAN REJECTION
          WEIGHT=WEIGHT*CXS*VREL*DIIN(IPLS,K)/SIGVCX(IRCX)
        ENDIF
C
        VXI=VXN
        VYI=VYN
        VZI=VZN

C
      ENDIF

C   CX = EXCHANGE OF IDENTITY (RELAXATION). NOTHING MORE TO BE DONE

      VELQ=VXI*VXI+VYI*VYI+VZI*VZI
      VEL=SQRT(VELQ)
      VN=1./VEL
      VELX=VXI*VN
      VELY=VYI*VN
      VELZ=VZI*VN
C
      RETURN
C
  999 CONTINUE
      WRITE (iunout,*)
     .  'PARAMETER ERROR IN SUBR. VELOCX. EXIT CALLED'
      CALL EIRENE_EXIT_OWN(1)
      END SUBROUTINE EIRENE_VELOCX

C  the following SUBROUTINE is for reinitialization of EIRENE
      SUBROUTINE EIRENE_VELOCX_REINIT
      IMPLICIT NONE
      IFIRST = 0
      return
      END SUBROUTINE EIRENE_VELOCX_REINIT

      END MODULE EIRMOD_VELOCX
