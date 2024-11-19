      MODULE EIRMOD_VELOEL

cym  04/2020 Module now contains several routines initially in volume-processes
cym  and called in veloel only : gaumeh,rstern,rtsaf,fi,fivec
cym  gaumeh & rstern had SAVE attributes
cym  there was a common shared between fivec & gaumeh (/CFI/)
cym  once this is done, it makes more sense to include fi&rtsaf here too


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

      PUBLIC :: EIRENE_VELOEL, EIRENE_VELOEL_REINIT

      INTEGER, SAVE :: IFIRST = 0

      REAL(DP) :: TEST, VRELX,
     .          VXN, VYN, VZN, VXI, VYI, VZI, VRELQ, VN, VREL,
     .          RLMS, RMSI, RMN, RMI, VRELY, EPS, CPH, CHI,
     .          RESULT, SPH, VRSX, VRSY, VRSZ, CEPS,
     .          SEPS, RS,
     .          VSX, VSY, VRYZ, VRELZ, VRQYZ, VSZ, PH, CCHI,
     .          BMAX, ER, ELMIN, ELMAX, B,
     .          VXDR, VYDR, VZDR,
     .          ZARGX, ZARGY, ZARGZ,
     .          VX, VY, VZ, ELAB,
     .          VR, CEL, VRQ

      REAL(DP) :: P(9)
      INTEGER :: IFLAG, IREAC, JM, ICOUNT

#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP THREADPRIVATE(IFIRST,
!$OMP&  TEST, VRELX,
!$OMP&  VXN, VYN, VZN, VXI, VYI, VZI, VRELQ, VN, VREL,
!$OMP&  RLMS, RMSI, RMN, RMI, VRELY, EPS, CPH, CHI,
!$OMP&  RESULT, SPH, VRSX, VRSY, VRSZ, CEPS,
!$OMP&  SEPS, RS,
!$OMP&  VSX, VSY, VRYZ, VRELZ, VRQYZ, VSZ, PH, CCHI,
!$OMP&  BMAX, ER, ELMIN, ELMAX, B,
!$OMP&  VXDR, VYDR, VZDR,
!$OMP&  ZARGX, ZARGY, ZARGZ,
!$OMP&  VX, VY, VZ, ELAB,
!$OMP&  VR, CEL, VRQ,
!$OMP&  P,
!$OMP& IFLAG, IREAC, JM, ICOUNT)
#endif

cym ccccccccccccccccccccc variables from RSTERN ccccccccccccccccccccccccccc
cym this comes here because of SAVE attribute

      REAL(DP) :: FIW, RV,  BQ, TOL, DFI,
     .        FITEST, rup, rlw, vw, rw, r0

#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP THREADPRIVATE(FIW,RV,BQ,TOL,DFI,
!$OMP&      FITEST, rup, rlw, vw, rw, r0)
#endif

cym ccccccccccccccccccccc variables from GAUMEH ccccccccccccccccccccccccccc
cym this comes here because of SAVE attibute

cym common to be treated fivec/gaumeh
cym      REAL(DP) :: AR(128), AFI(128)
cym      INTEGER :: NFI
cym      COMMON /CFI/ AR,AFI,NFI
cym common to be treated

      REAL(DP) :: AR(128), AFI(128)
      INTEGER :: NFI

#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP THREADPRIVATE(AR,AFI,NFI)
#endif

      REAL(DP) :: X5A(5),X5B(5),X10A(10),X10B(10),X20A(20),X20B(20)
C
      REAL(DP) :: XG10A(10),XG10B(10)
      REAL(DP) :: SUM, F, X

cym      DATA IFIRST1 /0/
cym CHECK THIS
      integer :: ifirst1=0
C
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP THREADPRIVATE(X5A,X5B,X10A,X10B,X20A,X20B,
!$OMP&       XG10A,XG10B,
!$OMP&       SUM, F, X,
!$OMP&       IFIRST1)
#endif

cym      integer :: icll
cym !$OMP THREADPRIVATE(icll)

cym      SAVE

      CONTAINS

c         061205: au_to_cm2 --> ccona
!         100107: SUBROUTINE VELOEL_REINIT added for reinitialization of EIRENE
CDR  Aug. 2015  : PROGRAMMING AND NOTATION SYNCHRONIZED WITH VELOCX.F
!DR
!DR
CDR  5.8.15 : ARGUMENTS ADDED TO VECUSR
cdr  4.9.15 : npbgkp(ipls,1).ne.0, rather than nchrgp(ipls)=0 to identify bgk collisions
cdr           There may be also neutral atom -test ion  bgk collisions, in which case
cdr           the virtual bgk background species may be an ion ?
cdr           (not contributing to electron density, however).
cdr           example He  He+  collisions, when both are test species.
cdr           Or, e.g. elastic component in H + p
cdr  Jan. 17: Added option: isotropic in COM frame, iflag=0, when modcol(5,0,..)=0
cdr                         exchange of identity in lab system when modcol(5,0,..)=-1
cdr                         (this was default for bgk collisions so far, with OLD iflag=0)
cdr March 18: Slight changes in notation, to sync with a new proprietary version of this
cdr           routine which carries out cross-section integration directly
cdr           from interaction potentials without intermediate fits. veloel_TEST.ff
cdr Nov.  19: More clearly separate NFLAG, IFLAG options
C
      SUBROUTINE EIRENE_VELOEL(K,VXO,VYO,VZO,VLO,IOLD,NOLD,VELQ,NFLAG,
     .                         IREL,RMASS)
C
C  THIS SUBROUTINE CARRIES OUT AN ELASTIC COLLISION OF A TEST PARTICLE
C  WITH A BULK PARTICLE.
C  IT RETURNS THE POST-COLLISION VELOCITY VECTOR.
C
C  NFLAG= 1: SAMPLING FROM MONOENERGETIC DISTRIBUTION
C            OF ION SPEED IN 3D, X,Y,Z DIRECTION
C            (I.E., DELTA FUNCTION IN ENERGY SPACE)
C            E=M/2 V_M^2 =3/2 KT, IN REST FRAME OF IPLS
C            USE WEIGHT CORRECTION OR REJECTION
C            to be generalized to E=ESIGCX(IRCX,1)
C  NFLAG= 2: SAMPLING FROM SHIFTED MAXWELLIAN
C            "FMAXW" AT TI AND V-DRIFT IN CELL K
C  NFLAG= 3: SAMPLING FROM SHIFTED MAXWELLIAN + WEIGHT CORRECTION
C            FACTOR = SIGMA*VREL*FMAXW/<SIGMA*VREL>
C            OR ALTERNATIVELY: REJECTION
C
C  1ST STEP: FIND COLLISION PARTNER FROM BULK ION SPECIES "IPLS":
C            (VXN,VYN,VZN)
C  2ND STEP: FIND CROSS-SECTION AS FUNCTION OF ELAB (AND
C            CARRY OUT WEIGHT CORRECTION, UNLESS REJECTION HAS BEEN USED
C            IN STEP 1
C  3RD STEP: FIND IMPACT PARAMETER B
C  4TH STEP: FIND NEW VELOCITY VECTOR
C
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
C  IREL: LABEL FOR EL-REACTION, E.G., FOR SIGVEL(IREL)
C        NOT NEEDED FOR NFLAG=2, THEN SET E.G.: IREL=1
C

cym moving modules & local variables to declaration section
cym this is because of SAVE attribute

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: RMASS
      REAL(DP), INTENT(IN) :: VXO, VYO, VZO, VLO
      REAL(DP), INTENT(OUT) :: VELQ
      INTEGER, INTENT(IN) :: K, IOLD, NOLD, NFLAG, IREL
      INTEGER :: IRL, J
      REAL(DP) :: EIRENE_CROSS
      EXTERNAL :: EIRENE_CROSS, EIRENE_FGAUSS, EIRENE_VECUSR,
     .            EIRENE_MASJ1R, EIRENE_MASJ4, EIRENE_MASJ5,
     .            EIRENE_LEER, EIRENE_EXIT_OWN

#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP THREADPRIVATE( IRL, J )
#endif

CYM/HJL Variables moved to module scope

!  PARAMETERS P_A_B FOR INTERACTION POTENTIALS for A on B collisions
!  ARE NOW READ FROM FILE AMJUEL,
!  NOT HARD-WIRED IN THIS ROUTINE OR (EVEN OLDER VERSIONS)
!  IN FUNTIONS FI, FIVEC
!     DATA P_HE_HE/2.55,2.35,0.90,1.9842,1.3990 ,2.6345  ,0.,-1.9125,0./
!     DATA P_H2_H /2.70,3.00,1.00,2.8355,2.18038,3.490687,0.,-2.025 ,0./
!     DATA P_HE_H /2.00,2.20,0.85,1.4556,0.99699,1.99515 ,0.,-1.5   ,0./
!     DATA P_NE_H /2.28,2.68,0.85,1.8709,1.3870 ,2.4402  ,0.,-1.71  ,0./
!     DATA P_AR_H /4.04,2.50,0.86,2.4756,1.7892 ,3.2737  ,0.,-3.03  ,0./
!     DATA P_KR_H /4.45,2.50,0.80,2.7779,2.0077 ,3.7406  ,0.,-3.3375,0./
!     DATA P_XE_H /6.75,3.80,1.08,3.2882,2.6884 ,3.8436  ,0.,-5.0625,0./
c  special relations for MORSE-Potential functions V(r):
c                iflag=2, projectile A on target B
c                these parameters P(5), P(6), P(8) are also now
c                read from file AMJUEL rather than initialising them here
c       P_A_B(5)=P(4)*(1.-LOG(2.)/P(2))         (=R0, ROOT OF V)
c       P_A_B(6)=P(4)*(1.+LOG(2.)/P(2)/P(3))    (=RW, INFLECTION OF V)
c       P_A_B(8)=-3.*P(1)/4.                    (=V(RW) )
      SAVE
C
c initialize arrays for "on the fly" rejection efficiency estimates
C NFLAG=1 AND NFLAG=3 OPTIONS
      IF (IFIRST.EQ.0) THEN
        IFIRST=1
        DO IRL=1,NRELI
          IFLREL(IRL)=0
          NEMEAN(IRL)=0
          XEMEAN(IRL)=0.D0
        ENDDO
      ENDIF
C
      IFLAG=-2
      IF (IFLREL(IREL).EQ.0.AND.NFLAG.NE.2) THEN
        IFLREL(IREL)=-1
C  PREPARE REJECTION SAMPLING OF INCIDENT ION VELOCITY
C  IS CROSS-SECTION AVAILABLE?
        IREAC=MODCOL(5,1,IREL)
        IF (IREAC.EQ.0) GOTO 1
C CURRENTLY: HARD-WIRED SEARCH RANGE
        elmin=log(0.01_dp)
        elmax=log(1.e3_dp)
        SGEVMX(IREL)=-1.D60
        JM=1
        do j=1,1000
c  elab: here ln(E), with E from 0.01 to 1e3 eV
          elab=elmin+(j-1)/999._dp*(elmax-elmin)

c  find cross-section at ENERGY ELAB from a fit or table.
          CEL=EIRENE_CROSS(ELAB,IREAC,IREL,FACREL(IREL,1),'VELOEL 1')
c
          vrq=exp(elab-defel(IREL))
          vr=sqrt(vrq)
          if (cel*vr.gt.SGEVMX(IREL)) then
            JM=J
            SGEVMX(IREL)=cel*vr
          endif
        enddo

        CALL EIRENE_LEER(1)
        WRITE (iunout,*) 'FIRST CALL TO VELOEL FOR IREL= ',IREL
        WRITE (iunout,*) 'PREPARE REJECTION TECHNIQUE '
        WRITE (iunout,*) 'FIND MAX. "SGCVMX" OF SIGMA(VEL) * VEL '
        CALL EIRENE_MASJ1R('JM, SGEVMX      ',JM, SGEVMX(IREL))
        IF (JM.NE.1.AND.JM.NE.1000) THEN
          elab=elmin+(JM-1)/999.*(elmax-elmin)
          ELAB=EXP(ELAB)
          WRITE (iunout,*) 'TRUE MAXIMUM FOUND AT ELAB(EV) = ',ELAB
          IFLREL(IREL)=1
        ELSE
          WRITE (iunout,*) 'NO TRUE MAXIMUM FOUND, USE WEIGHTING '
        ENDIF
        CALL EIRENE_LEER(1)
      ENDIF
    1 CONTINUE

c  preparations for process IREL done. Start sampling procedure here.
C
C  INITIALIZE COUNTER FOR REJECTION SAMPLING OF INCIDENT BULK PARTICLE
C
      ICOUNT=1
C
C  NEXT: STEP 1
C
C    set parameters for random sampling in cell icell=K
C
      IF (K.GT.0.AND.K.LE.NRAD) THEN  ! K is the grid cell number.
                                      ! Use local bulk medium parameters
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
C  DRIFTING, MONO-ENERGETIC ISOTROPIC DISTRIBUTION
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
C  DRIFTING MAXWELLIAN DISTRIBUTION (FOR MAXWELL-1/r^4-POTENTIAL:
C  SIGMA*V = CONST(T), BUT INDEPENDENT OF V)
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
        VRELQ=MAX((VXN-VX)**2+(VYN-VY)**2+(VZN-VZ)**2,EPS30)
        VREL=SQRT(VRELQ)
        ELAB=LOG(VRELQ)+DEFEL(IREL)
        IREAC=MODCOL(5,1,IREL)
        IF (IREAC.EQ.0) GOTO 995
        CEL=EIRENE_CROSS(ELAB,IREAC,IREL,FACREL(IREL,1),'VELOEL 2')
C
c.............................................................
cdr  test output only
c       elb=exp(elab)
c       if (elb.le.0.01) then
c         write (iunout,*) 'elb veloel ',elab,elb
c       endif
cdr
c.....................................................................

C
C       IF (NLREJC) THEN   !  REJECTION IS NOW DEFAULT OPTION
C
        IF (IFLREL(IREL).GT.0) THEN
          TEST=RANF_EIRENE()*SGEVMX(IREL)
          IF (TEST.GT.CEL*VREL) THEN
C  REJECT
            ICOUNT=ICOUNT+1
            IF (ICOUNT.LT.500) GOTO 123  ! fetch a new bulk ion velocity
c  rejection loop failed, too many attempts.
            WRITE (iunout,*)
     .        'ICOUNT TOO LARGE ( > 500) IN VELOEL. ACCEPT SAMPLE '
            WRITE (iunout,*) 'NPANU, IREAC, IREL, ELAB ',
     .                        NPANU, IREAC, IREL, ELAB
          ELSE
C  ACCEPT
            XEMEAN(IREL)=XEMEAN(IREL)+ICOUNT
            NEMEAN(IREL)=NEMEAN(IREL)+1
          ENDIF
C       ELSEIF (NLWEIGHT) THEN

        ELSE
C  FOR SOME REASON SGEVMX COULD NOT BE FOUND, or rejection is too inefficient.
C  SO USE WEIGHTING RATHER THAN REJECTION
          WEIGHT=WEIGHT*CEL*VREL*DIIN(IPLS,K)/SIGVEL(IREL)
        ENDIF
C
        VXI=VXN
        VYI=VYN
        VZI=VZN

C
      ENDIF
C
C  STEP 1 FINISHED, INCIDENT BULK ION (IPLS) VELOCITY IS SET: VXI,VYI,VZI
C
C  FIND TYPE OF COLLISION: IFLAG
C    IFLAG=-1       :  NOTHING TO BE DONE, POST-COLLISION VELOCITY IS ALREADY SAMPLED
C                      FROM MAXWELLIAN
C                      E.G. (relaxation)-BGK APPROXIMATION.
C    IFLAG=0        :  ISOTROPIC, IN CENTER OF MASS
C    IFLAG=IFTFLG>0 :  FIT FORM OF INTERACTION POTENTIAL (ELASTIC COLLISION IREL)
C
C
      IF (NPBGKP(IPLS,1).NE.0.OR.MODCOL(5,0,IREL).EQ.-1) THEN

C  ELASTIC TEST PARTICLE COLLISION IN BGK APPROXIMATION (E.G.: NEUTRAL-NEUTRAL)
C  ENFORCE: EXCHANGE OF IDENTITY COLLISION IN LAB FRAME
C  (AS IT IS ALSO DEFAULT FOR CX COLLISIONS).
C  FOR EQUAL MASSES OF BOTH COLLISION PARTNERS THIS CORRESPONDS
C  TO SCATTERING ANGLE = PI IN COM FRAME
        IFLAG=-1

      ELSEIF (MODCOL(5,0,IREL).EQ.0) THEN
C  NEITHER BGK-TYPE, NOR INTERACTION POTENTIAL PARAMETERS ARE GIVEN
C  USE ISOTROPIC SCATTERING IN CENTER OF MASS FRAME
        IFLAG= 0
C  SOME STUFF HERE FOR COM ISOTROPIC COLLISIONS.....

      ELSEIF (MODCOL(5,0,IREL).GT.0) THEN
C  INTERACTION POTENTIAL IS GIVEN, get fit coefficients of interaction potential
        IREAC=MODCOL(5,0,IREL)
        IFLAG=IFTFLG(IREAC,0)   !cdr  is iftflg correctly set for
                                !     repulsive potential?
cdr                                check slreac. There default is set to iftflg=2 (Morse)
CDR REACDAT should not be used during trajectory generation, see comments in xstel.f
c   to be done.
        P(1:9)=REACDAT(IREAC)%POT%POLY%DBLPOL(1:9,1)
      ENDIF
C

C
C  NEXT: STEP 2, FIND PRE-COLLISION DATA AND, IF NEEDED, ALSO THE TOTAL CROSS-SECTION
C                IN CASE IFLAG=-1 (RELAXATION COLLISION), NOTHING OF ALL THAT IS NEEDED.
C
      IF (IFLAG.NE.-1) THEN
        RMN=RMASS
        RMI=RMASSP(IPLS)
        RMSI=1./(RMN+RMI)
C  REDUCED MASS
        RLMS=RMN*RMI*RMSI
C  RELATIVE VELOCITY AND RELATED DATA
        VRELX=VX-VXI
        VRELY=VY-VYI
        VRELZ=VZ-VZI
        VRQYZ=         VRELY**2+VRELZ**2
        VRELQ=VRELX**2+VRQYZ
        ER=RLMS*VRELQ*CVELI2
c
        VREL=SQRT(VRELQ)
        VRYZ=SQRT(VRQYZ+EPS60)
C  IMPACT PARAMETER --> SCATTERING ANGLE --> NEW VELOCITY
C  CENTER OF MASS VELOCITY
        VSX=(RMI*VXI+RMN*VX)*RMSI
        VSY=(RMI*VYI+RMN*VY)*RMSI
        VSZ=(RMI*VZI+RMN*VZ)*RMSI
C
C  THAT IS ALL WE NEED FOR ISOTROPIC SCATTERING IN COM SYSTEM
        IF (IFLAG.GT.0) THEN
C
C  TOTAL CROSS-SECTION (ONLY IF NOT ALREADY COMPUTED EARLIER AT THIS CALL)
C  FOR RANDOM SAMPLING OF IMPACT PARAMETER
          IF (NFLAG.NE.3) THEN
C  IN PREVIOUS CALL TO CROSS WE HAD ALREADY SAME PARAMETER VRELQ=|VI-V0|**2
            ELAB=LOG(VRELQ)+DEFEL(IREL)
            IREAC=MODCOL(5,1,IREL)
            IF (IREAC.EQ.0) GOTO 995
            CEL=EIRENE_CROSS(ELAB,IREAC,IREL,FACREL(IREL,1),'VELOEL 3')
          ENDIF
C
C  STEP 2 FINISHED, CROSS-SECTION CEL IS FOUND

C  NEXT: STEP 3 ,  FIND IMPACT PARAMETER B (--> SCATTERING ANGLE --> NEW VELOCITY)

          BMAX=SQRT(CEL*PIAI)/0.52917E-8
          B= SQRT(RANF_EIRENE( ))*BMAX
        END IF

      ENDIF
C
C  STEP 3 FINISHED, COLLISION DATA IN CENTER OF MASS FRAME ARE SET
C         IN CASE IFLAG.GT.0: ALSO IMPACT PARAMETER IS FOUND , IN UNITS: BOHR RADIA
C  NEXT: STEP 4
C
      IF (IFLAG.EQ.-1) THEN
C
C  THIS PART: ONLY (BGK-TYPE) RELAXATION TO A MAXWELLIAN, I.E., POST-COLLISION
C             TEST PARTICLE IS SAMPLED FROM (WEIGHTED OR UNWEIGHTED)
C             BULK POPULATION (E.G.: EL. BGK COLLISION; EXCHANGE OF IDENTITY)
C
        VELQ=VXI*VXI+VYI*VYI+VZI*VZI
        VEL=SQRT(VELQ)
        VN=1./VEL
        VELX=VXI*VN
        VELY=VYI*VN
        VELZ=VZI*VN
C
        GOTO 1000

      ELSEIF (IFLAG.EQ.0) THEN

C   ISOTROPIC SCATTERING ANGLE PH IN CENTER OF MASS FRAME
        CPH=1.0-2.0*RANF_EIRENE( )  ! [-1 1]
        SPH=DSQRT(1.0-CPH*CPH)      ! [ 0 1]
C
      ELSEIF (IFLAG.GT.0) THEN
C
C  THIS PART: FIND DEFLECTION ANGLE
C             BINARY COLLISION KINETICS
C
C  COLLISION PARAMETERS IFLAG, ER AND B ARE DEFINED NOW.
C
C  FIND DISTANCE OF CLOSEST APPROACH: RSTERN
C
        RS=EIRENE_RSTERN(ER,B,IFLAG,P)
C
C  SCATTERING INTEGRAL TO FIND DEFLECTION ANGLE CHI
C
        CALL EIRENE_GAUMEH (RS,ER,B,IFLAG,P,10,1,RESULT)
        CHI=PIA-2.*B/RS*RESULT
        CCHI=COS(CHI)
C
C  CONVERT FROM DEFLECTION ANGLE CHI TO OBSERVABLE SCATTERING ANGLE PH, [0,...,PI]
        PH=ACOS(CCHI)
        CPH=COS(PH)
        SPH=SQRT(1.0-CPH*CPH)
      ENDIF

C  SCATTERING ANGLE PH, COS(PH)=CPH AND SIN(PH)=SPH ARE SET
C
C  POLAR ANGLE
      EPS=PI2A*RANF_EIRENE( )
C
C  CONVERT POST-COLLISON VELOCITY BACK INTO LAB FRAME
C
      CEPS=COS(EPS)
      SEPS=SIN(EPS)
      VRSX=VRELX*CPH+SPH*SEPS*VRYZ
      VRSY=VRELY*CPH+SPH*(VREL*VRELZ*CEPS-VRELX*VRELY*SEPS)/VRYZ
      VRSZ=VRELZ*CPH-SPH*(VREL*VRELY*CEPS+VRELX*VRELZ*SEPS)/VRYZ
C
      VELX=VSX+RLMS/RMN*VRSX
      VELY=VSY+RLMS/RMN*VRSY
      VELZ=VSZ+RLMS/RMN*VRSZ
      VELQ=VELX*VELX+VELY*VELY+VELZ*VELZ
      VEL=SQRT(VELQ)
      VELX=VELX/VEL
      VELY=VELY/VEL
      VELZ=VELZ/VEL

C
C  STEP 4 FINISHED, POST-COLLISION VELOCITY IS SET IN LAB FRAME
C  NEXT: RETURN
C
 1000 CONTINUE
C
      RETURN
C
  995 CONTINUE
      WRITE (iunout,*)
     . 'ERROR IN VELOEL, NO ELASTIC CROSS-SECTION DATA AVAILABLE'
      CALL EIRENE_MASJ5 ('ITYP,IATM,IMOL,IION,IPLS                ',
     .                    ITYP,IATM,IMOL,IION,IPLS)
      CALL EIRENE_MASJ4 ('NFLAG, IFLAG, IREL, IDREAC      ',
     .                    NFLAG, IFLAG, IREL, IDREAC)
      CALL EIRENE_EXIT_OWN(1)
  999 CONTINUE
      WRITE (iunout,*)
     .  'PARAMETER ERROR IN SUBR. VELOCX. EXIT CALLED'
      CALL EIRENE_EXIT_OWN(1)
      END SUBROUTINE EIRENE_VELOEL


C
C
      FUNCTION EIRENE_RSTERN (ER,B,IFLAGI,P)
C
C     IFLAG=1: H+ + H, PURELY REPULSIVE POTENTIAL
C     IFLAG=2: MORSE POTENTIAL, He+ + He, H+ + Noble Gases, H+ + H2
C
C               PARAMETERS: P(1):
C                           P(2):
C                           P(3):
C
C  RSTERN(ER,B) IS THE LARGEST ROOT OF THE EQUATION:
C
C     FI(R):=1.-V(R)/ER-(B/R)**2=0.
C
C   HERE: ER COLLISION ENERGY (EV)
C         B  IMPACT PARAMETER
C         V  INTERACTION POTENTIAL (EV)
C
C     DIMENSION RV0(3),RVM(3),RVW(3),VM(3),VW(3),VSW(3)

      USE EIRMOD_PRECISION
      USE EIRMOD_COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: ER, B, P(*)
      INTEGER, INTENT(IN) :: IFLAGI
      REAL(DP) :: EIRENE_RSTERN
      EXTERNAL :: EIRENE_EXIT_OWN

cym     REAL(DP) :: EIRENE_RTSAF,EIRENE_FI

cym ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cym local variables moved to declaration section of module
cym IFLAG -> IFLAGI to avoid confusion with IFLAG from veloel
cym ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     DATA IFIRST /0/
C     DATA RV0/0.,0.99699,2.18039/
C     DATA RVM/0.,1.4556 ,2.835539/
C     DATA RVW/0.,1.99515,3.490687/
C     DATA VM /0.,-2.     ,-2.7/
C     DATA VW /0.,-1.5    ,-2.025/
C     DATA VSW/0.,1.284688,1.4283/
C
C
      IF (IFLAGI.EQ.1) THEN
C  FIND UPPER AND LOWER BOUND FOR NUMERICAL ROOT FINDER RTSAF
C  FOR IFLAGI=1, RLW=B IS A LOWER BOUND AND RUP=INF IS AN UPPER BOUND.
C  FIRST TRY TO FIND BETTER BOUNDS:
        RUP=B
C
   10   RLW=RUP
        RUP=RLW*2.
        FITEST=EIRENE_FI(RUP,ER,B,IFLAGI,P,DFI)
        IF (FITEST.LT.0.D0) GOTO 10
C
        TOL=1.D-7*(RUP-RLW)
        EIRENE_RSTERN=EIRENE_RTSAF(RLW,RUP,TOL,ER,B,IFLAGI,P)
C
      ELSEIF (IFLAGI.EQ.2) THEN
C  GENERALISED MORSE POTENTIAL
C  FIND UPPER AND LOWER BOUND FOR NUMERICAL ROOT FINDER RTSAFE
C  FOR IFLAGI=2 AND IFLAGI=3, THE FOLLOWING INTERVAL CAN BE USED:
C  R0 IS THE ROOT OF THE INTERACTION POTENTIAL V(R): V(R0)=0.
C  RM IS THE RADIUS OF THE MINIMUM OF V(R)           V(RM)= VM = -EPS
C  RW IS THE RADIUS OF THE POINT OF INFLECTION       V(RM)= VW
C        RM=P(4)
         R0=P(5)
         RW=P(6)
         VW=P(8)
C
        BQ=B*B
        RV=R0
C
C  CASE 1: B < R0,
C
        IF (RV.GT.B) THEN
          RLW=B
          RUP=RV
        ENDIF
C
C  CASE 2: R0 < B,
C
        IF (RV.LT.B) THEN
          RLW=RV
          RUP=B
C  TRY TO IMPROVE LOWER BOUND RLW
          FIW=1.-VW/ER-BQ/(RW**2)
          IF (FIW.LT.0.D0) RLW=RW
        ENDIF
C
        TOL=1.D-7*(RUP-RLW)
        EIRENE_RSTERN=EIRENE_RTSAF(RLW,RUP,TOL,ER,B,IFLAGI,P)
      ELSE
        GOTO 990
      ENDIF
C
      RETURN
  990 CONTINUE
      WRITE (iunout,*) 'ERROR IN FUNCTION EIRENE_RSTERN '
      CALL EIRENE_EXIT_OWN(1)
      END FUNCTION EIRENE_RSTERN


C
C
      SUBROUTINE EIRENE_GAUMEH(RS,ER,B,IFLAGI,P,N,IGAUS,RESULT)
C
C  IGAUS=1
C  GAUSS MEHLER QUADRATURE, N=5,10,20, W(X)=1./SQRT(1+X)/SQRT(1-X)
C                           INTEGRATION FROM A=-1 TO B=1
C                           TRANSFORMED TO A=0,B=1, AND
C                           W(X)=1./SQRT(X)/SQRT(1-X)
C

      USE EIRMOD_PRECISION
      USE EIRMOD_COMPRT, ONLY: IUNOUT

      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: RS, ER, B, P(*)
      REAL(DP), INTENT(OUT) :: RESULT
      INTEGER, INTENT(IN) :: IFLAGI, N, IGAUS
      INTEGER :: I, IFI
      EXTERNAL :: EIRENE_EXIT_OWN

cym ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cym local variables moved to declaration section of module
cym IFLAG -> IFLAGI to avoid confusion with IFLAG from veloel
cym ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      REAL(DP) :: X5(5),X10(10),X20(20)
      REAL(DP) :: W5, W10, W20
      REAL(DP) :: XG10(10),WG10(10)

      DATA X5 /0.95105654E+00_DP, 0.58778542E+00_DP, 0.31391647E-06_DP,
     .       - 0.58778363E+00_DP,-0.95105600E+00_DP/,
     .     W5 /0.62831837E+00/
      DATA X10/9.876884E-01_DP, 8.910065E-01_DP, 7.071068E-01_DP,
     .         4.539905E-01_DP, 1.564344E-01_DP,-1.564345E-01_DP,
     .        -4.539905E-01_DP,-7.071068E-01_DP,-8.910065E-01_DP,
     .        -9.876884E-01_DP/,
     .     W10/3.141593E-01_DP/
      DATA X20/9.969173E-01_DP, 9.723699E-01_DP, 9.238795E-01_DP,
     .         8.526402E-01_DP, 7.604060E-01_DP, 6.494480E-01_DP,
     .         5.224985E-01_DP, 3.826834E-01_DP, 2.334453E-01_DP,
     .         7.845905E-02_DP,-7.845914E-02_DP,-2.334454E-01_DP,
     .        -3.826835E-01_DP,-5.224986E-01_DP,-6.494481E-01_DP,
     .        -7.604060E-01_DP,-8.526402E-01_DP,-9.238796E-01_DP,
     .        -9.723700E-01_DP,-9.969174E-01_DP/,
     .     W20/1.570796E-01_DP/
C
C  IGAUS=2
C  GAUSS MEHLER QUADRATURE, N=10, W(X)=1./SQRT(1-X)
C                           INTEGRATION FROM A=0 TO B=1
      DATA XG10/
     .    .013695585480651072296_DP,
     .    .070758123420104485088_DP,
     .    .16782834791297652828_DP,
     .    .29588270759990962383_DP,
     .    .44298868539955669668_DP,
     .    .59543571523425259033_DP,
     .    .73901490631777356542_DP,
     .    .86034375925702257262_DP,
     .    .94811360601967719799_DP,
     .    .99414369156320387287_DP/
      DATA WG10/
     .    .035228014278304236584_DP,
     .    .081202859600773882686_DP,
     .    .12534409666821812627_DP,
     .    .16655348315340949082_DP,
     .    .20386023963448093937_DP,
     .    .23638906392303662412_DP,
     .    .26337727689835336360_DP,
     .    .28419221863676482994_DP,
     .    .29834597294520589473_DP,
     .    .30550677426145261144_DP/

      IF (IFIRST1.EQ.0) THEN
C SET ROOTS AND WEIGHTS FOR GAUSS QUADRATURE RULES
        IFIRST1=1
        DO 5 I=1,5
          X=0.5*(1.+X5(I))
          X5A(I)=SQRT(X*(1.-X))
          X5B(I)=1./X
    5   CONTINUE
        DO 10 I=1,10
          X=0.5*(1.+X10(I))
          X10A(I)=SQRT(X*(1.-X))
          X10B(I)=1./X
   10   CONTINUE
        DO 11 I=1,10
          X=XG10(I)
          XG10A(I)=SQRT(1.-X)
          XG10B(I)=1./X
   11   CONTINUE
        DO 20 I=1,20
          X=0.5*(1.+X20(I))
          X20A(I)=SQRT(X*(1.-X))
          X20B(I)=1./X
   20   CONTINUE
      ENDIF
C
      NFI=N
      IF (N.EQ.5) THEN
        DO 50 IFI=1,5
          AR(IFI)=RS*X5B(IFI)
   50   CONTINUE
C
        CALL EIRENE_FIVEC(ER,B,IFLAGI,P)
        SUM=0.D0
        DO 51 IFI=1,5
          F=AFI(IFI)
          IF (F.GT.0.D0)
     .    SUM=SUM+X5A(IFI)/SQRT(F)
   51   CONTINUE
        RESULT=W5*SUM
C
      ELSEIF (N.EQ.10) THEN
        IF (IGAUS.EQ.1) THEN
          DO 100 IFI=1,10
            AR(IFI)=RS*X10B(IFI)
  100     CONTINUE
C
          CALL EIRENE_FIVEC(ER,B,IFLAGI,P)
          SUM=0.D0
          DO 101 IFI=1,10
            F=AFI(IFI)
            IF (F.GT.0.D0)
     .      SUM=SUM+X10A(IFI)/SQRT(F)
  101     CONTINUE
          RESULT=W10*SUM
        ELSEIF (IGAUS.EQ.2) THEN
          DO 110 IFI=1,10
            AR(IFI)=RS*XG10B(IFI)
  110     CONTINUE
C
          CALL EIRENE_FIVEC(ER,B,IFLAGI,P)
          SUM=0.D0
          DO 111 IFI=1,10
            F=AFI(IFI)
            IF (F.GT.0.D0)
     .      SUM=SUM+XG10A(IFI)/SQRT(F)*WG10(IFI)
  111     CONTINUE
          RESULT=SUM
        ENDIF
C
      ELSEIF (N.EQ.20) THEN
        DO 200 IFI=1,20
          AR(IFI)=RS*X20B(IFI)
  200   CONTINUE
C
        CALL EIRENE_FIVEC(ER,B,IFLAGI,P)
        SUM=0.D0
        DO 201 IFI=1,20
          F=AFI(IFI)
          IF (F.GT.0.D0)
     .    SUM=SUM+X20A(IFI)/SQRT(F)
  201   CONTINUE
        RESULT=W20*SUM
C
      ELSE
        WRITE (iunout,*) 'ERROR IN GAUMEH, WRONG PARAMETER N'
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
      RETURN
      END SUBROUTINE EIRENE_GAUMEH

C
C
      FUNCTION EIRENE_RTSAF(X1,X2,XACC,ER,B,IFLAGI,P)

C  TAKEN FROM: NUMERICAL RECIPES, W.H.PRESS ET AL.,
C              CAMBRIDGE UNIV. PRESS, 1989, P258
C  MODIFIED, TO FIND LARGEST ROOT, IN CASE MORE THAN ONE ROOTS
C  AND TO SPEED UP

      USE EIRMOD_PRECISION
      USE EIRMOD_COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP) :: EIRENE_RTSAF
      REAL(DP), INTENT(IN) :: X1, X2, XACC, ER, B, P(*)
      INTEGER, INTENT(IN) :: IFLAGI
      INTEGER, PARAMETER :: MAXIT=100
CYM/HJL Variables moved to module scope
      REAL(DP) :: DF, XDIST, F, TEMP, DX, DXOLD, XH,
     .            XL, POI
      INTEGER :: J

      XL=X1
      XH=X2
C
      F=EIRENE_FI(XH,ER,B,IFLAGI,P,DF)
      XDIST=(XH-XL)
      DX=MIN(XDIST,F/DF)
      DXOLD=DX
      EIRENE_RTSAF=XH-DX
      POI = EIRENE_RTSAF
C
      F=EIRENE_FI(POI,ER,B,IFLAGI,P,DF)

      IF(F.LT.0.D0) THEN
        XL=EIRENE_RTSAF
      ELSE
        XH=EIRENE_RTSAF
      ENDIF
C
      DO 11 J=1,MAXIT
        IF(((EIRENE_RTSAF-XH)*DF-F)*((EIRENE_RTSAF-XL)*DF-F).GT.0.
     *      .OR. ABS(F+F).GT.ABS(DXOLD*DF) ) THEN
          DXOLD=DX
          DX=0.5*(XH-XL)
          EIRENE_RTSAF=XL+DX
          POI=EIRENE_RTSAF
        ELSE
          DXOLD=DX
          DX=F/DF
          TEMP=EIRENE_RTSAF
          EIRENE_RTSAF=EIRENE_RTSAF-DX
          POI=EIRENE_RTSAF
        ENDIF
C
        IF(ABS(DX).LT.XACC) RETURN
C
        F=EIRENE_FI(POI,ER,B,IFLAGI,P,DF)
        IF(F.LT.0.D0) THEN
          XL=EIRENE_RTSAF
        ELSE
          XH=EIRENE_RTSAF
        ENDIF
   11 CONTINUE
      WRITE (iunout,*) 'RTSAF EXCEEDING MAXIMUM ITERATIONS'
      RETURN
      END FUNCTION EIRENE_RTSAF

C
C
      FUNCTION EIRENE_FI(R,ER,B,IFLAGI,P,DFI)

C  EVALUATE EFFECTIVE POTENTIAL FUNCTION FI AT R
C  EVALUATE DFI(R)/DR AT R
C  RETURN FI=FI(R), DFI=DFI(R)/DR
C     --------------
C  IFLAGI=1: H+ + H
C  IFLAGI=2: H+ + NOBLE GASES, H+ + H2, HE+ + HE
C

      USE EIRMOD_PRECISION
      USE EIRMOD_COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP) :: EIRENE_FI
      REAL(DP), INTENT(IN) :: R, ER, B, P(*)
      REAL(DP), INTENT(OUT) :: DFI
      INTEGER, INTENT(IN) :: IFLAGI

      REAL(DP) :: DSS, U, DU, SS, RI, RIQ, G1, REFF, RR,
     .            EMRR, G2, EX2, B2, V, DV, EX, RLOW, R2, R3,
     .            G, RMI, EPS
      EXTERNAL :: EIRENE_EXIT_OWN

      B2=B*B
C
      IF(IFLAGI.EQ.1) THEN
C  INTERACTION POTENTIAL V(R): H+ + H
C  R IN A0, V IN EV
        RLOW=1.D-7
        R2=R*R
        R3=R2*R
C  FIND V=V(R) and DV=DV(R)/DR
        IF (R.GT.160.D0) THEN
          V=0.D0
          DV=0.D0
        ELSEIF (R.LT.RLOW) THEN
          R2=RLOW*RLOW
          R3=R2*RLOW
          EX=EXP(-RLOW)
          EX2=EX*EX
          SS=(1.+RLOW+R2/3.)*EX-1.
          RI=1./RLOW
          V=27.211*((RI-(1.+RI)*EX2-(1.+RLOW)*EX)/SS+RI)
          DV=0.D0
        ELSE
          EX=EXP(-R)
          EX2=EX*EX
          SS=(1.+R+R2/3.)*EX-1.
          RI=1./R
          V=27.211*((RI-(1.+RI)*EX2-(1.+R)*EX)/SS+RI)
          RIQ=RI*RI
          DSS=-R/3.*(1.+R)*EX
          U=RI-(1.+RI)*EX2-(1.+R)*EX
          DU=-RIQ+(RIQ+2.*RI+2.)*EX2+R*EX
          DV=27.211*((DU*SS-U*DSS)/(SS*SS)-RIQ)
        ENDIF
C  FIND FI=FI(R) AND DFI=DFI(R)/DR
        EIRENE_FI=1.-V/ER-B2/R2
        DFI=-DV/ER+2.*B2/R3
C
      ELSEIF (IFLAGI.EQ.2) THEN
C  INTERACTION POTENTIAL V(R): H+ + NOBLE GASES (MORSE LIKE POTENTIAL)
C     R IN A0, V IN EV
        EPS=P(1)
        G1=P(2)
        G2=P(3)
        RMI=1./P(4)
        R2=R*R
        R3=R2*R
C
        RR=R*RMI
        EMRR=1.-RR
C
C       G2=1.00+(1.0-G2)*MAX(0.D0,-EMRR)/EMRR
C       REFF=-G1*G2*EMRR
        IF (RR.LT.1.0) THEN
          G=G1
        ELSE
          G=G1*G2
        ENDIF
        REFF=-G*EMRR
C
        IF (REFF.GT.160.D0) THEN
          V=0.D0
          DV=0.D0
        ELSE
          EX=EXP(-REFF)
          EX2=EX*EX
          V=EPS*(EX2-(EX+EX))
          DV=-2.*EPS*RMI*G*(EX2-EX)
        ENDIF
        EIRENE_FI=1.-V/ER-B2/R2
        DFI=-DV/ER+2.*B2/R3
C
      ELSE
        WRITE (iunout,*) 'ERROR IN FUNCTION EIRENE_FI. IFLAG INVALID.'
        WRITE (iunout,*) 'IFLAG = ',IFLAGI
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
C
      RETURN
      END FUNCTION EIRENE_FI

C
C
      SUBROUTINE EIRENE_FIVEC(ER,B,IFLAGI,P)

C  VECTORIZED VERSION OF FUNCTION FI FOR GAUSS MEHLER QUADRATURE
C  EVALUATE EFFECTIVE POTENTIAL FUNCTION AT AR(I),I=1,NFI
C  NOTE: NFI.LE.128 IS NOT CHECKED, BUT USED
C  RETURN FI(AR(I)) IN THE ARRAY AFI(I),I=1,NFI
C     --------------
C  IFLAGI=1: H+ + H
C  IFLAGI=2: H+ + NOBLE GASES, H+ + H2, HE+ + HE
C

      USE EIRMOD_PRECISION
      USE EIRMOD_COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: ER, B, P(*)
      INTEGER, INTENT(IN) :: IFLAGI
      REAL(DP) :: RI, G1, SS, EX, EX2, RR, REFF, EMRR, G2, V,
     .          R, R2, B2, RLOW, G, RMI, EPS
      INTEGER :: IFI
      EXTERNAL :: EIRENE_EXIT_OWN

cym cccccccccccc transformed into regular module variables cccccc
cym      REAL(DP) :: AR(128), AFI(128)
cym      INTEGER :: NFI
cym      COMMON /CFI/ AR,AFI,NFI
cym      SAVE /CFI/    !CDR DO WE NEED THIS ? cym good question ;-)
cym ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      B2=B*B
C
      IF(IFLAGI.EQ.1) THEN
C  INTERACTION POTENTIAL V(R): H+ + H
C  R IN A0, V IN EV
        RLOW=1.D-7
        DO 1 IFI=1,NFI
          R=AR(IFI)
          R2=R*R
C  FIND V=V(R)
          IF (R.GT.160.D0) THEN
            V=0.D0
          ELSEIF (R.LT.RLOW) THEN
            R2=RLOW*RLOW
            EX=EXP(-RLOW)
            EX2=EX*EX
            SS=(1.+RLOW+R2/3.)*EX-1.
            RI=1./RLOW
            V=27.211*((RI-(1.+RI)*EX2-(1.+RLOW)*EX)/SS+RI)
          ELSE
            EX=EXP(-R)
            EX2=EX*EX
            SS=(1.+R+R2/3.)*EX-1.
            RI=1./R
            V=27.211*((RI-(1.+RI)*EX2-(1.+R)*EX)/SS+RI)
          ENDIF
          AFI(IFI)=1.-V/ER-B2/R2
    1   CONTINUE
C
      ELSEIF(IFLAGI.EQ.2) THEN
C  INTERACTION POTENTIAL V(R): H+ + NOBLE GASES, (MORSE LIKE POTENTIAL)
C     R IN A0, V IN EV
        EPS=P(1)
        G1=P(2)
        G2=P(3)
        RMI=1./P(4)
        DO 2 IFI=1,NFI
          R=AR(IFI)
          R2=R*R
          RR=R*RMI
          EMRR=1.-RR
          IF (RR.LT.1) THEN
            G=G1
          ELSE
            G=G1*G2
          ENDIF
          REFF=-G*EMRR
C
          IF (REFF.GT.160.D0) THEN
            V=0.D0
          ELSE
            EX=EXP(-REFF)
            EX2=EX*EX
            V=EPS*(EX2-(EX+EX))
          ENDIF
          AFI(IFI)=1.-V/ER-B2/R2
    2   CONTINUE
C
      ELSE
        WRITE (iunout,*)
     .  'ERROR IN FUNCTION EIRENE_FIVEC. IFLAG INVALID.'
        WRITE (iunout,*) 'IFLAG = ',IFLAGI
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
C
      RETURN
      END SUBROUTINE EIRENE_FIVEC


C  the following SUBROUTINE is for reinitialization of EIRENE
      SUBROUTINE EIRENE_VELOEL_REINIT
      IMPLICIT NONE
      IFIRST = 0
      RETURN
      END SUBROUTINE EIRENE_VELOEL_REINIT

      END MODULE EIRMOD_VELOEL
