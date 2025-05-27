      module eirmod_profiles

      private

      public :: eirene_profe, eirene_profn, eirene_profr, eirene_profs

      interface eirene_profr
        module procedure eirene_profr_1d
        module procedure eirene_profr_2d
      end interface
      
      contains
CDR  june 17: COMMENTS, MINOR NOTATIONAL CHANGES
C
C
C
      SUBROUTINE EIRENE_PROFE (PRO,PRO0,RIN0,A1,E,SEP,PROVAC)
c
c  input profile parameters, see subr. input,  p0,...p5
c  p0: pro0
c  p1: rin0
c  p2: A1
c  p3: not used
c  p4: E
c  p5: SEP
C
C  EXPONENTIAL PROFILE, PLUS CONTINUATION EXPONENTIAL DECAY BEYOND "SEP"
C
C  PRO(X) =  PRO0   *EXP(-(X-RIN0)/A1)  FOR RINP0 <= X <= SEP
C  PRO(X) =  PROSEP *EXP(-(X-SEP)/A1)   FOR SEP   <= X <= RAA
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CCONA
      USE EIRMOD_CGRID
      USE EIRMOD_COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: PRO0, RIN0, A1, E, SEP, PROVAC
      REAL(DP), INTENT(OUT) :: PRO(:)
      REAL(DP) :: EP1R, ELLR, YR, EIRENE_AREAA, ARCA, RHO0, ZR,
     .            AR, DIST, PROSEP, RHOSEP
      INTEGER :: J, JM1, EIRENE_LEARCA, NLOCAL
      EXTERNAL :: EIRENE_LEARCA, EIRENE_AREAA

      RHOSEP=SEP
      RHO0=RIN0
      IF (LEVGEO.LE.2) THEN
        IF (RRA.GT.RAA) THEN
          NLOCAL=NR1STM
          PRO(NR1STM)=PROVAC
        ELSE
          NLOCAL=NR1ST
        ENDIF
C
C FIND modified RADIAL SURFACE LABELING COORDINATES RHO      AT "SEP"  :RHOSEP
C                                                        AND AT "RIN0" :RHO0
C
        IF (LEVGEO.EQ.2) THEN
          JM1=EIRENE_LEARCA(RHOSEP,RSURF,1,NLOCAL,1,'PROFE       ')
          AR=EIRENE_AREAA (RHOSEP,JM1,ARCA,YR,EP1R,ELLR)
          RHOSEP=SQRT(AR*PIAI)
          JM1=EIRENE_LEARCA(RHO0,RSURF,1,NLOCAL,1,'PROFE       ')
          AR=EIRENE_AREAA (RHO0,JM1,ARCA,YR,EP1R,ELLR)
          RHO0=SQRT(AR*PIAI)
        ENDIF

      ELSE
        WRITE (iunout,*)
     .  'WARNING: PROFE CALLED WITH LEVGEO.GT.2 '
        WRITE (iunout,*) 'NO PLASMA PARAMETERS RETURNED'
        RETURN
      ENDIF
C
C  VALUE OF PRO AT RHOSEP: PROSEP
      DIST=RHOSEP-RHO0
      PROSEP=PRO0*EXP(-DIST/(A1+EPS60))
C
      DO 20 J=1,NLOCAL-1
C  1ST EXPONENTIAL
        IF (RHOZNE(J).LE.RHOSEP) THEN
          ZR=RHOZNE(J)-RHO0
          PRO(J)=PRO0*EXP(-ZR/(A1+EPS60))
        ELSE
C  2ND EXPONENTIAL
          ZR=RHOZNE(J)-RHOSEP
          PRO(J)=PROSEP*EXP(-ZR/(E+EPS60))
        ENDIF
   20 CONTINUE
      PRO(NR1ST)=0.
      RETURN
      END SUBROUTINE EIRENE_PROFE
      
*****************************************************************
      
C sept.05:  activate this radial profile type also for levgeo=3
C           input radius SEP is relative to flux surface labelling grid RHOSRF
C           rra gt. raa indicates: one more vacuum cell in radial direction.
cdr jun 17: comments

C
      SUBROUTINE EIRENE_PROFN (PRO,PRO0,PROS,P,Q,E,SEP,PROVAC)
c  p1: pro0: value at x=rsurf(1)
c  p2: proS: value at x=sep
c  p3: profile parameter p, for rsurf(1) <= x <= sep
c  p4: profile parameter q, for rsurf(1) <= x <= sep
c  p5: E:    exponential decay length for  sep <= x
c  p6: sep:  x-position of sep:
C
C  PARABOLIC PROFILE, PLUS EXPONENTIAL DECAY BEYOND RHOSEP (SEP)
C
      USE EIRMOD_PRECISION, ONLY: DP
      USE EIRMOD_CCONA, ONLY: EPS60, PIAI
      USE EIRMOD_CGRID, ONLY: LEVGEO, NR1ST, NR1STM, RHOSRF, RHOZNE,
     >                        RAA, RRA, RSURF
      USE EIRMOD_COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: PRO0, PROS, P, Q, E, SEP, PROVAC
      REAL(DP), INTENT(OUT) :: PRO(:)
      REAL(DP) :: EP1R, ELLR, YR, EIRENE_AREAA, ARCA, RHOSEP, RORP,
     .            FACT, ZR, RDI, AR
      INTEGER  :: J, JM1, EIRENE_LEARCA, NLOCAL
      EXTERNAL :: EIRENE_AREAA, EIRENE_LEARCA

      IF (LEVGEO.LE.3) THEN
        IF (RRA.GT.RAA) THEN
          NLOCAL=NR1STM
          PRO(NR1STM)=PROVAC
        ELSE
          NLOCAL=NR1ST
        ENDIF
C
C FIND RADIAL SURFACE LABELING COORDINATE AT "SEP"
C
        select case (LEVGEO)
        case (1,3)
          RHOSEP=SEP
        case (2)
          JM1=EIRENE_LEARCA(SEP,RSURF,1,NLOCAL,1,'PROFN       ')
          AR=EIRENE_AREAA (SEP,JM1,ARCA,YR,EP1R,ELLR)
          RHOSEP=SQRT(AR*PIAI)
        end select
      ELSE
        WRITE (iunout,*)
     .  'WARNING: SUBR. PROFN CALLED WITH LEVGEO.GT.3 '
        WRITE (iunout,*) 'NO PLASMA PARAMETERS RETURNED'
        RETURN
      ENDIF
C
      RDI=1./(RHOSEP-RHOSRF(1))
      DO 20 J=1,NLOCAL-1
        IF (RHOZNE(J).LE.RHOSEP) THEN
          ZR=RHOZNE(J)-RHOSRF(1)
          RORP=ZR*RDI
          IF (RORP.LE.0.D0) THEN
            PRO(J)=PRO0
          ELSE
            FACT=(1.0-RORP**P)**Q
            PRO(J)=(PRO0-PROS)*FACT+PROS
          ENDIF
        ELSE
          PRO(J)=PROS*EXP((RHOSEP-RHOZNE(J))/(E+EPS60))
        ENDIF
   20 CONTINUE
      PRO(NR1ST)=0.
      RETURN
      END SUBROUTINE EIRENE_PROFN

********************************************************************
      
C     
C
      SUBROUTINE EIRENE_PROFR_1D (PRO,IINDEX,NSPZI,NSPZ1,NDAT)
C
C  READ ENTIRE PROFILE FROM TARGET DATA STRUCTURE PLASMA_BCKGRND  
C      (EIRMOD_CSPEI)
cdr  PRO is a 1D array, nspz1=1 necessarily (unused), 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CSPEI
      USE EIRMOD_COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IINDEX, NSPZI, NSPZ1, NDAT
      REAL(DP), INTENT(OUT) :: PRO(:)
      if (nspz1.ne.1) then
        write (iunout,*) 'error in PROFR'
        write (iunout,*) 'PRO: incorrect dimension in calling program'
        write (iunout,*) 'nspz1= ',NSPZ1
      endif

      PRO(1:NDAT) = PLASMA_BCKGRND(IINDEX+1,1:NDAT)
      RETURN
      END SUBROUTINE EIRENE_PROFR_1D

C
C
      SUBROUTINE EIRENE_PROFR_2D (PRO,IINDEX,NSPZI,NSPZ1,NDAT)
C
C  READ ENTIRE PROFILE FROM TARGET DATA STRUCTURE PLASMA_BCKGRND
C      (EIRMOD_CSPEI)
cdr  same as PROFR_1D, but PRO is a 2D array, 
C  NSPZ1: first dimension of PRO array as in calling program
C  NSPZI: fill the first NSPZI fields 1:NZPZI. NZPZI LE NSPZ1 necessarily.
C  NSPZ1,IINDEX:  = 1,     0        for TEIN
C                 = NPLSTI,1        for TIIN
C                 = NPLS,  NPLSTI   for DIIN
C  etc...
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CSPEI
      USE EIRMOD_COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IINDEX, NSPZI, NSPZ1, NDAT
!     REAL(DP), INTENT(OUT) :: PRO(NSPZ1,NDAT)
      REAL(DP), INTENT(OUT) :: PRO(:,:)
      if (nspzi.gt.nspz1 .or. nspz1.le.0) then
        write (iunout,*) 'error in PROFR'
        write (iunout,*) 'PRO: incorrect dimension in calling program'
        write (iunout,*) 'nspz1, nspzi= ',NSPZ1, NSPZI
      endif

      PRO(1:NSPZI,1:NDAT) = PLASMA_BCKGRND(IINDEX+1:IINDEX+NSPZI,1:NDAT)
      RETURN
      END SUBROUTINE EIRENE_PROFR_2D

**************************************************************************
      
C
      SUBROUTINE EIRENE_PROFS(PRO,PRO0,PRO1,SEP,PROVAC)
C
C  2-STEP PROFILE, VALUES PRO1, PRO2, 1:NR1ST, STEP AT POSITION SEP
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CCONA
      USE EIRMOD_CGRID
      USE EIRMOD_COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: PRO0, PRO1, SEP, PROVAC
      REAL(DP), INTENT(OUT) :: PRO(:)
      REAL(DP) :: EIRENE_AREAA, ARCA, AR, ELLR, YR, EP1R, RHOSEP
      INTEGER :: EIRENE_LEARCA, J, NLOCAL, JM1
      EXTERNAL :: EIRENE_AREAA, EIRENE_LEARCA

      RHOSEP=SEP
      IF (LEVGEO.LE.3) THEN
        IF (RRA.GT.RAA) THEN
          NLOCAL=NR1STM
          PRO(NR1STM)=PROVAC
        ELSE
          NLOCAL=NR1ST
        ENDIF
C
C FIND X SURFACE OR RADIAL SURFACE LABELING COORDINATE AT "SEP"
C
        select case (LEVGEO)
        case (1,3)
          RHOSEP=SEP
        case (2)
          JM1=EIRENE_LEARCA(SEP,RSURF,1,NLOCAL,1,'PROFS       ')
          AR=EIRENE_AREAA (SEP,JM1,ARCA,YR,EP1R,ELLR)
          RHOSEP=SQRT(AR*PIAI)
        end select
      ELSE
        WRITE (iunout,*)
     .  'WARNING: PROFS CALLED WITH LEVGEO.GT.3 '
        WRITE (iunout,*) 'CONSTANT PLASMA PARAMETERS RETURNED'
        NLOCAL=NR1ST
        DO 25 J=1,NLOCAL-1
          PRO(J)=PRO0
   25   CONTINUE
        PRO(NR1ST)=0.
        RETURN
      ENDIF
C
      DO 20 J=1,NLOCAL-1
        IF (RHOZNE(J).GT.RHOSEP) GOTO 15
        PRO(J)=PRO0
        GOTO 20
   15   PRO(J)=PRO1
   20 CONTINUE
      PRO(NR1ST)=0.
      RETURN
      END SUBROUTINE EIRENE_PROFS
      
      end module eirmod_profiles
