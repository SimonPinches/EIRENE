      MODULE EIRMOD_REFUSR
c
c--------------------------------------------------------------------------
c  REFUSR submodule for handling ammonia production
c
c         Author  : Nathan Bartlett
c         Contact : nbb2@illinois.edu
c--------------------------------------------------------------------------
c
c
c--------------------------------------------------------------------------
c  Import modules
c--------------------------------------------------------------------------
c
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT
      USE EIRMOD_CTRCEI
      USE EIRMOD_CCONA
      IMPLICIT NONE
      PRIVATE

      PUBLIC :: EIRENE_REFUSR, EIRENE_REFUSR_INIT,
     .          EIRENE_RF2USR,
     .          EIRENE_SPTUSR, EIRENE_SPTUSR_INIT,
     .          EIRENE_DEALLOC_REFUSR

      REAL(DP) :: AW         ! Initial distribution for N2 on W
      REAL(DP) :: BW         ! Mid-point temperature for N2 on W
      REAL(DP) :: M1W        ! Dependence on surface temp (%/eV) for W
      REAL(DP) :: ASS        ! Initial distribution for N2 on SS
      REAL(DP) :: BSS        ! Mid-point temperature for N2 on SS
      REAL(DP) :: M1SS       ! Dependence on surface temp (%/eV) for SS
      REAL(DP) :: TNORM      ! Normalization of temperature to 300K

      LOGICAL, ALLOCATABLE :: IS_ND(:), IS_ND2(:)
      LOGICAL, ALLOCATABLE :: IS_N(:)
      INTEGER :: N_ATOM, N2_MOL, ND3_MOL

      CONTAINS
c
c--------------------------------------------------------------------------
c  Subroutine REFUSR_INIT, set initial coefficients
c--------------------------------------------------------------------------
c    
      SUBROUTINE EIRENE_REFUSR_INIT
      IMPLICIT NONE
      INTEGER :: I          !Index for do loop
      CHARACTER*(9) :: HLP_FRM
c
c--------------------------------------------------------------------------
c  Set the starting distribution coefficients for N2 vs ND3
c       and the coefficients of dependence on temperature
c
c  We use two linear functions (one for W and one for Fe/SS) such that
c  at high temperature the reflected molecule is N2, and
c  at low temperature the probability of reflection as ND3 increases
c--------------------------------------------------------------------------
c
      IF (ALLOCATED(IS_N)) RETURN
      IF (TRCREF) WRITE (IUNOUT,*) '***REFUSR_INIT ENTRY TAG***'
c
      AW   = .5         !For W Surface
      BW   = 800.*EVKEL !mid-point temperature converted from K to eV
      ASS  = .5         !For SS Surface
      BSS  = 500.*EVKEL !mid-point temperature converted from K to eV
c
      M1W  = 11.6045 !Slope of .1 / 100K
      M1SS = 11.6045 !Slope of .1 / 100K
c
      TNORM = 300.*EVKEL !300K
c
c--------------------------------------------------------------------------
c  Determine the MOL species index in a general format
c--------------------------------------------------------------------------
c
      N2_MOL = 0
      ND3_MOL = 0
      ALLOCATE(IS_ND(NMOL))
      ALLOCATE(IS_ND2(NMOL))
      IS_ND = .FALSE.
      IS_ND2 = .FALSE.
c
      DO I = 1, NMOL
c  N2  ?
        IF (NCHARM(I).EQ.14.AND.NPRT(I+NATM).EQ.2.AND.
     &      N2_MOL.EQ.0) N2_MOL = I
c  ND  ?
        IS_ND(I)  = NCHARM(I).EQ.8.AND.NPRT(I+NATM).EQ.2
c  ND2 ?
        IS_ND2(I) = NCHARM(I).EQ.9.AND.NPRT(I+NATM).EQ.3
c  ND3 ?
        IF (NCHARM(I).EQ.10.AND.NPRT(I+NATM).EQ.4.AND.
     &      ND3_MOL.EQ.0) ND3_MOL = I
      ENDDO
c
c--------------------------------------------------------------------------
c  Determine the ATM species index in a general format
c--------------------------------------------------------------------------
c
      N_ATOM = 0
      ALLOCATE(IS_N(NATM))
      IS_N = .FALSE.
c
      DO I = 1, NATM
        IF (NCHARA(I).EQ.7.AND.NPRT(I).EQ.1) THEN
          IF (N_ATOM.EQ.0) N_ATOM = I
          IS_N(I) = .TRUE.
        ENDIF
      ENDDO
c
      IF (TRCREF) THEN
        WRITE(IUNOUT,'(A,I3)') 'N_ATOM  :: ', N_ATOM
        WRITE(IUNOUT,'(A,I3)') 'N2_MOL  :: ', N2_MOL
        WRITE(IUNOUT,'(A,I3)') 'ND3_MOL :: ', ND3_MOL
        WRITE(hlp_frm,'(A,I3,A)') '(A,',NATM,'L1)'
        WRITE(IUNOUT,hlp_frm)  'IS_N    :: ', IS_N
        WRITE(hlp_frm,'(A,I3,A)') '(A,',NMOL,'L1)'
        WRITE(IUNOUT,hlp_frm)  'IS_ND   :: ', IS_ND
        WRITE(IUNOUT,hlp_frm)  'IS_ND2  :: ', IS_ND2
        CALL EIRENE_LEER(1)
      END IF
c
      IF (N_ATOM.GT.0.AND.N2_MOL.GT.0.AND.ND3_MOL.GT.0) THEN
        WRITE(IUNOUT,*) 'REFUSR REFLECTION COEFFICIENTS FOR W:'
        CALL EIRENE_MASR3('A, B, M1                ', AW, BW, M1W)
        WRITE(IUNOUT,*) 'REFUSR REFLECTION COEFFICIENTS FOR SS:'
        CALL EIRENE_MASR3('A, B, M1                ', ASS, BSS, M1SS)
        CALL EIRENE_LEER(1)
      END IF
c
      RETURN
      END SUBROUTINE EIRENE_REFUSR_INIT

      SUBROUTINE EIRENE_SPTUSR_INIT
      IMPLICIT NONE
      RETURN
      END SUBROUTINE EIRENE_SPTUSR_INIT

      SUBROUTINE EIRENE_SPTUSR
      IMPLICIT NONE
      RETURN
      END SUBROUTINE EIRENE_SPTUSR
c
c--------------------------------------------------------------------------
c  Begin the REFUSR subroutine called from REFLEC.f
c--------------------------------------------------------------------------
c
      SUBROUTINE EIRENE_REFUSR (XMW,XCW,XMP,XCP,IGASF,IGAST,ZCOS,ZSIN,
     .                          EXPI,RPROB,E0TERM,ITYP,MSURF,ISPZO,IRET)
      USE EIRMOD_RANF, ONLY: RANF_EIRENE
      IMPLICIT NONE
      REAL(DP), INTENT(IN)   :: XMW        !Wall mass weight
      REAL(DP), INTENT(IN)   :: XCW        !Wall material nucl. charge number
      REAL(DP), INTENT(IN)   :: XMP        !projectile mass weight
      REAL(DP), INTENT(IN)   :: XCP        !projectile nucl. charge number
      REAL(DP), INTENT(IN)   :: ZCOS       !??
      REAL(DP), INTENT(IN)   :: ZSIN       !??
      REAL(DP), INTENT(IN)   :: EXPI       !??
      REAL(DP), INTENT(IN)   :: RPROB      !??
      REAL(DP), INTENT(IN)   :: E0TERM     !Thermal particle energy in eV
      INTEGER, INTENT(IN)    :: IGASF      !Species index for fast particle reflection
      INTEGER, INTENT(INOUT) :: IGAST      !Species index for thermal re-emission
      INTEGER, INTENT(IN)    :: ISPZO      !Species index
      INTEGER, INTENT(IN)    :: ITYP       !1 = Atom, 2 = Mol
      INTEGER, INTENT(IN)    :: MSURF      !Wall surface number
      INTEGER, INTENT(OUT) :: IRET
      REAL(DP) :: TW         !Wall temperature in eV
      REAL(DP) :: FW         !Calculated distribution for N2 on W
      REAL(DP) :: FSS        !Calculated distribution for N2 on SS
      REAL(DP) :: RF         !For random distribution

c     IF (TRCREF) WRITE (IUNOUT,*) '***REFUSR ENTRY TAG***'
      IRET = 0
      IF (E0TERM.LT.0._DP) THEN
        TW = -E0TERM
      ELSE IF (E0TERM.GT.0_DP) THEN
        TW = E0TERM/1.5_DP
      ELSE
        TW = TNORM
      ENDIF
c
c--------------------------------------------------------------------------
c  Begin new model for treating Nitrogen atoms and reflecting as either
c       N2 or ND3
c--------------------------------------------------------------------------
c
      IF (IS_N(ISPZO)) THEN   !Is it an N Atom ?

c       IF (TRCREF) WRITE (IUNOUT,*) 'N ATOM Hitting Surface'

        IF (ND3_MOL.EQ.0.AND.N2_MOL.EQ.0) THEN

c         IF (TRCREF) WRITE (IUNOUT,*) 'NO N CHEMISTRY: RETURN N ATOM' !No chemistry available
          IGAST = N_ATOM
          IRET = 5

        ELSEIF (ND3_MOL.EQ.0.AND.N2_MOL.NE.0) THEN

          IF (NINT(XCW).EQ.74) THEN     !Is it on a W surface

c           IF (TRCREF) WRITE (IUNOUT,*) 'NO AMMONIA MODEL: REFLECT N2'
            IGAST = -N2_MOL
            IRET = 2

          ELSEIF (NINT(XCW).EQ.26) THEN !Is it on Fe or SS surface

c           IF (TRCREF) WRITE (IUNOUT,*) 'NO AMMONIA MODEL: REFLECT N2'
            IGAST = -N2_MOL
            IRET = 2

          ELSE

c           IF (TRCREF) WRITE (IUNOUT,*) 'NOT ON CATALYST: RETURN N' !Not on Catalyst
            IGAST = N_ATOM
            IRET = 5

          ENDIF

        ELSE IF (ND3_MOL.NE.0.AND.N2_MOL.NE.0) THEN

          RF = RANF_EIRENE( )
c         IF (TRCREF) THEN
c           WRITE (IUNOUT,*) 'SURFACE ATOM NUMBER: ', XCW
c           WRITE (IUNOUT,*) 'SURFACE WEIGHT     : ', XMW
c           WRITE (IUNOUT,*) 'SURFACE TEMP       : ', TW
c           WRITE (IUNOUT,*) 'RANDOM NUMBER      : ', RF
c         ENDIF

          IF (NINT(XCW).EQ.74) THEN                !Is it on a W surface

c           IF (TRCREF) WRITE (IUNOUT,*) 'SURFACE: W'
            FW = AW + M1W*(BW-TW)                  !Calculate distribution
            FW = MIN(1._DP,MAX(0._DP,FW))
c           IF (TRCREF) WRITE (IUNOUT,*) 'CALC DISTRIBUTION: ', FW

            IF ((RF.LT.FW).AND.(FW.GT.0)) THEN     !Determine species
c             IF (TRCREF) WRITE (IUNOUT,*) 'REFLECT AS ND3'
              IGAST = -ND3_MOL
              IRET = 2
            ELSE
c             IF (TRCREF) WRITE (IUNOUT,*) 'REFLECT AS N2'
              IGAST = -N2_MOL
              IRET = 2
            ENDIF

          ELSEIF (NINT(XCW).EQ.26) THEN            !Is it on Fe or SS surface

c           IF (TRCREF) WRITE (IUNOUT,*) 'SURFACE: SS (or Fe)'
            FSS = ASS + M1SS*(BW-TW)               !Calculate distribution
            FSS = MIN(1._DP,MAX(0._DP,FSS))
c           IF (TRCREF) WRITE (IUNOUT,*) 'CALC DISTRIBUTION: ', FSS

            IF ((RF.LT.FSS).AND.(FSS.GT.0)) THEN   !Determine species
c             IF (TRCREF) WRITE (IUNOUT,*) 'REFLECT ND3'
              IGAST = -ND3_MOL
              IRET = 2
            ELSE
c             IF (TRCREF) WRITE (IUNOUT,*) 'REFLECT N2'
              IGAST = -N2_MOL
              IRET = 2
            ENDIF

          ELSE

c           IF (TRCREF) WRITE (IUNOUT,*) 'NOT ON CATALYST: RETURN N' !Not on Catalyst
            IGAST = N_ATOM
            IRET = 5

          ENDIF

        ENDIF

      ELSE

c       IF (TRCREF) WRITE (IUNOUT,*) 'NOT N ATOM Hitting Surface' !Not an N Atom
        IRET = 5

      ENDIF

      RETURN
      END SUBROUTINE EIRENE_REFUSR

      SUBROUTINE EIRENE_RF2USR (IMOL,MOL_DEFAULT)
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: IMOL        !Molecule index
      INTEGER, INTENT(IN)    :: MOL_DEFAULT !Default molecule reflection
c
c     IF (TRCREF) WRITE (IUNOUT,*) '***RF2USR ENTRY TAG***'
c
c--------------------------------------------------------------------------
c  If ND or ND2, then reflect back as ND3,
c  otherwise, will reflect as default species set in input
c--------------------------------------------------------------------------
c
c     IF (TRCREF) WRITE (IUNOUT,*) 'IMOL',IMOL

      IF ((IS_ND(IMOL)).OR.(IS_ND2(IMOL))) THEN
c       IF (TRCREF) WRITE (IUNOUT,*) 'REFLECT ND3 FROM SURFACE'
        IMOL = ND3_MOL
        RETURN
      ELSE
c       IF (TRCREF) WRITE (IUNOUT,*) 'REFLECT SOMETHING ELSE'
        IMOL = MOL_DEFAULT
        RETURN
      ENDIF
      END SUBROUTINE EIRENE_RF2USR

      SUBROUTINE EIRENE_DEALLOC_REFUSR
      IMPLICIT NONE

      IF (.NOT.ALLOCATED(IS_N)) RETURN

      DEALLOCATE(IS_N)
      DEALLOCATE(IS_ND)
      DEALLOCATE(IS_ND2)

      RETURN
      END SUBROUTINE EIRENE_DEALLOC_REFUSR

      END MODULE EIRMOD_REFUSR
