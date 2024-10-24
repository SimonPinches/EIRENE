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
      USE EIRMOD_WNEUTRALS
      IMPLICIT NONE
      PRIVATE

      PUBLIC :: EIRENE_REFUSR, EIRENE_REFUSR_INIT,
     .          EIRENE_RF2USR,
     .          EIRENE_SPTUSR_INIT,
     .          EIRENE_DEALLOC_REFUSR, EIRENE_WRITE_CONBE

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

      REAL(DP), ALLOCATABLE, SAVE, PUBLIC :: CONBE(:)
      REAL(DP), ALLOCATABLE, SAVE, PUBLIC :: CONBEPRIVIOS(:)
      REAL(DP), ALLOCATABLE, SAVE, PUBLIC :: SPUMPADD(:)
      REAL(DP), ALLOCATABLE, SAVE, PUBLIC :: SPUMPSAVE(:)
C     REAL(DP), ALLOCATABLE, SAVE, PUBLIC :: SPUMPPRIVIOS(:)
      REAL(DP), ALLOCATABLE, SAVE, PUBLIC :: SPTTOTADD(:)
      REAL(DP), ALLOCATABLE, SAVE, PUBLIC :: SPTTOTSAVE(:)
      REAL(DP), ALLOCATABLE, SAVE, PUBLIC :: DENSBE(:)
      REAL(DP), ALLOCATABLE, SAVE, PUBLIC :: TSTEP_ARR(:)

      INTEGER, SAVE :: IB2STEP, BEGINSURF, ENDSURF,
     .                 NOACUM_FLAG, CONBEOUTDELAY
      REAL(DP), SAVE :: SPUMPBE, SPTTOTBE, AVCOF, SBE, DENSBELOC
      CHARACTER*12, SAVE :: HLP_FRM

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
      EXTERNAL :: EIRENE_LEER, EIRENE_MASR3
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
c
c--------------------------------------------------------------------------
c  Begin the REFUSR subroutine called from REFLEC.f
c--------------------------------------------------------------------------
c
      SUBROUTINE EIRENE_REFUSR (XMW,XCW,XMP,XCP,IGASF,IGAST,ZCOS,ZSIN,
     .                          EXPI,RPROB,E0TERM,ITYP,MSURF,ISPZO,IRET)
      USE EIRMOD_RANF, ONLY: RANF_EIRENE
      IMPLICIT NONE
      REAL(DP), INTENT(IN)   :: XMW   !Wall mass weight
      REAL(DP), INTENT(IN)   :: XCW   !Wall material nucl. charge number
      REAL(DP), INTENT(IN)   :: XMP   !projectile mass weight
      REAL(DP), INTENT(IN)   :: XCP   !projectile nucl. charge number
      REAL(DP), INTENT(IN)   :: ZCOS  !??
      REAL(DP), INTENT(IN)   :: ZSIN  !??
      REAL(DP), INTENT(IN)   :: EXPI  !??
      REAL(DP), INTENT(IN)   :: RPROB !??
      REAL(DP), INTENT(IN)   :: E0TERM  !Thermal particle energy in eV
      INTEGER, INTENT(IN)    :: IGASF !Species index for
                                      !fast particle reflection
      INTEGER, INTENT(INOUT) :: IGAST !Species index for
                                      !thermal re-emission
      INTEGER, INTENT(IN)    :: ISPZO !Species index
      INTEGER, INTENT(IN)    :: ITYP  !1 = Atom, 2 = Mol
      INTEGER, INTENT(IN)    :: MSURF !Wall surface number
      INTEGER, INTENT(OUT)   :: IRET
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

          IF (NINT(XCW).EQ.74) THEN              !Is it on a W surface

c           IF (TRCREF) WRITE (IUNOUT,*) 'SURFACE: W'
            FW = AW + M1W*(BW-TW)                !Calculate distribution
            FW = MIN(1._DP,MAX(0._DP,FW))
c           IF (TRCREF) WRITE (IUNOUT,*) 'CALC DISTRIBUTION: ', FW

            IF ((RF.LT.FW).AND.(FW.GT.0)) THEN   !Determine species
c             IF (TRCREF) WRITE (IUNOUT,*) 'REFLECT AS ND3'
              IGAST = -ND3_MOL
              IRET = 2
            ELSE
c             IF (TRCREF) WRITE (IUNOUT,*) 'REFLECT AS N2'
              IGAST = -N2_MOL
              IRET = 2
            ENDIF

          ELSEIF (NINT(XCW).EQ.26) THEN   !Is it on Fe or SS surface

c           IF (TRCREF) WRITE (IUNOUT,*) 'SURFACE: SS (or Fe)'
            FSS = ASS + M1SS*(BW-TW)      !Calculate distribution
            FSS = MIN(1._DP,MAX(0._DP,FSS))
c           IF (TRCREF) WRITE (IUNOUT,*) 'CALC DISTRIBUTION: ', FSS

            IF ((RF.LT.FSS).AND.(FSS.GT.0)) THEN !Determine species
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

c--------------------------------------------------------------------------
c  SPTUSR submodule for dynamic follow-up of surface composition
c
c         Author  : Sergey Makarov
c
c--------------------------------------------------------------------------

      SUBROUTINE EIRENE_SPTUSR_INIT
      USE EIRMOD_PRECISION, ONLY: DP
      USE EIRMOD_PARMMOD, ONLY: NLIMPS
      USE EIRMOD_CADGEO, ONLY: NLIMI
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CPES, ONLY: MY_PE
      USE EIRMOD_MPI
      IMPLICIT NONE
      REAL(DP) :: TSTEP
      INTEGER :: TSTEP_FLAG, IER
      LOGICAL :: file_exists
      EXTERNAL :: EIRENE_CHECK_EXIT, EIRENE_EXIT_OWN
#if ( defined(USE_MPI) && !defined(GFORTRAN) && !defined(CRAY) )
      EXTERNAL :: MPI_BCAST
#endif

      IF (.NOT.ALLOCATED(DENSBE)) THEN
        ALLOCATE(DENSBE(0:NLIMPS))
        ALLOCATE(TSTEP_ARR(0:NLIMPS))
        DENSBE=0._DP
        TSTEP_ARR=0._DP
        IF (MY_PE == 0) THEN
          OPEN (UNIT=36,FILE="../input_for_refusr.dat",ERR=936)
          READ(36,*)
          READ(36,*) BEGINSURF, ENDSURF, CONBEOUTDELAY,
     .               SBE, TSTEP, AVCOF, TSTEP_FLAG, NOACUM_FLAG
          IF (BEGINSURF.LT.0) BEGINSURF=NLIMI+ABS(BEGINSURF)
          IF (ENDSURF.LT.0) ENDSURF=NLIMI+ABS(ENDSURF)
          WRITE(HLP_FRM,'(A,I4,A)') '(',ENDSURF-BEGINSURF+1,'ES14.7)'
          IF (TSTEP_FLAG == 1) THEN
            READ(36,*)
            READ(36,HLP_FRM) TSTEP_ARR(BEGINSURF:ENDSURF)
          ELSE
            TSTEP_ARR(BEGINSURF:ENDSURF)=TSTEP
          ENDIF
          CLOSE (UNIT=36)
          WRITE (IUNOUT,*) 'BEGINSURF', BEGINSURF
          WRITE (IUNOUT,*) 'ENDSURF', ENDSURF
          WRITE (IUNOUT,*) 'CONBEOUTDELAY', CONBEOUTDELAY
          WRITE (IUNOUT,*) 'SBE', SBE
          WRITE (IUNOUT,*) 'TSTEP', TSTEP
          WRITE (IUNOUT,*) 'AVCOF', AVCOF
          WRITE (IUNOUT,*) 'TSTEP_FLAG', TSTEP_FLAG
          WRITE (IUNOUT,*) 'NOACUM_FLAG', NOACUM_FLAG
          WRITE (IUNOUT,*)  TSTEP_ARR(BEGINSURF:ENDSURF)
          IB2STEP=0
          IF (BEGINSURF.LT.1 .OR. BEGINSURF.GT.NLIMPS .OR.
     .        ENDSURF.LT.BEGINSURF .OR. ENDSURF.GT.NLIMPS .OR.
     .        CONBEOUTDELAY.LT.1 .OR. SBE.LE.0._DP .OR.
     .        TSTEP.LE.0._DP .OR.
     .        MINVAL(TSTEP_ARR(BEGINSURF:ENDSURF)).LE.0._DP .OR.
     .        AVCOF.LT.0._DP .OR. AVCOF.GT.1._DP) THEN
            WRITE(IUNOUT,*) 'FAULTY INPUT in input_for_refusr.dat !'
            CALL EIRENE_EXIT_OWN(1)
          ENDIF
        ENDIF
        WRITE (IUNOUT,*) 'START DENSBE'
      ENDIF
      IF (.NOT.ALLOCATED(CONBE)) THEN
        ALLOCATE(CONBE(0:NLIMPS))
        CONBE(0:NLIMI)=0._DP
        CONBE(NLIMI+1:NLIMPS)=0.01_DP
        WRITE (IUNOUT,*) 'START CONBE'
      ENDIF
      IF (.NOT.ALLOCATED(SPUMPADD)) THEN
        ALLOCATE(SPUMPADD(0:NLIMPS))
        SPUMPADD=0._DP
        WRITE (IUNOUT,*) 'START SPUMPADD'
      ENDIF
      IF (.NOT.ALLOCATED(SPUMPSAVE)) THEN
        ALLOCATE(SPUMPSAVE(0:NLIMPS))
        SPUMPSAVE=0._DP
        WRITE (IUNOUT,*) 'START SPUMPSAVE'
       ENDIF
       IF (.NOT.ALLOCATED(SPTTOTADD)) THEN
        ALLOCATE(SPTTOTADD(0:NLIMPS))
        SPTTOTADD=0._DP
        WRITE (IUNOUT,*) 'START SPTTOTADD'
      ENDIF
      IF (.NOT.ALLOCATED(SPTTOTSAVE)) THEN
        ALLOCATE(SPTTOTSAVE(0:NLIMPS))
        SPTTOTSAVE=0._DP
        WRITE (IUNOUT,*) 'START SPTTOTSAVE'
        IF (MY_PE == 0) THEN
          INQUIRE(FILE="../conbe.dat", EXIST=file_exists)
          IF (file_exists) THEN
            OPEN (UNIT=36,FILE="../conbe.dat",position="append")
            BACKSPACE 36
            BACKSPACE 36
            BACKSPACE 36
            BACKSPACE 36
            READ(36,HLP_FRM) SPUMPSAVE(BEGINSURF:ENDSURF)
            READ(36,HLP_FRM) SPTTOTSAVE(BEGINSURF:ENDSURF)
            READ(36,*)
            READ(36,HLP_FRM) DENSBE(BEGINSURF:ENDSURF)
            CLOSE (UNIT=36)
          ENDIF
        ENDIF
      ENDIF

      CALL EIRENE_WRITE_CONBE('    ')
      CALL EIRENE_CHECK_EXIT
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
      CALL MPI_BCAST(CONBE,NLIMPS+1,MPI_REAL8,0,MPI_COMM_WORLD,ier)

      RETURN

  936 WRITE(IUNOUT,*) 'ERROR OPENING input_for_refusr.dat FILE !'
      CALL EIRENE_EXIT_OWN(1)
      END SUBROUTINE EIRENE_SPTUSR_INIT

      SUBROUTINE EIRENE_WRITE_CONBE(EDITION)
      USE EIRMOD_CPES, ONLY: MY_PE
      USE EIRMOD_CCOUPL, ONLY: NFLA
      IMPLICIT NONE
      CHARACTER*4, INTENT(IN) :: EDITION
      INTEGER :: K, ISPZ, ILIMPSUSER
      CHARACTER*14 :: FILENAME

      IF (.NOT.ALLOCATED(CONBE)) RETURN
      FILENAME = 'conbe.dat'
      if (edition.ne.'    ') filename = trim(filename)//'.'//edition

      IF (MY_PE == 0) THEN
        IF (edition.eq.'    ') IB2STEP=IB2STEP+1
C       WRITE (IUNOUT,*) 'IB2STEP', IB2STEP
        DO ILIMPSUSER=BEGINSURF, ENDSURF
        
C Instantaneous Be deposition flux calculation
          SPUMPADD(ILIMPSUSER)=0._DP
          DO K = 1, NATMI
            IF (NMASSA(K).EQ.9 .AND. NCHARA(K).EQ.4) THEN
              SPUMPADD(ILIMPSUSER)=SPUMPADD(ILIMPSUSER)+
     .         wlpump(K,ILIMPSUSER)
            END IF
          END DO
          DO K = 1, NFLA
            ISPZ = NATMI+NMOLI+NIONI+K
            IF (NMASSP(K).EQ.9 .AND. NCHARP(K).EQ.4) THEN
              SPUMPADD(ILIMPSUSER)=SPUMPADD(ILIMPSUSER)+
     .         wlpump(ISPZ,ILIMPSUSER)
            END IF
          END DO
          
C Time-averaged Be deposition flux calculation
          SPUMPADD(ILIMPSUSER)=SPUMPADD(ILIMPSUSER)/wlarea(ILIMPSUSER)
          SPUMPSAVE(ILIMPSUSER)=SPUMPSAVE(ILIMPSUSER)*
     .     AVCOF+SPUMPADD(ILIMPSUSER)*(1._DP-AVCOF)
          SPUMPBE=SPUMPSAVE(ILIMPSUSER)

C Instantaneous Be sputtering flux calculation
          SPTTOTADD(ILIMPSUSER)=wldspt(ILIMPSUSER,0)/wlarea(ILIMPSUSER)
          
C Time-averaged Be sputtering flux calculation
          SPTTOTSAVE(ILIMPSUSER)=SPTTOTSAVE(ILIMPSUSER)*
     .     AVCOF+SPTTOTADD(ILIMPSUSER)*(1._DP-AVCOF)
          SPTTOTBE=SPTTOTSAVE(ILIMPSUSER)

          DENSBELOC=DENSBE(ILIMPSUSER)
          
C Time-stepping
          DENSBE(ILIMPSUSER)=(SPUMPBE-SPTTOTBE)*
     .     TSTEP_ARR(ILIMPSUSER)+DENSBELOC
          DENSBE(ILIMPSUSER)=MAX(0._DP,DENSBE(ILIMPSUSER))

C Be concentration calculation
          IF (DENSBE(ILIMPSUSER)*SBE.LT.1._DP) THEN
            CONBE(ILIMPSUSER)=DENSBE(ILIMPSUSER)*SBE
          ELSE
            CONBE(ILIMPSUSER)=1._DP
            IF (NOACUM_FLAG == 1)  THEN
              DENSBE(ILIMPSUSER)=1._DP/SBE
            ENDIF
          ENDIF
        END DO
        
C Output
        IF (MOD(IB2STEP+CONBEOUTDELAY-1, CONBEOUTDELAY).EQ.0 .or.
     .      edition.ne.'    ') THEN
          write (iunout,*) 'Writing ',trim(filename)
          OPEN (UNIT=36,FILE=trim(FILENAME),position="append")
          WRITE (36,HLP_FRM) SPUMPADD(BEGINSURF:ENDSURF)
          WRITE (36,HLP_FRM) SPTTOTADD(BEGINSURF:ENDSURF)
          WRITE (36,HLP_FRM) SPUMPSAVE(BEGINSURF:ENDSURF)
          WRITE (36,HLP_FRM) SPTTOTSAVE(BEGINSURF:ENDSURF)
          WRITE (36,HLP_FRM) CONBE(BEGINSURF:ENDSURF)
          WRITE (36,HLP_FRM) DENSBE(BEGINSURF:ENDSURF)
          CLOSE (UNIT=36)
        ENDIF
      ENDIF

      RETURN
      END SUBROUTINE EIRENE_WRITE_CONBE

      END MODULE EIRMOD_REFUSR
