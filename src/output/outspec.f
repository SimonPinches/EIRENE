cdr  25.08.15:  formated spectrum printout improved
cdr  26.09.14:  commments, units added
cdr  oct.2014:  parameter istr (stratum number) in argument list


!pb  17.05.10:  write spectrum if the integral is nonzero
!               this change is necessary because spectra for bulk ions are sampled 
!               using negative weights
!pb  25.10.06:  format specifications corrected
 
      SUBROUTINE EIRENE_OUTSPEC(ISTR)
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_CCONA
      USE EIRMOD_CTRCEI
      USE EIRMOD_CTEXT
      USE EIRMOD_CSDVI
 
      IMPLICIT NONE
      INTEGER , INTENT(IN) :: ISTR
      INTEGER :: IADTYP(0:4)
      INTEGER :: IOUT, ISPC, I, IT, IE, IEND, IINI
      REAL(DP) :: EN,EN1,EN2
      CHARACTER(10) :: TEXTYP(0:4)
      CHARACTER(8) :: UNITINT(1:3),UNITOUT
 
C  SPECTRA
 
cdr   IOUT = 20+ifoff
cdr   OPEN (UNIT=IOUT,FILE='spectra.out')

      IOUT=IUNOUT
 
      TEXTYP(0) = 'PHOTONS   '
      TEXTYP(1) = 'ATOMS     '
      TEXTYP(2) = 'MOLECULES '
      TEXTYP(3) = 'TEST IONS '
      TEXTYP(4) = 'BULK IONS '

      UNITINT(1)= '(AMP)   '
      UNITINT(2)= '(WATT)  '
      UNITOUT   = '(?)     '

      IADTYP(0:4) = (/ 0, NSPH, NSPA, NSPAM, NSPAMI /)
 
      DO ISPC=1,NADSPC
        I = ESTIML(ISPC)%ISPCSRF
        IT = ESTIML(ISPC)%ISPCTYP
 
        WRITE (IOUT,*)
        IF (ISTR.GT.0) WRITE (IOUT,*) 'STRATUM NUMBER: ISTRA = ',istr
        IF (ISTR.EQ.0) WRITE (IOUT,*) 'SUM OVER STRATA'
        WRITE (IOUT,*)
 
        IF (ESTIML(ISPC)%ISRFCLL == 0)  THEN
c  surface-averaged spectra
          IF (I > NLIM) THEN
            WRITE (IOUT,'(A,A,I6)') ' SPECTRUM CALCULATED FOR',
     .                     ' NON-DEFAULT STANDARD SURFACE ',I-NLIM
          ELSE
            WRITE (IOUT,'(A,A,I6)') ' SPECTRUM CALCULATED FOR',
     .                     ' ADDITIONAL SURFACE ',I
          END IF
          IF (ESTIML(ISPC)%IDIREC > 0) THEN
            WRITE (iunout,'(A,3(ES12.4,A1))')
     .      ' IN DIRECTION (',ESTIML(ISPC)%SPCVX,',',
     .      ESTIML(ISPC)%SPCVY,',',ESTIML(ISPC)%SPCVZ,')'
          END IF
          IF (IT == 1) THEN
            WRITE (IOUT,'(A,A)') ' TYPE OF SPECTRUM : ',
     .                'INCIDENT PARTICLE FLUX IN AMP/BIN(EV)   '
            UNITOUT=UNITINT(1)
      
          ELSE IF (IT == 2) THEN
            WRITE (IOUT,'(A,A)') ' TYPE OF SPECTRUM : ',
     .                'INCIDENT ENERGY FLUX IN WATT/BIN(EV)    '
            UNITOUT=UNITINT(2)
          END IF
c  "cell based spectra"
        ELSE IF (ESTIML(ISPC)%ISRFCLL == 1)  THEN
          WRITE (IOUT,'(A,A,I6)') ' SPECTRUM CALCULATED FOR',
     .                   ' SCORING CELL ',I
          IF (ESTIML(ISPC)%IDIREC > 0) THEN
            WRITE (iunout,'(A,3(ES12.4,A1))')
     .      ' IN DIRECTION (',ESTIML(ISPC)%SPCVX,',',
     .      ESTIML(ISPC)%SPCVY,',',ESTIML(ISPC)%SPCVZ,')'
          END IF
          IF (IT == 1) THEN
            WRITE (iunout,'(A20,A)') ' TYPE OF SPECTRUM : ',
     .        'SPECTRAL PARTICLE DENSITY IN #/CM**3/BIN(EV)   '
          ELSEIF (IT == 2) THEN
            WRITE (iunout,'(A20,A)') ' TYPE OF SPECTRUM : ',
     .        'SPECTRAL ENERGY DENSITY IN EV/CM**3/BIN(EV)    '
          ELSEIF (IT == 3) THEN
            WRITE (iunout,'(A20,A)') ' TYPE OF SPECTRUM : ',
     .        'SPECTRAL MOMENTUM DENSITY IN (G*CM/S)/CM**3/BIN(EV)    '
          END IF
c  directional spectra in cell 
        ELSE IF (ESTIML(ISPC)%ISRFCLL == 2)  THEN
          WRITE (IOUT,'(A,A,I6)') ' SPECTRUM CALCULATED FOR',
     .                   ' GEOMETRICAL CELL ',I
          IF (ESTIML(ISPC)%IDIREC > 0) THEN
            WRITE (iunout,'(A,3(ES12.4,A1))')
     .      ' IN DIRECTION (',ESTIML(ISPC)%SPCVX,',',
     .      ESTIML(ISPC)%SPCVY,',',ESTIML(ISPC)%SPCVZ,')'
          END IF
          IF (IT == 1) THEN
            WRITE (iunout,'(A20,A)') ' TYPE OF SPECTRUM : ',
     .        'SPECTRAL PARTICLE DENSITY IN #/CM**3/BIN(EV)   '
          ELSEIF (IT == 2) THEN
            WRITE (iunout,'(A20,A)') ' TYPE OF SPECTRUM : ',
     .        'SPECTRAL ENERGY DENSITY IN EV/CM**3/BIN(EV)    '
          ELSEIF (IT == 3) THEN
            WRITE (iunout,'(A20,A)') ' TYPE OF SPECTRUM : ',
     .        'SPECTRAL MOMENTUM DENSITY IN (G*CM/S)/CM**3/BIN(EV)    '
          END IF
        END IF
 
        WRITE (IOUT,'(A20,A9)') ' TYPE OF PARTICLE : ',
     .         TEXTYP(ESTIML(ISPC)%IPRTYP)
        IF (ESTIML(ISPC)%IPRSP == 0) THEN
          WRITE (IOUT,'(A10,10X,A16)') ' SPECIES :',
     .                'SUM OVER SPECIES'
        ELSE
          WRITE (IOUT,'(A10,10X,A8)') ' SPECIES :',
     .          TEXTS(IADTYP(ESTIML(ISPC)%IPRTYP)+
     .          ESTIML(ISPC)%IPRSP)
        END IF

 
        IF (ESTIML(ISPC)%LOG) THEN
          WRITE (IOUT,'(A15,5X,ES12.4)') ' MINIMAL ENERGY ',
     .           10._DP**ESTIML(ISPC)%SPCMIN
          WRITE (IOUT,'(A15,5X,ES12.4)') ' MAXIMAL ENERGY ',
     .           10._DP**ESTIML(ISPC)%SPCMAX
          WRITE (IOUT,'(A)') ' LOGARITHMIC SPACING'
        ELSE
        WRITE (IOUT,'(A20,5X,ES12.4)') ' MINIMAL ENERGY (EV) ',
     .         ESTIML(ISPC)%SPCMIN
        WRITE (IOUT,'(A20,5X,ES12.4)') ' MAXIMAL ENERGY (EV) ',
     .         ESTIML(ISPC)%SPCMAX
          WRITE (IOUT,'(A)') ' LINEAR SPACING'
        END IF
        WRITE (IOUT,'(A20,4x,I6)') ' NUMBER OF BINS     ',
     .         ESTIML(ISPC)%NSPC
        CALL EIRENE_LEER(1)
C  HEADER DONE.

C  FORMATTED PRINTOUT OF SPECTRA STARTS HERE

        WRITE (IOUT,*)
        IF (ABS(ESTIML(ISPC)%SPCS) > EPS60) THEN
          IINI=0
          IEND=ESTIML(ISPC)%NSPC+1

          IF (NSIGI_SPC == 0) THEN
C  STANDARD DEVIATION IS NOT AVAILABLE
            WRITE (IOUT,'(A6,3A12)') '  BIN ',
     .                  '  B-LEFT    ','  B-RIGHT   ','  FLUX/BIN  '
            DO IE=IINI,IEND
C  DEAL WITH ENERGY BIN NO. IE
C  central energy bin value
              EN = ESTIML(ISPC)%SPCMIN +
     .            (IE-0.5_DP)*ESTIML(ISPC)%SPCDEL
c  LOWER energy bin value
              EN1= ESTIML(ISPC)%SPCMIN +
     .             (IE-1)*ESTIML(ISPC)%SPCDEL
c  UPPER energy bin value
              EN2= ESTIML(ISPC)%SPCMIN +
     .             (IE  )*ESTIML(ISPC)%SPCDEL
              IF (ESTIML(ISPC)%LOG) THEN
                EN = 10._DP**EN
                EN1= 10._DP**EN1
                EN2= 10._DP**EN2
              END IF

              IF (IE.EQ.IINI) EN1=0.0_DP ! even for log energy binning
              IF (IE.EQ.IEND) THEN
                WRITE (IOUT,'(I6,1ES12.4,A12,1ES12.4)') IE,EN1,
     .                                  ' INF       ',
     .                 ESTIML(ISPC)%SPC(IE)
              ELSE
                WRITE (IOUT,'(I6,3ES12.4)') IE,EN1,EN2,
     .                 ESTIML(ISPC)%SPC(IE)
              ENDIF
c  first and last bin: all the fluxes outside specified spectral range
              IF (IE.EQ.IINI.OR.IE.EQ.IEND-1)
     .          WRITE (IOUT,*) '.......................................'   
            END DO
          ELSE
c
C  STANDARD DEVIATION IS AVAILABLE
C
            WRITE (IOUT,'(A6,3A12,A14)') '  BIN ',
     .                  '  B-LEFT    ','  B-RIGHT   ','  FLUX/BIN  ',
     .                  '  STD.DEV. (%)'
            DO IE=IINI,IEND
C  DEAL WITH ENERGY BIN NO. IE
C  central energy bin value
              EN = ESTIML(ISPC)%SPCMIN +
     .            (IE-0.5_DP)*ESTIML(ISPC)%SPCDEL
c  LOWER energy bin value
              EN1= ESTIML(ISPC)%SPCMIN +
     .             (IE-1)*ESTIML(ISPC)%SPCDEL
c  UPPER energy bin value
              EN2= ESTIML(ISPC)%SPCMIN +
     .             (IE  )*ESTIML(ISPC)%SPCDEL
              IF (ESTIML(ISPC)%LOG) THEN
                EN =10._DP**EN
                EN1=10._DP**EN1
                EN2=10._DP**EN2
              END IF

              IF (IE.EQ.IINI) EN1=0.0_DP ! even for log energy binning
              IF (IE.EQ.IEND) THEN
                WRITE (IOUT,'(I6,1ES12.4,A12,2ES12.4)') IE,EN1,
     .                                  ' INF       ',
     .                 ESTIML(ISPC)%SPC(IE),
     .                 ESTIML(ISPC)%SGM(IE)
              ELSE
                WRITE (IOUT,'(I6,4ES12.4)') IE,EN1,EN2,
     .                 ESTIML(ISPC)%SPC(IE),
     .                 ESTIML(ISPC)%SGM(IE)
              ENDIF
c  first and last bin: all the fluxes outside specified spectral range
              IF (IE.EQ.IINI.OR.IE.EQ.IEND-1)
     .          WRITE (IOUT,*) '.......................................'   
            END DO 

          END IF
        ELSE
          WRITE (IOUT,'(A)') ' SPECTRUM IDENTICALLY 0 '
        END IF
C
C  PRINTOUT OF ENERGY INTEGRAL OVER SPECTRA
        WRITE (IOUT,*)
        WRITE (IOUT,'(A,A,ES12.4)') ' INTEGRAL OF SPECTRUM ',
     .             UNITOUT,ESTIML(ISPC)%SPCS 
        IF (NSIGI_SPC > 0)
     .    WRITE (IOUT,'(A,A,ES12.4)') ' STANDARD DEVIATION   ',
     .                  ' %      ',ESTIML(ISPC)%SGMS 
      END DO
 
      RETURN
      END SUBROUTINE EIRENE_OUTSPEC
