cdr  26.09.14:  commments, units added
cdr  oct.2014:  parameter istr (stratum number) in argument list

!pb  25.10.06:  format specifications corrected
!pb  17.05.10:  write spectrum if the integral is nonzero
!               this change is necessary because spectra for bulk ions are sampled 
!               using negative weights
 
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
      REAL(DP) :: EN
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
      UNITOUT   = '(???)   '

      IADTYP(0:4) = (/ 0, NSPH, NSPA, NSPAM, NSPAMI /)
 
      DO ISPC=1,NADSPC
        I = ESTIML(ISPC)%PSPC%ISPCSRF
        IT = ESTIML(ISPC)%PSPC%ISPCTYP
 
        WRITE (IOUT,*)
        IF (ISTR.GT.0) WRITE (IOUT,*) 'STRATUM NUMBER: ISTRA = ',istr
        IF (ISTR.EQ.0) WRITE (IOUT,*) 'SUM OVER STRATA'
        WRITE (IOUT,*)
 
        IF (ESTIML(ISPC)%PSPC%ISRFCLL == 0)  THEN
c  surface averaged spectra
          IF (I > NLIM) THEN
            WRITE (IOUT,'(A,A,I6)') ' SPECTRUM CALCULATED FOR',
     .                     ' NONDEFAULT STANDARD SURFACE ',I-NLIM
          ELSE
            WRITE (IOUT,'(A,A,I6)') ' SPECTRUM CALCULATED FOR',
     .                     ' ADDITIONAL SURFACE ',I
          END IF
          IF (ESTIML(ISPC)%PSPC%IDIREC > 0) THEN
            WRITE (iunout,'(A,3(ES12.4,A1))')
     .      ' IN DIRECTION (',ESTIML(ISPC)%PSPC%SPCVX,',',
     .      ESTIML(ISPC)%PSPC%SPCVY,',',ESTIML(ISPC)%PSPC%SPCVZ,')'
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
        ELSE IF (ESTIML(ISPC)%PSPC%ISRFCLL == 1)  THEN
          WRITE (IOUT,'(A,A,I6)') ' SPECTRUM CALCULATED FOR',
     .                   ' SCORING CELL ',I
          IF (ESTIML(ISPC)%PSPC%IDIREC > 0) THEN
            WRITE (iunout,'(A,3(ES12.4,A1))')
     .      ' IN DIRECTION (',ESTIML(ISPC)%PSPC%SPCVX,',',
     .      ESTIML(ISPC)%PSPC%SPCVY,',',ESTIML(ISPC)%PSPC%SPCVZ,')'
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
        ELSE IF (ESTIML(ISPC)%PSPC%ISRFCLL == 2)  THEN
          WRITE (IOUT,'(A,A)') ' SPECTRUM CALCULATED FOR',
     .                   ' GEOMETRICAL CELL ',I
          IF (ESTIML(ISPC)%PSPC%IDIREC > 0) THEN
            WRITE (iunout,'(A,3(ES12.4,A1))')
     .      ' IN DIRECTION (',ESTIML(ISPC)%PSPC%SPCVX,',',
     .      ESTIML(ISPC)%PSPC%SPCVY,',',ESTIML(ISPC)%PSPC%SPCVZ,')'
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
     .         TEXTYP(ESTIML(ISPC)%PSPC%IPRTYP)
        IF (ESTIML(ISPC)%PSPC%IPRSP == 0) THEN
          WRITE (IOUT,'(A10,10X,A16)') ' SPECIES :',
     .                'SUM OVER SPECIES'
        ELSE
          WRITE (IOUT,'(A10,10X,A8)') ' SPECIES :',
     .          TEXTS(IADTYP(ESTIML(ISPC)%PSPC%IPRTYP)+
     .          ESTIML(ISPC)%PSPC%IPRSP)
        END IF

 
        WRITE (IOUT,'(A19,5X,ES12.4)') ' MINIMAL ENERGY (EV) ',
     .         ESTIML(ISPC)%PSPC%SPCMIN
        WRITE (IOUT,'(A19,5X,ES12.4)') ' MAXIMAL ENERGY (EV) ',
     .         ESTIML(ISPC)%PSPC%SPCMAX
        WRITE (IOUT,'(A20,4x,I6)') ' NUMBER OF BINS     ',
     .         ESTIML(ISPC)%PSPC%NSPC
C  HEADER DONE.

C  FORMATTED PRINTOUT OF SPECTRA STARTS HERE

        WRITE (IOUT,*)
        IF (ABS(ESTIML(ISPC)%PSPC%SPCS) > EPS60) THEN
          IINI=0
          IEND=ESTIML(ISPC)%PSPC%NSPC+1
          IF (NSIGI_SPC == 0) THEN
C  STANDARD DEVIATION IS NOT AVAILABLE
            DO IE=IINI,IEND
              EN = ESTIML(ISPC)%PSPC%SPCMIN +
     .             (IE-0.5)*ESTIML(ISPC)%PSPC%SPCDEL
              WRITE (IOUT,'(I6,2ES12.4)') IE,EN,
     .               ESTIML(ISPC)%PSPC%SPC(IE)
c  first and last bin: all the fluxes outside specified spectral range
              IF (IE.EQ.IINI.OR.IE.EQ.IEND-1)
     .          WRITE (IOUT,*) '.......................................'   
            END DO
          ELSE
C  STANDARD DEVIATION IS AVAILABLE
C     
c  first bin: all the fluxes below specified spectral range
            EN = ESTIML(ISPC)%PSPC%SPCMIN  
            WRITE (IOUT,'(I6,A4,3ES12.4)') IINI,' <= ',EN,
     .             ESTIML(ISPC)%PSPC%SPC(IINI),
     .             ESTIML(ISPC)%PSPC%SGM(IINI)
            WRITE (IOUT,*) '.......................................'          
            DO IE=IINI+1,IEND-1
              EN = ESTIML(ISPC)%PSPC%SPCMIN +
     .             (IE-0.5)*ESTIML(ISPC)%PSPC%SPCDEL
              WRITE (IOUT,'(I6,A4,3ES12.4)') IE,'    ',EN,
     .               ESTIML(ISPC)%PSPC%SPC(IE),
     .               ESTIML(ISPC)%PSPC%SGM(IE)
            END DO
c  last bin: all the fluxes above specified spectral range
            WRITE (IOUT,*) '.......................................' 
            EN = ESTIML(ISPC)%PSPC%SPCMIN +
     .             (IEND-1)*ESTIML(ISPC)%PSPC%SPCDEL
            WRITE (IOUT,'(I6,A4,3ES12.4)') IEND,' >= ',EN,
     .             ESTIML(ISPC)%PSPC%SPC(IEND),
     .             ESTIML(ISPC)%PSPC%SGM(IEND)  

          END IF
        ELSE
          WRITE (IOUT,'(A)') ' SPECTRUM IDENTICAL 0 '
        END IF
C
C  PRINTOUT OF ENERGY INTEGRAL OVER SPECTRA
        WRITE (IOUT,*)
        WRITE (IOUT,'(A,A,ES12.4)') ' INTEGRAL OF SPECTRUM ',
     .             UNITOUT,ESTIML(ISPC)%PSPC%SPCS 
        IF (NSIGI_SPC > 0)
     .    WRITE (IOUT,'(A,A,ES12.4)') ' STANDARD DEVIATION   ',
     .                  ' %      ',ESTIML(ISPC)%PSPC%SGMS 
      END DO
 
      RETURN
      END SUBROUTINE EIRENE_OUTSPEC
