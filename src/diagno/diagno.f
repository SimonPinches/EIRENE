cdr Oct 17  :
cdr from W.Zholobenko: add         He emission lines, new options NCHTAL=5
cdr                    analogous to H emission lines,             NCHTAL=2
cdr  itp (select type) of component relevant for LOS (not in use yet)
cdr  Aug. 16:  re LOS option:
cdr            the option described in the manual regarding
cdr            use of emin1, emax1 to identify a particular spectroscopic
cdr            line  (by upper and lower quantum number in H-atom)
cdr            is apparently not available.  Lost from an earier version? to be checked.
C
C  CALLED IN POSTPROCESSING PHASE:
C  CALCULATE A NUMBER OF (ICHORI=1,NCHORI) LINE INTEGRALS ALONG LINES-OF-SIGHT
C  USING THE INPUT DATA OF INPUT BLOCK 12, AND THE VOLUMETRIC INPUT AND OUTPUT
C  TALLIES FROM THE EIRENE RUN.
C  THE LINE INTEGRALS CAN BE TOTAL EMISSIVITIES, OR THEY CAN HAVE
C  HAVE A SPECTRAL RESOLUTION (ENERGY-RESOLUTION), DEPENDING ON THE VALUE
C  OF INPUT FLAG NCHENI(ICHORD).
C  THIS ROUTINE SETS THE ENERGY ARRAYS AND THEN CALLS SUBR. SGNAL FOR
C  EACH LINE OF SIGHT ICHORI, THE SELECTED STRATUM ISTR AND THE SELECTED
C  INDEX OF THE TEST PARTICLE SPECIES ISP
C  AFTER ALL LINES OF SIGHT ARE DONE, THE OUTPUT ROUTINE OUTSIG IS CALLED
C  FOR PRINTING AND PLOTTING.
C
      SUBROUTINE EIRENE_DIAGNO
cdr main program for side on (line of sight) diagnostics, in post processing phase

c  step 1:  prepare arrays for energy (or spectral) resolution (binning)
c  step 2:  call sgnal, to carry out line of sight integration, all bins.
c  step 3:  call outsig, to print (if prspec) and plot (if plspec)
c           energy/spectrally resolved "side-on" data, line integrated.

c  nb    :  during step 2, also spatially resolved (along the line of sight)
c           information can be extracted. This is controlled by the flags
c           plargl, prargl

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CLOGAU
      USE EIRMOD_CPLOT
      USE EIRMOD_COMSIG
      USE EIRMOD_CTRCEI
      USE EIRMOD_COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP) :: ENSAVE(NCHOR,NCHEN)
      REAL(DP) :: EN, EQUOT, FMXENM, ALEMX, ALEMN
      INTEGER :: NSPTP(NCHOR) ! should come via comsig. not ready.
      INTEGER :: J, ISTR, ISP, ITP, NCHNI, ICHORI
      LOGICAL :: PLSAVE,L_CHOR(NCHOR)
C
C  INITIALISE LINE INTEGRATION ROUTINE
C
      NSPTP(1:NCHOR) = 0  !  NOT READY, NOT USED. TYPE OF RELEVANT COMPONENT, see nspspc...
      PLSAVE=PLHST
      PLHST=PLCHOR
C
      NSPNEW(1)=1
      NCHNI=IABS(NCHENI)    !Number of energy/wavelength grid points
      FMXENM=DBLE(NCHNI-1)

      IF (NCHORI.GT.0) THEN
      CALL EIRENE_LEER(2)
      CALL EIRENE_HEADNG
     .     ('DIAGNOSTICS MODULE (LINE-OF-SIGHT-INTEGRATION)',46)
      ENDIF

      DO 100 ICHORI=1,NCHORI

        CALL EIRENE_LEER(1)
        WRITE (IUNOUT,*) 'CHORD NO. ',ICHORI
C
C  SET ENERGY ARRAY IN CASE OF SPECTRALLY RESOLVED LINE INTEGRAL
C
        IF (((NCHTAL(ICHORI).EQ.2).OR.(NCHTAL(ICHORI).EQ.5))
     >    .AND.NCHENI.NE.1) THEN
          WRITE (IUNOUT,*) 'AUTOMATIC CORRECTION  IN DIAGNO'
          WRITE (IUNOUT,*) 'SET NCHENI = 1 (DEFAULT) '
          WRITE (IUNOUT,*) 'BECAUSE NCHTAL(ICHORI)=2 OR '
          WRITE (IUNOUT,*) 'OR NCHTAL(ICHORI)=5 ENCOUNTERED '
          NCHENI=1
        ENDIF

        IF (NCHENI.LT.-1) THEN
C  LOG. ENERGY SCALE
          ALEMN=LOG10(EMIN1(ICHORI))
          ALEMX=LOG10(EMAX1(ICHORI))
          EQUOT=(ALEMX-ALEMN)/FMXENM
          DO 10 J=1,NCHNI
            EN=ALEMN+(J-1)*EQUOT
            ENERGY(J)=10.**EN
            ENSAVE(ICHORI,J)=ENERGY(J)
   10     CONTINUE
        ELSEIF (NCHENI.GT.1) THEN
C  LIN. ENERGY SCALE
          EQUOT=(EMAX1(ICHORI)-EMIN1(ICHORI))/FMXENM
          DO 20 J=1,NCHNI
            ENERGY(J)=EMIN1(ICHORI)+(J-1)*EQUOT
            ENSAVE(ICHORI,J)=ENERGY(J)
   20     CONTINUE
        ELSEIF (NCHENI.EQ.1) THEN
C  NO ENERGY DEPENDENCE
          NCHNI=1
          ENERGY(1)=EMIN1(ICHORI)
          ENSAVE(ICHORI,1)=ENERGY(1)
        ELSEIF (NCHENI.EQ.0) THEN
C  NO CALL TO SUBR. SGNAL
          WRITE (IUNOUT,*) 'NO LOS SIGNALS COMPUTED, BECAUSE NCHENI = 0'
          L_CHOR(ICHORI)=.FALSE.
          GOTO 100
        ENDIF
C
C  CARRY OUT LINE INTEGRATION
C
        ISTR=NSPSTR(ICHORI)  ! Stratum index
        ISP =NSPSPZ(ICHORI)  ! Species index, or no. of contribution. Meaning depends on NCHTAL
        ITP =NSPTP (ICHORI)  ! Type index  (not in use)
C  TENTATIVELY ASSUME: THIS LINE OF SIGHT IS ACTIVE
        L_CHOR(ICHORI)=.TRUE.
        CALL EIRENE_SGNAL(ICHORI,ISTR,ISP,ITP,L_CHOR(ICHORI))
C
  100 CONTINUE
C
C  OUTPUT
C
      IF (NCHNI.NE.0) CALL EIRENE_OUTSIG(ENSAVE,L_CHOR)
C
      PLHST=PLSAVE
      RETURN
      END
