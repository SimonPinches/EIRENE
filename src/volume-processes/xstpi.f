! 30.08.06: data structure for reaction data redefined
! 12.10.06: modcol revised
! 22.11.06: flag for shift of first parameter to rate_coeff introduced
! 01.02.07: do not evaluate rates in vacuum region for IPL (use lgvac(..IPL)
! 20.01.14:  H.4 option for pi rate coefficients (e.g. CR rates: p + H-minus)
c            additional argument PLS, also in calling routines xsecta,xsectm,xsecti
C            additional argument  CHRDF0, also in calling routines
c 23.02.14:  additional argument IPL (was ISP, now ISP is incident test particle)
c            this fixes bug in printout texts(isp)
C            Option 4.3C: now ready,  lgvac for IPL, not for electrons.  corrected !
c 25.03.14:  further corrections. TII instead of plsti(:) array
cdr oct.14:  remove ctrcei, clogau, cspei
cdr oct.14:  synchronize with xstcx started
cdr march 15:  IN VERSION ..-new: nplrpl introduced: flag for bulk particle energy loss/gain
cdr            (was nelrpl, but that is also the flag for electron energy loss/gain)
CDR            eelec in parameter list, before chrdf0
cdr Aug.16  :  minor synchronisation with xstei.f, argument nend removed from prep_poly

C
C
      SUBROUTINE EIRENE_XSTPI (RMASS,IRPI,ISP,IPL,
     .                  IFRST,ISCND,ITHRD,IFRTH,
     .                  EBULK,EHEAVY,EELEC,CHRDF0,ISCDE,IESTM,
     .                  KK,FACTKK,PLS)
C
C       SET UP TABLES (E.G. OF REACTION RATE) FOR PI PROCESSES
C
C   MEANING OF INPUT VARIABLES: SEE XSTCX

C   KK:      COMMON IDENTIFIER FOR PROCESS, USED FOR CROSS-SECTION, RATES,
C                                           STORAGE SAVING MODE ETC...
C   FACTKK:  COMMON SCALING FACTOR FOR PROCESS KK
C   NREAPI(IRPI) = KK DURING MC RUN. THIS ESTABLISHES LINK BETWEEN IRPI AND KK, MUST BE UNIQUE


C  RETURNS:
C    MODCOL(4,...)
C    TABPI3(IRPI,NCELL,...)  1/s per incident test particle
C    EPLPI3(IRPI,NCELL,...) eV/s per incident test particle
C    DEFPI(IRPI)
C    EEFPI(IRPI)
C    IESTPI(IRPI,...)
C

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CCONA
      USE EIRMOD_CGRID
      USE EIRMOD_CZT1
      USE EIRMOD_COMXS
      use EIRMOD_ctrcei, only: trcamd

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: RMASS, EBULK, EHEAVY, EELEC, FACTKK,CHRDF0
      REAL(DP), INTENT(IN) :: PLS(NSTORDR)
      INTEGER, INTENT(IN) :: IRPI, ISP, IPL, IFRST, ISCND, ITHRD, IFRTH,
     .                       ISCDE, IESTM, KK
      REAL(DP) :: CF(9)
      REAL(DP) :: ADD, ADDL, RMTEST, RMBULK, FCTKKL, P2N, TMASS,
     .            ADDT, ADDTL, PMASS,
     .            CHRDIF, COU, ACCMAS, XLFTMAS,
     .            ACCINI, ACCINP, ACCMSM, ACCMSI, ACCMSA, ACCINA,
     .            ACCINM, ACCMSP, ACCINV,
     .            EFLAG, EIRENE_FEHVPI3,
     .            EIRENE_FEELPI1,
     .            EIRENE_RATE_COEFF,
     .            EIRENE_ENERGY_RATE_COEFF,
     .            EI, EA, EN, TB, TII,
     .            FP1(6),FP2(6)
      INTEGER :: NSEPI4, NSEPI5, NEND, J, IO, IA,
     .           ITYP1, ISPZ1, INUM1,
     .           IML, IM, MODC, IIO, IPLTI, IP, IAT,
     .           ICOUNT, IAA, IMM, III, IPP, KREAD, IERR, IMIN, IMAX,
     .           IRAD
      INTEGER, EXTERNAL :: EIRENE_IDEZ
      type(poly_data), pointer :: rp
      type(fit_forms), pointer :: rt

      SAVE

      NREAPI(IRPI) = KK  ! needed for storage saving mode
C
C  SET NON-DEFAULT ION IMPACT COLLISION PROCESS NO. IRPI
C
      IF (IPL.LE.0.OR.IPL.GT.NPLSI) GOTO 990
      IF (MASSP(KK).LE.0.OR.MASST(KK).LE.0) GOTO 992
      RMBULK=RMASSP(IPL)
      RMTEST=RMASS
      IPLTI=MPLSTI(IPL)

      XLFTMAS=RMASS+RMASSP(IPL)

c  In case of more than one secondary (either test or bulk particle)
c  we distribute the total kinetic energy release (KER, eV) per reaction
c  over the secondaries, proportional to their inverse masses.
c  this is stricly correct for reactions with two secondaries.

C ACCUMULATED MASS OF SECONDARIES: ACCMAS (AMU)
C ACCUMULATED INVERSE MASS OF SECONDARIES: ACCINV (1/AMU)
      ACCMAS=0.D0
      ACCMSA=0.D0
      ACCMSM=0.D0
      ACCMSI=0.D0
      ACCMSP=0.D0
      ACCINV=0.D0
      ACCINA=0.D0
      ACCINM=0.D0
      ACCINI=0.D0
      ACCINP=0.D0
C
C  LOOP OVER UP TO 4 SECONDARIES
      DO ICOUNT=1,4
        IF (ICOUNT == 1) THEN
C  SECONDARY INDEX, FIRST SECONDARY
          ITYP1=EIRENE_IDEZ(IFRST,1,3)
          INUM1=EIRENE_IDEZ(IFRST,2,3)
          ISPZ1=EIRENE_IDEZ(IFRST,3,3)
        ELSE IF (ICOUNT == 2) THEN
C  SECONDARY INDEX, SECOND SECONDARY (if any)
          IF (ISCND == 0) EXIT
          ITYP1=EIRENE_IDEZ(ISCND,1,3)
          INUM1=EIRENE_IDEZ(ISCND,2,3)
          ISPZ1=EIRENE_IDEZ(ISCND,3,3)
        ELSE IF (ICOUNT == 3) THEN
C  SECONDARY INDEX, THIRD SECONDARY (if any)
          IF (ITHRD == 0) EXIT
          ITYP1=EIRENE_IDEZ(ITHRD,1,3)
          INUM1=EIRENE_IDEZ(ITHRD,2,3)
          ISPZ1=EIRENE_IDEZ(ITHRD,3,3)
        ELSE IF (ICOUNT == 4) THEN
C  SECONDARY INDEX, FOURTH SECONDARY (if any)
          IF (IFRTH == 0) EXIT
          ITYP1=EIRENE_IDEZ(IFRTH,1,3)
          INUM1=EIRENE_IDEZ(IFRTH,2,3)
          ISPZ1=EIRENE_IDEZ(IFRTH,3,3)
        END IF

        IF ((ISPZ1 < 1) .OR. (ISPZ1 > MAXSPC(ITYP1))) GOTO 994

!  ACCMAS: accumulated mass of all secondaries (all types)
!  ACCINV: accumulated inverse mass of all secondaries (all types)

!  ACCMSA: accumulated mass of ATOMIC secondaries (type ITYP=1)
!  ACCINA: accumulated inverse mass of ATOMIC secondaries (type ITYP=1)
!  analogously for molecule, test ion and bulk secondaries

        IF (ITYP1.EQ.1) THEN
          IAT=ISPZ1
          IAA=NSPH+IAT
          PATPI(IRPI,IAT)=PATPI(IRPI,IAT)+INUM1
          P2NP(IRPI,IAA)=P2NP(IRPI,IAA)+INUM1
          ACCMAS=ACCMAS+INUM1*RMASSA(IAT)
          ACCMSA=ACCMSA+INUM1*RMASSA(IAT)
          ACCINV=ACCINV+INUM1/RMASSA(IAT)
          ACCINA=ACCINA+INUM1/RMASSA(IAT)
          EATPI(IRPI,IAT,1)=RMASSA(IAT)
          EATPI(IRPI,IAT,2)=1./RMASSA(IAT)
        ELSE IF (ITYP1.EQ.2) THEN
          IML=ISPZ1
          IMM=NSPA+IML
          PMLPI(IRPI,IML)=PMLPI(IRPI,IML)+INUM1
          P2NP(IRPI,IMM)=P2NP(IRPI,IMM)+INUM1
          ACCMAS=ACCMAS+INUM1*RMASSM(IML)
          ACCMSM=ACCMSM+INUM1*RMASSM(IML)
          ACCINV=ACCINV+INUM1/RMASSM(IML)
          ACCINM=ACCINM+INUM1/RMASSM(IML)
          EMLPI(IRPI,IML,1)=RMASSM(IML)
          EMLPI(IRPI,IML,2)=1./RMASSM(IML)
        ELSE IF (ITYP1.EQ.3) THEN
          IIO=ISPZ1
          III=NSPAM+IIO
          PIOPI(IRPI,IIO)=PIOPI(IRPI,IIO)+INUM1
          P2NP(IRPI,III)=P2NP(IRPI,III)+INUM1
          ACCMAS=ACCMAS+INUM1*RMASSI(IIO)
          ACCMSI=ACCMSI+INUM1*RMASSI(IIO)
          ACCINV=ACCINV+INUM1/RMASSI(IIO)
          ACCINI=ACCINI+INUM1/RMASSI(IIO)
          EIOPI(IRPI,IIO,1)=RMASSI(IIO)
          EIOPI(IRPI,IIO,2)=1./RMASSI(IIO)
        ELSE IF (ITYP1.EQ.4) THEN
          IPP=ISPZ1
          PPLPI(IRPI,IPP)=PPLPI(IRPI,IPP)+INUM1
C  NO P2NP HERE; P2NP ONLY COUNTS TEST PARTICLE SECONDARIES
          ACCMAS=ACCMAS+INUM1*RMASSP(IPP)
          ACCMSP=ACCMSP+INUM1*RMASSP(IPP)
          ACCINV=ACCINV+INUM1/RMASSP(IPP)
          ACCINP=ACCINP+INUM1/RMASSP(IPP)
          EPLPI(IRPI,IPP,1)=RMASSP(IPP)
          EPLPI(IRPI,IPP,2)=1./RMASSP(IPP)
        END IF
      END DO

      IF (ABS(ACCMAS-XLFTMAS).GT.1.D-10) THEN
        WRITE (IUNOUT,*) 'MESSAGE FROM XSTPI.F: '
        WRITE (IUNOUT,*) 'FOR INCIDENT TEST SPECIES ',TEXTS(ISP)
        WRITE (IUNOUT,*) 'FOR INCIDENT BULK SPECIES ',TEXTS(NSPAMI+IPL)
        WRITE (iunout,*)
     .    'MASS CONSERVATION VIOLATED FOR REACT. KK= ',KK
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
C
C  FRACTIONS OF KINETIC ENERGY RELEASE; TO SECONDARY ATOMS
      DO IAT=1,NATMI
        EATPI(IRPI,IAT,1)=EATPI(IRPI,IAT,1)/ACCMAS
        EATPI(IRPI,IAT,2)=EATPI(IRPI,IAT,2)/ACCINV
      ENDDO
      EATPI(IRPI,0,1)=ACCMSA/ACCMAS
      EATPI(IRPI,0,2)=ACCINA/ACCINV
C  FRACTIONS OF KINETIC ENERGY RELEASE; TO SECONDARY MOLECULES
      DO IML=1,NMOLI
        EMLPI(IRPI,IML,1)=EMLPI(IRPI,IML,1)/ACCMAS
        EMLPI(IRPI,IML,2)=EMLPI(IRPI,IML,2)/ACCINV
      ENDDO
      EMLPI(IRPI,0,1)=ACCMSM/ACCMAS
      EMLPI(IRPI,0,2)=ACCINM/ACCINV
C  FRACTIONS OF KINETIC ENERGY RELEASE; TO SECONDARY TEST IONS
      DO IIO=1,NIONI
        EIOPI(IRPI,IIO,1)=EIOPI(IRPI,IIO,1)/ACCMAS
        EIOPI(IRPI,IIO,2)=EIOPI(IRPI,IIO,2)/ACCINV
      ENDDO
      EIOPI(IRPI,0,1)=ACCMSI/ACCMAS
      EIOPI(IRPI,0,2)=ACCINI/ACCINV
C  FRACTIONS OF KINETIC ENERGY RELEASE; TO SECONDARY BULK IONS
      DO IPP=1,NPLSI
        EPLPI(IRPI,IPP,1)=EPLPI(IRPI,IPP,1)/ACCMAS
        EPLPI(IRPI,IPP,2)=EPLPI(IRPI,IPP,2)/ACCINV
      ENDDO
      EPLPI(IRPI,0,1)=ACCMSP/ACCMAS
      EPLPI(IRPI,0,2)=ACCINP/ACCINV
C
      CHRDIF=CHRDF0

      CHRDIF = CHRDIF-NCHRGP(IPL)
      DO 133 IIO=1,NIONI
        CHRDIF=CHRDIF+PIOPI(IRPI,IIO)*NCHRGI(IIO)
  133 CONTINUE
      DO 134 IP=1,NPLSI
        CHRDIF=CHRDIF+PPLPI(IRPI,IP)*NCHRGP(IP)
  134 CONTINUE
      PELPI(IRPI)=PELPI(IRPI)+CHRDIF
C
C
C  TARGET MASS IN <SIGMA*V> FORMULA: MAXW. BULK PARTICLE
C  (= PROJECTILE MASS IN CROSS-SECTION MEASUREMENT: TARGET AT REST)
      PMASS=MASSP(KK)*PMASSA
C  PROJECTILE MASS IN <SIGMA*V> FORMULA: MONOENERG. TEST PARTICLE
C  (= TARGET PARTICLE IN CROSS-SECTION MEASUREMENT; TARGET AT REST)
      TMASS=MASST(KK)*PMASSA
C
      ADDT=PMASS/RMASSP(IPL)
      ADDTL=LOG(ADDT)
      ADDPI(IRPI,IPL) = ADDTL
C
C CROSS-SECTION (E-LAB) AVAILABLE ?
      IF (EIRENE_IDEZ(MODCLF(KK),2,5).EQ.1) THEN
        MODCOL(4,1,IRPI)=KK
C  TENTATIVLEY ASSUME: SIGMA * V_EFF MODEL FOR RATE COEFFICIENT
        MODCOL(4,2,IRPI)=3
      ENDIF

C..................................................................
C 2. RATE COEFFICIENT  (CM**3/S) * TARGET DENSITY (CM**-3)
C..................................................................

      MODC=EIRENE_IDEZ(MODCLF(KK),3,5)

      IF (MODC.GE.1.AND.MODC.LE.2) THEN

        MODCOL(4,2,IRPI)=MODC
C  2.B)
        IF (MODC.EQ.1) NEND=1   ! rate coeff for (FIXED E0, e.g. E0=0.0, TI)
C  2.C)
        IF (MODC.EQ.2) NEND=NSTORDT ! rate coeff vs. (E0, TI) NEND=9 HERE
C   STORAGE SAVING MODE ?
        IF (NSTORDR >= NRAD) THEN
C   NO, NSTORDT=9 HERE

C  2.B) RATE COEFFICIENT(TI, FIXED E0, E.G. E0=0)
          IF (MODC.EQ.1) THEN
C           NEND=1
            DO 145 J=1,NSBOX
              IF (LGVAC(J,IPL)) CYCLE
              TII=TIINL(IPLTI,J)+ADDTL
              COU = EIRENE_RATE_COEFF(KK,J,TII,0._DP,.TRUE.,0)
              TABPI3(IRPI,J,1)=COU*DIIN(IPL,J)*FACTKK
  145       CONTINUE
          ELSEIF (MODC.EQ.2) THEN
C           NEND=9
C  2.C) RATE COEFFICIENT(TI,EBEAM)
C           NEND=9
            FCTKKL=LOG(FACTKK)
            rt => reacdat(kk)%rtc
            fp1(1:3) = rt%fp1l
            fp1(4:6) = rt%fp1r
            fp2(1:3) = rt%fp2b
            fp2(4:6) = rt%fp2t
            DO J=1,NSBOX
              IF (LGVAC(J,IPL)) CYCLE
              TII=TIINL(IPLTI,J)+ADDTL
              tii = max(-2.3_dp,tii) ! this is another cut-off, at TIIN <=0.1 eV rather than at TVAC = 0.02 ev
c old
c old         CALL EIRENE_PREP_RTCS (KK,3,TII,CF)
c old
              rp => reacdat(KK)%rtc%poly
              call EIRENE_dbl_poly (rp%dblpol,tii,0._dp,cou,cf,
     .             rt%rc1min, rt%rc1max, fp1, rt%jfex1mn, rt%jfex1mx,
     .             rt%rc2min, rt%rc2max, fp2, rt%jfex2mn, rt%jfex2mx,
     .             trcamd)

              TABPI3(IRPI,J,1:9) = CF(1:9)
              TABPI3(IRPI,J,1)=TABPI3(IRPI,J,1)+DIINL(IPL,J)+FCTKKL
            END DO
          END IF  ! MODC=1,2
        ELSE ! NOT SUFFICIENT STORAGE ON TABPI3
C  STORAGE SAVE MODE NOT READY FOR THIS OPTION ??
          GOTO 995

        ENDIF
      ELSEIF (EIRENE_IDEZ(MODCLF(KK),3,5).EQ.3) THEN
C  2.D) RATE COEFFICIENT(TI=TE, NE=NI ?, E0 FIXED, E.G. E0=0.)
C       IF (MODC.EQ.3) NEND=1  rate coeff vs. (N, T), NEND NOT NEEDED

        MODCOL(4,2,IRPI)=1 !  indicate: rate coefficient as fct. of local plasma conditions only
        FCTKKL=LOG(FACTKK)
        IF (NSTORDR >= NRAD) THEN

          DO J=1,NSBOX
            IF (LGVAC(J,IPL)) CYCLE
            COU = EIRENE_RATE_COEFF(KK,J,TEINL(J),PLS(J),.FALSE.,1)
            TB = COU + FCTKKL
            IF (IFTFLG(KK,2) < 100) TB = TB + DIINL(IPL,J)
            TB=MAX(-100._DP,TB)
            TABPI3(IRPI,J,1)=EXP(TB)
          END DO
        ELSE  ! ??
C  WHAT DO WE DO IN CASE NSTORDR < NRAD  ?
          write (iunout,*) 'storage save mode not available yet for PI'
          write (iunout,*) 'in case modc=3  (n,T-dependence).'
          GOTO 995
        ENDIF

      ELSE
C  NO RATE COEFFICIENT. IS THERE A CROSS-SECTION AT LEAST?
        IF (MODCOL(4,2,IRPI).NE.3) GOTO 996
      ENDIF

      FACRPI(IRPI,1) = FACTKK
      FACRPI(IRPI,2) = LOG(FACTKK)

      DEFPI(IRPI)=LOG(CVELI2*PMASS)
      EEFPI(IRPI)=LOG(CVELI2*TMASS)
C
C  3. BULK PARTICLE MOMENTUM LOSS RATE
C
C  4. BULK PARTICLE ENERGY LOSS RATE
C  4.1. HEAVY BULK PARTICLE ENERGY LOSS RATE
C
C  SET ENERGY LOSS RATE OF IMPACTING ION
C
      NSEPI4=EIRENE_IDEZ(ISCDE,4,5)
      IF (NSEPI4.EQ.0) THEN
C  4.1A)  ENERGY LOSS RATE OF IMP. BULK PARTICLE = CONST.*RATE COEFF.
C        SAMPLE COLLIDING ION FROM DRIFTING MONOENERGETIC ISOTROPIC DISTRIBUTION
c        WITH WEIGHTING/REJECTION
        IF (EBULK.LE.0.D0) THEN
          IF (NSTORDR >= NRAD) THEN
            EPLPI3(IRPI,1:NSBOX,1)=1.5*TIIN(IPLTI,1:NSBOX)
            IF (LEDRIFT) EPLPI3(IRPI,1:NSBOX,1)=
     .                   EPLPI3(IRPI,1:NSBOX,1)+EDRIFT(IPL,1:NSBOX)
            NELRPI(IRPI) = -3
          ELSE
            NELRPI(IRPI) = -3
          END IF
        ELSE ! EBULK GT.0
          IF (NSTORDR >= NRAD) THEN
            EPLPI3(IRPI,1:NSBOX,1)=EBULK
            IF (LEDRIFT) EPLPI3(IRPI,1:NSBOX,1)=
     .                   EPLPI3(IRPI,1:NSBOX,1)+EDRIFT(IPL,1:NSBOX)
            NELRPI(IRPI) = -2
          ELSE
            NELRPI(IRPI) = -2
            EPLPI3(IRPI,1,1)=EBULK
          END IF
C       ELSE
CDR   ERROR: EBULK < 0 IS NOT FORESEEN
        ENDIF
        MODCOL(4,4,IRPI)=3
      ELSEIF (NSEPI4.EQ.1) THEN
C  4.1B) ENERGY LOSS RATE OF IMP. ION = (1.5*TI+EDRIFT)* RATE COEFF.
C        SAMPLE COLLIDING ION FROM DRIFTING MAXWELLIAN
        IF (EBULK.LE.0.D0) THEN
          IF (NSTORDR >= NRAD) THEN
            EPLPI3(IRPI,1:NSBOX,1)=1.5*TIIN(IPLTI,1:NSBOX)
            IF (LEDRIFT) EPLPI3(IRPI,1:NSBOX,1)=
     .                   EPLPI3(IRPI,1:NSBOX,1)+EDRIFT(IPL,1:NSBOX)
            NELRPI(IRPI) = -3
          ELSE
            NELRPI(IRPI) = -3
          END IF
        ELSE ! EBULK GT.0
          WRITE (iunout,*) 'WARNING FROM SUBR. XSTPI: IRPI ', IRPI
          WRITE (iunout,*) 'MODIFIED TREATMENT OF BULK ION IMPACT '
          WRITE (iunout,*) 'SAMPLE FROM MAXWELLIAN WITH T = ',EBULK/1.5
          WRITE (iunout,*) 'RATHER THAN WITH T = TIIN '
          WRITE (iunout,*) 'NOT FULLY IMPLEMENTED (VELOPI) '
          CALL EIRENE_LEER(1)
          IF (NSTORDR >= NRAD) THEN
            EPLPI3(IRPI,1:NSBOX,1)=EBULK
            IF (LEDRIFT) EPLPI3(IRPI,1:NSBOX,1)=
     .                   EPLPI3(IRPI,1:NSBOX,1)+EDRIFT(IPL,1:NSBOX)
            NELRPI(IRPI) = -2
          ELSE
            NELRPI(IRPI) = -2
            EPLPI3(IRPI,1,1)=EBULK
          END IF
C       ELSE
CDR   ERROR: EBULK < 0 IS NOT FORESEEN
        ENDIF
        MODCOL(4,4,IRPI)=1
C     ELSEIF (NSEPI4.EQ.2) THEN
C  use i-integral expressions. to be written
      ELSEIF (NSEPI4.EQ.3) THEN
C  4.1C)  ENERGY LOSS RATE OF IMP. ION = EN.-WEIGHTED RATE
C       SAMPLE COLLIDING ION FROM DRIFTING MAXWELLIAN, WITH WEIGHTING/REJECTION
        KREAD=EBULK
        IF (KREAD.EQ.0) THEN
c  data for mean ion energy loss are not available
c  use collision estimator for energy balance
          IF (EIRENE_IDEZ(IESTM,3,3).NE.1) THEN
            WRITE (iunout,*)
     .        'COLLISION ESTIMATOR ENFORCED FOR ION ENERGY '
            WRITE (iunout,*) 'IN PI COLLISION IRPI= ',IRPI
            WRITE (iunout,*) 'BECAUSE NO ENERGY-WEIGHTED RATE AVAILABLE'
          ENDIF
          IESTPI(IRPI,3)=1
          MODCOL(4,4,IRPI)=2
        ELSE
C  ION ENERGY-AVERAGED RATE AVAILABLE AS REACTION NO. "KREAD"
        NELRPI(IRPI) = KREAD
        MODC=EIRENE_IDEZ(MODCLF(KREAD),5,5)
        IF (MODC.GE.1.AND.MODC.LE.2) THEN
          MODCOL(4,4,IRPI)=MODC
          IF (MODC.EQ.1) NEND=1
          IF (MODC.EQ.2) NEND=NSTORDT
C  STORAGE SAVING MODE ?
          IF (NSTORDR >= NRAD) THEN
C  NO
C           NSTORDT=9 HERE

            IF (MODC.EQ.1) THEN
C             NEND=1
C  ENERGY RATE COEFFICIENT(TI, EBEAM=0)
              ADD=FACTKK/ADDT
              DO 254 J=1,NSBOX
                IF (LGVAC(J,IPL)) CYCLE
                TII=TIINL(IPLTI,J)+ADDTL
                EPLPI3(IRPI,J,1)=EIRENE_ENERGY_RATE_COEFF
     .                          (KREAD,J,TII,
     .                           0._DP,.FALSE.,0)*DIIN(IPL,J)*ADD
  254         CONTINUE
            ELSEIF (MODC.EQ.2) THEN
C             NEND=9
C  ENERGY RATE COEFFICIENT(TI,EBEAM)
              ADDL=LOG(FACTKK)-ADDTL
              rt => reacdat(kread)%rtcew
              fp1(1:3) = rt%fp1l
              fp1(4:6) = rt%fp1r
              fp2(1:3) = rt%fp2b
              fp2(4:6) = rt%fp2t
              DO 257 J=1,NSBOX
                IF (LGVAC(J,IPL)) CYCLE
                TII=TIINL(IPLTI,J)+ADDTL
                tii = max(-2.3_dp,tii)
c old
c old           CALL EIRENE_PREP_RTCS (KREAD,5,TII,CF)
c old
                rp => reacdat(KREAD)%rtcew%poly
                call EIRENE_dbl_poly (rp%dblpol,tii,0._dp,cou,cf,
     .               rt%rc1min, rt%rc1max, fp1, rt%jfex1mn, rt%jfex1mx,
     .               rt%rc2min, rt%rc2max, fp2, rt%jfex2mn, rt%jfex2mx,
     .               trcamd)

                EPLPI3(IRPI,J,1:9) = CF(1:9)
                EPLPI3(IRPI,J,1) = EPLPI3(IRPI,J,1)+DIINL(IPL,J)+ADDL
  257         CONTINUE
            ENDIF

          ELSE  ! STORAGE SAVING MODE, no pre-defined tallies eplpi3
            IF (MODC.EQ.1) THEN
              ADD=FACTKK/ADDT
              EPLPI3(IRPI,1,1)=ADD   !  ????
            ELSEIF (MODC.EQ.2) THEN
              ADDL=LOG(FACTKK)-ADDTL
              FACRPI(IRPI,1) = EXP(ADDL)
              FACRPI(IRPI,2) = ADDL
            END IF
          END IF
        ENDIF
        ENDIF
      ELSE
        WRITE (iunout,*) 'NSEPI4 ILL-DEFINED IN XSTPI '
        WRITE (iunout,*) 'check parameter ISCDE for process irpi ',irpi
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
C
C  4.2. BULK ELECTRON ENERGY LOSS RATE
C
C  SET NET ENERGY LOSS RATE OF ELECTRON (IF ANY INVOLVED)
      NSEPI5=EIRENE_IDEZ(ISCDE,5,5)
      IF (NSEPI5.EQ.0) THEN
C  4.A1) ENERGY LOSS RATE OF IMP. ELECTRON = CONST.*RATE COEFF.
        IF (NSTORDR >= NRAD) THEN
          DO J=1,NSBOX
            EELPI1(IRPI,J)=EELEC
          ENDDO
          NELRPI(IRPI)=0
        ELSE
          NELRPI(IRPI)=0
          EELPI1(IRPI,1)=EELEC
        END IF
        MODCOL(4,4,IRPI)=1
      ELSE
        WRITE (iunout,*) 'NSEPI5 ILL-DEFINED IN XSTPI '
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
C
C  4.3. HEAVY PARTICLE ENERGY GAIN RATE
C
      EFLAG=EIRENE_IDEZ(ISCDE,3,5)
      IF (EFLAG.EQ.0) THEN
C  4.3A)  RATE = CONST.*RATE COEFF.
        IF (NSTORDR >= NRAD) THEN
          DO 201 J=1,NSBOX
            EHVPI3(IRPI,J,1)=EHEAVY
  201     CONTINUE
          NHVRPI(IRPI)=0
        ELSE
          NHVRPI(IRPI)=0
          EHVPI3(IRPI,1,1)=EHEAVY
        END IF

C     ELSEIF (EFLAG.EQ.1) THEN
C        NOT A VALID OPTION

      ELSEIF (EFLAG.EQ.3) THEN
C  4.3C)  SECONDARY HEAVY ENERGY GAIN RATE = EN.-WEIGHTED RATE(TI)
        KREAD=EHEAVY
        MODC=EIRENE_IDEZ(MODCLF(KREAD),5,5)
        IF (MODC.EQ.1) THEN
          IF (NSTORDR >= NRAD) THEN
            DO 202 J=1,NSBOX
              IF (LGVAC(J,IPL)) CYCLE
                TII=TIINL(IPLTI,J)+ADDTL
                EHVPI3(IRPI,J,1)=EIRENE_ENERGY_RATE_COEFF(KREAD,J,
     .                           TII,0._DP,.TRUE.,0)*
     .          DIIN(IPL,J)*FACTKK/(TABPI3(IRPI,J,1)+EPS60)
  202       CONTINUE
            NHVRPI(IRPI)=KREAD
          ELSE
            NHVRPI(IRPI)=KREAD
          END IF
        ELSE
          WRITE (iunout,*) 'INVALID OPTION IN XSTPI: MODC, EFLAG '
          WRITE (iunout,*) MODC, EFLAG
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
        FACRPI(IRPI,1)=FACTKK
        FACRPI(IRPI,2)=LOG(FACTKK)
      ELSE
        IERR=2
        GOTO 997
      ENDIF
C
C  ESTIMATOR FOR CONTRIBUTION TO COLLISION RATES FROM THIS REACTION
      IESTPI(IRPI,1)=EIRENE_IDEZ(IESTM,1,3)
      IESTPI(IRPI,2)=EIRENE_IDEZ(IESTM,2,3)
      IF (IESTPI(IRPI,3).EQ.0) IESTPI(IRPI,3)=EIRENE_IDEZ(IESTM,3,3)
C
      IF (IESTPI(IRPI,1).NE.0) THEN
        WRITE (iunout,*)
     .    'WARNING: COLL.EST NOT AVAILABLE FOR PART. BALANCE '
        WRITE (iunout,*) 'IRPI = ',IRPI
        WRITE (iunout,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
        IESTPI(IRPI,1)=0
        CALL EIRENE_LEER(1)
      ENDIF
      IF (IESTPI(IRPI,2).NE.0) THEN
        WRITE (iunout,*)
     .    'WARNING: COLL.EST NOT AVAILABLE FOR MOM. BALANCE '
        WRITE (iunout,*) 'IRPI = ',IRPI
        WRITE (iunout,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
        IESTPI(IRPI,2)=0
        CALL EIRENE_LEER(1)
      ENDIF

      IF (IESTPI(IRPI,3).NE.0) THEN
        WRITE (iunout,*)
     .    'WARNING: COLL.EST NOT AVAILABLE FOR EN. BALANCE '
        WRITE (iunout,*) 'IRPI = ',IRPI
        WRITE (iunout,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
        CALL EIRENE_LEER(1)
        IESTPI(IRPI,3)=0
      ENDIF
      RETURN
C
C-----------------------------------------------------------------------
C
      ENTRY EIRENE_XSTPI_1(IRPI)
C
C  SET TOTAL NUMBER OF SECONDARIES BY TYPE OF SECONDARY: P..PI(IRPI,0)
C  AND CUMULATE TOTAL: P2NP
C  CONVERT SECONDARY SPECIES DISTRIBUTION P2NP(IRPI)  INTO
C  CUMULATIVE DISTRIBUTION (NOT YET NORMALIZED, THIS IS DONE BELOW).


C  ATOM SECONDARIES
      DO 510 IAT=1,NATMI
        IA=NSPH+IAT
        PATPI(IRPI,0)=PATPI(IRPI,0)+
     +                      PATPI(IRPI,IAT)
        P2NP(IRPI,IA)=P2NP(IRPI,IA-1)+
     +                      P2NP(IRPI,IA)
  510 CONTINUE
C  MOLECULE SECONDARIES
      DO 520 IML=1,NMOLI
        IM=NSPA+IML
        PMLPI(IRPI,0)=PMLPI(IRPI,0)+
     +                      PMLPI(IRPI,IML)
        P2NP(IRPI,IM)=P2NP(IRPI,IM-1)+
     +                      P2NP(IRPI,IM)
  520 CONTINUE
C  TEST ION SECONDARIES
      DO 530 IIO=1,NIONI
        IO=NSPAM+IIO
        PIOPI(IRPI,0)=PIOPI(IRPI,0)+
     +                      PIOPI(IRPI,IIO)
        P2NP(IRPI,IO)=P2NP(IRPI,IO-1)+
     +                      P2NP(IRPI,IO)
  530 CONTINUE
C  BULK SECONDARIES (NOT ON P2NP)
      DO 540 IPP=1,NPLSI
        PPLPI(IRPI,0)=PPLPI(IRPI,0)+
     +                      PPLPI(IRPI,IPP)
  540 CONTINUE
C
C  TOTAL NUMBER OF TEST PARTICLE SECONDARIES
      P2NPI(IRPI)=PATPI(IRPI,0)+PMLPI(IRPI,0)+
     .            PIOPI(IRPI,0)

C  FINALY: NORMALIZE SECONDARY TEST PARTICLE SPECIES DISTRIBUTION P2NP
C          SUCH THAT IT BECOMES A CUMULATIVE SAMPLING DISTRIBUTION
C          FOR TEST PARTICLE SECONDARIES
C          NORMALIZATION DOES NOT EXTEND OVER SECONDARY BULK PARTICLES
      P2N=P2NP(IRPI,NSPAMI)
CDR  no photon secondaries here, otherwise loop would start with ispz1=1
      DO 550 ISPZ1=NSPH+1,NSPAMI
        IF (P2N.GT.0.D0)
     .  P2NP(IRPI,ISPZ1)=P2NP(IRPI,ISPZ1)/P2N
  550 CONTINUE
C
      RETURN
C
C-----------------------------------------------------------------------
C

      ENTRY EIRENE_XSTPI_2(IRPI,IPL)
C
      CALL EIRENE_LEER(2)
      WRITE (iunout,*) 'GENERAL ION IMPACT REACTION NO. IRPI= ', IRPI
      CALL EIRENE_LEER(1)
      WRITE (iunout,*) 'HEAVY PARTICLE COLLISION WITH BULK IONS IPLS:'
      WRITE (iunout,*) 'IPLS= ',TEXTS(NSPAMI+IPL)
      CALL EIRENE_LEER(1)

      WRITE (iunout,*) 'BACKGROUND SECONDARIES:'

C  ARE SECONDARY ELECTRONS INVOLVED?
      IF (PELPI(IRPI).NE.0.0) THEN
        EI=1.D30
        EA=-1.D30
        imin=0
        imax=0
        DO 875 IRAD=1,NSBOX
          IF (LGVAC(IRAD,IPL)) GOTO 875
          IF (NSTORDR >= NRAD) THEN
            EN=EELPI1(IRPI,IRAD)
          ELSE
            EN=EIRENE_FEELPI1(IRPI,IRAD)
          END IF
          if (en < ei) imin=irad
          if (en > ea) imax=irad
          EI=MIN(EI,EN)
          EA=MAX(EA,EN)
  875   CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10.OR.EI.EQ.1.D30) THEN
          WRITE (iunout,*) 'ELECTRONS: PELPI, CONSTANT ENERGY: EEL'
          WRITE (iunout,'(1X,A8,2(1PE12.4))') 'EL      ',PELPI(IRPI),EI
        ELSEIF (EI.NE.1.D30) THEN
          WRITE (iunout,*)
     .      'ELECTRONS: PELPI, ENERGY RANGE: EEL_MIN,EEL_MAX'
          WRITE (iunout,'(1X,A8,3(1PE12.4))') 'EL      ',
     .                   PELPI(IRPI),EI,EA
        ENDIF
      ENDIF
c     write (iunout,*) ' imin = ', imin, ' imax = ',imax
C
      EI=1.D30
      EA=-1.D30
      DO 876 IRAD=1,NSBOX
        IF (LGVAC(IRAD,IPL)) GOTO 876
        IF (NSTORDR >= NRAD) THEN
          EN=EHVPI3(IRPI,IRAD,1)
        ELSE
          EN=EIRENE_FEHVPI3(IRPI,IRAD)
        END IF
        EI=MIN(EI,EN)
        EA=MAX(EA,EN)
  876 CONTINUE

      IF (PPLPI(IRPI,0).GT.0.D0) THEN
        WRITE (iunout,*) 'BULK IONS: PPLPI, INCIDENT BULK SUBTRACTED '
        DO 874 IPP=1,NPLSI
          IP=NSPAMI+IPP
          IF (IPP.EQ.IPL) THEN
C  SUBTRACT ONE, BECAUSE INCIDENT BULK IS LOST
            WRITE (iunout,'(1X,A8,1PE12.4)') TEXTS(IP),PPLPI(IRPI,IPP)-1
          ELSEIF (PPLPI(IRPI,IPP).NE.0.D0) THEN
            WRITE (iunout,'(1X,A8,1PE12.4)') TEXTS(IP),PPLPI(IRPI,IPP)
          ENDIF
  874   CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10.OR.EI.EQ.1.D30) THEN
          WRITE (iunout,*) 'ENERGY: EPLPI '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4)') EPLPI(IRPI,0,1),
     .                                 ' * E0 + ',EPLPI(IRPI,0,2)*EI
C  PROBABLY INCORRECT: COM IS NOT EQ. E0 IN CASE OF PI, ONLY IN CASE OF EI
        ELSEIF (EI.NE.1.D30) THEN
          WRITE (iunout,*) 'ENERGY: EPLPI '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4,A10)') EPLPI(IRPI,0,1),
     .                                 ' * E0 + ',EPLPI(IRPI,0,2),
     .                                 ' * EHEAVY '
C  PROBABLY INCORRECT: COM IS NOT EQ. E0 IN CASE OF PI, ONLY IN CASE OF EI
          WRITE (iunout,*) 'ENERGY RANGE: EHEAVY_MIN, EHEAVY_MAX'
          WRITE (iunout,'(1X,2(1PE12.4))') EI,EA
        ENDIF
      ELSE
        WRITE (iunout,*) 'BULK IONS: NONE'
      ENDIF
      CALL EIRENE_LEER(1)
C
      WRITE (iunout,*) 'TEST PARTICLE SECONDARIES:'
      IF (P2NPI(IRPI).EQ.0.D0) THEN
        WRITE (iunout,*) 'NONE'
        CALL EIRENE_LEER(1)
        GOTO 880
      ENDIF
C
      IF (PATPI(IRPI,0).GT.0.D0) THEN
        WRITE (iunout,*) 'ATOMS    : PATPI '
        DO 871 IAT=1,NATMI
          IA=NSPH+IAT
          IF (PATPI(IRPI,IAT).NE.0.D0)
     .    WRITE (iunout,'(1X,A8,1PE12.4)') TEXTS(IA),PATPI(IRPI,IAT)
  871   CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
          WRITE (iunout,*) 'ENERGY: EATPI '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4)') EATPI(IRPI,0,1),
     .                                 ' * E0 + ',EATPI(IRPI,0,2)*EI
        ELSE
          WRITE (iunout,*) 'ENERGY: EATPI '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4,A10)') EATPI(IRPI,0,1),
     .                                 ' * E0 + ',EATPI(IRPI,0,2),
     .                                 ' * EHEAVY'
          WRITE (iunout,*) 'ENERGY RANGE: EHEAVY_MIN, EHEAVY_MAX'
          WRITE (iunout,'(1X,2(1PE12.4))') EI,EA
        ENDIF
      ENDIF
      IF (PMLPI(IRPI,0).GT.0.D0) THEN
        WRITE (iunout,*) 'MOLECULES: PMLPI '
        DO 872 IML=1,NMOLI
          IM=NSPA+IML
          IF (PMLPI(IRPI,IML).NE.0.D0)
     .    WRITE (iunout,'(1X,A8,1PE12.4)') TEXTS(IM),PMLPI(IRPI,IML)
  872   CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
          WRITE (iunout,*) 'ENERGY: EMLPI '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4)') EMLPI(IRPI,0,1),
     .                                 ' * E0 + ',EMLPI(IRPI,0,2)*EI
        ELSE
          WRITE (iunout,*) 'ENERGY: EMLPI '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4,A10)') EMLPI(IRPI,0,1),
     .                                 ' * E0 + ',EMLPI(IRPI,0,2),
     .                                 ' * EHEAVY'
          WRITE (iunout,*) 'ENERGY RANGE: EHEAVY_MIN, EHEAVY_MAX'
          WRITE (iunout,'(1X,2(1PE12.4))') EI,EA
        ENDIF
      ENDIF
      IF (PIOPI(IRPI,0).GT.0.D0) THEN
        WRITE (iunout,*) 'TEST IONS: PIOPI '
        DO 873 IIO=1,NIONI
          IO=NSPAM+IIO
          IF (PIOPI(IRPI,IIO).NE.0.D0)
     .    WRITE (iunout,'(1X,A8,1PE12.4)') TEXTS(IO),PIOPI(IRPI,IIO)
  873   CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
          WRITE (iunout,*) 'ENERGY: EIOPI '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4)') EIOPI(IRPI,0,1),
     .                                 ' * E0 + ',EIOPI(IRPI,0,2)*EI
        ELSE
          WRITE (iunout,*) 'ENERGY: EIOPI '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4,A10)') EIOPI(IRPI,0,1),
     .                                 ' * E0 + ',EIOPI(IRPI,0,2),
     .                                 ' * EHEAVY'
          WRITE (iunout,*) 'ENERGY RANGE: EHEAVY_MIN, EHEAVY_MAX'
          WRITE (iunout,'(1X,2(1PE12.4))') EI,EA
        ENDIF
      ENDIF

  880 CONTINUE

      CALL EIRENE_LEER(1)

      IF (IESTPI(IRPI,1).NE.0)
     .   WRITE (IUNOUT,*) 'COLLISION ESTIMATOR FOR PART. BALANCE '
      IF (IESTPI(IRPI,2).NE.0)
     .   WRITE (IUNOUT,*) 'COLLISION ESTIMATOR FOR MOM. BALANCE '
      IF (IESTPI(IRPI,3).NE.0)
     .   WRITE (IUNOUT,*) 'COLLISION ESTIMATOR FOR EN. BALANCE '
      CALL EIRENE_LEER(1)

      WRITE (IUNOUT,*) 'COLLISION MODEL: '
      WRITE (iunout,*) 'PROCESS NO. KK ',NREAPI(IRPI)
      WRITE (IUNOUT,*) 'MODCOL(1:4) ',
     .                  MODCOL(4,1,IRPI),MODCOL(4,2,IRPI),
     .                  MODCOL(4,3,IRPI),MODCOL(4,4,IRPI)
      WRITE (IUNOUT,'(1X,A15,1(1PE12.4))') 'SCALING FACTOR ',
     .                  FACRPI(IRPI,1)
      CALL EIRENE_LEER(1)


      RETURN
C
C
C-----------------------------------------------------------------------
C
  990 CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTPI: EXIT CALLED '
      WRITE (iunout,*) 'INVALID SPECIES INDEX FOR PI ',IRPI
      CALL EIRENE_EXIT_OWN(1)
  991 CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTPI: EXIT CALLED '
      WRITE (iunout,*) 'CHARGE CONSERVATION VIOLATED '
      WRITE (iunout,*) 'IRPI, TEST SPECIES, BULK SPECIES ',IRPI,
     .                  TEXTS(ISP),TEXTS(NSPAMI+IPL)
      CALL EIRENE_EXIT_OWN(1)
  992 CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTPI: EXIT CALLED '
      WRITE (iunout,*)
     .  'MASS NUMBERS OF INTERACTING PARTICLES INCONSISTENT'
      WRITE (iunout,*) 'KK ',KK
      CALL EIRENE_EXIT_OWN(1)
  993 CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTPI: EXIT CALLED '
      WRITE (iunout,*)
     .  'EBULK_ION .LE.0, BUT MONOENERGETIC DISTRIBUTION?'
      WRITE (iunout,*) 'CHECK ENERGY FLAG ISCDEA'
      WRITE (iunout,*) 'KK,ISCDEA ',KK,ISCDEA
      CALL EIRENE_EXIT_OWN(1)
  994 CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTPI: EXIT CALLED '
      WRITE (iunout,*)
     .  'SPECIES INDEX OF SECONDARY PARTICLE OUT OF RANGE'
      WRITE (iunout,*) 'KK ',KK
      CALL EIRENE_EXIT_OWN(1)
  995 CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTPI: EXIT CALLED '
      WRITE (iunout,*)
     .  'STORAGE SAVING MODE NOT READY; KK, IRPI'
      WRITE (iunout,*) 'KK, IRPI ',KK,IRPI
      CALL EIRENE_EXIT_OWN(1)
  996 CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTPI: INVALID DATA OPTION'
      WRITE (iunout,*) 'IRPI, MODC ',IRPI, MODC
      CALL EIRENE_EXIT_OWN(1)
  997 CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTPI: ISCDE FLAG'
      WRITE (iunout,*) IRPI
      CALL EIRENE_EXIT_OWN(1)
  999 CONTINUE
      WRITE (iunout,*) 'INSUFFICIENT STORAGE FOR PI: NRPI=',NRPI
      CALL EIRENE_EXIT_OWN(1)
      RETURN
C
      END
