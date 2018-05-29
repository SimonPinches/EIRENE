cdr modc=3 fuer cx rate coeff angefangen: um multi-step cx auch vs. t und n zu kriegen,
cdr aber dann die Frage:  te=Ti,  ne=ni ? Und E0 immer sehr klein? Korrektes te,ti,ne,ni 
cdr koennen zellweise kommen, z.b. aus CRM modell. 
cdr dann bleibt es bei einem 9-parameter fit (fuer E0 abhaengigkeit)  pro Zelle.


C 24.11.05: chrdf0 introduced (called from xsecta, xsectm, xsecti)
C 08.08.06: error exit 991 introduced: charge conservation violation
! 30.08.06: data structure for reaction data redefined
! 12.10.06: modcol revised
! 22.11.06: flag for shift of first parameter to rate_coeff introduced
! 08.01.07: pls = 0.dp, twice, preset.
! 01.02.07: do not evaluate rates in vacuum region for IPL (use lgvac(..IPL)
! 21.09.09: dr: a few more comments
! 20.01.14:  H.4 option for cx rate coefficients (e.g. CR rates: p + H-minus)
c            additional argument PLS, also in calling routines xsecta,xsectm,xsecti
C            remove plsti(nstordt), now: TII 
c 25.03.15:  rename nelrcx  to nplrcx, in order to enable 
c            consistency in notation with PI processes: not ready
cdr   sept.16: calls to prep_rtcs removed. prep_rtcs is now redundant
cdr   jan 17 : added nuclear charge number conservation test,
cdr            such that first secondary always corresponds to incident bulk particle,
cdr            just with charge state changed by an increment one (+1 or -1).
cdr   may 18 : this jan 17 fix was too narrow, e.g. He + He++ CX, charge state
cdr            increment=2.  Now fixed.
C

      SUBROUTINE EIRENE_XSTCX(RMASS,IRCX,ISP,IPL,
     .                        ISCD1,ISCD2,
     .                        EBULK, CHRDF0,ISCDE,IESTM,
     .                        KK,FACTKK,PLS)

c  set NON-DEFAULT cx collision cross-sections and rates  
c  IPL{n+} + ISP -->  IPL1{(n-m)+} + ISP2{m+}
c  defaults for CX type processes:  exchange of identity

c  carry out some consistency checks
c  first  secondary == previous bulk particle
c  second secondary == previous test particle

c   rmass: incident test particle mass (needed for check of mass conservation)
c   ircx:  counter for CX reaction in this run
c   isp:   incident test species index
c   ipl:   incident bulk ion species index (0 < ipl <= nplsi)
c
c   KK :   reaction number in modclf (input) array
c   FACTKK: scaling factor for this collision process (cross-section and rates)
c   PLS:   precomputed log of electron density

C  RETURNS:
C    MODCOL(3,...)
C    TABCX3(IRCX,NCELL,...)  1/s per incident test particle
C    EPLCX3(IRCX,NCELL,...) eV/s per incident test particle
C    DEFCX(IRCX)
C    EEFCX(IRCX)
C    IESTCX(IRCX,...)
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

      REAL(DP), INTENT(IN) :: RMASS, EBULK, FACTKK, CHRDF0
      REAL(DP), INTENT(IN) :: PLS(NSTORDR)
      INTEGER, INTENT(IN) :: IRCX, ISP, IPL, ISCD1, ISCD2, 
     .                       ISCDE, IESTM, KK
      REAL(DP) :: CF(9)
      REAL(DP) :: ADD, ADDL, RMTEST, RMBULK, FCTKKL, CHRDIF,
     .            ADDT, ADDTL, TMASS, PMASS, COU, 
     .            EIRENE_RATE_COEFF,
     .            EIRENE_ENERGY_RATE_COEFF, TB, TII,
     .            FP1(6),FP2(6)
      INTEGER :: ITYP1, ITYP2, ISPZ1, ISPZ2, KREAD,
     .           J, NEND, MODC, NSECX4, IPL2, IIO2, IPLTI,
     .           NCBULK, NCGBLK
      INTEGER, EXTERNAL :: EIRENE_IDEZ
      CHARACTER(8) :: TEXTS1, TEXTS2
      type(poly_data), pointer :: rp
      type(fit_forms), pointer :: rt

      SAVE

      NREACX(IRCX) = KK  ! needed for storage saving mode
C
C  SET NON-DEFAULT CHARGE EXCHANGE COLLISION PROCESS NO. IRCX
C
      IF (IPL.LE.0.OR.IPL.GT.NPLSI) GOTO 990
      IF (MASSP(KK).LE.0.OR.MASST(KK).LE.0) GOTO 992
      RMBULK=RMASSP(IPL)
      NCBULK=NCHARP(IPL)
      NCGBLK=NCHRGP(IPL)
      RMTEST=RMASS
      IPLTI=MPLSTI(IPL)
C
C  1ST SECONDARY INDEX, PREVIOUS BULK MASS
      N1STX(IRCX,1)=EIRENE_IDEZ(ISCD1,1,3)    !TYPE
      N1STX(IRCX,2)=EIRENE_IDEZ(ISCD1,3,3)    !SPECIES WITHIN TYPE CLASS
      N1STX(IRCX,3)=0
      IF (N1STX(IRCX,1).LT.4) N1STX(IRCX,3)=1 !DEFAULT: 1 "FIRST" TEST SECONDARY, IF ANY

      IF ((N1STX(IRCX,2) < 1) .OR. 
     .    (N1STX(IRCX,2) > MAXSPC(N1STX(IRCX,1)))) GOTO 994

C   CHECK MASS AND NUCLEAR CHARGE NUMBER CONSERVATION, FIRST SECONDARY, PREV. BULK
      IF (N1STX(IRCX,1).EQ.1) THEN
        IF (RMBULK.NE.RMASSA(N1STX(IRCX,2))) GOTO 992
        IF (NCBULK.NE.NCHARA(N1STX(IRCX,2))) GOTO 992 
cdr     IF (ABS(NCGBLK-0).ne.1) GOTO 992  ! allow also for double CX
      ELSEIF (N1STX(IRCX,1).EQ.2) THEN
        IF (RMBULK.NE.RMASSM(N1STX(IRCX,2))) GOTO 992
        IF (NCBULK.NE.NCHARM(N1STX(IRCX,2))) GOTO 992
cdr     IF (ABS(NCGBLK-0).ne.1) GOTO 992  ! allow also for double CX
      ELSEIF (N1STX(IRCX,1).EQ.3) THEN
        IF (RMBULK.NE.RMASSI(N1STX(IRCX,2))) GOTO 992
        IF (NCBULK.NE.NCHARI(N1STX(IRCX,2))) GOTO 992
cdr     IF (ABS(NCGBLK-NCHRGI(N1STX(IRCX,2))).ne.1) GOTO 992
      ELSEIF (N1STX(IRCX,1).EQ.4) THEN
        IF (RMBULK.NE.RMASSP(N1STX(IRCX,2))) GOTO 992
        IF (NCBULK.NE.NCHARP(N1STX(IRCX,2))) GOTO 992
cdr     IF (ABS(NCGBLK-NCHRGP(N1STX(IRCX,2))).ne.1) GOTO 992
      ENDIF
C
C  2ND SECONDARY INDEX, PREVIOUS TEST PARTICLE MASS
      N2NDX(IRCX,1)=EIRENE_IDEZ(ISCD2,1,3)    !TYPE
      N2NDX(IRCX,2)=EIRENE_IDEZ(ISCD2,3,3)    !SPECIES WITHIN TYPE CLASS
      N2NDX(IRCX,3)=N1STX(IRCX,3)             !CUMULATED NO. OF SECONDARIES
      IF (N2NDX(IRCX,1).LT.4) N2NDX(IRCX,3)=N2NDX(IRCX,3)+1 !DEFAULT: 1 "SECOND" TEST SECONDARY, IF ANY
C
      IF ((N2NDX(IRCX,2) < 1) .OR. 
     .    (N2NDX(IRCX,2) > MAXSPC(N2NDX(IRCX,1)))) GOTO 994
C   CHECK MASS CONSERVATION, SECOND SECONDARY
      IF (N2NDX(IRCX,1).EQ.1) THEN
        IF (RMTEST.NE.RMASSA(N2NDX(IRCX,2))) GOTO 992
      ELSEIF (N2NDX(IRCX,1).EQ.2) THEN
        IF (RMTEST.NE.RMASSM(N2NDX(IRCX,2))) GOTO 992
      ELSEIF (N2NDX(IRCX,1).EQ.3) THEN
        IF (RMTEST.NE.RMASSI(N2NDX(IRCX,2))) GOTO 992
      ELSEIF (N2NDX(IRCX,1).EQ.4) THEN
        IF (RMTEST.NE.RMASSP(N2NDX(IRCX,2))) GOTO 992
      ENDIF

C  CHECK CHARGE CONSERVATION
      CHRDIF=CHRDF0-NCHRGP(IPL)
      IF (N1STX(IRCX,1).EQ.3) THEN
        IIO2=N1STX(IRCX,2)
        CHRDIF=CHRDIF+NCHRGI(IIO2)
      ENDIF
      IF (N1STX(IRCX,1).EQ.4) THEN
        IPL2=N1STX(IRCX,2)
        CHRDIF=CHRDIF+NCHRGP(IPL2)
      ENDIF
      IF (N2NDX(IRCX,1).EQ.3) THEN
        IIO2=N2NDX(IRCX,2)
        CHRDIF=CHRDIF+NCHRGI(IIO2)
      ENDIF
      IF (N2NDX(IRCX,1).EQ.4) THEN
        IPL2=N2NDX(IRCX,2)
        CHRDIF=CHRDIF+NCHRGP(IPL2)
      ENDIF
      IF (CHRDIF.NE.0) GOTO 991
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
      ADDCX(IRCX,IPL) = ADDTL
C
C CROSS-SECTION (E-LAB) AVAILABLE ?
      IF (EIRENE_IDEZ(MODCLF(KK),2,5).EQ.1) THEN
        MODCOL(3,1,IRCX)=KK
C  TENTATIVLEY ASSUME: SIGMA * V_EFF MODEL FOR RATE COEFFICIENT
        MODCOL(3,2,IRCX)=3
      ENDIF

C..................................................................
C 2. RATE COEFFICIENT  (CM**3/S) * TARGET DENSITY (CM**-3)
C..................................................................

      MODC=EIRENE_IDEZ(MODCLF(KK),3,5)

C  2.B)
      IF (MODC.EQ.1) NEND=1   ! rate coeff for (FIXED e0, e.g. E=0, TI)
C  2.C)
      IF (MODC.EQ.2) NEND=NSTORDT ! rate coeff vs. (E, TI)
          
C  2.B) RATE COEFFICIENT(TI, FIXED E0, E.G. E0=0)
      IF (EIRENE_IDEZ(MODCLF(KK),3,5).EQ.1) THEN
C       NEND=1
        IF (NSTORDR >= NRAD) THEN
          DO 245 J=1,NSBOX
            IF (LGVAC(J,IPL)) CYCLE
              TII=TIINL(IPLTI,J)+ADDTL
              COU = EIRENE_RATE_COEFF(KK,J,TII,0._DP,.TRUE.,0)
              TABCX3(IRCX,J,1)=COU*DIIN(IPL,J)*FACTKK
245       CONTINUE          
        ELSE ! NOT SUFFICIENT STORAGE ON TABCX3
C  STORAGE SAVE MODE NOT READY FOR THIS OPTION ??
        ENDIF
        MODCOL(3,2,IRCX)=1 
 
      ELSEIF (EIRENE_IDEZ(MODCLF(KK),3,5).EQ.2) THEN
C  2.C) RATE COEFFICIENT(TI,EBEAM)
C       NEND=9
        IF (NSTORDR >= NRAD) THEN
          FCTKKL=LOG(FACTKK)
          rt => reacdat(kk)%rtc
          fp1(1:3) = rt%fp1l
          fp1(4:6) = rt%fp1r
          fp2(1:3) = rt%fp2b
          fp2(4:6) = rt%fp2t  
          DO J=1,NSBOX
            IF (LGVAC(J,IPL)) CYCLE
              TII=TIINL(IPLTI,J)+ADDTL
cdr  safety cut off at TI= 0.1 eV. (TVAC=0.02)
              tii = max(-2.3_dp,tii)
c  evaluate 2 parametric fit, 
c  collaps this to a one parameter fit CF for EB dependence, evaluated at TII. 
              rp => reacdat(KK)%rtc%poly
              call EIRENE_dbl_poly (rp%dblpol,tii,0._dp,cou,cf,
     .               rt%rc1min, rt%rc1max, fp1, rt%jfex1mn, rt%jfex1mx,
     .               rt%rc2min, rt%rc2max, fp2, rt%jfex2mn, rt%jfex2mx,
     .               trcamd) 
              TABCX3(IRCX,J,1:9) = CF(1:9)
              TABCX3(IRCX,J,1)=TABCX3(IRCX,J,1)+DIINL(IPL,J)+FCTKKL
          END DO
        ELSE ! NOT SUFFICIENT STORAGE ON TABCX3
C  STORAGE SAVE MODE NOT READY FOR THIS OPTION ??

        ENDIF
        MODCOL(3,2,IRCX)=2

      ELSEIF (EIRENE_IDEZ(MODCLF(KK),3,5).EQ.3) THEN
C  2.D) RATE COEFFICIENT(TI=TE, NE=NI ?, E0 FIXED, E.G. E0=0.)
C       IF (MODC.EQ.3) NEND=1  rate coeff vs. (N, T), NEND NOT NEEDED
        FCTKKL=LOG(FACTKK)
        IF (NSTORDR >= NRAD) THEN                
          DO J=1,NSBOX
            IF (LGVAC(J,IPL)) CYCLE
            COU = EIRENE_RATE_COEFF(KK,J,TEINL(J),PLS(J),.FALSE.,1)
            TB = COU + FCTKKL
            IF (IFTFLG(KK,2) < 100) TB = TB + DIINL(IPL,J)
            TB=MAX(-100._DP,TB)
            TABCX3(IRCX,J,1)=EXP(TB)
          END DO
C         JEREACX(IRCX) = 9
        ELSE  ! ??
C  WHAT DO WE DO IN CASE NSTORDR < NRAD  ?
          write (iunout,*) 'storage save mode not available yet for CX'
          write (iunout,*) 'in case modc=3  (n,T-dependence).'
          write (iunout,*) 'exit called '
          call eirene_exit_own(1) 
        ENDIF
        MODCOL(3,2,IRCX)=1 !  indicate: rate coefficient as fct. of local plasma conditions only
      ELSE
C  NO RATE COEFFICIENT. IS THERE A CROSS-SECTION AT LEAST?
        IF (MODCOL(3,2,IRCX).NE.3) GOTO 996
      ENDIF
 
      FACRCX(IRCX,1) = FACTKK
      FACRCX(IRCX,2) = LOG(FACTKK)

      DEFCX(IRCX)=LOG(CVELI2*PMASS)
      EEFCX(IRCX)=LOG(CVELI2*TMASS)
C
C  3. BULK PARTICLE MOMENTUM LOSS RATE
C
C  4. BULK PARTICLE ENERGY LOSS RATE
C  4.1. HEAVY BULK PARTICLE ENERGY LOSS RATE
C
C  SET ENERGY LOSS RATE OF IMPACTING ION
C
      NSECX4=EIRENE_IDEZ(ISCDE,4,5)
      IF (NSECX4.EQ.0) THEN
C  4.1A) ENERGY LOSS RATE OF IMP. BULK PARTICLE = CONST.*RATECOEFF.
C        SAMPLE COLLIDING ION FROM DRIFTING MONOENERGETIC ISOTROPIC DISTRIBUTION
c        WITH WEIGHTING/REJECTION
        IF (EBULK.LE.0.D0) THEN
          IF (NSTORDR >= NRAD) THEN
            DO J=1,NSBOX
              EPLCX3(IRCX,J,1)=1.5*TIIN(IPLTI,J)+EDRIFT(IPL,J)
            ENDDO
            NELRCX(IRCX) = -3
          ELSE
            NELRCX(IRCX) = -3
          END IF
        ELSEIF (EBULK.GT.0.D0) THEN ! EBULK GT.0
          WRITE (iunout,*) 'WARNING FROM SUBR. XSTCX: IRCX ', IRCX
          WRITE (iunout,*) 'MODIFIED TREATMENT OF CHARGE EXCHANGE '
          WRITE (iunout,*) 'SAMPLE FROM MAXWELLIAN WITH T = ',EBULK/1.5
          WRITE (iunout,*) 'RATHER THAN WITH T = TIIN '
          WRITE (iunout,*) 'NOT FULLY IMPLEMENTED (VELOCX) '    
          CALL EIRENE_LEER(1)
          IF (NSTORDR >= NRAD) THEN
            DO 251 J=1,NSBOX
              EPLCX3(IRCX,J,1)=EBULK+EDRIFT(IPL,J)
251         CONTINUE
            NELRCX(IRCX) = -2
          ELSE
            NELRCX(IRCX) = -2
            EPLCX3(IRCX,1,1)=EBULK
          END IF
C       ELSE
CDR   ERROR: EBULK < 0 IS NOT FORESEEN
        ENDIF
        MODCOL(3,4,IRCX)=3
      ELSEIF (NSECX4.EQ.1) THEN
C  4.1B) ENERGY LOSS RATE OF IMP. ION = (1.5*TI+EDRIFT)* RATECOEFF.
C       SAMPLE COLLIDING ION FROM DRIFTING MAXWELLIAN
        IF (EBULK.LE.0.D0) THEN
          IF (NSTORDR >= NRAD) THEN
            DO 252 J=1,NSBOX
              EPLCX3(IRCX,J,1)=1.5*TIIN(IPLTI,J)+EDRIFT(IPL,J)
252         CONTINUE
            NELRCX(IRCX) = -3
          ELSE
            NELRCX(IRCX) = -3
          END IF

        ELSEIF (EBULK.GT.0.D0) THEN ! EBULK GT.0  ! EBULK GT.0

          WRITE (iunout,*) 'WARNING FROM SUBR. XSTCX: IRCX ', IRCX
          WRITE (iunout,*) 'MODIFIED TREATMENT OF CHARGE EXCHANGE '
          WRITE (iunout,*) 'SAMPLE FROM MAXWELLIAN WITH T = ',EBULK/1.5
          WRITE (iunout,*) 'RATHER THAN WITH T = TIIN '
          WRITE (iunout,*) 'NOT FULLY IMPLEMENTED (VELOCX) '  
          CALL EIRENE_LEER(1)
          IF (NSTORDR >= NRAD) THEN
            DO 2511 J=1,NSBOX
              EPLCX3(IRCX,J,1)=EBULK+EDRIFT(IPL,J)
2511        CONTINUE
            NELRCX(IRCX) = -2
          ELSE
            NELRCX(IRCX) = -2
            EPLCX3(IRCX,1,1)=EBULK
          END IF
C       ELSE
CDR   ERROR: EBULK < 0 IS NOT FORESEEN
        ENDIF
        MODCOL(3,4,IRCX)=1
C     ELSEIF (NSECX4.EQ.2) THEN
C  use i-integral expressions. to be written
      ELSEIF (NSECX4.EQ.3) THEN
C  4.1C)  ENERGY LOSS RATE OF IMP. ION = EN.-WEIGHTED RATE
C       SAMPLE COLLIDING ION FROM DRIFTING MAXWELLIAN, WITH WEIGHTING/REJECTION
        KREAD=EBULK
        IF (KREAD.EQ.0) THEN
c  data for mean ion energy loss are not available
c  use collision estimator for energy balance
          IF (EIRENE_IDEZ(IESTM,3,3).NE.1) THEN
            WRITE (iunout,*)
     .        'COLLISION ESTIMATOR ENFORCED FOR ION ENERGY '
            WRITE (iunout,*) 'IN CX COLLISION IRCX= ',IRCX
            WRITE (iunout,*) 'BECAUSE NO ENERGY-WEIGHTED RATE AVAILABLE'
          ENDIF
          IESTCX(IRCX,3)=1
          MODCOL(3,4,IRCX)=2
        ELSE
C  ION ENERGY-AVERAGED RATE AVAILABLE AS REACTION NO. "KREAD"
        NELRCX(IRCX) = KREAD
        MODC=EIRENE_IDEZ(MODCLF(KREAD),5,5)
        IF (MODC.GE.1.AND.MODC.LE.2) THEN
          MODCOL(3,4,IRCX)=MODC
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
                EPLCX3(IRCX,J,1)=EIRENE_ENERGY_RATE_COEFF
     .                          (KREAD,J,TII,
     .                           0._DP,.FALSE.,0)*DIIN(IPL,J)*ADD
254           CONTINUE
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


                EPLCX3(IRCX,J,1:9) = CF(1:9)
                EPLCX3(IRCX,J,1) = EPLCX3(IRCX,J,1)+DIINL(IPL,J)+ADDL
257           CONTINUE
            ENDIF

          ELSE  ! STORAGE SAVING MODE, no pre-defined tallies eplcx3
            IF (MODC.EQ.1) THEN
              ADD=FACTKK/ADDT
              EPLCX3(IRCX,1,1)=ADD   !  ????
            ELSEIF (MODC.EQ.2) THEN
              ADDL=LOG(FACTKK)-ADDTL
              FACRCX(IRCX,1) = EXP(ADDL)
              FACRCX(IRCX,2) = ADDL
            END IF
          END IF
        ENDIF
        ENDIF
      ELSE
        WRITE (iunout,*) 'NSECX4 ILL-DEFINED IN XSTCX '
        WRITE (iunout,*) 'check parameter ISCDE for process ircx ',ircx
        CALL EIRENE_EXIT_OWN(1)
      ENDIF

C  ESTIMATOR FOR CONTRIBUTION TO COLLISION RATES FROM THIS REACTION
      IESTCX(IRCX,1)=EIRENE_IDEZ(IESTM,1,3)
      IESTCX(IRCX,2)=EIRENE_IDEZ(IESTM,2,3)
      IF (IESTCX(IRCX,3).EQ.0) IESTCX(IRCX,3)=EIRENE_IDEZ(IESTM,3,3)
C
      ITYP1=N1STX(IRCX,1)
      ITYP2=N2NDX(IRCX,1)

      IF (IESTCX(IRCX,1).NE.0.AND.(ITYP1.NE.1.OR.ITYP2.NE.4)) THEN
        CALL EIRENE_LEER(1)
        WRITE (iunout,*)
     .    'WARNING: COLL.EST NOT AVAILABLE FOR PART. BALANCE '
        WRITE (iunout,*) 'IRCX = ',IRCX
        WRITE (iunout,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
        IESTCX(IRCX,1)=0
      ENDIF
      IF (IESTCX(IRCX,2).NE.0.AND.(ITYP1.NE.1.OR.ITYP2.NE.4)) THEN
        CALL EIRENE_LEER(1)
        WRITE (iunout,*)
     .    'WARNING: COLL.EST NOT AVAILABLE FOR MOM. BALANCE '
        WRITE (iunout,*) 'IRCX = ',IRCX
        WRITE (iunout,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
        IESTCX(IRCX,2)=0
      ENDIF

      IF (IESTCX(IRCX,3).NE.0.AND.(ITYP1.NE.1.OR.ITYP2.NE.4)) THEN
        CALL EIRENE_LEER(1)
        WRITE (iunout,*)
     .    'WARNING: COLL.EST NOT AVAILABLE FOR EN. BALANCE '
        WRITE (iunout,*) 'IRCX = ',IRCX
        WRITE (iunout,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
        IESTCX(IRCX,3)=0
      ENDIF
      RETURN
C
C-----------------------------------------------------------------------
C

      ENTRY EIRENE_XSTCX_2(IRCX,IPL)
C
      CALL EIRENE_LEER(2)
      WRITE (iunout,*) 'CHARGE EXCHANGE REACTION NO. IRCX= ',IRCX
      CALL EIRENE_LEER(1)
      WRITE (iunout,*) 'CHARGE EXCHANGE WITH BULK IONS IPLS:'
      WRITE (iunout,*) 'IPLS= ',TEXTS(NSPAMI+IPL)
      CALL EIRENE_LEER(1)

      WRITE (iunout,*) '1ST AND 2ND NEXT GEN. SPECIES I2ND1, I2ND2:'

      ITYP1=N1STX(IRCX,1)
      ITYP2=N2NDX(IRCX,1)
      ISPZ1=N1STX(IRCX,2)
      ISPZ2=N2NDX(IRCX,2)
      IF (ITYP1.EQ.1) TEXTS1=TEXTS(NSPH+ISPZ1)
      IF (ITYP1.EQ.2) TEXTS1=TEXTS(NSPA+ISPZ1)
      IF (ITYP1.EQ.3) TEXTS1=TEXTS(NSPAM+ISPZ1)
      IF (ITYP1.EQ.4) TEXTS1=TEXTS(NSPAMI+ISPZ1)
      IF (ITYP2.EQ.1) TEXTS2=TEXTS(NSPH+ISPZ2)
      IF (ITYP2.EQ.2) TEXTS2=TEXTS(NSPA+ISPZ2)
      IF (ITYP2.EQ.3) TEXTS2=TEXTS(NSPAM+ISPZ2)
      IF (ITYP2.EQ.4) TEXTS2=TEXTS(NSPAMI+ISPZ2)
      WRITE (iunout,*) 'I2ND1= ',TEXTS1, 'I2ND2= ',TEXTS2

      CALL EIRENE_LEER(1)

      IF (IESTCX(IRCX,1).NE.0)
     .   WRITE (IUNOUT,*) 'COLLISION ESTIMATOR FOR PART. BALANCE '
      IF (IESTCX(IRCX,2).NE.0)
     .   WRITE (IUNOUT,*) 'COLLISION ESTIMATOR FOR MOM. BALANCE '
      IF (IESTCX(IRCX,3).NE.0)
     .   WRITE (IUNOUT,*) 'COLLISION ESTIMATOR FOR EN. BALANCE '
      CALL EIRENE_LEER(1)

      WRITE (IUNOUT,*) 'COLLISION MODEL: '
      WRITE (IUNOUT,*) 'PROCESS NO. KK ',NREACX(IRCX)
      WRITE (IUNOUT,*) 'MODCOL ',MODCOL(3,1,IRCX),MODCOL(3,2,IRCX),
     .                           MODCOL(3,3,IRCX),MODCOL(3,4,IRCX)
      WRITE (IUNOUT,'(1X,A15,1(1PE12.4))') 'SCALING FACTOR ',
     .                  FACRCX(IRCX,1) 
      CALL EIRENE_LEER(1)


      RETURN
C
C
C-----------------------------------------------------------------------
C
990   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTCX: EXIT CALLED '
      WRITE (iunout,*) 'INVALID SPECIES INDEX FOR CX ',IRCX
      CALL EIRENE_EXIT_OWN(1)
991   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTCX: EXIT CALLED '
      WRITE (iunout,*) 'CHARGE CONSERVATION VIOLATED '
      WRITE (iunout,*) 'IRCX, TEST-SPECIES, BULK SPECIES ',IRCX,
     .                  TEXTS(ISP),TEXTS(NSPAMI+IPL)
      CALL EIRENE_EXIT_OWN(1)
992   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTCX: EXIT CALLED '
      WRITE (iunout,*)
     .  'INTERACTING PARTICLES INCONSISTENT (MASS OR CHARGE)'
      WRITE (iunout,*) 'KK ',KK
      WRITE (iunout,*) 'IRCX, TEST-SPECIES, BULK SPECIES ',IRCX,
     .                  TEXTS(ISP),TEXTS(NSPAMI+IPL)
      CALL EIRENE_EXIT_OWN(1)
993   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTCX: EXIT CALLED '
      WRITE (iunout,*)
     .  'EBULK_ION .LE.0, BUT MONOENERGETIC DISTRIBUTION?'
      WRITE (iunout,*) 'CHECK ENERGY FLAG ISCDEA'
      WRITE (iunout,*) 'KK,ISCDEA ',KK,ISCDEA
      CALL EIRENE_EXIT_OWN(1)
994   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTCX: EXIT CALLED '
      WRITE (iunout,*)
     .  'SPECIES INDEX OF SECONDARY PARTICLE OUT OF RANGE'
      WRITE (iunout,*) 'KK ',KK
      CALL EIRENE_EXIT_OWN(1)
996   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTCX: EXIT CALLED '
      WRITE (iunout,*) 'NO CROSS-SECTION AVAILABLE FOR NON-DEFAULT CX'
      WRITE (iunout,*) 'KK ',KK
      WRITE (iunout,*)
     .  'EITHER PROVIDE CROSS-SECTION OR USE DIFFERENT '
      WRITE (iunout,*) 'POST-COLLISION SAMPLING FLAG ISCDEA'
      CALL EIRENE_EXIT_OWN(1)
      END
