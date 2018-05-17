!pb  30.08.06: data structure for reaction data redefined
!pb  12.10.06: modcol revised
!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced
cdr  05.01.07:  write(6,...) --> write(iunout,...) in one place
!    01.02.07: do not evaluate rates in vacuum region for IPL (use lgvac(..IPL)

cdr  20.04.14: bug fix: + edrift(...) was missing in eplel3, in case nseel4=0 and ebulk>0
cdr    oct.14: bug fix: use kread rather than kk in eplel3.
cdr    oct.14: remove pls array, synconize with xstcx started
cdr    aug.16: nend is always =1 or =9, remove redundant arguments in prep_poly
cdr   sept.16: calls to prep_rtcs removed. prep_rtcs is now redundant
cdr   jan .17: modcol(5,0,irel):  flag for differential cross section model, rather than =kk.
!              modcol(5,0,irel)=-1  : bgk (relaxation) collision, scattering angle =Pi in COM
!              modcol(5,0,irel)=0   : isotropic in COM, assume: the cross section
!                                     and rate coefficients are "diffusion" cross section,
!                                     and rate coefficients, respectively.
!              modcol(5,0,irel)=1,2,...: interaction potential is given via fit parameters
cdr     currently still: modcol(5,0,irel)=kk, and veloel uses reacdat(kk) directly.

cdr     Reaction identifyer KK is defined twice, within same routine veloel.
cdr     This risky exception can be removed by: modcol(5,0,irel)=iftflg(kk,0),
cdr     and by providing the potential p(1:9,irel) here, rather than in veloel.  
cdr  nov. 17:  added: parameter pls (as in xstcx,xstpi,...)                 
C
C
      SUBROUTINE EIRENE_XSTEL(IREL,ISP,IPL,
     .                 EBULK,ISCDE,IESTM,
     .                 KK,FACTKK,PLS)
C
C       SET UP TABLES (E.G. OF REACTION RATE ) FOR EL PROCESSES
C
C   MEANING OF INPUT VARIABLES: SEE XSTCX

C   KK:      COMMON IDENTIFIER FOR PROCESS, USED FOR POTENTIAL, CROSS SECTION, RATES,
C                                           STORAGE SAVING MODE ETC...
C   FACTKK:  COMMON SCALING FACTOR FOR PROCESS KK
C   NREAEL(IREL) = KK DURING MC RUN. THIS ESTABLISHES LINK BETWEEN IREL AND KK, MUST BE UNIQUE


C  RETURNS:
C    MODCOL(5,...)
C    TABEL3(IREL,NCELL,...)  1/s per incident test particle
C    EPLEL3(IREL,NCELL,...) eV/s per incident test particle
C    DEFEL(IREL)
C    EEFEL(IREL)
C    IESTEL(IREL,...)
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

      REAL(DP), INTENT(IN) :: EBULK, FACTKK
      REAL(DP), INTENT(IN) :: PLS(NSTORDR)
      INTEGER, INTENT(IN) :: IREL, ISP, IPL,
     .                       ISCDE, IESTM, KK
      REAL(DP) :: CF(9)
      REAL(DP) :: ADD, ADDL, ADDT, FCTKKL, ADDTL, PMASS, TMASS, COU,
     .            EIRENE_RATE_COEFF,
     .            EIRENE_ENERGY_RATE_COEFF, 
     .            TB, TII,
     .            FP1(6),FP2(6)
      INTEGER :: NSEEL4, NEND, J, KREAD, MODC,  IPLTI,
     .           IBGK,ISPZB,ITYPB
      INTEGER, EXTERNAL :: EIRENE_IDEZ
      type(poly_data), pointer :: rp
      type(fit_forms), pointer :: rt

      SAVE

      NREAEL(IREL) = KK  ! needed for storage saving mode
C
C  SET NON-DEFAULT ELASTIC COLLISION PROCESS NO. IREL
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
      ADDEL(IREL,IPL) = ADDTL
      
      IPLTI = MPLSTI(IPL)

C..................................................................
C 0. INTERACTION POTENTIAL, DIFFERENTIAL CROSS SECTION INFORMATION, ETC....
C..................................................................
      IF (EIRENE_IDEZ(MODCLF(KK),1,5).EQ.1) THEN
cdr  use total cross section and rate coefficients for transport.
cdr  differential cross section or interaction potential for collision kinetics
        MODCOL(5,0,IREL)=KK  !  fit parameters for interaction potential
cdr                          !  this should become = iftflg(kk,0),
cdr
cdr                          !  set here: pot(1:9,irel)=reacdat(kk):.....
      ELSEIF (EIRENE_IDEZ(MODCLF(KK),1,5).EQ.0) THEN      
cdr  use diffusion cross section and diffusion rate coeff. for transport
        modcol(5,0,irel)=0   !  isotropic scattering IN COM
        if (NPBGKP(IPL,1).eq.0) then
          WRITE (IUNOUT,*) 'WARNING FROM XSTEL: '
          WRITE (IUNOUT,*) 'KK, IREL ',KK,IREL
          WRITE (IUNOUT,*) 'NO SCATTERING ANGLE INFORMATION PROVIDED'
          WRITE (IUNOUT,*) 'BUT ALSO NO BGK RELAXATION.'
          WRITE (IUNOUT,*) 'USE ISOTROPIC SCATTERING'    
        endif

cdr  or
cdr  use 0.5*(diffusion cross section) and 0.5*(diffusion rate coeff.) for transport 
c       modcol(5,0,irel) =-1, scattering angle =PI IN COM (=exchange of identity in LAB)
      ENDIF
C
C...................................................................
C 1. CROSS SECTION (E-LAB) (CM**2) , AVAILABLE ?
C...................................................................

      IF (EIRENE_IDEZ(MODCLF(KK),2,5).EQ.1) THEN
        MODCOL(5,1,IREL)=KK
C  TENTATIVLEY ASSUME: SIGMA * V_EFF MODEL FOR RATE COEFFICIENT
        MODCOL(5,2,IREL)=3
      ENDIF

C..................................................................
C 2. RATE COEFFICIENT  (CM**3/S) * TARGET DENSITY (CM**-3)
C..................................................................

      MODC=EIRENE_IDEZ(MODCLF(KK),3,5)

      IF (MODC.GE.1.AND.MODC.LE.2) THEN

        MODCOL(5,2,IREL)=MODC
C  2.B)
        IF (MODC.EQ.1) NEND=1       ! rate coeff for (FIXED E0, e.g. E0=0.0, TI)
C  2.C)
        IF (MODC.EQ.2) NEND=NSTORDT ! rate coeff vs. (E0, TI) NEND=9 HERE
C   STORAGE SAVING MODE ?
        IF (NSTORDR >= NRAD) THEN
C   NO, NSTORDT=9 HERE
          
C  2.B) RATE COEFFICIENT(TI, FIXED E0, E.G. E0=0)
          IF (MODC.EQ.1) THEN
C           NEND=1
            DO 245 J=1,NSBOX
              IF (LGVAC(J,IPL)) CYCLE
              TII=TIINL(IPLTI,J)+ADDTL
              COU = EIRENE_RATE_COEFF(KK,J,TII,0._DP,.TRUE.,0)
              TABEL3(IREL,J,1)=COU*DIIN(IPL,J)*FACTKK
245         CONTINUE
          ELSEIF (MODC.EQ.2) THEN
C           NEND=9
C  2.C) RATE COEFFICIENT(TI,EBEAM)
C       NEND=9
          FCTKKL=LOG(FACTKK)
          rt => reacdat(kk)%rtc
          fp1(1:3) = rt%fp1l
          fp1(4:6) = rt%fp1r
          fp2(1:3) = rt%fp2b
          fp2(4:6) = rt%fp2t
          DO J=1,NSBOX
            IF (LGVAC(J,IPL)) CYCLE
              TII=TIINL(IPLTI,J)+ADDTL
! this is another cut off, at TIIN <=0.1 eV rather than at TVAC = 0.02 ev
              tii = max(-2.3_dp,tii) 
c old
c old         CALL EIRENE_PREP_RTCS (KK,3,TII,CF)
c old
              rp => reacdat(KK)%rtc%poly
              call EIRENE_dbl_poly (rp%dblpol,tii,0._dp,cou,cf,
     .               rt%rc1min, rt%rc1max, fp1, rt%jfex1mn, rt%jfex1mx,
     .               rt%rc2min, rt%rc2max, fp2, rt%jfex2mn, rt%jfex2mx,
     .               trcamd)

              TABEL3(IREL,J,1:9) = CF(1:9)
              TABEL3(IREL,J,1)=TABEL3(IREL,J,1)+DIINL(IPL,J)+FCTKKL
            END DO
          END IF  ! MODC=1,2
        ELSE ! NOT SUFFICIENT STORAGE ON TABEL3
C  STORAGE SAVE MODE NOT READY FOR THIS OPTION ??
          GOTO 995

        ENDIF
      ELSEIF (EIRENE_IDEZ(MODCLF(KK),3,5).EQ.3) THEN
C  2.D) RATE COEFFICIENT(TI=TE, NE=NI ?, E0 FIXED, E.G. E0=0.)
C       IF (MODC.EQ.3) NEND=1  rate coeff vs. (N, T), NEND NOT NEEDED

        MODCOL(5,2,IREL)=1 !  indicate: rate coefficient as fct. of local plasma conditions only
        FCTKKL=LOG(FACTKK)
        IF (NSTORDR >= NRAD) THEN 
                
          DO J=1,NSBOX
            IF (LGVAC(J,IPL)) CYCLE
            COU = EIRENE_RATE_COEFF(KK,J,TEINL(J),PLS(J),.FALSE.,1)
            TB = COU + FCTKKL
            IF (IFTFLG(KK,2) < 100) TB = TB + DIINL(IPL,J)
            TB=MAX(-100._DP,TB)
            TABEL3(IREL,J,1)=EXP(TB)
          END DO
C         JEREAEL(IREL) = 9
        ELSE  ! ??

C  WHAT DO WE DO IN CASE NSTORDR < NRAD  ?
          write (iunout,*) 'storage save mode not available yet for EL'
          write (iunout,*) 'in case modc=3  (n,T-dependence).'
          GOTO 995 
        ENDIF

      ELSE
C  NO RATE COEFFICIENT. IS THERE A CROSS-SECTION AT LEAST?
        IF (MODCOL(5,2,IREL).NE.3) GOTO 993
      ENDIF

      FACREL(IREL,1) = FACTKK
      FACREL(IREL,2) = LOG(FACTKK)
 
      DEFEL(IREL)=LOG(CVELI2*PMASS)
      EEFEL(IREL)=LOG(CVELI2*TMASS)
C
C  3. BULK PARTICLE MOMENTUM LOSS RATE
C
C  4. BULK PARTICLE ENERGY LOSS RATE
C  4.1. HEAVY BULK PARTICLE ENERGY LOSS RATE
C
C  SET ENERGY LOSS RATE OF IMPACTING ION
C
      NSEEL4=EIRENE_IDEZ(ISCDE,4,5)
      IF (NSEEL4.EQ.0) THEN
C  4.1A)  ENERGY LOSS RATE OF IMP. BULK PARTICLE = CONST.*RATECOEFF.
C        SAMPLE COLLIDING ION FROM DRIFTING MONOENERGETIC ISOTROPIC DISTRIBUTION
c        WITH WEIGHTING/REJECTION
        IF (EBULK.LE.0.D0) THEN
          IF (NSTORDR >= NRAD) THEN
            DO J=1,NSBOX
              EPLEL3(IREL,J,1)=1.5*TIIN(IPLTI,J)+EDRIFT(IPL,J)
            ENDDO
            NELREL(IREL) = -3
          ELSE
            NELREL(IREL) = -3
          END IF
        ELSE ! EBULK GT.0
          IF (NSTORDR >= NRAD) THEN
            DO 251 J=1,NSBOX
              EPLEL3(IREL,J,1)=EBULK+EDRIFT(IPL,J)
251         CONTINUE
            NELREL(IREL) = -1
          ELSE
            NELREL(IREL) = -1
            EPLEL3(IREL,1,1)=EBULK
          END IF
C       ELSE
CDR   ERROR: EBULK < 0 IS NOT FORESEEN
        ENDIF
        MODCOL(5,4,IREL)=3
      ELSEIF (NSEEL4.EQ.1) THEN
C  4.1B) ENERGY LOSS RATE OF IMP. ION = (1.5*TI+EDRIFT)* RATECOEFF.
C       SAMPLE COLLIDING ION FROM DRIFTING MAXWELLIAN
        IF (EBULK.LE.0.D0) THEN
          IF (NSTORDR >= NRAD) THEN
            DO 252 J=1,NSBOX
              EPLEL3(IREL,J,1)=1.5*TIIN(IPLTI,J)+EDRIFT(IPL,J)
252         CONTINUE
            NELREL(IREL) = -3
          ELSE
            NELREL(IREL) = -3
          END IF
        ELSE ! EBULK GT.0
          WRITE (iunout,*) 'WARNING FROM SUBR. XSTEL: IREL ', IREL
          WRITE (iunout,*) 'MODIFIED TREATMENT OF ELASTIC COLLISIONS '
          WRITE (iunout,*) 'SAMPLE FROM MAXWELLIAN WITH T = ',EBULK/1.5
          WRITE (iunout,*) 'RATHER THAN WITH T = TIIN '
          WRITE (iunout,*) 'NOT FULLY IMPLEMENTED (VELOEL) '  
          CALL EIRENE_LEER(1)
          IF (NSTORDR >= NRAD) THEN
            DO 2511 J=1,NSBOX
              EPLEL3(IREL,J,1)=EBULK+EDRIFT(IPL,J)
2511        CONTINUE
            NELREL(IREL) = -2
          ELSE
            NELREL(IREL) = -2
            EPLEL3(IREL,1,1)=EBULK
          END IF
C       ELSE
CDR   ERROR: EBULK < 0 IS NOT FORESEEN
        ENDIF
        MODCOL(5,4,IREL)=1
C     ELSEIF (NSEEL4.EQ.2) THEN
C  use i-integral expressions. to be written
      ELSEIF (NSEEL4.EQ.3) THEN
C  4.1C)  ENERGY LOSS RATE OF IMP. ION = EN.-WEIGHTED RATE
C       SAMPLE COLLIDING ION FROM DRIFTING MAXWELLIAN, WITH WEIGHTING/REJECTION
        KREAD=EBULK
        IF (KREAD.EQ.0) THEN
c  data for mean ion energy loss are not available
c  use collision estimator for energy balance
          IF (EIRENE_IDEZ(IESTM,3,3).NE.1) THEN
            WRITE (iunout,*)
     .        'COLLISION ESTIMATOR ENFORCED FOR ION ENERGY '
            WRITE (iunout,*) 'IN ELASTIC COLLISION IREL= ',IREL
            WRITE (iunout,*) 'BECAUSE NO ENERGY-WEIGHTED RATE AVAILABLE'
          ENDIF
          IESTEL(IREL,3)=1
          MODCOL(5,4,IREL)=2
        ELSE
C  ION ENERGY-AVERAGED RATE AVAILABLE AS REACTION NO. "KREAD"
        NELREL(IREL) = KREAD
        MODC=EIRENE_IDEZ(MODCLF(KREAD),5,5)
        IF (MODC.GE.1.AND.MODC.LE.2) THEN
          MODCOL(5,4,IREL)=MODC
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
                EPLEL3(IREL,J,1)=EIRENE_ENERGY_RATE_COEFF
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

                EPLEL3(IREL,J,1:9) = CF(1:9)
                EPLEL3(IREL,J,1) = EPLEL3(IREL,J,1)+DIINL(IPL,J)+ADDL
257           CONTINUE
            ENDIF

          ELSE  ! STORAGE SAVING MODE, no pre-defined tallies eplel3
            IF (MODC.EQ.1) THEN
              ADD=FACTKK/ADDT
              EPLEL3(IREL,1,1)=ADD   !  ????
            ELSEIF (MODC.EQ.2) THEN
              ADDL=LOG(FACTKK)-ADDTL
              FACREL(IREL,1) = EXP(ADDL)
              FACREL(IREL,2) = ADDL
            END IF
          END IF
        ENDIF
        ENDIF
      ELSE
        WRITE (iunout,*) 'NSEEL4 ILL-DEFINED IN XSTEL '
        WRITE (iunout,*) 'check parameter ISCDE for process irel ',irel
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
C
C  4.2. BULK ELECTRON ENERGY LOSS RATE  ! NOT APPLICABLE
C
C
C  4.3. HEAVY PARTICLE ENERGY GAIN RATE
C
C
C  ESTIMATOR FOR CONTRIBUTION TO COLLISION RATES FROM THIS REACTION
      IESTEL(IREL,1)=EIRENE_IDEZ(IESTM,1,3)
      IESTEL(IREL,2)=EIRENE_IDEZ(IESTM,2,3)
      IF (IESTEL(IREL,3).EQ.0) IESTEL(IREL,3)=EIRENE_IDEZ(IESTM,3,3)
C

      IF (IESTEL(IREL,2).EQ.0.AND.NPBGKP(IPL,1).EQ.0) THEN
        CALL EIRENE_LEER(1)
        WRITE (iunout,*)
     .    'WARNING: TR.L.EST NOT AVAILABLE FOR MOM. BALANCE '
        WRITE (iunout,*) 'IREL = ',IREL
        WRITE (iunout,*) 'AUTOMATICALLY RESET TO COLLISION ESTIMATOR '
        IESTEL(IREL,2)=1
      ENDIF
      IF (IESTEL(IREL,3).EQ.0.AND.NPBGKP(IPL,1).EQ.0) THEN
        CALL EIRENE_LEER(1)
        WRITE (iunout,*)
     .    'WARNING: TR.L.EST NOT AVAILABLE FOR EN. BALANCE '
        WRITE (iunout,*) 'IREL = ',IREL
        WRITE (iunout,*) 'AUTOMATICALLY RESET TO COLLISION ESTIMATOR '
        IESTEL(IREL,3)=1
      ENDIF
      RETURN
C
C-----------------------------------------------------------------------
C

      ENTRY EIRENE_XSTEL_2(IREL,IPL)
C
      CALL EIRENE_LEER(2)
      WRITE (iunout,*) 'ELASTIC COLLISION NO. IREL= ',IREL
      CALL EIRENE_LEER(1)
      WRITE (iunout,*) 'ELASTIC COLLISION WITH BULK IONS IPLS:'
      WRITE (iunout,*) 'IPLS= ',TEXTS(NSPAMI+IPL)
      CALL EIRENE_LEER(1)
C
      IF (NPBGKP(IPL,1).NE.0) THEN
        IBGK=NPBGKP(IPL,1)
        WRITE (iunout,*) 'THIS IS ALSO BGK COLLISION NO. IBGK= ',IBGK
        MODCOL(5,0,IREL)=-1
        IF (NPBGKP(IPL,2).EQ.0)
     .      WRITE (iunout,*) 'SELF COLLISION      ' 
        IF (NPBGKP(IPL,2).NE.0) THEN
          ITYPB=EIRENE_IDEZ(NPBGKP(IPL,2),1,3)
          ISPZB=EIRENE_IDEZ(NPBGKP(IPL,2),3,3)
          IF (ITYPB.EQ.1)
     .      WRITE (iunout,*) 'CROSS COLLISION WITH ATOM     ',ISPZB
          IF (ITYPB.EQ.2)
     .      WRITE (iunout,*) 'CROSS COLLISION WITH MOLECULE ',ISPZB
          IF (ITYPB.EQ.3)
     .      WRITE (iunout,*) 'CROSS COLLISION WITH TEST ION ',ISPZB
        ENDIF
      ENDIF

      CALL EIRENE_LEER(1)

      IF (IESTEL(IREL,1).NE.0)
     .   WRITE (IUNOUT,*) 'COLLISION ESTIMATOR FOR PART. BALANCE '
      IF (IESTEL(IREL,2).NE.0)
     .   WRITE (IUNOUT,*) 'COLLISION ESTIMATOR FOR MOM. BALANCE '
      IF (IESTEL(IREL,3).NE.0)
     .   WRITE (IUNOUT,*) 'COLLISION ESTIMATOR FOR EN. BALANCE '
      CALL EIRENE_LEER(1)

      WRITE (IUNOUT,*) 'COLLISION MODEL: '
      WRITE (iunout,*) 'PROCESS NO. KK ',NREAEL(IREL)
      WRITE (IUNOUT,*) 'MODCOL(0)   ',MODCOL(5,0,IREL)
      WRITE (IUNOUT,*) 'MODCOL(1:4) ',
     .                  MODCOL(5,1,IREL),MODCOL(5,2,IREL),
     .                  MODCOL(5,3,IREL),MODCOL(5,4,IREL)
      WRITE (IUNOUT,'(1X,A15,1(1PE12.4))') 'SCALING FACTOR ',
     .                  FACREL(IREL,1)
      CALL EIRENE_LEER(1)


      RETURN
C
993   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTEL, SPECIES ISP: '
      WRITE (iunout,*) ISP,IREL
      CALL EIRENE_EXIT_OWN(1)
995   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTEL: EXIT CALLED '
      WRITE (iunout,*)
     .  'STORAGE SAVING MODE NOT READY; KK, IREL'
      WRITE (iunout,*) 'KK, IREL ',KK,IREL
      CALL EIRENE_EXIT_OWN(1)
      END
