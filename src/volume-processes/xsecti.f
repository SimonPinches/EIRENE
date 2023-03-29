C  27.6.05  irds --> irei
c 24.11.05 use nprt(ispz) to check if iion is a molecular ion
c          otherwise he+ atomic ions could be confused with d2+ molecular
c          ions and then get assigned the wrong default collision model
c 24.11.05 chrdf0 in parameterlist for call to xstcx
c          (was ok already for call to xstei)
! 30.08.06: data structure for reaction data redefined
! 12.10.06: modcol revised
! 22.11.06: flag for shift of first parameter to rate_coeff introduced
!           setting of modcol corrected
! 02.03.07: remove escd2* arrays
! 22.03.07: PI reactions revised
! 25.03.07: 3rd and 4th secondary introduced
! 2013    : DSUB (RESCALING OF DENSITY IN H.4 FITS) REMOVED, NOW DONE IN RATE_COEFF.F
! 2013    : DENSITY LIMIT 1E8 SET FOR POLYNOMIAL FITS (ARRAY PLS).
! 23.02.14: call to xstcx: additional arguments: pls  (for H.4 option)
! 23.02.14: call to xstpi: additional arguments: III, pls (for H.4 option)
cdr  oct.14:  pls made allocatable
cdr  oct.14:  eelei1 set in storage save mode, for default models (was missing)
cdr  oct.14:  further synchronization with xsectm,xsecta,
cdr           remaining relevant differences in default models only.
cdr  aug.15:  ibgk_sp:  no of bgk species. to be distinguished from ibgk: no of bgk reaction.
cdr  apr.16:  accmas and accinv set explicitly also for reaction -9
cdr           (was missing, but accidentally correct)

!pb  APR  16:  pplds  -> pplei
!pb  APR  16:  patds  -> patei,  eatds -> eatei
!pb  APR  16:  pelds  -> pelei,  eelds -> eelei
!pb  MAY  16:  tabds1 -> tabds1
!pb  JUL  16:  ehvds1 -> ehvds1
cdr  sept 16:  nidsi  -> nieii
cdr  aug.18:   process kk=-9, default KER changed from 0.5 to 0.8  (=0.4 per particle)
cdr            to better sync with HYDHEL original data.
cdr  sept 18:  rationalization for nhvrei flags for default reaction:  same as -KK
C
      SUBROUTINE EIRENE_XSECTI
C
C       SET UP TABLES (E.G. OF REACTION RATES) FOR TEST IONS
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CCONA
      USE EIRMOD_CGRID
      USE EIRMOD_CZT1
      USE EIRMOD_CTRCEI
      USE EIRMOD_CTEXT
      USE EIRMOD_COMXS
      USE EIRMOD_CSPEI

      IMPLICIT NONE

      REAL(DP), ALLOCATABLE :: PLS(:)
      REAL(DP) :: FACTKK, CHRDF0, RMASS, DEIMIN,
     .          EHEAVY, EELEC, EBULK, COU, EIRENE_RATE_COEFF,
     .          ACCMAS, ACCINV,DE_10

      INTEGER :: ICOUNT, IA1, IP2, IPLS, ITEST, IIO, IION, IDSC1,
     .           NRC, J, IPLS1, IPLS2, IATM, KK, IATM1, IATM2, ITYPB,
     .           ISPZB, III, IDSC, IREL, IBGK_SP,
     .           IIEL, IIEI, IREI, IESTM, IFRST, ISCND, ISCDE, IPL,
     .           IICX, IRCX, ITHRD, IFRTH, IRPI, IIPI
      INTEGER, EXTERNAL :: EIRENE_IDEZ
      ALLOCATE (PLS(NSTORDR))

cdr  PLS: ELECTRON DENSITY PARAMETER in CR MODELS
cdr      (NOT TO BE CONFUSED WITH THE DENSITY FACTOR BETWEEN RATES AND RATE COEFF.)
cdr: set hard-wired lower density for H.4, H.10 type fits from AMJUEL: 1e8 cm**-3
cdr: at this lower limit density the fits are produced such
cdr: that they collapse to the corona limit values.
      DEIMIN=LOG(1.D8)
      IF (NSTORDR >= NRAD) THEN
        DO 10 J=1,NSBOX
          PLS(J)=MAX(DEIMIN,DEINL(J))
   10   CONTINUE
      END IF

C
C
C  SET TEST IONIC SPECIES ATOMIC AND MOLECULAR DATA;
C
C  STORE "DEFAULT DISSOCIATION MODEL" DATA
C  IN EACH CELL.
C  FOR HYDROGENIC MOLECULE IONS ONLY
C  FOR ALL OTHER SPECIES: INFINITE MFP, I.E. NO DEFAULT COLLISIONS
C
C
      DO 100 IION=1,NIONI
        IDSC1=0
        LGIEI(IION,0)=0
C
        DO NRC=1,NRCI(IION)
          KK=IREACI(IION,NRC)
          IF (ISWR(KK).LE.0.OR.ISWR(KK).GT.6) GOTO 994
        ENDDO
C
C   FIRST: DEAL WITH EI (ELECTRON IMPACT) COLLISIONS
C
        IF (NRCI(IION).EQ.0.AND.NCHARI(IION).EQ.2.AND.
     .      NPRT(NSPAM+IION).GT.1) THEN
C  APPLY THE DEFAULT MODEL FOR H2+ DISSOCIATION
C  USE NPRT.GT.1 TO DISTINGUISH ATOMIC FROM MOLECULAR IONS

C  FIRST: FIND SECONDARY SPECIES INDICES:
          IATM1=0
          IATM2=0
          IPLS1=0
          IPLS2=0
C  H2+:
          IF (NMASSI(IION).EQ.2) THEN
            DO 21 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.1) THEN
                IATM1=IATM
                IATM2=IATM
              ENDIF
   21       CONTINUE
            DO 23 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.1.AND.NCHARP(IPLS).EQ.1) THEN
                IPLS1=IPLS
                IPLS2=IPLS
              ENDIF
   23       CONTINUE
C  HD+:
          ELSEIF (NMASSI(IION).EQ.3) THEN
            DO 31 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.1) THEN
                IATM1=IATM
              ELSEIF (NMASSA(IATM).EQ.2) THEN
                IATM2=IATM
              ENDIF
   31       CONTINUE
            DO 33 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.1.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS1=IPLS
              ELSEIF (NMASSP(IPLS).EQ.2.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS2=IPLS
              ENDIF
   33       CONTINUE
C  D2+:
          ELSEIF (NMASSI(IION).EQ.4) THEN
C  TEST: D2+ OR HT+, USE TEXTS(IION)
            IF (INDEX(TEXTS(NSPAM+IION),'D').NE.0) THEN
C  D2+ TEST ION IDENTIFIED
              DO 41 IATM=1,NATMI
                IF (NMASSA(IATM).EQ.2) THEN
                  IATM1=IATM
                  IATM2=IATM
                ENDIF
   41         CONTINUE
              DO 43 IPLS=1,NPLSI
                IF (NMASSP(IPLS).EQ.2.AND.NCHRGP(IPLS).EQ.1) THEN
                  IPLS1=IPLS
                  IPLS2=IPLS
                ENDIF
   43         CONTINUE
            ELSEIF (INDEX(TEXTS(NSPAM+IION),'H').NE.0.OR.
     .              INDEX(TEXTS(NSPAM+IION),'T').NE.0) THEN
C  HT+ TEST ION IDENTIFIED
              DO 46 IATM=1,NATMI
                IF (NMASSA(IATM).EQ.1) THEN
                  IATM1=IATM
                ELSEIF (NMASSA(IATM).EQ.3) THEN
                  IATM2=IATM
                ENDIF
   46         CONTINUE
              DO 47 IPLS=1,NPLSI
                IF (NMASSP(IPLS).EQ.1.AND.NCHRGP(IPLS).EQ.1) THEN
                  IPLS1=IPLS
                ELSEIF (NMASSP(IPLS).EQ.3.AND.NCHRGP(IPLS).EQ.1) THEN
                  IPLS2=IPLS
                ENDIF
   47         CONTINUE
            ELSE
              CALL EIRENE_LEER(2)
              WRITE (iunout,*) 'TEST ION NO ',IION,
     .                         ' COULD NOT BE IDENTIFIED'
              WRITE (iunout,*) 'NO DEFAULT A&M DATA ASSIGNED'
              CALL EIRENE_LEER(2)
            ENDIF
C  DT+:
          ELSEIF (NMASSI(IION).EQ.5) THEN
            DO 51 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.2) THEN
                IATM1=IATM
              ELSEIF (NMASSA(IATM).EQ.3) THEN
                IATM2=IATM
              ENDIF
   51       CONTINUE
            DO 53 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.2.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS1=IPLS
              ELSEIF (NMASSP(IPLS).EQ.3.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS2=IPLS
              ENDIF
   53       CONTINUE
C  T2+:
          ELSEIF (NMASSI(IION).EQ.6) THEN
            DO 61 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.3) THEN
                IATM1=IATM
                IATM2=IATM
              ENDIF
   61       CONTINUE
            DO 63 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.3.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS1=IPLS
                IPLS2=IPLS
              ENDIF
   63       CONTINUE
          ENDIF
          ITEST=IATM1*IATM2*IPLS1*IPLS2
          IF (ITEST.EQ.0) GOTO 76
C
C  SET DEFAULT MODEL: 3 ELECTRON IMPACT PROCESSES, LABELED -8, -9 AND -10.
C
C  FIRST PROCESS (MAY BE SPLIT INTO 1A AND 1B)  H2+  -->  H + H+ :  DEFAULT PROCESS NO. KK=-8
          KK=-8
          IF (IATM1.NE.IATM2) THEN
            FACTKK=0.5
            ICOUNT=1
          ELSE
            FACTKK=1.D0
            ICOUNT=2
          ENDIF
C
          IA1=IATM1
          IP2=IPLS2
 7000     ACCMAS=0.D0
          ACCINV=0.D0
          IDSC1=IDSC1+1
          NREII=NREII+1
          IREI=NREII
          LGIEI(IION,IDSC1)=IREI
          PATEI(IREI,IA1)=PATEI(IREI,IA1)+1.
          PPLEI(IREI,IP2)=PPLEI(IREI,IP2)+1.
          ACCMAS=ACCMAS+RMASSA(IA1)
          ACCMAS=ACCMAS+RMASSP(IP2)
          ACCINV=ACCINV+1./RMASSA(IA1)
          ACCINV=ACCINV+1./RMASSP(IP2)
          P2ND(IREI,NSPH+IA1)=P2ND(IREI,NSPH+IA1)+1.

          EATEI(IREI,IA1,1)=RMASSA(IA1)/ACCMAS
          EATEI(IREI,IA1,2)=1./RMASSA(IA1)/ACCINV
          EATEI(IREI,0,    1)=EATEI(IREI,IA1,1)
          EATEI(IREI,0,    2)=EATEI(IREI,IA1,2)

          EPLEI(IREI,IP2,1)=RMASSP(IP2)/ACCMAS
          EPLEI(IREI,IP2,2)=1./RMASSP(IP2)/ACCINV
          EPLEI(IREI,0,    1)=RMASSP(IP2)/ACCMAS
          EPLEI(IREI,0,    2)=1./RMASSP(IP2)/ACCINV

          PELEI(IREI)=0.
          MODCOL(1,2,IREI)=1
          MODCOL(1,4,IREI)=1

          IF (NSTORDR >= NRAD) THEN
            DO 73 J=1,NSBOX
              COU = EIRENE_RATE_COEFF(-8,J,TEINL(J),0._DP,.TRUE.,0)
              TABEI1(IREI,J)=COU*DEIN(J)*FACTKK
   73       CONTINUE
            EELEI1(IREI,1:NSBOX)=-10.5
C           EPOTEI(IREI)=1.9   !  default for EDPOTI for this diss. excit. reaction)
C  TRANSFERRED KINETIC ENERGY: 8.6 EV
            EHVEI1(IREI,1:NSBOX)=8.6
            NREAEI(IREI) = -8
            NELREI(IREI) = -8
            NHVREI(IREI) = -8
          ELSE ! storage save mode
            EELEI1(IREI,1)=-10.5
            EHVEI1(IREI,1)=8.6
            NREAEI(IREI) = -8
            NELREI(IREI) = -8
            NHVREI(IREI) = -8
          END IF

          FACREI(IREI,1) = FACTKK
          FACREI(IREI,2) = LOG(FACTKK)
          IF (ICOUNT.EQ.1) THEN
            IA1=IATM2
            IP2=IPLS1
            ICOUNT=2
            GOTO 7000
          ENDIF

C  SECOND PROCESS,  H2+ --> H+ +  H+ + e :  DEFAULT PROCESS NO. KK=-9
cdr   KER = 0.8, ETH = -15.5, I.E. KER=0.4 PER PARTICLE
          KK=-9
          ACCMAS=0.D0
          ACCINV=0.D0
          IDSC1=IDSC1+1
          NREII=NREII+1
          IREI=NREII
          LGIEI(IION,IDSC1)=IREI
          PPLEI(IREI,IPLS1)=PPLEI(IREI,IPLS1)+1.
          PPLEI(IREI,IPLS2)=PPLEI(IREI,IPLS2)+1.
          ACCMAS=ACCMAS+RMASSP(IPLS1)
          ACCMAS=ACCMAS+RMASSP(IPLS2)
          ACCINV=ACCINV+1./RMASSP(IPLS1)
          ACCINV=ACCINV+1./RMASSP(IPLS2)

          EPLEI(IREI,IPLS1,1)=RMASSP(IPLS1)/ACCMAS
          EPLEI(IREI,IPLS2,1)=RMASSP(IPLS2)/ACCMAS     ! if ipls1=ipls2: eplei(ipls,1): only 1/2
          EPLEI(IREI,IPLS1,2)=1./RMASSP(IPLS1)/ACCINV
          EPLEI(IREI,IPLS2,2)=1./RMASSP(IPLS2)/ACCINV  ! if ipls1=ipls2: eplei(ipls,2): only 1/2
          EPLEI(IREI,0,    1)=EPLEI(IREI,IPLS1,1)+EPLEI(IREI,IPLS2,1)  ! if ipls1=ipls2: eplei(0,1): total, correct
          EPLEI(IREI,0,    2)=EPLEI(IREI,IPLS1,2)+EPLEI(IREI,IPLS2,2)  ! if ipls1=ipls2: eplei(0,2): total, correct

          PELEI(IREI)=1.

          MODCOL(1,2,IREI)=1
          MODCOL(1,4,IREI)=1
C
          IF (NSTORDR >= NRAD) THEN
            DO 71 J=1,NSBOX
              COU = EIRENE_RATE_COEFF(-9,J,TEINL(J),0._DP,.TRUE.,0)
              TABEI1(IREI,J)=COU*DEIN(J)
   71       CONTINUE
C  NO RADIATION LOSS INCLUDED
            EELEI1(IREI,1:NSBOX)=-15.5
C           EPOTEI(IREI) = 14.7!  default for EDPOTI for this diss ionis. reaction)
C  TRANSFERRED KINETIC ENERGY: 0.8  (=0.4 EV PER PARTICLE)
            EHVEI1(IREI,1:NSBOX)=0.8
            NREAEI(IREI) = -9
            NELREI(IREI) = -9   !  electron energy loss model
            NHVREI(IREI) = -9   !  KER model for process KK=-9
          ELSE  ! storage save mode
            EELEI1(IREI,1)=-15.5
            EHVEI1(IREI,1)=0.8
            NREAEI(IREI) = -9
            NELREI(IREI) = -9
            NHVREI(IREI) = -9
          END IF
          FACREI(IREI,1) = 1._DP
          FACREI(IREI,2) = 0._DP

C  THIRD PROCESS   H2+   -->   H + H:  DEFAULT PROCESS NO. KK=-10, diss.rec
          KK=-10
          ACCMAS=0.D0
          ACCINV=0.D0
          IDSC1=IDSC1+1
          NREII=NREII+1
          IREI=NREII
          LGIEI(IION,IDSC1)=IREI
          PATEI(IREI,IATM1)=PATEI(IREI,IATM1)+1.
          PATEI(IREI,IATM2)=PATEI(IREI,IATM2)+1.
          ACCMAS=ACCMAS+RMASSA(IATM1)
          ACCMAS=ACCMAS+RMASSA(IATM2)
          ACCINV=ACCINV+1./RMASSA(IATM1)
          ACCINV=ACCINV+1./RMASSA(IATM2)
          P2ND(IREI,NSPH+IATM1)=P2ND(IREI,NSPH+IATM1)+1.
          P2ND(IREI,NSPH+IATM2)=P2ND(IREI,NSPH+IATM2)+1.

          EATEI(IREI,IATM1,1)=RMASSA(IATM1)/ACCMAS
          EATEI(IREI,IATM2,1)=RMASSA(IATM2)/ACCMAS   ! problem is iatm1=iatm2. Then eatds(..iatm,..) is per particle.
          EATEI(IREI,IATM1,2)=1./RMASSA(IATM1)/ACCINV
          EATEI(IREI,IATM2,2)=1./RMASSA(IATM2)/ACCINV

cdr  eatei, emlei, eioei(...,.,...) are secondary energy per secondary particle iat, iml, iio.
cdr  eatei, emlei, eioei(...,0,...) is total per species
cdr  this is different now for eplei.
cdr  eplei  is secondary energy per secondary bulk species, i.e. twice, if two secondaries of same species.
cdr  eatei, ...is used for sampling energy, veloel. eplei is used for scoring bulk energy tallies.

          EATEI(IREI,0,    1)=EATEI(IREI,IATM1,1)+EATEI(IREI,IATM2,1)
          EATEI(IREI,0,    2)=EATEI(IREI,IATM1,2)+EATEI(IREI,IATM2,2)

          PELEI(IREI)=-1.
          MODCOL(1,2,IREI)=1
          MODCOL(1,4,IREI)=1
C
          IF (NSTORDR >= NRAD) THEN
            DO 72 J=1,NSBOX
              COU = EIRENE_RATE_COEFF(-10,J,TEINL(J),0._DP,.TRUE.,0)
              TABEI1(IREI,J)=COU*DEIN(J)
   72       CONTINUE
C  NO RADIATION LOSS INCLUDED
C  FOR THE FACTOR -0.896... SEE: EIRENE MANUAL, INPUT BLOCK 4, EXAMPLES
            DE_10=8.964355004318D-01
            EELEI1(IREI,1:NSBOX)=-DE_10*TEIN(1:NSBOX)
C  TRANSFERRED KINETIC ENERGY: = INGOING ELECTRON ENERGY
            EHVEI1(IREI,1:NSBOX)=DE_10*TEIN(1:NSBOX)
            NREAEI(IREI) = -10
            NELREI(IREI) = -10
            NHVREI(IREI) = -10
          ELSE ! storage save mode
cdr
cdr here: eelds1 und ehvds1 set in felee1 ?  there 0.88 times tein, i.e. not constant.
cdr
            NREAEI(IREI) = -10
            NELREI(IREI) = -10
            NHVREI(IREI) = -10
          END IF
          FACREI(IREI,1) = 1._DP
          FACREI(IREI,2) = 0._DP
C
   76     CONTINUE
C
          NIEII(IION)=IDSC1
C
C
C  NON-DEFAULT ELEC. IMP. COLLISION MODEL SPECIFIED IN INPUT BLOCK 4
C
        ELSEIF (NRCI(IION).GT.0) THEN
          DO 90 NRC=1,NRCI(IION)
            KK=IREACI(IION,NRC)
            IF (ISWR(KK).NE.1) GOTO 90
C
            FACTKK=FREACI(IION,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
            CHRDF0=-NCHRGI(IION)
            IIO=NSPAM+IION
            RMASS=RMASSI(IION)
            IFRST=ISCD1I(IION,NRC)
            ISCND=ISCD2I(IION,NRC)
            ITHRD=ISCD3I(IION,NRC)
            IFRTH=ISCD4I(IION,NRC)
            ISCDE=ISCDEI(IION,NRC)
            IESTM=IESTMI(IION,NRC)
            EHEAVY=ESCD1I(IION,NRC)
            EELEC=EELECI(IION,NRC)
            IDSC1=IDSC1+1
            NREII=NREII+1
            IREI=NREII
            LGIEI(IION,IDSC1)=IREI
            CALL EIRENE_XSTEI(RMASS,IREI,IIO,
     .                 IFRST,ISCND,ITHRD,IFRTH,EHEAVY,CHRDF0,
     .                 ISCDE,EELEC,IESTM,
     .                 KK,FACTKK,PLS)
   90     CONTINUE
          NIEII(IION)=IDSC1
        ENDIF
C
        NIEIIM(IION)=NIEII(IION)-1
        LGIEI(IION,0)=NIEII(IION)
C
        DO IIEI=1,NIEII(IION)
          IREI=LGIEI(IION,IIEI)
          CALL EIRENE_XSTEI_1(IREI)
        ENDDO


  100 CONTINUE

C   SECONDLY: DEAL WITH CX (CHARGE EXCHANGE) COLLISIONS

C  TENTATIVELY ASSUME: NO CHARGE EXCHANGE BETWEEN IION AND ANY IPLS
      DO 200 IION=1,NIONI
        IDSC=0
        LGICX(IION,0,0)=0
        LGICX(IION,0,1)=0
C
C  THERE ARE CURRENTLY NO DEFAULT CX RATES FOR TEST IONS
C
        IF (NRCI(IION).EQ.0) THEN
          NICXI(IION)=0
C
C  NON-DEFAULT CX MODEL:
C
        ELSEIF (NRCI(IION).GT.0) THEN
          DO 130 NRC=1,NRCI(IION)
            KK=IREACI(IION,NRC)
            IF (ISWR(KK).NE.3) CYCLE
C  make sure that incident particle is a bulk particle
            IF (EIRENE_IDEZ(IBULKI(IION,NRC),1,3).NE.4) THEN
C  WRONG TYPE OF INCIDENT BULK SPECIES
              WRITE (IUNOUT,*)
     .        'INPUT ERROR FOR CX PROCESS, IION,KK ',IION,KK
              CALL EIRENE_EXIT_OWN(1)
            ENDIF
C  CX PROCESS IDENTIFIED

            FACTKK=FREACI(IION,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
            CHRDF0=-NCHRGI(IION)
C  BULK PARTICLE INDEX
            IPLS=EIRENE_IDEZ(IBULKI(IION,NRC),3,3)
            IDSC=IDSC+1
            NRCXI=NRCXI+1
            IRCX=NRCXI
            LGICX(IION,IDSC,0)=IRCX
            LGICX(IION,IDSC,1)=IPLS

            IIO=NSPAM+IION
            IPL=IPLS
            RMASS=RMASSI(IION)
            IFRST=ISCD1I(IION,NRC)
            ISCND=ISCD2I(IION,NRC)
            ISCDE=ISCDEI(IION,NRC)
            IESTM=IESTMI(IION,NRC)
            EBULK=EBULKI(IION,NRC)
            CALL EIRENE_XSTCX(RMASS,IRCX,IIO,IPL,
     .                        IFRST,ISCND,EBULK,CHRDF0,
     .                        ISCDE,IESTM,KK,FACTKK,PLS)
C
  130     CONTINUE
C
          NICXI(IION)=IDSC
C  NO CX MODEL DEFINED
        ELSE
          NICXI(IION)=0
        ENDIF
C
        NICXIM(IION)=NICXI(IION)-1

        LGICX(IION,0,0)=0
        DO IICX=1,NICXI(IION)
          LGICX(IION,0,0)=LGICX(IION,0,0)+LGICX(IION,IICX,0)
        ENDDO
C
  200 CONTINUE
C
C   ELASTIC COLLISIONS
C
      DO 300 IION=1,NIONI
        IDSC=0
        LGIEL(IION,0,0)=0
        LGIEL(IION,0,1)=0
C
C  DEFAULT EL MODEL: NOT AVAILABLE
C
        IF (NRCI(IION).EQ.0) THEN
          NIELI(IION)=0
C
C  NON-DEFAULT EL MODEL:  240--
C
        ELSEIF (NRCI(IION).GT.0) THEN
          DO 230 NRC=1,NRCI(IION)
            KK=IREACI(IION,NRC)
            IF (ISWR(KK).NE.5) GOTO 230
C
            FACTKK=FREACI(IION,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
C  BULK PARTICLE INDEX
            IPLS=EIRENE_IDEZ(IBULKI(IION,NRC),3,3)
            IF (IPLS.LE.0.OR.IPLS.GT.NPLSI) GOTO 991
            IF (MASSP(KK).LE.0.OR.MASST(KK).LE.0) GOTO 992
            IDSC=IDSC+1
            NRELI=NRELI+1
            IREL=NRELI
            LGIEL(IION,IDSC,0)=IREL
            LGIEL(IION,IDSC,1)=IPLS
C
C  SPECIAL TREATMENT: BGK COLLISIONS AMONGST TEST PARTICLES
C  FOR THIS REACTION KK
            IF (IBGKI(IION,NRC).NE.0) THEN
              IF (NPBGKI(IION).EQ.0) THEN
C  IION HAS NOT YET BEEN LABELLED AS BGK SPECIES.
C  DO THIS HERE: IION IS BGK SPECIES NO. IBGK_SP, AND HAS 3 ADDITIONAL BGK TALLIES IN UPTBGK
                NRBGI=NRBGI+3
                IBGK_SP=NRBGI/3
                NPBGKI(IION)=IBGK_SP
              ENDIF
              IF (NPBGKP(IPLS,1).EQ.0) THEN
                NPBGKP(IPLS,1)=NPBGKI(IION)
cdr this is too restrictive?
              ELSE
                GOTO 999
              ENDIF
C  SELF- OR CROSS-COLLISION?
              ITYPB=EIRENE_IDEZ(IBGKI(IION,NRC),1,3)
              ISPZB=EIRENE_IDEZ(IBGKI(IION,NRC),3,3)
              IF (ITYPB.NE.3.OR.ISPZB.NE.IION) THEN
C  CROSS-COLLISION ! SET THE SECOND TEST PARTICLE SPECIES
C                    INVOLVED IN THIS PROCESS
                IF (NPBGKP(IPLS,2).EQ.0) THEN
                  NPBGKP(IPLS,2)=IBGKI(IION,NRC)
                ELSE
                  GOTO 999
                ENDIF
              ENDIF
            ENDIF
C  BGK COLLISION PARAMETERS DONE
C
            III=NSPAM+IION
            IPL=IPLS
            ISCDE=ISCDEI(IION,NRC)
            IESTM=IESTMI(IION,NRC)
            EBULK=EBULKI(IION,NRC)
            CALL EIRENE_XSTEL(IREL,III,IPL,EBULK,
     .                 ISCDE,IESTM,
     .                 KK,FACTKK,PLS)
C
  230     CONTINUE

          NIELI(IION)=IDSC
        ENDIF
C
        NIELIM(IION)=NIELI(IION)-1
C
        LGIEL(IION,0,0)=0.
        DO 280 IIEL=1,NIELI(IION)
          LGIEL(IION,0,0)=LGIEL(IION,0,0)+LGIEL(IION,IIEL,0)
  280   CONTINUE
C
C
  300 CONTINUE
C
C   GENERAL HEAVY PARTICLE IMPACT COLLISIONS
C

      DO IION=1,NIONI
        IDSC=0
        LGIPI(IION,0,0)=0
        LGIPI(IION,0,1)=0
C
C  NO DEFAULT MODEL
C
        IF (NRCI(IION).EQ.0) THEN
          NIPII(IION)=0
C
C  NON-DEFAULT ION IMPACT MODEL:  130--190
C
        ELSEIF (NRCI(IION).GT.0) THEN
          DO NRC=1,NRCI(IION)
            KK=IREACI(IION,NRC)
            IF (ISWR(KK).NE.4) CYCLE
            FACTKK=FREACI(IION,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
            IF (MASSP(KK).LE.0.OR.MASST(KK).LE.0) GOTO 992
C  INCIDENT BULK PARTICLE INDEX
            IPLS=EIRENE_IDEZ(IBULKI(IION,NRC),3,3)
            CHRDF0=-NCHRGI(IION)
            IF (IPLS.LE.0.OR.IPLS.GT.NPLSI) GOTO 990
            IDSC=IDSC+1
            NRPII=NRPII+1
            IF (NRPII.GT.NRPI) GOTO 998
            IRPI=NRPII
            LGIPI(IION,IDSC,0)=IRPI
            LGIPI(IION,IDSC,1)=IPLS

            III=NSPAM+IION
            IPL=IPLS
            RMASS=RMASSI(IION)
            IFRST=ISCD1I(IION,NRC)
            ISCND=ISCD2I(IION,NRC)
            ITHRD=ISCD3I(IION,NRC)
            IFRTH=ISCD4I(IION,NRC)
            ISCDE=ISCDEI(IION,NRC)
            IESTM=IESTMI(IION,NRC)
            EBULK=EBULKI(IION,NRC)
            EHEAVY=ESCD1I(IION,NRC)
            EELEC=EELECI(IION,NRC)
            CALL EIRENE_XSTPI (RMASS,IRPI,III,IPL,
     .                  IFRST,ISCND,ITHRD,IFRTH,
     .                  EBULK,EHEAVY,EELEC,CHRDF0,ISCDE,IESTM,
     .                  KK,FACTKK,PLS)
          END DO
C
          NIPII(IION)=IDSC
C  NO MODEL DEFINED
        ELSE
          NIPII(IION)=0
        ENDIF
C
        NIPIIM(IION)=NIPII(IION)-1
C
        LGIPI(IION,0,0)=0
        DO IIPI=1,NIPII(IION)
          LGIPI(IION,0,0)=LGIPI(IION,0,0)+LGIPI(IION,IIPI,0)
        END DO

        DO IIPI=1,NIPII(IION)
          IRPI=LGIPI(IION,IIPI,0)
          CALL EIRENE_XSTPI_1(IRPI)
        END DO
      END DO
C
C
      DO 1000 IION=1,NIONI
C
        IF (TRCAMD) THEN
          CALL EIRENE_MASBOX
     .  ('TEST ION SPECIES IION = '//TEXTS(NSPAM+IION))
          CALL EIRENE_LEER(1)
C
          IF (LGIEI(IION,0).EQ.0) THEN

            CALL EIRENE_LEER(1)
            WRITE (iunout,*) 'NO ELECTRON IMPACT COLLISIONS'
            CALL EIRENE_LEER(1)
          ELSE
            DO 870 IIEI=1,NIEII(IION)
              IREI=LGIEI(IION,IIEI)
              CALL EIRENE_XSTEI_2(IREI)
  870       CONTINUE
          ENDIF
C
C
          CALL EIRENE_LEER(2)
          IF (LGICX(IION,0,0).EQ.0) THEN
            CALL EIRENE_LEER(1)
            WRITE (iunout,*) 'NO CHARGE EXCHANGE WITH BULK IONS'
            CALL EIRENE_LEER(1)
          ELSE
            DO 890 IICX=1,NICXI(IION)
              IRCX=LGICX(IION,IICX,0)
              IPL =LGICX(IION,IICX,1)
              CALL EIRENE_XSTCX_2(IRCX,IPL)
  890       CONTINUE
          ENDIF
C
C
          CALL EIRENE_LEER(2)
          IF (LGIEL(IION,0,0).EQ.0) THEN
            CALL EIRENE_LEER(1)
            WRITE (iunout,*) 'NO ELASTIC COLLISIONS WITH BULK IONS'
            CALL EIRENE_LEER(1)
          ELSE
            DO 895 IIEL=1,NIELI(IION)
              IREL=LGIEL(IION,IIEL,0)
              IPL =LGIEL(IION,IIEL,1)
              CALL EIRENE_XSTEL_2(IREL,IPL)
  895       CONTINUE
          ENDIF
C
          CALL EIRENE_LEER(2)
          IF (LGIPI(IION,0,0).EQ.0) THEN
            CALL EIRENE_LEER(1)
            WRITE (iunout,*) 'NO GENERAL ION IMPACT COLLISIONS '
            CALL EIRENE_LEER(1)
          ELSE
            DO 885 IIPI=1,NIPII(IION)
              IRPI=LGIPI(IION,IIPI,0)
              IPL =LGIPI(IION,IIPI,1)
              CALL EIRENE_XSTPI_2(IRPI,IPL)
  885       CONTINUE
          ENDIF

        ENDIF
C
 1000 CONTINUE

      DEALLOCATE (PLS)
C
      RETURN
C
  990 CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTI: EXIT CALLED'
      WRITE (iunout,*) 'INVALID SPECIES INDEX FOR ION IMPACT COLLISION'
      CALL EIRENE_EXIT_OWN(1)
  991 CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTI: EXIT CALLED'
      WRITE (iunout,*) 'INVALID SPECIES INDEX FOR ELASTIC COLLISION'
      CALL EIRENE_EXIT_OWN(1)
  992 CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTI: EXIT CALLED'
      WRITE (iunout,*)
     .  'MASS NUMBERS OF INTERACTING PARTICLES INCONSISTENT'
      WRITE (iunout,*) 'KK,IION,IPLS ',KK,IION,IPLS
      CALL EIRENE_EXIT_OWN(1)
  994 CONTINUE
      WRITE (iunout,*) 'ERROR DETECTED IN XSECTI.'
      WRITE (iunout,*) 'REACTION NO. KK= ',KK, 'NOT READ FROM FILE '
      WRITE (iunout,*) 'IION = ',IION
      WRITE (iunout,*) 'ISWR(KK) = ',ISWR(KK)
      WRITE (iunout,*) 'EXIT CALLED'
      CALL EIRENE_EXIT_OWN(1)
  998 CONTINUE
      WRITE (iunout,*) 'INSUFFICIENT STORAGE FOR PI: NRPI=',NRPI
      CALL EIRENE_EXIT_OWN(1)
  999 CONTINUE
      WRITE (iunout,*) 'SPECIES CONFLICT FOR BGK COLLISIONS. IION,IREL'
      WRITE (iunout,*) IION,IREL,IPLS
      CALL EIRENE_EXIT_OWN(1)
      RETURN
C
      END SUBROUTINE EIRENE_XSECTI
