C 27.6.05  irds --> irei
c 24.11.05 use nprt(ispz) to check if imol is really a molecule.
c          otherwise He atoms could be confused with d2 molecules,
c          if they accidentally are specified in the molecule block 4b

c 24.11.05: chrdf0 in parameterlist for call to xstcx
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
! 23.02.14: call to xstpi: additional arguments: IML, pls (for H.4 option)
! oct.2014: call to xstpi: additional argument: chrdf0
cdr  oct.14:  clogau removed
cdr  oct.14:  PLS made allocatable,
cdr  oct.14:  further synchronization with xsecta,xsecti
cdr           remaining relevant differences in default models only.
cdr  aug.15:  ibgk_sp:  no of bgk species. to be distinguished from ibgk: no of bgk reaction.

!pb  APR  16:  pplds  -> pplei
!pb  APR  16:  patds  -> patei, eatds -> eatei
!pb  APR  16:  piods  -> pioei, eiods -> eioei
!pb  APR  16:  pelds  -> pelei, eelds -> eelei
!pb  MAY  16:  tabds1 -> tabds1
!pb  JUL  16:  ehvds1 -> ehvds1
cdr  Sept 16:  nmdsi  -> nmeii
cdr  May 18:  The fluid limit (critical cx Knudsen number) is now set from NGENM(imol) flag,
cdr           rather than from the former fldlmm(imol,kk) flag (which is removed now).
cdr           default: FDLMCX=0.0 (from initialisation phase) means:
CDR           no fluid limit cut-off at CX collisions.
cdr  sept 18: nhvrei rationalization for default reactions (==-KK)
C

      SUBROUTINE EIRENE_XSECTM
C
C       SET UP TABLES (E.G. OF REACTION RATES) FOR MOLECULAR SPECIES
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
      REAL(DP) :: FACTKK, DEIMIN, EELEC, CHRDF0,
     .            RMASS, EBULK, EHEAVY, COU, EIRENE_RATE_COEFF,
     .            ACCMAS, ACCINV

      INTEGER :: ITEST, IATM, IPLS, IION, IA1, IP2, ION, ICOUNT,
     .           IION3, IDSC1, NRC, KK, J, IMOL, IPLS1, IPLS2, IPLS3,
     .           IATM1, IATM2, ITYPB, ISPZB, IMEL, IDSC, IREL, IBGK_SP,
     .           IMEI, IMCX, ISCND, ISCDE, IESTM, IML, IFRST,
     .           IRCX, IREI, IPL, IMPI, IRPI, ITHRD, IFRTH,
     .           MFL
      INTEGER, EXTERNAL :: EIRENE_IDEZ

      ALLOCATE (PLS(NSTORDR))

cdr  PLS:  ELECTRON DENSITY PARAMETER in CR MODELS
cdr       (NOT TO BE CONFUSED WITH THE DENSITY FACTOR BETWEEN RATES AND RATE COEFF.)
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
C   ELECTRON IMPACT COLLISIONS:
C

C
C
      DO 100 IMOL=1,NMOLI
        IDSC1=0
        LGMEI(IMOL,0)=0
C
        DO NRC=1,NRCM(IMOL)
          KK=IREACM(IMOL,NRC)
          IF (ISWR(KK).LE.0.OR.ISWR(KK).GT.6) GOTO 994
        ENDDO
C
C  CHECK IF THIS REALLY IS A MOLECULE: USE NPRT(ISPZ).GT.1?

        IF (NPRT(NSPA+IMOL).LE.1) THEN
          WRITE (IUNOUT,*) 'SEVERE INPUT ERROR DETECTED IN XSECTM: '
          WRITE (IUNOUT,*) 'IMOL= ',IMOL,' CARRIES ONLY ONE FLUX UNIT'
          WRITE (IUNOUT,*) 'EXIT CALLED FROM XSECTM '
          CALL EIRENE_EXIT_OWN(1)
        ENDIF

C  YES, "IMOL" IS A MOLECULE !
C
        IF (NRCM(IMOL).EQ.0.AND.NCHARM(IMOL).EQ.2) THEN

C  APPLY THE DEFAULT MODEL FOR H2 DISSOCIATION AND IONIZATION
C  TO ALL HYDROGENIC MOLECULES

C   ELECTRON IMPACT COLLISIONS:
C
C  STORE "DEFAULT DISSOCIATION MODEL" DATA
C  IN EACH CELL. THIS DEFAULT MODEL (JANEV/LANGER, PPPL) MAY BE USED
C  FOR HYDROGENIC MOLECULES ONLY.
C
          EELEC=-EIONH2
C
C  FIRST: FIND SECONDARY SPECIES INDICES:
          IATM1=0
          IATM2=0
          IPLS1=0
          IPLS2=0
          IPLS3=0
          IION3=0
C  H2:
          IF (NMASSM(IMOL).EQ.2) THEN
            DO 21 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.1) THEN
                IATM1=IATM
                IATM2=IATM
              ENDIF
   21       CONTINUE
            DO 23 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.1.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS1=IPLS
                IPLS2=IPLS
              ENDIF
   23       CONTINUE
            DO 25 IION=1,NIONI
              IF (NMASSI(IION).EQ.2.AND.NCHARI(IION).EQ.2) THEN
                IION3=IION
              ENDIF
   25       CONTINUE
C  HD:
          ELSEIF (NMASSM(IMOL).EQ.3) THEN
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
            DO 35 IION=1,NIONI
              IF (NMASSI(IION).EQ.3.AND.NCHARI(IION).EQ.2) THEN
                IION3=IION
              ENDIF
   35       CONTINUE
C  D2:  (OR HT ?)
          ELSEIF (NMASSM(IMOL).EQ.4) THEN
C  TEST: D2 OR HT, USE TEXTS(IMOL)
            IF (INDEX(TEXTS(NSPA+IMOL),'D').NE.0) THEN
C  D2 MOLECULE IDENTIFIED
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
              DO 45 IION=1,NIONI
                IF (NMASSI(IION).EQ.4.AND.NCHARI(IION).EQ.2.AND.
     .              INDEX(TEXTS(NSPAM+IION),'D').NE.0) THEN
                  IION3=IION
                ENDIF
   45         CONTINUE
            ELSEIF (INDEX(TEXTS(NSPA+IMOL),'H').NE.0.OR.
     .              INDEX(TEXTS(NSPA+IMOL),'T').NE.0) THEN
C  HT MOLECULE IDENTIFIED
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
              DO 48 IION=1,NIONI
                IF (NMASSI(IION).EQ.4.AND.NCHARI(IION).EQ.2.AND.
     .             (INDEX(TEXTS(NSPAM+IION),'H').NE.0.OR.
     .              INDEX(TEXTS(NSPAM+IION),'T').NE.0)) THEN
                  IION3=IION
                ENDIF
   48         CONTINUE
            ELSE
              CALL EIRENE_LEER(2)
              WRITE (iunout,*) 'MOLECULE NO ',IMOL,
     .                         ' COULD NOT BE IDENTIFIED'
              WRITE (iunout,*) 'NO DEFAULT A&M DATA ASSIGNED'
              CALL EIRENE_LEER(2)
            ENDIF
C  DT:
          ELSEIF (NMASSM(IMOL).EQ.5) THEN
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
            DO 55 IION=1,NIONI
              IF (NMASSI(IION).EQ.5.AND.NCHARI(IION).EQ.2) THEN
                IION3=IION
              ENDIF
   55       CONTINUE
C  T2:
          ELSEIF (NMASSM(IMOL).EQ.6) THEN
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
            DO 65 IION=1,NIONI
              IF (NMASSI(IION).EQ.6.AND.NCHARI(IION).EQ.2) THEN
                IION3=IION
              ENDIF
   65       CONTINUE
          ENDIF

          ITEST=IATM1*IATM2*IPLS1*IPLS2*IION3
          IF (ITEST.EQ.0) GOTO 76
C
C  SET DEFAULT MODEL: 3 ELECTRON IMPACT PROCESSES, DEFAULT PROCESSES KK=-5,-6, -7
C
C  FIRST PROCESS, KK=-5   H2 --> H + H:  DEFAULT PROCESS NO KK=-5
          KK=-5
          ACCMAS=0.D0
          ACCINV=0.D0
          IDSC1=IDSC1+1
          NREII=NREII+1
          IREI=NREII
          LGMEI(IMOL,IDSC1)=IREI
          PATEI(IREI,IATM1)=PATEI(IREI,IATM1)+1.
          PATEI(IREI,IATM2)=PATEI(IREI,IATM2)+1.
          ACCMAS=ACCMAS+RMASSA(IATM1)
          ACCMAS=ACCMAS+RMASSA(IATM2)
          ACCINV=ACCINV+1./RMASSA(IATM1)
          ACCINV=ACCINV+1./RMASSA(IATM2)
          P2ND(IREI,NSPH+IATM1)=P2ND(IREI,NSPH+IATM1)+1.
          P2ND(IREI,NSPH+IATM2)=P2ND(IREI,NSPH+IATM2)+1.
          EATEI(IREI,IATM1,1)=RMASSA(IATM1)/ACCMAS
          EATEI(IREI,IATM2,1)=RMASSA(IATM2)/ACCMAS
          EATEI(IREI,IATM1,2)=1./RMASSA(IATM1)/ACCINV
          EATEI(IREI,IATM2,2)=1./RMASSA(IATM2)/ACCINV
          EATEI(IREI,0,    1)=EATEI(IREI,IATM1,1)+EATEI(IREI,IATM2,1)
          EATEI(IREI,0,    2)=EATEI(IREI,IATM1,2)+EATEI(IREI,IATM2,2)
          PELEI(IREI)=0.
          MODCOL(1,2,IREI)=1
          MODCOL(1,4,IREI)=1

          IF (NSTORDR >= NRAD) THEN
            DO 70 J=1,NSBOX
              COU = EIRENE_RATE_COEFF(-5,J,TEINL(J),0._DP,.TRUE.,0)
              TABEI1(IREI,J)=COU*DEIN(J)
   70       CONTINUE
            EELEI1(IREI,1:NSBOX)=-10.5
C  TRANSFERRED KINETIC ENERGY: 6 EV
            EHVEI1(IREI,1:NSBOX)=6.0
C           EPOTEI(IREI)=4.5   !  default for EDPOTM for this dissoc. reaction)
            NREAEI(IREI)=-5
            NELREI(IREI)=-5  ! FLAG FOR FEELEI1, FOR DEFAULT REACTION -5:
            NHVREI(IREI)=-5
          ELSE
            EELEI1(IREI,1)=-10.5
            NREAEI(IREI)=-5
            NELREI(IREI)=-5  ! FLAG FOR FEELEI1, FOR DEFAULT REACTION -5:
            NHVREI(IREI)=-5
          END IF
          FACREI(IREI,1) = 1._DP
          FACREI(IREI,2) = 0._DP

C  SECOND PROCESS (MAY BE SPLIT INTO 2A AND 2B)  H2 -->  H  + H+: DEFAULT PROCESS NO. KK=-6
          KK=-6
          IF (IATM1.NE.IATM2) THEN
c   e.g. DT -->  0.5 (D + T+)  + 0.5 (D+ + T)
            FACTKK=0.5
            ICOUNT=1
          ELSE
            FACTKK=1.D0
            ICOUNT=2
          ENDIF
C
          IA1=IATM1
          IP2=IPLS2
c   in case iatm1 ne iatm2:  this next segment is executed twice.
c   Split reaction kk=-6 into two ei processes irei and irei+1, with factkk=0.5 each.
c   Accumulate totals....
   73     ACCMAS=0.D0
          ACCINV=0.D0
          IDSC1=IDSC1+1
          NREII=NREII+1
          IREI=NREII
          LGMEI(IMOL,IDSC1)=IREI
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

          EPLEI(IREI,IP2,  1)=RMASSP(IP2)/ACCMAS
          EPLEI(IREI,IP2,  2)=1./RMASSP(IP2)/ACCINV
          EPLEI(IREI,0,    1)=EPLEI(IREI,IP2,1)
          EPLEI(IREI,0,    2)=EPLEI(IREI,IP2,2)

          PELEI(IREI)=1.0

          MODCOL(1,2,IREI)=1
          MODCOL(1,4,IREI)=1

          IF (NSTORDR >= NRAD) THEN
            DO 71 J=1,NSBOX
              COU = EIRENE_RATE_COEFF(-6,J,TEINL(J),0._DP,.TRUE.,0)
              TABEI1(IREI,J)=COU*DEIN(J)*FACTKK
   71       CONTINUE
            EELEI1(IREI,1:NSBOX)=-25.0
C           EPOTEI(IREI)=15.00   !  default for EDPOTM for this diss ionis. reaction)
C  TRANSFERRED KINETIC ENERGY: 10 EV
            EHVEI1(IREI,1:NSBOX)=10.0
            NREAEI(IREI) = -6
            NELREI(IREI) = -6  ! FLAG FOR FEELEI1, FOR DEFAULT REACTION -6:
            NHVREI(IREI) = -6
          ELSE
            EELEI1(IREI,1)=-25.0
C           EHVEI1(IREI,1)= 10.0 SET IN ....?
            NREAEI(IREI) = -6
            NELREI(IREI) = -6  ! FLAG FOR FEELEI1, FOR DEFAULT REACTION -6:
            NHVREI(IREI) = -6
          END IF
          FACREI(IREI,1) = FACTKK
          FACREI(IREI,2) = LOG(FACTKK)
          IF (ICOUNT.EQ.1) THEN
            IA1=IATM2
            IP2=IPLS1
            ICOUNT=2
            GOTO 73
          ENDIF
C
C  THIRD PROCESS  H2 --> H2+:  DEFAULT PROCESS NO. KK=-7
          KK=-7
          IDSC1=IDSC1+1
          NREII=NREII+1
          IREI=NREII
          LGMEI(IMOL,IDSC1)=IREI
          ION=NSPAM+IION3
          PIOEI(IREI,IION3)=PIOEI(IREI,IION3)+1.
          P2ND(IREI,ION)=P2ND(IREI,ION)+1.
          EIOEI(IREI,IION3,1)=1.
          EIOEI(IREI,IION3,2)=0.
          EIOEI(IREI,0,1)=1.
          EIOEI(IREI,0,2)=0.
          PELEI(IREI)=1.0

          MODCOL(1,2,IREI)=1
          MODCOL(1,4,IREI)=1
C
          IF (NSTORDR >= NRAD) THEN
            DO 72 J=1,NSBOX
              COU = EIRENE_RATE_COEFF(-7,J,TEINL(J),0._DP,.TRUE.,0)
              TABEI1(IREI,J)=COU*DEIN(J)
   72       CONTINUE
C  NO RADIATION LOSS INCLUDED
            EELEI1(IREI,1:NSBOX)=EELEC  ! =-EIONH2 = -15.45 EV
C           EPOTEI(IREI)=15.45   !  default for EDPOTM for this ionis. reaction)
C           EHVEI1(IREI,1:NSBOX)=0.0
C  PROBABLY NOT NEEDED, ONLY IN STORAGE SAVING MODE
            NREAEI(IREI) = -7  ! FLAG FOR FTABEI1, FOR DEFAULT REACTION -7
C  PROBABLY NOT NEEDED, ONLY IN STORAGE SAVING MODE
            NELREI(IREI) = -7  ! FLAG FOR FEELEI1, FOR DEFAULT REACTION -7:
            NHVREI(IREI) = -7
          ELSE  ! storage save mode
            EELEI1(IREI,1)=EELEC   ! =-EIONH2 = -15.45 EV
C           EHVEI1(IREI,1)= 0.0 SET IN ....?
            NREAEI(IREI) = -7  ! FLAG FOR FTABEI1, FOR DEFAULT REACTION -7
            NELREI(IREI) = -7  ! FLAG FOR FEELEI1, FOR DEFAULT REACTION -7:
            NHVREI(IREI) = -7
          END IF
          FACREI(IREI,1) = 1._DP
          FACREI(IREI,2) = 0._DP
C
   76     CONTINUE

          NMEII(IMOL)=IDSC1
C
C  NON-DEFAULT ELEC IMP. COLLISION MODEL SPECIFIED IN INPUT BLOCK 4
C
        ELSEIF (NRCM(IMOL).GT.0) THEN
          DO 90 NRC=1,NRCM(IMOL)
            KK=IREACM(IMOL,NRC)
            IF (ISWR(KK).NE.1) GOTO 90
C
C  EI PROCESS IDENTIFIED

            FACTKK=FREACM(IMOL,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
            CHRDF0=0.D0
            IML=NSPA+IMOL
            RMASS=RMASSM(IMOL)
            IFRST=ISCD1M(IMOL,NRC)
            ISCND=ISCD2M(IMOL,NRC)
            ITHRD=ISCD3M(IMOL,NRC)
            IFRTH=ISCD4M(IMOL,NRC)
            ISCDE=ISCDEM(IMOL,NRC)
            IESTM=IESTMM(IMOL,NRC)
            EHEAVY=ESCD1M(IMOL,NRC)
            EELEC=EELECM(IMOL,NRC)
            IDSC1=IDSC1+1
            NREII=NREII+1
            IREI=NREII
            LGMEI(IMOL,IDSC1)=IREI
            CALL EIRENE_XSTEI(RMASS,IREI,IML,
     .                 IFRST,ISCND,ITHRD,IFRTH,EHEAVY,CHRDF0,
     .                 ISCDE,EELEC,IESTM,KK,FACTKK,PLS)
   90     CONTINUE
          NMEII(IMOL)=IDSC1
        ENDIF
C
        NMEIIM(IMOL)=NMEII(IMOL)-1
        LGMEI(IMOL,0)=NMEII(IMOL)
C
        DO IMEI=1,NMEII(IMOL)
          IREI=LGMEI(IMOL,IMEI)
          CALL EIRENE_XSTEI_1(IREI)
        ENDDO


  100 CONTINUE
C
C
C   CHARGE EXCHANGE:
C
C  TENTATIVELY ASSUME: NO CHARGE EXCHANGE BETWEEN IATM AND ANY IPLS
      DO 200 IMOL=1,NMOLI
        IDSC=0
        LGMCX(IMOL,0,0)=0
        LGMCX(IMOL,0,1)=0
C
C  THERE ARE CURRENTLY NO DEFAULT CX RATES FOR MOLECULES
C
        IF (NRCM(IMOL).EQ.0) THEN
          NMCXI(IMOL)=0
C
C  NON-DEFAULT CX MODEL:
        ELSEIF (NRCM(IMOL).GT.0) THEN
          DO 130 NRC=1,NRCM(IMOL)
            KK=IREACM(IMOL,NRC)
            IF (ISWR(KK).NE.3) CYCLE
C  make sure that incident particle is a bulk particle
            IF (EIRENE_IDEZ(IBULKM(IMOL,NRC),1,3).NE.4) THEN
C  WRONG TYPE OF INCIDENT BULK SPECIES
              WRITE (IUNOUT,*)
     .        'INPUT ERROR FOR CX PROCESS, IMOL,KK ',IMOL,KK
              CALL EIRENE_EXIT_OWN(1)
            ENDIF
C  CX PROCESS IDENTIFIED

            FACTKK=FREACM(IMOL,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
            CHRDF0=0.D0
C  BULK PARTICLE INDEX
            IPLS=EIRENE_IDEZ(IBULKM(IMOL,NRC),3,3)
            IDSC=IDSC+1
            NRCXI=NRCXI+1
            IRCX=NRCXI
            LGMCX(IMOL,IDSC,0)=IRCX
            LGMCX(IMOL,IDSC,1)=IPLS
c
            if (ngenm(imol).lt.0) then  !  in range -1,...-infinity
c  set cx fluid limit FDLM (critical Knudsen number Kn_c = mfp_cx/delta
c  delta: typical length (could be cell size, or gradient length...)
c  use the integer input flag ngenm (generation limit).
              MFL=-(ngenm(imol)+1)  !  now MFL in range 0 to +infinity
c  ngena=-10001 produces Kn_c=1.0. Larger abs(ngenm) --> smaller Kn_c
              FDLMCX(IRCX)=1.0E4/(MFL+eps30)
            endif

            IML=NSPA+IMOL
            IPL=IPLS
            RMASS=RMASSM(IMOL)
            IFRST=ISCD1M(IMOL,NRC)
            ISCND=ISCD2M(IMOL,NRC)
            ISCDE=ISCDEM(IMOL,NRC)
            IESTM=IESTMM(IMOL,NRC)
            EBULK=EBULKM(IMOL,NRC)
            CALL EIRENE_XSTCX(RMASS,IRCX,IML,IPL,
     .                 IFRST,ISCND,EBULK,CHRDF0,
     .                 ISCDE,IESTM,KK,FACTKK,PLS)
C
  130     CONTINUE
C
          NMCXI(IMOL)=IDSC
C  NO CX MODEL DEFINED
        ELSE
          NMCXI(IMOL)=0
        ENDIF
C
        NMCXIM(IMOL)=NMCXI(IMOL)-1
C
        LGMCX(IMOL,0,0)=0
        DO IMCX=1,NMCXI(IMOL)
          LGMCX(IMOL,0,0)=LGMCX(IMOL,0,0)+LGMCX(IMOL,IMCX,0)
        ENDDO
C
  200 CONTINUE
C
C   ELASTIC COLLISIONS
C
      DO 300 IMOL=1,NMOLI
        IDSC=0
        LGMEL(IMOL,0,0)=0
        LGMEL(IMOL,0,1)=0
C
C  DEFAULT EL MODEL: NOT AVAILABLE
C
        IF (NRCM(IMOL).EQ.0) THEN
          NMELI(IMOL)=0
C
C  NON-DEFAULT EL MODEL:  240--
C
        ELSEIF (NRCM(IMOL).GT.0) THEN
          DO 230 NRC=1,NRCM(IMOL)
            KK=IREACM(IMOL,NRC)
            IF (ISWR(KK).NE.5) CYCLE
C  make sure that incident particle is a bulk particle
            IF (EIRENE_IDEZ(IBULKM(IMOL,NRC),1,3).NE.4) THEN
C  WRONG TYPE OF INCIDENT BULK SPECIES
              WRITE (IUNOUT,*)
     .        'INPUT ERROR FOR EL PROCESS, IMOL,KK ',IMOL,KK
              CALL EIRENE_EXIT_OWN(1)
            ENDIF
C  EL PROCESS IDENTIFIED
C
            FACTKK=FREACM(IMOL,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
C  BULK PARTICLE INDEX
            IPLS=EIRENE_IDEZ(IBULKM(IMOL,NRC),3,3)
            IF (IPLS.LE.0.OR.IPLS.GT.NPLSI) GOTO 991
            IF (MASSP(KK).LE.0.OR.MASST(KK).LE.0) GOTO 992
            IDSC=IDSC+1
            NRELI=NRELI+1
            IREL=NRELI
            LGMEL(IMOL,IDSC,0)=IREL
            LGMEL(IMOL,IDSC,1)=IPLS
C
C  SPECIAL TREATMENT: BGK COLLISIONS AMONGST TEST PARTICLES
            IF (IBGKM(IMOL,NRC).NE.0) THEN
              IF (NPBGKM(IMOL).EQ.0) THEN
C  IMOL HAS NOT YET BEEN LABELLED AS BGK SPECIES.
C  DO THIS HERE: IMOL IS BGK SPECIES NO. IBGK_SP, AND HAS 3 ADDITIONAL BGK TALLIES IN UPTBGK
                NRBGI=NRBGI+3
                IBGK_SP=NRBGI/3
                NPBGKM(IMOL)=IBGK_SP
              ENDIF
              IF (NPBGKP(IPLS,1).EQ.0) THEN
                NPBGKP(IPLS,1)=NPBGKM(IMOL)
              ELSE
                GOTO 999
              ENDIF
C  SELF- OR CROSS-COLLISION?
              ITYPB=EIRENE_IDEZ(IBGKM(IMOL,NRC),1,3)
              ISPZB=EIRENE_IDEZ(IBGKM(IMOL,NRC),3,3)
              IF (ITYPB.NE.2.OR.ISPZB.NE.IMOL) THEN
C  CROSS-COLLISION !
                IF (NPBGKP(IPLS,2).EQ.0) THEN
                  NPBGKP(IPLS,2)=IBGKM(IMOL,NRC)
                ELSE
                  GOTO 999
                ENDIF
              ENDIF
            ENDIF
C  BGK COLLISION PARAMETERS DONE
C
            IML=NSPA+IMOL
            IPL=IPLS
            ISCDE=ISCDEM(IMOL,NRC)
            IESTM=IESTMM(IMOL,NRC)
            EBULK=EBULKM(IMOL,NRC)
            CALL EIRENE_XSTEL(IREL,IML,IPL,EBULK,
     .                        ISCDE,IESTM,
     .                        KK,FACTKK,PLS)
C
  230     CONTINUE
C
          NMELI(IMOL)=IDSC
        ENDIF
C
        NMELIM(IMOL)=NMELI(IMOL)-1
C
        LGMEL(IMOL,0,0)=0.
        DO 280 IMEL=1,NMELI(IMOL)
          LGMEL(IMOL,0,0)=LGMEL(IMOL,0,0)+LGMEL(IMOL,IMEL,0)
  280   CONTINUE
C
C
  300 CONTINUE
C
C   GENERAL HEAVY PARTICLE IMPACT COLLISIONS
C

      DO IMOL=1,NMOLI
        IDSC=0
        LGMPI(IMOL,0,0)=0
        LGMPI(IMOL,0,1)=0
C
C  NO DEFAULT MODEL
C
        IF (NRCM(IMOL).EQ.0) THEN
          NMPII(IMOL)=0
C
C  NON-DEFAULT ION IMPACT MODEL:  130--190
C
        ELSEIF (NRCM(IMOL).GT.0) THEN
          DO NRC=1,NRCM(IMOL)
            KK=IREACM(IMOL,NRC)
            IF (ISWR(KK).NE.4) CYCLE
C  make sure that incident particle is a bulk particle
            IF (EIRENE_IDEZ(IBULKM(IMOL,NRC),1,3).NE.4) THEN
C  WRONG TYPE OF INCIDENT BULK SPECIES
              WRITE (IUNOUT,*)
     .        'INPUT ERROR FOR PI PROCESS, IMOL,KK ',IMOL,KK
              CALL EIRENE_EXIT_OWN(1)
            ENDIF
C  PI PROCESS IDENTIFIED

            FACTKK=FREACM(IMOL,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
            IF (MASSP(KK).LE.0.OR.MASST(KK).LE.0) GOTO 992
C  BULK PARTICLE INDEX
            IPLS=EIRENE_IDEZ(IBULKM(IMOL,NRC),3,3)
            IF (IPLS.LE.0.OR.IPLS.GT.NPLSI) GOTO 990
            IDSC=IDSC+1
            NRPII=NRPII+1
            IF (NRPII.GT.NRPI) GOTO 998
            IRPI=NRPII
            LGMPI(IMOL,IDSC,0)=IRPI
            LGMPI(IMOL,IDSC,1)=IPLS

            IML=NSPA+IMOL
            IPL=IPLS
            RMASS=RMASSM(IMOL)
            IFRST=ISCD1M(IMOL,NRC)
            ISCND=ISCD2M(IMOL,NRC)
            ITHRD=ISCD3M(IMOL,NRC)
            IFRTH=ISCD4M(IMOL,NRC)
            ISCDE=ISCDEM(IMOL,NRC)
            IESTM=IESTMM(IMOL,NRC)
            EBULK=EBULKM(IMOL,NRC)
            EHEAVY=ESCD1M(IMOL,NRC)
            EELEC=EELECM(IMOL,NRC)
            CALL EIRENE_XSTPI (RMASS,IRPI,IML,IPL,
     .                  IFRST,ISCND,ITHRD,IFRTH,
     .                  EBULK,EHEAVY,EELEC,CHRDF0,ISCDE,IESTM,
     .                  KK,FACTKK,PLS)
          END DO
C
          NMPII(IMOL)=IDSC
C  NO MODEL DEFINED
        ELSE
          NMPII(IMOL)=0
        ENDIF
C
        NMPIIM(IMOL)=NMPII(IMOL)-1
C
        LGMPI(IMOL,0,0)=0
        DO IMPI=1,NMPII(IMOL)
          LGMPI(IMOL,0,0)=LGMPI(IMOL,0,0)+LGMPI(IMOL,IMPI,0)
        END DO

        DO IMPI=1,NMPII(IMOL)
          IRPI=LGMPI(IMOL,IMPI,0)
          CALL EIRENE_XSTPI_1(IRPI)
        END DO
      END DO
C
C
      DO 1000 IMOL=1,NMOLI
C
        IF (TRCAMD) THEN
          CALL EIRENE_MASBOX
     .    ('MOLECULAR SPECIES IMOL = '//TEXTS(NSPA+IMOL))
          CALL EIRENE_LEER(1)
C
          IF (LGMEI(IMOL,0).EQ.0) THEN

            CALL EIRENE_LEER(1)
            WRITE (iunout,*) 'NO ELECTRON IMPACT COLLISIONS '
            CALL EIRENE_LEER(1)
          ELSE
            DO 870 IMEI=1,NMEII(IMOL)
              IREI=LGMEI(IMOL,IMEI)
              CALL EIRENE_XSTEI_2(IREI)
  870       CONTINUE
          ENDIF
C
C
          CALL EIRENE_LEER(2)
          IF (LGMCX(IMOL,0,0).EQ.0) THEN
            CALL EIRENE_LEER(1)
            WRITE (iunout,*) 'NO CHARGE EXCHANGE WITH BULK IONS'
            CALL EIRENE_LEER(1)
          ELSE
            DO 890 IMCX=1,NMCXI(IMOL)
              IRCX=LGMCX(IMOL,IMCX,0)
              IPL =LGMCX(IMOL,IMCX,1)
              CALL EIRENE_XSTCX_2(IRCX,IPL)
  890       CONTINUE
          ENDIF
C
C
          CALL EIRENE_LEER(2)
          IF (LGMEL(IMOL,0,0).EQ.0) THEN
            CALL EIRENE_LEER(1)
            WRITE (iunout,*) 'NO ELASTIC COLLISIONS WITH BULK IONS'
            CALL EIRENE_LEER(1)
          ELSE
            DO 895 IMEL=1,NMELI(IMOL)
              IREL=LGMEL(IMOL,IMEL,0)
              IPL =LGMEL(IMOL,IMEL,1)
              CALL EIRENE_XSTEL_2(IREL,IPL)
  895       CONTINUE
          ENDIF
C
          CALL EIRENE_LEER(2)
          IF (LGMPI(IMOL,0,0).EQ.0) THEN
            CALL EIRENE_LEER(1)
            WRITE (iunout,*) 'NO GENERAL ION IMPACT COLLISIONS '
            CALL EIRENE_LEER(1)
          ELSE
            DO 885 IMPI=1,NMPII(IMOL)
              IRPI=LGMPI(IMOL,IMPI,0)
              IPL =LGMPI(IMOL,IMPI,1)
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
      WRITE (iunout,*) 'ERROR IN XSECTM: EXIT CALLED'
      WRITE (iunout,*) 'INVALID SPECIES INDEX FOR ION IMPACT COLLISION'
      CALL EIRENE_EXIT_OWN(1)
  991 CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTM: EXIT CALLED'
      WRITE (iunout,*) 'INVALID SPECIES INDEX FOR ELASTIC COLLISION'
      CALL EIRENE_EXIT_OWN(1)
  992 CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTM: EXIT CALLED'
      WRITE (iunout,*)
     .  'MASS NUMBERS OF INTERACTING PARTICLES INCONSISTENT'
      WRITE (iunout,*) 'KK,IMOL,IPLS ',KK,IMOL,IPLS
      CALL EIRENE_EXIT_OWN(1)
  994 CONTINUE
      WRITE (iunout,*) 'ERROR DETECTED IN XSECTM.'
      WRITE (iunout,*) 'REACTION NO. KK= ',KK, 'NOT READ FROM FILE'
      WRITE (iunout,*) 'IMOL = ',IMOL
      WRITE (iunout,*) 'ISWR(KK) = ',ISWR(KK)
      WRITE (iunout,*) 'EXIT CALLED'
      CALL EIRENE_EXIT_OWN(1)
  996 CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTM: EXIT CALLED'
      WRITE (iunout,*) 'NO COLLISION DATA AVAILABLE FOR THE CHOICE'
      WRITE (iunout,*) 'OF POST-COLLISION SAMPLING FLAG ISCDEA'
      WRITE (iunout,*) 'OR OTHER COLLISION DATA INCONSISTENCY'
      CALL EIRENE_EXIT_OWN(1)
  998 CONTINUE
      WRITE (iunout,*) 'INSUFFICIENT STORAGE FOR PI: NRPI=',NRPI
      CALL EIRENE_EXIT_OWN(1)
  999 CONTINUE
      WRITE (iunout,*) 'SPECIES CONFLICT FOR BGK COLLISIONS. IMOL,IREL '
      WRITE (iunout,*) IMOL,IREL,IPLS
      CALL EIRENE_EXIT_OWN(1)
      RETURN
C
      END
