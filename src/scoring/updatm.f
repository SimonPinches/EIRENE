C  27.6.05 updphot: iadd removed
C  21.01.06: photon background for test atoms: removed
C  18.04.06: test ions and atoms: syncronized
C            bug fix: V0_para  --> parmom_0 for elastic momentum source
C                                  contribution from atoms.
C  10.01.07: parallel momentum exchange tallies MAPL, MMPL, MIPL
C            included as default EIRENE tallies.
C            Before these tallies have been updated in problem
C            specific section UPTCOP, as COPV tallies.
C  12.02.07: Add user supplied B field and plasma flow option indpro=8
C            to evaluation of parallel momentum sources,
C            Do not use BVIN, PARMOM arrays in this case, because they
C            may not have been initialized in subr. PLASMA_DERIV
C            for these options.
C            Use vsig_parp und val_parp instead.
C  25.04.07 update of tallies because of PI reactions revised
C  07.08.07 collision estimators vollstaendig fuer atom, mol und iion.
C           entries: atm, mol, ion voll syncronisiert.
C  28.8.07: esigpi(...,4) --> PL, esigpi(...,5)--> EL
c  oct.14:  some intermediate scoring of additional tally ADDV removed, back to development branch 
c  06.08.15 arguments added to vecusr
c  24.08.15 comments and documention wrt. BGK collision treatment
cdr dec.15: tracklength estimators for heavy test particle post collision energies 
cdr         in PI processes added. For A, M, I incident test particles.
cdr dec.15: further corrections, lea --> leio, and other logical flags for turning on-off estimators

cdr nov.15: tracklength estimators for eapl,empl,eipl: species ipl resolved.
cdr apr.16: bug fix J.Lore re index in lgiel. This part of code is still unused,
cdr          so no effect on any result.  Few further comments corrected

!pb APR  16: ipplds -> ipplei, pplds -> pplei
!pb APR  16: ipatds -> ipatei, patds -> patei
!pb APR  16: ipmlds -> ipmlei, pmlds -> pmlei
!pb APR  16: ipiods -> ipioei, piods -> pioei
!pb APR  16: pelds -> pelei
cdr sept 16: nmdsi -> nmeii, nidsi -> nieii
cdr dec. 16: some more comments re sign convention for momentum sources 

 
C
      SUBROUTINE EIRENE_UPDATM (XSTOR2,XSTORV2,IFLAG)
C
C ESTIMATORS ARE UPDATED FOR EACH TRACK TAKING T/VEL SEC.
C T (CM) IS STORED ON CLPD ARRAY FOR ONE OR MORE CELLS, THAT HAVE
C BEEN CROSSED WITHOUT COLLISION.

C
C  NCOU:  NUMBER OF PIECES OF TRACK IN DIFFERENT CELLS SCORED IN THIS PRESENT CALL (BUT FIXED NRCELL)
C     I:  INDIVIDUAL TRACK, I=1,NCOU
C  IRDO:  TRACK IS IN (FINE) GEOMETRY CELL IRDO (=NRCELL+NUPC(I)*NR1P2+NBLCKA)
C  IRD:   ESTIMATORS ARE UPDATED IN (COARSE) SCORING CELL IRD  (=NCLTAL(IRDO))
C
C  IFLAG:  CURRENTLY ONLY USED FOR PHOTON TALLIES, TO AVOID CANCELLATION OF TERMS

C  IFLAG=1:  
C  IFLAG=2:  
C  IFLAG=3:  
C  IFLAG=4:  CALLED FROM WITHIN STATIC LOOP  (PATH LENGTH SET TO MFP), OR CALLED AT POINT OF COLLISION
C  IFLAG=5: 

C  SPECIAL TREATMENT OF "BGK" COLLISIONS (= ELASTIC COLLISIONS WITH VIRTUEL BACKGROUND SPECIES) 
C
C  A) NPBGK..(ITEST) :  IF GT 0, THE CORRESPONDING PARTICLE (IATM, IMOL OR IION) IS A SO CALLED "BGK" SPECIES
C                             IF, ADDITIONALLY, LBGKV = T, THEN ADDITIONAL BGK TALLIES ARE SCORED VIA A CALL TO UPTBGK
C  B) SIGBGK         : TOTAL RATE OF BGK TYPE COLLISIONS. INCIDENT ATOM AND ITS ENERGY IS NOT LOST
C  C) NPBGKP (IPLS,1):  IREL ELASTIC COLLISION CONTRIBUTIONS WITH BULK COLLISION PARTNERS WITH NPBGKP(IPLS,1)>0
C          ARE NOT INCLUDED IN SOURCE/SINK TALLIES.

C          IN CASE OF EAPL THIS IS IMPORTANT, IN ORDER NOT TO MIX ENERGY SOURCES FOR REAL BACKGROUND
C          IONS WITH ENERGY SOURCES FOR VIRTUEL BACKGROUND "IONS"  (MISSING SPECIES INDEX) 
C          BUT:  CURRENTLY MISSING IN EAAT: CONTRIBUTIONS OF ENERGY EXCHANGE DUE TO BGK COLLISIONS 
C          (BOTH SOURCE (DUE TO C) AND SINK (DUE TO B)
  

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_CUPD
      USE EIRMOD_CGRID
      USE EIRMOD_CSPEZ
      USE EIRMOD_CGEOM
      USE EIRMOD_COMPRT
      USE EIRMOD_CSDVI
      USE EIRMOD_COMXS
      USE EIRMOD_CZT1
      USE EIRMOD_CCONA
      USE EIRMOD_PHOTON
      USE EIRMOD_CINIT
 
      IMPLICIT NONE
C
      REAL(DP), INTENT(IN OUT) :: XSTOR2(MSTOR1,MSTOR2,N2ND+N3RD),
     .                            XSTORV2(NSTORV,N2ND+N3RD)
      INTEGER, INTENT(IN) :: IFLAG
      REAL(DP) :: WTRSIG, DIST, WTR, WTRE0, WV, VELQ, CNDYNPH, WTRV,
     .            V0_PARB, PARMOM_0, P, BX, BY, BZ, BF, VION
      REAL(DP) :: VSIG_PARB(NPLS), VAL_PARB(NPLS), VX(NPLS), VY(NPLS),
     .            VZ(NPLS),XC,YC,ZC
      INTEGER :: IRD,  I, IRDO, INUM,
     .           IPL, IAT, IA,
     .           IM,  IIO, IP, IML, II, NPBGK,
     .           IBGK, IPLV
C SECONDARY SPECIES IDENTIFIERS
      INTEGER ::  IAT1,IAT2,IML1,IML2,IIO1,IIO2,IPH1,IPH2,IPL1,IPL2
C EL PROCESSES
      INTEGER ::      IAEL,IREL
      INTEGER ::      IMEL
      INTEGER ::      IIEL
C CX PROCESSES
      INTEGER ::      IACX,IRCX
      INTEGER ::      IMCX
      INTEGER ::      IICX
C OT PROCESSES
      INTEGER ::      IAOT,IROT,UPDF
C PI PROCESSES
      INTEGER ::      IAPI,IRPI
      INTEGER ::      IMPI
      INTEGER ::      IIPI
C EI PROCESSES
      INTEGER ::      IAEI,IREI
      INTEGER ::      IMEI
      INTEGER ::      IIEI
 
      REAL(DP) :: EIRENE_VDION


C
C  ESTIMATORS FOR ATOMS
C
C
      WV=WEIGHT/VEL
      NPBGK=NPBGKA(IATM)
C
      IF (NADVI.GT.0) CALL EIRENE_UPTUSR(XSTOR2,XSTORV2,WV,IFLAG)
      IF (NCPVI.GT.0) CALL EIRENE_UPTCOP(XSTOR2,XSTORV2,WV,IFLAG)
      IF ((NPBGK.GT.0).AND.LBGKV)
     .   CALL EIRENE_UPTBGK(XSTOR2,XSTORV2,WV,NPBGK,IFLAG)
 
      IF (IUPDTE == 2) RETURN
 
C
      VELQ=VEL*VEL
C
      DO 51 I=1,NCOU
        DIST=CLPD(I)
        WTR=WV*DIST
        WTRE0=WTR*E0
        WTRV=WTR*VEL*CNDYNA(IATM)
        IRDO=NRCELL+NUPC(I)*NR1P2+NBLCKA
        IRD=NCLTAL(IRDO)

C  FOR STANDARD DEVIATION: INDICATE CELLS THAT HAVE BEEN MET BY THE PRESENT MC HISTORY
        IF (IMETCL(IRD) == 0) THEN
          NCLMT = NCLMT+1
          ICLMT(NCLMT) = IRD
          IMETCL(IRD) = NCLMT
        END IF
C
C  PARTICLE, MOMENTUM AND ENERGY DENSITY ESTIMATORS
C
        IF (LEDENA) EDENA(IATM,IRD)=EDENA(IATM,IRD)+WTRE0
        IF (LPDENA) PDENA(IATM,IRD)=PDENA(IATM,IRD)+WTR
        IF (LEDENA.OR.LPDENA) LMETSP(NSPH+IATM)=.TRUE.
 
        IF (LVXDENA) VXDENA(IATM,IRD)=VXDENA(IATM,IRD)+WTRV*VELX
        IF (LVYDENA) VYDENA(IATM,IRD)=VYDENA(IATM,IRD)+WTRV*VELY
        IF (LVZDENA) VZDENA(IATM,IRD)=VZDENA(IATM,IRD)+WTRV*VELZ
        IF (LVXDENA.OR.LVYDENA.OR.LVZDENA) LMETSP(NSPH+IATM)=.TRUE.
C
C  ESTIMATORS FOR SOURCES AND SINKS
C  NEGATIVE SIGN MEANS: LOSS FOR PARTICLES
C  POSITIVE SIGN MEANS: GAIN FOR PARTICLES
C
        IF (LGVAC(IRDO,0)) GOTO 51
C
        if (ncou.gt.1) then
          XSTOR(:,:) = XSTOR2(:,:,I)
          XSTORV(:)  = XSTORV2(:,I)
        endif
C
C  PRE COLLISION RATES, ASSUME: TEST PARTICLES (AND THEIR ENERGY) ARE LOST
C
        WTRSIG=WTR*(SIGTOT-SIGBGK)
        IF (LPAAT) PAAT(IATM,IRD)=PAAT(IATM,IRD)-WTRSIG
        IF (LEAAT) EAAT(IRD)     =EAAT(IRD)     -WTRSIG*E0
C
C..........................................................................
C
C  CHARGE EXCHANGE CONTRIBUTION
C
        IF (LGACX(IATM,0,0).EQ.0) GOTO 43
C  DEFAULT TRACKLENGTH ESTIMATOR
        DO 44  IACX=1,NACXI(IATM)
          IRCX=LGACX(IATM,IACX,0)
          IPLS=LGACX(IATM,IACX,1)
          LOGPLS(IPLS,ISTRA)=.TRUE.
C
          WTRSIG=WTR*SIGVCX(IRCX)
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
C  COMPENSATE PRE COLLISION RATES HERE
C
          IF (IESTCX(IRCX,1).NE.0) THEN
            IF (LPAAT) PAAT(IATM,IRD)=PAAT(IATM,IRD)+WTRSIG
          ELSE
C
C  PRE COLLISION RATES, BULK IONS
C
            IF (LPAPL) THEN
              PAPL(IPLS,IRD)=PAPL(IPLS,IRD)-WTRSIG
              LMETSP(NSPAMI+IPLS)=.TRUE.
            END IF
C
C  POST COLLISION RATES, ALL SECONDARIES (TEST AND BULK PARTICLES)
C  FIRST SECONDARY: PREVIOUS BULK ION IPL
            IF (N1STX(IRCX,1).EQ.1) THEN
              IAT1=N1STX(IRCX,2)
              LOGATM(IAT1,ISTRA)=.TRUE.
              IF (LPAAT) THEN
                PAAT(IAT1,IRD)= PAAT(IAT1,IRD)+WTRSIG
                LMETSP(NSPH+IAT1)=.TRUE.
              END IF
            ELSEIF (N1STX(IRCX,1).EQ.2) THEN
              IML1=N1STX(IRCX,2)
              LOGMOL(IML1,ISTRA)=.TRUE.
              IF (LPAML) THEN
                PAML(IML1,IRD)= PAML(IML1,IRD)+WTRSIG
                LMETSP(NSPA+IML1)=.TRUE.
              END IF
            ELSEIF (N1STX(IRCX,1).EQ.3) THEN
              IIO1=N1STX(IRCX,2)
              LOGION(IIO1,ISTRA)=.TRUE.
              IF (LPAIO) THEN
                PAIO(IIO1,IRD)= PAIO(IIO1,IRD)+WTRSIG
                LMETSP(NSPAM+IIO1)=.TRUE.
              END IF
            ELSEIF (N1STX(IRCX,1).EQ.4) THEN
              IPL1=N1STX(IRCX,2)
              LOGPLS(IPL1,ISTRA)=.TRUE.
              IF (LPAPL) THEN
                PAPL(IPL1,IRD)= PAPL(IPL1,IRD)+WTRSIG
                LMETSP(NSPAMI+IPL1)=.TRUE.
              END IF
            ENDIF
C  SECOND SECONDARY: PREVIOUS ATOM IATM
            IF (N2NDX(IRCX,1).EQ.1) THEN
              IAT2=N2NDX(IRCX,2)
              LOGATM(IAT2,ISTRA)=.TRUE.
              IF (LPAAT) THEN
                PAAT(IAT2,IRD)= PAAT(IAT2,IRD)+WTRSIG
                LMETSP(NSPH+IAT2)=.TRUE.
              END IF
            ELSEIF (N2NDX(IRCX,1).EQ.2) THEN
              IML2=N2NDX(IRCX,2)
              LOGMOL(IML2,ISTRA)=.TRUE.
              IF (LPAML) THEN
                PAML(IML2,IRD)= PAML(IML2,IRD)+WTRSIG
                LMETSP(NSPA+IML2)=.TRUE.
              END IF
            ELSEIF (N2NDX(IRCX,1).EQ.3) THEN
              IIO2=N2NDX(IRCX,2)
              LOGION(IIO2,ISTRA)=.TRUE.
              IF (LPAIO) THEN
                PAIO(IIO2,IRD)= PAIO(IIO2,IRD)+WTRSIG
                LMETSP(NSPAM+IIO2)=.TRUE.
              END IF
            ELSEIF (N2NDX(IRCX,1).EQ.4) THEN
              IPL2=N2NDX(IRCX,2)
              LOGPLS(IPL2,ISTRA)=.TRUE.
              IF (LPAPL) THEN
                PAPL(IPL2,IRD)= PAPL(IPL2,IRD)+WTRSIG
                LMETSP(NSPAMI+IPL2)=.TRUE.
              END IF
            ENDIF
          ENDIF
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
C  COMPENSATE PRE COLLISION RATES HERE
C
          IF (LEA) THEN
            IF (IESTCX(IRCX,3).NE.0) THEN
              IF (LEAAT) EAAT(IRD) = EAAT(IRD) + WTRSIG*E0
            ELSE
C
C  PRE COLLISION RATES, BULK IONS
C
              IF (LEAPL) THEN
                EAPL(IPLS,IRD) = EAPL(IPLS,IRD) - WTRSIG*ESIGCX(IRCX,1)
                LMETSP(NSPAMI+IPLS)=.TRUE.
              END IF
C
C  POST COLLISION RATES, ALL SECONDARIES (TEST AND BULK PARTICLES)
C  FIRST SECONDARY: PREVIOUS BULK ION IPL
              IF (N1STX(IRCX,1).EQ.1) THEN
                IAT1=N1STX(IRCX,2)
                LOGATM(IAT1,ISTRA)=.TRUE.
                IF (LEAAT) EAAT(IRD) = EAAT(IRD) + WTRSIG*ESIGCX(IRCX,1)
              ELSEIF (N1STX(IRCX,1).EQ.2) THEN
                IML1=N1STX(IRCX,2)
                LOGMOL(IML1,ISTRA)=.TRUE.
                IF (LEAML) EAML(IRD) = EAML(IRD) + WTRSIG*ESIGCX(IRCX,1)
              ELSEIF (N1STX(IRCX,1).EQ.3) THEN
                IIO1=N1STX(IRCX,2)
                LOGION(IIO1,ISTRA)=.TRUE.
                IF (LEAIO) EAIO(IRD) = EAIO(IRD) + WTRSIG*ESIGCX(IRCX,1)
              ELSEIF (N1STX(IRCX,1).EQ.4) THEN
                IPL1=N1STX(IRCX,2)
                LOGPLS(IPL1,ISTRA)=.TRUE.
                IF (LEAPL) THEN
                  EAPL(IPL1,IRD) = EAPL(IPL1,IRD)+ WTRSIG*ESIGCX(IRCX,1)
                  LMETSP(NSPAMI+IPL1) = .TRUE.
                END IF
              ENDIF
C  SECOND SECONDARY: PREVIOUS ATOM IATM
              IF (N2NDX(IRCX,1).EQ.1) THEN
                IAT2=N2NDX(IRCX,2)
                LOGATM(IAT2,ISTRA)=.TRUE.
                IF (LEAAT) EAAT(IRD) = EAAT(IRD) + WTRSIG*E0
              ELSEIF (N2NDX(IRCX,1).EQ.2) THEN
                IML2=N2NDX(IRCX,2)
                LOGMOL(IML2,ISTRA)=.TRUE.
                IF (LEAML) EAML(IRD) = EAML(IRD) + WTRSIG*E0
              ELSEIF (N2NDX(IRCX,1).EQ.3) THEN
                IIO2=N2NDX(IRCX,2)
                LOGION(IIO2,ISTRA)=.TRUE.
                IF (LEAIO) EAIO(IRD) = EAIO(IRD) + WTRSIG*E0
              ELSEIF (N2NDX(IRCX,1).EQ.4) THEN
                IPL2=N2NDX(IRCX,2)
                LOGPLS(IPL2,ISTRA)=.TRUE.
                IF (LEAPL) THEN 
                  EAPL(IPL2,IRD) = EAPL(IPL2,IRD) + WTRSIG*E0
                  LMETSP(NSPAMI+IPL2) = .TRUE.
                END IF
              ENDIF
            ENDIF
          ENDIF
C
44      CONTINUE
43      CONTINUE
C
C  ELASTIC NEUTRAL BULK-ION COLLISION CONTRIBUTION
C
        IF (LGAEL(IATM,0,0).EQ.0) GOTO 60
C  DEFAULT TRACKLENGTH ESTIMATOR
        DO 61  IAEL=1,NAELI(IATM)
          IREL=LGAEL(IATM,IAEL,0)
          IPLS=LGAEL(IATM,IAEL,1)
C  DO NOT UPDATE BGK SOURCE RATE TALLIES HERE
          IBGK=NPBGKP(IPLS,1)
C  ELASTIC REACTION IAEL, BETWEEN SPECIES IATM/IPLS: 
C  IPLS IS A BGK VIRTUAL BACKGROUND SPECIES. 
C  CONTRIBUTIONS TO PAAT, PAPL ARE IDENTICALLY ZERO 
C  EAAT AND EAPL SHOULD NOT BE UPDATED HERE, BECAUSE THEY ARE SUMMED OVER SPECIES.
          IF (IBGK.NE.0) GOTO 61

          LOGPLS(IPLS,ISTRA)=.TRUE.
C
          WTRSIG=WTR*SIGVEL(IREL)
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
C  COMPENSATE PRE COLLISION RATES HERE
C
          IF (IESTEL(IREL,1).NE.0) THEN
            IF (LPAAT) PAAT(IATM,IRD)=PAAT(IATM,IRD)+WTRSIG
          ELSE
C  UPDATE TRACKLENGTH ESTIMATOR
C           IF (LPAPL) THEN
C             PAPL(IPLS,IRD)=PAPL(IPLS,IRD)-WTRSIG
C             PAPL(IPLS,IRD)=PAPL(IPLS,IRD)+WTRSIG
C             LMETSP(NSPAMI+IPLS)=.TRUE.
C           END IF
            IF (LPAAT) THEN
              PAAT(IATM,IRD)=PAAT(IATM,IRD)+WTRSIG
              LMETSP(NSPH+IATM)=.TRUE.
            END IF
          ENDIF
C
          IF (LEA) THEN
            IF (IESTEL(IREL,3).NE.0) THEN
 
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
C  COMPENSATE SUBTRACTED PRE COLLISION RATES HERE
 
              IF (LEAAT) EAAT(IRD)=EAAT(IRD)+WTRSIG*E0
            ELSE
C
C  DEFAULT TRACKLENGTH ESTIMATOR ("PERFECT IDENTITY EXCHANGE" APPROXIMATION)
C     THIS APPROXIMATION CORRESPONDS TO VELOCITY INDEPENDENT COLLISION RATES, 
C     SUCH AS E.G. THOSE USED FOR BGK APPROXIMATIONS TO NEUTRAL-NEUTRAL COLLISIONS
C  AVERAGE ENERGY OF POST COLLISION ATOM IS THAT OF PRE COLLISION BULK
C
C  PRE COLLISION RATES, BULK IONS
C
              IF (LEAPL) THEN
                EAPL(IPLS,IRD)=EAPL(IPLS,IRD)-WTRSIG*ESIGEL(IREL,1)
                LMETSP(NSPAMI+IPLS)=.TRUE.
              END IF
C
C  FIRST SECONDARY: = INCIDENT BULK ION. REMAINS SAME PARTICLE BY DEFAULT
              IF (LEAPL)  THEN
                EAPL(IPLS,IRD)=EAPL(IPLS,IRD)+WTRSIG*E0
                LMETSP(NSPAMI+IPLS)=.TRUE.
              END IF
C  SECOND SECONDARY: = INCIDENT ATOM. REMAINS SAME PARTICLE BY DEFAULT
              IF (LEAAT) EAAT(IRD)=EAAT(IRD)+WTRSIG*ESIGEL(IREL,1)
            ENDIF
          ENDIF
C
61      CONTINUE
60      CONTINUE
C
C.............................................................
C  ELECTRON IMPACT COLLISION CONTRIBUTION:  EL + IATM --> ....
C.............................................................
C
        IF (LGAEI(IATM,0).EQ.0) GOTO 57
C
        DO 55 IAEI=1,NAEII(IATM)
          IREI=LGAEI(IATM,IAEI)
          IF (SIGVEI(IREI).LE.0.D0) GOTO 55
C
          WTRSIG=WTR*SIGVEI(IREI)
C
C  EI PROCESS NO. IREI:
C
C  COLLISION ESTIMATOR FOR PARTICLE BALANCE IN SUBR. COLLIDE ?
C  COMPENSATE PRE COLLISION RATES HERE
C
          IF (IESTEI(IREI,1).NE.0) THEN
            IF (LPAAT) PAAT(IATM,IRD)=PAAT(IATM,IRD)+WTRSIG
          ELSE
C
C  TRACKLENGTH ESTIMATOR FOR PARTICLE BALANCE
C
C  ELECTRONS: DO NOT SEPARATE PRE AND POST COLLISION. UPDATE NET RATES
C
            IF (LPAEL) PAEL(IRD)=PAEL(IRD)+WTRSIG*PELEI(IREI)
C
C  POST COLLISION CONTRIBUTIONS
            DO IA=1,IPATEI(IREI,0)
              IAT=IPATEI(IREI,IA)
              LOGATM(IAT,ISTRA)=.TRUE.
              IF (LPAAT) THEN
                PAAT(IAT,IRD)=PAAT(IAT,IRD)+PATEI(IREI,IAT)*WTRSIG
                LMETSP(NSPH+IAT)=.TRUE.
              END IF
            END DO
 
            DO IM=1,IPMLEI(IREI,0)
              IML=IPMLEI(IREI,IM)
              LOGMOL(IML,ISTRA)=.TRUE.
              IF (LPAML) THEN
                PAML(IML,IRD)=PAML(IML,IRD)+PMLEI(IREI,IML)*WTRSIG
                LMETSP(NSPA+IML)=.TRUE.
              END IF
            END DO
 
            DO II=1,IPIOEI(IREI,0)
              IIO=IPIOEI(IREI,II)
              LOGION(IIO,ISTRA)=.TRUE.
              IF (LPAIO) THEN
                PAIO(IIO,IRD)=PAIO(IIO,IRD)+PIOEI(IREI,IIO)*WTRSIG
                LMETSP(NSPAM+IIO)=.TRUE.
              END IF
            END DO
 
            DO IP=1,IPPLEI(IREI,0)
              IPL=IPPLEI(IREI,IP)
              LOGPLS(IPL,ISTRA)=.TRUE.
              IF (LPAPL) THEN
                PAPL(IPL,IRD)=PAPL(IPL,IRD)+PPLEI(IREI,IPL)*WTRSIG
                LMETSP(NSPAMI+IPL)=.TRUE.
              END IF
            END DO
 
          ENDIF
C
C  PARTICLE BALANCE ESTIMATORS DONE.
C  NOW DEAL WITH ENERGY BALANCE ESTIMATORS
C  (STILL: EI PROCESSES)
C
          IF (IESTEI(IREI,3).EQ.0) THEN
            IF (LEAEL) EAEL(IRD)=EAEL(IRD)+WTRSIG*ESIGEI(IREI,5)
          ENDIF
 
          IF (LEA) THEN
            IF (IESTEI(IREI,3).NE.0) THEN
C
C  COLLISION ESTIMATOR
C  COMPENSATE PRE COLLISION CONTRIBUTION
C
              IF (LEAAT) EAAT(IRD)=EAAT(IRD)+WTRSIG*E0
C
            ELSE
C
              IF (LEAAT) EAAT(IRD)=EAAT(IRD)+WTRSIG*ESIGEI(IREI,1)
              IF (LEAML) EAML(IRD)=EAML(IRD)+WTRSIG*ESIGEI(IREI,2)
              IF (LEAIO) EAIO(IRD)=EAIO(IRD)+WTRSIG*ESIGEI(IREI,3)
              IF (LEAPL) THEN
                DO IP=1,IPPLEI(IREI,0)
                  IPL=IPPLEI(IREI,IP)
                  LOGPLS(IPL,ISTRA)=.TRUE.
cdr  this is incorrect. esigei is sum over ipl species.
cdr  it only happens to be correct if the post collision bulk species are the same (ipl),
cdr  because then esigei is the total for this species.
cdr  Must be fragmented into individual ipl contributions
                  EAPL(IPL,IRD)=EAPL(IPL,IRD)+WTRSIG*ESIGEI(IREI,4)
                  LMETSP(NSPAMI+IPL)=.TRUE.
                END DO
              END IF 
C
            ENDIF
          ENDIF
55      CONTINUE
 
57      CONTINUE
C
C........................................................
C  PLASMA ION IMPACT CONTRIBUTION: IPLS + IATM --> ......
C........................................................
C
        IF (LGAPI(IATM,0,0).EQ.0) GOTO 59
 
        DO 58  IAPI=1,NAPII(IATM)
          IRPI=LGAPI(IATM,IAPI,0)
          IPLS=LGAPI(IATM,IAPI,1)
          LOGPLS(IPLS,ISTRA)=.TRUE.
 
          WTRSIG=WTR*SIGVPI(IRPI)
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
C  COMPENSATE PRE COLLISION RATES HERE
C
          IF (IESTPI(IRPI,1).NE.0) THEN
            IF (LPAAT) PAAT(IATM,IRD)=PAAT(IATM,IRD)+WTRSIG
          ELSE
C
C  TRACKLENGTH ESTIMATOR FOR PARTICLE BALANCE
C
C
C  PRE COLLISION BULK ION CONTRIBUTION, ASSUME: INCIDENT ION IS LOST
C
            IF (LPAPL) THEN
              PAPL(IPLS,IRD)=PAPL(IPLS,IRD)-WTRSIG
              LMETSP(NSPAMI+IPLS)=.TRUE.
            END IF
C
C  ELECTRONS: HERE: ONLY POST COLLISION CONTRIBUTIONS
C
            IF (LPAEL) PAEL(IRD)=PAEL(IRD)+WTRSIG*PELPI(IRPI)
C
            DO IA=1,IPATPI(IRPI,0)
              IAT=IPATPI(IRPI,IA)
              LOGATM(IAT,ISTRA)=.TRUE.
              IF (LPAAT) THEN
                PAAT(IAT,IRD)= PAAT(IAT,IRD)+WTRSIG*PATPI(IRPI,IAT)
                LMETSP(NSPH+IAT)=.TRUE.
              END IF
            ENDDO
C
            DO IM=1,IPMLPI(IRPI,0)
              IML=IPMLPI(IRPI,IM)
              LOGMOL(IML,ISTRA)=.TRUE.
              IF (LPAML) THEN
                PAML(IML,IRD)= PAML(IML,IRD)+WTRSIG*PMLPI(IRPI,IML)
                LMETSP(NSPA+IML)=.TRUE.
              END IF
            ENDDO
C
            DO II=1,IPIOPI(IRPI,0)
              IIO=IPIOPI(IRPI,II)
              LOGION(IIO,ISTRA)=.TRUE.
              IF (LPAIO) THEN
                PAIO(IIO,IRD)= PAIO(IIO,IRD)+WTRSIG*PIOPI(IRPI,IIO)
                LMETSP(NSPAM+IIO)=.TRUE.
              END IF
            ENDDO
 
            DO IP=1,IPPLPI(IRPI,0)
              IPL=IPPLPI(IRPI,IP)
              LOGPLS(IPL,ISTRA)=.TRUE.
              IF (LPAPL) THEN
                PAPL(IPL,IRD)= PAPL(IPL,IRD)+WTRSIG*PPLPI(IRPI,IPL)
                LMETSP(NSPAMI+IPL)=.TRUE.
              END IF
            ENDDO
 
          ENDIF
C
C  PARTICLE BALANCE ESTIMATORS DONE.
C  NOW DEAL WITH ENERGY BALANCE ESTIMATORS
C  (STILL: PI PROCESSES)
C
          IF (IESTPI(IRPI,3).EQ.0) THEN
            IF (LEAEL) EAEL(IRD)=EAEL(IRD)+WTRSIG*ESIGPI(IRPI,5)
          ENDIF
 
          IF (LEA) THEN
            IF (IESTPI(IRPI,3).NE.0) THEN
C
C  COLLISION ESTIMATOR
C  COMPENSATE PRE COLLISION CONTRIBUTION
C
              IF (LEAAT) EAAT(IRD)=EAAT(IRD)+WTRSIG*E0
C
            ELSE
C
              IF (LEAAT) EAAT(IRD)=EAAT(IRD)+WTRSIG*ESIGPI(IRPI,1)
              IF (LEAML) EAML(IRD)=EAML(IRD)+WTRSIG*ESIGPI(IRPI,2)
              IF (LEAIO) EAIO(IRD)=EAIO(IRD)+WTRSIG*ESIGPI(IRPI,3)
              IF (LEAPL) THEN
                DO IP=1,IPPLPI(IRPI,0)
                  IPL=IPPLPI(IRPI,IP)
                  LOGPLS(IPL,ISTRA)=.TRUE.
cdr  this is incorrect. esigpi is sum over ipl species.
cdr  it only happens to be correct if the post collision bulk species are all the same (ipl),
cdr  because then esigpi is the total for this species.
cdr  Must be fragmented into individual ipl contributions
                  EAPL(IPL,IRD)=EAPL(IPL,IRD)+WTRSIG*ESIGPI(IRPI,4)
                  LMETSP(NSPAMI+IPL)=.TRUE.
                ENDDO
              END IF
            ENDIF
          ENDIF
58      CONTINUE
 
59      CONTINUE
C
C.........................................................................
C
C   PARALLEL MOMENTUM EXCHANGE RATE: DYN/CM**3,  CONTRIBUTIONS FROM ATOMS
C
C   CONTRIBUTIONS FROM CX, EI, PI, EL
C   PI: TO BE WRITTEN
C
C.........................................................................
C
C
C  PROJECTIONS, FIND PARALLEL COMPONENTS
C
C
        IF (LMAPL) THEN
 
          CALL EIRENE_BFIELD (IRDO, X0, Y0, Z0, BX, BY, BZ, BF,.FALSE.)

          DO IPL=1,NPLSI
            IF (INDPRO(4) == 8) THEN
              XC=0.
              YC=0.
              ZC=0.
              CALL EIRENE_VECUSR (2,IRDO,XC,YC,ZC,
     .             VX(IPL),VY(IPL),VZ(IPL),IPL,.FALSE.)
            ELSE
              IPLV=MPLSV(IPL)
              VX(IPL)=VXIN(IPLV,IRDO)
              VY(IPL)=VYIN(IPLV,IRDO)
              VZ(IPL)=VZIN(IPLV,IRDO)
            END IF
          END DO

c  set parameters for parallel momentum of incident bulk particle
c  val_parb   : parallel velocity component, incl. sign, relavive to B
c  vsig_parb  : parallel momentum, modulus (always positive)  
          IF ((INDPRO(4) == 8) .AND. (INDPRO(5) == 8)) THEN
            vion=EIRENE_vdion(irdo)
            VAL_PARB(1:NPLSI) =VION
            VSIG_PARB(1:NPLSI)=CNDYNP(1:NPLSI)*vion*
     .                         SIGN(1._DP,VION)
            
          ELSE IF ((INDPRO(4) == 8) .OR. (INDPRO(5) == 8)) THEN
C  PARMOM AND BVIN NOT KNOWN FROM PLASMA_DERIV
            DO IPL=1,NPLSI
              VAL_PARB(IPL) =(VX(IPL)*BX+VY(IPL)*BY+VZ(IPL)*BZ)
              VSIG_PARB(IPL)=CNDYNP(IPL)*VAL_PARB(IPL)*
     .                        SIGN(1._DP,VAL_PARB(IPL))
            END DO
          ELSE
            VAL_PARB(1:NPLSI) =BVIN(MPLSV(1:NPLSI),IRDO)
            VSIG_PARB(1:NPLSI)=PARMOM(1:NPLSI,IRDO)         
          END IF
c
c  set parameters for parallel momentum of incident neutral particle
c  v0_parb   : parallel velocity component, incl. sign, relavive to B
c  parmom_0  : parallel momentum   
          V0_PARB=VEL*(VELX*BX+VELY*BY+VELZ*BZ)
          PARMOM_0=CNDYNA(IATM)*V0_PARB
C                       *SIGN(1._DP,VAL_PARB(IPL))  !this sign factor is applied below
C  WITH THIS FACTOR: NO MATTER HOW THE SIGN OF PARALLEL MOMENTUM IS DEFINED:
C     THE PLASMA MOMENTUM IS TAKEN POSITIVE (PARMOM=|PARMOM|), AND  
C     |PARMOM_0| IS ADDED TO IPL MOMENTUM (SOURCE),  IF THE NEUTRAL V_PAR
C                     HAS THE SAME SIGN AS THE IPL PLASMA ION V_PAR.
c     |PARMOM_0| IS SUBTRACTED IF IT HAS OPPOSITE SIGN 
        

C  CHARGE EXCHANGE CONTRIBUTION FROM ATOMS
C
          IF (LGACX(IATM,0,0).EQ.0) GOTO 159
          DO 156 IACX=1,NACXI(IATM)
            IRCX=LGACX(IATM,IACX,0)
            IPLS=LGACX(IATM,IACX,1)
            IF (LGVAC(IRDO,IPLS)) GOTO 156
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
            IF (IESTCX(IRCX,2).NE.0) GOTO 156
C
C  PRESENTLY: PARALLEL COMPONENT OF VSIGCX(IRCX) IS NOT AVAILABLE
C             FROM FUNCTION FPATHA
C
            WTRSIG=WTR*SIGVCX(IRCX)
C  PREVIOUS BULK ION IPLS, NOW LOST.  REMOVE MODULUS OF PARALLEL MOMENTUM
            MAPL(IPLS,IRD)=MAPL(IPLS,IRD)-WTRSIG*VSIG_PARB(IPLS)
            LMETSP(NSPAMI+IPLS)=.TRUE.
C  NEW BULK ION IPL
            IF (N1STX(IRCX,1).EQ.4) THEN
              IPL=N1STX(IRCX,2)
              MAPL(IPL,IRD)=MAPL(IPL,IRD)+WTRSIG*VSIG_PARB(IPL)
              LMETSP(NSPAMI+IPL)=.TRUE.
            ENDIF
            IF (N2NDX(IRCX,1).EQ.4) THEN
              IPL=N2NDX(IRCX,2)
              MAPL(IPL,IRD)=MAPL(IPL,IRD)+WTRSIG*PARMOM_0*
     .                      SIGN(1._DP,VAL_PARB(IPL))
              LMETSP(NSPAMI+IPL)=.TRUE.
            ENDIF
156       CONTINUE
159       CONTINUE
C
C  ELECTRON IMPACT CONTRIBUTION
C
          DO 161 IAEI=1,NAEII(IATM)
            IREI=LGAEI(IATM,IAEI)
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
            IF (IESTEI(IREI,2).NE.0) GOTO 161
 
            IF (PPLEI(IREI,0).GT.0) THEN
              DO 162 IPL=1,NPLSI
                P=PPLEI(IREI,IPL)
                IF (P.GT.0) THEN
                  WTRSIG=WTR*SIGVEI(IREI)*P
C  NEW BULK ION IPL
                  MAPL(IPL,IRD)=MAPL(IPL,IRD)+WTRSIG*PARMOM_0*
     .                          SIGN(1._DP,VAL_PARB(IPL))
                  LMETSP(NSPAMI+IPL)=.TRUE.
                ENDIF
162           CONTINUE
            ENDIF
161       CONTINUE
C
C  ION IMPACT CONTRIBUTION: NOT READY
C
C
C  ELASTIC CONTRIBUTION FROM ATOMS
C
          IF (LGAEL(IATM,0,0).EQ.0) GOTO 180
C  DEFAULT TRACKLENGTH ESTIMATOR ("PERFECT IDENTITY EXCHANGE" APPROXIMATION)
          DO 181 IAEL=1,NAELI(IATM)
            IREL=LGAEL(IATM,IAEL,0)
            IPLS=LGAEL(IATM,IAEL,1)
            IBGK=NPBGKP(IPLS,1)
C
            IF (IBGK.NE.0) GOTO 181
C  THIS SPECIES IS A BGK VIRTUAL BACKGROUND SPECIES. 
C  MAPL NEEDS NOT BE UPDATED HERE, ALTHOUGH IT WOULD NOT CAUSE PROBLEMS, BECAUSE 
C       DISTINCT FROM EAPL THIS TALLY DOES HAVE A BULK SPECIES INDEX.
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
            IF (IESTEL(IREL,2).NE.0) GOTO 181
C
            WTRSIG=WTR*SIGVEL(IREL)
C
            MAPL(IPLS,IRD)=MAPL(IPLS,IRD)-WTRSIG*VSIG_PARB(IPLS)
            LMETSP(NSPAMI+IPLS)=.TRUE.
            IPL2=IPLS
            MAPL(IPL2,IRD)=MAPL(IPL2,IRD)+WTRSIG*PARMOM_0*
     .                     SIGN(1._DP,VAL_PARB(IPL2))
            LMETSP(NSPAMI+IPL2)=.TRUE.
181       CONTINUE
180       CONTINUE
C
        END IF
51    CONTINUE
      RETURN

      END
