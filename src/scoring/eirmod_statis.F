      MODULE EIRMOD_STATIS


C  may 06: bug fix: SDC initialized to zero
C  march 12: optimize calculations for surface tallies
cdr jan  18: comments
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CUPD
      USE EIRMOD_CGRID
      USE EIRMOD_CSDVI
      USE EIRMOD_COUTAU
      USE EIRMOD_CSPEI

      IMPLICIT NONE
      PRIVATE

      PUBLIC :: EIRENE_STATS0, EIRENE_STATS1,
     .          EIRENE_STATS2, EIRENE_STATS3,
cym will be removed once the parallel zone encompasses the whole code
     .          IIND
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      REAL(DP), ALLOCATABLE, SAVE :: VECTOR(:), VECTRC(:,:),
     .          SD(:),   SDC(:,:)
      INTEGER, ALLOCATABLE, SAVE :: IADD(:),    IGFF(:),
     .                              IADDW(:),   IGFFW(:),
     .                              IADDC(:,:), IGFFC(:,:),
     .                              IND(:,:),   IIND(:),    INDSS(:,:)

      REAL(DP) :: D1, DS1, D2S, DS2, DSA, DD22, DD11, D2, DD12,
     .          ZFLUXQ, DS, SD2S, SD2, SG2, SG, DA, D, DD,
     .          D2S11, D2S22, D2S12, SG12, SG1,
     .          SAV, SD1S, SD1, XNM
      INTEGER :: ISCO2, NR1, NP2, NT3, INP, IGF, IC,
     .           I, IRU, IIN, J, IR, IGS,
     .           ICO,  IS, ITL2, ISCO1, ITL1, IGE,
     .           NSYM, IGI, ITL, ISCO, NSYH, J1, J2, IP, IG, IT
      INTEGER, SAVE :: NSB, NRW

!$OMP THREADPRIVATE (vector,vectrc,iadd,igff,iaddw,igffw,
!$OMP& iaddc,igffc,ind,IIND,indss,SD,SDC,
!$OMP& D1, DS1, D2S, DS2, DSA, DD22, DD11, D2, DD12,
!$OMP& ZFLUXQ, DS, SD2S, SD2, SG2, SG, DA, D, DD,
!$OMP& D2S11, D2S22, D2S12, SG12, SG1,SAV, SD1S, SD1, XNM,
!$OMP& ISCO2, NR1, NP2, NT3, INP, IGF, IC, I, IRU, IIN, J, IR, IGS,
!$OMP& ICO,  IS, ITL2, ISCO1, ITL1, IGE,
!$OMP& NSYM, IGI, ITL, ISCO, NSYH, J1, J2, IP, IG, IT,nsb,nrw)

      CONTAINS

C
C
      SUBROUTINE EIRENE_STATS0

      IMPLICIT NONE

      INTEGER :: IIN, IR, J
      EXTERNAL :: EIRENE_INDTAL

      IMETCL = 0
      NCLMT = 0
      NCLMTS = 0
      LMETSP = .FALSE.

      IMETWL = 0
      NWLMT = 0
      NWLMTS = 0
      LMETSPW = .FALSE.

      IF (NSIGI.EQ.0) RETURN
C
      IF (.NOT.ALLOCATED(IADD)) THEN
        AllOCATE (IADD(NSD))
        AllOCATE (IGFF(NSD))
        AllOCATE (IADDW(NSDW))
        AllOCATE (IGFFW(NSDW))
        AllOCATE (IADDC(2,NCV))
        AllOCATE (IGFFC(2,NCV))
        AllOCATE (IND(NRTAL,8))
        AllOCATE (IIND(NRTAL))
        AllOCATE (INDSS(NRTAL,8))
        AllOCATE (VECTOR(MAX(NRTAL,NLMPGS)))
        AllOCATE (VECTRC(2,NRTAL))
        AllOCATE (SD(0:(MAX(NRTAL,NLMPGS))))
        AllOCATE (SDC(2,0:NRTAL))
        SD=0._DP
        SD2=0._DP
        SDC=0._DP
      END IF

C  FILL IIND, INDSS ARRAYS FOR THOSE "AVERAGE" CELLS, TO WHICH "REAL"
C  CELL IR ALSO CONTRIBUTES
C  IIND: HOW MANY CELLS
C  INDSS: WHICH CELLS
      CALL EIRENE_INDTAL(IND,NRTAL,NR1TAL,NP2TAL,NT3TAL,NBMLT)
      DO IR=1,NSBOX_TAL
        IIND(IR)=0
        IIN=0
        DO J=1,8
          IF (IND(IR,J).NE.0) THEN
            IIND(IR)=IIND(IR)+1
            IIN=IIN+1
            INDSS(IR,IIN)=J
          ENDIF
        ENDDO
      ENDDO
C
      DO 101 J=1,NSIGVI
        IGFF(J)=NFIRST(IIH(J))
        IADD(J)=NADDV(IIH(J))
  101 CONTINUE
C
      DO 102 J=1,NSIGCI
        IGFFC(1,J)=NFIRST(IIHC(1,J))
        IGFFC(2,J)=NFIRST(IIHC(2,J))
        IADDC(1,J)=NADDV(IIHC(1,J))
        IADDC(2,J)=NADDV(IIHC(2,J))
  102 CONTINUE
C
      DO 108 J=1,NSIGSI
        IGFFW(J)=NFRSTW(IIHW(J))
        IADDW(J)=NADDW(IIHW(J))
  108 CONTINUE
C
      RETURN
      END SUBROUTINE EIRENE_STATS0
C
C
C
      SUBROUTINE EIRENE_STATS1(NBIN,NRIN,NPIN,NTIN,NSIN,LP,LT)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NBIN, NRIN, NPIN, NTIN, NSIN
      LOGICAL, INTENT(IN) :: LP, LT
C
cym      REAL(DP) :: SAV, SD1S, SD1, SD2S
cym      INTEGER :: NR1, NP2, NT3, I, J, IRU, IGS, IC, IIN, INP, IGF, ITL,
cym     .           ISCO, IGE, IGI, NSYM, NSYH, J1, J2, IR, IP, IG, IT,
cym     .           ICO, ISCO2, ITL2, ISCO1, ITL1

      NSB=NBIN
      NR1=NRIN
      NP2=NPIN
      NT3=NTIN
      NRW=NSIN

      IF ((NSIGVI > 0) .OR. (NSIGCI > 0)) THEN
C  THE PRESENT HISTORY HAS TOUCHED NCLMT CELLS.
C  IT CONTRIBUTES ADDITIONALLY TO FURTHER CELLS (CONTAINING AVERAGES)
C  THESE ADDITIONAL CELL NUMBERS FOR AVERAGES ARE STORED HERE ON ICLMT ARRAY
C  FILL ICLMT(NCLMT+1),....ICLMT(NCLMTS)
        NCLMTS = NCLMT
        DO I=1,NCLMT
          IR = ICLMT(I)
          DO IIN=2,IIND(IR)
            J=INDSS(IR,IIN)
            IRU=IND(IR,J)
            IF (IMETCL(IRU) == 0) THEN
              NCLMTS = NCLMTS+1
              IMETCL(IRU) = NCLMTS
              ICLMT(NCLMTS) = IRU
            END IF
          END DO
        END DO
      END IF
C
C
      IF (NSIGVI.EQ.0) GOTO 1020
C
C  ARE THERE CONTRIBUTIONS TO THE REQUESTED TALLY FROM THIS HISTORY?
      DO 1012 IC=1,NSIGVI
        INP=IADD(IC)
        IGF=IGFF(IC)
        IGS=IGH(IC)
        ITL=IIH(IC)
        IF (NSPAN(ITL) == 0) THEN
          ISCO = 1
        ELSE
          ISCO = 0
          IF (IGS == 0) THEN
            IF ( ANY(LMETSP(NSPAN(ITL):NSPEN(ITL))) ) ISCO = 1
          ELSE
            IF (LMETSP(NSPAN(ITL)+IGS-1)) ISCO = 1
          END IF
        END IF
        IF (ISCO == 0) GOTO 1012
C
        IF (.NOT.LP.AND..NOT.LT) GOTO 1005
C  USE SYMMETRY IN POLOIDAL/Y AND/OR TOROIDAL/Z COORDINATE
        IF (IGS.LE.0) THEN
          IGI=1
          IGE=IGF
        ELSE
          IGI=IGS
          IGE=IGS
        ENDIF
        IF (LP) THEN
          NSYM=NP2
          NSYH=(NSYM-1)/2
          DO 1003 IG=IGI,IGE
            DO IR=1,NR1
              DO IT=1,NT3
                DO IP=1,NSYH
                  J1=IR+((IT-1)*NP2+IP-1)*NR1
                  J2=IR+((IT-1)*NP2+NSYM-IP-1)*NR1
                  SAV=(ESTIMV(INP+IG,J1)+ESTIMV(INP+IG,J2))*0.5
                  ESTIMV(INP+IG,J1)=SAV
                  ESTIMV(INP+IG,J2)=SAV
                END DO
              END DO
            END DO
 1003     CONTINUE
        ENDIF
        IF (LT) THEN
          NSYM=NT3
          NSYH=(NSYM-1)/2
          DO 1004 IG=IGI,IGE
            DO IR=1,NR1
              DO IP=1,NP2
                DO IT=1,NSYH
                  J1=IR+((IT-1)*NP2+IP-1)*NR1
                  J2=IR+((NSYM-IT-1)*NP2+IP-1)*NR1
                  SAV=(ESTIMV(INP+IG,J1)+ESTIMV(INP+IG,J2))*0.5
                  ESTIMV(INP+IG,J1)=SAV
                  ESTIMV(INP+IG,J2)=SAV
                END DO
              END DO
            END DO
 1004     CONTINUE
        ENDIF
 1005   CONTINUE
C
C  FILL ARRAY VECTOR, EITHER FOR SUM OVER SPECIES OR FOR INDIVIDUAL SPECIES
C  VECTOR IS FILLED ONLY FOR THOSE CELLS,
C  WHICH HAVE BEEN TOUCHED BY THIS HISTORY
C  VECTOR CONTAINS THE SUM FROM ALL HISTORIES IN THESE CELLS UP TO THE
C  PRESENT HISTORY
        IF (IGS.NE.0) THEN
           DO ICO = 1,NCLMT
             IR = ICLMT(ICO)
             VECTOR(ICO)=ESTIMV(INP+IGS,IR)
           END DO
        ELSE
          DO 1014 ICO=1,NCLMT
            VECTOR(ICO)=0.
 1014     CONTINUE
          DO 1015 IS=1,IGF
           DO ICO=1,NCLMT
            IR = ICLMT(ICO)
            VECTOR(ICO)=VECTOR(ICO)+ESTIMV(INP+IS,IR)
           END DO
 1015     CONTINUE
        ENDIF

        SD1S = 0.D0
C  FILL ARRAY SD WITH THE INDIVIDUAL CONTRIBUTION FROM THIS HISTORY,
C  IN EACH CELL THAT HAS BEEN TOUCHED BY THIS HISTORY
        DO ICO = 1,NCLMT
          IR = ICLMT(ICO)
          SD1 = VECTOR(ICO)-SDVIA(IC,IR)
          SD1S=SD1S+SD1
          SDVIA(IC,IR)=VECTOR(ICO)
          SD(IR) = SD1
C  SUM SD ALSO FOR OTHER "CELL", TO WHICH CELL IR CONTRIBUTES
C  I.E., AVERAGES OVER COORDINATES OR OVER THE ENTIRE COMPUTATIONAL DOMAIN
          DO IIN=2,IIND(IR)
            J=INDSS(IR,IIN)
            IRU=IND(IR,J)
            SD(IRU)=SD(IRU)+SD1
          END DO
        END DO

        DO ICO = 1,NCLMTS
          IR = ICLMT(ICO)
          SD1=SD(IR)
          SIGMA(IC,IR)=SIGMA(IC,IR)+SD1*SD1
          SD(IR)=0._DP
        END DO
        SGMS(IC)=SGMS(IC)+SD1S*SD1S
 1012 CONTINUE
c  updating of statistical variance (once per flight) for volume tallies: done

C
C  update statistical variance for surface tallies, once after each flight
 1020 CONTINUE
      IF (NSIGSI.EQ.0) GOTO 1030
      DO 1022 IC=1,NSIGSI
        INP=IADDW(IC)
        IGF=IGFFW(IC)
        IGS=IGHW(IC)
        ITL=IIHW(IC)
c  variance for tally itl, species igs.  in case igs=0: sum over species
        IF (NSPANW(ITL) == 0) THEN
c  tally itl has no (first) species index
          ISCO = 1
        ELSE
          ISCO = 0
          IF (IGS == 0) THEN
c  sum over species
            IF ( ANY(LMETSPW(NSPANW(ITL):NSPENW(ITL))) ) ISCO = 1
          ELSE
c  only for species igs
            IF (LMETSPW(NSPANW(ITL)+IGS-1)) ISCO = 1
          END IF
        END IF
        IF (ISCO == 0) GOTO 1022  ! do not score tally ITL,
c                                   because no relevant species in this trajectory

c  fill 'vector'
c  vector is cumulated contribution after present flight no. n
c
        IF (IGS.NE.0) THEN
c
          DO 1023 ICO=1,NWLMT
c  tally is for single species igs
            IR = IWLMT(ICO)
            VECTOR(ICO)=ESTIMS(INP+IGS,IR)
 1023     CONTINUE
        ELSE
c  tally is for sum over species
          DO 1024 ICO=1,NWLMT
            VECTOR(ICO)=0.
 1024     CONTINUE
          DO 1025 IS=1,IGF
           DO ICO=1,NWLMT
            IR = IWLMT(ICO)
            VECTOR(ICO)=VECTOR(ICO)+ESTIMS(INP+IS,IR)
           END DO
 1025     CONTINUE
        ENDIF

c  next:
C  contribution from current flight no. n only: sd1 =vector-sdviaw
c  sdvia  is cumulated contribution after previous flight no. n-1 (previous call)

        SD1S=0.
        DO 1021 ICO=1,NWLMT
          IR = IWLMT(ICO)
c
          SD1=VECTOR(ICO)-SDVIAW(IC,IR)
          SD1S=SD1S+SD1
          SIGMAW(IC,IR)=SIGMAW(IC,IR)+SD1*SD1
          SDVIAW(IC,IR)=VECTOR(ICO)
 1021   CONTINUE
        SGMWS(IC)=SGMWS(IC)+SD1S*SD1S
 1022 CONTINUE
c  sigma now is cumulated squared contribution after flight no. n
C
 1030 CONTINUE
C

c  surface tallies done.  next: covariances
      IF (NSIGCI.EQ.0) GOTO 1050
C
      DO 1032 IC=1,NSIGCI
        ITL1=IIHC(1,IC)
        ITL2=IIHC(2,IC)
        IF (NSPAN(ITL1) == 0) THEN
          ISCO1 = 1
        ELSE
          ISCO1 = 0
          IF (IGHC(1,IC) == 0) THEN
            IF ( ANY(LMETSP(NSPAN(ITL1):NSPEN(ITL1))) ) ISCO1 = 1
          ELSE
            IF (LMETSP(NSPAN(ITL1)+IGHC(1,IC)-1)) ISCO1 = 1
          END IF
        END IF
        IF (NSPAN(ITL2) == 0) THEN
          ISCO2 = 1
        ELSE
          ISCO2 = 0
          IF (IGHC(2,IC) == 0) THEN
            IF ( ANY(LMETSP(NSPAN(ITL2):NSPEN(ITL2))) ) ISCO2 = 1
          ELSE
            IF (LMETSP(NSPAN(ITL2)+IGHC(2,IC)-1)) ISCO2 = 1
          END IF
        END IF
        IF (ISCO1+ISCO2 == 0) GOTO 1032
C
        DO 1037 I=1,2
          INP=IADDC(I,IC)
          IGF=IGFFC(I,IC)
          IGS=IGHC(I,IC)
C
          IF (.NOT.LP.AND..NOT.LT) GOTO 1035
          IF (IGS.LE.0) THEN
            IGI=1
            IGE=IGF
          ELSE
            IGI=IGS
            IGE=IGS
          ENDIF
          IF (LP) THEN
            NSYM=NP2
            NSYH=(NSYM-1)/2
            DO 1033 IG=IGI,IGE
             DO IR=1,NR1
              DO IT=1,NT3
               DO IP=1,NSYH
                  J1=IR+((IT-1)*NP2+IP-1)*NR1
                  J2=IR+((IT-1)*NP2+NSYM-IP-1)*NR1
                  SAV=(ESTIMV(INP+IG,J1)+ESTIMV(INP+IG,J2))*0.5
                  ESTIMV(INP+IG,J1)=SAV
                  ESTIMV(INP+IG,J2)=SAV
               END DO
              END DO
             END DO
 1033       CONTINUE
          ENDIF
          IF (LT) THEN
            NSYM=NT3
            NSYH=(NSYM-1)/2
            DO 1034 IG=IGI,IGE
             DO IR=1,NR1
              DO IP=1,NP2
               DO IT=1,NSYH
                  J1=IR+((IT-1)*NP2+IP-1)*NR1
                  J2=IR+((NSYM-IT-1)*NP2+IP-1)*NR1
                  SAV=(ESTIMV(INP+IG,J1)+ESTIMV(INP+IG,J2))*0.5
                  ESTIMV(INP+IG,J1)=SAV
                  ESTIMV(INP+IG,J2)=SAV
               END DO
              END DO
             END DO
 1034       CONTINUE
          ENDIF
 1035     CONTINUE
C
          IF (IGS.NE.0) THEN
            DO ICO = 1,NCLMT
              IR = ICLMT(ICO)
              VECTRC(I,ICO)=ESTIMV(INP+IGS,IR)
            END DO
          ELSE
            DO 1044 ICO=1,NCLMT
              VECTRC(I,ICO)=0.
 1044       CONTINUE
            DO 1045 IS=1,IGF
             DO ICO=1,NCLMT
              IR = ICLMT(ICO)
              VECTRC(I,ICO)=VECTRC(I,ICO)+ESTIMV(INP+IS,IR)
             END DO
 1045       CONTINUE
          ENDIF
 1037   CONTINUE
C
C
        SD1S = 0.D0
        SD2S = 0.D0
        DO ICO = 1,NCLMT
          IR = ICLMT(ICO)
          SD1 = VECTRC(1,ICO)-SDVIAC(1,IC,IR)
          SD2 = VECTRC(2,ICO)-SDVIAC(2,IC,IR)
          SD1S=SD1S+SD1
          SD2S=SD2S+SD2
          SDVIAC(1,IC,IR)=VECTRC(1,ICO)
          SDVIAC(2,IC,IR)=VECTRC(2,ICO)
          SDC(1,IR) = SD1
          SDC(2,IR) = SD2
          DO IIN=2,IIND(IR)
            J=INDSS(IR,IIN)
            IRU=IND(IR,J)
            SDC(1,IRU)=SDC(1,IRU)+SD1
            SDC(2,IRU)=SDC(2,IRU)+SD2
          END DO
        END DO
C
        DO ICO = 1,NCLMTS
          IR = ICLMT(ICO)
          SD1=SDC(1,IR)
          SD2=SDC(2,IR)
          SIGMAC(0,IC,IR)=SIGMAC(0,IC,IR)+SD1*SD2
          SIGMAC(1,IC,IR)=SIGMAC(1,IC,IR)+SD1*SD1
          SIGMAC(2,IC,IR)=SIGMAC(2,IC,IR)+SD2*SD2
          SDC(1:2,IR)=0._DP
        END DO
        SGMCS(0,IC)=SGMCS(0,IC)+SD1S*SD2S
        SGMCS(1,IC)=SGMCS(1,IC)+SD1S*SD1S
        SGMCS(2,IC)=SGMCS(2,IC)+SD2S*SD2S
 1032 CONTINUE
C
C
 1050 CONTINUE
      RETURN
      END SUBROUTINE EIRENE_STATS1
C
C
C
c  scale statistical variance. called after all flights from a given stratum istra

      SUBROUTINE EIRENE_STATS2(XN,FSIG,ZFLUX)
      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: XN, FSIG, ZFLUX

cym      REAL(DP) :: D1, DS1, D2S, DS2, DSA, DD22, DD11, D2, DD12,
cym     .          ZFLUXQ, DS, SG2, SG, DA, D, DD,
cym     .          D2S11, D2S22, D2S12, SG12, SG1,
cym     .          XNM, SD1
cym      INTEGER :: IC, IR, IS, IIN, I

C
C  1. FALL  ALLE BEITRAEGE GLEICHES VORZEICHEN: SIG ZWISCHEN 0 UND 1
C           (=1, FALLS NUR EIN BEITRAG UNGLEICH 0, ODER (KUENSTLICH
C            ERZWUNGEN) FALLS GAR KEIN BEITRAG UNGLEICH NULL)
C  2. FALL  NEGATIVE UND POSITIVE BEITRAGE KOMMEN VOR:
C           LT. FORMEL SIND AUCH WERTE GROESSER 1 MOEGLICH.
C
      XNM=XN-1.
      IF (XNM.LE.0.D0) RETURN
      ZFLUXQ=ZFLUX*ZFLUX
C
      IF (NSIGVI.EQ.0) GOTO 2200
C
      DO 2112 IC=1,NSIGVI
        INP=IADD(IC)
        IGF=IGFF(IC)
        IGS=IGH(IC)
        IF (IGS.NE.0) THEN
          DO 2113 IR=1,NSB
            VECTOR(IR)=ESTIMV(INP+IGS,IR)
 2113     CONTINUE
        ELSE
          DO 2114 IR=1,NSB
            VECTOR(IR)=0.
 2114     CONTINUE
          DO 2115 IS=1,IGF
           DO IR=1,NSB
            VECTOR(IR)=VECTOR(IR)+ESTIMV(INP+IS,IR)
           END DO
 2115     CONTINUE
        ENDIF
C  tally is now on vector
        DS=0.
        DO 2011 IR=1,NSB
          SD1=VECTOR(IR)
          DS=DS+SD1
          DO 2016 IIN=1,IIND(IR)
            J=INDSS(IR,IIN)
            IRU=IND(IR,J)
            SD(IRU)=SD(IRU)+SD1
 2016     CONTINUE
 2011   CONTINUE
C
        DO 2111 IR=1,NSB
          D=SD(IR)
          DD=D*D
          DA=ABS(D)
          SG2=MAX(0._DP,SIGMA(IC,IR)-DD/XN)
C RELATIVE STANDARD DEVIATION FOR CURRENT STRATUM
          SG=SQRT(SG2)/(DA+EPS60)
          SIGMA(IC,IR)=SG*FSIG
C CUMULATED VARIANCE FOR SUM OVER STRATA
          STV(IC,IR)=STV(IC,IR)+SG2*ZFLUXQ/XNM/XN
          EE(IC,IR)=EE(IC,IR)+D*ZFLUX/XN
          SD(IR)=0._DP
 2111   CONTINUE
        D2S=DS*DS
        DSA=ABS(DS)
        SG2=MAX(0._DP,SGMS(IC)-D2S/XN)
        SG=SQRT(SG2)/(DSA+EPS60)
        SGMS(IC)=SG*FSIG
C
        STVS(IC)=STVS(IC)+SG2*ZFLUXQ/XNM/XN
        EES(IC)=EES(IC)+DS*ZFLUX/XN
 2112 CONTINUE


C
 2200 CONTINUE
      IF (NSIGSI.EQ.0) GOTO 2300
c  put requested tally on 'vector'
      DO 2212 IC=1,NSIGSI
        INP=IADDW(IC)
        IGF=IGFFW(IC)
        IGS=IGHW(IC)
        DS=0.
        IF (IGS.NE.0) THEN
c  tally for individual species
          DO 2213 IR=1,NRW
            VECTOR(IR)=ESTIMS(INP+IGS,IR)
 2213     CONTINUE
        ELSE
c  tally for sum over species
          DO 2214 IR=1,NRW
            VECTOR(IR)=0.
 2214     CONTINUE
          DO 2215 IS=1,IGF
           DO IR=1,NRW
            VECTOR(IR)=VECTOR(IR)+ESTIMS(INP+IS,IR)
           END DO
 2215     CONTINUE
        ENDIF

        DO 2211 IR=1,NRW
          D=VECTOR(IR)
          DS=DS+D
          DD=D*D
          DA=ABS(D)
          SG2=MAX(0._DP,SIGMAW(IC,IR)-DD/XN)
C RELATIVE STANDARD DEVIATION FOR CURRENT STRATUM
          SG=SQRT(SG2)/(DA+EPS60)
          SIGMAW(IC,IR)=SG*FSIG
C CUMULATED VARIANCE FOR SUM OVER STRATA
          STVW(IC,IR)=STVW(IC,IR)+SG2*ZFLUXQ/XNM/XN
          FF(IC,IR)=FF(IC,IR)+D*ZFLUX/XN
 2211   CONTINUE
        D2S=DS*DS
        DSA=ABS(DS)
        SG2=MAX(0._DP,SGMWS(IC)-D2S/XN)
        SG=SQRT(SG2)/(DSA+EPS60)
        SGMWS(IC)=SG*FSIG
C
        STVWS(IC)=STVWS(IC)+SG2*ZFLUXQ/XNM/XN
        FFS(IC)=FFS(IC)+DS*ZFLUX/XN
 2212 CONTINUE
C
 2300 CONTINUE
C

c  now deal with covariances
      IF (NSIGCI.EQ.0) GOTO 2400
C
      DO 2312 IC=1,NSIGCI
        DO 2317 I=1,2
          INP=IADDC(I,IC)
          IGF=IGFFC(I,IC)
          IGS=IGHC(I,IC)
          IF (IGS.NE.0) THEN
            DO 2313 IR=1,NSB
              VECTRC(I,IR)=ESTIMV(INP+IGS,IR)
 2313       CONTINUE
          ELSE
            DO 2314 IR=1,NSB
              VECTRC(I,IR)=0.
 2314       CONTINUE
            DO 2315 IS=1,IGF
             DO IR=1,NSB
              VECTRC(I,IR)=VECTRC(I,IR)+ESTIMV(INP+IS,IR)
             END DO
 2315       CONTINUE
          ENDIF
 2317   CONTINUE
C
        DS1=0.
        DS2=0.
        DO 2311 IR=1,NSB
          SD1=VECTRC(1,IR)
          SD2=VECTRC(2,IR)
          DS1=DS1+SD1
          DS2=DS2+SD2
          DO 2316 IIN=1,IIND(IR)
            J=INDSS(IR,IIN)
            IRU=IND(IR,J)
            SDC(1,IRU)=SDC(1,IRU)+SD1
            SDC(2,IRU)=SDC(2,IRU)+SD2
 2316     CONTINUE
 2311   CONTINUE
        DO 2411 IR=1,NSB
          D1=SDC(1,IR)
          D2=SDC(2,IR)
          DD12=D1*D2
          DD11=D1*D1
          DD22=D2*D2
          SG12=         SIGMAC(0,IC,IR)-DD12/XN
          SG1 =MAX(0._DP,SIGMAC(1,IC,IR)-DD11/XN)
          SG2 =MAX(0._DP,SIGMAC(2,IC,IR)-DD22/XN)
C ABSOLUTE STANDARD DEVIATION AND COVARIANCES
          SIGMAC(0,IC,IR)=SG12/XNM/XN
          SIGMAC(1,IC,IR)=SQRT(SG1/XNM/XN)
          SIGMAC(2,IC,IR)=SQRT(SG2/XNM/XN)
          SDC(1:2,IR)=0._DP
 2411   CONTINUE
        D2S12=DS1*DS2
        D2S11=DS1*DS1
        D2S22=DS2*DS2
        SG12=         SGMCS(0,IC)-D2S12/XN
        SG1 =MAX(0._DP,SGMCS(1,IC)-D2S11/XN)
        SG2 =MAX(0._DP,SGMCS(2,IC)-D2S22/XN)
        SGMCS(0,IC)=SG12/XNM/XN
        SGMCS(1,IC)=SQRT(SG1/XNM/XN)
        SGMCS(2,IC)=SQRT(SG2/XNM/XN)
 2312 CONTINUE
C
 2400 RETURN
      END SUBROUTINE EIRENE_STATS2
C
C
C
      SUBROUTINE EIRENE_STATS3
      IMPLICIT NONE

      IF (ALLOCATED(IADD)) THEN
         DEALLOCATE (IADD)
         DEALLOCATE (IGFF)
         DEALLOCATE (IADDW)
         DEALLOCATE (IGFFW)
         DEALLOCATE (IADDC)
         DEALLOCATE (IGFFC)
         DEALLOCATE (IND)
         DEALLOCATE (IIND)
         DEALLOCATE (INDSS)
         DEALLOCATE (VECTOR)
         DEALLOCATE (VECTRC)
         DEALLOCATE (SD)
         DEALLOCATE (SDC)
      END IF

      RETURN
      END SUBROUTINE EIRENE_STATS3

      END MODULE EIRMOD_STATIS
