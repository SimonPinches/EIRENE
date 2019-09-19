C  statistical variance of spectra
C  NOTE:
C  distinct from the other variances (volume tallies, surface tallies, bgk and cop tallies)
c  here in case of spectra tallies the variances are contained in the same structure (ESTIML)
c  as the tallies themselves.
c  nomenclature, however has been synchronized (oct. 2014)
c  e.g.  ESTIML(ISPC)%SGM  <--> sigma, sigmaw
c        ESTIML(ISPC)%SDV  <--> sdvia, sdviaw
c  etc.
c and similarly for the intermediate storage structure SMESTL
c  e.g.  SMESTL(ISPC)%GG   <--> ee, ff
c  etc.
C
      SUBROUTINE EIRENE_STATIS_SPC
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_CCONA
      USE EIRMOD_CGRID
      USE EIRMOD_CSDVI
      USE EIRMOD_COUTAU
      USE EIRMOD_COMSOU

      IMPLICIT NONE

      CALL EIRENE_STATS0_SPC
      END SUBROUTINE EIRENE_STATIS_SPC

C
      SUBROUTINE EIRENE_STATS0_SPC
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_CCONA
      USE EIRMOD_CGRID
      USE EIRMOD_CSDVI
      USE EIRMOD_COUTAU
      USE EIRMOD_COMSOU
      IMPLICIT NONE

      REAL(DP), ALLOCATABLE :: SD(:)
      INTEGER :: ISPC

      SAVE
      IF (NADSPC > 0) THEN
        DO ISPC=1,NADSPC
          ESTIML(ISPC)%IMETSP = 0
        END DO
      END IF
C
      END SUBROUTINE EIRENE_STATS0_SPC

C
      SUBROUTINE EIRENE_STATS1_SPC(NBIN,NRIN,NPIN,NTIN,NSIN,LP,LT)
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_CCONA
      USE EIRMOD_CGRID
      USE EIRMOD_CSDVI
      USE EIRMOD_COUTAU
      USE EIRMOD_COMSOU
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NBIN, NRIN, NPIN, NTIN, NSIN
      LOGICAL, INTENT(IN) :: LP, LT

      REAL(DP), ALLOCATABLE :: SD(:)
      REAL(DP) :: SD1, SD1S
      INTEGER :: NSB, NP2, NR1, NT3, I, ISPC, NSPECI, NSPECE

      NSB=NBIN
      NR1=NRIN
      NP2=NPIN
      NT3=NTIN
C
C
      IF (NADSPC.EQ.0) RETURN
C
C
C  STATISTICS FOR SPECTRA

      DO ISPC=1,NADSPC
c  vector = ESTIML(ISPC)%SPC(I) is cumulated contribution after present (n-th) flight

c  sdvia = ESTIML(ISPC)%SDV(I) is cumulated contribution after previous flight no. n-1 (previous call)
C  contribution from current flight no. n only: sd1 =vector-sdviaw
C
        IF (ESTIML(ISPC)%IMETSP > 0) THEN
          SD1S=0.
c  size of tally ISPC, ADD BIN 0 AND BIN NSPC+1 for low and high end of spectrum
          NSPECI=0
          NSPECE=ESTIML(ISPC)%NSPC+1
          ALLOCATE (SD(NSPECI:NSPECE))
          SD = 0.D0
          DO I = NSPECI,NSPECE
            SD1=ESTIML(ISPC)%SPC(I)-ESTIML(ISPC)%SDV(I)
            SD1S=SD1S+SD1
            ESTIML(ISPC)%SDV(I)=ESTIML(ISPC)%SPC(I)
            SD(I) = SD1
          END DO
c  now  sdvia = ESTIML(ISPC)%SDV(I) is cumulated contribution after present flight no. n

          DO I = NSPECI,NSPECE
            SD1=SD(I)
            ESTIML(ISPC)%SGM(I)=ESTIML(ISPC)%SGM(I)+SD1*SD1
          END DO
c  sigma = ESTIML(ISPC)%SGM(I)  now is cumulated squared contribution after flight no. n
          ESTIML(ISPC)%SGMS=ESTIML(ISPC)%SGMS+SD1S*SD1S
          DEALLOCATE (SD)
        END IF
      END DO
C
C
      END SUBROUTINE EIRENE_STATS1_SPC


C  next entry:
c  scale statistical variance. called after all flights from a given stratum istra
      SUBROUTINE EIRENE_STATS2_SPC(XN,FSIG,ZFLUX)
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_CCONA
      USE EIRMOD_CGRID
      USE EIRMOD_CSDVI
      USE EIRMOD_COUTAU
      USE EIRMOD_COMSOU
      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: XN, FSIG, ZFLUX

      REAL(DP), ALLOCATABLE :: SD(:)
      REAL(DP) :: XNM, DS, ZFLUXQ, D2S, SG, SG2, DSA, D, DA, DD
      INTEGER I, ISPC, NSPECI, NSPECE

      SAVE
C
C  1. FALL  ALLE BEITRAEGE GLEICHES VORZEICHEN: SIG ZWISCHEN 0 UND 1
C           (=1, FALLS NUR EIN BEITRAG UNGLEICH 0, ODER (KUENSTLICH
C            ERZWUNGEN) FALLS GAR KEIN BEITRAG UNGLEICH NULL)
C  2. FALL  NEGATIVE UND POSITIVE BEITRAGE KOMMEN VOR:
C           LT. FORMEL SIND AUCH WERTE GROESSER 1 MOEGLICH.
C
      XNM=XN-1.
      IF (XNM.LE.0.) RETURN
      ZFLUXQ=ZFLUX*ZFLUX
C
      IF (NADSPC.EQ.0) RETURN
C
C  STATISTICS FOR SPECTRA
      DO ISPC=1,NADSPC
C
c  size of tally ISPC, ADD BIN 0 AND BIN NSPC+1 for low and high end of spectrum
        NSPECI=0
        NSPECE=ESTIML(ISPC)%NSPC+1
        ALLOCATE (SD(NSPECI:NSPECE))
C ESTIML(ISPC)%SPC cumulated tally score after all flights from present stratum istra
        SD=ESTIML(ISPC)%SPC
        DS=SUM(SD)

        DO I=NSPECI,NSPECE
          D=SD(I)
          DD=D*D
          DA=ABS(D)
          SG2=MAX(0._DP,ESTIML(ISPC)%SGM(I)-DD/XN)
C RELATIVE STANDARD DEVIATION
          SG=SQRT(SG2)/(DA+EPS60)
          ESTIML(ISPC)%SGM(I)=SG*FSIG

C CUMULATED VARIANCE FOR SUM OVER STRATA
! STV, CORRESPONDS TO STV, STVW
          IF ((NSMSTRA > 0 ) .AND. (NSTRAI > 1)) THEN
            SMESTL(ISPC)%STV(I)=SMESTL(ISPC)%STV(I)+
     .                               SG2*ZFLUXQ/XNM/XN
! GG, CORRESPONDS TO EE, FF
            SMESTL(ISPC)%GG(I)=SMESTL(ISPC)%GG(I)+D*ZFLUX/XN
          END IF
        END DO

c variance of summed (total) contribution: STVS, GGS
        D2S=DS*DS
        DSA=ABS(DS)
        SG2=MAX(0._DP,ESTIML(ISPC)%SGMS-D2S/XN)
        SG=SQRT(SG2)/(DSA+EPS60)
        ESTIML(ISPC)%SGMS=SG*FSIG
C
        IF ((NSMSTRA > 0 ) .AND. (NSTRAI > 1)) THEN
          SMESTL(ISPC)%STVS=SMESTL(ISPC)%STVS+
     .                           SG2*ZFLUXQ/XNM/XN
          SMESTL(ISPC)%GGS =SMESTL(ISPC)%GGS+DS*ZFLUX/XN
        END IF

        DEALLOCATE (SD)
      END DO
C
      RETURN
      END SUBROUTINE EIRENE_STATS2_SPC
