C  statistical variance of spectra
C  NOTE:
C  distinct from the other variances (volume tallies, surface tallies, bgk and cop tallies)
c  here in case of spectra tallies the variances are contained in the same structure (ESTIML)
c  as the tallies themselves.
c  nomenclature, however has been syncronized (oct. 2014)
c  e.g.  ESTIML(ISPC)%PSPC%SGM  <--> sigma, sigmaw
c        ESTIML(ISPC)%PSPC%SDV  <--> sdvia, sdviaw
c  etc.
c and similarly for the intermediate storage structure SMESTL
c  e.g.  SMESTL(ISPC)%PSPC%GG   <--> ee, ff
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
      USE EIRMOD_CSDVI_COP
      USE EIRMOD_COUTAU
      USE EIRMOD_COMSOU
 
      IMPLICIT NONE
 
      REAL(DP), INTENT(IN) :: XN, FSIG, ZFLUX
      INTEGER, INTENT(IN) :: NBIN, NRIN, NPIN, NTIN, NSIN
      LOGICAL, INTENT(IN) :: LP, LT
 
      REAL(DP), ALLOCATABLE :: SD(:)
      REAL(DP) :: XNM, DS, ZFLUXQ, D2S, SG,
     .            DSA, DD, D, SG2, DA, SD1, SD1S
      INTEGER :: NSB, NP2, NR1, NT3, I, ISPC,NSPECI,NSPECE
C
      SAVE
C
      ENTRY EIRENE_STATS0_SPC
 
      IF (NADSPC > 0) THEN
        DO ISPC=1,NADSPC
          ESTIML(ISPC)%PSPC%IMETSP = 0
        END DO
      END IF
C
      RETURN
 
C
      ENTRY EIRENE_STATS1_SPC(NBIN,NRIN,NPIN,NTIN,NSIN,LP,LT)
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
c  vector = ESTIML(ISPC)%PSPC%SPC(I) is cummulated contribution after present (n-th) flight

c  sdvia = ESTIML(ISPC)%PSPC%SDV(I) is cummulated contribution after previous flight no. n-1 (previous call)
C  contribution from current flight no. n only: sd1 =vector-sdviaw
C
        IF (ESTIML(ISPC)%PSPC%IMETSP > 0) THEN
          SD1S=0.
c  size of tally ISPC, ADD BIN 0 AND BIN NSPC+1 for low and high end of spectrum
          NSPECI=0
          NSPECE=ESTIML(ISPC)%PSPC%NSPC+1
          ALLOCATE (SD(NSPECI:NSPECE))
          SD = 0.D0
          DO I = NSPECI,NSPECE
            SD1=ESTIML(ISPC)%PSPC%SPC(I)-ESTIML(ISPC)%PSPC%SDV(I)
            SD1S=SD1S+SD1
            ESTIML(ISPC)%PSPC%SDV(I)=ESTIML(ISPC)%PSPC%SPC(I)
            SD(I) = SD1
          END DO
c  now  sdvia = ESTIML(ISPC)%PSPC%SDV(I) is cummulated contribution after present flight no. n

          DO I = NSPECI,NSPECE
            SD1=SD(I)
            ESTIML(ISPC)%PSPC%SGM(I)=ESTIML(ISPC)%PSPC%SGM(I)+SD1*SD1
          END DO
c  sigma = ESTIML(ISPC)%PSPC%SGM(I)  now is cummulated squared contribution after flight no. n
          ESTIML(ISPC)%PSPC%SGMS=ESTIML(ISPC)%PSPC%SGMS+SD1S*SD1S
          DEALLOCATE (SD)
        END IF
      END DO
C
C
      RETURN


C  next entry:
c  scale statistical variance. called after all flights from a given stratum istra
      ENTRY EIRENE_STATS2_SPC(XN,FSIG,ZFLUX)
C
C  1. FALL  ALLE BEITRAEGE GLEICHES VORZEICHEN: SIG ZWISCHEN 0 UND 1
C           (=1, FALLS NUR EIN BEITRAG UNGLEICH 0, ODER (KUENSTLICH
C            ERZWUNGEN) FALLS GAR KEIN BEITRAG UNGLEICH NULL)
C  2. FALL  NEGATIVE UND POSITIVE BEITRAGE KOMMEN VOR:
C           LT. FORMEL SIND AUCH WERTE GROESSER 1  MOEGLICH.
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
        NSPECE=ESTIML(ISPC)%PSPC%NSPC+1
        ALLOCATE (SD(NSPECI:NSPECE))
C ESTIML(ISPC)%PSPC%SPC cummulated tally score after all flights from present stratum istra 
        SD=ESTIML(ISPC)%PSPC%SPC
        DS=SUM(SD)
 
        DO I=NSPECI,NSPECE
          D=SD(I)
          DD=D*D
          DA=ABS(D)
          SG2=MAX(0._DP,ESTIML(ISPC)%PSPC%SGM(I)-DD/XN)
C RELATIV STANDARD DEVIATION
          SG=SQRT(SG2)/(DA+EPS60)
          ESTIML(ISPC)%PSPC%SGM(I)=SG*FSIG

C CUMMULATED VARIANCE FOR SUM OVER STRATA
! STV, CORRESPONDS TO STV, STVW
          IF ((NSMSTRA > 0 ) .AND. (NSTRAI > 1)) THEN
            SMESTL(ISPC)%PSPC%STV(I)=SMESTL(ISPC)%PSPC%STV(I)+
     .                               SG2*ZFLUXQ/XNM/XN
! GG, CORRESPONDS TO EE, FF
            SMESTL(ISPC)%PSPC%GG(I)=SMESTL(ISPC)%PSPC%GG(I)+D*ZFLUX/XN
          END IF
        END DO

c variance of summed (total) contribution: STVS, GGS
        D2S=DS*DS
        DSA=ABS(DS)
        SG2=MAX(0._DP,ESTIML(ISPC)%PSPC%SGMS-D2S/XN)
        SG=SQRT(SG2)/(DSA+EPS60)
        ESTIML(ISPC)%PSPC%SGMS=SG*FSIG
C
        IF ((NSMSTRA > 0 ) .AND. (NSTRAI > 1)) THEN
          SMESTL(ISPC)%PSPC%STVS=SMESTL(ISPC)%PSPC%STVS+
     .                           SG2*ZFLUXQ/XNM/XN
          SMESTL(ISPC)%PSPC%GGS =SMESTL(ISPC)%PSPC%GGS+DS*ZFLUX/XN
        END IF
 
        DEALLOCATE (SD)
      END DO
C
2200  CONTINUE
      RETURN
      END
