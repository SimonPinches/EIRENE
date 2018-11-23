CDR:  evaluate standard deviation for specific tallies needed for coupling
CDR   which would not be available otherwise.
cdr   in early 2014 the sum over atomic (a) molecular (m) and test ion (i) components
cdr   for particle, momentum and energy sources was removed here, so this routine is currently
cdr   empty.
cdr   These standard deviations are now, together with other linear combinations of default
cdr   tallies, obtained "on the fly" by scoring per history, (in upfcop).




!PB  17.11.05  USAGE OF SIGMA_COP CHANGED
!PB            SIGMA_COP(      1:  NPLSI) : STATISTICS FOR MOMENTUM SOURCES
!PB            SIGMA_COP(NPLSI+1:2*NPLSI) : STATISTICS FOR PARTICLE SOURCES
!PB            SIGMA_COP(2*NPLSI+1)       : STATISTICS FOR ELECTRON ENERGY SOURCES
!PB            SIGMA_COP(2*NPLSI+2)       : STATISTICS FOR ION ENERGY SOURCES
cdr:  2015: all old preprogrammed sigma_cop removed.  COPV is a default tally,
cdr         and hence has its default variance options.
C
C
      SUBROUTINE EIRENE_STATIS_COP

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_CCONA
      USE EIRMOD_CGRID
      USE EIRMOD_CSDVI
      USE EIRMOD_CSDVI_COP
      USE EIRMOD_COUTAU

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: XN, FSIG, ZFLUX
      INTEGER, INTENT(IN) :: NBIN, NRIN, NPIN, NTIN, NSIN
      LOGICAL, INTENT(IN) :: LP, LT

      INTEGER, ALLOCATABLE :: IND(:,:), IIND(:), INDSS(:,:)
      REAL(DP) :: SD(0:NRTAL), SDD(0:NRTAL)
      REAL(DP) :: XNM, SD2, DS, ZFLUXQ, SD2S, SDS, SDI, SDE, D2S, SG,
     .          DSA, DD, D, SG2, DA, SD1, SD1S
      INTEGER :: IPLS, NSB, IR, ICO, IPL, NP2, NR1, ICPV, NT3, J, IIN,
     .           IRU, I
C
      SAVE
C
      ENTRY EIRENE_STATS0_COP
C
      IF (.NOT.ALLOCATED(IND)) THEN
        AllOCATE (IND(NRTAL,8))
        AllOCATE (IIND(NRTAL))
        AllOCATE (INDSS(NRTAL,8))
      END IF

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

      RETURN

C
      ENTRY EIRENE_STATS1_COP(NBIN,NRIN,NPIN,NTIN,NSIN,LP,LT)
      NSB=NBIN
      NR1=NRIN
      NP2=NPIN
      NT3=NTIN
C
C
      IF (NCPVI.EQ.0) RETURN
C
      IF (NCLMTS < NCLMT) NCLMTS = NCLMT
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
C
C
 1020 CONTINUE
      RETURN
C
      ENTRY EIRENE_STATS2_COP(XN,FSIG,ZFLUX)
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
      IF (NCPVI.EQ.0) GOTO 2200
C
 2200 CONTINUE
      RETURN
      END
