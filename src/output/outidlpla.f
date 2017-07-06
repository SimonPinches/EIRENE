cdr  tally 22 (electric potential) added, and a few comments, started...
C
      SUBROUTINE EIRENE_OUTIDLPLA
C
C  This routine prints background tallies for plotting in IDL tool.
C
C  PRINT INPUT TALLIES ONTO OUTPUT FILE IUNOUT
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_COMUSR
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CGRID
      USE EIRMOD_CSPEZ
      USE EIRMOD_CTRCEI
      USE EIRMOD_CGEOM
      USE EIRMOD_CTEXT
      USE EIRMOD_COUTAU
      USE EIRMOD_CSPEI
      USE EIRMOD_CINIT
 
      IMPLICIT NONE
 
      REAL(DP), ALLOCATABLE :: HELPP(:,:),HELPW(:,:),TALAV(:),TALTOT(:)
      REAL(DP) :: TALTYP(NTALI)
      REAL(DP) :: HELPI, TOTAL
      INTEGER :: IR, IP, IT, I, NBLCKA, IB, ITAL, NXM, NYM, NZM,
     .           K, NFTI, NFTE, MXSPZ, IOUT
      CHARACTER(50) :: FNAME, FORMA, FORME, FORME2
C
C  TYPE OF TALLY: TALTYP=0: #              (#-UNITS)
C                 TALTYP=1: #-DENSITY      (#-UNITS/CM**3)
C                 TALTYP=2: VOLUME         (CM**3)
C                 TALTYP=3: DIMENSIONLESS  (1)
C                 TALTYP=4: UNKNOWN        (?)
c
c  number of volumetric input tallies: ntali = 22
      TALTYP(1)=0
      TALTYP(2)=0
      TALTYP(3)=1
      TALTYP(4)=1
      TALTYP(5)=0
      TALTYP(6)=0
      TALTYP(7)=0
      TALTYP(8)=3
      TALTYP(9)=3
      TALTYP(10)=3
      TALTYP(11)=3
      TALTYP(12)=4
      TALTYP(13)=0
      TALTYP(14)=2
      TALTYP(15)=3
      TALTYP(16)=3
      TALTYP(17)=3
      TALTYP(18)=0
      TALTYP(19)=0
      TALTYP(20)=0
      TALTYP(21)=0
      TALTYP(22)=0

      MXSPZ = MAXVAL(NFSTPI(1:NTALI))
      
      ALLOCATE (HELPP(NRAD,MXSPZ))
      ALLOCATE (HELPW(NRAD,MXSPZ))
      ALLOCATE (TALTOT(MXSPZ))
      ALLOCATE (TALAV(MXSPZ))

      IOUT = 3 + IFOFF

      FORMA=REPEAT(' ',50)
      FORMA='(6X,   A25)'
      WRITE (FORMA(5:7),'(I3)') MXSPZ
      
      FORME=REPEAT(' ',50)
      FORME='(I10,   ES25.7)'
      WRITE (FORME(6:8),'(I3)') MXSPZ 
      
      FORME2=REPEAT(' ',50)
      FORME2='(6X,   ES25.7)'
      WRITE (FORME2(5:7),'(I3)') MXSPZ 
C     
C  PRINT INPUT VOLUME AVERAGED TALLIES
C
      NXM=MAX(1,NR1STM)
      NYM=MAX(1,NP2NDM)
      NZM=MAX(1,NT3RDM)

      DO 100 ITAL=1,NTALI

        NFTI=1
        NFTE=NFSTPI(ITAL)

        DO 119 K=NFTI,NFTE

          SELECT CASE (ITAL)
          CASE (1)
            HELPP(1:NSBOX,K) = TEIN(1:NSBOX)
          CASE (2)
            HELPP(1:NSBOX,K) = TIIN(MPLSTI(K),1:NSBOX)
          CASE (3)
            HELPP(1:NSBOX,K) = DEIN(1:NSBOX)
          CASE (4)
            HELPP(1:NSBOX,K) = DIIN(K,1:NSBOX)
          CASE (5)
            HELPP(1:NSBOX,K) = VXIN(MPLSV(K),1:NSBOX)
          CASE (6)
            HELPP(1:NSBOX,K) = VYIN(MPLSV(K),1:NSBOX)
          CASE (7)
            HELPP(1:NSBOX,K) = VZIN(MPLSV(K),1:NSBOX)
          CASE (8)
            HELPP(1:NSBOX,K) = BXIN(1:NSBOX)
          CASE (9)
            HELPP(1:NSBOX,K) = BYIN(1:NSBOX)
          CASE (10)
            HELPP(1:NSBOX,K) = BZIN(1:NSBOX)
          CASE (11)
            HELPP(1:NSBOX,K) = BFIN(1:NSBOX)
          CASE (12)
            HELPP(1:NSBOX,K) = ADIN(K,1:NSBOX)
          CASE (13)
            HELPP(1:NSBOX,K) = EDRIFT(K,1:NSBOX)
          CASE (14)
            HELPP(1:NSBOX,K) = VOL(1:NSBOX)
          CASE (15)
            HELPP(1:NSBOX,K) = WGHT(K,1:NSBOX)
          CASE (16)
            HELPP(1:NSBOX,K) = BXPERP(1:NSBOX)
          CASE (17)
            HELPP(1:NSBOX,K) = BYPERP(1:NSBOX)
          CASE (18)
            HELPP(1:NSBOX,K) = EXIN(1:NSBOX)
          CASE (19)
            HELPP(1:NSBOX,K) = EYIN(1:NSBOX)
          CASE (20)
            HELPP(1:NSBOX,K) = EZIN(1:NSBOX)
          CASE (21)
            HELPP(1:NSBOX,K) = EFIN(1:NSBOX)
          CASE (22)
            HELPP(1:NSBOX,K) = POT(1:NSBOX)
          CASE DEFAULT
            WRITE (iunout,*) ' WRONG TALLY NUMBER (OUTIDLPLA), ITAL = '
     .                        ,ITAL
            WRITE (iunout,*) ' NO OUTPUT PERFORMED '
            CALL EIRENE_LEER(1)
            GOTO 100
          END SELECT
C
          TOTAL=0.D0
          DO 121 IB=1,NBMLT
          NBLCKA=NSTRD*(IB-1)
          DO 121 IR=1,NXM
          DO 121 IP=1,NYM
          DO 121 IT=1,NZM
            I=IR + ((IP-1)+(IT-1)*NP2T3)*NR1P2 + NBLCKA
            IF (ITAL.EQ.1) THEN
C  ELECTR. TEMPERATURE: NE*VOLUME WEIGHTED AVERAGES
              HELPW(I,K)=DEIN(I)*VOL(I)
            ELSEIF (ITAL.EQ.2) THEN
C  ION TEMPERTURE: NI(K)*VOLUME WEIGHTED AVERAGES
              HELPW(I,K)=DIIN(K,I)*VOL(I)
            ELSEIF (ITAL.EQ.3.OR.ITAL.EQ.4) THEN
C  PARTICLE DENSITY PROFILES: VOLUME WEIGHTED AVERAGES
              HELPW(I,K)=VOL(I)
            ELSEIF (ITAL.EQ.5.OR.ITAL.EQ.6.OR.ITAL.EQ.7) THEN
C  ION DRIFT VELOCITY: NI(K)*VOLUME WEIGHTED AVERAGES
              HELPW(I,K)=DIIN(K,I)*VOL(I)
            ELSEIF (ITAL.GE.8.AND.ITAL.LE.11) THEN
C  B-FIELD UNIT VECTOR, B-FIELD STRENGTH   " 1 - WEIGHTED" AVERAGES
              HELPW(I,K)=1.D0
              IF (NSTGRD(I).GT.0) HELPW(I,K)=0.D0
            ELSEIF (ITAL.EQ.16.OR.ITAL.LE.17) THEN
C  B_PERP-FIELD: " 1 - WEIGHTED" AVERAGES
              HELPW(I,K)=1.D0
              IF (NSTGRD(I).GT.0) HELPW(I,K)=0.D0
            ELSEIF (ITAL.EQ.12.OR.ITAL.EQ.14.OR.ITAL.EQ.15) THEN
C  ADDITIONAL TALLY, CELL VOLUME ,WEIGHT FUNCTION " 1 - WEIGHTED" AVERAGES
              HELPW(I,K)=1.D0
              IF (NSTGRD(I).GT.0) HELPW(I,K)=0.D0
            ELSEIF (ITAL.EQ.13) THEN
C  ION DRIFT ENERGY
              HELPW(I,K)=DIIN(K,I)*VOL(I)
            ELSEIF (ITAL.GE.18.AND.ITAL.LE.21) THEN
C  E-FIELD UNIT VECTOR, E-FIELD STRENGTH   
              HELPW(I,K)=1.D0
            ELSEIF (ITAL.EQ.22) THEN
C  ELECTRIC POTENTIAL   
              HELPW(I,K)=1.D0
            ENDIF
            TOTAL=TOTAL+HELPW(I,K)
121       CONTINUE
C
C  SAME LOOP AGAIN, OVER ADDITIONAL CELL REGION
          DO 122 I=NSURF+1,NSURF+NRADD
            IF (ITAL.EQ.1) THEN
C  ELECTR. TEMPERATURE: NE*VOLUME WEIGHTED AVERAGES
              HELPW(I,K)=DEIN(I)*VOL(I)
            ELSEIF (ITAL.EQ.2) THEN
C  ION TEMPERTURE: NI(K)*VOLUME WEIGHTED AVERAGES
              HELPW(I,K)=DIIN(K,I)*VOL(I)
            ELSEIF (ITAL.EQ.3.OR.ITAL.EQ.4) THEN
C  PARTICLE DENSITY PROFILES: VOLUME WEIGHTED AVERAGES
              HELPW(I,K)=VOL(I)
            ELSEIF (ITAL.EQ.5.OR.ITAL.EQ.6.OR.ITAL.EQ.7) THEN
C  ION DRIFT VELOCITY: NI(K)*VOLUME WEIGHTED AVERAGES
              HELPW(I,K)=DIIN(K,I)*VOL(I)
            ELSEIF (ITAL.GE.8.AND.ITAL.LE.11) THEN
C  B-FIELD UNIT VECTOR, B-FIELD STRENGTH   " 1 - WEIGHTED" AVERAGES
              HELPW(I,K)=1.D0
              IF (NSTGRD(I).GT.0) HELPW(I,K)=0.D0
            ELSEIF (ITAL.EQ.16.OR.ITAL.LE.17) THEN
C  B_PERP-FIELD: " 1 - WEIGHTED" AVERAGES
              HELPW(I,K)=1.D0
              IF (NSTGRD(I).GT.0) HELPW(I,K)=0.D0
            ELSEIF (ITAL.EQ.12.OR.ITAL.EQ.14.OR.ITAL.EQ.15) THEN
C  ADDITIONAL TALLY (NO.12)
C  CELL VOLUME      (NO.14)
C  WEIGHT FUNCTION  (NO.15)
              HELPW(I,K)=1.D0
              IF (NSTGRD(I).GT.0) HELPW(I,K)=0.D0
            ELSEIF (ITAL.EQ.13) THEN
C  ION DRIFT ENERGY: NI(K)*VOLUME WEIGHTED AVERAGES
              HELPW(I,K)=DIIN(K,I)*VOL(I)
            ELSEIF (ITAL.GE.18.AND.ITAL.LE.21) THEN
C  E-FIELD UNIT VECTOR, E-FIELD STRENGTH   
              HELPW(I,K)=1.D0
            ELSEIF (ITAL.EQ.22) THEN
C  ELECTRIC POTENTIAL   
              HELPW(I,K)=1.D0
            ENDIF
            TOTAL=TOTAL+HELPW(I,K)
122       CONTINUE
C
          IF (ITAL.NE.NTALO) THEN
            CALL EIRENE_INTTAL (HELPP(1,K),HELPW(1,K),1,1,NSBOX,HELPI,
     .                          NR1ST,NP2ND,NT3RD,NBMLT)
          ELSEIF (ITAL.EQ.NTALO) THEN
            CALL EIRENE_INTVOL (HELPP(1,K),      1,1,NSBOX,HELPI,
     .                          NR1ST,NP2ND,NT3RD,NBMLT)
          ENDIF

          TALTOT(K)=HELPI
          TALAV(K)=HELPI/(TOTAL+EPS60)

119     CONTINUE
C
        FNAME = 'intal_  '
        WRITE (FNAME(7:8),'(i0)') ITAL

        OPEN (UNIT=IOUT,FILE=FNAME,FORM='FORMATTED',
     .        ACCESS='SEQUENTIAL')

        WRITE (IOUT,'(A)') TXTPLS(1,ITAL)
        WRITE (IOUT,'(A,I10)') 'NCELLS:   ',NSBOX
        WRITE (IOUT,'(A,I10)') 'NSPECIES: ',NFSTPI(ITAL)
        WRITE (IOUT,'(A)') 'SPECIES'
        WRITE (IOUT,FORMA) (TRIM(TXTPSP(K,ITAL)), K=NFTI, NFTE)
        WRITE (IOUT,'(A)') 'UNITS'
        WRITE (IOUT,FORMA) (TRIM(TXTPUN(K,ITAL)), K=NFTI, NFTE)

        IF (TALTYP(ITAL).EQ.0) THEN
          WRITE (IOUT,'(A)')
     .        'WEIGHTED MEAN VALUE ("UNITS")               '
          WRITE (IOUT,'(A)') 'MEAN'
          WRITE (IOUT,FORME2) (TALAV(K), K=NFTI, NFTE)
        ELSEIF (TALTYP(ITAL).EQ.1) THEN
          WRITE (IOUT,'(A)')
     .        'WEIGHTED TOTAL ("UNITS*CM**3"), MEAN ("UNITS")'
          WRITE (IOUT,'(A)') 'TOTAL'
          WRITE (IOUT,FORME2) (TALTOT(K), K=NFTI, NFTE)
          WRITE (IOUT,'(A)') 'MEAN'
          WRITE (IOUT,FORME2) (TALAV(K), K=NFTI, NFTE)
        ELSEIF (TALTYP(ITAL).EQ.2) THEN
          WRITE (IOUT,'(A)')
     .        'ARITHMETIC TOTAL ("UNITS")                         '
          WRITE (IOUT,'(A)') 'TOTAL'
          WRITE (IOUT,FORME2) (TALTOT(K), K=NFTI, NFTE)
        ELSEIF (TALTYP(ITAL).EQ.3) THEN
          WRITE (IOUT,'(A)')
     .        'ARITHMETIC MEAN ("UNITS")                          '
          WRITE (IOUT,'(A)') 'MEAN'
          WRITE (IOUT,FORME2) (TALAV(K), K=NFTI, NFTE)
        ELSEIF (TALTYP(ITAL).EQ.4) THEN
          WRITE (IOUT,'(A)')
     .        'ARITHMETIC TOTAL ("UNITS"), MEAN ("UNITS")'
          WRITE (IOUT,'(A)') 'TOTAL'
          WRITE (IOUT,FORME2) (TALTOT(K), K=NFTI, NFTE)
          WRITE (IOUT,'(A)') 'MEAN'
          WRITE (IOUT,FORME2) (TALAV(K), K=NFTI, NFTE)
        ENDIF

        WRITE (IOUT,'(A)') '==========================================='
          DO I=1, NSBOX
            IF (ANY(ABS(HELPP(I,NFTI:NFTE)) > EPS30)) 
     .        WRITE (IOUT,FORME) I,(HELPP(I,K), K=NFTI, NFTE) 
          END DO
        WRITE (IOUT,'(A)') '==========================================='

        CLOSE (UNIT=IOUT)
C
100   CONTINUE
 
      DEALLOCATE (HELPP)
      DEALLOCATE (HELPW)
      DEALLOCATE (TALTOT)
      DEALLOCATE (TALAV)
C
      RETURN
      END
