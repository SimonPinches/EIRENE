C  TO BE DONE: ADD PARAMETERS DUMMY AND ESBPARM, AS IN EIRENE_RDTRIM
C
C  read (try to read) IFLR = NHD6  (=12) pre-defined target projectile combinations.
C
C
      SUBROUTINE EIRENE_REFDAT(TMM,TCC,WMM,WCC)

c  input:  nhd6=12
c          iun =21+ifoff  (input stream for TRIM.dat)

C  This is the old (and default) model for reading TRIM conditional quantile tables, for reflection.

c  more recent input of TRIM files:  subr. RDTRIM: there: read individual trim files: A_on_B
c  as selected in input block 6.
C
C  THIS SUBROUTINE READS REFLECTION DATA PRODUCED BY BCA MONTE CARLO CODES,
C  DATA ARE STORED IN CONDITIONAL QUANTILE FORMAT (E.G. TRIM)
C
C  distinct from rdtrim.f this routines reads one single file containing many
C  (in this present version: NHD6=12)
C  fixed target-projectile combinations, hard-wired in the following order
C
C    IFILE=1  H ON FE
C    IFILE=2  D ON FE
C    IFILE=3  H ON C
C    IFILE=4  D ON C
C    IFILE=5  HE ON FE
C    IFILE=6  HE ON C
C    IFILE=7  T ON FE
C    IFILE=8  T ON C
C  W targets added at some later stage to the trim.dat file.
C  nhd6 increased from 8 to 12
C    IFILE=9  D ON W
C    IFILE=10 HE ON W
C    IFILE=11 H ON W
C    IFILE=12 T ON W
C
C
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CREF
      USE EIRMOD_CSPEI
      USE EIRMOD_CINIT

      IMPLICIT NONE
C
      REAL(DP), INTENT(OUT) :: TMM(*), TCC(*), WMM(*), WCC(*)
      REAL(DP) :: TML(12), TCL(12), WML(12), WCL(12),
     .          FELD(1092)
      REAL(DP) :: PID180
      INTEGER :: I1, I2, I3, I4, I5, NRECL, IUN, IFILE, I, J, IFLR
C
      IFLR=NHD6
      IF (NHD6.GT.12) THEN
        WRITE (iunout,*) 'STORAGE ERROR. NHD6 MUST BE LESS OR EQUAL 12 '
        WRITE (iunout,*) 'NHD6= ',NHD6
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
C  HARD-WIRED FORMAT:  12*7*5*5*5,  AND NHD6=12 such files.
      NRECL=1092

      TML(1)=1._DP
      TCL(1)=1._DP
      WML(1)=56._DP
      WCL(1)=26._DP
      TML(2)=2._DP
      TCL(2)=1._DP
      WML(2)=56._DP
      WCL(2)=26._DP
      TML(3)=1._DP
      TCL(3)=1._DP
      WML(3)=12._DP
      WCL(3)=6._DP
      TML(4)=2._DP
      TCL(4)=1._DP
      WML(4)=12._DP
      WCL(4)=6._DP
      TML(5)=4._DP
      TCL(5)=2._DP
      WML(5)=56._DP
      WCL(5)=26._DP
      TML(6)=4._DP
      TCL(6)=2._DP
      WML(6)=12._DP
      WCL(6)=6._DP
      TML(7)=3._DP
      TCL(7)=1._DP
      WML(7)=56._DP
      WCL(7)=26._DP
      TML(8)=3._DP
      TCL(8)=1._DP
      WML(8)=12._DP
      WCL(8)=6._DP
      TML(9)=2._DP
      TCL(9)=1._DP
      WML(9)=184._DP
      WCL(9)=74._DP
      TML(10)=4._DP
      TCL(10)=2._DP
      WML(10)=184._DP
      WCL(10)=74._DP
      TML(11)=1._DP
      TCL(11)=1._DP
      WML(11)=184._DP
      WCL(11)=74._DP
      TML(12)=3._DP
      TCL(12)=1._DP
      WML(12)=184._DP
      WCL(12)=74._DP

      DO 10 J=1,NHD6
        TMM(J)=TML(J)
        TCC(J)=TCL(J)
        WMM(J)=WML(J)
        WCC(J)=WCL(J)
   10 CONTINUE

c  INE=12 incident energies
c  (not to be confused with the nhd6=12 projectile-target cases)
      INE=12
      INEM=INE-1

      ENAR(1)=1._DP
      ENAR(2)=2._DP
      ENAR(3)=5._DP
      ENAR(4)=10._DP
      ENAR(5)=20._DP
      ENAR(6)=50._DP
      ENAR(7)=100._DP
      ENAR(8)=200._DP
      ENAR(9)=500._DP
      ENAR(10)=1000._DP
      ENAR(11)=2000._DP
      ENAR(12)=5000._DP
      DO 11 I=1,INEM
        DENAR(I)=1./(ENAR(I+1)-ENAR(I))
   11 CONTINUE

      INW=7
      INWM=INW-1
      PID180=ATAN(1._DP)/45._DP
      WIAR(1)=COS(0._DP)
      WIAR(2)=COS(30._DP*PID180)
      WIAR(3)=COS(45._DP*PID180)
      WIAR(4)=COS(60._DP*PID180)
      WIAR(5)=COS(70._DP*PID180)
      WIAR(6)=COS(80._DP*PID180)
      WIAR(7)=COS(85._DP*PID180)
      DO 12 I=1,INWM
        DWIAR(I)=1._DP/(WIAR(I+1)-WIAR(I))
   12 CONTINUE
C
      INR=5
      INRM=INR-1
      RAAR(1)=0.1_DP
      RAAR(2)=0.3_DP
      RAAR(3)=0.5_DP
      RAAR(4)=0.7_DP
      RAAR(5)=0.9_DP
      DO 15 I=1,INRM
        DRAAR(I)=1./(RAAR(I+1)-RAAR(I))
   15 CONTINUE
C
      IF (INE*INW*IFLR.GT.NH0 .OR.
     .    INE*INW*INR*IFLR.GT.NH1  .OR.
     .    INE*INW*INR*INR*IFLR.GT.NH2  .OR.
     .    INE*INW*INR*INR*INR*IFLR.GT.NH3) THEN
        WRITE (iunout,*)
     .    'ERROR IN PARAMETER STATEMENT FOR REFLECTION DATA'
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
C
      DO IFILE=1, NDBNAMES
        IF (INDEX(DBHANDLE(IFILE),'TRIM') /= 0) EXIT
      END DO

      IF (IFILE > NDBNAMES) THEN
        WRITE (IUNOUT,*) ' NO DATABASE NAME FOR TRIM DEFINED '
        WRITE (IUNOUT,*) ' CALCULATION ABANDONED '
        CALL EIRENE_EXIT_OWN(1)
      END IF

      IUN=21+ifoff
      OPEN (UNIT=IUN,FILE=DBFNAME(IFILE))
      REWIND IUN
C

C  IFLR=12, HARD-WIRED
      DO 1 IFILE=1,IFLR
        DO 2 I1=1,INE
          I=0
!read nrecl=1092 data blocks, one for each of the INE=12 incident energies
!  7 + 7*5 + 7*5*5 + 7*5*5+5 = 1092
          READ (IUN,661) (FELD(J),J=1,NRECL)
c  7: reflection coefficients for each of the 7 incident angles
          DO I2=1,INW
            I=I+1
            HFTR0(I1,I2,IFILE)=FELD(I)
          END DO
c  7*5: energy quantiles
          DO I2=1,INW
            DO I3=1,INR
              I=I+1
              HFTR1(I1,I2,I3,IFILE)=FELD(I)
            END DO
          END DO
c  7*5*5:  polar angle quantiles
          DO I2=1,INW
            DO I3=1,INR
              DO I4=1,INR
                I=I+1
                HFTR2(I1,I2,I3,I4,IFILE)=FELD(I)
              END DO
            END DO
          END DO
c  7*5*5*5:  azimuthal angle quantiles
          DO I2=1,INW
            DO I3=1,INR
              DO I4=1,INR
                DO I5=1,INR
                  I=I+1
                  HFTR3(I1,I2,I3,I4,I5,IFILE)=FELD(I)
                END DO
              END DO
            END DO
          END DO
    2   CONTINUE
    1 CONTINUE
C

  661 FORMAT (4E20.12)


      RETURN
      END
