cdr  aug. 2016: added x coordinate, only for printing 1D profiles on file IFILE
cdr             to be done: exclude levgeo .gt. 3 from this 1D output format.

cdr  jan. 2017: synchronize with prtvol started (goal: remove special case prtvol)
cdr  june 2017: bug fix, printed tally output, calls to eirene_write_tally
cdr  jun. 2017: fully synchronized with prtvol, except printed text
C
C  INPUT:  T1,T2,T3:         TALLY TEXT, SPECIES AND UNITS, RESP.
C          PROF:             TALLY DATA, ON 1d ARRAX PROF(1:NRAD)
C          X:                X-COORDINATE: ONLY FOR 1D STANDARD GRIDS
C                            CURRENTLY: ZONE CENTERED
C          NR,NP,NT,NB,NTT:  GRID STRUCTURE FOR 2D OR 3D CASES
C          IFLAG:            SEE BELOW
C          IFILE:            WRITE OUT TALLY "PROF" ONTO STREAM: FORT.IFILE
C
      SUBROUTINE
     .  EIRENE_PRTTAL(T1,T2,T3,PROF,X,NR,NP,NT,NB,NTT,IFLAG,IFILE)
C
C  PRINT VOLUME AVERAGED TALLIES
C  IFLAG=-1:  ONLY HEADER IS PRINTED
C  IFLAG= 0:  ONLY MEAN VALUES IN EACH BLOCK
C  IFLAG= 1:  ADDITIONALLY: 1D PROFILES (THESE MAY BE AVERAGES)
C  IFLAG= 2:  ADDITIONALLY: 2D PROFILES (THESE MAY BE AVERAGES)
C  IFLAG= 3:  ADDITIONALLY: 3D PROFILES
C  IFLAG> 3:  ONLY FULL PROFILES, NO AVERAGES
C
C  IFILE> 0:  WRITE FULL TALLY ONTO STREAM FORT.IFILE
C
      USE EIRMOD_PRECISION
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: T1, T2, T3
      REAL(DP), INTENT(IN) :: PROF(*),X(*)
      INTEGER, INTENT(IN) :: NR, NP, NT, NB, NTT, IFLAG, IFILE
      INTEGER, PARAMETER :: NSTREAM=15
      REAL(DP) :: H(6)
      INTEGER :: K(6), ISTREAM(NSTREAM)
      INTEGER :: JR, JP, JT, IJ, N1DEL, N2DEL, IADD, JA, IA, IB, I,
     .           IC, IT, IP, NRM, NS, NTM, NPM, IRAD, IST, NCOL, IR,
     .           NTTS
      CHARACTER(1) :: TL(72)

      DATA TL/72*'='/
C  BLOCK A FEW RESERVED OUTPUT STREAMS.
      DATA ISTREAM/6,50,20,21,29,30,31,32,33,10,11,12,13,14,15/
      SAVE

      CALL EIRENE_LEER(3)
      WRITE (iunout,'(72A1)') TL
      WRITE (iunout,'(72A1)') TL
      WRITE (iunout,'(A9,A72)') 'TALLY:   ',T1
      WRITE (iunout,'(A9,A24)') 'SPECIES: ',T2
      WRITE (iunout,'(A9,A24)') 'UNITS:   ',T3
      WRITE (iunout,'(72A1)') TL
      WRITE (iunout,'(72A1)') TL
      CALL EIRENE_LEER(1)

C..................................................................
C  WRITE ONTO STREAM "IFILE"

      IF (IFILE.GT.0) THEN
        DO IST=1,NSTREAM
          IF (IFILE.EQ.ISTREAM(IST)) GOTO 11111
        ENDDO
        OPEN (UNIT=IFILE,POSITION='APPEND')
        WRITE (IFILE,*) TL
        WRITE (IFILE,*) TL
        WRITE (IFILE,*) T1
        WRITE (IFILE,*) T2
        WRITE (IFILE,*) T3
        WRITE (IFILE,*) TL
        WRITE (IFILE,*) TL
        WRITE (IFILE,*) NR,NP,NT,NB,NTT

        IF (NP.GT.1.OR.NT.GT.1.OR.NB.GT.1) THEN
C  THIS WAS A RUN WITH AN AT LEAST 2D GRID STRUCTURE
          DO IRAD=1,NTT,5
            WRITE (IFILE,*) (PROF(IR),IR=IRAD,MIN(IRAD+4,NTT))
          ENDDO

        ELSE
C  THIS WAS A 1D RUN, AT LEAST THERE IS ONLY A 1D GRID STRUCTURE
C  TO BE DONE: ALSO REQUIRED: LEVGEO.LE.3, OTHERWISE: TRIANGLES, TETRAHEDRA, HERE
          NTTS=NR-1
          DO IRAD=1,NTTS,1
            WRITE (IFILE,5) IRAD,X(IRAD),PROF(IRAD)
          ENDDO
C   NR: AVERAGED VALUE
          WRITE (IFILE,'(72A1)') TL
          WRITE (IFILE,55) PROF(NR)
          IF (NTT.GT.NR) THEN
C  ADDITIONAL CELL REGION
            WRITE (IFILE,'(72A1)') TL
            WRITE (IFILE,56)
            DO IRAD=NR+1,NTT
              WRITE (IFILE,57) IRAD-NR, PROF(IRAD)
            ENDDO
          ENDIF
        ENDIF
        WRITE (IFILE,'(72A1)') TL
        CLOSE (UNIT=IFILE)
      ENDIF
C
11111 IF (IFLAG.LT.0) RETURN

C...................................................................
C  WRITING ONTO STREAM "IFILE" COMPLETED FOR THIS TALLY

C  NCOL: NUMBER OF PRINTED DATA PER LINE, .LE.6
      NCOL=5
C
      NS=NR*NP*NT*NB
      NRM=MAX(1,NR-1)
      NPM=MAX(1,NP-1)
      NTM=MAX(1,NT-1)
      N1DEL=0
      IF (NP.GT.1.OR.NT.GT.1) N1DEL=NR
      N2DEL=0
      IF (NT.GT.1) N2DEL=NP
C
C LOOP OVER STANDARD MESH BLOCKS
C
      IF (NS.EQ.0) GOTO 50000
      DO 10000 IB=1,NB
      IF (NB.GT.1) THEN
        WRITE (iunout,*) TL
        WRITE (iunout,777) IB
        WRITE (iunout,*) TL
      ENDIF
      IADD=(IB-1)*NR*NP*NT
C
      IF (IFLAG.LE.2) GOTO 1000
      IF (NR.LE.1.OR.NP.LE.1.OR.NT.LE.1) GOTO 1000
C
C  3 D PROFILES
C
      WRITE (iunout,*) TL
      WRITE (iunout,81)
      WRITE (iunout,*) TL
      CALL EIRENE_LEER(1)
      DO 1 JT=1,NTM
        WRITE (iunout,77) JT
        CALL EIRENE_LEER(1)
        DO 11 JP=1,NPM
          IJ=1
          IR=0
          WRITE (iunout,7) JP
  110     DO 111 JR=IJ,NRM
            IC=JR+((JP-1)+(JT-1)*N2DEL)*N1DEL+IADD
            IR=IR+1
            IJ=IJ+1
            K(IR)=JR
            H(IR)=PROF(IC)
            IF (IR.GE.NCOL) GOTO 112
  111     CONTINUE
  112     CONTINUE
          IF (IR.GT.0) CALL EIRENE_WRITE_TALLY (K,H,IR,IUNOUT)
          IR=0
          IF (IJ.LE.NRM) GOTO 110
C  NEXT SEGMENT
          CALL EIRENE_LEER(2)
   11   CONTINUE
        WRITE (iunout,*) TL
    1 CONTINUE
      WRITE (iunout,*) TL
      IF (IFLAG.GT.3) GOTO 10000
C
C  2 D PROFILES
C
 1000 CONTINUE
C
      IF (IFLAG.LE.1) GOTO 2000
      IF (NR.LE.1.AND.NP.LE.1) GOTO 2000
      IF (NP.LE.1.AND.NT.LE.1) GOTO 2000
      IF (NR.LE.1.AND.NT.LE.1) GOTO 2000
C
C  POLOIDAL AND TOROIDAL PROFILE, RADIALLY AVERAGED
C
      IF (NTM.GT.1.AND.NPM.GT.1) THEN
        WRITE (iunout,82)
        IF (NR.GT.1) WRITE (iunout,881)
        IF (NR.GT.1) CALL EIRENE_LEER(1)
        DO 2 JT=1,NTM
          WRITE (iunout,77) JT
          IJ=1
          IP=0
  220     DO 222 JP=IJ,NPM
            IC=NR+((JP-1)+(JT-1)*N2DEL)*N1DEL+IADD
            IP=IP+1
            IJ=IJ+1
            K(IP)=JP
            H(IP)=PROF(IC)
            IF (IP.GE.NCOL) GOTO 223
  222     CONTINUE
  223     CONTINUE
          IF (IP.GT.0) CALL EIRENE_WRITE_TALLY (K,H,IP,IUNOUT)
C         WRITE (iunout,64) (K(I),H(I),I=1,IP)
          IP=0
          IF (IJ.LE.NPM) GOTO 220
C  NEXT SEGMENT
          CALL EIRENE_LEER(2)
    2   CONTINUE
        WRITE (iunout,*) TL
      ENDIF
C
C  RADIAL AND POLOIDAL PROFILE, TOROIDALLY AVERAGED
C
      IF (NRM.GT.1.AND.NPM.GT.1) THEN
        WRITE (iunout,81)
        IF (NT.GT.1) WRITE (iunout,883)
        IF (NT.GT.1) CALL EIRENE_LEER(1)
        DO 3 JP=1,NPM
          WRITE (iunout,7) JP
          IJ=1
          IR=0
  330     DO 333 JR=IJ,NRM
            IC=JR+((JP-1)+(NT-1)*N2DEL)*N1DEL+IADD
            IR=IR+1
            IJ=IJ+1
            K(IR)=JR
            H(IR)=PROF(IC)
            IF (IR.GE.NCOL) GOTO 334
  333     CONTINUE
  334     CONTINUE
          IF (IR.GT.0) CALL EIRENE_WRITE_TALLY (K,H,IR,IUNOUT)
C         WRITE (iunout,64) (K(I),H(I),I=1,IR)
          IR=0
          IF (IJ.LE.NRM) GOTO 330
C  NEXT SEGMENT
          CALL EIRENE_LEER(2)
    3   CONTINUE
        WRITE (iunout,*) TL
      ENDIF
C
C  RADIAL AND TOROIDAL PROFILE, POLOIDALLY AVERAGED
C
      IF (NRM.GT.1.AND.NTM.GT.1) THEN
        WRITE (iunout,81)
        IF (NP.GT.1) WRITE (iunout,882)
        IF (NP.GT.1) CALL EIRENE_LEER(1)
        DO 4 JT=1,NTM
          WRITE (iunout,77) JT
          IJ=1
          IR=0
  440     DO 444 JR=IJ,NRM
            IC=JR+((NP-1)+(JT-1)*N2DEL)*N1DEL+IADD
            IR=IR+1
            IJ=IJ+1
            K(IR)=JR
            H(IR)=PROF(IC)
            IF (IR.GE.NCOL) GOTO 445
  444     CONTINUE
  445     CONTINUE
          IF (IR.GT.0) CALL EIRENE_WRITE_TALLY (K,H,IR,IUNOUT)
C         WRITE (iunout,64) (K(I),H(I),I=1,IR)
          IR=0
          IF (IJ.LE.NRM) GOTO 440
C  NEXT SEGMENT
          CALL EIRENE_LEER(2)
    4   CONTINUE
        WRITE (iunout,*) TL
      ENDIF
      IF (IFLAG.GT.3) GOTO 10000
C
C  1 D PROFILES
C
 2000 CONTINUE
C
      IF (IFLAG.LE.0) GOTO 3000
C  RADIAL PROFILE, POLOIDALLY AND TOROIDALLY AVERAGED
C
      IF (NRM.GT.1) THEN
        WRITE (iunout,81)
        IF (NP.GT.1.AND.NT.EQ.1) WRITE (iunout,882)
        IF (NP.EQ.1.AND.NT.GT.1) WRITE (iunout,883)
        IF (NP.GT.1.AND.NT.GT.1) WRITE (iunout,8883)
        IJ=1
        IR=0
 1110   DO 1111 JR=IJ,NRM
          IC=JR+((NP-1)+(NT-1)*N2DEL)*N1DEL+IADD
          IR=IR+1
          IJ=IJ+1
          K(IR)=JR
          H(IR)=PROF(IC)
          IF (IR.GE.NCOL) GOTO 1112
 1111   CONTINUE
 1112   CONTINUE
        IF (IR.GT.0) CALL EIRENE_WRITE_TALLY (K,H,IR,IUNOUT)
C       WRITE (iunout,64) (K(I),H(I),I=1,IR)
        IR=0
        IF (IJ.LE.NRM) GOTO 1110
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) TL
      ENDIF
C
C  POLOIDAL PROFILE, RADIALLY AND TOROIDALLY AVERAGED
C
      IF (NPM.GT.1) THEN
        WRITE (iunout,82)
        IF (NR.GT.1.AND.NT.EQ.1) WRITE (iunout,881)
        IF (NR.EQ.1.AND.NT.GT.1) WRITE (iunout,883)
        IF (NR.GT.1.AND.NT.GT.1) WRITE (iunout,8882)
        IJ=1
        IP=0
 1220   DO 1222 JP=IJ,NPM
          IC=NR+((JP-1)+(NT-1)*N2DEL)*N1DEL+IADD
          IP=IP+1
          IJ=IJ+1
          K(IP)=JP
          H(IP)=PROF(IC)
          IF (IP.GE.NCOL) GOTO 1223
 1222   CONTINUE
 1223   CONTINUE
        IF (IP.GT.0) CALL EIRENE_WRITE_TALLY (K,H,IP,IUNOUT)
C       WRITE (iunout,64) (K(I),H(I),I=1,IP)
        IP=0
        IF (IJ.LE.NPM) GOTO 1220
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) TL
      ENDIF
C
C  TOROIDAL PROFILE, RADIALLY AND POLOIDALLY AVERAGED
C
      IF (NTM.GT.1) THEN
        IJ=1
        IT=0
        WRITE (iunout,83)
        IF (NR.GT.1.AND.NP.EQ.1) WRITE (iunout,881)
        IF (NR.EQ.1.AND.NP.GT.1) WRITE (iunout,882)
        IF (NR.GT.1.AND.NP.GT.1) WRITE (iunout,8881)
 1330   DO 1333 JT=IJ,NTM
          IC=NR+((NP-1)+(JT-1)*N2DEL)*N1DEL+IADD
          IT=IT+1
          IJ=IJ+1
          K(IT)=JT
          H(IT)=PROF(IC)
          IF (IT.GE.NCOL) GOTO 1334
 1333   CONTINUE
 1334   CONTINUE
        IF (IT.GT.0) CALL EIRENE_WRITE_TALLY (K,H,IT,IUNOUT)
C       WRITE (iunout,64) (K(I),H(I),I=1,IT)
        IT=0
        IF (IJ.LE.NTM) GOTO 1330
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) TL
      ENDIF
      IF (IFLAG.GT.3) GOTO 10000
C
 3000 CONTINUE
      IC=NR+((NP-1)+(NT-1)*N2DEL)*N1DEL+IADD
      WRITE (iunout,8888) PROF(IC)
      WRITE (iunout,*) TL
      CALL EIRENE_LEER(4)
C
10000 CONTINUE
C
50000 CONTINUE
C  ADDITIONAL CELLS
      IF (NTT.GT.NS) THEN
        WRITE (iunout,*) TL
        WRITE (iunout,7777)
        WRITE (iunout,*) TL
      ENDIF
      IJ=NS+1
      IA=0
  550 DO 555 JA=IJ,NTT
        IC=JA
        IA=IA+1
        IJ=IJ+1
        K(IA)=JA-NS
        H(IA)=PROF(IC)
        IF (IA.GE.NCOL) GOTO 556
  555 CONTINUE
  556 CONTINUE
      IF (IA.GT.0) CALL EIRENE_WRITE_TALLY (K,H,IA,IUNOUT)
C     WRITE (iunout,64) (K(IC),H(IC),IC=1,IA)
      IA=0
      IF (IJ.LE.NTT) GOTO 550
      CALL EIRENE_LEER(2)
C
    5 FORMAT (1X,I6,2X,2(1PE12.4,2X))
   55 FORMAT (1X,'AVERAGE VALUE ',1PE12.4)
   56 FORMAT (1X,'ADDITIONAL CELLS ')
   57 FORMAT (1X,I6,2X,1PE12.4)

    7 FORMAT (1X,'Y- OR POLOIDAL SEGMENT NUMBER ',I4)
   77 FORMAT (1X,'Z- OR TOROIDAL SEGMENT NUMBER ',I4)
  777 FORMAT (1X,'STANDARD MESH BLOCK NUMBER ',I4)
 7777 FORMAT (1X,'ADDITIONAL CELLS ')
   81 FORMAT (1X,'X- OR RADIAL PROFILE ')
   82 FORMAT (1X,'Y- OR POLOIDAL PROFILE ')
   83 FORMAT (1X,'Z- OR TOROIDAL PROFILE ')
  881 FORMAT (1X,'X- OR RADIAL AVERAGE ')
  882 FORMAT (1X,'Y- OR POLOIDAL AVERAGE ')
  883 FORMAT (1X,'Z- OR TOROIDAL AVERAGE ')
 8881 FORMAT (1X,'X- OR RAD. AND Y- OR POL. AVERAGE ',1PE12.4)
 8882 FORMAT (1X,'X- OR RAD. AND Z- OR TOR. AVERAGE ',1PE12.4)
 8883 FORMAT (1X,'Y- OR POL. AND Z- OR TOR. AVERAGE ',1PE12.4)
 8888 FORMAT (1X,'BLOCK AVERAGE ',1PE12.4)

      RETURN


      CONTAINS

      SUBROUTINE EIRENE_WRITE_TALLY(K,H,NR,IUNOUT)

c  write formatted line of tally values H, with their cell number K
c  1 .le.NR.le.6 is already verified in calling program.

      REAL(DP), INTENT(IN) :: H(6)
      INTEGER, INTENT(IN)  :: K(6),NR,IUNOUT

      IF (K(NR).LT.1E4) THEN
        WRITE (iunout,64) (K(I),H(I),I=1,NR)
      ELSEIF (K(NR).LT.1E5) THEN
        WRITE (iunout,65) (K(I),H(I),I=1,NR)
      ELSE
        WRITE (iunout,66) (K(I),H(I),I=1,NR)
      ENDIF

   64 FORMAT (1X,6(I4,2X,1PE12.4,2X))
   65 FORMAT (1X,6(I5,2X,1PE12.4,2X))
   66 FORMAT (1X,6(I6,2X,1PE12.4,2X))

      END SUBROUTINE EIRENE_WRITE_TALLY


      END
