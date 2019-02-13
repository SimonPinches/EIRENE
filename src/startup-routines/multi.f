CDR   June 17:  comments, and fix re option indpro=4 (unused so far)
C     May  05:  "no multip on averaging cells" corrected for 3D grids
C
      SUBROUTINE EIRENE_MULTI
cdr  entry multig:  copy grid data NBMLT times
cdr  entry multip:  indpro<=3: copy 1D profiles NP2ND*NT3RD*NBMLT times
cdr                 indpro>=4: copy    profiles            *NBMLT times
cdr                 indpro =4: check this option: tbd.

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CINIT
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM

      IMPLICIT NONE

      INTEGER :: I, J, K
C
C  GEOMETRY DATA
C
      ENTRY EIRENE_MULTIG
C
C  ZONE VOLUMES, KNOWN IN ZONE 1 TO NSTRD
      DO 130 J=2,NBMLT
        DO 120 I=1,NSTRD
          VOL(I+(J-1)*NSTRD)=VOL(I)*VOLCOR(J)
  120   CONTINUE
  130 CONTINUE
      DO 140 I=1,NSTRD
        VOL(I)=VOL(I)*VOLCOR(1)
  140 CONTINUE
C
      RETURN
C
C  PLASMA DATA
C
      ENTRY EIRENE_MULTIP
C
C  INDPRO.LT.4: ONLY RADIAL PLASMA PROFILES ARE GIVEN
C  RADIAL PLASMA PROFILES, KNOWN IN ZONES 1 TO NR1ST
C  NBLCKS=NP2ND*NT3RD*NBMLT
C  COPY THESE 1D (RADIAL) PROFILES, NBLCKS TIMES
c  EXCEPTION: B-FIELD DATA, INDPRO(5). SEE BELOW
C
      DO 210 J=2,NBLCKS
C  RADIAL "BLOCK" NO J
C  IS THIS A SPACE FOR AVERAGING: THEN DO NOT COPY
        IF (NP2ND.GT.1.AND.MOD(J,NP2ND).EQ.0) GOTO 210
        IF (NT3RD.GT.1.AND.NP2ND.LE.1.AND.MOD(J,NT3RD).EQ.0) GOTO 210
        IF (NT3RD.GT.1.AND.NP2ND.GT.1.AND.J.GT.NP2ND*(NT3RD-1)) GOTO 210
        IF (INDPRO(1).LT.4) THEN
          DO 201 I=1,NR1ST
            TEIN(I+(J-1)*NR1ST)=TEIN(I)
  201     CONTINUE
        ENDIF
        IF (INDPRO(2).LT.4) THEN
          DO 202 K=1,NPLSTI
            DO I=1,NR1ST
              TIIN(K,I+(J-1)*NR1ST)=TIIN(K,I)
            END DO
  202     CONTINUE
        ENDIF
        IF (INDPRO(3).LT.4) THEN
          DO 204 K=1,NPLSI
            DO I=1,NR1ST
              DIIN(K,I+(J-1)*NR1ST)=DIIN(K,I)
            END DO
  204     CONTINUE
        ENDIF
        IF (INDPRO(4).LT.4) THEN
          DO 205 K=1,NPLSV
            DO I=1,NR1ST
              VXIN(K,I+(J-1)*NR1ST)=VXIN(K,I)
              VYIN(K,I+(J-1)*NR1ST)=VYIN(K,I)
              VZIN(K,I+(J-1)*NR1ST)=VZIN(K,I)
            END DO
  205     CONTINUE
        ENDIF
        IF (INDPRO(5).LT.4) THEN
C  BFIELD DATA, INDPRO(5), ARE ALREADY SET ON 1:NSURF, SET IN PLASMA.F
        ENDIF
        IF (INDPRO(6).LT.4) THEN
          DO 207 K=1,NAINI
            DO I=1,NR1ST
              ADIN(K,I+(J-1)*NR1ST)=ADIN(K,I)
            END DO
  207     CONTINUE
        ENDIF
  210 CONTINUE
C
C  INDPRO.GT.4: ONLY NSTRD=NR1ST*NP2ND*NT3RD PLASMA DATA GIVEN
      DO 310 J=2,NBMLT
        IF (INDPRO(1).GT.4) THEN
          DO 301 I=1,NSTRD
            TEIN(I+(J-1)*NSTRD)=TEIN(I)
  301     CONTINUE
        ENDIF
        IF (INDPRO(2).GT.4) THEN
          DO 302 K=1,NPLSTI
            DO I=1,NSTRD
              TIIN(K,I+(J-1)*NSTRD)=TIIN(K,I)
            END DO
  302     CONTINUE
        ENDIF
        IF (INDPRO(3).GT.4) THEN
          DO 304 K=1,NPLSI
            DO I=1,NSTRD
              DIIN(K,I+(J-1)*NSTRD)=DIIN(K,I)
            END DO
  304     CONTINUE
        ENDIF
        IF (INDPRO(4).GT.4) THEN
          DO 305 K=1,NPLSV
            DO I=1,NSTRD
              VXIN(K,I+(J-1)*NSTRD)=VXIN(K,I)
              VYIN(K,I+(J-1)*NSTRD)=VYIN(K,I)
              VZIN(K,I+(J-1)*NSTRD)=VZIN(K,I)
            END DO
  305     CONTINUE
        ENDIF
        IF (INDPRO(5).GT.4) THEN
          DO 306 I=1,NSTRD
            BXIN(I+(J-1)*NSTRD)=BXIN(I)
            BYIN(I+(J-1)*NSTRD)=BYIN(I)
            BZIN(I+(J-1)*NSTRD)=BZIN(I)
            BFIN(I+(J-1)*NSTRD)=BFIN(I)
  306     CONTINUE
        ENDIF
        IF (INDPRO(6).GT.4) THEN
          DO 307 K=1,NAINI
            DO I=1,NSTRD
              ADIN(K,I+(J-1)*NSTRD)=ADIN(K,I)
            END DO
  307     CONTINUE
        ENDIF
  310 CONTINUE
C
      RETURN
      END
