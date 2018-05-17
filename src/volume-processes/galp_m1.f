      SUBROUTINE EIRENE_GALPD_M(A,NA,NG,B,NB,NBI,IW,IER)
C
C     COPY OF ORIGINAL ROUTINE EIRENE_GALPD FOR USE OF MULTIPLE RIGHT HAND SIDES
C     NEW ARGUMENT: nb:  NO. OF RIGHT HAND SIDE (INHOMOGENEOUS) VECTORS
C
      USE EIRMOD_PRECISION
C
C***********************************************************************
C*   GAUSS-ALGORITHMUS ZUR LOESUNG LINEARER GLEICHUNGS-SYSTEME MIT     *
C*   PIVOTIERUNG.                                                      *
C*   GENAUIGKEIT:   DOUBLE-PRECISION                   (01.07.1991)    *
C***********************************************************************
C    A(NA,NA): KOEFFIZIENTEN-MATRIX DES GLEICHUNGS-SYSTEMS
C    NA      : DIMENSION VON A WIE IM AUFRUFENDEN PROGRAMM ANGEGEBEN
C    NG      : ANZAHL DER UNBEKANNTEN   (NG <= NA)
C    B(NB,NG): ELEMENTE DER RECHTEN SEITE DES GLEICHUNGS-SYSTEMS
C    NB      : DIMENSION ANZAHL DER RECHTEN SEITEN (NB >= 1)
C    NBI     : ANZAHL DER RECHTEN SEITEN (NB >= 1)
c    B WIRD MODIFIZIERT UND ENTHAELT BEIM OUTPUT DIE NBI LOESUNGSVEKTOREN
C    IW(NG)  : INTEGER-HILFS-ARRAY FUER EINE MOEGLICHE PROGRAMM-
C              INTERNE UMNUMERIERUNG DER GLEICHUNGEN
C    IER     : ERROR-INDEX (IER = 1: MATRIX SINGULAER)
C***********************************************************************
C
      IMPLICIT REAL(DP) (A-H,O-Z)
      DIMENSION A(NA,NA),B(NB,NG),IW(NG)
      DIMENSION HB(NB), R(NB,NG)
      DATA ZERO /1.E-71_DP/
      IER=0
C
C     ******************************************************************
C     DER FALL:   NG = 2
C     ******************************************************************
C
      IF(NG.EQ.2) THEN
                  DO IB = 1, NBI
                    H=A(1,1)*A(2,2)-A(2,1)*A(1,2)
                    AI=B(IB,1)*A(2,2)-B(IB,2)*A(1,2)
                    AK=A(1,1)*B(IB,2)-A(2,1)*B(IB,1)
                    B(IB,1)=AI/H
                    B(IB,2)=AK/H
                  END DO
                  RETURN
                  ENDIF
C
C     ******************************************************************
C     NUMMERN DER UNBEKANNTEN AUF IW ABSPEICHERN.
C     ******************************************************************
C
      DO 1 K=1,NG
    1    IW(K)=K
C
C     ******************************************************************
C     DIE A-MATRIX AUF DREIECKS-FORM BRINGEN.
C     ******************************************************************
C
      DO 10 I=1,NG
C
C        ===============================================================
C        Pivot-Element suchen  (  Zeile IZ,  Spalte KS )
C        ===============================================================
C
         AP=0
         IZ=0
         KS=0
         DO 3 M=I,NG
            DO 2 N=I,NG
               AMN=ABS(A(M,N))
               IF(AMN.GT.AP) THEN
                             AP=AMN
 
                             IZ=M
 
                             KS=N
 
                             ENDIF
    2          CONTINUE
    3      CONTINUE
C
C        ============================
C        Zeilen umordnen, wenn IZ > I
C        ============================
C
         IF(IZ.GT.I) THEN
                     DO 4 N=I,NG
                        H=A(I,N)
                        A(I,N)=A(IZ,N)
    4                   A(IZ,N)=H
                     HB=B(1:NBI,I)
                     B(1:NBI,I)=B(1:NBI,IZ)
                     B(1:NBI,IZ)=HB(1:NBI)
                     ENDIF
C
C        ===============================================
C        Spalten umordnen und Unbekannte neu numerieren
C        ===============================================
C
         IF(KS.GT.I) THEN
                     IH=IW(I)
                     IW(I)=IW(KS)
                     IW(KS)=IH
                     DO 5 M=1,NG
                        AK=A(M,KS)
                        A(M,KS)=A(M,I)
    5                   A(M,I)=AK
                     ENDIF
C
C        ===============================================================
C        Total-Pivotierung. Spalte i,  Zeile 1 ... i-1, i+1 ... NG
C        zu Null machen.
C        ===============================================================
C
         AP=A(I,I)
         IF(ABS(AP).LT.ZERO) GOTO 15
         AP=1/AP
         DO 8 M=1,NG
            IF(M.EQ.I) GOTO 8
            IF(ABS(A(M,I)).GT.ZERO) THEN
                                    Q=A(M,I)*AP
                                    DO 7 N=I,NG
    7                                  A(M,N)=A(M,N)-A(I,N)*Q
                                    B(1:NBI,M)=B(1:NBI,M)-B(1:NBI,I)*Q
                                    ENDIF
    8       CONTINUE
   10    CONTINUE
C
C     ******************************************************************
C     WENN A(N,N) ^= 0, KOENNEN DIE UNBEKANNTEN BERECHNET WERDEN.
C     -:  SIE WERDEN ZUNAECHST AUF A(N,I), I=1,...,N GESETZT.
C     -:  DANN DIE BERECHNETEN UNBEKANNTEN IN DER RICHTIGEN ANORDNUNG
C         AUF B(I) SCHREIBEN UND AN DAS AUFRUFENDE PROGRAMM ZURUECKGEBEN
C     ******************************************************************
C
      DO 12 M=1,NG
   12    R(1:NBI,M)=B(1:NBI,M)/A(M,M)
      DO 14 M=1,NG
         II=IW(M)
   14    B(1:NBI,II)=R(1:NBI,M)
C
       RETURN
 
C     ******************************************************************
C     MATRIX IST SINGULAER.
C       -: IER = 1 SETZTEN
C       -: RETURN
C     ******************************************************************
C
   15 IER=1
      RETURN
      END
