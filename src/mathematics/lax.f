cdr  Old CR-model linear matrix equation solver: LAX, GALPD.
cdr  Single (inhomogeneous) right hand side vector only.
cdr  Added here because the internal H2-colrad that we are working on still uses
cdr  these versions.
cdr
cdr  March 2020:
cdr  LAX   --> LAX_M
cdr  GALPD --> GALPD_M
cdr  linear matrix equation solver A x= b, simultaneously for multiple right hand side terms b.
cdr  based on (and calling) routine GALPD.F: Gauss elimination, full pivoting.
cdr  Then: after transformation to upper triangular form: closed form solution.
cdr
cdr  This routine was/is used in CR codes: H_colrad, He_colrad, H2_colrad
cdr  to solve for the population coefficients.
cdr  So far: no benefit is taken from possible block triagonal structures.
cdr  See: Eigenvalue problems for special unsymmetric matrices.

cdr tbd: for CX coupling to cr models, we would need to solve two systems,
cdr      which differ only in the diagonal of the matrix.
cdr      Note:  if D is diagonal, then still A and D do not commute.
cdr      Note:  we do not need A**-1 nor (A+D)**-1 explicitly. If both A and A+D
cdr             could be brought to upper triangular form simultaneously, if would suffice.


      subroutine EIRENE_LAX_M(A,N1,N,B,nb,nbi,eps,ifl,is,vw,ip,icon)
C
C     COPY OF ORIGINAL ROUTINE EIRENE_LAX, extended FOR USE OF MULTIPLE RIGHT HAND SIDES
C     NEW ARGUMENT: NB:  STORAGE FOR NO. OF RIGHT HAND SIDE (INHOMOGENEOUS) VECTORS
C                        AS IN CALLING PROGRAM
C     NEW ARGUMENT: NBI: NO. OF RIGHT HAND SIDE (INHOMOGENEOUS) VECTORS. .LE.NB
C
      USE EIRMOD_PRECISION
      USE EIRMOD_COMPRT, ONLY: IUNOUT

      integer n1,n,nb,nbi,ifl,is,icon
      real(dp) a(n1,n1),B(nb,*),vw(*)
      real(dp) eps
      integer ip(*)
      integer iw(100), ier
      external eirene_exit_own, eirene_galpd_m

      if (n1.gt.100) then
        write (iunout,*) 'error in lax'
        call eirene_exit_own(1)
      endif
      call EIRENE_galpd_m(a,n1,n,b,nb,nbi,iw,ier)
      if (ier.eq.1) then
        write (iunout,*) 'error in lax, matrix is singular'
      endif
      return
      end subroutine EIRENE_LAX_M
c

c
      SUBROUTINE EIRENE_GALPD_M(A,NA,NG,B,NB,NBI,IW,IER)
C
C     COPY OF ORIGINAL ROUTINE EIRENE_GALPD FOR USE OF MULTIPLE RIGHT HAND SIDES
C     NEW ARGUMENT: nb: NO. OF RIGHT HAND SIDE (INHOMOGENEOUS) VECTORS
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
C    IW(NG)  : INTEGER HILFS-ARRAY FUER EINE MOEGLICHE PROGRAMM-
C              INTERNE UMNUMERIERUNG DER GLEICHUNGEN
C    IER     : ERROR-INDEX (IER = 1: MATRIX SINGULAER)
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER NA, NG, NB, NBI, IER
      INTEGER IW(NG)
      REAL(DP) :: A(NA,NA), B(NB,NG)
      REAL(DP) :: HB(NB), R(NB,NG)
      REAL(DP) :: AI, AK, AMN, AP, H, Q, ZERO
      INTEGER IB, I, K, M, N, IZ, KS, IH, II
      DATA ZERO /1.E-71_DP/
      IER=0
C
C     ******************************************************************
C     DER FALL: NG = 2
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
        IW(K)=K
    1 CONTINUE
C
C     ******************************************************************
C     DIE A-MATRIX AUF DREIECKS-FORM BRINGEN.
C     ******************************************************************
C
      DO 10 I=1,NG
C
C        ===============================================================
C        Pivot-Element suchen  ( Zeile IZ,  Spalte KS )
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
                        A(IZ,N)=H
    4                CONTINUE
                     HB(1:NBI)=B(1:NBI,I)
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
                        A(M,I)=AK
    5                CONTINUE
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
                                       A(M,N)=A(M,N)-A(I,N)*Q
    7                               CONTINUE
                                    B(1:NBI,M)=B(1:NBI,M)-B(1:NBI,I)*Q
                                    ENDIF
    8       CONTINUE
   10    CONTINUE
C
C     ******************************************************************
C     WENN A(N,N) ^= 0, KOENNEN DIE UNBEKANNTEN BERECHNET WERDEN.
C     -: SIE WERDEN ZUNAECHST AUF A(N,I), I=1,...,N GESETZT.
C     -: DANN DIE BERECHNETEN UNBEKANNTEN IN DER RICHTIGEN ANORDNUNG
C        AUF B(I) SCHREIBEN UND AN DAS AUFRUFENDE PROGRAMM ZURUECKGEBEN
C     ******************************************************************
C
      DO 12 M=1,NG
         R(1:NBI,M)=B(1:NBI,M)/A(M,M)
   12 CONTINUE
      DO 14 M=1,NG
         II=IW(M)
         B(1:NBI,II)=R(1:NBI,M)
   14 CONTINUE
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
      END SUBROUTINE EIRENE_GALPD_M

cdr  old versions of these routines, for only one right hand side vector
c
      subroutine LAX(A,N1,N,B,eps,ifl,is,vw,ip,icon)
      integer n1,n,iw,ifl,is,ip,icon,ier
      double precision a(n1,n1),B(*),vw(*),eps
      dimension ip(*)
      dimension iw(100)
      external galpd

      if (n1.gt.100) then
        write (6,*) 'error in lax'
      stop
      endif
      call galpd(a,n1,n,b,iw,ier)
      if (ier.eq.1) then
        write (6,*) 'error in lax, matrix is singular'
      endif
      return
      end subroutine lax
c

      SUBROUTINE GALPD(A,NA,NG,B,IW,IER)
C
C***********************************************************************
C*   GAUSS-ALGORITHMUS ZUR LOESUNG LINEARER GLEICHUNGS-SYSTEME MIT     *
C*   PIVOTIERUNG.                                                      *
C*   GENAUIGKEIT:   DOUBLE-PRECISION                   (01.07.1991)    *
C***********************************************************************
C    A(NA,NA): KOEFFIZIENTEN-MATRIX DES GLEICHUNGS-SYSTEMS
C    NA      : DIMENSION VON A WIE IM AUFRUFENDEN PROGRAMM ANGEGEBEN
C    NG      : ANZAHL DER UNBEKANNTEN   (NG <= NA)
C    B(NG)   : ELEMENTE DER RECHTEN SEITE DES GLEICHUNGS-SYSTEMS
c    B WIRD MODIFIZIERT UND ENTHAELT BEIM OUTPUT DIE NB LOESUNGSVEKTOREN
C    IW(NG)  : INTEGER-HILFS-ARRAY FUER EINE MOEGLICHE PROGRAMM-
C              INTERNE UMNUMERIERUNG DER GLEICHUNGEN
C    IER     : ERROR-INDEX (IER = 1: MATRIX SINGULAER)
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER NA,NG
      INTEGER IW(NG),IER
      DOUBLE PRECISION A(NA,NA),B(NG)
      INTEGER I,IH,II,IZ,K,KS,M,N
      DOUBLE PRECISION AI,AK,AMN,AP,H,Q,ZERO
      DATA ZERO /1.D-71/
      IER=0
C
C     ******************************************************************
C     DER FALL:   NG = 2
C     ******************************************************************
C
      IF(NG.EQ.2) THEN
                  H=A(1,1)*A(2,2)-A(2,1)*A(1,2)
                  AI=B(1)*A(2,2)-B(2)*A(1,2)
                  AK=A(1,1)*B(2)-A(2,1)*B(1)
                  B(1)=AI/H
                  B(2)=AK/H
                  RETURN
                  ENDIF
C
C     ******************************************************************
C     NUMMERN DER UNBEKANNTEN AUF IW ABSPEICHERN.
C     ******************************************************************
C
      DO K=1,NG
         IW(K)=K
      END DO
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
                     DO N=I,NG
                        H=A(I,N)
                        A(I,N)=A(IZ,N)
                        A(IZ,N)=H
                     END DO
                     H=B(I)
                     B(I)=B(IZ)
                     B(IZ)=H
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
                     DO M=1,NG
                        AK=A(M,KS)
                        A(M,KS)=A(M,I)
                        A(M,I)=AK
                     END DO
                     ENDIF
C
C        ===============================================================
C        Total-Pivotierung. Spalte i,  Zeile 1 ... i-1, i+1 ... NG
C        zu Null machen.
C        ===============================================================
C
         AP=A(I,I)
         IF(DABS(AP).LT.ZERO) GOTO 15
         AP=1/AP
         DO 8 M=1,NG

            IF(M.EQ.I) GOTO 8
cdr  schliesse diagonale aus. aber verwende A(I,I)
            IF(DABS(A(M,I)).GT.ZERO) THEN
                                     Q=A(M,I)*AP
                                     DO N=I,NG
                                        A(M,N)=A(M,N)-A(I,N)*Q
                                     END DO
                                     B(M)=B(M)-B(I)*Q
                                     ENDIF
    8       CONTINUE
   10    CONTINUE
C
C     ******************************************************************
C     WENN A(N,N) ^= 0, KOENNEN DIE UNBEKANNTEN BERECHNET WERDEN.
C     -:  SIE WERDEN ZUNAECHST AUF A(I,I), I=1,...,N GESETZT.
C     -:  DANN DIE BERECHNETEN UNBEKANNTEN IN DER RICHTIGEN ANORDNUNG
C         AUF B(I) SCHREIBEN UND AN DAS AUFRUFENDE PROGRAMM ZURUECKGEBEN
C     ******************************************************************
C
cdr  loesung auf A(M,M),
      DO M=1,NG
         A(M,M)=B(M)/A(M,M)
      END DO
cdr  wg spaltenvertauchung: ruecksortieren mit IW index array
      DO M=1,NG
         II=IW(M)
         B(II)=A(M,M)
      END DO
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
      END SUBROUTINE GALPD
