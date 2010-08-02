C
C
      SUBROUTINE EIRENE_SYMET(ESTIM,NTAL,NRAD,NR1ST,NP2ND,NT3RD,
     .                 LP,LT)
C
C  SYMMETRISE TALLIES
C
      USE EIRMOD_PRECISION
      IMPLICIT NONE
 
      REAL(DP), INTENT(INOUT) :: ESTIM(NTAL,NRAD)
      INTEGER, INTENT(IN) :: NTAL, NRAD, NR1ST, NP2ND, NT3RD
      LOGICAL, INTENT(IN) :: LP,LT
      REAL(DP) :: SAV
      INTEGER :: IR, IP, IT, IND, I1, INDEX1, INDEX2, J1, J2, ITAL,
     .           NSYM, NSYH
 
      IF (LP) THEN
        NSYM=NP2ND
        NSYH=(NSYM-1)/2
        DO 30 ITAL=1,NTAL
            DO 10 IR=1,NR1ST
            DO 10 IT=1,NT3RD
            DO 10 IP=1,NSYH
                  J1=IR+((IT-1)*NP2ND+IP-1)*NR1ST
                  J2=IR+((IT-1)*NP2ND+NSYM-IP-1)*NR1ST
                  SAV=(ESTIM(ITAL,J1)+ESTIM(ITAL,J2))*0.5
                  ESTIM(ITAL,J1)=SAV
                  ESTIM(ITAL,J2)=SAV
10          CONTINUE
30      CONTINUE
      ENDIF
      IF (LT) THEN
        NSYM=NT3RD
        NSYH=(NSYM-1)/2
        DO 130 ITAL=1,NTAL
            DO 110 IR=1,NR1ST
            DO 110 IP=1,NP2ND
            DO 110 IT=1,NSYH
                  J1=IR+((IT-1)*NP2ND+IP-1)*NR1ST
                  J2=IR+((NSYM-IT-1)*NP2ND+IP-1)*NR1ST
                  SAV=(ESTIM(ITAL,J1)+ESTIM(ITAL,J2))*0.5
                  ESTIM(ITAL,J1)=SAV
                  ESTIM(ITAL,J2)=SAV
110         CONTINUE
130     CONTINUE
      ENDIF
      RETURN
      END
