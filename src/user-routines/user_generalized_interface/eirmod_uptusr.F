! 20.11.09: use NSBOX as number of accounted cells instead of NTRII in order to
!           generalize the routine

      MODULE EIRMOD_UPTUSR

      PUBLIC

      CONTAINS

C
C
      SUBROUTINE EIRENE_UPTUSR(XSTOR2,XSTORV2,WV,IFLAG)
C
C  USER-SUPPLIED TRACKLENGTH ESTIMATOR, VOLUME-AVERAGED
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CESTIM
      USE EIRMOD_COMUSR
      USE EIRMOD_COMSOU
      USE EIRMOD_COMPRT
      USE EIRMOD_CUPD
      USE EIRMOD_COMXS
      USE EIRMOD_CSPEZ
      USE EIRMOD_CGRID
      USE EIRMOD_CLOGAU
      USE EIRMOD_CCONA
      USE EIRMOD_CPOLYG
      USE EIRMOD_CZT1
      USE EIRMOD_CTRIG
      USE EIRMOD_CGEOM

      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: XSTOR2(MSTOR1,MSTOR2,N2ND+N3RD),
     .                        XSTORV2(NSTORV,N2ND+N3RD), WV
      INTEGER, INTENT(IN) :: IFLAG
      REAL(DP), ALLOCATABLE :: PPPL_COP(:,:), MPPL_COP(:,:),
     .                         EPPL_COP(:,:), EPEL_COP(:)
      REAL(DP) :: RECTOT, SUMN, SUMM(0:2), SUMEI, SUMEE, RECADD, EEADD,
     .            EIRENE_FTABRC1, EIRENE_FEELRC1, PIADD_VEC, EIADD
      INTEGER :: IPL, JPLS, ISR, ISTEP, IIRC, IU, IPLSTI, IRRC, IN, INC,
     .           ICO, KK
      INTEGER, SAVE :: ISTROLD=-1


      IF (ISTRA /= ISTROLD) THEN
C
C
C  ADD CONTRIBUTIONS FROM VOLUME RECOMBINATION SOURCE
C
        ALLOCATE (PPPL_COP(NPLS,NRAD))
        ALLOCATE (MPPL_COP(NPLS,NRAD))
        ALLOCATE (EPPL_COP(NPLS,NRAD))
        ALLOCATE (EPEL_COP(NRAD))
        PPPL_COP =0.D0
        MPPL_COP = 0.D0
        EPPL_COP = 0.D0
        EPEL_COP = 0.D0

        IF (NLVOL(ISTRA)) THEN
C
          RECTOT = 0._DP
          DO 7473 IPLS=1,NPLSI
            JPLS = NSPEZ(ISTRA)
            IF ((JPLS > 0) .AND. (JPLS <= NPLSI) .AND.
     .          (IPLS /= JPLS)) CYCLE
            IPLSTI= MPLSTI(IPLS)
            DO ISR=1, NSRFSI(ISTRA)
              ISTEP = SORIND(ISR,ISTRA)
              DO 7472 IIRC=1,NPRCI(IPLS)
                IRRC=LGPRC(IPLS,IIRC)
                IF ((ISTEP > 0) .AND. (ISTEP /= IRRC)) CYCLE
                SUMN=0.0
                SUMM(0:2)=0.0
                SUMEI=0.0
                SUMEE=0.0
!pb                DO IN=1,NTRII
                DO IN=1,NSBOX
                  INC=NCLTAL(IN)
                  IF (NSTORDR >= NRAD) THEN
                    RECADD=-TABRC1(IRRC,IN)*DIIN(IPLS,IN)*ELCHA
                    EEADD=  EELRC1(IRRC,IN)*DIIN(IPLS,IN)*ELCHA
                  ELSE
                    RECADD=-EIRENE_FTABRC1(IRRC,IN)*DIIN(IPLS,IN)*ELCHA
                    EEADD=  EIRENE_FEELRC1(IRRC,IN)*DIIN(IPLS,IN)*ELCHA
                  END IF
                  PPPL_COP(IPLS,INC)=PPPL_COP(IPLS,INC)+RECADD
                  SUMN=SUMN+RECADD*VOL(IN)
                  PIADD_VEC(0:2)=0._DP
                  IF (LPMOM_VEC) THEN
                    DO ICO = 0,2
                      KK=ICO*NPLS
                      PIADD_VEC(ICO)=PMOM_VEC(KK+IPLS,IN)*RECADD
                      MPPL_COP(KK+IPLS,INC)=MPPL_COP(KK+IPLS,INC)+
     .                                      PIADD_VEC(ICO)
                      SUMM_VEC(ICO)=SUMM_VEC(ICO)+
     .                              PIADD_VEC(ICO)*VOL(IN)
                    END DO
                  END IF
                  EIADD=1.5*TIIN(IPLSTI,IN)*RECADD
                  IF (LEDRIFT) EIADD=EIADD+EDRIFT(IPLS,IN)*RECADD
                  EPPL_COP(IPLS,INC)=EPPL_COP(IPLS,INC)+EIADD
                  SUMEI=SUMEI+EIADD*VOL(IN)
                  EPEL_COP(INC)=EPEL_COP(INC)+EEADD
                  SUMEE=SUMEE+EEADD*VOL(IN)
                END DO
                RECTOT = RECTOT + SUMN
                WRITE (iunout,*) 'PARTIAL: IRRC ',IRRC
                CALL EIRENE_MASR6
     .           ('SUMN, SUMM_VEC, SUMEI, SUMEE                    ',
     .             SUMN, SUMM_VEC(0), SUMM_VEC(1), SUMM_VEC(2),
     .             SUMEI, SUMEE)
 7472         CONTINUE
            END DO
 7473     CONTINUE
        END IF

        IF (LCOPV) THEN
          IU = UBOUND(COPV,1)

          IF (IU >= NPLSI)
     .      COPV(1:NPLSI,:) = PPPL_COP(1:NPLSI,:)
          IF (IU >= 2*NPLSI)
     .      COPV(NPLSI+1:2*NPLSI,:) = MPPL_COP(1:NPLSI,:)
          IF (IU >= 3*NPLSI)
     .      COPV(2*NPLSI+1:3*NPLSI,:) = EPPL_COP(1:NPLSI,:)
          IF (IU >= 3*NPLSI+1)
     .      COPV(3*NPLSI+1,:) = EPEL_COP(:)
        END IF

        DEALLOCATE (PPPL_COP)
        DEALLOCATE (MPPL_COP)
        DEALLOCATE (EPPL_COP)
        DEALLOCATE (EPEL_COP)

        ISTROLD = ISTRA

      END IF

      RETURN
      END SUBROUTINE EIRENE_UPTUSR

      END MODULE EIRMOD_UPTUSR
