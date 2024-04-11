      SUBROUTINE EIRENE_SCAL_SURF_TALLIES (ISTR)
cdr May 19:    scoring of surface-averaged spectra:
cdr          part. balance scaling with FATM, FMOL,... is done in SCALE_TALLIES
cdr Mar 23:  remove ADDS tally from universal scaling with FLXFAC,
cdr          and apply the input flags IADSE from block 10D instead
cdr          Not ready, but ZWW, ZW are already here?

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CESTIM
      USE EIRMOD_COUTAU
      USE EIRMOD_CCONA

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ISTR
      INTEGER :: I, IATM, IMOL, IION, IPHOT, IPLS, ISPC, D
      REAL(DP) :: DEL, DELI, ELEFT, ERIGHT

cdr  surface flux tallies
      ESTIMS=ESTIMS*FLXFAC(ISTR)

C
C  these next tallies are not surface tallies.
C  We scale them here anyway.
C
      ETOTA(ISTR)=ETOTA(ISTR)*FLXFAC(ISTR)
      ETOTM(ISTR)=ETOTM(ISTR)*FLXFAC(ISTR)
      ETOTI(ISTR)=ETOTI(ISTR)*FLXFAC(ISTR)
      ETOTP(ISTR)=ETOTP(ISTR)*FLXFAC(ISTR)
      ETOTPH(ISTR)=ETOTPH(ISTR)*FLXFAC(ISTR)
      PTRASH(ISTR)=PTRASH(ISTR)*FLXFAC(ISTR)
      ETRASH(ISTR)=ETRASH(ISTR)*FLXFAC(ISTR)
      DO 610 IATM=0,NATM
        WTOTA(IATM,ISTR)=WTOTA(IATM,ISTR)*FLXFAC(ISTR)
  610 CONTINUE
      DO 615 IMOL=0,NMOL
        WTOTM(IMOL,ISTR)=WTOTM(IMOL,ISTR)*FLXFAC(ISTR)
  615 CONTINUE
      DO 620 IION=0,NION
        EELFI(IION,ISTR)=EELFI(IION,ISTR)*FLXFAC(ISTR)
        WTOTI(IION,ISTR)=WTOTI(IION,ISTR)*FLXFAC(ISTR)
  620 CONTINUE
      DO 622 IPHOT=0,NPHOT
        WTOTPH(IPHOT,ISTR)=WTOTPH(IPHOT,ISTR)*FLXFAC(ISTR)
  622 CONTINUE
      DO 625 IPLS=0,NPLS
        WTOTP(IPLS,ISTR)=WTOTP(IPLS,ISTR)*FLXFAC(ISTR)
  625 CONTINUE
      WTOTE(ISTR)=WTOTE(ISTR)*FLXFAC(ISTR)

C   SCALE AND INTEGRATE SURFACE SPECTRA

      DO ISPC=1,NADSPC
        IF (ESTIML(ISPC)%ISRFCLL == 0) THEN
C  SURFACE FLUX SPECTRA
          ESTIML(ISPC)%SPCS = 0._DP
          DO I = 0, ESTIML(ISPC)%NSPC+1
            IF (I.EQ.0) THEN

              ELEFT = MIN(ESTIML(ISPC)%SPCMIN,
     .                   ESTIML(ISPC)%ESP_MIN)
              ERIGHT= ESTIML(ISPC)%SPCMIN
            ELSE IF (I.EQ.ESTIML(ISPC)%NSPC+1) THEN
              ELEFT = ESTIML(ISPC)%SPCMAX
              ERIGHT= MAX(ESTIML(ISPC)%SPCMAX,
     .                   ESTIML(ISPC)%ESP_MAX)
            ELSE
              ELEFT = ESTIML(ISPC)%SPCMIN+
     .               ESTIML(ISPC)%SPCDEL*(I-1)
              ERIGHT= ESTIML(ISPC)%SPCMIN+
     .               ESTIML(ISPC)%SPCDEL*I
            END IF
            IF (ESTIML(ISPC)%LOG) THEN
              DEL = 10._DP**ERIGHT-10._DP**ELEFT+EPS60
            ELSE
              DEL = ERIGHT-ELEFT+EPS60
            END IF
cdr 1/DE, DE= energy increment for bin no. I.
            DELI = 1._DP/(DEL+EPS60)

C  SCALE: FROM SCORING TALLY UNITS PER ENERGY BIN --> TALLY UNITS PER EV
            ESTIML(ISPC)%SPC(I) =
     .      ESTIML(ISPC)%SPC(I)*FLXFAC(ISTR)*DELI
chk Also scale the Legendre expansion coefficients if applicable
            IF (ESTIML(ISPC)%ISPCOPT==2) THEN
              DO D=1,ESTIML(ISPC)%ISPLDEG
                ESTIML(ISPC)%SPCAN(D,I) =
     .            ESTIML(ISPC)%SPCAN(D,I)*FLXFAC(ISTR)*DELI
              END DO
            END IF
C  INTEGRATE  --> TALLY UNITS
cdr  Test tbd: in case of total (not directional) spectrum,
cdr            i.e. for IDIREC=0, this
cdr            integral must coincide with the particle outflux POT..
cdr            or the energy outflux EOT..,
cdr            at the selected surface, depending on ISPTYP=1,
cdr            or =2, respectively.
cdr            I believe this test must also work
cdr            in the same way for directional spectra.
cdr            The spectral resolution is along a line of sight,
cdr            but cumulative wrt. to the orthogonal
cdr            directions.
            ESTIML(ISPC)%SPCS = ESTIML(ISPC)%SPCS +
     .                          ESTIML(ISPC)%SPC(I)*DEL
          END DO
        END IF
      END DO

      RETURN

      END SUBROUTINE EIRENE_SCAL_SURF_TALLIES
