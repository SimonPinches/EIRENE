cdr jan. 2018: added: outgoing current tallies: ind=1
cdr sept.2014: only comments....

c  SCORE "INCIDENT" and "EMITTED" SURFACE FLUX TALLIES, FOR SURFACE MSURF, OR SURFACE SEGMENT MSURFG.
c  INCIDENT: update tallies POT_A_B(iout,msurf) and EOT_A_B(iout,msurf)
c  EMITTED : update tallies PRF_A_B(iout,msurf) and ERF_A_B(iout,msurf)
c  A code-letter for incident type of particle: A, M, I, P, PH
c  B code-letter for emitted type of particle :  AT, ML, IO, PL, PHT
c  iout:  species index for emitted particle

      SUBROUTINE EIRENE_UPDATE_SURFACE (ITOLD,WGHTSG,IND)

C       PARTICLE FLUXES (WEIGHT=WGHTSG),     score POT.., PRF...
C       ENERGY FLUXES   (E0*WGHTSG),         score EOT.., ERF...

c  input:
c
c  ind=1:  score incident currents
c    itold:  type of incident particle
c    ispez  (iatm, imol, iion, ipls, iphot): of incident particle
c    msurf:  surface index
c    msurfg:  sub-segment of surface MSURF, for spatial resolution on surface
c    E0:     energy (eV) of incident particle
c    WGHTSG:   stat. weight of incident particle WEIGHT* PROB* SIGN.

c  ind=2:  score reemitted currents
c    ityp :  type of emitted particle
c    ispez  (iatm, imol, iion, ipls, iphot): of emitted particle
c    msurf:  surface index
c    msurfg:  sub-segment of surface MSURF, for spatial resolution on surface
c    E0:     energy (eV) of reemitted particle
c    WGHTSG:   stat. weight of reemitted particle* PROB* SIGN.

c  output:
c  lmetspw(ispz):  species ispz is emitted, emitted flux tally is scored.

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CESTIM
      USE EIRMOD_COMPRT
      USE EIRMOD_COMUSR
      USE EIRMOD_CSPEZ
      USE EIRMOD_CSDVI

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ITOLD, IND
      REAL(dp), INTENT(IN) :: WGHTSG
      REAL(dp) ::             EWGHTSG

      IF (MSURF .LE. 0) RETURN

      EWGHTSG=E0*WGHTSG

      select case(IND)

      CASE(1)
C  incident fluxes
      IF (ITOLD.EQ.0) THEN
C  INCIDENT PHOTONS
        IF (LEOTPHT) THEN
!$OMP ATOMIC
          EOTPHT(IPHOT,MSURF)=EOTPHT(IPHOT,MSURF)+EWGHTSG
        ENDIF
        IF (LPOTPHT) THEN
!$OMP ATOMIC
          POTPHT(IPHOT,MSURF)=POTPHT(IPHOT,MSURF)+WGHTSG
        ENDIF

        IF (MSURFG.GT.0) THEN
          IF (LEOTPHT) THEN
!$OMP ATOMIC
            EOTPHT(IPHOT,MSURFG)=EOTPHT(IPHOT,MSURFG)+EWGHTSG
          ENDIF
          IF (LPOTPHT) THEN
!$OMP ATOMIC
            POTPHT(IPHOT,MSURFG)=POTPHT(IPHOT,MSURFG)+WGHTSG
          ENDIF
        ENDIF
        IF (LEOTPHT .OR. LPOTPHT) LMETSPW(IPHOT) = .TRUE.
      ELSEIF (ITOLD.EQ.1) THEN
C  INCIDENT ATOMS
        IF (LEOTAT) THEN
!$OMP ATOMIC
          EOTAT(IATM,MSURF)=EOTAT(IATM,MSURF)+EWGHTSG
        ENDIF
        IF (LPOTAT) THEN
!$OMP ATOMIC
          POTAT(IATM,MSURF)=POTAT(IATM,MSURF)+WGHTSG
        ENDIF

        IF (MSURFG.GT.0) THEN
          IF (LEOTAT) THEN
!$OMP ATOMIC
            EOTAT(IATM,MSURFG)=EOTAT(IATM,MSURFG)+EWGHTSG
          ENDIF
          IF (LPOTAT) THEN
!$OMP ATOMIC
            POTAT(IATM,MSURFG)=POTAT(IATM,MSURFG)+WGHTSG
          ENDIF
        ENDIF
        IF (LEOTAT .OR. LPOTAT) LMETSPW(NSPH+IATM) = .TRUE.
      ELSEIF (ITOLD.EQ.2) THEN
C  INCIDENT MOLECULES
        IF (LEOTML) THEN
!$OMP ATOMIC
          EOTML(IMOL,MSURF)=EOTML(IMOL,MSURF)+EWGHTSG
        ENDIF
        IF (LPOTML) THEN
!$OMP ATOMIC
          POTML(IMOL,MSURF)=POTML(IMOL,MSURF)+WGHTSG
        ENDIF

        IF (MSURFG.GT.0) THEN
          IF (LEOTML) THEN
!$OMP ATOMIC
             EOTML(IMOL,MSURFG)=EOTML(IMOL,MSURFG)+EWGHTSG
          ENDIF
          IF (LPOTML) THEN
!$OMP ATOMIC
             POTML(IMOL,MSURFG)=POTML(IMOL,MSURFG)+WGHTSG
          ENDIF
        ENDIF
        IF (LEOTML .OR. LPOTML) LMETSPW(NSPA+IMOL) = .TRUE.
      ELSEIF (ITOLD.EQ.3) THEN
C  INCIDENT TEST IONS
        IF (LEOTIO) THEN
!$OMP ATOMIC
          EOTIO(IION,MSURF)=EOTIO(IION,MSURF)+EWGHTSG
        ENDIF
        IF (LPOTIO) THEN
!$OMP ATOMIC
          POTIO(IION,MSURF)=POTIO(IION,MSURF)+WGHTSG
        ENDIF

        IF (MSURFG.GT.0) THEN
          IF (LEOTIO) THEN
!$OMP ATOMIC
            EOTIO(IION,MSURFG)=EOTIO(IION,MSURFG)+EWGHTSG
          ENDIF
          IF (LPOTIO) THEN
!$OMP ATOMIC
            POTIO(IION,MSURFG)=POTIO(IION,MSURFG)+WGHTSG
          ENDIF
        ENDIF
        IF (LEOTIO .OR. LPOTIO) LMETSPW(NSPAM+IION) = .TRUE.
      ELSEIF (ITOLD.EQ.4) THEN
C  INCIDENT BULK IONS
        IF (LEOTPL) THEN
!$OMP ATOMIC
          EOTPL(IPLS,MSURF)=EOTPL(IPLS,MSURF)+EWGHTSG
        ENDIF
        IF (LPOTPL) THEN
!$OMP ATOMIC
          POTPL(IPLS,MSURF)=POTPL(IPLS,MSURF)+WGHTSG
        ENDIF

        IF (MSURFG.GT.0) THEN
          IF (LEOTPL) THEN
!$OMP ATOMIC
            EOTPL(IPLS,MSURFG)=EOTPL(IPLS,MSURFG)+EWGHTSG
          ENDIF
          IF (LPOTPL) THEN
!$OMP ATOMIC
            POTPL(IPLS,MSURFG)=POTPL(IPLS,MSURFG)+WGHTSG
          ENDIF
        ENDIF
        IF (LEOTPL .OR. LPOTPL) LMETSPW(NSPAMI+IPLS) = .TRUE.
      ELSE
        goto 999
      ENDIF  !ITOLD


      CASE(2)
C  reemitted fluxes
c  a photon is reemitted. currently only foreseen for incident photons
      IF (ITYP.EQ.0) THEN
        LOGPHOT(IPHOT,ISTRA)=.TRUE.
        IF (ITOLD.EQ.0) THEN
c... from an incident photon
          IF (LPRFPHPHT) THEN
!$OMP ATOMIC
            PRFPHPHT(IPHOT,MSURF)=PRFPHPHT(IPHOT,MSURF)+WGHTSG
          ENDIF
          IF (LERFPHPHT) THEN
!$OMP ATOMIC
            ERFPHPHT(IPHOT,MSURF)=ERFPHPHT(IPHOT,MSURF)+EWGHTSG
          ENDIF
          IF (MSURFG.GT.0) THEN
            IF (LPRFPHPHT) THEN
!$OMP ATOMIC
              PRFPHPHT(IPHOT,MSURFG)=PRFPHPHT(IPHOT,MSURFG)+WGHTSG
            ENDIF
            IF (LERFPHPHT) THEN
!$OMP ATOMIC
              ERFPHPHT(IPHOT,MSURFG)=ERFPHPHT(IPHOT,MSURFG)+EWGHTSG
            ENDIF
          ENDIF
          IF (LPRFPHPHT .OR. LERFPHPHT) LMETSPW(IPHOT) = .TRUE.
        ELSE
          goto 999
        ENDIF

c  an atom is reemitted...
      ELSEIF (ITYP.EQ.1) THEN
        LOGATM(IATM,ISTRA)=.TRUE.
        IF (ITOLD.EQ.1) THEN
c... from an incident atom
          IF (LPRFAAT) THEN
!$OMP ATOMIC
            PRFAAT(IATM,MSURF)=PRFAAT(IATM,MSURF)+WGHTSG
          ENDIF
          IF (LERFAAT) THEN
!$OMP ATOMIC
            ERFAAT(IATM,MSURF)=ERFAAT(IATM,MSURF)+EWGHTSG
          ENDIF
          IF (MSURFG.GT.0) THEN
            IF (LPRFAAT) THEN
!$OMP ATOMIC
              PRFAAT(IATM,MSURFG)=PRFAAT(IATM,MSURFG)+WGHTSG
            ENDIF
            IF (LERFAAT) THEN
!$OMP ATOMIC
              ERFAAT(IATM,MSURFG)=ERFAAT(IATM,MSURFG)+EWGHTSG
            ENDIF
          ENDIF
          IF (LPRFAAT .OR. LERFAAT) LMETSPW(NSPH+IATM) = .TRUE.
        ELSEIF (ITOLD.EQ.2) THEN
c... from an incident molecule
          IF (LPRFMAT) THEN
!$OMP ATOMIC
            PRFMAT(IATM,MSURF)=PRFMAT(IATM,MSURF)+WGHTSG
          ENDIF
          IF (LERFMAT) THEN
!$OMP ATOMIC
            ERFMAT(IATM,MSURF)=ERFMAT(IATM,MSURF)+EWGHTSG
          ENDIF
          IF (MSURFG.GT.0) THEN
            IF (LPRFMAT) THEN
!$OMP ATOMIC
              PRFMAT(IATM,MSURFG)=PRFMAT(IATM,MSURFG)+WGHTSG
            ENDIF
            IF (LERFMAT) THEN
!$OMP ATOMIC
              ERFMAT(IATM,MSURFG)=ERFMAT(IATM,MSURFG)+EWGHTSG
            ENDIF
          ENDIF
          IF (LPRFMAT .OR. LERFMAT) LMETSPW(NSPH+IATM) = .TRUE.
        ELSEIF (ITOLD.EQ.3) THEN
c... from an incident test ion
          IF (LPRFIAT) THEN
!$OMP ATOMIC
            PRFIAT(IATM,MSURF)=PRFIAT(IATM,MSURF)+WGHTSG
          ENDIF
          IF (LERFIAT) THEN
!$OMP ATOMIC
            ERFIAT(IATM,MSURF)=ERFIAT(IATM,MSURF)+EWGHTSG
          ENDIF
          IF (MSURFG.GT.0) THEN
            IF (LPRFIAT) THEN
!$OMP ATOMIC
              PRFIAT(IATM,MSURFG)=PRFIAT(IATM,MSURFG)+WGHTSG
            ENDIF
            IF (LERFIAT) THEN
!$OMP ATOMIC
              ERFIAT(IATM,MSURFG)=ERFIAT(IATM,MSURFG)+EWGHTSG
            ENDIF
          ENDIF
          IF (LPRFIAT .OR. LERFIAT) LMETSPW(NSPH+IATM) = .TRUE.
        ELSEIF (ITOLD.EQ.4) THEN
c... from an incident bulk ion
          IF (LPRFPAT) THEN
!$OMP ATOMIC
            PRFPAT(IATM,MSURF)=PRFPAT(IATM,MSURF)+WGHTSG
          ENDIF
          IF (LERFPAT) THEN
!$OMP ATOMIC
            ERFPAT(IATM,MSURF)=ERFPAT(IATM,MSURF)+EWGHTSG
          ENDIF
          IF (MSURFG.GT.0) THEN
            IF (LPRFPAT) THEN
!$OMP ATOMIC
              PRFPAT(IATM,MSURFG)=PRFPAT(IATM,MSURFG)+WGHTSG
            ENDIF
            IF (LERFPAT) THEN
!$OMP ATOMIC
              ERFPAT(IATM,MSURFG)=ERFPAT(IATM,MSURFG)+EWGHTSG
            ENDIF
          ENDIF
          IF (LPRFPAT .OR. LERFPAT) LMETSPW(NSPH+IATM) = .TRUE.
        ELSE
          goto 999
        ENDIF

c  a molecule is reemitted
      ELSEIF (ITYP.EQ.2) THEN
        LOGMOL(IMOL,ISTRA)=.TRUE.
        IF (ITOLD.EQ.1) THEN
c... from an incident atom
          IF (LPRFAML) THEN
!$OMP ATOMIC
             PRFAML(IMOL,MSURF)=PRFAML(IMOL,MSURF)+WGHTSG
          ENDIF
          IF (LERFAML) THEN
!$OMP ATOMIC
            ERFAML(IMOL,MSURF)=ERFAML(IMOL,MSURF)+EWGHTSG
          ENDIF
          IF (MSURFG.GT.0) THEN
            IF (LPRFAML) THEN
!$OMP ATOMIC
              PRFAML(IMOL,MSURFG)=PRFAML(IMOL,MSURFG)+WGHTSG
            ENDIF
            IF (LERFAML) THEN
!$OMP ATOMIC
              ERFAML(IMOL,MSURFG)=ERFAML(IMOL,MSURFG)+EWGHTSG
            ENDIF
          ENDIF
          IF (LPRFAML .OR. LERFAML) LMETSPW(NSPA+IMOL) = .TRUE.
        ELSEIF (ITOLD.EQ.2) THEN
c... from an incident molecule
          IF (LPRFMML) THEN
!$OMP ATOMIC
            PRFMML(IMOL,MSURF)=PRFMML(IMOL,MSURF)+WGHTSG
          ENDIF
          IF (LERFMML) THEN
!$OMP ATOMIC
            ERFMML(IMOL,MSURF)=ERFMML(IMOL,MSURF)+EWGHTSG
          ENDIF
          IF (MSURFG.GT.0) THEN
            IF (LPRFMML) THEN
!$OMP ATOMIC
              PRFMML(IMOL,MSURFG)=PRFMML(IMOL,MSURFG)+WGHTSG
            ENDIF
            IF (LERFMML) THEN
!$OMP ATOMIC
              ERFMML(IMOL,MSURFG)=ERFMML(IMOL,MSURFG)+EWGHTSG
            ENDIF
          ENDIF
          IF (LPRFMML .OR. LERFMML) LMETSPW(NSPA+IMOL) = .TRUE.
        ELSEIF (ITOLD.EQ.3) THEN
c... from an incident test ion
          IF (LPRFIML) THEN
!$OMP ATOMIC
            PRFIML(IMOL,MSURF)=PRFIML(IMOL,MSURF)+WGHTSG
          ENDIF
          IF (LERFIML) THEN
!$OMP ATOMIC
            ERFIML(IMOL,MSURF)=ERFIML(IMOL,MSURF)+EWGHTSG
          ENDIF
          IF (MSURFG.GT.0) THEN
            IF (LPRFIML) THEN
!$OMP ATOMIC
              PRFIML(IMOL,MSURFG)=PRFIML(IMOL,MSURFG)+WGHTSG
            ENDIF
            IF (LERFIML) THEN
!$OMP ATOMIC
              ERFIML(IMOL,MSURFG)=ERFIML(IMOL,MSURFG)+EWGHTSG
            ENDIF
          ENDIF
          IF (LPRFIML .OR. LERFIML) LMETSPW(NSPA+IMOL) = .TRUE.
        ELSEIF (ITOLD.EQ.4) THEN
c... from an incident bulk ion
          IF (LPRFPML) THEN
!$OMP ATOMIC
            PRFPML(IMOL,MSURF)=PRFPML(IMOL,MSURF)+WGHTSG
          ENDIF
          IF (LERFPML) THEN
!$OMP ATOMIC
            ERFPML(IMOL,MSURF)=ERFPML(IMOL,MSURF)+EWGHTSG
          ENDIF
          IF (MSURFG.GT.0) THEN
            IF (LPRFPML) THEN
!$OMP ATOMIC
              PRFPML(IMOL,MSURFG)=PRFPML(IMOL,MSURFG)+WGHTSG
            ENDIF
            IF (LERFPML) THEN
!$OMP ATOMIC
              ERFPML(IMOL,MSURFG)=ERFPML(IMOL,MSURFG)+EWGHTSG
            ENDIF
          ENDIF
          IF (LPRFPML .OR. LERFPML) LMETSPW(NSPA+IMOL) = .TRUE.
        ELSE
          goto 999
        ENDIF

c  a test ion is reemitted
      ELSEIF (ITYP.EQ.3) THEN
        LOGION(IION,ISTRA)=.TRUE.
        IF (ITOLD.EQ.1) THEN
          IF (LPRFAIO) THEN
!$OMP ATOMIC
            PRFAIO(IION,MSURF)=PRFAIO(IION,MSURF)+WGHTSG
          ENDIF
          IF (LERFAIO) THEN
!$OMP ATOMIC
            ERFAIO(IION,MSURF)=ERFAIO(IION,MSURF)+EWGHTSG
          ENDIF
          IF (MSURFG.GT.0) THEN
            IF (LPRFAIO) THEN
!$OMP ATOMIC
              PRFAIO(IION,MSURFG)=PRFAIO(IION,MSURFG)+WGHTSG
            ENDIF
            IF (LERFAIO) THEN
!$OMP ATOMIC
              ERFAIO(IION,MSURFG)=ERFAIO(IION,MSURFG)+EWGHTSG
            ENDIF
          ENDIF
          IF (LPRFAIO .OR. LERFAIO) LMETSPW(NSPAM+IION) = .TRUE.
        ELSEIF (ITOLD.EQ.2) THEN
          IF (LPRFMIO) THEN
!$OMP ATOMIC
            PRFMIO(IION,MSURF)=PRFMIO(IION,MSURF)+WGHTSG
          ENDIF
          IF (LERFMIO) THEN
!$OMP ATOMIC
            ERFMIO(IION,MSURF)=ERFMIO(IION,MSURF)+EWGHTSG
          ENDIF
          IF (MSURFG.GT.0) THEN
            IF (LPRFMIO) THEN
!$OMP ATOMIC
              PRFMIO(IION,MSURFG)=PRFMIO(IION,MSURFG)+WGHTSG
            ENDIF
            IF (LERFMIO) THEN
!$OMP ATOMIC
              ERFMIO(IION,MSURFG)=ERFMIO(IION,MSURFG)+EWGHTSG
            ENDIF
          ENDIF
          IF (LPRFMIO .OR. LERFMIO) LMETSPW(NSPAM+IION) = .TRUE.
        ELSEIF (ITOLD.EQ.3) THEN
          IF (LPRFIIO) THEN
!$OMP ATOMIC
            PRFIIO(IION,MSURF)=PRFIIO(IION,MSURF)+WGHTSG
          ENDIF
          IF (LERFIIO) THEN
!$OMP ATOMIC
            ERFIIO(IION,MSURF)=ERFIIO(IION,MSURF)+EWGHTSG
          ENDIF
          IF (MSURFG.GT.0) THEN
            IF (LPRFIIO) THEN
!$OMP ATOMIC
              PRFIIO(IION,MSURFG)=PRFIIO(IION,MSURFG)+WGHTSG
            ENDIF
            IF (LERFIIO) THEN
!$OMP ATOMIC
              ERFIIO(IION,MSURFG)=ERFIIO(IION,MSURFG)+EWGHTSG
            ENDIF
          ENDIF
          IF (LPRFIIO .OR. LERFIIO) LMETSPW(NSPAM+IION) = .TRUE.
        ELSEIF (ITOLD.EQ.4) THEN
          IF (LPRFPIO) THEN
!$OMP ATOMIC
            PRFPIO(IION,MSURF)=PRFPIO(IION,MSURF)+WGHTSG
          ENDIF
          IF (LERFPIO) THEN
!$OMP ATOMIC
            ERFPIO(IION,MSURF)=ERFPIO(IION,MSURF)+EWGHTSG
          ENDIF
          IF (MSURFG.GT.0) THEN
            IF (LPRFPIO) THEN
!$OMP ATOMIC
              PRFPIO(IION,MSURFG)=PRFPIO(IION,MSURFG)+WGHTSG
            ENDIF
            IF (LERFPIO) THEN
!$OMP ATOMIC
              ERFPIO(IION,MSURFG)=ERFPIO(IION,MSURFG)+EWGHTSG
            ENDIF
          ENDIF
          IF (LPRFPIO .OR. LERFPIO) LMETSPW(NSPAM+IION) = .TRUE.
        ELSE
          goto 999
        ENDIF  ! itold
      ELSE
c  no tallies for emitted bulk particles
        goto 999
      ENDIF  ! ityp


      case default  ! ind
        goto 999
      end select ! ind

      RETURN

  999 CONTINUE
      WRITE (IUNOUT,*) 'ERROR EXIT IN UPDATE_SURFACE'
      WRITE (IUNOUT,*) 'IND,ITOLD,ITYP ',IND,ITOLD,ITYP
      CALL EIRENE_EXIT_OWN(1)

      END SUBROUTINE EIRENE_UPDATE_SURFACE
