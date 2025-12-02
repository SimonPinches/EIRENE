cdr jan. 2018: added: outgoing current tallies: ind=1
cdr sept.2014: only comments....

c  SCORE "INCIDENT" and "EMITTED" SURFACE FLUX TALLIES, FOR SURFACE MSURF, OR SURFACE SEGMENT MSURFG.
c  INCIDENT: update tallies POT_A_B(iout,msurf) and EOT_A_B(iout,msurf)
c  EMITTED : update tallies PRF_A_B(iout,msurf) and ERF_A_B(iout,msurf)
c  A code-letter for incident type of particle: A, M, I, P, PH
c  B code-letter for emitted type of particle : AT, ML, IO, PL, PHT
c  iout:  species index for emitted particle

      SUBROUTINE EIRENE_UPDATE_SURFACE (ITOLD,IOLD,WGHTSG,IND)

C       PARTICLE FLUXES (WEIGHT=WGHTSG),     score POT.., PRF...
C       ENERGY FLUXES   (E0*WGHTSG),         score EOT.., ERF...

c  input:
c
c  ind=1: score incident currents
c    itold:  type of incident particle
c    iold:   index of incident particle
c    ispez  (iatm, imol, iion, ipls, iphot): of incident particle
c    msurf:  surface index
c    msurfg: sub-segment of surface MSURF, for spatial resolution on surface
c    E0:     energy (eV) of incident particle
c    WGHTSG: stat. weight of incident particle WEIGHT* PROB* SIGN.

c  ind=2: score reemitted currents
c    itold:  type of incident particle
c    iold:   index of incident particle
c    ityp :  type of emitted particle
c    ispez  (iatm, imol, iion, ipls, iphot): of emitted particle
c    msurf:  surface index
c    msurfg: sub-segment of surface MSURF, for spatial resolution on surface
c    E0:     energy (eV) of reemitted particle
c    WGHTSG: stat. weight of reemitted particle* PROB* SIGN.

c  output:
c  lmetspw(ispz): species ispz is emitted, emitted flux tally is scored.

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CESTIM
      USE EIRMOD_COMPRT
      USE EIRMOD_COMUSR
      USE EIRMOD_CSPEZ
      USE EIRMOD_CSDVI
      USE EIRMOD_CLOGAU

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ITOLD, IOLD, IND
      REAL(dp), INTENT(IN) :: WGHTSG
      REAL(dp) ::             EWGHTSG
      EXTERNAL :: EIRENE_EXIT_OWN

      IF (MSURF .LE. 0) RETURN

      EWGHTSG=E0*WGHTSG

      select case(IND)

      CASE(1)
C  incident fluxes
      IF (ITOLD.EQ.0) THEN
C  INCIDENT PHOTONS
        IF (LEOTPHT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
          EOTPHT(IPHOT,MSURF)=EOTPHT(IPHOT,MSURF)+EWGHTSG
        ENDIF
        IF (LPOTPHT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
          POTPHT(IPHOT,MSURF)=POTPHT(IPHOT,MSURF)+WGHTSG
        ENDIF

        IF (MSURFG.GT.0) THEN
          IF (LEOTPHT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            EOTPHT(IPHOT,MSURFG)=EOTPHT(IPHOT,MSURFG)+EWGHTSG
          ENDIF
          IF (LPOTPHT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            POTPHT(IPHOT,MSURFG)=POTPHT(IPHOT,MSURFG)+WGHTSG
          ENDIF
        ENDIF
        IF (LEOTPHT .OR. LPOTPHT) LMETSPW(IPHOT) = .TRUE.
      ELSEIF (ITOLD.EQ.1) THEN
C  INCIDENT ATOMS
        IF (LEOTAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
          EOTAT(IATM,MSURF)=EOTAT(IATM,MSURF)+EWGHTSG
        ENDIF
        IF (LPOTAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
          POTAT(IATM,MSURF)=POTAT(IATM,MSURF)+WGHTSG
        ENDIF

        IF (MSURFG.GT.0) THEN
          IF (LEOTAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            EOTAT(IATM,MSURFG)=EOTAT(IATM,MSURFG)+EWGHTSG
          ENDIF
          IF (LPOTAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            POTAT(IATM,MSURFG)=POTAT(IATM,MSURFG)+WGHTSG
          ENDIF
        ENDIF
        IF (LEOTAT .OR. LPOTAT) LMETSPW(NSPH+IATM) = .TRUE.
      ELSEIF (ITOLD.EQ.2) THEN
C  INCIDENT MOLECULES
        IF (LEOTML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
          EOTML(IMOL,MSURF)=EOTML(IMOL,MSURF)+EWGHTSG
        ENDIF
        IF (LPOTML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
          POTML(IMOL,MSURF)=POTML(IMOL,MSURF)+WGHTSG
        ENDIF

        IF (MSURFG.GT.0) THEN
          IF (LEOTML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            EOTML(IMOL,MSURFG)=EOTML(IMOL,MSURFG)+EWGHTSG
          ENDIF
          IF (LPOTML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            POTML(IMOL,MSURFG)=POTML(IMOL,MSURFG)+WGHTSG
          ENDIF
        ENDIF
        IF (LEOTML .OR. LPOTML) LMETSPW(NSPA+IMOL) = .TRUE.
      ELSEIF (ITOLD.EQ.3) THEN
C  INCIDENT TEST IONS
        IF (LEOTIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
          EOTIO(IION,MSURF)=EOTIO(IION,MSURF)+EWGHTSG
        ENDIF
        IF (LPOTIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
          POTIO(IION,MSURF)=POTIO(IION,MSURF)+WGHTSG
        ENDIF

        IF (MSURFG.GT.0) THEN
          IF (LEOTIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            EOTIO(IION,MSURFG)=EOTIO(IION,MSURFG)+EWGHTSG
          ENDIF
          IF (LPOTIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            POTIO(IION,MSURFG)=POTIO(IION,MSURFG)+WGHTSG
          ENDIF
        ENDIF
        IF (LEOTIO .OR. LPOTIO) LMETSPW(NSPAM+IION) = .TRUE.
      ELSEIF (ITOLD.EQ.4) THEN
C  INCIDENT BULK IONS
        IF (LEOTPL) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
          EOTPL(IPLS,MSURF)=EOTPL(IPLS,MSURF)+EWGHTSG
        ENDIF
        IF (LPOTPL) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
          POTPL(IPLS,MSURF)=POTPL(IPLS,MSURF)+WGHTSG
        ENDIF

        IF (MSURFG.GT.0) THEN
          IF (LEOTPL) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            EOTPL(IPLS,MSURFG)=EOTPL(IPLS,MSURFG)+EWGHTSG
          ENDIF
          IF (LPOTPL) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
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
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            PRFPHPHT(IPHOT,MSURF)=PRFPHPHT(IPHOT,MSURF)+WGHTSG
            IF (NLSPCSCL_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
              PRFPHPHT2(1:NPHOT,0:NPHOT) => PRFPHPHT(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFPHPHT2(IPHOT,IOLD)=PRFPHPHT2(IPHOT,IOLD)+WGHTSG
              LMETSPW2(1:NPHOT,0:NPHOT) => LMETSPW(NTWS_PHI+1:NTWS_PHPH)
              LMETSPW2(IPHOT,0) = .TRUE.
              LMETSPW2(IPHOT,IOLD) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFPHPHT(IPHOT+IOLD*NPHOT,MSURF)=
     .         PRFPHPHT(IPHOT+IOLD*NPHOT,MSURF)+WGHTSG
              LMETSPW(NTWS_PHI+IPHOT) = .TRUE.
              LMETSPW(NTWS_PHI+IPHOT+IOLD*NPHOT) = .TRUE.
#endif
            END IF
          ENDIF
          IF (LERFPHPHT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            ERFPHPHT(IPHOT,MSURF)=ERFPHPHT(IPHOT,MSURF)+EWGHTSG
          ENDIF
          IF (MSURFG.GT.0) THEN
            IF (LPRFPHPHT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFPHPHT(IPHOT,MSURFG)=PRFPHPHT(IPHOT,MSURFG)+WGHTSG
              IF (NLSPCSCL_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFPHPHT2(1:NPHOT,0:NPHOT) => PRFPHPHT(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFPHPHT2(IPHOT,IOLD)=PRFPHPHT2(IPHOT,IOLD)+WGHTSG
                LMETSPW2(1:NPHOT,0:NPHOT)=>LMETSPW(NTWS_PHI+1:NTWS_PHPH)
                LMETSPW2(IPHOT,0) = .TRUE.
                LMETSPW2(IPHOT,IOLD) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFPHPHT(IPHOT+IOLD*NPHOT,MSURFG)=
     .           PRFPHPHT(IPHOT+IOLD*NPHOT,MSURFG)+WGHTSG
                LMETSPW(NTWS_PHI+IPHOT) = .TRUE.
                LMETSPW(NTWS_PHI+IPHOT+IOLD*NPHOT) = .TRUE.
#endif
              END IF
            ENDIF
            IF (LERFPHPHT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              ERFPHPHT(IPHOT,MSURFG)=ERFPHPHT(IPHOT,MSURFG)+EWGHTSG
            ENDIF
          ENDIF
          IF (LPRFPHPHT .OR. LERFPHPHT) LMETSPW(IPHOT) = .TRUE.
        ELSE
          goto 999
        ENDIF

c  an atom IATM is reemitted...
      ELSEIF (ITYP.EQ.1) THEN
        LOGATM(IATM,ISTRA)=.TRUE.
        IF (ITOLD.EQ.1) THEN
c... from an incident atom IOLD
          IF (LPRFAAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            PRFAAT(IATM,MSURF)=PRFAAT(IATM,MSURF)+WGHTSG
            IF (NLSPCSCL_ATM) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
              PRFAAT2(1:NATM,0:NATM) => PRFAAT(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFAAT2(IATM,IOLD)=PRFAAT2(IATM,IOLD)+WGHTSG
              LMETSPW2(1:NATM,0:NATM)=>LMETSPW(NSPZTOTWS+1:NTWS_AA)
              LMETSPW2(IATM,0) = .TRUE.
              LMETSPW2(IATM,IOLD) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFAAT(IATM+IOLD*NATM,MSURF)=
     .         PRFAAT(IATM+IOLD*NATM,MSURF)+WGHTSG
              LMETSPW(NSPZTOTWS+IATM) = .TRUE.
              LMETSPW(NSPZTOTWS+IATM+IOLD*NATM) = .TRUE.
#endif
            END IF
          ENDIF
          IF (LERFAAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            ERFAAT(IATM,MSURF)=ERFAAT(IATM,MSURF)+EWGHTSG
          ENDIF
          IF (MSURFG.GT.0) THEN
            IF (LPRFAAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFAAT(IATM,MSURFG)=PRFAAT(IATM,MSURFG)+WGHTSG
              IF (NLSPCSCL_ATM) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFAAT2(1:NATM,0:NATM) => PRFAAT(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFAAT2(IATM,IOLD)=PRFAAT2(IATM,IOLD)+WGHTSG
                LMETSPW2(1:NATM,0:NATM)=>LMETSPW(NSPZTOTWS+1:NTWS_AA)
                LMETSPW2(IATM,0) = .TRUE.
                LMETSPW2(IATM,IOLD) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFAAT(IATM+IOLD*NATM,MSURFG)=
     .           PRFAAT(IATM+IOLD*NATM,MSURFG)+WGHTSG
                LMETSPW(NSPZTOTWS+IATM) = .TRUE.
                LMETSPW(NSPZTOTWS+IATM+IOLD*NATM) = .TRUE.
#endif
              END IF
            ENDIF
            IF (LERFAAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              ERFAAT(IATM,MSURFG)=ERFAAT(IATM,MSURFG)+EWGHTSG
            ENDIF
          ENDIF
          IF (LPRFAAT .OR. LERFAAT) LMETSPW(NSPH+IATM) = .TRUE.
        ELSEIF (ITOLD.EQ.2) THEN
c... from an incident molecule IOLD
          IF (LPRFMAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            PRFMAT(IATM,MSURF)=PRFMAT(IATM,MSURF)+WGHTSG
            IF (NLSPCSCL_MOL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
              PRFMAT2(1:NATM,0:NMOL) => PRFMAT(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFMAT2(IATM,IOLD)=PRFMAT2(IATM,IOLD)+WGHTSG
              LMETSPW2(1:NATM,0:NMOL)=>LMETSPW(NTWS_AA+1:NTWS_AM)
              LMETSPW2(IATM,0) = .TRUE.
              LMETSPW2(IATM,IOLD) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFMAT(IATM+IOLD*NATM,MSURF)=
     .         PRFMAT(IATM+IOLD*NATM,MSURF)+WGHTSG
              LMETSPW(NTWS_AA+IATM) = .TRUE.
              LMETSPW(NTWS_AA+IATM+IOLD*NATM) = .TRUE.
#endif
            END IF
          ENDIF
          IF (LERFMAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            ERFMAT(IATM,MSURF)=ERFMAT(IATM,MSURF)+EWGHTSG
          ENDIF
          IF (MSURFG.GT.0) THEN
            IF (LPRFMAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFMAT(IATM,MSURFG)=PRFMAT(IATM,MSURFG)+WGHTSG
              IF (NLSPCSCL_MOL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFMAT2(1:NATM,0:NMOL) => PRFMAT(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFMAT2(IATM,IOLD)=PRFMAT2(IATM,IOLD)+WGHTSG
                LMETSPW2(1:NATM,0:NMOL)=>LMETSPW(NTWS_AA+1:NTWS_AM)
                LMETSPW2(IATM,0) = .TRUE.
                LMETSPW2(IATM,IOLD) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFMAT(IATM+IOLD*NATM,MSURFG)=
     .           PRFMAT(IATM+IOLD*NATM,MSURFG)+WGHTSG
                LMETSPW(NTWS_AA+IATM) = .TRUE.
                LMETSPW(NTWS_AA+IATM+IOLD*NATM) = .TRUE.
#endif
              END IF
            ENDIF
            IF (LERFMAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              ERFMAT(IATM,MSURFG)=ERFMAT(IATM,MSURFG)+EWGHTSG
            ENDIF
          ENDIF
          IF (LPRFMAT .OR. LERFMAT) LMETSPW(NSPH+IATM) = .TRUE.
        ELSEIF (ITOLD.EQ.3) THEN
c... from an incident test ion IOLD
          IF (LPRFIAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            PRFIAT(IATM,MSURF)=PRFIAT(IATM,MSURF)+WGHTSG
            IF (NLSPCSCL_ION) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
              PRFIAT2(1:NATM,0:NION) => PRFIAT(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFIAT2(IATM,IOLD)=PRFIAT2(IATM,IOLD)+WGHTSG
              LMETSPW2(1:NATM,0:NION)=>LMETSPW(NTWS_AM+1:NTWS_AI)
              LMETSPW2(IATM,0) = .TRUE.
              LMETSPW2(IATM,IOLD) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFIAT(IATM+IOLD*NATM,MSURF)=
     .         PRFIAT(IATM+IOLD*NATM,MSURF)+WGHTSG
              LMETSPW(NTWS_AM+IATM) = .TRUE.
              LMETSPW(NTWS_AM+IATM+IOLD*NATM) = .TRUE.
#endif
            END IF
          ENDIF
          IF (LERFIAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            ERFIAT(IATM,MSURF)=ERFIAT(IATM,MSURF)+EWGHTSG
          ENDIF
          IF (MSURFG.GT.0) THEN
            IF (LPRFIAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFIAT(IATM,MSURFG)=PRFIAT(IATM,MSURFG)+WGHTSG
              IF (NLSPCSCL_ION) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFIAT2(1:NATM,0:NION) => PRFIAT(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFIAT2(IATM,IOLD)=PRFIAT2(IATM,IOLD)+WGHTSG
                LMETSPW2(1:NATM,0:NION)=>LMETSPW(NTWS_AM+1:NTWS_AI)
                LMETSPW2(IATM,0) = .TRUE.
                LMETSPW2(IATM,IOLD) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFIAT(IATM+IOLD*NATM,MSURFG)=
     .           PRFIAT(IATM+IOLD*NATM,MSURFG)+WGHTSG
                LMETSPW(NTWS_AM+IATM) = .TRUE.
                LMETSPW(NTWS_AM+IATM+IOLD*NATM) = .TRUE.
#endif
              END IF
            ENDIF
            IF (LERFIAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              ERFIAT(IATM,MSURFG)=ERFIAT(IATM,MSURFG)+EWGHTSG
            ENDIF
          ENDIF
          IF (LPRFIAT .OR. LERFIAT) LMETSPW(NSPH+IATM) = .TRUE.
        ELSEIF (ITOLD.EQ.4) THEN
c... from an incident bulk ion
          IF (LPRFPAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            PRFPAT(IATM,MSURF)=PRFPAT(IATM,MSURF)+WGHTSG
          ENDIF
          IF (LERFPAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            ERFPAT(IATM,MSURF)=ERFPAT(IATM,MSURF)+EWGHTSG
          ENDIF
          IF (MSURFG.GT.0) THEN
            IF (LPRFPAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFPAT(IATM,MSURFG)=PRFPAT(IATM,MSURFG)+WGHTSG
            ENDIF
            IF (LERFPAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              ERFPAT(IATM,MSURFG)=ERFPAT(IATM,MSURFG)+EWGHTSG
            ENDIF
          ENDIF
          IF (LPRFPAT .OR. LERFPAT) LMETSPW(NSPH+IATM) = .TRUE.
        ELSE
          goto 999
        ENDIF

c  a molecule IMOL is reemitted
      ELSEIF (ITYP.EQ.2) THEN
        LOGMOL(IMOL,ISTRA)=.TRUE.
        IF (ITOLD.EQ.1) THEN
c... from an incident atom IOLD
          IF (LPRFAML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            PRFAML(IMOL,MSURF)=PRFAML(IMOL,MSURF)+WGHTSG
            IF (NLSPCSCL_ATM) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
              PRFAML2(1:NMOL,0:NATM) => PRFAML(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFAML2(IMOL,IOLD)=PRFAML2(IMOL,IOLD)+WGHTSG
              LMETSPW2(1:NMOL,0:NATM)=>LMETSPW(NTWS_APH+1:NTWS_MA)
              LMETSPW2(IMOL,0) = .TRUE.
              LMETSPW2(IMOL,IOLD) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFAML(IMOL+IOLD*NMOL,MSURF)=
     .         PRFAML(IMOL+IOLD*NMOL,MSURF)+WGHTSG
              LMETSPW(NTWS_APH+IMOL) = .TRUE.
              LMETSPW(NTWS_APH+IMOL+IOLD*NMOL) = .TRUE.
#endif
            END IF
          ENDIF
          IF (LERFAML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            ERFAML(IMOL,MSURF)=ERFAML(IMOL,MSURF)+EWGHTSG
          ENDIF
          IF (MSURFG.GT.0) THEN
            IF (LPRFAML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFAML(IMOL,MSURFG)=PRFAML(IMOL,MSURFG)+WGHTSG
              IF (NLSPCSCL_ATM) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFAML2(1:NMOL,0:NATM) => PRFAML(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFAML2(IMOL,IOLD)=PRFAML2(IMOL,IOLD)+WGHTSG
                LMETSPW2(1:NMOL,0:NATM)=>LMETSPW(NTWS_APH+1:NTWS_MA)
                LMETSPW2(IMOL,0) = .TRUE.
                LMETSPW2(IMOL,IOLD) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFAML(IMOL+IOLD*NMOL,MSURFG)=
     .           PRFAML(IMOL+IOLD*NMOL,MSURFG)+WGHTSG
                LMETSPW(NTWS_APH+IMOL) = .TRUE.
                LMETSPW(NTWS_APH+IMOL+IOLD*NMOL) = .TRUE.
#endif
              END IF
            ENDIF
            IF (LERFAML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              ERFAML(IMOL,MSURFG)=ERFAML(IMOL,MSURFG)+EWGHTSG
            ENDIF
          ENDIF
          IF (LPRFAML .OR. LERFAML) LMETSPW(NSPA+IMOL) = .TRUE.
        ELSEIF (ITOLD.EQ.2) THEN
c... from an incident molecule IOLD
          IF (LPRFMML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            PRFMML(IMOL,MSURF)=PRFMML(IMOL,MSURF)+WGHTSG
            IF (NLSPCSCL_MOL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
              PRFMML2(1:NMOL,0:NMOL) => PRFMML(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFMML2(IMOL,IOLD)=PRFMML2(IMOL,IOLD)+WGHTSG
              LMETSPW2(1:NMOL,0:NMOL)=>LMETSPW(NTWS_MA+1:NTWS_MM)
              LMETSPW2(IMOL,0) = .TRUE.
              LMETSPW2(IMOL,IOLD) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFMML(IMOL+IOLD*NMOL,MSURF)=
     .         PRFMML(IMOL+IOLD*NMOL,MSURF)+WGHTSG
              LMETSPW(NTWS_MA+IMOL) = .TRUE.
              LMETSPW(NTWS_MA+IMOL+IOLD*NMOL) = .TRUE.
#endif
            END IF
          ENDIF
          IF (LERFMML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            ERFMML(IMOL,MSURF)=ERFMML(IMOL,MSURF)+EWGHTSG
          ENDIF
          IF (MSURFG.GT.0) THEN
            IF (LPRFMML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFMML(IMOL,MSURFG)=PRFMML(IMOL,MSURFG)+WGHTSG
              IF (NLSPCSCL_MOL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFMML2(1:NMOL,0:NMOL) => PRFMML(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFMML2(IMOL,IOLD)=PRFMML2(IMOL,IOLD)+WGHTSG
                LMETSPW2(1:NMOL,0:NMOL)=>LMETSPW(NTWS_MA+1:NTWS_MM)
                LMETSPW2(IMOL,0) = .TRUE.
                LMETSPW2(IMOL,IOLD) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFMML(IMOL+IOLD*NMOL,MSURFG)=
     .           PRFMML(IMOL+IOLD*NMOL,MSURFG)+WGHTSG
                LMETSPW(NTWS_MA+IMOL) = .TRUE.
                LMETSPW(NTWS_MA+IMOL+IOLD*NMOL) = .TRUE.
#endif
              END IF
            ENDIF
            IF (LERFMML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              ERFMML(IMOL,MSURFG)=ERFMML(IMOL,MSURFG)+EWGHTSG
            ENDIF
          ENDIF
          IF (LPRFMML .OR. LERFMML) LMETSPW(NSPA+IMOL) = .TRUE.
        ELSEIF (ITOLD.EQ.3) THEN
c... from an incident test ion IOLD
          IF (LPRFIML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            PRFIML(IMOL,MSURF)=PRFIML(IMOL,MSURF)+WGHTSG
            IF (NLSPCSCL_ION) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
              PRFIML2(1:NMOL,0:NION) => PRFIML(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFIML2(IMOL,IOLD)=PRFIML2(IMOL,IOLD)+WGHTSG
              LMETSPW2(1:NMOL,0:NION)=>LMETSPW(NTWS_MM+1:NTWS_MI)
              LMETSPW2(IMOL,0) = .TRUE.
              LMETSPW2(IMOL,IOLD) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFIML(IMOL+IOLD*NMOL,MSURF)=
     .         PRFIML(IMOL+IOLD*NMOL,MSURF)+WGHTSG
              LMETSPW(NTWS_MM+IMOL) = .TRUE.
              LMETSPW(NTWS_MM+IMOL+IOLD*NMOL) = .TRUE.
#endif
            END IF
          ENDIF
          IF (LERFIML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            ERFIML(IMOL,MSURF)=ERFIML(IMOL,MSURF)+EWGHTSG
          ENDIF
          IF (MSURFG.GT.0) THEN
            IF (LPRFIML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFIML(IMOL,MSURFG)=PRFIML(IMOL,MSURFG)+WGHTSG
              IF (NLSPCSCL_ION) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFIML2(1:NMOL,0:NION) => PRFIML(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFIML2(IMOL,IOLD)=PRFIML2(IMOL,IOLD)+WGHTSG
                LMETSPW2(1:NMOL,0:NION)=>LMETSPW(NTWS_MM+1:NTWS_MI)
                LMETSPW2(IMOL,0) = .TRUE.
                LMETSPW2(IMOL,IOLD) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFIML(IMOL+IOLD*NMOL,MSURFG)=
     .           PRFIML(IMOL+IOLD*NMOL,MSURFG)+WGHTSG
                LMETSPW(NTWS_MM+IMOL) = .TRUE.
                LMETSPW(NTWS_MM+IMOL+IOLD*NMOL) = .TRUE.
#endif
              END IF
            ENDIF
            IF (LERFIML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              ERFIML(IMOL,MSURFG)=ERFIML(IMOL,MSURFG)+EWGHTSG
            ENDIF
          ENDIF
          IF (LPRFIML .OR. LERFIML) LMETSPW(NSPA+IMOL) = .TRUE.
        ELSEIF (ITOLD.EQ.4) THEN
c... from an incident bulk ion IOLD
          IF (LPRFPML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            PRFPML(IMOL,MSURF)=PRFPML(IMOL,MSURF)+WGHTSG
          ENDIF
          IF (LERFPML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            ERFPML(IMOL,MSURF)=ERFPML(IMOL,MSURF)+EWGHTSG
          ENDIF
          IF (MSURFG.GT.0) THEN
            IF (LPRFPML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFPML(IMOL,MSURFG)=PRFPML(IMOL,MSURFG)+WGHTSG
            ENDIF
            IF (LERFPML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              ERFPML(IMOL,MSURFG)=ERFPML(IMOL,MSURFG)+EWGHTSG
            ENDIF
          ENDIF
          IF (LPRFPML .OR. LERFPML) LMETSPW(NSPA+IMOL) = .TRUE.
        ELSE
          goto 999
        ENDIF

c  a test ion IION is reemitted
      ELSEIF (ITYP.EQ.3) THEN
        LOGION(IION,ISTRA)=.TRUE.
        IF (ITOLD.EQ.1) THEN
          IF (LPRFAIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            PRFAIO(IION,MSURF)=PRFAIO(IION,MSURF)+WGHTSG
            IF (NLSPCSCL_ATM) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
              PRFAIO2(1:NION,0:NATM) => PRFAIO(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFAIO2(IION,IOLD)=PRFAIO2(IION,IOLD)+WGHTSG
              LMETSPW2(1:NION,0:NATM)=>LMETSPW(NTWS_MPH+1:NTWS_IA)
              LMETSPW2(IION,0) = .TRUE.
              LMETSPW2(IION,IOLD) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFAIO(IION+IOLD*NION,MSURF)=
     .         PRFAIO(IION+IOLD*NION,MSURF)+WGHTSG
              LMETSPW(NTWS_MPH+IION) = .TRUE.
              LMETSPW(NTWS_MPH+IION+IOLD*NION) = .TRUE.
#endif
            END IF
          ENDIF
          IF (LERFAIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            ERFAIO(IION,MSURF)=ERFAIO(IION,MSURF)+EWGHTSG
          ENDIF
          IF (MSURFG.GT.0) THEN
            IF (LPRFAIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFAIO(IION,MSURFG)=PRFAIO(IION,MSURFG)+WGHTSG
              IF (NLSPCSCL_ION) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFAIO2(1:NION,0:NATM) => PRFAIO(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFAIO2(IION,IOLD)=PRFAIO2(IION,IOLD)+WGHTSG
                LMETSPW2(1:NION,0:NATM)=>LMETSPW(NTWS_MPH+1:NTWS_IA)
                LMETSPW2(IION,0) = .TRUE.
                LMETSPW2(IION,IOLD) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFAIO(IION+IOLD*NION,MSURFG)=
     .           PRFAIO(IION+IOLD*NION,MSURFG)+WGHTSG
                LMETSPW(NTWS_MPH+IION) = .TRUE.
                LMETSPW(NTWS_MPH+IION+IOLD*NION) = .TRUE.
#endif
              END IF
            ENDIF
            IF (LERFAIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              ERFAIO(IION,MSURFG)=ERFAIO(IION,MSURFG)+EWGHTSG
            ENDIF
          ENDIF
          IF (LPRFAIO .OR. LERFAIO) LMETSPW(NSPAM+IION) = .TRUE.
        ELSEIF (ITOLD.EQ.2) THEN
          IF (LPRFMIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            PRFMIO(IION,MSURF)=PRFMIO(IION,MSURF)+WGHTSG
            IF (NLSPCSCL_MOL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
              PRFMIO2(1:NION,0:NMOL) => PRFMIO(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFMIO2(IION,IOLD)=PRFMIO2(IION,IOLD)+WGHTSG
              LMETSPW2(1:NION,0:NMOL)=>LMETSPW(NTWS_IA+1:NTWS_IM)
              LMETSPW2(IION,0) = .TRUE.
              LMETSPW2(IION,IOLD) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFMIO(IION+IOLD*NION,MSURF)=
     .         PRFMIO(IION+IOLD*NION,MSURF)+WGHTSG
              LMETSPW(NTWS_IA+IION) = .TRUE.
              LMETSPW(NTWS_IA+IION+IOLD*NION) = .TRUE.
#endif
            END IF
          ENDIF
          IF (LERFMIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            ERFMIO(IION,MSURF)=ERFMIO(IION,MSURF)+EWGHTSG
          ENDIF
          IF (MSURFG.GT.0) THEN
            IF (LPRFMIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFMIO(IION,MSURFG)=PRFMIO(IION,MSURFG)+WGHTSG
              IF (NLSPCSCL_MOL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFMIO2(1:NION,0:NMOL) => PRFMIO(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFMIO2(IION,IOLD)=PRFMIO2(IION,IOLD)+WGHTSG
                LMETSPW2(1:NION,0:NMOL)=>LMETSPW(NTWS_IA+1:NTWS_IM)
                LMETSPW2(IION,0) = .TRUE.
                LMETSPW2(IION,IOLD) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFMIO(IION+IOLD*NION,MSURFG)=
     .           PRFMIO(IION+IOLD*NION,MSURFG)+WGHTSG
                LMETSPW(NTWS_IA+IION) = .TRUE.
                LMETSPW(NTWS_IA+IION+IOLD*NION) = .TRUE.
#endif
              END IF
            ENDIF
            IF (LERFMIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              ERFMIO(IION,MSURFG)=ERFMIO(IION,MSURFG)+EWGHTSG
            ENDIF
          ENDIF
          IF (LPRFMIO .OR. LERFMIO) LMETSPW(NSPAM+IION) = .TRUE.
        ELSEIF (ITOLD.EQ.3) THEN
          IF (LPRFIIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            PRFIIO(IION,MSURF)=PRFIIO(IION,MSURF)+WGHTSG
            IF (NLSPCSCL_ION) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
              PRFIIO2(1:NION,0:NION) => PRFIIO(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFIIO2(IION,IOLD)=PRFIIO2(IION,IOLD)+WGHTSG
              LMETSPW2(1:NION,0:NION)=>LMETSPW(NTWS_IM+1:NTWS_II)
              LMETSPW2(IION,0) = .TRUE.
              LMETSPW2(IION,IOLD) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFIIO(IION+IOLD*NION,MSURF)=
     .         PRFIIO(IION+IOLD*NION,MSURF)+WGHTSG
              LMETSPW(NTWS_IM+IION) = .TRUE.
              LMETSPW(NTWS_IM+IION+IOLD*NION) = .TRUE.
#endif
            END IF
          ENDIF
          IF (LERFIIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            ERFIIO(IION,MSURF)=ERFIIO(IION,MSURF)+EWGHTSG
          ENDIF
          IF (MSURFG.GT.0) THEN
            IF (LPRFIIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFIIO(IION,MSURFG)=PRFIIO(IION,MSURFG)+WGHTSG
              IF (NLSPCSCL_ION) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFIIO2(1:NION,0:NION) => PRFIIO(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFIIO2(IION,IOLD)=PRFIIO2(IION,IOLD)+WGHTSG
                LMETSPW2(1:NION,0:NION)=>LMETSPW(NTWS_IM+1:NTWS_II)
                LMETSPW2(IION,0) = .TRUE.
                LMETSPW2(IION,IOLD) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFIIO(IION+IOLD*NION,MSURFG)=
     .           PRFIIO(IION+IOLD*NION,MSURFG)+WGHTSG
                LMETSPW(NTWS_IM+IION) = .TRUE.
                LMETSPW(NTWS_IM+IION+IOLD*NION) = .TRUE.
#endif
              END IF
            ENDIF
            IF (LERFIIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              ERFIIO(IION,MSURFG)=ERFIIO(IION,MSURFG)+EWGHTSG
            ENDIF
          ENDIF
          IF (LPRFIIO .OR. LERFIIO) LMETSPW(NSPAM+IION) = .TRUE.
        ELSEIF (ITOLD.EQ.4) THEN
          IF (LPRFPIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            PRFPIO(IION,MSURF)=PRFPIO(IION,MSURF)+WGHTSG
          ENDIF
          IF (LERFPIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            ERFPIO(IION,MSURF)=ERFPIO(IION,MSURF)+EWGHTSG
          ENDIF
          IF (MSURFG.GT.0) THEN
            IF (LPRFPIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFPIO(IION,MSURFG)=PRFPIO(IION,MSURFG)+WGHTSG
            ENDIF
            IF (LERFPIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
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
