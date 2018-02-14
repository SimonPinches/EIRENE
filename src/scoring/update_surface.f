cdr jan. 2018: added: outgoing current tallies: ind=1
cdr sept.2014: only comments....

c  SCORE "EMITTED" SURFACE FLUX TALLIES, FOR SURFACE MSURF, OR SURFACE SEGMENT MSURFG.
c  update tallies PRF_A_B(iout,msurf) and ERF_A_B(iout,msurf)
c  A code-letter for incident type of particle: A, M, I, P, PH
c  B code-letter for emitted type of particle :  AT, ML, IO, PL, PHT
c  iout:  species index for emitted particle
 
      SUBROUTINE EIRENE_UPDATE_SURFACE (ITOLD,WGHTSG,IND)
 
C       PARTICLE FLUXES (WEIGHT=WGHTSG),     score PRF...
C       ENERGY FLUXES   (E0*WGHTSG),         score ERF...

c  input:
c  itold:  type of incident particle
c  ind=1:  score incident currents
c    msurf:  surface index
c    msurfg:  sub-segment of surface MSURF, for spatial resolution on surface
c    E0:     energy (eV) of incident particle
c    WGHTSG:   stat. weight of incident particle WEIGHT* PROB* SIGN.

c  ind=2:  score re-emitted currents
c    ityp :  type of emitted particle
c    ispez  (iatm, imol, iion, ipls, iphot): of emitted particle
c    msurf:  surface index
c    msurfg:  sub-segment of surface MSURF, for spatial resolution on surface
c    E0:     energy (eV) of re-emitted particle
c    WGHTSG:   stat. weight of re-emitted particle* PROB* SIGN.

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
        IF (LEOTPHT) EOTPHT(IPHOT,MSURF)=EOTPHT(IPHOT,MSURF)+EWGHTSG
        IF (LPOTPHT) POTPHT(IPHOT,MSURF)=POTPHT(IPHOT,MSURF)+WGHTSG

        IF (MSURFG.GT.0) THEN
          IF (LEOTPHT) 
     .      EOTPHT(IPHOT,MSURFG)=EOTPHT(IPHOT,MSURFG)+EWGHTSG
          IF (LPOTPHT) 
     .      POTPHT(IPHOT,MSURFG)=POTPHT(IPHOT,MSURFG)+WGHTSG
        ENDIF
        IF (LEOTPHT .OR. LPOTPHT) LMETSPW(IPHOT) = .TRUE.
      ELSEIF (ITOLD.EQ.1) THEN
C  INCIDENT ATOMS
        IF (LEOTAT) EOTAT(IATM,MSURF)=EOTAT(IATM,MSURF)+EWGHTSG
        IF (LPOTAT) POTAT(IATM,MSURF)=POTAT(IATM,MSURF)+WGHTSG

        IF (MSURFG.GT.0) THEN
          IF (LEOTAT) EOTAT(IATM,MSURFG)=EOTAT(IATM,MSURFG)+EWGHTSG
          IF (LPOTAT) POTAT(IATM,MSURFG)=POTAT(IATM,MSURFG)+WGHTSG
        ENDIF
        IF (LEOTAT .OR. LPOTAT) LMETSPW(NSPH+IATM) = .TRUE.
      ELSEIF (ITOLD.EQ.2) THEN
C  INCIDENT MOLECULES
        IF (LEOTML) EOTML(IMOL,MSURF)=EOTML(IMOL,MSURF)+EWGHTSG
        IF (LPOTML) POTML(IMOL,MSURF)=POTML(IMOL,MSURF)+WGHTSG

        IF (MSURFG.GT.0) THEN
          IF (LEOTML) EOTML(IMOL,MSURFG)=EOTML(IMOL,MSURFG)+EWGHTSG
          IF (LPOTML) POTML(IMOL,MSURFG)=POTML(IMOL,MSURFG)+WGHTSG
        ENDIF
        IF (LEOTML .OR. LPOTML) LMETSPW(NSPA+IMOL) = .TRUE.
      ELSEIF (ITOLD.EQ.3) THEN
C  INCIDENT TEST IONS
        IF (LEOTIO) EOTIO(IION,MSURF)=EOTIO(IION,MSURF)+EWGHTSG
        IF (LPOTIO) POTIO(IION,MSURF)=POTIO(IION,MSURF)+WGHTSG

        IF (MSURFG.GT.0) THEN
          IF (LEOTIO) EOTIO(IION,MSURFG)=EOTIO(IION,MSURFG)+EWGHTSG
          IF (LPOTIO) POTIO(IION,MSURFG)=POTIO(IION,MSURFG)+WGHTSG
        ENDIF
        IF (LEOTIO .OR. LPOTIO) LMETSPW(NSPAM+IION) = .TRUE.
      ELSEIF (ITOLD.EQ.4) THEN
C  INCIDENT BULK IONS
        IF (LEOTPL) EOTPL(IPLS,MSURF)=EOTPL(IPLS,MSURF)+EWGHTSG
        IF (LPOTPL) POTPL(IPLS,MSURF)=POTPL(IPLS,MSURF)+WGHTSG

        IF (MSURFG.GT.0) THEN
          IF (LEOTPL) EOTPL(IPLS,MSURFG)=EOTPL(IPLS,MSURFG)+EWGHTSG
          IF (LPOTPL) POTPL(IPLS,MSURFG)=POTPL(IPLS,MSURFG)+WGHTSG
        ENDIF
        IF (LEOTPL .OR. LPOTPL) LMETSPW(NSPAMI+IPLS) = .TRUE.
      ELSE
        goto 999
      ENDIF  !ITOLD


      CASE(2)
C  re-emitted fluxes
c  a photon is re-emitted. currently only foreseen for incident photons 
      IF (ITYP.EQ.0) THEN
        LOGPHOT(IPHOT,ISTRA)=.TRUE.
        IF (ITOLD.EQ.0) THEN
c... from an incident photon
          IF (LPRFPHPHT)
     .      PRFPHPHT(IPHOT,MSURF)=PRFPHPHT(IPHOT,MSURF)+WGHTSG
          IF (LERFPHPHT)
     .      ERFPHPHT(IPHOT,MSURF)=ERFPHPHT(IPHOT,MSURF)+EWGHTSG

          IF (MSURFG.GT.0) THEN
            IF (LPRFPHPHT)
     .        PRFPHPHT(IPHOT,MSURFG)=PRFPHPHT(IPHOT,MSURFG)+WGHTSG
            IF (LERFPHPHT)
     .        ERFPHPHT(IPHOT,MSURFG)=ERFPHPHT(IPHOT,MSURFG)+EWGHTSG
          ENDIF
          IF (LPRFPHPHT .OR. LERFPHPHT) LMETSPW(IPHOT) = .TRUE.
        ELSE
          goto 999
        ENDIF

c  an atom is re-emitted...
      ELSEIF (ITYP.EQ.1) THEN
        LOGATM(IATM,ISTRA)=.TRUE.
        IF (ITOLD.EQ.1) THEN
c... from an incident atom
          IF (LPRFAAT) PRFAAT(IATM,MSURF)=PRFAAT(IATM,MSURF)+WGHTSG
          IF (LERFAAT)
     .      ERFAAT(IATM,MSURF)=ERFAAT(IATM,MSURF)+EWGHTSG
          IF (MSURFG.GT.0) THEN
            IF (LPRFAAT)
     .        PRFAAT(IATM,MSURFG)=PRFAAT(IATM,MSURFG)+WGHTSG
            IF (LERFAAT)
     .        ERFAAT(IATM,MSURFG)=ERFAAT(IATM,MSURFG)+EWGHTSG
          ENDIF
          IF (LPRFAAT .OR. LERFAAT) LMETSPW(NSPH+IATM) = .TRUE.
        ELSEIF (ITOLD.EQ.2) THEN
c... from an incident molecule
          IF (LPRFMAT) PRFMAT(IATM,MSURF)=PRFMAT(IATM,MSURF)+WGHTSG
          IF (LERFMAT)
     .      ERFMAT(IATM,MSURF)=ERFMAT(IATM,MSURF)+EWGHTSG
          IF (MSURFG.GT.0) THEN
            IF (LPRFMAT)
     .        PRFMAT(IATM,MSURFG)=PRFMAT(IATM,MSURFG)+WGHTSG
            IF (LERFMAT)
     .        ERFMAT(IATM,MSURFG)=ERFMAT(IATM,MSURFG)+EWGHTSG
          ENDIF
          IF (LPRFMAT .OR. LERFMAT) LMETSPW(NSPH+IATM) = .TRUE.
        ELSEIF (ITOLD.EQ.3) THEN
c... from an incident test ion
          IF (LPRFIAT) PRFIAT(IATM,MSURF)=PRFIAT(IATM,MSURF)+WGHTSG
          IF (LERFIAT)
     .      ERFIAT(IATM,MSURF)=ERFIAT(IATM,MSURF)+EWGHTSG
          IF (MSURFG.GT.0) THEN
            IF (LPRFIAT)
     .        PRFIAT(IATM,MSURFG)=PRFIAT(IATM,MSURFG)+WGHTSG
            IF (LERFIAT)
     .        ERFIAT(IATM,MSURFG)=ERFIAT(IATM,MSURFG)+EWGHTSG
          ENDIF
          IF (LPRFIAT .OR. LERFIAT) LMETSPW(NSPH+IATM) = .TRUE.
        ELSEIF (ITOLD.EQ.4) THEN
c... from an incident bulk ion
          IF (LPRFPAT)
     .      PRFPAT(IATM,MSURF)=PRFPAT(IATM,MSURF)+WGHTSG
          IF (LERFPAT)
     .      ERFPAT(IATM,MSURF)=ERFPAT(IATM,MSURF)+EWGHTSG
          IF (MSURFG.GT.0) THEN
            IF (LPRFPAT)
     .        PRFPAT(IATM,MSURFG)=PRFPAT(IATM,MSURFG)+WGHTSG
            IF (LERFPAT)
     .        ERFPAT(IATM,MSURFG)=ERFPAT(IATM,MSURFG)+EWGHTSG
          ENDIF
          IF (LPRFPAT .OR. LERFPAT) LMETSPW(NSPH+IATM) = .TRUE.
        ELSE
          goto 999
        ENDIF

c  a molecule is re-emitted
      ELSEIF (ITYP.EQ.2) THEN
        LOGMOL(IMOL,ISTRA)=.TRUE.
        IF (ITOLD.EQ.1) THEN
c... from an incident atom
          IF (LPRFAML) PRFAML(IMOL,MSURF)=PRFAML(IMOL,MSURF)+WGHTSG
          IF (LERFAML)
     .      ERFAML(IMOL,MSURF)=ERFAML(IMOL,MSURF)+EWGHTSG
          IF (MSURFG.GT.0) THEN
            IF (LPRFAML)
     .        PRFAML(IMOL,MSURFG)=PRFAML(IMOL,MSURFG)+WGHTSG
            IF (LERFAML)
     .        ERFAML(IMOL,MSURFG)=ERFAML(IMOL,MSURFG)+EWGHTSG
          ENDIF
          IF (LPRFAML .OR. LERFAML) LMETSPW(NSPA+IMOL) = .TRUE.
        ELSEIF (ITOLD.EQ.2) THEN
c... from an incident molecule
          IF (LPRFMML) PRFMML(IMOL,MSURF)=PRFMML(IMOL,MSURF)+WGHTSG
          IF (LERFMML)
     .      ERFMML(IMOL,MSURF)=ERFMML(IMOL,MSURF)+EWGHTSG
          IF (MSURFG.GT.0) THEN
            IF (LPRFMML)
     .        PRFMML(IMOL,MSURFG)=PRFMML(IMOL,MSURFG)+WGHTSG
            IF (LERFMML)
     .        ERFMML(IMOL,MSURFG)=ERFMML(IMOL,MSURFG)+EWGHTSG
          ENDIF
          IF (LPRFMML .OR. LERFMML) LMETSPW(NSPA+IMOL) = .TRUE.
        ELSEIF (ITOLD.EQ.3) THEN
c... from an incident test ion
          IF (LPRFIML) PRFIML(IMOL,MSURF)=PRFIML(IMOL,MSURF)+WGHTSG
          IF (LERFIML)
     .      ERFIML(IMOL,MSURF)=ERFIML(IMOL,MSURF)+EWGHTSG
          IF (MSURFG.GT.0) THEN
            IF (LPRFIML)
     .        PRFIML(IMOL,MSURFG)=PRFIML(IMOL,MSURFG)+WGHTSG
            IF (LERFIML)
     .        ERFIML(IMOL,MSURFG)=ERFIML(IMOL,MSURFG)+EWGHTSG
          ENDIF
          IF (LPRFIML .OR. LERFIML) LMETSPW(NSPA+IMOL) = .TRUE.
        ELSEIF (ITOLD.EQ.4) THEN
c... from an incident bulk ion
          IF (LPRFPML)
     .      PRFPML(IMOL,MSURF)=PRFPML(IMOL,MSURF)+WGHTSG
          IF (LERFPML)
     .      ERFPML(IMOL,MSURF)=ERFPML(IMOL,MSURF)+EWGHTSG
          IF (MSURFG.GT.0) THEN
            IF (LPRFPML)
     .        PRFPML(IMOL,MSURFG)=PRFPML(IMOL,MSURFG)+WGHTSG
            IF (LERFPML)
     .        ERFPML(IMOL,MSURFG)=ERFPML(IMOL,MSURFG)+EWGHTSG
          ENDIF
          IF (LPRFPML .OR. LERFPML) LMETSPW(NSPA+IMOL) = .TRUE.
        ELSE
          goto 999
        ENDIF

c  a test-ion is reemitted
      ELSEIF (ITYP.EQ.3) THEN
        LOGION(IION,ISTRA)=.TRUE.
        IF (ITOLD.EQ.1) THEN
          IF (LPRFAIO) PRFAIO(IION,MSURF)=PRFAIO(IION,MSURF)+WGHTSG
          IF (LERFAIO)
     .      ERFAIO(IION,MSURF)=ERFAIO(IION,MSURF)+EWGHTSG
          IF (MSURFG.GT.0) THEN
            IF (LPRFAIO)
     .        PRFAIO(IION,MSURFG)=PRFAIO(IION,MSURFG)+WGHTSG
            IF (LERFAIO)
     .        ERFAIO(IION,MSURFG)=ERFAIO(IION,MSURFG)+EWGHTSG
          ENDIF
          IF (LPRFAIO .OR. LERFAIO) LMETSPW(NSPAM+IION) = .TRUE.
        ELSEIF (ITOLD.EQ.2) THEN
          IF (LPRFMIO) PRFMIO(IION,MSURF)=PRFMIO(IION,MSURF)+WGHTSG
          IF (LERFMIO)
     .      ERFMIO(IION,MSURF)=ERFMIO(IION,MSURF)+EWGHTSG
          IF (MSURFG.GT.0) THEN
            IF (LPRFMIO)
     .        PRFMIO(IION,MSURFG)=PRFMIO(IION,MSURFG)+WGHTSG
            IF (LERFMIO)
     .        ERFMIO(IION,MSURFG)=ERFMIO(IION,MSURFG)+EWGHTSG
          ENDIF
          IF (LPRFMIO .OR. LERFMIO) LMETSPW(NSPAM+IION) = .TRUE.
        ELSEIF (ITOLD.EQ.3) THEN
          IF (LPRFIIO) PRFIIO(IION,MSURF)=PRFIIO(IION,MSURF)+WGHTSG
          IF (LERFIIO)
     .      ERFIIO(IION,MSURF)=ERFIIO(IION,MSURF)+EWGHTSG
          IF (MSURFG.GT.0) THEN
            IF (LPRFIIO)
     .        PRFIIO(IION,MSURFG)=PRFIIO(IION,MSURFG)+WGHTSG
            IF (LERFIIO)
     .        ERFIIO(IION,MSURFG)=ERFIIO(IION,MSURFG)+EWGHTSG
          ENDIF
          IF (LPRFIIO .OR. LERFIIO) LMETSPW(NSPAM+IION) = .TRUE.
        ELSEIF (ITOLD.EQ.4) THEN
          IF (LPRFPIO)
     .      PRFPIO(IION,MSURF)=PRFPIO(IION,MSURF)+WGHTSG
          IF (LERFPIO)
     .      ERFPIO(IION,MSURF)=ERFPIO(IION,MSURF)+EWGHTSG
          IF (MSURFG.GT.0) THEN
            IF (LPRFPIO)
     .        PRFPIO(IION,MSURFG)=PRFPIO(IION,MSURFG)+WGHTSG
            IF (LERFPIO)
     .        ERFPIO(IION,MSURFG)=ERFPIO(IION,MSURFG)+EWGHTSG
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

999   CONTINUE
      WRITE (IUNOUT,*) 'ERROR EXIT IN UPDATE_SURFACE'
      WRITE (IUNOUT,*) 'IND,ITOLD,ITYP ',IND,ITOLD,ITYP
      CALL EIRENE_EXIT_OWN(1)

      END SUBROUTINE EIRENE_UPDATE_SURFACE
