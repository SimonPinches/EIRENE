!pb  24.04.07:  allow for logarithmic equidistant energy bins
cdr  29.09.14:  only comments
cdr             All current calls are with either isc=0 or isc=1.
cdr             isc=2: is apparently unused but ready?
cdr  16.05.19:  log spectra disabled for directional spectra (with negative signs of EB possible)
cdr             IDIREC=1 option is not available here for surface spectra (forgotten?).
cdr             but already programed in OUTSPEC. "Failed Safe" added for now.
cdr             Lots of identical code twice in ISC=0 and ISC > 0. Can be reduced.
cdr  29.05.19:  Bug fix: binning for directional volumetrix spectra
cdr             use E0_par times sign(vel,chord), 
cdr             to bin in (parallel) energy units with sign.
cdr             More convenient units: parallel velocity (parallel to chord)
cdr             with sign (tbd)

      SUBROUTINE EIRENE_UPDATE_SPECTRUM (WT,IND,ISC)
C  update contributions to surface- or volume/line-averaged energy spectra
c  wt:  particle weight, (or wt=wpr, conditional particle weight)
c
c  cell crossing   : (conditional) tracklength estimator for cell based spectra
c  surface crossing: here tracklength estim. collapses to a collision estim.

c  isc:    =0: update surface-averaged spectra,
c       ind:  =1: particle incident on surface
c       ind:  =2: particle reemitted from surface

c  isc:  =1,2: else (update cell based spectra)
c       isc =1:  score in coarse (scoring) grid
c       isc =2:  score in fine (geometry)  grid
c       ind:  not in use  (often: ind = iflag in calling programs,
c                          IFLAG is a flag used for special (non-standard)
c                          options for volume averged tally estimators)
c  ityp:  type of particle

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CESTIM
      USE EIRMOD_COMPRT
      USE EIRMOD_CUPD
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      USE EIRMOD_CCONA
      USE EIRMOD_COMUSR
      USE EIRMOD_CZT1

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IND, ISC
      REAL(DP), INTENT(IN) :: WT
      INTEGER :: ISPC, I, IS, IC, IRDO, IRD
      REAL(DP) :: ADD, WV, DIST, WTR, SPCVX, SPCVY, SPCVZ, CDYN, 
     .            EB, SIG

      TYPE(EIRENE_SPECTRUM), POINTER :: P

c  currently: surface based spectra only from particles incident onto surface (ind=1)
c             no spectra of emitted particles (ind=2) yet

      IF ((ISC == 0) .AND. (IND .NE. 1)) RETURN

C  set "type" specific parameters:  IS, CDYN
      SELECT CASE (ITYP)
      CASE (0)
        IS = IPHOT
        CDYN = 1._DP
      CASE (1)
        IS = IATM
        CDYN = CNDYNA(IATM)
      CASE (2)
        IS = IMOL
        CDYN = CNDYNM(IMOL)
      CASE (3)
        IS = IION
        CDYN = CNDYNI(IION)
      CASE (4)
        IS = IPLS
        CDYN = CNDYNP(IPLS)
      END SELECT

      IF (ISC == 0) THEN
! SURFACE-AVERAGED SPECTRA

cdr  IDIREC=1 option not written. See ISC > 0 for required coding.

        DO ISPC=1,NADSPC
          P => ESTIML(ISPC)
          IF ((P%ISRFCLL == ISC) .AND.
     .        (P%ISPCSRF == MSURF) .AND.
     .        (P%IPRTYP == ITYP) .AND.
     .        ((P%IPRSP == IS) .OR. (P%IPRSP == 0))) THEN

            SELECT CASE(ESTIML(ISPC)%ISPCTYP)
            CASE (1)
              ADD = WT  ! particle flux per bin
            CASE (2)
              ADD = WT*E0 ! energy-weighted flux (power flux) per bin
            CASE (3)
              ADD = WTR*VEL*CDYN
            CASE DEFAULT
              ADD = 0._DP ! no scoring
            END SELECT

            EB = E0

            IF (ESTIML(ISPC)%IDIREC > 0) THEN
              WRITE (IUNOUT,*) 'error in update_spectrum (block 10F)'
              WRITE (IUNOUT,*) 'Unfinished option, spectrum ISPC:',ISPC
              write (iunout,*) 'Directional surface spectrum ?'
              call eirene_exit_own(1)
            ENDIF

            IF (ESTIML(ISPC)%LOG) EB=LOG10(EB)

c  find the bin I  (energy units)
            IF (EB < ESTIML(ISPC)%SPCMIN) THEN
              I = 0
            ELSEIF (EB >= ESTIML(ISPC)%SPCMAX) THEN
              I = ESTIML(ISPC)%NSPC + 1
            ELSE
              I = INT((EB - ESTIML(ISPC)%SPCMIN) *
     .                 ESTIML(ISPC)%SPCDELI + 1)
            END IF

cdr  score SPC(I), ESP_MIN and ESP_MAX
            ESTIML(ISPC)%SPC(I) = ESTIML(ISPC)%SPC(I) + ADD
            ESTIML(ISPC)%ESP_MIN= MIN(ESTIML(ISPC)%ESP_MIN,EB)
            ESTIML(ISPC)%ESP_MAX= MAX(ESTIML(ISPC)%ESP_MAX,EB)
            ESTIML(ISPC)%IMETSP = 1
          END IF
        END DO

      ELSE

! CELL BASED SPECTRA

cdr  meaning of isc = 1,2  see subr. input, flag ISRFCLL
cdr  meaning of ind:   not in use for cell based spectra    ??  iflag in calling program ??

        WV=WEIGHT/VEL
        DO IC=1,NCOU
          DIST=CLPD(IC)
          WTR=WV*DIST
          IRDO=NRCELL+NUPC(IC)*NR1P2+NBLCKA
          IRD=NCLTAL(IRDO)


          DO ISPC=1,NADSPC
            P => ESTIML(ISPC)
            IF ((P%ISRFCLL > 0) .AND.
     .          (((P%ISRFCLL == 1).AND.(P%ISPCSRF == IRD)) .OR.      ! scoring cell, coarse grid
     .           ((P%ISRFCLL == 2).AND.(P%ISPCSRF == IRDO))) .AND.   ! geometry cell, fine grid
     .          (P%IPRTYP == ITYP) .AND.
     .          ((P%IPRSP == IS) .OR. (P%IPRSP == 0))) THEN

              SELECT CASE(ESTIML(ISPC)%ISPCTYP)
              CASE (1)
                ADD = WTR       ! particle density per bin
              CASE (2)
                ADD = WTR*E0    ! energy density per bin
              CASE (3)
                ADD = WTR*VEL*CDYN
              CASE DEFAULT
                ADD = 0._DP
              END SELECT
cdr
cdr  binning according to value of parameter EB
              EB = E0
              IF (ESTIML(ISPC)%IDIREC > 0) THEN
                SPCVX = ESTIML(ISPC)%SPCVX
                SPCVY = ESTIML(ISPC)%SPCVY
                SPCVZ = ESTIML(ISPC)%SPCVZ
cdr   both signs possible for EB now.
cdr             EB = EB * (SPCVX*VELX+SPCVY*VELY+SPCVZ*VELZ)
cdr   bug fix, with Sergej Makarov, May 29th 2019
                SIG=SIGN(1._DP, SPCVX*VELX+SPCVY*VELY+SPCVZ*VELZ)
                EB = SIG*EB * (SPCVX*VELX+SPCVY*VELY+SPCVZ*VELZ)**2
              END IF
cdr   log scale makes sense only for IDIREC=0 option. Use ABS(EB) for safety
              IF (ESTIML(ISPC)%LOG) EB=LOG10(abs(EB))


c  find the bin I  (energy units)
              IF (EB < ESTIML(ISPC)%SPCMIN) THEN
                I = 0
              ELSEIF (EB >= ESTIML(ISPC)%SPCMAX) THEN
                I = ESTIML(ISPC)%NSPC + 1
              ELSE
                I = INT((EB - ESTIML(ISPC)%SPCMIN) *
     .                   ESTIML(ISPC)%SPCDELI + 1)
              END IF

cdr  score SPC(I), ESP_MIN and ESP_MAX
              ESTIML(ISPC)%SPC(I) = ESTIML(ISPC)%SPC(I) + ADD
              ESTIML(ISPC)%ESP_MIN= MIN(ESTIML(ISPC)%ESP_MIN,EB)
              ESTIML(ISPC)%ESP_MAX= MAX(ESTIML(ISPC)%ESP_MAX,EB)
              ESTIML(ISPC)%IMETSP = 1
            END IF
          END DO

        END DO

      END IF

      RETURN
      END SUBROUTINE EIRENE_UPDATE_SPECTRUM

csw 21oct08
      SUBROUTINE EIRENE_update_spectrum_reinit
      IMPLICIT NONE

      return
csw
      END SUBROUTINE EIRENE_update_spectrum_reinit
