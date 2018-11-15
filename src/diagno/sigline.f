!pb  100107 ENTRY SIGHA_REINIT added
CDR  parameter PEN introduced, to identify hydrogen line by central energy
Cdr Aug.16:  The identification of particular lines
cdr          by upper and lower energy level (input flags EMIN1,EMAX1 in block 12)
cdr          is not functional in this version, distinct from the manual description
cdr          currently lines can only be identified by their central energy PEN (EMIN1)
cdr          and input parameter EMAX1 is not used at all.
cpb Feb 17:  refresh ADDV tallies (volumetric line emissivities) 
c            not only for new stratum, but also when
c            PEN parameter is different from that from previous call,
c            i.e. a new line is requested for same stratum flag.
cdr Jan 18:  parameter ICHORI added
c   june 18: renamed from sigha (hydrogen only) to sigline (generalized,
c            any transition line)
c            
C
      SUBROUTINE EIRENE_SIGLINE(INIT,JJJ,ZDS,PEN,PSIG,
     .                          DUMMY2,ARGST,ICHORI)
CDR  this routine evaluates ("side on") emissivities of certain transition lines,
cdr  integrated along a line of sight (PSIG) and also the integrand resolved along
cdr  line of sight (ARGST).
cdr new version:  the lines, components and contributions are specified in input block 12.

cdr old version (up to May 2018):
c    Currently there are up to 6 contributions to each particular pre-programmed
c    transition (depending on population coefficient data stored 
c    in file AMJUEL, section H.11 and H.12) 
c  aug.16: available transitions in H-atom:
c          ly-alpha  (2 - 1)
c          ly-beta   (3 - 1)
c          ba-alpha  (3 - 2)
c          ba-beta   (4 - 2)
c          ba-gamma  (5 - 2)
c          ba-delta  (6 - 2)
c    for each of these lines there are separate components from
c    1) coupling to H
c    2) coupling to H+
c    3) coupling to H2
c    4) coupling to H2+
c    5) coupling to H-
c    6) coupling to H3+
c    0) total, sum over these 6 components.
c
c    for each of this components there may be several contributions,
c    e.g. component 1) may have contributions from H, D and T atoms.
c
c
C
C  INPUT:
C          INIT: FLAG FOR INITIALISATION (DO NOT CHANGE!)
C          NCELL (COMPRT): INDEX IN TALLY ARRAYS FOR CURRENT ZONE
C          JJJ:    INDEX OF SEGMENT ALONG CHORD
C          ZDS:    LENGTH OF SEGMENT NO. JJJ
C          PEN:    CENTRAL ENERGY OF LINE (EV)
C  OUTPUT: PSIG:  LINE INTEGRAL OF EMISSION,I=0,6 COMPONENTS
C          ARGST: CONTRIB. FROM CELL NCELL AND CHORD SEGMENT JJJ TO:
C          THE H LINE FLUX PSIG(I),I=0,6 COMPONENTS
C          FROM ATOMS, MOLECULES, TEST IONS, BULK IONS AND NEGATIVE IONS
C          THE INTEGRAND ARGST IS SUCH THAT INTEGR.(ARGST*DL) = PSIG
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CESTIM
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      USE EIRMOD_COMPRT
      USE EIRMOD_COMUSR
      USE EIRMOD_COMSIG
 
      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: INIT, JJJ, ICHORI
      REAL(DP), INTENT(IN) :: ZDS, DUMMY2, PEN
      REAL(DP), INTENT(IN OUT) :: PSIG(0:), ARGST(0:,:)
      REAL(DP) :: PENOLD
      INTEGER, SAVE :: LNO
      INTEGER :: JCOMP, IADV
      INTEGER :: ISTOLD, ISP, NCELC, ICELL, ITROLD
      LOGICAL :: LARGST
      DATA ISTOLD/-1/
      DATA ITROLD/-1/
      DATA PENOLD/-1._DP/
C
      SAVE
C
c     WRITE (IUNOUT,*) 'SIGLINE,INIT,PEN,ISTRA ',
c    .                  INIT,PEN,ISTRA,ISTOLD,IITER,ITROLD

      LARGST = SIZE(ARGST,2) >= NSBOX

      IF (INIT.EQ.0) THEN
        PSIG=0.
        IF (LARGST) ARGST=0.

C  INITIALISE ATOMIC H-LINE ARRAYS FOR CURRENT STRATUM ?
        IF ((ISTRA .NE. ISTOLD) .OR. (IITER .NE. ITROLD) .OR.
     .      (PEN .NE. PENOLD) ) then
c  new, unified routine for line emissivities, replacing: Ly_alpha, Ba_alpha, Ba_beta, etc.

c  identify the selected emission line LNO from the input flags.
          CALL EIRENE_FIND_EMIS_LINE (ISTRA,ICHORI,PEN,LNO)

        endif   ! additional tallies ADDV are now filled, for new LINE, and for present stratum

        ISTOLD=ISTRA
        ITROLD=IITER
        PENOLD=PEN
        RETURN
      ENDIF
C
C  LINE INTEGRAL: PHOTONS/SEC/CM**2
C
C
      ncelc=ncltal(ncell)

      IF (LNO == 0) THEN
! NO MATCHING EMISSION LINE FOUND 
        PSIG(0) = 0._DP
        IF (LARGST) ARGST(0,JJJ) = 0._DP
      ELSE
! USE DATA PROVIDED FOR EMISSION LINE LNO
        DO JCOMP = 1, EMIS_LINES(LNO)%NUM_COMPO
          IADV = EMIS_LINES(LNO)%COMPO(JCOMP)%IADV
          PSIG(JCOMP) = PSIG(JCOMP) + ZDS*ADDV(IADV,NCELC)
          IF (LARGST) ARGST(JCOMP,JJJ) = ADDV(IADV,NCELC)
        END DO

c  sum over components of line LNO
        IADV = EMIS_LINES(LNO)%IADV_TOTAL

        PSIG(0) = PSIG(0) + ZDS*ADDV(IADV,NCELC)
        IF (LARGST) ARGST(0,JJJ) = ADDV(IADV,NCELC)
      END IF
C
      RETURN
 
C     Following lines added for reinitialisation of eirene (DMH)
 
      ENTRY EIRENE_SIGLINE_REINIT
      ISTOLD = -1
      ITROLD = -1
      PENOLD = -1._DP
      RETURN

      END
