      MODULE EIRMOD_SIGHE
      USE EIRMOD_PRECISION
      
      IMPLICIT NONE
      PRIVATE

      PUBLIC :: EIRENE_SIGHE, EIRENE_SIGHE_REINIT

      REAL(DP), SAVE :: PENOLD=-1._DP
      INTEGER, SAVE :: ISTOLD=-1, ITROLD=-1

      CONTAINS
     
CDR  from W.Z, (HELIUM) DERIVED FROM FROM SIGHA (HYDROGEN)
c
c
      SUBROUTINE EIRENE_SIGHE(INIT,JJJ,ZDS,PEN,PSIG,DUMMY2,ARGST)

CDR  this routine evaluates ("side on") helium atom ("He") emissivities,
cdr  integrated along a line of sight (PSIG) and also the integrand resolved along
cdr  line of sight (ARGST).
c    Currently there are up to 2 contributions to each particular preprogrammed
c    transition (depending on population coefficient data stored
c    in file AMJUEL, section H.11 and H.12
c  aug.16: available transitions in H-atom:
c          ba-alpha  (3S - 2P, singlet,  728 nm)
c          ba-alpha  (3S - 2P, triplet,  706 nm)
c          ba-alpha  (3P - 2S, singlet,  501 nm)
c          ba-alpha  (3D - 2P, singlet,  667 nm)
c          ba-beta   (4D - 2P, singlet,  492 nm)
c    for each of these lines there are separate contributions from
c    1) coupling to He
c    2) coupling to He+
c    0) total, sum over these 6 contributions
c
c
C
C  INPUT:
C          INIT: FLAG FOR INITIALISATION (DO NOT CHANGE!)
C          NCELL (COMPRT): INDEX IN TALLY ARRAYS FOR CURRENT ZONE
C          JJJ:    INDEX OF SEGMENT ALONG CHORD
C          ZDS:    LENGTH OF SEGMENT NO. JJJ
C          PEN:    CENTRAL wavelength OF LINE (nm)
C  OUTPUT: CONTRIB. FROM CELL NCELL AND CHORD SEGMENT JJJ TO:
C          THE H LINE FLUX PSIG(I),I=0,5 CONTRIBUTIONS
C          FROM ATOMS (only Ground state, MS unresolved), and BULK IONS
C          THE INTEGRAND ARGST IS SUCH THAT INTEGR.(ARGST*DL) = PSIG
C

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CESTIM
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      USE EIRMOD_COMPRT
      USE EIRMOD_COMUSR

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: INIT, JJJ
      REAL(DP), INTENT(IN) :: ZDS, DUMMY2, PEN
      REAL(DP), INTENT(IN OUT) :: PSIG(0:), ARGST(0:,:)
      INTEGER :: ISP, NCELC, ICELL
      LOGICAL :: LARGST
      CHARACTER(9) :: REAC1, REAC2
C
      SAVE
C
c     WRITE (IUNOUT,*) 'SIGHA,INIT,PEN,ISTRA ',
c    .                  INIT,PEN,ISTRA,ISTOLD,IITER,ITROLD

      LARGST = SIZE(ARGST,2) >= NSBOX

      IF (INIT.EQ.0) THEN
        PSIG=0.
        IF (LARGST) ARGST=0.

C  INITIALISE He-LINE ARRAYS FOR CURRENT STRATUM ?
        IF ((ISTRA .NE. ISTOLD) .OR. (IITER .NE. ITROLD) .OR.
     .      (PEN .NE. PENOLD) ) then
          if (PEN.EQ.728._DP) THEN
            write (iunout,*) ' He: 31S->21P '
            REAC1='2.2a     '   ! Amjuel population coeff. R1
            REAC2='2.3.2a   '   ! Amjuel population coeff. R0
            CALL EIRENE_HE_EMIS(ISTRA,REAC1,REAC2,
     .          NADVI+1,NADVI+2,NADVI+3)
          elseif (PEN.EQ.706._DP) THEN
            write (iunout,*) ' He: 33S->23P '
            REAC1='2.2b     '
            REAC2='2.3.2b   '
            CALL EIRENE_HE_EMIS(ISTRA,REAC1,REAC2,
     .          NADVI+1,NADVI+2,NADVI+3)
           elseif (PEN.EQ.501._DP) THEN
            write (iunout,*) ' He: 31P->21S '
            REAC1='2.2c     '
            REAC2='2.3.2c   '
            CALL EIRENE_HE_EMIS(ISTRA,REAC1,REAC2,
     .          NADVI+1,NADVI+2,NADVI+3)
          elseif (PEN.EQ.667._DP) THEN
            write (iunout,*) ' He: 31D->21P '
            REAC1='2.2d     '
            REAC2='2.3.2d   '
            CALL EIRENE_HE_EMIS(ISTRA,REAC1,REAC2,
     .          NADVI+1,NADVI+2,NADVI+3)
          elseif (PEN.EQ.492._DP) THEN
            write (iunout,*) ' He: 41D->21P '
            REAC1='2.2e     '
            REAC2='2.3.2e   '
            CALL EIRENE_HE_EMIS(ISTRA,REAC1,REAC2,
     .          NADVI+1,NADVI+2,NADVI+3)
          else
            WRITE (IUNOUT,*) 'NO LINE DEFINITION FOUND FOR PEN=',PEN
            WRITE (IUNOUT,*) 'SIGNAL IS SET TO 0'
            ADDV(NADVI+1:NADVI+3,:) = 0._DP
          endif
        endif   ! NEW INTERNAL ITERATION, OR NEW LINE, OR NEW STRATUM

        ISTOLD=ISTRA
        ITROLD=IITER
        PENOLD=PEN
        RETURN
      ENDIF
C
C  LINE INTEGRAL: PHOTONS/SEC/CM**2
C
!WZ:  This error message might be obsolete.
      IF (NSPZ.LT.1) THEN
        WRITE (iunout,*) 'ERROR EXIT FROM SIGHE '
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
C
      ncelc=ncltal(ncell)
      PSIG(1)=PSIG(1)+ZDS*ADDV(NADVI+1,NCELC)
      PSIG(2)=PSIG(2)+ZDS*ADDV(NADVI+2,NCELC)
      PSIG(0)=PSIG(0)+ZDS*ADDV(NADVI+3,NCELC)

      IF (LARGST) THEN
        ARGST(1,JJJ)=ADDV(NADVI+1,NCELC)
        ARGST(2,JJJ)=ADDV(NADVI+2,NCELC)
        ARGST(0,JJJ)=ADDV(NADVI+3,NCELC)
      END IF
C
      END

C     Following lines added for reinitialisation of eirene (DMH)

      SUBROUTINE EIRENE_SIGHE_REINIT
      ISTOLD = -1
      ITROLD = -1
      PENOLD = -1._DP
      END

      END MODULE EIRMOD_SIGHE
