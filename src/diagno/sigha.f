!pb  100107 ENTRY SIGHA_REINIT added
CDR  parameter PEN introduced, to identify hydrogen line by central energy
Cdr Aug.16:  The idenitifcation of particular lines 
cdr          by upper and lower energy level (input flags EMIN1,EMAX1 in block 12)
cdr          is not functional in this version, distinct from the manual description
cdr          currently lines can only be identified by their central energy PEN (EMIN1)
C
      SUBROUTINE EIRENE_SIGHA(INIT,JJJ,ZDS,PEN,PSIG,DUMMY2,ARGST)
CDR  this routine evaluates ("side on") hydrogen atom ("HA") emissivities,
cdr  integrated along a line of side (PSIG) and also the integrant resolved along 
cdr  line of side (ARGST)
c    currently there are up to 6 contributions to each particular preprogrammed
c    transition (depending on population coefficient data stored 
c    in file AMJUEL, section H.11 and H.12 
c  aug.16: available transitions in H-atom:
c          ly-alpha  (2 - 1)
c          ly-beta   (3 - 1)
c          ba-alpha  (3 - 2)
c          ba-beta   (4 - 2)
c          ba-gamma  (5 - 2)
c          ba-delta  (6 - 2)
c    for each of these lines there are separate contributions from
c    1) coupling to H
c    2) coupling to H+
c    3) coupling to H2
c    4) coupling to H2+
c    5) coupling to H-
c    6) coupling to H3+
c    0) total, sum over these 6 contributions
c
c
C
C  INPUT:
C          INIT: FLAG FOR INITIALISATION (DO NOT CHANGE!)
C          NCELL (COMPRT): INDEX IN TALLY ARRAYS FOR CURRENT ZONE
C          JJJ:    INDEX OF SEGMENT ALONG CHORD
C          ZDS:    LENGTH OF SEGMENT NO. JJJ
C          PEN:    CENTRAL ENERGY OF LINE (EV)
C  OUTPUT: PSIG:  LINE INTEGRAL OF EMISSION,I=0,6 CONTRIBUTIONS
C          ARGST: CONTRIB. FROM CELL NCELL AND CHORD SEGMENT JJJ TO:
C          THE H LINE FLUX PSIG(I),I=0,6 CONTRIBUTIONS
C          FROM ATOMS, MOLECULES, TEST IONS, BULK IONS AND NEGATIV IONS
C          THE INTEGRANT ARGST IS SUCH THAT INTEGR.(ARGST*DL) = PSIG
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
      REAL(DP), INTENT(IN OUT) :: PSIG(0:NSPZ+10), ARGST(0:NSPZ+10,NRAD)
      REAL(DP) :: PENOLD
      INTEGER :: ISTOLD, ISP, NCELC, ICELL, ITROLD
      DATA ISTOLD/-1/
      DATA ITROLD/-1/
      DATA PENOLD/-1._DP/
C
      SAVE
C
c     WRITE (IUNOUT,*) 'SIGHA,INIT,PEN,ISTRA ',
c    .                  INIT,PEN,ISTRA,ISTOLD,IITER,ITROLD
      IF (INIT.EQ.0) THEN
        DO 100 ISP=0,NSPZ+10
          PSIG(ISP)=0.
          DO 100 ICELL=1,NSBOX
            ARGST(ISP,ICELL)=0.
100     CONTINUE
C  INITIALISE ATOMIC H-LINE ARRAYS FOR CURRENT STRATUM ?
        IF ((ISTRA .NE. ISTOLD) .OR. (IITER .NE. ITROLD) .OR.
     .      (PEN .NE. PENOLD) ) then
          if (PEN.EQ.12.089_DP) THEN
            write (iunout,*) ' ly_beta '
            CALL EIRENE_Ly_beta 
     .          (ISTRA,NADVI+1,NADVI+2,NADVI+3,NADVI+4,NADVI+5,NADVI+6,
     .                 NADVI+7)
          elseif (PEN.EQ.10.2375_DP) THEN
            write (iunout,*) ' ly_alpha '
            CALL EIRENE_Ly_alpha 
     .          (ISTRA,NADVI+1,NADVI+2,NADVI+3,NADVI+4,NADVI+5,NADVI+6,
     .                 NADVI+7)
          elseif (PEN.EQ.3.0222_DP) THEN
            write (iunout,*) ' ba_delta '
            CALL EIRENE_Ba_delta
     .          (ISTRA,NADVI+1,NADVI+2,NADVI+3,NADVI+4,NADVI+5,NADVI+6,
     .                 NADVI+7)
           elseif (PEN.EQ.2.8560_DP) THEN
            write (iunout,*) ' ba_gamma '
            CALL EIRENE_Ba_gamma
     .          (ISTRA,NADVI+1,NADVI+2,NADVI+3,NADVI+4,NADVI+5,NADVI+6,
     .                 NADVI+7)
          elseif (PEN.EQ.2.5500_DP) THEN
            write (iunout,*) ' ba_beta '
            CALL EIRENE_Ba_beta
     .          (ISTRA,NADVI+1,NADVI+2,NADVI+3,NADVI+4,NADVI+5,NADVI+6,
     .                 NADVI+7)
          elseif (PEN.EQ.1.8889_DP) THEN 
            write (iunout,*) ' ba_alpha '
            CALL EIRENE_Ba_alpha
     .          (ISTRA,NADVI+1,NADVI+2,NADVI+3,NADVI+4,NADVI+5,NADVI+6,
     .                 NADVI+7)
          else
            WRITE (IUNOUT,*) 'NO LINE DEFINITION FOUND FOR PEN=',PEN
            WRITE (IUNOUT,*) 'SIGNAL IS SET TO 0'
            ADDV(NADVI+1:NADVI+7,:) = 0._DP
          endif
        endif
        ISTOLD=ISTRA
        ITROLD=IITER
        PENOLD=PEN
        RETURN
      ENDIF
C
C  LINE INTEGRAL: PHOTONS/SEC/CM**2
C
      IF (NSPZ+2.LT.6) THEN
        WRITE (iunout,*) 'ERROR EXIT FROM SIGHA '
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
C
      ncelc=ncltal(ncell)
      PSIG(1)=PSIG(1)+ZDS*ADDV(NADVI+1,NCELC)
      PSIG(2)=PSIG(2)+ZDS*ADDV(NADVI+2,NCELC)
      PSIG(3)=PSIG(3)+ZDS*ADDV(NADVI+3,NCELC)
      PSIG(4)=PSIG(4)+ZDS*ADDV(NADVI+4,NCELC)
      PSIG(5)=PSIG(5)+ZDS*ADDV(NADVI+5,NCELC)
      PSIG(6)=PSIG(6)+ZDS*ADDV(NADVI+6,NCELC)
      PSIG(0)=PSIG(0)+ZDS*ADDV(NADVI+7,NCELC)
      ARGST(1,JJJ)=ADDV(NADVI+1,NCELC)
      ARGST(2,JJJ)=ADDV(NADVI+2,NCELC)
      ARGST(3,JJJ)=ADDV(NADVI+3,NCELC)
      ARGST(4,JJJ)=ADDV(NADVI+4,NCELC)
      ARGST(5,JJJ)=ADDV(NADVI+5,NCELC)
      ARGST(6,JJJ)=ADDV(NADVI+6,NCELC)
      ARGST(0,JJJ)=ADDV(NADVI+7,NCELC)
C
      RETURN
 
C     Following lines added for reinitialisation of eirene (DMH)
 
      ENTRY EIRENE_SIGHA_REINIT
      ISTOLD = -1
      ITROLD = -1
      RETURN
      END
