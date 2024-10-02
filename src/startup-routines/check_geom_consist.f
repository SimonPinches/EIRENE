
      SUBROUTINE EIRENE_CHECK_GEOM_CONSIST
C
C   MODIFICATION OF INPUT DUE TO EITHER INCONSISTENCIES OR DUE
C   TO COUPLED NEUTRAL-PLASMA (OR NEUTRAL-NEUTRAL) CALCULATIONS
C   SOME FURTHER CONSTANTS ARE SET.      STATEM. NO. 2000 --> 4000

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CLOGAU
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      USE EIRMOD_CINIT
      USE EIRMOD_COMPRT, ONLY : IUNOUT

      IMPLICIT NONE

      INTEGER :: I

C
C
C   GEOMETRY, GRIDS
C
      IF (NLSLB) LEVGEO=1
      IF (NLCRC) LEVGEO=2
      IF (NLELL) LEVGEO=2
      IF (NLTRI) LEVGEO=2
      IF (NLPLG) LEVGEO=3
      IF (NLFEM) LEVGEO=4
      IF (NLTET) LEVGEO=5
      IF (NLGEN) LEVGEO=10

      IF (NLFEM) THEN
        NLPOL=.FALSE.
        IF (INDPRO(1).NE.3.AND.INDPRO(1).LT.5) THEN
          WRITE (iunout,*) ' PROFILE OPTION ',INDPRO(1),
     .                     ' NOT FORESEEN FOR'
          WRITE (iunout,*) ' FINITE ELEMENT MESH'
          WRITE (iunout,*)
     .      ' INDPRO(1) IS SET TO 3 <=> CONSTANT PROFILE'
          INDPRO(1)=3
        ENDIF
        IF (INDPRO(2).NE.3.AND.INDPRO(2).LT.5) THEN
          WRITE (iunout,*) ' PROFILE OPTION ',INDPRO(2),
     .                     ' NOT FORESEEN FOR'
          WRITE (iunout,*) ' FINITE ELEMENT MESH'
          WRITE (iunout,*)
     .      ' INDPRO(2) IS SET TO 3 <=> CONSTANT PROFILE'
          INDPRO(2)=3
        ENDIF
        IF (INDPRO(3).NE.3.AND.INDPRO(3).LT.5) THEN
          WRITE (iunout,*) ' PROFILE OPTION ',INDPRO(3),
     .                     ' NOT FORESEEN FOR'
          WRITE (iunout,*) ' FINITE ELEMENT MESH'
          WRITE (iunout,*)
     .      ' INDPRO(3) IS SET TO 3 <=> CONSTANT PROFILE'
          INDPRO(3)=3
        ENDIF
        IF (INDPRO(4).NE.3.AND.INDPRO(4).LT.5) THEN
          WRITE (iunout,*) ' PROFILE OPTION ',INDPRO(4),
     .                     ' NOT FORESEEN FOR'
          WRITE (iunout,*) ' FINITE ELEMENT MESH'
          WRITE (iunout,*)
     .      ' INDPRO(4) IS SET TO 3 <=> CONSTANT PROFILE'
          INDPRO(4)=3
        ENDIF
      ENDIF

      IF (NLTET) THEN
        NLPOL=.FALSE.
        NLTOR=.FALSE.
        IF (INDPRO(1).NE.3.AND.INDPRO(1).LT.5) THEN
          WRITE (iunout,*) ' PROFILE OPTION ',INDPRO(1),
     .                     ' NOT FORESEEN FOR'
          WRITE (iunout,*) ' TETRAHEDRON MESH'
          WRITE (iunout,*)
     .      ' INDPRO(1) IS SET TO 3 <=> CONSTANT PROFILE'
          INDPRO(1)=3
        ENDIF
        IF (INDPRO(2).NE.3.AND.INDPRO(2).LT.5) THEN
          WRITE (iunout,*) ' PROFILE OPTION ',INDPRO(2),
     .                     ' NOT FORESEEN FOR'
          WRITE (iunout,*) ' TETRAHEDRON MESH'
          WRITE (iunout,*)
     .      ' INDPRO(2) IS SET TO 3 <=> CONSTANT PROFILE'
          INDPRO(2)=3
        ENDIF
        IF (INDPRO(3).NE.3.AND.INDPRO(3).LT.5) THEN
          WRITE (iunout,*) ' PROFILE OPTION ',INDPRO(3),
     .                     ' NOT FORESEEN FOR'
          WRITE (iunout,*) ' TETRAHEDRON MESH'
          WRITE (iunout,*)
     .      ' INDPRO(3) IS SET TO 3 <=> CONSTANT PROFILE'
          INDPRO(3)=3
        ENDIF
        IF (INDPRO(4).NE.3.AND.INDPRO(4).LT.5) THEN
          WRITE (iunout,*) ' PROFILE OPTION ',INDPRO(4),
     .                     ' NOT FORESEEN FOR'
          WRITE (iunout,*) ' TETRAHEDRON MESH'
          WRITE (iunout,*)
     .      ' INDPRO(4) IS SET TO 3 <=> CONSTANT PROFILE'
          INDPRO(4)=3
        ENDIF
      ENDIF
C
      NP2ND=MIN0(N2ND,NP2ND)
      NT3RD=MIN0(N3RD,NT3RD)
      IF (NLTOR.AND.NLTRA) NTTRA=NT3RD
      NR1P2=0
      IF (NLPOL.OR.NLTOR) NR1P2=NR1ST
      NP2T3=0
      IF (NLTOR) NP2T3=NP2ND
      IF (.NOT.NLADD) NRADD=0
      IF (.NOT.NLMLT) NBMLT=1

C  DEFAULT: NO COARSE-GRAINING OR REFINING OF GRID
C           THESE SETTINGS MAY BE ALTERED LATER, E.G. IN IF0COP.
      DO I=1,NRAD
        NCLTAL(I)=I
      END DO
      NSURF=NR1ST*NP2ND*NT3RD*NBMLT
      NSBOX=NSURF+NRADD

      NR1TAL = NR1ST
      NP2TAL = NP2ND
      NT3TAL = NT3RD
      NRADD_TAL = NRADD
      NSURF_TAL = NSURF
      NSBOX_TAL = NSBOX

      RETURN
      END SUBROUTINE EIRENE_CHECK_GEOM_CONSIST
