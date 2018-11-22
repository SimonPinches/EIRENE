C
C
C
      SUBROUTINE EIRENE_GEOUSR
C
C   PREPARE/MODIFY GEOMETRICAL DATA, OUTSIDE STANDARD INPUT OPTIONS
C   CALLED FROM SUBR. INPUT, AFTER READING FORMATED INPUT FILE
C
C   ALSO BEST HERE: CASE SPECIFIC SPEED-UP OF GEOMETRICAL WORK
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD

      USE EIRMOD_COMPRT

      IMPLICIT NONE

      CHARACTER(80) :: ZEILE
      integer, save :: ifirst=0

cdr Carry out geousr only once. It may be called from interfacing: infcop.
cdr Then avoid second (default) call from input.f

      if (ifirst /= 0) return
      ifirst=1

      READ (IUNIN,'(A80)') ZEILE    !*** 15
      READ (IUNIN,'(A80)') ZEILE    !  geometry comment

      CALL EIRENE_UPPERCASE(ZEILE)

      IF (INDEX(ZEILE,'GENERAL') /= 0) THEN
        CALL EIRENE_GEOUSR_GENERAL
      ELSEIF (INDEX(ZEILE,'BIASED_GARCHING') /= 0) THEN
        CALL EIRENE_GEOUSR_BIASED_GARCHING
      ELSE
c  no card: "geometry comment" found in input block 14.
        BACKSPACE(IUNIN)
        CALL EIRENE_GEOUSR_BIASED
      END IF

      RETURN
      END


C
C
