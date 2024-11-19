C  june 05:  new option: K=3 (wg. iteravitve mode)
C  june 05:  name: calling routine
C
      SUBROUTINE EIRENE_PLNXTB(K,name)
c  prepare a next frame (a new picture) for eirene plotting options

c  PLNXTB is called from eirene routine "NAME.f"
C
C  K=1 CALL VOR EINEM BILD AUS EIGENER GRSOFTWARE
C  K=2 CALL VOR EINEM BILD MIT GRBLD  (KURVEF,....)
C  K=3 CALL NACH EINER VOLLEN ITERATION. RESET IFRST1,IFRST2
C                                        CLOSE PICTURES
C
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: K
      CHARACTER(LEN=*), INTENT(IN) :: NAME
      INTEGER, SAVE :: IFRST1=0, IFRST2=0
      EXTERNAL :: GRNXTF, GRSCLC

      WRITE (iunout,*) 'PLNXTB CALLED FROM ',NAME
      WRITE (iunout,*) 'K,IFRST1,IFRST2 ',K,IFRST1,IFRST2

      select case (K)
      case (1)
        IF (IFRST1.EQ.1) THEN
          CALL GRNXTF
        ELSE
          IFRST1=1
        ENDIF
      case (2)
        IF (IFRST1.EQ.0) THEN
          CALL GRSCLC (0.,0.0,39.5,28.7)
          CALL GRSCLC (3.,3.5,39.5,28.7)
          IFRST1=1
          IFRST2=1
        ELSE IF (IFRST2.EQ.0) THEN
          CALL GRNXTF
          CALL GRSCLC (0.,0.,39.5,28.7)
          CALL GRSCLC (3.,3.5,39.5,28.7)
          IFRST2=1
        ELSE
          CALL GRSCLC (3.,3.5,39.5,28.7)
        ENDIF
      case (3)
        IF (IFRST1.EQ.1) CALL GRNXTF
        IFRST1=0
        IF (IFRST2.EQ.1) THEN
          CALL GRSCLC (0.,0.,39.5,28.7)
          CALL GRSCLC (3.,3.5,39.5,28.7)
        ENDIF
        IFRST2=0
      end select
      RETURN
      END SUBROUTINE EIRENE_PLNXTB
