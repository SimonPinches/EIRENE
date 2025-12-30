cdr  Mar. 24: safety added: MASURF=0, for surface tallies
C
      SUBROUTINE EIRENE_TORCOL(IRET)

C  PARTICLE ON TOROIDAL PERIODICITY SURFACE MTSURF, FLIGHT INTO CELL NNTCLL
C  RETURN:  IRET = 1, NO SURFACE TALLIES, FLIGHT CONTINUES

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CUPD
      USE EIRMOD_CGRID
      USE EIRMOD_COMPRT
      USE EIRMOD_CLGIN
      USE EIRMOD_CFPLK
      USE EIRMOD_LEARC1, ONLY: EIRENE_LEARC1
      USE EIRMOD_PLT2D, ONLY: EIRENE_CHCTRC

      IMPLICIT NONE

      INTEGER, INTENT(OUT) :: IRET
      REAL(DP) :: VELX_OLD, SINROT,COSROT
      INTEGER :: ISTS, NN
      EXTERNAL :: EIRENE_EXIT_OWN
C
      IRET = 0
      TIME=TIME+ZT/VEL
C
      IF (NINCZ.EQ.0) THEN
        WRITE (iunout,*) 'ERROR IN TORCOL, NINCZ ?  '
        WRITE (iunout,*) 'NINCZ ',NINCZ
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
C
C
      IF (ITYP.LE.2) THEN
C  ANGLE FOR ROTATION IS 2*ALPHA
        COSROT=COSAL
        SINROT=SINAL*FLOAT(NINCZ)
        VELX_OLD=VELX
        VELX= VELX_OLD*COSROT +VELZ*SINROT
        VELZ=-VELX_OLD*SINROT +VELZ*COSROT
CNR     ADVANCE TO NEXT TOROIDAL CELL
CNR     For neutrals, use X01 computed in EIRENE_TIMET
CNR     which is already advanced
        X0=X01-RMTOR
CNR   For test ions, they follow B and VELX and VELZ are not modified,
CNR   advance as with other coordinates.
      ELSEIF (ITYP.EQ.3) THEN
        X0=X0+ZT*VELX
      ENDIF
CNR   Now advance Y and Z
      Y0=Y0+ZT*VELY
      Z01=-Z01
      Z0=Z01
C
      NLSRFX=.FALSE.
      NLSRFZ=.FALSE.
      NLSRFZ=.TRUE.
      IPERID=NNTCLL

!PB
!      IF ((NNTCLL==1) .AND. (MTSURF == 1)) PHI=ZSURF(1)
!      IF ((NNTCLL==NTTRAM) .AND. (MTSURF == NTTRA)) PHI=ZSURF(NTTRA)
      PHI = ZSURF(MTSURF)
C
      IF (LEVGEO.EQ.3.OR.(LEVGEO.EQ.2.AND.NLPOL)) THEN
cdr find index ipolg on polygonal ring  no. NRCELL
        IF (NRCELL.GT.0)
     .  NN=EIRENE_LEARC1(X0,Y0,Z0,IPOLG,NRCELL,NRCELL,.FALSE.,.FALSE.,
     .                           NPANU,'TORCOL 1')
      ENDIF
C
      IF (NLTRC) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP CRITICAL
#endif
        CALL EIRENE_CHCTRC(X0,Y0,Z0,16,11)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP END CRITICAL
#endif
      ENDIF
C
      IF (.NOT.NLTOR) THEN
C
        MSURF=0
        MASURF=0
        colflag = .true.
        IRET = 1
        RETURN
      ELSE
C
C  NLTOR=TRUE:
C
        MSURF=0
        MASURF=0
        IF (NRCELL.GT.0) NTCELL=NNTCLL
        ISTS=INMP3I(IRCELL,IPCELL,MTSURF)
        IF (ISTS.EQ.0) THEN
          IRET = 1
          RETURN
        END IF
        ZT=0.
        IRET = 0
        RETURN
      ENDIF
C
      END SUBROUTINE EIRENE_TORCOL
