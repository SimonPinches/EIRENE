C
      SUBROUTINE SRFCHK(VX,VY,VZ,SG,*)
c  RETURN 1   TRY ONCE AGAIN, WITH FRESH NCELL NUMBERS, SG
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMPRT
      USE EIRMOD_CCONA
      USE EIRMOD_CGEOM
      USE EIRMOD_CTRIG
      USE EIRMOD_CTETRA
      USE EIRMOD_CPOLYG
      USE EIRMOD_CGRID

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: VX,VY,VZ
      REAL(DP) :: SG, SH, PUX,PUY,PN, XOLD,YOLD
      INTEGER :: NRCELL_OLD,NPCELL_OLD,NTCELL_OLD, NTEST, ICO,
     .           IDUM, IFPB,
     .           EIRENE_LEARC1, EIRENE_LEARC2

      IF (NLSRFX) THEN

c  particle is exactly on one of the radial grid surfaces (MRSURF)
c  radial cell no. NRCELL may be wrong
c  check orientation of parallel motion relative to radial coordinate

        NRCELL_OLD=NRCELL

        select case (levgeo)
        case(1)
          SG=SIGN(1._DP,VX)
          IF (SG.LT.0) THEN
            NRCELL=MRSURF-1
          ELSEIF (SG.GT.0) THEN
            NRCELL=MRSURF
          ENDIF
        case(2)
          PUX= X0-EP1(MRSURF)
          PUY= Y0/ELL(MRSURF)/ELL(MRSURF)
          PN=SQRT(PUX*PUX+PUY*PUY+EPS60)
          PUX=PUX/PN
          PUY=PUY/PN
          SG=VX*PUX+VY*PUY
          IF (ABS(SG) .LT. EPS6) THEN
            NLSRFX=.FALSE.
            SH=SIGN(1._DP,SG)*CELDIA(NCELL)*1.D-2
            X0 = X0 + SH*PUX
            Y0 = Y0 + SH*PUY
          END IF
          IF (SG.LT.0) THEN
            NRCELL=NGHPLS(1,MRSURF,NPCELL)
          ELSEIF (SG.GT.0) THEN
            NRCELL=NGHPLS(3,MRSURF,NPCELL)
          ENDIF
        case (3)
          IFPB = 1
          XOLD = X0
          YOLD = Y0
          IDUM = NPCELL
          SG=VX*PLNX(MRSURF,NPCELL)+VY*PLNY(MRSURF,NPCELL)
          DO
            IF (ABS(SG) .LT. EPS6) THEN
              NLSRFX=.FALSE.
              SH=SIGN(1._DP,SG)*CELDIA(NCELL)*1.D-2
              X0 = XOLD + SH*PLNX(MRSURF,NPCELL)*IFPB
              Y0 = YOLD + SH*PLNY(MRSURF,NPCELL)*IFPB
            END IF
            IF (SG.LT.0) THEN
              NRCELL=NGHPLS(1,MRSURF,NPCELL)
            ELSEIF (SG.GT.0) THEN
              NRCELL=NGHPLS(3,MRSURF,NPCELL)
            ELSE
              NRCELL=EIRENE_LEARC1(X0,Y0,Z0,IDUM,MRSURF-1,MRSURF,
     .                             NLSRFX,NLSRFY,NPANU,'SRFCHK      ')
            ENDIF
            IF (NPCELL == IDUM) EXIT
            IFPB = -1
          END DO
        case (4)
          SG=VX*PTRIX(IPOLG,MRSURF)+
     .       VY*PTRIY(IPOLG,MRSURF)
          IF (ABS(SG) .LT. EPS6) THEN
            SH=SIGN(1._DP,SG)*CELDIA(NCELL)*1.D-2
            X0 = X0  +SH*PTRIX(IPOLG,MRSURF)
            Y0 = Y0  +SH*PTRIY(IPOLG,MRSURF)
            WRITE (IUNOUT,*) 'ON SURFACE IN FOLION, NPANU = ',NPANU
            WRITE (IUNOUT,*) 'AND MOVING PARALLEL TO SURFACE'
            WRITE (IUNOUT,*) 'PUSH INTO SUSPECTED NEXT CELL, SH = ',SH
            NLSRFX=.FALSE.
            IF (SG.GT.0.0_DP) THEN
c             NTEST=EIRENE_LEARC1(X0,Y0,Z0,IPOLG,1,NR1STM,
c    .                            NLSRFX,NLSRFY,NPANU,'FOLION      ')
c             if (ntest.ne.nrcell)
c    .           write (iunout,*) 'sg,ntest,nchbar ',
C    .                             SG,NTEST,NCHBAR(IPOLG,MRSURF)
              NRCELL=NCHBAR(IPOLG,MRSURF)
              IPOLG=NSEITE(IPOLG,MRSURF)
              MRSURF=NRCELL
            ENDIF
          ELSEIF (SG.GT.0.0_DP) THEN  !  SG IS GT EPS6
            NTEST=NCHBAR(IPOLG,MRSURF)
            IF (NTEST.EQ.0) THEN
C  NO NEIGHBOR. PUSH BACK INTO OLD CELL.
              SH=-CELDIA(NCELL)*1.D-2
              WRITE (IUNOUT,*) 'ON SURFACE IN FOLION, NPANU = ',NPANU
              WRITE (IUNOUT,*) 'PUSH BACK INTO OLD CELL: SH = ',SH
              WRITE (iunout,*) 'NRCELL = ',NRCELL
              NLSRFX=.FALSE.
c  strictly: particle should be pushed towards COM.
              X0 = X0  +SH*PTRIX(IPOLG,MRSURF)
              Y0 = Y0  +SH*PTRIY(IPOLG,MRSURF)
            ELSE
c  neighbor found. continue in neighbor cell.
              NRCELL=NTEST
              IPOLG=NSEITE(IPOLG,MRSURF)
              MRSURF=NRCELL
            ENDIF
          ELSEIF (SG.LT.0.0_DP) THEN ! SG IS LT.- EPS6
C  CONTINUE FLIGHT IN ORIGINAL CELL.
C  NOTHING TO BE DONE
          ENDIF
        case (5)
          SG=VX*PTETX(IPOLG,MRSURF)+
     .       VY*PTETY(IPOLG,MRSURF)+
     .       VZ*PTETZ(IPOLG,MRSURF)
          IF (ABS(SG) .LT. EPS6) THEN
C  TO BE WRITTEN
            WRITE (iunout,*) 'PARALLEL TO SURFACE IN FOLION ',NPANU
            WRITE (IUNOUT,*) 'CORRECTION FOR LEVGEO=5: TO BE DONE'
            CALL EIRENE_EXIT_OWN(1)
          ELSEIF (SG.GT.0) THEN
            NRCELL=NTBAR(IPOLG,MRSURF)
            IPOLG=NTSEITE(IPOLG,MRSURF)
            MRSURF=NRCELL
          ELSEIF (SG.LT.0) THEN
C  NOTHING TO BE DONE
          ENDIF
        case (10)
!PB EXPLICITLY ALLOW FOR LEVGEO=10
!PB NOTHING TO BE DONE
        case default
          write (iunout,*) 'levgeo in SRFCHK  ', levgeo
          write (iunout,*) 'option not ready, exit called'
          call EIRENE_exit_own(1)
        end select

        IF (NRCELL.NE.NRCELL_OLD) THEN
          ico=ico+1
          if (ico.le.1) return 1  ! GO BACK AND TRY AGAIN WITH NEW CELL NUMBER
        ENDIF


      ELSEIF (NLSRFY) THEN


c  particle is on one of the poloidal grid surfaces (MPSURF)
C  POLOIDAL CELL NO. NPCELL MAY BE WRONG
C  CHECK ORIENTATION OF PARALLEL MOTION RELATIV TO POLOIDAL COORDINATE
C
        NPCELL_OLD=NPCELL
        select case (LEVGEO)
        case (1)
          SG=SIGN(1._DP,VY)
          IF (SG.LT.0) THEN
            NPCELL=MPSURF-1
          ELSEIF (SG.GT.0) THEN
            NPCELL=MPSURF
          ENDIF
        case (2:3)
          SG=VX*PPLNX(NRCELL,MPSURF)+VY*PPLNY(NRCELL,MPSURF)
          IF (SG.LT.0) THEN
            npcell=nghpls(4,nrcell,mpsurf)
            ipolg=npcell
C  ACCOUNT FOR CUTS, PERIODICITY, ETC.
C           mpsurf is correct
          ELSEIF (SG.GT.0) THEN
            npcell=nghpls(2,nrcell,mpsurf)
            ipolg=npcell
C  ACCOUNT FOR CUTS, PERIODICITY, ETC.
            mpsurf=npcell
          ENDIF
        end select
        IF (NPCELL.NE.NPCELL_OLD) THEN
          ico=ico+1
          if (ico.le.1) return 1  ! GO BACK AND TRY AGAIN WITH NEW CELL NUMBER
        ENDIF


      ELSEIF (NLSRFZ) THEN


c  particle is on one of the toroidal grid surfaces (MTSURF)
C  TOROIDAL CELL NO. NTCELL MAY BE WRONG
C  CHECK ORIENTATION OF PARALLEL MOTION RELATIV TO POLOIDAL COORDINATE
C
        NTCELL_OLD=NTCELL
C  VLZPAR IS THE RELEVANT VELOCITY COMPONENT, BOTH FOR
C  NLTRZ AND NLTRT OPTION
        SG=SIGN(1._DP,VZ)
        IF (SG.LT.0) THEN
          NTCELL=MTSURF-1
        ELSEIF (SG.GT.0) THEN
          NTCELL=MTSURF
        ENDIF
        IF (NTCELL.NE.NTCELL_OLD) THEN
          ico=ico+1
          if (ico.le.1) return 1 ! GO BACK AND TRY AGAIN WITH NEW CELL NUMBER
        ENDIF

      ENDIF


      RETURN
      END
