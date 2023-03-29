C
cdr  march 2019: bug fix: counter ICO was not transferred --> infinite loops possible
cdr  still not working properly, in case of trace ions: the discontinuity
cdr  of the B-field at cell boundaries (and projection of velocity onto
cdr  the B-field) may lead ambiguity wrt next cell and SG value
cdr  Much better performance results with INDPRO(5)=106 (smooth B-field)
C
      SUBROUTINE EIRENE_SRFCHK(VX,VY,VZ,SG,ICO,EPSLIM,IRET)
c  input: sg:     cosine of tracjetory with surface normal.
c         epslim: limiting value for SG, to identify motion parallel to surface.
c
C  newly added in March 2019: code snippet removed from FOLION.F
C                              similar code snippet still to be removed from LOCATE.F
CDR  CHECK FOR CORRECT ORIENTATION OF FLIGHT "SG" AND CELL NUMBERS,
CDR  IF A TEST PARTICLE STARTS FROM STITTING EXACTLY ON A SURFACE
CDR  |SG| LE EPSLIM IS TAKEN TO MEAN: MOTION PARALLEL TO A SURFACE

c  RETURN  IRET=1   TRY ONCE AGAIN, WITH FRESH NCELL NUMBERS, AND FRESH VALUE OF SG
c  ICO        =0 CALLED ON FIRST TRY, ICO SET IN CALLING PROGRAM
C  ICO        =1 CALLED ON SECOND TRY, WITH NEW GUESSES FOR SG, NEW CELL NUMBERS
C  ICO        >1 GIVE UP, TRY TO CONTINUE TRAJECTORY ANYWAY.
      USE EIRMOD_PRECISION, ONLY: DP
      USE EIRMOD_COMPRT, ONLY: IPOLG, IUNOUT, MPSURF, MRSURF, MTSURF,
     >                         NCELL, NLSRFX, NLSRFY, NLSRFZ, NPANU, 
     >                         NPCELL, NRCELL, NTCELL, X0, Y0, Z0 
      USE EIRMOD_CCONA, ONLY: EPS60
      USE EIRMOD_CGEOM, ONLY: CELDIA, NGHPLS
      USE EIRMOD_CTRIG, ONLY: NCHBAR, NSEITE, PTRIX, PTRIY
      USE EIRMOD_CTETRA, ONLY: NTBAR, NTSEITE, PTETX, PTETY, PTETZ
      USE EIRMOD_CPOLYG, ONLY: PLNX, PLNY, PPLNX, PPLNY
      USE EIRMOD_CGRID, ONLY: ELL, EP1, LEVGEO
      USE EIRMOD_LEARC1, ONLY: EIRENE_LEARC1

      IMPLICIT NONE

      INTEGER, INTENT(OUT) :: IRET
      REAL(DP), INTENT(IN) :: VX,VY,VZ
      REAL(DP), INTENT(INOUT) :: EPSLIM
      REAL(DP), INTENT(OUT) :: SG
      INTEGER, INTENT(INOUT) :: ICO
      REAL(DP) :: SH, PUX, PUY, PN, XOLD,YOLD
      INTEGER :: NRCELL_OLD,NPCELL_OLD,NTCELL_OLD, NTEST,
     .           IDUM, IFPB

      IRET = 0
      
      IF (NLSRFX) THEN

c  particle is exactly on one of the radial grid surfaces (MRSURF)
c  radial cell no. NRCELL may be wrong.
c  Check orientation of particle motion relative to radial coordinate

        NRCELL_OLD=NRCELL

        select case (levgeo)
        case(1)
cdr  celdia not yet defined. Tbd.
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
          IF (ABS(SG) .LT. EPSLIM) THEN
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
            IF (ABS(SG) .LT. EPSLIM) THEN
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
          IF (ABS(SG) .LT. EPSLIM) THEN
            SH=SIGN(1._DP,SG)*CELDIA(NCELL)*1.D-2
            X0 = X0  +SH*PTRIX(IPOLG,MRSURF)
            Y0 = Y0  +SH*PTRIY(IPOLG,MRSURF)
cdiag       WRITE (IUNOUT,*) 'ON SURFACE IN SRFCHK, NPANU = ',NPANU
cdiag       WRITE (IUNOUT,*) 'AND MOVING PARALLEL TO SURFACE'
cdiag       WRITE (IUNOUT,*) 'PUSH INTO SUSPECTED NEXT CELL, SH = ',SH
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
          ELSEIF (SG.GT.0.0_DP) THEN  ! SG IS GT EPSLIM (=EPS6 FOR ICO LE 1)
            NTEST=NCHBAR(IPOLG,MRSURF)
            IF (NTEST.EQ.0) THEN
C  NO NEIGHBOR. PUSH BACK INTO OLD CELL.
              SH=-CELDIA(NCELL)*1.D-2
              WRITE (IUNOUT,*) 'ON SURFACE IN SRFCHK, NPANU = ',NPANU
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
          ELSEIF (SG.LT.0.0_DP) THEN ! SG IS LT.- EPSLIM (=EPS6 FOR ICO LE 1)
C  CONTINUE FLIGHT IN ORIGINAL CELL.
C  NOTHING TO BE DONE
          ENDIF
        case (5)
          SG=VX*PTETX(IPOLG,MRSURF)+
     .       VY*PTETY(IPOLG,MRSURF)+
     .       VZ*PTETZ(IPOLG,MRSURF)
          IF (ABS(SG) .LT. EPSLIM) THEN
C  TO BE WRITTEN
            WRITE (iunout,*) 'PARALLEL TO SURFACE IN SRFCHK ',NPANU
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
          if (ico.le.1) then
! GO BACK AND TRY AGAIN WITH NEW CELL NUMBER
            iret = 1
            return
          elseif (ico.le.2) then
            epslim=epslim*10.0_DP
            iret = 1
            return  ! GO BACK AND TRY AGAIN WITH NEW CELL NUMBER
          end if
          ENDIF

      ELSEIF (NLSRFY) THEN


c  particle is on one of the poloidal grid surfaces (MPSURF)
C  POLOIDAL CELL NO. NPCELL MAY BE WRONG
C  CHECK ORIENTATION OF PARTICLE MOTION RELATIVE TO POLOIDAL COORDINATE
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
          if (ico.le.1) then      ! GO BACK AND TRY AGAIN WITH NEW CELL NUMBER
            iret = 1
            return
          elseif (ico.le.2) then
            epslim=epslim*10.0_DP
            iret = 1
            return  ! GO BACK AND TRY AGAIN WITH NEW CELL NUMBER
          end if
        ENDIF


      ELSEIF (NLSRFZ) THEN


c  particle is on one of the toroidal grid surfaces (MTSURF)
C  TOROIDAL CELL NO. NTCELL MAY BE WRONG
C  CHECK ORIENTATION OF PARALLEL MOTION RELATIVE TO POLOIDAL COORDINATE
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
          if (ico.le.1) then     ! GO BACK AND TRY AGAIN WITH NEW CELL NUMBER
            iret = 1
            return
          elseif (ico.le.2) then
            epslim=epslim*10.0_DP
            iret = 1
            return  ! GO BACK AND TRY AGAIN WITH NEW CELL NUMBER
          end if
        ENDIF

      ENDIF

      IF (ICO.LE.1) THEN
        IRET=0
        RETURN
      ENDIF

  999 CONTINUE
cdiag WRITE (IUNOUT,*) 'WARNING: ICO GE.3, INFINITE LOOP IN SRFCHK ?'
cdiag WRITE (IUNOUT,*) 'PARTICLE MOTION NEARLY WITHIN A SURFACE'
cdiag WRITE (IUNOUT,*) 'ICO, SG ', ICO, SG
cdiag WRITE (IUNOUT,*) 'TRY TO CONTINUE TRACK ANYWAY', NPANU, ICO

      IRET = 0
      RETURN
      END SUBROUTINE EIRENE_SRFCHK
