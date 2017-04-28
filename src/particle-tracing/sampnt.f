!pb 22.03.07:  LEVGEO=6 --> LEVGEO=10
c   jet-2005, patch 1:  new arguments shwl and efwl in parameter list
c   06.08.15:  arguments added to vecusr
c   28.04.17:  minor cleanup, more info printed at error exit.
C
C       ..............................
C       .                            .
C       .  SOURCE SAMPLING ROUTINES  .
C       .                            .
C       ..............................
C
C     SUBROUTINE SAMPNT
C     SUBROUTINE SAMLNE
C     SUBROUTINE SAMSRF
C     SUBROUTINE SAMVOL
C
      SUBROUTINE EIRENE_SAMPNT (NLPT,TIWL,TEWL,DIWL,VXWL,VYWL,VZWL,
     .                        EFWL,SHWL,WEISPZ)

cdr  point source. identify the starting point coordinates (from IPOINT)
c  input:
c    NLPT = IPOINT (no. of sub-stratum, i.e. "point" source) 
c    ISTRA         (no. of stratum), via Common
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CINIT
      USE EIRMOD_CGRID
      USE EIRMOD_COMPRT
      USE EIRMOD_COMSOU
      USE EIRMOD_CTRIG
 
      IMPLICIT NONE
 
      REAL(DP), INTENT(OUT) :: TEWL, SHWL, TIWL(*), DIWL(*), EFWL(*),
     .                         VXWL(*), VYWL(*), VZWL(*), WEISPZ(*)
      INTEGER, INTENT(IN) :: NLPT
      REAL(DP) :: X01, CNORM, WINK
      INTEGER :: NT, EIRENE_LEARCA, EIRENE_LEARC1, EIRENE_LEARC2, 
     .           EIRENE_LEARCT, EIRENE_LEAUSR, IPOINT, JSPZ,
     .           IAUSR, IBUSR, IRUSR, IPUSR, ITUSR, IPLSTI, IPLSV
C
      IPOINT=NLPT
      IF (SORLIM(IPOINT,ISTRA).LT.0.D0) THEN
        CALL EIRENE_SM1USR(NLPT,X0,Y0,Z0,
     .              SORAD1(NLPT,ISTRA),SORAD2(NLPT,ISTRA),
     .              SORAD3(NLPT,ISTRA),SORAD4(NLPT,ISTRA),
     .              SORAD5(NLPT,ISTRA),SORAD6(NLPT,ISTRA),
     .              IRUSR,IPUSR,ITUSR,IAUSR,IBUSR,
     .              TIWL,TEWL,DIWL,VXWL,VYWL,VZWL,EFWL,SHWL,WEISPZ)
      ELSE
        X0=SORAD1(IPOINT,ISTRA)
        Y0=SORAD2(IPOINT,ISTRA)
        Z0=SORAD3(IPOINT,ISTRA)
      ENDIF
C
      IF (NLTRA) PHI=SORAD3(IPOINT,ISTRA)*DEGRAD
      MRSURF=0
      MPSURF=0
      MTSURF=0
      MASURF=0
      NLSRFX=.FALSE.
      NLSRFY=.FALSE.
      NLSRFZ=.FALSE.
      NLSRFA=.FALSE.

c  set 1st grid (radial) grid point
      NRCELL=1
      NACELL=0
      IPOLG=1
      IF (NRSOR(IPOINT,ISTRA).GT.0.OR.NASOR(IPOINT,ISTRA).GT.0) THEN
        NRCELL=NRSOR(IPOINT,ISTRA)
        NACELL=NASOR(IPOINT,ISTRA)
        IPOLG =NISOR(IPOINT,ISTRA)
      ELSEIF (NRSOR(IPOINT,ISTRA).EQ.0.AND.
     .        NASOR(IPOINT,ISTRA).EQ.0) THEN
C  find nrcell, ipolg automatically. 
        IF (LEVGEO.LE.4) THEN
          NRCELL=EIRENE_LEARC1(X0,Y0,Z0,IPOLG,1,NR1STM,.FALSE.,.FALSE.,
     .                         NPANU,'SAMPNT      ')
        ELSEIF (LEVGEO.EQ.5) THEN
          NRCELL=EIRENE_LEARCT(X0,Y0,Z0)
        ELSEIF (LEVGEO.EQ.10) THEN
          nrcell=EIRENE_leausr(x0,y0,z0)
        END IF
        IF (NRCELL.GT.0.AND.NRCELL.LT.NR1ST) then
          NACELL=0
        ELSE
          GOTO 991 
        ENDIF
      ELSE
        GOTO 991
      ENDIF

c  set 3rd grid (toroidal) grid point
      NTCELL=1
      IPERID=1
      IF (NLTOR.AND.NACELL.EQ.0) THEN
        IF (NLTRZ) THEN
          IF (NTSOR(NLPT,ISTRA).GT.0) THEN
C  NTCELL IS EXPLICITLY DEFINED BY INPUT VARIABLE NTSOR
            NTCELL=NTSOR(NLPT,ISTRA)
          ELSEIF (NTSOR(NLPT,ISTRA).EQ.0) THEN
C  NTCELL IS COMPUTED IN STANDARD MESH
            NTCELL=EIRENE_LEARCA(Z0,ZSURF,1,NT3RD,1,'SAMPNT      ')
          ELSE
            GOTO 991
          ENDIF
        ELSEIF (NLTRA) THEN
C  NTSOR NOT AVAILABLE FOR NLTRA OPTION
C  FIND Z0,NTCELL FROM X01,PHI
          NTCELL=EIRENE_LEARCA(PHI,ZSURF,1,NT3RD,1,'SAMPNT      ')
          IF (NTCELL.LE.0.OR.NTCELL.GT.NT3RDM) THEN
            WRITE (iunout,*) 'NTCELL OUT OF RANGE IN SAMPNT '
            WRITE (iunout,*) 'PHI,ZHALF ',PHI,ZHALF
            CALL EIRENE_EXIT_OWN(1)
          ENDIF
          X01=X0+RMTOR
          CALL EIRENE_FZRTRI(X0,Z0,NTCELL,X01,PHI,NTCELL)
          IPERID=NTCELL
        ELSEIF (NLTRT) THEN
          WRITE (iunout,*) 'NLTRT: TO BE WRITTEN IN SAMPNT '
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
      ELSEIF (.NOT.NLTOR.OR.NACELL.GT.0) THEN
        IF (NLTRA) THEN
C  FIND Z0, NT,  FROM X0,PHI
          NT=EIRENE_LEARCA(PHI,ZSURF,1,NTTRA,1,'SAMPNT      ')
          IF (NT.LE.0.OR.NT.GT.NTTRAM) THEN
            WRITE (iunout,*) 'NT OUT OF RANGE IN SAMPNT '
            WRITE (iunout,*) 'PHI,ZFULL ',PHI,ZFULL
            CALL EIRENE_EXIT_OWN(1)
          ENDIF
          X01=X0+RMTOR
          CALL EIRENE_FZRTRI(X0,Z0,NT,X01,PHI,NT)
          IPERID=NT
        ELSEIF (NLTRT) THEN
          WRITE (iunout,*) 'NLTRT: TO BE WRITTEN IN SAMPNT '
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
      ENDIF

c  set 2rd grid (poloidal) grid point
      NPCELL=1
      IF (NLPOL.AND.NACELL.EQ.0) THEN
        IF (NPSOR(NLPT,ISTRA).GT.0) THEN
C  NPCELL IS EXPLICITLY DEFINED BY INPUT VARIABLE NPSOR
          NPCELL=NPSOR(NLPT,ISTRA)
        ELSEIF (NPSOR(NLPT,ISTRA).EQ.0) THEN
C  NPCELL IS COMPUTED IN STANDARD MESH
          IF (LEVGEO.EQ.1) THEN
            NPCELL=EIRENE_LEARCA(Y0,PSURF,1,NP2ND,1,'SAMPNT')
          ELSEIF (LEVGEO.EQ.2) THEN
            IF (NLCRC) THEN
              WINK=MOD(ATAN2(Y0,X0)+PI2A-PSURF(1),PI2A)+PSURF(1)
              NPCELL=EIRENE_LEARCA(WINK,PSURF,1,NP2ND,1,'SAMPNT')
            ELSE
              NPCELL=EIRENE_LEARC2(X0,Y0,NRCELL,NPANU,'SAMPNT')
            ENDIF
          ELSEIF (LEVGEO.EQ.3) THEN
            NPCELL=IPOLG
          ELSE
            WRITE (iunout,*) 'ERROR EXIT FROM SAMPNT. NLPOL ',LEVGEO
            CALL EIRENE_EXIT_OWN(1)
          ENDIF
        ELSE
          GOTO 991
        ENDIF
      ENDIF
C
      NBLOCK=NBSOR(IPOINT,ISTRA)
      NBLOCK=MAX0(1,NBLOCK)
      NBLOCK=MIN0(NBLOCK,NBMLT)
C
      IF (NRCELL.GT.0) NACELL=0
      IF (NACELL.GT.0) NBLOCK=NBMLTP
      NBLCKA=NSTRD*(NBLOCK-1)+NACELL
      NCELL=NRCELL+((NPCELL-1)+(NTCELL-1)*NP2T3)*NR1P2+NBLCKA
C 
c  set local plasma parameters in cell of point source
c
      TEWL=TEIN(NCELL)
      DO 13 IPLS=1,NPLSI
        IPLSTI=MPLSTI(IPLS)
        IPLSV=MPLSV(IPLS)
        TIWL(IPLS)=TIIN(IPLSTI,NCELL)
        IF (INDPRO(4) == 8) THEN
          CALL EIRENE_VECUSR (2,NCELL,X0,Y0,Z0,
     .         VXWL(IPLS),VYWL(IPLS),VZWL(IPLS),IPLS,.TRUE.)
        ELSE
          VXWL(IPLS)=VXIN(IPLSV,NCELL)
          VYWL(IPLS)=VYIN(IPLSV,NCELL)
          VZWL(IPLS)=VZIN(IPLSV,NCELL)
        END IF
        DIWL(IPLS)=DIIN(IPLS,NCELL)
13    CONTINUE
C
      DO 20 JSPZ=1,NSPZ
        WEISPZ(JSPZ)=-1.
20    CONTINUE
      IF (NSPEZ(ISTRA).LE.0) THEN
C  ANALOG SPECIES SAMPLING DISTRIBUTION NOT AVAILABLE FOR POINT SOURCE
        GOTO 992
      ENDIF
C
C  set reference directional unit vector 
C
      CRTX=SORAD4(IPOINT,ISTRA)
      CRTY=SORAD5(IPOINT,ISTRA)
      CRTZ=SORAD6(IPOINT,ISTRA)
      CNORM=SQRT(CRTX**2+CRTY**2+CRTZ**2)+EPS60
      CRTX=CRTX/CNORM
      CRTY=CRTY/CNORM
      CRTZ=CRTZ/CNORM
C
      RETURN
991   CONTINUE
      WRITE (iunout,*) 'ERROR IN SAMPNT , IPOINT, ISTRA ', IPOINT,ISTRA
      WRITE (iunout,*) 'NRCELL,NACELL,IPOLG ',NRCELL,NACELL,IPOLG
      WRITE (iunout,*) 'NPCELL              ',NPCELL
      WRITE (iunout,*) 'NTCELL,IPERID       ',NTCELL,IPERID
      CALL EIRENE_EXIT_OWN(1)
992   CONTINUE
      WRITE (iunout,*) 'ERROR IN SAMPNT, NSPEZ OUT OF RANGE           '
      CALL EIRENE_EXIT_OWN(1)
      END
