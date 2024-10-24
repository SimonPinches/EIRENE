      MODULE EIRMOD_PLT2D
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CLOGAU
      USE EIRMOD_CPLOT
      USE EIRMOD_CTRCEI
      USE EIRMOD_CTETRA
      USE EIRMOD_COMPRT
      USE EIRMOD_CFPLK
      USE EIRMOD_COMSOU
      USE EIRMOD_CCONA
      USE EIRMOD_CUPD
      IMPLICIT NONE
      PRIVATE

      PUBLIC :: EIRENE_PLT2D, EIRENE_CHCTRC, EIRENE_PLT2D_REINIT

      INTEGER,PARAMETER :: NTXHST=21
      CONTAINS

      subroutine EIRENE_plt2d
      IMPLICIT NONE

      RETURN
      END SUBROUTINE EIRENE_PLT2D

c--------------------------------------------------------------

C
C  PRINT AND PLOT PARTICLE HISTORIES IN GEOMETRY PLOT
C
      SUBROUTINE EIRENE_CHCTRC(XPLO,YPLO,ZPLO,IFLAG,ISYM)

cdr  this routine is identical to the subroutine chctrc(...) inside eirene routine plt2d.
cdr  It is kept as separate routine here, in case no further eirene default plotting routines are used,
cdr  to still be able to provide printed trajectory output.

      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: XPLO, YPLO, ZPLO
      INTEGER, INTENT(IN) :: IFLAG, ISYM
      INTEGER :: ISTR, I
      LOGICAL :: LWR
      CHARACTER(20) :: TXTHST(NTXHST)
      EXTERNAL :: EIRENE_LEER, EIRENE_MASJ1, EIRENE_MASJ1R,
     .            EIRENE_MASJ3, EIRENE_MASJ4, EIRENE_MASR1,
     .            EIRENE_MASR2, EIRENE_MASR3, EIRENE_MASR5,
     .            EIRENE_MASR6

      DATA TXTHST
     .           /'LOCATE(1)           ',
     .            'ELECTR. IMPACT(2)   ',
     .            'HEAVY PAR. IMPACT(3)',
     .            'PHOTON IMPACT(4)    ',
     .            'ELASTIC COLL.(5)    ',
     .            'CHARGE EXCHANGE(6)  ',
     .            'FOKKER-PLANCK(7)    ',
     .            'SURFACE(8)          ',
     .            'SPLITTING(9)        ',
     .            'RUSSIAN ROULETTE(10)',
     .            'PERIODICITY(11)     ',
     .            'RESTART:A. SPLT.(12)',
     .            'SAVE:COND. EXP.(13) ',
     .            'RESTART:COND EXP(14)',
     .            'TIME LIMIT(15)      ',
     .            'GENERATION LIMIT(16)',
     .            'FLUID LIMIT(17)     ',
     .            'ERROR DETECTED      ',     ! SYMBOL FOR PARTICLE
                                              ! TRACING ERROR.
c  next symbols/text: only for printout, not on plot.
     .            'INT.GRID SURFACE(19)',
     .            'DIFFUSION STEP(20)  ',
     .            'STATIC LOOP(21)     '/
C
C  WRITE TRACK DATA
C
      LWR = .TRUE.
      DO I = 1, 8
        IF (-ISYM == ISYPLT(I)) LWR = .FALSE.
      END DO
      IF (TRCHST .AND. LWR) THEN
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) TXTHST(ISYM)
        IF (ISPZ.GT.0.AND.ISPZ.LE.NSPZ) THEN
          IF (ISYM.EQ.1.OR..NOT.NLTRC)
     .      CALL EIRENE_MASJ1('NPANU   ',NPANU)
          WRITE (iunout,'(1X,A8)') TEXTS(ISPZ)
        ELSE
          WRITE (iunout,'(1X,A14,1X,I6)') 'LINE OF SIGHT ',NPANU
        ENDIF
        CALL EIRENE_MASJ4 ('ITIME,IFPATH,IUPDTE,ICOL        ',
     .               ITIME,IFPATH,IUPDTE,ICOL)
        CALL EIRENE_MASR3 ('X0,Y0,Z0                ',XPLO,YPLO,ZPLO)
C  FOR TRACE IONS: VELOCITY IS EITHER CARTESIAN (LCART) OR THE REDUCED (GC) VELOCITY
        IF (ITYP.EQ.3) THEN
          IF (LCART) THEN
            CALL EIRENE_MASR5
     .           ('VELX,VELY,VELZ,VEL,E0                   ',
     .             VELX,VELY,VELZ,VEL,E0)
          ELSE
            CALL EIRENE_MASR6
     .           ('VLXPAR,VLYPAR,VLZPAR,VELPAR,E0PAR,E0             ',
     .             VLXPAR,VLYPAR,VLZPAR,VELPAR,E0PAR,E0)
          ENDIF
        ELSE
C  FOR NEUTRALS OR PHOTONS: VELOCITY IS ALWAYS GIVEN BY THE CARTESIAN COMPONENTS
          CALL EIRENE_MASR5
     .         ('VELX,VELY,VELZ,VEL,E0                   ',
     .           VELX,VELY,VELZ,VEL,E0)
        ENDIF
        CALL EIRENE_MASR2 ('WEIGHT,TIME     ',WEIGHT,TIME)
        CALL EIRENE_MASR1 ('XGENER  ',XGENER)
        CALL EIRENE_MASJ3
     .  ('NRCELL,NACELL,NBLOCK    ',NRCELL,NACELL,NBLOCK)
        IF (NLPLG) CALL EIRENE_MASJ1('IPOLG   ',IPOLG)
        IF (NLTRA) THEN
          CALL EIRENE_MASJ1('IPERID  ',IPERID)
          CALL EIRENE_MASJ1R ('NNTCLL,PHI      ',NNTCLL,PHI/DEGRAD)
        ENDIF
        IF (NLPOL) THEN
          CALL EIRENE_MASJ1 ('NPCELL  ',NPCELL)
        ENDIF
        IF (NLTOR) THEN
          CALL EIRENE_MASJ1 ('NTCELL  ',NTCELL)
        ENDIF
        IF (NLTET.AND.(IPOLGN>0).AND.(MRSURF>0)) THEN
          CALL EIRENE_MASJ1 ('INMTIT  ',INMTIT(IPOLGN,MRSURF))
        END IF
C  CALLED FROM DIAGNO, WITH ISTRA=0?
        IF (ISTRA.EQ.0) THEN
          ISTR=1
        ELSE
          ISTR=ISTRA
        ENDIF
C  SURFACE EVENT
        IF ((ISYM.GE.8.AND.ISYM.LE.11.OR.ISYM.EQ.19).OR.
     .      (ISYM.EQ.1.AND.NLSRF(ISTR))) THEN
          IF (NLSRFX) THEN
            CALL EIRENE_MASJ1 ('MRSURF  ',MRSURF)
          ELSEIF (NLSRFY) THEN
            CALL EIRENE_MASJ1 ('MPSURF  ',MPSURF)
          ELSEIF (NLSRFZ) THEN
            CALL EIRENE_MASJ1 ('MTSURF  ',MTSURF)
          ELSEIF (NLSRFA) THEN
            CALL EIRENE_MASJ1 ('MASURF  ',MASURF)
          ENDIF
C         CALL MASJ4 ('MRSURF,MPSURF,MTSURF,MASURF     ',
C    .                 MRSURF,MPSURF,MTSURF,MASURF)
          CALL EIRENE_MASR1 ('SCOS    ',SCOS)
        ENDIF
      ENDIF

      RETURN
      END SUBROUTINE EIRENE_CHCTRC

C     following SUBROUTINE is for reinitialization of EIRENE (DMH)
      SUBROUTINE EIRENE_PLT2D_REINIT
      IMPLICIT NONE
      RETURN
      END SUBROUTINE EIRENE_PLT2D_REINIT

      END MODULE EIRMOD_PLT2D
