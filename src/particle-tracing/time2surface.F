      SUBROUTINE EIRENE_TIME_TO_STANDARD_SURFACE
     .           (TL, TF, TT, TS, ZDT1, ZT, ZTST,
     .            VELXS, VELYS, VELZS, VELS, ISRFCL, IRET, NEVNT)

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMPRT, ONLY : ITIME, ITYP, VELX, VELY, VELZ, VEL,
     .                          LGPART, MRSURF, IPOLGN, IUNOUT, NPANU,
     .                          NCELL, NRCELL, NPCELL, NTCELL, NACELL
      USE EIRMOD_CLOGAU
      USE EIRMOD_CCONA
      USE EIRMOD_CUPD, ONLY : NCOU, NUPC, CLPD, NCOUNT, NCOUNP
      USE EIRMOD_COMUSR, ONLY: NEVNTMAX, LVNMXTR
      use eirmod_timer
      use eirmod_timep
      use eirmod_cfplk, only : lcart, velpar, vlxpar, vlypar, vlzpar
      use eirmod_openmp

      IMPLICIT NONE
      REAL(DP), INTENT(INOUT) :: TL, TF, TT, TS, ZDT1, ZT, ZTST,
     .                           VELS, VELXS, VELYS, VELZS
      INTEGER, INTENT(INOUT) :: ISRFCL
      INTEGER, INTENT(OUT) :: IRET
      INTEGER, INTENT(IN) :: NEVNT
      REAL(DP) :: T
      EXTERNAL :: EIRENE_TIMET, EIRENE_LEER, EIRENE_MASAGE,
     .            EIRENE_MASJ4, EIRENE_MASR2
C
C  TS:   DISTANCE TO NEXT RADIAL SURFACE OF STANDARD MESH
C  ZDT1: DISTANCE TRAVELLED IN CURRENT RADIAL CELL
C  ZT:   ACCUMULATED DISTANCE, UNTIL THIS SEGMENT
C
C  USE PARALLEL VELOCITY, I.E., COMPUTE PARALLEL DISTANCES IN GRID
C  THUS ZT,TS,ZTST,ZDT1,CLPD ETC. ARE PARALLEL DISTANCES
C  I.E., LCART=F AT THIS POINT
C
      IRET = 0
      IF (ITIME.EQ.1) THEN
c  switch to gc velocity
        IF ((ITYP == 3) .AND. LCART) THEN
          VELXS=VELX
          VELYS=VELY
          VELZS=VELZ
          VELS =VEL
          VELX=VLXPAR
          VELY=VLYPAR
          VELZ=VLZPAR
          VEL =VELPAR
          LCART=.FALSE.
        ENDIF

        IF (NLRAD) THEN
          CALL EIRENE_TIMER(TS)
          IF (.NOT.LGPART) THEN
            IRET = 1
            IF (NEVNTMAX.GT.0.AND.NEVNT.GT.NEVNTMAX) GOTO 9913
            GOTO 9911
          END IF
C
          T=TS/TL-1.0_DP
          IF (ABS(T).LE.EPS10.AND.TL.NE.1.E30_DP) THEN
            IRET = 2
            GOTO 992
          END IF
          IF (TL.LT.TS.OR.TT.LT.TS.OR.TF.LT.TS) THEN
            MRSURF=0
            IPOLGN=0
C  CHECK FOR INTERSECTION WITH ADDITIONAL SURFACE
            IF (TL.LE.TT.AND.TL.LE.TF) THEN
              ZDT1=TL-ZT
              TL=ZT+ZDT1
              ZTST=TL
              ISRFCL=1
C  INTERSECTION WITH TIME SURFACE. TIME LIMIT REACHED ?
            ELSEIF (TT.LT.TL.AND.TL.LE.TF) THEN
              ZDT1=TT-ZT
              TT=ZT+ZDT1
              ZTST=TT
              ISRFCL=2
C  Fokker-Planck collision, DIFFUSIVE STEP
            ELSEIF (TF.LT.TL.AND.TF.LE.TT) THEN
              ZDT1=TF-ZT
              TF=ZT+ZDT1
              ZTST=TF
              ISRFCL=4
              IF (ITYP /= 3) THEN
                IRET = 3
                GOTO 994
              ENDIF
            ENDIF
          ELSE
C  INTERSECTION A  WITH 1-ST (RADIAL) GRID SURFACE
            ZDT1=TS-ZT
            ZTST=TS
            ISRFCL=0
          ENDIF
        ENDIF
C
        NCOU=1
        NUPC(1)=0
        CLPD(1)=ZDT1
        NCOUNT(1)=1
        NCOUNP(1)=1

C  CHECK SUB-GRIDS.  FOR OPTIONAL Y,Z RESOLUTION, ON 1D BACKGROUND MEDIUM
c  (can be switched off: nlpol, nltor, nltra)
C  sub-cells have the same background parameters as the parent (x-or-radial) cell,
C  i.e. same collision rates, same mean free path.

C  3RD Z (OR TOROIDAL) SUB-GRID, ALSO:  TOROIDAL PERIODICITY SURFACES
C  SUBDIVIDE GIVEN TRACK INTO Z (OR TOROIDAL) SMALLER SEGMENTS
        IF (NLTOR.OR.NLTRA) THEN
          CALL EIRENE_TIMET (ZDT1)
          TS=ZT+ZDT1
          ZTST=TS
        ENDIF
C  2ND (OR POLOIDAL) SUB-GRID
        IF (NLPOL) THEN
          CALL EIRENE_TIMEP(ZDT1)
          TS=ZT+ZDT1
          ZTST=TS
        ENDIF
C
        IF (ZDT1.LE.0.D0) THEN
          IRET = 3
          GOTO 990
        END IF

c  switch to full velocity but gc velocity is not saved
        IF ((ITYP == 3) .AND. .NOT.LCART) THEN
          VELX=VELXS
          VELY=VELYS
          VELZ=VELZS
          VEL =VELS
          LCART=.TRUE.
        ENDIF
C
      ELSEIF (ITIME.NE.1) THEN
C
        IF (NLTOR.OR.NLTRA) THEN
          CALL EIRENE_TIMET (ZDT1)
          TS=ZT+ZDT1
          ZTST=TS
        ENDIF

      ENDIF
      RETURN
C
  990 CONTINUE
cym
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP CRITICAL
#endif
      CALL EIRENE_LEER(1)
      CALL EIRENE_MASAGE('ERROR IN TIME_TO_STANDARD_SURFACE, '//
     .                   'ZDT1 OR NCELL OUT OF RANGE')
      CALL EIRENE_MASAGE('PARTICLE IS KILLED')

      write(iunout,*) 'ERROR for NPANU,thread =',NPANU,
     .                EIRENE_ITHREAD
      WRITE (iunout,*) 'ERROR NPANU,NCELL,ZDT1,ZTST,TL,TS '
      WRITE (iunout,'(I8,1X,I6,1P,4(1X,1E14.7))')
     .                  NPANU,NCELL,ZDT1,ZTST,TL,TS
      CALL EIRENE_MASJ4('NRCELL,NPCELL,NTCELL,NACELL     ',
     .                   NRCELL,NPCELL,NTCELL,NACELL)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP END CRITICAL
#endif
cym
      RETURN
C
  994 CONTINUE
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP CRITICAL
#endif
      CALL EIRENE_LEER(1)
      CALL EIRENE_MASAGE('ERROR IN TIME_TO_STANDARD_SURFACE, '//
     .                   'FOKKER-PLANCK COLLISION INDICATED BUT '//
     .                   'PARTICLE IS NO ION')
      CALL EIRENE_MASAGE('PARTICLE IS KILLED')

      write(iunout,*) 'ERROR for NPANU,thread =',NPANU,
     .                EIRENE_ITHREAD
      WRITE (iunout,*) 'ERROR NPANU,NCELL,ZDT1,ZTST,TL,TS,TT,TF '
      WRITE (iunout,'(I8,1X,I6,1P,6(1X,1E14.7))')
     .                  NPANU,NCELL,ZDT1,ZTST,TL,TS,TT,TF
      CALL EIRENE_MASJ4('NRCELL,NPCELL,NTCELL,NACELL     ',
     .                   NRCELL,NPCELL,NTCELL,NACELL)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP END CRITICAL
#endif
      RETURN
C
 9911 CONTINUE
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP CRITICAL
#endif
      CALL EIRENE_LEER(1)
      CALL EIRENE_MASAGE('ERROR IN TIME_TO_STANDARD_SURFACE, '//
     .                   'NO INTERSECTION FOUND')
      CALL EIRENE_MASAGE('PARTICLE IS KILLED')
      WRITE (iunout,*) 'NPANU,NCELL,NRCELL,NPCELL,NTCELL '
      WRITE (iunout,*)  NPANU,NCELL,NRCELL,NPCELL,NTCELL
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP END CRITICAL
#endif
      RETURN
C
 9913 CONTINUE
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP CRITICAL
#endif
      IF (LVNMXTR) THEN
        CALL EIRENE_LEER(1)
        CALL EIRENE_MASAGE('TIME_TO_STANDARD_SURFACE: '//
     .                     'STOPPING THE LONG TRACE')
        CALL EIRENE_MASAGE('MAXIMUM NUMBER OF EVENTS PER '//
     .                     'TRAJECTORY REACHED')
        CALL EIRENE_MASAGE('PARTICLE IS KILLED')
        WRITE (iunout,*) 'NPANU,NCELL,NRCELL,NPCELL,NTCELL '
        WRITE (iunout,*)  NPANU,NCELL,NRCELL,NPCELL,NTCELL
      ENDIF
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP END CRITICAL
#endif
      RETURN
C
  992 CONTINUE
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP CRITICAL
#endif
      CALL EIRENE_LEER(1)
      CALL EIRENE_MASAGE('ERROR IN TIME_TO_STANDARD_SURFACE, '//
     .                   ' SURFACE CONFLICT')
      CALL EIRENE_MASR2('TL,TS           ',TL,TS)
      WRITE (iunout,*) 'NPANU ',NPANU
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP END CRITICAL
#endif
      ZT=TL
      RETURN

      END SUBROUTINE EIRENE_TIME_TO_STANDARD_SURFACE
