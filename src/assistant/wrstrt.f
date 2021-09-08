!pb  27.11.06: open and close statements for fort.10 moved here
Cdr  15.10.14: variances of spectral tallies: names synchronized with other variance tallies,
c              range of spectra corrected: 0 -- nspc+1, rather than 1 -- nspc
c    Jan.  16: remove redundant PSGM
C.........................................................................................

cdr  SUBROUTINE WRSTRT:
cdr  write MC estimated tallies, per stratum, onto fort.10
cdr    (volume-averaged, surface-averaged, spectra, and their standard deviations)

cdr  SUBROUTINE RSTRT:
cdr  read MC estimated tallies, per stratum, onto fort.10
cdr    (volume-averaged, surface-averaged, spectra, and their standard deviations)
cdr     e.g. for printout, plotting, etc.. of results from specified strata

cdr  on input:  IG     :  number of stratum ISTRA
cdr             IG=0   :  sum over strata
cdr             TRCFLE :  print diagnostics
cpb  Dec. 2017: remove type SPECT_ARRAY, not needed in Fortran 2003

      SUBROUTINE EIRENE_WRSTRT(IG,NSTRAI,IESTM1,IESTM2,IESTM3,
     .                  TALLYV,TALLYS,TALLYL,
     .                  ISDVI1,STAT1,ISDVI2,STAT2,
     .                  ISDVC1,SIGC,ISDVC2,SIGCS,
     .                  ISPCI,TRCFLE)

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD, ONLY: EIRENE_SPECTRUM, IFOFF
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      IMPLICIT NONE

      TYPE(EIRENE_SPECTRUM), INTENT(INOUT) :: TALLYL(*)
      REAL(DP), INTENT(INOUT) :: TALLYV(*), TALLYS(*),
     .                           STAT1(*)
      REAL(DP), INTENT(INOUT) :: STAT2(*), SIGC(*), SIGCS(*)
      INTEGER, INTENT(IN) :: IG, NSTRAI, IESTM1, IESTM2, ISDVI1, ISDVI2,
     .                       ISDVC1, ISDVC2, 
     .                       IESTM3, ISPCI
      LOGICAL, INTENT(IN) :: TRCFLE

      INTEGER :: IMAX11, IMAX12, IMAX21, IMAX22, IMAX23, IMAX24, IMAX2,
     .           NRECL, IRC, ISTRA,
     .           JINI, J, JEND, IMAX, ISPC, IMAXS, NSPECI,NSPECE

C
C  WRITE DATA FOR SINGLE STRATA OR FROM SUM OVER STRATA ON TEMP. FILE FORT.10
C
      NRECL=1500
      IMAX11=IESTM1/NRECL+1
      IMAX12=IESTM2/NRECL+1
      IMAX21=ISDVI1/NRECL+1
      IMAX22=ISDVI2/NRECL+1
      IMAX23=ISDVC1/NRECL+1
      IMAX24=ISDVC2/NRECL+1
      IMAX2=IMAX21+IMAX22+IMAX23+IMAX24
 
      IMAXS=0
      DO ISPC=1,IESTM3
        IMAXS=IMAXS+1
C  SPECTRUM BINS RANGE FROM 0 TO NSPC+1
        IMAXS=IMAXS+(1+TALLYL(ISPC)%NSPC+1)/NRECL+1
        IF (ISPCI.NE.0) THEN
          IMAXS=IMAXS+4*((1+TALLYL(ISPC)%NSPC+1)/NRECL+1)
        END IF
      END DO
      IMAX=IMAX11+IMAX12+IMAX2+IMAXS
      ISTRA=IG
      IRC=ISTRA*IMAX+1
      IF (TRCFLE.AND.IG.NE.0) WRITE (iunout,*) 'WRITE STRATUM NO. ',IG
      IF (TRCFLE.AND.IG.EQ.0) WRITE (iunout,*) 'WRITE SUM OVER STRATA'
C
      OPEN (UNIT=10+ifoff,ACCESS='DIRECT',FORM='UNFORMATTED',
     .      RECL=8*NRECL,STATUS='UNKNOWN')

      JINI=1
      IF (TRCFLE) WRITE (iunout,*) 'ESTIMV'
    1 JEND=MIN0(JINI-1+NRECL,IESTM1)
      IF (TRCFLE.AND.(JINI.EQ.1.OR.JEND.EQ.IESTM1)) THEN
        WRITE (iunout,*) 'WRITE 10 IRC,JINI,JEND ',
     .                             IRC,JINI,JEND
      ENDIF
      WRITE (10+ifoff,REC=IRC) (TALLYV(J),J=JINI,JEND)
      IF (JEND.EQ.IESTM1) GOTO 12
      JINI=JEND+1
      IRC=IRC+1
      GOTO 1

   12 CONTINUE
      IF (TRCFLE) WRITE (iunout,*) 'ESTIMS'
      IRC=IRC+1
      JINI=1
   11 JEND=MIN0(JINI-1+NRECL,IESTM2)
      IF (TRCFLE.AND.(JINI.EQ.1.OR.JEND.EQ.IESTM2)) THEN
        WRITE (iunout,*) 'WRITE 10 IRC,JINI,JEND ',
     .                             IRC,JINI,JEND
      ENDIF
      WRITE (10+ifoff,REC=IRC) (TALLYS(J),J=JINI,JEND)
      IF (JEND.EQ.IESTM2) GOTO 2
      JINI=JEND+1
      IRC=IRC+1
      GOTO 11
C
    2 CONTINUE
      IF (TRCFLE) WRITE (iunout,*) 'STATIS 1'
      IRC=IRC+1
      JINI=1
    3 JEND=MIN0(JINI-1+NRECL,ISDVI1)
      IF (TRCFLE.AND.(JINI.EQ.1.OR.JEND.EQ.ISDVI1)) THEN
        WRITE (iunout,*) 'WRITE 10 IRC,JINI,JEND ',
     .                             IRC,JINI,JEND
      ENDIF
      WRITE (10+ifoff,REC=IRC) (STAT1(J),J=JINI,JEND)
      IF (JEND.EQ.ISDVI1) GOTO 21
      JINI=JEND+1
      IRC=IRC+1
      GOTO 3
C
   21 CONTINUE
      IF (TRCFLE) WRITE (iunout,*) 'STATIS 2'
      IRC=IRC+1
      JINI=1
   22 JEND=MIN0(JINI-1+NRECL,ISDVI2)
      IF (TRCFLE.AND.(JINI.EQ.1.OR.JEND.EQ.ISDVI2)) THEN
        WRITE (iunout,*) 'WRITE 10 IRC,JINI,JEND ',
     .                             IRC,JINI,JEND
      ENDIF
      WRITE (10+ifoff,REC=IRC) (STAT2(J),J=JINI,JEND)
      IF (JEND.EQ.ISDVI2) GOTO 23
      JINI=JEND+1
      IRC=IRC+1
      GOTO 22
C
   23 CONTINUE
      IF (TRCFLE) WRITE (iunout,*) 'STATIS 3'
      IRC=IRC+1
      JINI=1
   24 JEND=MIN0(JINI-1+NRECL,ISDVC1)
      IF (TRCFLE.AND.(JINI.EQ.1.OR.JEND.EQ.ISDVC1)) THEN
        WRITE (iunout,*) 'WRITE 10 IRC,JINI,JEND ',
     .                             IRC,JINI,JEND
      ENDIF
      WRITE (10+ifoff,REC=IRC) (SIGC(J),J=JINI,JEND)
      IF (JEND.EQ.ISDVC1) GOTO 25
      JINI=JEND+1
      IRC=IRC+1
      GOTO 24
C
   25 CONTINUE
      IF (TRCFLE) WRITE (iunout,*) 'STATIS 4'
      IRC=IRC+1
      JINI=1
   26 JEND=MIN0(JINI-1+NRECL,ISDVC2)
      IF (TRCFLE.AND.(JINI.EQ.1.OR.JEND.EQ.ISDVC2)) THEN
        WRITE (iunout,*) 'WRITE 10 IRC,JINI,JEND ',
     .                             IRC,JINI,JEND
      ENDIF
      WRITE (10+ifoff,REC=IRC) (SIGCS(J),J=JINI,JEND)
      IF (JEND.EQ.ISDVC2) GOTO 62
      JINI=JEND+1
      IRC=IRC+1
      GOTO 26
C
   62 CONTINUE

      IF (TRCFLE) WRITE (iunout,*) 'SPECTRA'
      DO ISPC=1,IESTM3
C  SET RANGE OF SPECTRUM ISPC, ADD BIN 0 AND NSPC+1 FOR LOW AND HIGH END OF SPECTRUM
        NSPECI=0
        NSPECE=TALLYL(ISPC)%NSPC+1
        IRC=IRC+1
        WRITE (10+ifoff,REC=IRC) TALLYL(ISPC)%SPCMIN,
     .                     TALLYL(ISPC)%SPCMAX,
     .                     TALLYL(ISPC)%SPCDEL,
     .                     TALLYL(ISPC)%SPCDELI,
     .                     TALLYL(ISPC)%SPCS,
     .                     TALLYL(ISPC)%SGMS,
     .                     TALLYL(ISPC)%STVS,
     .                     TALLYL(ISPC)%GGS,
     .                     TALLYL(ISPC)%NSPC,
     .                     TALLYL(ISPC)%ISPCTYP,
     .                     TALLYL(ISPC)%ISPCSRF,
     .                     TALLYL(ISPC)%IPRTYP,
     .                     TALLYL(ISPC)%IPRSP,
     .                     TALLYL(ISPC)%IMETSP
        DO JINI=NSPECI,NSPECE,NRECL
          IRC=IRC+1
          JEND=MIN(NSPECE, JINI+NRECL-1)
          WRITE (10+ifoff,REC=IRC)
     .      (TALLYL(ISPC)%SPC(J),J=JINI,JEND)
        END DO
        IF (ISPCI.NE.0) THEN
          DO JINI=NSPECI,NSPECE,NRECL
            IRC=IRC+1
            JEND=MIN(NSPECE, JINI+NRECL-1)
            WRITE (10+ifoff,REC=IRC)
     .        (TALLYL(ISPC)%SGM(J),J=JINI,JEND)
          END DO
          DO JINI=NSPECI,NSPECE,NRECL
            IRC=IRC+1
            JEND=MIN(NSPECE, JINI+NRECL-1)
            WRITE (10+ifoff,REC=IRC)
     .        (TALLYL(ISPC)%SDV(J),J=JINI,JEND)
          END DO
          DO JINI=NSPECI,NSPECE,NRECL
            IRC=IRC+1
            JEND=MIN(NSPECE, JINI+NRECL-1)
            WRITE (10+ifoff,REC=IRC)
     .        (TALLYL(ISPC)%STV(J),J=JINI,JEND)
          END DO
          DO JINI=NSPECI,NSPECE,NRECL
            IRC=IRC+1
            JEND=MIN(NSPECE, JINI+NRECL-1)
            WRITE (10+ifoff,REC=IRC)
     .        (TALLYL(ISPC)%GG(J),J=JINI,JEND)
          END DO
        END IF
      END DO

      CLOSE (UNIT=10+ifoff)
C
      RETURN
      END SUBROUTINE EIRENE_WRSTRT
C
      SUBROUTINE EIRENE_RSTRT(IG,NSTRAI,IESTM1,IESTM2,IESTM3,
     .            TALLYV,TALLYS,TALLYL,
     .            ISDVI1,STAT1,ISDVI2,STAT2,
     .            ISDVC1,SIGC,ISDVC2,SIGCS,
     .            ISPCI,TRCFLE)
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD, ONLY: EIRENE_SPECTRUM, IFOFF
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      IMPLICIT NONE

      TYPE(EIRENE_SPECTRUM), INTENT(INOUT) :: TALLYL(*)
      REAL(DP), INTENT(INOUT) :: TALLYV(*), TALLYS(*),
     .                         STAT1(*)
      REAL(DP), INTENT(INOUT) :: STAT2(*), SIGC(*), SIGCS(*)
      INTEGER, INTENT(IN) :: IG, NSTRAI, IESTM1, IESTM2, ISDVI1, ISDVI2,
     .                       ISDVC1, ISDVC2,  
     .                       IESTM3, ISPCI
      LOGICAL, INTENT(IN) :: TRCFLE

      INTEGER :: IMAX11, IMAX12, IMAX21, IMAX22, IMAX23, IMAX24, IMAX2,
     .           NRECL, IRC, ISTRA,
     .           JINI, J, JEND, IMAX, ISPC, IMAXS, NSPECI,NSPECE
C
C  READ DATA FOR SINGLE STRATA OR SUM OVER STRATA FROM TEMP. FILE FORT.10
C
      NRECL=1500
      IMAX11=IESTM1/NRECL+1
      IMAX12=IESTM2/NRECL+1
      IMAX21=ISDVI1/NRECL+1
      IMAX22=ISDVI2/NRECL+1
      IMAX23=ISDVC1/NRECL+1
      IMAX24=ISDVC2/NRECL+1
      IMAX2=IMAX21+IMAX22+IMAX23+IMAX24
 
      IMAXS=0
      DO ISPC=1,IESTM3
        IMAXS=IMAXS+1
C  SPECTRUM BINS RANGE FROM 0 TO NSPC+1
        IMAXS=IMAXS+(1+TALLYL(ISPC)%NSPC+1)/NRECL+1
        IF (ISPCI.NE.0) THEN
          IMAXS=IMAXS+4*((1+TALLYL(ISPC)%NSPC+1)/NRECL+1)
        END IF
      END DO
      IMAX=IMAX11+IMAX12+IMAX2+IMAXS
      ISTRA=IG
      IRC=ISTRA*IMAX+1
      IF (TRCFLE.AND.IG.NE.0) WRITE (iunout,*) 'READ STRATUM NO. ',IG
      IF (TRCFLE.AND.IG.EQ.0) WRITE (iunout,*) 'READ SUM OVER STRATA'

      OPEN (UNIT=10+ifoff,ACCESS='DIRECT',FORM='UNFORMATTED',
     .      RECL=8*NRECL,STATUS='OLD')

C
      JINI=1
      IF (TRCFLE) WRITE (iunout,*) 'ESTIMV'
   10 JEND=MIN0(JINI-1+NRECL,IESTM1)
      IF (TRCFLE.AND.(JINI.EQ.1.OR.JEND.EQ.IESTM1)) THEN
        WRITE (iunout,*) 'READ 10 IRC,JINI,JEND ',
     .                            IRC,JINI,JEND
      ENDIF
      READ (10+ifoff,REC=IRC) (TALLYV(J),J=JINI,JEND)
      IF (JEND.EQ.IESTM1) GOTO 15
      JINI=JEND+1
      IRC=IRC+1
      GOTO 10
C
   15 CONTINUE
      IF (TRCFLE) WRITE (iunout,*) 'ESTIMS'
      JINI=1
      IRC=IRC+1
   16 JEND=MIN0(JINI-1+NRECL,IESTM2)
      IF (TRCFLE.AND.(JINI.EQ.1.OR.JEND.EQ.IESTM2)) THEN
        WRITE (iunout,*) 'READ 10 IRC,JINI,JEND ',
     .                            IRC,JINI,JEND
      ENDIF
      READ (10+ifoff,REC=IRC) (TALLYS(J),J=JINI,JEND)
      IF (JEND.EQ.IESTM2) GOTO 20
      JINI=JEND+1
      IRC=IRC+1
      GOTO 16
C
   20 CONTINUE
      IF (TRCFLE) WRITE (iunout,*) 'STATIS 1'
      IRC=IRC+1
      JINI=1
   30 JEND=MIN0(JINI-1+NRECL,ISDVI1)
      IF (TRCFLE.AND.(JINI.EQ.1.OR.JEND.EQ.ISDVI1)) THEN
        WRITE (iunout,*) 'READ 10 IRC,JINI,JEND ',
     .                            IRC,JINI,JEND
      ENDIF
      READ (10+ifoff,REC=IRC) (STAT1(J),J=JINI,JEND)
      IF (JEND.EQ.ISDVI1) GOTO 31
      JINI=JEND+1
      IRC=IRC+1
      GOTO 30
C
   31 CONTINUE
      IF (TRCFLE) WRITE (iunout,*) 'STATIS 2'
      IRC=IRC+1
      JINI=1
   32 JEND=MIN0(JINI-1+NRECL,ISDVI2)
      IF (TRCFLE.AND.(JINI.EQ.1.OR.JEND.EQ.ISDVI2)) THEN
        WRITE (iunout,*) 'READ 10 IRC,JINI,JEND ',
     .                                      IRC,JINI,JEND
      ENDIF
      READ (10+ifoff,REC=IRC) (STAT2(J),J=JINI,JEND)
      IF (JEND.EQ.ISDVI2) GOTO 33
      JINI=JEND+1
      IRC=IRC+1
      GOTO 32
C
   33 CONTINUE
      IF (TRCFLE) WRITE (iunout,*) 'STATIS 3'
      IRC=IRC+1
      JINI=1
   34 JEND=MIN0(JINI-1+NRECL,ISDVC1)
      IF (TRCFLE.AND.(JINI.EQ.1.OR.JEND.EQ.ISDVC1)) THEN
        WRITE (iunout,*) 'READ 10 IRC,JINI,JEND ',
     .                            IRC,JINI,JEND
      ENDIF
      READ (10+ifoff,REC=IRC) (SIGC(J),J=JINI,JEND)
      IF (JEND.EQ.ISDVC1) GOTO 35
      JINI=JEND+1
      IRC=IRC+1
      GOTO 34
C
   35 CONTINUE
      IF (TRCFLE) WRITE (iunout,*) 'STATIS 4'
      IRC=IRC+1
      JINI=1
   36 JEND=MIN0(JINI-1+NRECL,ISDVC2)
      IF (TRCFLE.AND.(JINI.EQ.1.OR.JEND.EQ.ISDVC2)) THEN
        WRITE (iunout,*) 'READ 10 IRC,JINI,JEND ',
     .                            IRC,JINI,JEND
      ENDIF
      READ (10+ifoff,REC=IRC) (SIGCS(J),J=JINI,JEND)
      IF (JEND.EQ.ISDVC2) GOTO 66
      JINI=JEND+1
      IRC=IRC+1
      GOTO 36
C
   66 CONTINUE

      IF (TRCFLE) WRITE (iunout,*) 'SPECTRA'
      DO ISPC=1,IESTM3
        IRC=IRC+1
C  SET RANGE OF SPECTRUM ISPC, ADD BIN 0 AND NSPC+1 FOR LOW AND HIGH END OF SPECTRUM
        NSPECI=0
        NSPECE=TALLYL(ISPC)%NSPC+1
        READ (10+ifoff,REC=IRC) TALLYL(ISPC)%SPCMIN,
     .                    TALLYL(ISPC)%SPCMAX,
     .                    TALLYL(ISPC)%SPCDEL,
     .                    TALLYL(ISPC)%SPCDELI,
     .                    TALLYL(ISPC)%SPCS,
     .                    TALLYL(ISPC)%SGMS,
     .                    TALLYL(ISPC)%STVS,
     .                    TALLYL(ISPC)%GGS,
     .                    TALLYL(ISPC)%NSPC,
     .                    TALLYL(ISPC)%ISPCTYP,
     .                    TALLYL(ISPC)%ISPCSRF,
     .                    TALLYL(ISPC)%IPRTYP,
     .                    TALLYL(ISPC)%IPRSP,
     .                    TALLYL(ISPC)%IMETSP
        DO JINI=NSPECI,NSPECE,NRECL
          IRC=IRC+1
          JEND=MIN(NSPECE, JINI+NRECL-1)
          READ (10+ifoff,REC=IRC)
     .      (TALLYL(ISPC)%SPC(J),J=JINI,JEND)
        END DO
        IF (ISPCI.NE.0) THEN
          DO JINI=NSPECI,NSPECE,NRECL
            IRC=IRC+1
            JEND=MIN(NSPECE, JINI+NRECL-1)
            READ (10+ifoff,REC=IRC)
     .        (TALLYL(ISPC)%SGM(J),J=JINI,JEND)
          END DO
          DO JINI=NSPECI,NSPECE,NRECL
            IRC=IRC+1
            JEND=MIN(NSPECE, JINI+NRECL-1)
            READ (10+ifoff,REC=IRC)
     .        (TALLYL(ISPC)%SDV(J),J=JINI,JEND)
          END DO
          DO JINI=NSPECI,NSPECE,NRECL
            IRC=IRC+1
            JEND=MIN(NSPECE, JINI+NRECL-1)
            READ (10+ifoff,REC=IRC)
     .        (TALLYL(ISPC)%STV(J),J=JINI,JEND)
          END DO
          DO JINI=NSPECI,NSPECE,NRECL
            IRC=IRC+1
            JEND=MIN(NSPECE, JINI+NRECL-1)
            READ (10+ifoff,REC=IRC)
     .        (TALLYL(ISPC)%GG(J),J=JINI,JEND)
          END DO
        END IF
      END DO

      CLOSE (UNIT=10+ifoff)
C
      RETURN
      END SUBROUTINE EIRENE_RSTRT
