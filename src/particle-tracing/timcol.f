C  SEPT 05:  IN CASE OF TEST IONS, VEL IS THE PARALLEL VELOCITY ONLY
C            THIS HAS STILL TO BE TAKEN INTO ACCOUNT WHEN STORING AND SAMPLING
C            THE CENSUS ARRAY
!pb 08.11.06: definition of CENSUS arrays changed
!             RPART (NPRNL,1:NPARTT) --> RPART (1:NPARTT,NPRNL)
!             IPART (NPRNL,1:MPARTT) --> IPART (1:MPARTT,NPRNL)
!             RPARTC(NPRNL,1:NPARTT) --> RPARTC(1:NPARTT,NPRNL)
!             IPARTC(NPRNL,1:MPARTT) --> IPARTC(1:MPARTT,NPRNL)
cdr Jan 2016 : comments,  and: stop scoring census not only after total number
cdr            of allowed census scores is reached,
cdr            but instead do so also for each stratum, and for the scores per stratum limit.

!pb   SUBROUTINE EIRENE_TIMCOL (PR,*,*)
      SUBROUTINE EIRENE_TIMCOL (PR,IRET)
C
C  COLLISION WITH "TIME SURFACE", FIND NEW COORDINATES
C  UPDATE (TIME-) SURFACE TALLIES
C  UPDATE USER-SUPPLIED SNAPSHOT-ESTIMATED TALLIES (CALL UPNUSR)
C  PUT PARTICLE ONTO CENSUS ARRAYS
C  AND EITHER STOP HISTORY OR CONTINUE

C  RETURN 1: IRET = 1, CONTINUE FLIGHT
C  RETURN 2: IRET = 2, STOP FLIGHT
C
      USE EIRMOD_PRECISION, ONLY: DP
      USE EIRMOD_PARMMOD, ONLY: MPARTT, NLIM, NPARTT, NPRNL
      USE EIRMOD_COMUSR, ONLY: ISPEZ, NSNVI
      USE EIRMOD_CESTIM, ONLY: LEOTPHT, LEOTAT, LEOTIO, LEOTML,
     >                         LPOTPHT, LPOTAT, LPOTIO, LPOTML,
     >                         LSPUMP,
     >                         EOTPHT, EOTAT, EOTIO, EOTML,
     >                         POTPHT, POTAT, POTIO, POTML,
     >                         SPUMP
      USE EIRMOD_CCONA, ONLY: PI2A
      USE EIRMOD_CLOGAU, ONLY: NLMOVIE, NLTRA
      USE EIRMOD_CUPD, ONLY: NNTCLL, X00, X01, Y00, Z00, Z01
      USE EIRMOD_CGRID, ONLY: RMTOR
      USE EIRMOD_COMPRT, ONLY: IATM, IION, IMOL, IPHOT, IPLS, ISPZ,
     >                         ITYP, IPERID, IPOLG, IPOLGN, IPSTT,
     >                         ISTRA, E0, LGLAST, MSURF, MSURFG,
     >                         MASURF, MRSURF, MPSURF, MTSURF, NLSRFX,
     >                         NLSRFY, NLSRFZ, NLTRC, NPANU, PHI,
     >                         RPSTT, TIME, TT, VEL, VELX, VELY, VELZ,
     >                         WEIGHT, X0, Y0, Z0,
     >                         IUNOUT     
      USE EIRMOD_COMNNL, ONLY: IPART, IPRNLI, IPRNLS, ITMSTP, NPRNLS,
     >                         NTMSTP, RPART, TIME0
      USE EIRMOD_CLGIN, ONLY: NSTSI
      USE EIRMOD_CSDVI, ONLY: LMETSPW
      USE EIRMOD_PLT2D, ONLY: EIRENE_CHCTRC

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: PR
      INTEGER, INTENT(OUT) :: IRET
      INTEGER  :: IND
      REAL(DP) :: DIST, WGHTSG
C
      IRET = 0
      X0=X0+VELX*TT
      Y0=Y0+VELY*TT
      Z0=Z0+VELZ*TT
      TIME=TIME+TT/VEL
      MSURF=0
      NLSRFX=.FALSE.
      NLSRFY=.FALSE.
      NLSRFZ=.FALSE.
      MRSURF=0
      MPSURF=0
      MTSURF=0
      MASURF=0
C
      IPOLG=IPOLGN
      IPERID=NNTCLL
      IF (NLTRA) THEN
        PHI=MOD(PHI-ATAN2(Z01,X01)+ATAN2(Z0,(RMTOR+X0)),PI2A)
        X01=X0+RMTOR
      ENDIF
      X00=X0
      Y00=Y0
      Z00=Z0
      Z01=Z0
      WEIGHT=WEIGHT*PR
C
      IF (NLTRC) CALL EIRENE_CHCTRC(X0,Y0,Z0,16,15)
C
C  UPDATE SNAPSHOT ESTIMATORS
      IF (NSNVI.GT.0) CALL EIRENE_UPNUSR
C
cdpc
CDR:  this must be generalized, towards a more general horizon
CDR   rather than fixed horizon at 100 meters in x-y plane
      dist=sqrt(x0**2+y0**2)
      if(dist.gt.1e4) then
        write(iunout,*) 'timcol: ERROR!  dist = ',dist,
     1   ' (particle more than 100 m from the origin)'
        write(iunout,*) 'npanu,x0,y0,z0,velx,vely,velz ',
     1   npanu,x0,y0,z0,velx,vely,velz
        weight=0.
        goto 112
      endif
cdpc
C
C  TOTAL NO. OF SCORES ON CENSUS
      IPRNLI=IPRNLI+1
C  NO. OF SCORES ON CENSUS FOR PRESENT STRATUM ISTRA
      IPRNLS=IPRNLS+1
      IF (IPRNLS.GE.NPRNLS(ISTRA).AND..NOT.NLMOVIE) THEN
C  THIS IS THE LAST SCORE FOR THIS STRATUM TO BE STORED
        LGLAST=.TRUE.
      ENDIF
C
C   CENSUS ARRAYS:
C   SAVE LOCATION, WEIGHT AND OTHER PARAMETERS
C   STOP SCORING ON CENSUS AFTER NPRNL SCORES TOTAL

CDR STOP ALSO AFTER NPRNLS SCORES FOR STRATUM ISTRA ??
      if (iprnli <= nprnl.and.iprnls <= nprnls(istra)) then
cdr   if (iprnli <= nprnl) then

        RPART(1:NPARTT,IPRNLI)=RPSTT(1:NPARTT)
        IPART(1:MPARTT,IPRNLI)=IPSTT(1:MPARTT)
      end if

C  DO NOT SCORE ON CENSUS ANYMORE FOR THIS STRATUM
      if (iprnls > nprnls(istra)) iprnls = nprnls(istra)
      if (iprnli > nprnl)         iprnli = nprnl

C
  112 continue

C  DECIDE: CONTINUE OR STOP TRAJECTORY
      IF (NTMSTP.GE.0.AND.ITMSTP.GE.NTMSTP) THEN
C
C  DO NOT CONTINUE THIS TRACK
C  UPDATE PARTICLE EFFLUX ONTO TIME SURFACE MSURF=NLIM+NSTSI
C  UPDATE ENERGY FLUX ONTO TIME SURFACE MSURF=NLIM+NSTSI
C  THEN STOP HISTORY
C
        MSURF=NLIM+NSTSI
cdr  to replace cdr out ini -- cdr out end code below with call to update_surface.
cdr  Still to be tested first...
cdr     ITYP_OLD=ITYP
        MSURFG=0
        WGHTSG=WEIGHT
        IND=1
cdr     CALL EIRENE_UPDATE_SURFACE (ITYP_OLD,WGHTSG,IND)
cdr out ini
        IF (ITYP.EQ.0) THEN
          IF (LEOTPHT) EOTPHT(IPHOT,MSURF)=EOTPHT(IPHOT,MSURF)+E0*WEIGHT
          IF (LPOTPHT) POTPHT(IPHOT,MSURF)=POTPHT(IPHOT,MSURF)+WEIGHT
        ELSEIF (ITYP.EQ.1) THEN
          IF (LEOTAT) EOTAT(IATM,MSURF)=EOTAT(IATM,MSURF)+E0*WEIGHT
          IF (LPOTAT) POTAT(IATM,MSURF)=POTAT(IATM,MSURF)+WEIGHT
        ELSEIF (ITYP.EQ.2) THEN
          IF (LEOTML) EOTML(IMOL,MSURF)=EOTML(IMOL,MSURF)+E0*WEIGHT
          IF (LPOTML) POTML(IMOL,MSURF)=POTML(IMOL,MSURF)+WEIGHT
        ELSEIF (ITYP.EQ.3) THEN
          IF (LEOTIO) EOTIO(IION,MSURF)=EOTIO(IION,MSURF)+E0*WEIGHT
          IF (LPOTIO) POTIO(IION,MSURF)=POTIO(IION,MSURF)+WEIGHT
        ENDIF
cdr out end
        ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
c spatial resolution on time-surface is not available. MSURFG ?
        IF (LSPUMP) SPUMP(ISPZ,MSURF)=SPUMP(ISPZ,MSURF)+WEIGHT
        IF (LSPUMP) LMETSPW(ISPZ)    = .TRUE.
        IRET = 2
        RETURN 

      ELSE
C  OTHERWISE: RESTORE WEIGHT = WEIGHT/PR, TIME, AND CONTINUE ANOTHER TIME STEP
        WEIGHT=WEIGHT/PR
        ITMSTP=ITMSTP+1
        TIME=TIME0
        IRET = 1
        RETURN
      ENDIF
      IRET = 0
      RETURN
      END
