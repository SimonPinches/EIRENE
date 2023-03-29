C
c  written by P. Boerner, for FZJ proprietary IDL plotting tool.
c  not intended for 3rd party use.
c  last modified: jan 2017
cdr correction: 3 digits rather than 2 digits for I0 in outtal file name

      SUBROUTINE EIRENE_OUTIDLTAL

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_CADGEO
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CGRID
      USE EIRMOD_CSPEZ
      USE EIRMOD_CTRCEI
      USE EIRMOD_CGEOM
      USE EIRMOD_CSDVI
      USE EIRMOD_COMPRT
      USE EIRMOD_COMNNL
      USE EIRMOD_COMSOU
      USE EIRMOD_CTEXT
      USE EIRMOD_COUTAU
      USE EIRMOD_CSPEI

      IMPLICIT NONE
C
      REAL(DP), ALLOCATABLE :: VECTOR(:,:),TALAV(:),TALTOT(:)
      REAL(DP) :: OUTAUI
      INTEGER :: NFTI, NFTE, K, ITAL, I, ISTR, MXSPZ, IOUT, IN, KK
      LOGICAL :: LFIRST
C
      CHARACTER(50) :: FNAME, FORMA, FORME, FORME2
C
!      IF (NSBOX_TAL /= NSBOX) THEN
!         WRITE (IUNOUT,*) ' ERROR IN OUTIDLTAL '
!         WRITE (IUNOUT,*) ' NSBOX_TAL /= NSBOX '
!         WRITE (IUNOUT,*) ' THIS CASE IS NOT YET FORESEEN '
!         WRITE (IUNOUT,*) ' NO DATA WRITTEN '
!         RETURN
!      END IF

      MXSPZ = MAXVAL(NFSTVI(1:NTALV))
      MXSPZ = MAX(MXSPZ, NATM, NMOL, NION, NPHOT, NADV, NALV)
      ALLOCATE (VECTOR(NRAD,MXSPZ))
      ALLOCATE (TALTOT(MXSPZ))
      ALLOCATE (TALAV(MXSPZ))

      IOUT = 3 + IFOFF

      FORMA=REPEAT(' ',50)
      FORMA='(6X,   A25)'
      WRITE (FORMA(5:7),'(I3)') MXSPZ

      FORME=REPEAT(' ',50)
      FORME='(I10,   ES25.7)'
      WRITE (FORME(6:8),'(I3)') MXSPZ

      FORME2=REPEAT(' ',50)
      FORME2='(6X,   ES25.7)'
      WRITE (FORME2(5:7),'(I3)') MXSPZ

      LFIRST=.TRUE.
C
C ISTRA IS THE STRATUM NUMBER. ISTRA=0 STANDS FOR: SUM OVER STRATA

      DO ISTR = 0, NSTRAI
        ISTRA=ISTR

        IF ((ISTRA == 0) .AND. (NSMSTRA /= 1)) CYCLE
C
        IF (XMCP(ISTRA).LT.1.) CYCLE
C
C
        IF (ISTRA.EQ.IESTR) THEN
C  NOTHING TO BE DONE
        ELSEIF (NFILEN.EQ.1.OR.NFILEN.EQ.2) THEN
          IESTR=ISTRA
          CALL EIRENE_RSTRT(ISTRA,NSTRAI,NESTM1,NESTM2,NADSPC,
     .             ESTIMV,ESTIMS,ESTIML,
     .             NSDVI1,SDVI1,NSDVI2,SDVI2,
     .             NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .             NSIGI_SPC,TRCFLE)
          IF (NLSYMP(ISTRA).OR.NLSYMT(ISTRA)) THEN
            CALL EIRENE_SYMET(ESTIMV,NVOLTL,NRTAL,NR1TAL,NP2TAL,NT3TAL,
     .               NLSYMP(ISTRA),NLSYMT(ISTRA))
          ENDIF
        ELSEIF ((NFILEN.EQ.6.OR.NFILEN.EQ.7).AND.ISTRA.EQ.0) THEN
          IESTR=ISTRA
          CALL EIRENE_RSTRT(ISTRA,NSTRAI,NESTM1,NESTM2,NADSPC,
     .             ESTIMV,ESTIMS,ESTIML,
     .             NSDVI1,SDVI1,NSDVI2,SDVI2,
     .             NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .             NSIGI_SPC,TRCFLE)
          IF (NLSYMP(ISTRA).OR.NLSYMT(ISTRA)) THEN
            CALL EIRENE_SYMET(ESTIMV,NVOLTL,NRTAL,NR1TAL,NP2TAL,NT3TAL,
     .               NLSYMP(ISTRA),NLSYMT(ISTRA))
          ENDIF
        ELSE
          WRITE (iunout,*) 'ERROR IN OUTEIR: DATA FOR STRATUM ISTRA= ',
     .                     ISTRA
          WRITE (iunout,*) 'ARE NOT AVAILABLE. PRINTOUT ABANDONED'
          CYCLE
        ENDIF
C
C
C  PRINT VOLUME-AVERAGED TALLIES

        DO 100 ITAL = 1, NTALV

          IF (.NOT.LIVTALV(ITAL)) CYCLE
C
          NFTI=1
          NFTE=NFSTVI(ITAL)

          IF (ITAL == NTALA) THEN
            NFTE = SUM(VERIFY(TXTSPC(1:NADV,NTALA),' '))
          END IF

          KK = 0
          DO 119 K=NFTI,NFTE
            IF (NEXTVI(ITAL) > 0) THEN
              KK = KK + 1
              IF ((K > NEXTVI(ITAL) .AND. 
     .             (MOD(K,NEXTVI(ITAL)) == 1))) KK = KK + 1
            ELSE
              KK = K
            END IF
            CALL EIRENE_FETCH_OUTAU (OUTAUI,ITAL,KK,ISTRA,IUNOUT)
C
            IF (NSBOX_TAL /= NSBOX) THEN
              DO I=1,NSBOX
                IN = NCLTAL(I)
                IF (IN > 0) THEN
                  VECTOR(I,K)=ESTIMV(NADDV(ITAL)+K,IN)
                ELSE
                  VECTOR(I,K)=0._DP
                END IF
              END DO
            ELSE
              DO 110 I=1,NSBOX_TAL
                VECTOR(I,K)=ESTIMV(NADDV(ITAL)+K,I)
  110         CONTINUE
          END IF

            TALTOT(K)=OUTAUI
            TALAV(K)=TALTOT(K)/VOLTOT
  119     CONTINUE
C
          FNAME = 'outtal_   '
          WRITE (FNAME(8:10),'(i0)') ITAL

          IF (LFIRST) THEN
            OPEN (UNIT=IOUT,FILE=FNAME,FORM='FORMATTED',
     .            ACCESS='SEQUENTIAL')
          ELSE
            OPEN (UNIT=IOUT,FILE=FNAME,FORM='FORMATTED',
     .            ACCESS='SEQUENTIAL',POSITION='APPEND')
          END IF

          WRITE (IOUT,'(A)')
     .      '+++++++++++++++++++++++++++++++++++++++++++++++++'
          WRITE (IOUT,'(A,I6)') 'ISTRA = ',ISTRA
          WRITE (IOUT,'(A)')
     .      '+++++++++++++++++++++++++++++++++++++++++++++++++'

          WRITE (IOUT,'(A)') TXTTAL(1,ITAL)
          WRITE (IOUT,'(A,I10)') 'NCELLS:   ',NSBOX
!pb          WRITE (IOUT,'(A,I10)') 'NSPECIES: ',NFSTVI(ITAL)
          WRITE (IOUT,'(A,I10)') 'NSPECIES: ',NFTE
          WRITE (IOUT,'(A)') 'SPECIES'
          WRITE (IOUT,FORMA) (TRIM(TXTSPC(K,ITAL)), K=NFTI, NFTE)
          WRITE (IOUT,'(A)') 'UNITS'
          WRITE (IOUT,FORMA) (TRIM(TXTUNT(K,ITAL)), K=NFTI, NFTE)

          WRITE (IOUT,'(A)')
     .        'TOTAL ("UNITS*CM**3), AND MEAN VALUE ("UNITS") '
          WRITE (IOUT,'(A)') 'TOTAL'
          WRITE (IOUT,FORME2) (TALTOT(K), K=NFTI, NFTE)
          WRITE (IOUT,'(A)') 'MEAN'
          WRITE (IOUT,FORME2) (TALAV(K), K=NFTI, NFTE)

          WRITE (IOUT,'(A)')
     .      '================================================='
          DO I=1, NSBOX
            IF (ANY(ABS(VECTOR(I,NFTI:NFTE)) > EPS30))
     .        WRITE (IOUT,FORME) I,(VECTOR(I,K), K=NFTI, NFTE)
          END DO
          WRITE (IOUT,'(A)')
     .      '================================================='

          CLOSE (UNIT=IOUT)
C
  100   CONTINUE

        LFIRST = .FALSE.

      END DO

      DEALLOCATE (VECTOR)
      DEALLOCATE (TALTOT)
      DEALLOCATE (TALAV)

      RETURN
      END SUBROUTINE EIRENE_OUTIDLTAL
