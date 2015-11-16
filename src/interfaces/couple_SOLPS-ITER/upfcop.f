      SUBROUTINE EIRENE_UPFCOP

!  update sources portions on couple tally COPV after completion of 
!  trajektory

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CESTIM
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT
      USE EIRMOD_CSDVI
      USE EIRMOD_CGEOM
      USE EIRMOD_CZT1
      USE EIRMOD_CCONA
      USE EIRMOD_COMSOU
      USE EIRMOD_COUTAU
      USE EIRMOD_CCOUPL

      IMPLICIT NONE
      
      INTEGER :: ICP, ICP2, ICP3, ICO, IR, IPL, NMTSP, IS, IRD
      INTEGER, SAVE :: IFIRST=0, ird1

      REAL(DP), ALLOCATABLE, SAVE :: UAH(:,:)
      REAL(DP) :: EKIN, wtt, zvoliw, seit, sni, smo, seii

      IF (IFIRST == 0) THEN
         ALLOCATE (UAH(NPLS,NRTAL))
         ird1 = 0
         DO IPL = 1, NPLSI
           IS = MPLSV(IPL)
           DO IR = 1, NRAD
             IRD = NCLTAL(IR)
             IF (IRD > 0) UAH(IPL,IRD) = BVIN(IS,IR)
             if ((ird1==0) .and. (ird==1)) ird1 = ir
           END DO
         END DO
         IFIRST = 1
      END IF

      ICP = NPLSI
      ICP2 = 2*NPLSI
      ICP3 = 3*NPLSI
      NMTSP=NPHOTI+NATMI+NMOLI+NIONI+NPLSI+NADVI+NALVI+NCLVI

      IF (NCPVI < ICP3+4) THEN
         IF (IFIRST == 0) THEN
            WRITE (IUNOUT,*) 'COUPLE TALLY COPV TOO SMALL '
            WRITE (IUNOUT,*) 'NCPVI NEEDS TO BE AT LEAST ',ICP3+4
            WRITE (IUNOUT,*) 'COPV IS NOT UPDATED '
            IFIRST = 1
         END IF
         RETURN
      END IF

!  particle source (sni)
      DO IPL = 1,NPLSI
        IF (LMETSP(NSPAN(14)+IPL-1) .OR.
     .      LMETSP(NSPAN(20)+IPL-1) .OR.
     .      LMETSP(NSPAN(26)+IPL-1) ) THEN

           DO ICO = 1,NCLMT
             IR = ICLMT(ICO)

             COPV(ICP+IPL,IR) = 0._DP
             IF (LPAPL) COPV(ICP+IPL,IR)=COPV(ICP+IPL,IR)+PAPL(IPL,IR)
             IF (LPIPL) COPV(ICP+IPL,IR)=COPV(ICP+IPL,IR)+PIPL(IPL,IR)
             IF (LPMPL) COPV(ICP+IPL,IR)=COPV(ICP+IPL,IR)+PMPL(IPL,IR)
             IF (LPAPL.OR.LPMPL.OR.LPIPL) LMETSP(NMTSP+ICP+IPL)=.TRUE.
           END DO
        END IF
      END DO

!  momentum source (smo)
      DO IPL = 1,NPLSI
        IF (LMETSP(NSPAN(97)+IPL-1) .OR.
     .      LMETSP(NSPAN(98)+IPL-1) .OR.
     .      LMETSP(NSPAN(99)+IPL-1) ) THEN

           DO ICO = 1,NCLMT
             IR = ICLMT(ICO)

             COPV(ICP2+IPL,IR) = 0._DP
             IF (LMAPL) COPV(ICP2+IPL,IR)=COPV(ICP2+IPL,IR)+MAPL(IPL,IR)
             IF (LMIPL) COPV(ICP2+IPL,IR)=COPV(ICP2+IPL,IR)+MIPL(IPL,IR)
             IF (LMMPL) COPV(ICP2+IPL,IR)=COPV(ICP2+IPL,IR)+MMPL(IPL,IR)
             IF (LMAPL.OR.LMMPL.OR.LMIPL) LMETSP(NMTSP+ICP2+IPL)=.TRUE.

           END DO
        END IF
      END DO

      DO ICO = 1,NCLMT
        IR = ICLMT(ICO)

!  electron energy source (see)
        COPV(ICP3+1,IR) = 0._DP
        IF (LEAEL) COPV(ICP3+1,IR)=COPV(ICP3+1,IR)+EAEL(IR)
        IF (LEIEL) COPV(ICP3+1,IR)=COPV(ICP3+1,IR)+EIEL(IR)
        IF (LEMEL) COPV(ICP3+1,IR)=COPV(ICP3+1,IR)+EMEL(IR)
        IF (LEAEL.OR.LEMEL.OR.LEIEL) LMETSP(NMTSP+ICP3+1)=.TRUE.
      END DO

      DO IPL = 1,NPLSI
        IF ((LMETSP(NSPAN(38)+IPL-1) .OR.
     .      LMETSP(NSPAN(44)+IPL-1) .OR.
     .      LMETSP(NSPAN(50)+IPL-1) ) .AND. (IFLB(IPL) > 0)) THEN

          DO ICO = 1,NCLMT
            IR = ICLMT(ICO)

!  ion energy source (sei)
            COPV(ICP3+2,IR) = 0._DP
            IF (LEAPL) COPV(ICP3+2,IR)=COPV(ICP3+2,IR)+EAPL(IPL,IR)
            IF (LEIPL) COPV(ICP3+2,IR)=COPV(ICP3+2,IR)+EIPL(IPL,IR)
            IF (LEMPL) COPV(ICP3+2,IR)=COPV(ICP3+2,IR)+EMPL(IPL,IR)
            IF (LEAPL.OR.LEMPL.OR.LEIPL) LMETSP(NMTSP+ICP3+2)=.TRUE.

          END DO
        END IF
      END DO


!  internal ion energy source
      
      DO ICO = 1,NCLMT
        IR = ICLMT(ICO)

        COPV(ICP3+3,IR) = COPV(ICP3+2,IR)   ! SEI_TOTAL
        
        DO IPL = 1,NPLSI

!pb          EKIN = 0.5_DP * RMASSP(IPL) * UAH(IPL,IR)**2
          EKIN = cvrssp(IPL) * UAH(IPL,IR)**2                  ! in eV
          COPV(ICP3+3,IR) = COPV(ICP3+3,IR)
     .         - UAH(IPL,IR) * COPV(ICP2+IPL,IR)*              ! UA*SMO
     .           cveli2/amua*2._DP * SIGN(1._DP,UAH(IPL,IR))
     .         + EKIN * COPV(ICP+IPL,IR)                       ! EKIN*SNI
          LMETSP(NMTSP+ICP3+3)=.TRUE. 
        END DO

      END DO

!  Detlevs tests

      if (.false.) then
! taken out again
      if (laddv) then

        DO ICO = 1,NCLMT
          IR = ICLMT(ICO)
          
          COPV(ICP3+4,IR) = ADDV(NADVI,IR) ! SEI_TOTAL
        
          DO IPL = 1,NPLSI

!pb          EKIN = 0.5_DP * RMASSP(IPL) * UAH(IPL,IR)**2
            EKIN = cvrssp(IPL) * UAH(IPL,IR)**2                ! in eV
            COPV(ICP3+4,IR) = COPV(ICP3+4,IR)
     .         - UAH(IPL,IR) * COPV(ICP2+IPL,IR)*              ! UA*SMO
     .           cveli2/amua*2._DP * SIGN(1._DP,UAH(IPL,IR))
     .         + EKIN * COPV(ICP+IPL,IR)                       ! EKIN*SNI
            LMETSP(NMTSP+ICP3+4)=.TRUE. 
          END DO

        END DO
      end if
      end if

      RETURN

      ENTRY EIRENE_RESET_UPFCOP

      IFIRST = 0
      IF (ALLOCATED(UAH)) DEALLOCATE (UAH)

      RETURN

      END SUBROUTINE EIRENE_UPFCOP
