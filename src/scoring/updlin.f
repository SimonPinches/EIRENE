      SUBROUTINE EIRENE_UPDLIN

!  update tallies (currently: COPV) after completion of 
!  trajectory. Use linear algebraic expressions of default tallies
!  
!  score per history --> automatically variances per history are available
!                        distinct from aposteriori evaluation of linear combinations

!  current version:
!    1)   total particle source             (sni=papl+pmpl+pipl)  
!    2)   total parallel momentum source    (smo=mapl+mmpl+mipl) 
!    3)   total electr. energy source       (see=eael+emel+eiel)  
!    4)   total ion energy source           (sei=eapl+empl+eipl) 
!    5)   internal energy source            (sei_int=sei-u*smo+ek*sni)

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

      IMPLICIT NONE
      
      INTEGER :: ICP, ICP2, ICP3, 
     .           ICO, IR, IPL, NMTSP, IS, IRD
      INTEGER, SAVE :: IFIRST=0

      REAL(DP), ALLOCATABLE, SAVE :: UAH(:,:),EKIN(:,:)
      REAL(DP) :: wtt, zvoliw, seit, sni, smo, seii

      IF (IFIRST == 0) THEN
         ALLOCATE (UAH(NPLS,NRTAL))
         ALLOCATE (EKIN(NPLS,NRTAL))
         
         DO IPL = 1, NPLSI
           IS = MPLSV(IPL)
           DO IR = 1, NRAD
             IRD = NCLTAL(IR)
cdr  ir is fine grid for background medium, geometry, etc...
cdr  ird is coarse grid for scoring
             IF (IRD > 0) THEN
               UAH(IPL,IRD) = BVIN(IS,IR)
               EKIN(IPL,IRD)= cvrssp(IPL) * UAH(IPL,IR)**2       ! eV
             ENDIF
           END DO
         END DO
         IFIRST = 1
      END IF

      ICP = NPLSI     ! ...+1:  summed ipls part. source, a+m+i+ph
      ICP2 = 2*NPLSI  ! ...+1:  summed ipls parallel mom. source, a+m+i+ph
      ICP3 = 3*NPLSI


      NMTSP=NPHOTI+NATMI+NMOLI+NIONI+NPLSI+NADVI+NALVI+NCLVI

      IF (NCPVI < ICP3+4) THEN
         IF (IFIRST == 0) THEN
            WRITE (IUNOUT,*) 'UPDLIN: COUPLE TALLY COPV IS TOO SMALL '
            WRITE (IUNOUT,*) 'NCPVI NEEDS TO BE AT LEAST ',ICP3+4
            WRITE (IUNOUT,*) 'COPV IS NOT UPDATED IN UPDLIN '
            IFIRST = 1
         END IF
         RETURN
      END IF

!  particle source (sni), ipls (=ipl) resolved, copv(icp+1:icp+nplsi)
      DO IPL = 1,NPLSI
        IF (LMETSP(NSPAN(14)+IPL-1) .OR.
     .      LMETSP(NSPAN(20)+IPL-1) .OR.
     .      LMETSP(NSPAN(26)+IPL-1) ) THEN

           DO ICO = 1,NCLMT
CDR  the present trajectory has visited NCLMT (coarse) scoring cells
             IR = ICLMT(ICO)

             COPV(ICP+IPL,IR) = 0._DP
             IF (LPAPL) COPV(ICP+IPL,IR)=COPV(ICP+IPL,IR)+PAPL(IPL,IR)
             IF (LPIPL) COPV(ICP+IPL,IR)=COPV(ICP+IPL,IR)+PIPL(IPL,IR)
             IF (LPMPL) COPV(ICP+IPL,IR)=COPV(ICP+IPL,IR)+PMPL(IPL,IR)
             IF (LPAPL.OR.LPMPL.OR.LPIPL) LMETSP(NMTSP+ICP+IPL)=.TRUE.
           END DO
        END IF
      END DO

!  parallel momentum source (smo), ipls (=ipl) resolved, copv(icp2+1:icp2+nplsi)
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

!  electron energy source (see),  no species index here, copv(icp3+1)
 
      DO ICO = 1,NCLMT
        IR = ICLMT(ICO)

        COPV(ICP3+1,IR) = 0._DP
        IF (LEAEL) COPV(ICP3+1,IR)=COPV(ICP3+1,IR)+EAEL(IR)
        IF (LEIEL) COPV(ICP3+1,IR)=COPV(ICP3+1,IR)+EIEL(IR)
        IF (LEMEL) COPV(ICP3+1,IR)=COPV(ICP3+1,IR)+EMEL(IR)
        IF (LEAEL.OR.LEMEL.OR.LEIEL) LMETSP(NMTSP+ICP3+1)=.TRUE.

!  ion energy source (sei)
        COPV(ICP3+2,IR) = 0._DP
        IF (LEAPL) COPV(ICP3+2,IR)=COPV(ICP3+2,IR)+EAPL(IR)
        IF (LEIPL) COPV(ICP3+2,IR)=COPV(ICP3+2,IR)+EIPL(IR)
        IF (LEMPL) COPV(ICP3+2,IR)=COPV(ICP3+2,IR)+EMPL(IR)
        IF (LEAPL.OR.LEMPL.OR.LEIPL) LMETSP(NMTSP+ICP3+2)=.TRUE.

      END DO

!  copv(...+2,ir) above was:  total energy source
!  copv(...+3,ir) below is :  internal ion energy source
      
      DO ICO = 1,NCLMT
        IR = ICLMT(ICO)

        COPV(ICP3+3,IR) = COPV(ICP3+2,IR)   
! SEI_TOTAL, now correct to find SEI_INTERNAL.....
        
        DO IPL = 1,NPLSI
          COPV(ICP3+3,IR) = COPV(ICP3+3,IR)
     .         - UAH(IPL,IR) * COPV(ICP2+IPL,IR)*              ! UA*SMO
     .           cveli2/amua*2._DP * SIGN(1._DP,UAH(IPL,IR))
     .         + EKIN(IPL,IR) * COPV(ICP+IPL,IR)               ! EKIN*SNI
          LMETSP(NMTSP+ICP3+3)=.TRUE. 
        END DO

      END DO


      RETURN

      ENTRY EIRENE_RESET_UPDLIN

      IFIRST = 0
      IF (ALLOCATED(UAH)) DEALLOCATE (UAH)
      IF (ALLOCATED(EKIN)) DEALLOCATE (EKIN)

      RETURN

      END SUBROUTINE EIRENE_UPDLIN
