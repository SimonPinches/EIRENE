      MODULE EIRMOD_UPDLIN
cdr  2020: comments, explanation added

cdr  Score linear algebraic expressions of eirene volumetric tallies.
cdr  Update the score after each completed history.
cdr  Currently: scored on COPV. A more logical tally would be ALGV,
cdr             as this routine provides an "algebraic expression of tallies".

cdr             Distinct from post processing routine ALGEBRA  this present routine
cdr             it is limited to linear algebraic expressions.
cdr             But for such the default variances (standard deviations) options
cdr             are available, distinct from the general algebraic ALGV tallies.
cdr             To achieve this, when UPDLIN is called after trajectory no. IHIST,
cdr             COPV contains the cumulated (summed) scores from trajectories 1 to IHIST.
cdr             Subr. STATIS can then infer the score per single trajectory IHIST,
cdr             and provide the squared score per history needed for updating the variance
cdr
cdr
cdr  Called after each completed trajectory.
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      PRIVATE

      PUBLIC :: EIRENE_UPDLIN, EIRENE_RESET_UPDLIN,
     .          EIRENE_PREPARE_UPDLIN

      INTEGER, SAVE :: IFIRST=0

      REAL(DP), ALLOCATABLE, SAVE :: UAH(:,:), EKIN(:,:)

cym UAH/EKIN should be shared - large memory usage
cymtest !$OMP THREADPRIVATE(IFIRST,UAH,EKIN)

      CONTAINS

cdr  Nov. 2015

cdr  internal energy sources:  make EAPL, EMPL, EIPL, EPPL also IPLS species-dependent
cdr  check for storage (COPV) and return, if not enough storage.
cdr  UPDLIN should be made a default eirene option
cdr  for linear combination of tallies (with their stat. variances)

      SUBROUTINE EIRENE_PREPARE_UPDLIN

cpb  prepare array UAL and EKIN for use by EIRENE_UPDLIN
cpb  this needs to be done before the first thread enters EIRENE_UPDLIN
cpb  for scoring on COPV

      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CGEOM
      USE EIRMOD_CZT1

      IMPLICIT NONE

      INTEGER :: IR, IPL, IPLV, IRD

      IF (IFIRST == 0) THEN
         ALLOCATE (UAH(NPLS,NRTAL))
         ALLOCATE (EKIN(NPLS,NRTAL))
         UAH = 0.0_dp
         EKIN = 0.0_dp

         DO IPL = 1, NPLSI
           IPLV = MPLSV(IPL)
           DO IR = 1, NRAD
             IRD = NCLTAL(IR)
cdr  ir is fine grid for background medium, geometry, etc...
cdr  ird is coarse grid for scoring
             IF (IRD > 0) THEN
               UAH(IPL,IRD) = 0._DP
               IF (LBVIN) UAH(IPL,IRD) = BVIN(IPLV,IR)
               EKIN(IPL,IRD)= cvrssp(IPL) * UAH(IPL,IRD)**2       ! eV
             ENDIF
           END DO
         END DO
         IFIRST = 1
      END IF
      RETURN

      END SUBROUTINE EIRENE_PREPARE_UPDLIN


      SUBROUTINE EIRENE_UPDLIN

!  update tallies (currently on: COPV) after completion of
!  trajectory. Use linear algebraic expressions of default tallies
!
!  score per history --> automatically variances per history are available
!                        distinct from aposteriori evaluation of linear combinations

!  current version:
!    1) total particle source          (sni=papl+pmpl+pipl       , ICP+1  ,ICP2)
!    2) total parallel momentum source (smo=mapl+mmpl+mipl       , ICP2+1 ,ICP3)
!    3) total ion energy source        (sei_tot=eapl+empl+eipl   , ICP3+1 ,ICP4)
!    4) internal ion energy source     (sei_int=sei-u*smo+ek*sni , ICP4+1 ,ICP5)
!    5) total electr. energy source    (see=eael+emel+eiel       , ICP5+1)

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

      INTEGER :: ICP, ICP2, ICP3, ICP4, ICP5,
     .           ICO, IR, IPL, NMTSP

      ICP = NPLSI     ! ...+1: summed ipls part. source, a+m+i+ph
      ICP2 = 2*NPLSI  ! ...+1: summed ipls parallel mom.
                      !                          source, a+m+i+ph
      ICP3 = 3*NPLSI  ! ...+1: summed ipls ion energy
                      !                   source, total, a+m+i+ph
      ICP4 = 4*NPLSI  ! ...+1: summed ipls ion energy
                      !                source, internal, a+m+i+ph
      ICP5 = 5*NPLSI  ! ...+1: summed electr. energy
                      !                         sources, a+m+i+ph
      NMTSP=NPHOTI+NATMI+NMOLI+NIONI+NPLSI+NADVI+NALVI+NCLVI

      IF (NCPVI < ICP5+1) THEN
         IF (IFIRST == 0) THEN
            WRITE (IUNOUT,*) 'UPDLIN: COUPLE TALLY COPV IS TOO SMALL '
            WRITE (IUNOUT,*) 'NCPVI NEEDS TO BE AT LEAST ',ICP5+1
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

#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP CRITICAL
#endif
             COPV(ICP+IPL,IR) = 0._DP
             IF (LPAPL) COPV(ICP+IPL,IR)=COPV(ICP+IPL,IR)+PAPL(IPL,IR)
             IF (LPIPL) COPV(ICP+IPL,IR)=COPV(ICP+IPL,IR)+PIPL(IPL,IR)
             IF (LPMPL) COPV(ICP+IPL,IR)=COPV(ICP+IPL,IR)+PMPL(IPL,IR)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP END CRITICAL
#endif
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

#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP CRITICAL
#endif
             COPV(ICP2+IPL,IR) = 0._DP
             IF (LMAPL) COPV(ICP2+IPL,IR)=
     .                  COPV(ICP2+IPL,IR)+MAPL(IPL,IR)
             IF (LMIPL) COPV(ICP2+IPL,IR)=
     .                  COPV(ICP2+IPL,IR)+MIPL(IPL,IR)
             IF (LMMPL) COPV(ICP2+IPL,IR)=
     .                  COPV(ICP2+IPL,IR)+MMPL(IPL,IR)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP END CRITICAL
#endif
             IF (LMAPL.OR.LMMPL.OR.LMIPL)
     .           LMETSP(NMTSP+ICP2+IPL)=.TRUE.

           END DO
        END IF
      END DO

!  electron energy source (see),  no species index here, copv(icp5+1)

      DO ICO = 1,NCLMT
        IR = ICLMT(ICO)

#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP CRITICAL
#endif
        COPV(ICP5+1,IR) = 0._DP
        IF (LEAEL) COPV(ICP5+1,IR)=COPV(ICP5+1,IR)+EAEL(IR)
        IF (LEIEL) COPV(ICP5+1,IR)=COPV(ICP5+1,IR)+EIEL(IR)
        IF (LEMEL) COPV(ICP5+1,IR)=COPV(ICP5+1,IR)+EMEL(IR)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP END CRITICAL
#endif
        IF (LEAEL.OR.LEMEL.OR.LEIEL) LMETSP(NMTSP+ICP5+1)=.TRUE.
      END DO


!  total ion energy source (sei), ipls (=ipl) resolved, copv(icp3+1:icp3+nplsi)
      DO IPL = 1,NPLSI
        IF (LMETSP(NSPAN(38)+IPL-1) .OR.
     .      LMETSP(NSPAN(44)+IPL-1) .OR.
     .      LMETSP(NSPAN(50)+IPL-1) ) THEN

          DO ICO = 1,NCLMT
            IR = ICLMT(ICO)


#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP CRITICAL
#endif
            COPV(ICP3+IPL,IR) = 0._DP
            IF (LEAPL) COPV(ICP3+IPL,IR)=COPV(ICP3+IPL,IR)+EAPL(IPL,IR)
            IF (LEIPL) COPV(ICP3+IPL,IR)=COPV(ICP3+IPL,IR)+EIPL(IPL,IR)
            IF (LEMPL) COPV(ICP3+IPL,IR)=COPV(ICP3+IPL,IR)+EMPL(IPL,IR)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP END CRITICAL
#endif
            IF (LEAPL.OR.LEMPL.OR.LEIPL) LMETSP(NMTSP+ICP3+IPL)=.TRUE.

          END DO
        END IF
      END DO

! SEI_TOTAL, now correct to find SEI_INTERNAL.....
!  internal ion energy source, ipls (=ipl) resolved, copv(icp4+1:icp4+nplsi)
      DO IPL = 1,NPLSI
        IF (LMETSP(NSPAN(38)+IPL-1) .OR.
     .      LMETSP(NSPAN(44)+IPL-1) .OR.
     .      LMETSP(NSPAN(50)+IPL-1) ) THEN

          DO ICO = 1,NCLMT
            IR = ICLMT(ICO)

#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP CRITICAL
#endif
            COPV(ICP4+IPL,IR) = 0._DP
            COPV(ICP4+IPL,IR) = COPV(ICP4+IPL,IR)
     .          - UAH(IPL,IR) * COPV(ICP2+IPL,IR)*      ! UA*SMO
     .           cveli2/amua*2._DP * SIGN(1._DP,UAH(IPL,IR))
     .          + EKIN(IPL,IR) * COPV(ICP+IPL,IR)       ! EKIN*SNI
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP END CRITICAL
#endif
            IF (LEAPL.OR.LEMPL.OR.LEIPL) LMETSP(NMTSP+ICP4+IPL)=.TRUE.
          END DO
        ENDIF

      END DO


      RETURN
      END SUBROUTINE EIRENE_UPDLIN

      SUBROUTINE EIRENE_RESET_UPDLIN

      IFIRST = 0
      IF (ALLOCATED(UAH)) DEALLOCATE (UAH)
      IF (ALLOCATED(EKIN)) DEALLOCATE (EKIN)

      RETURN

      END SUBROUTINE EIRENE_RESET_UPDLIN

      END MODULE EIRMOD_UPDLIN
