c  score sputtered fluxes
c  modified in spring 2014: old version: resolved wrt. incident type
cdr:  Sept. 2014: input flag 'IND' added,
cdr   Jan 18    : ind=0 --> ind=1, and ind=1 --> ind=2, to synchronize with
cdr               other surface scoring (e.g. update_surface)
c
cdr       ind=1: FOR TOTAL TALLIES (due to sputtering by "incident" particles)
cdr              ALSO IN CASE SPUTTERED PARTICLES ARE NOT FOLLOWED
cdr       ind=2: additionally: reemitted sputtered fluxes,
cdr              resolved wrt. species of emitted particle




      SUBROUTINE EIRENE_UPDATE_SPTFLX (ITOLD,IOLD,WGH,IND)

cdr  score sputtered fluxes on
cdr            sputtered flux tallies SP_A_OT(MSURF) (_A_: type of incident particle), IND=1, IND=2)
cdr        and sputtered (reemitted) particle species-resolved tallies: SP_A_B(MSURF)  IND=2

C  present version: resolved with respect to incidence species type (for nlscl option)
C  and also resolved wrt.  emitted species type and species
c  INPUT:
c  incident type :               itold (parameter list)
c  weight of sputtered particle: wgh   (parameter list)

c  ind=1:  type and species index of sputtered particle may not be not known:
c          Update only total sputtered fluxes resolved by outgoing flux particle type.
c          These total sputter tallies may contain sputtered (emitted) fluxes which are not
c          identified eirene test particles, i.e. these totals may be larger
c          than the sum over emitted species
c          of species-resolved sputtered tallies.
c  ind=2:  Type and species index of sputtered particle is known: update both: total and species-resolved fluxes
c          only in case IND=2:
c          sputtered (emitted) particle type:    ityp (common)
c          sputtered (emitted) particle species: iphot,iatm,imol,iion,ipls (common)

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CESTIM
      USE EIRMOD_COMPRT
      USE EIRMOD_COMUSR
      USE EIRMOD_CSPEZ
      USE EIRMOD_CSDVI
      USE EIRMOD_CLOGAU

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ITOLD,IOLD,IND
      REAL(DP), INTENT(IN) :: WGH

      IF (MSURF .LE. 0) RETURN

C  ALWAYS UPDATE TOTAL FLUXES, NOT RESOLVED WRT. EMITTED PARTICLE TYPE OR SPECIES.
c
C  THIS IS NEEDED IN CASE THE SPUTTERED (emitted) SPECIES IS NOT AN EIRENE TEST SPECIES IN THIS RUN
C  ONLY THE SPUTTERED FLUX IS SCORED

      IF (LSPTTOT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
        SPTTOT(MSURF) = SPTTOT(MSURF) + WGH
      ENDIF
      IF (MSURFG.GT.0) THEN
         IF (LSPTTOT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
           SPTTOT(MSURFG) = SPTTOT(MSURFG) + WGH
         ENDIF
      ENDIF
C  AT THIS PLACE: NEW ITYP AND NEW SPECIES INDEX NOT NECESSARILY KNOWN
      SELECT CASE (ITOLD)
      CASE (0)
        if (lsptphtot) then
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
          sptphtot(msurf) = sptphtot(msurf) + wgh
        endif
      CASE (1)
        if (lsptatot) then
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
           sptatot(msurf)  = sptatot(msurf)  + wgh
        endif
      CASE (2)
        if (lsptmtot) then
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
           sptmtot(msurf)  = sptmtot(msurf)  + wgh
        endif
      CASE (3)
        if (lsptitot) then
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
           sptitot(msurf)  = sptitot(msurf)  + wgh
        endif
      CASE (4)
         if (lsptpltot) then
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            sptpltot(msurf) = sptpltot(msurf) + wgh
         endif
      END SELECT

      IF (IND.EQ.1) RETURN

C  FROM HERE ON: INCIDENT TYPE (ITOLD), SPECIES (IOLD) AND
C  SPUTTERED TYPE (ITYP) AND SPECIES (IPHOT,IATM,....) RESOLVED FLUXES,
C  ITYP AND ISPEZ ARE SET TO SPUTTERED (emitted) PARTICLE SPECIES

      SELECT CASE (ITYP)  ! SAME FOR ALL INCIDENT TYPES ITOLD
      CASE (0)
        LOGPHOT(IPHOT,ISTRA)=.TRUE.
      CASE (1)
        LOGATM(IATM,ISTRA)=.TRUE.
      CASE (2)
        LOGMOL(IMOL,ISTRA)=.TRUE.
      CASE (3)
        LOGION(IION,ISTRA)=.TRUE.
      CASE (4)
        LOGPLS(IPLS,ISTRA)=.TRUE.
      END SELECT

      SELECT CASE (ITOLD)  ! incident type

! INCIDENT PARTICLE IS PHOTON
      CASE (0)

         SELECT CASE (ITYP)  ! emitted type
!      OUTGOING PARTICLE IS PHOTON
         CASE(0)
           IF (LSPTPHPHT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
             SPTPHPHT(IPHOT,MSURF)=SPTPHPHT(IPHOT,MSURF)+WGH
             IF (NLSPCSPT_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
               SPTPHPHT2(1:NPHOT,0:NPHOT) => SPTPHPHT(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPHPHT2(IPHOT,IOLD)=SPTPHPHT2(IPHOT,IOLD)+WGH
               LMETSPW2(1:NPHOT,0:NPHOT)=>LMETSPW(NTWS_PHI+1:NTWS_PHPH)
               LMETSPW2(IPHOT,0)=.TRUE.
               LMETSPW2(IPHOT,IOLD)=.TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPHPHT(IPHOT+IOLD*NPHOT,MSURF)=
     .          SPTPHPHT(IPHOT+IOLD*NPHOT,MSURF)+WGH
               LMETSPW(NTWS_PHI+IPHOT)=.TRUE.
               LMETSPW(NTWS_PHI+IPHOT+IOLD*NPHOT)=.TRUE.
#endif
             ENDIF
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPHPHT(IPHOT,MSURFG)=SPTPHPHT(IPHOT,MSURFG)+WGH
               IF (NLSPCSPT_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                 SPTPHPHT2(1:NPHOT,0:NPHOT) => SPTPHPHT(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTPHPHT2(IPHOT,IOLD)=SPTPHPHT2(IPHOT,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTPHPHT(IPHOT+IOLD*NPHOT,MSURFG)=
     .            SPTPHPHT(IPHOT+IOLD*NPHOT,MSURFG)+WGH
#endif
               ENDIF
             ENDIF
             LMETSPW(IPHOT) = .TRUE.
           END IF
!     OUTGOING PARTICLE IS ATOM
         CASE(1)
           IF (LSPTPHAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
             SPTPHAT(IATM,MSURF)=SPTPHAT(IATM,MSURF)+WGH
             IF (NLSPCSPT_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
               SPTPHAT2(1:NATM,0:NPHOT) => SPTPHAT(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPHAT2(IATM,IOLD)=SPTPHAT2(IATM,IOLD)+WGH
               LMETSPW2(1:NATM,0:NPHOT)=>LMETSPW(NTWS_IPH+1:NTWS_PHA)
               LMETSPW2(IATM,0)=.TRUE.
               LMETSPW2(IATM,IOLD)=.TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPHAT(IATM+IOLD*NATM,MSURF)=
     .          SPTPHAT(IATM+IOLD*NATM,MSURF)+WGH
               LMETSPW(NTWS_IPH+IATM)=.TRUE.
               LMETSPW(NTWS_IPH+IATM+IOLD*NATM)=.TRUE.
#endif
             ENDIF
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPHAT(IATM,MSURFG)=SPTPHAT(IATM,MSURFG)+WGH
               IF (NLSPCSPT_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                 SPTPHAT2(1:NATM,0:NPHOT) => SPTPHAT(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTPHAT2(IATM,IOLD)=SPTPHAT2(IATM,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTPHAT(IATM+IOLD*NATM,MSURFG)=
     .            SPTPHAT(IATM+IOLD*NATM,MSURFG)+WGH
#endif
               ENDIF
             ENDIF
             LMETSPW(NSPH+IATM) = .TRUE.
           END IF
!     OUTGOING PARTICLE IS MOLECULE
         CASE(2)
           IF (LSPTPHML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
             SPTPHML(IMOL,MSURF)=SPTPHML(IMOL,MSURF)+WGH
             IF (NLSPCSPT_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
               SPTPHML2(1:NMOL,0:NPHOT) => SPTPHML(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPHML2(IMOL,IOLD)=SPTPHML2(IMOL,IOLD)+WGH
               LMETSPW2(1:NMOL,0:NPHOT)=>LMETSPW(NTWS_PHA+1:NTWS_PHM)
               LMETSPW2(IMOL,0)=.TRUE.
               LMETSPW2(IMOL,IOLD)=.TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPHML(IMOL+IOLD*NMOL,MSURF)=
     .          SPTPHML(IMOL+IOLD*NMOL,MSURF)+WGH
               LMETSPW(NTWS_PHA+IMOL)=.TRUE.
               LMETSPW(NTWS_PHA+IMOL+IOLD*NMOL)=.TRUE.
#endif
             ENDIF
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPHML(IMOL,MSURFG)=SPTPHML(IMOL,MSURFG)+WGH
               IF (NLSPCSPT_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                 SPTPHML2(1:NMOL,0:NPHOT) => SPTPHML(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTPHML2(IMOL,IOLD)=SPTPHML2(IMOL,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTPHML(IMOL+IOLD*NMOL,MSURFG)=
     .            SPTPHML(IMOL+IOLD*NMOL,MSURFG)+WGH
#endif
               ENDIF
             ENDIF
             LMETSPW(NSPA+IMOL) = .TRUE.
           END IF
!     OUTGOING PARTICLE IS TEST ION
         CASE(3)
           IF (LSPTPHIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
             SPTPHIO(IION,MSURF)=SPTPHIO(IION,MSURF)+WGH
             IF (NLSPCSPT_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
               SPTPHIO2(1:NION,0:NPHOT) => SPTPHIO(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPHIO2(IION,IOLD)=SPTPHIO2(IION,IOLD)+WGH
               LMETSPW2(1:NION,0:NPHOT)=>LMETSPW(NTWS_PHM+1:NTWS_PHI)
               LMETSPW2(IION,0)=.TRUE.
               LMETSPW2(IION,IOLD)=.TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPHIO(IION+IOLD*NION,MSURF)=
     .          SPTPHIO(IION+IOLD*NION,MSURF)+WGH
               LMETSPW(NTWS_PHM+IION)=.TRUE.
               LMETSPW(NTWS_PHM+IION+IOLD*NION)=.TRUE.
#endif
             ENDIF
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPHIO(IION,MSURFG)=SPTPHIO(IION,MSURFG)+WGH
               IF (NLSPCSPT_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                 SPTPHIO2(1:NION,0:NPHOT) => SPTPHIO(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTPHIO2(IION,IOLD)=SPTPHIO2(IION,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTPHIO(IION+IOLD*NION,MSURFG)=
     .            SPTPHIO(IION+IOLD*NION,MSURFG)+WGH
#endif
               ENDIF
             ENDIF
             LMETSPW(NSPAM+IION) = .TRUE.
           END IF
!     OUTGOING PARTICLE IS BULK ION
         CASE(4)
           IF (LSPTPHPL) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
             SPTPHPL(IPLS,MSURF)=SPTPHPL(IPLS,MSURF)+WGH
             IF (NLSPCSPT_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
               SPTPHPL2(1:NPLS,0:NPHOT) => SPTPHPL(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPHPL2(IPLS,IOLD)=SPTPHPL2(IPLS,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPHPL(IPLS+IOLD*NPLS,MSURF)=
     .          SPTPHPL(IPLS+IOLD*NPLS,MSURF)+WGH
#endif
             ENDIF
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPHPL(IPLS,MSURFG)=SPTPHPL(IPLS,MSURFG)+WGH
               IF (NLSPCSPT_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                 SPTPHPL2(1:NPLS,0:NPHOT) => SPTPHPL(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTPHPL2(IPLS,IOLD)=SPTPHPL2(IPLS,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTPHPL(IPLS+IOLD*NPLS,MSURFG)=
     .            SPTPHPL(IPLS+IOLD*NPLS,MSURFG)+WGH
#endif
               ENDIF
             ENDIF
             LMETSPW(NSPAMI+IPLS) = .TRUE.
           END IF
!     OUTGOING PARTICLE IS OF UNKNOWN TYPE
         CASE DEFAULT
           WRITE (IUNOUT,*) ' ERROR IN EIRENE_UPDATE_SPTFLX '
           WRITE (IUNOUT,*)
     .           ' PARTICLE OF UNKNOWN TYPE SPUTTERED '
           WRITE (IUNOUT,*) ' ITYP = ', ITYP
         END SELECT


! INCIDENT PARTICLE IS ATOM
      CASE (1)

         SELECT CASE (ITYP)
!      OUTGOING PARTICLE IS PHOTON
         CASE(0)
           IF (LSPTAPHT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
             SPTAPHT(IPHOT,MSURF)=SPTAPHT(IPHOT,MSURF)+WGH
             IF (NLSPCSPT_ATM) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
               SPTAPHT2(1:NPHOT,0:NATM) => SPTAPHT(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTAPHT2(IPHOT,IOLD)=SPTAPHT2(IPHOT,IOLD)+WGH
               LMETSPW2(1:NPHOT,0:NATM)=>LMETSPW(NTWS_AI+1:NTWS_APH)
               LMETSPW2(IPHOT,0)=.TRUE.
               LMETSPW2(IPHOT,IOLD)=.TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTAPHT(IPHOT+IOLD*NPHOT,MSURF)=
     .          SPTAPHT(IPHOT+IOLD*NPHOT,MSURF)+WGH
               LMETSPW(NTWS_AI+IPHOT)=.TRUE.
               LMETSPW(NTWS_AI+IPHOT+IOLD*NPHOT)=.TRUE.
#endif
             ENDIF
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTAPHT(IPHOT,MSURFG)=SPTAPHT(IPHOT,MSURFG)+WGH
               IF (NLSPCSPT_ATM) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                 SPTAPHT2(1:NPHOT,0:NATM) => SPTAPHT(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTAPHT2(IPHOT,IOLD)=SPTAPHT2(IPHOT,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTAPHT(IPHOT+IOLD*NPHOT,MSURFG)=
     .            SPTAPHT(IPHOT+IOLD*NPHOT,MSURFG)+WGH
#endif
               ENDIF
             ENDIF
             LMETSPW(IPHOT) = .TRUE.
           END IF
!     OUTGOING PARTICLE IS ATOM
         CASE(1)
           IF (LSPTAAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
             SPTAAT(IATM,MSURF)=SPTAAT(IATM,MSURF)+WGH
             IF (NLSPCSPT_ATM) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
               SPTAAT2(1:NATM,0:NATM) => SPTAAT(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTAAT2(IATM,IOLD)=SPTAAT2(IATM,IOLD)+WGH
               LMETSPW2(1:NATM,0:NATM)=>LMETSPW(NSPZTOTWS+1:NTWS_AA)
               LMETSPW2(IATM,0)=.TRUE.
               LMETSPW2(IATM,IOLD)=.TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTAAT(IATM+IOLD*NATM,MSURF)=
     .          SPTAAT(IATM+IOLD*NATM,MSURF)+WGH
               LMETSPW(NSPZTOTWS+IATM)=.TRUE.
               LMETSPW(NSPZTOTWS+IATM+IOLD*NATM)=.TRUE.
#endif
             ENDIF
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTAAT(IATM,MSURFG)=SPTAAT(IATM,MSURFG)+WGH
               IF (NLSPCSPT_ATM) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                 SPTAAT2(1:NATM,0:NATM) => SPTAAT(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTAAT2(IATM,IOLD)=SPTAAT2(IATM,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTAAT(IATM+IOLD*NATM,MSURFG)=
     .            SPTAAT(IATM+IOLD*NATM,MSURFG)+WGH
#endif
               ENDIF
             ENDIF
             LMETSPW(NSPH+IATM) = .TRUE.
           END IF
!     OUTGOING PARTICLE IS MOLECULE
         CASE(2)
           IF (LSPTAML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
             SPTAML(IMOL,MSURF)=SPTAML(IMOL,MSURF)+WGH
             IF (NLSPCSPT_ATM) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
               SPTAML2(1:NMOL,0:NATM) => SPTAML(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTAML2(IMOL,IOLD)=SPTAML2(IMOL,IOLD)+WGH
               LMETSPW2(1:NMOL,0:NATM)=>LMETSPW(NTWS_AA+1:NTWS_AM)
               LMETSPW2(IMOL,0)=.TRUE.
               LMETSPW2(IMOL,IOLD)=.TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTAML(IMOL+IOLD*NMOL,MSURF)=
     .          SPTAML(IMOL+IOLD*NMOL,MSURF)+WGH
               LMETSPW(NTWS_AA+IMOL)=.TRUE.
               LMETSPW(NTWS_AA+IMOL+IOLD*NMOL)=.TRUE.
#endif
             ENDIF
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTAML(IMOL,MSURFG)=SPTAML(IMOL,MSURFG)+WGH
               IF (NLSPCSPT_ATM) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                 SPTAML2(1:NMOL,0:NATM) => SPTAML(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTAML2(IMOL,IOLD)=SPTAML2(IMOL,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTAML(IMOL+IOLD*NMOL,MSURFG)=
     .            SPTAML(IMOL+IOLD*NMOL,MSURFG)+WGH
#endif
               ENDIF
             ENDIF
             LMETSPW(NSPA+IMOL) = .TRUE.
           END IF
!     OUTGOING PARTICLE IS TEST ION
         CASE(3)
           IF (LSPTAIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
             SPTAIO(IION,MSURF)=SPTAIO(IION,MSURF)+WGH
             IF (NLSPCSPT_ATM) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
               SPTAIO2(1:NION,0:NATM) => SPTAIO(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTAIO2(IION,IOLD)=SPTAIO2(IION,IOLD)+WGH
               LMETSPW2(1:NION,0:NATM)=>LMETSPW(NTWS_AM+1:NTWS_AI)
               LMETSPW2(IION,0)=.TRUE.
               LMETSPW2(IION,IOLD)=.TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTAIO(IION+IOLD*NION,MSURF)=
     .          SPTAIO(IION+IOLD*NION,MSURF)+WGH
               LMETSPW(NTWS_AM+IION)=.TRUE.
               LMETSPW(NTWS_AM+IION+IOLD*NION)=.TRUE.
#endif
             ENDIF
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTAIO(IION,MSURFG)=SPTAIO(IION,MSURFG)+WGH
               IF (NLSPCSPT_ATM) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                 SPTAIO2(1:NION,0:NATM) => SPTAIO(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTAIO2(IION,IOLD)=SPTAIO2(IION,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTAIO(IION+IOLD*NION,MSURFG)=
     .            SPTAIO(IION+IOLD*NION,MSURFG)+WGH
#endif
               ENDIF
             ENDIF
             LMETSPW(NSPAM+IION) = .TRUE.
           END IF
!     OUTGOING PARTICLE IS BULK ION
         CASE(4)
           IF (LSPTAPL) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
             SPTAPL(IPLS,MSURF)=SPTAPL(IPLS,MSURF)+WGH
             IF (NLSPCSPT_ATM) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
               SPTAPL2(1:NPLS,0:NATM) => SPTAPL(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTAPL2(IPLS,IOLD)=SPTAPL2(IPLS,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTAPL(IPLS+IOLD*NPLS,MSURF)=
     .          SPTAPL(IPLS+IOLD*NPLS,MSURF)+WGH
#endif
             ENDIF
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTAPL(IPLS,MSURFG)=SPTAPL(IPLS,MSURFG)+WGH
               IF (NLSPCSPT_ATM) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                 SPTAPL2(1:NPLS,0:NATM) => SPTAPL(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTAPL2(IPLS,IOLD)=SPTAPL2(IPLS,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTAPL(IPLS+IOLD*NPLS,MSURFG)=
     .            SPTAPL(IPLS+IOLD*NPLS,MSURFG)+WGH
#endif
               ENDIF
             ENDIF
             LMETSPW(NSPAMI+IPLS) = .TRUE.
           END IF
!     OUTGOING PARTICLE IS OF UNKNOWN TYPE
         CASE DEFAULT
           WRITE (IUNOUT,*) ' ERROR IN EIRENE_UPDATE_SPTFLX '
           WRITE (IUNOUT,*) ' PARTICLE OF UNKNOWN TYPE SPUTTERED '
           WRITE (IUNOUT,*) ' ITYP = ', ITYP
         END SELECT


! INCIDENT PARTICLE IS MOLECULE
      CASE (2)

         SELECT CASE (ITYP)
! OUTGOING PARTICLE IS PHOTON
         CASE(0)
           IF (LSPTMPHT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
             SPTMPHT(IPHOT,MSURF)=SPTMPHT(IPHOT,MSURF)+WGH
             IF (NLSPCSPT_MOL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
               SPTMPHT2(1:NPHOT,0:NMOL) => SPTMPHT(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTMPHT2(IPHOT,IOLD)=SPTMPHT2(IPHOT,IOLD)+WGH
               LMETSPW2(1:NPHOT,0:NMOL)=>LMETSPW(NTWS_MI+1:NTWS_MPH)
               LMETSPW2(IPHOT,0)=.TRUE.
               LMETSPW2(IPHOT,IOLD)=.TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTMPHT(IPHOT+IOLD*NPHOT,MSURF)=
     .          SPTMPHT(IPHOT+IOLD*NPHOT,MSURF)+WGH
               LMETSPW(NTWS_MI+IPHOT)=.TRUE.
               LMETSPW(NTWS_MI+IPHOT+IOLD*NPHOT)=.TRUE.
#endif
             ENDIF
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTMPHT(IPHOT,MSURFG)=SPTMPHT(IPHOT,MSURFG)+WGH
               IF (NLSPCSPT_MOL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                 SPTMPHT2(1:NPHOT,0:NMOL) => SPTMPHT(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTMPHT2(IPHOT,IOLD)=SPTMPHT2(IPHOT,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTMPHT(IPHOT+IOLD*NPHOT,MSURFG)=
     .            SPTMPHT(IPHOT+IOLD*NPHOT,MSURFG)+WGH
#endif
               ENDIF
             ENDIF
             LMETSPW(IPHOT) = .TRUE.
           END IF
! OUTGOING PARTICLE IS ATOM
         CASE(1)
           IF (LSPTMAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
             SPTMAT(IATM,MSURF)=SPTMAT(IATM,MSURF)+WGH
             IF (NLSPCSPT_MOL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
               SPTMAT2(1:NATM,0:NMOL) => SPTMAT(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTMAT2(IATM,IOLD)=SPTMAT2(IATM,IOLD)+WGH
               LMETSPW2(1:NATM,0:NMOL)=>LMETSPW(NTWS_APH+1:NTWS_MA)
               LMETSPW2(IATM,0)=.TRUE.
               LMETSPW2(IATM,IOLD)=.TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTMAT(IATM+IOLD*NATM,MSURF)=
     .          SPTMAT(IATM+IOLD*NATM,MSURF)+WGH
               LMETSPW(NTWS_APH+IATM)=.TRUE.
               LMETSPW(NTWS_APH+IATM+IOLD*NATM)=.TRUE.
#endif
             ENDIF
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTMAT(IATM,MSURFG)=SPTMAT(IATM,MSURFG)+WGH
               IF (NLSPCSPT_MOL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                 SPTMAT2(1:NATM,0:NMOL) => SPTMAT(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTMAT2(IATM,IOLD)=SPTMAT2(IATM,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTMAT(IATM+IOLD*NATM,MSURFG)=
     .            SPTMAT(IATM+IOLD*NATM,MSURFG)+WGH
#endif
               ENDIF
             ENDIF
             LMETSPW(NSPH+IATM) = .TRUE.
           END IF
!     OUTGOING PARTICLE IS MOLECULE
         CASE(2)
           IF (LSPTMML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
             SPTMML(IMOL,MSURF)=SPTMML(IMOL,MSURF)+WGH
             IF (NLSPCSPT_MOL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
               SPTMML2(1:NMOL,0:NMOL) => SPTMML(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTMML2(IMOL,IOLD)=SPTMML2(IMOL,IOLD)+WGH
               LMETSPW2(1:NMOL,0:NMOL)=>LMETSPW(NTWS_MA+1:NTWS_MM)
               LMETSPW2(IMOL,0)=.TRUE.
               LMETSPW2(IMOL,IOLD)=.TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTMML(IMOL+IOLD*NMOL,MSURF)=
     .          SPTMML(IMOL+IOLD*NMOL,MSURF)+WGH
               LMETSPW(NTWS_MA+IMOL)=.TRUE.
               LMETSPW(NTWS_MA+IMOL+IOLD*NMOL)=.TRUE.
#endif
             ENDIF
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTMML(IMOL,MSURFG)=SPTMML(IMOL,MSURFG)+WGH
               IF (NLSPCSPT_MOL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                 SPTMML2(1:NMOL,0:NMOL) => SPTMML(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTMML2(IMOL,IOLD)=SPTMML2(IMOL,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTMML(IMOL+IOLD*NMOL,MSURFG)=
     .            SPTMML(IMOL+IOLD*NMOL,MSURFG)+WGH
#endif
               ENDIF
             ENDIF
             LMETSPW(NSPA+IMOL) = .TRUE.
           END IF
!     OUTGOING PARTICLE IS TEST ION
         CASE(3)
           IF (LSPTMIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
             SPTMIO(IION,MSURF)=SPTMIO(IION,MSURF)+WGH
             IF (NLSPCSPT_MOL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
               SPTMIO2(1:NION,0:NMOL) => SPTMIO(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTMIO2(IION,IOLD)=SPTMIO2(IION,IOLD)+WGH
               LMETSPW2(1:NION,0:NMOL)=>LMETSPW(NTWS_MM+1:NTWS_MI)
               LMETSPW2(IION,0)=.TRUE.
               LMETSPW2(IION,IOLD)=.TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTMIO(IION+IOLD*NMOL,MSURF)=
     .          SPTMIO(IION+IOLD*NMOL,MSURF)+WGH
               LMETSPW(NTWS_MM+IION)=.TRUE.
               LMETSPW(NTWS_MM+IION+IOLD*NION)=.TRUE.
#endif
             ENDIF
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTMIO(IION,MSURFG)=SPTMIO(IION,MSURFG)+WGH
               IF (NLSPCSPT_MOL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                 SPTMIO2(1:NION,0:NMOL) => SPTMIO(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTMIO2(IION,IOLD)=SPTMIO2(IION,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTMIO(IION+IOLD*NION,MSURFG)=
     .            SPTMIO(IION+IOLD*NION,MSURFG)+WGH
#endif
               ENDIF
             ENDIF
             LMETSPW(NSPAM+IION) = .TRUE.
           END IF
!     OUTGOING PARTICLE IS BULK ION
         CASE(4)
           IF (LSPTMPL) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
             SPTMPL(IPLS,MSURF)=SPTMPL(IPLS,MSURF)+WGH
             IF (NLSPCSPT_MOL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
               SPTMPL2(1:NPLS,0:NMOL) => SPTMPL(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTMPL2(IPLS,IOLD)=SPTMPL2(IPLS,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTMPL(IPLS+IOLD*NPLS,MSURF)=
     .          SPTMPL(IPLS+IOLD*NPLS,MSURF)+WGH
#endif
             ENDIF
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTMPL(IPLS,MSURFG)=SPTMPL(IPLS,MSURFG)+WGH
               IF (NLSPCSPT_MOL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                 SPTMPL2(1:NPLS,0:NMOL) => SPTMPL(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTMPL2(IPLS,IOLD)=SPTMPL2(IPLS,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTMPL(IPLS+IOLD*NPLS,MSURFG)=
     .            SPTMPL(IPLS+IOLD*NPLS,MSURFG)+WGH
#endif
               ENDIF
             ENDIF
             LMETSPW(NSPAMI+IPLS) = .TRUE.
           END IF
!     OUTGOING PARTICLE IS OF UNKNOWN TYPE
         CASE DEFAULT
           WRITE (IUNOUT,*) ' ERROR IN EIRENE_UPDATE_SPTFLX '
           WRITE (IUNOUT,*)
     .           ' PARTICLE OF UNKNOWN TYPE SPUTTERED '
           WRITE (IUNOUT,*) ' ITYP = ', ITYP
         END SELECT


! INCIDENT PARTICLE IS TEST ION
      CASE (3)

         SELECT CASE (ITYP)
! OUTGOING PARTICLE IS PHOTON
         CASE(0)
           IF (LSPTIPHT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
             SPTIPHT(IPHOT,MSURF)=SPTIPHT(IPHOT,MSURF)+WGH
             IF (NLSPCSPT_ION) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
               SPTIPHT2(1:NPHOT,0:NION) => SPTIPHT(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTIPHT2(IPHOT,IOLD)=SPTIPHT2(IPHOT,IOLD)+WGH
               LMETSPW2(1:NPHOT,0:NION)=>LMETSPW(NTWS_II+1:NTWS_IPH)
               LMETSPW2(IPHOT,0)=.TRUE.
               LMETSPW2(IPHOT,IOLD)=.TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTIPHT(IPHOT+IOLD*NPHOT,MSURF)=
     .          SPTIPHT(IPHOT+IOLD*NPHOT,MSURF)+WGH
               LMETSPW(NTWS_II+IPHOT)=.TRUE.
               LMETSPW(NTWS_II+IPHOT+IOLD*NPHOT)=.TRUE.
#endif
             ENDIF
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTIPHT(IPHOT,MSURFG)=SPTIPHT(IPHOT,MSURFG)+WGH
               IF (NLSPCSPT_ION) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                 SPTIPHT2(1:NPHOT,0:NION) => SPTIPHT(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTIPHT2(IPHOT,IOLD)=SPTIPHT2(IPHOT,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTIPHT(IPHOT+IOLD*NPHOT,MSURFG)=
     .            SPTIPHT(IPHOT+IOLD*NPHOT,MSURFG)+WGH
#endif
               ENDIF
             ENDIF
             LMETSPW(IPHOT) = .TRUE.
           END IF
! OUTGOING PARTICLE IS ATOM
         CASE(1)
           IF (LSPTIAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
             SPTIAT(IATM,MSURF)=SPTIAT(IATM,MSURF)+WGH
             IF (NLSPCSPT_ION) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
               SPTIAT2(1:NATM,0:NION) => SPTIAT(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTIAT2(IATM,IOLD)=SPTIAT2(IATM,IOLD)+WGH
               LMETSPW2(1:NATM,0:NION)=>LMETSPW(NTWS_MPH+1:NTWS_IA)
               LMETSPW2(IATM,0)=.TRUE.
               LMETSPW2(IATM,IOLD)=.TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTIAT(IATM+IOLD*NATM,MSURF)=
     .          SPTIAT(IATM+IOLD*NATM,MSURF)+WGH
               LMETSPW(NTWS_MPH+IATM)=.TRUE.
               LMETSPW(NTWS_MPH+IATM+IOLD*NATM)=.TRUE.
#endif
             ENDIF
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTIAT(IATM,MSURFG)=SPTIAT(IATM,MSURFG)+WGH
               IF (NLSPCSPT_ION) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                 SPTIAT2(1:NATM,0:NION) => SPTIAT(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTIAT2(IATM,IOLD)=SPTIAT2(IATM,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTIAT(IATM+IOLD*NATM,MSURFG)=
     .            SPTIAT(IATM+IOLD*NATM,MSURFG)+WGH
#endif
               ENDIF
             ENDIF
             LMETSPW(NSPH+IATM) = .TRUE.
           END IF
!     OUTGOING PARTICLE IS MOLECULE
         CASE(2)
           IF (LSPTIML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
             SPTIML(IMOL,MSURF)=SPTIML(IMOL,MSURF)+WGH
             IF (NLSPCSPT_ION) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
               SPTIML2(1:NMOL,0:NION) => SPTIML(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTIML2(IMOL,IOLD)=SPTIML2(IMOL,IOLD)+WGH
               LMETSPW2(1:NMOL,0:NION)=>LMETSPW(NTWS_IA+1:NTWS_IM)
               LMETSPW2(IMOL,0)=.TRUE.
               LMETSPW2(IMOL,IOLD)=.TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTIML(IMOL+IOLD*NMOL,MSURF)=
     .          SPTIML(IMOL+IOLD*NMOL,MSURF)+WGH
               LMETSPW(NTWS_IA+IMOL)=.TRUE.
               LMETSPW(NTWS_IA+IMOL+IOLD*NMOL)=.TRUE.
#endif
             ENDIF
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTIML(IMOL,MSURFG)=SPTIML(IMOL,MSURFG)+WGH
               IF (NLSPCSPT_ION) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                 SPTIML2(1:NMOL,0:NION) => SPTIML(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTIML2(IMOL,IOLD)=SPTIML2(IMOL,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTIML(IMOL+IOLD*NMOL,MSURFG)=
     .            SPTIML(IMOL+IOLD*NMOL,MSURFG)+WGH
#endif
               ENDIF
             ENDIF
             LMETSPW(NSPA+IMOL) = .TRUE.
           END IF
!     OUTGOING PARTICLE IS TEST ION
         CASE(3)
           IF (LSPTIIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
             SPTIIO(IION,MSURF)=SPTIIO(IION,MSURF)+WGH
             IF (NLSPCSPT_ION) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
               SPTIIO2(1:NION,0:NION) => SPTIIO(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTIIO2(IION,IOLD)=SPTIIO2(IION,IOLD)+WGH
               LMETSPW2(1:NION,0:NMOL)=>LMETSPW(NTWS_IM+1:NTWS_II)
               LMETSPW2(IION,0)=.TRUE.
               LMETSPW2(IION,IOLD)=.TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTIIO(IION+IOLD*NMOL,MSURF)=
     .          SPTIIO(IION+IOLD*NMOL,MSURF)+WGH
               LMETSPW(NTWS_MM+IION)=.TRUE.
               LMETSPW(NTWS_MM+IION+IOLD*NION)=.TRUE.
#endif
             ENDIF
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTIIO(IION,MSURFG)=SPTIIO(IION,MSURFG)+WGH
               IF (NLSPCSPT_ION) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                 SPTIIO2(1:NION,0:NION) => SPTIIO(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTIIO2(IION,IOLD)=SPTIIO2(IION,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTIIO(IION+IOLD*NION,MSURFG)=
     .            SPTIIO(IION+IOLD*NION,MSURFG)+WGH
#endif
               ENDIF
             ENDIF
             LMETSPW(NSPAM+IION) = .TRUE.
           END IF
!     OUTGOING PARTICLE IS BULK ION
         CASE(4)
           IF (LSPTIPL) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
             SPTIPL(IPLS,MSURF)=SPTIPL(IPLS,MSURF)+WGH
             IF (NLSPCSPT_ION) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
               SPTIPL2(1:NPLS,0:NION) => SPTIPL(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTIPL2(IPLS,IOLD)=SPTIPL2(IPLS,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTIPL(IPLS+IOLD*NPLS,MSURF)=
     .          SPTIPL(IPLS+IOLD*NPLS,MSURF)+WGH
#endif
             ENDIF
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTIPL(IPLS,MSURFG)=SPTIPL(IPLS,MSURFG)+WGH
               IF (NLSPCSPT_ION) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                 SPTIPL2(1:NPLS,0:NION) => SPTIPL(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTIPL2(IPLS,IOLD)=SPTIPL2(IPLS,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTIPL(IPLS+IOLD*NPLS,MSURFG)=
     .            SPTIPL(IPLS+IOLD*NPLS,MSURFG)+WGH
#endif
               ENDIF
             ENDIF
             LMETSPW(NSPAMI+IPLS) = .TRUE.
           END IF
!     OUTGOING PARTICLE IS OF UNKNOWN TYPE
         CASE DEFAULT
           WRITE (IUNOUT,*) ' ERROR IN EIRENE_UPDATE_SPTFLX '
           WRITE (IUNOUT,*)
     .           ' PARTICLE OF UNKNOWN TYPE SPUTTERED '
           WRITE (IUNOUT,*) ' ITYP = ', ITYP
         END SELECT


! INCIDENT PARTICLE IS BULK ION
      CASE (4)

         SELECT CASE (ITYP)
! OUTGOING PARTICLE IS PHOTON
         CASE(0)
           IF (LSPTPPHT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
             SPTPPHT(IPHOT,MSURF)=SPTPPHT(IPHOT,MSURF)+WGH
             IF (NLSPCSPT_PLS) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
               SPTPPHT2(1:NPHOT,0:NPLS) => SPTPPHT(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPPHT2(IPHOT,IOLD)=SPTPPHT2(IPHOT,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPPHT(IPHOT+IOLD*NPHOT,MSURF)=
     .          SPTPPHT(IPHOT+IOLD*NPHOT,MSURF)+WGH
#endif
             ENDIF
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPPHT(IPHOT,MSURFG)=SPTPPHT(IPHOT,MSURFG)+WGH
               IF (NLSPCSPT_PLS) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                 SPTPPHT2(1:NPHOT,0:NPLS) => SPTPPHT(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTPPHT2(IPHOT,IOLD)=SPTPPHT2(IPHOT,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTPPHT(IPHOT+IOLD*NPHOT,MSURFG)=
     .            SPTPPHT(IPHOT+IOLD*NPHOT,MSURFG)+WGH
#endif
               ENDIF
             ENDIF
             LMETSPW(IPHOT) = .TRUE.
           END IF
! OUTGOING PARTICLE IS ATOM
         CASE(1)
           IF (LSPTPAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
             SPTPAT(IATM,MSURF)=SPTPAT(IATM,MSURF)+WGH
             IF (NLSPCSPT_PLS) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
               SPTPAT2(1:NATM,0:NPLS) => SPTPAT(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPAT2(IATM,IOLD)=SPTPAT2(IATM,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPAT(IATM+IOLD*NATM,MSURF)=
     .          SPTPAT(IATM+IOLD*NATM,MSURF)+WGH
#endif
             ENDIF
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPAT(IATM,MSURFG)=SPTPAT(IATM,MSURFG)+WGH
               IF (NLSPCSPT_PLS) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                 SPTPAT2(1:NATM,0:NPLS) => SPTPAT(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTPAT2(IATM,IOLD)=SPTPAT2(IATM,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTPAT(IATM+IOLD*NATM,MSURFG)=
     .            SPTPAT(IATM+IOLD*NATM,MSURFG)+WGH
#endif
               ENDIF
             ENDIF
             LMETSPW(NSPH+IATM) = .TRUE.
           END IF
!     OUTGOING PARTICLE IS MOLECULE
         CASE(2)
           IF (LSPTPML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
             SPTPML(IMOL,MSURF)=SPTPML(IMOL,MSURF)+WGH
             IF (NLSPCSPT_PLS) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
               SPTPML2(1:NMOL,0:NPLS) => SPTPML(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPML2(IMOL,IOLD)=SPTPML2(IMOL,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPML(IMOL+IOLD*NMOL,MSURF)=
     .          SPTPML(IMOL+IOLD*NMOL,MSURF)+WGH
#endif
             ENDIF
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPML(IMOL,MSURFG)=SPTPML(IMOL,MSURFG)+WGH
               IF (NLSPCSPT_PLS) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                 SPTPML2(1:NMOL,0:NPLS) => SPTPML(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTPML2(IMOL,IOLD)=SPTPML2(IMOL,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTPML(IMOL+IOLD*NMOL,MSURFG)=
     .            SPTPML(IMOL+IOLD*NMOL,MSURFG)+WGH
#endif
               ENDIF
             ENDIF
             LMETSPW(NSPA+IMOL) = .TRUE.
           END IF
!     OUTGOING PARTICLE IS TEST ION
         CASE(3)
           IF (LSPTPIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
             SPTPIO(IION,MSURF)=SPTPIO(IION,MSURF)+WGH
             IF (NLSPCSPT_PLS) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
               SPTPIO2(1:NION,0:NPLS) => SPTPIO(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPIO2(IION,IOLD)=SPTPIO2(IION,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPIO(IION+IOLD*NMOL,MSURF)=
     .          SPTPIO(IION+IOLD*NMOL,MSURF)+WGH
#endif
             ENDIF
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPIO(IION,MSURFG)=SPTPIO(IION,MSURFG)+WGH
               IF (NLSPCSPT_PLS) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                 SPTPIO2(1:NION,0:NPLS) => SPTPIO(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTPIO2(IION,IOLD)=SPTPIO2(IION,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTPIO(IION+IOLD*NION,MSURFG)=
     .            SPTPIO(IION+IOLD*NION,MSURFG)+WGH
#endif
               ENDIF
             ENDIF
             LMETSPW(NSPAM+IION) = .TRUE.
           END IF
!     OUTGOING PARTICLE IS BULK ION
         CASE(4)
           IF (LSPTPPL) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
             SPTPPL(IPLS,MSURF)=SPTPPL(IPLS,MSURF)+WGH
             IF (NLSPCSPT_PLS) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
               SPTPPL2(1:NPLS,0:NPLS) => SPTPPL(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPPL2(IPLS,IOLD)=SPTPPL2(IPLS,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPPL(IPLS+IOLD*NPLS,MSURF)=
     .          SPTPPL(IPLS+IOLD*NPLS,MSURF)+WGH
#endif
             ENDIF
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPPL(IPLS,MSURFG)=SPTPPL(IPLS,MSURFG)+WGH
               IF (NLSPCSPT_PLS) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                 SPTPPL2(1:NPLS,0:NPLS) => SPTPPL(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTPPL2(IPLS,IOLD)=SPTPPL2(IPLS,IOLD)+WGH
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                 SPTPPL(IPLS+IOLD*NPLS,MSURFG)=
     .            SPTPPL(IPLS+IOLD*NPLS,MSURFG)+WGH
#endif
               ENDIF
             ENDIF
             LMETSPW(NSPAMI+IPLS) = .TRUE.
           END IF
!     OUTGOING PARTICLE IS OF UNKNOWN TYPE
         CASE DEFAULT
           WRITE (IUNOUT,*) ' ERROR IN EIRENE_UPDATE_SPTFLX '
           WRITE (IUNOUT,*)
     .           ' PARTICLE OF UNKNOWN TYPE SPUTTERED '
           WRITE (IUNOUT,*) ' ITYP = ', ITYP
         END SELECT


! INCIDENT PARTICLE IS OF UNKNOWN TYP
      CASE DEFAULT
        WRITE (IUNOUT,*) ' ERROR IN EIRENE_UPDATE_SPTFLX '
        WRITE (IUNOUT,*) ' UNKNOWN TYPE OF INCIDENT PARTICLE '
        WRITE (IUNOUT,*) ' ITOLD = ', ITOLD
      END SELECT

      RETURN
      END SUBROUTINE EIRENE_UPDATE_SPTFLX
