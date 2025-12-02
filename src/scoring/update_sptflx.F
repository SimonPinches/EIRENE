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




      SUBROUTINE EIRENE_UPDATE_SPTFLX (ITOLD,WGH,IND)

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

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ITOLD,IND
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

C  FROM HERE ON: INCIDENT TYPE (ITOLD) AND
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
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPHPHT(IPHOT,MSURFG)=SPTPHPHT(IPHOT,MSURFG)+WGH
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
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPHAT(IATM,MSURFG)=SPTPHAT(IATM,MSURFG)+WGH
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
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPHML(IMOL,MSURFG)=SPTPHML(IMOL,MSURFG)+WGH
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
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPHIO(IION,MSURFG)=SPTPHIO(IION,MSURFG)+WGH
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
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPHPL(IPLS,MSURFG)=SPTPHPL(IPLS,MSURFG)+WGH
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
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTAPHT(IPHOT,MSURFG)=SPTAPHT(IPHOT,MSURFG)+WGH
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
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTAAT(IATM,MSURFG)=SPTAAT(IATM,MSURFG)+WGH
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
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTAML(IMOL,MSURFG)=SPTAML(IMOL,MSURFG)+WGH
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
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTAIO(IION,MSURFG)=SPTAIO(IION,MSURFG)+WGH
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
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTAPL(IPLS,MSURFG)=SPTAPL(IPLS,MSURFG)+WGH
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
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTMPHT(IPHOT,MSURFG)=SPTMPHT(IPHOT,MSURFG)+WGH
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
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTMAT(IATM,MSURFG)=SPTMAT(IATM,MSURFG)+WGH
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
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTMML(IMOL,MSURFG)=SPTMML(IMOL,MSURFG)+WGH
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
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTMIO(IION,MSURFG)=SPTMIO(IION,MSURFG)+WGH
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
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTMPL(IPLS,MSURFG)=SPTMPL(IPLS,MSURFG)+WGH
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
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTIPHT(IPHOT,MSURFG)=SPTIPHT(IPHOT,MSURFG)+WGH
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
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTIAT(IATM,MSURFG)=SPTIAT(IATM,MSURFG)+WGH
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
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTIML(IMOL,MSURFG)=SPTIML(IMOL,MSURFG)+WGH
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
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTIIO(IION,MSURFG)=SPTIIO(IION,MSURFG)+WGH
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
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTIPL(IPLS,MSURFG)=SPTIPL(IPLS,MSURFG)+WGH
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
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPPHT(IPHOT,MSURFG)=SPTPPHT(IPHOT,MSURFG)+WGH
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
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPAT(IATM,MSURFG)=SPTPAT(IATM,MSURFG)+WGH
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
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPML(IMOL,MSURFG)=SPTPML(IMOL,MSURFG)+WGH
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
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPIO(IION,MSURFG)=SPTPIO(IION,MSURFG)+WGH
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
             IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
               SPTPPL(IPLS,MSURFG)=SPTPPL(IPLS,MSURFG)+WGH
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
