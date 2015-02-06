c  score sputtered fluxes
c  modified in spring 2014: old version: resolved wrt. incident type
cdr:  Sept. 2014: input flag 'IND' added, FOR TOTAL TALLIES ALSO IN CASE SPUTTERED PARTICLES ARE NOT FOLLOWED

C  present version: resolved with respect to incidence species type (for nlscl option)
C  and also resolved wrt.  emitted species type and species 
c  INPUT:
c  incident type :                      itold  (parameter list)
c  ind=0:  species index of sputtered particle not known: update only total sputtered fluxes
c  ind=1:  species index of sputtered particle is known:  update total and species resolved fluxes
c  sputtered type:                      ityp   (common)
c  sputtered species:                   iphot,iatm,imol,iion,ipls  (common)
c  weight of sputtered particle:        wgh    (parameter list)

      SUBROUTINE EIRENE_UPDATE_SPTFLX (ITOLD, WGH,IND)
 
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

C  UPDATE TOTAL FLUXES, NOT RESOLVED WRT. EMITTED PARTICLE SPECIES
C  THIS IS NEEDED IN CASE SPUTTERED SPECIES IS NOT AN EIRENE TEST SPECIES IN THIS RUN
      
      IF (LSPTTOT) SPTTOT(MSURF) = SPTTOT(MSURF) + WGH
      IF (MSURFG.GT.0) THEN
         IF (LSPTTOT) SPTTOT(MSURFG) = SPTTOT(MSURFG) + WGH
      ENDIF
!pb  additional sputter tallies introduced
!pb      IF (IND.EQ.0) RETURN

C  FROM HERE ON:  INCIDENT TYPE AND SPUTTERED TYPE AND SPECIES RESOLVED FLUXES, 
C  ITPY AND ISPEZ ARE SET TO SPUTTERED PARTICLE SPECIES
 
      SELECT CASE (ITYP)
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
             
      SELECT CASE (ITOLD)

! INCIDENT PARTICLE IS PHOTON
      CASE (0)

         if (lsptphtot) sptphtot(msurf) = sptphtot(msurf) + wgh
         if (ind == 0) return

         SELECT CASE (ITYP)
!      OUTGOING PARTICLE IS PHOTON
         CASE(0)
            IF (LSPTPHPHT) THEN
               SPTPHPHT(IPHOT,MSURF)=SPTPHPHT(IPHOT,MSURF)+WGH
               IF (MSURFG.GT.0) 
     .              SPTPHPHT(IPHOT,MSURFG)=SPTPHPHT(IPHOT,MSURFG)+WGH
               LMETSPW(IPHOT) = .TRUE.
            END IF
!     OUTGOING PARTICLE IS ATOM
         CASE(1)
            IF (LSPTPHAT) THEN
               SPTPHAT(IATM,MSURF)=SPTPHAT(IATM,MSURF)+WGH
               IF (MSURFG.GT.0) 
     .              SPTPHAT(IATM,MSURFG)=SPTPHAT(IATM,MSURFG)+WGH
               LMETSPW(NSPH+IATM) = .TRUE.
            END IF
!     OUTGOING PARTICLE IS MOLECULE
         CASE(2)
            IF (LSPTPHML) THEN
               SPTPHML(IMOL,MSURF)=SPTPHML(IMOL,MSURF)+WGH
               IF (MSURFG.GT.0) 
     .              SPTPHML(IMOL,MSURFG)=SPTPHML(IMOL,MSURFG)+WGH
               LMETSPW(NSPA+IMOL) = .TRUE.
            END IF
!     OUTGOING PARTICLE IS TEST ION
         CASE(3)
            IF (LSPTPHIO) THEN
               SPTPHIO(IION,MSURF)=SPTPHIO(IION,MSURF)+WGH
               IF (MSURFG.GT.0) 
     .              SPTPHIO(IION,MSURFG)=SPTPHIO(IION,MSURFG)+WGH
               LMETSPW(NSPAM+IION) = .TRUE.
            END IF
!     OUTGOING PARTICLE IS BULK ION
         CASE(4)
            IF (LSPTPHPL) THEN
               SPTPHPL(IPLS,MSURF)=SPTPHPL(IPLS,MSURF)+WGH
               IF (MSURFG.GT.0) 
     .              SPTPHPL(IPLS,MSURFG)=SPTPHPL(IPLS,MSURFG)+WGH
               LMETSPW(NSPAMI+IPLS) = .TRUE.
            END IF
!     OUTGOING PARTICLE IS OF UNKNOWN TYPE
         CASE DEFAULT
            WRITE (IUNOUT,*) ' ERROR IN EIRENE_UPDATE_SPTFLX '
            WRITE (IUNOUT,*) 
     .           ' PARTICLE OF UNKNOWN TYPE SPUTERED '
            WRITE (IUNOUT,*) ' ITYP = ', ITYP
         END SELECT


! INCIDENT PARTICLE IS ATOM
      CASE (1)

         if (lsptatot) sptatot(msurf) = sptatot(msurf) + wgh
         if (ind == 0) return

         SELECT CASE (ITYP)
!      OUTGOING PARTICLE IS PHOTON
         CASE(0)
            IF (LSPTAPHT) THEN
               SPTAPHT(IPHOT,MSURF)=SPTAPHT(IPHOT,MSURF)+WGH
               IF (MSURFG.GT.0) 
     .              SPTAPHT(IPHOT,MSURFG)=SPTAPHT(IPHOT,MSURFG)+WGH
               LMETSPW(IPHOT) = .TRUE.
            END IF
!     OUTGOING PARTICLE IS ATOM
         CASE(1)
            IF (LSPTAAT) THEN
               SPTAAT(IATM,MSURF)=SPTAAT(IATM,MSURF)+WGH
               IF (MSURFG.GT.0) 
     .              SPTAAT(IATM,MSURFG)=SPTAAT(IATM,MSURFG)+WGH
               LMETSPW(NSPH+IATM) = .TRUE.
            END IF
!     OUTGOING PARTICLE IS MOLECULE
         CASE(2)
            IF (LSPTAML) THEN
               SPTAML(IMOL,MSURF)=SPTAML(IMOL,MSURF)+WGH
               IF (MSURFG.GT.0) 
     .              SPTAML(IMOL,MSURFG)=SPTAML(IMOL,MSURFG)+WGH
               LMETSPW(NSPA+IMOL) = .TRUE.
            END IF
!     OUTGOING PARTICLE IS TEST ION
         CASE(3)
            IF (LSPTAIO) THEN
               SPTAIO(IION,MSURF)=SPTAIO(IION,MSURF)+WGH
               IF (MSURFG.GT.0) 
     .              SPTAIO(IION,MSURFG)=SPTAIO(IION,MSURFG)+WGH
               LMETSPW(NSPAM+IION) = .TRUE.
            END IF
!     OUTGOING PARTICLE IS BULK ION
         CASE(4)
            IF (LSPTAPL) THEN
               SPTAPL(IPLS,MSURF)=SPTAPL(IPLS,MSURF)+WGH
               IF (MSURFG.GT.0) 
     .              SPTAPL(IPLS,MSURFG)=SPTAPL(IPLS,MSURFG)+WGH
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

         if (lsptmtot) sptmtot(msurf) = sptmtot(msurf) + wgh
         if (ind == 0) return

         SELECT CASE (ITYP)
! OUTGOING PARTICLE IS PHOTON
         CASE(0)
            IF (LSPTMPHT) THEN
               SPTMPHT(IPHOT,MSURF)=SPTMPHT(IPHOT,MSURF)+WGH
               IF (MSURFG.GT.0) 
     .              SPTMPHT(IPHOT,MSURFG)=SPTMPHT(IPHOT,MSURFG)+WGH
               LMETSPW(IPHOT) = .TRUE.
            END IF
! OUTGOING PARTICLE IS ATOM
         CASE(1)
            IF (LSPTMAT) THEN
               SPTMAT(IATM,MSURF)=SPTMAT(IATM,MSURF)+WGH
               IF (MSURFG.GT.0) 
     .              SPTMAT(IATM,MSURFG)=SPTMAT(IATM,MSURFG)+WGH
               LMETSPW(NSPH+IATM) = .TRUE.
            END IF
!     OUTGOING PARTICLE IS MOLECULE
         CASE(2)
            IF (LSPTMML) THEN
               SPTMML(IMOL,MSURF)=SPTMML(IMOL,MSURF)+WGH
               IF (MSURFG.GT.0) 
     .              SPTMML(IMOL,MSURFG)=SPTMML(IMOL,MSURFG)+WGH
               LMETSPW(NSPA+IMOL) = .TRUE.
            END IF
!     OUTGOING PARTICLE IS ATOMTEST ION
         CASE(3)
            IF (LSPTMIO) THEN
               SPTMIO(IION,MSURF)=SPTMIO(IION,MSURF)+WGH
               IF (MSURFG.GT.0) 
     .              SPTMIO(IION,MSURFG)=SPTMIO(IION,MSURFG)+WGH
               LMETSPW(NSPAM+IION) = .TRUE.
            END IF
!     OUTGOING PARTICLE IS BULK ION
         CASE(4)
            IF (LSPTMPL) THEN
               SPTMPL(IPLS,MSURF)=SPTMPL(IPLS,MSURF)+WGH
               IF (MSURFG.GT.0) 
     .              SPTMPL(IPLS,MSURFG)=SPTMPL(IPLS,MSURFG)+WGH
               LMETSPW(NSPAMI+IPLS) = .TRUE.
            END IF
!     OUTGOING PARTICLE IS OF UNKNOWN TYPE
         CASE DEFAULT
            WRITE (IUNOUT,*) ' ERROR IN EIRENE_UPDATE_SPTFLX '
            WRITE (IUNOUT,*) 
     .           ' PARTICLE OF UNKNOWN TYPE SPUTERED '
            WRITE (IUNOUT,*) ' ITYP = ', ITYP
         END SELECT


! INCIDENT PARTICLE IS TEST ION
      CASE (3)

         if (lsptitot) sptitot(msurf) = sptitot(msurf) + wgh
         if (ind == 0) return

         SELECT CASE (ITYP)
! OUTGOING PARTICLE IS PHOTON
         CASE(0)
            IF (LSPTIPHT) THEN
               SPTIPHT(IPHOT,MSURF)=SPTIPHT(IPHOT,MSURF)+WGH
               IF (MSURFG.GT.0) 
     .              SPTIPHT(IPHOT,MSURFG)=SPTIPHT(IPHOT,MSURFG)+WGH
               LMETSPW(IPHOT) = .TRUE.
            END IF
! OUTGOING PARTICLE IS ATOM
         CASE(1)
            IF (LSPTIAT) THEN
               SPTIAT(IATM,MSURF)=SPTIAT(IATM,MSURF)+WGH
               IF (MSURFG.GT.0) 
     .              SPTIAT(IATM,MSURFG)=SPTIAT(IATM,MSURFG)+WGH
               LMETSPW(NSPH+IATM) = .TRUE.
            END IF
!     OUTGOING PARTICLE IS MOLECULE
         CASE(2)
            IF (LSPTIML) THEN
               SPTIML(IMOL,MSURF)=SPTIML(IMOL,MSURF)+WGH
               IF (MSURFG.GT.0) 
     .              SPTIML(IMOL,MSURFG)=SPTIML(IMOL,MSURFG)+WGH
               LMETSPW(NSPA+IMOL) = .TRUE.
            END IF
!     OUTGOING PARTICLE IS TEST ION
         CASE(3)
            IF (LSPTIIO) THEN
               SPTIIO(IION,MSURF)=SPTIIO(IION,MSURF)+WGH
               IF (MSURFG.GT.0) 
     .              SPTIIO(IION,MSURFG)=SPTIIO(IION,MSURFG)+WGH
               LMETSPW(NSPAM+IION) = .TRUE.
            END IF
!     OUTGOING PARTICLE IS BULK ION
         CASE(4)
            IF (LSPTIPL) THEN
               SPTIPL(IPLS,MSURF)=SPTIPL(IPLS,MSURF)+WGH
               IF (MSURFG.GT.0) 
     .              SPTIPL(IPLS,MSURFG)=SPTIPL(IPLS,MSURFG)+WGH
               LMETSPW(NSPAMI+IPLS) = .TRUE.
            END IF
!     OUTGOING PARTICLE IS OF UNKNOWN TYPE
         CASE DEFAULT
            WRITE (IUNOUT,*) ' ERROR IN EIRENE_UPDATE_SPTFLX '
            WRITE (IUNOUT,*) 
     .           ' PARTICLE OF UNKNOWN TYPE SPUTERED '
            WRITE (IUNOUT,*) ' ITYP = ', ITYP
         END SELECT


! INCIDENT PARTICLE IS BULK ION
      CASE (4)

         if (lsptpltot) sptpltot(msurf) = sptpltot(msurf) + wgh
         if (ind == 0) return

         SELECT CASE (ITYP)
! OUTGOING PARTICLE IS PHOTON
         CASE(0)
            IF (LSPTPPHT) THEN
               SPTPPHT(IPHOT,MSURF)=SPTPPHT(IPHOT,MSURF)+WGH
               IF (MSURFG.GT.0) 
     .              SPTPPHT(IPHOT,MSURFG)=SPTPPHT(IPHOT,MSURFG)+WGH
               LMETSPW(IPHOT) = .TRUE.
            END IF
! OUTGOING PARTICLE IS ATOM
         CASE(1)
            IF (LSPTPAT) THEN
               SPTPAT(IATM,MSURF)=SPTPAT(IATM,MSURF)+WGH
               IF (MSURFG.GT.0) 
     .              SPTPAT(IATM,MSURFG)=SPTPAT(IATM,MSURFG)+WGH
               LMETSPW(NSPH+IATM) = .TRUE.
            END IF
!     OUTGOING PARTICLE IS MOLECULE
         CASE(2)
            IF (LSPTPML) THEN
               SPTPML(IMOL,MSURF)=SPTPML(IMOL,MSURF)+WGH
               IF (MSURFG.GT.0) 
     .              SPTPML(IMOL,MSURFG)=SPTPML(IMOL,MSURFG)+WGH
               LMETSPW(NSPA+IMOL) = .TRUE.
            END IF
!     OUTGOING PARTICLE IS TEST ION
         CASE(3)
            IF (LSPTPIO) THEN
               SPTPIO(IION,MSURF)=SPTPIO(IION,MSURF)+WGH
               IF (MSURFG.GT.0) 
     .              SPTPIO(IION,MSURFG)=SPTPIO(IION,MSURFG)+WGH
               LMETSPW(NSPAM+IION) = .TRUE.
            END IF
!     OUTGOING PARTICLE IS BULK ION
         CASE(4)
            IF (LSPTPPL) THEN
               SPTPPL(IPLS,MSURF)=SPTPPL(IPLS,MSURF)+WGH
               IF (MSURFG.GT.0) 
     .              SPTPPL(IPLS,MSURFG)=SPTPPL(IPLS,MSURFG)+WGH
               LMETSPW(NSPAMI+IPLS) = .TRUE.
            END IF
!     OUTGOING PARTICLE IS OF UNKNOWN TYPE
         CASE DEFAULT
            WRITE (IUNOUT,*) ' ERROR IN EIRENE_UPDATE_SPTFLX '
            WRITE (IUNOUT,*) 
     .           ' PARTICLE OF UNKNOWN TYPE SPUTERED '
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
