cdr    dec. 15:  added species index ipls, for volumetric energy source tallies for bulk ions
cdr              eapl,empl,eipl,ephpl
cdr   24.09.14:  scaling of new sputter tallies with fatm, fmol,fion,nphot: corrected
c  spring 2014:  new sputter tallies introduced: emitted species-resolved
C  15.02.05 :    double printout: fatm2,....taken out. use only getscl4, not getscl
C   6. 7.05 :    call ph_integrate for photon-background tallies taken out.
C                no more additional photon background tallies active
C  15.12.05 :    rescaling connected to spump surface tally


      SUBROUTINE EIRENE_SCALE_TALLIES (ISTRA)
C
C  RESCALE TRACKLENGTH-ESTIMATED VOLUME-AVERAGED TALLIES
C  WITH PARTICLE BALANCE CORRECTION FACTORS FATM,FMOL,FION,FPHOT
C  TO ENFORCE PERFECT GLOBAL PARTICLE BALANCE
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY : IUNOUT
      USE EIRMOD_CLOGAU
      USE EIRMOD_COUTAU
      USE EIRMOD_CGRID
      USE EIRMOD_CCONA
      USE EIRMOD_CESTIM
      USE EIRMOD_CSPEZ

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ISTRA
      REAL(DP), ALLOCATABLE :: FATM(:), FMOL(:), FION(:), FPHOT(:)
      REAL(DP) :: FADD, SUMP, SUMPRF
      INTEGER :: IATM, IMOL, IION, IPLS, IPHOT, IADV, ICLV, IBGV, ICPV,
     .           JATM, JMOL, JION, JPHOT, J, ISPC, 
     .           IAD, EIRENE_INDIRECT_ADDRESS
C
      ALLOCATE(FATM(0:NATMI))
      ALLOCATE(FMOL(0:NMOLI))
      ALLOCATE(FION(0:NIONI))
      ALLOCATE(FPHOT(0:NPHOTI))
      CALL EIRENE_GETSCL4 (ISTRA,FATM(0),FMOL(0),FION(0),FPHOT(0))
      IF (NLSPCSCL) THEN
        CALL EIRENE_GETSCLS (ISTRA,NATMI,NMOLI,NIONI,NPHOTI,
     .                       FATM,FMOL,FION,FPHOT)
      ELSE
        FATM(1:NATMI) = FATM(0)
        FMOL(1:NMOLI) = FMOL(0)
        FION(1:NIONI) = FION(0)
        FPHOT(1:NPHOTI) = FPHOT(0)
      END IF
C
      IF (.NOT.NLSCL) THEN

        CALL EIRENE_LEER(1)
        WRITE (iunout,*) 'NO RESCALING DONE (NLSCL=FALSE)'
        CALL EIRENE_LEER(2)

      ELSEIF (NLSCL) THEN

        IF (NLSPCSCL) THEN
          FASCL(:,ISTRA) = FATM
          FMSCL(:,ISTRA) = FMOL
          FISCL(:,ISTRA) = FION
          FPHSCL(:,ISTRA)= FPHOT
        ELSE
          FASCL(:,ISTRA) = FATM(0)
          FMSCL(:,ISTRA) = FMOL(0)
          FISCL(:,ISTRA) = FION(0)
          FPHSCL(:,ISTRA)= FPHOT(0)
        END IF
C
C  CARRY OUT SCALING OF VOLUME AND SURFACE TALLIES, RESP.
C
C  ATOM TALLIES
C
c  volumetric atomic tallies
        DO 2101 IATM=1,NATMI
          if (.not.logatm(iatm,istra)) cycle
          DO 111 J=1,NSBOX_TAL
            IF (LPDENA) PDENA(IATM,J)=PDENA(IATM,J)*FATM(IATM)
            IF (LEDENA) EDENA(IATM,J)=EDENA(IATM,J)*FATM(IATM)
            IF (LPAAT) THEN
              IF (NLSPCSCL_ATM) THEN
                PAAT(IATM,J)=0._DP
                DO JATM = 1,NATMI
                  IAD = EIRENE_INDIRECT_ADDRESS(IATM,JATM,NATM)
                  PAAT(IATM,J)=PAAT(IATM,J)+PAAT(IAD,J)*FATM(JATM)
                END DO
              ELSE  
                PAAT(IATM,J)=PAAT(IATM,J)*FATM(0)
              END IF
            END IF
            IF (LPMAT) THEN
              IF (NLSPCSCL_MOL) THEN
                PMAT(IATM,J)=0._DP
                DO JMOL = 1,NMOLI
                  IAD = EIRENE_INDIRECT_ADDRESS(IATM,JMOL,NATM)
                  PMAT(IATM,J)=PMAT(IATM,J)+PMAT(IAD,J)*FMOL(JMOL)
                END DO
              ELSE  
                PMAT(IATM,J)=PMAT(IATM,J)*FMOL(0)
              END IF
            END IF
            IF (LPIAT) THEN
              IF (NLSPCSCL_ION) THEN
                PIAT(IATM,J)=0._DP
                DO JION = 1,NIONI
                  IAD = EIRENE_INDIRECT_ADDRESS(IATM,JION,NATM)
                  PIAT(IATM,J)=PIAT(IATM,J)+PIAT(IAD,J)*FION(JION)
                END DO
              ELSE  
                PIAT(IATM,J)=PIAT(IATM,J)*FION(0)
              END IF
            END IF
            IF (LPPHAT) THEN
              IF (NLSPCSCL_PHOT) THEN
                PPHAT(IATM,J)=0._DP
                DO JPHOT = 1,NPHOTI
                  IAD = EIRENE_INDIRECT_ADDRESS(IATM,JPHOT,NATM)
                  PPHAT(IATM,J)=PPHAT(IATM,J)+PPHAT(IAD,J)*FPHOT(JPHOT)
                END DO
              ELSE  
                PPHAT(IATM,J)=PPHAT(IATM,J)*FPHOT(0)
              END IF
            END IF
            IF (LPGENA) PGENA(IATM,J)=PGENA(IATM,J)*FATM(IATM)
            IF (LEGENA) EGENA(IATM,J)=EGENA(IATM,J)*FATM(IATM)
            IF (LVGENA) VGENA(IATM,J)=VGENA(IATM,J)*FATM(IATM)
            IF (LVXDENA) VXDENA(IATM,J)=VXDENA(IATM,J)*FATM(IATM)
            IF (LVYDENA) VYDENA(IATM,J)=VYDENA(IATM,J)*FATM(IATM)
            IF (LVZDENA) VZDENA(IATM,J)=VZDENA(IATM,J)*FATM(IATM)
            IF (LRAEL) RAEL(IATM,J)=RAEL(IATM,J)*FATM(IATM)
  111     CONTINUE
c  atomic surface tallies
          DO 310 J=1,NLMPGS
            IF (LPOTAT) POTAT(IATM,J)=POTAT(IATM,J)*FATM(IATM)
            IF (LPRFAAT) THEN
              IF (NLSPCSCL_ATM) THEN
                PRFAAT(IATM,J)=0._DP
                DO JATM = 1,NATMI
                  IAD = EIRENE_INDIRECT_ADDRESS(IATM,JATM,NATM)
                  PRFAAT(IATM,J)=PRFAAT(IATM,J)+PRFAAT(IAD,J)*FATM(JATM)
                END DO
              ELSE
                PRFAAT(IATM,J)=PRFAAT(IATM,J)*FATM(0)
              END IF
            END IF
            IF (LPRFMAT) THEN
              IF (NLSPCSCL_MOL) THEN
                PRFMAT(IATM,J)=0._DP
                DO JMOL = 1,NMOLI
                  IAD = EIRENE_INDIRECT_ADDRESS(IATM,JMOL,NATM)
                  PRFMAT(IATM,J)=PRFMAT(IATM,J)+PRFMAT(IAD,J)*FMOL(JMOL)
                END DO
              ELSE
              PRFMAT(IATM,J)=PRFMAT(IATM,J)*FMOL(0)
              END IF
            END IF
            IF (LPRFIAT) THEN
              IF (NLSPCSCL_ION) THEN
                PRFIAT(IATM,J)=0._DP
                DO JION = 1,NIONI
                  IAD = EIRENE_INDIRECT_ADDRESS(IATM,JION,NATM)
                  PRFIAT(IATM,J)=PRFIAT(IATM,J)+PRFIAT(IAD,J)*FION(JION)
                END DO
              ELSE
                PRFIAT(IATM,J)=PRFIAT(IATM,J)*FION(0)
              END IF
            END IF
            IF (LPRFPHAT) THEN
              IF (NLSPCSCL_PHOT) THEN
                PRFPHAT(IATM,J)=0._DP
                DO JPHOT = 1,NPHOT
                  IAD = EIRENE_INDIRECT_ADDRESS(IATM,JPHOT,NATM)
                  PRFPHAT(IATM,J)=PRFPHAT(IATM,J)+
     .                            PRFPHAT(IAD,J)*FPHOT(JPHOT)
                END DO
              ELSE
                PRFPHAT(IATM,J)=PRFPHAT(IATM,J)*FPHOT(0)
              END IF
            END IF
            IF (LEOTAT) EOTAT(IATM,J)=EOTAT(IATM,J)*FATM(IATM)
            IF (LERFAAT) ERFAAT(IATM,J)=ERFAAT(IATM,J)*FATM(IATM)
            IF (LERFMAT) ERFMAT(IATM,J)=ERFMAT(IATM,J)*FMOL(0)
            IF (LERFIAT) ERFIAT(IATM,J)=ERFIAT(IATM,J)*FION(0)
            IF (LERFPHAT) ERFPHAT(IATM,J)=ERFPHAT(IATM,J)*FPHOT(0)
cdr sputter tallies
            IF (LSPTAAT) SPTAAT(IATM,J)=SPTAAT(IATM,J)*FATM(IATM)
            IF (LSPTMAT) SPTMAT(IATM,J)=SPTMAT(IATM,J)*FMOL(0)
            IF (LSPTIAT) SPTIAT(IATM,J)=SPTIAT(IATM,J)*FION(0)
            IF (LSPTPHAT) SPTPHAT(IATM,J)=SPTPHAT(IATM,J)*FPHOT(0)
cdr  ?? scaling with bulk flux?
            IF (LSPTPAT) SPTPAT(IATM,J)=SPTPAT(IATM,J)*FATM(IATM)
            IF (LSPUMP) SPUMP(NSPH+IATM,J)=SPUMP(NSPH+IATM,J)*FATM(IATM)
  310     CONTINUE

 2101   CONTINUE

        if (lsptatot) sptatot = sptatot * fatm(0)

c  integrated atomic tallies, both volumetric and surface-averaged
        DO 2111 IATM=1,NATMI
          if (.not.logatm(iatm,istra)) cycle
          PDENAI(IATM,ISTRA)=PDENAI(IATM,ISTRA)*FATM(IATM)
          EDENAI(IATM,ISTRA)=EDENAI(IATM,ISTRA)*FATM(IATM)
          PGENAI(IATM,ISTRA)=PGENAI(IATM,ISTRA)*FATM(IATM)
          EGENAI(IATM,ISTRA)=EGENAI(IATM,ISTRA)*FATM(IATM)
          VGENAI(IATM,ISTRA)=VGENAI(IATM,ISTRA)*FATM(IATM)
          VXDENAI(IATM,ISTRA)=VXDENAI(IATM,ISTRA)*FATM(IATM)
          VYDENAI(IATM,ISTRA)=VYDENAI(IATM,ISTRA)*FATM(IATM)
          VZDENAI(IATM,ISTRA)=VZDENAI(IATM,ISTRA)*FATM(IATM)

          POTATI(IATM,ISTRA)=POTATI(IATM,ISTRA)*FATM(IATM)
          EOTATI(IATM,ISTRA)=EOTATI(IATM,ISTRA)*FATM(IATM)
          ERFAAI(IATM,ISTRA)=ERFAAI(IATM,ISTRA)*FATM(IATM)
          ERFMAI(IATM,ISTRA)=ERFMAI(IATM,ISTRA)*FMOL(0)
          ERFIAI(IATM,ISTRA)=ERFIAI(IATM,ISTRA)*FION(0)
          ERFPHAI(IATM,ISTRA)=ERFPHAI(IATM,ISTRA)*FPHOT(0)
          SPTAATI(IATM,ISTRA)=SPTAATI(IATM,ISTRA)*FATM(IATM)
          SPTMATI(IATM,ISTRA)=SPTMATI(IATM,ISTRA)*FMOL(0)
          SPTIATI(IATM,ISTRA)=SPTIATI(IATM,ISTRA)*FION(0)
          SPTPHATI(IATM,ISTRA)=SPTPHATI(IATM,ISTRA)*FPHOT(0)
cdr  ?? scaling with bulk flux ??
          SPTPATI(IATM,ISTRA)=SPTPATI(IATM,ISTRA)*FATM(IATM)
          SPUMPI(NSPH+IATM,ISTRA)=SPUMPI(NSPH+IATM,ISTRA)*FATM(IATM)
          RAELI(IATM,ISTRA)=RAELI(IATM,ISTRA)*FATM(IATM)

          IF (NLSPCSCL_ATM) THEN
            PAATI(IATM,ISTRA)=0._DP
            PRFAAI(IATM,ISTRA)=0._DP
            DO JATM=1,NATMI
              IAD = EIRENE_INDIRECT_ADDRESS(IATM,JATM,NATM)
              PAATI(IATM,ISTRA)=PAATI(IATM,ISTRA)+
     .                          PAATI(IAD,ISTRA)*FATM(JATM)
              PRFAAI(IATM,ISTRA)=PRFAAI(IATM,ISTRA)+
     .                          PRFAAI(IAD,ISTRA)*FATM(JATM)
            END DO
          ELSE
            PAATI(IATM,ISTRA)=PAATI(IATM,ISTRA)*FATM(0)
            PRFAAI(IATM,ISTRA)=PRFAAI(IATM,ISTRA)*FATM(0)
          END IF

          IF (NLSPCSCL_MOL) THEN
            PMATI(IATM,ISTRA)=0._DP
            PRFMAI(IATM,ISTRA)=0._DP
            DO JMOL=1,NMOLI
              IAD = EIRENE_INDIRECT_ADDRESS(IATM,JMOL,NATM)
              PMATI(IATM,ISTRA)=PMATI(IATM,ISTRA)+
     .                          PMATI(IAD,ISTRA)*FMOL(JMOL)
              PRFMAI(IATM,ISTRA)=PRFMAI(IATM,ISTRA)+
     .                          PRFMAI(IAD,ISTRA)*FMOL(JMOL)
            END DO
          ELSE
            PMATI(IATM,ISTRA)=PMATI(IATM,ISTRA)*FMOL(0)
            PRFMAI(IATM,ISTRA)=PRFMAI(IATM,ISTRA)*FMOL(0)
          END IF

          IF (NLSPCSCL_ION) THEN
            PIATI(IATM,ISTRA)=0._DP
            PRFIAI(IATM,ISTRA)=0._DP
            DO JION=1,NIONI
              IAD = EIRENE_INDIRECT_ADDRESS(IATM,JION,NATM)
              PIATI(IATM,ISTRA)=PIATI(IATM,ISTRA)+
     .                          PIATI(IAD,ISTRA)*FION(JION)
              PRFIAI(IATM,ISTRA)=PRFIAI(IATM,ISTRA)+
     .                          PRFIAI(IAD,ISTRA)*FION(JION)
            END DO
          ELSE
            PIATI(IATM,ISTRA)=PIATI(IATM,ISTRA)*FION(0)
            PRFIAI(IATM,ISTRA)=PRFIAI(IATM,ISTRA)*FION(0)
          END IF

          IF (NLSPCSCL_PHOT) THEN
            PPHATI(IATM,ISTRA)=0._DP
            PRFPHAI(IATM,ISTRA)=0._DP
            DO JPHOT=1,NPHOTI
              IAD = EIRENE_INDIRECT_ADDRESS(IATM,JPHOT,NATM)
              PPHATI(IATM,ISTRA)=PPHATI(IATM,ISTRA)+
     .                          PPHATI(IAD,ISTRA)*FPHOT(JPHOT)
              PRFPHAI(IATM,ISTRA)=PRFPHAI(IATM,ISTRA)+
     .                          PRFPHAI(IAD,ISTRA)*FPHOT(JPHOT)
            END DO
          ELSE
          PPHATI(IATM,ISTRA)=PPHATI(IATM,ISTRA)*FPHOT(0)
          PRFPHAI(IATM,ISTRA)=PRFPHAI(IATM,ISTRA)*FPHOT(0)
          END IF

 2111   CONTINUE
 
       IF (LOGATM(0,ISTRA)) THEN
          PDENAI(0,ISTRA)=SUM(PDENAI(1:NATMI,ISTRA))
          EDENAI(0,ISTRA)=SUM(EDENAI(1:NATMI,ISTRA))
          PGENAI(0,ISTRA)=SUM(PGENAI(1:NATMI,ISTRA))
          EGENAI(0,ISTRA)=SUM(EGENAI(1:NATMI,ISTRA))
          VGENAI(0,ISTRA)=SUM(VGENAI(1:NATMI,ISTRA))
          VXDENAI(0,ISTRA)=SUM(VXDENAI(1:NATMI,ISTRA))
          VYDENAI(0,ISTRA)=SUM(VYDENAI(1:NATMI,ISTRA))
          VZDENAI(0,ISTRA)=SUM(VZDENAI(1:NATMI,ISTRA))

          POTATI(0,ISTRA)=SUM(POTATI(1:NATMI,ISTRA))
          EOTATI(0,ISTRA)=SUM(EOTATI(1:NATMI,ISTRA))
          ERFAAI(0,ISTRA)=SUM(ERFAAI(1:NATMI,ISTRA))
          ERFMAI(0,ISTRA)=SUM(ERFMAI(1:NATMI,ISTRA))
          ERFIAI(0,ISTRA)=SUM(ERFIAI(1:NATMI,ISTRA))
          ERFPHAI(0,ISTRA)=SUM(ERFPHAI(1:NATMI,ISTRA))
          SPTAATI(0,ISTRA)=SUM(SPTAATI(1:NATMI,ISTRA))
          SPTMATI(0,ISTRA)=SUM(SPTMATI(1:NATMI,ISTRA))
          SPTIATI(0,ISTRA)=SUM(SPTIATI(1:NATMI,ISTRA))
          SPTPHATI(0,ISTRA)=SUM(SPTPHATI(1:NATMI,ISTRA))
cdr  ?? scaling with bulk flux ??
          SPTPATI(0,ISTRA)=SUM(SPTPATI(1:NATMI,ISTRA))
          RAELI(0,ISTRA)=SUM(RAELI(1:NATMI,ISTRA))

          IF (NLSPCSCL_ATM) THEN
            PAATI(0,ISTRA)=0._DP
            PRFAAI(0,ISTRA)=0._DP
            DO JATM = 1,NATMI
              SUMP = 0._DP
              SUMPRF = 0._DP
              DO IATM = 1,NATMI
                IAD = EIRENE_INDIRECT_ADDRESS(IATM,JATM,NATM)
                SUMP = SUMP + PAATI(IAD,ISTRA)
                SUMPRF = SUMPRF + PRFAAI(IAD,ISTRA)
              END DO
              PAATI(0,ISTRA)=PAATI(0,ISTRA)+SUMP*FATM(JATM)
              PRFAAI(0,ISTRA)=PRFAAI(0,ISTRA)+SUMPRF*FATM(JATM)
            END DO
          ELSE
            PAATI(0,ISTRA)=SUM(PAATI(1:NATMI,ISTRA))
            PRFAAI(0,ISTRA)=SUM(PRFAAI(1:NATMI,ISTRA))
          END IF

          IF (NLSPCSCL_MOL) THEN
            PMATI(0,ISTRA)=0._DP
            PRFMAI(0,ISTRA)=0._DP
            DO IMOL = 1,NMOLI
              SUMP = 0._DP
              SUMPRF = 0._DP
              DO IATM = 1,NATMI
                IAD = EIRENE_INDIRECT_ADDRESS(IATM,IMOL,NATM)
                SUMP = SUMP + PMATI(IAD,ISTRA)
                SUMPRF = SUMPRF + PRFMAI(IAD,ISTRA)
              END DO
              PMATI(0,ISTRA)=PMATI(0,ISTRA)+SUMP*FMOL(IMOL)
              PRFMAI(0,ISTRA)=PRFMAI(0,ISTRA)+SUMPRF*FMOL(IMOL)
            END DO
          ELSE
            PMATI(0,ISTRA)=SUM(PMATI(1:NATMI,ISTRA))
            PRFMAI(0,ISTRA)=SUM(PRFMAI(1:NATMI,ISTRA))
          END IF

          IF (NLSPCSCL_ION) THEN
            PIATI(0,ISTRA)=0._DP
            PRFIAI(0,ISTRA)=0._DP
            DO IION = 1,NIONI
              SUMP = 0._DP
              SUMPRF = 0._DP
              DO IATM = 1,NATMI
                IAD = EIRENE_INDIRECT_ADDRESS(IATM,IION,NATM)
                SUMP = SUMP + PIATI(IAD,ISTRA)
                SUMPRF = SUMPRF + PRFIAI(IAD,ISTRA)
              END DO
              PIATI(0,ISTRA)=PIATI(0,ISTRA)+SUMP*FION(IION)
              PRFIAI(0,ISTRA)=PRFIAI(0,ISTRA)+SUMPRF*FION(IION)
            END DO
          ELSE
            PIATI(0,ISTRA)=SUM(PIATI(1:NATMI,ISTRA))
            PRFIAI(0,ISTRA)=SUM(PRFIAI(1:NATMI,ISTRA))
          END IF

          IF (NLSPCSCL_PHOT) THEN
            PPHATI(0,ISTRA)=0._DP
            PRFPHAI(0,ISTRA)=0._DP
            DO IPHOT = 1,NPHOTI
              SUMP = 0._DP
              SUMPRF = 0._DP
              DO IATM = 1,NATMI
                IAD = EIRENE_INDIRECT_ADDRESS(IATM,IPHOT,NATM)
                SUMP = SUMP + PPHATI(IAD,ISTRA)
                SUMPRF = SUMPRF + PRFPHAI(IAD,ISTRA)
              END DO
              PPHATI(0,ISTRA)=PPHATI(0,ISTRA)+SUMP*FPHOT(IPHOT)
              PRFPHAI(0,ISTRA)=PRFPHAI(0,ISTRA)+SUMPRF*FPHOT(IPHOT)
            END DO
          ELSE
            PMATI(0,ISTRA)=SUM(PMATI(1:NATMI,ISTRA))
            PRFMAI(0,ISTRA)=SUM(PRFMAI(1:NATMI,ISTRA))
          END IF
             
        END IF

        sptatti(istra) = sptatti(istra)*fatm(0)

        DO 2112 J=1,NSBOX_TAL
          IF (LEAAT) EAAT(J)=EAAT(J)*FATM(0)
          IF (LEMAT) EMAT(J)=EMAT(J)*FMOL(0)
          IF (LEIAT) EIAT(J)=EIAT(J)*FION(0)
          IF (LEPHAT) EPHAT(J)=EPHAT(J)*FPHOT(0)
 2112   CONTINUE
        EAATI(ISTRA)=EAATI(ISTRA)*FATM(0)
        EMATI(ISTRA)=EMATI(ISTRA)*FMOL(0)
        EIATI(ISTRA)=EIATI(ISTRA)*FION(0)
        EPHATI(ISTRA)=EPHATI(ISTRA)*FPHOT(0)
C
C  MOLECULE TALLIES
C
        DO 2115 IMOL=1,NMOLI
          if (.not.logmol(imol,istra)) cycle
          DO 115 J=1,NSBOX_TAL
            IF (LPDENM) PDENM(IMOL,J)=PDENM(IMOL,J)*FMOL(IMOL)
            IF (LEDENM) EDENM(IMOL,J)=EDENM(IMOL,J)*FMOL(IMOL)
            IF (LPAML) THEN
              IF (NLSPCSCL_ATM) THEN
                PAML(IMOL,J)=0._DP
                DO IATM = 1,NATMI
                  IAD = EIRENE_INDIRECT_ADDRESS(IMOL,IATM,NMOL)
                  PAML(IMOL,J)=PAML(IMOL,J)+PAML(IAD,J)*FATM(IATM)
                END DO
              ELSE  
                PAML(IMOL,J)=PAML(IMOL,J)*FATM(0)
              END IF
            END IF
            IF (LPMML) THEN
              IF (NLSPCSCL_MOL) THEN
                PMML(IMOL,J)=0._DP
                DO JMOL = 1,NMOLI
                  IAD = EIRENE_INDIRECT_ADDRESS(IMOL,JMOL,NMOL)
                  PMML(IMOL,J)=PMML(IMOL,J)+PMML(IAD,J)*FMOL(JMOL)
                END DO
              ELSE  
                PMML(IMOL,J)=PMML(IMOL,J)*FMOL(0)
              END IF
            END IF
            IF (LPIML) THEN
              IF (NLSPCSCL_ION) THEN
                PIML(IMOL,J)=0._DP
                DO IION = 1,NIONI
                  IAD = EIRENE_INDIRECT_ADDRESS(IMOL,IION,NMOL)
                  PIML(IMOL,J)=PIML(IMOL,J)+PIML(IAD,J)*FION(IION)
                END DO
              ELSE  
                PIML(IMOL,J)=PIML(IMOL,J)*FION(0)
              END IF
            END IF
            IF (LPPHML) THEN 
              IF (NLSPCSCL_PHOT) THEN
                PPHML(IMOL,J)=0._DP
                DO IPHOT = 1,NPHOTI
                  IAD = EIRENE_INDIRECT_ADDRESS(IMOL,IPHOT,NMOL)
                  PPHML(IMOL,J)=PPHML(IMOL,J)+PPHML(IAD,J)*FPHOT(IPHOT)
                END DO
              ELSE  
                PPHML(IMOL,J)=PPHML(IMOL,J)*FPHOT(0)
              END IF
            END IF
            IF (LPGENM) PGENM(IMOL,J)=PGENM(IMOL,J)*FMOL(IMOL)
            IF (LEGENM) EGENM(IMOL,J)=EGENM(IMOL,J)*FMOL(IMOL)
            IF (LVGENM) VGENM(IMOL,J)=VGENM(IMOL,J)*FMOL(IMOL)
            IF (LVXDENM) VXDENM(IMOL,J)=VXDENM(IMOL,J)*FMOL(IMOL)
            IF (LVYDENM) VYDENM(IMOL,J)=VYDENM(IMOL,J)*FMOL(IMOL)
            IF (LVZDENM) VZDENM(IMOL,J)=VZDENM(IMOL,J)*FMOL(IMOL)
            IF (LRMEL) RMEL(IMOL,J)=RMEL(IMOL,J)*FMOL(IMOL)
  115     CONTINUE
          DO 315 J=1,NLMPGS
            IF (LPOTML) POTML(IMOL,J)=POTML(IMOL,J)*FMOL(IMOL)
            IF (LPRFAML) THEN
              IF (NLSPCSCL_ATM) THEN
                PRFAML(IMOL,J)=0._DP
                DO IATM = 1,NATMI
                  IAD = EIRENE_INDIRECT_ADDRESS(IMOL,IATM,NMOL)
                  PRFAML(IMOL,J)=PRFAML(IMOL,J)+PRFAML(IAD,J)*FATM(IATM)
                END DO
              ELSE
                PRFAML(IMOL,J)=PRFAML(IMOL,J)*FATM(0)
              END IF
            END IF
            IF (LPRFMML) THEN
              IF (NLSPCSCL_MOL) THEN
                PRFMML(IMOL,J)=0._DP
                DO JMOL = 1,NMOLI
                  IAD = EIRENE_INDIRECT_ADDRESS(IMOL,JMOL,NMOL)
                  PRFMML(IMOL,J)=PRFMML(IMOL,J)+PRFMML(IAD,J)*FMOL(IMOL)
                END DO
              ELSE
                PRFMML(IMOL,J)=PRFMML(IMOL,J)*FMOL(0)
              END IF
            END IF
            IF (LPRFIML) THEN
              IF (NLSPCSCL_ION) THEN
                PRFIML(IMOL,J)=0._DP
                DO IION = 1,NIONI
                  IAD = EIRENE_INDIRECT_ADDRESS(IMOL,IION,NMOL)
                  PRFIML(IMOL,J)=PRFIML(IMOL,J)+PRFIML(IAD,J)*FION(IION)
                END DO
              ELSE
                PRFIML(IMOL,J)=PRFIML(IMOL,J)*FION(0)
              END IF
            END IF
            IF (LPRFPHML) THEN
              IF (NLSPCSCL_PHOT) THEN
                PRFPHML(IMOL,J)=0._DP
                DO IPHOT = 1,NPHOTI
                  IAD = EIRENE_INDIRECT_ADDRESS(IMOL,IPHOT,NMOL)
                  PRFPHML(IMOL,J)=PRFPHML(IMOL,J)+
     .                            PRFPHML(IAD,J)*FPHOT(IPHOT)
                END DO
              ELSE
                PRFPHML(IMOL,J)=PRFPHML(IMOL,J)*FPHOT(0)
              END IF
            END IF
            IF (LEOTML) EOTML(IMOL,J)=EOTML(IMOL,J)*FMOL(IMOL)
            IF (LERFAML) ERFAML(IMOL,J)=ERFAML(IMOL,J)*FATM(0)
            IF (LERFMML) ERFMML(IMOL,J)=ERFMML(IMOL,J)*FMOL(IMOL)
            IF (LERFIML) ERFIML(IMOL,J)=ERFIML(IMOL,J)*FION(0)
            IF (LERFPHML) ERFPHML(IMOL,J)=ERFPHML(IMOL,J)*FPHOT(0)
            IF (LSPTAML) SPTAML(IMOL,J)=SPTAML(IMOL,J)*FATM(0)
            IF (LSPTMML) SPTMML(IMOL,J)=SPTMML(IMOL,J)*FMOL(IMOL)
            IF (LSPTIML) SPTIML(IMOL,J)=SPTIML(IMOL,J)*FION(0)
            IF (LSPTPHML) SPTPHML(IMOL,J)=SPTPHML(IMOL,J)*FPHOT(0)
cdr  ?? scaling with bulk flux ??
            IF (LSPTPML) SPTPML(IMOL,J)=SPTPML(IMOL,J)*FMOL(IMOL)
            IF (LSPUMP) SPUMP(NSPA+IMOL,J)=SPUMP(NSPA+IMOL,J)*FMOL(IMOL)
  315     CONTINUE
 2115   CONTINUE

        if (lsptmtot) sptmtot = sptmtot*fmol(0)

        DO 2116 IMOL=1,NMOLI
          if (.not.logmol(imol,istra)) cycle
          PDENMI(IMOL,ISTRA)=PDENMI(IMOL,ISTRA)*FMOL(IMOL)
          EDENMI(IMOL,ISTRA)=EDENMI(IMOL,ISTRA)*FMOL(IMOL)
          PGENMI(IMOL,ISTRA)=PGENMI(IMOL,ISTRA)*FMOL(IMOL)
          EGENMI(IMOL,ISTRA)=EGENMI(IMOL,ISTRA)*FMOL(IMOL)
          VGENMI(IMOL,ISTRA)=VGENMI(IMOL,ISTRA)*FMOL(IMOL)
          VXDENMI(IMOL,ISTRA)=VXDENMI(IMOL,ISTRA)*FMOL(IMOL)
          VYDENMI(IMOL,ISTRA)=VYDENMI(IMOL,ISTRA)*FMOL(IMOL)
          VZDENMI(IMOL,ISTRA)=VZDENMI(IMOL,ISTRA)*FMOL(IMOL)

          POTMLI(IMOL,ISTRA)=POTMLI(IMOL,ISTRA)*FMOL(IMOL)
          EOTMLI(IMOL,ISTRA)=EOTMLI(IMOL,ISTRA)*FMOL(IMOL)
          ERFAMI(IMOL,ISTRA)=ERFAMI(IMOL,ISTRA)*FATM(0)
          ERFMMI(IMOL,ISTRA)=ERFMMI(IMOL,ISTRA)*FMOL(IMOL)
          ERFIMI(IMOL,ISTRA)=ERFIMI(IMOL,ISTRA)*FION(0)
          ERFPHMI(IMOL,ISTRA)=ERFPHMI(IMOL,ISTRA)*FPHOT(0)
          SPTAMLI(IMOL,ISTRA)=SPTAMLI(IMOL,ISTRA)*FATM(0)
          SPTMMLI(IMOL,ISTRA)=SPTMMLI(IMOL,ISTRA)*FMOL(IMOL)
          SPTIMLI(IMOL,ISTRA)=SPTIMLI(IMOL,ISTRA)*FION(0)
          SPTPHMLI(IMOL,ISTRA)=SPTPHMLI(IMOL,ISTRA)*FPHOT(0)
cdr  ?? scaling with bulk flux ??
          SPTPMLI(IMOL,ISTRA)=SPTPMLI(IMOL,ISTRA)*FMOL(IMOL)
          SPUMPI(NSPA+IMOL,ISTRA)=SPUMPI(NSPA+IMOL,ISTRA)*FMOL(IMOL)
          RMELI(IMOL,ISTRA)=RMELI(IMOL,ISTRA)*FMOL(IMOL)

          IF (NLSPCSCL_ATM) THEN
            PAMLI(IMOL,ISTRA)=0._DP
            PRFAMI(IMOL,ISTRA)=0._DP
            DO IATM=1,NATMI
              IAD = EIRENE_INDIRECT_ADDRESS(IMOL,IATM,NMOL)
              PAMLI(IMOL,ISTRA)=PAMLI(IMOL,ISTRA)+
     .                          PAMLI(IAD,ISTRA)*FATM(IATM)
              PRFAMI(IMOL,ISTRA)=PRFAMI(IMOL,ISTRA)+
     .                          PRFAMI(IAD,ISTRA)*FATM(IATM)
            END DO
          ELSE
            PAMLI(IMOL,ISTRA)=PAMLI(IMOL,ISTRA)*FATM(0)
            PRFAMI(IMOL,ISTRA)=PRFAMI(IMOL,ISTRA)*FATM(0)
          END IF

          IF (NLSPCSCL_MOL) THEN
            PMMLI(IMOL,ISTRA)=0._DP
            PRFMMI(IMOL,ISTRA)=0._DP
            DO JMOL=1,NMOLI
              IAD = EIRENE_INDIRECT_ADDRESS(IMOL,JMOL,NMOL)
              PMMLI(IMOL,ISTRA)=PMMLI(IMOL,ISTRA)+
     .                          PMMLI(IAD,ISTRA)*FMOL(IMOL)
              PRFMMI(IMOL,ISTRA)=PRFMMI(IMOL,ISTRA)+
     .                          PRFMMI(IAD,ISTRA)*FMOL(IMOL)
            END DO
          ELSE
            PMMLI(IMOL,ISTRA)=PMMLI(IMOL,ISTRA)*FMOL(0)
            PRFMMI(IMOL,ISTRA)=PRFMMI(IMOL,ISTRA)*FMOL(0)
          END IF

          IF (NLSPCSCL_ION) THEN
            PIMLI(IMOL,ISTRA)=0._DP
            PRFIMI(IMOL,ISTRA)=0._DP
            DO IION=1,NIONI
              IAD = EIRENE_INDIRECT_ADDRESS(IMOL,IION,NMOL)
              PIMLI(IMOL,ISTRA)=PIMLI(IMOL,ISTRA)+
     .                          PIMLI(IAD,ISTRA)*FION(IION)
              PRFIMI(IMOL,ISTRA)=PRFIMI(IION,ISTRA)+
     .                          PRFIMI(IAD,ISTRA)*FION(IION)
            END DO
          ELSE
            PIMLI(IMOL,ISTRA)=PIMLI(IMOL,ISTRA)*FION(0)
            PRFIMI(IMOL,ISTRA)=PRFIMI(IMOL,ISTRA)*FION(0)
          END IF

          IF (NLSPCSCL_PHOT) THEN
            PPHMLI(IMOL,ISTRA)=0._DP
            PRFPHMI(IMOL,ISTRA)=0._DP
            DO IPHOT=1,NPHOTI
              IAD = EIRENE_INDIRECT_ADDRESS(IMOL,IPHOT,NMOL)
              PPHMLI(IMOL,ISTRA)=PPHMLI(IMOL,ISTRA)+
     .                          PPHMLI(IAD,ISTRA)*FPHOT(IPHOT)
              PRFPHMI(IMOL,ISTRA)=PRFPHMI(IION,ISTRA)+
     .                          PRFPHMI(IAD,ISTRA)*FPHOT(IPHOT)
            END DO
          ELSE
            PPHMLI(IMOL,ISTRA)=PPHMLI(IMOL,ISTRA)*FPHOT(0)
            PRFPHMI(IMOL,ISTRA)=PRFPHMI(IMOL,ISTRA)*FPHOT(0)
          END IF

 2116   CONTINUE

        if (logmol(0,istra)) then
          PDENMI(0,ISTRA)=SUM(PDENMI(1:NMOLI,ISTRA))
          EDENMI(0,ISTRA)=SUM(EDENMI(1:NMOLI,ISTRA))
          PGENMI(0,ISTRA)=SUM(PGENMI(1:NMOLI,ISTRA))
          EGENMI(0,ISTRA)=SUM(EGENMI(1:NMOLI,ISTRA))
          VGENMI(0,ISTRA)=SUM(VGENMI(1:NMOLI,ISTRA))
          VXDENMI(0,ISTRA)=SUM(VXDENMI(1:NMOLI,ISTRA))
          VYDENMI(0,ISTRA)=SUM(VYDENMI(1:NMOLI,ISTRA))
          VZDENMI(0,ISTRA)=SUM(VZDENMI(1:NMOLI,ISTRA))

          POTMLI(0,ISTRA)=SUM(POTMLI(1:NMOLI,ISTRA))
          EOTMLI(0,ISTRA)=SUM(EOTMLI(1:NMOLI,ISTRA))
          ERFAMI(0,ISTRA)=SUM(ERFAMI(1:NMOLI,ISTRA))
          ERFMMI(0,ISTRA)=SUM(ERFMMI(1:NMOLI,ISTRA))
          ERFIMI(0,ISTRA)=SUM(ERFIMI(1:NMOLI,ISTRA))
          ERFPHMI(0,ISTRA)=SUM(ERFPHMI(1:NMOLI,ISTRA))
          SPTAMLI(0,ISTRA)=SUM(SPTAMLI(1:NMOLI,ISTRA))
          SPTMMLI(0,ISTRA)=SUM(SPTMMLI(1:NMOLI,ISTRA))
          SPTIMLI(0,ISTRA)=SUM(SPTIMLI(1:NMOLI,ISTRA))
          SPTPHMLI(0,ISTRA)=SUM(SPTPHMLI(1:NMOLI,ISTRA))
cdr  ?? scaling with bulk flux ??
          SPTPMLI(0,ISTRA)=SUM(SPTPMLI(1:NMOLI,ISTRA))
          RMELI(0,ISTRA)=SUM(RMELI(1:NMOLI,ISTRA))

          IF (NLSPCSCL_ATM) THEN
            PAMLI(0,ISTRA)=0._DP
            PRFAMI(0,ISTRA)=0._DP
            DO IATM = 1,NATMI
              SUMP = 0._DP
              SUMPRF = 0._DP
              DO IMOL = 1,NMOLI
                IAD = EIRENE_INDIRECT_ADDRESS(IMOL,IATM,NMOL)
                SUMP = SUMP + PAMLI(IAD,ISTRA)
                SUMPRF = SUMPRF + PRFAMI(IAD,ISTRA)
              END DO
              PAMLI(0,ISTRA)=PAMLI(0,ISTRA)+SUMP*FATM(IATM)
              PRFAMI(0,ISTRA)=PRFAMI(0,ISTRA)+SUMPRF*FATM(IATM)
            END DO
          ELSE
            PAMLI(0,ISTRA)=SUM(PAMLI(1:NMOLI,ISTRA))
            PRFAMI(0,ISTRA)=SUM(PRFAMI(1:NMOLI,ISTRA))
          END IF

          IF (NLSPCSCL_MOL) THEN
            PMMLI(0,ISTRA)=0._DP
            PRFMMI(0,ISTRA)=0._DP
            DO JMOL = 1,NMOLI
              SUMP = 0._DP
              SUMPRF = 0._DP
              DO IMOL = 1,NMOLI
                IAD = EIRENE_INDIRECT_ADDRESS(IMOL,JATM,NMOL)
                SUMP = SUMP + PMMLI(IAD,ISTRA)
                SUMPRF = SUMPRF + PRFMMI(IAD,ISTRA)
              END DO
              PMMLI(0,ISTRA)=PMMLI(0,ISTRA)+SUMP*FMOL(JMOL)
              PRFMMI(0,ISTRA)=PRFMMI(0,ISTRA)+SUMPRF*FMOL(JMOL)
            END DO
          ELSE
            PMMLI(0,ISTRA)=SUM(PMMLI(1:NMOLI,ISTRA))
            PRFMMI(0,ISTRA)=SUM(PRFMMI(1:NMOLI,ISTRA))
          END IF

          IF (NLSPCSCL_ION) THEN
            PIMLI(0,ISTRA)=0._DP
            PRFIMI(0,ISTRA)=0._DP
            DO IION = 1,NIONI
              SUMP = 0._DP
              SUMPRF = 0._DP
              DO IMOL = 1,NMOLI
                IAD = EIRENE_INDIRECT_ADDRESS(IMOL,IION,NMOL)
                SUMP = SUMP + PIMLI(IAD,ISTRA)
                SUMPRF = SUMPRF + PRFIMI(IAD,ISTRA)
              END DO
              PIMLI(0,ISTRA)=PIMLI(0,ISTRA)+SUMP*FION(IION)
              PRFIMI(0,ISTRA)=PRFIMI(0,ISTRA)+SUMPRF*FION(IION)
            END DO
          ELSE
            PIMLI(0,ISTRA)=SUM(PIMLI(1:NMOLI,ISTRA))
            PRFIMI(0,ISTRA)=SUM(PRFIMI(1:NMOLI,ISTRA))
          END IF

          IF (NLSPCSCL_PHOT) THEN
            PPHMLI(0,ISTRA)=0._DP
            PRFPHMI(0,ISTRA)=0._DP
            DO IPHOT = 1,NPHOTI
              SUMP = 0._DP
              SUMPRF = 0._DP
              DO IMOL = 1,NMOLI
                IAD = EIRENE_INDIRECT_ADDRESS(IMOL,IPHOT,NMOL)
                SUMP = SUMP + PPHMLI(IAD,ISTRA)
                SUMPRF = SUMPRF + PRFPHMI(IAD,ISTRA)
              END DO
              PPHMLI(0,ISTRA)=PPHMLI(0,ISTRA)+SUMP*FPHOT(IPHOT)
              PRFPHMI(0,ISTRA)=PRFPHMI(0,ISTRA)+SUMPRF*FPHOT(IPHOT)
            END DO
          ELSE
            PPHMLI(0,ISTRA)=SUM(PPHMLI(1:NMOLI,ISTRA))
            PRFPHMI(0,ISTRA)=SUM(PRFPHMI(1:NMOLI,ISTRA))
          END IF

        end if
  
        sptmtti(istra) = sptmtti(istra)*fmol(0)

        DO 2117 J=1,NSBOX_TAL
          IF (LEAML) EAML(J)=EAML(J)*FATM(0)
          IF (LEMML) EMML(J)=EMML(J)*FMOL(0)
          IF (LEIML) EIML(J)=EIML(J)*FION(0)
          IF (LEPHML) EPHML(J)=EPHML(J)*FPHOT(0)
 2117   CONTINUE
        EAMLI(ISTRA)=EAMLI(ISTRA)*FATM(0)
        EMMLI(ISTRA)=EMMLI(ISTRA)*FMOL(0)
        EIMLI(ISTRA)=EIMLI(ISTRA)*FION(0)
        EPHMLI(ISTRA)=EPHMLI(ISTRA)*FPHOT(0)
C
C  TEST ION TALLIES
C
        DO 420 IION=1,NIONI
          if (.not.logion(iion,istra)) cycle
          DO 421 J=1,NSBOX_TAL
            IF (LPDENI) PDENI(IION,J)=PDENI(IION,J)*FION(IION)
            IF (LEDENI) EDENI(IION,J)=EDENI(IION,J)*FION(IION)
            IF (LPAIO) THEN
              IF (NLSPCSCL_ATM) THEN
                PAIO(IION,J)=0._DP
                DO IATM = 1,NATMI
                  IAD = EIRENE_INDIRECT_ADDRESS(IION,IATM,NION)
                  PAIO(IION,J)=PAIO(IION,J)+PAIO(IAD,J)*FATM(IATM)
                END DO
              ELSE  
                PAIO(IION,J)=PAIO(IION,J)*FATM(0)
              END IF
            END IF
            IF (LPMIO) THEN
              IF (NLSPCSCL_MOL) THEN
                PMIO(IION,J)=0._DP
                DO IMOL = 1,NMOLI
                  IAD = EIRENE_INDIRECT_ADDRESS(IION,IMOL,NION)
                  PMIO(IION,J)=PMIO(IION,J)+PMIO(IAD,J)*FMOL(IMOL)
                END DO
              ELSE  
                PMIO(IION,J)=PMIO(IION,J)*FMOL(0)
              END IF
            END IF
            IF (LPIIO) THEN
              IF (NLSPCSCL_ION) THEN
                PIIO(IION,J)=0._DP
                DO JION = 1,NIONI
                  IAD = EIRENE_INDIRECT_ADDRESS(IION,JION,NION)
                  PIIO(IION,J)=PIIO(IION,J)+PIIO(IAD,J)*FION(IION)
                END DO
              ELSE  
                PIIO(IION,J)=PIIO(IION,J)*FION(0)
              END IF
            END IF
            IF (LPPHIO) THEN
              IF (NLSPCSCL_PHOT) THEN
                PPHIO(IION,J)=0._DP
                DO IPHOT = 1,NPHOTI
                  IAD = EIRENE_INDIRECT_ADDRESS(IION,IPHOT,NION)
                  PPHIO(IION,J)=PPHIO(IION,J)+PPHIO(IAD,J)*FPHOT(IPHOT)
                END DO
              ELSE  
                PPHIO(IION,J)=PPHIO(IION,J)*FPHOT(0)
              END IF
            END IF
            IF (LPGENI) PGENI(IION,J)=PGENI(IION,J)*FION(IION)
            IF (LEGENI) EGENI(IION,J)=EGENI(IION,J)*FION(IION)
            IF (LVGENI) VGENI(IION,J)=VGENI(IION,J)*FION(IION)
            IF (LVXDENI) VXDENI(IION,J)=VXDENI(IION,J)*FION(IION)
            IF (LVYDENI) VYDENI(IION,J)=VYDENI(IION,J)*FION(IION)
            IF (LVZDENI) VZDENI(IION,J)=VZDENI(IION,J)*FION(IION)
            IF (LRIEL) RIEL(IION,J)=RIEL(IION,J)*FION(IION)
  421     CONTINUE
          DO 422 J=1,NLMPGS
            IF (LPOTIO) POTIO(IION,J)=POTIO(IION,J)*FION(IION)
            IF (LPRFAIO) THEN
              IF (NLSPCSCL_ATM) THEN
                PRFAIO(IION,J)=0._DP
                DO IATM = 1,NATMI
                  IAD = EIRENE_INDIRECT_ADDRESS(IION,IATM,NION)
                  PRFAIO(IION,J)=PRFAIO(IION,J)+PRFAIO(IAD,J)*FATM(IATM)
                END DO
              ELSE
                PRFAIO(IION,J)=PRFAIO(IION,J)*FATM(0)
              END IF
            END IF
            IF (LPRFMIO) THEN
              IF (NLSPCSCL_MOL) THEN
                PRFMIO(IION,J)=0._DP
                DO IMOL = 1,NMOLI
                  IAD = EIRENE_INDIRECT_ADDRESS(IION,IMOL,NION)
                  PRFMIO(IION,J)=PRFMIO(IION,J)+PRFMIO(IAD,J)*FMOL(IMOL)
                END DO
              ELSE
                PRFMIO(IION,J)=PRFMIO(IION,J)*FMOL(0)
              END IF
            END IF
            IF (LPRFIIO) THEN
              IF (NLSPCSCL_ION) THEN
                PRFIIO(IION,J)=0._DP
                DO JION = 1,NIONI
                  IAD = EIRENE_INDIRECT_ADDRESS(IION,JION,NION)
                  PRFIIO(IION,J)=PRFIIO(IION,J)+PRFIIO(IAD,J)*FION(IION)
                END DO
              ELSE
                PRFIIO(IION,J)=PRFIIO(IION,J)*FION(0)
              END IF
            END IF
            IF (LPRFPHIO) THEN
              IF (NLSPCSCL_PHOT) THEN
                PRFPHIO(IION,J)=0._DP
                DO IPHOT = 1,NPHOTI
                  IAD = EIRENE_INDIRECT_ADDRESS(IION,IPHOT,NION)
                  PRFPHIO(IION,J)=PRFPHIO(IION,J)+
     .                            PRFPHIO(IAD,J)*FPHOT(IPHOT)
                END DO
              ELSE
                PRFPHIO(IION,J)=PRFPHIO(IION,J)*FPHOT(0)
              END IF
            END IF
            IF (LEOTIO) EOTIO(IION,J)=EOTIO(IION,J)*FION(IION)
            IF (LERFAIO) ERFAIO(IION,J)=ERFAIO(IION,J)*FATM(0)
            IF (LERFMIO) ERFMIO(IION,J)=ERFMIO(IION,J)*FMOL(0)
            IF (LERFIIO) ERFIIO(IION,J)=ERFIIO(IION,J)*FION(IION)
            IF (LERFPHIO) ERFPHIO(IION,J)=ERFPHIO(IION,J)*FPHOT(0)
            IF (LSPTAIO) SPTAIO(IION,J)=SPTAIO(IION,J)*FATM(0)
            IF (LSPTMIO) SPTMIO(IION,J)=SPTMIO(IION,J)*FMOL(0)
            IF (LSPTIIO) SPTIIO(IION,J)=SPTIIO(IION,J)*FION(IION)
            IF (LSPTPHIO) SPTPHIO(IION,J)=SPTPHIO(IION,J)*FPHOT(0)
cdr  ?? scaling with bulk flux ??
            IF (LSPTPIO) SPTPIO(IION,J)=SPTPIO(IION,J)*FION(IION)
            IF (LSPUMP) SPUMP(NSPAM+IION,J)=
     .                             SPUMP(NSPAM+IION,J)*FION(IION)
  422     CONTINUE
  420   CONTINUE

        if (lsptitot) sptitot = sptitot*fion(0)

        DO 431 IION=1,NIONI
          if (.not.logion(iion,istra)) cycle
          PDENII(IION,ISTRA)=PDENII(IION,ISTRA)*FION(IION)
          EDENII(IION,ISTRA)=EDENII(IION,ISTRA)*FION(IION)
          PGENII(IION,ISTRA)=PGENII(IION,ISTRA)*FION(IION)
          EGENII(IION,ISTRA)=EGENII(IION,ISTRA)*FION(IION)
          VGENII(IION,ISTRA)=VGENII(IION,ISTRA)*FION(IION)
          VXDENII(IION,ISTRA)=VXDENII(IION,ISTRA)*FION(IION)
          VYDENII(IION,ISTRA)=VYDENII(IION,ISTRA)*FION(IION)
          VZDENII(IION,ISTRA)=VZDENII(IION,ISTRA)*FION(IION)

          POTIOI(IION,ISTRA)=POTIOI(IION,ISTRA)*FION(IION)
          EOTIOI(IION,ISTRA)=EOTIOI(IION,ISTRA)*FION(IION)
          ERFAII(IION,ISTRA)=ERFAII(IION,ISTRA)*FATM(0)
          ERFMII(IION,ISTRA)=ERFMII(IION,ISTRA)*FMOL(0)
          ERFIII(IION,ISTRA)=ERFIII(IION,ISTRA)*FION(IION)
          ERFPHII(IION,ISTRA)=ERFPHII(IION,ISTRA)*FPHOT(0)
          SPTAIOI(IION,ISTRA)=SPTAIOI(IION,ISTRA)*FATM(0)
          SPTMIOI(IION,ISTRA)=SPTMIOI(IION,ISTRA)*FMOL(0)
          SPTIIOI(IION,ISTRA)=SPTIIOI(IION,ISTRA)*FION(IION)
          SPTPHIOI(IION,ISTRA)=SPTPHIOI(IION,ISTRA)*FPHOT(0)
cdr  ?? scaling with bulk flux ??
          SPTPIOI(IION,ISTRA)=SPTPIOI(IION,ISTRA)*FION(IION)
          SPUMPI(NSPAM+IION,ISTRA)=SPUMPI(NSPAM+IION,ISTRA)*FION(IION)
          RIELI(IION,ISTRA)=RIELI(IION,ISTRA)*FION(IION)

          IF (NLSPCSCL_ATM) THEN
            PAIOI(IION,ISTRA)=0._DP
            PRFAII(IION,ISTRA)=0._DP
            DO IATM=1,NATMI
              IAD = EIRENE_INDIRECT_ADDRESS(IION,IATM,NION)
              PAIOI(IION,ISTRA)=PAIOI(IION,ISTRA)+
     .                          PAIOI(IAD,ISTRA)*FATM(IATM)
              PRFAII(IION,ISTRA)=PRFAII(IION,ISTRA)+
     .                          PRFAII(IAD,ISTRA)*FATM(IATM)
            END DO
          ELSE
            PAIOI(IION,ISTRA)=PAIOI(IION,ISTRA)*FATM(0)
            PRFAII(IION,ISTRA)=PRFAII(IION,ISTRA)*FATM(0)
          END IF

          IF (NLSPCSCL_MOL) THEN
            PMIOI(IION,ISTRA)=0._DP
            PRFMII(IION,ISTRA)=0._DP
            DO IMOL=1,NMOLI
              IAD = EIRENE_INDIRECT_ADDRESS(IION,IMOL,NION)
              PMIOI(IION,ISTRA)=PMIOI(IION,ISTRA)+
     .                          PMIOI(IAD,ISTRA)*FMOL(IMOL)
              PRFMII(IION,ISTRA)=PRFMII(IION,ISTRA)+
     .                          PRFMII(IAD,ISTRA)*FMOL(IMOL)
            END DO
          ELSE
            PMIOI(IION,ISTRA)=PMIOI(IION,ISTRA)*FMOL(0)
            PRFMII(IION,ISTRA)=PRFMII(IION,ISTRA)*FMOL(0)
          END IF

          IF (NLSPCSCL_ION) THEN
            PIIOI(IION,ISTRA)=0._DP
            PRFIII(IION,ISTRA)=0._DP
            DO JION=1,NIONI
              IAD = EIRENE_INDIRECT_ADDRESS(IION,JION,NION)
              PIIOI(IION,ISTRA)=PIIOI(IION,ISTRA)+
     .                          PIIOI(IAD,ISTRA)*FION(IION)
              PRFIII(IION,ISTRA)=PRFIII(IION,ISTRA)+
     .                          PRFIII(IAD,ISTRA)*FION(IION)
            END DO
          ELSE
            PIIOI(IION,ISTRA)=PIIOI(IION,ISTRA)*FION(0)
            PRFIII(IION,ISTRA)=PRFIII(IION,ISTRA)*FION(0)
          END IF

          IF (NLSPCSCL_PHOT) THEN
            PPHIOI(IION,ISTRA)=0._DP
            PRFPHII(IION,ISTRA)=0._DP
            DO IPHOT=1,NPHOTI
              IAD = EIRENE_INDIRECT_ADDRESS(IION,IPHOT,NION)
              PPHIOI(IION,ISTRA)=PPHIOI(IION,ISTRA)+
     .                          PPHIOI(IAD,ISTRA)*FPHOT(IPHOT)
              PRFPHII(IION,ISTRA)=PRFPHII(IION,ISTRA)+
     .                          PRFPHII(IAD,ISTRA)*FPHOT(IPHOT)
            END DO
          ELSE
            PPHIOI(IION,ISTRA)=PPHIOI(IION,ISTRA)*FPHOT(0)
            PRFPHII(IION,ISTRA)=PRFPHII(IION,ISTRA)*FPHOT(0)
          END IF

  431   CONTINUE

        if (logion(0,istra)) then
          PDENII(0,ISTRA)=SUM(PDENII(1:NIONI,ISTRA))
          EDENII(0,ISTRA)=SUM(EDENII(1:NIONI,ISTRA))
          PGENII(0,ISTRA)=SUM(PGENII(1:NIONI,ISTRA))
          EGENII(0,ISTRA)=SUM(EGENII(1:NIONI,ISTRA))
          VGENII(0,ISTRA)=SUM(VGENII(1:NIONI,ISTRA))
          VXDENII(0,ISTRA)=SUM(VXDENII(1:NIONI,ISTRA))
          VYDENII(0,ISTRA)=SUM(VYDENII(1:NIONI,ISTRA))
          VZDENII(0,ISTRA)=SUM(VZDENII(1:NIONI,ISTRA))

          POTIOI(0,ISTRA)=SUM(POTIOI(1:NIONI,ISTRA))
          EOTIOI(0,ISTRA)=SUM(EOTIOI(1:NIONI,ISTRA))
          ERFAII(0,ISTRA)=SUM(ERFAII(1:NIONI,ISTRA))
          ERFMII(0,ISTRA)=SUM(ERFMII(1:NIONI,ISTRA))
          ERFIII(0,ISTRA)=SUM(ERFIII(1:NIONI,ISTRA))
          ERFPHII(0,ISTRA)=SUM(ERFPHII(1:NIONI,ISTRA))
          SPTAIOI(0,ISTRA)=SUM(SPTAIOI(1:NIONI,ISTRA))
          SPTMIOI(0,ISTRA)=SUM(SPTMIOI(1:NIONI,ISTRA))
          SPTIIOI(0,ISTRA)=SUM(SPTIIOI(1:NIONI,ISTRA))
          SPTPHIOI(0,ISTRA)=SUM(SPTPHIOI(1:NIONI,ISTRA))
cdr  ?? scaling with bulk flux ??
          SPTPIOI(0,ISTRA)=SUM(SPTPIOI(1:NIONI,ISTRA))
          RIELI(0,ISTRA)=SUM(RIELI(1:NIONI,ISTRA))

          IF (NLSPCSCL_ATM) THEN
            PAIOI(0,ISTRA)=0._DP
            PRFAII(0,ISTRA)=0._DP
            DO IATM = 1,NATMI
              SUMP = 0._DP
              SUMPRF = 0._DP
              DO IION = 1,NIONI
                IAD = EIRENE_INDIRECT_ADDRESS(IION,IATM,NION)
                SUMP = SUMP + PAIOI(IAD,ISTRA)
                SUMPRF = SUMPRF + PRFAII(IAD,ISTRA)
              END DO
              PAIOI(0,ISTRA)=PAIOI(0,ISTRA)+SUMP*FATM(IATM)
              PRFAII(0,ISTRA)=PRFAII(0,ISTRA)+SUMPRF*FATM(IATM)
            END DO
          ELSE
            PAIOI(0,ISTRA)=SUM(PAIOI(1:NIONI,ISTRA))
            PRFAII(0,ISTRA)=SUM(PRFAII(1:NIONI,ISTRA))
          END IF

          IF (NLSPCSCL_MOL) THEN
            PMIOI(0,ISTRA)=0._DP
            PRFMII(0,ISTRA)=0._DP
            DO IMOL = 1,NMOLI
              SUMP = 0._DP
              SUMPRF = 0._DP
              DO IION = 1,NIONI
                IAD = EIRENE_INDIRECT_ADDRESS(IION,IMOL,NION)
                SUMP = SUMP + PMIOI(IAD,ISTRA)
                SUMPRF = SUMPRF + PRFMII(IAD,ISTRA)
              END DO
              PMIOI(0,ISTRA)=PMIOI(0,ISTRA)+SUMP*FMOL(IMOL)
              PRFMII(0,ISTRA)=PRFMII(0,ISTRA)+SUMPRF*FMOL(IMOL)
            END DO
          ELSE
            PMIOI(0,ISTRA)=SUM(PMIOI(1:NIONI,ISTRA))
            PRFMII(0,ISTRA)=SUM(PRFMII(1:NIONI,ISTRA))
          END IF

          IF (NLSPCSCL_ION) THEN
            PIIOI(0,ISTRA)=0._DP
            PRFIII(0,ISTRA)=0._DP
            DO JION = 1,NIONI
              SUMP = 0._DP
              SUMPRF = 0._DP
              DO IION = 1,NIONI
                IAD = EIRENE_INDIRECT_ADDRESS(IION,JION,NION)
                SUMP = SUMP + PIIOI(IAD,ISTRA)
                SUMPRF = SUMPRF + PRFIII(IAD,ISTRA)
              END DO
              PIIOI(0,ISTRA)=PIIOI(0,ISTRA)+SUMP*FION(IION)
              PRFIII(0,ISTRA)=PRFIII(0,ISTRA)+SUMPRF*FION(IION)
            END DO
          ELSE
            PIIOI(0,ISTRA)=SUM(PIIOI(1:NIONI,ISTRA))
            PRFIII(0,ISTRA)=SUM(PRFIII(1:NIONI,ISTRA))
          END IF

          IF (NLSPCSCL_PHOT) THEN
            PPHIOI(0,ISTRA)=0._DP
            PRFPHII(0,ISTRA)=0._DP
            DO IPHOT = 1,NPHOTI
              SUMP = 0._DP
              SUMPRF = 0._DP
              DO IION = 1,NIONI
                IAD = EIRENE_INDIRECT_ADDRESS(IION,IPHOT,NION)
                SUMP = SUMP + PPHIOI(IAD,ISTRA)
                SUMPRF = SUMPRF + PRFPHII(IAD,ISTRA)
              END DO
              PPHIOI(0,ISTRA)=PPHIOI(0,ISTRA)+SUMP*FPHOT(IPHOT)
              PRFPHII(0,ISTRA)=PRFPHII(0,ISTRA)+SUMPRF*FPHOT(IPHOT)
            END DO
          ELSE
            PPHIOI(0,ISTRA)=SUM(PPHIOI(1:NIONI,ISTRA))
            PRFPHII(0,ISTRA)=SUM(PRFPHII(1:NIONI,ISTRA))
          END IF

        end if

        sptitti(istra) = sptitti(istra)*fion(0)

        DO 432 J=1,NSBOX_TAL
          IF (LEAIO) EAIO(J)=EAIO(J)*FATM(0)
          IF (LEMIO) EMIO(J)=EMIO(J)*FMOL(0)
          IF (LEIIO) EIIO(J)=EIIO(J)*FION(0)
          IF (LEPHIO) EPHIO(J)=EPHIO(J)*FPHOT(0)
  432   CONTINUE
        EAIOI(ISTRA)=EAIOI(ISTRA)*FATM(0)
        EMIOI(ISTRA)=EMIOI(ISTRA)*FMOL(0)
        EIIOI(ISTRA)=EIIOI(ISTRA)*FION(0)
        EPHIOI(ISTRA)=EPHIOI(ISTRA)*FPHOT(0)
C
C  PHOTON TALLIES
C
        DO IPHOT=1,NPHOTI
          if (.not.logphot(iphot,istra)) cycle
          DO J=1,NSBOX_TAL
            IF (LPDENPH) PDENPH(IPHOT,J)=PDENPH(IPHOT,J)*FPHOT(IPHOT)
            IF (LEDENPH) EDENPH(IPHOT,J)=EDENPH(IPHOT,J)*FPHOT(IPHOT)
            IF (LPAPHT) THEN
              IF (NLSPCSCL_ATM) THEN
                PAPHT(IPHOT,J)=0._DP
                DO IATM = 1,NATMI
                  IAD = EIRENE_INDIRECT_ADDRESS(IPHOT,IATM,NPHOT)
                  PAPHT(IPHOT,J)=PAPHT(IPHOT,J)+PAPHT(IAD,J)*FATM(IATM)
                END DO
              ELSE  
                PAPHT (IPHOT,J)=PAPHT (IPHOT,J)*FATM(0)
              END IF
            END IF
            IF (LPMPHT) THEN
              IF (NLSPCSCL_MOL) THEN
                PMPHT(IPHOT,J)=0._DP
                DO IMOL = 1,NMOLI
                  IAD = EIRENE_INDIRECT_ADDRESS(IPHOT,IMOL,NPHOT)
                  PMPHT(IPHOT,J)=PMPHT(IPHOT,J)+PMPHT(IAD,J)*FMOL(IMOL)
                END DO
              ELSE  
                PMPHT (IPHOT,J)=PMPHT (IPHOT,J)*FMOL(0)
              END IF
            END IF
            IF (LPIPHT) THEN
              IF (NLSPCSCL_ION) THEN
                PIPHT(IPHOT,J)=0._DP
                DO IION = 1,NIONI
                  IAD = EIRENE_INDIRECT_ADDRESS(IPHOT,IION,NPHOT)
                  PIPHT(IPHOT,J)=PIPHT(IPHOT,J)+PIPHT(IAD,J)*FION(IION)
                END DO
              ELSE  
                PIPHT (IPHOT,J)=PIPHT (IPHOT,J)*FION(0)
              END IF
            END IF
            IF (LPPHPHT) THEN
              IF (NLSPCSCL_PHOT) THEN
                PPHPHT(IPHOT,J)=0._DP
                DO JPHOT = 1,NPHOTI
                  IAD = EIRENE_INDIRECT_ADDRESS(IPHOT,JPHOT,NPHOT)
                  PPHPHT(IPHOT,J)=PPHPHT(IPHOT,J)+
     .                            PPHPHT(IAD,J)*FPHOT(JPHOT)
                END DO
              ELSE  
                PPHPHT(IPHOT,J)=PPHPHT(IPHOT,J)*FPHOT(0)
              END IF
            END IF
            IF (LPGENPH) PGENPH(IPHOT,J)=PGENPH(IPHOT,J)*FPHOT(IPHOT)
            IF (LEGENPH) EGENPH(IPHOT,J)=EGENPH(IPHOT,J)*FPHOT(IPHOT)
            IF (LVGENPH) VGENPH(IPHOT,J)=VGENPH(IPHOT,J)*FPHOT(IPHOT)
            IF (LVXDENPH) VXDENPH(IPHOT,J)=VXDENPH(IPHOT,J)*FPHOT(IPHOT)
            IF (LVYDENPH) VYDENPH(IPHOT,J)=VYDENPH(IPHOT,J)*FPHOT(IPHOT)
            IF (LVZDENPH) VZDENPH(IPHOT,J)=VZDENPH(IPHOT,J)*FPHOT(IPHOT)
          END DO
C  SURFACE-AVERAGED TALLIES
          DO J=1,NLMPGS
            IF (LPOTPHT)  POTPHT (IPHOT,J)=POTPHT (IPHOT,J)*FPHOT(IPHOT)
            IF (LPRFAPHT) THEN 
              IF (NLSPCSCL_ATM) THEN
                PRFAPHT(IPHOT,J)=0._DP
                DO IATM = 1,NATMI
                  IAD = EIRENE_INDIRECT_ADDRESS(IPHOT,IATM,NPHOT)
                  PRFAPHT(IPHOT,J)=PRFAPHT(IPHOT,J)+
     .                             PRFAPHT(IAD,J)*FATM(IATM)
                END DO
              ELSE
                PRFAPHT(IPHOT,J)=PRFAPHT(IPHOT,J)*FATM(0)
              END IF
            END IF
            IF (LPRFMPHT) THEN
              IF (NLSPCSCL_MOL) THEN
                PRFMPHT(IPHOT,J)=0._DP
                DO IMOL = 1,NMOLI
                  IAD = EIRENE_INDIRECT_ADDRESS(IPHOT,IMOL,NPHOT)
                  PRFMPHT(IPHOT,J)=PRFMPHT(IPHOT,J)+
     .                             PRFMPHT(IAD,J)*FMOL(IMOL)
                END DO
              ELSE
                PRFMPHT(IPHOT,J)=PRFMPHT(IPHOT,J)*FMOL(0)
              END IF
            END IF
            IF (LPRFIPHT) THEN
              IF (NLSPCSCL_ION) THEN
                PRFIPHT(IPHOT,J)=0._DP
                DO IION = 1,NIONI
                  IAD = EIRENE_INDIRECT_ADDRESS(IPHOT,IION,NPHOT)
                  PRFIPHT(IPHOT,J)=PRFIPHT(IPHOT,J)+
     .                             PRFIPHT(IAD,J)*FION(IION)
                END DO
              ELSE
                PRFIPHT(IPHOT,J)=PRFIPHT(IPHOT,J)*FION(0)
              END IF
            END IF
            IF (LPRFPHPHT) THEN
              IF (NLSPCSCL_PHOT) THEN
                PRFPHPHT(IPHOT,J)=0._DP
                DO JPHOT = 1,NPHOTI
                  IAD = EIRENE_INDIRECT_ADDRESS(IPHOT,JPHOT,NPHOT)
                  PRFPHPHT(IPHOT,J)=PRFPHPHT(IPHOT,J)+
     .                              PRFPHPHT(IAD,J)*FPHOT(JPHOT)
                END DO
              ELSE
                PRFPHPHT(IPHOT,J)=PRFPHPHT(IPHOT,J)*FPHOT(0)
              END IF
            END IF
            IF (LPRFPPHT) PRFPPHT(IPHOT,J)=PRFPPHT(IPHOT,J)*FPHOT(IPHOT)
            IF (LEOTPHT) EOTPHT (IPHOT,J)=EOTPHT (IPHOT,J)*FPHOT(IPHOT)
            IF (LERFAPHT) ERFAPHT(IPHOT,J)=ERFAPHT(IPHOT,J)*FATM(0)
            IF (LERFMPHT) ERFMPHT(IPHOT,J)=ERFMPHT(IPHOT,J)*FMOL(0)
            IF (LERFIPHT) ERFIPHT(IPHOT,J)=ERFIPHT(IPHOT,J)*FION(0)
            IF (LERFPHPHT) ERFPHPHT(IPHOT,J)=
     .                                    ERFPHPHT(IPHOT,J)*FPHOT(IPHOT)
            IF (LERFPPHT) ERFPPHT(IPHOT,J)=ERFPPHT(IPHOT,J)*FPHOT(IPHOT)
            IF (LSPTAPHT) SPTAPHT(IPHOT,J)=SPTAPHT(IPHOT,J)*FATM(0)
            IF (LSPTMPHT) SPTMPHT(IPHOT,J)=SPTMPHT(IPHOT,J)*FMOL(0)
            IF (LSPTIPHT) SPTIPHT(IPHOT,J)=SPTIPHT(IPHOT,J)*FION(0)
            IF (LSPTPHPHT) SPTPHPHT(IPHOT,J)=
     .                                    SPTPHPHT(IPHOT,J)*FPHOT(IPHOT)
cdr  ?? scaling with bulk flux ??
            IF (LSPTPPHT) SPTPPHT(IPHOT,J)=SPTPPHT(IPHOT,J)*FPHOT(IPHOT)
            IF (LSPUMP) SPUMP(IPHOT,J)=SPUMP(IPHOT,J)*FPHOT(IPHOT)
          END DO
        END DO

        if (lsptphtot) sptphtot = sptphtot*fphot(0)

        DO IPHOT=1,NPHOTI
          if (.not.logphot(iphot,istra)) cycle
          PDENPHI(IPHOT,ISTRA)=PDENPHI(IPHOT,ISTRA)*FPHOT(IPHOT)
          EDENPHI(IPHOT,ISTRA)=EDENPHI(IPHOT,ISTRA)*FPHOT(IPHOT)
          PGENPHI(IPHOT,ISTRA)=PGENPHI(IPHOT,ISTRA)*FPHOT(IPHOT)
          EGENPHI(IPHOT,ISTRA)=EGENPHI(IPHOT,ISTRA)*FPHOT(IPHOT)
          VGENPHI(IPHOT,ISTRA)=VGENPHI(IPHOT,ISTRA)*FPHOT(IPHOT)
          VXDENPHI(IPHOT,ISTRA)=VXDENPHI(IPHOT,ISTRA)*FPHOT(IPHOT)
          VYDENPHI(IPHOT,ISTRA)=VYDENPHI(IPHOT,ISTRA)*FPHOT(IPHOT)
          VZDENPHI(IPHOT,ISTRA)=VZDENPHI(IPHOT,ISTRA)*FPHOT(IPHOT)

          POTPHTI(IPHOT,ISTRA)=POTPHTI(IPHOT,ISTRA)*FPHOT(IPHOT)
          EOTPHTI(IPHOT,ISTRA)=EOTPHTI(IPHOT,ISTRA)*FPHOT(IPHOT)
          ERFAPHTI(IPHOT,ISTRA)=ERFAPHTI(IPHOT,ISTRA)*FATM(0)
          ERFMPHTI(IPHOT,ISTRA)=ERFMPHTI(IPHOT,ISTRA)*FMOL(0)
          ERFIPHTI(IPHOT,ISTRA)=ERFIPHTI(IPHOT,ISTRA)*FION(0)
          ERFPHPHTI(IPHOT,ISTRA)=ERFPHPHTI(IPHOT,ISTRA)*FPHOT(IPHOT)
          SPTAPHTI(IPHOT,ISTRA)=SPTAPHTI(IPHOT,ISTRA)*FATM(0)
          SPTMPHTI(IPHOT,ISTRA)=SPTMPHTI(IPHOT,ISTRA)*FMOL(0)
          SPTIPHTI(IPHOT,ISTRA)=SPTIPHTI(IPHOT,ISTRA)*FION(0)
          SPTPHPHTI(IPHOT,ISTRA)=SPTPHPHTI(IPHOT,ISTRA)*FPHOT(IPHOT)
cdr  ?? scaling with bulk flux ??
          SPTPPHTI(IPHOT,ISTRA)=SPTPPHTI(IPHOT,ISTRA)*FPHOT(IPHOT)
          SPUMPI(IPHOT,ISTRA)=SPUMPI(IPHOT,ISTRA)*FPHOT(IPHOT)

          IF (NLSPCSCL_ATM) THEN
            PAPHTI(IPHOT,ISTRA)=0._DP
            PRFAPHTI(IPHOT,ISTRA)=0._DP
            DO IATM=1,NATMI
              IAD = EIRENE_INDIRECT_ADDRESS(IPHOT,IATM,NPHOT)
              PAPHTI(IPHOT,ISTRA)=PAPHTI(IPHOT,ISTRA)+
     .                          PAPHTI(IAD,ISTRA)*FATM(IATM)
              PRFAPHTI(IPHOT,ISTRA)=PRFAPHTI(IPHOT,ISTRA)+
     .                          PRFAPHTI(IAD,ISTRA)*FATM(IATM)
            END DO
          ELSE
            PAPHTI (IPHOT,ISTRA)=PAPHTI (IPHOT,ISTRA)*FATM(0)
            PRFAPHTI(IPHOT,ISTRA)=PRFAPHTI(IPHOT,ISTRA)*FATM(0)
          END IF

          IF (NLSPCSCL_MOL) THEN
            PMPHTI(IPHOT,ISTRA)=0._DP
            PRFMPHTI(IPHOT,ISTRA)=0._DP
            DO IMOL=1,NMOLI
              IAD = EIRENE_INDIRECT_ADDRESS(IPHOT,IMOL,NPHOT)
              PMPHTI(IPHOT,ISTRA)=PMPHTI(IPHOT,ISTRA)+
     .                          PMPHTI(IAD,ISTRA)*FMOL(IMOL)
              PRFMPHTI(IPHOT,ISTRA)=PRFMPHTI(IPHOT,ISTRA)+
     .                          PRFMPHTI(IAD,ISTRA)*FMOL(IMOL)
            END DO
          ELSE
            PMPHTI (IPHOT,ISTRA)=PMPHTI (IPHOT,ISTRA)*FMOL(0)
            PRFMPHTI(IPHOT,ISTRA)=PRFMPHTI(IPHOT,ISTRA)*FMOL(0)
          END IF

          IF (NLSPCSCL_ION) THEN
            PIPHTI(IPHOT,ISTRA)=0._DP
            PRFIPHTI(IPHOT,ISTRA)=0._DP
            DO IION=1,NIONI
              IAD = EIRENE_INDIRECT_ADDRESS(IPHOT,IION,NPHOT)
              PIPHTI(IPHOT,ISTRA)=PIPHTI(IPHOT,ISTRA)+
     .                          PIPHTI(IAD,ISTRA)*FION(IION)
              PRFIPHTI(IPHOT,ISTRA)=PRFIPHTI(IPHOT,ISTRA)+
     .                          PRFIPHTI(IAD,ISTRA)*FION(IION)
            END DO
          ELSE
            PIPHTI (IPHOT,ISTRA)=PIPHTI (IPHOT,ISTRA)*FION(0)
            PRFIPHTI(IPHOT,ISTRA)=PRFIPHTI(IPHOT,ISTRA)*FION(0)
          END IF

          IF (NLSPCSCL_PHOT) THEN
            PPHPHTI(IPHOT,ISTRA)=0._DP
            PRFPHPHTI(IPHOT,ISTRA)=0._DP
            DO JPHOT=1,NPHOTI
              IAD = EIRENE_INDIRECT_ADDRESS(IPHOT,JPHOT,NPHOT)
              PPHPHTI(IPHOT,ISTRA)=PPHPHTI(IPHOT,ISTRA)+
     .                          PPHPHTI(IAD,ISTRA)*FPHOT(IPHOT)
              PRFPHPHTI(IPHOT,ISTRA)=PRFPHPHTI(IPHOT,ISTRA)+
     .                          PRFPHPHTI(IAD,ISTRA)*FPHOT(IPHOT)
            END DO
          ELSE
            PPHPHTI(IPHOT,ISTRA)=PPHPHTI(IPHOT,ISTRA)*FPHOT(0)
            PRFPHPHTI(IPHOT,ISTRA)=PRFPHPHTI(IPHOT,ISTRA)*FPHOT(0)
          END IF

        END DO
        if (logphot(0,istra)) then
          PDENPHI(0,ISTRA)=SUM(PDENPHI(1:NPHOTI,ISTRA))
          EDENPHI(0,ISTRA)=SUM(EDENPHI(1:NPHOTI,ISTRA))
          PGENPHI(0,ISTRA)=SUM(PGENPHI(1:NPHOTI,ISTRA))
          EGENPHI(0,ISTRA)=SUM(EGENPHI(1:NPHOTI,ISTRA))
          VGENPHI(0,ISTRA)=SUM(VGENPHI(1:NPHOTI,ISTRA))
          VXDENPHI(0,ISTRA)=SUM(VXDENPHI(1:NPHOTI,ISTRA))
          VYDENPHI(0,ISTRA)=SUM(VYDENPHI(1:NPHOTI,ISTRA))
          VZDENPHI(0,ISTRA)=SUM(VZDENPHI(1:NPHOTI,ISTRA))

          POTPHTI(0,ISTRA)=SUM(POTPHTI(1:NPHOTI,ISTRA))
          EOTPHTI(0,ISTRA)=SUM(EOTPHTI(1:NPHOTI,ISTRA))
          ERFAPHTI(0,ISTRA)=SUM(ERFAPHTI(1:NPHOTI,ISTRA))
          ERFMPHTI(0,ISTRA)=SUM(ERFMPHTI(1:NPHOTI,ISTRA))
          ERFIPHTI(0,ISTRA)=SUM(ERFIPHTI(1:NPHOTI,ISTRA))
          ERFPHPHTI(0,ISTRA)=SUM(ERFPHPHTI(1:NPHOTI,ISTRA))
          SPTAPHTI(0,ISTRA)=SUM(SPTAPHTI(1:NPHOTI,ISTRA))
          SPTMPHTI(0,ISTRA)=SUM(SPTMPHTI(1:NPHOTI,ISTRA))
          SPTIPHTI(0,ISTRA)=SUM(SPTIPHTI(1:NPHOTI,ISTRA))
          SPTPHPHTI(0,ISTRA)=SUM(SPTPHPHTI(1:NPHOTI,ISTRA))
cdr  ?? scaling with bulk flux ??
          SPTPPHTI(0,ISTRA)=SUM(SPTPPHTI(1:NPHOTI,ISTRA))

          IF (NLSPCSCL_ATM) THEN
            PAPHTI(0,ISTRA)=0._DP
            PRFAPHTI(0,ISTRA)=0._DP
            DO IATM = 1,NATMI
              SUMP = 0._DP
              SUMPRF = 0._DP
              DO IPHOT = 1,NPHOTI
                IAD = EIRENE_INDIRECT_ADDRESS(IPHOT,IATM,NPHOT)
                SUMP = SUMP + PAPHTI(IAD,ISTRA)
                SUMPRF = SUMPRF + PRFAPHTI(IAD,ISTRA)
              END DO
              PAPHTI(0,ISTRA)=PAPHTI(0,ISTRA)+SUMP*FATM(IATM)
              PRFAPHTI(0,ISTRA)=PRFAPHTI(0,ISTRA)+SUMPRF*FATM(IATM)
            END DO
          ELSE
            PAPHTI(0,ISTRA)=SUM(PAPHTI(1:NPHOTI,ISTRA))
            PRFAPHTI(0,ISTRA)=SUM(PRFAPHTI(1:NPHOTI,ISTRA))
          END IF

          IF (NLSPCSCL_MOL) THEN
            PMPHTI(0,ISTRA)=0._DP
            PRFMPHTI(0,ISTRA)=0._DP
            DO IMOL = 1,NMOLI
              SUMP = 0._DP
              SUMPRF = 0._DP
              DO IPHOT = 1,NPHOTI
                IAD = EIRENE_INDIRECT_ADDRESS(IPHOT,IMOL,NPHOT)
                SUMP = SUMP + PMPHTI(IAD,ISTRA)
                SUMPRF = SUMPRF + PRFMPHTI(IAD,ISTRA)
              END DO
              PMPHTI(0,ISTRA)=PMPHTI(0,ISTRA)+SUMP*FMOL(IMOL)
              PRFMPHTI(0,ISTRA)=PRFMPHTI(0,ISTRA)+SUMPRF*FMOL(IMOL)
            END DO
          ELSE
            PMPHTI(0,ISTRA)=SUM(PMPHTI(1:NPHOTI,ISTRA))
            PRFMPHTI(0,ISTRA)=SUM(PRFMPHTI(1:NPHOTI,ISTRA))
          END IF
          IF (NLSPCSCL_ION) THEN
            PIPHTI(0,ISTRA)=0._DP
            PRFIPHTI(0,ISTRA)=0._DP
            DO IION = 1,NIONI
              SUMP = 0._DP
              SUMPRF = 0._DP
              DO IPHOT = 1,NPHOTI
                IAD = EIRENE_INDIRECT_ADDRESS(IPHOT,IION,NPHOT)
                SUMP = SUMP + PIPHTI(IAD,ISTRA)
                SUMPRF = SUMPRF + PRFIPHTI(IAD,ISTRA)
              END DO
              PIPHTI(0,ISTRA)=PIPHTI(0,ISTRA)+SUMP*FION(IION)
              PRFIPHTI(0,ISTRA)=PRFIPHTI(0,ISTRA)+SUMPRF*FION(IION)
            END DO
          ELSE
            PIPHTI(0,ISTRA)=SUM(PIPHTI(1:NPHOTI,ISTRA))
            PRFIPHTI(0,ISTRA)=SUM(PRFIPHTI(1:NPHOTI,ISTRA))
          END IF
          IF (NLSPCSCL_PHOT) THEN
            PPHPHTI(0,ISTRA)=0._DP
            PRFPHPHTI(0,ISTRA)=0._DP
            DO JPHOT = 1,NPHOTI
              SUMP = 0._DP
              SUMPRF = 0._DP
              DO IPHOT = 1,NPHOTI
                IAD = EIRENE_INDIRECT_ADDRESS(IPHOT,JPHOT,NPHOT)
                SUMP = SUMP + PPHPHTI(IAD,ISTRA)
                SUMPRF = SUMPRF + PRFPHPHTI(IAD,ISTRA)
              END DO
              PPHPHTI(0,ISTRA)=PPHPHTI(0,ISTRA)+SUMP*FPHOT(JPHOT)
              PRFPHPHTI(0,ISTRA)=PRFPHPHTI(0,ISTRA)+SUMPRF*FPHOT(JPHOT)
            END DO
          ELSE
            PPHPHTI(0,ISTRA)=SUM(PPHPHTI(1:NPHOTI,ISTRA))
            PRFPHPHTI(0,ISTRA)=SUM(PRFPHPHTI(1:NPHOTI,ISTRA))
          END IF

        END IF

        SPUMPI(0,ISTRA) = SUM(SPUMPI(1:NSPAM+NIONI,ISTRA))

        sptphtti(istra) = sptphtti(istra)*fphot(0)

        DO J=1,NSBOX_TAL
          IF (LEAPHT) EAPHT(J)=EAPHT(J)*FATM(0)
          IF (LEMPHT) EMPHT(J)=EMPHT(J)*FMOL(0)
          IF (LEIPHT) EIPHT(J)=EIPHT(J)*FION(0)
          IF (LEPHPHT) EPHPHT(J)=EPHPHT(J)*FPHOT(0)
        END DO
        EAPHTI(ISTRA)=EAPHTI(ISTRA)*FATM(0)
        EMPHTI(ISTRA)=EMPHTI(ISTRA)*FMOL(0)
        EIPHTI(ISTRA)=EIPHTI(ISTRA)*FION(0)
        EPHPHTI(ISTRA)=EPHPHTI(ISTRA)*FPHOT(0)

csw integrate/scale internal photon tallies
cdr     IF (NPHOTI > 0)
cdr  .    CALL PH_INTEGRATE(ISTRA,FLXFAC(ISTRA)/ELCHA,
cdr  .                    fphot,fatm,fmol,fion)
C
C  ADDITIONAL TRACKLENGTH-ESTIMATED TALLIES
C
        IF (LADDV) THEN
          DO 423 IADV=1,NADVI
            IF (IADRC(IADV).EQ.0) THEN
              FADD=FPHOT(IADVS(IADV))
            ELSEIF (IADRC(IADV).EQ.1) THEN
              FADD=FATM(IADVS(IADV))
            ELSEIF (IADRC(IADV).EQ.2) THEN
              FADD=FMOL(IADVS(IADV))
            ELSEIF (IADRC(IADV).EQ.3) THEN
              FADD=FION(IADVS(IADV))
            ELSE
              GOTO 423
            ENDIF
            DO 424 J=1,NSBOX_TAL
              ADDV(IADV,J)=ADDV(IADV,J)*FADD
  424       CONTINUE
            ADDVI(IADV,ISTRA)=ADDVI(IADV,ISTRA)*FADD
  423     CONTINUE
        END IF
C
C  ADDITIONAL COLLISION-ESTIMATED TALLIES
C
        IF (LCOLV) THEN
          DO 426 ICLV=1,NCLVI
            IF (ICLRC(ICLV).EQ.0) THEN
              FADD=FPHOT(ICLVS(ICLV))
            ELSEIF (ICLRC(ICLV).EQ.1) THEN
              FADD=FATM(ICLVS(ICLV))
            ELSEIF (ICLRC(ICLV).EQ.2) THEN
              FADD=FMOL(ICLVS(ICLV))
            ELSEIF (ICLRC(ICLV).EQ.3) THEN
              FADD=FION(ICLVS(ICLV))
            ELSE
              GOTO 426
            ENDIF
            DO 427 J=1,NSBOX_TAL
              COLV(ICLV,J)=COLV(ICLV,J)*FADD
  427       CONTINUE
            COLVI(ICLV,ISTRA)=COLVI(ICLV,ISTRA)*FADD
  426     CONTINUE
        END IF
C
C  SNAPSHOT TALLIES
C
C  TO BE WRITTEN
C
C  TALLIES FOR COUPLING TO FLUID PLASMA CODE
C
        IF (LCOPV) THEN
          DO 435 ICPV=1,NCPVI
            IF (ICPRC(ICPV).EQ.1) THEN
              FADD=FATM(ICPVS(ICPV))
            ELSEIF (ICPRC(ICPV).EQ.2) THEN
              FADD=FMOL(ICPVS(ICPV))
            ELSEIF (ICPRC(ICPV).EQ.3) THEN
              FADD=FION(ICPVS(ICPV))
            ELSEIF (ICPRC(ICPV).EQ.0) THEN
              FADD=FPHOT(ICPVS(ICPV))
            ELSE
              GOTO 435
            ENDIF
            DO 436 J=1,NSBOX_TAL
              COPV(ICPV,J)=COPV(ICPV,J)*FADD
  436       CONTINUE
            COPVI(ICPV,ISTRA)=COPVI(ICPV,ISTRA)*FADD
  435     CONTINUE
        END IF
C
C  TALLIES FOR BGK SELF-COLLISION ITERATIONS
C
        IF (LBGKV) THEN
          DO 437 IBGV=1,NBGVI
            IF (IBGRC(IBGV).EQ.1) THEN
              FADD=FATM(IBGVS(IBGV))
            ELSEIF (IBGRC(IBGV).EQ.2) THEN
              FADD=FMOL(IBGVS(IBGV))
            ELSEIF (IBGRC(IBGV).EQ.3) THEN
              FADD=FION(IBGVS(IBGV))
            ELSEIF (IBGRC(IBGV).EQ.0) THEN
              FADD=FPHOT(IBGVS(IBGV))
            ELSE
              GOTO 437
            ENDIF
            DO 438 J=1,NSBOX_TAL
              BGKV(IBGV,J)=BGKV(IBGV,J)*FADD
  438       CONTINUE
            BGKVI(IBGV,ISTRA)=BGKVI(IBGV,ISTRA)*FADD
  437     CONTINUE
        END IF
C
C  BULK ION TALLIES
C
        DO IPLS=1,NPLSI
          if (.not.logpls(ipls,istra)) cycle
          DO J=1,NSBOX_TAL
            IF (LPAPL) THEN
              IF (NLSPCSCL_ATM) THEN
                PAPL(IPLS,J)=0._DP
                DO IATM = 1,NATMI
                  IAD = EIRENE_INDIRECT_ADDRESS(IPLS,IATM,NPLS)
                  PAPL(IPLS,J)=PAPL(IPLS,J)+PAPL(IAD,J)*FATM(IATM)
                END DO
              ELSE  
                PAPL(IPLS,J)=PAPL(IPLS,J)*FATM(0)
              END IF
            END IF
            IF (LPMPL) THEN
              IF (NLSPCSCL_MOL) THEN
                PMPL(IPLS,J)=0._DP
                DO IMOL = 1,NMOLI
                  IAD = EIRENE_INDIRECT_ADDRESS(IPLS,IMOL,NPLS)
                  PMPL(IPLS,J)=PMPL(IPLS,J)+PMPL(IAD,J)*FMOL(IMOL)
                END DO
              ELSE  
                PMPL(IPLS,J)=PMPL(IPLS,J)*FMOL(0)
              END IF
            END IF
            IF (LPIPL) THEN
              IF (NLSPCSCL_ION) THEN
                PIPL(IPLS,J)=0._DP
                DO IION = 1,NIONI
                  IAD = EIRENE_INDIRECT_ADDRESS(IPLS,IION,NPLS)
                  PIPL(IPLS,J)=PIPL(IPLS,J)+PIPL(IAD,J)*FION(IION)
                END DO
              ELSE  
                PIPL(IPLS,J)=PIPL(IPLS,J)*FION(0)
              END IF
            END IF
            IF (LPPHPL) THEN
              IF (NLSPCSCL_PHOT) THEN
                PPHPL(IPLS,J)=0._DP
                DO IPHOT = 1,NPHOTI
                  IAD = EIRENE_INDIRECT_ADDRESS(IPLS,IPHOT,NPLS)
                  PPHPL(IPLS,J)=PPHPL(IPLS,J)+PPHPL(IAD,J)*FPHOT(IPHOT)
                END DO
              ELSE  
                PPHPL(IPLS,J)=PPHPL(IPLS,J)*FPHOT(0)
              END IF
            END IF
            IF (LEAPL) EAPL(IPLS,J)=EAPL(IPLS,J)*FATM(0)
            IF (LEMPL) EMPL(IPLS,J)=EMPL(IPLS,J)*FMOL(0)
            IF (LEIPL) EIPL(IPLS,J)=EIPL(IPLS,J)*FION(0)
            IF (LEPHPL) EPHPL(IPLS,J)=EPHPL(IPLS,J)*FPHOT(0)
            IF (LMAPL) MAPL(IPLS,J)=MAPL(IPLS,J)*FATM(0)
            IF (LMMPL) MMPL(IPLS,J)=MMPL(IPLS,J)*FMOL(0)
            IF (LMIPL) MIPL(IPLS,J)=MIPL(IPLS,J)*FION(0)
            IF (LMPHPL) MPHPL(IPLS,J)=MPHPL(IPLS,J)*FPHOT(0)
          END DO
        END DO
        DO IPLS=1,NPLSI
          IF (NLSPCSCL_ATM) THEN
            PAPLI(IPLS,ISTRA)=0._DP
            DO IATM=1,NATMI
              IAD = EIRENE_INDIRECT_ADDRESS(IPLS,IATM,NPLS)
              PAPLI(IPLS,ISTRA)=PAPLI(IPLS,ISTRA)+
     .                          PAPLI(IAD,ISTRA)*FATM(IATM)
            END DO
          ELSE
            PAPLI(IPLS,ISTRA)=PAPLI(IPLS,ISTRA)*FATM(0)
          END IF
          IF (NLSPCSCL_MOL) THEN
            PMPLI(IPLS,ISTRA)=0._DP
            DO IMOL=1,NMOLI
              IAD = EIRENE_INDIRECT_ADDRESS(IPLS,IMOL,NPLS)
              PMPLI(IPLS,ISTRA)=PMPLI(IPLS,ISTRA)+
     .                          PMPLI(IAD,ISTRA)*FMOL(IMOL)
            END DO
          ELSE
            PMPLI(IPLS,ISTRA)=PMPLI(IPLS,ISTRA)*FMOL(0)
          END IF
          IF (NLSPCSCL_ION) THEN
            PIPLI(IPLS,ISTRA)=0._DP
            DO IION=1,NIONI
              IAD = EIRENE_INDIRECT_ADDRESS(IPLS,IION,NPLS)
              PIPLI(IPLS,ISTRA)=PIPLI(IPLS,ISTRA)+
     .                          PIPLI(IAD,ISTRA)*FION(IION)
            END DO
          ELSE
            PIPLI(IPLS,ISTRA)=PIPLI(IPLS,ISTRA)*FION(0)
          END IF
          IF (NLSPCSCL_PHOT) THEN
            PPHPLI(IPLS,ISTRA)=0._DP
            DO IPHOT=1,NPHOTI
              IAD = EIRENE_INDIRECT_ADDRESS(IPLS,IPHOT,NPLS)
              PPHPLI(IPLS,ISTRA)=PPHPLI(IPLS,ISTRA)+
     .                          PPHPLI(IAD,ISTRA)*FPHOT(IPHOT)
            END DO
          ELSE
            PPHPLI(IPLS,ISTRA)=PPHPLI(IPLS,ISTRA)*FPHOT(0)
          END IF
          EAPLI(IPLS,ISTRA)=EAPLI(IPLS,ISTRA)*FATM(0)
          EMPLI(IPLS,ISTRA)=EMPLI(IPLS,ISTRA)*FMOL(0)
          EIPLI(IPLS,ISTRA)=EIPLI(IPLS,ISTRA)*FION(0)
          EPHPLI(IPLS,ISTRA)=EPHPLI(IPLS,ISTRA)*FPHOT(0)
          MAPLI(IPLS,ISTRA)=MAPLI(IPLS,ISTRA)*FATM(0)
          MMPLI(IPLS,ISTRA)=MMPLI(IPLS,ISTRA)*FMOL(0)
          MIPLI(IPLS,ISTRA)=MIPLI(IPLS,ISTRA)*FION(0)
          MPHPLI(IPLS,ISTRA)=MPHPLI(IPLS,ISTRA)*FPHOT(0)
        END DO

        PAPLI(0,ISTRA)=SUM(PAPLI(1:NPLSI,ISTRA))
        PMPLI(0,ISTRA)=SUM(PMPLI(1:NPLSI,ISTRA))
        PIPLI(0,ISTRA)=SUM(PIPLI(1:NPLSI,ISTRA))
        PPHPLI(0,ISTRA)=SUM(PPHPLI(1:NPLSI,ISTRA))
        EAPLI(0,ISTRA)=SUM(EAPLI(1:NPLSI,ISTRA))
        EMPLI(0,ISTRA)=SUM(EMPLI(1:NPLSI,ISTRA))
        EIPLI(0,ISTRA)=SUM(EIPLI(1:NPLSI,ISTRA))
        EPHPLI(0,ISTRA)=SUM(EPHPLI(1:NPLSI,ISTRA))
        MAPLI(0,ISTRA)=SUM(MAPLI(1:NPLSI,ISTRA))
        MMPLI(0,ISTRA)=SUM(MMPLI(1:NPLSI,ISTRA))
        MIPLI(0,ISTRA)=SUM(MIPLI(1:NPLSI,ISTRA))
        MPHPLI(0,ISTRA)=SUM(MPHPLI(1:NPLSI,ISTRA))
C
C  ELECTRON TALLIES
C
        DO 551 J=1,NSBOX_TAL
          IF (LPAEL) PAEL(J)=PAEL(J)*FATM(0)
          IF (LPMEL) PMEL(J)=PMEL(J)*FMOL(0)
          IF (LPIEL) PIEL(J)=PIEL(J)*FION(0)
          IF (LPPHEL) PPHEL(J)=PPHEL(J)*FPHOT(0)
  551   CONTINUE
        PAELI(ISTRA)=PAELI(ISTRA)*FATM(0)
        PMELI(ISTRA)=PMELI(ISTRA)*FMOL(0)
        PIELI(ISTRA)=PIELI(ISTRA)*FION(0)
        PPHELI(ISTRA)=PPHELI(ISTRA)*FPHOT(0)
        DO 552 J=1,NSBOX_TAL
          IF (LEAEL) EAEL(J)=EAEL(J)*FATM(0)
          IF (LEMEL) EMEL(J)=EMEL(J)*FMOL(0)
          IF (LEIEL) EIEL(J)=EIEL(J)*FION(0)
          IF (LEPHEL) EPHEL(J)=EPHEL(J)*FPHOT(0)
  552   CONTINUE
        EAELI(ISTRA)=EAELI(ISTRA)*FATM(0)
        EMELI(ISTRA)=EMELI(ISTRA)*FMOL(0)
        EIELI(ISTRA)=EIELI(ISTRA)*FION(0)
        EPHELI(ISTRA)=EPHELI(ISTRA)*FPHOT(0)
C
C  SPECTRUM TALLIES
C
        DO ISPC=1,NADSPC
          SELECT CASE (ESTIML(ISPC)%IPRTYP)
          CASE (0)
            FADD = FPHOT(ESTIML(ISPC)%IPRSP)
          CASE (1)
            FADD = FATM(ESTIML(ISPC)%IPRSP)
          CASE (2)
            FADD = FMOL(ESTIML(ISPC)%IPRSP)
          CASE (3)
            FADD = FION(ESTIML(ISPC)%IPRSP)
          CASE DEFAULT
            FADD = 1._DP
          END SELECT
          ESTIML(ISPC)%SPC = ESTIML(ISPC)%SPC * FADD
          ESTIML(ISPC)%SPCS = ESTIML(ISPC)%SPCS * FADD
        END DO
C
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) 'RESCALING OF TRACKLENGTH TALLIES COMPLETED'
        WRITE (iunout,*) 'RESCALING FACTORS:'
        CALL EIRENE_MASR4 ('FATM,FMOL,FION,FPHOT            ',
     .               FATM(0), FMOL(0), FION(0), FPHOT(0))
        IF (NLSPCSCL) THEN
          IF (NLSPCSCL_ATM)
     .       CALL EIRENE_MASRR1 ('FATM       ', FATM(1:NATMI), NATMI, 6)
          IF (NLSPCSCL_MOL)
     .       CALL EIRENE_MASRR1 ('FMOL       ', FMOL(1:NMOLI), NMOLI, 6)
          IF (NLSPCSCL_ION)
     .       CALL EIRENE_MASRR1 ('FION       ', FION(1:NIONI), NIONI, 6)
          IF (NLSPCSCL_PHOT)
     .       CALL EIRENE_MASRR1 ('FPHOT      ',FPHOT(1:NPHOTI),NPHOTI,6)
        END IF
        CALL EIRENE_LEER(2)
C
      ENDIF

      DEALLOCATE(FATM,FMOL,FION,FPHOT)

      RETURN
      END SUBROUTINE EIRENE_SCALE_TALLIES
