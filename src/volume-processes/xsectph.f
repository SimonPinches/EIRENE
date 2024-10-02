Cdr  analogue to  xsecta, xsectm,.... printout to be sync.  
cdr  And: call xstph for non-default photonic reaction options. 
cdr  Unfinished code, for photons.
c
Cdr  in xsectp, there only call to xstrc.f.
cdr  At present: it seems to be just the other way round.
C
      SUBROUTINE EIRENE_XSECTPH
C
C  TABLE FOR CROSS-SECTION AND REACTION RATES FOR PHOTONS
C
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMXS
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CTRCEI
      USE EIRMOD_PHOTON
      IMPLICIT NONE
c
cdr   PHOTON COLLISIONS, PH - type (separate from OT processes, which are
cdr                                 of H.11, H.12 type, popul. ratios)
c
      integer :: kk,iphot,idsc,nrc,ipl0,ipl1,ipl2,ityp1,ityp2,ifnd,
     .           updf,mode,idph

      IDPH=0

      DO IPHOT=1,NPHOTI
        IDSC=0
        PHV_LGPHOT(IPHOT,0,0)=0
        PHV_LGPHOT(IPHOT,0,1)=0
C
C   AT PRESENT NO DEFAULT MODEL
C
        IF (NRCPH(IPHOT).EQ.0) THEN
          PHV_NPHOTI(IPHOT)=0
C
C  NON-DEFAULT "PH" MODEL:
C
        ELSEIF(NRCPH(IPHOT) > 0) THEN
          DO NRC=1,NRCPH(IPHOT)
            KK=IREACPH(IPHOT,NRC)
            IF (ISWR(KK).NE.7) CYCLE
            IDSC=IDSC+1
            IDPH=IDPH+1
            NREAPH(IDPH) = KK
            CALL EIRENE_PH_XSECTPH (IPHOT,NRC,IDSC)
CDR  HERE SHOULD BE CALL TO XSTPH, GENERAL ROUTINE FOR PH PROCESSES
          ENDDO
          PHV_NPHOTI(IPHOT)=IDSC
C  NO "PH" MODEL DEFINED
        ELSE
          PHV_NPHOTI(IPHOT)=0
        ENDIF

CDR     PHV_NPHOTIM(IPHOT)=PHV_NPHOTI(IPHOT)-1

        PHV_LGPHOT(IPHOT,0,0)=IDSC

      ENDDO
csw
csw output:
csw
      DO IPHOT=1,NPHOTI
C
        IF (TRCAMD) THEN
          CALL EIRENE_MASBOX ('PHOTON SPECIES IPHOT = '//TEXTS(IPHOT))
          CALL EIRENE_LEER(1)
C
          IF(PHV_NPHOTI(iphot).eq.0) then
            CALL EIRENE_LEER(1)
            WRITE (iunout,*) 'NO "PH"-REACTION FOR THIS PHOTON'
            CALL EIRENE_LEER(1)
          ELSE
            DO IDSC=1,PHV_NPHOTI(IPHOT)
              CALL EIRENE_LEER(1)
              WRITE (iunout,*) '(OTHER) REACTION NO. IRPH= ',IDSC
              CALL EIRENE_LEER(1)

                  ipl0=PHV_LGPHOT(iphot,idsc,1)
                  ifnd=PHV_LGPHOT(iphot,idsc,2)
                  kk=PHV_LGPHOT(iphot,idsc,3)
                  updf=PHV_LGPHOT(iphot,idsc,4)
                  mode=PHV_LGPHOT(iphot,idsc,5)

                  ityp1=PHV_N1STOTph(iphot,idsc,1)
                  ipl1= PHV_N1STOTph(iphot,idsc,2)
                  ityp2=PHV_N2NDOTph(iphot,idsc,1)
                  ipl2= PHV_N2NDOTph(iphot,idsc,2)

                  write (iunout,*) 'irph,ipl0,il,kk,updf,mode'
                  write (iunout,*)  idsc,ipl0,ifnd,kk,updf,mode
                  write (iunout,*) 'ityp1,ipl1,ityp2,ipl2'
                  write (iunout,*)  ityp1,ipl1,ityp2,ipl2
                  call EIRENE_leer(1)
               enddo
            endif
         endif
      enddo

      RETURN
      END SUBROUTINE EIRENE_XSECTPH
