cdr These are the SHORT versions of WRPLAM.F, RPLAM.F
cdr They write and read (RPLAM_SHRT, below) background data onto/from fort.13
cdr Distinct from WRPLAM_long here only the virtual background tallies are written/read,
cdr (tallies T, n, V for ipls=npls_fix+1,nplsi),
cdr and, new, march 2020 (optional) only some related atomic data.

cdr Currently still not the primary source sampling information.
cdr That would be needed
cdr if also some primary source data would be iterated via a MODUSR routine.
cdr
c
cdr  jan 16: started to cleanup, comments
cdr          remove redundant parameter iflg
cdr  jan 18: comments
cdr  2019,2020, more comments
cdr  feb 2020: revived within ISFN activities: retain selected A&M rates
cdr            involving the NPLS_FIX+1:NPLS virtual field species, lgvac, etc...

c
C
      SUBROUTINE EIRENE_WRPLAM_SHRT(TRCFLE,CALLEDFROM)

cdr Only the input tallies of the last (virtual) plasma species: npls_fix+1,...nplsi
cdr are written/read using I/O stream fort.13. 
cdr The other ones 1,...,npls_fix are directly transferred 
cdr from external plasma code/ external data set, or re-computed

cdr Only cross-sections, collision rates, other parameters,
cdr involving some virt. species need to be retained.

c  Here: NLSRT13=true : wrplam_short and rplam_short are called from WRPLAM, RPLAM, resp.
c

cdr this is the SHORT version of WRPLAM.F
cdr It writes and reads (RPLAM_SHRT) background data onto/from fort.13
cdr Distinct from WRPLAM_long here only the background tallies are written/read,
cdr (tallies T, n, V for ipls=1,nplsi), but not the atomic data,
cdr nor the primary source sampling information.
cdr for BGK type nonlinear iterations, with velocity-independent rates,
cdr this may be sufficient.
cdr Better: add here also the rates and other atomic data needed to streamline
cdr         nonlinear iterations.



      USE EIRMOD_PRECISION, ONLY: DP
      USE EIRMOD_PARMMOD, ONLY: IFOFF, NRAD, NPLS
      USE EIRMOD_CINIT, ONLY: FORT
      USE EIRMOD_COMUSR, ONLY: NPLSI, TIIN, DIIN, VXIN, VYIN, VZIN,
     .                         NPLS_FIX
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_COMUSR, ONLY: LGVAC
      USE EIRMOD_COMNNL, ONLY: NFLA_VIRT, NREA_VIRT,
     .                         NFLA_IPLS, NFLA_ISWR, NFLA_IR
      USE EIRMOD_COMXS, ONLY: TABEI1, TABCX3, TABPI3, TABEL3,
     .                        EELEI1, EPLCX3, EELPI3, EPLEL3

      IMPLICIT NONE
      LOGICAL,INTENT(IN) :: TRCFLE
      CHARACTER(*) CALLEDFROM
      INTEGER :: IREA, ISWR, IREI, IRCX, IRPI, IREL, IFL, IPLS, ND, IO
      REAL(DP), ALLOCATABLE, SAVE :: TIAR(:,:), DIAR(:,:),
     R                         VXAR(:,:), VYAR(:,:), VZAR(:,:)

      write (iunout,*) 'WRPLAM_SHRT called from ',calledfrom

cdr  We need to verify somehow, that NLSHRT13 is the same still as it
cdr  was when writing FORT.13

      OPEN (UNIT=13+ifoff,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 13+ifoff
#ifdef CHECKBIN
      OPEN (UNIT=113,ACCESS='SEQUENTIAL',FORM='FORMATTED')
      REWIND 113
#endif
cdr  only write plasma background data for species, which are not already
cdr  transferred via common BRAEIR, i.e. only: npls_fix+1,....nplsi
cdr  I.e. the virtual background species for nonlinear iterations
cdr  have to come last in the list of all background species.
      IF(NPLS_FIX.LT.NPLSI) THEN
        ND = NPLSI - NPLS_FIX
        IF (.NOT.ALLOCATED(TIAR)) THEN
          ALLOCATE (TIAR(ND,NRAD))
          ALLOCATE (DIAR(ND,NRAD))
          ALLOCATE (VXAR(ND,NRAD))
          ALLOCATE (VYAR(ND,NRAD))
          ALLOCATE (VZAR(ND,NRAD))
        END IF
        TIAR(1:ND,1:NRAD) = TIIN(NPLS_FIX+1:NPLSI,1:NRAD)
        DIAR(1:ND,1:NRAD) = DIIN(NPLS_FIX+1:NPLSI,1:NRAD)
        VXAR(1:ND,1:NRAD) = VXIN(NPLS_FIX+1:NPLSI,1:NRAD)
        VYAR(1:ND,1:NRAD) = VYIN(NPLS_FIX+1:NPLSI,1:NRAD)
        VZAR(1:ND,1:NRAD) = VZIN(NPLS_FIX+1:NPLSI,1:NRAD)
        WRITE (13+ifoff,IOSTAT=IO)
     w           TIAR, DIAR, VXAR, VYAR, VZAR
#ifdef CHECKBIN
        WRITE (113,*)
     w           TIAR, DIAR,  VXAR, VYAR, VZAR
#endif

        IF (TRCFLE) THEN
          CALL EIRENE_LEER(1)
          WRITE (iunout,*)
     w   'WRPLAM_SHRT: VIRT. BACKGROUND IS WRITTEN ON '//FORT//'13'
          WRITE (iunout,*) 'SPECIES NPLS_FIX+1, NPLSI ',NPLS_FIX+1,NPLSI
          CALL EIRENE_LEER(1)
        ENDIF
      END IF

cdr  next: write atomic/molecular data set in MODUSR, for collision processes
cdr        involving the npls_fix+1:npls virtual background species
cdr        and their possible dependencies on parameters distinct from PLSTLS.
cdr        There are NFLA_VIRT such data arrays.
cdr  All other A&M data will be recomputed in SETAMD.f

c  nfla_virt is not a pointer to some target
      IF(NFLA_VIRT == 0) THEN
        WRITE(iunout,'(4a)')
     w         "WRPLAM_SHRT: NFLA_VIRT = 0. ",
     w         "NO COLLISION DATA WILL BE STORED IN ", FORT, "13"

        CLOSE (UNIT=13+ifoff)
        RETURN
      END IF

cdr  now we have: NFLA_VIRT
      write (13+ifoff) nfla_virt, nrea_virt
      write (13+ifoff) nfla_ipls(1:nfla_virt), nfla_iswr(1:nrea_virt),
     w                 nfla_ir(1:nrea_virt)
#ifdef CHECKBIN
      write (113,*) nfla_virt, nrea_virt
      write (113,*) nfla_ipls(1:nfla_virt), nfla_iswr(1:nrea_virt),
     w              nfla_ir(1:nrea_virt)
#endif

      DO IFL=1,NFLA_VIRT
        ipls=nfla_ipls(ifl)
        if (ipls.ge.npls_fix+1. .and. ipls.le.npls) then
c  data for virt. species ipls
          write (13+ifoff)
     w      lgvac(1:nrad,ipls)
#ifdef CHECKBIN
          write (113,*) lgvac(1:nrad,ipls)
#endif
        endif
      END DO

cdr
c   Here we should add a safety:
c   The rates stored below should refer to virtual bulk species reactants
c   ipls= npls_fix+1,...,npls
c   In some cases, however, these virtual bulk species may be an intermediate
c   reactant only, i.e.: not stored.
c   Then only the final bulk reactant in a multi-step
c   process is retained. Example: photon excitation driven
c   additional (EI) electron impact ionisation rates for atoms.
c   The photon absorption rates, tabph1,..., however, relate to
c   the absorber, which must be a virtual background atom species.
c
cdr
      do irea = 1, nrea_virt
        iswr=nfla_iswr(irea)
        select case (iswr)
c  data for reaction irei, ircx, irpi or irel
        case (1)
          irei=nfla_ir(irea)
          write (13+ifoff)
     w      tabei1(irei,1:nrad),eelei1(irei,1:nrad)
#ifdef CHECKBIN
          write (113,*) tabei1(irei,1:nrad),eelei1(irei,1:nrad)
#endif
        case (3)
          ircx=nfla_ir(irea)
          write (13+ifoff)
     w      tabcx3(ircx,1:nrad,1:9),eplcx3(ircx,1:nrad,1:9)
#ifdef CHECKBIN
          write(113,*) tabcx3(ircx,1:nrad,1:9),eplcx3(ircx,1:nrad,1:9)
#endif
        case (4)
          irpi=nfla_ir(irea)
          write (13+ifoff)
     w      tabpi3(irei,1:nrad,1:9),eelpi3(irpi,1:nrad,1:9)
#ifdef CHECKBIN
          write(113,*) tabpi3(irei,1:nrad,1:9),eelpi3(irpi,1:nrad,1:9)
#endif
        case (5)
          irel=nfla_ir(irea)
          write (13+ifoff)
     w      tabel3(irel,1:nrad,1:9),eplel3(irel,1:nrad,1:9)
#ifdef CHECKBIN
          write(113,*) tabel3(irel,1:nrad,1:9),eplel3(irel,1:nrad,1:9)
#endif
cdr  case 6,7: to be done
cdr  for photons we may also need "irrc" data
        end select
      ENDDO

cdr  some diagnostics, if trcfle.....

      CLOSE (UNIT=13+ifoff)
#ifdef CHECKBIN
      CLOSE (UNIT=113)
#endif

      RETURN
      END SUBROUTINE EIRENE_WRPLAM_SHRT
C .......................................................................

      SUBROUTINE EIRENE_RPLAM_SHRT(TRCFLE,CALLEDFROM)
      USE EIRMOD_PRECISION, ONLY: DP
      USE EIRMOD_PARMMOD, ONLY: IFOFF, NRAD, NPLS
      USE EIRMOD_CINIT, ONLY: FORT, FORT_LC
      USE EIRMOD_COMUSR, ONLY: NPLSI, TIIN, DIIN, VXIN, VYIN, VZIN,
     .                         NPLS_FIX
!pb   USE EIRMOD_CCOUPL, ONLY: NFLA
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_COMNNL, ONLY: NFLA_VIRT, NREA_VIRT,
     .                         NFLA_IPLS, NFLA_ISWR,
     .                         NFLA_IR, LG_STORE, TAB_STORE, E_STORE

      IMPLICIT NONE

      LOGICAL,INTENT(IN) :: TRCFLE
      CHARACTER(*) CALLEDFROM
      CHARACTER*8 FILENAME
      CHARACTER*3 FILENUMBER
      INTEGER :: IO
      INTEGER :: IREA, ISWR, IREI, IRCX, IRPI, IREL, IFL, IPLS
      LOGICAL EX

C ........................................................................

      if (13+ifoff.ge.100) then
      WRITE(FILENUMBER,'(I3)') 13+ifoff
      else
        WRITE(FILENUMBER,'(I2)') 13+ifoff
      end if
      FILENAME = FORT_LC//trim(FILENUMBER)
      INQUIRE(file=trim(FILENAME),exist=EX)
      IF (EX) THEN
#ifdef CRAY
        OPEN (UNIT=13+ifoff,FORM='UNFORMATTED',IOSTAT=IO)
#else
        OPEN (UNIT=13+ifoff,ACCESS='SEQUENTIAL',FORM='UNFORMATTED',
     o        STATUS='OLD',IOSTAT=IO)
#endif
        write (iunout,*) 'RPLAM_SHRT called from ',calledfrom
        IF(IO.NE.0) THEN
          CALL EIRENE_LEER(1)
          WRITE(iunout,*)
     w     'ERROR IN RPLAM_SHRT: CANNOT READ ', FORT, '13'
          CALL EIRENE_LEER(1)
          CLOSE (UNIT=13+ifoff)
          GOTO 200
        END IF
      ELSE
        GOTO 200
      END IF

      REWIND 13+ifoff

      IF(NPLS_FIX.LT.NPLSI) THEN
cdr  only read plasma background data for species, which are not already
cdr  transferred via common BRAEIR, i.e. only: npls_fix+1,....nplsi
cdr  I.e. the virtual background species for nonlinear iterations
cdr  have to come last in the list of all background species.  
        REWIND 13+ifoff
        READ (13+ifoff,IOSTAT=IO)
     R       TIIN(NPLS_FIX+1:NPLSI,1:NRAD),
     R       DIIN(NPLS_FIX+1:NPLSI,1:NRAD),
     R       VXIN(NPLS_FIX+1:NPLSI,1:NRAD),
     R       VYIN(NPLS_FIX+1:NPLSI,1:NRAD),
     R       VZIN(NPLS_FIX+1:NPLSI,1:NRAD)
        IF(IO.NE.0) GOTO 200
        IF (TRCFLE) THEN
          CALL EIRENE_LEER(1)
          WRITE (iunout,*)
     w   'RPLAM_SHRT: VIRT. BACKGROUND IS READ FROM ', FORT, '13'
          WRITE (iunout,*) 'SPECIES NPLS_FIX+1, NPLSI ',NPLS_FIX+1,NPLSI
          CALL EIRENE_LEER(1)
        ENDIF
      END IF

      NFLA_VIRT = 0
      NREA_VIRT = 0
      READ (13+IFOFF,IOSTAT=IO) NFLA_VIRT, NREA_VIRT
      IF(IO /= 0) THEN
        WRITE(iunout,'(4a)')
     w         "RPLAM_SHRT: NFLA_VIRT COULD NOT BE READ. ",
     w         "NO COLLISION DATA WILL BE READ FROM ", FORT, "13"
        CALL EIRENE_LEER(1)
        CLOSE (UNIT=13+ifoff)
        RETURN
      END IF

cdr  now we have: NFLA_VIRT
cdr  next: read the NFLA_VIRT collisional process data related to these virtual species
      IF (NFLA_VIRT == 0) THEN
        WRITE(iunout,'(4a)')
     w         "RPLAM_SHRT: NFLA_VIRT = 0. ",
     w         "NO COLLISION DATA WILL BE READ FROM ", FORT, "13"
        CALL EIRENE_LEER(1)
        CLOSE (UNIT=13+ifoff)
      END IF
      IF (NREA_VIRT == 0) THEN
        WRITE(iunout,'(4a)')
     w         "RPLAM_SHRT: NREA_VIRT = 0. ",
     w         "NO COLLISION DATA WILL BE READ FROM ", FORT, "13"
        CALL EIRENE_LEER(1)
        CLOSE (UNIT=13+ifoff)
      END IF

      IF (ALLOCATED(NFLA_ISWR) .AND. (SIZE(NFLA_ISWR) < NREA_VIRT)) THEN
        DEALLOCATE (NFLA_ISWR)
        DEALLOCATE (NFLA_IR)
        DEALLOCATE (NFLA_IPLS)
      END IF

      IF (.NOT.ALLOCATED(NFLA_ISWR)) THEN
#ifdef F2003
        ALLOCATE (NFLA_ISWR(NREA_VIRT), SOURCE = 0)
        ALLOCATE (NFLA_IR(NREA_VIRT), SOURCE = 0)
        ALLOCATE (NFLA_IPLS(NPLS), SOURCE = 0)
#else
        ALLOCATE (NFLA_ISWR(NREA_VIRT))
        ALLOCATE (NFLA_IR(NREA_VIRT))
        ALLOCATE (NFLA_IPLS(NPLS))
        NFLA_ISWR = 0
        NFLA_IR = 0
        NFLA_IPLS = 0
#endif
      END IF

      IF (ALLOCATED(TAB_STORE) .AND. SIZE(TAB_STORE,1) < NREA_VIRT) THEN
        DEALLOCATE (TAB_STORE)
        DEALLOCATE (E_STORE)
      END IF

      IF (ALLOCATED(LG_STORE) .AND. SIZE(LG_STORE,2) < NFLA_VIRT) THEN
        DEALLOCATE(LG_STORE)
      END IF

      IF (.NOT.ALLOCATED(TAB_STORE)) THEN
#ifdef F2003
        ALLOCATE (TAB_STORE(NREA_VIRT,NRAD,9), SOURCE=0._DP)
        ALLOCATE (E_STORE(NREA_VIRT,NRAD,9), SOURCE=0._DP)
#else
        ALLOCATE (TAB_STORE(NREA_VIRT,NRAD,9))
        ALLOCATE (E_STORE(NREA_VIRT,NRAD,9))
        TAB_STORE = 0._DP
        E_STORE = 0._DP
#endif
      END IF

      IF (.NOT.ALLOCATED(LG_STORE)) THEN
#ifdef F2003
        ALLOCATE (LG_STORE(NRAD,NFLA_VIRT), SOURCE=.FALSE.)
#else
        ALLOCATE (LG_STORE(NRAD,NFLA_VIRT))
        LG_STORE = .FALSE.
#endif
      END IF

      READ (13+IFOFF) NFLA_IPLS(1:NFLA_VIRT), NFLA_ISWR(1:NREA_VIRT),
     R                NFLA_IR(1:NREA_VIRT)

      DO IFL = 1, NFLA_VIRT
        ipls=nfla_ipls(ifl)
        if (ipls.ge.npls_fix+1 .and. ipls.le.npls) then
c  data for virt. species ipls
          read (13+ifoff) lg_store(1:nrad,ifl)
        endif
      END DO

      do irea = 1, nrea_virt
        iswr=nfla_iswr(irea)
        select case (iswr)
        case (1)
cdr restore some EI rates for iterative procedure (nonlinearities)
          irei=nfla_ir(irea)
          read (13+ifoff)
     w      tab_store(irea,1:nrad,1),e_store(irea,1:nrad,1)
        case (3)
          ircx=nfla_ir(irea)
          read (13+ifoff)
     w      tab_store(irea,1:nrad,1:9),e_store(irea,1:nrad,1:9)
        case (4)
          irpi=nfla_ir(irea)
          read (13+ifoff)
     w      tab_store(irei,1:nrad,1:9),e_store(irea,1:nrad,1:9)
        case (5)
cdr restore some EL rates for iterative procedure (nonlinearities)
          irel=nfla_ir(irea)
          read (13+ifoff)
     w      tab_store(irea,1:nrad,1:9),e_store(irea,1:nrad,1:9)
cdr  case 6,7: to be done
        end select
      ENDDO

      CLOSE (UNIT=13+ifoff)
      RETURN

  200 CONTINUE

      WRITE(iunout,'(4a)')
     w 'RPLAM_SHRT: CANNOT READ ', FORT, '13: ',
     w 'ZERO VIRTUAL BACKGROUND IS ASSIGNED'
      WRITE (iunout,*) 'SPECIES NPLS_FIX+1, NPLSI ',NPLS_FIX+1,NPLSI
      CALL EIRENE_LEER(1)

      TIIN(NPLS_FIX+1:NPLSI,1:NRAD)=0._DP
      DIIN(NPLS_FIX+1:NPLSI,1:NRAD)=0._DP
      VXIN(NPLS_FIX+1:NPLSI,1:NRAD)=0._DP
      VYIN(NPLS_FIX+1:NPLSI,1:NRAD)=0._DP
      VZIN(NPLS_FIX+1:NPLSI,1:NRAD)=0._DP

      IF (EX) CLOSE (UNIT=13+ifoff)

      RETURN

      END SUBROUTINE EIRENE_RPLAM_SHRT
