cdr  feb, 16., 2015, added: naint=22, modcol=2 option, EB=1.5 Ti
cdr  aug,  4., 2015, added: naint=24, modcol=2 option, EB=1.5 Ti
cdr  aug,  4., 2015, added: naint=26, modcol=2 option, EB=1.5 Ti
cdr  nov.      2015: noted: modcol=3: take sigma(E) * sqrt(E), to be done
cdr  aug.      2016: set e0 low energy cut-off, as on fpath routines, for H.3 rates
cdr                  also: lgvac(i,ipl), lgvac(i,npls+1) is used, not finished.
!pb  apr       2016: eelds -> eelei
!pb  may       2016: tabds1 -> tabei1
cdr  Nov.      2016: final ds --> ei notational unifications
CDR  July      2017: RC reactions connected. trcamd in parameter list
c                    for function eirene_sngl_poly
cdr  Sept      2018  modcol(...,4,...), (energy-weighted rates) (rather than (..,3,..)
cdr                  naint=21 and =29: done
cdr  Feb       2020: sync code for EI, PI, CX, and EL processes. add energy weighted rates
cdr                  for modcol=1, na = 23, 25, 27
cdr                  e.g. also now for EL processes (because of bgk balances)
cdr  Feb       2021  ND2 for sngl_poly. Not used yet. ND2=9 so far.
cdr  Feb       2022  ADIN reset to zero only for tallies IAIN with NA=20,...29,
cdr                  but not for other ADIN tallies.
cdr                  Otherwise other ADIN tallies, e.g. from INFCOP, are lost here.

CDR:  A&M Data diagnostics routine, added in Jan. 2014
C  PUT SELECTED EIRENE ATOMIC DATA FIELDS ONTO ADIN ARRAY FOR OUTPUT.
C  ADIN CONTAINS RATE COEFFICIENTS (VOL/TIME) IN ATOMIC UNITS
c
c  modcol=1: rate coefficients only dependent on local background data, not on test particle parameters
c            tabcx3(...,1),tabel3(...1),tabpi3(...1),tabei1(...)
c            are rates, density of impacting bulk ion included,
c            so here we divide again by ne or ni.
c
c  modcol=2: rate coefficients depend on Eb = energy of impacting test particle
c            tabcx3(...,1:nend),tabel3(...,1:nend),tabpi3(...,1:nend)
c               are ln(rate)
c               with ln(rate)= sum_i=1^nend  ln^i(Eb) tab..3(...,i)
c            ADIN is evaluated with Eb = 3/2 T, T =T(ipls)
c
c  modcol=3: rate coefficients depend on Eb = energy of impacting test particle
c            tabcx3(...,1:nend),tabel3(...,1:nend),tabpi3(...,1:nend)
c               are ln(rate)
c               with ln(rate)= sum_i=1^nend  ln^i(Eb) tab..3(...,i)
c            ADIN is evaluated with Eb = ????, T =T(ipls)
c
c
c
C
C  ATOMIC UNITS FOR REACTION RATE COEFFICIENTS: A0^2 V0 = 0.612E-08 CM^3-S
C  TO CONVERT THE ADDITIONAL TALLIES ADIN INTO UNITS OF CM**3/S, MULTIPLY ADIN BY 0.612 e-08
c
c  done for naint=20,22,24,26 and modcol=1

c  naint=20:   Tabei1(irei,....) electron impact collision rate, 1/s --> cm^3/s      ! done
c  naint=21:   eelei1(irei,....) electron cooling rate           eV/s --> cm^3 eV/s  ! done

c  naint=22:   Tabcx3(ircx,..,1) charge exchange collision rate, 1/s --> cm^3/s      ! done
c  naint=23:   eplcx3(ircx,..,1) cx        energy-weighted rate, eV/s --> cm^3 eV/s  ! not ready

c  naint=24:   Tabel3(irel,..,1) elastic         collision rate, 1/s --> cm^3/s      ! done
c  naint=25:   eplel3(irel,..,1) elastic   energy-weighted rate, eV/s --> cm^3 eV/s  ! not ready

c  naint=26:   Tabpi3(irpi,..,1) heavy particle imp.  coll.rate, 1/s --> cm^3/s      ! done
c  naint=27:   eplpi3(irpi,..,1) ditto,    energy-weighted rate, eV/s --> cm^3 eV/s  ! not ready

c  naint=28:   Tabrc1(irrc,....) electron-ion volume recomb.rate, 1/s --> cm^3/s     ! done
c  naint=29:   eelrc1(irrc,....) ditto,    energy-weighted rate, eV/s --> cm^3 eV/s  ! done

c

      SUBROUTINE EIRENE_AMDIAG
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CTRCEI, ONLY: TRCAMD
      USE EIRMOD_CCONA
      USE EIRMOD_CGRID
      USE EIRMOD_COMXS
      USE EIRMOD_CTEXT

      IMPLICIT NONE

      REAL(DP) :: AU, ELB, FP(6), RCMIN, RCMAX,
     .            TBCX3(9), TBPI3(9), TBEL3(9),
     .            EIRENE_SNGL_POLY, EARRH,
     .            RMASSS,EBFAC,RATE,
     .            TII, TBCX, TBEI, TBPI, TBEL, TBRC, 
     .            ELEI, EPPI, EPCX, EPEL, ELRC, 
cdr  functions for 'on the fly' evaluation of A&M data
     .          EIRENE_FEELEI1, EIRENE_FEELPI3,
     .          EIRENE_FEELRC1,
     .          EIRENE_FEPLCX3, EIRENE_FEPLEL3,
     .          EIRENE_FTABCX3, EIRENE_FTABPI3,
     .          EIRENE_FTABEI1, EIRENE_FTABRC1,
     .          EIRENE_FTABEL3,
     .          EIRENE_RATE_COEFF
      REAL(DP),PARAMETER :: EMIN = 0.1003_DP  ! hard coded cut-off for EBEAM parameter in H.3 fits
      REAL(DP),PARAMETER :: EMINL=-2.3000_DP  ! hard coded cut-off for EBEAM parameter in H.3 fits
      INTEGER :: NS,NA,IAIN,MM,KK,ND2,
     .           irei,ircx,irpi,irel,irrc,
     .           iat,iml,iio,ipl,isp,iplti,
     .           icell,iapi,impi,iipi,iacx,imcx,iicx,
     .           iael,imel,iiel,iaei,imei,iiei,iprc
      CHARACTER(4) :: CNO, CN1
      LOGICAL :: LEXP



      AU=0.6120D-08

      IF (.NOT.LADIN) THEN
        WRITE (IUNOUT,*) ' INPUT TALLY ADIN NOT AVAILABLE',
     .                   ' FOR STORING OF RATES IN AMDIAG '
        CALL EIRENE_LEER(1)
        RETURN
      END IF

      DO 190 IAIN=1,NAINI
        NS=NAINS(IAIN)    !  ns stands for ircx,irei,irel,irpi,irpi,..., internal number of process
        NA=NAINT(IAIN)    !  na stands for tally:  TAB..3(...),  EPL..3(...)

        IF ((NA < 20) .OR. (NA > 29)) CYCLE
! reinitialize ADIN to 0 in order to avoid residual values from prior
! iterations, e.g. in case of changed LGVAC
        if (any(adin(iain,:).ne.0.0)) then
          write (iunout,*) 'AMDIAG, Reset ADIN tally no. ',iain
c         call eirene_exit_own(1)
        endif
        ADIN(IAIN,:) = 0._DP

        IF (NSTORDR < NRAD) THEN
          WRITE (IUNOUT,*) 'AMDIAG NOT READY FOR STORAGE SAVING MODE'
          WRITE (IUNOUT,*) 'NS, NA ',NS,NA
          CALL EIRENE_LEER(1)
          CYCLE
        ENDIF

        MM = 0
        KK = 0

c  currently:  only tabei1, tabcx3, tabel3 and tabpi3 are available,
c              and only for modcol(..,2,ns)=1,2, also: some modcol(..,4,ns) started
c              modcol=1: rates depend only on background parameters,
c                        not on test particle parameters
c              modcol=2: rates depend also on test particle energy.
c                        Use E_test=1.5 kT_background
c  to be done:  (e.g. for Beams)
c              modcol=3: use sigma(E_test) * sqrt(E_test), ignore thermal background parameters
c
c.....................................................................

C  ELECTRON IMPACT RATE COEFFICIENT NO. IREI
        IF (NA.EQ.20.OR.NA.EQ.21) THEN
          irei=ns
          kk=NREAEI(irei)
c find collision partners corresponding to process irei: ISP
          IPL=0  ! ELECTRONS
c  first: try atoms
          LATEI: do iat=1,natmi
          do iaei=1,NAEII(iat)
            if (IREI.eq.LGAEI(IAT,IAEI)) then
              ISP=NSPH+IAT
              RMASSS=RMASSA(IAT)
              GOTO 170
            endif
          enddo
          enddo LATEI
c  irei is a process for atom iat, colliding with electron

c  next: try molecules
          LMLEI: do iml=1,nmoli
          do imei=1,NMEII(iml)
            if (IREI.eq.LGMEI(IML,IMEI)) then
              ISP=NSPA+IML
              RMASSS=RMASSM(IML)
              goto 170
            endif
          enddo
          enddo LMLEI
c  irei is a process for molecule iml, colliding with electron

c  next: try test ions
          LIOEI: do iio=1,nioni
          do iiei=1,NIEII(iio)
            if (IREI.eq.LGIEI(IIO,IIEI)) then
              ISP=NSPAM+IIO
              RMASSS=RMASSI(IIO)
              goto 170
            endif
          enddo
          enddo LIOEI
c  irei is a process for test ion iio, colliding with electron
c
c  no interacting particle species found
          TXTPLS(IAIN,NTALN) =
     .      'ELECTRON IMPACT REACTION RATE COEFFICIENT IREI ='//CNO
     .      //' KK='//CN1
          TXTPSP(IAIN,NTALN) = 'unidentified species    '
          TXTPUN(IAIN,NTALN) = ' '
          goto 3000
  170   CONTINUE
      ENDIF 
 
        IF (NA.EQ.20) THEN
          mm=modcol(1,2,irei)
          KK=NREAEI(irei)
          WRITE (CNO,'(I4)') IREI
          WRITE (CN1,'(I4)') KK
          TXTPLS(IAIN,NTALN) =
     .      'ELECTRON IMPACT REACTION RATE COEFFICIENT IREI ='//CNO
     .      //' KK='//CN1
          TXTPSP(IAIN,NTALN) = TEXTS(ISP)// ' on ELECTRONS'
          TXTPUN(IAIN,NTALN) = 'A.U. (0.612E-8 cm3/s)   '

          if (mm.eq.1) then
            DO 1720 ICELL=1,NSBOX
              if (lgvac(icell,npls+1)) cycle
              IF (NSTORDR >= NRAD) THEN
                TBEI=TABEI1(irei,ICELL)
              ELSE
                TBEI=EIRENE_FTABEI1(IREI,ICELL)
              END IF
              ADIN(IAIN,ICELL)=TBEI/(DEIN(ICELL)+EPS30)/AU
 1720       CONTINUE

            GOTO 5000 !done
          else  !  mm= MODCOL(1,2,irei)=2, not ready
            goto 3000
          endif

        ELSEIF (NA.EQ.21) THEN
C  ELECTRON IMPACT ENERGY LOSS RATE COEFFICIENT NO. IREI
          mm=modcol(1,4,irei)
          kk=nelrei(irei)
c
          WRITE (CNO,'(I4)') IREI
          WRITE (CN1,'(I4)') KK
          TXTPLS(IAIN,NTALN) =
     .      'ELECTRON IMPACT ENERGY LOSS RATE COEFFICIENT IREI ='//CNO
     .      //' KK='//CN1
          TXTPSP(IAIN,NTALN) = TEXTS(ISP)// ' on ELECTRONS'
          TXTPUN(IAIN,NTALN) = 'eV A.U. (0.612E-8 cm3/s)'
          if (mm.eq.1) then
            DO 1721 ICELL=1,NSBOX
              if (lgvac(icell,npls+1)) cycle
              IF (NSTORDR >= NRAD) THEN
                TBEI=TABEI1(irei,ICELL)
                ELEI=EELEI1(irei,ICELL)
              ELSE
                TBEI=EIRENE_FTABEI1(IREI,ICELL)
                ELEI=EIRENE_FEELEI1(IREI,ICELL)
              END IF
              RATE=TBEI/(DEIN(ICELL)+EPS30)/AU
              ADIN(IAIN,ICELL)=ELEI*RATE
 1721       CONTINUE

            goto 5000
          else  !  mm= MODCOL(1,4,irei)=2, not ready
            goto 3000
          endif

c.............................................................................

C  CHARGE EXCHANGE RATE COEFFICIENT NO. IRCX
        ELSEIF (NA.EQ.22.OR.NA.EQ.23) THEN
          ircx=ns
          kk=NREACX(ircx)
c find collision partners corresponding to process ircx: IPL and ISP
          IPL=0
c  first: try atoms
          LATCX: do iat=1,natmi
          do iacx=1,NACXI(iat)
            if (IRCX.eq.LGACX(IAT,IACX,0)) then
              IPL =LGACX(IAT,IACX,1)
              ISP=NSPH+IAT
              RMASSS=RMASSA(IAT)
              GOTO 172
            endif
          enddo
          enddo LATCX
c  ircx is a process for atom iat, colliding with bulk ipl

c  next: try molecules
          LMLCX: do iml=1,nmoli
          do imcx=1,NMCXI(iml)
            if (IRCX.eq.LGMCX(IML,IMCX,0)) then
              IPL =LGMCX(IML,IMCX,1)
              ISP=NSPA+IML
              RMASSS=RMASSM(IML)
              goto 172
            endif
          enddo
          enddo LMLCX
c  ircx is a process for molecule iml, colliding with bulk ipl

c  next: try test ions
          LIOCX: do iio=1,nioni
          do iicx=1,NICXI(iio)
            if (IRCX.eq.LGICX(IIO,IICX,0)) then
              IPL =LGICX(IIO,IICX,1)
              ISP=NSPAM+IIO
              RMASSS=RMASSI(IIO)
              goto 172
            endif
          enddo
          enddo LIOCX
c  ircx is a process for test ion iio, colliding with bulk ipl
c
c  no interacting particle species found
          TXTPLS(IAIN,NTALN) =
     .      'CHARGE EXCHANGE REACTION RATE COEFFICIENT IRCX ='//CNO
     .      //' KK='//CN1
          TXTPSP(IAIN,NTALN) = 'unidentified species    '
          TXTPUN(IAIN,NTALN) = ' '
          goto 3000
  172     CONTINUE
        ENDIF

        IF (NA.EQ.22) THEN
          mm=modcol(3,2,ircx)
          kk=NREACX(ircx)
          WRITE (CNO,'(I4)') IRCX
          WRITE (CN1,'(I4)') KK
          TXTPLS(IAIN,NTALN) =
     .      'CHARGE EXCHANGE REACTION RATE COEFFICIENT IRCX ='//CNO
     .      //' KK='//CN1
          TXTPSP(IAIN,NTALN) = TEXTS(ISP)//' on '//TEXTS(NSPAMI+IPL)
          TXTPUN(IAIN,NTALN) = 'A.U. (0.612E-8 cm3/s)   '

          if (mm.eq.1) then
            DO 1722 ICELL=1,NSBOX
              if (lgvac(icell,ipl)) cycle
              IF (NSTORDR >= NRAD) THEN
                TBCX=TABCX3(IRCX,ICELL,1)
              ELSE
                TBCX=EIRENE_FTABCX3(IRCX,ICELL)
              END IF
              ADIN(IAIN,ICELL)=TBCX/(diin(ipl,icell)+eps30)/AU
 1722       CONTINUE

            GOTO 5000 !done

          ELSEIF (MM.EQ.2) THEN
C  USE EB (ENERGY OF TEST PARTICLE) = 1.5 TI
            IPLTI = MPLSTI(IPL)
            FP = 0._DP
            RCMIN = -HUGE(1._DP)
            RCMAX = HUGE(1._DP)
            EARRH = 0._DP
            ND2   = 9
c   TEST PARTICLE VELOCITY NOT KNOWN HERE, TAKE Tn = Ti, and apply mass scaling
c      MASST(KK)=  TARGET MASS FOR CROSS-SECTION, BEAM MASS FOR BEAM MAXWELLIAN RATE COEFF.
            EBFAC= MASST(KK)*PMASSA/RMASSS
            DO ICELL=1,NSBOX
              if (lgvac(icell,ipl)) cycle
c  in fpath we use: ELB=MAX(-2.3_DP,LOG(PVELQ(IPLSV))+EEFCX(IRCX))
              ELB=log(max(EMIN,1.5_DP*TIIN(iplti,icell)*EBFAC))
              IF (NSTORDR >= NRAD) THEN
                TBCX3(1:NSTORDT) = TABCX3(IRCX,ICELL,1:NSTORDT)
                LEXP = .TRUE.
                RATE = EIRENE_SNGL_POLY(TBCX3,ELB,RCMIN,RCMAX,FP,0,0,
     .                                  EARRH,TRCAMD,LEXP)
              ELSE
! CALCULATE RATE COEFFICIENT ON THE FLY
CDR  THIS SHOULD BE DONE IN FTABCX3.  NOT READY
                KK=NREACX(IRCX)
                TII=TIINL(IPLTI,ICELL)+ADDCX(IRCX,IPL)
                RATE = EIRENE_RATE_COEFF(KK,ICELL,TII,ELB,.FALSE.,0)
     .                 + DIINL(IPL,ICELL) + FACRCX(IRCX,2)
                rate=exp(rate)
              ENDIF
              ADIN(IAIN,ICELL)=rate/(diin(ipl,icell)+eps30)/AU
            enddo
            goto 5000  !done

          ELSE ! mm=modcol(3,2,ircx).gt.2: NOT READY
            GOTO 3000
          ENDIF

        ELSEIF (NA.EQ.23) THEN
C  BULK ION IMPACT ENERGY LOSS RATE COEFFICIENT NO. IRCX
          mm=modcol(3,4,ircx)
          kk=NELRCX(IRCX)
c
          WRITE (CNO,'(I4)') IRCX
          WRITE (CN1,'(I4)') KK
          TXTPLS(IAIN,NTALN) =
     .      ' BULK ION IMPACT ENERGY LOSS RATE COEFFICIENT IRCX ='//CNO
     .      //' KK='//CN1
          TXTPSP(IAIN,NTALN) = TEXTS(ISP)//' on '//TEXTS(NSPAMI+IPL)
          TXTPUN(IAIN,NTALN) = 'eV A.U. (0.612E-8 cm3/s)'
          if (mm.eq.1) then
          DO 1723 ICELL=1,NSBOX
            if (lgvac(icell,ipl)) cycle
            IF (NSTORDR >= NRAD) THEN
              EPCX=EPLCX3(IRCX,ICELL,1)
              RATE=TABCX3(IRCX,ICELL,1)/(diin(ipl,icell)+EPS30)/AU
            ELSE
              EPCX=EIRENE_FEPLCX3(IRCX,ICELL)
              RATE=EIRENE_FTABCX3(IRCX,ICELL)/(diin(ipl,icell)+EPS30)/AU
            END IF
            ADIN(IAIN,ICELL)=EPCX*RATE
 1723     CONTINUE

            goto 5000
          else  !  mm= MODCOL(3,4,ircx)>=2, not ready
            GOTO 3000
          endif
c........................................................................

C  ELASTIC COLLISION RATE COEFFICIENT NO. IREL
        ELSEIF (NA.EQ.24.or.NA.EQ.25) THEN
          irel=ns
          kk=NREAEL(irel)
c find collision partners corresponding to process irel: IPL and ISP
          IPL=0
c  first: try atoms
          LATEL: do iat=1,natmi
          do iael=1,NAELI(iat)
            if (IREL.eq.LGAEL(IAT,IAEL,0)) then
              IPL =LGAEL(IAT,IAEL,1)
              ISP=NSPH+IAT
              RMASSS=RMASSA(IAT)
              GOTO 174
            endif
          enddo
          enddo LATEL
c  irel is a process for atom iat, colliding with bulk ipl

c  next: try molecules
          LMLEL: do iml=1,nmoli
          do imel=1,NMELI(iml)
            if (IREL.eq.LGMEL(IML,IMEL,0)) then
              IPL =LGMEL(IML,IMEL,1)
              ISP=NSPA+IML
              RMASSS=RMASSM(IML)
              goto 174
            endif
          enddo
          enddo LMLEL
c  irel is a process for molecule iml, colliding with bulk ipl

c  next: try test ions
          LIOEL: do iio=1,nioni
          do iiel=1,NIELI(iio)
            if (IREL.eq.LGIEL(IIO,IIEL,0)) then
              IPL =LGIEL(IIO,IIEL,1)
              ISP=NSPAM+IIO
              RMASSS=RMASSI(IIO)
              goto 174
            endif
          enddo
          enddo LIOEL
c  irel is a process for test ion iio, colliding with bulk ipl

c  no interacting particle species found
          TXTPLS(IAIN,NTALN) =
     .      'ELASTIC REACTION RATE COEFFICIENT IREL ='//CNO
     .      //' KK='//CN1
          TXTPSP(IAIN,NTALN) = 'unidentified species    '
          TXTPUN(IAIN,NTALN) = ' '
          GOTO 3000
  174     CONTINUE
        ENDIF

        IF (NA.EQ.24) THEN
          mm=modcol(5,2,irel)
          kk=NREAEL(irel)
          WRITE (CNO,'(I4)') IREL
          WRITE (CN1,'(I4)') KK
          TXTPLS(IAIN,NTALN) =
     .      'ELASTIC REACTION RATE COEFFICIENT IREL ='//CNO
     .      //' KK='//CN1
          TXTPSP(IAIN,NTALN) = TEXTS(ISP)// ' on '// TEXTS(NSPAMI+IPL)
          TXTPUN(IAIN,NTALN) = 'A.U. (0.612E-8 cm3/s)   '

          if (mm.eq.1) then
            IPLTI = MPLSTI(IPL)
            DO 1724 ICELL=1,NSBOX
              if (lgvac(icell,ipl)) cycle
              IF (NSTORDR >= NRAD) THEN
                TBEL=TABEL3(IREL,ICELL,1)
              ELSE
cdr  here should be call to ftabel3,  to be done
                KK=NREAEL(IREL)
                TII=TIINL(IPLTI,ICELL)+ADDEL(IREL,IPL)
                TBEL = EIRENE_RATE_COEFF(KK,ICELL,TII,0._DP,.TRUE.,0)*
     .                 DIIN(IPL,ICELL)
              END IF
              ADIN(IAIN,ICELL)=TBEL/(diin(ipl,icell)+eps30)/AU
 1724       CONTINUE
            GOTO 5000  !done

          ELSEIF (mm.EQ.2) THEN
C  USE EB (ENERGY OF TEST PARTICLE) = 1.5 TI
            IPLTI = MPLSTI(IPL)
            FP = 0._DP
            RCMIN = -HUGE(1._DP)
            RCMAX = HUGE(1._DP)
            EARRH = 0._DP
            ND2   = 9
c   TEST PARTICLE VELOCITY NOT KNOWN HERE, TAKE Tn = Ti, and apply mass scaling
c      MASST(KK)=  TARGET MASS FOR CROSS-SECTION, BEAM MASS FOR BEAM MAXWELLIAN RATE COEFF.
            EBFAC= MASST(KK)*PMASSA/RMASSS
            DO ICELL=1,NSBOX
              if (lgvac(icell,ipl)) cycle
c  in fpath we use: ELB=MAX(-2.3_DP,LOG(PVELQ(IPLSV))+EEFEL(IREL))
              ELB=log(max(EMIN,1.5_DP*TIIN(iplti,icell)*EBFAC))
              IF (NSTORDR >= NRAD) THEN
                TBEL3(1:NSTORDT) = TABEL3(IREL,ICELL,1:NSTORDT)
                LEXP = .TRUE.
                RATE = EIRENE_SNGL_POLY(TBEL3,ELB,RCMIN,RCMAX,FP,0,0,
     .                                  EARRH,TRCAMD,LEXP)
              ELSE
cdr  here should be call to ftabel3,  to be done
! CALCULATE RATE COEFFICIENT ON THE FLY
                KK=NREAEL(IREL)
                TII=TIINL(IPLTI,ICELL)+ADDEL(IREL,IPL)
                RATE = EIRENE_RATE_COEFF(KK,ICELL,TII,ELB,.FALSE.,0)
     .                 + DIINL(IPL,ICELL) + FACREL(IREL,2)
              END IF
              ADIN(IAIN,ICELL)=rate/(diin(ipl,icell)+eps30)/AU
            enddo
            goto 5000  !done

          ELSE ! mm=modcol(5,2,irel).gt.2:  NOT READY
            GOTO 3000
          ENDIF

        ELSEIF (NA.EQ.25) THEN
C  BULK ION IMPACT ENERGY LOSS RATE COEFFICIENT NO. IREL
          mm=modcol(5,4,irel)
          kk=NELREL(irel)
c
          WRITE (CNO,'(I4)') IREL
          WRITE (CN1,'(I4)') KK
          TXTPLS(IAIN,NTALN) =
     .      ' BULK ION IMPACT ENERGY LOSS RATE COEFFICIENT IREL ='//CNO
     .      //' KK='//CN1
          TXTPSP(IAIN,NTALN) = TEXTS(ISP)//' on '//TEXTS(NSPAMI+IPL)
          TXTPUN(IAIN,NTALN) = 'eV A.U. (0.612E-8 cm3/s)'
          if (mm.eq.1) then
            irel=ns
            DO 1725 ICELL=1,NSBOX
              if (lgvac(icell,ipl)) cycle
              IF (NSTORDR >= NRAD) THEN
                EPEL=EPLEL3(IREL,ICELL,1)
                RATE=TABEL3(IREL,ICELL,1)/(diin(ipl,icell)+EPS30)/AU
              ELSE
                EPEL=EIRENE_FEPLEL3(IREL,ICELL)
                RATE=EIRENE_FTABEL3(IREL,ICELL)/
     .               (diin(ipl,icell)+EPS30)/AU
              ENDIF
              ADIN(IAIN,ICELL)=EPEL*RATE
 1725       CONTINUE

            goto 5000
          else !  mm=modcol(5,4,irel)=2, not ready
            GOTO 3000
          endif

c.............................................................................

C  GENERAL ION IMPACT COLLISION RATE COEFFICIENT NO. IRPI
      ELSEIF (NA.EQ.26.OR.NA.EQ.27) THEN
        irpi=ns
        kk=nreapi(irpi)
c find collision partners corresponding to process irpi: IPL and ISP
        IPL=0
c  first: try atoms
        LATPI: do iat=1,natmi
        do iapi=1,NAPII(iat)
          if (IRPI.eq.LGAPI(IAT,IAPI,0)) then
            IPL =LGAPI(IAT,IAPI,1)
            ISP=NSPH+IAT
            RMASSS=RMASSA(IAT)
            GOTO 176
          endif
        enddo
        enddo LATPI
c  irpi is a process for atom iat, colliding with bulk ipl

c  next: try molecules
        LMLPI: do iml=1,nmoli
        do impi=1,NMPII(iml)
          if (IRPI.eq.LGMPI(IML,IMPI,0)) then
            IPL =LGMPI(IML,IMPI,1)
            ISP=NSPA+IML
            RMASSS=RMASSM(IML)
            goto 176
          endif
        enddo
        enddo LMLPI
c  irpi is a process for molecule iml, colliding with bulk ipl

c  next: try test ions
        LIOPI: do iio=1,nioni
        do iipi=1,NIPII(iio)
          if (IRPI.eq.LGIPI(IIO,IIPI,0)) then
            IPL =LGIPI(IIO,IIPI,1)
            ISP=NSPAM+IIO
            RMASSS=RMASSI(IIO)
            goto 176
          endif
        enddo
        enddo LIOPI
c  irpi is a process for test ion iio, colliding with bulk ipl

c  no interacting particle species found
        TXTPLS(IAIN,NTALN) =
     .    'HEAVY PARTICLE REACTION RATE COEFFICIENT IRPI ='//CNO
     .    //' KK='//CN1
        TXTPSP(IAIN,NTALN) = 'unidentified species    '
        TXTPUN(IAIN,NTALN) = ' '
        GOTO 3000
 176    CONTINUE
      ENDIF

      IF (NA.EQ.26) THEN
        mm=modcol(4,2,irpi)
        kk=nreapi(irpi)
        WRITE (CNO,'(I4)') IRPI
        WRITE (CN1,'(I4)') KK
        TXTPLS(IAIN,NTALN) =
     .    'BULK ION IMPACT REACTION RATE COEFFICIENT IRPI ='//CNO
     .    //' KK='//CN1
        TXTPSP(IAIN,NTALN) = TEXTS(ISP)// ' on '// TEXTS(NSPAMI+IPL)
        TXTPUN(IAIN,NTALN) = 'A.U. (0.612E-8 cm3/s)   '

        if (mm.eq.1) then
          DO 1726 ICELL=1,NSBOX
            if (lgvac(icell,ipl)) cycle
            IF (NSTORDR >= NRAD) THEN
              TBPI=TABPI3(IRPI,ICELL,1)
            ELSE
              TBPI=EIRENE_FTABPI3(IRPI,ICELL)
            ENDIF  
            ADIN(IAIN,ICELL)=TBPI/(diin(ipl,icell)+eps30)/AU
 1726     CONTINUE
          GOTO 5000  !DONE !

        ELSEIF (MM.EQ.2) THEN
C  USE EB (ENERGY OF TEST PARTICLE) = 1.5 TI
          IPLTI = MPLSTI(IPL)
          FP = 0._DP
          RCMIN = -HUGE(1._DP)
          RCMAX = HUGE(1._DP)
          EARRH = 0._DP
          ND2   = 9
c   TEST PARTICLE VELOCITY NOT KNOWN HERE, TAKE T_TEST = T-IPLS, and apply mass scaling
c      MASST(KK)=  TARGET MASS FOR CROSS-SECTION, BEAM MASS FOR BEAM MAXWELLIAN RATE COEFF.
          EBFAC= MASST(KK)*PMASSA/RMASSS
          DO ICELL=1,NSBOX
            if (lgvac(icell,ipl)) cycle
c  in fpath we use: ELB=MAX(-2.3_DP,LOG(PVELQ(IPLSV))+EEFPI(IRPI))
            ELB=log(max(EMIN,1.5_DP*TIIN(iplti,icell)*EBFAC))
            IF (NSTORDR >= NRAD) THEN
              TBPI3(1:NSTORDT) = TABPI3(IRPI,ICELL,1:NSTORDT)
              LEXP = .TRUE.
              RATE = EIRENE_SNGL_POLY(TBPI3,ELB,RCMIN,RCMAX,FP,0,0,
     .                                EARRH,TRCAMD,LEXP)
            ELSE
              TII=TIINL(IPLTI,ICELL)+ADDPI(IRPI,IPL)
! CALCULATE RATE COEFFICIENT "ON THE FLY"
              KK=NREAPI(IRPI)
              RATE = EIRENE_RATE_COEFF(KK,ICELL,TII,ELB,.FALSE.,0)
     .             + DIINL(IPL,ICELL) + FACRPI(IRPI,2)
            END IF
            ADIN(IAIN,ICELL)=rate/(diin(ipl,icell)+eps30)/AU
          enddo
          goto 5000  !done

        ELSE ! mm= modcol(4,2,irpi).gt.2:  NOT READY
            GOTO 3000
          ENDIF

        ELSEIF (NA.EQ.27) THEN
C  BULK ION IMPACT ENERGY LOSS RATE COEFFICIENT NO. IRCX
          mm=modcol(4,4,irpi)
          kk=NELRPI(irpi)
c
          WRITE (CNO,'(I4)') IRPI
          WRITE (CN1,'(I4)') KK
          TXTPLS(IAIN,NTALN) =
     .      ' BULK ION IMPACT ENERGY LOSS RATE COEFFICIENT IRPI ='//CNO
     .      //' KK='//CN1
          TXTPSP(IAIN,NTALN) = TEXTS(ISP)//' on '//TEXTS(NSPAMI+IPL)
          TXTPUN(IAIN,NTALN) = 'eV A.U. (0.612E-8 cm3/s)'
          if (mm.eq.1) then
          IF (NSTORDR >= NRAD) THEN
            DO 1727 ICELL=1,NSBOX
              if (lgvac(icell,ipl)) cycle
              RATE=TABPI3(IRPI,ICELL,1)/(diin(ipl,icell)+EPS30)/AU
              ADIN(IAIN,ICELL)=EPLPI3(irpi,ICELL,1)*RATE
 1727       CONTINUE
          ELSE
            WRITE (IUNOUT,*) 'BULK ION IMPACT LOSS RATE NOT READY'
            WRITE (IUNOUT,*) 'FOR CALCULATION ON THE FLY'
            ADIN(IAIN,1:NSBOX)=0._DP
          END IF 

          goto 5000
        else  !  mm=modcol(4,4,irpi)>=2, not ready
          GOTO 3000
        endif

c.................................................................

C  RECOMBINATION RATE COEFFICIENT NO. IRRC
        ELSEIF (NA.EQ.28.OR.NA.EQ.29) THEN
          irrc=ns
          mm=modcol(6,2,irrc)
          kk=NREARC(irrc)
c find collision partners corresponding to process irrc: IPL AND ISP
          IPL=0  ! ELECTRONS

c  here: only try bulk ions:
          LPLRC: do ipl=1,nplsi
          do iprc=1,NPRCI(ipl)
            if (IRRC.eq.LGPRC(IPL,IPRC)) then
              ISP=NSPAMI+IPL
              GOTO 178
            endif
          enddo
          enddo LPLRC
c  irrc is a process for bulk ion ipl, colliding with electron
c
c  no interacting particle species found
          TXTPLS(IAIN,NTALN) =
     .      'RECOMBINATION REACTION RATE COEFFICIENT IRRC ='//CNO
     .      //' KK='//CN1
          TXTPSP(IAIN,NTALN) = 'unidentified species    '
          TXTPUN(IAIN,NTALN) = ' '
          goto 3000
        ENDIF
  178   CONTINUE

        IF (NA.EQ.28) THEN
          mm=modcol(6,2,irrc)
          kk=NREARC(irrc)
          WRITE (CNO,'(I4)') IRRC
          WRITE (CN1,'(I4)') NREARC(IRRC)
          TXTPLS(IAIN,NTALN) =
     .      'RECOMBINATION REACTION RATE COEFFICIENT IRRC ='//CNO
     .      //' KK='//CN1
          TXTPSP(IAIN,NTALN) = TEXTS(ISP)// ' on ELECTRONS'
          TXTPUN(IAIN,NTALN) = 'A.U. (0.612E-8 cm3/s)   '

          if (mm.eq.1) then
            DO 1728 ICELL=1,NSBOX
              if (lgvac(icell,npls+1)) cycle
              IF (NSTORDR >= NRAD) THEN
                TBRC=TABRC1(irrc,ICELL)
              ELSE
                TBRC=EIRENE_FTABRC1(irrc,ICELL)
              END IF
              ADIN(IAIN,ICELL)=TBRC/(DEIN(ICELL)+EPS30)/AU
 1728       CONTINUE

            goto 5000
          else  !  mm= MODCOL(1,2,irei)=2, not ready
            goto 3000
          endif

        ELSEIF (NA.EQ.29) THEN
          mm=modcol(6,4,irrc)
          kk=nelrrc(irrc)
c  recombination energy loss rate coefficient no. irrc
c
          WRITE (CNO,'(I4)') IRRC
          WRITE (CN1,'(I4)') NELRRC(IRRC)
          TXTPLS(IAIN,NTALN) =
     .      'RECOMBINATION ENERGY LOSS RATE COEFFICIENT IRRC ='//CNO
     .      //' KK='//CN1
          TXTPSP(IAIN,NTALN) = TEXTS(ISP)// ' on ELECTRONS'
          TXTPUN(IAIN,NTALN) = 'eV A.U. (0.612E-8 cm3/s)'
          irrc=ns
          if (mm.eq.1) then
            DO 1729 ICELL=1,NSBOX
              if (lgvac(icell,npls+1)) cycle
c             RATE=TABRC1(irrc,ICELL)/(DEIN(ICELL)+EPS30)/AU
cdr  distinct from eelei1:  here eelrc1 already contains tabrc1 as factor
              RATE=1.0/(DEIN(ICELL)+EPS30)/AU
              IF (NSTORDR >= NRAD) THEN
                ELRC=EELRC1(irrc,ICELL)
              ELSE
                ELRC=EIRENE_FEELRC1(irrc,ICELL)
              END IF
              ADIN(IAIN,ICELL)=ELRC*RATE
 1729       CONTINUE
            goto 5000
          else  !  mm= MODCOL(6,4,irrc)=2, not ready
            goto 3000
          endif

        ENDIF

 3000   CONTINUE  !  unfinished option, or error

        if (mm.ne.0) then
          call eirene_leer(1)
          WRITE (iunout,*) 'ERROR IN AMDIAG, OPTION NOT READY '
          write (iunout,'(1X,A72)') TXTPLS(IAIN,NTALN)
          write (iunout,'(1X,A24)') TXTPSP(IAIN,NTALN)
          WRITE (iunout,*) 'IAIN, NS,NA      ', IAIN,NS,NA
          WRITE (iunout,*) 'PROCESS NO. KK, MODCOL(.,.,.)   ', KK,MM
          GOTO 190
        else   !mm = 0,  reaction kk has not been assigned to any particle
          call eirene_leer(1)
          WRITE (iunout,*) 'ERROR IN AMDIAG, ',
     .                     'PROCESS KK NOT ASSIGNED TO ANY PARTICLE '
          write (iunout,'(1X,A72)') TXTPLS(IAIN,NTALN)
          write (iunout,'(1X,A24)') TXTPSP(IAIN,NTALN)
          WRITE (iunout,*) 'IAIN, NS,NA      ', IAIN,NS,NA
          WRITE (iunout,*) 'PROCESS NO. KK, MODCOL(.,.,.)   ', KK,MM
          GOTO 190
        endif

 5000 CONTINUE
      CALL eirene_leer(1)
      WRITE (iunout,*) 'AMDIAG: ADDITIONAL INPUT TALLY ADIN(IAIN) SET'
      write (iunout,'(1X,A72)') TXTPLS(IAIN,NTALN)
      write (iunout,'(1X,A24)') TXTPSP(IAIN,NTALN)
      WRITE (iunout,*) 'IAIN, NS,NA      ', IAIN, NS,NA
      WRITE (iunout,*) 'PROCESS NO. KK, MODCOL(.,.,.)   ', KK,MM

  190 CONTINUE

      RETURN
      END SUBROUTINE EIRENE_AMDIAG
