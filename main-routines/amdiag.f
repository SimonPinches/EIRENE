cdr  feb, 16., 2015, added: naint=22, modcol=2 option, EB=1.5 Ti

CDR:  A&M Data diagnostics routine, added in Jan. 2014
C PUT SELECTED EIRENE ATOMIC DATA FIELDS ONTO ADIN-ARRAY FOR OUTPUT
C  ADIN CONTAINES RATE COEFFICIENTS (VOL/TIME) IN ATOMIC UNITS
c 
c  modcol=1: only dependent on local background data, not on test particle parameters
c            tabcx3(...,1),tabel3(...1),tabpi3(...1),tabds1(...) 
c            are rates, density of impacting bulk ion included, 
c            so here we divide again by ne or ni.
c
c  modcol=2: depends on Eb = energy of impacting test particle
c            tabcx3(...,1,nend),tabel3(...,1,nend),tabpi3(...,1,nend) 
c               are ln(rate)
c               with ln(rate)= sum_i=1^nend  ln^i(Eb) tab..3(...,i)
c  
c
C
C  RATE COEFFCIENTS IN ATOMIC UNITS FOR REACTION RATE COEFFICIENTS: A0^2 V0 = 0.612E-08 CM^3-S
C  TO CONVERT THESE ADDITIONAL TALLIES ADIN INTO UNITS OF 1/S, DIVIDE ADIN BY 0.612 e-08
c
c  done for naint=20,22,24,26 and modcol=1

c  naint=20:   Tabds1(irei,....) electron impact collision rate, 1/s --> cm^3/s      ! done 
c  naint=21:   eelds1(irei,....) electron cooling rate           eV/s --> cm^3 eV/s  ! not ready

c  naint=22:   Tabcx3(ircx,..,1) charge exchange collision rate, 1/s --> cm^3/s      ! done 
c  naint=23:   eplcx3(ircx,..,1) cx        energy weighted rate, eV/s --> cm^3 eV/s  ! not ready

c  naint=24:   Tabel3(irel,..,1) elastic         collision rate, 1/s --> cm^3/s      ! done 
c  naint=25:   eplel3(irel,..,1) elastic   energy weighted rate, eV/s --> cm^3 eV/s  ! not ready


c  naint=26:   Tabpi3(irpi,..,1) heavy particle imp.  coll.rate, 1/s --> cm^3/s      ! done 
c  naint=27:   eplpi3(irpi,..,1) ditto,    energy weighted rate, eV/s --> cm^3 eV/s  ! not ready


c  to be done:  nomenclature nmdsi(iml) --> nmeii(iml)
c                            nidsi(iio) --> nieii(iio)

      SUBROUTINE EIRENE_AMDIAG
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CCONA
      USE EIRMOD_CGRID
      USE EIRMOD_COMXS
      USE EIRMOD_CTEXT
      USE EIRMOD_CCOUPL

      IMPLICIT NONE

      REAL(DP) :: AU, ELB, EXPO, FP, RCMIN,RCMAX,
     .            TBCX3(9),
     .            EIRENE_SNGL_POLY
      INTEGER :: NS,NA,IAIN,MM,KK,
     .           irei,ircx,irpi,irel,irrc,
     .           iat,iml,iio,ipl,isp,iplti,
     .           icell,iapi,impi,iipi,iacx,imcx,iicx,
     .           iael,imel,iiel,iaei,imei,iiei
      CHARACTER(4) :: CNO, CN1



      AU=0.6120D-08

      IF (NSTORDR < NRAD) RETURN

      DO 190 IAIN=1,NAINI
        NS=NAINS(IAIN)    !  ns stands for ircx,irei,irel,irpi,irpi,..., internal number of process
        NA=NAINT(IAIN)    !  na stands for tally:  TAB..3(...),  EPL..3(...)

c  currently:  only tabcx3, tabel3 and tabpi3 are available, and only for modcol(..,2,ns)=1
c
        IF (NA.EQ.20.OR.NA.EQ.21) THEN       
c  electron impact rate coefficient no. irei
          irei=ns
          mm=modcol(1,2,irei)
          kk=NREAEI(irei)
c find collision partners corresponding to process irei: IPL ANBD ISP
          IPL=0  ! ELECTRONS
          
c  first: try atoms
          LATEI: do iat=1,natmi
          do iaei=1,NAEII(iat)
            if (IREI.eq.LGAEI(IAT,IAEI)) then
              ISP=NSPH+IAT
              GOTO 170
            endif
          enddo
          enddo LATEI
c  irei is a process for atom iat, colliding with electron

c  next: try molecules
          LMLEI: do iml=1,nmoli
          do imei=1,NMDSI(iml)
            if (IREI.eq.LGMEI(IML,IMEI)) then
              ISP=NSPA+IML
              goto 170
            endif
          enddo
          enddo LMLEI
c  irei is a process for molecule iml, colliding with electron

c  next: try test ions
          LIOEI: do iio=1,nioni
          do iiei=1,NIDSI(iio)
            if (IREI.eq.LGIEI(IIO,IIEI)) then
              ISP=NSPAM+IIO
              goto 170
            endif
          enddo
          enddo LIOEI
c  irei is a process for test ion iio, colliding with electron          
c
c  no interacting particle species found 
          TXTPLS(IAIN,NTALN) = 
     .      'ELECTRON IMPACT REACTION  IREI ='//CNO
     .      //' KK='//CN1            
          TXTPSP(IAIN,NTALN) = 'un-identified species on ELECTRONS'     
          TXTPUN(IAIN,NTALN) = ' '
          goto 171       
        ENDIF
170     CONTINUE

        IF (NA.EQ.20) THEN
          mm=modcol(1,2,irei)
          kk=NREAEI(irei)
          WRITE (CNO,'(I4)') IREI
          WRITE (CN1,'(I4)') NREAEI(IREI)
          TXTPLS(IAIN,NTALN) = 
     .      'ELECTRON IMPACT REACTION RATE COEFFICIENT IREI ='//CNO
     .      //' KK='//CN1            
          TXTPSP(IAIN,NTALN) = TEXTS(ISP)// ' on ELECTRONS'     
          TXTPUN(IAIN,NTALN) = 'A.U. (0.612 E-8 cm3/s)'
 
          if (mm.eq.1) then
            DO 1720 ICELL=1,NSBOX
              ADIN(IAIN,ICELL)=TABDS1(irei,ICELL)/(DEIN(ICELL)+EPS30)/AU
1720        CONTINUE

            goto 180
          else  !  mm= MODCOL(1,2,irei)=2, not ready            
            goto 171
          endif

        ELSEIF (NA.EQ.21) THEN
          mm=modcol(1,3,irei)
          kk=nelrei(irei)
c  electron impact energy loss rate coefficient no. irei
c  not ready
          WRITE (CNO,'(I4)') IREI
          WRITE (CN1,'(I4)') NREAEI(IREI)
          TXTPLS(IAIN,NTALN) = 
     .      'ELECTRON IMPACT ENERGY LOSS RATE COEFFICIENT IREI ='//CNO
     .      //' KK='//CN1            
          TXTPSP(IAIN,NTALN) = TEXTS(ISP)// ' on ELECTRONS'     
          TXTPUN(IAIN,NTALN) = 'eV times A.U. (0.612 E-8 cm3/s)'
          irei=ns
          if (mm.eq.1) then
            DO 1721 ICELL=1,NSBOX
              ADIN(IAIN,ICELL)=EELDS1(irei,ICELL)
1721        CONTINUE
          else  !  mm= MODCOL(1,2,irei)=2, not ready
            goto 171
          endif     
          GOTO 171
        endif

        IF (NA.EQ.22.OR.NA.EQ.23) THEN          
c  charge exchange rate coefficient no. ircx   
          ircx=ns
          mm=modcol(3,2,ircx)
c find collision partners corresponding to process ircx
          IPL=0
c  first: try atoms
          LATCX: do iat=1,natmi
          do iacx=1,NACXI(iat)
            if (IRCX.eq.LGACX(IAT,IACX,0)) then
              IPL =LGACX(IAT,IACX,1)
              ISP=NSPH+IAT
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
              goto 172
            endif
          enddo
          enddo LIOCX
c  ircx is a process for test ion iio, colliding with bulk ipl
c
c  no interacting particle species found 
          TXTPLS(IAIN,NTALN) = 
     .      'CHARGE EXCHANGE REACTION  IRCX ='//CNO
     .      //' KK='//CN1            
          TXTPSP(IAIN,NTALN) = 'un-identified colliding species'     
          TXTPUN(IAIN,NTALN) = ' '
          goto 171                
172       CONTINUE      
        ENDIF

        IF (NA.EQ.22) THEN
          mm=modcol(3,2,ircx)
          kk=NREACX(ircx)
          WRITE (CNO,'(I4)') IRCX
          WRITE (CN1,'(I4)') NREACX(IRCX)
          TXTPLS(IAIN,NTALN) = 
     .      'CHARGE EXCHANGE REACTION RATE COEFFICIENT IRCX ='//CNO
     .      //' KK='//CN1
          TXTPSP(IAIN,NTALN) = TEXTS(ISP)//' on '//TEXTS(NSPAMI+IPL)     
          TXTPUN(IAIN,NTALN) = 'A.U. (0.612 E-8 cm3/s)'

          if (mm.eq.1) then            
            DO 1722 ICELL=1,NSBOX            
              ADIN(IAIN,ICELL)=
     .        TABCX3(IRCX,ICELL,1)/(diin(ipl,icell)+eps30)/AU
1722        CONTINUE
            
            GOTO 180 !done
          ELSEIF (MM.EQ.2) THEN
C  USE EB (ENERGY OF TEST PARTICLE) = 1.5 TI
            IPLTI = MPLSTI(IPL)
            FP = 0._DP
            RCMIN = -HUGE(1._DP)
            RCMAX = HUGE(1._DP)
c           ELB=MAX(-2.3_DP,LOG(PVELQ(IPLSV))+EEFCX(IRCX))
            DO ICELL=1,NSBOX
              ELB=log(1.5*TIIN(iplti,icell))
              TBCX3(1:NSTORDT) = TABCX3(IRCX,ICELL,1:NSTORDT)             
              EXPO = EIRENE_SNGL_POLY(TBCX3,ELB,RCMIN,RCMAX,FP,0,0)
              ADIN(IAIN,ICELL)=
     .        exp(expo)/(diin(ipl,icell)+eps30)/AU
            enddo
            goto 180  !done
          ELSE ! MODCOL.gt.2:  NOT READY            
            GOTO 171
          ENDIF


        ELSEIF (NA.EQ.23) THEN
          mm=modcol(3,3,ircx)
          kk=NELRCX(IRCX)
c  not ready
          DO 1723 ICELL=1,NSBOX            
            ADIN(IAIN,ICELL)=EPLCX3(IRCX,ICELL,1)
1723      CONTINUE
          
          GOTO 171


        ELSEIF (NA.EQ.24.or.NA.EQ.25) THEN          
c  elastic collision rate coefficient no. irel
          irel=ns 
          mm=modcol(5,2,irel)
          kk=NREAEL(irel)
c find collision partners corresponding to process irel
          IPL=0
c  first: try atoms
          LATEL: do iat=1,natmi
          do iael=1,NAELI(iat)
            if (IREL.eq.LGAEL(IAT,IAEL,0)) then
              IPL =LGAEL(IAT,IAEL,1)
              ISP=NSPH+IAT
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
              goto 174
            endif
          enddo
          enddo LIOEL
c  irel is a process for test ion iio, colliding with bulk ipl

c  no interacting particle species found 
          TXTPLS(IAIN,NTALN) = 
     .      'ELASTIC REACTION  IREL ='//CNO
     .      //' KK='//CN1            
          TXTPSP(IAIN,NTALN) = 'un-identified colliding species'     
          TXTPUN(IAIN,NTALN) = ' '         
          GOTO 171
174       CONTINUE
        ENDIF

        IF (NA.EQ.24) THEN
          mm=modcol(5,2,irel)
          kk=NREAEL(irel)
          WRITE (CNO,'(I4)') IREL
          WRITE (CN1,'(I4)') NREAEL(IREL)
          TXTPLS(IAIN,NTALN) = 
     .      'ELASTIC REACTION RATE COEFFICIENT IREL ='//CNO
     .      //' KK='//CN1
          TXTPSP(IAIN,NTALN) = TEXTS(ISP)// ' on '// TEXTS(NSPAMI+IPL)     
          TXTPUN(IAIN,NTALN) = 'A.U. (0.612 E-8 cm3/s)'
          if (mm.eq.1) then
            DO 1724 ICELL=1,NSBOX            
              ADIN(IAIN,ICELL)=
     .        TABEL3(IREL,ICELL,1)/(diin(ipl,icell)+eps30)/AU
1724        CONTINUE
            GOTO 180
          ELSE ! MODCOL(...)=2:  NOT READY             
            
            GOTO 171
          ENDIF

        ELSEIF (NA.EQ.25) THEN
          mm=modcol(5,3,irel)
          kk=NELREL(irel)
c  not ready
          irel=ns
          DO 1725 ICELL=1,NSBOX            
            ADIN(IAIN,ICELL)=EPLEL3(NS,ICELL,1)
1725      CONTINUE
          
          GOTO 171

        ELSEIF (NA.EQ.26.OR.NA.EQ.27) THEN
          
c  general ion impact collision rate coefficient no. irpi
          IRPI=NS
          mm=modcol(4,2,irpi)
          kk=nreapi(irpi)
c find collision partners corresponding to process irpi
          IPL=0
c  first: try atoms
          LATPI: do iat=1,natmi
          do iapi=1,NAPII(iat)
            if (IRPI.eq.LGAPI(IAT,IAPI,0)) then
              IPL =LGAPI(IAT,IAPI,1)
              ISP=NSPH+IAT
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
              goto 176
            endif
          enddo
          enddo LIOPI
c  irpi is a process for test ion iio, colliding with bulk ipl
          
c  no interacting particle species found 
          TXTPLS(IAIN,NTALN) = 
     .      'HEAVY PARTICLE REACTION  IRPI ='//CNO
     .      //' KK='//CN1            
          TXTPSP(IAIN,NTALN) = 'un-identified colliding species'     
          TXTPUN(IAIN,NTALN) = ' '         
          GOTO 171
176       CONTINUE
        ENDIF

        IF (NA.EQ.26) THEN
          mm=modcol(4,2,irpi)
          kk=nreapi(irpi)
          WRITE (CNO,'(I4)') IRPI
          WRITE (CN1,'(I4)') NREAPI(IRPI)
          TXTPLS(IAIN,NTALN) = 
     .      'BULK ION IMPACT REACTION RATE COEFFICIENT IRPI ='//CNO
     .      //' KK='//CN1
          TXTPSP(IAIN,NTALN) = TEXTS(ISP)// ' on '// TEXTS(NSPAMI+IPL)     
          TXTPUN(IAIN,NTALN) = 'A.U. (0.612 E-8 cm3/s)'

          if (mm.eq.1) then   
            DO 1726 ICELL=1,NSBOX            
              ADIN(IAIN,ICELL)=
     .        TABPI3(NS,ICELL,1)/(diin(ipl,icell)+eps30)/AU
1726        CONTINUE
            GOTO 180  !DONE !

          ELSE   ! MODCOL(.....) =2:  NOT READY           
            
            GOTO 171
          ENDIF
        ELSEIF (NA.EQ.27) THEN
          mm=modcol(4,3,irpi)
          kk=NELRPI(irpi)
          irpi=ns
          DO 1727 ICELL=1,NSBOX            
            ADIN(IAIN,ICELL)=EPLPI3(irpi,ICELL,1)
1727      CONTINUE
          GOTO 171

        ELSEIF (NA.EQ.28) THEN
c  volume recombination rate coefficent no. irrc
c  not ready
          irrc=ns
          DO 1728 ICELL=1,NSBOX          
            ADIN(IAIN,ICELL)=TABRC1(irrc,ICELL)
1728      CONTINUE
          mm=0
          kk=0
          GOTO 171

        ELSEIF (NA.EQ.29) THEN
c  not ready
          irrc=ns
          DO 1729 ICELL=1,NSBOX
            ADIN(IAIN,ICELL)=EELRC1(irrc,ICELL)
1729      CONTINUE
          mm=0
          kk=0
          GOTO 171
        ENDIF

171     CONTINUE

        if (mm.ne.0) then
          call eirene_leer(1)
          WRITE (iunout,*) 'ERROR IN AMDIAG, OPTION NOT READY '
          write (iunout,'(A72)') txtpls(IAIN,NTALN)
          write (iunout,'(A72)') TXTPSP(IAIN,NTALN)
          WRITE (iunout,*) 'IAIN, NS,NA      ', IAIN,NS,NA
          WRITE (iunout,*) 'process no. KK, MODCOL(...)   ', KK,MM
          GOTO 190
        else   !mm = 0,  reaction kk has not been assgined to any particle
          call eirene_leer(1)
          WRITE (iunout,*) 'ERROR IN AMDIAG,', 
     .                     'PROCESS KK NOT ASSIGNED TO ANY PARTICLE '
          write (iunout,'(A72)') txtpls(IAIN,NTALN)
          write (iunout,'(A72)') TXTPSP(IAIN,NTALN)
          WRITE (iunout,*) 'IAIN, NS,NA      ', IAIN,NS,NA
          WRITE (iunout,*) 'process no. KK, MODCOL(...)   ', KK,MM
          GOTO 190
        endif 

180   CONTINUE
      CALL eirene_leer(1)
      WRITE (iunout,*) 'AMDIAG: ADDITIONAL INPUT TALLY ADIN(IAIN) SET'
      write (iunout,'(A72)') txtpls(IAIN,NTALN)
      write (iunout,'(A72)') TXTPSP(IAIN,NTALN)
      WRITE (iunout,*) 'IAIN, NS,NA      ', IAIN, NS,NA
      WRITE (iunout,*) 'process no. KK, MODCOL(...)   ', KK,MM
       

190   CONTINUE

      RETURN
      END

