c april05:  *sqrt(ze) moved from here (for cx spectra) into sigcx
c april06:  restriction to iphot.eq.isp in case of los-radiances
c
cdr aug.16:  to be done: psig: allocatable, psig(0,nspi), NSPI depends on NCHTAL option
c            option NCHTAL=4 is unfinished. print warning and return 
cdr nov.16:  avoid reading strata, in case of single stratum runs (NSTRAI=1)
c            set default ncheni=1 for nchtal=2 already in calling routine,
c            to avoid that chords are erroneously turned off there. 
cpb jul.17:  request from aug.16: psig and ARGST depends on NCHTAL done
cdr       :  commit PART 1: allocatable storage: ARGST, VPLOT, AA, XNTG
cdr       :  PART 2: automatic detection of 1st dimension ND: tb committed later 
cdr Oct 17  :
cdr from W.Zholobenko: add         He emission lines, new options NCHTAL=5       
cdr                    analogous to H emission lines,             NCHTAL=2 
cdr       : added ITP (type of relevant component)           
C
C
      SUBROUTINE EIRENE_SGNAL(ICHORI,IISTR,ISP,ITP,LCHOR)
C
C  THIS SUBROUTINE CALCULATES LINE INTEGRATED SIGNALS, USING THE EIRENE
C  VOLUME AVERAGED TALLIES AND THE PLASMA BACKGROUND DATA.
C  THERE MAY BE A CONTRIBUTION DIRECTLY FROM A PRIMARY SOURCE,
C  DUE TO DIRECT EMISSION FROM THE SOURCE INTO THE LINE OF SIGHT,
C  AS WELL AS A SECONDARY SOURCE (POST COLLISION) CONTRIBUTION, DUE TO
C  SCATTERING INTO THE LINE OF SIGHT

C  STEP 1)  FETCH THE APPROPRIATE STRATUM DATA (OR: SUM OVER STRATA) IISTR
C  STEP 2)  PREPARE DIRECT CONTRIBUTION FROM PRIMARY SOURCE (IF ANY)
C           (PROBABLY NOT READY)
C  STEP 3)  INTEGRATE ALONG LINE OF SIGHT, CALL LININT, AND LOOP OVER ENERGY/WAVELENGTH
C  STEP 4)  PROCESS LINE INTEGRALS: CURVE FITTING, SCALING, ETC..
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CPLOT
      USE EIRMOD_COMSIG
      USE EIRMOD_CGRID
      USE EIRMOD_CTRCEI
      USE EIRMOD_CSDVI
      USE EIRMOD_CSDVI_COP
      USE EIRMOD_CSDVI_BGK
      USE EIRMOD_COMPRT
      USE EIRMOD_COMSOU
      USE EIRMOD_COUTAU
      USE EIRMOD_COMXS
      USE EIRMOD_CSPEI
      USE EIRMOD_CGEOM
      USE EIRMOD_CTEXT
      USE EIRMOD_CUPD
 
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: ICHORI, ISP, ITP
      INTEGER :: IISTR
      REAL(DP), ALLOCATABLE, SAVE :: PSIG(:)
      REAL(DP) :: C1(3),C2(3),
     .          BUFFER(NCHOR,NCHEN),ESTART(NCHOR),ENDFIT(NCHOR),
     .          FP1(6), FP2(6), DUM(9)
      REAL(DP) :: ZE1, ZE2, ZSCALE, ZZ, EIRENE_SLOPE, STEIG, PMI, PMA,
     .            XMI, XMAX, XMIN, ZSI, TIMAX, ZE, SUMM, ADD, FAC32,
     .            TEF, DEF, DE, TE, ZDS, RATE, CHKSUM
     .           ,summt,addt, XMA, 
     .            RC1MIN, RC1MAX, RC2MIN, RC2MAX 
      REAL(DP) :: EIRENE_FTABRC1
      INTEGER :: I1, I2, IN, I, IS, NAC2, NBC2, ICHRD, IPVOT, NCHNI,
     .           IFIRST, JEN, NSPI, ISTR, ICOUNT, KK, IR, IZ,
     .           KREC, IRRC, MAXREC, IFLAG, ISPC,
     .           ICELL, 
     .           JFEX1MN, JFEX1MX, JFEX2MN, JFEX2MX
      INTEGER, SAVE :: MX_COMPO, ND
      LOGICAL :: NLVL(0:NSTRAI),LCHOR
      CHARACTER(8) :: FILNAM
      CHARACTER(4) :: H123
      CHARACTER(9) :: REAC
      CHARACTER(3) :: CRC
      CHARACTER(2) :: ELNAME
      TYPE(CELL_INFO), POINTER :: FIRST, CUR

      INTERFACE
        SUBROUTINE EIRENE_LININT
     .           (IFIRST,ICHORI,C1,C2,ICHRD,IPVOT,NBC2,NAC2,PEN,
     .            PSIG,TIMAX,ISP,NSPI,JEN,NCHNI)
        USE EIRMOD_PRECISION
        USE EIRMOD_PARMMOD
        INTEGER, INTENT(IN) :: IFIRST,ICHORI, ICHRD,IPVOT,NBC2,NAC2,ISP,
     .                         NSPI, JEN, NCHNI
        REAL(DP), INTENT(IN) :: C1(3),C2(3),PEN
        REAL(DP), INTENT(IN OUT) :: PSIG(0:)
        REAL(DP), INTENT(IN OUT) :: TIMAX
        END SUBROUTINE EIRENE_LININT
      END INTERFACE
C
      ISTRA=IISTR
      NCHNI=IABS(NCHENI)

cdr   IF ((NCHTAL(ICHORI).EQ.2).OR.(NCHTAL(ICHORI).EQ.5)) NCHNI=1  
cdr   this is now already done in calling routine diagno.f

cdr:  aug. 2016
cdr:  to be written: use emin1, emax1 to identify upper and lower state of a transition,
cdr                  as already described in manual, but apparently not programmed.
C
 
C
      IF (NCHTAL(ICHORI) .NE. 3) THEN
c  no volumetric tallies at all are needed in case NCHTAL=3
c
c  in single stratum runs no dump files for strata are needed.
c  all tallies are still in active storage
      IF (NSTRAI.EQ.1) THEN
        IISTR=1
        ISTRA=1
        IESTR=1
      ENDIF 
  
      IF (ISTRA.EQ.IESTR) THEN
C  NOTHING TO BE DONE

      ELSEIF (NFILEN.EQ.1.OR.NFILEN.EQ.2) THEN
        IESTR=ISTRA
        CALL EIRENE_RSTRT(ISTRA,NSTRAI,NESTM1,NESTM2,NADSPC,
     .             ESTIMV,ESTIMS,ESTIML,
     .             NSDVI1,SDVI1,NSDVI2,SDVI2,
     .             NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .             NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .             NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .             NSIGI_SPC,TRCFLE)
        IF (NLSYMP(ISTRA).OR.NLSYMT(ISTRA)) THEN
          CALL EIRENE_SYMET(ESTIMV,NVOLTL,NRAD,NR1ST,NP2ND,NT3RD,
     .               NLSYMP(ISTRA),NLSYMT(ISTRA))
        ENDIF
      ELSEIF ((NFILEN.EQ.6.OR.NFILEN.EQ.7).AND.ISTRA.EQ.0) THEN
        IESTR=ISTRA
        CALL EIRENE_RSTRT(ISTRA,NSTRAI,NESTM1,NESTM2,NADSPC,
     .             ESTIMV,ESTIMS,ESTIML,
     .             NSDVI1,SDVI1,NSDVI2,SDVI2,
     .             NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .             NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .             NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .             NSIGI_SPC,TRCFLE)
        IF (NLSYMP(ISTRA).OR.NLSYMT(ISTRA)) THEN
          CALL EIRENE_SYMET(ESTIMV,NVOLTL,NRAD,NR1ST,NP2ND,NT3RD,
     .               NLSYMP(ISTRA),NLSYMT(ISTRA))
        ENDIF
      ELSE
        WRITE (iunout,*) 'ERROR IN DIAGNO: DATA FOR STRATUM ISTRA= ',
     .                   ISTRA
        WRITE (iunout,*)
     .    'ARE NOT AVAILABLE. LINE INTEGRATION ABANDONNED'
        RETURN
      ENDIF

      END IF

C  STEP 1 DONE

 
cdr  intermediate part, 
cdr  proprietary option only. nchtal=4 option is not ready. 
cdr  complain and then return to calling program  
!  INTEGRATE SPECTRA COLLECTED IN THE CELLS ALONG THE LINE OF SIGHT
 
      IF ((NCHTAL(ICHORI) == 4) .AND. NLSTCHR(ICHORI)) THEN
 
        if (.not.associated(traj(ichori)%trj%cells)) then
           write (iunout,*) 'WARNING FROM MODULE: DIAGNO '
           write (iunout,*) 'error in SGNAL, NCHTAL=4 '
           write (iunout,*) 'no proper spectra found for NCHORI '
           call EIRENE_MASJ1('NCHORI= ',NCHORI)
           write (iunout,*) 'OPTION NOT READY, RETURN '
           return
        end if
C
        first => traj(ichori)%trj%cells
        cur => first
C
CDR : this next part is not generally valid, nor ready to use.
cdr : an attempt had been made, apparently, to use velocity resolved 
cdr : neutral distributions (the velocity component along the line of sight),
cdr : then to turn that into a doppler broadened line shape of the Ba-alpha line

        write (iunout,*) 'WARNING FROM MODULE: DIAGNO '
        write (iunout,*) 'error in proprietary section NCHTAL=4 '
        write (iunout,*) 'OPTION NOT READY, RETURN '
        return  

C 
C  RADIATIVE TRANSITION RATES (1/S)
C  BALMER ALPHA
        FAC32=4.410E7
 
        FILNAM='AMJUEL  '
        H123='H.12'
        CRC='OT '
        FP1 = 0._DP
        RC1MIN = -HUGE(1._DP)
        RC1MAX =  HUGE(1._DP)
        JFEX1MN = 0
        JFEX1MX = 0
        FP2 = 0._DP
        RC2MIN = -HUGE(1._DP)
        RC2MAX =  HUGE(1._DP)
        JFEX2MN = 0
        JFEX2MX = 0
C        
        ELNAME = 'H     '
        IZ=0


C
C  H(n=3)/H(n=1)
        REAC='2.1.5a   '
        REACDAT(NREACI+1)%LOTH = .FALSE.
        CALL EIRENE_SLREAC(NREACI+1,FILNAM,H123,REAC,CRC,
     .              RC1MIN, RC1MAX, FP1, JFEX1MN, JFEX1MX,
     .              RC2MIN, RC2MAX, FP2, JFEX2MN, JFEX2MX,
     .              ELNAME,IZ)
 
        chksum = 0._dp
        do
          ispc = cur%no_spect
          icell = cur%no_cell
          zds = cur%flight
 
          IF (NSTGRD(ICELL) > 0) GOTO 500
 
          IF (LGVAC(ICELL,NPLS+1)) GOTO 500
 
          TE=TEIN(ICELL)
          DE=DEIN(ICELL)
 
          DEF=LOG(DE*1.D-8)
          TEF=LOG(TE)
 
          call EIRENE_dbl_poly(REACDAT(NREACI+1)%OTH%POLY%DBLPOL,
     .                         TEF, DEF, RATE, DUM, 
     .                         RC1MIN, RC1MAX, FP1, JFEX1MN, JFEX1MX,
     .                         RC2MIN, RC2MAX, FP2, JFEX2MN, JFEX2MX,
     .                         TRCAMD)

          rate = exp(rate)
 
          fuffer(ichori,1:ncheni) = fuffer(ichori,1:ncheni) +
     .                              zds * rate * FAC32 /(4.*PIA) *
     .                              estiml(ispc)%pspc%spc(1:ncheni)
 
          chksum = chksum + zds * rate * FAC32 /(4.*PIA) *
     .                      estiml(ispc)%pspc%spcs
 
 500      continue
          cur => cur%nextc
!pb associated with two arguments tests if both arguments point to the same target
          if (associated(cur,first)) exit
        end do
 
        write (iunout,*) ' chord ', ichori
        write (iunout,*) ' integral along chord calculated from '
        write (iunout,*) ' spectra in the cells '
        write (iunout,*) ' integral ',chksum
        return
 
      END IF  !  END OF UNFINISHED NCHTAL=4 OPTION

C  STEP 2
C
C  PREPARE DIRECT (PRIMARY) EMISSION FROM SOURCE, INTO LINE OF SIGHT
C  THE SCATTERED (SECONDARY) CONTRIBUTION WILL BE DONE 
C  IN SUBR. SIGCX, SIGRAD, SIGHA, SIGUSR, ETC.
 
      NLVL=.FALSE.
      DO 10 ISTR=1,NSTRAI
C  PRIMARY SOURCE CONTRIBUTION, REFER TO INPUT BLOCK 7
C  (PRIMARY) VOLUME RECOMBINATION SOURCE:
C  HOWEVER, THIS STRATUM MIGHT NOT NECESSARILY HAVE BEEN ACTIVE?
        NLVL(ISTR)=NLVOL(ISTR).AND.NLPLS(ISTR)
10    CONTINUE
      NLVL(0)=ANY(NLVL(1:NSTRAI))
C.................................................................
      IF (NCHTAL(ICHORI).EQ.1) THEN
C.................................................................
C  FOR CX SIGNAL:  TO BE WRITTEN
C     CALL ZEROA2(RECADD,NATM,NRAD)
C     IF (NLVL) THEN
C       WRITE (iunout,*) 'WARNING:'
C       WRITE (iunout,*) 'VOLUME RECOMBINATION CONTRIBUTION TO SIGNAL:'
C       WRITE (iunout,*) 'FIRST GENERATION CONTRIBUTION: TO BE WRITTEN'
C  RECADD: #/S/CM**3
C  RECADD*ELCHA*VOL: AMP/CELL
C       DO 100 IPLS=1,NPLSI
C         DO 100 KREC=1,NPRCI(IPLS)
C           IATM=NATPRC(KREC)
C           IF (IATM.LE.0.OR.IATM.GT.NATMI) GOTO 100
C           DO 101 IR=1,NSBOX
C             RECADD(IATM,IR)=RECADD(IATM,IR)+
C    .                        TABRC1(KREC,IR)*DIIN(IPLS,IR)*ELCHA
101         CONTINUE
100     CONTINUE
        write (iunout,*) 'sgnal, cx: ichord,istra,sum ',
     .                      ichori,istra
        write (iunout,*) 'volumetric emission to be written'
C     ENDIF
C................................................................. 
C  FOR RADIANCE OF LINE ISP=IPHOT, IN STRATUM ISTR
      ELSEIF (NCHTAL(ICHORI).EQ.3) THEN
C.................................................................
        MAXREC=SUM(NPRCI(1:NPLSI))
        ALLOCATE(RECADD(MAXREC,NRAD))
        ALLOCATE(INTADD(3,MAXREC))
        RECADD=0._DP
        INTADD=0
        ICOUNT=0
        NCTSIG=0
        IF (NLVL(ISTRA)) THEN
C  RECADD: #/S/CM**3
C  RECADD*ELCHA*VOL: AMP/CELL
          IF (ISTRA.EQ.0) THEN
            WRITE (iunout,*) 'ISTRA=0 NOT READY IN SUBR. SGNAL '
            DEALLOCATE(RECADD)
            DEALLOCATE(INTADD)
            LCHOR=.FALSE.
            RETURN
          ENDIF
          IF (NSPEZ(ISTRA).LE.0.OR.NSPEZ(ISTRA).GT.NPLSI) THEN
            WRITE (iunout,*) 'NSPEZ(ISTRA) OUT OF RANGE IN SUBR. SGNAL '
            WRITE (iunout,*) 'ISTRA, NSPEZ ', ISTRA, NSPEZ(ISTRA)
            DEALLOCATE(RECADD)
            DEALLOCATE(INTADD)
            LCHOR=.FALSE.
            RETURN
          ENDIF
          IPLS=NSPEZ(ISTRA)
          IFLAG=0
          DO 130 KREC=1,NPRCI(IPLS)
            IRRC=LGPRC(IPLS,KREC)
            IPHOT=NPHPRC(IRRC)
            IF (IPHOT.EQ.0) IPHOT=NPHPRC_2(IRRC)
            IF (IPHOT.LE.0.OR.IPHOT.GT.NPHOTI) GOTO 130
            IF (IPHOT.NE.ISP) THEN
              IF (NPRCI(IPLS).EQ.1) THEN
                WRITE (IUNOUT,*) 'IN SGNAL: IPHOT.NE.ISP '
                WRITE (IUNOUT,*) 'ICHORD, IRRC, IPHOT, ISP ',
     .                            ICHORI, IRRC,IPHOT,ISP
                WRITE (IUNOUT,*) 'EMISSION FROM THIS LINE IS SKIPPED'
              ELSE
C  THIS PARTICULAR IRRC DOES NOT FIT TO SPECIES ISP FOR PRESENT CHORD
C  SEARCH FOR NEXT IRRC FOR THIS SAME SOURCE PARTICLE IPLS=NSPEZ(ISTRA)
              ENDIF
              GOTO 130
            ENDIF
            KK=NREARC(IRRC)
            ICOUNT=ICOUNT+1
            INTADD(1,ICOUNT)=KK
            INTADD(2,ICOUNT)=IPLS
            INTADD(3,ICOUNT)=IPHOT
            SUMM=0._DP
            summt=0.
            DO 131 IR=1,NSBOX
              ADD=0._DP
              ADDt=0._DP
              IF (NSTGRD(IR).EQ.0.AND..NOT.LGVAC(IR,IPLS)) THEN
                IF (NSTORDR >= NRAD) THEN
                  ADD=TABRC1(IRRC,IR)*DIIN(IPLS,IR)
                ELSE
                  ADD=EIRENE_FTABRC1(IRRC,IR)*DIIN(IPLS,IR)
                END IF
              else
                IF (NSTORDR >= NRAD) THEN
                  ADDt=TABRC1(IRRC,IR)*DIIN(IPLS,IR)
                ELSE
                  ADDt=EIRENE_FTABRC1(IRRC,IR)*DIIN(IPLS,IR)
                END IF
              END IF
              RECADD(ICOUNT,IR)=ADD
              SUMM=SUMM+ADD*VOL(IR)*ELCHA
              SUMMt=SUMMt+ADDt*VOL(IR)*ELCHA
131         CONTINUE
            IFLAG=1
            write (iunout,*) 'sgnal, rad: ichord,istra,icount,',
     .                                                     'vol-source',
     .                            ichori,istra,icount,summ,summt
            write (iunout,*) 'ipls, irrc ',texts(nspami+ipls),' ',irrc
            write (iunout,*) 'iphot, isp ',texts(iphot),' ',texts(isp)
130       CONTINUE
          IF (IFLAG.EQ.0) THEN
C  NO EMISSION FOUND FOR THIS STRATUM, turn off this chord
            LCHOR=.FALSE.
            write (iunout,*) 'sgnal, rad: ichord, istra, is turned off',
     .                               ichori, istra
            write (iunout,*) 'no photon emission rate found '
            DEALLOCATE(RECADD)
            DEALLOCATE(INTADD)
            return
          ENDIF
          NCTSIG=ICOUNT
        ELSE
c  nlvl is not true:
          LCHOR=.FALSE.
          write (iunout,*) 'sgnal, rad: ichord, istra, is turned off',
     .                               ichori, istra
          write (iunout,*) 'specified stratum is not a volume source '
          DEALLOCATE(RECADD)
          DEALLOCATE(INTADD)
          return
        ENDIF
C.................................................................
      ELSEIF ((NCHTAL(ICHORI).EQ.2).OR.(NCHTAL(ICHORI).EQ.5)) THEN
C.................................................................
        write (iunout,*) 'sgnal, emis: ichord,istra ',
     .                            ichori,istra
        write (iunout,*) 'volumetric line emission '
        write (iunout,*) 'contribution no. isp ',isp
      ENDIF

C   STEP 2 DONE

C
C     CALCULATE SIGNAL STRENGTHS, LOOP OVER ENERGY/WAVELENGTH, IF ANY
C
      IPVOT=IPIVOT(ICHORI)
      C1(1)=XPIVOT(ICHORI)
      C1(2)=YPIVOT(ICHORI)
      C1(3)=ZPIVOT(ICHORI)
C
      ICHRD=ICHORD(ICHORI)
      C2(1)=XCHORD(ICHORI)
      C2(2)=YCHORD(ICHORI)
      C2(3)=ZCHORD(ICHORI)
C
      NBC2=NSPBLC(ICHORI)
      NAC2=NSPADD(ICHORI)
C
C  ENERGY LOOP (IF ANY)
C
      PMA=-1.E30
      PMI=1.E30
      XMA=-1.E30
      XMI=1.E30

      IF (.NOT.ALLOCATED(PSIG)) THEN
        ND = 0
        IF (ANY(NCHTAL == 1)) ND = MAX(ND, NATMI)
        IF (ANY(NCHTAL == 2)) THEN
          MX_COMPO = 0
          IF (ALLOCATED(EMIS_LINES)) THEN
            DO I=1, NO_LINES
              MX_COMPO = MAX(MX_COMPO, EMIS_LINES(I)%NO_COMPO)
            END DO
          END IF
          ND = MAX(ND, MX_COMPO)
        END IF
        IF (ANY(NCHTAL == 3)) ND = MAX(ND, NPHOTI)
        IF (ANY(NCHTAL == 10)) ND = MAX(ND, NSPZ)
        ALLOCATE (PSIG(0:ND))
      END IF

      IF (NCHTAL(ICHORI).EQ.1)  NSPI=NATMI  ! post collision CX atomic species
      IF (NCHTAL(ICHORI).EQ.2)  NSPI=MX_COMPO ! use maximum number of components to spectral line emissivities (transitions) in one single LOS evaluation
      IF (NCHTAL(ICHORI).EQ.3)  NSPI=NPHOTI ! one spectrally resolved radiance per LOS and per photon species ("transition")
      IF (NCHTAL(ICHORI).EQ.5)  NSPI=10  
      IF (NCHTAL(ICHORI).EQ.10) NSPI=NSPZ   ! 3rd party specified LOS integrals.
      PSIG = 0._DP
      IFIRST=0
      DO 231 JEN=1,NCHNI
        ZE=ENERGY(JEN)
        CALL EIRENE_LININT
     .  (IFIRST,ICHORI,C1,C2,ICHRD,IPVOT,NBC2,NAC2,ZE,
     .               PSIG,TIMAX,ISP,NSPI,JEN,NCHNI)
        IFIRST=1
        IF (ISP.GT.0.AND.ISP.LE.NSPI) THEN
C  SINGLE SPECIES INDEX ISP
          BUFFER(ICHORI,JEN)=PSIG(ISP)
        ELSEIF (ISP.EQ.0) THEN
C  SUM OVER SPECIES INDEX
          ZSI=0.
          DO 239 IS=1,NSPI
            ZSI=ZSI+PSIG(IS)
239       CONTINUE
          BUFFER(ICHORI,JEN)=ZSI
          WRITE (80,'(I6,3ES12.4)') ICHORI,C2
          WRITE (80,'(6ES12.4)') PSIG(0:NSPI)
        ELSE
          WRITE (iunout,*) 'ERROR IN SUBR. SGNAL: ISP= ',ISP
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
C
C  PROCESS DATA FROM LINE INTEGRAL ROUTINES INTO REQUESTED DATA & UNITS
C
        IF (NCHTAL(ICHORI).EQ.1) THEN
C  LINE INTEGRAL: CX ATOMS/SEC/CM**2/EV/STERAD
C  THE NUMERICAL FACTOR 1./11.137 ARISES FROM A TRANSFORMATION
C  OF A MAXWELLIAN VELOCITY DISTRIBUTION TO A MAXW. ENERGY DISTR.
C  1./11.137=0.5*(1./PI)**1.5, IN SIGCX
          FUFFER(ICHORI,JEN)=BUFFER(ICHORI,JEN)/11.137
        ELSEIF ((NCHTAL(ICHORI).EQ.2).OR.(NCHTAL(ICHORI).EQ.5)) THEN
C  LINE INTEGRAL: PHOTONS/SEC/CM**2/STERAD (EMISSIVITY), JEN=1 HERE.
          FUFFER(ICHORI,JEN)=BUFFER(ICHORI,JEN)/(4.*PIA)
        ELSEIF (NCHTAL(ICHORI).EQ.3) THEN
C  LINE INTEGRAL: PHOTONS/SEC/CM**2/EV/STERAD (SPECTRAL RADIANCE)
          FUFFER(ICHORI,JEN)=BUFFER(ICHORI,JEN)/(4.*PIA)
        ELSEIF (NCHTAL(ICHORI).EQ.10) THEN
C  LINE INTEGRAL: USER SUPPLIED INTEGRAND ALONG LINE OF SIGHT
          FUFFER(ICHORI,JEN)=BUFFER(ICHORI,JEN)
        ENDIF
231   CONTINUE
C
C  ENERGY LOOP FINISHED
C
C
C  FIT A STRAIGHT LINE TO SPECTRUM, BETWEEN ESTART AND ENDFIT
C  TO BE DONE ONLY FOR NCHTAL=1
C
      IF ((NCHTAL(ICHORI).NE.1) .OR. (NCHENI == 1)) GOTO 300
 
      ESTART(ICHORI)=NSPINI(ICHORI)*TIMAX
      ENDFIT(ICHORI)=NSPEND(ICHORI)*TIMAX
      TINP(ICHORI)=TIMAX
      IF (TRCSIG) THEN
        WRITE (iunout,*) 'SLOPE FITTING RANGE: E1--E2, TIMAX'
        WRITE (iunout,*) ESTART(ICHORI),'--',ENDFIT(ICHORI),'  ',TIMAX
      ENDIF
C
C  FIND MAX. VALUE OF SIGNAL
      CALL EIRENE_MAXMN2(BUFFER,NCHOR,ICHORI,ICHORI,1,NCHNI,XMIN,XMAX)
C  SCALE RESULT
      IF (XMAX.GT.0.) GOTO 235
      CALL EIRENE_MASAGE
     .  ('NO SLOPE IN SIGNAL, BECAUSE MAX(BUFFER).LE.0   ')
      PLSPEC=.FALSE.
      GOTO 300
235   ZSCALE=1./XMAX
      DO 233 I=1,NCHNI
        BUFFER(ICHORI,I)=BUFFER(ICHORI,I)*ZSCALE
        ZZ=MAX(1.E-10_DP,BUFFER(ICHORI,I))
233     BUFFER(ICHORI,I)=LOG(ZZ)
C
C  CURVE FITTING
C
C  MIN. ENERGY FOR FIT
      ZE1=ESTART(ICHORI)
C  MAX. ENERGY
      ZE2=ENDFIT(ICHORI)
C
C  FIND ELEMENTS
      DO 241 JEN=2,NCHNI
        I1=JEN
        IF (ENERGY(I1).GE.ZE1) GO TO 242
241   CONTINUE
242   CONTINUE
C
      DO 243 JEN=I1,NCHNI
        I2=JEN
        IF (ENERGY(I2).GE.ZE2) GO TO 244
243   CONTINUE
244   CONTINUE
C
C   NUMBER OF POINTS FOR FITTING
      IN=I2-I1+1
      IF (IN.LT.2) THEN
        CALL EIRENE_MASAGE
     .  ('WRONG ENERGY RANGE FOR CURVE FITTING IN SIGNAL ')
        CALL EIRENE_MASJ1 ('CHORD.NO. I=           ',ICHORI)
        WRITE (iunout,*) 'I1,I2,IN ',I1,I2,IN
        TILINE(ICHORI)=0.
      ELSE
C  FITTED RESULT
        STEIG=EIRENE_SLOPE(IN,ICHORI,I1,ENERGY,BUFFER,NCHOR,NCHNI)
        IF (STEIG.GE.0.) THEN
          TILINE(ICHORI)=0.
        ELSE
          TILINE(ICHORI)=-1./STEIG
        ENDIF
      ENDIF
C
      IF (TRCSIG) THEN
        WRITE (iunout,*) 'ICHORI,TILINE(ICHORI) ',ICHORI,TILINE(ICHORI)
      ENDIF
 
300   CONTINUE
C
      IF (NCHTAL(ICHORI).EQ.3) THEN
        DEALLOCATE(RECADD)
        DEALLOCATE(INTADD)
      ENDIF
      RETURN
      END
