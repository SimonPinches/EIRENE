      subroutine eirene_emissivity(ist)

      use eirmod_precision
      use eirmod_parmmod
      use eirmod_comsig
      use eirmod_ccona
      use eirmod_comusr
      use eirmod_comsou
      use eirmod_cgeom
      use eirmod_cgrid
      use eirmod_ctext
      use eirmod_coutau
      use eirmod_cspei
      use eirmod_ctrcei
      USE EIRMOD_CESTIM
      USE EIRMOD_CSDVI
      USE EIRMOD_CSDVI_BGK
      USE EIRMOD_CSDVI_COP
      USE EIRMOD_COMPRT

      implicit none

      integer, intent(in) :: ist
      integer :: l1, l2, i, j, k, iads, iadv, isp, itp, iratio, irc,
     .           irc_rat, ncelc
      real(dp) :: density, sigadd, add, ratio, powalf, powalfs, fac,
     .            ry, dee00, facte, DE, TE, TEF, DEF, rate, erate,
     .            EIRENE_OTHER_RATE_COEFF
      REAL(DP) :: DUMMY(NRTAL)
      REAL(DP), ALLOCATABLE :: OUTAU(:)
      logical :: lwrite
      CHARACTER(6) :: CISTRA

      CALL EIRENE_LEER(2)
      CALL EIRENE_FTCRI(IST,CISTRA)
      IF (IST.GT.0) CALL EIRENE_MASBOX
     .   ('SUBR. EMISSIVITY CALLED, FOR STRATUM NO. '//CISTRA)
      IF (IST.EQ.0) CALL EIRENE_MASBOX
     .   ('SUBR. EMISSIVITY CALLED, FOR SUM OVER STRATA')
      CALL EIRENE_LEER(1)
      WRITE (iunout,*) ' AFTER INTEGRATION OVER COMPUTATIONAL DOMAIN'

      do i = 1, no_lines
        WRITE (iunout,*) ' FLUX (AMP) AND POWER (WATT) BY ' //
     .                 EMIS_LINES(I)%LINE_NAME // ':'

        l1 = emis_lines(i)%l1
        l2 = emis_lines(i)%l2
        fac = emis_lines(i)%fac
        iads = emis_lines(i)%iadv_total

C  ENERGY FACTOR FOR POWER LOSS (W)
        RY=13.605
!pb     DEE00=RY*(1./(2.*2.)-1./(3.*3.))
        DEE00=RY*(1./REAL(L2*L2,KIND(1._DP))
     .        -1./REAL(L1*L1,KIND(1._DP)))
        FACTE=DEE00*ELCHA

        powalfs = 0._dp
        do j = 1, emis_lines(i)%no_compo
          iadv = emis_lines(i)%compo(j)%iadv
          sigadd = 0._dp
          powalf = 0._dp
            
          do k = 1, emis_lines(i)%compo(j)%no_contrib
            isp = emis_lines(i)%compo(j)%contrib(k)%isp
            itp = emis_lines(i)%compo(j)%contrib(k)%itp
            iratio = emis_lines(i)%compo(j)%contrib(k)%iratio
            irc = emis_lines(i)%compo(j)%contrib(k)%irc
            irc_rat = emis_lines(i)%compo(j)%contrib(k)%irc_rat

            lwrite = .true.

            DO NCELL=1,NSBOX
C
C  LOCAL BACKGROUND DATA ARE IN CELL NCELL
C  LOCAL TEST PARTICLE DATA ARE IN (PERHAPS COARSER) SCORING CELL NCELC
C  ACCUMULATE THE EMISSIVITIES ALSO ON THE "SCORING" GRID.
C
              NCELC=NCLTAL(NCELL)
C
              IF (NSTGRD(NCELL) > 0) CYCLE
              IF (LGVAC(NCELL,NPLS+1)) CYCLE

              TE=TEIN(NCELL)
              DE=DEIN(NCELL)
              
              DEF=LOG(DE*1.D-8)
              TEF=LOG(TE)

              select case (itp)
                case (0)
                  density = pdenph(isp,ncelc)
                case (1)
                  density = pdena(isp,ncelc)
                case (2)
                  density = pdenm(isp,ncelc)
                case (3)
                  density = pdeni(isp,ncelc)
                case (4)
                  density = diin(isp,ncell)
                case default
                  density = 0._dp
                  if (lwrite) then
                    write (iunout,*) ' ERROR IN EMISSIVITY' 
                    write (iunout,*) 
     .                ' WRONG PARTICLE TYPE SPECIFIED FOR'
                    write (iunout,*) ' line ',i,emis_lines(i)%line_name
                    write (iunout,*) ' component ',j,
     .                 emis_lines(i)%compo(j)%compo_name
                    write (iunout,*) ' contribution ',k
                    lwrite = .false.
                  end if
              end select
             
              rate = EIRENE_OTHER_RATE_COEFF(IRC,TEF,DEF,.TRUE.,0,ERATE)
              add = rate*density

              if (iratio > 0) then 
                ratio = EIRENE_OTHER_RATE_COEFF(IRC_RAT,TEF,DEF,
     .                                           .TRUE.,0,ERATE)
                add = add*ratio
              end if

              sigadd = sigadd + add * fac * vol(ncell)
              addv(iadv,ncelc) = addv(iadv,ncelc) + sigadd
              addv(iads,ncelc) = addv(iads,ncelc) + sigadd

              powalf = powalf + sigadd
            end do               ! ncell
      
          end do ! k contributions of component j of line i

          addv(iadv,1:nsbox_tal) = addv(iadv,1:nsbox_tal) 
     .                             / voltal(1:nsbox_tal)

          powalf = powalf * facte
          powalfs = powalfs + powalf

          WRITE (iunout,'(A,2ES16.7)') ' COUPL. TO ' // 
     .                     TRIM(EMIS_LINES(I)%COMPO(J)%COMPO_NAME)
     .                    ,POWALF/FACTE*ELCHA,POWALF

          DUMMY(1:NSBOX_TAL) = ADDV(IADV,1:NSBOX_TAL)
          CALL EIRENE_INTTAL
     .         (DUMMY,VOLTAL,1,1,NSBOX_TAL,ADDVI(IADV,IST),
     .          NR1TAL,NP2TAL,NT3TAL,NBMLT)
          ADDV(IADV,1:NSBOX_TAL) = DUMMY(1:NSBOX_TAL)

          TXTTAL(IADV,NTALA) =REPEAT(' ',72)
          TXTTAL(IADV,NTALA) =TRIM(EMIS_LINES(I)%LINE_NAME) // ', ' //
     .                    ' SOURCE RATE '
          TXTSPC(IADV,NTALA) =TRIM(EMIS_LINES(I)%COMPO(J)%COMPO_NAME)
          TXTUNT(IADV,NTALA) ='PHOTONS/S/CM**3         '

        end do ! j components of line i

        addv(iads,1:nsbox_tal) = addv(iads,1:nsbox_tal) * fac 
     .                           / voltal(1:nsbox_tal)

        WRITE (iunout,'(A,2ES16.7)') 
     ,                  ' TOTAL FLUX (AMP) AND POWER (WATT) ' 
     .                  ,POWALFS/FACTE*ELCHA,POWALFS
        CALL EIRENE_LEER(2)

        DUMMY(1:NSBOX_TAL) = ADDV(IADS,1:NSBOX_TAL)
        CALL EIRENE_INTTAL
     .       (DUMMY,VOLTAL,1,1,NSBOX_TAL,ADDVI(IADS,IST),
     .        NR1TAL,NP2TAL,NT3TAL,NBMLT)
        ADDV(IADS,1:NSBOX_TAL) = DUMMY(1:NSBOX_TAL)

        TXTTAL(IADS,NTALA) =TXTTAL(IADV,NTALA)
        TXTSPC(IADS,NTALA) ='SUM_OVER_ALL            '
        TXTUNT(IADS,NTALA) ='PHOTONS/S/CM**3         '

      end do ! line i

C
C  WRITE ON STREAM 11 DATA FOR STRATUM NO. IST
      IF (NFILEN.EQ.1.OR.NFILEN.EQ.2) THEN
        IESTR=IST
        CALL EIRENE_WRSTRT(IST,NSTRAI,NESTM1,NESTM2,NADSPC,
     .              ESTIMV,ESTIMS,ESTIML,
     .              NSDVI1,SDVI1,NSDVI2,SDVI2,
     .              NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .              NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .              NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .              NSIGI_SPC,TRCFLE)
C
        IRC=2
        ALLOCATE (OUTAU(NOUTAU))
        CALL EIRENE_WRITE_COUTAU (OUTAU, IUNOUT)
        WRITE (11+ifoff,REC=IRC) OUTAU
        DEALLOCATE (OUTAU)
        IF (TRCFLE)   WRITE (iunout,*) 'WRITE 11  IRC= ',IRC

C  WRITE ON STREAM 11 ONLY DATA FOR SUM OVER STRATA
      ELSEIF ((NFILEN.EQ.6.OR.NFILEN.EQ.7).AND.IST.EQ.0) THEN
        IESTR=IST
        CALL EIRENE_WRSTRT(IST,NSTRAI,NESTM1,NESTM2,NADSPC,
     .              ESTIMV,ESTIMS,ESTIML,
     .              NSDVI1,SDVI1,NSDVI2,SDVI2,
     .              NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .              NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .              NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .              NSIGI_SPC,TRCFLE)
C
        IRC=2
        ALLOCATE (OUTAU(NOUTAU))
        CALL EIRENE_WRITE_COUTAU (OUTAU, IUNOUT)
        WRITE (11+ifoff,REC=IRC) OUTAU
        DEALLOCATE (OUTAU)
        IF (TRCFLE)   WRITE (iunout,*) 'WRITE 11  IRC= ',IRC
      ENDIF
C
      RETURN

      end subroutine eirene_emissivity
