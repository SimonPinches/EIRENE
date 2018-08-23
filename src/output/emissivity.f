cdr  comments

cdr  may 18: some comments tried......NOT FINISHED



      subroutine eirene_emissivity(istr, lstart, lend)

cdr  probably something to fill ADDV tallies with emissivities, stratum ISTR
cdr  for lines lstart to lend ?? Contained parts of old routines Ba_alpha,....,Ly-Beta.
cdr  write the newly defined tallies ADDV onto stream fort.11, stratum ISTR



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

      integer, intent(in) :: istr, lstart, lend
      integer :: i, j, k, iline, jcomp, kcontr,
     .           iads, iadv, isp(3), itp(3), iratio, irc,
     .           irc_rat(2), ncelc, ndens, idens
      real(dp) :: density(3), sigadd, add, ratio, powalf, powalfs, 
     .            einstein, trans_en, DE, TE, TEF, DEF, popcf, 
     .            EIRENE_OTHER_RATE_COEFF, 
     .            ratio2
      REAL(DP) :: DUMMY(NRTAL)
      REAL(DP), ALLOCATABLE :: OUTAU(:)
      logical :: lwrite
      CHARACTER(6) :: CISTRA

      CALL EIRENE_LEER(2)
      CALL EIRENE_FTCRI(ISTR,CISTRA)
      IF (ISTR.GT.0) CALL EIRENE_MASBOX
     .   ('SUBR. EMISSIVITY CALLED, FOR STRATUM NO. '//CISTRA)
      IF (ISTR.EQ.0) CALL EIRENE_MASBOX
     .   ('SUBR. EMISSIVITY CALLED, FOR SUM OVER STRATA')
      CALL EIRENE_LEER(1)

      WRITE (iunout,*) ' AFTER INTEGRATION OVER COMPUTATIONAL DOMAIN'

      do i = lstart, lend
        ILINE=I
        WRITE (iunout,*) 'LINE no. ',ILINE       
        WRITE (iunout,*) ' FLUX (AMP) AND POWER (WATT) BY ' //
     .                 EMIS_LINES(I)%LINE_NAME // ':'
        write (iunout,'(A,ES12.4)') 'EINSTEIN COEFFICIENT',
     .                               emis_lines(i)%einstein
        write (iunout,'(A,ES12.4/1x)') 'TRANSITION ENERGY   ',
     .                               emis_lines(i)%trans_en

        einstein = emis_lines(i)%einstein
C  ENERGY FACTOR FOR POWER LOSS (W)
        trans_en = emis_lines(i)%trans_en * elcha

cdr initialize sum over components
        iads = emis_lines(i)%iadv_total
        addv(iads,:) = 0._dp
        powalfs = 0._dp

cdr run over components
        do j = 1, emis_lines(i)%num_compo
          iadv = emis_lines(i)%compo(j)%iadv
          addv(iadv,:) = 0._dp
          sigadd = 0._dp
          powalf = 0._dp
            
          do k = 1, emis_lines(i)%compo(j)%num_contrib
            isp = emis_lines(i)%compo(j)%contrib(k)%isp
            itp = emis_lines(i)%compo(j)%contrib(k)%itp
            iratio = emis_lines(i)%compo(j)%contrib(k)%iratio
            irc = emis_lines(i)%compo(j)%contrib(k)%irc
            irc_rat = emis_lines(i)%compo(j)%contrib(k)%irc_rat

            ndens = count(itp >= 0)
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
              
              DEF=LOG(DE)
              TEF=max(-2.30,LOG(TE)) ! cut off at 0.1 eV

              do idens = 1, ndens
                select case (itp(idens))
                  case (0)
                    density(idens) = pdenph(isp(idens),ncelc)
                  case (1)
                    density(idens) = pdena(isp(idens),ncelc)
                  case (2)
                    density(idens) = pdenm(isp(idens),ncelc)
                  case (3)
                    density(idens) = pdeni(isp(idens),ncelc)
                  case (4)
                    density(idens) = diin(isp(idens),ncell)
                  case (5)
                    density(idens) = dein(ncell)
                  case default
                    density(idens) = 0._dp
                    if (lwrite) then
                      write (iunout,*) ' ERROR IN EMISSIVITY' 
                      write (iunout,*) 
     .                  ' WRONG PARTICLE TYPE SPECIFIED FOR'
                      write (iunout,*) ' line ',i,
     .                   emis_lines(i)%line_name
                      write (iunout,*) ' component ',j,
     .                   emis_lines(i)%compo(j)%compo_name
                      write (iunout,*) ' contribution ',k
                      lwrite = .false.
                    end if
                end select
              end do
  
c  density is the "true" parent density         
c  density(1) is taken as "intermediate" parent density. Fetch reduced population coefficent
c  and density ratios ratio="density"/"density(1)" will be applied below, 
c  to turn density(1)it into "density"
              popcf= EIRENE_OTHER_RATE_COEFF(IRC,NCELL,TEF,DEF,.TRUE.,1)
              add = popcf*density(1)

c  density ratio, if true parent density is not available (or in QSS mode)
c  then: ratio converts from density(1) to density 
              if (iratio > 0) then 

                ratio = EIRENE_OTHER_RATE_COEFF(IRC_RAT(1),NCELL,
     .                                          TEF,DEF,.TRUE.,1)
                add = add*ratio
c  second conversion to yet another parent density
c  e.g: density    = H3+.   = [H2+] * [H2/ne] *ratio2
c       density(1) = H2
c       ratio      = H2+/H2(Te,ne) (CR equilibrium)
                if (iratio == 2) then
                  ratio2 = EIRENE_OTHER_RATE_COEFF(IRC_RAT(2),NCELL,
     .                                             TEF,DEF,.TRUE.,1)
                  add = add * density(2) / density(3) *ratio2
                end if
              end if

cdr so far: add is scored on the fine grid cell "ncell".
cdr         add volume weighted contribution to coarse cell "ncelc"
              sigadd = add * einstein * vol(ncell)

              addv(iadv,ncelc) = addv(iadv,ncelc) + sigadd
              addv(iads,ncelc) = addv(iads,ncelc) + sigadd

              powalf = powalf + sigadd
            end do               ! ncell
      
          end do ! k contributions (summed) of component j of line iline

cdr addv was volume weighted (extensive) sum. now divide by coarse cell volume
cdr      to turn it into an intensive score:  [...] per cm**3  
          addv(iadv,1:nsbox_tal) = addv(iadv,1:nsbox_tal) 
     .                             / voltal(1:nsbox_tal)

          powalf = powalf * trans_en
          powalfs = powalfs + powalf

          WRITE (iunout,'(A50,2ES16.7)') ' COUPL. TO ' // 
     .                     TRIM(EMIS_LINES(I)%COMPO(J)%COMPO_NAME)
     .                    ,POWALF/TRANS_EN*ELCHA,POWALF

          DUMMY(1:NSBOX_TAL) = ADDV(IADV,1:NSBOX_TAL)
          CALL EIRENE_INTTAL
     .         (DUMMY,VOLTAL,1,1,NSBOX_TAL,ADDVI(IADV,ISTR),
     .          NR1TAL,NP2TAL,NT3TAL,NBMLT)
          ADDV(IADV,1:NSBOX_TAL) = DUMMY(1:NSBOX_TAL)

          TXTTAL(IADV,NTALA) =REPEAT(' ',72)
          TXTTAL(IADV,NTALA) =TRIM(EMIS_LINES(I)%LINE_NAME) // ', ' //
     .                    ' SOURCE RATE '
          TXTSPC(IADV,NTALA) =TRIM(EMIS_LINES(I)%COMPO(J)%COMPO_NAME)
          TXTUNT(IADV,NTALA) ='PHOTONS/S/CM**3         '

          WRITE (iunout,*) ' TALLY ADDV(IADV) prepared. IADV=',IADV 

        end do ! j components of line ILINE are done

cdr  now sum over compontents: on tally ADDV(IADS)
        call eirene_leer(1)
        addv(iads,1:nsbox_tal) = addv(iads,1:nsbox_tal) 
     .                           / voltal(1:nsbox_tal)

        WRITE (iunout,'(A50,2ES16.7)') 
     ,                  ' TOTAL FLUX (AMP) AND POWER (WATT) ' 
     .                  ,POWALFS/TRANS_EN*ELCHA,POWALFS

        DUMMY(1:NSBOX_TAL) = ADDV(IADS,1:NSBOX_TAL)
        CALL EIRENE_INTTAL
     .       (DUMMY,VOLTAL,1,1,NSBOX_TAL,ADDVI(IADS,ISTR),
     .        NR1TAL,NP2TAL,NT3TAL,NBMLT)
        ADDV(IADS,1:NSBOX_TAL) = DUMMY(1:NSBOX_TAL)

        TXTTAL(IADS,NTALA) =REPEAT(' ',72)
        TXTTAL(IADS,NTALA) ='SUM OVER COMPONENTS  '
        TXTSPC(IADS,NTALA) ='  '
        TXTUNT(IADS,NTALA) ='PHOTONS/S/CM**3         '
        WRITE (iunout,*) ' TALLY ADDV(IADV) prepared. IADV=',IADS 

        CALL EIRENE_LEER(2)

      end do ! line no. ILINE

C
C  WRITE ON STREAM 11 DATA FOR STRATUM NO. ISTR
      IF (NFILEN.EQ.1.OR.NFILEN.EQ.2) THEN
        IESTR=ISTR
        CALL EIRENE_WRSTRT(ISTR,NSTRAI,NESTM1,NESTM2,NADSPC,
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
      ELSEIF ((NFILEN.EQ.6.OR.NFILEN.EQ.7).AND.ISTR.EQ.0) THEN
        IESTR=ISTR
        CALL EIRENE_WRSTRT(ISTR,NSTRAI,NESTM1,NESTM2,NADSPC,
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
