cdr  aug 19 learc not needed
c
C ===== SOURCE: plausr.f
c
c modified: s.wiesen@fz-juelich.de
c
c reads target plasma information from casename.zplasma
c sets step-functions
c
      SUBROUTINE EIRENE_PLAUSR
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CSTEP
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      USE EIRMOD_CINIT
      USE EIRMOD_COMSOU
      USE EIRMOD_CTRIG
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CCONA
      USE EIRMOD_LEARC1, ONLY: EIRENE_LEARC1

      IMPLICIT NONE

      REAL(DP) :: FACTOR, EIRENE_STEP
      INTEGER :: NLINES, ITRI, ISIDE, I, NBIN, ISTRA, ISRFS, ISOR, JJJ,
     .           ITEC1, ITEC2, ITEC3, ISTEP, INDSRF, IS1, IERROR, IPLS,
     .           ITARG, IG, ISP
      INTEGER :: EIRENE_IDEZ
      INTEGER, ALLOCATABLE :: KSTEP(:), INOSRC(:), IPLAN(:), IPLEN(:)
      REAL(DP) :: FLX, TE, TI, DE, MC, FE, FI, FSH, VP, FEL, DUM
      REAL(DP) :: DELR, FL, MCC, FFEL, CS, vx,vy,vz,di,usrval
      real(dp) :: xref, yref, bzref, facbz, x, y, rad, bx, by, bz, bf,
     .            errbx, errby, errbz, errbf, bxmax, bymax, bzmax,
     .            bfmax, dfdx, dfdy, dfdz, xref2, yref2, bzref2, facbz2,
     .            xref3, yref3, bzref3, facbz3
      integer :: nref, icell, nplcll, ipolg, nref2, nref3
      CHARACTER(256) :: line,sstr,filename
      character(2) :: cstr2
      integer, allocatable :: itritmp(:),isidetmp(:)
      integer, allocatable, save :: inmass(:), inchar(:), inchrg(:)
      real(dp), allocatable :: tetmp(:),fetmp(:),fshtmp(:),
     .              flxtmp(:,:),titmp(:,:),ditmp(:,:),
     .              vxtmp(:,:),vytmp(:,:),vztmp(:,:),
     .              fitmp(:,:),feltmp(:,:),
     .              vptmp(:,:),mctmp(:,:),
     .              usrvaltmp(:,:)
      REAL(DP), ALLOCATABLE :: PS(:), PS_CORNER(:), copy(:)

      integer, parameter :: fp=31
      integer :: ll,ier,j,ind, jj

      CALL EIRENE_ALLOC_CSTEP
      ALLOCATE (KSTEP(NSTEP))
      ALLOCATE (INOSRC(NSTEP))
      ALLOCATE (IPLAN(NSTEP))
      ALLOCATE (IPLEN(NSTEP))
      allocate (inmass(npls))
      allocate (inchar(npls))
      allocate (inchrg(npls))
      KSTEP = 0
      INOSRC = 0
      IPLAN = 0
      IPLEN = 0
      inmass = 0
      inchar = 0
      inchrg = 0

      ll=len_trim(casename)
      filename=casename(1:ll) // '.plasma'
      open (unit=fp+ifoff,file=filename,access='sequential',
     .     form='formatted')

c firstly, read misc target data
      write(sstr,'(a20)') '*** MISC TARGET DATA'
      CALL EIRENE_locstr_usr(fp+ifoff,sstr,ier)
      if(ier /=0) then
         write(iunout,*) 'PLAUSR: ',trim(sstr),' not found'
         close(fp+ifoff)
         CALL EIRENE_exit_own(1)
      endif
      read(fp+ifoff,'(a)') line
      do while (line(1:1) == '*')
         read(fp+ifoff,'(a)') line
      enddo
      READ (line,*) NLINES
      if(nlines <= 0) then
         call eirene_leer(2)
         write(iunout,*) ' WARNING!!! '
         write(iunout,*) ' NO TARGET DATA FOUND IN FILE ',trim(filename)
         write(iunout,*) ' TRY TO CONTINUE WITHOUT THIS DATA'
         call eirene_leer(2)
         close(fp+ifoff)
         DEALLOCATE (KSTEP)
         DEALLOCATE (INOSRC)
         DEALLOCATE (IPLAN)
         DEALLOCATE (IPLEN)
         deallocate (inmass)
         deallocate (inchar)
         deallocate (inchrg)
         return
      endif

      allocate(itritmp(nlines))
      allocate(isidetmp(nlines))
      allocate(tetmp(nlines))
      allocate(fetmp(nlines))
      allocate(fshtmp(nlines))
      do j=1,nlines
!pb         read(fp+ifoff,'(3i7,3(1x,e14.7))')
         read(fp+ifoff,*)
     .        jj, itritmp(j), isidetmp(j),
     .        tetmp(j), fetmp(j), fshtmp(j)
      enddo

      allocate(flxtmp(nplsi,nlines))
      allocate(titmp(nplsi,nlines))
      allocate(ditmp(nplsi,nlines))
      allocate(vxtmp(nplsi,nlines))
      allocate(vytmp(nplsi,nlines))
      allocate(vztmp(nplsi,nlines))
      allocate(fitmp(nplsi,nlines))
      allocate(feltmp(nplsi,nlines))
      allocate(vptmp(nplsi,nlines))
      allocate(mctmp(nplsi,nlines))
      allocate(usrvaltmp(nplsi,nlines))

      flxtmp = 0._dp
      titmp = 0._dp
      ditmp = 0._dp
      vxtmp = 0._dp
      vytmp = 0._dp
      vztmp = 0._dp
      fitmp = 0._dp
      feltmp = 0._dp
      vptmp = 0._dp
      mctmp = 0._dp
      usrvaltmp = 0._dp

      DO IPLS = 1, NPLSI
c     search target tag in .plasma file
         write(cstr2,'(i2.2)') ipls
         write(sstr,'(a23)')
     .          '*** ION #'//cstr2//' TARGET DATA'
         CALL EIRENE_locstr_usr(fp+ifoff,sstr,ier)
         if(ier /=0) then
            write(iunout,*) 'PLAUSR: tag ',trim(sstr),' not found'
!            close(fp+ifoff)
!            call EIRENE_exit_own(1)
            write(iunout,*) 'using zero plasma and fluxes'
            write(iunout,*)
            cycle
         endif
         read(fp+ifoff,'(a)') line
         do while (line(1:1) == '*')
            read(fp+ifoff,'(a)') line
         enddo

         read(line,*) inmass(ipls),inchar(ipls),inchrg(ipls)

         if ((nmassp(ipls) /= inmass(ipls)) .or.
     .       (ncharp(ipls) /= inchar(ipls)) .or.
     .       (nchrgp(ipls) /= inchrg(ipls))) then
            write(iunout,*) 'ERROR in plasma file found in PLAUSR'
            write(iunout,*) 'mass or charge of ion target data '
            write(iunout,*) 'does not match species data from '
            write(iunout,*) 'Eirene input file '
            write(iunout,*) 'species index ', ipls
            write(iunout,*) 'nmassp, ncharp, nchrgp ',
     .                       nmassp(ipls),ncharp(ipls),nchrgp(ipls)
            write(iunout,*) 'inmass, inchar, inchrg ',
     .                       inmass(ipls),inchar(ipls),inchrg(ipls)
            close(fp+ifoff)
            call EIRENE_exit_own(1)
         end if

c     read step functions
         READ (fp+ifoff,*) NLINES
         DO I=1, NLINES
!pb               read(fp+ifoff,'(i7,11(1x,e14.7))')
            read(fp+ifoff,*)
     .           ind,
     .           flxtmp(ipls,i),titmp(ipls,i),ditmp(ipls,i),
     .           vxtmp(ipls,i),vytmp(ipls,i),vztmp(ipls,i),
     .           fitmp(ipls,i),feltmp(ipls,i),
     .           vptmp(ipls,i),mctmp(ipls,i),
     .           usrvaltmp(ipls,i)
         END DO  ! i
      END DO ! ipls

c now read species-dependent target data
      DO ISTRA=1,NSTRAI
         IF (.NOT.NLSRF(ISTRA)) CYCLE
         DO ISRFS=1,NSRFSI(ISTRA)

c     get step function index ISTEP
            ISOR=SORLIM(ISRFS,ISTRA)
            ITEC1=EIRENE_IDEZ(ISOR,1,4)
            ITEC2=EIRENE_IDEZ(ISOR,2,4)
            ITEC3=EIRENE_IDEZ(ISOR,3,4)
            IF ((ITEC1 /= 4).AND.(ITEC2 /= 4).AND.(ITEC3 /= 4)) CYCLE
            ISTEP=SORIND(ISRFS,ISTRA)
            IF (ISTEP.EQ.0) THEN
               WRITE (iunout,*) 'ERROR IN PRIMARY SOURCE DATA '
               WRITE (iunout,*) 'STEPFUNCTION REQUESTED FOR SOURCE ',
     .                          'SURFACE '
               WRITE (iunout,*) 'NO. ',INSOR(ISRFS,ISTRA),
     .                          ' BUT SORIND.EQ.0.'
               CALL EIRENE_EXIT_own(1)
            ELSEIF (ISTEP.GT.NSTEP) THEN
               CALL EIRENE_MASPRM
     .              ('NSTEP',5,NSTEP,'ISTEP',5,ISTEP,IERROR)
               CALL EIRENE_EXIT_own(1)
            ENDIF

c     get surface index INDSRF
            INDSRF=INSOR(ISRFS,ISTRA)
            IF (INDSRF < 0) INDSRF=NLIM+ABS(INDSRF)

c     get species index/indices IPLAN(ISTEP) --> IPLEN(ISTEP)
            IF (NSPEZ(ISTRA) <= 0) THEN
               IPLAN(ISTEP)=1
               IPLEN(ISTEP)=NPLSI
               ipls = 0
            ELSE
               IPLAN(ISTEP)=NSPEZ(ISTRA)
               IPLEN(ISTEP)=NSPEZ(ISTRA)
               ipls=nspez(istra)
            END IF

            DO I=1, NLINES

               itri = itritmp(i)
               iside= isidetmp(i)
               te   = tetmp(i)
               fe   = fetmp(i)
               fsh  = fshtmp(i)

               IF (INMTI(ISIDE,ITRI) == INDSRF) THEN
                  IF (KSTEP(ISTEP) == 0) RRSTEP(ISTEP,1) = 0._DP
                  KSTEP(ISTEP) = KSTEP(ISTEP) + 1
                  INOSRC(ISTEP) = ISTRA
                  IS1 = ISIDE + 1
                  IF (IS1.GT.3) IS1=1
                  IRSTEP(ISTEP,KSTEP(ISTEP))=ITRI
                  IPSTEP(ISTEP,KSTEP(ISTEP))=ISIDE
                  ITSTEP(ISTEP,KSTEP(ISTEP))=1
                  IASTEP(ISTEP,KSTEP(ISTEP))=0
                  IBSTEP(ISTEP,KSTEP(ISTEP))=1
                  DELR =  SQRT(
     .                 (XTRIAN(NECKE(ISIDE,ITRI))
     .                 -XTRIAN(NECKE(IS1,ITRI)))**2+
     .                 (YTRIAN(NECKE(ISIDE,ITRI))
     .                 -YTRIAN(NECKE(IS1,ITRI)))**2)
                  RRSTEP(ISTEP,KSTEP(ISTEP)+1)=
     .                 RRSTEP(ISTEP,KSTEP(ISTEP)) + DELR
                  TESTEP(ISTEP,KSTEP(ISTEP)) = TE
                  FESTEP(ISTEP,KSTEP(ISTEP)) = FE
C     IF NO SHEATH POTENTIAL SPECIFIED, DERIVE IT FROM ELECTRON ENERGY
C     FLUX BY SUBTRACTING THE KINETIC CONTRIBUTION 2.0*TE
                  IF (FSH.EQ.0..AND.FE.GE.2.0) FSH=FE-2.0
                  SHSTEP(ISTEP,KSTEP(ISTEP)) = FSH
                  DO IPLS=IPLAN(ISTEP), IPLEN(ISTEP)

                     flx = flxtmp(ipls,i)
                     ti = titmp(ipls,i)
                     di = ditmp(ipls,i)
                     vx = vxtmp(ipls,i)
                     vy = vytmp(ipls,i)
                     vz = vztmp(ipls,i)
                     fi = fitmp(ipls,i)
                     fel = feltmp(ipls,i)
                     vp = vptmp(ipls,i)
                     mc = mctmp(ipls,i)
                     usrval =usrvaltmp(ipls,i)

                     TISTEP(IPLS,ISTEP,KSTEP(ISTEP)) = TI ! eV
                     DISTEP(IPLS,ISTEP,KSTEP(ISTEP)) = DI ! 1/cm**3
                     FISTEP(IPLS,ISTEP,KSTEP(ISTEP)) = FI !   1
                     VPSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = abs(VP) ! cm/s
C     VP OVERRULES MC, IF VP IS GIVEN and MC=0
                     MCC=0.0
                     IF (VP.NE.0.0) THEN
                        CS=SQRT((TI+TE)/RMASSP(IPLS))*CVEL2A
                        IF (CS.GT.0.) MCC=VP/CS
                     ENDIF
                     jjj=kstep(istep)
                     IF (MC.EQ.0.) MC=MCC
!pb                     MCSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = abs(MC) ! 1
C     THIS NEXT VECTOR IS V-PARALLEL, IN CARTESIAN COORDINATES
                     IF (ABS(VX) > EPS10) THEN
                       VXSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = VX
                       VYSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = VY
                       VZSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = VZ
                     ELSE
                       VXSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = VXIN(IPLS,ITRI)
                       VYSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = VYIN(IPLS,ITRI)
                       VZSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = VZIN(IPLS,ITRI)
                     END IF
                     FLSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = ABS(FLX)/DELR
C     IF NO ION KINETIC ENERGY FLUX IS SPECIFIED, DERIVE IT FROM  FI,MC,VP
                     ffel=0.
                     IF (FI.gt.0..or.mc.gt.0.) then
                        ffel=(FI*TI+0.5*MC*MC*(TE+TI))*abs(FLX)
                     endif
                     IF (FEL.EQ.0.) FEL=FFEL
                     ELSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = FEL
csw
!pb                     usrstep(ipls,istep,kstep(istep)) = usrval
csw
c                 enddo ipls,iplan
                  ENDDO

c              endif inmti
               ENDIF

c           enddo nlines
            ENDDO

c        enddo isrfs
         ENDDO

         IF (KSTEP(ISTEP) > 0) THEN
            NBIN=KSTEP(ISTEP)+1
            FL=EIRENE_STEP(IPLAN(ISTEP),IPLEN(ISTEP),NBIN,ISTEP)
            FLUX(INOSRC(ISTEP))=FL
         END IF

         CALL EIRENE_LEER(1)
         WRITE (IUNOUT,*) 'TARGET DATA: TARGET NO. ITARG=ISTRA= ',ISTRA
         WRITE (IUNOUT,*)
     .' IG,  ARC,     P-FLUX,   E-FLUX,     TE,       TI,    SHEATH/TE'
         ITARG = ISTRA
         DO 6100 IG=1,KSTEP(ISTEP)
           WRITE (IUNOUT,'(1X,I3,1P,6E11.3)')
     .             IG,RRSTEP(ITARG,IG),FLSTEP(0,ITARG,IG),
     .             ELSTEP(0,ITARG,IG),
     .             TESTEP(ITARG,IG),TISTEP(1,ITARG,IG),
     .             SHSTEP(ITARG,IG)
 6100    CONTINUE
         WRITE (IUNOUT,'(1X,I3,1P,1E11.3)') KSTEP(ISTEP)+1,
     .                                      RRSTEP(ITARG,KSTEP(ISTEP)+1)
C
         WRITE (IUNOUT,*) 'PARTICLE FLUX(IPLS), IPLS=1,NPLSI '
         WRITE (IUNOUT,'(1X,1P,6E12.4)') (FLTOT(ISP,ITARG),ISP=1,NPLSI)
         CALL EIRENE_LEER(2)
C

c     enddo istra
      ENDDO

      DEALLOCATE (KSTEP)
      DEALLOCATE (INOSRC)
      DEALLOCATE (IPLAN)
      DEALLOCATE (IPLEN)
      deallocate(itritmp)
      deallocate(isidetmp)
      deallocate(tetmp)
      deallocate(fetmp)
      deallocate(fshtmp)

      deallocate(flxtmp)
      deallocate(titmp)
      deallocate(ditmp)
      deallocate(vxtmp)
      deallocate(vytmp)
      deallocate(vztmp)
      deallocate(fitmp)
      deallocate(feltmp)
      deallocate(vptmp)
      deallocate(mctmp)
      deallocate(usrvaltmp)

      close(fp+ifoff)


   99 RETURN
      END SUBROUTINE EIRENE_PLAUSR
