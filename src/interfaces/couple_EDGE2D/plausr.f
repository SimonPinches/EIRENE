c
c modified: s.wiesen@fz-juelich.de
c
C+---------------------------------------------------------------+
C| Purpose:                                                      |
C| --------                                                      |
C| reads target plasma information from casename.zplasma         |
C| sets step-functions                                           |
C| reads also neutral particle fluxes from 'eirene.chemFluxDep'  |
C| and sets FLXOUT for fluxdependency of chemical sputtering.    |
C+---------------------------------------------------------------+
C| Modifications:                                                |
C| --------------                                                |
C| 16/07/2010   D.Harting    Added reading of neutral fluxes from|
C|                           file eirene.chemFluxDep. Added also |
C|                           two variables to eirene_user        |
C|                           namelist.                           |
C| 23/11/2010   D.Harting    If EDGE2D is used with density      |
C|                           control by puff+recycling, the      |
C|                           number of puffeing surfaces and thus|
C|                           the number of additional surfaces   |
C|                           (NLIM) may vary. Before, the neutral|
C|                           flux file eirene.chemFluxDep from   |
C|                           previous run was checked to have the|
C|                           same number of add. surfaces as the |
C|                           actual run. This forced a stopping  |
C|                           of the code. Now the actual triangle|
C|                           number and its side is checked, to  |
C|                           asure that the right neutral flux   |
C|                           is used.                            |
C| 24/11/2010   D.Harting    Added for backward compatibility    |
C|                           a version number to the neutral flux|
C|                           file eirene.chemFluxDep. If the     |
C|                           Version number in the file and in   |
C|                           code are not matching, the file is  |
C|                           not read and zero neutral flux is   |
C|                           assumed in the actual run. At the   |
C|                           end of the eirene run, a new neutral|
C|                           flux file with the current version  |
C|                           number is generated.                |
C| 22/03/2011   D.Harting    Do not exit anymore if one of the   |
C|                           geometrical parameters in           |
C|                           eirene.chemFluxDep are not matching.|
C|                           Just ignore the neutral flux from   |
C|                           the previous run and continue.      |
C+---------------------------------------------------------------+
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
      IMPLICIT NONE
      REAL(DP) :: FACTOR, EIRENE_STEP
      INTEGER :: NLINES, ITRI, ISIDE, I, NBIN, ISTRA, ISRFS, ISOR, JJJ,
     .           ITEC1, ITEC2, ITEC3, ISTEP, INDSRF, IS1, IERROR, IPLS
      INTEGER :: EIRENE_IDEZ
      INTEGER, ALLOCATABLE :: KSTEP(:), INOSRC(:), IPLAN(:), IPLEN(:)
      REAL(DP) :: FLX, TE, TI, DE, MC, FE, FI, FSH, VP, FEL, DUM
      REAL(DP) :: DELR, FL, MCC, FFEL, CS, vx,vy,vz,di,usrval
      real(dp) :: xref, yref, bzref, facbz, x, y, rad, bx, by, bz, bf,
     .            errbx, errby, errbz, errbf, bxmax, bymax, bzmax, 
     .            bfmax, dfdx, dfdy, dfdz, xref2, yref2, bzref2, facbz2,
     .            xref3, yref3, bzref3, facbz3
      integer :: nref, icell, EIRENE_learc1, nplcll, ipolg, nref2, nref3
      CHARACTER(256) :: line,sstr,filename
      character(2) :: cstr2
      integer, allocatable :: indextmp(:),itritmp(:),isidetmp(:)
      integer, allocatable :: isegtmp(:)
      real(dp), allocatable :: tetmp(:),fetmp(:),fshtmp(:)
      REAL(DP), ALLOCATABLE :: PSI(:), PSI_CORNER(:), copy(:)

      integer, parameter :: fp=31
      integer :: ll,ier,j,iseg,ind, fp2

      integer, save :: eirene_nbirth,eirene_njetto
      character(len=256), save :: eirene_fbirth,eirene_ftransfer,
     &     eirene_fstoreneutflux
      real(dp) :: eirene_phi_offsets(9)

      integer :: ntr, NLIM_tmp, NSTS_tmp, NGITT_tmp, NGSTAL_tmp, 
     &     NATM_tmp, IPLS_tmp, NTR_tmp, NLMPGS_tmp
      REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: hydIonFLX, EIRENE_wall_area
      REAL(DP),ALLOCATABLE,DIMENSION(:)   :: hydNeutFLX
      INTEGER, ALLOCATABLE,DIMENSION(:,:) :: hydNeutFLX_info
      logical :: lex,dbg_out
      real(dp):: leng, twopi, tmp
      integer :: eirene_wallFluxModel ! calculation of wall fluxes for chemical sputtering
c                = 0: no wall fluxes are used (old edge2d model)
c                = 1: only ion fluxes are used
c                = 2: ion fluxes and neutral fluxes from last eirene iteration are used
c                = 3: ion and neutral fluxes are used, and EIRENE is iterated to give 
c                     converged neutral fluxes.
      logical :: eirene_use_elstepdat_bug
      real(dp) :: neutralFluxFileVersion

      namelist /eirene_user/eirene_nbirth,eirene_njetto,
     .                      eirene_fbirth,eirene_ftransfer,
     .                      eirene_phi_offsets,
     .                      eirene_fstoreneutflux,
     .                      eirene_wallFluxModel,
     .                      eirene_use_elstepdat_bug



      interface
        subroutine EIRENE_cell_to_corner (f, fcorner)
          use eirmod_precision
          implicit none
          real(dp), intent(in) :: f(:)
          real(dp), intent(out) :: fcorner(:)
        end subroutine EIRENE_cell_to_corner

        subroutine EIRENE_df_dxyz (fcorner, icell, x, y, z, 
     .                      dfdx, dfdy, dfdz) 
          use eirmod_precision
          implicit none
          real(dp), intent(in) :: fcorner(:), x, y, z
          real(dp), intent(out) :: dfdx, dfdy, dfdz
          integer, intent(in) :: icell
        end subroutine EIRENE_df_dxyz 
      end interface

      CALL EIRENE_ALLOC_CSTEP
      ALLOCATE (KSTEP(NSTEP))
      ALLOCATE (INOSRC(NSTEP))
      ALLOCATE (IPLAN(NSTEP))
      ALLOCATE (IPLEN(NSTEP))
      KSTEP = 0
      INOSRC = 0
      IPLAN = 0
      IPLEN = 0

      dbg_out=.true.

c     set default namelist values
      eirene_ftransfer = 'eirene.transfer'
      eirene_njetto=0
cdmh
      eirene_fstoreneutflux = 'eirene.chemFluxDep'
      eirene_wallFluxModel = 1
      NeutralFluxFileVersion = 1.0
cdmh      

c     get namelist config
      open(unit=9998,file='eirene_user.namelist')
      read(9998,eirene_user)
      close(9998)

c begin dmh added 21.06.2010
      ll=len_trim(casename)
      filename=casename(1:ll) // '.zplasma'
      open (unit=fp+ifoff,file=filename,access='sequential',
     .     form='formatted')

c first get number of triangles from misc plasma data:
        write(sstr,'(a20)') 
     .          '*** MISC PLASMA DATA'
        CALL EIRENE_locstr_usr(fp+ifoff,sstr,ier)
        if(ier /=0) then
           write(*,*) 'PLAUSR: ',sstr,' not found'
           close(fp+ifoff)
           call EIRENE_exit_own(1)
        endif
        do j=1,2
           read(fp+ifoff,'(a)') line
        enddo

        read(fp+ifoff,'(i7)') ntr
        if (ntr /= nr1st-1) then
           write (*,*) 'PLAUSR:', sstr
           write (*,*) ' wrong number of triangles in plasma file'
           write (*,*) ' check for correct number in file ',filename
           call EIRENE_exit_own(1)
        endif
c allocate temporary array to store plasma flux (hydrogen isotope)
        allocate(hydIonFLX(3,ntr))
        hydIonFLX(:,:) = 0.0
c allocate temporary array to store wall area of elements
        allocate(EIRENE_wall_area(3,ntr))
        EIRENE_wall_area(:,:) = 0.0
c end dmh added 21.06.2010

c firstly, read misc target data
      write(sstr,'(a20)') '*** MISC TARGET DATA'
      CALL EIRENE_locstr_usr(fp+ifoff,sstr,ier)
      if(ier /=0) then
         write(*,*) 'PLAUSR: ',sstr,' not found'
         close(fp+ifoff)
         CALL EIRENE_exit_own(1)
      endif
      do j=1,2
         read(fp+ifoff,'(a)') line
      enddo
      READ (fp+ifoff,*) NLINES
cswx 24sep07
      if(nlines == 0) then
        close(fp+ifoff)
        return
      endif
cswx
      if(nlines <= 0) then
         write(*,*) ' PLAUSR: error, nlines=',nlines
         close(fp+ifoff)
         call EIRENE_exit_own(1)
      endif

      allocate(indextmp(nlines))
      allocate(itritmp(nlines))
      allocate(isidetmp(nlines))
      allocate(tetmp(nlines))
      allocate(fetmp(nlines))
      allocate(fshtmp(nlines))
      allocate(isegtmp(nlines))
      do j=1,nlines
         read(fp+ifoff,'(3i7,3(1x,e14.7))') 
     .        indextmp(j),
     .        itritmp(j),
     .        isidetmp(j),
     .        tetmp(j),
     .        fetmp(j),
     .        fshtmp(j)
         isegtmp(indextmp(j)) = j
      enddo


c now read species dependent target data
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
               WRITE (*,*) 'ERROR IN PRIMARY SOURCE DATA '
               WRITE (*,*) 'STEPFUNCTION REQUESTED FOR SOURCE SURFACE '
               WRITE (*,*) 'NO. ',INSOR(ISRFS,ISTRA),' BUT SORIND.EQ.0.'
               CALL EIRENE_EXIT_own(1)
            ELSEIF (ISTEP.GT.NSTEP) THEN
               CALL EIRENE_MASPRM('NSTEP',5,NSTEP,'ISTEP',5,ISTEP,IERROR)
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
c     fudge species index for atomic impurity flux (get it from NEMODS index K)
            IPLS_tmp = EIRENE_IDEZ(NEMODS(ISTRA),4,4)
            IF (IPLS_tmp.gt.1) then 
               ipls=IPLS_tmp
            ENDIF
c     search target tag in .zplasma file
           IPLS_tmp = IPLS 
           write(cstr2,'(i2.2)') ipls
           write(sstr,'(a23)') 
     .          '*** ION #'//cstr2//' TARGET DATA'
           CALL EIRENE_locstr_usr(fp+ifoff,sstr,ier)
           if(ier /=0) then
              write(*,*) 'PROUSR: tag ',sstr,'not found'
              close(fp+ifoff)
              call EIRENE_exit_own(1)
           endif
           do j=1,2
              read(fp+ifoff,'(a)') line
           enddo

c     read step functions
            READ (fp+ifoff,*) NLINES
            DO I=1, NLINES
               read(fp+ifoff,'(i7,11(1x,e14.7))')
     .              ind,
     .              flx,ti,di,
     .              vx,vy,vz,
     .              fi,fel,
     .              vp,mc,
     .              usrval

               iseg = isegtmp(ind)
               itri = itritmp(iseg)
               iside= isidetmp(iseg)
               te   = tetmp(iseg)
               fe   = fetmp(iseg)
               fsh  = fshtmp(iseg)

c begin added dmh 21.06.2010
c     store ion flux in [A] of main plasma
               if(ipls_tmp.eq.1) then
                  hydIonFLX(iside,itri) = flx
               endif
c end added dmh 21.06.2010

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
                     VXSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = VXIN(IPLS,ITRI)
                     VYSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = VYIN(IPLS,ITRI)
                     VZSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = VZIN(IPLS,ITRI)
c                     VXSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = VX
c                     VYSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = VY
c                     VZSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = VZ
                     FLSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = ABS(FLX)/DELR
C     IF NO ION KINETIC ENERGY FLUX IS SPECIFIED, DERIVE IT FROM  FI,MC,VP
                     ffel=0.
                     IF (FI.gt.0..or.mc.gt.0.) then
                        ffel=(FI*TI+0.5*MC*MC*(TE+TI))*abs(FLX)
                     endif
                     IF (FEL.EQ.0.) FEL=FFEL
                     IF (eirene_use_elstepdat_bug) THEN
                        ELSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = FEL
                     ELSE
                        ELSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = FEL/DELR
                     ENDIF
csw
                     usrstep(ipls,istep,kstep(istep)) = usrval
csw
c                 enddo ipls,iplan
                  ENDDO

c              endif inmti
               ENDIF

c           enddo nlines
            ENDDO 

c        enddo isrfs
         ENDDO 

c     enddo istra
      ENDDO



      DO ISTEP = 1, NSTEP
         IF (KSTEP(ISTEP) > 0) THEN
            NBIN=KSTEP(ISTEP)+1
            FL=EIRENE_STEP(IPLAN(ISTEP),IPLEN(ISTEP),NBIN,ISTEP)
            FLUX(INOSRC(ISTEP))=FL
         END IF
      END DO
 
c begin added dmh 21.06.2010 for flux dependency of chemical sputtering
c     read neutral flux [A] to target and walls from last EIRENE run

      allocate(hydNeutFLX(NLMPGS))
      hydNeutFLX(:) = 0.0         
      allocate(hydNeutFLX_info(NLMPGS,3))
      hydNeutFLX_info = 0
      NLIM_tmp   = NLIM
      NSTS_tmp   = NSTS
      NGITT_tmp  = NGITT
      NGSTAL_tmp = NGSTAL
      NATM_tmp   = NATM
      NLMPGS_tmp = NLMPGS
      NTR_tmp    = NTR
      fp2 = 4999
      inquire(file=trim(eirene_fstoreneutflux),exist=lex)
      if (lex) then
          open(unit=fp2,file=trim(eirene_fstoreneutflux),
     &        access='sequential')
         read(fp2,'(a)') line
         read(line,'(a28)') sstr

c        check if netral flux file is compatible with actual code version
         if (index(sstr,"* Neutral flux file version:").ne.1) then
            WRITE(IUNOUT,*) "* EIRENE_PLAUSR: Warning"
            write(IUNOUT,*) "Found obsolete neutral flux file: ",
     &           trim(eirene_fstoreneutflux)
            write(IUNOUT,*) "Ignoring neutral flux from previous run"
            WRITE(IUNOUT,*) "* EIRENE_PLAUSR: Warning end"
            lex =.false.
         else
            read(line,'(a28,f14.6)') sstr,tmp
            if (tmp.ne.NeutralFluxFileVersion) then
               WRITE(IUNOUT,*) "* EIRENE_PLAUSR: Warning"
               write(IUNOUT,'(a,a,a,f14.6)') 
     &              "Found obsolete neutral flux file: ",
     &              trim(eirene_fstoreneutflux),"; version:",tmp
               write(IUNOUT,'(a,f14.6)') 
     &              "but actual neutral flux file should be version:",
     &              NeutralFluxFileVersion
               write(IUNOUT,*) "Ignoring neutral flux from previous run"
               WRITE(IUNOUT,*) "* EIRENE_PLAUSR: Warning end"
               lex =.false.
            else
c              neutral flux file is compatible with actual code version, so read it
               read(fp2,'(a)') line
               read(fp2,'(7i8)') NLIM_tmp, NSTS_tmp, NGITT_tmp,
     &              NGSTAL_tmp, NATM_tmp, NLMPGS_tmp, NTR_tmp
               if ((NLMPGS_tmp-NLIM_tmp-NSTS_tmp.ne.NLMPGS-NLIM-NSTS)
     &              .or.(NSTS_tmp.ne.NSTS).or.(NTR_tmp.ne.NTR).or.
     &              (NGITT_tmp.ne.NGITT).or.(NGSTAL_tmp.ne.NGSTAL)) then
                  WRITE(IUNOUT,*) "* EIRENE_PLAUSR: Warning"
                  WRITE(IUNOUT,*) "Neutral flux file ",
     &                 trim(eirene_fstoreneutflux),
     &                 " does not fit to simulation!"
                  WRITE(IUNOUT,*)"NLIM = ",NLIM,"; NLIM_tmp = ",NLIM_tmp
                  WRITE(IUNOUT,*)"NSTS = ",NSTS,"; NSTS_tmp = ",NSTS_tmp
                  WRITE(IUNOUT,*)"NGITT = ",NGITT,
     &                 "; NGITT_tmp = ",NGITT_tmp
                  WRITE(IUNOUT,*)"NGSTAL = ",NGSTAL,
     &                 "; NGSTAL_tmp = ",NGSTAL_tmp
                  WRITE(IUNOUT,*)"NATM = ",NATM,"; NATM_tmp = ",NATM_tmp
                  WRITE(IUNOUT,*)"NTRII= ",NTR,"; NTRII_tmp = ",NTR_tmp
                  WRITE(IUNOUT,*)"NLMPGS_tmp-NLIM_tmp-NSTS_tmp = ",
     &                 NLMPGS_tmp-NLIM_tmp-NSTS_tmp,
     &                 "NLMPGS-NLIM-NSTS = ",NLMPGS-NLIM-NSTS
                  write(IUNOUT,*) 
     &                 "Ignoring neutral flux from previous run"
                  WRITE(IUNOUT,*) "* EIRENE_PLAUSR: Warning end"
                  lex =.false.
               endif

               if (lex) then
                  deallocate(hydNeutFLX, hydNeutFLX_info)
                  allocate(hydNeutFLX(NLMPGS_tmp))
                  hydNeutFLX(:) = 0.0            
                  allocate(hydNeutFLX_info(NLMPGS_tmp,3))
                  hydNeutFLX_info = 0
                  
                  read(fp2,'(a)') line
                  read(fp2,'(a)') line
                  do i=1,NLMPGS_tmp
                     read(fp2,'(i6,1x,e14.6,1x,i6,1x,i6,1x,i6)') 
     &                    j, hydNeutFLX(i),
     &                    hydNeutFLX_info(i,1), hydNeutFLX_info(i,2),
     &                    hydNeutFLX_info(i,3)
                     if (i.ne.j) then
                        WRITE(IUNOUT,*) "* EIRENE_PLAUSR:"
                        write(IUNOUT,*) "Error reading from file: ",
     &                       trim(eirene_fstoreneutflux)
                        write(IUNOUT,*) "i = ",i,";  in file j = ",j
                        CALL EIRENE_EXIT_own(1)
                     endif
                  enddo         !i
               endif            !(lex)
            endif               !(tmp.ne.NeutralFluxFileVersion)
         endif                  !(index(sstr,"* Neutral flux file version:").ne.1)
         close(fp2)
      endif 
      
c     calculate area of wall surface elements
      twopi = 2.d0*dabs(dacos(-1.0))
      do itri=1,ntr
         do iside=1,3
            IS1 = ISIDE + 1
            IF (IS1.GT.3) IS1=1
            leng  = SQRT(
     .           (XTRIAN(NECKE(ISIDE,ITRI))
     .           -XTRIAN(NECKE(IS1,ITRI)))**2+
     .           (YTRIAN(NECKE(ISIDE,ITRI))
     .           -YTRIAN(NECKE(IS1,ITRI)))**2)
            EIRENE_wall_area(iside,itri) = twopi*leng
     &           *( XTRIAN(NECKE(ISIDE,ITRI)) 
     &           +  XTRIAN(NECKE(IS1,ITRI)) )/2.D0
         enddo                  ! iside
      enddo                     ! itri              

c     calculate target flux for chemical sputtering of Roth-formula
      FLXOUT(:) = 0.0

      if (eirene_wallFluxModel.ge.1) then
c     use hydrogen ion flux [A] to wall
         do itri=1,ntr
            do iside=1,3
               if ( (hydIonFLX(iside,itri).ne.0)
     &              .and.(INMTI(iside,itri).ne.0) )then
                  FLXOUT(NLIM+NSTS + INSPAT(iside,itri)) = 
     &                 FLXOUT(NLIM+NSTS + INSPAT(iside,itri))
     &                 + DABS(hydIonFLX(iside,itri))
               endif
            enddo               ! iside
         enddo                  ! itri
      endif

      if ((eirene_wallFluxModel.ge.2).and.lex) then
c     use hydrogen neutral flux [A] to wall
         do i=NLIM_tmp+NSTS_tmp+1,NLMPGS_tmp
            j = i-NLIM_tmp-NSTS_tmp+NLIM+NSTS            
            if (hydNeutFLX(i).ne.0) then
               if (hydNeutFLX_info(i,1).ne.0) then
                  lex=.false.
                  do itri=1,ntr
                     do iside=1,3
                        if ((INSPAT(iside,itri).eq. j -(NLIM+NSTS))
     &                       .and.(INSPAT(iside,itri).ne.0) )then
                           if (lex) then
                              WRITE(IUNOUT,*) "* EIRENE_PLAUSR:"
                              write(IUNOUT,*) "* Edge twice found"
                              CALL EIRENE_EXIT_own(1)
                           endif
                           lex=.true.
                           
                           if  ((itri .eq. hydNeutFLX_info(i,1)).or.
     &                          (iside .eq. hydNeutFLX_info(i,2)).or.
     &                          (INMTI(iside,itri)-NLIM-NSTS.eq.
     &                          hydNeutFLX_info(i,3)-NLIM_tmp-NSTS_tmp))
     &                          then
                              FLXOUT(j) = FLXOUT(j)
     &                             + DABS(hydNeutFLX(i))
                           else
                              WRITE(IUNOUT,*) "* EIRENE_PLAUSR:"
                              write(IUNOUT,*) 
     &                             "Neutral flux from previous run ",
     &                             "is not associated with the same ",
     &                             "triangle in this run"
                              write(IUNOUT,*) "itri_old = ",
     &                             hydNeutFLX_info(i,1),
     &                             "; itri_new = ",itri
                              write(IUNOUT,*) "iside_old = ",
     &                             hydNeutFLX_info(i,2),
     &                             "; iside_new = ",iside
                              write(IUNOUT,*) "isurf_old = ",
     &                             hydNeutFLX_info(i,3),
     &                             "; isurf_new = ",itri
                              CALL EIRENE_EXIT_own(1)
                           endif ! same triangle associated
                        endif   ! triangle found
                     enddo      ! iside
                  enddo         ! itri
               else
                  WRITE(IUNOUT,*) "* EIRENE_PLAUSR:"
                  write(IUNOUT,*) 
     &                 "* No triangle associated with neutral flux"
                  CALL EIRENE_EXIT_own(1)
               endif            ! hydNeutFLX_info(i).ne.0
            endif               ! hydNeutFLX(i).ne.0
         enddo                  ! i
      endif

c     Debug output of boundary
      if (dbg_out) then
      open(unit=fp2, file="eirene.chemSput_dbgOut",access='sequential')
      write(fp2,*) NLMPGS
      write(fp2,'(a,a,a)') 
     &     "* index, nsurf,           R1,           R2,",
     &     "           Z1,           Z2,       FLXOUT,",
     &     "      NeutFLX,       IonFLX,         area"
      do i=1,NLMPGS
         lex=.false.
         do itri=1,ntr
            do iside=1,3
               if ((INSPAT(iside,itri).eq. i -(NLIM+NSTS))
     &              .and.(INSPAT(iside,itri).ne.0) )then
                  if (lex) write(fp2,*)"*Edge twice found"
                  lex=.true.
                  IS1 = ISIDE + 1
                  IF (IS1.GT.3) IS1=1                 
                  write(fp2,'(I8,I7,8(x,e13.6))')
     &                 i,INMTI(iside,itri),
     &                 XTRIAN(NECKE(ISIDE,ITRI)),
     &                 XTRIAN(NECKE(IS1,ITRI)),
     &                 YTRIAN(NECKE(ISIDE,ITRI)),
     &                 YTRIAN(NECKE(IS1,ITRI)), 
     &                 FLXOUT(i)/
     &                 (1.6022D-19*EIRENE_wall_area(iside,itri)),
     &                 hydNeutFLX(i)
     &                 /(1.6022D-19*EIRENE_wall_area(iside,itri)),
     &                 hydIonFLX(iside,itri)
     &                 /(1.6022D-19*EIRENE_wall_area(iside,itri)),
     &                 EIRENE_wall_area(iside,itri)
               endif
            enddo
         enddo
         if (.not.lex) then
            write(fp2,'(I8,I7,8(x,e13.6))')
     &           i,0,0.D0,0.D0,0.D0,0.D0,FLXOUT(i),hydNeutFLX(i),0.D0
         endif
      enddo
      close(fp2)

      
      open(unit=fp2, file="eirene.chemSput_dbgOut2",access='sequential')
      write(fp2,*) ntr
      write(fp2,'(a,a)') 
     &     "*   itri, iside, nsurf,           R1,           R2,",
     &     "           Z1,           Z2"
      do itri=1,ntr
         do iside=1,3
            IS1 = ISIDE + 1
            IF (IS1.GT.3) IS1=1                 
            write(fp2,'(I8,x,I6,x,I6,4(x,e13.6))')
     &                 itri,iside,INMTI(iside,itri),
     &                 XTRIAN(NECKE(ISIDE,ITRI)),
     &                 XTRIAN(NECKE(IS1,ITRI)),
     &                 YTRIAN(NECKE(ISIDE,ITRI)),
     &                 YTRIAN(NECKE(IS1,ITRI))
         enddo
      enddo
      close(fp2)
      endif

c     FLXOUT is needed in #/(cm^2 s) (covert from A to #/(cm^2 s))
      do itri=1,ntr
         do iside=1,3
            if (INMTI(iside,itri).ne.0) then
               FLXOUT(NLIM+NSTS + INSPAT(iside,itri)) = 
     &              FLXOUT(NLIM+NSTS + INSPAT(iside,itri))
     &              /(1.6022D-19*EIRENE_wall_area(iside,itri))
            endif
         enddo                  ! iside
      enddo                     ! itri


c end added dmh 21.06.2010 for flux dependency of chemical sputtering

c dmh begin added output of step function data
      if (dbg_out) then
      open(unit=fp2, file="eirene.stepdat_dbg",access='sequential')
      write(fp2,"(A,A,A)") "ipls istep i  IRSTEP(ITRI)  IPSTEP(ISIDE)",
     &     "  X1   Y1  DR  FLSTEP  ELSTEP  SHSTEP  TISTEP",
     &     "  TESTEP"
      do ipls=1,npls
         do istep=1,nstep
            do i=1,KSTEP(ISTEP)
               ITRI = IRSTEP(ISTEP,I)
               ISIDE= IPSTEP(ISTEP,I)
               IS1 = ISIDE+1
               if(IS1.GT.3) IS1=1
               write(fp2,"(2(I3,x),I4,x,I6,x,I1,8(x,e13.6))")
     &              ipls, istep, i, ITRI, ISIDE, 
     &              XTRIAN(NECKE(ISIDE,ITRI)),YTRIAN(NECKE(ISIDE,ITRI)),
     &              (RRSTEP(istep,i+1)-RRSTEP(istep,i)),
     &              FLSTEP(IPLS,ISTEP,i), ELSTEP(IPLS,ISTEP,i), 
     &              SHSTEP(ISTEP,i), TISTEP(IPLS,ISTEP,i), 
     &              TESTEP(ISTEP,i)
               write(fp2,"(2(I3,x),I4,x,I6,x,I1,8(x,e13.6))")
     &              ipls, istep, i, ITRI, ISIDE, 
     &              XTRIAN(NECKE(IS1,ITRI)), YTRIAN(NECKE(IS1,ITRI)),
     &              (RRSTEP(istep,i+1)-RRSTEP(istep,i)),
     &              FLSTEP(IPLS,ISTEP,i), ELSTEP(IPLS,ISTEP,i), 
     &              SHSTEP(ISTEP,i), TISTEP(IPLS,ISTEP,i), 
     &              TESTEP(ISTEP,i)
               write(fp2,*)
            enddo
         enddo
      enddo
      close(fp2)
      endif
c dmh begin added output of step function data


c     cleanup
      DEALLOCATE (KSTEP)
      DEALLOCATE (INOSRC)
      DEALLOCATE (IPLAN)
      DEALLOCATE (IPLEN)
      deallocate(indextmp)
      deallocate(itritmp)
      deallocate(isidetmp)
      deallocate(tetmp)
      deallocate(fetmp)
      deallocate(fshtmp)
      deallocate(isegtmp)
      deallocate(hydIonFLX)
      deallocate(hydNeutFLX)
      deallocate(hydNeutFLX_info)
      deallocate(EIRENE_wall_area)

      close(fp+ifoff)


      return
!pb      return
 4711 continue

      if (nain < 1) then
         write (iunout,*) ' no additional input tally available '
         write (iunout,*) ' for PSI-function'
         write (iunout,*) ' calculation of B field from PSI function'
         write (iunout,*) ' abandonned '
         return
      end if

      ALLOCATE (PSI(NRAD))
      ALLOCATE (PSI_CORNER(NRAD))
      ALLOCATE (COPY(NRAD))
      PSI=0.
      PSI_CORNER=0.
      COPY=0._dp

      if (naini >= 8) then
        call EIRENE_prousr(copy,1+5*npls,0._dp,0._dp,
     .        0._dp,0._dp,0._dp,0._dp,0._dp,nsbox)
        adin(7,1:nsbox) = copy(1:nsbox)

        call EIRENE_prousr(copy,2+5*npls,0._dp,0._dp,
     .       0._dp,0._dp,0._dp,0._dp,0._dp,nsbox)
        adin(8,1:nsbox) = copy(1:nsbox)
      end if 

      call EIRENE_prousr(psi,5+5*npls,0._dp,0._dp,0._dp,
     .     0._dp,0._dp,0._dp,0._dp,nsbox)

      psi(1:nsbox) = psi(1:nsbox) * 1.e4_dp

      adin(1,1:nsbox) = psi(1:nsbox)

      call EIRENE_cell_to_corner (psi, psi_corner)

      xref = 200._dp
      yref = 0._dp
      NREF=EIRENE_LEARC1(XREF,YREF,0._DP,IPOLG,1,NR1STM,
     .     .FALSE.,.FALSE.,1,'PLAUSR      ')
      BZREF=BZIN(NREF)*BFIN(NREF)
      FACBZ=xcom(nref)*BZREF

      write (iunout,*) ' xref,  yref,  bzref,  facbz ',
     .                   xref,  yref,  bzref,  facbz

      xref2= 380._dp
      yref2 = 0._dp
      NREF2=EIRENE_LEARC1(XREF2,YREF2,0._DP,IPOLG,1,NR1STM,
     .     .FALSE.,.FALSE.,1,'PLAUSR      ')
      BZREF2=BZIN(NREF2)*BFIN(NREF2)
      FACBZ2=xcom(nref2)*BZREF2

      write (iunout,*) ' xref2, yref2, bzref2, facbz2 ',
     .                   xref2, yref2, bzref2, facbz2

      xref3= 250._dp
      yref3 = 150._dp
      NREF3=EIRENE_LEARC1(XREF3,YREF3,0._DP,IPOLG,1,NR1STM,
     .     .FALSE.,.FALSE.,1,'PLAUSR      ')
      BZREF3=BZIN(NREF3)*BFIN(NREF3)
      FACBZ3=xcom(nref3)*BZREF3

      write (iunout,*) ' xref3, yref3, bzref3, facbz3 ',
     .                   xref3, yref3, bzref3, facbz3

      errbx = 0._dp
      errby = 0._dp
      errbz = 0._dp
      errbf = 0._dp
      bxmax = 0._dp
      bymax = 0._dp
      bzmax = 0._dp
      bfmax = 0._dp
      nplcll = 0

      do icell = 1, ntrii
        x = xcom(icell)
        y = ycom(icell)
        rad = x
        call EIRENE_df_dxyz (psi_corner, icell, x, y, 
     &       0._dp, dfdx, dfdy, dfdz) 

        bx = -dfdy / rad
        by =  dfdx / rad
        bz = facbz / rad
        bf = sqrt(bx*bx + by*by + bz*bz)

        if (naini >= 6) then
        adin(3,icell) = bx
        adin(4,icell) = by
        adin(5,icell) = bz
        adin(6,icell) = bf

        else

        bx = bx / bf
        by = by / bf
        bz = bz / bf
        if (ixtri(icell)+iytri(icell) == 0) then
! outside plasma region: set b-field
          bxin(icell) = bx
          byin(icell) = by
          bzin(icell) = bz
          bfin(icell) = bf
        else
          nplcll = nplcll + 1
          bxmax = max(bxmax, abs(bxin(icell)))
          bymax = max(bymax, abs(byin(icell)))
          bzmax = max(bzmax, abs(bzin(icell)))
          bfmax = max(bfmax, abs(bfin(icell)))

          errbx = errbx + abs(bx-bxin(icell))
          errby = errby + abs(by-byin(icell))
          errbz = errbz + abs(bz-bzin(icell))
          errbf = errbf + abs(bf-bfin(icell))
          if (icell == nref) then
            write (iunout,*) ' reference cell ',nref
            write (iunout,*) ' bxin, byin, bzin, bfin ',
     .             bxin(icell), byin(icell), bzin(icell), bfin(icell)
            write (iunout,*) ' bx,   by,   bz,   bf   ',
     .             bx, by, bz, bf
          elseif (icell == nref2) then
            write (iunout,*) ' reference cell ',nref2
            write (iunout,*) ' bxin, byin, bzin, bfin ',
     .             bxin(icell), byin(icell), bzin(icell), bfin(icell)
            write (iunout,*) ' bx,   by,   bz,   bf   ',
     .             bx, by, bz, bf
          elseif (icell == nref3) then
            write (iunout,*) ' reference cell ',nref3
            write (iunout,*) ' bxin, byin, bzin, bfin ',
     .             bxin(icell), byin(icell), bzin(icell), bfin(icell)
            write (iunout,*) ' bx,   by,   bz,   bf   ',
     .             bx, by, bz, bf
          end if
        end if
        end if

      end do

      adin(2,1:nsbox) = bzin(1:nsbox) * bfin(1:nsbox)

      if (.false.) then
      errbx = errbx / nplcll
      errby = errby / nplcll
      errbz = errbz / nplcll
      errbf = errbf / nplcll

      write (iunout,*) ' no of cell in plasma region ', nplcll
      write (iunout,*) ' errbx, errby, errbz, errbf ',
     .                   errbx, errby, errbz, errbf

      write (iunout,*) ' bxmax, bymax, bzmax, bfmax ',
     .                   bxmax, bymax, bzmax, bfmax
      end if

 99   RETURN
      END
