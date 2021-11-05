Cdr  Purpose: find storage parameters NPARM for a number of allocatable
c             arrays.
c             Set default for NPARM             (for example NATM=0, no atomic species in this run)
c             Read input file, to find NPARMI,  (for example NATMI)
c             Then set NPARM=MAX(NPARM,NPARMI)  (for example NATM=MAX(NATM,NATMI))
c
c    Later these parameters may be modified, but NPARM >= NPARMI always.
C
!pb  11.12.06:  allow letters 'f' or 't' in case name of fem or tetrahedron
!pb             calculation
!pb  27.12.06:  bug fix: increase NSTS in case of time dependent mode
!pb  15.01.07:  additional line in input block 4 defining HYDKIN  model
!pb  02.03.07:  NUMSEC=4 introduced
!pb  20.03.07:  include input block written by HYDKIN model
!pb  22.03.07:  input for NLFEM and NLTET corrected.
!dr  16.01.14:  default NOPTIM changed from 1 to NRAD (automatically), some printout rearranged
!cd  29.10.14:  reading external file for block 4&5: allow comment lines at the beginning of file
!               (same in find-param)
!cd  2.2.15:    nfr (number of TRIM A_on_B files), now same name as in input.f
cdr  Jan 2016:  storage for second dimension only if nlpol=true
cdr             to be tested: storage for nplg, if nlpol=false?
cdr             storage for third dimension only if nltor=true
cdr             to be tested:  storage for nltra, if nltor=false?
cdr             to be done: check for comment lines *... syncronized with input.f?
!pb  June 16:   default for NPLSTI changed from 1 to NPLS
!pb  MAY  16:   nrds -> nrei
cdr  March 17:  NPTRGT printed. May have been changed in call to if0parm, block 14.
CDR  May 2017:  try to fix NSTRAI, NSRFSI, consistent with input.f
cdr             same thing: NCPVI, NCPV  (and eliminate old parameters NCOP, NCOPI)
cdr  July 17 :  lmulti, lmulvi:  automatic options for multiple ion temperatures,
cdr                              multiple ion velocities in case of BGK non-lin. colisions
cdr  July 17 :  initialize 2D CFD code coupling parameters NDX,....
c               move nrad=... after call to if0prm, because of emc3 coupling
C
      SUBROUTINE EIRENE_FIND_PARAM_JSON
C
C   SET DIMENSIONS (PARAMETERS) FOR ALLOCATABLE ARRAYS (TALLIES, GRIDS, ETC..) AND SET SOME FURTHER DEFAULT VALUES
C
      USE EIRMOD_PRECISION, ONLY: DP
      USE EIRMOD_PARMMOD, ONLY: N1ST, N2ND, N3RD, NADD, NADS, NADV,
     .                          NAIN, NALS, NALV, NATM, NBGK, NCLV,
     .                          NCOORD, NCOUP_INPUT, NCPV, NCV,
     .                          NGEOM_USR, NGSTAL, NION, NKNOT, NLIM,
     .                          NLIMPS, NMOL, NOPTIM, NOPTM1, NPHOT,
     .                          NPLG, NPLS, NPPART, NPRNL, NPTRGT, NRAD,
     .                          NREAC, NRTAL, NSD, NSDW, NSMSTRA, NSNV,
     .                          NSRFS, NSTEP, NSTORAM, NSTRA, NSTS,
     .                          NTETRA, NTOR, NTRI, NCHEN, NCHOR, NPLT,
     .                          NSRPR, NVLPR, NHD6, NPLSTI, NPLSV,
     .                          NREAC_LINES, IFOFF, NBMAX, NGITT,
     .                          NKNOTS, NTRIS
      USE EIRMOD_COMUSR, ONLY: NMODE, NTIME,
     .                         NATMI, NMOLI, NPLSI
      USE EIRMOD_COMSOU, ONLY: NSTRAI
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_JSON, ONLY: jtrees, blks, itree_num,
     .                       ldef_time_horizon,nlfem_in, nlplg_in,
     .                       nlpol_in, np2nd_in, nr1st_in, nt3rd_in,
     .                       noptim_in, nrtal_in, nsmstra_in,
     .                       eirene_init_input_blocks 
      
      use json_module, ck => json_ck

      IMPLICIT NONE

      INTERFACE
        SUBROUTINE EIRENE_IF0PRM_JSON(json,p)
        use json_module
        class(json_core),intent(inout) :: json
        type(json_value), pointer, intent(in) :: p
        END SUBROUTINE EIRENE_IF0PRM_JSON
      END INTERFACE
      
      type(json_value),pointer :: p

!pb   integer :: nreac_add
!      integer :: nstsi, ll, ic, i, ipls, ispz, j
      integer :: nstsi, j
      logical :: found0, ldefstor, nlerg
!      logical :: found0, ldefstor, lext, nlerg, lhyddef, ladapt
      logical :: lmulti, lmulvi ! multiple ion temperatures (per species) multiple ion velocities (per species)
!      CHARACTER(420) :: FILENAME, ULINE, ZEILE
!      CHARACTER(12) :: CHR
!      CHARACTER(1000) :: HLINE
!      character(len=:), allocatable :: name
      character(kind=CK,len=:), allocatable :: header

C
C  SET DEFAULT VALUES FOR STORAGE PARAMETERS
C
      CALL EIRENE_INIT_PARAMS
!pb   NREAC_ADD=0
      LDEFSTOR=.false.
      
      lmulti = .false.
      lmulvi = .false.
  
      CALL EIRENE_LEER(3)

      WRITE (IUNOUT,*) 'PRINTOUT FROM EIRENE PRE-PROCESSING:'
      WRITE (IUNOUT,*) 'BROWSE INPUT FOR STORAGE NEEDS (FIND_PARAM.F)'
      WRITE (IUNOUT,*) 'FILE = ', 'eirene.input.json'
      CALL EIRENE_LEER(1)
C
c  start reading json file

      call eirene_init_input_blocks ('eirene.input.json', iunout)

C  read and write header

      j = itree_num(0)        ! index of the jtree that holds the data of the input block   
      call jtrees(j)%get(blks(0)%p,'TXTRUN',header,found0)
      write (iunout,'(A)') header
      CALL EIRENE_LEER(1)

!  browse block 1 "GENERAL_DATA" for parameters
      j = itree_num(1)
      call eirene_browse_block_1(jtrees(j),blks(1)%p)

!  browse block 2 "STANDARD_MESH" for parameters
      j = itree_num(2)
      call eirene_browse_block_2(jtrees(j),blks(2)%p)

      WRITE (iunout,*) '*** 3. DATA FOR BREP SURFACES'

!  browse block 3a "NONDEF_STD_SURFACES" for parameters
      j = itree_num(3)
      call eirene_browse_block_3a(jtrees(j),blks(3)%p)

!  browse block 3b "ADDITIONAL_SURFACES" for parameters
      j = itree_num(4)
      call eirene_browse_block_3b(jtrees(j),blks(4)%p)

!  browse block 4 and 5 "PHYSICS_MODEL" for parameters
      j = itree_num(5)
      call eirene_browse_block_4(jtrees(j),blks(5)%p)
      call eirene_browse_block_5(jtrees(j),blks(5)%p)

!  browse block 6 "REFLECTION_MODELS" for parameters
      j = itree_num(6)
      call eirene_browse_block_6(jtrees(j),blks(6)%p)

!  browse block 7 "SOURCES" for parameters
      j = itree_num(7)
      call eirene_browse_block_7(jtrees(j),blks(7)%p)

      WRITE (iunout,*) '*** 8. ADDITIONAL DATA FOR SPECIFIC ZONES'

!  browse block 8 "SPECIFIC_ZONES" for parameters
!  nothing to be done
      j = itree_num(8)

!  browse block 9 "STATISTICS" for parameters
      j = itree_num(9)
      call eirene_browse_block_9(jtrees(j),blks(9)%p)

!  browse block 10 "ADDITIONAL_TALLIES" for parameters
      j = itree_num(10)
      call eirene_browse_block_10(jtrees(j),blks(10)%p)

!  browse block 11 "OUTPUT" for parameters
      j = itree_num(11)
      call eirene_browse_block_11(jtrees(j),blks(11)%p)

!  browse block 12 "DIAGNOSTICS" for parameters
      j = itree_num(12)
      call eirene_browse_block_12(jtrees(j),blks(12)%p)

!  browse block 13 "TIMEDEPENDENT_MODE" for parameters
      j = itree_num(13)
      call eirene_browse_block_13(jtrees(j),blks(13)%p)

!  browse block 14 "INTERFACING" for parameters
      if (nmode == 0) then
        j = itree_num(14)
        call eirene_browse_block_14(jtrees(j),blks(14)%p)
      else 
! to be prepared
!       CALL EIRENE_IF0PRM(IUNIN)
        j = itree_num(14)
        call eirene_if0prm_json(jtrees(j),blks(14)%p)
      end if
     
      if ((NTIME.GE.1.AND.NPRNL > 0).OR.NLERG) THEN
        NSTSI=NSTSI+1
        NSTRAI=NSTRAI+1
      ENDIF
      NSTS = MAX(NSTS,NSTSI)
      NSTRA = MAX(NSTRA,NSTRAI)
      NLIMPS = NLIM + NSTS

C  SWITCH OFF SUM OVER STRATA IF THERE IS ONLY ONE STRATUM TO BE CALCULATED
!pb      IF (NSTRAI == 1) NSMSTRA = 0
 
cdr  due to these changes there, also some derived storage parmeters may have changed....
      NRAD=MAX(N1ST*N2ND*N3RD,NTRI*N3RD,NTETRA)+NADD+1 ! as in parmmod

      
!  output paramters for check 

      CALL EIRENE_LEER(1)
      WRITE (IUNOUT,'(A)') ' AUTOMATED STORAGE SETTING (FIND_PARAM.F)'
      CALL EIRENE_LEER(1)
C
cdr  grid size
      WRITE (iunout,'(a14,i8)') 'N1ST        = ',N1ST
      WRITE (iunout,'(a14,i8)') 'N2ND        = ',N2ND
      WRITE (iunout,'(a14,i8)') 'N3RD        = ',N3RD
      WRITE (iunout,'(a14,i8)') 'NADD        = ',NADD
      WRITE (iunout,'(a14,i8)') 'NTOR        = ',NTOR
      WRITE (iunout,'(a14,i8)') 'NRTAL       = ',NRTAL
      WRITE (iunout,'(a14,i8)') 'NLIM        = ',NLIM
      WRITE (iunout,'(a14,i8)') 'NSTS        = ',NSTS
      WRITE (iunout,'(a14,i8)') 'NPLG        = ',NPLG
      WRITE (iunout,'(a14,i8)') 'NPPART      = ',NPPART
      WRITE (iunout,'(a14,i8)') 'NKNOT       = ',NKNOT
      WRITE (iunout,'(a14,i8)') 'NTRI        = ',NTRI
      WRITE (iunout,'(a14,i8)') 'NTETRA      = ',NTETRA
      WRITE (iunout,'(a14,i8)') 'NCOORD      = ',NCOORD

      WRITE (iunout,*) ' '
      WRITE (iunout,'(a14,i8)') 'NRAD        = ',NRAD
cdr  primary source
      CALL EIRENE_LEER(1)
      WRITE (iunout,'(a14,i8)') 'NSTRA       = ',NSTRA
      WRITE (iunout,'(a14,i8)') 'NSRFS       = ',NSRFS
      WRITE (iunout,'(a14,i8)') 'NSTEP       = ',NSTEP
      WRITE (iunout,'(a14,i8)') 'NPTRGT      = ',NPTRGT
cdr species
      CALL EIRENE_LEER(1)
      WRITE (iunout,'(a14,i8)') 'NATM        = ',NATM
      WRITE (iunout,'(a14,i8)') 'NMOL        = ',NMOL
      WRITE (iunout,'(a14,i8)') 'NION        = ',NION
      WRITE (iunout,'(a14,i8)') 'NPHOT       = ',NPHOT
      WRITE (iunout,'(a14,i8)') 'NPLS        = ',NPLS

      CALL EIRENE_LEER(1)
      WRITE (iunout,'(a14,i8)') 'NADV        = ',NADV    
      WRITE (iunout,'(a14,i8)') 'NCLV        = ',NCLV
      WRITE (iunout,'(a14,i8)') 'NSNV        = ',NSNV
      WRITE (iunout,'(a14,i8)') 'NALV        = ',NALV

      WRITE (iunout,'(a14,i8)') 'NADS        = ',NADS
      WRITE (iunout,'(a14,i8)') 'NALS        = ',NALS

      WRITE (iunout,'(a14,i8)') 'NAIN        = ',NAIN
      WRITE (iunout,'(a14,i8)') 'NCPV        = ',NCPV
      WRITE (iunout,'(a14,i8)') 'NBGK        = ',NBGK
      WRITE (iunout,'(a14,i8)') 'NSD         = ',NSD
      WRITE (iunout,'(a14,i8)') 'NSDW        = ',NSDW
      WRITE (iunout,'(a14,i8)') 'NCV         = ',NCV

      CALL EIRENE_LEER(1)
      WRITE (IUNOUT,*) 'MAX. NO. OF ATOMIC/MOLECULAR "REACTIONS"'
      WRITE (iunout,'(a14,i8)') 'NREAC       = ',NREAC
      WRITE (IUNOUT,*) 'NREC,NREI,NRCX,NREL,NRPI: DETERMINED LATER'
C     WRITE (iunout,'(a14,i8)') 'NREC        = ',NREC
C     WRITE (iunout,'(a14,i8)') 'NREI        = ',NREI
C     WRITE (iunout,'(a14,i8)') 'NRCX        = ',NRCX
C     WRITE (iunout,'(a14,i8)') 'NREL        = ',NREL
C     WRITE (iunout,'(a14,i8)') 'NRPI        = ',NRPI

      CALL EIRENE_LEER(1)
      WRITE (IUNOUT,*) 'SETTING OF STORAGE OPTIMIZATION OPTIONS'
      IF (.NOT.LDEFSTOR) THEN
C  ADJUST SOME DEFAULT STORAGE OPTIMIZATION SETTING
        NOPTIM=N1ST*N2ND*N3RD+NADD     ! = NRAD ??
C       NOPTM1=   ??
      ENDIF
C  TRY TO BE INTELLIGENT:  AUTOMATIC STORAGE REDUCTION....
C  SWITCH OFF SUM OVER STRATA IF THERE IS ONLY ONE STRATUM TO BE CALCULATED
c  TO BE TESTED, E.G. CHECK VARIANCES, AND VARIANCES SUM OVER STRATA
      IF (NSTRAI == 1.AND.NSMSTRA.NE.0) THEN
         WRITE (iunout,*) 'NSMSTRA SET = 0, BECAUSE NSTRAI=1 '
         NSMSTRA = 0
      ENDIF
      CALL EIRENE_LEER(1)

C  OPTIONAL STORAGE/PERFORMANCE HANDLING FLAGS
      WRITE (iunout,'(a14,i8)') 'NOPTIM      = ',NOPTIM
      WRITE (iunout,'(a14,i8)') 'NOPTM1      = ',NOPTM1
      WRITE (iunout,'(a14,i8)') 'NGEOM_USR   = ',NGEOM_USR
      WRITE (iunout,'(a14,i8)') 'NCOUP_INPUT = ',NCOUP_INPUT
      WRITE (iunout,'(a14,i8)') 'NSMSTRA     = ',NSMSTRA
      WRITE (iunout,'(a14,i8)') 'NSTORAM     = ',NSTORAM
      WRITE (iunout,'(a14,i8)') 'NGSTAL      = ',NGSTAL
C 
      CALL EIRENE_LEER(1)
      WRITE (IUNOUT,*) 'SETTING OF CENSUS STORAGE FOR T-DEP. MODE'
cdr  time dependent options: census array size
      WRITE (iunout,'(a14,i8)') 'NPRNL       = ',NPRNL
C
      CALL EIRENE_LEER(2)
      
      return

      contains

!******************************************************************************


      subroutine eirene_browse_block_1 (json,p)

      implicit none
      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(inout) :: p
      logical :: found

      WRITE (iunout,*) '*** 1. DATA FOR OPERATING MODE'

      call json%get(p,'NMODE',nmode,found)
      call json%get(p,'NTIME',ntime,found)
      call json%get(p,'NOPTIM',noptim,found)
      call json%get(p,'NOPTM1',noptm1,found)
      call json%get(p,'NGEOM_USR',ngeom_usr,found)
      call json%get(p,'NCOUP_INPUT',ncoup_input,found)
      call json%get(p,'NSMSTRA',nsmstra,found)
      call json%get(p,'NSTORAM',nstoram,found)
      call json%get(p,'NGSTAL',ngstal,found)
      call json%get(p,'NRTAL',nrtal,found)
      ldefstor = .true.
      nrtal_in = nrtal
      noptim_in = noptim
      nsmstra_in = nsmstra
!pb   call json%get(p,'NREAC_ADD',nreac_add,found)

      call json%get(p,'NLERG',nlerg,found)
     
C  NSTORAM IS REDEFINED, FINALLY EITHER =0  (A&M STORAGE SAVE MODE) 
C                                    OR =9  (FULL A&M STORAGE MODE, =DEFAULT) 
      NSTORAM = MIN(NSTORAM,9)
      IF (NSTORAM < 9) NSTORAM = 0
      NOPTM1 = MAX(NOPTM1,1)
C
      return

      end subroutine eirene_browse_block_1


!******************************************************************************


      subroutine eirene_browse_block_2(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p

      INTEGER, allocatable :: indgrd(:)
      integer :: nr1st, nrplg, npplg, nrknot, ncoor, ntrii, ntet,
     .           np2nd, nt3rd, nttra, nbmlt, nradd, ll
      logical :: nlplg, nlfem, nltet, nlrad, nlpol, nltor, nlmlt, nladd,
     .           nltra
      logical :: found
      character(420) :: casename
      character(kind=CK,len=:), allocatable :: case

      WRITE (iunout,*) '*** 2. DATA FOR VOXEL GRID GENERATION'

      call json%get(p,'INDGRD',indgrd,found)
      call json%get(p,'NLRAD',nlrad,found)
      call json%get(p,'NLPOL',nlpol,found)
      call json%get(p,'NLTOR',nltor,found)
      call json%get(p,'NLMLT',nlmlt,found)
      call json%get(p,'NLADD',nladd,found)
      nlpol_in = nlpol

C
C  RADIAL MESH
C
      call json%get(p,'RADIAL_GRID.NLPLG',nlplg,found)
      call json%get(p,'RADIAL_GRID.NLFEM',nlfem,found)
      call json%get(p,'RADIAL_GRID.NLTET',nltet,found)
      nlplg_in = nlplg
      nlfem_in = nlfem

      call json%get(p,'RADIAL_GRID.NR1ST',nr1st,found)
      nr1st_in = nr1st
      n1st = max(n1st,nr1st)

      if (nlplg) then
        call json%get(p,'RADIAL_GRID.NRPLG',nrplg,found)
        call json%get(p,'RADIAL_GRID.NPPLG',npplg,found)
        nplg = max(nplg,nrplg)
        nppart = max(nppart,npplg)
      end if
      
      call json%get(p,'RADIAL_GRID.NRKNOT',nrknot,found)
      call json%get(p,'RADIAL_GRID.NCOOR',ncoor,found)
      nknot = max(nknot,nrknot)
      ncoord = max(ncoord,ncoor)

      ntrii = 0
      ntet = 0
      if (indgrd(1) <= 5) then
        if (nlfem .or. nltet) then
          call json%get(p,'RADIAL_GRID.CASE',case,found)

          IF (.NOT.FOUND) THEN
            WRITE (IUNOUT,*) ' ERROR IN GEOMETRY SPECIFICATION '
            WRITE (IUNOUT,*)
     .           ' TRIANGLE OR TETRAHEDRON GRID SWITCHED ON '
            WRITE (IUNOUT,*) ' BUT NO CASENAME SPECIFIED '
            CALL EIRENE_EXIT_OWN(1)
          END IF

          LL = LEN_TRIM(CASE)
          CASENAME=REPEAT(' ',420)
          IF (LL > LEN(CASENAME)) THEN
            WRITE (IUNOUT,*) 'CASE NAME TOO LONG'
            CALL EIRENE_EXIT_OWN(1)
          END IF
          CASENAME(1:LL) = CASE(1:LL)
          
          IF (NLFEM) THEN
            CALL EIRENE_FIND_TRIANG_DIM(NR1ST, NTRI, NTRII, NKNOT,
     .                                  NGITT, CASENAME, IFOFF)
          ENDIF

          IF (NLTET) THEN
            CALL EIRENE_FIND_TET_DIM(NR1ST, NTETRA, NTET, NCOORD,
     .                               NGITT, CASENAME, IFOFF)
          ENDIF       
        end if
      end if
      NTRIS=NTRII
      NKNOTS=NKNOT
      N1ST = MAX(N1ST,NR1ST)

C
C  POLOIDAL MESH
C
      call json%get(p,'POLOIDAL_GRID.NP2ND',np2nd, found)
      np2nd_in = np2nd
cdr  storage for 2nd coordinate grid only, if NLPOL=T
      if (nlpol) n2nd = max(n2nd,np2nd)

C
C  TOROIDAL MESH
C
      call json%get(p,'TOROIDAL_GRID.NLTRA',nltra,found)
      call json%get(p,'TOROIDAL_GRID.NT3RD',nt3rd,found)
      call json%get(p,'TOROIDAL_GRID.NTTRA',nttra,found)
      nt3rd_in = nt3rd
      if (nltor.and.nltra) nttra = nt3rd
cdr  storage for 3nd coordinate grid only, if NLTOR=T
      if (nltor) n3rd = max(n3rd,nt3rd)
CDR STORAGE FOR TOROIDAL EFFECTS, EVEN IF
CDR NO TOROIDAL RESOLUTION IS USED
      NTOR = MAX(NTOR,NTTRA)  
C
C  MESH MULTIPLICATION
C
      call json%get(p,'MESH_MULTIPLICATION.NBMLT',nbmlt,found)
      NBMLT=MAX0(1,NBMLT)
      NBMAX=NBMLT
C
C  ADDITIONAL CELLS OUTSIDE STANDARD MESH
C
      if (nladd) then
        call json%get(p,'ADD_CELLS.NRADD',nradd,found)
        NADD = MAX(NADD,NRADD)
      end if

      end subroutine eirene_browse_block_2

!******************************************************************************


      subroutine eirene_browse_block_3a(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p

      logical :: found

      WRITE (iunout,*) '*** 3A. DATA FOR NON-DEFAULT STANDARD SURFACES'

      call json%get(p,'NSTSI',nstsi,found)
      NSTS = MAX(NSTS,NSTSI)
        
      end subroutine eirene_browse_block_3a


!******************************************************************************


      subroutine eirene_browse_block_3b(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p

      integer :: nlimi
      logical :: found

      WRITE (iunout,*) '*** 3B. DATA FOR ADDITIONAL SURFACES'

      call json%get(p,'NLIMI',nlimi,found)
      NLIM = MAX(NLIM,NLIMI)
        
      end subroutine eirene_browse_block_3b


!******************************************************************************


      subroutine eirene_browse_block_4(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p

      type(json_value), pointer :: species, preac, pspecs, pspc 
      integer :: nreaci, natmi, nmoli, nioni, nphoti
      logical :: found
      character(50) :: s1, s2, sss

      WRITE (iunout,*) '*** 4. DATA FOR SPECIES SPECIFICATION AND'
      WRITE (iunout,*) '       ATOMIC PHYSICS MODULE'

      s1 = 'SPECIES_SPEC'
      s2 = '.REACTIONS'

      call json%get(p,s1//s2//'.NREACI',nreaci,found)
      NREAC = MAX(NREAC,NREACI)
      

      call json%get_child(p,'SPECIES_SPEC',species,found)
      call json%get_child(species,'REACTIONS',preac,found)
      call json%get_child(preac,'REAC_SPECS',pspecs,found)
      NREAC_LINES=json%count(pspecs)

!  ATOMS

      WRITE (iunout,*)
     .  '*** 4A. NEUTRAL ATOMS SPECIES CARDS, NATMI SPECIES'

      call json%get_child(species,'ATOMS',pspc,found)
      call json%get(pspc,'NATMI',natmi,found)
      NATM = MAX(NATM,NATMI)

      call eirene_read_spc_block (json,pspc,'A',natmi,lmulti,lmulvi)

!  MOLECULES

      WRITE (iunout,*)
     .  '*** 4B. NEUTRAL MOLECULE SPECIES CARDS, NMOLI SPECIES'

      call json%get_child(species,'MOLECULES',pspc,found)
      call json%get(pspc,'NMOLI',nmoli,found)
      NMOL = MAX(NMOL,NMOLI)

      call eirene_read_spc_block (json,pspc,'M',nmoli,lmulti,lmulvi)

!  TEST IONS

      WRITE (iunout,*)
     . '*** 4C. TEST IONS SPECIES CARDS, NIONI SPECIES'

      call json%get_child(species,'TEST_IONS',pspc,found)
      call json%get(pspc,'NIONI',nioni,found)
      NION = MAX(NION,NIONI)

      call eirene_read_spc_block (json,pspc,'I',nioni,lmulti,lmulvi)

!  PHOTONS

      WRITE (iunout,*) '*4D.   PHOTONS SPECIES CARDS, NPHOTI SPECIES '

      call json%get_child(species,'PHOTONS',pspc,found)
      call json%get(pspc,'NPHOTI',nphoti,found)
      NPHOT = MAX(NPHOT,NPHOTI)

      call eirene_read_spc_block (json,pspc,'PH',nphoti,lmulti,lmulvi)

      nullify(species)
      nullify(preac)
      nullify(pspecs)
      nullify(pspc)
       
      end subroutine eirene_browse_block_4

!******************************************************************************

      subroutine eirene_read_spc_block (json,p,ch,nspc,logt,logv)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p
      character(*), intent(in) :: ch
      integer, intent(in) :: nspc
      logical, intent(inout) :: logt, logv

      type(json_value), pointer :: pspecies, spchild, preac, pir
      integer :: nrc, ibgk, i, j
      logical found
      character(kind=CK,len=:), allocatable :: spname

      call json%get_child(p,'SPECIES',pspecies,found)
      call json%get_child(pspecies,spchild)
      if (.not.associated(spchild)) return
      
      do i = 1, nspc
        call json%get(spchild,'NRC'//ch,nrc,found)
        if (.not.found) cycle

        call json%get_child(spchild,'REACTIONS',preac,found)
        if (.not.found) cycle

        call json%get_child(preac,pir)
        do j = 1, nrc
          if (.not.associated(pir)) exit
cdr:  try to identify if there are so called BKG collisions, input flag IBGK::
cdr:  to be generalized: there may be other reactions, which require multiple Ti profiles
          call json%get(pir,'IBGK'//ch,ibgk,found)
          logt = logt .or. (ibgk /= 0)
          logv = logv .or. (ibgk /= 0)
          call json%get_next(pir,pir)
        end do 
        call json%get_next(spchild,spchild)
        if (.not.associated(spchild)) exit
      end do

      nullify (pspecies)
      nullify (spchild)
      nullify (preac)
      nullify (pir)

      end subroutine eirene_read_spc_block 


!******************************************************************************


      subroutine eirene_browse_block_5(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p

      type(json_value), pointer :: pback, pbulk, pspecies, 
     .                             pspc, pdm, pplas
      integer :: nreaci, natmi, nmoli, nioni, nphoti, i, ico
      INTEGER, allocatable :: indpro(:)
      character(kind=CK,len=:), allocatable :: spname
      logical :: found

      WRITE (iunout,*) '*** 5. DATA FOR PLASMA BACKGROUND'

!  BULK IONS

      WRITE (iunout,*) '*** 5A. BULK ION SPECIES CARDS, NPLSI SPECIES'

      call json%get_child(p,'BACKGROUND',pback,found)
      call json%get_child(pback,'BULK_IONS',pbulk,found)
      call json%get(pbulk,'NPLSI',nplsi,found)
      NPLS = MAX(NPLS,NPLSI)

      call json%get_child(pbulk,'SPECIES',pspecies,found)
      call json%get_child(pspecies,pspc)

      ico = 0
      do i = 1, nplsi
        if (associated(pspc)) then
          call json%get_child(pspc,'DENS_MODEL',pdm,found)
          if (found) ico = ico + 1 
        end if
        call json%get_next(pspc,pspc)
      end do

!  enhance number of reactions by 1 for use by density model
      if (ico > 0) nreac = nreac + 1

!  PLASMA

      WRITE (IUNOUT,*) '*** 5B. PLASMA BACKGROUND DATA'

      call json%get_child(pback,'PLASMA',pplas,found)
      call json%get(pplas,'INDPRO',indpro,found)


cdr to be done: syncronisation of options for Ti and Vi.
cdr these next 2 lines for Ti(ipls)  
!pb   NPLSTI = 1
!pb   IF ((INDPRO(2) < 0) .OR. (MOD(INDPRO(2),100) > 9)) NPLSTI=NPLS

      NPLSTI = NPLS  ! now same as for VI.  good
cdr   IF (MOD(ABS(INDPRO(2)),100) > 9) NPLSTI = 1  this should be here, to syncronize with V
      IF (INDPRO(2)<0) NPLSTI=1  !dr  different still from VI logic.

      IF ((NPLS > 1) .AND. (NPLSTI == 1)) THEN
        WRITE (IUNOUT,*) 'WARNING FROM FIND_PARAM'
        WRITE (IUNOUT,*) 'TIIN PROVIDED FOR ONE SPECIES ONLY'
        WRITE (IUNOUT,*) 'DUE TO INDPRO(2) < 0'  !dr  or:  > 10 ???
        IF (LMULTI) THEN
          WRITE (IUNOUT,*) 'DIMENSION OF TIIN OVERWRITTEN'
          WRITE (IUNOUT,*) 'BECAUSE BGK REACTIONS ARE PRESENT'
          NPLSTI = NPLS
        END IF
        WRITE (IUNOUT,*) ' NPLSTI = ',NPLSTI
      END IF


cdr these next 2 lines for Vi(ipls)
      NPLSV = NPLS
      IF (MOD(ABS(INDPRO(4)),100) > 9) NPLSV = 1

      IF ((NPLS > 1) .AND. (NPLSV == 1)) THEN
        WRITE (IUNOUT,*) 'WARNING FROM FIND_PARAM'
        WRITE (IUNOUT,*) 'V_IN PROVIDED FOR ONE SPECIES ONLY',
     .                   'DUE TO INDPRO(4) > 10'  !dr above, for Ti, we say:  < 0
        IF (LMULVI) THEN
          WRITE (IUNOUT,*) 'DIMENSION OF V_IN ARRAYS OVERWRITTEN',
     .                     'BECAUSE BGK REACTIONS ARE  PRESENT'
          NPLSV = NPLS
        END IF
        WRITE (IUNOUT,*) ' NPLSV = ',NPLSV
      END IF
      
      deallocate(indpro)
      nullify(pback)
      nullify(pbulk)
      nullify(pspecies)
      nullify(pspc)
      nullify(pdm)
      nullify(pplas)

      end subroutine eirene_browse_block_5


!******************************************************************************


      subroutine eirene_browse_block_6(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p

      type(json_value), pointer :: ptom
      character(kind=CK,len=:), allocatable :: path
      logical :: found, nltrim

      WRITE (iunout,*) '*** 6. GENERAL DATA FOR REFLECTION MODEL'

      call json%get(p,'NLTRIM',nltrim,found)

      if (nltrim) then
         call json%get_child(p,'PROJ_ON_MAT',ptom,found)
      else
         found = .false.
      end if

      if (found) then
        NHD6 = json%count(ptom)
      else
        NHD6 = 12
      end if

      nullify(ptom)
        
      end subroutine eirene_browse_block_6


!******************************************************************************


      subroutine eirene_browse_block_7(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p

      type(json_value), pointer :: psrc, pstra, pst, psub, psubs
      integer :: j, noc, nsrfsi, isor, ist, i
      real(dp) :: sorlim, sorind
      integer, allocatable :: indsrc(:)
      logical :: found, fsub
      character(kind=CK,len=:),allocatable :: txt

      WRITE (iunout,*) '*** 7. DATA FOR PRIMARY SOURCES, NSTRAI STRATA'

      call json%get(p,'NSTRAI',nstrai,found)
      NSTRA = MAX(NSTRA,NSTRAI)
      
      call json%get(p,'INDSRC',indsrc,found)
      IF (ANY(INDSRC == 6)) THEN
        DO J=NSTRAI,1,-1
          IF (INDSRC(J) == 6) THEN
            NSTEP = J
            EXIT
          END IF
        END DO
      END IF

      call json%get_child(p,'STRATA',pstra,found)

      if (found) then
!  noc stratum blocks in input file 
!  noc < nstrai if any(indsrc==6)
         noc = json%count(pstra)
         call json%get_child(pstra,pst)
         do j = 1, noc
           if (associated(pst)) then
             call json%get(pst,'TXTSOU',txt,found)
             WRITE (IUNOUT,'(A1,A)') ' ',trim(txt)
             deallocate(txt)
!  check for substrata
             call json%get(pst,'NSRFSI',nsrfsi,found)
             NSRFS = MAX(NSRFS,NSRFSI)

             call json%get_child(pst,'SUBSTRATA',psubs,fsub)
             
             call json%get_child(psubs,psub)
             do i = 1, nsrfsi
                if (associated(psub)) then
!  find maximum step function index
                  call json%get(psub,'SORLIM',sorlim,fsub)
                  call json%get(psub,'SORIND',sorind,fsub)
                  isor = int(sorlim)
                  do while (isor > 0)
                    ist =  mod(isor,10)
                    if ((ist == 4).OR.(ist==5)) 
     .                nstep = max(nstep,int(sorind))
                    isor = isor / 10
                  end do
                  call json%get_next(psub,psub)          
                end if
             end do
             
             call json%get_next(pst,pst)          
           end if
         end do
      end if
        
      deallocate (indsrc)
      nullify(pstra)
      nullify(pst)
      nullify(psub)
      nullify(psubs)

      end subroutine eirene_browse_block_7


!******************************************************************************


      subroutine eirene_browse_block_9(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p

      integer :: nsigvi, nsigsi, nsigci
      logical :: found

      WRITE (iunout,*) '*** 9. DATA FOR STATISTIC AND NONANALOG MODEL'

      WRITE (iunout,*) '       CARDS FOR STANDARD DEVIATION'
      
      if (associated(p)) then
        call json%get(p,'NSIGVI',nsigvi,found)
        call json%get(p,'NSIGSI',nsigsi,found)
        call json%get(p,'NSIGCI',nsigci,found)
      else
        nsigvi = 0
        nsigsi = 0
        nsigci = 0
      end if

      NSD = MAX(NSD,NSIGVI)
      NSDW = MAX(NSDW,NSIGSI)
      NCV = MAX(NCV,NSIGCI)
        
      end subroutine eirene_browse_block_9


!******************************************************************************


      subroutine eirene_browse_block_10(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p

      integer :: nadvi, nclvi, nalvi, nadsi, nalsi, nadspc
      logical :: found

      WRITE (iunout,*)
     . '*** 10. DATA FOR ADDITIONAL TALLIES, COLLISION'
      WRITE (iunout,*)
     . '        ESTIMATORS AND ALGEBRAIC EXPRESSIONS'

      if (associated(p)) then
        call json%get(p,'NADVI',nadvi,found)
        call json%get(p,'NCLVI',nclvi,found)
        call json%get(p,'NALVI',nalvi,found)
        call json%get(p,'NADSI',nadsi,found)
        call json%get(p,'NALSI',nalsi,found)
        call json%get(p,'NADSPC',nadspc,found)
      else
        nadvi = 0
        nclvi = 0
        nalvi = 0
        nadsi = 0
        nalsi = 0
        nadspc = 0
      end if

      NADV = MAX(NADV,NADVI)
      NCLV = MAX(NCLV,NCLVI)
      NALV = MAX(NALV,NALVI)
      NADS = MAX(NADS,NADSI)
      NALS = MAX(NALS,NALSI)

      WRITE (iunout,*) '*** 10A. DATA FOR ADDITIONAL TALLIES'
      WRITE (iunout,*) '*** 10B. DATA FOR COLLISION ESTIMATORS'
      WRITE (iunout,*) '*** 10C. DATA FOR ALGEBRAIC EXPRESSIONS'
      WRITE (iunout,*) '*** 10D. DATA FOR ADDITIONAL SURFACE TALLIES'
      WRITE (iunout,*) '*** 10E. DATA FOR ALGEBRAIC SURFACE TALLIES'
      WRITE (iunout,*) '*** 10F. DATA FOR SPECTRA'
        
      end subroutine eirene_browse_block_10


!******************************************************************************


      subroutine eirene_browse_block_11(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p

      type(json_value), pointer :: ptals, ptal
      integer :: nvolpr, nsurpr, nvolpl, i, nsp
      logical :: found, lrpscut

      WRITE (iunout,*)
     .  '*** 11. DATA FOR NUMERICAL AND GRAPHICAL OUTPUT'

      call json%get(p,'NVOLPR',nvolpr,found)
      call json%get(p,'NSURPR',nsurpr,found)

      NVLPR=NVOLPR
C  ERGODIC OPTION NEEDS PRINTOUT OF VOLUME, AND ONE, TWO OR THREE FURTHER TALLIES AT LEAST
      IF (NLERG) NVLPR=MAX(4,NVLPR)

      NSRPR=NSURPR
C  ERGODIC OPTION NEEDS PRINTOUT AT LEAST FROM TIME-HORIZON
      IF (NLERG) NSRPR=MAX(1,NSRPR)

      call json%get(p,'LRPSCUT',lrpscut,found)
      call json%get(p,'NVOLPL',nvolpl,found)

      NPLT = 1
      if (nvolpl > 0) then
        call json%get_child(p,'PLOT_TALLIES',ptals,found)
        call json%get_child(ptals,ptal)

        do i = 1, nvolpl
          if (associated(ptal)) then
            call json%get(ptal,'NSP',nsp,found)
            NPLT = MAX(NPLT, NSP)
            call json%get_next(ptal, ptal)
          end if
        end do

        nullify(ptals)
        nullify(ptal)
      end if
       
      end subroutine eirene_browse_block_11


!******************************************************************************


      subroutine eirene_browse_block_12(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p

      type(json_value), pointer :: plines, pline, pcomps, pcomp,
     .                             pcnts, pcnt
      character(kind=CK,len=:), allocatable :: ckey
      integer :: num_lines, mod_addv, i, j, k, num_compo, num_contrib, 
     .           lines, iratio, nadv_add, nchori, ncheni
      logical :: lkey, found, nlemis

      WRITE (iunout,*) '*** 12. DATA FOR DIAGNOSTIC MODULE'

      NADV_ADD = 0
      NLEMIS = .FALSE.
      NUM_LINES = 0
      NCHORI = 0
      NCHENI = 0

      if (associated(p)) then
        call json%get(p,'KEYWORD',ckey,lkey)

        if (lkey) then
          if (index(ckey,'DEFINE_LINES') == 0) lkey = .false.
        end if
      
        if (lkey) then
!  read emission lines
          call json%get(p,'.NUM_LINES',num_lines,found)
          call json%get(p,'.MOD_ADDV',mod_addv, found)
          
          call json%get_child(p,'LINES',plines,found)
          call json%get_child(plines,pline)

!  lines
          do i = 1, num_lines
            if (associated(pline)) then
              call json%get(pline,'NUM_COMPO',num_compo,found)
              IF (MOD_ADDV == 0) THEN
                NADV_ADD = MAX(NADV_ADD,NUM_COMPO+1)
              ELSE 
                NADV_ADD = NADV_ADD + NUM_COMPO + 1
              END IF

              if (.false.) then  ! read components and contributions, not needed
! components
              call json%get_child(pline,'COMPONENTS',pcomps,found)
              call json%get_child(pcomps,pcomp)
              do j = 1, num_compo
                if (associated(pcomp)) then
                  call json%get(pcomp,'NUM_CONTRIB',num_contrib,found)

! contributions
                  call json%get_child(pcomp,'CONTRIBUTIONS',pcnts,found)
                  call json%get_child(pcnts,pcnt)
                  do k = 1, num_contrib
                    if (associated(pcnt)) then
                      call json%get(pcnt,'IRATIO',iratio,found)
                      call json%get_next(pcnt,pcnt)
                    end if
                  end do

                  call json%get_next(pcomp,pcomp)
                end if
              end do
              end if  ! if .false.

              call json%get_next(pline,pline)
            end if
          end do

! ADD 1 FOR TOTAL  --  already done
!         IF (MOD_ADDV == 0) NADV_ADD = NADV_ADD + 1
          NADV = NADV + NADV_ADD 
!         NREAC = NREAC + LINES
         
        end if

        WRITE (iunout,*) '*** 12A. DATA FOR LINES OF SIGHT'

        call json%get(p,'NCHORI',nchori,found)
        call json%get(p,'NCHENI',ncheni,found)

        nullify(plines) 
        nullify(pline)
        nullify(pcomps) 
        nullify(pcomp)
        nullify(pcnts) 
        nullify(pcnt)
 
      end if
      
      NCHOR = MAX(NCHOR,NCHORI)
      NCHEN = MAX(NCHEN,NCHENI)
      NLEMIS = NLEMIS .OR. (NCHOR > 0)
      IF (NLEMIS.AND.(NUM_LINES == 0)) THEN
! USE DEFAULT LINES FOR EMISSIVITY
        MOD_ADDV = 0
        NADV=NADV+7
        NUM_LINES = 6
        NUM_COMPO = 6
! USE MAXIMUM AS NCHAR AND NCHRG ARE NOT YET AVAILABLE
        NUM_CONTRIB = NATMI + NMOLI + 2*NMOLI + 2*NMOLI + 2*NMOLI +
     .                NPLSI
        NREAC = NREAC + NUM_LINES*NUM_COMPO + 3           
      END IF
    
      end subroutine eirene_browse_block_12


!******************************************************************************


      subroutine eirene_browse_block_13(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p

      integer :: nprnli, nprmul, nsnvi
      logical :: found

      WRITE (iunout,*)
     .  '*** 13. DATA FOR ITERATIVE AND TIME DEP. OPTION'

      nprnli = 0
      nprmul = 1
      NSNVI = 0
      
      if (associated(p)) then
        call json%get(p,'NPRNLI',nprnli,found)
        call json%get(p,'DEFAULT_TIME_HORIZON',ldef_time_horizon,found)
        call json%get(p,'NPRMUL',nprmul,found)
        if (nprnli > 0) call json%get(p,'NSNVI',nsnvi,found)
      end if

      IF (NPRMUL > 1) NPRNLI = NPRNLI * NPRMUL
      NPRNL = MAX(NPRNL,NPRNLI)

      NSNV = MAX(NSNV,NSNVI)

      IF (NLERG.AND.NPRNLI.LE.0) THEN
C  NO TIME HORIZON DEFINED, DESPITE NLERG=.TRUE.
C  THEREFORE: SET A DEFAULT TIME HORIZON HERE
        IF (NTIME.EQ.0) NTIME=1
        NPRNLI=100
      ENDIF
      NPRNL = MAX(NPRNL,NPRNLI)
 
      end subroutine eirene_browse_block_13


!******************************************************************************


      subroutine eirene_browse_block_14(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p

      integer :: naini, ncopii, ncopie, ncpvi
      logical :: found

      WRITE (iunout,*) '*** 14. DATA FOR INTERFACING ROUTINE "INFCOP"'

      IF (NMODE.EQ.0) THEN
        call json%get(p,'NAINI',naini,found)
        call json%get(p,'NCOPII',ncopii,found)
        call json%get(p,'NCOPIE',ncopie,found)
        NCPVI=NCOPIE
      ELSE
        NAINI=0
        NCPVI=0
!pb     CALL EIRENE_IF0PRM(IUNIN)
        j = itree_num(14)
        call eirene_if0prm_json(json,p)
      ENDIF
      
cdr  some parameters may have gotten changed in IF0PRM, case-specific
      NAIN = MAX(NAIN,NAINI)
      NCPV = MAX(NCPV,NCPVI)
        
      end subroutine eirene_browse_block_14

 
!******************************************************************************


      END SUBROUTINE EIRENE_FIND_PARAM_JSON
