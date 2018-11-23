C ===== SOURCE: prousr.f
CDK USER
C
C   USER SUPPLIED SUBROUTINES
C
C           ************
C           *  EDGE2D  *  (fem-interface)
C           ************
C
C s.wiesen@fz-juelich.de (mar06)
C
C get plasma values for bulk species IPLS (module comprt, set by subr. plasma)
c
c reads casename.zplasma
c
      SUBROUTINE EIRENE_PROUSR (PRO,INDX,P0,P1,P2,P3,P4,P5,PROVAC,N)
C
C  EXAMPLE:
C     P0 : CENTRAL VALUE
C     P1 : STARTING RADIUS FOR POLYNOMIAL
C     P2 : SWITCH FOR PHASE 1: OH-PHASE
C                           2: NI-PHASE
C     P3 : FACTOR FOR TI: TI=K*TE
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      USE EIRMOD_CINIT
      USE EIRMOD_CCONA
      USE EIRMOD_COMPRT

      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: P0, P1, P2, P3, P4, P5, PROVAC
      REAL(DP), INTENT(OUT) :: PRO(*)
      INTEGER, INTENT(IN) :: INDX, N
      CHARACTER(256) :: FILENAME, line,sstr
      INTEGER :: NTR, I, J, LL,ier, ipl
      REAL(DP), ALLOCATABLE, SAVE :: PLAS_INDEP(:,:)
      REAL(DP), ALLOCATABLE, SAVE :: PLAS_DEP(:,:,:)
      REAL(DP) :: bnorm
      INTEGER, ALLOCATABLE, SAVE :: INDAR_indep(:)
      INTEGER, ALLOCATABLE, SAVE :: INDAR_dep(:)
      integer, allocatable, save :: inmass(:), inchar(:), inchrg(:)
      integer, parameter :: fp=31
      character(2) :: cstr2
      integer :: idum
      real(dp) :: rdum
csw
c     reset?
      if(indx < 0) then
         if(allocated(plas_indep)) deallocate(plas_indep)
         if(allocated(plas_dep)) deallocate(plas_dep)
         if(allocated(indar_indep)) deallocate(indar_indep)
         if(allocated(indar_dep)) deallocate(indar_dep)
         if(allocated(inmass)) deallocate(inmass)
         if(allocated(inchar)) deallocate(inchar)
         if(allocated(inchrg)) deallocate(inchrg)
c just in case...:
         close(fp+ifoff)
         return
      endif
csw

c read in plasma data from fort.31 ?
      if (.not.allocated(plas_indep)) then
        allocate(plas_indep(11,nrad))
        allocate(plas_dep(5,nrad,npls))
        plas_indep=0.
        plas_dep=0.
        allocate (indar_indep(11))
        allocate (indar_dep(5))
        indar_indep=0
        indar_dep=0
        allocate (inmass(npls))
        allocate (inchar(npls))
        allocate (inchrg(npls))
        inmass = 0
        inchar = 0
        inchrg = 0

        ll=len_trim(casename)
        filename=casename(1:ll) // '.plasma'
        open (unit=fp+ifoff,file=filename,access='sequential',
     .        form='formatted')


c misc plasma data:
        write(sstr,'(a20)')
     .          '*** MISC PLASMA DATA'
        CALL EIRENE_locstr_usr(fp+ifoff,sstr,ier)
        if(ier /=0) then
           write(*,*) 'PROUSR: ',trim(sstr),' not found'
           close(fp+ifoff)
           call EIRENE_exit_own(1)
        endif
        read(fp+ifoff,'(a)') line
        do while (line(1:1) == '*')
           read(fp+ifoff,'(a)') line
        enddo

        read(line,*) ntr
        if (ntr /= nr1st-1) then
           write (*,*) 'PROUSR:', trim(sstr)
           write (*,*) ' wrong number of triangles in plasma file'
           write (*,*) ' check for correct number in file ',filename
           call EIRENE_exit_own(1)
        endif

        do j=1,ntr
c     .       i,te,ne,bx,by,bz,ex,ey,ez,pot,psi
!pb           read(fp+ifoff,'(i7,7(1x,e14.7))') idum,plas_indep(1:10,j)
           read(fp+ifoff,*) idum,plas_indep(1:10,j)

           bnorm=sqrt(plas_indep(3,j)**2 + plas_indep(4,j)**2 +
     .                plas_indep(5,j)**2)
           if (bnorm < eps12) then
              plas_indep(3,j) = 0._dp
              plas_indep(4,j) = 0._dp
              plas_indep(5,j) = 1._dp
              plas_indep(11,j) = 0._dp
           else
              plas_indep(3:5,j) = plas_indep(3:5,j)/bnorm
              plas_indep(11,j) =  bnorm
           endif
        enddo
! default BZ to 1
        plas_indep(5,ntr+1:nrad) = 1._dp
! default BF to 1
        plas_indep(11,ntr+1:nrad) = 1._dp

c loop over species:
        do i=1,nplsi
           write(cstr2,'(i2.2)') i
           write(sstr,'(a23)')
     .          '*** ION #'//cstr2//' PLASMA DATA'
           CALL EIRENE_locstr_usr(fp+ifoff,sstr,ier)
           if(ier /=0) then
              write(*,*) 'PROUSR: tag ',trim(sstr),' not found'
!pb              close(fp+ifoff)
!pb              call EIRENE_exit_own(1)
              write(*,*) 'using zero plasma and fluxes'
              write (iunout,*)
              cycle
           endif

           read(fp+ifoff,'(a)') line
           do while (line(1:1) == '*')
             read(fp+ifoff,'(a)') line
           enddo

           read(line,*) inmass(i),inchar(i),inchrg(i)

           read(fp+ifoff,*) ntr
           if (ntr /= nr1st-1) then
              write (*,*) 'PROUSR:', trim(sstr)
              write (*,*) ' wrong number of triangles in plasma file'
              write (*,*) ' check for correct number in file ',filename
              call EIRENE_exit_own(1)
           endif

           do j=1,ntr
c     .             ti,ni,vx,vy,vz
!pb              read(fp+ifoff,'(i7,5(1x,e14.7))') idum,plas_dep(1:5,j,i)
              read(fp+ifoff,*) idum,plas_dep(1:5,j,i)
           enddo
        enddo

        close(fp+ifoff)
      endif

c
c fort.31 read in now
c
c sort in values

      if ((indx == 1)        .or. (indx == 1+1*npls) .or.
     .    (indx == 1+2*npls) .or. (indx == 1+3*npls) .or.
     .    (indx == 1+4*npls)) then

        ipl = 0
        do i=1,nplsi

           if (nmassp(ipls) /= inmass(i)) cycle
           if (ncharp(ipls) /= inchar(i)) cycle
           if (nchrgp(ipls) /= inchrg(i)) cycle

           ipl = i
           exit
        end do

        if (ipl == 0) then
          write (iunout,*) ' NO MATCHING SPECIES FOUND IN PLASMA FILE'
          write (iunout,*) ' FOR IPLS = ',ipls, ' INDEX = ',indx
          write (iunout,*) ' PROFILE IS SET TO 0 '
          write (iunout,*)
          pro(1:n) = 0._dp
          return
        end if
      end if
c
      if(indx == 0) then
c te
         pro(1:n) = plas_indep(1,1:n)
         indar_indep(1) = indar_indep(1)+1

      elseif (indx == 1) then
c ti
         pro(1:n) = plas_dep(1,1:n,ipl)
         indar_dep(1) = indar_dep(1)+1

      elseif (indx == 1+1*npls) then
c ni
         pro(1:n) = plas_dep(2,1:n,ipl)
         indar_dep(2) = indar_dep(2) + 1

      elseif (indx == 1+2*npls) then
c vx
         pro(1:n) = plas_dep(3,1:n,ipl)
         indar_dep(3) = indar_dep(3) + 1

      elseif (indx == 1+3*npls) then
c vy
         pro(1:n) = plas_dep(4,1:n,ipl)
         indar_dep(4) = indar_dep(4) + 1

      elseif (indx == 1+4*npls) then
c vz
         pro(1:n) = plas_dep(5,1:n,ipl)
         indar_dep(5) = indar_dep(5) + 1

      elseif (indx == 1+1*npls+nplsti+3*nplsv) then
c bx
         pro(1:n) = plas_indep(3,1:n)
         indar_indep(3) = indar_indep(3) + 1

      elseif (indx == 2+1*npls+nplsti+3*nplsv) then
c by
         pro(1:n) = plas_indep(4,1:n)
         indar_indep(4) = indar_indep(4) + 1

      elseif (indx == 3+1*npls+nplsti+3*nplsv) then
c bz
         pro(1:n) = plas_indep(5,1:n)
         indar_indep(5) = indar_indep(5) + 1

      elseif (indx == 4+1*npls+nplsti+3*nplsv) then
c bf
         pro(1:n) = plas_indep(11,1:n)
         indar_indep(11) = indar_indep(11) + 1

      elseif (indx == 5+5*npls) then
! psi
         pro(1:n) = plas_indep(10,1:n)

      else
         write (iunout,*) ' prousr: no data provided for index ',indx
         write (iunout,*)
         pro(1:n) = 0._dp
      endif

      return
      end
