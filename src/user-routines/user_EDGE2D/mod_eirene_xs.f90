module mod_eirene_xs
  implicit none

  public :: eirene_xs_read_species, eirene_xs_dealloc_species
  public :: eirene_xs_init, eirene_xs_kill, eirene_xs_linkdb
  public :: eirene_xs_getrc, eirene_xs_isadas

  integer, parameter :: fp1=4977,fp2=4978,fpdb=4979
  integer, public, parameter :: MAXNRC=10

  real*8, save, public :: rca(MAXNRC),rcm(MAXNRC),rcp(MAXNRC)
  integer, save, public :: nnrca,nnrcm,nnrcp
  logical, save, public :: linit
  logical, save, public :: lreadspecies=.false.

  private

  ! eirene.species
  integer, save :: nfuel,nimps
  integer, save, allocatable :: ihfuel(:),itypfuel(:),ispcfuel(:),iplsfuel(:)
  integer, save, allocatable :: izzfuel(:), izzimps(:),immfuel(:), immimps(:)
  real*8, save, allocatable :: factfuel(:)
  integer, save, allocatable :: izimps(:),itypimps(:),ispcimps(:),iplsimps(:)
  real*8, save, allocatable :: factimps(:)
  integer, save, allocatable :: ih2iatm(:),ih2imol(:),ih2iion(:),ih2ipls(:)
  integer, save, allocatable :: iz2iatm(:),iz2imol(:),iz2iion(:),iz2ipls(:)

  integer, parameter :: MAXIZZ = 10
  integer, save, allocatable, public :: ih2izz(:),iz2izz(:), izz2ih(:), izz2iz(:)

  !cc


  integer :: natm,nmol,npls

  integer, save, allocatable :: ireaca(:,:),ireacm(:,:),ireacp(:,:)
  integer, save, allocatable :: nmassa(:),nchara(:),nmassm(:),ncharm(:),nmassp(:),ncharp(:)
  integer, save, allocatable :: nrca(:),nrcm(:),nrcp(:),massp(:),masst(:)
  character(len=8), save, allocatable :: texta(:),textm(:),textp(:)  

  integer, save, allocatable :: iftflg(:,:), iadas(:),iswr(:),modclf(:)
  real*8, save, allocatable :: creac(:,:,:)

  real*8, parameter :: tvac=0.02, dvac=1.e2, eps60=1.e-60
  real*8, parameter :: PMASSA=1.0073, AMUA=1.6606D-24, ELCHA=1.6022D-19
  real*8, save :: amuakg,cvel2a,cvelaa,cveli2

  logical, parameter :: leirxs_debug=.false.
  integer, parameter,public :: ifeirxs=598

contains
  logical function eirene_xs_isadas(ir) result(lres)
    implicit none
    integer, intent(in) :: ir
    if(iadas(ir) /= 0) then
      lres = .true.
    else
      lres=.false.
    endif
    return
  end function eirene_xs_isadas

  subroutine eirene_xs_init
    implicit none
    integer :: idum,iatm,imol,ipls,kk,k,ier,numsec,i
    real*8 :: rdum
    character(len=256) :: sstr


    ! constants
    AMUAKG=AMUA*1.D-3
    CVEL2A=SQRT(1.D4*ELCHA/AMUAKG)
    CVELAA=sqrt(2.d0)*CVEL2A
    CVELI2=1.d0/CVELAA/CVELAA

    open (unit=fp1,file='eirene.input',access='sequential',form='formatted')

    ! read atoms, block4
    write(sstr,'(a7)') '** 4a '
    call locstr(fp1,sstr,ier)    
    if(ier /=0) call serror(sstr)
    call skipasterisk(fp1,ier)
    if(ier /=0) call serror('skipasterisks')

    read(fp1,*) natm
    allocate(texta(natm))
    allocate(nmassa(natm))
    allocate(nchara(natm))
    allocate(nrca(0:natm))
    nrca=0
    allocate(ireaca(natm,MAXNRC))

    do iatm=1,natm
       read (fp1,'(I2,1X,A8,12(I3),1X,A10,1X,I2)') i,texta(iatm),nmassa(iatm),nchara(iatm),idum,idum,idum,idum,numsec,nrca(iatm),idum,idum,idum
       if(numsec == 3) call serror('numsec==3')
       if(nrca(iatm) > MAXNRC) call serror('increase MAXNRC')

       if(nrca(iatm) < 0) call serror(' nrc<0 not supported yet')

       do k=1,nrca(iatm)
          READ (fp1,'(12I6)') IREACA(IATM,K), idum, idum, idum, idum, idum, idum
          read (fp1,'(6e12.4)') rdum,rdum,rdum,rdum,rdum,rdum
          
          kk = ireaca(iatm,k)
          call read_reac(kk)
       enddo
    enddo


    ! read molecules, block4
    write(sstr,'(a7)') '** 4b '
    call locstr(fp1,sstr,ier)    
    if(ier /=0) call serror(sstr)
    call skipasterisk(fp1,ier)
    if(ier /=0) call serror('skipasterisks')

    read(fp1,*) nmol
    allocate(textm(nmol))
    allocate(nmassm(nmol))
    allocate(ncharm(nmol))
    allocate(nrcm(0:nmol))
    nrcm=0
    allocate(ireacm(nmol,MAXNRC))

    do imol=1,nmol
       read (fp1,'(I2,1X,A8,12(I3),1X,A10,1X,I2)') i,textm(imol),nmassm(imol),ncharm(imol),idum,idum,idum,idum,numsec,nrcm(imol),idum,idum,idum
       if(numsec == 3) call serror('numsec==3')
       if(nrcm(imol) > MAXNRC) call serror('increase MAXNRC')

       if(nrcm(imol) < 0) call serror(' nrc<0 not supported yet')

       do k=1,nrcm(imol)
          READ (fp1,'(12I6)') IREACm(Imol,K), idum, idum, idum, idum, idum, idum
          read (fp1,'(6e12.4)') rdum,rdum,rdum,rdum,rdum,rdum
          
          kk = ireacm(imol,k)
          call read_reac(kk)
       enddo
    enddo

    ! read bulk particles, block5
    write(sstr,'(a7)') '*** 5. '
    call locstr(fp1,sstr,ier)    
    if(ier /=0) call serror(sstr)
    call skipasterisk(fp1,ier)
    if(ier /=0) call serror('skipasterisks')

    read(fp1,*) npls
    allocate(textp(npls))
    allocate(nmassp(npls))
    allocate(ncharp(npls))
    allocate(nrcp(0:npls))
    nrcp=0
    allocate(ireacp(npls,MAXNRC))

    do ipls=1,npls
       read (fp1,'(I2,1X,A8,12(I3),1X,A10,1X,I2)') i,textp(ipls),nmassp(ipls),ncharp(ipls),idum,idum,idum,idum,numsec,nrcp(ipls),idum,idum,idum
       if(numsec == 3) call serror('numsec==3')
       if(nrcp(ipls) > MAXNRC) call serror('increase MAXNRC')

       if(nrcp(ipls) < 0) call serror(' nrc<0 not supported yet')

       do k=1,nrcp(ipls)
          READ (fp1,'(12I6)') IREACp(Ipls,K), idum, idum, idum, idum, idum, idum
          read (fp1,'(6e12.4)') rdum,rdum,rdum,rdum,rdum,rdum
          
          kk = ireacp(ipls,k)
          call read_reac(kk)          
       enddo
    enddo

    close(fp1)
    linit=.true.
  end subroutine eirene_xs_init


  subroutine eirene_xs_kill
    implicit none
    if(allocated(ireaca)) deallocate(ireaca)
    if(allocated(ireacm)) deallocate(ireacm)
    if(allocated(ireacp)) deallocate(ireacp)
    if(allocated(nmassa)) deallocate(nmassa)
    if(allocated(nchara)) deallocate(nchara)
    if(allocated(nmassm)) deallocate(nmassm)
    if(allocated(ncharm)) deallocate(ncharm)
    if(allocated(nmassp)) deallocate(nmassp)
    if(allocated(ncharp)) deallocate(ncharp)
    if(allocated(nrca)) deallocate(nrca)
    if(allocated(nrcm)) deallocate(nrcm)
    if(allocated(nrcp)) deallocate(nrcp)
    if(allocated(texta)) deallocate(texta)
    if(allocated(textm)) deallocate(textm)
    if(allocated(textp)) deallocate(textp)
    if(allocated(iftflg)) deallocate(iftflg)
    if(allocated(creac)) deallocate(creac)
    if(allocated(iadas)) deallocate(iadas)
    if(allocated(iswr)) deallocate(iswr)
    if(allocated(modclf)) deallocate(modclf)

    call eirene_xs_dealloc_species
    linit=.false.
  end subroutine eirene_xs_kill

  subroutine eirene_xs_linkdb(ierr)
    implicit none
    integer,intent(out) :: ierr
    character*256 :: path
    character(len=80) :: HOME='/u/sim/'
    character(len=80) :: AM_DB='edge2d/data/eirene/AMdata/'
    character(len=80) :: SF_DB='edge2d/data/eirene/Surfacedata/'
    character(len=80) :: ADD_DB='edge2d/data/eirene/Additionals/'
    character(len=80) :: ADAS_DB='edge2d/data/adas'
    logical :: lex    
    
    inquire(file='eirene.input', exist=lex)
    if(.not.lex) then
       ierr = 1
       return
    endif
    

    !c     A&M data
    write(path,'(a)') trim(home)//trim(am_db)
    call lfile (trim(path)//'amjuel.tex', 'AMJUEL')
    call lfile (trim(path)//'h2vibr.tex', 'H2VIBR')
    call lfile (trim(path)//'hydhel.tex', 'HYDHEL')
    call lfile (trim(path)//'methane.tex', 'METHANE')
    call lfile (trim(path)//'spectral.tex', 'SPECTR')         
    write(path,'(a)') trim(home)//trim(adas_db)
    call lfile (trim(path), 'ADAS')

    !c     Surface data
    write(path,'(a)') trim(home)//trim(sf_db)
    call lfile (trim(path)//'SPUTER', 'SPUTER')
    call lfile (trim(path)//'TRIM/trim.dat', 'fort.521')
         
    !c     additional data (photons, etc..)
    write(path,'(a)') trim(home)//trim(add_db)
    call lfile (trim(path)//'PHOTON', 'PHOTON')
    call lfile (trim(path)//'POLARI', 'POLARI')
    call lfile (trim(path)//'graphite_ext.dat', 'graphite_ext.dat')
    call lfile (trim(path)//'mo_ext.dat', 'mo_ext.dat')         

    ierr = 0
  end subroutine eirene_xs_linkdb
  

  subroutine lfile(src,dst)
    implicit none
    character(len=*), intent(in) :: src,dst
    integer :: ier
    logical :: lex
    integer, external :: link,unlink,symlnk,getcwd
    
    if(leirxs_debug) then
       write(ifeirxs,'(a,a,a)') 'eirene_link_file: ',trim(src),' --> ',trim(dst)
    endif
    
    inquire(file=trim(src), exist=lex)
    if(.not.lex) then
       write(*,*) 'mod_eirene_xs: error: src does not exist:'
       write(*,*) src
       stop
    endif
    
    
    ier = unlink(trim(dst))
    ier = symlnk(trim(src),trim(dst))
    if(ier /= 0) then
       write(*,*) 'mod_eirene_xs: link error ', ier
       stop
    endif
  end subroutine lfile
  

  subroutine read_reac(kk)
    implicit none
    integer, intent(in) :: kk
    character(8) :: filnam
    character(4) :: h123
    character(9) :: reac
    character(3) :: crc
    integer :: ir,idum,nreac,ier,i
    real*8 :: rdum
    logical :: found
    character(256) :: sstr,line
    integer, external :: unlink,symlnk
    integer :: mp,mt

    ier = unlink(trim('temp'))
    ier = symlnk(trim('eirene.input'),trim('temp'))
    open (unit=fp2,file='temp',access='sequential',form='formatted')
    write(sstr,'(a7)') '*** 4. '
    call locstr(fp2,sstr,ier)    
    if(ier /=0) call serror(sstr)
    call skipasterisk(fp2,ier)
    if(ier /=0) call serror('skipasterisks')
    read(fp2,*) nreac

    if(.not.allocated(creac)) then
       allocate(iftflg(nreac,0:5))
       allocate(creac(9,-1:9,-10:nreac))
       allocate(iadas(nreac))
       allocate(iswr(nreac))
       allocate(modclf(nreac))
       allocate(massp(nreac))
       allocate(masst(nreac))
       massp=0
       masst=0
       iadas=0
       iswr=0
    endif

    do
       read(fp2,'(a)', iostat=ier) line
       if(ier /= 0) call serror('read_reac')

       if(index(line,'*') == 1) then
          write(*,*) 'read_reac: kk=',kk,' not found'
          call serror('read_reac')
       endif

       READ (line,'(I3,1X,A6,1X,A4,A9,A3,2I3,3E12.4)') IR,FILNAM,H123,REAC,CRC,mp,mt,rdum,rdum,rdum
       if( ir == kk) then
          call slreac_xs (ir,filnam,h123,reac,crc)
          massp(kk) = mp
          masst(kk) = mt
          if(leirxs_debug) then
             write(ifeirxs,'(a,I3,1X,A6,1X,A4,A9,A3,2I3,3E12.4)') ' read:',IR,FILNAM,H123,REAC,CRC
             write(ifeirxs,*) 'ISWR=',iswr(kk)
          endif
          exit
       endif

       ! skip 9 const coeffs ! (todo: CHECK IFTFLG == 10 --> just 1 coeff ?!)
       IF (INDEX(FILNAM,'CONST').NE.0) THEN
          read(fp2,'(6e12.4)',iostat=ier) rdum,rdum,rdum,rdum,rdum,rdum,rdum,rdum,rdum
          if(ier /= 0) call serror('read_reac')
       endif
       
    enddo    
    close(fp2)
    ier = unlink(trim('temp'))
  end subroutine read_reac


  subroutine eirene_xs_getrc(ih, iswrr, te,de, ti,di, v0,vi)
    ! ih: nimbus species, ih > 0: fuel
    !                     ih < 0: impurity
    ! (see eirene.species)
    ! 
    ! iswrr: collision type
    !  1 : EI/DS (el.impact)
    !  2 :
    !  3 : CX (charge exch.)
    !  4 : II (ion impact)
    !  5 : EL (elastic)
    !  6 : RC (recombination)
    !  7 : OT (other)
    !
    ! results:
    ! rca: atomic ratecoeffs <sigma v>   (nnrca)
    ! rcm: molec ratecoeffs <sigma v>    (nnrcm)
    ! rcp: bulk-ion ratecoeffs <sigma v> (nnrcp)
    !
    implicit none
    integer, intent(in) :: ih,iswrr
    real*8, intent(in) :: te,de,ti,di,v0(3),vi(3)
    logical :: limpure
    integer :: iz,iatm,imol,ipls,i,kk,nchar
    real*8 :: res

    ! impurity?
    if(ih < 0) then
       limpure = .true.
       iz = -ih
    else
       limpure = .false.
    endif

    ! get ityp/iatm/imol...
    if(.not.limpure) then
       iatm = ih2iatm(ih)
       imol = ih2imol(ih)
       ipls = ih2ipls(ih)
    else
       iatm = iz2iatm(iz)
       imol = iz2imol(iz)
       ipls = iz2ipls(iz)
    endif

    nnrca = 0
    nnrcm = 0
    nnrcp = 0
    rca=0.
    rcm=0.
    rcp=0.

    ! atoms
    do i=1,nrca(iatm)       
       kk = ireaca(iatm,i)
       if(iswrr /= iswr(kk)) cycle
       nchar = nchara(iatm)
       select case(iswrr)
       case(1)
          call eirene_xs_EI(kk, nchar, te,de, res)
          nnrca = nnrca+1
          rca(nnrca) = res
       case(2)
       case(3)
          call eirene_xs_CX(kk, ipls, ti,di, v0,vi, res)
          nnrca = nnrca+1
          rca(nnrca) = res
       case(4)
       case(5)
       case(6)
       case(7)
       end select
    enddo

    ! molecules
    do i=1,nrcm(imol)
       kk = ireacm(imol,i)
       if(iswrr /= iswr(kk)) cycle
       nchar = ncharm(imol)
       select case(iswrr)
       case(1)
          call eirene_xs_EI(kk, nchar, te,de, res)
          nnrcm = nnrcm+1
          rcm(nnrcm) = res
       case(2)
       case(3)
       case(4)
       case(5)
       case(6)
       case(7)
       end select
    enddo

    ! bulk-ions
    do i=1,nrcp(ipls)
       kk = ireacp(ipls,i)
       if(iswrr /= iswr(kk)) cycle
       nchar = ncharp(ipls)
       select case(iswrr)
       case(1)
       case(2)
       case(3)
       case(4)
       case(5)
       case(6)
          call eirene_xs_RC(kk,ipls, te,de, res)
          nnrcp = nnrcp+1
          rcp(nnrcp) = res
       case(7)
       end select
    enddo

  end subroutine eirene_xs_getrc

  subroutine slreac_xs(ir,filnam,h123,reac,crc)
    implicit none
    INTEGER,      INTENT(IN) :: IR
    CHARACTER(8), INTENT(IN) :: FILNAM
    CHARACTER(4), INTENT(IN) :: H123
    CHARACTER(LEN=*), INTENT(IN) :: REAC
    CHARACTER(3), INTENT(IN) :: CRC
    CHARACTER(6),parameter :: AMJUEL='AMJUEL', HYDHEL='HYDHEL', H2VIBR='H2VIBR'
    CHARACTER(7),parameter :: METHANE='METHANE'
    CHARACTER(11) :: REACSTR
    REAL*8 :: CONST
    INTEGER :: I, IND, J, K, IH, I0P1, I0, IC, IREAC, ISW, INDFF,IFLG, INC, IANF,ier
    CHARACTER(80) :: ZEILE
    CHARACTER(2) :: CHR
    CHARACTER(3) :: CHRL, CHRR
    LOGICAL :: LCONST,LGEMIN,LGEMAX

    if(leirxs_debug) write(ifeirxs,*) ir,filnam,h123,reac,crc

    LGEMIN=.FALSE.
    LGEMAX=.FALSE.
    ISWR(IR)=0
    modclf(ir) = 0
    CONST=0.
    CHR='l0'
    I0=0

    IF (INDEX(CRC,'EI').NE.0.OR.INDEX(CRC,'DS').NE.0) ISWR(ir)=1
    IF (INDEX(CRC,'CX').NE.0) ISWR(ir)=3
    IF (INDEX(CRC,'II').NE.0.OR.INDEX(CRC,'PI').NE.0) ISWR(ir)=4
    IF (INDEX(CRC,'EL').NE.0) ISWR(ir)=5
    IF (INDEX(CRC,'RC').NE.0) ISWR(ir)=6
    IF (INDEX(CRC,'OT').NE.0) ISWR(ir)=7

    IF (INDEX(FILNAM,'AMJUEL').NE.0) THEN
       LCONST=.FALSE.
       OPEN (UNIT=fpdb,FILE='AMJUEL')
    ELSEIF (INDEX(FILNAM,'METHAN').NE.0) THEN
       LCONST=.FALSE.
       OPEN (UNIT=fpdb,FILE='METHANE')
    ELSEIF (INDEX(FILNAM,'HYDHEL').NE.0) THEN
       LCONST=.FALSE.
       OPEN (UNIT=fpdb,FILE='HYDHEL')
    ELSEIF (INDEX(FILNAM,'H2VIBR').NE.0) THEN
       LCONST=.FALSE.
       OPEN (UNIT=fpdb,FILE='H2VIBR')
    ELSEIF (INDEX(FILNAM,'CONST').NE.0) THEN
       LCONST=.TRUE.
    ELSEIF (INDEX(FILNAM,'ADAS').ne.0) THEN
!cswxx       read(filnam(5:6),'(i2)') iadas(ir)
       read(fpdb,'(a)',iostat=ier) zeile
       iadas(ir) = ir
       lconst=.false.
       return
    else
       call serror('slreac_xs, filnam='//filnam)
    ENDIF

    IF (H123(4:4).EQ.' ') THEN
       READ (H123(3:3),'(I1)') ISW
    ELSE
       READ (H123(3:4),'(I2)') ISW
    ENDIF

    REACSTR=REPEAT(' ',11)
    IANF=VERIFY(REAC,' ')
    IF (IANF > 0) THEN
       IREAC=INDEX(REAC(IANF:),' ')-1
       IF (IREAC.LT.0) IREAC=LEN(REAC(IANF:))
       REACSTR(2:IREAC+1)=REAC(IANF:IREAC+IANF-1)
       !C  ADD ONE MORE BLANK, IF POSSIBLE
       IREAC=IREAC+2
    ELSE
       IF (.NOT.LCONST) THEN
          write(*,*) ' NO REACTION SPECIFIED IN REACTION CARD ',IR
          close(fpdb)
          call serror('slreac_xs(1)')
       END IF
    END IF


!C  H.0
    IF (ISW.EQ.0) THEN
       CHR='p0'
       CHRL='pl0'
       CHRR='pr0'
       I0=-1
       MODCLF(IR)=MODCLF(IR)+1
       IFLG=0
!C  DEFAULT POTENTIAL: GENERALISED MORSE
       IFTFLG(IR,IFLG)=2
!C  H.1
    ELSEIF (ISW.EQ.1) THEN
       CHR='a0'
       CHRL='al0'
       CHRR='ar0'
       I0=0
       MODCLF(IR)=MODCLF(IR)+10
       IFLG=1
!C  DEFAULT CROSS SECTION: 8TH ORDER POLYNOM OF LN(SIGMA)
       IFTFLG(IR,IFLG)=0
!C  H.2
    ELSEIF (ISW.EQ.2) THEN
       CHR='b0'
       CHRL='bl0'
       CHRR='br0'
       I0=1
       MODCLF(IR)=MODCLF(IR)+100
       IFLG=2
!C  DEFAULT RATE COEFFICIENT: 8TH ORDER POLYNOM OF LN(<SIGMA V>) FOR E0=0.
       IFTFLG(IR,IFLG)=0
!C  H.3
    ELSEIF (ISW.EQ.3) THEN
       MODCLF(IR)=MODCLF(IR)+200
       I0=1
       IFLG=2
       IFTFLG(IR,IFLG)=0
!C  H.4
    ELSEIF (ISW.EQ.4) THEN
       MODCLF(IR)=MODCLF(IR)+300
       I0=1
       IFLG=2
       IFTFLG(IR,IFLG)=0
!C  H.5
    ELSEIF (ISW.EQ.5) THEN
       CHR='e0'
       CHRL='el0'
       CHRR='er0'
       I0=1
       MODCLF(IR)=MODCLF(IR)+1000
       IFLG=3
       IFTFLG(IR,IFLG)=0
!C  H.6
    ELSEIF (ISW.EQ.6) THEN
       MODCLF(IR)=MODCLF(IR)+2000
       I0=1
       IFLG=3
       IFTFLG(IR,IFLG)=0
!C  H.7
    ELSEIF (ISW.EQ.7) THEN
       MODCLF(IR)=MODCLF(IR)+3000
       I0=1
       IFLG=3
       IFTFLG(IR,IFLG)=0
!C  H.8
    ELSEIF (ISW.EQ.8) THEN
       CHR='h0'
       CHRL='hl0'
       CHRR='hr0'
       I0=1
       MODCLF(IR)=MODCLF(IR)+10000
       IFLG=4
       IFTFLG(IR,IFLG)=0
!C  H.9
    ELSEIF (ISW.EQ.9) THEN
       MODCLF(IR)=MODCLF(IR)+20000
       I0=1
       IFLG=4
       IFTFLG(IR,IFLG)=0
!C  H.10
    ELSEIF (ISW.EQ.10) THEN
       MODCLF(IR)=MODCLF(IR)+30000
       I0=1
       IFLG=4
       IFTFLG(IR,IFLG)=0
!C  H.11
    ELSEIF (ISW.EQ.11) THEN
       CHR='k0'
       CHRL='kl0'
       CHRR='kr0'
       I0=1
       IFLG=5
       IFTFLG(IR,IFLG)=0
!C  H.12
    ELSEIF (ISW.EQ.12) THEN
       I0=1
       IFLG=5
       IFTFLG(IR,IFLG)=0
    ENDIF

! check adas
    if(iadas(ir) /= 0) return         

! constant?
    IF (LCONST) THEN
       IND=INDEX(REACSTR,'FT')
       IF (IND /= 0) THEN
          READ (REACSTR(IND+2:),*) IFTFLG(IR,IFLG)
       END IF

       IF (MOD(IFTFLG(IR,IFLG),100) == 10) THEN
!C
!C  READ ONLY ONE FIT COEFFICIENT FROM INPUT FILE
          READ (fp2,'(6e12.4)') CREAC(1,I0,IR)
       ELSE
!C
!C  READ 9 FIT COEFFICIENTS FROM INPUT FILE
          READ (fp2,'(6e12.4)') (CREAC(IC,I0,IR),IC=1,9)
       END IF
       RETURN


!C
!C  READ FROM DATA FILE
!C
    ELSEIF (.NOT.LCONST) THEN

       call locstr(fpdb,'##BEGIN DATA HERE##',ier)
       if(ier /=0) call serror('slreac_xs(2)')
       
       do
          READ (fpdb,'(A80)',iostat=ier) ZEILE
          if(ier /= 0) call serror('slreac_xs(3)')

          IF (INDEX(ZEILE,H123).ne.0) exit
       enddo
!C
       do 
          READ (fpdb,'(A80)',iostat=ier) ZEILE
          if(ier /=0) call serror('slreac_xs(4)')

          IF (INDEX(ZEILE,'H.').NE.0) call serror('slreac_xs(5)')
          IF (INDEX(ZEILE,'Reaction ').EQ.0.or. INDEX(ZEILE,REACSTR(1:ireac)).EQ.0) then
             cycle
          else
             exit
          endif
       enddo
    ENDIF
!C
!C  SINGLE PARAM. FIT, ISW=0,1,2,5,8,11
    IF (ISW.EQ.0.OR.ISW.EQ.1.OR.ISW.EQ.2.OR.ISW.EQ.5.OR.ISW.EQ.8.OR.ISW.EQ.11) THEN
       IF (.NOT.LCONST) THEN
          do
             READ (fpdb,'(A80)',iostat=ier) ZEILE
             if(ier /=0) call serror('slreac_xs(6)')

             INDFF=INDEX(ZEILE,'fit-flag')
             IF (INDEX(ZEILE,CHR)+INDFF.EQ.0) cycle

             IF (INDFF > 0) THEN
                READ (ZEILE((INDFF+8):80),*) IFTFLG(IR,IFLG)
                cycle
             ENDIF
             exit
          enddo

          IF (MOD(IFTFLG(IR,IFLG),100) == 10) THEN
             IND=INDEX(ZEILE,CHR(1:1))
             CREAC(:,I0,IR)=0.
             READ (ZEILE((IND+2):80),'(E20.12)') CREAC(1,I0,IR)
          ELSE
             DO J=0,2
                IND=0
                DO I=1,3
                   IND=IND+INDEX(ZEILE((IND+1):80),CHR(1:1))
                   READ (ZEILE((IND+2):80),'(E20.12)') CREAC(J*3+I,I0,IR)
                enddo
                READ (fpdb,'(A80)',iostat=ier) ZEILE
                if(ier/=0) call serror('slreac_xs(7)')
             enddo
          END IF
!C
!C  READ ASYMPTOTICS, IF AVAILABLE
!C  I0P1=1 FOR CROSS SECTION
!C  I0P1=2 FOR (WEIGHTED) RATE
!          I0P1=I0+1
!          IF (ISW.EQ.0) then
!             ! NO ASYMPTOTICS FOR POTENTIALS
!             CLOSE (UNIT=fpdb)
!             RETURN
!          endif
!c
!          IF (INDEX(ZEILE,CHRL).NE.0.AND.IFEXMN(IR,I0P1).EQ.0) THEN
!             IND=0
!             DO I=1,3
!                INC=INDEX(ZEILE((IND+1):80),CHR(1:1))
!                IF (INC.GT.0) THEN
!                   IND=IND+INDEX(ZEILE((IND+1):80),CHR(1:1))
!                   READ (ZEILE((IND+3):80),'(E20.12)') FPARM(IR,I,I0P1)
!                ENDIF
!             enddo
!             LGEMIN=.true.
!             READ (fpdb,'(A80)',iostat=ier) ZEILE
!             if(ier/=0) call serror('slreac_xs')
!          ENDIF
!c
!          IF (INDEX(ZEILE,CHRR).NE.0.AND.IFEXMX(IR,I0P1).EQ.0) THEN
!             IND=0
!             DO I=4,6
!                INC=INDEX(ZEILE((IND+1):80),CHR(1:1))
!                IF (INC.GT.0) THEN
!                   IND=IND+INDEX(ZEILE((IND+1):80),CHR(1:1))
!                   READ (ZEILE((IND+3):80),'(E20.12)') FPARM(IR,I,I0P1)
!                ENDIF
!             enddo
!             LGEMAX=.true.
!             READ (fpdb,'(A80)',iostat=ier) ZEILE
!             if(ier/=0) call serror('slreac_xs')
!          ENDIF
!c
!          if (lgemin.and.ifexmn(ir,I0P1).eq.0) then
!             IND=INDEX(ZEILE,'=')
!             READ (ZEILE((IND+2):80),'(E12.5)') rcmn(IR,I0P1)
!             rcmn(ir,I0P1)=log(rcmn(ir,I0P1))
!             ifexmn(ir,I0P1)=5
!             READ (fpdb,'(A80)',iostat=ier) ZEILE
!             if(ier/=0) call serror('slreac_xs')
!          endif
!c
!          if (lgemax.and.ifexmx(ir,I0P1).eq.0) then
!             IND=INDEX(ZEILE,'=')
!             READ (ZEILE((IND+2):80),'(E12.5)') rcmx(IR,I0P1)
!             rcmx(ir,I0P1)=log(rcmx(ir,I0P1))
!             ifexmx(ir,I0P1)=5
!             READ (fpdb,'(A80)',iostat=ier) ZEILE
!             if(ier/=0) call serror('slreac_xs')
!          endif
!
!C
!C  ANY OTHER ASYMPTOTICS INFO ON FILE?  SEARCH FOR Tmin, or Emin
!          IF ((INDEX(ZEILE,'Tmin').NE.0.and.I0P1==2).or.
!             .        (INDEX(ZEILE,'Emin').NE.0.and.I0P1==1)) then
!             IND=INDEX(ZEILE,'n')
!             READ (ZEILE((IND+2):80),'(E9.2)') rcmn(IR,I0P1)
!             rcmn(ir,I0P1)=log(rcmn(ir,I0P1))
!C  extrapolation from subr. CROSS
!             if (I0P1.eq.1.and.iswr(ir).eq.1) ifexmn(ir,1)=1
!             if (I0P1.eq.1.and.iswr(ir).eq.3) ifexmn(ir,1)=-1
!             if (I0P1.eq.1.and.iswr(ir).eq.5) ifexmn(ir,1)=-1
!C  extrapolation from subr. CDEF
!C   ??      if (I0PT.eq.2) ifexmn(ir,1)=-1
!             READ (fpdb,'(A80)',iostat=ier) ZEILE
!             if(ier/=0) call serror('slreac_xs')
!          ENDIF

!C      ELSEIF (LCONST) THEN
!C      NOTHING TO BE DONE
       ENDIF


!C
!C  TWO PARAM. FIT, ISW=3,4,6,7,9,10,12
    ELSEIF (ISW.EQ.3.OR.ISW.EQ.4.OR.ISW.EQ.6.OR.ISW.EQ.7.OR.ISW.EQ.9.OR.ISW.EQ.10.OR.ISW.EQ.12) THEN
       l1: DO J=0,2

          do
             READ (fpdb,'(A80)',iostat=ier) ZEILE
             if(ier/=0) call serror('slreac_xs(8)')
             INDFF=INDEX(ZEILE,'fit-flag')
             IF (INDEX(ZEILE,'Index')+INDFF.EQ.0) cycle
             IF (INDFF > 0) THEN
                READ (ZEILE((INDFF+8):80),*) IFTFLG(IR,IFLG)
                cycle
             ENDIF
             exit
          enddo

          READ (fpdb,'(1X)')
          IF (MOD(IFTFLG(IR,IFLG),100) == 10) THEN
             CREAC(:,:,IR)=0.
             READ (fpdb,*) IH,CREAC(1,1,IR)
             EXIT l1
          ELSE
             DO  I=1,9
                READ (fpdb,*) IH,(CREAC(I,K,IR),K=J*3+1,J*3+3)
             enddo
          END IF
       enddo l1
!C   NO ASYMPTOTICS AVAILABLE YET
!C
    ENDIF

    CLOSE (UNIT=fpdb)
    RETURN
  end subroutine slreac_xs

  subroutine serror(sstr)
    implicit none
    character(len=*), intent(in) :: sstr
    write(*,*) 'ERROR eirene_xs: ',trim(adjustl(sstr))
    close(fp1)
    close(fp2)
    close(fpdb)
    call eirene_xs_kill
    stop
  end subroutine serror

  subroutine eirene_xs_EI(kk,ncharz,te,de,res)
    use mod_adas
    implicit none
    integer, intent(in) :: kk,ncharz
    real*8, intent(in) :: te,de
    real*8, intent(out) :: res
    integer, parameter :: nsbox=1, irei=1
    real*8 :: coun(0:9,nsbox), factkk, cf(9,0:9), fctkkl
    real*8 :: ztei,ztne,teinl(1:nsbox),deinl(1:nsbox), dein(1:nsbox), tein(1:nsbox)
    real*8 :: tabds1(1:irei,1:nsbox), dsub,deimin, pls(1:nsbox),tb
    integer :: j,i
    !csw adas:
    real*8 :: tee
    integer, parameter :: izmax=28
    integer :: ier
    real :: svi,sa4(izmax),rta4(izmax),pt04,pta4(izmax)
    character(len=80) :: adasuser='ADAS'

    ztei=max(tvac,min(te,dble(1.e10)))
    teinl(1)=log(ztei)
    ztne=max(dvac,min(de,dble(1.e20)))
    deinl(1)=log(ztne)
    dein(1) = de
    tein(1) = te

    DSUB=LOG(1.D8)
    DEIMIN=LOG(1.D8)
    DO  J=1,NSBOX
       PLS(J)=MAX(DEIMIN,DEINL(J))-DSUB
    enddo

    factkk=1.
    tabds1=0.


    !C  1.) CROSS SECTION(TE) : NOT NEEDED
    !C

    !C  2.) RATE COEFFICIENT (CM**3/S) * ELECTRON DENSITY (CM**-3)
    !C
    !C  2.A) RATE COEFFICIENT = CONST.
    !C     TO BE WRITTEN

    !C  2.B) RATE COEFFICIENT(TE)
    IF (MYIDEZ(MODCLF(KK),3,5).EQ.1) THEN
       IF (MOD(IFTFLG(KK,2),100) == 10) THEN
          !C  RATE:  (1/S)
          COUN(1,1:NSBOX)=CREAC(1,1,KK)
       ELSE
          !C  RATE COEFFICIENT: (CM^3/S)
          CALL myCDEF (TEINL,1,1,KK,COUN,NSBOX,CF,.TRUE.,.FALSE.,.TRUE.)
       END IF
       
       IF (IFTFLG(KK,2) < 100) THEN
          DO J=1,NSBOX
             TABDS1(IREI,J)=COUN(1,J)*DEIN(J)*FACTKK
          END DO
       ELSE
          DO J=1,NSBOX
             TABDS1(IREI,J)=COUN(1,J)*FACTKK
          END DO
       END IF

       !C     ELSEIF (MYIDEZ(MODCLF(KK),3,5).EQ.2) THEN
       !C  2.C) RATE COEFFICIENT(TE,EBEAM)
       !C  TO BE WRITTEN
       !C       MODCOL(1,2,ISP,1)=2

    ELSEIF (MYIDEZ(MODCLF(KK),3,5).EQ.3) THEN
       !C  2.D) RATE COEFFICIENT(TE,NE)
       !csw check adas?
       if(iadas(kk) == 0) then
          CALL myCDEF (TEINL,1,9,KK,COUN,NSBOX,CF,.FALSE.,.FALSE.,.TRUE.)
          DO  J=1,NSBOX
             TABDS1(IREI,J)=COUN(9,J)
          enddo
          DO  I=8,1,-1
             DO  J=1,NSBOX
                TABDS1(IREI,J)=TABDS1(IREI,J)*PLS(J)+COUN(I,J)
             enddo
          enddo
          FCTKKL=LOG(FACTKK)
          IF (IFTFLG(KK,2) < 100) THEN
             DO J=1,NSBOX
                TB=MAX(dble(-100.),TABDS1(IREI,J)+FCTKKL+DEINL(J))
                TABDS1(IREI,J)=EXP(TB)
             enddo
          ELSE
             DO J=1,NSBOX
                TB=MAX(dble(-100.),TABDS1(IREI,J)+FCTKKL)
                TABDS1(IREI,J)=EXP(TB)
             END DO
          END IF
       else
          do j =1,nsbox
             tee = max(tein(j),tvac)
             CALL ADAS_eirene( ncharz, iadas(kk),89, adasuser,izmax,0, 0.00E+00 , 60       , ifeirxs , 0, sngl(tee)  , sngl(tee), sngl(dein(j))  , 1,(/ 0.E0 /) ,(/ 1.E0 /) ,SVI , SA4(1)   , RTA4(1)  , PT04   , PTA4(1), IER )

             IF( IER.NE.0 ) then
                write(*,*) ' eirene_xs_EI: adas error, ier=',ier
                stop
             endif
             if(dein(j) > 0.) then
                tabds1(irei,j) = svi*dein(j)
             else
                tabds1(irei,j) = 0.
             endif
          enddo
       endif

    ENDIF

    res=tabds1(1,1)/(dein(1)+1.e-20)

    ! NO MOMENTUM LOSSES !
    ! NO ENERGY LOSSES !


    RETURN


  end subroutine eirene_xs_EI

  subroutine eirene_xs_RC(kk,ipls,te,de, res)
    use mod_adas
    implicit none
    integer, intent(in) :: kk,ipls
    real*8, intent(in) :: te,de
    real*8, intent(out) :: res
    integer, parameter :: nsbox=1, irrc=1
    real*8 :: coun(0:9,nsbox), factkk, cf(9,0:9)
    real*8 :: ztei,ztne,teinl(1:nsbox),deinl(1:nsbox), dein(1:nsbox), tein(1:nsbox)
    real*8 :: tabrc1(1:irrc,1:nsbox), dsub,deimin, pls(1:nsbox)
    integer :: j
    !csw adas:
    integer, parameter :: izmax=28
    integer :: ier
    real :: svi,sa4(izmax),rta4(izmax),pt04,pta4(izmax)
    character(len=80) :: adasuser='ADAS'

    ztei=max(tvac,min(te,dble(1.e10)))
    teinl(1)=log(ztei)
    ztne=max(dvac,min(de,dble(1.e20)))
    deinl(1)=log(ztne)
    dein(1) = de
    tein(1) = te

    DSUB=LOG(1.D8)
    DEIMIN=LOG(1.D8)
    DO  J=1,NSBOX
       PLS(J)=MAX(DEIMIN,DEINL(J))-DSUB
    enddo

    factkk=1.
    tabrc1=0.


    !C  1.) CROSS SECTION(TE)
    !C           NOT NEEDED

    !C  2.  RATE COEFFICIENT (CM**3/S) * DENSITY (CM**-3) --> RATE (1/S)
    !C

    !C  2.A) RATE COEFFICIENT = CONST.
    !C           TO BE WRITTEN

    !C  2.B) RATE COEFFICIENT(TE)
    IF (MYIDEZ(MODCLF(KK),3,5).EQ.1) THEN
       IF (MOD(IFTFLG(KK,2),100) == 10) THEN
          !C  RATE:  (1/S)
          COUN(1,1:NSBOX)=CREAC(1,1,KK)
       ELSE
          !C  RATE COEFFICIENT: (CM^3/S)
          CALL myCDEF (TEINL,1,1,KK,COUN,NSBOX,CF,.TRUE.,.FALSE.,.TRUE.)
       END IF
       IF (IFTFLG(KK,2) < 100) THEN
          DO J=1,NSBOX
             TABRC1(IRRC,J)=TABRC1(IRRC,J)+COUN(1,J)*DEIN(J)*FACTKK
          END DO
       ELSE
          DO J=1,NSBOX
             TABRC1(IRRC,J)=TABRC1(IRRC,J)+COUN(1,J)*FACTKK
          END DO
       ENDIF

    !ELSEIF (MYIDEZ(MODCLF(KK),3,5).EQ.2) THEN
       !C  2.C) RATE COEFFICIENT(TE,EBEAM): IRRELEVANT

    ELSEIF (MYIDEZ(MODCLF(KK),3,5).EQ.3) THEN
       !C  2.D) RATE COEFFICIENT(TE,NE)
       !csw check adas?
       if(iadas(kk) == 0) then
          CALL myCDEFN(TEINL,PLS,KK,COUN,NSBOX,CF,.TRUE.,.FALSE.,.TRUE.)
          IF (IFTFLG(KK,2) < 100) THEN
             DO J=1,NSBOX
                TABRC1(IRRC,J)=COUN(1,J)*DEIN(J)*FACTKK
             enddo
          ELSE
             DO J=1,NSBOX
                TABRC1(IRRC,J)=COUN(1,J)*FACTKK
             END DO
          ENDIF
       else
          !csw adas branch                                                       
          do j =1,nsbox
             CALL ADAS_eirene( ncharp(ipls), iadas(kk), 89, adasuser, izmax, 0, 0.00E+00 , 60       , ifeirxs , 0,  sngl(tein(j)) , sngl(tein(j)),  sngl(dein(j)) ,1 ,(/ 0.E0 /) ,(/ 1.E0 /), SVI , SA4(1)  , RTA4(1)  , PT04, PTA4(1), IER )
             
             IF( IER.NE.0 ) then
                write(*,*) ' eirene_xs_RC: adas error, ier=',ier
                stop
             endif
             if(dein(j) > 0.) then
                tabrc1(irrc,j) = rta4(1)*dein(j)
             else
                tabrc1(irrc,j) = 0.
             endif
          enddo
                     
       endif

    ENDIF

    res=tabrc1(1,1)/(dein(1)+1.e-20)

    ! NO MOMENTUM LOSSES !
    ! NO ENERGY LOSSES !

  end subroutine eirene_xs_RC


  subroutine eirene_xs_CX(kk, ipl, ti,di, v0, vi,res)
    implicit none
    integer, intent(in) :: kk,ipl
    real*8, intent(in) :: ti,di,v0(3),vi(3)
    real*8, intent(out) :: res

    integer, parameter :: nsbox=1, ircx=1, nstordt=9

    real*8 :: factkk, pls(nsbox),coun(0:9,nsbox),cf(9,0:9)
    real*8 :: ztii,ztni,tiinl(1:nsbox),diinl(1:nsbox), diin(1:nsbox), tiin(1:nsbox)
    real*8 :: tabcx3(1:ircx,1:nsbox,0:9),defcx(1:ircx),eefcx(1:ircx),sigvcx(1:ircx)
    real*8 :: pmass,tmass,addt,addtl,addcx,elab,cxs,tbcx3(1:nstordt),vrelq,vrel,veffq,veff
    real*8 :: pvelq,fct2,zti,rmassp,fctkkl
    integer :: modcol32,modc,nend,j,ipls,iplsti,iplsv,k,ireac,jan,ii,if8,i
    real*8 :: elb,expo

    factkk = 1.

    ztii=max(tvac,min(ti,dble(1.e10)))
    tiinl(1)=log(ztii)
    ztni=max(dvac,min(di,dble(1.e20)))
    diinl(1)=log(ztni)
    diin(1) = di
    tiin(1) = ti

    pvelq= (v0(1)-vi(1))**2+(v0(2)-vi(2))**2+(v0(3)-vi(3))**2

    !C  FACTOR FOR ROOT MEAN SQUARE SPEED
    rmassp = nmassp(ipl)*pmassa
    FCT2=1./RMASSP*3.*CVEL2A*CVEL2A
    ZTI=FCT2*ZTII

    !C
    !C  TARGET MASS IN <SIGMA*V> FORMULA: MAXW. BULK PARTICLE
    !C  (= PROJECTILE MASS IN CROSS SECTION MEASUREMENT: TARGET AT REST)
    PMASS=MASSP(KK)*PMASSA
    !C  PROJECTILE MASS IN <SIGMA*V> FORMULA: MONOENERG. TEST PARTICLE
    !C  (= TARGET PARTICLE IN CROSS SECTION MEASUREMENT; TARGET AT REST)
    TMASS=MASST(KK)*PMASSA
    ADDT=PMASS/RMASSP
    ADDTL=LOG(ADDT)
    ADDCX = ADDTL

    !C
    !C CROSS SECTION (E-LAB)
    IF (myIDEZ(MODCLF(KK),2,5).EQ.1) THEN
       modcol32 = 3
       IF (FACTKK.NE.1.D0) then
          WRITE (*,*) 'FREAC NOT READY FOR CROSS SECTION IN XSTCX'
          stop
       endif
    ENDIF

    !C
    !C RATE COEFFICIENT
    MODC=myIDEZ(MODCLF(KK),3,5)
    IF (MODC.GE.1.AND.MODC.LE.2) THEN
       modcol32=modc
       !MODCOL(3,2,ISP,IPL)=MODC       
       IF (MODC.EQ.1) NEND=1
       IF (MODC.EQ.2) NEND=NSTORDT
       DO  J=1,NSBOX
          !PLS(J)=TIINL(IPLTI,J)+ADDTL
          PLS(J)=TIINL(J)+ADDTL
       enddo
       IF (MODC.EQ.1) THEN
          CALL myCDEF (PLS,1,NEND,KK,COUN,NSBOX,CF,.TRUE.,.FALSE.,.TRUE.)
          DO  J=1,NSBOX
             !TABCX3(IRCX,J,1)=COUN(1,J)*DIIN(IPL,J)*FACTKK
             TABCX3(IRCX,J,1)=COUN(1,J)*DIIN(J)*FACTKK
          enddo
       ELSEIF (MODC.EQ.2) THEN
          CALL myCDEF (PLS,1,NEND,KK,COUN,NSBOX,CF,.FALSE.,.FALSE.,.TRUE.)
          FCTKKL=LOG(FACTKK)
          DO  J=1,NSBOX
             !TABCX3(IRCX,J,1)=COUN(1,J)+DIINL(IPL,J)+FCTKKL
             TABCX3(IRCX,J,1)=COUN(1,J)+DIINL(J)+FCTKKL
          enddo
       ENDIF
       DO  I=2,NEND
          DO J=1,NSBOX
             TABCX3(IRCX,J,I)=COUN(I,J)
          enddo
       enddo
    ELSE
       !C  NO RATE COEFFICIENT. IS THERE A CROSS SECTION AT LEAST?
       IF (MODCOL32 .NE.3) then
          write(*,*) 'eirene_xs_cx error'
          stop
       endif
    ENDIF
    DEFCX(IRCX)=LOG(CVELI2*PMASS)
    EEFCX(IRCX)=LOG(CVELI2*TMASS)


    !C
    !C  CHARGE EXCHANGE RATE COEFFICIENT OF ATOMS IATM
    !C  WITH BULK IONS OF SPEZIES IPLS=1,NPLSI
    !C
    IPLS=ipl
    IPLSTI=ipl
    IPLSV=ipl
    k=1
    !IPLSTI=MPLSTI(IPLS)
    !IPLSV=MPLSV(IPLS)

    !C
    !C  1.) RATE COEFFICIENT
    !C
    IF (MODCOL32.EQ.1) THEN
       !c  MODEL 1:
       !c  MAXWELLIAN RATE, IGNORE NEUTRAL VELOCITY
       SIGVCX(IRCX)=TABCX3(IRCX,K,1)

    ELSEIF (MODCOL32.EQ.2) THEN
       !C  MODEL 2:
       !C  BEAM - MAXWELLIAN RATE
       !IF (TIIN(IPLSTI,K).LT.TVAC) THEN
       IF (TIIN(K).LT.TVAC) THEN
          !C     HERE: T_I IS SO LOW, THAT ALL ION ENERGY IS IN DRIFT MOTION.
          !C           HENCE: USE BEAM-BEAM RATE INSTEAD.
          VRELQ=PVELQ
          VREL=SQRT(VRELQ)
          ELAB=LOG(VRELQ)+DEFCX(IRCX)
          IREAC=kk
          CXS=eirene_xs_CROSS(ELAB,IREAC,IRCX,'FPATHA CX1 ')
          SIGVCX(IRCX)=CXS*VREL*Diin(k)
          !SIGVCX(IRCX)=CXS*VREL*DENIO(IPLS)
       ELSE
          TBCX3(1:NSTORDT) = TABCX3(IRCX,K,1:NSTORDT)
          JAN=NSTORDT
          IF (JAN < 9) THEN
             !PLS=TIINL(IPLSTI,K)+ADDCX(IRCX,IPLS)
             PLS(k)=TIINL(K)+ADDCX
             DO J=JAN+1,9
                TBCX3(J)=CREAC(9,J,KK)
                DO II=8,1,-1
                   TBCX3(J)=TBCX3(J)*PLS(k)+CREAC(II,J,KK)
                END DO
             END DO
          END IF
          !IF (JAN < 1) TBCX3(1)=TBCX3(1)+DIINL(IPLS,K)
          IF (JAN < 1) TBCX3(1)=TBCX3(1)+DIInl(K)
          !C  MINIMUM PROJECTILE ENERGY: 0.1 EV
          ELB=MAX(-2.3d0,LOG(PVELQ)+EEFCX(IRCX))
          EXPO=TBCX3(9)
          DO  II=1,8
             IF8=9-II
             EXPO=EXPO*ELB+TBCX3(IF8)
          enddo
          SIGVCX(IRCX)=EXP(EXPO)
       ENDIF

    ELSEIF (MODCOL32.EQ.3) THEN
       !C  MODEL 3:
       !C  BEAM - BEAM RATE, BUT WITH EFFECTIVE INTERACTION ENERGY
       VEFFQ=ZTI+PVELQ
       VEFF=SQRT(VEFFQ)
       ELAB=LOG(VEFFQ)+DEFCX(IRCX)
       !IREAC=MODCOL(3,1,NSPH+IATM,IPLS)
       IREAC=kk
       CXS=eirene_xs_CROSS(ELAB,IREAC,IRCX,'FPATHA CX2')
       SIGVCX(IRCX)=CXS*VEFF*DIIN(k)

    ELSEIF (MODCOL32.EQ.4) THEN
       !C  MODEL 4
       !C  BEAM - BEAM RATE, IGNORE THERMAL ION ENERGY
       VRELQ=PVELQ
       VREL=SQRT(VRELQ)
       ELAB=LOG(VRELQ)+DEFCX(IRCX)
       !IREAC=MODCOL(3,1,NSPH+IATM,IPLS)
       IREAC=kk
       CXS=eirene_xs_CROSS(ELAB,IREAC,IRCX,'FPATHA CX3')
       SIGVCX(IRCX)=CXS*VREL*DIIN(k)
    ELSE
       write(*,*) ' eirene_xs_cx error(2)'
       stop
    ENDIF

    res=SIGVCX(ircx)/(diin(1)+1.e-20)

    ! NO ENERGY LOSS RATE !

  end subroutine eirene_xs_CX


  subroutine eirene_xs_read_species(nhs,nzs)
    implicit none
    integer, intent(in) :: nhs,nzs
    integer, parameter :: fp=4999
    integer :: ih,ityp,ispc,ipls,idum,i,iz,izz,imm
    real*8 :: fact
    logical :: lex

    if(lreadspecies) then
       return
    endif

    inquire(file='eirene.species',exist=lex)
    if(.not.lex) then
       lreadspecies=.false.
       return
    endif

    !c read eirene.species (defines which eirene species belongs to which edge2d species)
    !c
    open(unit=fp,file='eirene.species',form='formatted')
    read(fp,'(3i6)') idum,idum,idum
    read(fp,'(i6)') nfuel
    allocate(ihfuel(nfuel))
    allocate(itypfuel(nfuel))
    allocate(ispcfuel(nfuel))
    allocate(iplsfuel(nfuel))
    allocate(izzfuel(nfuel))
    allocate(immfuel(nfuel))
    allocate(factfuel(nfuel))
    do i=1,nfuel
       read(fp,'(6(i6,1x),f6.2)') ih,ityp,ispc,ipls,izz,imm,fact
       ihfuel(i)   = ih
       itypfuel(i) = ityp
       ispcfuel(i) = ispc
       iplsfuel(i) = ipls
       izzfuel(i)  = izz
       immfuel(i)  = imm
       factfuel(i) = fact
    enddo
    read(fp,'(i6)') nimps
    allocate(izimps(nimps))
    allocate(itypimps(nimps))
    allocate(ispcimps(nimps))
    allocate(iplsimps(nimps))
    allocate(izzimps(nimps))
    allocate(immimps(nimps))
    allocate(factimps(nimps))
    do i=1,nimps
       read(fp,'(6(i6,1x),f6.2)') iz,ityp,ispc,ipls,izz,imm,fact
       izimps(i)   = iz
       itypimps(i) = ityp
       ispcimps(i) = ispc
       iplsimps(i) = ipls
       izzimps(i)  = izz
       immimps(i)  = imm
       factimps(i) = fact
    enddo
    close(fp)

    !c make fuel table
    allocate(ih2iatm(nhs))
    allocate(ih2imol(nhs))
    allocate(ih2iion(nhs))
    allocate(ih2ipls(nhs))
    allocate(ih2izz(nhs))
    allocate(izz2ih(MAXIZZ))
    izz2ih = 0
    do ih=1,nhs
       ih2iatm(ih) = 0
       ih2imol(ih) = 0
       ih2iion(ih) = 0
       ih2ipls(ih) = 0
       ih2izz(ih)  = 0
    enddo
    do ih=1,nhs
       do i=1,nfuel
          ispc = ispcfuel(i)
          ipls = iplsfuel(i)
          izz = izzfuel(i)

          if(itypfuel(i) .eq. 1) then
             if( ih2iatm(ih) .gt. 0) then
                if(leirxs_debug) then
                   write(ifeirxs,*) ' eirene_xs_read_species WARN:'
                   write(ifeirxs,*) '  ih2iatm(',ih,') already set'
                   write(ifeirxs,*) '  (use iatm=',ih2iatm(ih),')'
                endif
             else
                ih2iatm(ih) = ispc
             endif             
             if( ih2izz(ih) .gt. 0) then
                if(leirxs_debug) then
                   write(ifeirxs,*) ' eirene_xs_read_species WARN:'
                   write(ifeirxs,*) '  ih2izz(',ih,') already set'
                   write(ifeirxs,*) '  (use izz=',ih2izz(ih),')'
                endif
             else
                ih2izz(ih) = izz
                izz2ih(izz) = ih
             endif             
          endif

          if(itypfuel(i) .eq. 2) then
             if( ih2imol(ih) .gt. 0) then
                if(leirxs_debug) then
                   write(ifeirxs,*) ' eirene_xs_read_species WARN:'
                   write(ifeirxs,*) '  ih2imol(',ih,') already set'
                   write(ifeirxs,*) '  (use imol=',ih2imol(ih),')'
                endif
             else
                ih2imol(ih) = ispc
             endif
             if( ih2izz(ih) .gt. 0) then
                if(leirxs_debug) then
                   write(ifeirxs,*) ' eirene_xs_read_species WARN:'
                   write(ifeirxs,*) '  ih2izz(',ih,') already set'
                   write(ifeirxs,*) '  (use izz=',ih2izz(ih),')'
                endif
             else
                ih2izz(ih) = izz
                izz2ih(izz) = ih
             endif             
          endif

          if(itypfuel(i) .eq. 3) then
             if( ih2iion(ih) .gt. 0) then
                if(leirxs_debug) then
                   write(ifeirxs,*) ' eirene_xs_read_species WARN:'
                   write(ifeirxs,*) '  ih2iion(',ih,') already set'
                   write(ifeirxs,*) '  (use iion=',ih2iion(ih),')'
                endif
             else
                ih2iion(ih) = ispc
             endif
             if( ih2izz(ih) .gt. 0) then
                if(leirxs_debug) then
                   write(ifeirxs,*) ' eirene_xs_read_species WARN:'
                   write(ifeirxs,*) '  ih2izz(',ih,') already set'
                   write(ifeirxs,*) '  (use izz=',ih2izz(ih),')'
                endif
             else
                ih2izz(ih) = izz
                izz2ih(izz) = ih
             endif             
          endif

          if(ipls .gt. 0) then
             if( ih2ipls(ih) .gt. 0) then
                if(leirxs_debug) then
                   write(ifeirxs,*) ' eirene_xs_read_species WARN:'
                   write(ifeirxs,*) '  ih2ipls(',ih,') already set'
                   write(ifeirxs,*) '  (use ipls=',ih2ipls(ih),')'
                endif
             else
                ih2ipls(ih) = ipls
             endif
             if( ih2izz(ih) .gt. 0) then
                if(leirxs_debug) then
                   write(ifeirxs,*) ' eirene_xs_read_species WARN:'
                   write(ifeirxs,*) '  ih2izz(',ih,') already set'
                   write(ifeirxs,*) '  (use izz=',ih2izz(ih),')'
                endif
             else
                ih2izz(ih) = izz
                izz2ih(izz) = ih
             endif             
          endif
          
       enddo
    enddo


    !c make impurity table
    allocate(iz2iatm(nzs))
    allocate(iz2imol(nzs))
    allocate(iz2iion(nzs))
    allocate(iz2ipls(nzs))
    allocate(iz2izz(nzs))
    allocate(izz2iz(MAXIZZ))
    izz2iz = 0
    do iz=1,nzs
       iz2iatm(iz) = 0
       iz2imol(iz) = 0
       iz2iion(iz) = 0
       iz2ipls(iz) = 0
       iz2izz(iz) = 0
    enddo
    do iz=1,nzs
       do i=1,nimps
          ispc = ispcimps(i)
          ipls = iplsimps(i)
          izz = izzimps(i)

          if(itypimps(i) .eq. 1) then
             if( iz2iatm(iz) .gt. 0) then
                if(leirxs_debug) then
                   write(ifeirxs,*) ' eirene_xs_read_species WARN:'
                   write(ifeirxs,*) '  iz2iatm(',iz,') already set'
                   write(ifeirxs,*) '  (use iatm=',iz2iatm(iz),')'
                endif
             else
                iz2iatm(iz) = ispc
             endif
             if( iz2izz(iz) .gt. 0) then
                if(leirxs_debug) then
                   write(ifeirxs,*) ' eirene_xs_read_species WARN:'
                   write(ifeirxs,*) '  iz2izz(',iz,') already set'
                   write(ifeirxs,*) '  (use izz=',iz2izz(iz),')'
                endif
             else
                iz2izz(iz) = izz
                izz2iz(izz) = iz
             endif             
          endif
          
          if(itypimps(i) .eq. 2) then
             if( iz2imol(iz) .gt. 0) then
                if(leirxs_debug) then
                   write(ifeirxs,*) ' eirene_xs_read_species WARN:'
                   write(ifeirxs,*) '  iz2imol(',iz,') already set'
                   write(ifeirxs,*) '  (use imol=',iz2imol(iz),')'
                endif
             else
                iz2imol(iz) = ispc
             endif
             if( iz2izz(iz) .gt. 0) then
                if(leirxs_debug) then
                   write(ifeirxs,*) ' eirene_xs_read_species WARN:'
                   write(ifeirxs,*) '  iz2izz(',iz,') already set'
                   write(ifeirxs,*) '  (use izz=',iz2izz(iz),')'
                endif
             else
                iz2izz(iz) = izz
                izz2iz(izz) = iz
             endif             
          endif
               
          if(itypimps(i) .eq. 3) then
             if( iz2iion(iz) .gt. 0) then
                if(leirxs_debug) then
                   write(ifeirxs,*) ' eirene_xs_read_species WARN:'
                   write(ifeirxs,*) '  iz2iion(',iz,') already set'
                   write(ifeirxs,*) '  (use iion=',iz2iion(iz),')'
                endif
             else
                iz2iion(iz) = ispc
             endif
             if( iz2izz(iz) .gt. 0) then
                if(leirxs_debug) then
                   write(ifeirxs,*) ' eirene_xs_read_species WARN:'
                   write(ifeirxs,*) '  iz2izz(',iz,') already set'
                   write(ifeirxs,*) '  (use izz=',iz2izz(iz),')'
                endif
             else
                iz2izz(iz) = izz
                izz2iz(izz) = iz
             endif             
          endif
          
          if(ipls .gt. 0) then
             if( iz2ipls(iz) .gt. 0) then
                if(leirxs_debug) then
                   write(ifeirxs,*) ' eirene_xs_read_species WARN:'
                   write(ifeirxs,*) '  iz2ipls(',iz,') already set'
                   write(ifeirxs,*) '  (use ipls=',iz2ipls(iz),')'
                endif
             else
                iz2ipls(iz) = ipls
             endif
             if( iz2izz(iz) .gt. 0) then
                if(leirxs_debug) then
                   write(ifeirxs,*) ' eirene_xs_read_species WARN:'
                   write(ifeirxs,*) '  iz2izz(',iz,') already set'
                   write(ifeirxs,*) '  (use izz=',iz2izz(iz),')'
                endif
             else
                iz2izz(iz) = izz
                izz2iz(izz) = iz
             endif             
          endif
          
       enddo
    enddo
    lreadspecies = .true.
  end subroutine eirene_xs_read_species

  subroutine eirene_xs_dealloc_species
    implicit none
    if(allocated(ihfuel)) deallocate(ihfuel)
    if(allocated(itypfuel)) deallocate(itypfuel)
    if(allocated(ispcfuel)) deallocate(ispcfuel)
    if(allocated(iplsfuel)) deallocate(iplsfuel)
    if(allocated(factfuel)) deallocate(factfuel)
    if(allocated(izimps)) deallocate(izimps)
    if(allocated(itypimps)) deallocate(itypimps)
    if(allocated(ispcimps)) deallocate(ispcimps)
    if(allocated(iplsimps)) deallocate(iplsimps)
    if(allocated(factimps)) deallocate(factimps)
    if(allocated(ih2iatm)) deallocate(ih2iatm)
    if(allocated(ih2imol)) deallocate(ih2imol)
    if(allocated(ih2iion)) deallocate(ih2iion)
    if(allocated(ih2ipls)) deallocate(ih2ipls)
    if(allocated(iz2iatm)) deallocate(iz2iatm)
    if(allocated(iz2imol)) deallocate(iz2imol)
    if(allocated(iz2iion)) deallocate(iz2iion)
    if(allocated(iz2ipls)) deallocate(iz2ipls)

    if(allocated(izzfuel)) deallocate(izzfuel)
    if(allocated(izzimps)) deallocate(izzimps)
    if(allocated(immfuel)) deallocate(immfuel)
    if(allocated(immimps)) deallocate(immimps)
    if(allocated(ih2izz)) deallocate(ih2izz)
    if(allocated(iz2izz)) deallocate(iz2izz)
    if(allocated(izz2ih)) deallocate(izz2ih)
    if(allocated(izz2iz)) deallocate(izz2iz)
    
    lreadspecies = .false.
  end subroutine eirene_xs_dealloc_species


  subroutine locstr(fp,sstr,ier)
    implicit none
    integer, intent(in) :: fp
    integer, intent(out) :: ier
    character(*), intent(in) :: sstr
    character(256) :: line,mystr
    integer :: ierror,ilen
    
    ier=1
    rewind(fp)
    write(mystr,'(a)') adjustl(sstr)
    ilen = len_trim(mystr)
    do
       read(fp,'(a)',iostat=ierror) line
       if(ierror /= 0) return
       
       if( index(line,mystr(1:ilen)) /= 0) then
          ier = 0
          return
       endif
    enddo
  end subroutine locstr

  subroutine skipasterisk(fp,ier)
    implicit none
    integer, intent(in) :: fp
    integer, intent(out) :: ier
    character(256) :: line
    integer :: ierror

    ier=1
    do
       read(fp,'(a)',iostat=ierror) line
       if(ierror /= 0) return
          
       if(index(line,'*') /= 1) then
          ier=0
          backspace(fp)
          return
       endif
   enddo
  end subroutine skipasterisk

  subroutine mycdef(al,ji,je,k,cou,nte,cf,lexp,ltest,lsum)
    implicit none

    real*8, intent(in) :: al(*)
    real*8, intent(out) :: cou(0:9,*),cf(9,0:9)
    integer, intent(in) :: ji, je, k, nte
    logical, intent(in) :: lexp, ltest, lsum

    real*8 :: ctest
    integer :: icell, j, ii

    if(k <=0 ) then
       write(*,*) 'mycdef: k <=0'
       stop
    endif

    !c  data from a&m data files
    do  j=ji,je
       do  ii=1,9
          cf(ii,j)=creac(ii,j,k)
       enddo
    enddo

    if (ltest) then
       do  j=ji,je
          ctest=0.
          do  ii=1,9
             ctest=ctest+abs(cf(ii,j))
          enddo
          if (ctest.le.eps60) then
             write (*,*) 'error in subroutine cdef: zero fit coefficients'
             write (*,*) 'j,k = ',j,k,'  exit called!'
             stop
          endif
       enddo
    endif

    if (lsum) then
       do j=ji,je
          do icell=1,nte
             cou(j,icell)=cf(9,j)
          end do
       end do
       !c
       do  j=ji,je
          do  ii=8,1,-1
             do  icell=1,nte
                cou(j,icell)=cou(j,icell)*al(icell)+cf(ii,j)
             enddo
          enddo
       enddo
    endif
    !c
    if (lexp) then
       do  icell=1,nte
          cou(je,icell)=exp(max(dble(-100.),cou(je,icell)))
       enddo
    endif
    !c
    return
  end subroutine mycdef

  subroutine mycdefn(al,pl,k,cou,nte,cf,lexp,ltest,lsum)
    implicit none

    real*8, intent(in) :: al(*), pl(*)
    real*8, intent(out) :: cou(0:9,*), cf(9,0:9)
    integer, intent(in) :: k, nte
    logical, intent(in) :: lexp, ltest, lsum
    real*8 :: dummp(9)
    real*8 :: ccxm1, ccxm2, expo1, expo2, fpar1, fpar2, fpar3, s01,s02, ds12, ctest, extrap
    integer :: i, jj, kk, ifex, j, ii, icell
    
    if(k <=0 ) then
       write(*,*) 'mycdefn: k <=0'
       stop
    endif    

    !c  k>0:
    !c  data from array creac(9,0:9,k)
    do  j=1,9
       do  ii=1,9
          cf(ii,j)=creac(ii,j,k)
       enddo
    enddo

    if (ltest) then
       do  j=1,9
          ctest=0.
          do  ii=1,9
             ctest=ctest+abs(cf(ii,j))
          enddo
          if (ctest.le.eps60) then
             write (*,*)  'error in subroutine cdefn: zero fit coefficients'
             write (*,*) 'j,k = ',j,k,'  exit called!'
             stop
          endif
       enddo
    endif

    if (lsum) then
       do  icell=1,nte
          cou(1,icell)=0.
       enddo

       do  icell=1,nte
!          if (al(icell).lt.rcmn(k,2)) then             
!             !c  determine extrapolation coefficients for linear extrap. in ln(<s*v>)
!             s01=rcmn(k,2)
!             s02=log(2.)+rcmn(k,2)
!             ds12=s02-s01
!             expo1=0.
!             expo2=0.
!             do  j=1,9
!                jj=j-1
!                do  i=1,9
!                   ii=i-1
!                   expo1=expo1+s01**ii*pl(icell)**jj*cf(i,j)
!                   expo2=expo2+s02**ii*pl(icell)**jj*cf(i,j)
!                enddo
!             enddo
!             ccxm1=expo1
!             ccxm2=expo2
!             fpar1=ccxm1+(ccxm2-ccxm1)/ds12*(-s01)
!             fpar2=      (ccxm2-ccxm1)/ds12
!             fpar3=0.d0
!             !c
!             ifex=5
!             cou(1,icell)=extrap(al(icell),ifex,fpar1,fpar2,fpar3)
!             if (.not.lexp) cou(1,icell)=log(cou(1,icell))
!          else
             do  jj=9,1,-1
                dummp(jj)=cf(9,jj)
                do  kk=8,1,-1
                   dummp(jj)=dummp(jj)*al(icell)+cf(kk,jj)
                enddo
             enddo
             cou(1,icell)=dummp(9)
             do  jj=8,1,-1
                cou(1,icell)=cou(1,icell)*pl(icell)+dummp(jj)
             enddo
             if (lexp) cou(1,icell)=exp(max(dble(-100.),cou(1,icell)))
!          endif
         
       enddo
    endif
    return
  end subroutine mycdefn

  function myidez(i,j,izif)
    !c
    !c   i is an integer with  n dezimals, the function returns
    !c   myidez, the j.th dezimal
    !c   j must be .le. izif
    !c      if n .lt. izif and j .gt. n  , myidez=0
    !c      if n .gt. izif and j .eq. izif  , myidez=part of i on the left of
    !c      the dezimal j, including the j th dezimal
    !c      e.g.: i=1234,izif=3, j=1: myidez=4
    !c                           j=2: myidez=3
    !c                           j=3: myidez=12
    !c                           j=4: error ,j gt izif
    implicit none
    integer, intent(in) :: i, j, izif
    integer :: myidez, k, iz, iq, it, idif

    idif=izif-j+1
    if (idif.le.0) then
       write(*,*) 'error in function myidez            '
       stop       
    endif
    iz=i
    do  k=1,idif
       it=10**(izif-k)
       iq=iz/it
       iz=iz-iq*it
    enddo
    myidez=iq
    return
  end function myidez


  real*8 FUNCTION eirene_xs_CROSS(AL,K,IR,TEXT) result(res)
    !C
    !C  CROSS SECTION
    !C    AL=LN(ELAB), ELAB IN (EV)
    !C    RETURN CROSS SECTION IN CM**2
    !C
    !C  K>0 :  DATA FROM ARRAY CREAC(9,0:9,K)
    !C
    !C  K<0 :  DEFAULT MODEL
    !C
    !C  K=-1:  H + H+ --> H+ + H   CROSS SECTION, JANEV, 3.1.8
    !C         LINEAR EXTRAPOLATION AT LOW ENERGY END FOR LN(SIGMA)
    !C         IDENTICAL TO hydhel.tex, H.1, 3.1.8
    !C
    IMPLICIT NONE

    REAL*8, INTENT(IN) :: AL
    INTEGER, INTENT(IN) :: K, IR
    CHARACTER(LEN=*), INTENT(IN) :: TEXT
    REAL*8 :: B(8)
    REAL*8 :: RMIN, FP1, FP2, S01, S02, DS12, EXPO1, EXPO2, CROSS, CCXM1, CCXM2, EXPO, E, XI
    INTEGER :: IF8, II, I

    real*8,parameter :: CF1(9)=(/-3.274124E+01,-8.916457E-02,-3.016991E-02,9.205482E-03, 2.400267E-03,-1.927122E-03,3.654750E-04,-2.788866E-05, 7.422296E-07/)

    RMIN=-2.3025851E+00
    FP1=-3.2945896E+01
    FP2=-1.713112E-01
    !C
    IF (K.EQ.-1) THEN
       !C  ELAB BELOW MIMINUM ENERGY FOR FIT:
!       IF (AL.LT.RMIN) THEN
!          !C  USE ASYMPTOTIC EXPRESSION NO. IFMN=5
!          CROSS=myEXTRAP(AL,5,FP1,FP2,0._DP)
!          !C  ELAB ABOVE MAXIMUM ENERGY FOR FIT: NOT IN USE FOR K=-1
!          !C       ELSEIF (ELAB.GT.RMAX) THEN
!          !C  USE ASYMPTOTIC EXPRESSION NO. IFMX=5
!          !C         CROSS=myEXTRAP(AL,5,FP4,FP5,0.D0)
!          !C
!       ELSE
          EXPO=CF1(9)
          DO  II=1,8
             IF8=9-II
             EXPO=EXPO*AL+CF1(IF8)
          enddo
          CROSS=EXP(MAX(-100.d0,EXPO))
!       ENDIF
       !C
    ELSEIF (K.GT.0) THEN

       IF (IFTFLG(K,1) == 0) THEN

          !C  ELAB BELOW MINIMUM ENERGY FOR FIT:
!          IF (AL.LT.RCMN(K,1)) THEN
!             !C  USE ASYMPTOTIC EXPRESSION NO. IFEXMN(K)
!             IF (IFEXMN(K,1).LT.0) THEN
!                !C  DETERMINE EXTRAPOLATION COEFFICIENTS FOR LINEAR EXTRAP. IN LN(SIGMA)
!                S01=RCMN(K,1)
!                S02=LOG(2._DP)+RCMN(K,1)
!                DS12=S02-S01
!                EXPO1=CREAC(9,0,K)
!                EXPO2=CREAC(9,0,K)
!                DO  II=1,8
!                   IF8=9-II
!                   EXPO1=EXPO1*S01+CREAC(IF8,0,K)
!                   EXPO2=EXPO2*S02+CREAC(IF8,0,K)
!                enddo
!                CCXM1=EXPO1
!                CCXM2=EXPO2
!                FPARM(K,1,1)=CCXM1+(CCXM2-CCXM1)/DS12*(-S01)
!                FPARM(K,2,1)=      (CCXM2-CCXM1)/DS12
!                FPARM(K,3,1)=0.D0
!                !C
!                IFEXMN(K,1)=5
!             ENDIF
!             CROSS=myEXTRAP(AL,IFEXMN(K,1),FPARM(K,1,1),FPARM(K,2,1),FPARM(K,3,1))
!             !C  ELAB ABOVE MAXIMUM ENERGY FOR FIT:
!          ELSEIF (AL.GT.RCMX(K,1)) THEN
!             !C  USE ASYMPTOTIC EXPRESSION NO. IFEXMX(K,1)
!             CROSS=myEXTRAP(AL,IFEXMX(K,1),FPARM(K,4,1),FPARM(K,5,1),FPARM(K,6,1))
!          ELSE
             EXPO=CREAC(9,0,K)
             DO  II=1,8
                IF8=9-II
                EXPO=EXPO*AL+CREAC(IF8,0,K)
             enddo
             CROSS=EXP(MAX(-100.d0,EXPO))
!          ENDIF

       ELSE IF (IFTFLG(K,1) == 3) THEN
          !C  default extrapolation ifexmn=-1 not yet available
          !C  ELAB BELOW MINIMUM ENERGY FOR FIT:
!          IF (AL.LT.RCMN(K,1)) THEN
!             !C  USE ASYMPTOTIC EXPRESSION NO. IFEXMN(K)
!             CROSS=myEXTRAP(AL,IFEXMN(K,1),FPARM(K,1,1),FPARM(K,2,1), FPARM(K,3,1))
!             !C  ELAB ABOVE MAXIMUM ENERGY FOR FIT:
!          ELSEIF (AL.GT.RCMX(K,1)) THEN
!             !C  USE ASYMPTOTIC EXPRESSION NO. IFEXMX(K,1)
!             CROSS=myEXTRAP(AL,IFEXMX(K,1),FPARM(K,4,1),FPARM(K,5,1), FPARM(K,6,1))
!          ELSE
             E = EXP(AL)
             XI = CREAC(1,0,K)
             B(1:8) = CREAC(2:9,0,K)
             CROSS = B(1)*LOG(E/XI)
             DO I=1,7
                CROSS = CROSS + B(I+1)*(1.D0-XI/E)**I
             END DO
             CROSS = CROSS * 1.D-13/(XI*E)
!          ENDIF
       ELSE
          write(*,*) 'eirene_xs_cross:'
          WRITE (*,*) ' WRONG FITTING FLAG IN CROSS '
          WRITE (*,*) ' K = ',K,' IFTFLG = ',IFTFLG(K,1)
          WRITE (*,*) 'REACTION NO. ',IR
          stop
       END IF
    ELSE
       write(*,*) 'eirene_xs_cross:'
       WRITE (*,*) 'ERROR IN CROSS: K= ',K
       WRITE (*,*) 'CALLED FROM ',TEXT
       WRITE (*,*) 'REACTION NO. ',IR
       stop
    ENDIF

    res=cross
    RETURN
  END FUNCTION eirene_xs_CROSS
  
  real*8 FUNCTION myEXTRAP(ELAB,IFLAG,FP1,FP2,FP3) result(res)
    !C
    !C  NOTE:
    !C  INPUT:  ELAB IS LOG OF RELATIVE ENERGY, OR LOG OF TEMP
    !C  OUTPUT: EXTRAP IS NOT LOG, BUT THE TRUE VALUE
    !C
    !C  FUNCTION FOR EXTRAPOLATING SINGLE PARAMETER FITS BEYOND THEIR
    !C  RANGE OF VALIDITY
    !C  TYPE  IFLAG=1--4: JANEV ET AL. , SPRINGER, 1987, P13
    !C  TYPE  IFLAG=5  BACHMANN ET AL., IPP-REPORT, .....ELASTIC
    !C
    IMPLICIT NONE

    REAL*8, INTENT(IN) :: ELAB, FP1, FP2, FP3
    INTEGER, INTENT(IN) :: IFLAG
    REAL*8 :: X, EL, EXTRAP

    IF (IFLAG.EQ.1) THEN
       !C  NON ZERO THRESHOLD
       EXTRAP=0.
    ELSEIF (IFLAG.EQ.2) THEN
       !C  EXTRAPOLATION AT HIGH ENERGY END FOR REACTIONS WITH NON ZERO THRESHOLD
       !C  FP1 SHOULD BE = E_THRESHOLD (EV)
       EL=EXP(ELAB)
       X=EL/FP1
       EXTRAP=FP2*X**FP3*LOG(X)
    ELSEIF (IFLAG.EQ.3) THEN
       EXTRAP=EXP(FP1+FP2*ELAB)
    ELSEIF (IFLAG.EQ.4) THEN
       !C
       !C  OUT
       !C       EXTRAP=EXP((FP1+FP2*ELAB)**2)
       !C
    ELSEIF (IFLAG.EQ.5) THEN
       !C  LINEAR OR QUADRATIC EXTRAPOLATION IN LN(SIGMA)
       EXTRAP=EXP(FP1+FP2*ELAB+FP3*ELAB**2)
    ELSE
       write(*,*) 'eirene_xs_extrap:'
       WRITE (*,*) 'ERROR IN EXTRAP. EXIT CALLED '
       stop
    ENDIF
    res=extrap
    RETURN
  END FUNCTION myEXTRAP

end module mod_eirene_xs


