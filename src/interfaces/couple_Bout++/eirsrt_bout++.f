C
C
C  21 Oct. 2013 Samad Mekkaoui

C  BOUT++ <-> Eierene cycling subroutine (LAPD)

      SUBROUTINE EIRSRT(DT,nx,ny,nz,niei,teei,viz,src_transp,
     . msrc_transp,esrc_transp,eisrc_transp,at_den_transp,
     . mol_den_transp,
     . ndim,prtr,nstrab)
      USE ALLOC_BOUT
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CVARUSR
      USE EIRMOD_CSPEI 
      USE EIRMOD_CPES 
      USE eirmod_braspoi
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, only: iunout
      USE EIRMOD_COMSOU
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      USE EIRMOD_CSPEZ
      USE EIRMOD_CPLOT
      USE EIRMOD_CESTIM
      USE EIRMOD_COUTAU
      USE eirmod_comnnl
      
c     USE EIRMOD_COMPRT, ONLY: IUNOUT
      implicit none
      REAL(DP),INTENT(IN) :: DT
      
      REAL(DP) :: DTIMVO
      real(dp), parameter :: ech=1.602e-19
      INTEGER, INTENT(IN) :: ndim
      INTEGER, INTENT(INOUT) :: prtr(ndim),nstrab
      integer ip,it,ir,i,j,k,ncell,ipls,nx,ny,nz,istri
      double precision :: niei(nz,ny,nx),viz(nz,ny,nx)
      double precision :: teei(nz,ny,nx)
      double precision,intent(out)::src_transp(nz,ny,nx)
      double precision,intent(out)::msrc_transp(nz,ny,nx)
      double precision,intent(out)::esrc_transp(nz,ny,nx)
      double precision,intent(out)::eisrc_transp(nz,ny,nx)
      double precision,intent(out)::at_den_transp(nz,ny,nx)
      double precision,intent(out)::mol_den_transp(nz,ny,nx)
      logical,save :: ifirst
      data ifirst /.true./

      
      if(.not.allocated(ni_bout))allocate(ni_bout(nx,ny,nz))
       if(.not.allocated(te_bout))allocate(te_bout(nx,ny,nz))
        if(.not.allocated(viz_bout))allocate(viz_bout(nx,ny,nz))

      print*, nt3rd, np2nd,nr1st
      print*,'***************',nx,ny,nz
      nx_bout=nx
      ny_bout=ny
      nz_bout=nz
      do i=1,nx
        do j=1,ny
          do k=1,nz
         ni_bout(i,j,k)=niei(k,j,i)
         te_bout(i,j,k)=teei(k,j,i)       
         viz_bout(i,j,k)=viz(k,j,i)
          enddo
        enddo
      enddo
     

      print *, '==========w108========', dt, '==========' 
      if(ifirst) then       
      call eirene_eirene(DT,.FALSE.,.FALSE.,1,.FALSE.)
     
      
      if (ndim < nstra) then
        write (iunout,*) ' ERROR in COUPLING TO EIRENE '
        write (iunout,*) ' ARRAY PRTR FOR TRANSFERRING THE INDICES OF',
     .                   ' PROCESSORS IS TOO SMALL ! '
        CALL EIRENE_EXIT_OWN(1)
      END IF

      prtr(1:nstra) = npesta(1:nstra)
      write (iunout,*) prtr(2) ,npesta(2)
      write (iunout,*) prtr(1) ,npesta(1)  
      nstrab=nstra
      else
      IITER=1
      ITIMV=2 
      DTIMVN=DT
c     DTIMVO=DT
      DTIMV=DT
      
      call EIRENE_EIRENE_COUPLE (DT,.FALSE.,1,.false.)
      nstrab=nstrai 
       prtr(1:nstrai) = npesta(1:nstrai)
       print *, '==========w116========', dt, '=========='
       print*, 'eirenefinished02'
      endif 
      ifirst=.false.
      ntmstp = 1 
c      dtimvo=DT
c      NFILEJ = 3
       DTIMVN=DT
       DTIMV=DT
         if (my_pe==0) then
	do ip=1, ny
           do it = 1, nz
              do ir = 1, nx
          src_transp(it,ip,ir) = sum(src_bck_bout(ir,ip,it,:))/ech
           msrc_transp(it,ip,ir) = sum(msrc_bck_bout(ir,ip,it,:))/ech
            esrc_transp(it,ip,ir) = sum(esrc_bck_bout(ir,ip,it,:))
             eisrc_transp(it,ip,ir) = sum(eisrc_bck_bout(ir,ip,it,:))
             at_den_transp(it,ip,ir)=sum(at_den(ir,ip,it,:))
               mol_den_transp(it,ip,ir)=sum(mol_den(ir,ip,it,:))
	enddo
	   enddo
	      enddo
        

	 endIF       
       write (iunout,*) prtr(2) ,npesta(2) 
       write (iunout,*) prtr(1) ,npesta(1)    
       END subroutine eirsrt
