      subroutine EIRENE_collect_census
cdr:  Aug. 2015 comments added
c
c this routine is called for each processer my_pe
c it first defines the census array rpartw(i) and total flux peflux, for each processor.
c it then tries to combine these onto a single new census.
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMNNL
      USE EIRMOD_COUTAU
      USE EIRMOD_COMSOU
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CPES
 
      IMPLICIT NONE
 
      INCLUDE 'mpif.h'
      real(dp), allocatable :: rpselect(:), rand(:), rdistrib(:),
     .                         rscat(:), rbuf(:,:)
      real(dp) :: ra, weight, peflux, totflux, sumrpw, sclfac, add,
     .            totrpw
      real(dp), external :: ranf_eirene
      integer, allocatable :: iranpro(:), ibuf(:,:)
      integer :: ier, i, istr, ncoreal, itotal, il, im, iu, ipe,
     .           ityp, iphot, iatm, imol, iion
      integer :: icopro(0:nprs), idistrib(0:nprs), icosend(0:nprs)
 
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (FLXFAC,NSTRAI+1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
 
      RPARTW(0)=0.0
 
!  each processor prepares his census for transfer to processor 0
 
!  processor my_pe has accumulated iprnli scores on census

      PEFLUX=0._DP
      DO I=1,IPRNLI
        ISTR=IPART(8,I)
        ITYP=ISPEZI(IPART(9,I),-1)
        WEIGHT=RPART(9,I)
        IF (ITYP.EQ.0) THEN
          IPHOT=ISPEZI(IPART(I,9),0)
          ADD=WEIGHT*FLXFAC(ISTR)*NPRT(IPHOT)
        ELSEIF (ITYP.EQ.1) THEN
          IATM=ISPEZI(IPART(9,I),1)
          ADD=WEIGHT*FLXFAC(ISTR)*NPRT(NSPH+IATM)
        ELSEIF (ITYP.EQ.2) THEN
          IMOL=ISPEZI(IPART(9,I),2)
          ADD=WEIGHT*FLXFAC(ISTR)*NPRT(NSPA+IMOL)
        ELSEIF (ITYP.EQ.3) THEN
          IION=ISPEZI(IPART(9,I),3)
          ADD=WEIGHT*FLXFAC(ISTR)*NPRT(NSPAM+IION)
        ENDIF
        RPARTW(I)=RPARTW(I-1)+WEIGHT*FLXFAC(ISTR)
        PEFLUX = PEFLUX + ADD
      END DO

c  peflux is the total, fully scaled census "atomic" flux accumulated on my_pe
c  rpartw(i) is the cummulative, scaled, flux distribution on census accumulated on my_pe 
      write (iunout,*) ' collect census, from my_pe            ',my_pe
      write (iunout,*) ' scores on census from my_pe: iprnli   ',iprnli
      write (iunout,*) ' atomic flux on census from my_pe (Amp)',peflux
 

! transfer maximum possible rpartw to processor 0
 
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
 
      call mpi_allreduce(iprnli,itotal,1,MPI_INTEGER,
     .                   MPI_SUM,MPI_COMM_WORLD,ier)
 
      call mpi_allreduce(peflux,totflux,1,MPI_REAL8,
     .                   MPI_SUM,MPI_COMM_WORLD,ier)
c
c  cummulated scores, and atomic flux 
      write (iunout,*) ' tentative: itotal, totflux', itotal, totflux

      if (itotal <= nprnl) then
! THERE IS ENOUGH STORAGE for all scores from all processors.
!                          send all particles to processor 0
 
 
        allocate (rbuf(size(rpart,1),size(rpart,2)))
        allocate (ibuf(size(ipart,1),size(ipart,2)))       
        rbuf = 0._dp
        ibuf = 0
         
        call mpi_gather(iprnli,1,MPI_INTEGER,
     .                  icopro,1,MPI_INTEGER,0,
     .                  MPI_COMM_WORLD,ier)
 

        if (my_pe == 0) then
          icosend = icopro*npartt
          idistrib(0) = 0
          do ipe = 1, nprs
            idistrib(ipe) = idistrib(ipe-1) + icosend(ipe-1)
          end do
        end if
        call mpi_gatherv(rpart,iprnli*npartt,MPI_REAL8,
     .                   rbuf,icosend,idistrib,MPI_REAL8,
     .                   0,MPI_COMM_WORLD,ier)
        if (my_pe == 0) then
           rpart = rbuf
        end if        
 
        if (my_pe == 0) then
          icosend = icopro*mpartt
          idistrib(0) = 0
          do ipe = 1, nprs
            idistrib(ipe) = idistrib(ipe-1) + icosend(ipe-1)
          end do
        end if
        call mpi_gatherv(ipart,iprnli*mpartt,MPI_INTEGER,
     .                   ibuf,icosend,idistrib,MPI_INTEGER,
     .                   0,MPI_COMM_WORLD,ier)
        if (my_pe == 0) then
           ipart = ibuf
        end if        
  
        iprnli = itotal
        write (iunout,*) 'total no. of scores on census ', itotal
        write (iunout,*) 'total atomic flux on census   ', totflux
 
        deallocate (rbuf)
        deallocate (ibuf)

! THERE IS NOT ENOUGH STORAGE for all scores from all processors.

      else  ! here: itotal > nprnl:  carry out some condensation: 
!                                    sample exactly nprnl scores from the full set of itotal scores
 
        itotal = nprnl
 
        allocate (rpselect(-1:nprs))
        if (my_pe == 0) then
          rpselect(-1) = 0._dp
          rpselect(0) = RPARTW(iprnli)  !cdr  start with my_pe=0 
        end if
 
        CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
        call mpi_gather(rpartw(iprnli),1,MPI_REAL8,
     .                    rpselect(0:),1,MPI_REAL8,0,
     .                    MPI_COMM_WORLD,ier)
 
!pb        write (iunout,*) ' rpselect before summation '
!pb        write (iunout,'(i6,es12.4)') (ipe,rpselect(ipe),ipe=-1,nprs)
 
        if (my_pe == 0) then
          do ipe=1, nprs-1
            rpselect(ipe) = rpselect(ipe-1) + rpselect(ipe) !cdr rpselect(0) war schon gesetzt.
          end do
 
!pb          write (iunout,*) 'total cummulated flux on census ', totflux
!pb          write (iunout,*) 'rpselect '
!pb          write (iunout,'(i6,es12.4)') (ipe,rpselect(ipe),ipe=-1,nprs-1)
 
! now find random numbers for NPRNL particles
 
          allocate (rand(nprnl))
          allocate (iranpro(nprnl))
          icopro = 0
 
!  first step : find processor iu, from which to sample a census score
!  second step:  after that sample from that processor iu
          do i = 1, nprnl
            ra = ranf_eirene() * rpselect(nprs-1)
 
            IL=0
            IU=NPRS
 
            if ( ra <= rpselect(0) ) then
              iu = 0
            else
c  binary search amongst processors
              DO WHILE (IU-IL.gt.1)
                IM=(IU+IL)*0.5
                IF (RA.GE.rpselect(IM)) THEN
                  IL=IM
                ELSE
                  IU=IM
                ENDIF
              END DO
            end if
 
            icopro(iu) = icopro(iu) + 1
            iranpro(i) = iu
            rand(i) = ra - rpselect(iu-1)
c  processor iu identified, random number rand(i) set for sampling from census
c                           restricted to this processor 
          end do
 
! setup displacements for distribution of random numbers
          idistrib(0) = 0
          do ipe = 1, nprs-1
            idistrib(ipe) = idistrib(ipe-1) + icopro(ipe-1)
          end do
 
! sort random numbers
          allocate (rdistrib(nprnl))
          do i = 1, nprnl
            ipe = iranpro(i)
            idistrib(ipe) = idistrib(ipe) + 1
            rdistrib(idistrib(ipe)) = rand(i)
          end do
 
! reset displacements for distribution of random numbers
          idistrib(0) = 0
          do ipe = 1, nprs
            idistrib(ipe) = idistrib(ipe-1) + icopro(ipe-1)
          end do
 
        write (iunout,*) 'number of particles per processor '
        write (iunout,'(10i6)') (icopro(ipe),ipe=0,nprs)
 
        end if
 
! broadcast numbers of required particles per processor

       	if (.not.allocated(rdistrib)) allocate(rdistrib(nprnl))
 
        CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
        call mpi_scatter(icopro ,1,MPI_INTEGER,
     .                   ncoreal,1,MPI_INTEGER,0,
     .                   MPI_COMM_WORLD,ier)
 
! broadcast random numbers (position) for each processor
        allocate(rscat(ncoreal))
        call mpi_scatterv(rdistrib,icopro,idistrib,MPI_REAL8,
     .                       rscat,ncoreal,        MPI_REAL8,0,
     .                       MPI_COMM_WORLD,ier)
 
 
! on each processor look for the indices of the particles to be
! put into the global census arrays rpartc,ipartc. 
! Sampling with replacement, fill rpartc,ipartc per processor, 
! then gather theses into one single array on my_pe=0
 

        sumrpw = 0._dp
        do i = 1, ncoreal
 
          RA = RSCAT(I)
 
          IL=0
          IU=IPRNLI
 
c  binary search
          DO WHILE (IU-IL.gt.1)
            IM=(IU+IL)*0.5
            IF (RA.GE.RPARTW(IM)) THEN
              IL=IM
            ELSE
              IU=IM
            ENDIF
          end do
 
          rpartc(:,i) = rpart(:,iu)
          ipartc(:,i) = ipart(:,iu)
 
          ISTR=IPARTC(8,I)
          ITYP=ISPEZI(IPARTC(9,I),-1)
          WEIGHT=RPARTC(9,I)
          IF (ITYP.EQ.0) THEN
             IPHOT=ISPEZI(IPARTC(I,9),0)
             ADD=WEIGHT*FLXFAC(ISTR)*NPRT(IPHOT)
          ELSEIF (ITYP.EQ.1) THEN
             IATM=ISPEZI(IPARTC(9,I),1)
             ADD=WEIGHT*FLXFAC(ISTR)*NPRT(NSPH+IATM)
          ELSEIF (ITYP.EQ.2) THEN
             IMOL=ISPEZI(IPARTC(9,I),2)
             ADD=WEIGHT*FLXFAC(ISTR)*NPRT(NSPA+IMOL)
          ELSEIF (ITYP.EQ.3) THEN
             IION=ISPEZI(IPARTC(9,I),3)
             ADD=WEIGHT*FLXFAC(ISTR)*NPRT(NSPAM+IION)
          ENDIF
c   accumulated atomic flux from current processor
          sumrpw = sumrpw + add

 
!          write (iunout,*) i, ra, iu
!          write (iunout,'(15i6)') ipartc(:,i)
!          write (iunout,'(6es12.4,/,6es12.4)') rpartc(:,i)
 
 
        end do
 
! send chosen particles to processor 0
 
        if (my_pe == 0) then
          icosend = icopro*npartt
          idistrib(0) = 0
          do ipe = 1, nprs
            idistrib(ipe) = idistrib(ipe-1) + icosend(ipe-1)
          end do
        end if
        call mpi_gatherv(rpartc,ncoreal*npartt,MPI_REAL8,
     .                   rpart,icosend,idistrib,MPI_REAL8,
     .                   0,MPI_COMM_WORLD,ier)
 
        if (my_pe == 0) then
          icosend = icopro*mpartt
          idistrib(0) = 0
          do ipe = 1, nprs
            idistrib(ipe) = idistrib(ipe-1) + icosend(ipe-1)
          end do
        end if
        call mpi_gatherv(ipartc,ncoreal*mpartt,MPI_INTEGER,
     .                   ipart,icosend,idistrib,MPI_INTEGER,
     .                   0,MPI_COMM_WORLD,ier)
 
        call mpi_reduce (sumrpw,totrpw,1,MPI_REAL8,MPI_SUM,
     .                   0,MPI_COMM_WORLD,ier)
 
!pb        if (my_pe == 0) then
!pb           write (iunout,*) ' rpart collected from all '
!pb           do i=1, nprnl
!pb             write (iunout,'(i6,4es12.4)') i,rpart(1:3,i),weight
!              write (iunout,*) ' i, my_pe ',i,my_pe
!              write (iunout,'(15i6)') ipart(:,i)
!              write (iunout,'(6es12.4,/,6es12.4)') rpart(:,i)
!pb           end do
!pb        end if
 
        iprnli = itotal
 
! rescale
 
        if (my_pe == 0) then
          sclfac = totflux / totrpw
          write (iunout,*) ' totrpw ',totrpw
          write (iunout,*) ' sclfac ',sclfac
          do i=1,iprnli
            rpart(9,i) = rpart(9,i) * sclfac
          end do
        end if
 
        if (allocated(rscat)) deallocate(rscat)
        if (allocated(rdistrib)) deallocate(rdistrib)
        if (allocated(rand)) deallocate (rand)
        if (allocated(iranpro)) deallocate (iranpro)
        if (allocated(rpselect)) deallocate (rpselect)
 
      end if
 
      return
 
      end subroutine EIRENE_collect_census
