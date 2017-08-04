      subroutine EIRENE_collect_census

cpb July  17: bug fix, rpselect allocation from -1 , not from 0
cdr sept. 15: bug fix:   after re-sampling (with replacement) from census, the weight of
cdr                      sampled census particles is set to 1.0, rather than keeping the old weight.
cdr                      The census flux is regarded as "discrete distribution" for the index "i" of a particle,
cdr                      and the weight stored on census during particle tracing is the probability mass of index "i"

cdr:  Aug. 2015 comments added
c
c this routine is called for each processer my_pe
c it first defines the census array rpartw(i) and total flux peflux, for each processor.
c it then tries to combine these onto a single new census.
c If the combined census from all processors contains too many particles, then the
c reduction is done by re-sampling, just like in subr. locate for re-launch from census.

cdr  rpselect reduziert auf 0,nprs-1
cdr  addph,adda,addm,addi: type resolved census fluxes.
c

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
      real(dp) :: ra, weight, peflux, 
     .            totflux, sumrpw, sclfac, add, totrpw,
     .            addph, adda, addm, addi
C      real(dp) :: pefluxp(0:nprs-1), sumrpwp(0,nprs-1)
      real(dp), external :: ranf_eirene
      integer, allocatable :: iranpro(:), ibuf(:,:)
      integer :: ier, i, istr, ncoreal, itotal, il, im, iu, ipe,
     .           ityp, iphot, iatm, imol, iion
      integer :: icopro(0:nprs), idistrib(0:nprs), icosend(0:nprs)

      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (FLXFAC,NSTRAI+1,MPI_REAL8,0,MPI_COMM_WORLD,ier)

      RPARTW(0)=0.0

!  each processor prepares his census for later transfer to processor 0

!  processor my_pe has accumulated iprnli scores on census

!  distinct from tmstep (scoring on census)

      PEFLUX =0._DP
      ADDPH =0._DP
      ADDA  =0._DP
      ADDM  =0._DP
      ADDI  =0._DP

      DO I=1,IPRNLI
        ISTR=IPART(8,I)
        ITYP=ISPEZI(IPART(9,I),-1)
        WEIGHT=RPART(9,I)    !  THIS WEIGHT SHOULD ALREADY CONTAIN THE PARTICLE BALANCE 
!                               RESCALING FACTORS, DONE LATER IN TMSTEP.
        IF (ITYP.EQ.0) THEN
          IPHOT=ISPEZI(IPART(I,9),0)
          ADD=WEIGHT*FLXFAC(ISTR)*NPRT(IPHOT)
          ADDPH=ADDPH+ADD
        ELSEIF (ITYP.EQ.1) THEN
          IATM=ISPEZI(IPART(9,I),1)
          ADD=WEIGHT*FLXFAC(ISTR)*NPRT(NSPH+IATM)
          ADDA=ADDA+ADD
        ELSEIF (ITYP.EQ.2) THEN
          IMOL=ISPEZI(IPART(9,I),2)
          ADD=WEIGHT*FLXFAC(ISTR)*NPRT(NSPA+IMOL)
          ADDM=ADDM+ADD
        ELSEIF (ITYP.EQ.3) THEN
          IION=ISPEZI(IPART(9,I),3)
          ADD=WEIGHT*FLXFAC(ISTR)*NPRT(NSPAM+IION)
          ADDI=ADDI+ADD
        ENDIF
! cummulativ distribution of weight of particle no I, for sampling. not atomic flux
        RPARTW(I)=RPARTW(I-1)+WEIGHT*FLXFAC(ISTR)
! total flux von census, atomic flux (AMP)
        PEFLUX   = PEFLUX + ADD
      END DO

c  peflux is the total, fully scaled census "atomic" flux accumulated on my_pe
c  rpartw(i) is the cummulative, scaled, flux distribution on census accumulated on my_pe
      call eirene_leer(1)
      write (iunout,*) 'COLLECT CENSUS, from my_pe      ',my_pe
      write (iunout,*) 'scores on census: iprnli    ',iprnli
      write (iunout,*) 'atomic flux on census (Amp) ',peflux
      call eirene_masr4('at. flx.: addph, adda, addm, addi',
     .                   addph, adda, addm, addi)

      
c     pefluxp(my_pe)=peflux
      

! transfer maximum possible rpartw to processor 0

      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

      call mpi_allreduce(iprnli,itotal,1,MPI_INTEGER,
     .                   MPI_SUM,MPI_COMM_WORLD,ier)

      call mpi_allreduce(peflux,totflux,1,MPI_REAL8,
     .                   MPI_SUM,MPI_COMM_WORLD,ier)
c
c  cummulated number of scores, and atomic flux, summed from all PEs.
      if (my_pe.eq.0) THEN
        write (iunout,*) ' itotal, totflux', itotal, totflux
      ENDIF

      if (itotal <= nprnl) then
! THERE IS ENOUGH STORAGE for all scores from all processors.
!                          send all particles to processor 0


        allocate (rbuf(size(rpart,1),size(rpart,2)))
        allocate (ibuf(size(ipart,1),size(ipart,2)))
        rbuf = 0._dp
        ibuf = 0

! store the numbers of stored particles per processor in array icopro
        call mpi_gather(iprnli,1,MPI_INTEGER,
     .                  icopro,1,MPI_INTEGER,0,
     .                  MPI_COMM_WORLD,ier)


        if (my_pe == 0) then
! icosend is the number of real values gathered from the individual processors
          icosend = icopro*npartt
! idistrib gives the starting points for each processor on receiving buffer rbuf
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
! icosend is the number of integer values gathered from the individual processors
          icosend = icopro*mpartt
! idistrib gives the starting points for each processor on receiving buffer rbuf
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

        if (my_pe == 0) then
          write (iunout,*) 
     .           ' There is enough storage for collected census'
          write (iunout,*) ' No bootstrapping done in "collect_census"'
          write (iunout,*) ' total no. of scores on census ', itotal
          write (iunout,*) ' total atomic flux on census   ', totflux
        endif

        deallocate (rbuf)
        deallocate (ibuf)

! THERE IS NOT ENOUGH STORAGE for all scores from all processors.

      else  ! here: itotal > nprnl:  carry out some condensation:
!                                    sample exactly nprnl scores from the full set of itotal scores
        if (my_pe == 0) then
          write (iunout,*) 
     .           ' There is not enough storage for collected census '
          write (iunout,*) ' bootstrapping done in "collect_census"'
          write (iunout,*) ' itotal, nprnl = ',itotal,nprnl
        endif

        itotal = nprnl

        allocate (rpselect(-1:nprs))    !  reicht wohl:  0,nprs-1
!pb     allocate (rpselect(0:nprs-1))   !  reicht wohl:  0,nprs-1 !pb nein!

        if (my_pe == 0) then
          rpselect(-1) = 0._dp
!cdr  start with my_pe=0, flux to census (in amp, not atomic flux!)
          rpselect(0) = RPARTW(iprnli)  
        end if

        CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
! fetch the total flux from all the individual processors
! damit wird obiges rpselect(0) nochmal ueberschrieben
        call mpi_gather(rpartw(iprnli),1,MPI_REAL8,
     .                  rpselect(0:nprs-1),1,MPI_REAL8,0,   
     .                  MPI_COMM_WORLD,ier)

        if (my_pe == 0) then
! build accumulated flux distribution for all processores on my_pe=0
cdr rpselect(0) war schon gesetzt.
          do ipe=1, nprs-1
            rpselect(ipe) = rpselect(ipe-1) + rpselect(ipe)
          end do

          write (iunout,*) 'total cummulated flux on census (AMP) ' 
          write (iunout,*) 'rpselect '
          write (iunout,'(i6,es12.4)') (ipe,rpselect(ipe),ipe=0,nprs-1)

! now find random numbers to select NPRNL particles from the total of itotal scores

          allocate (rand(nprnl))
          allocate (iranpro(nprnl))
          icopro = 0

!  first step : for each of the nprnl new census scores,
!               find processor iu, from which to sample a census score
!  second step:  after that sample from that processor iu

!  now: first step:
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
                if (im.gt.nprs-1) then 
                  write (iunout,*) 'error in do while, collect census'
                  call eirene_exit_own(1)
                endif
                IF (RA.GE.rpselect(IM)) THEN
                  IL=IM
                ELSE
                  IU=IM
                ENDIF
              END DO
            end if


c  icopro(iu)       entries to be sampled from sub-census from processor iu
c  iranpro(i) = iu: random number i samples from sub-census from processor iu
c  rand(i) is the reduced random number, for sampling within sub-census iu only
            icopro(iu) = icopro(iu) + 1
            iranpro(i) = iu
            rand(i) = ra - rpselect(iu-1)
c
c  for each random number i the processor iu is identified (iu=iranpro(i)), 
c  random number rand(i) set for sampling from census
c                           restricted to this processor iu
          end do


          write (iunout,*) 'number of particles to be resampled ',
     .                     'per processor '
          do ipe=0,nprs-1
            write (iunout,'(2i10)') ipe, icopro(ipe)
          enddo
          write (iunout,'(A7,i10)') 'total  ', sum(icopro(0:nprs-1))


! setup displacements for distribution of random numbers, for each processor IPE
! idistrib(ipe) contains the number of random samples summed up until processor ipe-1, ipe=1,nprs-1 
! idistrib(ipe)+1 is the initial storage for resampled particles from processor ipe
          idistrib(0) = 0
          do ipe = 1, nprs-1
            idistrib(ipe) = idistrib(ipe-1) + icopro(ipe-1)
          end do
          do ipe=0,nprs-1
            write (iunout,*) 'ipe,idistrib(ipe),I ',ipe,idistrib(ipe)
          enddo

! assign random numbers to their corresponding processors, 
! store them in rdistrib, one such array for each processor ipe
          allocate (rdistrib(nprnl))
          do i = 1, nprnl
            ipe = iranpro(i)   !  this ipe is running from 0 to nprs-1
            idistrib(ipe) = idistrib(ipe) + 1
            rdistrib(idistrib(ipe)) = rand(i)
          end do

! reset displacements for distribution of random numbers
! idistrib is now cumulative number of random numbers per processor, ipe= 1,nprs
          idistrib(0) = 0
          do ipe = 1, nprs
            idistrib(ipe) = idistrib(ipe-1) + icopro(ipe-1)
          end do
          do ipe=1,nprs
            write (iunout,*) 'ipe,idistrib(ipe),II ',ipe,idistrib(ipe)
          enddo
      
c   now we know: for each random number i: i=1,nprnl
c         iu=iranpro(i)  : use this on processor iu
c         
c                   

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
        addph  = 0._dp
        adda   = 0._dp
        addm   = 0._dp
        addi   = 0._dp

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


c  partc, ipartc will later be used in tmstep to store census for 
c  re-sampling in locate at next time-step
c  here we abuse this storage to for the re-sampled census per stratum.
          rpartc(:,i) = rpart(:,iu)
          ipartc(:,i) = ipart(:,iu)

          ISTR=IPARTC(8,I)
          ITYP=ISPEZI(IPARTC(9,I),-1)
cdr> Sept. 2015
cdr  reset weight to one, because sampling according to weight is already accounting for rpartc.
c                       (same as in locate, except in case of one-by-one relaunch: then keep weight)
          RPARTC(9,I)=1.0
c
          WEIGHT=RPARTC(9,I)
          IF (ITYP.EQ.0) THEN
             IPHOT=ISPEZI(IPARTC(I,9),0)
             ADD=WEIGHT*FLXFAC(ISTR)*NPRT(IPHOT)
             addph=addph+add
          ELSEIF (ITYP.EQ.1) THEN
             IATM=ISPEZI(IPARTC(9,I),1)
             ADD=WEIGHT*FLXFAC(ISTR)*NPRT(NSPH+IATM)
             adda=adda+add
          ELSEIF (ITYP.EQ.2) THEN
             IMOL=ISPEZI(IPARTC(9,I),2)
             ADD=WEIGHT*FLXFAC(ISTR)*NPRT(NSPA+IMOL)
             addm=addm+add
          ELSEIF (ITYP.EQ.3) THEN
             IION=ISPEZI(IPARTC(9,I),3)
             ADD=WEIGHT*FLXFAC(ISTR)*NPRT(NSPAM+IION)
             addi=addi+add
          ENDIF

c   accumulated atomic flux from current processor
          sumrpw = sumrpw + add
        end do

cdr diagnose resampling procedure:
        write (iunout,*) ' resampled atomic flux from my_pe '
        write (iunout,*) ' total, sumrpw  ', sumrpw
        call eirene_masr4('at. flx.: addph, adda, addm, addi',
     .                     addph, adda, addm, addi)
        
c       sumrpwp(ipe)=sumrpw

!          write (iunout,*) i, ra, iu
!          write (iunout,'(15i6)') ipartc(:,i)
!          write (iunout,'(6es12.4,/,6es12.4)') rpartc(:,i)


    

! send re-sampled particles to processor 0
c  combine all the resampled census from all processors into a single one: rpart, ipart

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

! rescale to preserve census flux despite of resampling

        if (my_pe == 0) then
          sclfac = totflux / totrpw
          write (iunout,*) ' totrpw ',totrpw
          write (iunout,*) ' sclfac ',sclfac
          do i=1,iprnli
            rpart(9,i) = rpart(9,i) * sclfac
          end do
        end if
cdr  for resampling in locate at next timestep:
cdr  rpartc=rpart is done in tmstep.

        write (iunout,*) 'collect_census finished for MY_PE = ',my_pe
        call eirene_leer(2)

        if (allocated(rscat)) deallocate(rscat)
        if (allocated(rdistrib)) deallocate(rdistrib)
        if (allocated(rand)) deallocate (rand)
        if (allocated(iranpro)) deallocate (iranpro)
        if (allocated(rpselect)) deallocate (rpselect)

      end if

      return

      end subroutine EIRENE_collect_census
