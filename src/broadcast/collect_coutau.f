cdr  Nov.17: comments started...
!pb  060309  mpi_real8 --> mpi_double_precision


      subroutine eirene_collect_coutau
cdr this routine......
cdr .....
c
c    npesta(istra):  master processor ("group-leader") for each stratum ISTRA
c
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CPES
      USE EIRMOD_COUTAU
      USE EIRMOD_CSPEZ
      use EIRMOD_CESTIM
      use EIRMOD_COMUSR
      use EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CGRID
      USE EIRMOD_COMSOU
      USE EIRMOD_CSDVI
      USE EIRMOD_CSPEI
      USE EIRMOD_MPI
      IMPLICIT NONE

!      include 'mpif.h'
      REAL(DP), ALLOCATABLE :: OUTAU(:), help(:)
      integer :: ier, icolor, istr, icomgrp, ier1, i, my_pe_gr,
     .           mxdim, ns, ir
      logical, allocatable :: lhelp(:)

      if (nsteff < nprs) then
! collect from group leader pe
        icolor=MPI_UNDEFINED
        do istr=1,nstrai
          if (npts(istr) == 0) cycle
          if (my_pe == npesta(istr)) then
            icolor = 1
            write (0,*) ' my_pe, istr ',my_pe, istr
            exit
          end if
        end do
      else
! collect from all pes
        icolor = 1
      end if

      call mpi_barrier(mpi_comm_world,ier)

      call mpi_comm_split (mpi_comm_world,icolor,my_pe,icomgrp,ier)
      
      if (icolor == 1) then
        ALLOCATE (OUTAU(NOUTAU))
        CALL EIRENE_WRITE_COUTAU (OUTAU, IUNOUT)
        
        mxdim = max(noutau,nidv,nids,3*nsigci,nsigvi,nsigsi)
        allocate (help(mxdim))
        
        call mpi_comm_rank(mpi_comm_world,my_pe_gr,ier)

        CALL MPI_REDUCE(OUTAU,help,NOUTAU,
     .                  mpi_double_precision,mpi_sum,0,icomgrp,ier)

        if (my_pe == 0) CALL EIRENE_READ_COUTAU (help, IUNOUT)
        DEALLOCATE (OUTAU)

        mxdim = (max(NMOLI,NATMI,NIONI,NPHOTI,NPLSI)+1)*(NSTRA+1)
        allocate (lhelp(mxdim))

        call mpi_reduce(LOGMOL,lhelp,(NMOLI+1)*(NSTRA+1),
     .                  mpi_logical,mpi_LOR,0,icomgrp,ier)
     	if (my_pe == 0) 
     .    LOGMOL(0:nmoli,0:nstra) = 
     .      reshape(lhelp(1:(NMOLI+1)*(NSTRA+1)),(/nmoli+1,nstra+1/))

        call mpi_reduce(LOGATM,lhelp,(NATMI+1)*(NSTRA+1),
     .                  mpi_logical,mpi_LOR,0,icomgrp,ier)
     	if (my_pe == 0) 
     .    LOGATM(0:natmi,0:nstra) = 
     .      reshape(lhelp(1:(NATMI+1)*(NSTRA+1)),(/natmi+1,nstra+1/))

        call mpi_reduce(LOGION,lhelp,(NIONI+1)*(NSTRA+1),
     .                  mpi_logical,mpi_LOR,0,icomgrp,ier)
     	if (my_pe == 0) 
     .    LOGION(0:nioni,0:nstra) = 
     .      reshape(lhelp(1:(NIONI+1)*(NSTRA+1)),(/nIONi+1,nstra+1/))

        call mpi_reduce(LOGPHOT,lhelp,(NPHOTI+1)*(NSTRA+1),
     .                  mpi_logical,mpi_LOR,0,icomgrp,ier)
     	if (my_pe == 0) 
     .    LOGPHOT(0:nphoti,0:nstra) = 
     .      reshape(lhelp(1:(NPHOTI+1)*(NSTRA+1)),(/nPHOTi+1,nstra+1/))

        call mpi_reduce(LOGPLS,lhelp,(NPLSI+1)*(NSTRA+1),
     .                  mpi_logical,mpi_LOR,0,icomgrp,ier)
     	if (my_pe == 0) 
     .    LOGPLS(0:nplsi,0:nstra) = 
     .      reshape(lhelp(1:(NPLSI+1)*(NSTRA+1)),(/nPLSi+1,nstra+1/))

        deallocate(lhelp)

        if (nsmstra > 0) then
c  volume averaged output tallies

	do ir = 1, nrtal
          CALL MPI_REDUCE(SMESTV(1:nidv,ir),help,NIDV,
     .                  mpi_double_precision,mpi_sum,0,icomgrp,ier1)
          if (my_pe == 0) SMESTV(1:nidv,ir) = help(1:nidv)
        end do
c  surface averaged output tallies

	do ir = 1, nlmpgs
          CALL MPI_REDUCE(SMESTS(1:nids,ir),help,NIDS,
     .                  mpi_double_precision,mpi_sum,0,icomgrp,ier1)
          if (my_pe == 0) SMESTS(1:nids,ir) = help(1:nids)
        end do

c  spectrally resolved output tallies
        DO I=1,NADSPC
          ns = SMESTL(I)%PSPC%NSPC
          CALL MPI_REDUCE(SMESTL(I)%PSPC%SPC(0:ns+1),help,
     .                    SMESTL(I)%PSPC%NSPC+2,
     .                    MPI_DOUBLE_PRECISION,MPI_SUM,0,icomgrp,IER1)
          if (my_pe == 0) SMESTL(I)%PSPC%SPC(0:ns+1) = help(1:ns+2)

          CALL MPI_REDUCE(SMESTL(I)%PSPC%SPCS,help(1),
     .                    1,MPI_DOUBLE_PRECISION,MPI_SUM,0,icomgrp,IER1)
          if (my_pe == 0) SMESTL(I)%PSPC%SPCS = help(1)

c  variances of spectrally resolved output tallies
          if (nsigi_spc > 0) then
            call mpi_reduce(smestl(i)%pspc%GG(0:ns+1),help,
     .                      smestl(i)%pspc%nspc+2,
     .                      mpi_double_precision,mpi_sum,0,icomgrp,ier1)
            if (my_pe == 0) SMESTL(I)%PSPC%GG(0:ns+1) = help(1:ns+2)

            call mpi_reduce(smestl(i)%pspc%STV(0:ns+1),help,
     .                      smestl(i)%pspc%nspc+2,
     .                      mpi_double_precision,mpi_sum,0,icomgrp,ier1)
            if (my_pe == 0) SMESTL(I)%PSPC%STV(0:ns+1) = help(1:ns+2)

            call mpi_reduce(smestl(i)%pspc%stvs,
     .                      help(1),1,
     .                      mpi_double_precision,mpi_sum,0,icomgrp,ier1)
            if (my_pe == 0) SMESTL(I)%PSPC%STVS = help(1)
  
            call mpi_reduce(smestl(i)%pspc%ggs,
     .                      help(1),1,
     .                      mpi_double_precision,mpi_sum,0,icomgrp,ier1)
            if (my_pe == 0) SMESTL(I)%PSPC%GGS = help(1)
          end if
        END DO

        end if

        IF (NSIGCI > 0) THEN
c  covariances, volume tallies
          DO IR=1,NSBOX_TAL
            CALL MPI_REDUCE(reshape(STVC(0:2,1:NSIGCI,IR),(/3*nsigci/)),
     .                      help,3*NSIGCI,mpi_double_precision,
     .                      mpi_sum,0,icomgrp,ier1)
            if (my_pe == 0) STVC(0:2,1:NSIGCI,IR) = 
     .                      RESHAPE(help(1:3*nsigci),(/3,nsigci/))

          ENDDO
c  covariances, surface tallies
          CALL MPI_REDUCE(reshape(STVCS(0:2,1:NSIGCI),(/3*nsigci/)),
     .                    help,3*NSIGCI,
     .                    mpi_double_precision,mpi_sum,0,icomgrp,ier1)
          if (my_pe == 0) STVCS(0:2,1:NSIGCI) = 
     .	                  RESHAPE(help(1:3*nsigci),(/3,nsigci/))
        END IF

c  variances of volumetric output tallies

        IF (NSIGVI > 0) THEN
          DO IR=1,NSBOX_TAL
            CALL MPI_REDUCE(STV(1:NSIGVI,ir),help,
     .           NSIGVI,MPI_DOUBLE_PRECISION,MPI_SUM,0,ICOMGRP,IER1)
            if (my_pe == 0) STV(1:NSIGVI,ir) = help(1:nsigvi)

            CALL MPI_REDUCE(EE(1:NSIGVI,IR),help,
     .           NSIGVI,MPI_DOUBLE_PRECISION,MPI_SUM,0,ICOMGRP,IER1)
            if (my_pe == 0) EE(1:NSIGVI,ir) = help(1:nsigvi)
          ENDDO

          CALL MPI_REDUCE(STVS,help,NSIGVI,
     .                    MPI_DOUBLE_PRECISION,MPI_SUM,0,ICOMGRP,IER1)
          if (my_pe == 0) STVS(1:NSIGVI) = help(1:nsigvi)

          CALL MPI_REDUCE(EES,help,NSIGVI,
     .                    MPI_DOUBLE_PRECISION,MPI_SUM,0,ICOMGRP,IER1)
          if (my_pe == 0) EES(1:NSIGVI) = help(1:nsigvi)
        END IF


        IF (NSIGSI > 0) THEN
c  variances of surface averaged output tallies

          do ir=1,nlimps
            CALL MPI_REDUCE(STVW(1:NSIGSI,IR),help,NSIGSI,
     .                      MPI_DOUBLE_PRECISION,MPI_SUM,0,ICOMGRP,IER1)
            if (my_pe == 0) STVW(1:NSIGSI,IR) = help(1:nsigsi)
            
            CALL MPI_REDUCE(FF(1:NSIGSI,IR),help,NSIGSI,
     .                      MPI_DOUBLE_PRECISION,MPI_SUM,0,ICOMGRP,IER1)
            if (my_pe == 0) FF(1:NSIGSI,IR) = help(1:nsigsi)
     	  end do
          CALL MPI_REDUCE(STVWS,help,NSIGSI,
     .                    MPI_DOUBLE_PRECISION,MPI_SUM,0,ICOMGRP,IER1)
          if (my_pe == 0) STVWS(1:NSIGSI) = help(1:nsigsi)

          CALL MPI_REDUCE(FFS,help,NSIGSI,
     .                    MPI_DOUBLE_PRECISION,MPI_SUM,0,ICOMGRP,IER1)
          if (my_pe == 0) FFS(1:NSIGSI) = help(1:nsigsi)
        END IF
        
        deallocate (help)

!csw 08jul2011 missing comm_free fixed (invalid comms from icolor/=1 are actually not allocated, is this conform to MPI standard?)
        call mpi_comm_free(icomgrp,ier)

      end if ! icolor=1
      
      call mpi_barrier(mpi_comm_world,ier)
      
      return
      end subroutine eirene_collect_coutau
