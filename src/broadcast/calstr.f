!pb 300806  reduce commands for WTOTE and EELFI added
!           reduction of spectrum data corrected
!pb 181206  group management changed
!pb 110707  calls to mpi_reduce corrected
!pb 060309  mpi_real8 --> mpi_double_precision
!pb 090309  rewritten to use automatic arrays as output buffer in mpi_reduce
!pb 090309  loops reorganized
!pb 270309  typos corrected
!sw 091112  added support for csdvi_cop and csdvi_bgk

cdr Nov. 15:  comments needed. copv tallies: variances for coupling ??
cdr                            to be checked again after changes in 2013
cdr dec. 15:  eppli: now resolved wrt. species index ipls, added
cdr july 17:  comments re. call to user routine: calstr_usr.

      SUBROUTINE EIRENE_CALSTR
cdr
c
c  called from MCARLO.f, from within strata loop, at the end of each stratum,
c  if there are more processors than active strata.
c
c  Unclear: if more strata than processors: is it then excluded that still
c           there may be strata with more than one processor dealing with them?
c          
c  Purpose:
c   collect data from processors belonging to one particular stratum istra (COMPRT)
c   put merged data for output tallies for stratum istra then on:  my_pe_gr=0
c
c
cdr
c   input:
c       istra  (stratum number, from common COMPRT)
c       npesta(istra)  : number of master processor for stratum istra
c       npestr(istra)  : total no. of processors working on stratum istra
c
c   results:
c       npean, npeen:   the processors in the range npean,...,npeen work on stratum istra

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_CSPEZ
      USE EIRMOD_COMPRT
      USE EIRMOD_CPES
      USE EIRMOD_CSDVI
      USE EIRMOD_CSDVI_COP
      USE EIRMOD_CSDVI_BGK
      USE EIRMOD_COUTAU
      USE EIRMOD_CSTEP
      IMPLICIT NONE

C
      INCLUDE 'mpif.h'
      real(dp), allocatable :: help(:), helpest(:), helpv(:), dummyv(:)
      real(dp) :: helpa(0:natm), helpm(0:nmol), helpi(0:nion),
     .            helpp(0:npls), helpph(0:nphot),
     .            helps(nlmpgs+1), helpc
C     real(dp) :: dummyv(nrtal+1), dummys(nlmpgs+1)
      real(dp) :: dummys(nlmpgs+1)
      real(dp), allocatable :: dummyw(:), helpw(:)
      integer :: icomgrp(0:nstra)
      integer :: ier1, ier, ir, npean, npeen, i, mpicw, ispc, my_pe_gr,
     .           mxdim, ns, j
      logical, allocatable :: lhelp(:)
      logical :: lhelpa(0:natm),lhelpm(0:nmol), lhelpi(0:nion),
     .           lhelpp(0:npls), lhelpph(0:nphot)
c
c  range of processors working on stratum ISTRA
      npean = npesta(istra)
      npeen = npesta(istra)+npestr(istra)-1
c
      call mpi_comm_group (mpi_comm_world,mpicw,ier)
      call mpi_comm_split (mpi_comm_world,istra,my_pe-npesta(istra),
     .                     icomgrp(istra),ier)


      if(      count( procforstra(istra,0:nprs-1) ) >1
     .   .and.        procforstra(istra,my_pe)) then
CDR  more than one single processor was active on this stratum ISTRA,
CDR  and my_pe is one of them

        call mpi_barrier(icomgrp(istra),ier)
c
c  my_pe_gr=0 indicates: my_pe is the master processor for istra
c   
        my_pe_gr = my_pe-npesta(istra)

        mxdim = max(nvoltl,nsrftl,nsd,nsdw,
     .              nmoli+1,natmi+1,nioni+1,nphoti+1,nplsi+1)

        allocate (help(mxdim))
     	
        call mpi_reduce(WTOTM(0,istra),helpm,nmoli+1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	    if (my_pe_gr==0) WTOTM(0:nmoli,istra) = helpm(0:nmoli)

        call mpi_reduce(WTOTA(0,istra),helpa,natmi+1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	    if (my_pe_gr==0) WTOTA(0:natmi,istra) = helpa(0:natmi)

        call mpi_reduce(WTOTI(0,istra),helpi,nioni+1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	    if (my_pe_gr==0) WTOTI(0:nioni,istra) = helpi(00:nioni)

        call mpi_reduce(WTOTP(0,istra),helpp,nplsi+1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	    if (my_pe_gr==0) WTOTP(0:nplsi,istra) = helpp(0:nplsi)

        call mpi_reduce(WTOTE(istra),helpc,1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	    if (my_pe_gr==0) WTOTE(istra) = helpc
C
C
        call mpi_reduce(XMCP(istra),helpc,1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	    if (my_pe_gr==0) XMCP(istra) = helpc

csw 19feb2013, added XMCT
        call mpi_reduce(XMCT(istra),helpc,1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	    if (my_pe_gr==0) XMCT(istra) = helpc
csw

        call mpi_reduce(PTRASH(istra),helpc,1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	    if (my_pe_gr==0) PTRASH(istra) = helpc

        call mpi_reduce(ETRASH(istra),helpc,1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
        if (my_pe_gr==0) ETRASH(istra) = helpc

        call mpi_reduce(ETOTA(istra),helpc,1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	    if (my_pe_gr==0) ETOTA(istra) = helpc

        call mpi_reduce(ETOTM(istra),helpc,1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	    if (my_pe_gr==0) ETOTM(istra) = helpc

        call mpi_reduce(ETOTI(istra),helpc,1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	    if (my_pe_gr==0) ETOTI(istra) = helpc

        call mpi_reduce(ETOTP(istra),helpc,1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	    if (my_pe_gr==0) ETOTP(istra) = helpc

        call mpi_reduce(EELFI(0,istra),helpi,nioni+1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	    if (my_pe_gr==0) EELFI(0:nioni,istra) = helpi(0:nioni)

c  particle balance tallies:  from bulk (ipls) to species a,m,i,ph,pl
        call mpi_reduce(PPATI(0,istra),helpa,natmi+1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	    if (my_pe_gr==0) PPATI(0:natmi,istra) = helpa(0:natmi)

        call mpi_reduce(PPMLI(0,istra),helpm,nmoli+1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	    if (my_pe_gr==0) PPMLI(0:nmoli,istra) = helpm(0:nmoli)

        call mpi_reduce(PPIOI(0,istra),helpi,nioni+1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	    if (my_pe_gr==0) PPIOI(0:nioni,istra) = helpi(0:nioni)

        if (nphoti > 0) then
          call mpi_reduce(PPPHTI(0,istra),helpph,nphoti+1,
     .         mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	      if (my_pe_gr==0) PPPHTI(0:nphoti,istra) = helpph(0:nphoti)
        end if

        call mpi_reduce(PPPLI(0,istra),helpp,nplsi+1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	    if (my_pe_gr==0) PPPLI(0:nplsi,istra) = helpp(0:nplsi)

c  energy balance tallies:  from bulk (ipls) to species a,m,i,ph,pl
        call mpi_reduce(EPATI(istra),helpc,1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	    if (my_pe_gr==0) EPATI(istra) = helpc

        call mpi_reduce(EPMLI(istra),helpc,1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	    if (my_pe_gr==0) EPMLI(istra) = helpc

        call mpi_reduce(EPIOI(istra),helpc,1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	    if (my_pe_gr==0) EPIOI(istra) = helpc

        call mpi_reduce(EPPHTI(istra),helpc,1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	    if (my_pe_gr==0) EPPHTI(istra) = helpc

        call mpi_reduce(EPPLI(0,istra),helpp,nplsi+1,
     .       mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	    if (my_pe_gr==0) EPPLI(0:nplsi,istra) = helpp(0:nplsi)

C
C
c  all other volume averaged tallies: estimv
        allocate (helpv(nrtal+1), dummyv(nrtal+1))
        do ir=1,nvoltl
          dummyv(1:nrtal) = estimv(ir,1:nrtal)
          call mpi_reduce(dummyv,helpv,nrtal,
     .         mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	      if (my_pe_gr==0) estimv(ir,1:nrtal) = helpv(1:nrtal)
        end do

c  all surface averaged tallies: estims
	    do ir=1,nsrftl
          dummys(1:nlmpgs) = estims(ir,1:nlmpgs)
          call mpi_reduce(dummys,helps,nlmpgs,
     .         mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	      if (my_pe_gr==0) estims(ir,1:nlmpgs) = helps(1:nlmpgs)
        end do

c   energy resolved ("spectra") tallies
        do ispc=1,nadspc
          ns = estiml(ispc)%pspc%nspc
          allocate (helpest(ns+2))
          call mpi_reduce(estiml(ispc)%pspc%spc,helpest,
     .                    estiml(ispc)%pspc%nspc+2,
     .         mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
          if (my_pe_gr==0) estiml(ispc)%pspc%spc(0:ns+1)=helpest(1:ns+2)

c  standard deviation of energy resolved "spectra"
          if (nsigi_spc > 0) then
            call mpi_reduce(estiml(ispc)%pspc%sdv,helpest,
     .                      estiml(ispc)%pspc%nspc+2,
     .           mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
            if (my_pe_gr==0)
     .      estiml(ispc)%pspc%sdv(0:ns+1) = helpest(1:ns+2)

            call mpi_reduce(estiml(ispc)%pspc%sgm,helpest,
     .                      estiml(ispc)%pspc%nspc+2,
     .           mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
            if (my_pe_gr==0)
     .        estiml(ispc)%pspc%sgm(0:ns+1) = helpest(1:ns+2)

            call mpi_reduce(estiml(ispc)%pspc%sgms,helpest,1,
     .           mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
            if (my_pe_gr==0) estiml(ispc)%pspc%sgms = helpest(1)
          end if

          deallocate (helpest)
        end do   !nadspc

C  standard deviation of volume averaged tallies
        if (nsd > 0) then
          do ir=1,nsd
            dummyv(1:nrtal+1) = sdvi1(ir,1:nrtal+1)
            call mpi_reduce(dummyv,helpv,nrtal+1,
     .           mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	        if (my_pe_gr==0) sdvi1(ir,1:nrtal+1) = helpv(1:nrtal+1)
          end do
        end if

C  standard deviation of surface averaged tallies
	    if (nsdw > 0) then
          do ir=1,nsdw
            dummys(1:nlimps+1) = sdvi2(ir,1:nlimps+1)
            call mpi_reduce(dummys,helps,nlimps+1,
     .           mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	        if (my_pe_gr==0) sdvi2(ir,1:nlimps+1) = helps(1:nlimps+1)
          end do
	    end if

C  covariances between two volume averaged tallies
	    if (ncv > 0) then
          do i=0,2
            do j=1,ncv
              dummyv(1:nrtal) = sigmac(i,j,1:nrtal)
     	      call mpi_reduce(dummyv,helpv,nrtal,
     .             mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	          if (my_pe_gr==0) sigmac(i,j,1:nrtal) = helpv(1:nrtal)
            end do
          end do

          dummyv(1:ncv) = sgmcs(0,1:ncv)
          call mpi_reduce(dummyv,helpv,ncv,
     .         mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	      if (my_pe_gr==0) sgmcs(0,1:ncv) = helpv(1:ncv)

          dummyv(1:ncv) = sgmcs(1,1:ncv)
          call mpi_reduce(dummyv,helpv,ncv,
     .         mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
	      if (my_pe_gr==0) sgmcs(1,1:ncv) = helpv(1:ncv)

          dummyv(1:ncv) = sgmcs(2,1:ncv)
          call mpi_reduce(dummyv,helpv,ncv,
     .         mpi_double_precision,mpi_sum,0,icomgrp(istra),ier1)
     	  if (my_pe_gr==0) sgmcs(2,1:ncv) = helpv(1:ncv)
	    end if

csw 09nov2012 reduce csdvi_cop
        if (ncpv_stat > 0) then
          allocate(dummyw(max(nrtals,ncpv_stat)+1))
          allocate(helpw(max(nrtals,ncpv_stat)+1))

          do i=1,ncpv_stat
            dummyw(1:nrtals) = sigma_cop(i,1:nrtals)
            call mpi_reduce(dummyw,helpw,nrtals,
     .                      mpi_double_precision,mpi_sum,
     .                      0,icomgrp(istra),ier1)
            if(my_pe_gr==0) sigma_cop(i,1:nrtals) = helpw(1:nrtals)

            dummyw(1:nrtals) = stv_cop(i,1:nrtals)
            call mpi_reduce(dummyw,helpw,nrtals,
     .                      mpi_double_precision,mpi_sum,
     .                      0,icomgrp(istra),ier1)
            if(my_pe_gr==0) stv_cop(i,1:nrtals) = helpw(1:nrtals)

            dummyw(1:nrtals) = sdvia_cop(i,1:nrtals)
            call mpi_reduce(dummyw,helpw,nrtals,
     .                      mpi_double_precision,mpi_sum,
     .                      0,icomgrp(istra),ier1)
            if(my_pe_gr==0) sdvia_cop(i,1:nrtals) = helpw(1:nrtals)

            dummyw(1:nrtals) = ee_cop(i,1:nrtals)
            call mpi_reduce(dummyw,helpw,nrtals,
     .                      mpi_double_precision,mpi_sum,
     .                      0,icomgrp(istra),ier1)
            if(my_pe_gr==0) ee_cop(i,1:nrtals) = helpw(1:nrtals)
          enddo

          call mpi_reduce(sgms_cop(1:ncpv_stat),helpv,ncpv_stat,
     .                    mpi_double_precision,mpi_sum,
     .                    0,icomgrp(istra),ier1)
          if (my_pe_gr==0) sgms_cop(1:ncpv_stat) = helpv(1:ncpv_stat)

          call mpi_reduce(stvs_cop(1:ncpv_stat),helpv,ncpv_stat,
     .                    mpi_double_precision,mpi_sum,
     .                    0,icomgrp(istra),ier1)
          if (my_pe_gr==0) stvs_cop(1:ncpv_stat) = helpv(1:ncpv_stat)

          call mpi_reduce(ees_cop(1:ncpv_stat),helpv,ncpv_stat,
     .                    mpi_double_precision,mpi_sum,
     .                    0,icomgrp(istra),ier1)
          if (my_pe_gr==0) ees_cop(1:ncpv_stat) = helpv(1:ncpv_stat)

          deallocate(dummyw)
          deallocate(helpw)
        endif
csw

csw 09nov2012 reduce csdvi_bgk
        if (nbgv_stat > 0) then
          allocate(dummyw(max(nrtals,nbgv_stat)+1))
          allocate(helpw(max(nrtals,nbgv_stat)+1))

          do i=1,nbgv_stat
            dummyw(1:nrtals) = sigma_bgk(i,1:nrtals)
            call mpi_reduce(dummyw,helpw,nrtals,
     .                      mpi_double_precision,mpi_sum,
     .                      0,icomgrp(istra),ier1)
            if(my_pe_gr==0) sigma_bgk(i,1:nrtals) = helpw(1:nrtals)

            dummyw(1:nrtals) = stv_bgk(i,1:nrtals)
            call mpi_reduce(dummyw,helpw,nrtals,
     .                      mpi_double_precision,mpi_sum,
     .                      0,icomgrp(istra),ier1)
            if(my_pe_gr==0) stv_bgk(i,1:nrtals) = helpw(1:nrtals)

            dummyw(1:nrtals) = sdvia_bgk(i,1:nrtals)
            call mpi_reduce(dummyw,helpw,nrtals,
     .                      mpi_double_precision,mpi_sum,
     .                      0,icomgrp(istra),ier1)
            if(my_pe_gr==0) sdvia_bgk(i,1:nrtals) = helpw(1:nrtals)

            dummyw(1:nrtals) = ee_bgk(i,1:nrtals)
            call mpi_reduce(dummyw,helpw,nrtals,
     .                      mpi_double_precision,mpi_sum,
     .                      0,icomgrp(istra),ier1)
            if(my_pe_gr==0) ee_bgk(i,1:nrtals) = helpw(1:nrtals)
          enddo

          call mpi_reduce(sgms_bgk(1:nbgv_stat),helpv,nbgv_stat,
     .                    mpi_double_precision,mpi_sum,
     .                    0,icomgrp(istra),ier1)
          if (my_pe_gr==0) sgms_bgk(1:nbgv_stat) = helpv(1:nbgv_stat)

          call mpi_reduce(stvs_bgk(1:nbgv_stat),helpv,nbgv_stat,
     .                    mpi_double_precision,mpi_sum,
     .                    0,icomgrp(istra),ier1)
          if (my_pe_gr==0) stvs_bgk(1:nbgv_stat) = helpv(1:nbgv_stat)

          call mpi_reduce(ees_bgk(1:nbgv_stat),helpv,nbgv_stat,
     .                    mpi_double_precision,mpi_sum,
     .                    0,icomgrp(istra),ier1)
          if (my_pe_gr==0) ees_bgk(1:nbgv_stat) = helpv(1:nbgv_stat)

          deallocate(dummyw)
          deallocate(helpw)
        endif
csw

        deallocate(help)
	    mxdim = max (nmoli+1,natmi+1,nioni+1,nphoti+1,nplsi+1)

        allocate (lhelp(mxdim))

        call mpi_reduce(LOGMOL(0,ISTRA),lhelpm,NMOLI+1,
     .       mpi_logical,mpi_LOR,0,icomgrp(istra),ier1)
	    if (my_pe_gr==0) LOGMOL(0:nmoli,ISTRA) = lhelpm(0:nmoli)

        call mpi_reduce(LOGATM(0,ISTRA),lhelpa,NATMI+1,
     .       mpi_logical,mpi_LOR,0,icomgrp(istra),ier1)
        if (my_pe_gr==0) LOGATM(0:natmi,ISTRA) = lhelpa(0:natmi)

        call mpi_reduce(LOGION(0,ISTRA),lhelpi,NIONI+1,
     .       mpi_logical,mpi_LOR,0,icomgrp(istra),ier1)
	    if (my_pe_gr==0) LOGION(0:nioni,ISTRA) = lhelpi(0:nioni)

        if (nphoti > 0) then
          call mpi_reduce(LOGPHOT(0,ISTRA),lhelpph,NPHOTI+1,
     .         mpi_logical,mpi_LOR,0,icomgrp(istra),ier1)
	      if (my_pe_gr==0) LOGPHOT(0:nphoti,ISTRA) = lhelpph(0:nphoti)
        end if

        call mpi_reduce(LOGPLS(0,ISTRA),lhelpp,NPLSI+1,
     .       mpi_logical,mpi_LOR,0,icomgrp(istra),ier1)
	    if (my_pe_gr==0) LOGPLS(0:nplsi,ISTRA) = lhelpp(0:nplsi)
	
        deallocate(lhelp)
c
c  collect user or case specific information from all Pes that worked on
c  stratum no. ISTRA.  Depends on ...usr.f  or ...cop.f routines.
c  Strictly there should also be an analogue  call to eirene_calstr_cop.f 

        call mpi_barrier(icomgrp(istra),ier)
        call eirene_calstr_usr (my_pe_gr, icomgrp(istra))

      endif

      call mpi_comm_free (icomgrp(istra),ier)
      call mpi_group_free(mpicw,ier)
      call mpi_barrier(mpi_comm_world,ier)

      if (allocated(helpv)) deallocate (helpv)
      if (allocated(dummyv)) deallocate (dummyv)
      RETURN
      END
