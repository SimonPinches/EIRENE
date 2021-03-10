!     HJL Container for openmp persistent variables and openMP specific routines

      MODULE EIRMOD_OPENMP

      USE EIRMOD_COMXS
      USE EIRMOD_REFLEC
      USE EIRMOD_STATIS
      USE EIRMOD_CLAST
      USE EIRMOD_CFPLK, ONLY: FNUIAR
!HJL ADDED FOR ROUTINE MOVE FROM EIRMOD_MCARLO
      USE EIRMOD_COMSPL
      USE EIRMOD_COMSOU
      USE EIRMOD_COMPRT
      USE EIRMOD_CSDVI
      USE EIRMOD_CGRID
      USE EIRMOD_CUPD
      USE EIRMOD_PARMMOD
      USE EIRMOD_PRECISION
      USE EIRMOD_SPUTER, ONLY: ETH,Q,M2M1,ES,ETF
      
#ifdef USE_OPENMP     
      USE OMP_LIB
#endif
      
      IMPLICIT NONE

      PRIVATE

      PUBLIC  :: EIRENE_INIT_OPENMP, EIRENE_DEALLOCATE_OPENMP
      
      INTEGER, PUBLIC, SAVE :: EIRENE_ITHREAD, EIRENE_NTHREADS

!$OMP THREADPRIVATE(EIRENE_ITHREAD,EIRENE_NTHREADS)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ym shared variables used as buffer to initialize private pointer variables
c     ym copyin does not work for allocatable pointer arrrays      

      INTEGER :: BISTRA
      INTEGER,TARGET,ALLOCATABLE :: BISDVI(:)      
      INTEGER,TARGET,ALLOCATABLE :: BIPSTD(:),BICMSPL(:)
      REAL(DP),TARGET,ALLOCATABLE :: BRPST(:),BRCMSPL(:),
     .     BXSTOR(:,:),BXSTORV(:)
      LOGICAL,TARGET,ALLOCATABLE :: BLCMSOU(:,:)
      REAL(DP), DIMENSION(28,0:11) :: BETH,BQ,BM2M1,BETF
      REAL(DP), DIMENSION(28) :: BES

      CONTAINS

!     Initialisation routine that calls allocate and array copy routines
!     The preprocessor directives differentiate between the parallel region
!     being created in eirmod_mcarlo and being created by an external coupled code.
!      
!     There is also a parallel region in eirene_main that allows for compilation
!     of the standalone code while using the EXT_OPENMP flag, this is not seen
!     by a coupled code as it is before the coupling ENTRY point
      
      SUBROUTINE EIRENE_INIT_OPENMP
#ifndef USE_EXT_OPENMP      
!$OMP PARALLEL DEFAULT(SHARED)
#endif

#ifdef USE_OPENMP     
      EIRENE_ITHREAD  = OMP_GET_THREAD_NUM()
      EIRENE_NTHREADS = OMP_GET_NUM_THREADS()
!$OMP MASTER
      write(iunout,*)"Compiled with OpenMP - NTHREADS =",EIRENE_NTHREADS
!$OMP END MASTER
#else
      EIRENE_ITHREAD  = 0
      EIRENE_NTHREADS = 1
#endif
      
!$OMP BARRIER

#ifdef USE_EXT_OPENMP      
      IF(EIRENE_ITHREAD==0) THEN
         CALL EIRENE_INIT_BUFFERS_OPENMP()
      ELSE
#endif
#ifdef USE_OPENMP
         CALL EIRENE_ALLOCATE_OPENMP()
#endif
#ifdef USE_EXT_OPENMP
      ENDIF

!$OMP BARRIER

      IF(EIRENE_ITHREAD>0) THEN
         CALL EIRENE_COPY_BUFFERS_OPENMP()
      ENDIF
      
!$OMP BARRIER

      IF(EIRENE_ITHREAD==0) THEN
         CALL EIRENE_DEALLOCATE_BUFFERS_OPENMP()
      ENDIF
#else
!$OMP END PARALLEL
#endif

      END SUBROUTINE EIRENE_INIT_OPENMP

      
!!! Allocate THREAD_PRIVATE arrays for non master threads      
      SUBROUTINE EIRENE_ALLOCATE_OPENMP

      IF(EIRENE_ITHREAD==0) RETURN
      
      ALLOCATE (ISDVI(MSDVI))              
      ALLOCATE (LCMSOU(14,NSTRA))
      ALLOCATE (TIMINT(NRADS))
      ALLOCATE (TIMPOL(N1STS,N2NDPLGS))
      ALLOCATE (NTIM(NRADS))
      ALLOCATE (IIMPOL(N1STS,N2NDPLGS))
      ALLOCATE (IIMINT(NRADS))
      
      ALLOCATE (RPST(NPARTC))
      ALLOCATE (IPSTD(MPARTC+1))
      ALLOCATE (RCMSPL(NCMSPL))
      ALLOCATE (ICMSPL(MCMSPL))
      
      ALLOCATE (ALPD(N2ND))
      ALLOCATE (BLPD(N3RD))
      ALLOCATE (CLPD(N2ND+N3RD))
      
      ALLOCATE (JUPC(N2ND))
      ALLOCATE (KUPC(N3RD))
      ALLOCATE (NUPC(N2ND+N3RD))
      ALLOCATE (NCOUNP(N2ND+N3RD))
      ALLOCATE (NCOUNT(N2ND+N3RD))
      ALLOCATE (LUPC(N2ND))
      ALLOCATE (MUPC(N2ND))
      
!!      ALLOCATE (RCGRID(NCGRD))

      ALLOCATE(EREDUC(NSPZ,0:NLIMPS))
      ALLOCATE(FREDUC(NSPZ,0:NLIMPS))
      ALLOCATE(IREDUC(NSPZ,0:NLIMPS))
      
      AllOCATE (IIND(NRTAL))
      ALLOCATE (XSTOR(MSTOR1,MSTOR2))       
      ALLOCATE (XSTORV(NSTORV))

c     ym arrays from eirmod_clast - these are not pointers
      ALLOCATE (XCMEAN(NRCX))
      ALLOCATE (SGCVMX(NRCX))
      ALLOCATE (XEMEAN(NREL))
      ALLOCATE (SGEVMX(NREL))
      ALLOCATE (XPMEAN(NRPI))
      ALLOCATE (SGPVMX(NRPI))

      ALLOCATE (NCMEAN(NRCX))
      ALLOCATE (IFLRCX(NRCX))
      ALLOCATE (NEMEAN(NREL))
      ALLOCATE (IFLREL(NREL))
      ALLOCATE (NPMEAN(NRPI))
      ALLOCATE (IFLRPI(NRPI))
      
c     ym test iter
      ALLOCATE (RSPLST(NPARTC,MAXLEVEL))
      ALLOCATE (ISPLST(MPARTC,MAXLEVEL))

      RSPLST=0._dp
      ISPLST=0

c     ym make sure the clast variables do not take exotic values      
      call eirene_init_clast

      ALLOCATE (FNUIAR(NPLS))
      FNUIAR=0._dp
      
      NCLMT     => ISDVI(8)
      NCLMTS    => ISDVI(9)
      NWLMT     => ISDVI(10)         
      NWLMTS    => ISDVI(11)
      ICLMT     => ISDVI(12+2*NSD+2*NSDW+NCV+NRTAL :
     .     11+2*NSD+2*NSDW+NCV+2*NRTAL)
      IMETCL    => ISDVI(12+2*NSD+2*NSDW+NCV :
     .     11+2*NSD+2*NSDW+NCV+NRTAL)
      IMETWL    => ISDVI(12+2*NSD+2*NSDW+NCV+2*NRTAL :
     .     11+2*NSD+2*NSDW+NCV+2*NRTAL+NLIMPS)
      IWLMT    => ISDVI(12+2*NSD+2*NSDW+NCV+2*NRTAL+NLIMPS : MSDVI)
      
      ISPZ   => IPSTD( 9)
!HJL      NT3RD  => ICGRID( 8)
      MRSURF => IPSTD(10)
      MPSURF => IPSTD(11)
      MTSURF => IPSTD(12)
      MASURF => IPSTD(13)
      MSURF  => IPSTD(14)
      NLRAY  => LCMSOU(14,:)
      
      RPSTT => RPST
      
      X0     => RPST( 1)
      Y0     => RPST( 2)
      Z0     => RPST( 3)
      VEL    => RPST( 4)
      VELX   => RPST( 5)
      VELY   => RPST( 6)
      VELZ   => RPST( 7)
      E0     => RPST( 8)
      WEIGHT => RPST( 9)
      TIME   => RPST(10)
      PHI    => RPST(11)
      
      XGENER => RPST(12)
      
      IPST  => IPSTD(2:MPARTC+1)
      IPSTT => IPSTD(1:MPARTT)

      NPANU  => IPSTD(1)
      IPOLG  => IPSTD(2)
      IPERID => IPSTD(3)
      NCELL  => IPSTD(4)
      ITIME  => IPSTD(5)
      IFPATH => IPSTD(6)
      IUPDTE => IPSTD(7)
      !HJL Put this back
      ISTRA  => IPSTD( 8)
      ISPZ   => IPSTD(9)
      
      MSURFG => IPSTD(15)
      WMINV  => RCMSPL(1)
      WMINS  => RCMSPL(2)
      WMINC  => RCMSPL(3)
      WMINL  => RCMSPL(4)
      SPLPAR => RCMSPL(5)
      RNUMB  => RCMSPL(6:5+ N1ST+N2ND+N3RD+NLIM)
      PRMSPL => RCMSPL(6+   N1ST+N2ND+N3RD+NLIM : NCMSPL)
      
      MAXLEV => ICMSPL(1)
      NLEVEL => ICMSPL(2)
      MAXRAD => ICMSPL(3)
      MAXPOL => ICMSPL(4)
      MAXTOR => ICMSPL(5)
      MAXADD => ICMSPL(6)
      
      NODES  => ICMSPL(7:6+ MAXLEVEL)
      NSSPL  => ICMSPL(7  + MAXLEVEL:MCMSPL)
      
      SIGVCX => XSTOR(:,1)
      SIGVPI => XSTOR(:,2)
      SIGVEI => XSTOR(:,3)
      SIGVEL => XSTOR(:,4)
      SIGVPH => XSTOR(:,22)
      
      ESIGCX => XSTOR(:,5:6)
      ESIGPI => XSTOR(:,7:11)
      ESIGEI => XSTOR(:,12:16)
      ESIGEL => XSTOR(:,17:18)
      ESIGPH => XSTOR(:,23:24)
      
      VSIGCX => XSTOR(:,19)
      VSIGPI => XSTOR(:,20)
      VSIGEL => XSTOR(:,21)
      
      SIGCXT  => XSTORV(1)
      SIGPIT  => XSTORV(2)
      SIGEIT  => XSTORV(3)
      SIGELT  => XSTORV(4)
      SIGPHT  => XSTORV(5)
      SIGTOT  => XSTORV(6)
      SIGBGK  => XSTORV(7)
      ZMFPI   => XSTORV(8)
      
      END SUBROUTINE EIRENE_ALLOCATE_OPENMP


!!! Deallocate THREADPRIVATE arrays
      SUBROUTINE EIRENE_DEALLOCATE_OPENMP
      

      IF(EIRENE_ITHREAD == 0 ) RETURN
      
      DEALLOCATE(ISDVI)
      DEALLOCATE(TIMINT)
      DEALLOCATE(TIMPOL)
      DEALLOCATE(NTIM)
      DEALLOCATE(IIMPOL)
      DEALLOCATE(IIMINT)
      
      DEALLOCATE(RPST)
      DEALLOCATE(IPSTD)
      DEALLOCATE(RCMSPL)
      DEALLOCATE(ICMSPL)
      
      DEALLOCATE(ALPD)
      DEALLOCATE(BLPD)
      DEALLOCATE(CLPD)
      
      DEALLOCATE (JUPC)
      DEALLOCATE (KUPC)
      DEALLOCATE (LUPC)
      DEALLOCATE (MUPC)
      DEALLOCATE(NUPC)
      DEALLOCATE(NCOUNP)
      DEALLOCATE(NCOUNT)
      
      DEALLOCATE(EREDUC)
      DEALLOCATE(FREDUC)
      DEALLOCATE(IREDUC)
      
      DEALLOCATE(XSTOR)
      DEALLOCATE(XSTORV)
      DEALLOCATE(LCMSOU)
      
      DEALLOCATE(IIND)
      
      DEALLOCATE (XCMEAN)
      DEALLOCATE (SGCVMX)
      DEALLOCATE (XEMEAN)
      DEALLOCATE (SGEVMX)
      DEALLOCATE (XPMEAN)
      DEALLOCATE (SGPVMX)
      
      DEALLOCATE (NCMEAN)
      DEALLOCATE (IFLRCX)
      DEALLOCATE (NEMEAN)
      DEALLOCATE (IFLREL)
      DEALLOCATE (NPMEAN)
      DEALLOCATE (IFLRPI)
      
      DEALLOCATE (FNUIAR)
      
      DEALLOCATE (RSPLST)
      DEALLOCATE (ISPLST)

      END SUBROUTINE EIRENE_DEALLOCATE_OPENMP

      
!!! Allocate shared arrays and copy in the data from master
      SUBROUTINE EIRENE_INIT_BUFFERS_OPENMP

      ALLOCATE (BISDVI(MSDVI))              
      BISDVI=ISDVI
      ALLOCATE (BIPSTD(MPARTC+1))      
      BIPSTD=IPSTD
      ALLOCATE (BRPST(NPARTC))     
      BRPST=RPST
      ALLOCATE (BRCMSPL(NCMSPL))
      BRCMSPL=RCMSPL
      ALLOCATE (BICMSPL(MCMSPL))
      BICMSPL=ICMSPL
      ALLOCATE (BLCMSOU(14,NSTRA))
      BLCMSOU=LCMSOU
      ALLOCATE (BXSTOR(MSTOR1,MSTOR2))       
      BXSTOR=XSTOR
      ALLOCATE (BXSTORV(NSTORV))
      BXSTORV=XSTORV
      ! New for complete threading
      BISTRA=ISTRA
      BETH=ETH
      BQ=Q
      BM2M1=M2M1
      BETF=ETF
      BES=ES

      END SUBROUTINE EIRENE_INIT_BUFFERS_OPENMP


!!! Copy the shared buffers into the thread private arrays
      SUBROUTINE EIRENE_COPY_BUFFERS_OPENMP
      
      ISDVI=BISDVI             
      IPSTD=BIPSTD
      RPST=BRPST
      RCMSPL=BRCMSPL
      ICMSPL=BICMSPL
      LCMSOU=BLCMSOU
      XSTOR=BXSTOR
      XSTORV=BXSTORV
!      RCGRID=BRCGRID
!     New for complete threading
      ISTRA=BISTRA
      ETH=BETH
      Q=BQ
      M2M1=BM2M1
      ETF=BETF
      ES=BES
      
      END SUBROUTINE EIRENE_COPY_BUFFERS_OPENMP

      
 !!! Allocate shared arrays and copy in the data from master
      SUBROUTINE EIRENE_DEALLOCATE_BUFFERS_OPENMP

      DEALLOCATE(BISDVI)              
      DEALLOCATE(BIPSTD)      
      DEALLOCATE(BRPST)     
      DEALLOCATE(BRCMSPL)
      DEALLOCATE(BICMSPL)
      DEALLOCATE(BLCMSOU)
      DEALLOCATE(BXSTOR)       
      DEALLOCATE(BXSTORV)

      END SUBROUTINE EIRENE_DEALLOCATE_BUFFERS_OPENMP    
      
      
      END MODULE EIRMOD_OPENMP
