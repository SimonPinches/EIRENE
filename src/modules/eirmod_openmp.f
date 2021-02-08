!     HJL Container for openmp persistent variables and openMP specific routines

      MODULE EIRMOD_OPENMP

      USE OMP_LIB
      
      IMPLICIT NONE

      PRIVATE

      PUBLIC  :: EIRENE_INIT_OPENMP, EIRENE_DEALLOCATE_OPENMP
      
      INTEGER, PUBLIC, SAVE :: ITHREAD, NTHREADS

!$OMP THREADPRIVATE(ITHREAD,NTHREADS)

      CONTAINS

! Initialisation routine that also allocates THREADPRIVATE arrays
      SUBROUTINE EIRENE_INIT_OPENMP
      
!$OMP PARALLEL
      ITHREAD  = OMP_GET_THREAD_NUM()
      NTHREADS = OMP_GET_NUM_THREADS()

      CALL EIRENE_ALLOCATE_OPENMP()      
!$OMP END PARALLEL

      
      END SUBROUTINE EIRENE_INIT_OPENMP
      
      SUBROUTINE EIRENE_ALLOCATE_OPENMP

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

      IF(ITHREAD > 0) THEN

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
         
         ALLOCATE (RCGRID(NCGRD))

         ALLOCATE(EREDUC(NSPZ,0:NLIMPS))
         ALLOCATE(FREDUC(NSPZ,0:NLIMPS))
         ALLOCATE(IREDUC(NSPZ,0:NLIMPS))
       
         AllOCATE (IIND(NRTAL))
         ALLOCATE (XSTOR(MSTOR1,MSTOR2))       
         ALLOCATE (XSTORV(NSTORV))

cym arrays from eirmod_clast - these are not pointers
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
         
cym test iter
         ALLOCATE (RSPLST(NPARTC,MAXLEVEL))
         ALLOCATE (ISPLST(MPARTC,MAXLEVEL))

         RSPLST=0._dp
         ISPLST=0

cym make sure the clast variables do not take exotic values      
         call eirene_init_clast

         ALLOCATE (FNUIAR(NPLS))
         FNUIAR=0._dp
                   
         NCLMT     => ISDVI(8)
         NCLMTS    => ISDVI(9)
         NWLMT     => ISDVI(10)         
         NWLMTS    => ISDVI(11)
         ICLMT     => ISDVI(12+2*NSD+2*NSDW+NCV+NRTAL :
     .        11+2*NSD+2*NSDW+NCV+2*NRTAL)
         IMETCL    => ISDVI(12+2*NSD+2*NSDW+NCV :
     .        11+2*NSD+2*NSDW+NCV+NRTAL)
         IMETWL    => ISDVI(12+2*NSD+2*NSDW+NCV+2*NRTAL :
     .        11+2*NSD+2*NSDW+NCV+2*NRTAL+NLIMPS)
         IWLMT    => ISDVI(12+2*NSD+2*NSDW+NCV+2*NRTAL+NLIMPS : MSDVI)
         
         ISPZ   => IPSTD( 9)
         NT3RD  => ICGRID( 8)
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
cpg      ISTRA  => IPSTD( 8)
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
         
         EP1    => RCGRID(1+1*N1ST : 2*N1ST)

      END IF

      END SUBROUTINE EIRENE_ALLOCATE_OPENMP

      SUBROUTINE EIRENE_DEALLOCATE_OPENMP
     
      USE EIRMOD_COMXS
      USE EIRMOD_REFLEC
      USE EIRMOD_STATIS
      USE EIRMOD_CLAST
      USE EIRMOD_CFPLK, ONLY: FNUIAR  
      !HJL ADDED FOR ROUTINE MOVE
      USE EIRMOD_COMSPL
      USE EIRMOD_COMSOU
      USE EIRMOD_COMPRT
      USE EIRMOD_CSDVI
      USE EIRMOD_CGRID
      USE EIRMOD_CUPD
      USE EIRMOD_PARMMOD

      IF(ITHREAD > 0 ) THEN
      
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
         DEALLOCATE(RCGRID)
             
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
      END IF

      END SUBROUTINE EIRENE_DEALLOCATE_OPENMP

      
      END MODULE EIRMOD_OPENMP
