c
c
      subroutine EIRENE_talusr (ICOUNT,VECTOR,TALTOT,TALAV,
     .              TXTTL,TXTSP,TXTUN,ILAST,*)
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      USE EIRMOD_EIRBRA
      USE EIRMOD_COMPRT
      USE EIRMOD_COMUSR
      USE EIRMOD_CCOUPL
      USE EIRMOD_CPOLYG
      USE EIRMOD_CESTIM
      USE EIRMOD_CSDVI
      implicit NONE
      integer, intent(in) :: icount
      integer, intent(out) :: ilast
      real(dp), intent(in) :: vector(*), TALTOT, TALAV
      real(dp) :: dumout(0:ndxp,0:ndyp), dummy(0:ndxp,0:ndyp)
      real(dp) :: sigout(0:ndxp,0:ndyp)
      character(len=*) :: txttl,txtsp,txtun
      integer :: i, ix, iy, icp, icp2, icp3, in
 
      ilast = 1

      if ( .true. ) return 1

! pb quick and dirty for quantities on quadrangular mesh only
      if ((levgeo.ne.3).and.(levgeo.ne.4)) return
      if (istra == 0) return

      open (unit=57,file='copv.out',position='APPEND')

      icp=nplsi
      icp2=2*nplsi
      icp3=3*nplsi
    
      write (57,*) ' ISTRA ', istra
      
      write (57,'(A)') 'SNI(COPV)'

      do ipls = 1,nfla
        write (57,'(a,i6)') ' IPLS = ',ipls
        dumout = 0._dp
        sigout = 0._dp
        do ix = 1, ndxa
          do iy=1,ndya
            in=IY+(IX-1)*NR1TAL
            dumout(ix,iy) = copv(icp+ipls,in)
            sigout(ix,iy) = sigma(ipls,in)
          end do
        end do
        CALL EIRENE_INDMPI (dumout,DUMMY,NDX,NDY,1  ,NDXA,NDYA,1   ,
     .               NCUTB,NCUTL,NPOINT,NPPLG,1,1)
        do ix = 1, ndxa
          write (57,'(10es15.7)') (dumout(ix,iy),iy=1, ndya)
        end do
        CALL EIRENE_INDMPI (sigout,DUMMY,NDX,NDY,1  ,NDXA,NDYA,1   ,
     .               NCUTB,NCUTL,NPOINT,NPPLG,1,1)
        do ix = 1, ndxa
          write (57,'(10es15.7)') (sigout(ix,iy),iy=1, ndya)
        end do
      end do
      
      write (57,'(/1x,A)') 'SMO(COPV)'

      do ipls = 1,nfla
        write (57,'(a,i6)') ' IPLS = ',ipls
        dumout = 0._dp
        sigout = 0._dp
        do ix = 1, ndxa
          do iy=1,ndya
            in=IY+(IX-1)*NR1TAL
            dumout(ix,iy) = copv(icp2+ipls,in)
            sigout(ix,iy) = sigma(icp+ipls,in)
          end do
        end do
        CALL EIRENE_INDMPI (dumout,DUMMY,NDX,NDY,1  ,NDXA,NDYA,1   ,
     .               NCUTB,NCUTL,NPOINT,NPPLG,1,1)
        do ix = 1, ndxa
          write (57,'(10es15.7)') (dumout(ix,iy),iy=1, ndya)
        end do
        CALL EIRENE_INDMPI (sigout,DUMMY,NDX,NDY,1  ,NDXA,NDYA,1   ,
     .               NCUTB,NCUTL,NPOINT,NPPLG,1,1)
        do ix = 1, ndxa
          write (57,'(10es15.7)') (sigout(ix,iy),iy=1, ndya)
        end do
      end do
      
      write (57,'(/1x,A)') 'SEE(COPV)'

      dumout = 0._dp
      sigout = 0._dp
      do ix = 1, ndxa
        do iy=1,ndya
          in=IY+(IX-1)*NR1TAL
          dumout(ix,iy) = copv(icp3+1,in)
          sigout(ix,iy) = sigma(icp2+1,in)
        end do
      end do
      CALL EIRENE_INDMPI (dumout,DUMMY,NDX,NDY,1  ,NDXA,NDYA,1   ,
     .               NCUTB,NCUTL,NPOINT,NPPLG,1,1)
      do ix = 1, ndxa
        write (57,'(10es15.7)') (dumout(ix,iy),iy=1, ndya)
      end do
      CALL EIRENE_INDMPI (sigout,DUMMY,NDX,NDY,1  ,NDXA,NDYA,1   ,
     .               NCUTB,NCUTL,NPOINT,NPPLG,1,1)
      do ix = 1, ndxa
        write (57,'(10es15.7)') (sigout(ix,iy),iy=1, ndya)
      end do
      
      write (57,'(/1x,A)') 'SEI(COPV)'

      dumout = 0._dp
      sigout = 0._dp
      do ix = 1, ndxa
        do iy=1,ndya
          in=IY+(IX-1)*NR1TAL
          dumout(ix,iy) = copv(icp3+2,in)
          sigout(ix,iy) = sigma(icp2+2,in)
        end do
      end do
      CALL EIRENE_INDMPI (dumout,DUMMY,NDX,NDY,1  ,NDXA,NDYA,1   ,
     .               NCUTB,NCUTL,NPOINT,NPPLG,1,1)
      do ix = 1, ndxa
        write (57,'(10es15.7)') (dumout(ix,iy),iy=1, ndya)
      end do
      CALL EIRENE_INDMPI (sigout,DUMMY,NDX,NDY,1  ,NDXA,NDYA,1   ,
     .               NCUTB,NCUTL,NPOINT,NPPLG,1,1)
      do ix = 1, ndxa
        write (57,'(10es15.7)') (sigout(ix,iy),iy=1, ndya)
      end do

      close (unit=57)



      open (unit=56,file='sources.out',position='APPEND')
      
      write (56,*) ' ISTRA ', istra
      
      write (56,'(A)') 'SNI'
      do ipls = 1,nfla
        write (56,'(a,i6)') ' IPLS = ',ipls
        do ix = 1, ndxa
          write (56,'(10es15.7)') (sni(ix,iy,ipls,istra),iy=1, ndya)
        end do
      end do
      
      write (56,'(/1x,A)') 'SMO'
      do ipls = 1,nfla
        write (56,'(a,i6)') ' IPLS = ',ipls
        do ix = 1, ndxa
          write (56,'(10es15.7)') (smo(ix,iy,ipls,istra),iy=1, ndya)
        end do
      end do
      
      write (56,'(/1x,A)') 'SEE'
      do ix = 1, ndxa
        write (56,'(10es15.7)') (see(ix,iy,istra),iy=1, ndya)
      end do
      
      write (56,'(/1x,A)') 'SEI'
      do ix = 1, ndxa
        write (56,'(10es15.7)') (sei(ix,iy,istra),iy=1, ndya)
      end do
      
      close (unit=56)
      ilast=1
      return 1
      end
