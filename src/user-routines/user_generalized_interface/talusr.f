c
c
      subroutine EIRENE_talusr (ICOUNT,VECTOR,TALTOT,TALAV,
     .              TXTTL,TXTSP,TXTUN,ILAST,IRET)
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD

      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      USE EIRMOD_CINIT
      USE EIRMOD_CTEXT
      USE EIRMOD_CTRIG
      USE EIRMOD_COMSOU
      USE EIRMOD_COMPRT
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      implicit NONE
      integer, intent(in) :: icount
      integer, intent(out) :: ilast, iret
      real(dp), intent(inout) :: VECTOR(*), TALTOT, TALAV

      real(dp), allocatable, save :: algv_corner(:,:), dummy(:)
      character(len=*) :: txttl,txtsp,txtun
      CHARACTER(256) :: FILENAME
      integer, parameter :: fp1=33, fp2=34
      integer :: i, ll, ip, itr

      interface
        subroutine eirene_cell_to_corner (f,fcorner)
          use eirmod_precision
          real(dp), intent(in) :: f(:)
          real(dp), intent(out) :: fcorner(:)
        end subroutine eirene_cell_to_corner
      end interface

      TXTTL=' '
      TXTSP=' '
      TXTUN=' '

      ilast=1
      iret = 0

      if (istra > 0) then

      if (istra == 1) then
        ll=len_trim(casename)
        filename=casename(1:ll) // '.neut_cell'
        open (unit=fp1+ifoff,file=filename,access='sequential',
     .        form='formatted')

        write (fp1+ifoff,'(A)') '* EIRENE NEUTRAL FILE '
        write (fp1+ifoff,'(A)') '* '
        write (fp1+ifoff,'(A)') '* No. of strata  No. of fluid species'
        write (fp1+ifoff,'(i10,5x,i10)') nstrai, nplsi


        filename=casename(1:ll) // '.neut_vert'
        open (unit=fp2+ifoff,file=filename,access='sequential',
     .        form='formatted')

        write (fp2+ifoff,'(A)') '* EIRENE NEUTRAL FILE '
        write (fp2+ifoff,'(A)') '* '
        write (fp2+ifoff,'(A)') '* No. of strata  No. of fluid species'
        write (fp2+ifoff,'(i10,5x,i10)') nstrai, nplsi

        allocate (algv_corner(nrknot,nalv))
        allocate (dummy(ntri))

        if (nalv < 2*nplsi+2) then
          write (iunout,*) ' no sources available '
          call eirene_exit_own(1)
        end if
      end if


      write (fp1+ifoff,'(A)') '* '
      write (fp1+ifoff,'(A)') '* '
      write (fp1+ifoff,'(A,i2.2)') '*** STRATUM #', istra
      write (fp1+ifoff,'(A)') '* '
      write (fp1+ifoff,'(A)') '* '

      write (fp1+ifoff,'(A)') txtsou(istra)
      write (fp1+ifoff,'(es16.7)') flux(istra)

      write (fp1+ifoff,'(A)') '* '
      write (fp1+ifoff,'(A)') '* '
      write (fp1+ifoff,'(A)') '*** PARTICLE SOURES '
      write (fp1+ifoff,'(A)') '* '
      write (fp1+ifoff,'(A)') '* '

      write (fp1+ifoff,'(i10)') ntrii
      do itr = 1, ntrii
         write (fp1+ifoff,'(i10,2es16.7)')
     .          itr, (algv(ip,itr), ip=1,nplsi)
      end do

      write (fp1+ifoff,'(A)') '* '
      write (fp1+ifoff,'(A)') '* '
      write (fp1+ifoff,'(A)') '*** MOMENTUM SOURES '
      write (fp1+ifoff,'(A)') '* '
      write (fp1+ifoff,'(A)') '* '

      write (fp1+ifoff,'(i10)') ntrii
      do itr = 1, ntrii
         write (fp1+ifoff,'(i10,2es16.7)')
     .          itr, (algv(nplsi+ip,itr), ip=1,nplsi)
      end do

      write (fp1+ifoff,'(A)') '* '
      write (fp1+ifoff,'(A)') '* '
      write (fp1+ifoff,'(A)') '*** ENERGY SOURES '
      write (fp1+ifoff,'(A)') '* '
      write (fp1+ifoff,'(A)') '* '

      write (fp1+ifoff,'(i10)') ntrii
      do itr = 1, ntrii
         write (fp1+ifoff,'(i10,2es16.7)')
     .          itr, algv(2*nplsi+1,itr), algv(2*nplsi+2,itr)
      end do

      do ip = 1, 2*nplsi+2

        dummy(1:ntrii) = algv(ip,1:ntrii)
        call eirene_cell_to_corner (dummy, algv_corner(1:nrknot,ip))

      end do

      write (fp2+ifoff,'(A)') '* '
      write (fp2+ifoff,'(A)') '* '
      write (fp2+ifoff,'(A,i2.2)') '*** STRATUM #', istra
      write (fp2+ifoff,'(A)') '* '
      write (fp2+ifoff,'(A)') '* '

      write (fp2+ifoff,'(A)') txtsou(istra)
      write (fp2+ifoff,'(es16.7)') flux(istra)

      write (fp2+ifoff,'(A)') '* '
      write (fp2+ifoff,'(A)') '* '
      write (fp2+ifoff,'(A)') '*** PARTICLE SOURES '
      write (fp2+ifoff,'(A)') '* '
      write (fp2+ifoff,'(A)') '* '

      write (fp2+ifoff,'(i10)') nrknot
      do itr = 1, nrknot
        write (fp2+ifoff,'(i10,2es16.7)') itr,
     .         (algv_corner(itr,ip), ip=1,nplsi)
      end do

      write (fp2+ifoff,'(A)') '* '
      write (fp2+ifoff,'(A)') '* '
      write (fp2+ifoff,'(A)') '*** MOMENTUM SOURES '
      write (fp2+ifoff,'(A)') '* '
      write (fp2+ifoff,'(A)') '* '

      write (fp2+ifoff,'(i10)') nrknot
      do itr = 1, nrknot
        write (fp2+ifoff,'(i10,2es16.7)') itr,
     .         (algv_corner(itr,nplsi+ip), ip=1,nplsi)
      end do

      write (fp2+ifoff,'(A)') '* '
      write (fp2+ifoff,'(A)') '* '
      write (fp2+ifoff,'(A)') '*** ENERGY SOURES '
      write (fp2+ifoff,'(A)') '* '
      write (fp2+ifoff,'(A)') '* '

      write (fp2+ifoff,'(i10)') nrknot
      do itr = 1, nrknot
        write (fp2+ifoff,'(i10,2es16.7)')
     .         itr, algv_corner(itr,2*nplsi+1),
     .              algv_corner(itr,2*nplsi+2)
      end do

      if (istra == nstrai) then
        deallocate (dummy)
        deallocate (algv_corner)
        close (fp1+ifoff)
        close (fp2+ifoff)
      end if

      end if

      VECTOR(1:nrad) = 0.
      TALTOT = 0.
      TALAV = 0.

      iret = 1
      return
      end
