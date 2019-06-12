      subroutine eirene_outidlconf

      use eirmod_parmmod
      USE EIRMOD_CTEXT
      USE EIRMOD_CGRID
      USE EIRMOD_CPOLYG
      USE EIRMOD_CTRIG
      USE EIRMOD_CTETRA
      USE EIRMOD_CINIT

      implicit none
      integer :: iout,I1,I2,I3
      CHARACTER(10) :: CDATE, CTIME
      CHARACTER(LEN(CASENAME)) :: CNAME


      IOUT = 3 + IFOFF

      open (unit=iout,file='idl.conf',form='FORMATTED',
     .            ACCESS='SEQUENTIAL')

      write (iout,'(A)') ' PRODUCING PROGRAM: EIRENE'

      write (iout,'(//A,I6)') ' LEVGEO ',LEVGEO

      select case (levgeo)
      case (1,2)
!  slab or circ geometry
        write (iout,'(A,I10)') ' NR1ST ',NR1ST
        write (iout,'(A,I10)') ' NP2ND ',NP2ND
        write (iout,'(A,I10)') ' NT3RD ',NT3RD
      case (3)
!  quadrangle mesh (B2)
        write (iout,'(A,I10)') ' NR1ST ',NR1ST
        write (iout,'(A,I10)') ' NP2ND ',NP2ND
        write (iout,'(A,I10)') ' NT3RD ',NT3RD
        write (iout,'(A,I10)') ' NRPLG ',NRPLG
      case (4)
!  triangle mesh
        CNAME=CASENAME
        IF (LEN_TRIM(CNAME) == 0) CNAME='triang'
        write (iout,'(A,A)')  ' CASENAME ',CNAME
        write (iout,'(A,I10)') ' NRKNOT ',NRKNOT
        write (iout,'(A,I10)') ' NTRII  ',NTRII
      case (5)
!  tetrahedron mesh
        write (iout,'(A,A)')  ' CASENAME ',CASENAME
        write (iout,'(A,I10)') ' NCOORD ',NCOORD
        write (iout,'(A,I10)') ' NTET   ',NTET
      case default

      end select

      write (iout,'(//A)') txtrun

      CALL DATE_AND_TIME(CDATE,CTIME)
      READ(CDATE(1:4),*) I1
      READ(CDATE(5:6),*) I2
      READ(CDATE(7:8),*) I3
      WRITE (iout,'(//1X,A6,1X,2(I2,1X),I4)') 'DATE: ',I3,I2,I1
      READ(CTIME(1:2),*) I1
      READ(CTIME(3:4),*) I2
      READ(CTIME(5:6),*) I3
      WRITE (iout,'(1X,A6,1X,3(I2,1X))') 'TIME: ',I1,I2,I3

      close (unit=iout)

      return

      end subroutine eirene_outidlconf
