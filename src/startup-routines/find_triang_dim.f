      subroutine eirene_find_triang_dim(nr1st,ntri,ntrii,nknot,ngitt,
     .                                  casename,ifoff)

      use eirmod_precision
      USE EIRMOD_COMPRT, ONLY: IUNOUT

      implicit none

      integer, intent(inout) :: nr1st,ntri,ntrii,nknot,ngitt
      integer, intent(in) :: ifoff
      character(*), intent(in) :: casename 
      integer :: ll, nrknot, i, id, nb1, ns1, inm1, nb2, ns2, inm2,
     .           nb3, ns3, inm3, iunco, iunel, iunng
      character(:), allocatable :: filename
      character(80) :: zeile
      
      LL=LEN_TRIM(CASENAME)
      NTRII=NR1ST
      NTRI=MAX(NTRI,NR1ST)

      if (ll > 0) then
        allocate (character(ll+10) :: filename)

        FILENAME=CASENAME(1:LL) // '.npco_char'
        OPEN (NEWUNIT=iunco,FILE=FILENAME,ACCESS='SEQUENTIAL',
     .      FORM='FORMATTED')
      
        FILENAME=CASENAME(1:LL) // '.elemente'
        OPEN (NEWUNIT=iunel,FILE=FILENAME,ACCESS='SEQUENTIAL',
     .       FORM='FORMATTED')
      
        FILENAME=CASENAME(1:LL) // '.neighbors'
        OPEN (NEWUNIT=IUNNG,FILE=FILENAME,ACCESS='SEQUENTIAL',
     .      FORM='FORMATTED')
      else
        iunco = 33
        iunel = 34
        iunng = 35
      end if
     

      ZEILE='*   '
      DO WHILE (ZEILE(1:1) == '*')
        READ (IUNCO,'(A)') ZEILE
      END DO
      
      READ (ZEILE,*) NRKNOT
      WRITE(IUNOUT,*) 'NRKNOT = ',NRKNOT
      NKNOT=NRKNOT
      CLOSE (UNIT=IUNCO)
      
      ZEILE='*   '
      DO WHILE (ZEILE(1:1) == '*')
        READ (IUNEL,'(A)') ZEILE
      END DO
      
      READ (ZEILE,*) NTRII
      NTRI=NTRII+1
      NR1ST=NTRI
      WRITE(IUNOUT,*) 'NTRII =  ',NTRII
      CLOSE (UNIT=IUNEL)

      ZEILE='*   '
      DO WHILE (ZEILE(1:1) == '*')
        READ (IUNNG,'(A)') ZEILE
      END DO

      NGITT=1
      DO I=1,NTRII
        READ (IUNNG,*) ID, NB1, NS1, INM1,
     .                     NB2, NS2, INM2,
     .                     NB3, NS3, INM3
        IF (INM1 /= 0) NGITT = NGITT + 1
        IF (INM2 /= 0) NGITT = NGITT + 1
        IF (INM3 /= 0) NGITT = NGITT + 1
      END DO

      CLOSE (UNIT=IUNNG)

      if (allocated(filename)) deallocate (filename)

      end subroutine eirene_find_triang_dim
