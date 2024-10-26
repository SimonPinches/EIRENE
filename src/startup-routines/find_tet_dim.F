      subroutine eirene_find_tet_dim(nr1st,ntetra,ntet,ncoord,ngitt,
     .                                  casename,ifoff)

      use eirmod_precision

      implicit none

      integer, intent(inout) :: nr1st,ntet,ntetra, ncoord,ngitt
      integer, intent(in) :: ifoff
      character(*), intent(in) :: casename
      integer :: ll, ncoor, i, id, nb1, ns1, inm1, nb2, ns2, inm2,
     .           nb3, ns3, inm3, nb4, ns4, inm4
      character(:), allocatable :: filename
      character(80) :: zeile

      LL=LEN_TRIM(CASENAME)
      allocate (character(ll+10) :: filename)

      NTET=NR1ST
      NTETRA=MAX(NTETRA,NR1ST)

      FILENAME=TRIM(CASENAME(1:LL)) // '.npco_char'
      OPEN (UNIT=30+ifoff,FILE=FILENAME,ACCESS='SEQUENTIAL',
     .      FORM='FORMATTED')
      ZEILE='*   '
      DO WHILE (ZEILE(1:1) == '*')
        READ (30+ifoff,'(A)') ZEILE
      END DO

      READ (ZEILE,*) NCOOR
      CLOSE (UNIT=30+ifoff)
      NCOORD=NCOOR

      FILENAME=TRIM(CASENAME(1:LL)) // '.elemente'
      OPEN (UNIT=30+ifoff,FILE=FILENAME,ACCESS='SEQUENTIAL',
     .      FORM='FORMATTED')

      ZEILE='*   '
      DO WHILE (ZEILE(1:1) == '*')
        READ (30+ifoff,'(A)') ZEILE
      END DO

      READ (ZEILE,*) NTET
      CLOSE (UNIT=30+ifoff)
      NTETRA=NTET+1
      NR1ST=NTETRA

      FILENAME=TRIM(CASENAME(1:LL)) // '.neighbors'
      OPEN (UNIT=30+ifoff,FILE=FILENAME,ACCESS='SEQUENTIAL',
     .      FORM='FORMATTED')

      ZEILE='*   '
      DO WHILE (ZEILE(1:1) == '*')
        READ (30+ifoff,'(A100)') ZEILE
      END DO

      NGITT=1
      DO I=1,NTET
        READ (30+ifoff,*) ID, NB1, NS1, INM1,
     .                        NB2, NS2, INM2,
     .                        NB3, NS3, INM3,
     .                        NB4, NS4, INM4
        IF (INM1 /= 0) NGITT = NGITT + 1
        IF (INM2 /= 0) NGITT = NGITT + 1
        IF (INM3 /= 0) NGITT = NGITT + 1
        IF (INM4 /= 0) NGITT = NGITT + 1
      END DO

      CLOSE (UNIT=30+ifoff)

      deallocate (filename)

      return
      end subroutine eirene_find_tet_dim
