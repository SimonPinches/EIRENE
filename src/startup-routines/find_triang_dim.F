      subroutine eirene_find_triang_dim(nr1st,ntri,ntrii,nknot,ngitt,
     .                                  casename)

      use eirmod_precision
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_OPENFILE, ONLY: EIRENE_OPENFILE

      implicit none

      integer, intent(inout) :: nr1st,ntri,ntrii,nknot,ngitt
      character(*), intent(in) :: casename
      integer :: ll, nrknot, i, id, nb1, ns1, inm1, nb2, ns2, inm2,
     .           nb3, ns3, inm3, iunco, iunel, iunng
#if ( !defined(G95) && !defined(LEGACYCOMP) )
      character(:), allocatable :: filename
#else
      character*80 :: filename
#endif
      character(80) :: zeile

      LL=LEN_TRIM(CASENAME)
      NTRII=NR1ST
      NTRI=MAX(NTRI,NR1ST)

      if (ll > 0) then
#if ( !defined(G95) && !defined(LEGACYCOMP) )
        allocate (character(ll+10) :: filename)
#endif
        FILENAME=CASENAME(1:LL) // '.npco_char'
        IUNCO = -9999
        CALL EIRENE_OPENFILE (iunco,FILE=FILENAME,ACCESS='SEQUENTIAL',
     .      FORM='FORMATTED')

        FILENAME=CASENAME(1:LL) // '.elemente'
        IUNEL = -9999
        CALL EIRENE_OPENFILE (iunel,FILE=FILENAME,ACCESS='SEQUENTIAL',
     .       FORM='FORMATTED')

        FILENAME=CASENAME(1:LL) // '.neighbors'
        IUNNG = -9999
        CALL EIRENE_OPENFILE (IUNNG,FILE=FILENAME,ACCESS='SEQUENTIAL',
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
      WRITE(IUNOUT,*) 'NRKNOT      = ',NRKNOT
      NKNOT=NRKNOT
      CLOSE (UNIT=IUNCO)

      ZEILE='*   '
      DO WHILE (ZEILE(1:1) == '*')
        READ (IUNEL,'(A)') ZEILE
      END DO

      READ (ZEILE,*) NTRII
      NTRI=NTRII+1
      NR1ST=NTRI
      WRITE(IUNOUT,*) 'NTRII       =  ',NTRII
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
#if ( !defined(G95) && !defined(LEGACYCOMP) )
      if (allocated(filename)) deallocate (filename)
#endif
      return
      end subroutine eirene_find_triang_dim
