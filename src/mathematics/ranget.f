cdr   Routine ranget is used to provide a new seed derived
cdr   in a deterministic (reproducible) way from a current random number generatur status
cdr   The current seed (status) is iseed on Commen CMEM, a new seed is returned as integer ranget
c

cdr   April 17:
cdr   Used only in case of correlated sampling, (and of MPI ??)
cdr   two random generators:
cdr   1) Random number generator H1RN
cdr   see reference of F. James 1990 review paper: F. James, CPC, 60 (1990) 329, sec 3.3
cdr   2) old IBM congruential generator (1968) (loc.cit.)


      integer function ranget_eirene(ISE)
cdr  return a legal next seed for random number generator.
cdr  output:  ise (=ranget_eirene),  return integer seed

cdr  input:  NLOLDRAN:
cdr      T :  use H1RN, which is RANMAR, F. James, CPC, 60 (1990) 329, sec 3.3
c             A legal seed must be 0<=ise<=900.000.000
cdr      F :  use old IBM (1968) generator
c             A legal seed must be 1<=ise<=2147483647 (=2**31-1)
C            ISE:   old reference seed, 
C                   from which the current status of random generator is set 
C                   and from which new seed should result in a determinsitic way
c                   

      USE EIRMOD_PRECISION
      USE EIRMOD_CLOGAU
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      implicit none
      integer, intent(in) :: ise
      integer :: iseed,idumran
      common /cmem/ iseed
      real(dp) :: ran, ranf_eirene
      integer, external :: ranset_eirene

      IF (NLOLDRAN) THEN
c  1st generator: H1RN  (RANMAR)

         if (iseed.le.0.or.iseed.gt.900000000) then
c  no legal seed availble
           write (iunout,*) 'error in fct. ranget of random generator'
           write (iunout,*) 'exit called from subr. ranget'
           call eirene_exit_own(1)
         endif

         ranget_eirene= mod(iseed+1,900000000)+1

      ELSE

c  old seed already known from previous call to ranf, or ranset
c  (from common CMEM).
         if (iseed.le.0.or.iseed.gt.2147483647.or.iseed.ne.ISE) then
c  no legal seed availble
           write (iunout,*) 'error in fct. ranget of random generator'
           write (iunout,*) 'exit called from subr. ranget'
           call eirene_exit_own(1)
         endif
c
c  set a new seed iseed 
c  Return a "derived seed" for a fresh sequence for random number starting from there
c  call ranf with seed ISEED=ISE
         ran=ranf_eirene()   !  switch to a next seed, by vasting a call to ranf().
c                               now: new ISEED on CMEM
         ranget_eirene=2147483647-ISEED
         idumran=ranset_eirene(ISE) !  return to current seed for continuation
c        write (iunout,*) 'ranget  ',ranget_eirene,ISEED,ISE

      END IF

      return
      end
