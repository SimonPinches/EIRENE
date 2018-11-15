cdr   initialize random number generator, set the random number seed "iseed",
cdr   store that in Common CMEM, and initialize random generator.
c
      integer function ranset_eirene(ise)
cdr   input:
c       nloldran:  (CLOGAU)
c         T: use random number generator  H1RN, H1RNIN,...., recommended.
c         F: use random number generator (IBM 1968),     not recommended.
c       ise: input seed, must be an 2**32 bit integer, i.e. 0<=ise<=2**31-1
c            additionally: 0<= ise <= 900.000.000 for H1RN, H1RNIN generator

c       ref.: review paper  F. James, CPC, 60 (1990) 329,  for both generators
c
c  if in very first call ise <=0:  try a default seed.
c  if that happens in a later call:  error exit.

      USE EIRMOD_PRECISION
      USE EIRMOD_CLOGAU
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      implicit none
      integer, intent(in) :: ise
cdr   integer :: iseed1, iseed2, ! older version of h1rn with two seeds
      integer :: iseed
cdr  status of RANMAR (H1RN) generator is stored in Common RASET (after call to H1RNIN)
      common /cmem/ iseed
      integer :: ranset_eirene_reinit
      integer, save :: ifirst=0

      IF (NLOLDRAN) THEN
cdr there are various variants of initializer RMARIN, taking either
cdr one, two or 4 input seeds.
cdr currently we take one seed, 0<=ise<=900.000.000,
cdr then produce two smaller integer seeds iseed1,iseed2 from that,
cdr in SUBR. H1RNIN, and initialize the status of H1RN with the content of 
cdr Common RASET1 (full 102 words)

         if (ise <= 0) then
           if (ifirst == 0) then
cdr this default seed produces the 4 small Marsaglia-Zaman seeds. 
c   Loc.cit. F.James, CPC (1990)
             iseed = 54217137
             call h1rnin(iseed)
!   else  seed is already set, do nothing
           else
             write (iunout,*) 'RANSET.F: '
             write (iunout,*) 'no initialization of random generator'
             write (iunout,*) 'continue sequence without new seed'
           end if
         else  ! here: ise > 0
cdr  900000000 is the maximum value for the single input seed initialization. Loc.cit.
           iseed=mod(ise-1,900000000)+1
           call h1rnin(iseed)
         end if


cdr  initializing the old IBM (1968) generator
      ELSE

         if (ise <= 0) then
           if (ifirst == 0) then
             iseed = 9876543
!   else seed is already set, do nothing, some older seed is already available
           else
             write (iunout,*) 'RANSET.F: '
             write (iunout,*) 'no initialization of random generator'
             write (iunout,*) 'continue sequence without new seed'
           end if
         else  ! here: ise > 0
! constrain input seed to be 0 < iseed <= 2**31-1
           iseed=mod(ise-1,2147483647)+1
c          write (iunout,*) 'ranset: seed set to ',iseed
         end if

      END IF

      ifirst = 1
cdr  just return iseed, the current legal seed used, on which the next random
cdr  number returned from ranf() will be based. This seed is also saved in Common CMEM
      ranset_eirene=iseed

      return

C     The following ENTRY is for reinitialization of EIRENE

      ENTRY ranset_eirene_reinit
      ifirst = 0
      ranset_eirene_reinit = 0
      return
      end
