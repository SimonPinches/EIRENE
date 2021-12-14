      MODULE EIRMOD_RANF
      USE EIRMOD_PRECISION
      USE EIRMOD_CLOGAU
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_H1RNM
      IMPLICIT NONE
      PRIVATE

      PUBLIC :: ranf_eirene, ranf_eirene_reinit,
     .          ranset_eirene, ranset_eirene_reinit,
     .          ranget_eirene 

      integer, save :: ifirst_ranf=0, ifirst_ranset=0
cccccccccccccccccccccccccc
cym common /cmem/ replaced
cccccccccccccccccccccccccc
      integer :: iseed

!$omp threadprivate(ifirst_ranf,ifirst_ranset,iseed)
 
      CONTAINS

cdr April 2016:  looked at current default random number generator.
cdr               it seems to be a rather trivial congruential generator,
cdr               even without additive constant  (c=0.0)
cdr               very likely that this generator must be removed urgently !
cdr              maybe the CERN generator (nlranmar) H1rn is superior by far

cdr April 2017:  references found, see F. James,
c       ref.: review paper  F. James, CPC, 60 (1990) 329,  for both generators
c
cdr              H1RN is a variant of the RANMAR generator described there,
cdr              while in case NLOLDRAN=F the generator SURAND, IBM 1968, is activated


      function ranf_eirene ()

cdr  if nloldran = F:
C
C RANDOM NUMBER GENERATOR FROM
C  http://www.srcc.msu.su/num_anal/lib_na/cat/g/gsu1r.htm (in russian)

cdr:  web page is not accessible anymore,
C
C SOURCE:  Knuth, D.E. 1981, Seminumerical Algorithms, 2nd ed., vol. 2 of The Art
C          of Computer Programming (Reading, MA: Addison-Wesley)
C
cdr:  I cannot find this random generator in that reference.
c
cdr   More likely:
cdr:  quick and dirty home made?
cdr:  Tbd: Why did it get into EIRENE ?
cdr   April 17: Found it elsewhere: it is the SURAND generator,
cdr                                 for IBM for its system/360
cdr                                 It was state of the art in 1968.
cdr                                 Not recommended anymore.

C ISEED IS THE INTEGER FROM 1 TO 2147483647, AFTER FINISHING ITS VALUE IS
C (2**31) * R (N) AND CAN BE USED FOR THE FUTURE CALLS
C RETURNS ONE RANDOM NUMBER FROM 0 TO 1

cdr  if nloldran = T:  H1RN generator is used
C
      USE EIRMOD_PRECISION
      implicit none
cym      integer :: iseed
cym      common /cmem/ iseed
      real(dp) :: ra, dummy, ranf_eirene !, h1rn

cdr  parameters for SURAND  (IBM, 1968), aka: "GGL generator". 
cdr  generally quoted as: "should not be used anymore"
cdr  since 19-eighties
      INTEGER D2P32M
      DOUBLE PRECISION Z,D2P31M,D2PN31
      DATA  D2PN31/4.656612873077393D-10/,  !    = 1 / 2**31   = 1/m
     .      D2P31M/2147483647.D0/,          !    = 2**31-1 = m
     .      D2P32M/16807/                   !    = a = 7**5

      IF (NLOLDRAN) THEN
cdr April 2017:  h1rn is a variant of RANMAR. It has period 2**144, if properly used.
cdr              initialization is by a 32 bit integer
cdr here: initialization enforced with seeds 0<= iseed<=900.000.000 in subr. ranset.
cdr       for each such seed a different sequence of average length 10**30 is produced.
         ranf_eirene=h1rn(dummy)

      ELSE

cdr  IBM, 1968 generator, kept here for historic reasons and backward compatibility.
         if (ifirst_ranf == 0) then
!$OMP MASTER
            call eirene_leer(2)
            call eirene_headng
     .      ('ANCIENT (1968) RANDOM NUMBER GENERATOR ACTIVATED',48)
            call eirene_headng
     .      ('IS THAT INTENTIONAL? Check flag NLOLDRAN in RANF.f',50)
            call eirene_leer(2)
!$OMP END MASTER            
            ifirst_ranf=1
         end if

cdr  congruential generator, I_i+1 = a * I_i + c  (mod m),  c=0, a=16807, m=2**31-1
         Z=REAL(ISEED,KIND=DP)
         Z=DMOD(D2P32M*Z,D2P31M)

         RA=Z*D2PN31
c  save the seed for next random number.
         ISEED=INT(Z)

         ranf_eirene=ra
      END IF

c  done
      return
      end function ranf_eirene

C     The following FUNCTION is for reinitialization of EIRENE

      FUNCTION ranf_eirene_reinit()
      real(dp) :: ranf_eirene_reinit
cdr   indicate that random number generator is not initialized.
      ifirst_ranf = 0
      ranf_eirene_reinit = 0.D0
      return
      end FUNCTION ranf_eirene_reinit

cdr   initialize random number generator, set the random number seed "iseed",
cdr   store that in Common CMEM, and initialize random generator.
c
      INTEGER FUNCTION RANSET_EIRENE(ISE)
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

      implicit none
      integer, intent(in) :: ise
cdr   integer :: iseed1, iseed2, ! older version of h1rn with two seeds
cym      integer :: iseed
cdr  status of RANMAR (H1RN) generator is stored in Common RASET (after call to H1RNIN)
cym      common /cmem/ iseed

      IF (NLOLDRAN) THEN
cdr there are various variants of initializer RMARIN, taking either
cdr one, two or 4 input seeds.
cdr currently we take one seed, 0<=ise<=900.000.000,
cdr then produce two smaller integer seeds iseed1,iseed2 from that,
cdr in SUBR. H1RNIN, and initialize the status of H1RN with the contents of
cdr common RASET1 (full 102 words)

         if (ise <= 0) then
           if (ifirst_ranset == 0) then
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
           if (ifirst_ranset == 0) then
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

      ifirst_ranset = 1
cdr  just return iseed, the current legal seed used, on which the next random
cdr  number returned from ranf() will be based. This seed is also saved in common CMEM
      ranset_eirene=iseed

      return
      end function ranset_eirene

C     The following FUNCTION is for reinitialization of EIRENE

      FUNCTION ranset_eirene_reinit()
      implicit none
      integer :: ranset_eirene_reinit
      ifirst_ranset = 0
      ranset_eirene_reinit = 0
      return
      end FUNCTION ranset_eirene_reinit

cdr   Routine ranget is used to provide a new seed derived
cdr   in a deterministic (reproducible) way from a current random number generator status
cdr   The current seed (status) is
cdr   ISEED in Common CMEM for the SURAND generator,
cdr   and is coded in Commons RASET1, RASET2 for the RANMAR (= H1RN) generator
cdr   as array of length 100.
cdr   A new seed is returned as integer RANGET_EIRENE
c

cdr   April 17:
cdr   Used only in case of correlated sampling, (and for MPI, OPENMP parallelization)
cdr   Two random generators:
cdr   1) Random number generator H1RN
cdr   see reference of F. James 1990 review paper: F. James, CPC, 60 (1990) 329, sec 3.3
cdr   2) pre-historic IBM congruential generator SURAND (1968) (loc.cit.)


      INTEGER FUNCTION RANGET_EIRENE(ISE)
cdr  return a legal next seed for random number generator.
cdr  output:  ISE (=ranget_eirene),  return a legal integer seed

cdr  input:  NLOLDRAN:
cdr      T :  use H1RN, which is RANMAR, F. James, CPC, 60 (1990) 329, sec 3.3
c           A legal seed must be 0<=ISE<=900.000.000
cdr     F : use pre-historic IBM (1968) generator
c           A legal seed must be 1<=ISE<=2147483647 (=2**31-1)
C            ISE:   old reference seed,
C                   from which the current status of random generator is set
C                   and from which new seed should result in a deterministic way
c

      implicit none
      integer, intent(in) :: ise
      integer :: idumran
cym cccccccccccccccccccccccccccccccccc
cym      integer :: iseed,idumran
cym      common /cmem/ iseed
cym cccccccccccccccccccccccccccccccccc
      real(dp) :: ran

      IF (NLOLDRAN) THEN
c  1st generator: H1RN  (RANMAR)

         if (iseed.le.0.or.iseed.gt.900000000) then
c  no legal seed available
           write (iunout,*) 'error in fct. ranget of random generator'
           write (iunout,*) 'exit called from subr. ranget'
           call eirene_exit_own(1)
         endif

         ranget_eirene= mod(iseed+1,900000000)+1

      ELSE

c  old seed already known from previous call to ranf, or ranset
c  (from common CMEM).
         if (iseed.le.0.or.iseed.gt.2147483647.or.iseed.ne.ISE) then
c  no legal seed available
           write (iunout,*) 'error in fct. ranget of random generator'
           write (iunout,*) 'exit called from subr. ranget'
           call eirene_exit_own(1)
         endif
c
c  set a new seed ISEED
c  Return a "derived seed" for a fresh sequence for random number starting from there
c  call ranf with seed ISEED=ISE
         ran=ranf_eirene()   !  switch to a next seed, by wasting a call to ranf().
c                               now: new ISEED on CMEM
         ranget_eirene=2147483647-ISEED
         idumran=ranset_eirene(ISE) !  return to current seed for continuation
c        write (iunout,*) 'ranget  ',ranget_eirene,ISEED,ISE

      END IF

      return
      end function ranget_eirene

C     The following ENTRY is for reinitialization of EIRENE

      FUNCTION ranset_eirene_reinit()
      implicit none
      integer :: ranset_eirene_reinit
      ifirst_ranset = 0
      ranset_eirene_reinit = 0
      return
      end FUNCTION ranset_eirene_reinit

C     The following ENTRY is for reinitialization of EIRENE

      FUNCTION ranf_eirene_reinit()
      real(dp) :: ranf_eirene_reinit
cdr   indicate that random number generator is not initialized.
      ifirst_ranf = 0
      ranf_eirene_reinit = 0.D0
      return
      end FUNCTION ranf_eirene_reinit

      
      END MODULE EIRMOD_RANF
