cdr  april 2016:  looked at current default random number generator.
cdr               it seems to be a rather trivial congruential generator,
cdr               even without additive constant  (c=0.0)
cdr               very likely that this generator must be removed urgently !
cdr               maybe the original generator (nloldran) H1rn is superior by far 


      function ranf_eirene ()
 
C
C RANDOM NUMBER GENERATOR FROM
C  http://www.srcc.msu.su/num_anal/lib_na/cat/g/gsu1r.htm (in russian)

cdr:  web page is not accessible any more,  
C
C SOURCE:  Knuth, D.E. 1981, Seminumerical Algorithms, 2nd ed., vol. 2 of The Art
C          of Computer Programming (Reading, MA: Addison-Wesley)
C
cdr:  I cannot find this random generator in that reference. More likely:  
cdr:  quick and dirty home made?
cdr:  Tbd: Why did it get into EIRENE ?

C ISEED IS THE INTEGER FROM 1 TO  2147483646, AFTER FINISHING ITS VALUE IS
C (2**31) * R (N) AND CAN BE USED FOR THE FUTURE CALLS
C RETURNS ONE RANDOM NUMBER FROM 0 TO 1
C
      USE EIRMOD_PRECISION
      USE EIRMOD_CLOGAU
      implicit none
      integer :: iseed
      common /cmem/ iseed
      integer, save :: ifirst=0
      integer :: ise
      real(dp) :: ra, dummy, ranf_eirene, ranset_eirene,
     &            ranf_eirene_reinit, h1rn
 
      INTEGER D2P32M
      DOUBLE PRECISION Z,D2P31M,D2PN31,DMOD,DFLOAT
      DATA  D2PN31/4.656612873077393D-10/,  !    = 1 / 2**31   = 1/m
     .      D2P31M/2147483647.D0/,          !    = 2**31 = m
     .      D2P32M/16807/                   !    = a
 
      IF (NLOLDRAN) THEN
         ranf_eirene=h1rn(dummy)
      ELSE
         if (ifirst == 0) then
cdr  in very first call: set a fixed seed ISEED=9876543
            ise = -1
            dummy = ranset_eirene(ise)
            ifirst = 1

            call eirene_leer(2)
            call eirene_headng
     .      ('ROTTEN RANDOM NUMBER GENERATOR ACTIVATED',40)
            call eirene_headng
     .      ('IS THAT INTENTIONAL? Check flag  NLOLDRAN in RANF.f',51)
            call eirene_leer(2)

         end if
         
!pb      Z=DFLOAT(ISEED)
         Z=REAL(ISEED,KIND=DP)
         Z=DMOD(D2P32M*Z,D2P31M)  ! congruential generators, I_i+1 = a * I_i + c  (mod m),  c=0, a=16807, m=2**31
         RA=Z*D2PN31
         ISEED=Z
         
         ranf_eirene=ra
         
      END IF
      return
 
C     The following ENTRY is for reinitialization of EIRENE
 
      ENTRY ranf_eirene_reinit ()
      ifirst = 0
      ranf_eirene_reinit = 0.D0
      return
      end
