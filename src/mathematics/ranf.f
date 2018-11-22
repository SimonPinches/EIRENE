cdr  april 2016:  looked at current default random number generator.
cdr               it seems to be a rather trivial congruential generator,
cdr               even without additive constant  (c=0.0)
cdr               very likely that this generator must be removed urgently !
cdr               maybe the original generator (nloldran) H1rn is superior by far 

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
      USE EIRMOD_CLOGAU
      implicit none
      integer :: iseed
      common /cmem/ iseed
      integer, save :: ifirst=0
      real(dp) :: ra, dummy, ranf_eirene,
     .	            ranf_eirene_reinit, h1rn

cdr  parameters for SURAND  (IBM, 1968)
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
         if (ifirst == 0) then
            call eirene_leer(2)
            call eirene_headng
     .      ('ANCIENT (1968) RANDOM NUMBER GENERATOR ACTIVATED',48)
            call eirene_headng
     .      ('IS THAT INTENTIONAL? Check flag NLOLDRAN in RANF.f',50)
            call eirene_leer(2)
            ifirst=1
         end if

cdr  congruential generator, I_i+1 = a * I_i + c  (mod m),  c=0, a=16807, m=2**31-1
         Z=REAL(ISEED,KIND=DP)
         Z=DMOD(D2P32M*Z,D2P31M)

         RA=Z*D2PN31
c  save the seed for next random number.
         ISEED=Z

         ranf_eirene=ra
      END IF

c  done
      return

C     The following ENTRY is for reinitialization of EIRENE

      ENTRY ranf_eirene_reinit
cdr   indicate that random number generator is not initialized.
      ifirst = 0
      ranf_eirene_reinit = 0.D0
      return
      end
