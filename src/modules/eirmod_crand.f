      MODULE EIRMOD_CRAND
cdr  precomputed arrays of 3D random vectors, length up to 64 vectors

      USE EIRMOD_PRECISION

      IMPLICIT NONE

      PUBLIC
cdr  begin threadprivate here
      REAL(DP), DIMENSION(64), PUBLIC, SAVE ::
     R FM1, FM2, FM3,  ! Maxwellian flux, FMAXW.f, INIV1 FM1 is in direction of flux.
     R FG1, FG2, FG3,  ! std. Gaussian, Mean=0, Variance =1, FGAUSS.f, INIV2
     R FC1, FC2, FC3,  ! Cosine (Lambertian),  FCOSIN.f, INIV4, centered around (-1.,0.,0.)
     R FI1, FI2, FI3   ! Isotropic distribution, FISOTR, INIV3

      INTEGER, PUBLIC, SAVE ::
! counters for still unused random vectors. Arrays are refreshed if counter = 0
     I INIV1,  INIV2,  INIV3,  INIV4,

cdr  end threadprivate here
     I IRNDVC,  !  number of 3d vectors generated in each call
     I IRNDVH,
c
     I INTMAX, ISEEDR

      END MODULE EIRMOD_CRAND
