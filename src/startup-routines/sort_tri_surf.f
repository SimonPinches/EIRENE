      subroutine eirene_sort_tri_surf(orig_srf, srf)

      use eirmod_precision
      use eirmod_parmmod
      use eirmod_ctrig

      implicit none

      type(tri_surf),intent(in) :: orig_srf
      type(tri_surf), intent(inout) :: srf
      integer :: i, nt, it, is, is1, i1, i2, nsingle, nmulti, inew,
     .           ipoint
      integer, allocatable :: ivert(:,:)
      integer, allocatable, save :: nodes(:)
      logical, allocatable :: visit(:)

      if (.not.allocated(nodes)) allocate(nodes(nknots))
      nodes = 0

      nt = orig_srf%numtr

      if (.not.associated(srf%itrias)) allocate(srf%itrias(nt))
      if (.not.associated(srf%itrisi)) allocate(srf%itrisi(nt))
      if (.not.associated(srf%bglt)) allocate(srf%bglt(nt+1))
      srf%numtr = nt

      allocate (ivert(nt,2))
      allocate (visit(nt))
      ivert = 0
      visit = .true.

      do i = 1,nt
         it = orig_srf%itrias(i)
         is = orig_srf%itrisi(i)
         is1 = mod(is,3) + 1
         i1 = necke(is,it)
         i2 = necke(is1,it)
!  find vertices of triangle face
         ivert(i,1) = i1
         ivert(i,2) = i2
         nodes(i1) = nodes(i1) + 1
         nodes(i2) = nodes(i2) + 1
      end do

!  number of endpoints, indicating nonclosed contour
      nsingle = count(nodes == 1)
!  number of x-points, resulting in multiple parts of surface
      nmulti = count(nodes > 2)
      srf%lsplit = (nsingle > 2) .or. (nmulti > 0)

      inew = 0
      srf%bglt(1) = 0._dp

 100  continue

!  identify starting point
      if (nsingle > 0) then
!  find endpoint of open contour
         do i = 1, nknots
            if (nodes(i) == 1) then
               ipoint = i
               nsingle = nsingle -1
               exit
            end if
         end do
      else
!  use the first point on a closed contour
         do i = 1, nt
            if (visit(i)) then
               ipoint = ivert(i,1)
               exit
            end if
         end do
      end if

!     identify triangle face containing vertex ipoint
      outer_loop: do
        do i = 1, nt
           if (.not.visit(i)) cycle ! triangle face was already used
           if ((ivert(i,1) == ipoint) .or. (ivert(i,2) == ipoint)) then
              inew = inew + 1
              it = orig_srf%itrias(i)
              is = orig_srf%itrisi(i)
              srf%itrias(inew) = it
              srf%itrisi(inew) = is
              srf%bglt(inew+1) = srf%bglt(inew) +
     .           sqrt(vtrix(is,it)**2 + vtriy(is,it)**2)
              visit(i) = .false.
              i1 = ivert(i,1)
              i2 = ivert(i,2)
              nodes(i1) = nodes(i1) - 1
              nodes(i2) = nodes(i2) - 1
              ipoint = ivert(i,1) + ivert(i,2) - ipoint
              if (nodes(ipoint) == 0) then
! end of open countour reached
                 nsingle = nsingle - 1
                 exit outer_loop
              end if
              cycle outer_loop
           end if
        end do

      end do outer_loop

      if (any(visit)) goto 100

      deallocate(ivert)
      deallocate(visit)

      return
      end subroutine eirene_sort_tri_surf
