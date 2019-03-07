      subroutine eirene_outidlsrf
c
c     this subroutine prints information about surface properties for IDL tool
c
      USE EIRMOD_PARMMOD
      USE EIRMOD_CLGIN
      USE EIRMOD_CGRID
      USE EIRMOD_CADGEO

      IMPLICIT NONE

      integer :: iout, j, iad_cell, icos, ists

      OPEN (NEWUNIT=IOUT,FILE='srf_properties',FORM='FORMATTED',
     .      ACCESS='SEQUENTIAL')

      do ists = 1, nstsi
        j = nlim + ists
        iad_cell = 0
        if ((iliin(j) < 0) .and. (ilswch(j) /= 0)) then
          if (iswich(5,j) /= 0) then
            iad_cell = nsurf + ilacll(j)
          else if (iswich(6,j) /= 0) then
            icos = 1
            iad_cell = nsurf + ILBLCK(j)+ICOS*ISWICH(6,j)*ILACLL(j)
          end if
        end if
        if (levgeo == 5) then
          write (iout,'(3i9)') ists, iliin(j), iad_cell
        else
          write (iout,'(3i9)') -ists, iliin(j), iad_cell
        end if
      end do

      if (levgeo /= 5) then
        do j = 1, nlimi
          iad_cell = 0
          if ((iliin(j) < 0) .and. (ilswch(j) /= 0)) then
            if (iswich(5,j) /= 0) then
              iad_cell = ilacll(j)
            else if (iswich(6,j) /= 0) then
              icos = 1
              iad_cell = ILBLCK(j)+ICOS*ISWICH(6,j)*ILACLL(j)
            end if
          end if
          write (iout,'(3i9)') j, iliin(j), iad_cell
        end do
      end if

      close (iout)

      end  subroutine eirene_outidlsrf
