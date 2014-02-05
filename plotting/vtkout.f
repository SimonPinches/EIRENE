      subroutine eirene_vtkout_head

        use eirmod_precision
        use eirmod_parmmod
        use eirmod_cadgeo
        use eirmod_clgin
        use eirmod_cplot
        use eirmod_ctext
        use eirmod_comprt, only: ivtkout
  
        implicit none
  
        open (ivtkout,file='vtk.out')
        write (ivtkout,'(a)') 
     .    '<?xml version="1.0" encoding="ISO-8859-1" ?>'
  
        write (ivtkout,'(a)')  
        write (ivtkout,'(a)') '<input>'
  
        write (ivtkout,'(a,ss,6es12.4,a)') 
     .                 '<bound>', ch3x0-ch3mx, ch3x0+ch3mx,
     .                            ch3y0-ch3my, ch3y0+ch3my,
     .                            ch3z0-ch3mz, ch3z0+ch3mz,
     .                 ' </bound>'

!       camera
        write (ivtkout,'(a)') '<camera>'
        write (ivtkout,'(a)') '<background color="1,0,1" /> '
        write (ivtkout,'(a,ss,3(es12.4,a),a)')
     .                 '<focalpoint coordinate="',
     .                  ch3x0,',',ch3y0,',',ch3z0, '"/>'
        write (ivtkout,'(a,ss,3(es12.4,a),a)')
     .                 '<position coordinate="',
     .                  ch3x0+2*ch3mx,',',ch3y0,',',ch3z0+2*ch3mz,
     .                 '"/>'
        write (ivtkout,'(a)') '<viewup coordinate="0,0,0"/>'
        write (ivtkout,'(a)') '</camera>'
      end subroutine eirene_vtkout_head
      
      
      
      subroutine eirene_vtkout_tail
        use eirmod_comprt, only: ivtkout
        implicit none
        
        write (ivtkout,'(a)') '</geometry>'
        write (ivtkout,'(a)') '</input>'
      end subroutine eirene_vtkout_tail
      
      
      
      subroutine eirene_vtkout_surfaces
        use eirmod_precision
        use eirmod_parmmod
        use eirmod_cadgeo
        use eirmod_clgin
        use eirmod_cplot
        use eirmod_ctext
        use eirmod_comprt, only: ivtkout
        
        implicit none
        
        integer :: ilim,j
        integer :: rgbcol(3,7) = reshape( 
     .                         (/ 0, 0, 0,         ! black
     .                            1, 0, 0,         ! red
     .                            0, 0, 1,         ! blue
     .                            0, 1, 0,         ! green
     .                            1, 1, 0,         ! yellow
     .                            0, 1, 1,         ! cyan
     .                            1, 0, 1 /),      ! magenta
     .                            (/3,7/))
        
!       geometry
        write (ivtkout,'(a)') '<geometry resolution="10">'
c       output all quadric surfaces (no tri/quad/pentagons!)
        write (ivtkout,'(a)') '<surfaces>'
        do ilim=1,nlimi
          if (IGJUM0(ilim).NE.0 .or. rlb(ilim) .ge. 3.) cycle
          
          if (rlb(ilim) .lt. 3.) then
            write (ivtkout,'(3a)') '<!-- ',txtsfl(ilim),' -->'
            write (ivtkout,'(a,3(i6,a),a)')
     .                 '<quadric color="',
     .                 rgbcol(1,ilcol(ilim)),',',
     .                 rgbcol(2,ilcol(ilim)),',',
     .                 rgbcol(3,ilcol(ilim)),
     .                 ' " >'
            write (ivtkout,'(a,ss,10es12.4,a)') 
     .                 '<coefficients>',
     .                 a0lm(ilim),a1lm(ilim),a2lm(ilim),a3lm(ilim),
     .                            a4lm(ilim),a5lm(ilim),a6lm(ilim),
     .                            a7lm(ilim),a8lm(ilim),a9lm(ilim),
     .                 ' </coefficients>'
            if (rlb(ilim) < 0.) then
              do j=1,ilin(ilim)
                write (ivtkout,'(a,ss,4es12.4,a)') 
     .                 '<limquadric>',
     .                   alims(j,ilim),xlims(j,ilim),
     .                   ylims(j,ilim),zlims(j,ilim),
     .                 ' </limquadric>'
              end do
              do j=1,iscn(ilim)
                write (ivtkout,'(a,ss,10es12.4,a)') 
     .                 '<limquadric>',
     .                   alims0(j,ilim),xlims1(j,ilim),ylims1(j,ilim),
     .                   zlims1(j,ilim),xlims2(j,ilim),ylims2(j,ilim),
     .                   zlims2(j,ilim),xlims3(j,ilim),ylims3(j,ilim),
     .                   zlims3(j,ilim),
     .                 ' </limquadric>'
              end do
            else if ((rlb(ilim) > 0.) .and. (rlb(ilim) < 3.)) then
              write (ivtkout,'(a,ss,6es12.4,a)') 
     .                 '<bound>', xlims1(1,ilim),xlims2(1,ilim),
     .                            ylims1(1,ilim),ylims2(1,ilim),
     .                            zlims1(1,ilim),zlims2(1,ilim),
     .                 ' </bound>'
            end if
            write (ivtkout,'(a)') '</quadric>'
          end if
        end do
        write (ivtkout,'(a)') '</surfaces>'

c       output all plane surfaces (tri/quad/pentagons)
        write (ivtkout,'(a)') '<planes>'
        do ilim=1,nlimi
c         skip invis and non tri/quad/pentagons
          if (IGJUM0(ilim).NE.0 .or. rlb(ilim) < 3.) cycle
          
          write (ivtkout,'(3a)') '<!-- ',txtsfl(ilim),' -->'
          if (rlb(ilim) < 4.) then
            write (ivtkout,'(a)') '<tri>'
            write (ivtkout,'(a,3es12.4,a)')
     .                   '<point>',P1(:,ilim),'</point>'
            write (ivtkout,'(a,3es12.4,a)')
     .                   '<point>',P2(:,ilim),'</point>'
            write (ivtkout,'(a,3es12.4,a)')
     .                   '<point>',P3(:,ilim),'</point>'
            write (ivtkout,'(a)') '</tri>'
          else if (rlb(ilim) < 5.) then
        
          else if (rlb(ilim) < 6.) then
                     
          end if
      
        end do
        write (ivtkout,'(a)') '</planes>'
      end subroutine eirene_vtkout_surfaces
      
