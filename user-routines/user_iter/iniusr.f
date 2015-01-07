 
 
      SUBROUTINE EIRENE_iniUSR
csw 22jul2011 modifications of input from b2.5.neutrals.namelist
csw already read in into extrab25 module
      use eirmod_precision
      use eirmod_parmmod
      use eirmod_clgin
      use eirmod_extrab25
      IMPLICIT NONE
      integer :: k,i,l,j
      do k=1,bn_spcsrf 
        do i=bi_spcsrf(k),bj_spcsrf(k)
          l=bl_spcsrf(i)
          if(l.lt.0) l=nlim-l
          if(bsps_absr(k).ge.0.) then 
            do j=1,nspz 
              recyct(j,l)=1.-bsps_absr(k)
            end do 
          end if 
          if(bsps_trno(k).ge.0.) then 
            do j=1,nspz
              transp(j,1,l)=bsps_trno(k)
              transp(j,2,l)=bsps_trno(k)
            enddo
          end if 
          if(bsps_trni(k).ge.0.) then 
            do j=1,nspz
              transp(j,2,l)=bsps_trni(k)
            enddo
          end if 
          if(bsps_mtri(k).gt.0.) then 
            znml(l)=bsps_mtri(k)
          end if 
          if(abs(bsps_tmpr(k)).le.1.d10) then
            ewall(l)=bsps_tmpr(k)
          end if 
          if(bsps_spph(k).ge.0.) then 
            do j=1,nspz 
              recycs(j,l)=bsps_spph(k)
            end do 
          end if 
          if(bsps_spch(k).ge.0.) then
            do j=1,nspz 
              recycc(j,l)=bsps_spch(k)
            enddo
          end if 
          if(bsps_sgrp(k).ge.0) then 
            do j=1,nspz
              isrc(j,l)=bsps_sgrp(k)
            enddo
          end if 
        end do 
      end do 
csw
      RETURN
      END
