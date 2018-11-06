cdr  nov. 17:  more error exits for unwritten options.
CDR  Nov. 17:  no internal boundaries considered ?
c              Interpolation at grid boundaries ?
c              tbd: Compare with routine plotting/celint.f, remove dublicated code
c
c  (was formerly:  cell-to-corner)

      subroutine eirene_fem_cell_corner (f, fcorner) 

c  interpolate cell averaged scalar tallies onto cell vertices     
c  FOR EACH CELL VERTIX USE INVERSE DISTANCE TO NEIGHBORING CELL CELL-CENTERS (com)
c  FOR WEIGHTING
cdr To be done:
c  For nltra or nltrt options the z-coordinate is a polygon or circle, resp.
c  Interpolation of vectorial quantities must account for that curved z-coordinate,
c  even in 2D (x,y)  toroidally symmetric cases.
     
      use eirmod_precision
      use eirmod_parmmod
      use eirmod_ctrig
      use eirmod_cgrid
      use eirmod_cgeom
      use eirmod_ctetra
      use eirmod_ccona
      USE eirmod_CPOLYG
      USE eirmod_CLOGAU
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      
      implicit none

      real(dp), intent(in) :: f(:)
      real(dp), intent(out) :: fcorner(:)

      real(dp), allocatable :: volsum(:)
      real(dp) :: dist, ages, xc, yc, dist1, dist2, dist3, dist4
      integer :: i, j, ir, ipart, ip, ic, in, nrk, ic1, ic2, ic3, ic4,
     .           it
      TYPE(CELL_ELEM), POINTER :: CUR

c  2d cartesian x-y- grid      
      if ((levgeo == 1) 
     .    .and. nlrad .and. nlpol.and..not.nltor) then
c  ready for 2d x-y- slab grid.  

         nrk = indpoint(nr1st,np2nd)
         allocate(volsum(nrk))
         fcorner = 0._dp
         volsum = 0._dp

         IT = 1 
         DO IR=1,NR1STM
           XC = 0.5_DP * (RSURF(IR) + RSURF(IR+1)) 
           DO IP=1,NP2NDM
             YC = 0.5_DP * (PSURF(IP) + PSURF(IP+1)) 
             IN = IR + ((IP-1)+(IT-1)*NP2T3)*NR1P2

             IC1 = INDPOINT(IR,IP)
             dist1 = 1._DP/SQRT((XC-RSURF(IR))**2+
     .                          (YC-PSURF(IP))**2) 
             FCORNER(IC1) = FCORNER(IC1) + F(IN)*DIST1
             VOLSUM(IC1) = VOLSUM(IC1) + DIST1

             IC2 = INDPOINT(IR+1,IP)
             dist2 = 1._DP/SQRT((XC-RSURF(IR+1))**2+
     .                          (YC-PSURF(IP))**2) 
             FCORNER(IC2) = FCORNER(IC2) + F(IN)*DIST2
             VOLSUM(IC2) = VOLSUM(IC2) + DIST2

             IC3 = INDPOINT(IR+1,IP+1)
             dist3 = 1._DP/SQRT((XC-RSURF(IR+1))**2+
     .                          (YC-PSURF(IP+1))**2) 
             FCORNER(IC3) = FCORNER(IC3) + F(IN)*DIST3
             VOLSUM(IC3) = VOLSUM(IC3) + DIST3

             IC4 = INDPOINT(IR,IP+1)
             dist4 = 1._DP/SQRT((XC-RSURF(IR))**2+
     .                          (YC-PSURF(IP+1))**2) 
             FCORNER(IC4) = FCORNER(IC4) + F(IN)*DIST4
             VOLSUM(IC4) = VOLSUM(IC4) + DIST4
           END DO  ! ip
         END DO  ! ir

         fcorner(1:nrk) = fcorner(1:nrk)/(volsum(1:nrk)+eps60)
         deallocate (volsum)

c  2d cartesian x-z grid
      elseif ((levgeo == 1) 
     .    .and. nlrad .and..not.nlpol.and.nltor.and.nltrz) then 
c   TO BE DONE: 2d  x,z grid
         goto 999
c  2d polar x-phi grid
      elseif ((levgeo == 1) 
     .    .and. nlrad .and..not.nlpol.and.nltor.and.nltra) then 
c   TO BE DONE: 2d  x,phi grid
         goto 999
      elseif ((levgeo == 1) 
     .    .and. nlrad .and.nlpol.and.nltor) then 
c   TO BE DONE: 3d  x,y,z  or x,y,phi grid
         goto 999


c  2d  r-theta grid, cell vertices along a coordinate line are given as polygons
      elseif (((levgeo == 2) .and. nlpol.and..not.nlcrc) .or.
     .         (levgeo == 3) .and. nlpol) then

         fcorner = 0.
         DO IR=1,NR1ST
           DO IPART=1,NPPLG
             DO IP=NPOINT(1,IPART),NPOINT(2,IPART)
               IC=INDPOINT(IR,IP)
               AGES = 0._DP
               FCORNER(IC) = 0.
               CUR => COORCELL(IC)%PCELL
               DO WHILE (ASSOCIATED(CUR))
                 IN=CUR%NOCELL
                 IF (NSTGRD(IN) == 0) THEN
                   DIST=1._DP/SQRT((XCOM(IN)-XPOL(IR,IP))**2+
     .                  (YCOM(IN)-YPOL(IR,IP))**2)
!pb                   AGES = AGES + DIST*XSTGRD(IN)
                   AGES = AGES + DIST
                   FCORNER(IC) = FCORNER(IC) + F(IN)*DIST
                 END IF
                 CUR => CUR%NEXT_CELL
               END DO
               FCORNER(IC) = FCORNER(IC) / (AGES+EPS60)
             END DO  ! ip
           END DO  ! ipart 
         END DO  ! ir

c  2d  r-theta grid, cell vertices along a coordinate line are straight lines
      elseif (((levgeo == 2) .and. nlpol.and.nlcrc)) then
c   TO BE DONE: 2d  r,theta grid, but no polygons
         goto 999
      elseif ((levgeo == 3) .and. .not. nlpol) then
c   TO BE DONE: 1d  r grid of polygons
         goto 999     

c  2d grid of triangles         
      elseif (levgeo == 4) then

         allocate(volsum(nrknot))
         volsum = eps60
         fcorner(1:nrknot) = 0.
         do i=1,ntrii
           do j=1,3
             DIST=1._DP/SQRT((XTRIAN(NECKE(J,I))-XCOM(I))**2+
     .                       (YTRIAN(NECKE(J,I))-YCOM(I))**2)  
             volsum(necke(j,i)) = volsum(necke(j,i)) + dist
             fcorner(necke(j,i)) = fcorner(necke(j,i)) + dist*f(i)
           end do
         end do

         fcorner(1:nrknot) = fcorner(1:nrknot)/volsum(1:nrknot)
         deallocate (volsum)

c  3d grid of tetrahedra
      elseif (levgeo.eq.5) then
         
         allocate(volsum(ncoord))
         volsum = eps60
         fcorner(1:ncoord) = 0.
         do i=1,ntet
           do j=1,4
             dist=1._dp/sqrt((xtetra(nteck(j,i))-xtcen(i))**2 +
     .                       (ytetra(nteck(j,i))-ytcen(i))**2 +
     .                       (ztetra(nteck(j,i))-ztcen(i))**2) 
             volsum(nteck(j,i))=volsum(nteck(j,i))+dist
             fcorner(nteck(j,i)) = fcorner(nteck(j,i)) + dist*f(i)
           enddo
         enddo
         fcorner(1:ncoor) = fcorner(1:ncoor)/volsum(1:ncoor)
         deallocate (volsum)

      else
        goto 999
      endif 

      return

999   continue
      write (iunout,*) ' levgeo = ',levgeo,' to be written in',
     .               ' subroutine cell_to_corner '
      write (iunout,*) ' calculation abandonned '
         call  eirene_exit_own(1)

      return
      end subroutine  eirene_fem_cell_corner

         
