c  oct 18: calc. gradient of a volumetric function f (input or output tally)
c          using fem interpolation options.

      subroutine eirene_calc_grad (f, fdx, fdy, fdz, lfdx, lfdy, lfdz)
c  input:
c         f(:) a function given on computational mesh 1:nrtal, 2D or 3D.
c         lfdx,lfdy,lfdz: =.t. : storage for resp. derivates is allocated
c  output:
c         fdx,fdy,fdz, corresponding derivatives df/dx (dx: in cm)
c                      df/dy, and df/dz, resp.

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
      real(dp), intent(out) :: fdx(:), fdy(:), fdz(:)
      logical, intent(in) :: lfdx, lfdy, lfdz

      integer :: i, ir, ipart, ip, in, it, ierr
      real(dp) :: xc, yc, zc, dx, dy, dz
      external :: eirene_exit_own

      interface
        subroutine eirene_df_dxyz (f,in,XC,YC,ZC,DX,DY,DZ)
          use eirmod_precision
          real(dp), intent(in) :: f(:)
          real(dp), intent(in) :: xc, yc, zc
          real(dp), intent(out) :: dx, dy, dz
          integer, intent(in) :: in
        end subroutine eirene_df_dxyz
      end interface

      ierr=0
c  1d cartesian x -grid
      if ((levgeo == 1)
     .    .and. nlrad .and. (.not.nlpol).and.(.not.nltor)) then
c  TO BE DONE: 1d x grid
         ierr=11
         goto 999
      endif

c  1d r -grid
      if ((levgeo == 2 .or. levgeo == 3)
     .    .and. nlrad .and. (.not.nlpol).and.(.not.nltor)) then
c   TO BE DONE: 1d  r grid
         ierr=12
         goto 999
      endif

c  2d cartesian x-y- grid
      if ((levgeo == 1)
     .    .and. nlrad .and. nlpol.and..not.nltor) then
c  ready for 2d x-y- slab grid.
         IT = 1
         ZC = 0._DP
         DO IR=1,NR1STM
           XC = 0.5_DP * (RSURF(IR) + RSURF(IR+1))
           DO IP=1,NP2NDM
             YC = 0.5_DP * (PSURF(IP) + PSURF(IP+1))
             IN = IR + ((IP-1)+(IT-1)*NP2T3)*NR1P2

             CALL EIRENE_DF_DXYZ (F,IN,XC,YC,ZC,DX,DY,DZ)
             IF (LFDX) FDX(IN) = DX
             IF (LFDY) FDY(IN) = DY
cdr          IF (LFDZ) FDZ(IN) = DZ  ! no z-dependence
             IF (LFDZ) FDZ(IN) = 0.0
           END DO  ! ip
         END DO  ! ir

c  2d cartesian x-z grid
      elseif ((levgeo == 1)
     .    .and. nlrad .and..not.nlpol.and.nltor.and.nltrz) then
c   TO BE DONE: 2d  x,z grid
         ierr=1
         goto 999
c  2d polar x-phi grid
      elseif ((levgeo == 1)
     .    .and. nlrad .and..not.nlpol.and.nltor.and.nltra) then
c   TO BE DONE: 2d  x,phi grid
         ierr=2
         goto 999
c  3d cartesian x-y-z grid
      elseif ((levgeo == 1)
     .    .and. nlrad .and.nlpol.and.nltor.and.nltrz) then
         DO IR=1,NR1STM
           XC = 0.5_DP * (RSURF(IR) + RSURF(IR+1))
           DO IP=1,NP2NDM
             YC = 0.5_DP * (PSURF(IP) + PSURF(IP+1))
             DO IT=1,NT3RDM
               ZC = 0.5_DP * (ZSURF(IT) + ZSURF(IT+1))
               IN = IR + ((IP-1)+(IT-1)*NP2T3)*NR1P2

               CALL EIRENE_DF_DXYZ (F,IN,XC,YC,ZC,DX,DY,DZ)
               IF (LFDX) FDX(IN) = DX
               IF (LFDY) FDY(IN) = DY
               IF (LFDZ) FDZ(IN) = DZ
             END DO  ! it
           END DO  ! ip
         END DO  ! ir

c  3d semi-toroidal grid: x-y-phi, phi approximated by polygon
      elseif ((levgeo == 1)
     .    .and. nlrad .and.nlpol.and.nltor.and.nltra) then
c  TO BE DONE: 3d x,y,phi grid
         ierr=3
         goto 999

c  2d r-theta grid, cell vertices along a coordinate line are given as polygons
      elseif (((levgeo == 2) .and. nlpol.and..not.nlcrc) .or.
     .         (levgeo == 3) .and. nlpol) then

         IT = 1
         ZC = 0._DP
         DO IR=1,NR1STM
           DO IPART=1,NPPLG
             DO IP=NPOINT(1,IPART),NPOINT(2,IPART)-1

               IN = IR + ((IP-1)+(IT-1)*NP2T3)*NR1P2
               XC = XCOM(IN)
               YC = YCOM(IN)

               CALL EIRENE_DF_DXYZ (F,IN,XC,YC,ZC,DX,DY,DZ)
               IF (LFDX) FDX(IN) = DX
               IF (LFDY) FDY(IN) = DY
               IF (LFDZ) FDZ(IN) = DZ

             END DO  ! ip
           END DO  ! ipart
         END DO  ! ir

c  2d r-theta grid, cell vertices along a coordinate line are straight lines
      elseif (((levgeo == 2) .and. nlpol.and.nlcrc)) then
c  TO BE DONE: 2d r,theta grid, but no polygons
         ierr=4
         goto 999
      elseif ((levgeo == 3) .and. .not. nlpol) then
         ierr=5
c  TO BE DONE: 1d r grid of polygons
         goto 999

c  2d grid of triangles
      elseif (levgeo == 4) then

         do i=1,ntrii
           XC = XCOM(I)
           YC = YCOM(I)
           CALL EIRENE_DF_DXYZ (F,I,XC,YC,ZC,DX,DY,DZ)
           IF (LFDX) FDX(I) = DX
           IF (LFDY) FDY(I) = DY
           IF (LFDZ) FDZ(I) = DZ
         end do

c  3d grid of tetrahedra
      elseif (levgeo.eq.5) then

        do i=1,ntet
          XC = XTCEN(I)
          YC = YTCEN(I)
          ZC = ZTCEN(I)
          CALL EIRENE_DF_DXYZ (F,I,XC,YC,ZC,DX,DY,DZ)
          IF (LFDX) FDX(I) = DX
          IF (LFDY) FDY(I) = DY
          IF (LFDZ) FDZ(I) = DZ
        enddo

      else
        ierr=6
        goto 999
      endif

      return

  999 continue
      write (iunout,*) ' levgeo = ',levgeo,' to be written in',
     .                 ' subroutine calc_grad'
      write (iunout,*) ' ierr = ',ierr
      write (iunout,*) ' calculation abandoned'
      call  eirene_exit_own(1)

      return

      end subroutine eirene_calc_grad
