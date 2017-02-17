
subroutine eirene_xs_init_driver(np,nh,nz,lread)
  use mod_eirene_xs
  implicit none
  integer, intent(in) :: np,nh,nz
  logical, intent(out) :: lread
  integer :: ierr
  
  if(.not.linit) then
     call eirene_xs_linkdb(ierr)
     if(ierr /= 0) then
        ! database/eirene.input cannot be found
        lread=.false.
        return
     endif

     call eirene_xs_init
  endif
  if(.not.lreadspecies) then
     call eirene_xs_read_species(nh,nz)
  endif

  lread = lreadspecies.and.linit
end subroutine eirene_xs_init_driver

subroutine eirene_xs_dealloc_driver
  use mod_eirene_xs
  implicit none
  call eirene_xs_dealloc_species
  call eirene_xs_kill
end subroutine eirene_xs_dealloc_driver

subroutine GetEireneXS(k,ihh,aneut,eneut,vneut,aion,vion,driftc,svi,dsvi,scx,dscx,smi,dsmi,smd,dsmd,svr,dsvr,de,te,ti,ierr)
  use mod_eirene_xs
  implicit none
  integer,intent(in):: k,ihh      
  real*8, intent(in) :: aneut,eneut,vneut(3),aion,vion,driftc(3),de,ti,te
  real*8,intent(out):: svi,dsvi,scx,dscx,smi,dsmi,smd,dsmd,svr,dsvr
  integer, intent(out) :: ierr
  real*8 :: rsvi,rdsvi,rscx,rdscx,rsmi,rdsmi,rsmd,rdsmd,rsvr,rdsvr
  real*8 :: vi(3),v0(3),di,tep,tem,tip,tim,vneu
  real*8 :: svi1,svi2,smi1,smi2,smd1,smd2,scx1,scx2,svr1,svr2
  integer :: iswrr,ih
  real*8, parameter :: dlinrat=0.01  

  ! prepare data
  di=de
  vneu = sqrt(2.d0*eneut/aneut/1.04394d-12) 
  v0(1:3) = vneu*vneut(1:3)
  vi(1:3) = vion*driftc(1:3)
  tep = te*(1.d0+dlinrat)
  tem = te*(1.d0-dlinrat)
  tip = ti*(1.d0+dlinrat)
  tim = ti*(1.d0-dlinrat)
  ih=ihh

  if(ihh < 0) then
     ih = -izz2iz(-ihh)
  endif


  ! electron impact + te-linearisation
  iswrr = 1
  call eirene_xs_getrc(ih,iswrr, te,de,ti,di, v0,vi)
  rsvi = rca(1)
  rsmi = rcm(2)
  rsmd = rcm(1)
  call eirene_xs_getrc(ih,iswrr, tep,de,tip,di, v0,vi)
  svi1 = rca(1)
  smi1 = rcm(2)
  smd1 = rcm(1)
  call eirene_xs_getrc(ih,iswrr, tem,de,tim,di, v0,vi)
  svi2 = rca(1)
  smi2 = rcm(2)
  smd2 = rcm(1)
  rdsvi= (svi1-svi2)/2./te/dlinrat
  rdsmi= (smi1-smi2)/2./te/dlinrat
  rdsmd= (smd1-smd2)/2./te/dlinrat


  ! charge exchange + ti-linearisation
  iswrr=3
  call eirene_xs_getrc(ih,iswrr, te,de,ti,di, v0,vi)
  rscx = rca(1)
  call eirene_xs_getrc(ih,iswrr, tep,de,tip,di, v0,vi)
  scx1 = rca(1)
  call eirene_xs_getrc(ih,iswrr, tem,de,tim,di, v0,vi)
  scx2 = rca(1)
  rdscx= (scx1-scx2)/2./ti/dlinrat

  ! recombination + te-linearisation
  iswrr=6
  call eirene_xs_getrc(ih,iswrr, te,de,ti,di, v0,vi)
  rsvr = rcp(1)
  call eirene_xs_getrc(ih,iswrr, tep,de,tip,di, v0,vi)
  svr1 = rcp(1)
  call eirene_xs_getrc(ih,iswrr, tem,de,tim,di, v0,vi)
  svr2 = rcp(1)
  rdsvr= (svr1-svr2)/2./te/dlinrat

  
  svi=rsvi
  dsvi=rdsvi
  scx=rscx
  dscx=rdscx
  smi=rsmi
  dsmi=rdsmi
  smd=rsmd
  dsmd=rdsmd
  svr=rsvr
  dsvr=rdsvr

  ierr=0

end subroutine GetEireneXS
