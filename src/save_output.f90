subroutine save_output(nstep)

 use params
 use model_vars
 use spectral_vars
 
 implicit none
 
 integer, intent(in) :: nstep
 character(len=1)    :: ichst1
 character(len=2)    :: ichst2
 character(len=3)    :: ichst3
 character(len=4)    :: ichst
 character(len=3)    :: tt
 character(len=2)    :: tt1
 integer :: i1, j1, i2, ms, js, j_index2, ihour
 real :: zlon, zlat, zu, zv, zphi, zfac  
  
 ihour = nstep*dt/3600.0
 
 if (ihour < 10) then 
   write(ichst1,'(i1)') ihour
   ichst = '000'//ichst1
 elseif (ihour < 100) then
   write(ichst2,'(i2)') ihour
   ichst = '00'//ichst2
 elseif (ihour < 1000) then
   write(ichst3,'(i3)') ihour
   ichst = '0'//ichst3   
 else
   write(ichst,'(i4)') ihour
 endif      
    
 if (mm < 99) then    
   write(tt1,'(i2)') mm 
   tt = '0'//tt1
 else
   write(tt,'(i3)') mm
 endif    
  
! Back to physical space - vorticity

 call legt_i(vor_m,vor_mn(:,2),0)
 call fft_i(vor,vor_m) 
  
! Back to physical space - divergence  
  
 call legt_i(div_m,div_mn(:,2),0)
 call fft_i(div,div_m)  
  
! Back to physical space - geopotential
  
 call legt_i(phi_m,phi_mn(:,2),0)
 call fft_i(phi,phi_m)   
 
! Back to physical space - passive tracer
  
 call legt_i(qv_m,qv_mn(:,2),0)
 call fft_i(qv,qv_m)     
  
! Back to physical space - wind components

 call legt_i(u_m,u_mn,1)
 call fft_i(u,u_m) 
  
 call legt_i(v_m,v_mn,1)
 call fft_i(v,v_m) 
 
! Compute streamfunction and velocity potential  

 psi_mn(:) = (0.,0.)
 khi_mn(:) = (0.,0.) 
 do i1 = 0,mm
   ms=abs(i1)
   do i2 = ms,mm
     js = j_index2(mm,ms,i2)
     if (i2 > 0) then
       psi_mn(js) = (-a*a/(i2*(i2+1.0)))*vor_mn(js,2)
       khi_mn(js) = (-a*a/(i2*(i2+1.0)))*div_mn(js,2)       
     endif  
   enddo
 enddo 
 
! Back to physical space - khi and psi

 call legt_i(psi_m,psi_mn,0)
 call fft_i(psi,psi_m) 
  
 call legt_i(khi_m,khi_mn,0)
 call fft_i(khi,khi_m)  
 
! Write results in ASCII file for plotting purposes
  
 open (unit=20,file='../data_out/SSWE_model_T'//tt//'_step_'//ichst//'_expid_'//expid//'.dat',status='unknown')
  
 do j1=1,nlat 
   zfac = 1.0/(a*sqrt(1.0 - x(j1)*x(j1)))
   do i1=1,nlon
     zlon = float(i1-1)/float(nlon)*360.0
     zlat = asin(x(j1))*180.0/pi
     zu = u(i1,j1)*zfac   ! real U wind
     zv = v(i1,j1)*zfac   ! real V wind
     zphi = phi(i1,j1) + phis(i1,j1)
     write(20,*) zlon,zlat,vor(i1,j1),div(i1,j1),zphi,zu,zv,psi(i1,j1),khi(i1,j1),qv(i1,j1)
   enddo
   zu = u(1,j1)*zfac
   zv = v(1,j1)*zfac
   zphi = phi(1,j1) 
   write(20,*) 360.0,zlat,vor(1,j1),div(1,j1),zphi,zu,zv,psi(1,j1),khi(1,j1),qv(1,j1)
 enddo  
  
 close (unit=20)  
  
 return  
 
end subroutine save_output
