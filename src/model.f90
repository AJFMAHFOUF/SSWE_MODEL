subroutine model(xin,xout,dt1,npdt_max,ldfi,loutput)

 use params
 use model_vars
 use spectral_vars
 
 implicit none
 
 type (prog_var), intent(in)  :: xin
 type (prog_var), intent(out) :: xout
 real, intent(in)             :: dt1
 integer, intent(in)          :: npdt_max
 logical, intent(in)          :: ldfi, loutput
 
 complex, dimension(mmax) :: tend_vor_mn, tend_div_mn, tend_phi_mn, tend_qv_mn, filter, phi_diff
 integer :: nstep, i1, i2, ms, js, j_index2, nstep1
 real    :: t1, t2, zlap, zwindow, zsum, zsum2  
 
! Put initial conditions in spectral arrays
 
 vor_mn(:,1) = xin%vormn
 vor_mn(:,2) = xin%vormn
 vor_mn(:,3) = xin%vormn
 div_mn(:,1) = xin%divmn
 div_mn(:,2) = xin%divmn
 div_mn(:,3) = xin%divmn
 phi_mn(:,1) = xin%phimn
 phi_mn(:,2) = xin%phimn
 phi_mn(:,3) = xin%phimn
 
 zsum = 0.0

! Start temporal loop

 do nstep = 0,npdt_max-1
 
   call convert_vordiv2uv
 
   call compute_kinetic_energy
 
   call compute_vorticity_tendency(tend_vor_mn)
   call compute_divergence_tendency(nstep,dt1,tend_div_mn)
   call compute_geopotential_tendency(nstep,dt1,tend_phi_mn)
   !call compute_tracer_tendency(tend_qv_mn)
   tend_qv_mn(:) = (0.,0.)
   !tend_vor_mn(:) = (0.,0.)
   !tend_div_mn(:) = (0.,0.)
   !tend_phi_mn(:) = (0.,0.)
   
   nstep1 = nstep + 1
   
   zwindow = sin(nstep1*pi/(npdt_max+1.0))/(nstep1*pi/(npdt_max+1.0))*sin(nstep1*pi/npdt_max)/(nstep1*pi)
   zsum = zsum + zwindow
   
   if (nstep > 0) then
   
! Semi-implicit or explicit time scheme
   
     if (lsemimp) then 
       do i1 = 0,mm
         ms = abs(i1)       
         do i2 = ms,mm 
           js = j_index2(mm,ms,i2)
           zlap = i2*(i2 + 1.0)/(a*a)
           phi_mn(js,3) = (tend_phi_mn(js) - phi_bar*dt1*tend_div_mn(js))/(1.0 + phi_bar*zlap*dt1*dt1)  
           div_mn(js,3) = tend_div_mn(js) + dt1*zlap*phi_mn(js,3)                              
         enddo
       enddo
     else
       phi_mn(:,3) = phi_mn(:,1) + 2.0*dt1*tend_phi_mn(:) 
       div_mn(:,3) = div_mn(:,1) + 2.0*dt1*tend_div_mn(:) 
     endif
     
     vor_mn(:,3) = vor_mn(:,1) + 2.0*dt1*tend_vor_mn(:)
     qv_mn(:,3) = qv_mn(:,1) + 2.0*dt1*tend_qv_mn(:)
     
! Apply horizontal diffusion in spectral space

     if (.not.ldfi .and. dt1 > 0.0) then 
       call numerical_diffusion(vor_mn(:,3),dt1,1)
       call numerical_diffusion(div_mn(:,3),dt1,1)
      !call numerical_diffusion(qv_mn(:,3),dt1,0)

! Include orography for geopotential filtering

       phi_diff(:) = phi_mn(:,3) + phis_mn(:)
       call numerical_diffusion(phi_mn(:,3),dt1,0)
       phi_mn(:,3) = phi_diff(:) - phis_mn(:)     
     endif 
 
! Apply Robert Asselin Williams filter to remove 2*dt noise
     
     filter(:) = vor_mn(:,1) - 2.0*vor_mn(:,2) + vor_mn(:,3)
     vor_mn(:,2) = vor_mn(:,2) + nu*wk*filter(:)
     vor_mn(:,3) = vor_mn(:,3) - nu*(1.0-wk)*filter(:)
     
     filter(:) = div_mn(:,1) - 2.0*div_mn(:,2) + div_mn(:,3)
     div_mn(:,2) = div_mn(:,2) + nu*wk*filter(:)
     div_mn(:,3) = div_mn(:,3) - nu*(1.0-wk)*filter(:)
     
     filter(:) = phi_mn(:,1) - 2.0*phi_mn(:,2) + phi_mn(:,3)
     phi_mn(:,2) = phi_mn(:,2) + nu*wk*filter(:)
     phi_mn(:,3) = phi_mn(:,3) - nu*(1.0-wk)*filter(:)
     
     filter(:) = qv_mn(:,1) - 2.0*qv_mn(:,2) + qv_mn(:,3)
     qv_mn(:,2) = qv_mn(:,2) + nu*wk*filter(:)
     qv_mn(:,3) = qv_mn(:,3) - nu*(1.0-wk)*filter(:)
  
  
! Swap time steps    
   
     vor_mn(:,1) = vor_mn(:,2)
     vor_mn(:,2) = vor_mn(:,3)
     div_mn(:,1) = div_mn(:,2)
     div_mn(:,2) = div_mn(:,3)
     phi_mn(:,1) = phi_mn(:,2)
     phi_mn(:,2) = phi_mn(:,3)
     qv_mn(:,1)  = qv_mn(:,2)
     qv_mn(:,2)  = qv_mn(:,3)
     
     if (ldfi) then
       xout%vormn = xout%vormn + zwindow*vor_mn(:,2) 
       xout%divmn = xout%divmn + zwindow*div_mn(:,2) 
       xout%phimn = xout%phimn + zwindow*phi_mn(:,2)    
     endif
        
   else
   
     vor_mn(:,2) = vor_mn(:,1) + dt1*tend_vor_mn(:)
     div_mn(:,2) = div_mn(:,1) + dt1*tend_div_mn(:)
     phi_mn(:,2) = phi_mn(:,1) + dt1*tend_phi_mn(:)   
     qv_mn(:,2) = qv_mn(:,1) + dt1*tend_qv_mn(:)  
     
     if (ldfi) then
       xout%vormn = 0.5/float(npdt_max)*vor_mn(:,1) + zwindow*vor_mn(:,2) 
       xout%divmn = 0.5/float(npdt_max)*div_mn(:,1) + zwindow*div_mn(:,2) 
       xout%phimn = 0.5/float(npdt_max)*phi_mn(:,1) + zwindow*phi_mn(:,2)    
     endif
       
   endif     
  
! Write fields in physical space - spectral transforms in the subroutine
   
   if (mod(nstep,nfreq) == 0 .and. loutput) then
     call save_output(nstep)
     call compute_ke_spectrum(nstep)
   endif  
   
   if(.not.ldfi) then
     call legt_i(phi_m,phi_mn(:,2),0)
     call fft_i(phi,phi_m)  
     write (101,*) nstep,phi(64,48) 
   endif   
   
 enddo  
 zsum2 = 2.0*zsum + 1.0/float(npdt_max)
 
! Put final results in output arrays 

  if (.not.ldfi) then
    xout%vormn = vor_mn(:,2)
    xout%divmn = div_mn(:,2)
    xout%phimn = phi_mn(:,2)
  else
    xout%vormn = xout%vormn/zsum2
    xout%divmn = xout%divmn/zsum2
    xout%phimn = xout%phimn/zsum2  
  endif
 
 return

end subroutine model
