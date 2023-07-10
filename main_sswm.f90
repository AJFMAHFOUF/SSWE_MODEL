program main_sswm
!=================================================================================
!
! Global spectral shallow water model in Fortran 90
! with vorticity, divergence and geopotential as prognostic variables
! Explicit leapfrog scheme - Triangular truncation (defined in module params)
!
! Equations taken from Coiffier (2011) - pages 128-129
!
! Jean-François Mahfouf (22/06/2023)
!
! Use of FFT99 defined by Temperton (ECMWF) 
!
!==================================================================================
 use params
 use model_vars
 use spectral_vars
 
 implicit none
 
 complex, dimension(mmax) :: tend_vor_mn, tend_div_mn, tend_phi_mn, filter
 integer :: nstep, i1, i2, ms, js, j_index2
 real :: t0, t1, t2, t3, zlap

! Required initialisations and read fields in physical space

 call cpu_time(time=t0)
 call init

 if (lreaduv) call convert_uv2vordiv
 
! Time integration  
 
 call cpu_time(time=t1)
 do nstep = 0,12*4
 
   call convert_vordiv2uv(nstep)
 
   call compute_kinetic_energy
 
   call compute_vorticity_tendency(tend_vor_mn)
   call compute_divergence_tendency(nstep,tend_div_mn)
   call compute_geopotential_tendency(nstep,tend_phi_mn)
   
   if (nstep > 0) then
   
! Semi-implicit or explicit time scheme
   
     if (lsemimp) then 
       do i1 = 0,mm
         ms = abs(i1)       
         do i2 = ms,mm 
           js = j_index2(mm,ms,i2)
           zlap = i2*(i2 + 1.0)/(a*a)
           phi_mn(js,3) = (tend_phi_mn(js) - phi_bar*dt*tend_div_mn(js))/(1.0 + phi_bar*zlap*dt*dt)  
           div_mn(js,3) = tend_div_mn(js) + dt*zlap*phi_mn(js,3)                              
         enddo
       enddo
     else
       phi_mn(:,3) = phi_mn(:,1) + 2.0*dt*tend_phi_mn(:) 
       div_mn(:,3) = div_mn(:,1) + 2.0*dt*tend_div_mn(:) 
     endif
     
     vor_mn(:,3) = vor_mn(:,1) + 2.0*dt*tend_vor_mn(:)
 
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
  
! Swap time steps    
   
     vor_mn(:,1) = vor_mn(:,2)
     vor_mn(:,2) = vor_mn(:,3)
     div_mn(:,1) = div_mn(:,2)
     div_mn(:,2) = div_mn(:,3)
     phi_mn(:,1) = phi_mn(:,2)
     phi_mn(:,2) = phi_mn(:,3)
        
   else
     vor_mn(:,2) = vor_mn(:,1) + dt*tend_vor_mn(:)
     div_mn(:,2) = div_mn(:,1) + dt*tend_div_mn(:)
     phi_mn(:,2) = phi_mn(:,1) + dt*tend_phi_mn(:)   
   endif     
  
   !if (nstep == 1) call compute_ke_spectrum(nstep)
   
   if (mod(nstep,4) == 0) call save_output(nstep)
   
 enddo  
 
! Write fields in physical space - spectral transforms in the subroutine

 nstep = nstep - 1
 call cpu_time(time=t2)
 call compute_ke_spectrum(nstep)
 !call save_output(nstep) 
 call cpu_time(time=t3)
 
 print *,'temps execution modele =',t2-t1
 print *,'temps total avec initialisation et ecriture resultats =',t3-t0
 

end program main_sswm
