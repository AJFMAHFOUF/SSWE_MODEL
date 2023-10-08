subroutine model(xin,xout,dt1,npdt_max,loutput)

 use params
 use model_vars
 use spectral_vars
 
 implicit none
 
 type (prog_var), intent(in)  :: xin
 type (prog_var), intent(out) :: xout
 real, intent(in)             :: dt1
 integer, intent(in)          :: npdt_max
 logical, intent(in)          :: loutput
 
 complex, dimension(mmax) :: tend_vor_mn, filter
 integer :: nstep
 
! Put initial conditions in spectral arrays
 
 vor_mn(:,1) = xin%vormn
 vor_mn(:,2) = xin%vormn
 vor_mn(:,3) = xin%vormn

! Start temporal loop

 do nstep = 0,npdt_max-1
 
   call convert_vor2uv
 
   call compute_kinetic_energy
 
   call compute_vorticity_tendency(tend_vor_mn)
   
   if (nstep > 0) then
     
     vor_mn(:,3) = vor_mn(:,1) + 2.0*dt1*tend_vor_mn(:)
     
! Apply horizontal diffusion in spectral space

     call numerical_diffusion(vor_mn(:,3),dt1,1)
 
! Apply Robert Asselin Williams filter to remove 2*dt noise
     
     filter(:) = vor_mn(:,1) - 2.0*vor_mn(:,2) + vor_mn(:,3)
     vor_mn(:,2) = vor_mn(:,2) + nu*wk*filter(:)
     vor_mn(:,3) = vor_mn(:,3) - nu*(1.0-wk)*filter(:)
      
! Swap time steps    
   
     vor_mn(:,1) = vor_mn(:,2)
     vor_mn(:,2) = vor_mn(:,3)
        
   else
   
     vor_mn(:,2) = vor_mn(:,1) + dt1*tend_vor_mn(:)
       
   endif     
  
! Write fields in physical space - spectral transforms in the subroutine
   
   if (mod(nstep,nfreq) == 0 .and. loutput) then
     call save_output(nstep)
     call compute_ke_spectrum(nstep)
   endif  
   
 enddo  
 
! Put final results in output arrays 

   xout%vormn = vor_mn(:,2)
 
 return

end subroutine model
