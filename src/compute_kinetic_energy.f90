subroutine compute_kinetic_energy

 use params
 use model_vars
 use spectral_vars
 
 implicit none
 
 integer :: j1
 
! Physical space 
 
 do j1 = 1,nlat
   ke(:,j1) = 0.5*(u(:,j1)*u(:,j1) + v(:,j1)*v(:,j1))/(a*a*(1.0 - x(j1)*x(j1)))
 enddo   
 
! Spectral space

 call fft_d(ke,ke_m) 
 call legt_d(ke_m,ke_mn)
 
 return 

end subroutine compute_kinetic_energy
