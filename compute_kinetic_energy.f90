subroutine compute_kinetic_energy

 use params
 use model_vars
 use spectral_vars
 implicit none
 integer :: i1, j1, i2, js , ms, j_index2
 
! Physical space 
 
 do j1 = 1,nlat
   do i1  = 1,nlon
     ke(i1,j1) = 0.5*(u(i1,j1)*u(i1,j1) + v(i1,j1)*v(i1,j1))/(a*a*(1.0 - x(j1)*x(j1)))
   enddo
 enddo    
 
! Spectral space

 call fft_d(ke,ke_m) 
 call legt_d(ke_m,ke_mn)
 
 return 

end subroutine compute_kinetic_energy
