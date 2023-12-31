subroutine compute_geopotential_tendency(nstep,dt1,tend_phi_mn)

 use params
 use model_vars
 use spectral_vars
 
 implicit none
 
 integer, intent(in)                    :: nstep
 real, intent(in)                       :: dt1
 complex, dimension(mmax), intent(out)  :: tend_phi_mn
 real    :: d_legpol
 integer :: i1, i2, j1, j_index2, ms, js
 
! Physical space 
 
 do j1 = 1,nlat
   do i1  = 1,nlon
     uvar(i1,j1) = u(i1,j1)*phi(i1,j1)
     vvar(i1,j1) = v(i1,j1)*phi(i1,j1)
   enddo
 enddo  
 
! Fourier space for non linear terms

 call fft_d(uvar,uvar_m) 
 call fft_d(vvar,vvar_m) 

! Compute tendencies in spectral space

 tend_phi_mn(:) = (0.,0.)
 do i1 = 0,mm
   ms = abs(i1)
   do i2 = ms,mm 
     js = j_index2(mm,ms,i2)
     do j1 = 1,nlat       
       tend_phi_mn(js) = tend_phi_mn(js) - (j*i1*uvar_m(j1,i1)*pl_legendr(js,j1) - &
                         vvar_m(j1,i1)*d_legpol(j1,i1,i2))*w(j1)/(a*a*(1.0 - x(j1)*x(j1)))
     enddo
   enddo
 enddo                           
 tend_phi_mn(:) = 0.5*tend_phi_mn(:) 
 
! Spectral quantities for semi-implicit scheme
 
 if (nstep > 0 .and. lsemimp) then 
   do i1 = 0,mm
     ms = abs(i1)
     do i2 = ms,mm 
       js = j_index2(mm,ms,i2)        
       tend_phi_mn(js) = phi_mn(js,1) + 2.0*dt1*(tend_phi_mn(js) + phi_bar*(div_mn(js,2) - 0.5*div_mn(js,1))) 
     enddo
   enddo  
 endif
 
 return 

end subroutine compute_geopotential_tendency
