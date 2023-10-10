subroutine compute_balance_equation

 use params
 use model_vars
 use spectral_vars
 
 implicit none
 
 real, dimension (nlon,nlat)     :: linbal
 complex, dimension(nlat,-mm:mm) :: linbal_m
 complex, dimension(mmax)        :: linbal_mn
 complex                         :: vor1, vor2
 real                            :: beta1, beta2
 
 integer :: i1,i2,js,ms,j1,j_index2,jsp,jsm
 
! Physical space - compute div (f x nabla psi) - Eq.(7.5.3) from Daley (1991)
 
 do j1 = 1,nlat
   linbal(:,j1) = vor(:,j1)*f(j1) - 2.0*u(:,j1)*omega/(a*a)
 enddo   
 
! Spectral space

 call fft_d(linbal,linbal_m) 
 call legt_d(linbal_m,linbal_mn)
 
! Solve linear balance equation to obtain the geopotential - Derber and Bouttier (1999)

 phi_mn(1) = (50000.,0.)  ! mean geopotential
 do i1 = 0,mm
   ms=abs(i1)
   do i2 = ms,mm
     js = j_index2(mm,ms,i2)
     !if (ms+1 <= mm) then
     !  jsp = j_index2(mm,ms+1,i2)
     !  vor1 = vor_mn(jsp,2)
     !else
     !  vor1 = (0.,0.)
     !endif
     !if (ms-1 >= 0) then    
     !  jsm = j_index2(mm,ms-1,i2)
     !  vor2 = vor_mn(jsm,2)
     !else
     !  vor2 = (0.,0.)
     !endif
     if (i2 > 0) then
       !beta1 = -2.0*omega*a*a/(i2*(i2+1.0))*(1.0 - 1.0/(i2+1.0))*sqrt((i2+1.0+ms)*(i2+1.0-ms)/((2.*i2+3.)*(2.*i2+1.)))  
       !beta2 = -2.0*omega*a*a/(i2*(i2+1.0))*(1.0 + 1.0/i2)*sqrt((i2+ms)*(i2-ms)/((2.*i2+1.)*(2.*i2-1.)))
       !phi_mn(js) = beta1*vor1 + beta2*vor2
       phi_mn(js) = (-a*a/(i2*(i2+1.0)))*linbal_mn(js)     
     endif  
   enddo
 enddo  
 
 return 

end subroutine compute_balance_equation
