subroutine convert_uv2vordiv

 use params
 use model_vars
 use spectral_vars
 
 implicit none
 
 integer :: i1, i2, j1, ms, js, j_index2
 real    :: d_legpol

! Spectral coefficients for u and v components (transformed fields)

 call fft_d(u,u_m)
 call fft_d(v,v_m)

! Spectral coefficients for vorticity and divergence

 vor_mn(:,2) = (0.,0.)
 div_mn(:,2) = (0.,0.)
 do i1 = 0,mm
   ms = abs(i1)
   do i2 = ms,mm
     js = j_index2(mm,ms,i2)
     do j1 = 1,nlat        
       vor_mn(js,2) = vor_mn(js,2) + w(j1)*(j*i1*v_m(j1,i1)*pl_legendr(js,j1) + &
                      u_m(j1,i1)*d_legpol(j1,i1,i2))/(a*a*(1.0-x(j1)*x(j1)))
       div_mn(js,2) = div_mn(js,2) + w(j1)*(j*i1*u_m(j1,i1)*pl_legendr(js,j1) - &
                      v_m(j1,i1)*d_legpol(j1,i1,i2))/(a*a*(1.0-x(j1)*x(j1)))
     enddo  
   enddo
 enddo
 vor_mn(:,2) = 0.5*vor_mn(:,2)
 div_mn(:,2) = 0.5*div_mn(:,2)  
 
!  For initial step

 vor_mn(:,1) = vor_mn(:,2)
 vor_mn(:,3) = vor_mn(:,2)
 div_mn(:,1) = div_mn(:,2)
 div_mn(:,3) = div_mn(:,2)   
      
 return

end subroutine convert_uv2vordiv
