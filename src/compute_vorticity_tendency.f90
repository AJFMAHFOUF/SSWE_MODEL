subroutine compute_vorticity_tendency(tend_vor_mn)

 use params
 use model_vars
 use spectral_vars 
 
 implicit none
 
 real, dimension(nlon,nlat)             :: uvor, vvor
 complex, dimension(mmax), intent(out)  :: tend_vor_mn
 real    :: d_legpol
 integer :: i1, i2, j1, j_index2, js, ms
 
! Physical space  (absolute vorticity)
 
 do j1 = 1,nlat
   do i1  = 1,nlon
     uvor(i1,j1) = u(i1,j1)*(vor(i1,j1) + f(j1))
     vvor(i1,j1) = v(i1,j1)*(vor(i1,j1) + f(j1))
   enddo
 enddo    
 
! Fourier space for non linear terms

 call fft_d(uvor,uvor_m) 
 call fft_d(vvor,vvor_m) 

! Compute tendencies in spectral space

 tend_vor_mn(:) = (0.,0.)
 do i1 = 0,mm
   ms = abs(i1)
   do i2 = ms,mm 
     js = j_index2(mm,ms,i2)
     do j1 = 1,nlat    
       tend_vor_mn(js) = tend_vor_mn(js) - (j*i1*uvor_m(j1,i1)*pl_legendr(js,j1) - &
                         vvor_m(j1,i1)*d_legpol(j1,i1,i2))*w(j1)/(a*a*(1.0 - x(j1)*x(j1)))                 
     enddo
   enddo
 enddo                           
 tend_vor_mn(:) = 0.5*tend_vor_mn(:) 
 
 return 

end subroutine compute_vorticity_tendency
