subroutine convert_uv2vor

 use params
 use model_vars
 use spectral_vars
 
 implicit none
 
 integer :: i1, i2, j1, ms, js, j_index2
 real    :: zvor, ztheta0, zthetaw, zm, za, zlon, zlat
 real    :: d_legpol

! Spectral coefficients for u and v components (transformed fields)

 call fft_d(u,u_m)
 call fft_d(v,v_m)

! Spectral coefficients for vorticity 

 vor_mn(:,2) = (0.,0.)
 do i1 = 0,mm
   ms = abs(i1)
   do i2 = ms,mm
     js = j_index2(mm,ms,i2)
     do j1 = 1,nlat        
       vor_mn(js,2) = vor_mn(js,2) + w(j1)*(j*i1*v_m(j1,i1)*pl_legendr(js,j1) + &
                      u_m(j1,i1)*d_legpol(j1,i1,i2))/(a*a*(1.0-x(j1)*x(j1)))
     enddo  
   enddo
 enddo
 vor_mn(:,2) = 0.5*vor_mn(:,2)
 
!  For initial step
  
 if (.not.l_real_ic) then 
 
   call legt_i(vor_m,vor_mn(:,2),0)
   call fft_i(vor,vor_m)
 
! Analytical initial conditions proposed by Held and Phillips (1987) JAS
 
   ztheta0 = 0.25*pi
   zthetaw = 15.0*pi/180. 
   zm = 4.0
   za = 8.0E-5
   
   do j1=1,nlon
     do i1=1,nlat
       zlon = (j1 - 1.0)/float(nlon)*2.0*pi
       zlat = asin(x(i1))
       zvor = 0.5*za*sqrt(1.0 - x(i1)*x(i1))*exp(-((zlat -ztheta0)/zthetaw)**2)*cos(zm*zlon)
       vor(j1,i1) = vor(j1,i1) + zvor
     enddo
   enddo  

   call fft_d(vor,vor_m)    
   call legt_d(vor_m,vor_mn(:,2))
   
 endif  

 vor_mn(:,1) = vor_mn(:,2)
 vor_mn(:,3) = vor_mn(:,2)
      
 return

end subroutine convert_uv2vor
