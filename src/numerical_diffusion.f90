subroutine numerical_diffusion(sf_mn,dt1,itype)

 use params
 
 implicit none
 
 complex, dimension(mmax), intent(inout) :: sf_mn
 real, intent(in)                        :: dt1
 integer, intent(in)                     :: itype

! itype = 1 for vorticity and divergence and 0 for geopotential 
 
 integer :: i1, i2, ms, js, j_index2

! Apply horizontal diffusion to spectral coefficients - Nabla^4 

 do i1 = 0,mm
   ms = abs(i1)
   do i2 = ms,mm
     js = j_index2(mm,ms,i2)        
     sf_mn(js) = sf_mn(js)/(1.0 + 2.0*dt1*kdiff*((i2*(i2+1.0))**2 -4.0*itype)/a**4) 
   enddo  
 enddo
      
 return

end subroutine numerical_diffusion
