subroutine compute_divergence_tendency(nstep,tend_div_mn)

 use params
 use model_vars
 use spectral_vars
 
 implicit none
 
 integer, intent(in)                    :: nstep
 complex, dimension(mmax), intent(out)  :: tend_div_mn
 real    :: d_legpol, zlap
 integer :: i1, i2, j1, j_index2, js, ms
 
! One assumes that compute_vorticity_tendency has been called before
! => uvor_m and vvor_m fields have been already computed 
 
! Compute tendencies in spectral space

 tend_div_mn(:) = (0.,0.)
 do i1 = 0,mm
   ms = abs(i1)
   do i2 = ms,mm 
     js = j_index2(mm,ms,i2)
     do j1 = 1,nlat
       tend_div_mn(js) = tend_div_mn(js) + (j*i1*vvor_m(j1,i1)*pl_legendr(js,j1) + &
                         uvor_m(j1,i1)*d_legpol(j1,i1,i2))*w(j1)/(a*a*(1.0 - x(j1)*x(j1)))
     enddo
     tend_div_mn(js) = 0.5*tend_div_mn(js)  + i2*(i2+1.0)/a**2*(ke_mn(js) + phi_mn(js,2))
   enddo
 enddo                           
 
! Spectral quantities for semi-implicit scheme
 
 if (nstep > 0 .and. lsemimp) then
   do i1 = 0,mm
     ms = abs(i1)
     do i2 = ms,mm 
       js = j_index2(mm,ms,i2)
       zlap = i2*(i2 + 1.0)/(a*a)
       tend_div_mn(js) = tend_div_mn(js) - zlap*(ke_mn(js) + phi_mn(js,2))
       tend_div_mn(js) = div_mn(js,1) + 2.0*dt*(tend_div_mn(js) + zlap*(ke_mn(js) + 0.5*phi_mn(js,1)))
     enddo
   enddo  
 endif
 
 return 

end subroutine compute_divergence_tendency
