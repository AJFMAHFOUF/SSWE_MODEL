subroutine compute_tracer_tendency(tend_qv_mn)

 use params
 use model_vars
 use spectral_vars

 implicit none

 complex, dimension(mmax), intent(out)  :: tend_qv_mn
 real    :: d_legpol
 integer :: i1, i2, j1, j_index2, js, ms

! Physical space

 do j1 = 1,nlat
   do i1  = 1,nlon
     uvar(i1,j1) = u(i1,j1)*qv(i1,j1) !+ f(j1)) ! to mimic vorticity
     vvar(i1,j1) = v(i1,j1)*qv(i1,j1) !+ f(j1)) ! to mimic vorticity
   enddo
 enddo

! Fourier space for non linear terms

 call fft_d(uvar,uvar_m)
 call fft_d(vvar,vvar_m)

! Compute tendencies in spectral space

 tend_qv_mn(:) = (0.,0.)
 do i1 = 0,mm
   ms = abs(i1)
   do i2 = ms,mm
     js = j_index2(mm,ms,i2)
     do j1 = 1,nlat
       tend_qv_mn(js) = tend_qv_mn(js) - (j*i1*uvar_m(j1,i1)*pl_legendr(js,j1) - &
                         vvar_m(j1,i1)*d_legpol(j1,i1,i2))*w(j1)/(a*a*(1.0 - x(j1)*x(j1)))
     enddo
   enddo
 enddo
 tend_qv_mn(:) = 0.5*tend_qv_mn(:)

 return

end subroutine compute_tracer_tendency
